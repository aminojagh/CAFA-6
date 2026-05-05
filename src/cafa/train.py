from __future__ import annotations

import gzip
from collections.abc import Iterable
from collections.abc import Iterator
from pathlib import Path

from Bio import SwissProt
from tqdm.auto import tqdm

from .config import ProjectConfig
from .ontology import GeneOntology, canonicalize_go_id
from .sources import download_source, extract_tar_gz_member
from .types import ProteinTaxonRecord, ProteinTermRecord, SourceSnapshot, Subontology

_SWISSPROT_FLATFILE_MEMBER = "uniprot_sprot.dat.gz"
_NAMESPACE_TO_ASPECT: dict[str, str] = {
    "molecular_function": "F",
    "biological_process": "P",
    "cellular_component": "C",
}
_NAMESPACE_TO_SUBONTOLOGY: dict[str, Subontology] = {
    "molecular_function": "MF",
    "biological_process": "BP",
    "cellular_component": "CC",
}


def extract_train_taxonomy_records(
    config: ProjectConfig,
    swissprot_snapshot: SourceSnapshot,
) -> tuple[ProteinTaxonRecord, ...]:
    """Extract the recreated train taxonomy rows from Swiss-Prot.

    Parameters
    ----------
    config:
        Normalized project configuration. Only `train_taxon_ids` is used in
        this first train-extraction slice.
    swissprot_snapshot:
        Resolved Swiss-Prot release archive that contains
        `uniprot_sprot.dat.gz`.

    Returns
    -------
    tuple[ProteinTaxonRecord, ...]
        Deterministically ordered protein-to-taxon mappings for Swiss-Prot
        entries whose primary taxonomy ID is in `config.train_taxon_ids`.

    Notes
    -----
    - The primary UniProt accession is `record.accessions[0]`.
    - The primary taxon ID is `int(record.taxonomy_id[0])`.
    - This function intentionally reads the annotation-bearing Swiss-Prot
      flatfile, not FASTA, because the same source remains authoritative for
      later train terms and train sequences.
    - The raw `uniprot_sprot.dat.gz` member is extracted once under the cache
      tree and then reused on later runs. This speeds up notebook reruns
      without caching any logic-bearing transformation output.
    """

    allowed_taxon_ids = set(config.train_taxon_ids)
    if not allowed_taxon_ids:
        return ()

    archive_path = download_source(swissprot_snapshot)
    flatfile_gz_path = extract_tar_gz_member(archive_path, _SWISSPROT_FLATFILE_MEMBER)
    protein_to_taxon: dict[str, int] = {}
    for record in _iter_swissprot_records(flatfile_gz_path):
        protein_id = _primary_accession(record)
        taxon_id = _primary_taxon_id(record)
        if taxon_id not in allowed_taxon_ids:
            continue
        existing_taxon_id = protein_to_taxon.get(protein_id)
        if existing_taxon_id is not None and existing_taxon_id != taxon_id:
            raise ValueError(
                f"Conflicting taxon assignments for {protein_id}: "
                f"{existing_taxon_id} vs {taxon_id}."
            )
        protein_to_taxon[protein_id] = taxon_id

    return tuple(
        ProteinTaxonRecord(protein_id=protein_id, taxon_id=taxon_id)
        for protein_id, taxon_id in sorted(protein_to_taxon.items())
    )


def extract_train_term_records(
    config: ProjectConfig,
    swissprot_snapshot: SourceSnapshot,
    ontology: GeneOntology,
    taxonomy_rows: tuple[ProteinTaxonRecord, ...],
) -> tuple[ProteinTermRecord, ...]:
    """Extract recreated train term rows from the pinned Swiss-Prot flatfile.

    Parameters
    ----------
    config:
        Normalized project configuration. This slice uses `train_taxon_ids`,
        `subontologies`, and `evidence_codes`.
    swissprot_snapshot:
        Resolved Swiss-Prot release archive that contains
        `uniprot_sprot.dat.gz`.
    ontology:
        Downloaded and validated GO ontology used as the source of truth for
        canonical GO IDs and aspect/subontology derivation.
    taxonomy_rows:
        Recreated train-taxonomy rows. These define the allowed protein set for
        train-term extraction.

    Returns
    -------
    tuple[ProteinTermRecord, ...]
        Deterministically ordered direct protein-to-GO-term mappings filtered by
        the configured taxonomy gate, evidence-code whitelist, and selected
        subontologies.

    Notes
    -----
    - GO terms are canonicalized through the ontology before writing.
    - Aspect is derived from the ontology namespace, not from Swiss-Prot GO
      text, because ontology-derived fields have higher confidence.
    - Unknown, obsolete, or out-of-scope GO IDs are ignored.
    """

    allowed_protein_ids = {
        row.protein_id
        for row in taxonomy_rows
        if row.taxon_id in set(config.train_taxon_ids)
    }
    if not allowed_protein_ids:
        return ()

    allowed_evidence_codes = {code.upper() for code in config.evidence_codes}
    allowed_subontologies = set(config.subontologies)

    archive_path = download_source(swissprot_snapshot)
    flatfile_gz_path = extract_tar_gz_member(archive_path, _SWISSPROT_FLATFILE_MEMBER)
    extracted_rows: set[ProteinTermRecord] = set()

    for record in _iter_swissprot_records(flatfile_gz_path):
        protein_id = _primary_accession(record)
        if protein_id not in allowed_protein_ids:
            continue

        for go_term_id, evidence_code in _iter_go_term_and_evidence_pairs(record):
            if evidence_code not in allowed_evidence_codes:
                continue
            canonical_term_id = canonicalize_go_id(ontology, go_term_id)
            term = ontology.terms.get(canonical_term_id)
            if term is None:
                continue
            subontology = _NAMESPACE_TO_SUBONTOLOGY.get(term.namespace)
            if subontology not in allowed_subontologies:
                continue
            aspect = _NAMESPACE_TO_ASPECT.get(term.namespace)
            if aspect is None:
                continue
            extracted_rows.add(
                ProteinTermRecord(
                    protein_id=protein_id,
                    term_id=canonical_term_id,
                    aspect=aspect,
                )
            )

    return tuple(
        sorted(
            extracted_rows,
            key=lambda row: (row.protein_id, row.term_id, row.aspect),
        )
    )


def _iter_swissprot_records(flatfile_gz_path: str | Path) -> Iterator[SwissProt.Record]:
    """Yield parsed Swiss-Prot records from an extracted flatfile member."""

    total_records = _count_swissprot_records(flatfile_gz_path)
    with gzip.open(flatfile_gz_path, mode="rt", encoding="utf-8") as flatfile_handle:
        for record in tqdm(
            SwissProt.parse(flatfile_handle),
            total=total_records,
            desc="Swiss-Prot records",
            unit="record",
        ):
            yield record


def _count_swissprot_records(flatfile_gz_path: str | Path) -> int:
    """Count Swiss-Prot flatfile records using `//` record terminators."""

    count = 0
    with gzip.open(flatfile_gz_path, mode="rt", encoding="utf-8") as flatfile_handle:
        for line in flatfile_handle:
            if line.rstrip("\n") == "//":
                count += 1
    return count


def _iter_go_term_and_evidence_pairs(
    record: SwissProt.Record,
) -> Iterator[tuple[str, str]]:
    """Yield `(go_term_id, evidence_code)` pairs from Swiss-Prot GO cross-references."""

    for reference in _go_cross_references(record):
        if len(reference) < 4:
            continue
        go_term_id = str(reference[1]).strip()
        evidence_code = _go_evidence_code(reference)
        if not go_term_id or not evidence_code:
            continue
        yield go_term_id, evidence_code


def _go_cross_references(record: SwissProt.Record) -> Iterator[tuple[str, ...]]:
    """Yield GO cross-references from one Swiss-Prot record."""

    for reference in record.cross_references:
        if not reference:
            continue
        if str(reference[0]).strip().upper() != "GO":
            continue
        yield tuple(str(value) for value in reference)


def _go_evidence_code(reference: Iterable[str]) -> str:
    """Extract the GO evidence-code prefix from one GO cross-reference tuple."""

    values = tuple(reference)
    if len(values) < 4:
        return ""
    return values[3].split(":", 1)[0].strip().upper()


def _primary_accession(record: SwissProt.Record) -> str:
    if not record.accessions:
        raise ValueError("Swiss-Prot record is missing a primary accession.")
    return record.accessions[0]


def _primary_taxon_id(record: SwissProt.Record) -> int:
    if not record.taxonomy_id:
        primary_accession = record.accessions[0] if record.accessions else "<missing accession>"
        raise ValueError(
            f"Swiss-Prot record {primary_accession} is missing a taxonomy identifier."
        )
    return int(record.taxonomy_id[0])
