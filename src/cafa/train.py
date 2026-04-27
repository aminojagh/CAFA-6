from __future__ import annotations

import gzip
import tarfile
from collections.abc import Iterator
from pathlib import Path

from Bio import SwissProt

from .config import ProjectConfig
from .sources import download_source
from .types import ProteinTaxonRecord, SourceSnapshot

_SWISSPROT_FLATFILE_MEMBER = "uniprot_sprot.dat.gz"


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
    """

    allowed_taxon_ids = set(config.train_taxon_ids)
    if not allowed_taxon_ids:
        return ()

    archive_path = download_source(swissprot_snapshot)
    protein_to_taxon: dict[str, int] = {}
    for protein_id, taxon_id in _iter_primary_accession_and_taxon_pairs(archive_path):
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


def _iter_primary_accession_and_taxon_pairs(
    archive_path: str | Path,
    member_name: str = _SWISSPROT_FLATFILE_MEMBER,
) -> Iterator[tuple[str, int]]:
    """Yield `(primary_accession, primary_taxon_id)` from a Swiss-Prot archive."""

    with tarfile.open(archive_path, "r:gz") as archive:
        try:
            member = archive.getmember(member_name)
        except KeyError as exc:
            raise FileNotFoundError(
                f"Swiss-Prot archive does not contain {member_name}."
            ) from exc

        compressed_handle = archive.extractfile(member)
        if compressed_handle is None:
            raise FileNotFoundError(f"Unable to read {member_name} from Swiss-Prot archive.")

        with gzip.open(compressed_handle, mode="rt", encoding="utf-8") as flatfile_handle:
            for record in SwissProt.parse(flatfile_handle):
                yield _primary_accession(record), _primary_taxon_id(record)


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
