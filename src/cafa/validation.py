from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from collections.abc import Iterable
from typing import AbstractSet

from .config import ProjectConfig
from .io import (
    read_fasta_records,
    read_ia_values,
    read_train_taxonomy_rows,
    read_train_term_rows,
)
from .ontology import GeneOntology, read_go_obo
from .types import (
    ProteinId,
    ProteinTaxonRecord,
    ProteinTermRecord,
    SequenceRecord,
    Subontology,
    TaxonId,
    ValidationReport,
)

_ASPECT_TO_SUBONTOLOGY: dict[str, Subontology] = {
    "F": "MF",
    "P": "BP",
    "C": "CC",
}


def validate_train_taxonomy(
    recreated_path: str | Path,
    reference_path: str | Path,
    config: ProjectConfig,
) -> ValidationReport:
    """Validate recreated train taxonomy rows against the filtered reference."""

    recreated_rows = filter_reference_train_taxonomy_rows(
        read_train_taxonomy_rows(recreated_path),
        allowed_taxon_ids=set(config.train_taxon_ids),
    )
    reference_rows = filter_reference_train_taxonomy_rows(
        read_train_taxonomy_rows(reference_path),
        allowed_taxon_ids=set(config.train_taxon_ids),
    )

    recreated_mapping = {row.protein_id: row.taxon_id for row in recreated_rows}
    reference_mapping = {row.protein_id: row.taxon_id for row in reference_rows}

    return _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message="Protein-to-taxon membership mismatch.",
        left_mapping=recreated_mapping,
        right_mapping=reference_mapping,
        formatter=lambda protein_id, taxon_id: f"{protein_id}\t{taxon_id}",
    )


def validate_go_obo(
    recreated_path: str | Path,
    reference_path: str | Path,
) -> ValidationReport:
    """Validate GO release metadata, canonical term set, and parent edges."""

    recreated_ontology = read_go_obo(recreated_path)
    reference_ontology = read_go_obo(reference_path)

    if recreated_ontology.release != reference_ontology.release:
        return ValidationReport(
            left_path=Path(recreated_path),
            right_path=Path(reference_path),
            passed=False,
            message="GO release metadata mismatch.",
            sample_left_only=((recreated_ontology.release or "None"),),
            sample_right_only=((reference_ontology.release or "None"),),
        )

    recreated_terms = {
        term_id: tuple(sorted(term.parent_ids))
        for term_id, term in recreated_ontology.terms.items()
    }
    reference_terms = {
        term_id: tuple(sorted(term.parent_ids))
        for term_id, term in reference_ontology.terms.items()
    }

    return _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message="GO canonical term set or parent-edge structure mismatch.",
        left_mapping=recreated_terms,
        right_mapping=reference_terms,
        formatter=lambda term_id, parent_ids: (
            f"{term_id}\t{','.join(parent_ids)}" if parent_ids else term_id
        ),
    )


def validate_train_terms(
    recreated_path: str | Path,
    reference_path: str | Path,
    config: ProjectConfig,
    ontology: GeneOntology,
    taxonomy_rows: tuple[ProteinTaxonRecord, ...],
) -> ValidationReport:
    """Validate recreated train terms against the filtered reference."""

    allowed_protein_ids = {
        row.protein_id
        for row in filter_reference_train_taxonomy_rows(
            taxonomy_rows,
            allowed_taxon_ids=set(config.train_taxon_ids),
        )
    }
    recreated_rows = filter_reference_train_term_rows(
        read_train_term_rows(recreated_path, ontology=ontology),
        allowed_protein_ids=allowed_protein_ids,
        allowed_subontologies=set(config.subontologies),
    )
    filtered_reference_rows = filter_reference_train_term_rows(
        read_train_term_rows(reference_path, ontology=ontology),
        allowed_protein_ids=allowed_protein_ids,
        allowed_subontologies=set(config.subontologies),
    )

    recreated_mapping = _group_terms_by_protein(recreated_rows)
    reference_mapping = _group_terms_by_protein(filtered_reference_rows)

    return _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message="Protein-to-direct-GO-term mapping mismatch.",
        left_mapping=recreated_mapping,
        right_mapping=reference_mapping,
        formatter=lambda protein_id, term_records: (
            f"{protein_id}\t"
            + ",".join(f"{record.term_id}:{record.aspect}" for record in term_records)
        ),
    )


def validate_sequence_mapping(
    recreated_path: str | Path,
    reference_path: str | Path,
    artifact_name: str,
) -> ValidationReport:
    """Validate exact protein-to-sequence mapping between two FASTA artifacts."""

    recreated_records = read_fasta_records(recreated_path)
    reference_records = read_fasta_records(reference_path)

    recreated_mapping = {record.protein_id: record.sequence for record in recreated_records}
    reference_mapping = {record.protein_id: record.sequence for record in reference_records}

    return _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message=f"{artifact_name} protein-to-sequence mapping mismatch.",
        left_mapping=recreated_mapping,
        right_mapping=reference_mapping,
        formatter=lambda protein_id, sequence: f"{protein_id}\t{sequence}",
    )


def validate_ia_values(
    recreated_path: str | Path,
    reference_path: str | Path,
    ontology: GeneOntology | None = None,
    relative_tolerance: float = 1e-9,
    absolute_tolerance: float = 1e-12,
) -> ValidationReport:
    """Validate IA term coverage and numeric agreement within tolerance."""

    recreated_values = read_ia_values(recreated_path, ontology=ontology)
    reference_values = read_ia_values(reference_path, ontology=ontology)

    left_only_keys = sorted(set(recreated_values) - set(reference_values))
    right_only_keys = sorted(set(reference_values) - set(recreated_values))
    shared_mismatch_keys = sorted(
        term_id
        for term_id in set(recreated_values) & set(reference_values)
        if not _is_close(
            recreated_values[term_id],
            reference_values[term_id],
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
    )

    sample_left_only = tuple(
        f"{term_id}\t{recreated_values[term_id]}"
        for term_id in left_only_keys[:5]
    ) + tuple(
        f"{term_id}\t{recreated_values[term_id]}"
        for term_id in shared_mismatch_keys[:5]
    )
    sample_right_only = tuple(
        f"{term_id}\t{reference_values[term_id]}"
        for term_id in right_only_keys[:5]
    ) + tuple(
        f"{term_id}\t{reference_values[term_id]}"
        for term_id in shared_mismatch_keys[:5]
    )

    passed = not left_only_keys and not right_only_keys and not shared_mismatch_keys
    return ValidationReport(
        left_path=Path(recreated_path),
        right_path=Path(reference_path),
        passed=passed,
        message="" if passed else "IA term coverage or numeric values mismatch.",
        sample_left_only=sample_left_only,
        sample_right_only=sample_right_only,
    )


def filter_reference_train_taxonomy_rows(
    rows: tuple[ProteinTaxonRecord, ...],
    allowed_taxon_ids: AbstractSet[TaxonId],
) -> tuple[ProteinTaxonRecord, ...]:
    """Filter reference train taxonomy rows by allowed taxa."""

    return tuple(
        row
        for row in rows
        if row.taxon_id in allowed_taxon_ids
    )


def filter_reference_train_term_rows(
    rows: tuple[ProteinTermRecord, ...],
    allowed_protein_ids: AbstractSet[ProteinId],
    allowed_subontologies: AbstractSet[Subontology],
) -> tuple[ProteinTermRecord, ...]:
    """Filter reference train-term rows by train taxa and selected ontologies."""

    return tuple(
        row
        for row in rows
        if row.protein_id in allowed_protein_ids
        and _ASPECT_TO_SUBONTOLOGY[row.aspect] in allowed_subontologies
    )


def filter_reference_sequence_records(
    records: tuple[SequenceRecord, ...],
    allowed_protein_ids: AbstractSet[ProteinId],
) -> tuple[SequenceRecord, ...]:
    """Filter reference FASTA records by allowed protein IDs."""

    return tuple(
        record
        for record in records
        if record.protein_id in allowed_protein_ids
    )


def _mapping_comparison_report(
    recreated_path: str | Path,
    reference_path: str | Path,
    message: str,
    left_mapping: dict[ProteinId, object],
    right_mapping: dict[ProteinId, object],
    formatter,
) -> ValidationReport:
    left_only_keys = sorted(set(left_mapping) - set(right_mapping))
    right_only_keys = sorted(set(right_mapping) - set(left_mapping))
    shared_mismatch_keys = sorted(
        protein_id
        for protein_id in set(left_mapping) & set(right_mapping)
        if left_mapping[protein_id] != right_mapping[protein_id]
    )

    sample_left_only = tuple(
        formatter(protein_id, left_mapping[protein_id])
        for protein_id in left_only_keys[:5]
    ) + tuple(
        formatter(protein_id, left_mapping[protein_id])
        for protein_id in shared_mismatch_keys[:5]
    )
    sample_right_only = tuple(
        formatter(protein_id, right_mapping[protein_id])
        for protein_id in right_only_keys[:5]
    ) + tuple(
        formatter(protein_id, right_mapping[protein_id])
        for protein_id in shared_mismatch_keys[:5]
    )

    passed = not left_only_keys and not right_only_keys and not shared_mismatch_keys
    return ValidationReport(
        left_path=Path(recreated_path),
        right_path=Path(reference_path),
        passed=passed,
        message="" if passed else message,
        sample_left_only=sample_left_only,
        sample_right_only=sample_right_only,
    )


def _group_terms_by_protein(
    rows: Iterable[ProteinTermRecord],
) -> dict[ProteinId, tuple[ProteinTermRecord, ...]]:
    grouped: dict[ProteinId, list[ProteinTermRecord]] = defaultdict(list)
    for row in rows:
        grouped[row.protein_id].append(row)
    return {
        protein_id: tuple(
            sorted(records, key=lambda record: (record.term_id, record.aspect))
        )
        for protein_id, records in grouped.items()
    }


def _is_close(
    left: float,
    right: float,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> bool:
    difference = abs(left - right)
    tolerance = max(absolute_tolerance, relative_tolerance * max(abs(left), abs(right)))
    return difference <= tolerance
