from __future__ import annotations

from pathlib import Path
from collections.abc import Iterable
from collections.abc import Callable
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
    """Validate recreated train taxonomy rows against the filtered reference.

    Train-taxonomy validation is intentionally asymmetric. At this stage the
    recreated taxonomy artifact is treated as a candidate superset gate for the
    later protein-to-term extraction step, so extra recreated proteins are
    tolerated. The validation is therefore considered passed when:

    - every filtered reference protein is present in the recreated mapping, and
    - every shared protein has the same taxon assignment.

    In report terms, that means `right_only_count == 0` and
    `shared_mismatch_count == 0`. Any remaining `left_only_count` is reported as
    diagnostic information but does not fail this stage.
    """

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

    report = _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message="Protein-to-taxon membership mismatch.",
        left_mapping=recreated_mapping,
        right_mapping=reference_mapping,
        comparison_unit="protein_id",
        formatter=lambda protein_id, taxon_id: f"{protein_id}\t{taxon_id}",
        pass_rule=lambda left_only, right_only, shared_mismatch: (
            not right_only and not shared_mismatch
        ),
    )
    return report


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
            comparison_unit="go_release",
            left_only_count=1,
            right_only_count=1,
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
        comparison_unit="go_id",
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

    recreated_mapping = _map_term_pairs_to_aspect(recreated_rows)
    reference_mapping = _map_term_pairs_to_aspect(filtered_reference_rows)

    return _mapping_comparison_report(
        recreated_path=recreated_path,
        reference_path=reference_path,
        message="Protein-to-direct-GO-term pair mapping mismatch.",
        left_mapping=recreated_mapping,
        right_mapping=reference_mapping,
        comparison_unit="protein_go_pair",
        formatter=lambda key, aspect: (
            f"{key[0]}\t{key[1]}\t{aspect}"
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
        comparison_unit="protein_id",
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
        comparison_unit="go_id",
        left_only_count=len(left_only_keys),
        right_only_count=len(right_only_keys),
        shared_mismatch_count=len(shared_mismatch_keys),
        sample_left_only=sample_left_only,
        sample_right_only=sample_right_only,
        sample_shared_mismatch_left=tuple(
            f"{term_id}\t{recreated_values[term_id]}"
            for term_id in shared_mismatch_keys[:5]
        ),
        sample_shared_mismatch_right=tuple(
            f"{term_id}\t{reference_values[term_id]}"
            for term_id in shared_mismatch_keys[:5]
        ),
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
    left_mapping: dict[object, object],
    right_mapping: dict[object, object],
    comparison_unit: str,
    formatter,
    pass_rule: Callable[[list[object], list[object], list[object]], bool] | None = None,
) -> ValidationReport:
    left_only_keys = sorted(set(left_mapping) - set(right_mapping))
    right_only_keys = sorted(set(right_mapping) - set(left_mapping))
    shared_mismatch_keys = sorted(
        protein_id
        for protein_id in set(left_mapping) & set(right_mapping)
        if left_mapping[protein_id] != right_mapping[protein_id]
    )

    sample_left_only = tuple(
        formatter(key, left_mapping[key])
        for key in left_only_keys[:5]
    )
    sample_right_only = tuple(
        formatter(key, right_mapping[key])
        for key in right_only_keys[:5]
    )
    sample_shared_mismatch_left = tuple(
        formatter(key, left_mapping[key])
        for key in shared_mismatch_keys[:5]
    )
    sample_shared_mismatch_right = tuple(
        formatter(key, right_mapping[key])
        for key in shared_mismatch_keys[:5]
    )

    passed = (
        pass_rule(left_only_keys, right_only_keys, shared_mismatch_keys)
        if pass_rule is not None
        else (not left_only_keys and not right_only_keys and not shared_mismatch_keys)
    )
    return ValidationReport(
        left_path=Path(recreated_path),
        right_path=Path(reference_path),
        passed=passed,
        message="" if passed else message,
        comparison_unit=comparison_unit,
        left_only_count=len(left_only_keys),
        right_only_count=len(right_only_keys),
        shared_mismatch_count=len(shared_mismatch_keys),
        sample_left_only=sample_left_only,
        sample_right_only=sample_right_only,
        sample_shared_mismatch_left=sample_shared_mismatch_left,
        sample_shared_mismatch_right=sample_shared_mismatch_right,
    )


def _map_term_pairs_to_aspect(
    rows: Iterable[ProteinTermRecord],
) -> dict[tuple[ProteinId, str], str]:
    pair_to_aspect: dict[tuple[ProteinId, str], str] = {}
    for row in rows:
        key = (row.protein_id, row.term_id)
        existing_aspect = pair_to_aspect.get(key)
        if existing_aspect is not None and existing_aspect != row.aspect:
            raise ValueError(
                f"Conflicting aspects for protein-term pair {key[0]} / {key[1]}: "
                f"{existing_aspect} vs {row.aspect}."
            )
        pair_to_aspect[key] = row.aspect
    return pair_to_aspect


def _is_close(
    left: float,
    right: float,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> bool:
    difference = abs(left - right)
    tolerance = max(absolute_tolerance, relative_tolerance * max(abs(left), abs(right)))
    return difference <= tolerance
