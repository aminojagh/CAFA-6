"""Shared domain types for the CAFA project.

This module intentionally contains only immutable data containers and simple
type aliases. It is the narrow first implementation step because every later
module depends on these symbols.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, TypeAlias


ProteinId: TypeAlias = str
GoId: TypeAlias = str
TaxonId: TypeAlias = int
Subontology: TypeAlias = Literal["MF", "BP", "CC"]
AspectCode: TypeAlias = Literal["F", "P", "C"]


@dataclass(frozen=True, slots=True)
class SourceSnapshot:
    """Pinned public-source descriptor.

    Attributes
    ----------
    name:
        Stable internal name of the source.
    url:
        Remote URL used to download the source.
    local_path:
        Intended local cache path.
    description:
        Short human-readable description of the source payload.
    """

    name: str
    url: str
    local_path: Path
    description: str = ""


@dataclass(frozen=True, slots=True)
class ProteinTaxonRecord:
    """Association between one protein and one taxon.

    This is the same logical entity regardless of whether it later appears in a
    train-derived artifact, a test-derived artifact, or an intermediate table.
    """

    protein_id: ProteinId
    taxon_id: TaxonId


@dataclass(frozen=True, slots=True)
class ProteinTermRecord:
    """Association between one protein and one direct GO term.

    `term_id` is always intended to be the direct asserted GO annotation, not a
    propagated ancestor term.
    """

    protein_id: ProteinId
    term_id: GoId
    aspect: AspectCode


@dataclass(frozen=True, slots=True)
class TaxonRecord:
    """Taxon identifier and display name pair.

    This is the same logical entity regardless of whether it is later written to
    a test-superset taxon list or used as an intermediate taxonomy lookup.
    """

    taxon_id: TaxonId
    species_name: str


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """FASTA record container used by train and test artifact pipelines."""

    protein_id: ProteinId
    header: str
    sequence: str
    taxon_id: TaxonId | None = None


@dataclass(frozen=True, slots=True)
class ValidationMismatch:
    """One concrete mismatch discovered by artifact validation."""

    artifact_name: str
    kind: str
    message: str
    sample_left_only: tuple[str, ...] = ()
    sample_right_only: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Canonical validation result for one recreated artifact."""

    artifact_name: str
    left_path: Path
    right_path: Path
    passed: bool
    mode: str
    checked_properties: tuple[str, ...] = ()
    mismatches: tuple[ValidationMismatch, ...] = ()


@dataclass(frozen=True, slots=True)
class RecreatedLayout:
    """Canonical output layout for recreated CAFA-style artifacts.

    Attributes
    ----------
    root_dir:
        Root directory that contains all recreated artifacts.
    train_dir:
        Directory for recreated training artifacts.
    test_dir:
        Directory for recreated test-superset artifacts.
    pred_dir:
        Directory reserved for prediction-like artifacts.
    benchmark_dir:
        Directory reserved for ontology-specific benchmark outputs.
    ia_path:
        Canonical path for the recreated `IA.tsv` artifact.
    """

    root_dir: Path
    train_dir: Path
    test_dir: Path
    pred_dir: Path
    benchmark_dir: Path
    ia_path: Path


@dataclass(frozen=True, slots=True)
class BenchmarkBundle:
    """Ontology-specific benchmark artifacts grouped in memory.

    Attributes
    ----------
    ontology:
        The ontology branch this bundle belongs to: `MF`, `BP`, or `CC`.
    benchmark_ids:
        Protein IDs that qualify for scoring in this ontology.
    known_t0:
        Direct eligible GO annotations already known for each benchmark protein
        at `t0`. These are used to exclude prior knowledge from evaluation.
    ground_truth_delta:
        Direct eligible GO annotations newly accumulated in `(t0, t1]` for each
        benchmark protein. This is the actual ontology-specific ground truth to
        be predicted and scored.
    terms_of_interest:
        Canonical GO terms allowed to participate in scoring for this ontology.
    metadata:
        Extra bookkeeping fields, such as date parameters or source identifiers.
    """

    ontology: Subontology
    benchmark_ids: tuple[ProteinId, ...]
    known_t0: dict[ProteinId, tuple[GoId, ...]] = field(default_factory=dict)
    ground_truth_delta: dict[ProteinId, tuple[GoId, ...]] = field(default_factory=dict)
    terms_of_interest: tuple[GoId, ...] = ()
    metadata: dict[str, str] = field(default_factory=dict)

@dataclass(frozen=True, slots=True)
class CurvePoint:
    """One threshold-specific evaluation point."""

    threshold: float
    precision: float
    recall: float
    f1: float
    coverage: float
    weighted_precision: float | None = None
    weighted_recall: float | None = None
    weighted_f1: float | None = None


@dataclass(frozen=True, slots=True)
class EvaluationResult:
    """Summary of one ontology-specific evaluation run."""

    ontology: Subontology
    curve: tuple[CurvePoint, ...]
    best_f1: float
    best_threshold: float
    best_weighted_f1: float | None = None
    best_weighted_threshold: float | None = None
