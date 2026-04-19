"""Project configuration domain model.

This module is part of the `G1` bundle and therefore intentionally keeps a
small scope: it defines the normalized project configuration object, but it does
not yet implement parsing or normalization helpers.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path

from .types import Subontology, TaxonId


@dataclass(frozen=True, slots=True)
class ProjectConfig:
    """Normalized project-wide configuration.

    This type represents the post-normalization state expected by later tasks.
    Future normalization code must enforce the constraints already specified in
    `PLANS.md`.

    Attributes
    ----------
    go_release:
        GO release date in `YYYY-MM-DD` format.
    train_uniprot_release:
        UniProt release identifier used for Swiss-Prot extraction.
    submission_deadline:
        `t0` date used for test/benchmark split logic.
    evaluation_time:
        `t1` date used for test/benchmark split logic.
    train_taxon_ids:
        Sorted unique train taxa chosen by the user.
    test_taxon_ids:
        Sorted unique test taxa chosen by the user.
    subontologies:
        Selected ontology branches to include. Allowed values are `MF`, `BP`,
        and `CC`.
    evidence_codes:
        Sorted unique evidence-code whitelist used by annotation extraction
        tasks.
    similarity_backend:
        Baseline sequence-retrieval backend. Initial expected value is
        `diamond`.
    validation_mode:
        Artifact-validation mode. Initial expected value is `canonical`.
    project_root:
        Root path of the repository.
    cache_dir:
        Relative or absolute cache directory.
    recreated_data_dir:
        Relative or absolute output directory for recreated artifacts.
    artifacts_dir:
        Relative or absolute model/artifact directory.
    results_dir:
        Relative or absolute results directory.
    """

    go_release: str
    train_uniprot_release: str
    submission_deadline: date
    evaluation_time: date
    train_taxon_ids: tuple[TaxonId, ...]
    test_taxon_ids: tuple[TaxonId, ...]
    subontologies: tuple[Subontology, ...]
    evidence_codes: tuple[str, ...]
    similarity_backend: str
    validation_mode: str
    project_root: Path
    cache_dir: Path
    recreated_data_dir: Path
    artifacts_dir: Path
    results_dir: Path
