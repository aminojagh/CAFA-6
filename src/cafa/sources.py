from __future__ import annotations

from hashlib import sha256
from pathlib import Path

from .config import ProjectConfig
from .types import SourceSnapshot


class ResearchRequiredError(RuntimeError):
    """Raised when a requested source chain is still blocked by a research gate."""


def _resolve_from_project_root(config: ProjectConfig, path: Path) -> Path:
    """Resolve `path` relative to `config.project_root` when not absolute."""

    if path.is_absolute():
        return path
    return config.project_root / path


def resolve_go_obo_snapshot(config: ProjectConfig) -> SourceSnapshot:
    """Resolve the pinned GO OBO source snapshot for a config.

    Parameters
    ----------
    config:
        Normalized project configuration.

    Returns
    -------
    SourceSnapshot
        Metadata for the `go-basic.obo` release, including its download URL and
        intended local cache path.
    """

    local_path = _resolve_from_project_root(
        config,
        config.cache_dir / "go" / config.go_release / "go-basic.obo",
    )
    return SourceSnapshot(
        name="go-basic.obo",
        url=f"https://release.geneontology.org/{config.go_release}/ontology/go-basic.obo",
        local_path=local_path,
        description="Pinned GO ontology release in OBO format.",
    )


def resolve_uniprot_sprot_snapshot(config: ProjectConfig) -> SourceSnapshot:
    """Resolve the pinned UniProt Swiss-Prot archive for a config.

    Parameters
    ----------
    config:
        Normalized project configuration.

    Returns
    -------
    SourceSnapshot
        Metadata for the archive that contains the release-specific Swiss-Prot
        FASTA payload used by later train and test extraction work.
    """

    release = config.train_uniprot_release
    filename = f"uniprot_sprot-only{release}.tar.gz"
    local_path = _resolve_from_project_root(
        config,
        config.cache_dir / "uniprot" / release / filename,
    )
    return SourceSnapshot(
        name="uniprot_sprot",
        url=(
            "https://ftp.uniprot.org/pub/databases/uniprot/"
            f"previous_releases/release-{release}/knowledgebase/{filename}"
        ),
        local_path=local_path,
        description="Pinned UniProt Swiss-Prot archive for the selected release.",
    )


def resolve_annotation_source_chain(config: ProjectConfig) -> tuple[SourceSnapshot, ...]:
    """Resolve the authoritative annotation source chain.

    This remains intentionally blocked because the project has not yet settled
    the official public-source chain for train annotations, benchmark evidence,
    and IA reconstruction.
    """

    raise ResearchRequiredError(
        "Annotation source-chain resolution is still blocked by research gates "
        "in PLANS.md (benchmark evidence policy and IA source chain)."
    )


def sha256_file(path: Path) -> str:
    """Compute the SHA-256 digest of a local file.

    Parameters
    ----------
    path:
        Local filesystem path to an existing file.

    Returns
    -------
    str
        Lowercase hexadecimal SHA-256 digest.
    """

    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()
