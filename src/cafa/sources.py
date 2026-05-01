from __future__ import annotations

import shutil
import tarfile
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

from .config import ProjectConfig
from .types import SourceSnapshot
from .validation import validate_go_obo


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
        Metadata for the release-specific Swiss-Prot archive that contains the
        flatfile, FASTA, isoform FASTA, and XML payloads used by later train
        and test extraction work.
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


def download_source(snapshot: SourceSnapshot, chunk_size: int = 1024 * 1024) -> Path:
    """Download a source snapshot to its resolved local path.

    Existing files are reused as-is to keep notebook runs deterministic and to
    avoid unnecessary network traffic.
    """

    if snapshot.local_path.exists():
        return snapshot.local_path

    snapshot.local_path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(snapshot.url) as response, snapshot.local_path.open("wb") as handle:
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            handle.write(chunk)
    return snapshot.local_path


def extract_tar_gz_member(
    archive_path: str | Path,
    member_name: str,
    output_path: str | Path | None = None,
) -> Path:
    """Extract and cache one member from a `.tar.gz` archive.

    Parameters
    ----------
    archive_path:
        Local path to an existing `.tar.gz` archive.
    member_name:
        Archive member to extract.
    output_path:
        Optional explicit extraction target. When omitted, the extracted member
        is cached alongside the archive under a sibling directory whose name is
        derived from the archive stem.

    Returns
    -------
    Path
        Local path to the extracted member. Existing extracted files are reused
        as-is so only the raw unpacking work is cached.
    """

    archive_path = Path(archive_path)
    extracted_path = Path(output_path) if output_path is not None else _default_extracted_member_path(
        archive_path,
        member_name,
    )
    if extracted_path.exists():
        return extracted_path

    extracted_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:gz") as archive:
        try:
            member = archive.getmember(member_name)
        except KeyError as exc:
            raise FileNotFoundError(
                f"Archive {archive_path} does not contain {member_name}."
            ) from exc

        member_handle = archive.extractfile(member)
        if member_handle is None:
            raise FileNotFoundError(f"Unable to read {member_name} from {archive_path}.")

        with extracted_path.open("wb") as output_handle:
            shutil.copyfileobj(member_handle, output_handle)
    return extracted_path


def download_and_validate_go_obo(
    config: ProjectConfig,
    reference_path: str | Path,
) -> Path:
    """Download the pinned GO OBO file and validate it against the reference."""

    snapshot = resolve_go_obo_snapshot(config)
    downloaded_path = download_source(snapshot)
    report = validate_go_obo(downloaded_path, reference_path)
    if not report.passed:
        message = report.message or "Downloaded GO ontology validation failed."
        raise RuntimeError(message)
    return downloaded_path


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


def _default_extracted_member_path(archive_path: Path, member_name: str) -> Path:
    archive_stem = archive_path.name
    if archive_stem.endswith(".tar.gz"):
        archive_stem = archive_stem[: -len(".tar.gz")]
    return archive_path.parent / archive_stem / member_name
