from __future__ import annotations

from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import unittest
from unittest.mock import patch

from cafa.config import ProjectConfig
from cafa.sources import (
    ResearchRequiredError,
    download_and_validate_go_obo,
    download_source,
    resolve_annotation_source_chain,
    resolve_go_obo_snapshot,
    resolve_uniprot_sprot_snapshot,
    sha256_file,
)
from cafa.types import SourceSnapshot

TEST_OBO = """format-version: 1.2
ontology: go
name: mini-go
data-version: releases/2025-06-01

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function
"""


def build_config(project_root: Path) -> ProjectConfig:
    return ProjectConfig(
        go_release="2025-06-01",
        train_uniprot_release="2025_03",
        submission_deadline=date(2026, 1, 27),
        evaluation_time=date(2026, 3, 17),
        train_taxon_ids=(),
        test_taxon_ids=(),
        subontologies=("MF", "BP", "CC"),
        evidence_codes=("EXP", "IDA"),
        similarity_backend="diamond",
        validation_mode="canonical",
        project_root=project_root,
        cache_dir=Path(".cache/cafa"),
        recreated_data_dir=Path("recreated_comp_data"),
        artifacts_dir=Path("artifacts"),
        results_dir=Path("results"),
    )

class SourceResolutionTests(unittest.TestCase):
    def test_resolve_go_obo_snapshot(self) -> None:
        config = build_config(Path("/tmp/project"))
        snapshot = resolve_go_obo_snapshot(config)
        self.assertEqual(snapshot.name, "go-basic.obo")
        self.assertEqual(
            snapshot.url,
            "https://release.geneontology.org/2025-06-01/ontology/go-basic.obo",
        )
        self.assertEqual(
            snapshot.local_path,
            Path("/tmp/project/.cache/cafa/go/2025-06-01/go-basic.obo"),
        )

    def test_resolve_uniprot_sprot_snapshot(self) -> None:
        config = build_config(Path("/tmp/project"))
        snapshot = resolve_uniprot_sprot_snapshot(config)
        self.assertEqual(snapshot.name, "uniprot_sprot")
        self.assertEqual(
            snapshot.url,
            "https://ftp.uniprot.org/pub/databases/uniprot/"
            "previous_releases/release-2025_03/knowledgebase/"
            "uniprot_sprot-only2025_03.tar.gz",
        )
        self.assertEqual(
            snapshot.local_path,
            Path("/tmp/project/.cache/cafa/uniprot/2025_03/uniprot_sprot-only2025_03.tar.gz"),
        )

    def test_sha256_file(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.txt"
            payload = b"cafa-source-resolution\n"
            path.write_bytes(payload)
            self.assertEqual(sha256_file(path), hashlib.sha256(payload).hexdigest())

    def test_download_source_downloads_file_url(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_path = root / "source.txt"
            source_path.write_text("hello go\n", encoding="utf-8")
            target_path = root / "cache" / "downloaded.txt"
            snapshot = SourceSnapshot(
                name="sample",
                url=source_path.resolve().as_uri(),
                local_path=target_path,
                description="sample file",
            )

            downloaded_path = download_source(snapshot, chunk_size=4)

            self.assertEqual(downloaded_path, target_path)
            self.assertEqual(target_path.read_text(encoding="utf-8"), "hello go\n")

    def test_download_and_validate_go_obo_uses_downloaded_snapshot(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_path = root / "source_go.obo"
            source_path.write_text(TEST_OBO, encoding="utf-8")
            reference_path = root / "reference_go.obo"
            reference_path.write_text(TEST_OBO, encoding="utf-8")
            target_path = root / "cache" / "go-basic.obo"
            snapshot = SourceSnapshot(
                name="go-basic.obo",
                url=source_path.resolve().as_uri(),
                local_path=target_path,
                description="test go file",
            )

            with patch("cafa.sources.resolve_go_obo_snapshot", return_value=snapshot):
                downloaded_path = download_and_validate_go_obo(
                    build_config(root),
                    reference_path,
                )

            self.assertEqual(downloaded_path, target_path)
            self.assertEqual(target_path.read_text(encoding="utf-8"), TEST_OBO)

    def test_resolve_annotation_source_chain_raises(self) -> None:
        config = build_config(Path("/tmp/project"))
        with self.assertRaises(ResearchRequiredError):
            resolve_annotation_source_chain(config)


if __name__ == "__main__":
    unittest.main()
