from __future__ import annotations

from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import unittest

from cafa.config import ProjectConfig
from cafa.sources import (
    ResearchRequiredError,
    resolve_annotation_source_chain,
    resolve_go_obo_snapshot,
    resolve_uniprot_sprot_snapshot,
    sha256_file,
)


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

    def test_resolve_annotation_source_chain_raises(self) -> None:
        config = build_config(Path("/tmp/project"))
        with self.assertRaises(ResearchRequiredError):
            resolve_annotation_source_chain(config)


if __name__ == "__main__":
    unittest.main()