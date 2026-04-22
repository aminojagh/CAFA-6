from __future__ import annotations

from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from cafa.config import ProjectConfig
from cafa.ontology import read_go_obo
from cafa.types import ProteinTaxonRecord, ProteinTermRecord, SequenceRecord
from cafa.validation import (
    filter_reference_sequence_records,
    filter_reference_train_taxonomy_rows,
    filter_reference_train_term_rows,
    validate_go_obo,
    validate_ia_values,
    validate_sequence_mapping,
    validate_train_taxonomy,
    validate_train_terms,
)

TEST_OBO = """format-version: 1.2
ontology: go
name: mini-go
data-version: releases/2025-06-01

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0000002
name: child term
namespace: molecular_function
alt_id: GO:0009002
is_a: GO:0003674 ! molecular_function
"""


def build_config(project_root: Path) -> ProjectConfig:
    return ProjectConfig(
        go_release="2025-06-01",
        train_uniprot_release="2025_03",
        submission_deadline=date(2026, 1, 27),
        evaluation_time=date(2026, 3, 17),
        train_taxon_ids=(9606,),
        test_taxon_ids=(),
        subontologies=("MF",),
        evidence_codes=("EXP", "IDA"),
        similarity_backend="diamond",
        validation_mode="canonical",
        project_root=project_root,
        cache_dir=Path(".cache/cafa"),
        recreated_data_dir=Path("recreated_comp_data"),
        artifacts_dir=Path("artifacts"),
        results_dir=Path("results"),
    )


class FilterTests(unittest.TestCase):
    def test_filter_reference_train_taxonomy_rows(self) -> None:
        rows = (
            ProteinTaxonRecord(protein_id="P1", taxon_id=9606),
            ProteinTaxonRecord(protein_id="P2", taxon_id=10090),
        )
        self.assertEqual(
            filter_reference_train_taxonomy_rows(rows, allowed_taxon_ids={9606}),
            (ProteinTaxonRecord(protein_id="P1", taxon_id=9606),),
        )

    def test_filter_reference_train_term_rows(self) -> None:
        rows = (
            ProteinTermRecord(protein_id="P1", term_id="GO:0000002", aspect="F"),
            ProteinTermRecord(protein_id="P1", term_id="GO:0008150", aspect="P"),
            ProteinTermRecord(protein_id="P2", term_id="GO:0000002", aspect="F"),
        )
        self.assertEqual(
            filter_reference_train_term_rows(
                rows,
                allowed_protein_ids={"P1"},
                allowed_subontologies={"MF"},
            ),
            (ProteinTermRecord(protein_id="P1", term_id="GO:0000002", aspect="F"),),
        )

    def test_filter_reference_sequence_records(self) -> None:
        records = (
            SequenceRecord(protein_id="P1", header=">P1", sequence="AAAA"),
            SequenceRecord(protein_id="P2", header=">P2", sequence="BBBB"),
        )
        self.assertEqual(
            filter_reference_sequence_records(records, allowed_protein_ids={"P2"}),
            (SequenceRecord(protein_id="P2", header=">P2", sequence="BBBB"),),
        )


class ValidationTests(unittest.TestCase):
    def test_validate_go_obo_passes_for_structurally_equal_files(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated.obo"
            reference_path = root / "reference.obo"
            recreated_path.write_text(TEST_OBO, encoding="utf-8")
            reference_path.write_text(TEST_OBO, encoding="utf-8")

            report = validate_go_obo(recreated_path, reference_path)

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")

    def test_validate_go_obo_reports_release_mismatch(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated.obo"
            reference_path = root / "reference.obo"
            recreated_path.write_text(TEST_OBO.replace("2025-06-01", "2025-06-02"), encoding="utf-8")
            reference_path.write_text(TEST_OBO, encoding="utf-8")

            report = validate_go_obo(recreated_path, reference_path)

            self.assertFalse(report.passed)
            self.assertIn("release", report.message.lower())

    def test_validate_train_taxonomy_passes_after_reference_filtering(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated_train_taxonomy.tsv"
            reference_path = root / "reference_train_taxonomy.tsv"
            recreated_path.write_text("EntryID\ttaxon_id\nP1\t9606\n", encoding="utf-8")
            reference_path.write_text("P1\t9606\nP2\t10090\n", encoding="utf-8")

            report = validate_train_taxonomy(recreated_path, reference_path, build_config(root))

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")

    def test_validate_train_terms_passes_after_reference_filtering(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            obo_path = root / "mini.obo"
            obo_path.write_text(TEST_OBO, encoding="utf-8")
            ontology = read_go_obo(obo_path)

            recreated_path = root / "train_terms.tsv"
            reference_path = root / "reference_train_terms.tsv"
            recreated_taxonomy_rows = (ProteinTaxonRecord(protein_id="P1", taxon_id=9606),)

            recreated_path.write_text(
                "EntryID\tterm\taspect\nP1\tGO:0000002\tF\n",
                encoding="utf-8",
            )
            reference_path.write_text(
                "EntryID\tterm\taspect\nP1\tGO:0009002\tF\nP2\tGO:0000002\tF\nP1\tGO:0008150\tP\n",
                encoding="utf-8",
            )

            report = validate_train_terms(
                recreated_path,
                reference_path,
                build_config(root),
                ontology,
                taxonomy_rows=recreated_taxonomy_rows,
            )

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")

    def test_validate_sequence_mapping_reports_mismatch(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated.fasta"
            reference_path = root / "reference.fasta"
            recreated_path.write_text(">sp|P1|ONE\nAAAA\n", encoding="utf-8")
            reference_path.write_text(">sp|P1|ONE\nAAAT\n", encoding="utf-8")

            report = validate_sequence_mapping(
                recreated_path,
                reference_path,
                artifact_name="Train/train_sequences.fasta",
            )

            self.assertFalse(report.passed)
            self.assertIn("P1", report.sample_left_only[0])
            self.assertIn("mismatch", report.message.lower())

    def test_validate_ia_values_passes_within_tolerance(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated_ia.tsv"
            reference_path = root / "reference_ia.tsv"
            recreated_path.write_text("GO:0000001\t1.0000000001\n", encoding="utf-8")
            reference_path.write_text("GO:0000001\t1.0\n", encoding="utf-8")

            report = validate_ia_values(
                recreated_path,
                reference_path,
                relative_tolerance=1e-9,
                absolute_tolerance=1e-9,
            )

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")

    def test_validate_ia_values_reports_term_or_value_mismatch(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recreated_path = root / "recreated_ia.tsv"
            reference_path = root / "reference_ia.tsv"
            recreated_path.write_text("GO:0000001\t1.5\n", encoding="utf-8")
            reference_path.write_text("GO:0000001\t1.0\nGO:0000002\t2.0\n", encoding="utf-8")

            report = validate_ia_values(
                recreated_path,
                reference_path,
                relative_tolerance=1e-12,
                absolute_tolerance=1e-12,
            )

            self.assertFalse(report.passed)
            self.assertIn("mismatch", report.message.lower())
            self.assertTrue(report.sample_left_only)
            self.assertTrue(report.sample_right_only)


if __name__ == "__main__":
    unittest.main()
