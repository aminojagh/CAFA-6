from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from cafa.io import (
    benchmark_output_dir,
    build_recreated_layout,
    write_sequences,
    write_test_taxon_rows,
    write_train_taxonomy,
    write_train_terms,
)
from cafa.types import ProteinTaxonRecord, ProteinTermRecord, SequenceRecord, TaxonRecord


class RecreatedLayoutTests(unittest.TestCase):
    def test_build_recreated_layout_creates_expected_directories(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "recreated_comp_data"
            layout = build_recreated_layout(root)

            self.assertEqual(layout.root_dir, root)
            self.assertEqual(layout.train_dir, root / "Train")
            self.assertEqual(layout.test_dir, root / "Test")
            self.assertEqual(layout.pred_dir, root / "Pred")
            self.assertEqual(layout.benchmark_dir, root / "Benchmark")
            self.assertEqual(layout.ia_path, root / "IA.tsv")

            self.assertTrue(layout.train_dir.is_dir())
            self.assertTrue(layout.test_dir.is_dir())
            self.assertTrue(layout.pred_dir.is_dir())
            self.assertTrue(layout.benchmark_dir.is_dir())

    def test_benchmark_output_dir_creates_subontology_directory(self) -> None:
        with TemporaryDirectory() as tmpdir:
            layout = build_recreated_layout(Path(tmpdir) / "recreated_comp_data")
            output_dir = benchmark_output_dir(layout, "MF")
            self.assertEqual(output_dir, layout.benchmark_dir / "MF")
            self.assertTrue(output_dir.is_dir())


class WriterTests(unittest.TestCase):
    def test_write_train_taxonomy_sorts_rows_and_writes_header(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "Train" / "train_taxonomy.tsv"
            rows = (
                ProteinTaxonRecord(protein_id="Q9Z", taxon_id=2),
                ProteinTaxonRecord(protein_id="A0A", taxon_id=9606),
            )

            write_train_taxonomy(rows, output_path)

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "EntryID\ttaxon_id\nA0A\t9606\nQ9Z\t2\n",
            )

    def test_write_train_terms_sorts_rows_and_writes_header(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "Train" / "train_terms.tsv"
            rows = (
                ProteinTermRecord(protein_id="P2", term_id="GO:0000002", aspect="P"),
                ProteinTermRecord(protein_id="P1", term_id="GO:0000003", aspect="F"),
                ProteinTermRecord(protein_id="P1", term_id="GO:0000002", aspect="C"),
            )

            write_train_terms(rows, output_path)

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "EntryID\tterm\taspect\n"
                "P1\tGO:0000002\tC\n"
                "P1\tGO:0000003\tF\n"
                "P2\tGO:0000002\tP\n",
            )

    def test_write_sequences_preserves_header_and_wraps_at_sixty(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "Train" / "train_sequences.fasta"
            records = (
                SequenceRecord(
                    protein_id="P2",
                    header="sp|P2|SECOND second protein",
                    sequence="B" * 5,
                ),
                SequenceRecord(
                    protein_id="P1",
                    header=">sp|P1|FIRST first protein",
                    sequence="A" * 65,
                ),
            )

            write_sequences(records, output_path)

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                ">sp|P1|FIRST first protein\n"
                + ("A" * 60)
                + "\n"
                + ("A" * 5)
                + "\n"
                + ">sp|P2|SECOND second protein\n"
                + ("B" * 5)
                + "\n",
            )

    def test_write_test_taxon_rows_sorts_rows_and_writes_header(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "Test" / "testsuperset-taxon-list.tsv"
            rows = (
                TaxonRecord(taxon_id=9606, species_name="Homo sapiens"),
                TaxonRecord(taxon_id=10090, species_name="Mus musculus"),
            )

            write_test_taxon_rows(tuple(reversed(rows)), output_path)

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "taxon_id\tspecies_name\n9606\tHomo sapiens\n10090\tMus musculus\n",
            )


if __name__ == "__main__":
    unittest.main()
