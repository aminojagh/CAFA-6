from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from cafa.io import (
    benchmark_output_dir,
    build_recreated_layout,
    read_fasta_records,
    read_ia_values,
    read_test_taxon_rows,
    read_train_taxonomy_rows,
    read_train_term_rows,
    write_sequences,
    write_test_taxon_rows,
    write_train_taxonomy,
    write_train_terms,
)
from cafa.ontology import read_go_obo
from cafa.types import ProteinTaxonRecord, ProteinTermRecord, SequenceRecord, TaxonRecord

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


class ReaderTests(unittest.TestCase):
    def test_read_fasta_records_round_trips_writer_output(self) -> None:
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
            parsed = read_fasta_records(output_path)

            self.assertEqual(
                parsed,
                (
                    SequenceRecord(
                        protein_id="P1",
                        header=">sp|P1|FIRST first protein",
                        sequence="A" * 65,
                    ),
                    SequenceRecord(
                        protein_id="P2",
                        header=">sp|P2|SECOND second protein",
                        sequence="B" * 5,
                    ),
                ),
            )

    def test_read_train_taxonomy_rows_supports_reference_file_without_header(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "train_taxonomy.tsv"
            path.write_text("A0A\t9606\nQ9Z\t10090\n", encoding="utf-8")

            self.assertEqual(
                read_train_taxonomy_rows(path),
                (
                    ProteinTaxonRecord(protein_id="A0A", taxon_id=9606),
                    ProteinTaxonRecord(protein_id="Q9Z", taxon_id=10090),
                ),
            )

    def test_read_train_term_rows_can_canonicalize_with_ontology(self) -> None:
        with TemporaryDirectory() as tmpdir:
            obo_path = Path(tmpdir) / "mini.obo"
            obo_path.write_text(TEST_OBO, encoding="utf-8")
            ontology = read_go_obo(obo_path)

            path = Path(tmpdir) / "train_terms.tsv"
            path.write_text(
                "EntryID\tterm\taspect\nP1\tGO:0009002\tf\nP2\tGO:0008150\tP\n",
                encoding="utf-8",
            )

            self.assertEqual(
                read_train_term_rows(path, ontology=ontology),
                (
                    ProteinTermRecord(protein_id="P1", term_id="GO:0000002", aspect="F"),
                    ProteinTermRecord(protein_id="P2", term_id="GO:0008150", aspect="P"),
                ),
            )

    def test_read_test_taxon_rows_supports_reference_header_shape(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "testsuperset-taxon-list.tsv"
            path.write_text(
                "ID\tSpecies\n9606\tHomo sapiens\n10116\tRattus norvegicus\n",
                encoding="utf-8",
            )

            self.assertEqual(
                read_test_taxon_rows(path),
                (
                    TaxonRecord(taxon_id=9606, species_name="Homo sapiens"),
                    TaxonRecord(taxon_id=10116, species_name="Rattus norvegicus"),
                ),
            )

    def test_read_ia_values_reads_two_column_file(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "IA.tsv"
            path.write_text(
                "GO:0000001\t0.0\nGO:0000012\t6.038630248164372\n",
                encoding="utf-8",
            )

            self.assertEqual(
                read_ia_values(path),
                {
                    "GO:0000001": 0.0,
                    "GO:0000012": 6.038630248164372,
                },
            )


if __name__ == "__main__":
    unittest.main()
