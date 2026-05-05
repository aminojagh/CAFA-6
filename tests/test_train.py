from __future__ import annotations

import gzip
import io
import tarfile
import unittest
from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory

from cafa.config import ProjectConfig
from cafa.io import write_train_taxonomy, write_train_terms
from cafa.ontology import read_go_obo
from cafa.train import (
    _count_swissprot_records,
    extract_train_taxonomy_records,
    extract_train_term_records,
)
from cafa.types import ProteinTaxonRecord, ProteinTermRecord, SourceSnapshot
from cafa.validation import validate_train_taxonomy, validate_train_terms

TEST_SWISSPROT_FLATFILE = """ID   HUMAN_TWO               Reviewed;         5 AA.
AC   Q99999;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MHHHH
//
ID   HUMAN_ONE               Reviewed;         5 AA.
AC   P12345; Q54321;
OX   NCBI_TaxID=9606;
DR   GO; GO:0003674; F:mock function; EXP:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MAAAA
//
ID   HUMAN_THREE             Reviewed;         5 AA.
AC   P54321;
OX   NCBI_TaxID=9606;
DR   GO; GO:0008150; P:mock process; IEA:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MCCCC
//
ID   MOUSE_ONE               Reviewed;         5 AA.
AC   P67890;
OX   NCBI_TaxID=10090;
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MBBBB
//
"""

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
name: child function
namespace: molecular_function
alt_id: GO:0009002
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000100
name: process term
namespace: biological_process
is_a: GO:0008150 ! biological_process
"""

TEST_SWISSPROT_TERM_FLATFILE = """ID   HUMAN_ONE               Reviewed;         5 AA.
AC   P12345; Q54321;
OX   NCBI_TaxID=9606;
DR   GO; GO:0009002; F:mock function; EXP:UniProtKB-KW.
DR   GO; GO:0000100; P:mock process; IDA:UniProtKB-KW.
DR   GO; GO:9999999; F:unknown term; EXP:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MAAAA
//
ID   HUMAN_TWO               Reviewed;         5 AA.
AC   P54321;
OX   NCBI_TaxID=9606;
DR   GO; GO:0000100; P:mock process; IEA:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MCCCC
//
ID   HUMAN_THREE             Reviewed;         5 AA.
AC   P77777;
OX   NCBI_TaxID=9606;
DR   GO; GO:0009002; F:mock function; EXP:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MDDDD
//
ID   MOUSE_ONE               Reviewed;         5 AA.
AC   P67890;
OX   NCBI_TaxID=10090;
DR   GO; GO:0009002; F:mock function; EXP:UniProtKB-KW.
SQ   SEQUENCE   5 AA;  555 MW;  ABCDEF CRC64;
     MBBBB
//
"""

# TODO: In order to keep files clean, put each mock input in a separate file


def build_config(
    project_root: Path,
    train_taxon_ids: tuple[int, ...],
    *,
    subontologies: tuple[str, ...] = ("MF", "BP", "CC"),
    evidence_codes: tuple[str, ...] = ("EXP", "IDA"),
) -> ProjectConfig:
    return ProjectConfig(
        go_release="2025-06-01",
        train_uniprot_release="2025_03",
        submission_deadline=date(2026, 1, 27),
        evaluation_time=date(2026, 3, 17),
        train_taxon_ids=train_taxon_ids,
        test_taxon_ids=(),
        subontologies=subontologies,
        evidence_codes=evidence_codes,
        similarity_backend="diamond",
        validation_mode="canonical",
        project_root=project_root,
        cache_dir=Path(".cache/cafa"),
        recreated_data_dir=Path("recreated_comp_data"),
        artifacts_dir=Path("artifacts"),
        results_dir=Path("results"),
    )


def write_swissprot_archive(path: Path, flatfile_text: str, member_name: str = "uniprot_sprot.dat.gz") -> None:
    compressed_payload = io.BytesIO()
    with gzip.GzipFile(fileobj=compressed_payload, mode="wb") as gzip_handle:
        gzip_handle.write(flatfile_text.encode("utf-8"))

    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo(member_name)
        member_bytes = compressed_payload.getvalue()
        member.size = len(member_bytes)
        archive.addfile(member, io.BytesIO(member_bytes))


class TrainTaxonomyExtractionTests(unittest.TestCase):
    def test_count_swissprot_records_counts_record_terminators(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_FLATFILE)
            with tarfile.open(archive_path, "r:gz") as archive:
                member_handle = archive.extractfile("uniprot_sprot.dat.gz")
                self.assertIsNotNone(member_handle)
                flatfile_gz_path = root / "uniprot_sprot.dat.gz"
                flatfile_gz_path.write_bytes(member_handle.read())

            self.assertEqual(_count_swissprot_records(flatfile_gz_path), 4)

    def test_extract_train_taxonomy_records_filters_selected_taxa(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_FLATFILE)
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )

            rows = extract_train_taxonomy_records(build_config(root, (9606,)), snapshot)

            self.assertEqual(
                rows,
                (
                    # The primary accession is used even when secondary accessions exist.
                    # Output order is deterministic by protein ID.
                    ProteinTaxonRecord(protein_id="P12345", taxon_id=9606),
                    ProteinTaxonRecord(protein_id="P54321", taxon_id=9606),
                    ProteinTaxonRecord(protein_id="Q99999", taxon_id=9606),
                ),
            )

    def test_extract_train_taxonomy_records_returns_empty_for_empty_taxa(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_FLATFILE)
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )

            rows = extract_train_taxonomy_records(build_config(root, ()), snapshot)

            self.assertEqual(rows, ())

    def test_extract_train_taxonomy_records_raises_when_flatfile_member_is_missing(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(
                archive_path,
                TEST_SWISSPROT_FLATFILE,
                member_name="unexpected_name.dat.gz",
            )
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )

            with self.assertRaises(FileNotFoundError):
                extract_train_taxonomy_records(build_config(root, (9606,)), snapshot)

    def test_extracted_train_taxonomy_validates_against_filtered_reference(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_FLATFILE)
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )
            config = build_config(root, (9606,))

            rows = extract_train_taxonomy_records(config, snapshot)
            recreated_path = write_train_taxonomy(rows, root / "Train" / "train_taxonomy.tsv")
            reference_path = root / "reference_train_taxonomy.tsv"
            reference_path.write_text(
                "P12345\t9606\nP67890\t10090\nQ99999\t9606\n",
                encoding="utf-8",
            )

            report = validate_train_taxonomy(recreated_path, reference_path, config)

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")
            self.assertEqual(report.left_only_count, 1)
            self.assertEqual(report.right_only_count, 0)
            self.assertEqual(report.shared_mismatch_count, 0)


class TrainTermExtractionTests(unittest.TestCase):
    def test_extract_train_term_records_filters_by_taxonomy_evidence_and_ontology(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_TERM_FLATFILE)
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )
            obo_path = root / "mini.obo"
            obo_path.write_text(TEST_OBO, encoding="utf-8")
            ontology = read_go_obo(obo_path)
            config = build_config(root, (9606,), subontologies=("MF", "BP"))

            rows = extract_train_term_records(
                config,
                snapshot,
                ontology,
                taxonomy_rows=(
                    ProteinTaxonRecord(protein_id="P12345", taxon_id=9606),
                    ProteinTaxonRecord(protein_id="P54321", taxon_id=9606),
                ),
            )

            self.assertEqual(
                rows,
                (
                    ProteinTermRecord(protein_id="P12345", term_id="GO:0000002", aspect="F"),
                    ProteinTermRecord(protein_id="P12345", term_id="GO:0000100", aspect="P"),
                ),
            )

    def test_extracted_train_terms_validate_against_reference(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive_path = root / "uniprot_sprot-only2025_03.tar.gz"
            write_swissprot_archive(archive_path, TEST_SWISSPROT_TERM_FLATFILE)
            snapshot = SourceSnapshot(
                name="uniprot_sprot",
                url=archive_path.resolve().as_uri(),
                local_path=archive_path,
                description="synthetic Swiss-Prot archive",
            )
            obo_path = root / "mini.obo"
            obo_path.write_text(TEST_OBO, encoding="utf-8")
            ontology = read_go_obo(obo_path)
            config = build_config(root, (9606,), subontologies=("MF", "BP"))
            taxonomy_rows = (
                ProteinTaxonRecord(protein_id="P12345", taxon_id=9606),
                ProteinTaxonRecord(protein_id="P54321", taxon_id=9606),
            )

            rows = extract_train_term_records(
                config,
                snapshot,
                ontology,
                taxonomy_rows=taxonomy_rows,
            )
            recreated_path = write_train_terms(rows, root / "Train" / "train_terms.tsv")
            reference_path = root / "reference_train_terms.tsv"
            reference_path.write_text(
                "EntryID\tterm\taspect\n"
                "P12345\tGO:0009002\tF\n"
                "P12345\tGO:0000100\tP\n"
                "P67890\tGO:0009002\tF\n",
                encoding="utf-8",
            )

            report = validate_train_terms(
                recreated_path,
                reference_path,
                config,
                ontology,
                taxonomy_rows=taxonomy_rows,
            )

            self.assertTrue(report.passed)
            self.assertEqual(report.message, "")
            self.assertEqual(report.left_only_count, 0)
            self.assertEqual(report.right_only_count, 0)
            self.assertEqual(report.shared_mismatch_count, 0)


if __name__ == "__main__":
    unittest.main()
