from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO

from .ontology import GeneOntology, canonicalize_go_id
from .types import (
    GoId,
    ProteinTaxonRecord,
    ProteinTermRecord,
    RecreatedLayout,
    SequenceRecord,
    Subontology,
    TaxonRecord,
)


def build_recreated_layout(root: str | Path) -> RecreatedLayout:
    """Create and return the canonical recreated artifact layout.

    Parameters
    ----------
    root:
        Root directory under which recreated artifact subdirectories should be
        created.

    Returns
    -------
    RecreatedLayout
        Canonical layout with created directories for train, test, prediction,
        and benchmark outputs.
    """

    root_path = Path(root)
    train_dir = root_path / "Train"
    test_dir = root_path / "Test"
    pred_dir = root_path / "Pred"
    benchmark_dir = root_path / "Benchmark"

    for directory in (root_path, train_dir, test_dir, pred_dir, benchmark_dir):
        directory.mkdir(parents=True, exist_ok=True)

    return RecreatedLayout(
        root_dir=root_path,
        train_dir=train_dir,
        test_dir=test_dir,
        pred_dir=pred_dir,
        benchmark_dir=benchmark_dir,
        ia_path=root_path / "IA.tsv",
    )


def benchmark_output_dir(layout: RecreatedLayout, subontology: Subontology) -> Path:
    """Create and return the ontology-specific benchmark output directory."""

    output_dir = layout.benchmark_dir / subontology
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def write_train_taxonomy(
    rows: tuple[ProteinTaxonRecord, ...],
    output_path: str | Path,
) -> Path:
    """Write train taxonomy rows in deterministic TSV order."""

    header = "EntryID\ttaxon_id\n"
    ordered_rows = sorted(rows, key=lambda row: (row.protein_id, row.taxon_id))
    body = "".join(f"{row.protein_id}\t{row.taxon_id}\n" for row in ordered_rows)
    return _write_text(output_path, header + body)


def write_train_terms(
    rows: tuple[ProteinTermRecord, ...],
    output_path: str | Path,
) -> Path:
    """Write train term rows in deterministic TSV order."""

    header = "EntryID\tterm\taspect\n"
    ordered_rows = sorted(rows, key=lambda row: (row.protein_id, row.term_id, row.aspect))
    body = "".join(
        f"{row.protein_id}\t{row.term_id}\t{row.aspect}\n" for row in ordered_rows
    )
    return _write_text(output_path, header + body)


def write_sequences(
    records: tuple[SequenceRecord, ...],
    output_path: str | Path,
) -> Path:
    """Write FASTA records using verbatim headers and 60-character wrapping."""

    ordered_records = sorted(records, key=lambda record: record.protein_id)
    lines: list[str] = []
    for record in ordered_records:
        lines.append(_format_fasta_header(record.header))
        lines.extend(_wrap_sequence(record.sequence, width=60))
    text = "\n".join(lines)
    if text:
        text += "\n"
    return _write_text(output_path, text)


def write_test_taxon_rows(
    rows: tuple[TaxonRecord, ...],
    output_path: str | Path,
) -> Path:
    """Write test taxon rows in deterministic TSV order."""

    header = "taxon_id\tspecies_name\n"
    ordered_rows = sorted(rows, key=lambda row: (row.taxon_id, row.species_name))
    body = "".join(f"{row.taxon_id}\t{row.species_name}\n" for row in ordered_rows)
    return _write_text(output_path, header + body)


def read_fasta_records(path: str | Path) -> tuple[SequenceRecord, ...]:
    """Read FASTA records into canonical sequence containers."""

    records: list[SequenceRecord] = []
    with Path(path).open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header = f">{record.description}"
            records.append(
                SequenceRecord(
                    protein_id=_protein_id_from_fasta_header(header),
                    header=header,
                    sequence=str(record.seq),
                )
            )
    return tuple(records)


def read_train_taxonomy_rows(path: str | Path) -> tuple[ProteinTaxonRecord, ...]:
    """Read train taxonomy rows from a recreated or reference TSV artifact."""

    rows: list[ProteinTaxonRecord] = []
    for raw_row in _read_tsv_rows(path):
        if raw_row[0] in {"EntryID", "ID"}:
            continue
        rows.append(ProteinTaxonRecord(protein_id=raw_row[0], taxon_id=int(raw_row[1])))
    return tuple(rows)


def read_train_term_rows(
    path: str | Path,
    ontology: GeneOntology | None = None,
) -> tuple[ProteinTermRecord, ...]:
    """Read train term rows from a recreated or reference TSV artifact."""

    rows: list[ProteinTermRecord] = []
    for raw_row in _read_tsv_rows(path):
        if raw_row[0] in {"EntryID", "ID"}:
            continue
        term_id = raw_row[1]
        if ontology is not None:
            term_id = canonicalize_go_id(ontology, term_id)
        rows.append(
            ProteinTermRecord(
                protein_id=raw_row[0],
                term_id=term_id,
                aspect=raw_row[2].strip().upper(),
            )
        )
    return tuple(rows)


def read_test_taxon_rows(path: str | Path) -> tuple[TaxonRecord, ...]:
    """Read test taxon rows from a recreated or reference TSV artifact."""

    rows: list[TaxonRecord] = []
    for raw_row in _read_tsv_rows(path):
        if raw_row[0] in {"taxon_id", "ID"}:
            continue
        rows.append(TaxonRecord(taxon_id=int(raw_row[0]), species_name=raw_row[1]))
    return tuple(rows)


def read_ia_values(
    path: str | Path,
    ontology: GeneOntology | None = None,
) -> dict[GoId, float]:
    """Read IA values from a two-column TSV artifact."""

    values: dict[GoId, float] = {}
    for raw_row in _read_tsv_rows(path):
        if raw_row[0] in {"term", "go_id", "GO_ID"}:
            continue
        term_id = raw_row[0]
        if ontology is not None:
            term_id = canonicalize_go_id(ontology, term_id)
        values[term_id] = float(raw_row[1])
    return values


def _write_text(output_path: str | Path, text: str) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _format_fasta_header(header: str) -> str:
    if header.startswith(">"):
        return header
    return f">{header}"


def _wrap_sequence(sequence: str, width: int) -> list[str]:
    if not sequence:
        return [""]
    return [sequence[start : start + width] for start in range(0, len(sequence), width)]


def _read_tsv_rows(path: str | Path) -> list[list[str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


def _protein_id_from_fasta_header(header: str) -> str:
    bare_header = header[1:] if header.startswith(">") else header
    token = bare_header.split(maxsplit=1)[0]
    if token.startswith(("sp|", "tr|")):
        parts = token.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return token
