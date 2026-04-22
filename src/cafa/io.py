from __future__ import annotations

from pathlib import Path

from .types import (
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
