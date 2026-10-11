#!/usr/bin/env python
# coding: utf-8

"""The one reader of comparison-genome sequence files (design D12).

A Circular similarity ring compares the displayed reference with a comparison
genome. The comparison file is FASTA or a GenBank / DDBJ flat file; the format
is read from the content (``>`` starts FASTA, ``LOCUS`` starts a flat file,
which Biopython parses the same way for GenBank and DDBJ). One file is one
genome: every record is part of the LOSAT query set.

The CLI, the Python API and the Web worker helper call
:func:`read_comparison_sequence_file`. It returns records with only the ID and
the sequence, so a ring and its interactive span export do not depend on the
file format.
"""

from __future__ import annotations

from dataclasses import dataclass
import io
import os
import re
from typing import Literal

from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord

from gbdraw.exceptions import ValidationError

ComparisonSequenceFormat = Literal["fasta", "genbank"]

_FIELD = "comparison_sequence"


@dataclass(frozen=True)
class ComparisonSequenceFile:
    """The records of one comparison genome and its default ring label.

    ``record_label`` is the label the file names itself (GenBank / DDBJ: the
    first record's DEFINITION, else its organism), or ``None``; ``label`` falls
    back to :func:`comparison_file_stem` when it is ``None``.
    """

    path: str
    format: ComparisonSequenceFormat
    records: tuple[SeqRecord, ...]
    label: str
    record_label: str | None = None


def comparison_file_stem(path: str | os.PathLike[str]) -> str:
    """The file name without its last extension: ``genome.v2.fasta`` -> ``genome.v2``.

    The default label of a ring whose file names no label (D-03) and the stem
    of its raw LOSAT TSV name; the Web ``defaultConservationSeriesLabel`` uses
    the same rule (``tests/fixtures/comparison_ring_default_label_cases.json``).
    """

    return re.sub(r"\.[^.]+$", "", os.path.basename(os.fspath(path)))


def _default_label(name: str) -> str:
    return comparison_file_stem(name).strip() or name


def _unreadable(message: str, *, reason: str | None = None) -> ValidationError:
    diagnostic: dict[str, object] = {"code": "INPUT_UNREADABLE", "field": _FIELD}
    if reason is not None:
        diagnostic["reason"] = reason
    return ValidationError(message, diagnostic=diagnostic)


def _detect_format(text: str, path: str) -> ComparisonSequenceFormat:
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            return "fasta"
        if stripped.startswith("LOCUS"):
            return "genbank"
        break
    raise _unreadable(
        f"Comparison sequence file {os.path.basename(path)} is not FASTA (first line "
        "starting with '>') or a GenBank / DDBJ flat file (first line starting with "
        "'LOCUS')."
    )


def _sequence_text(record: SeqRecord) -> str:
    try:
        return str(record.seq)
    except UndefinedSequenceError:
        return ""


def _record_label(
    file_format: ComparisonSequenceFormat,
    records: tuple[SeqRecord, ...],
) -> str | None:
    """GenBank / DDBJ: the first record's DEFINITION, else its organism."""

    if file_format == "genbank" and records:
        first = records[0]
        definition = str(first.description or "").strip()
        if definition and definition != ".":
            return definition
        organism = str(first.annotations.get("organism") or "").strip()
        if organism and organism != ".":
            return organism
    return None


def read_comparison_sequence_file(path: str | os.PathLike[str]) -> ComparisonSequenceFile:
    """Read one comparison genome (FASTA, GenBank or DDBJ) from ``path``.

    An empty file has no records. Raises ``INPUT_UNREADABLE`` when the file
    cannot be read or parsed, and ``INPUT_UNREADABLE`` / ``SEQUENCE_MISSING``
    when a record has no sequence (an empty ``ORIGIN`` or a ``CONTIG``-only
    flat file).
    """

    text_path = os.fspath(path)
    name = os.path.basename(text_path)
    try:
        with open(text_path, "r", encoding="utf-8-sig") as handle:
            text = handle.read()
    except (OSError, UnicodeDecodeError) as exc:
        raise _unreadable(f"Could not read comparison sequence file {name}: {exc}") from exc
    if not text.strip():
        # An empty file names no genome; a LOSAT ring rejects it, a precomputed
        # ring's span export has no sequence to offer.
        return ComparisonSequenceFile(text_path, "fasta", (), _default_label(name))
    file_format = _detect_format(text, text_path)
    try:
        parsed = list(SeqIO.parse(io.StringIO(text), file_format))
    except ValueError as exc:
        raise _unreadable(
            f"Could not parse comparison sequence file {name} as {file_format}: {exc}"
        ) from exc
    if not parsed:
        raise _unreadable(f"Comparison sequence file {name} has no sequence record.")
    records: list[SeqRecord] = []
    for record in parsed:
        sequence = _sequence_text(record)
        if not sequence:
            raise _unreadable(
                f"Comparison sequence file {name}: record {record.id} has no sequence "
                "(empty ORIGIN or CONTIG only); use a file that contains the sequence.",
                reason="SEQUENCE_MISSING",
            )
        records.append(
            SeqRecord(
                Seq(sequence),
                id=str(record.id),
                name=str(record.id),
                description=str(record.description or ""),
            )
        )
    normalized = tuple(records)
    record_label = _record_label(file_format, tuple(parsed))
    return ComparisonSequenceFile(
        path=text_path,
        format=file_format,
        records=normalized,
        label=record_label or _default_label(name),
        record_label=record_label,
    )


def read_comparison_sequence_records(path: str | os.PathLike[str]) -> tuple[SeqRecord, ...]:
    """The records of :func:`read_comparison_sequence_file`."""

    return read_comparison_sequence_file(path).records


__all__ = [
    "ComparisonSequenceFile",
    "comparison_file_stem",
    "read_comparison_sequence_file",
    "read_comparison_sequence_records",
]
