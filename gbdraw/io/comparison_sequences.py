#!/usr/bin/env python
# coding: utf-8

"""The one reader of comparison-genome sequence files (design D12).

A Circular similarity ring compares the displayed reference with a comparison
genome. Each comparison file is FASTA or a GenBank / DDBJ flat file; the format
is read from the content (``>`` starts FASTA, ``LOCUS`` starts a flat file,
which Biopython parses the same way for GenBank and DDBJ). A genome is one or
more files (a chromosome and plasmids downloaded separately): every record of
every file, in the listed order, is part of the LOSAT query set
(:func:`comparison_genome`).

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
from typing import Literal, Sequence

from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord

from gbdraw.exceptions import ValidationError
from gbdraw.io.input_text import GZIP_ERRORS, read_sequence_text

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
        text = read_sequence_text(text_path, "".join)
    except (OSError, UnicodeDecodeError, *GZIP_ERRORS) as exc:
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


@dataclass(frozen=True)
class RepeatedRecordIds:
    """Record IDs that each name more than one record of one comparison file (OV-2954).

    The ring draws every record, but a search names a hit by its query record
    ID, so the hover popup and the span export take the first record with that
    ID for each of them. Two files of one genome that share an ID are refused
    instead (:func:`comparison_genome`).
    """

    source_name: str
    record_ids: tuple[str, ...]
    # The file's path as the request named it: its identity, which the shown ``source_name`` is not.
    source_path: str

    @property
    def message(self) -> str:
        """The warning every surface shows (the CLI logs it, the Web shows it after Generate)."""

        ids = ", ".join(self.record_ids)
        named = f"record ID {ids} names" if len(self.record_ids) == 1 else f"record IDs {ids} each name"
        return (
            f"Comparison sequence file {self.source_name}: {named} more than one record. The ring draws "
            "every record, but the hover popup and the span export show the first record with that ID."
        )


@dataclass(frozen=True)
class ComparisonGenome:
    """The comparison genome of one ring: its files, in the listed order.

    ``records`` are the records of every file, file after file, as if the
    files were one file; ``label`` and ``record_label`` are those of the first
    file (the ring label rule applied to the first file).
    """

    files: tuple[ComparisonSequenceFile, ...]

    @property
    def path(self) -> str:
        return self.files[0].path

    @property
    def records(self) -> tuple[SeqRecord, ...]:
        return tuple(record for file in self.files for record in file.records)

    @property
    def label(self) -> str:
        return self.files[0].label

    @property
    def record_label(self) -> str | None:
        return self.files[0].record_label

    @property
    def repeated_record_ids(self) -> tuple[RepeatedRecordIds, ...]:
        """Each file with a record ID that names more than one of its records, in file order."""

        repeated = []
        for file in self.files:
            counts: dict[str, int] = {}
            for record in file.records:
                counts[str(record.id)] = counts.get(str(record.id), 0) + 1
            ids = tuple(record_id for record_id, count in counts.items() if count > 1)
            if ids:
                repeated.append(RepeatedRecordIds(os.path.basename(file.path), ids, file.path))
        return tuple(repeated)


def comparison_genome(files: Sequence[ComparisonSequenceFile]) -> ComparisonGenome:
    """One comparison genome of ``files``, refused when two of its files share a record ID.

    A search reports a hit by the query record ID, so an ID that names records
    of two files would make the ring's hits and its span export ambiguous
    (``INPUT_INVALID`` / ``UNIQUE_IDS``), and an empty file among others is
    refused (``INPUT_UNREADABLE`` / ``SEQUENCE_MISSING``). One file is read as before; an ID
    that names two records of one file is drawn and warned about
    (:attr:`ComparisonGenome.repeated_record_ids`, Owner 2026-11-09 「Q2=C」).
    """

    if not files:
        raise ValueError("A comparison genome has at least one file.")
    if len(files) > 1:
        # One empty file is a ring without a query (a LOSAT ring rejects it); among
        # other files it would leave the genome without saying so.
        for file in files:
            if not file.records:
                raise _unreadable(
                    f"Comparison sequence file {os.path.basename(file.path)} has no sequence record.",
                    reason="SEQUENCE_MISSING",
                )
    owner: dict[str, int] = {}
    for index, file in enumerate(files):
        for record in file.records:
            record_id = str(record.id)
            first = owner.setdefault(record_id, index)
            if first != index:
                raise ValidationError(
                    f"The comparison genome has record {record_id} in {files[first].path} and "
                    f"{file.path}; the files of one genome must not share a record ID.",
                    diagnostic={
                        "code": "INPUT_INVALID",
                        "reason": "UNIQUE_IDS",
                        "field": _FIELD,
                        "recordId": record_id,
                    },
                )
    return ComparisonGenome(tuple(files))


__all__ = [
    "ComparisonGenome",
    "ComparisonSequenceFile",
    "RepeatedRecordIds",
    "comparison_file_stem",
    "comparison_genome",
    "read_comparison_sequence_file",
]
