"""The segments of a GFA 1 assembly graph, each a record with its sequence and no features (D-86, D-126).

The one reader of a GFA file: the CLI (``--gfa``), the Python API
(``read_gfa``, ``GfaInputSource``), the typed request and the Web's Worker
helper read GFA through :func:`gbdraw.io.genome.load_gfa`, which reads here.

Only S (segment) and L (link) lines are read. A segment's name is its record
ID, and segments keep their order in the file. A segment with a link from its
end back to its own start (``L s + s +`` or ``L s - s -``) is circular; every
other segment is linear. Where that link says the two ends overlap by n bases
(overlap ``nM``), the record drops its last n bases, so each base of the circle
is drawn once. Other links are not drawn. P (path) and W (walk) lines are
counted, not drawn: the input summary names them once per file. A GFA 2 file,
a segment without sequence (``*``) and a self-loop overlap in any other CIGAR
form are refused with the reason.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterator
from dataclasses import dataclass

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..exceptions import ValidationError
from .input_text import read_sequence_text

# The annotation of a GFA record that counts the P and W lines of its file, which are not drawn
# (``ResolvedRecordProvenance.graph_lines_not_drawn``; the input summary names them).
GFA_PATH_LINES_ANNOTATION = "gbdraw_gfa_path_lines"
# Record types only GFA 2 has; a file with one of them, or whose header says VN:Z:2, or whose S line has a length where GFA 1 has the sequence, is GFA 2.
_GFA2_RECORD_TYPES = frozenset("EFGOU")
_GFA2_VERSION = re.compile(r"\tVN:Z:2(\.|\t|$)")
_MATCH_OVERLAP = re.compile(r"([0-9]+)M")
_ORIENTATIONS = frozenset("+-")


def _refusal(message: str, reason: str, *, code: str = "INPUT_INVALID", **context: object) -> ValidationError:
    return ValidationError(message, diagnostic={"code": code, "reason": reason, **context})


@dataclass
class _Segment:
    name: str
    sequence: str
    # The overlap of its self-loop, once a link from its end back to its start is read.
    loop_overlap: int | None = None


def _loop_overlap(overlap: str, name: str, file_name: str) -> int:
    """The bases a self-loop's overlap repeats: n for ``nM``; none for ``*`` or no overlap; any other form refused."""

    if overlap in ("", "*"):
        return 0
    match = _MATCH_OVERLAP.fullmatch(overlap)
    if match is None:
        raise _refusal(
            f"GFA segment {name} in {file_name} links its end to its start with the overlap {overlap!r}: gbdraw reads "
            "a self-loop overlap of the form nM (n matching bases), * or none.",
            "GFA_SELF_LOOP_OVERLAP", recordName=name,
        )
    return int(match.group(1))


def _segment_lines(lines: Iterator[str], file_name: str) -> tuple[list[_Segment], int]:
    """The segments of a GFA 1 file in file order, each with its self-loop overlap, and the count of P and W lines."""

    segments: dict[str, _Segment] = {}
    loops: list[tuple[str, str, int]] = []
    path_lines = 0
    for line_number, line in enumerate(lines, start=1):
        fields = line.rstrip("\r\n").split("\t")
        kind = fields[0]
        if (
            kind == "H" and _GFA2_VERSION.search(line.rstrip("\r\n"))
            or kind in _GFA2_RECORD_TYPES
            or kind == "S" and len(fields) > 2 and fields[2].isdigit()  # GFA 2: S <sid> <slen> <sequence>
        ):
            raise _refusal(
                f"{file_name} is a GFA 2 file: gbdraw reads GFA 1 (S and L lines).", "GFA_VERSION", code="INPUT_UNREADABLE"
            )
        if kind in ("P", "W"):
            path_lines += 1
        elif kind == "S":
            if len(fields) < 3 or not fields[1] or not fields[2]:
                raise _refusal(
                    f"Line {line_number} of {file_name} is an S line without a segment name and sequence.",
                    "GFA_LINE", code="INPUT_UNREADABLE", row=line_number,
                )
            name, sequence = fields[1], fields[2]
            if sequence == "*":
                raise _refusal(
                    f"GFA segment {name} in {file_name} has no sequence in the file (*). "
                    "Write the graph with its segment sequences.",
                    "GFA_NO_SEQUENCE", recordName=name,
                )
            if name in segments:
                raise _refusal(
                    f"{file_name} has more than one segment named {name}. Give each segment a unique name.",
                    "UNIQUE_IDS", recordName=name,
                )
            segments[name] = _Segment(name, sequence)
        elif kind == "L":
            if len(fields) < 5:
                raise _refusal(
                    f"Line {line_number} of {file_name} is an L line without its two segments and their orientations.",
                    "GFA_LINE", code="INPUT_UNREADABLE", row=line_number,
                )
            source, source_side, target, target_side = fields[1:5]
            if source == target and source_side == target_side and source_side in _ORIENTATIONS:
                loops.append((source, fields[5] if len(fields) > 5 else "", line_number))
    for name, overlap, _line_number in loops:
        segment = segments.get(name)
        if segment is None:
            continue
        bases = _loop_overlap(overlap, name, file_name)
        if segment.loop_overlap is not None and segment.loop_overlap != bases:
            raise _refusal(
                f"GFA segment {name} in {file_name} has self-loops with different overlaps "
                f"({segment.loop_overlap} and {bases} bases).",
                "GFA_SELF_LOOP_OVERLAP", recordName=name,
            )
        if bases >= len(segment.sequence):
            raise _refusal(
                f"GFA segment {name} in {file_name} is {len(segment.sequence):,} bp, and its self-loop overlaps "
                f"{bases:,} bases: nothing would remain to draw.",
                "GFA_SELF_LOOP_OVERLAP", recordName=name,
            )
        segment.loop_overlap = bases
    return list(segments.values()), path_lines


def read_gfa_records(path: str) -> list[SeqRecord]:
    """The segments of the GFA 1 file ``path``, in file order: each a record with its sequence and no features.

    A circular segment (one with a self-loop) has the topology ``circular`` and
    drops the bases its self-loop overlap repeats; every other segment is
    ``linear``. Each record counts the file's P and W lines
    (:data:`GFA_PATH_LINES_ANNOTATION`). A file without an S line has no record.
    """

    file_name = os.path.basename(path)
    segments, path_lines = read_sequence_text(path, lambda lines: _segment_lines(lines, file_name))
    if not segments:
        raise _refusal(f"{file_name} has no segment (S line): it is not a GFA 1 file.", "GFA_NO_SEGMENTS", code="NO_RECORDS")
    records = []
    for segment in segments:
        bases = segment.loop_overlap or 0
        record = SeqRecord(
            Seq(segment.sequence[: len(segment.sequence) - bases]), id=segment.name, name=segment.name, description=""
        )
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear" if segment.loop_overlap is None else "circular"
        record.annotations[GFA_PATH_LINES_ANNOTATION] = path_lines
        records.append(record)
    return records


__all__ = ["GFA_PATH_LINES_ANNOTATION", "read_gfa_records"]
