"""The drawn units of a figure: which record, or which window of it, each drawn record shows (D-01).

A figure draws its records as an ordered list of units. A unit names the record it draws (its key is shared by the
windows of one record), the part of that record it shows in source bases, whether it is drawn reverse-complemented,
and the record's length and topology. The list belongs to no diagram mode: a Linear row of windows reads its captions
and marks from it (``gbdraw.layout.linear_row_windows``), and the Linear figure keys its record translations and the
records the Web maps back by it.

This module is the one reader of those facts on the drawing side. The record planner writes them on each drawn record
(``gbdraw.api.record_planning``: the window span and the record a window is one of); the transform a record is drawn
by gives its orientation, length and topology.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import WINDOW_OF_ANNOTATION, record_identity_key, record_window_span
from gbdraw.exceptions import ValidationError
from gbdraw.layout.record_coordinates import CirclePlace, RecordDisplayTransform
from gbdraw.layout.record_labels import record_map_span
from gbdraw.layout.record_transform import bind_record_transforms


@dataclass(frozen=True)
class RecordUnit:
    """One drawn unit: the record, or the window of a record, that one drawn record shows.

    ``record_index`` is the unit's place among the drawn records. ``record_key`` is the key of the record it draws,
    shared by the windows of one record (``record_identity_key``). ``span`` is the part of the record the unit shows
    in source bases, start before end in the record's stored order (``start > end`` runs across the origin of a
    circular record); ``reverse`` is true when the unit is drawn reverse-complemented. ``is_window`` is false for a
    record drawn whole. ``window_index`` is the unit's place among the windows of its record when the record draws
    several, else ``None``.
    """

    record_index: int
    record_key: str
    record_id: str
    span: tuple[int, int]
    reverse: bool
    source_length: int
    circular: bool
    is_window: bool
    window_index: int | None

    @property
    def length(self) -> int:
        """The bases the unit shows."""

        start, end = self.span
        return end - start + 1 if start <= end else self.source_length - start + 1 + end

    def edge_bases(self) -> tuple[int, int]:
        """The source bases at the unit's drawn start and end."""

        start, end = self.span
        return (end, start) if self.reverse else (start, end)


def _record_unit(
    record: SeqRecord, record_index: int, transform: RecordDisplayTransform, window_index: int | None
) -> RecordUnit:
    window = record_window_span(record)
    span = window
    if span is None:
        low, high = sorted(record_map_span(record, len(record.seq)))
        span = (max(1, low), min(transform.length, high))
    return RecordUnit(
        record_index=record_index,
        record_key=str(record_identity_key(record) or f"record-{record_index + 1}"),
        record_id=str(record.id),
        span=(int(span[0]), int(span[1])),
        reverse=transform.source_step == -1,
        source_length=int(transform.length),
        circular=bool(transform.is_circular),
        is_window=window is not None,
        window_index=window_index,
    )


def record_units(
    records: Sequence[SeqRecord], transforms: Sequence[RecordDisplayTransform] | None
) -> tuple[RecordUnit, ...]:
    """The drawn units of ``records``, one per record in drawing order, read through the transforms they are drawn
    by (``bind_record_transforms``)."""

    windows_seen: dict[str, int] = {}
    units: list[RecordUnit] = []
    for index, (record, transform) in enumerate(zip(records, bind_record_transforms(records, transforms), strict=True)):
        window_of = (record.annotations or {}).get(WINDOW_OF_ANNOTATION)
        window_index = None
        if window_of is not None:
            window_index = windows_seen.get(str(window_of), 0)
            windows_seen[str(window_of)] = window_index + 1
        units.append(_record_unit(record, index, transform, window_index))
    return tuple(units)


@dataclass(frozen=True)
class CircleUnits:
    """Records drawn one after another on one circle, in drawing order, with a gap between neighbours.

    ``placed`` pairs each drawn unit with where it lies on the circle (``CirclePlace``: its first display base
    and the circle's display length, records and gaps together, which the drawing reads as the record length).
    A gap is ``gap_bp`` display bases; half a gap lies before the first unit and half after the last, so the
    gaps are equal round the circle. ``hidden`` holds the units shorter than the shortest drawn length: they are
    not drawn, and the figure says how many there are.
    """

    placed: tuple[tuple[RecordUnit, CirclePlace], ...]
    hidden: tuple[RecordUnit, ...]
    gap_bp: int

    @property
    def display_length(self) -> int:
        return self.placed[0][1].length


def left_off_circle(unit: RecordUnit, min_length_bp: int) -> bool:
    """Whether one circle leaves ``unit`` out: it is shorter than ``min_length_bp``, the shortest record drawn."""

    return unit.length < min_length_bp


def circle_units(units: Sequence[RecordUnit], *, gap_degrees: float, min_length_bp: int) -> CircleUnits:
    """Lay ``units`` out on one circle: each unit's arc is its length at one scale for the whole circle, and
    ``gap_degrees`` separates neighbours. Units shorter than ``min_length_bp`` are left out and counted."""

    if not math.isfinite(gap_degrees) or gap_degrees < 0:
        raise ValidationError(
            "The gap between records on one circle must be zero or more degrees.",
            diagnostic={"code": "INPUT_INVALID", "reason": "NONNEGATIVE"},
        )
    drawn = [unit for unit in units if not left_off_circle(unit, min_length_bp)]
    hidden = tuple(unit for unit in units if left_off_circle(unit, min_length_bp))
    if not drawn:
        raise ValidationError(
            f"Every record is shorter than {min_length_bp:,} bp, the shortest record drawn on one circle; "
            "lower that length to draw them.",
            diagnostic={"code": "NO_RECORDS"},
        )
    gaps = len(drawn)
    if gaps * gap_degrees >= 360.0:
        raise ValidationError(
            f"{gaps} records with gaps of {gap_degrees:g} degrees leave no room on one circle; use a smaller gap.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    bases = sum(unit.length for unit in drawn)
    gap_bp = math.ceil(gap_degrees * bases / (360.0 - gaps * gap_degrees)) if gap_degrees > 0 else 0
    display_length = bases + gaps * gap_bp
    spans: list[tuple[int, int]] = []
    offset = gap_bp // 2
    for unit in drawn:
        spans.append((offset, unit.length))
        offset += unit.length + gap_bp
    longest = max(unit.length for unit in drawn)
    record_spans = tuple(spans)
    placed = tuple(
        (unit, CirclePlace(start, display_length, longest, position, record_spans))
        for position, (unit, (start, _)) in enumerate(zip(drawn, spans, strict=True))
    )
    return CircleUnits(placed, hidden, gap_bp)


__all__ = ["CircleUnits", "RecordUnit", "circle_units", "left_off_circle", "record_units"]
