"""What a Linear row of record windows shows beside its windows (PD-OI-122).

A row that holds two or more drawn units, at least one of them a window with a
region, is a *windowed row*. It shows one name at its left and, under it, a
numbered list with one line per window; each window carries its number at its
upper left; each window edge is a bar (the first or last base of a linear
record) or a dashed continuation of the axis (an edge inside a record, or a
stored end of a circular record); between two windows of one record a "//"
stands in the gap, with the length left out above it when the two windows run
in the same direction, are next to each other in the record and leave an
interval between them.

This module is the one owner of those values: which rows are windowed, the row
name, the numbers, the list lines, the edge and gap marks, the omitted-length
text, the sizes of the marks and the room they take beside each window. The
layout and the drawing read them from here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.config.models import GbdrawConfig
from gbdraw.core.record_metadata import record_identity_key, record_window_span
from gbdraw.core.text import drawn_text_box
from gbdraw.layout.linear import VerticalBand
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_labels import record_map_span

EdgeMark = Literal["bar", "dash"]

# The minus sign a reverse-complemented window's list line ends with.
REVERSE_MARK = " (−)"
# Between a record's own label and its window in a list line.
LIST_LABEL_SEPARATOR = " · "
# The list text colour of the approved caption design (D-26, mock C2).
LIST_TEXT_FILL = "#6b6b6b"


@dataclass(frozen=True)
class RowWindowUnit:
    """One drawn unit of a row as the row's captions and marks read it.

    ``span`` is the part of the record the unit shows in source bases, start before end in the record's stored
    order (``start > end`` runs across the origin of a circular record); ``reverse`` is true when the unit is drawn
    reverse-complemented; ``has_region`` is false for a record drawn whole.
    """

    record_index: int
    identity_key: str
    record_id: str
    label: str
    span: tuple[int, int]
    reverse: bool
    source_length: int
    circular: bool
    has_region: bool

    @property
    def length(self) -> int:
        start, end = self.span
        return end - start + 1 if start <= end else self.source_length - start + 1 + end

    def edge_bases(self) -> tuple[int, int]:
        """The source bases at the unit's drawn left and right ends."""

        start, end = self.span
        return (end, start) if self.reverse else (start, end)

    def edge_mark(self, base: int) -> EdgeMark:
        """A bar at the first or last base of a linear record; a dashed axis anywhere else."""

        return "bar" if not self.circular and base in (1, self.source_length) else "dash"


def row_window_unit(
    record: SeqRecord,
    record_index: int,
    transform: RecordDisplayTransform,
    label: str | None,
    axis_length: int,
) -> RowWindowUnit:
    """Read one drawn unit: its record identity, its window span (or the bases it draws when whole), its
    orientation and its record's length and topology (``transform``)."""

    window = record_window_span(record)
    if window is None:
        low, high = sorted(record_map_span(record, axis_length))
        window = (max(1, low), min(transform.length, high))
    return RowWindowUnit(
        record_index=record_index,
        identity_key=str(record_identity_key(record) or f"record-{record_index + 1}"),
        record_id=str(record.id),
        label=(label or "").strip(),
        span=(int(window[0]), int(window[1])),
        reverse=transform.source_step == -1,
        source_length=int(transform.length),
        circular=bool(transform.is_circular),
        has_region=record_window_span(record) is not None,
    )


def window_span_text(start: int, end: int) -> str:
    """A window's span with thousands separators, start first (``65,000-4,000`` runs across the origin)."""

    return f"{start:,}-{end:,}"


def omitted_length_text(length_bp: int) -> str:
    """The length a "//" leaves out, as the row draws it above the gap (``1.54 Mb``, ``12.3 kb``, ``850 bp``)."""

    # Rounded half up in integers, so a length such as 1,535,000 bp reads 1.54 Mb on every platform. The unit is
    # chosen after rounding: 999,950 bp reads 1.00 Mb, not 1,000.0 kb.
    if length_bp < 1_000:
        return f"{length_bp:,} bp"
    tenths = (length_bp + 50) // 100
    if tenths < 10_000:
        return f"{tenths // 10:,}.{tenths % 10} kb"
    hundredths = (length_bp + 5_000) // 10_000
    return f"{hundredths // 100:,}.{hundredths % 100:02d} Mb"


@dataclass(frozen=True)
class RowWindowLine:
    """One line of a windowed row's list: its number, the record's own label when it differs from the row
    name (else empty), and the window text (record ID, span, and the reverse mark)."""

    number: int
    record_index: int
    label: str
    text: str


@dataclass(frozen=True)
class RowGapMark:
    """The "//" between two adjacent windows of one record; ``omitted_bp`` is the length left out, or ``None``
    when the windows overlap, touch, run backwards or in opposite directions, or another window of the record
    lies in between."""

    left_record_index: int
    right_record_index: int
    omitted_bp: int | None

    @property
    def text(self) -> str:
        return "" if self.omitted_bp is None else omitted_length_text(self.omitted_bp)


@dataclass(frozen=True)
class LinearRowWindows:
    """The captions and marks of one windowed row. ``record_indices`` are its units in drawing order (left to
    right); ``edges`` are each unit's left and right marks in that order."""

    record_indices: tuple[int, ...]
    name_record_index: int | None
    lines: tuple[RowWindowLine, ...]
    edges: tuple[tuple[EdgeMark, EdgeMark], ...]
    gaps: tuple[RowGapMark, ...]

    @property
    def leading_record_index(self) -> int:
        return self.record_indices[0]

    @property
    def heading_record_index(self) -> int:
        """The record whose identity the row heading carries: the one whose label it shows, else the leading one."""

        return self.name_record_index if self.name_record_index is not None else self.leading_record_index

    def number_of(self, record_index: int) -> int:
        return self.record_indices.index(record_index) + 1

    def edges_of(self, record_index: int) -> tuple[EdgeMark, EdgeMark]:
        return self.edges[self.record_indices.index(record_index)]

    def gap_after(self, record_index: int) -> RowGapMark | None:
        return next((gap for gap in self.gaps if gap.left_record_index == record_index), None)


def is_windowed_row(units: Sequence[RowWindowUnit]) -> bool:
    """Two or more units, at least one a window with a region."""

    return len(units) >= 2 and any(unit.has_region for unit in units)


def _arc_contains(start: int, end: int, base: int) -> bool:
    return start <= base <= end if start <= end else base >= start or base <= end


def _source_interval_left_out(earlier: RowWindowUnit, later: RowWindowUnit) -> tuple[int, int] | None:
    """The source interval between ``earlier``'s end and ``later``'s start, read in the record's stored order,
    or ``None`` when the two overlap, touch or (on a linear record) run backwards. On a circular record the
    interval runs forward from ``earlier``'s end to ``later``'s start, through the origin when needed, for every
    pair: windows of a circular record never run backwards."""

    length = earlier.source_length
    first, last = earlier.span[1] + 1, later.span[0] - 1
    if not earlier.circular:
        return (first, last) if earlier.span[1] < later.span[0] and first <= last else None
    gap = (later.span[0] - earlier.span[1] - 1) % length
    if gap <= 0 or gap + earlier.length + later.length > length:
        return None
    return ((first - 1) % length + 1, (last - 1) % length + 1)


def _omitted_bp(left: RowWindowUnit, right: RowWindowUnit, others: Sequence[RowWindowUnit]) -> int | None:
    if left.reverse != right.reverse:
        return None
    earlier, later = (right, left) if left.reverse else (left, right)
    interval = _source_interval_left_out(earlier, later)
    if interval is None:
        return None
    start, end = interval
    # The windows are next to each other in their record: no other window of it lies in the interval left out.
    if any(
        _arc_contains(start, end, other.span[0])
        or _arc_contains(start, end, other.span[1])
        or _arc_contains(other.span[0], other.span[1], start)
        for other in others
    ):
        return None
    return end - start + 1 if start <= end else left.source_length - start + 1 + end


def plan_row_windows(units: Sequence[RowWindowUnit]) -> LinearRowWindows | None:
    """The captions and marks of a row whose units are ``units`` in drawing order, or ``None`` when the row is
    not windowed (it keeps the text of each record above it)."""

    if not is_windowed_row(units):
        return None
    named = next((unit for unit in units if unit.label), None)
    name = named.label if named is not None else ""
    lines = tuple(
        RowWindowLine(
            number=number,
            record_index=unit.record_index,
            label=unit.label if unit.label and unit.label != name else "",
            text=f"{unit.record_id}:{window_span_text(*unit.span)}{REVERSE_MARK if unit.reverse else ''}",
        )
        for number, unit in enumerate(units, 1)
    )
    edges = tuple(
        (unit.edge_mark(unit.edge_bases()[0]), unit.edge_mark(unit.edge_bases()[1])) for unit in units
    )
    gaps = tuple(
        RowGapMark(
            left.record_index,
            right.record_index,
            _omitted_bp(
                left,
                right,
                [
                    unit
                    for unit in units
                    if unit.identity_key == left.identity_key and unit not in (left, right)
                ],
            ),
        )
        for left, right in zip(units, units[1:])
        if left.identity_key == right.identity_key
    )
    return LinearRowWindows(
        record_indices=tuple(unit.record_index for unit in units),
        name_record_index=named.record_index if named is not None else None,
        lines=lines,
        edges=edges,
        gaps=gaps,
    )


@dataclass(frozen=True)
class RowWindowStyle:
    """The look of a windowed row's list and marks: sizes from the feature block height and the list font,
    colours from the axis and the definition text (proportions of the approved designs D-26 mock C2 and D-27
    mock M2).

    The list reads the definition's ``accession`` line style (the line that names a record today); without its
    own size it is drawn at three quarters of the definition size, in the grey of the approved design.
    """

    block_height: float
    axis_stroke_width: float
    axis_stroke_color: str
    list_font_size: float
    list_font_weight: str
    list_fill: str
    text_fill: str
    font_family: str

    @classmethod
    def of(cls, cfg: GbdrawConfig, length_param: str, default_cds_height: float) -> RowWindowStyle:
        definition = cfg.objects.definition.linear
        accession = definition.line_styles.accession
        definition_size = definition.font_size.for_length_param(length_param)
        return cls(
            # A feature block on one strand: the height the marks scale with.
            block_height=0.5 * float(default_cds_height),
            axis_stroke_width=float(cfg.objects.axis.linear.stroke_width.for_length_param(length_param)),
            axis_stroke_color=str(cfg.objects.axis.linear.stroke_color),
            list_font_size=float(accession.font_size) if accession.font_size is not None else 0.75 * definition_size,
            list_font_weight=str(accession.font_weight or definition.font_weight),
            list_fill=str(accession.fill) if accession.fill is not None else LIST_TEXT_FILL,
            text_fill=str(definition.fill),
            font_family=str(cfg.objects.text.font_family),
        )

    @property
    def bar_half(self) -> float:
        return 0.4 * self.block_height

    @property
    def bar_stroke(self) -> float:
        return max(1.0, 0.09 * self.block_height)

    @property
    def bar_offset(self) -> float:
        """How far outside the window edge the bar's centre stands."""
        return 0.5 * self.bar_stroke + 0.5

    @property
    def dash_length(self) -> float:
        return 0.5 * self.block_height

    @property
    def dash_pattern(self) -> tuple[float, float]:
        return (0.12 * self.block_height, 0.085 * self.block_height)

    @property
    def slash_half(self) -> float:
        return 0.28 * self.block_height

    @property
    def slash_slant(self) -> float:
        return 0.093 * self.block_height

    @property
    def slash_spacing(self) -> float:
        """Half the distance between the two strokes of "//"."""
        return 0.128 * self.block_height

    @property
    def slash_stroke(self) -> float:
        """"//" breaks the axis, so it is drawn with the axis's stroke (its colour is ``axis_stroke_color``)."""
        return self.axis_stroke_width

    @property
    def slash_width(self) -> float:
        return 2.0 * (self.slash_spacing + self.slash_slant) + self.slash_stroke

    @property
    def circle_radius(self) -> float:
        return 0.45 * self.list_font_size

    @property
    def circle_stroke(self) -> float:
        return max(0.75, 0.18 * self.circle_radius)

    @property
    def circle_gap(self) -> float:
        return 0.5 * self.circle_radius

    @property
    def number_font_size(self) -> float:
        return 1.45 * self.circle_radius

    @property
    def number_room(self) -> float:
        """The width a window's number takes at its upper left, with its gap to the window."""
        return 2.0 * self.circle_radius + self.circle_gap

    @property
    def pad(self) -> float:
        return 0.1 * self.block_height

    def edge_room(self, mark: EdgeMark) -> float:
        return self.dash_length if mark == "dash" else self.bar_stroke + 0.5

    def omitted_text_box(self, text: str) -> tuple[float, float]:
        return drawn_text_box(text, self.font_family, self.list_font_size, "bold") if text else (0.0, 0.0)

    def marks_band(self, feature_top_y: float) -> VerticalBand:
        """The vertical extent of a window's marks, number and gap text, from the top of its feature lanes
        (axis coordinates): numbers and gap texts are centred on that top."""

        text_half = 0.5 * drawn_text_box("0", self.font_family, self.list_font_size, "bold")[1]
        top = float(feature_top_y) - max(self.circle_radius + self.circle_stroke, text_half + 1.0)
        bottom = max(self.bar_half, self.slash_half + 0.5 * self.slash_stroke)
        return VerticalBand(min(top, -bottom), bottom)


def row_window_insets(plan: LinearRowWindows, style: RowWindowStyle) -> dict[int, tuple[float, float]]:
    """The room each unit of a windowed row keeps on its left and right for its marks, by record index: the gap
    between two units is at least the sum of the two sides that face it, and a gap with "//" fits the two edge
    marks, the "//" and the omitted length clear of the next window's number."""

    left_room = {
        index: max(style.edge_room(edges[0]), style.number_room)
        for index, edges in zip(plan.record_indices, plan.edges)
    }
    right_room = {index: style.edge_room(edges[1]) for index, edges in zip(plan.record_indices, plan.edges)}
    for gap in plan.gaps:
        marks = (
            right_room[gap.left_record_index]
            + style.slash_width
            + left_room[gap.right_record_index]
            + 2.0 * style.pad
        )
        # The length's backing box keeps a circle gap from the next window's number on each side.
        text = style.omitted_text_box(gap.text)[0] + 2.0 * (style.number_room + style.circle_gap + style.pad)
        right_room[gap.left_record_index] = max(marks, text) - left_room[gap.right_record_index]
    return {index: (left_room[index], right_room[index]) for index in plan.record_indices}


__all__ = [
    "EdgeMark",
    "LIST_LABEL_SEPARATOR",
    "LIST_TEXT_FILL",
    "LinearRowWindows",
    "REVERSE_MARK",
    "RowGapMark",
    "RowWindowLine",
    "RowWindowStyle",
    "RowWindowUnit",
    "is_windowed_row",
    "omitted_length_text",
    "plan_row_windows",
    "row_window_insets",
    "row_window_unit",
    "window_span_text",
]
