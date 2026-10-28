"""Pure placement contracts for Linear diagrams with multiple records per row."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Mapping, NewType, Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.core.numeric import compensated_sum
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear import AxisGapResolution, CollisionBand, resolve_axis_gap
from gbdraw.layout.record_placement import (
    parse_record_row_position,
    resolve_record_row_positions,
)


RecordKey = NewType("RecordKey", str)


@dataclass(frozen=True)
class LinearRecordMeasurement:
    """Canvas enclosure and local collision bands used by the layout solver."""

    record_index: int
    record_key: RecordKey
    sequence_length: int
    left_inset: float = 0.0
    right_inset: float = 0.0
    top_extent: float = 0.0
    bottom_extent: float = 0.0
    comparison_top_extent: float | None = None
    comparison_bottom_extent: float | None = None
    collision_bands: tuple[CollisionBand, ...] = ()

    def __post_init__(self) -> None:
        if self.record_index < 0:
            raise ValidationError("record_index must be non-negative.")
        if not str(self.record_key).strip():
            raise ValidationError("record_key must be a non-empty string.")
        if self.sequence_length <= 0:
            raise ValidationError("sequence_length must be positive.")
        for name in ("left_inset", "right_inset", "top_extent", "bottom_extent"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0:
                raise ValidationError(f"{name} must be a finite non-negative number.")
        for name, layout_extent_name in (
            ("comparison_top_extent", "top_extent"),
            ("comparison_bottom_extent", "bottom_extent"),
        ):
            raw_value = getattr(self, name)
            if raw_value is None:
                continue
            value = float(raw_value)
            if not math.isfinite(value) or value < 0:
                raise ValidationError(f"{name} must be a finite non-negative number.")
            if value > float(getattr(self, layout_extent_name)):
                raise ValidationError(f"{name} cannot exceed {layout_extent_name}.")
        object.__setattr__(self, "collision_bands", tuple(self.collision_bands))


def _x_for_position(x: float, px_per_bp: float, position: float) -> float:
    return x + (float(position) * px_per_bp)


# A comparison's ends and attachment offsets are kept to this many decimals of a px: far below the 3 that path data
# writes (``PATH_DECIMALS``), far above the float noise (about 1e-13 px) that two equal ribbon heights reached by
# different sums of canvas ys carry. Rows that move without changing a ribbon then leave its key equal (OV-2089).
COMPARISON_FRAME_DECIMALS = 9


def comparison_frame_px(value: float) -> float:
    """``value`` on the comparison frame grid (``COMPARISON_FRAME_DECIMALS``)."""

    return round(float(value), COMPARISON_FRAME_DECIMALS) + 0.0


@dataclass(frozen=True)
class LinearFeatureAttachment:
    """How far a drawn feature's attachment band lies from its record's comparison anchors, in px (positive is
    down): a ribbon end on that feature moves by this much from the anchor it would meet."""

    top_offset_px: float
    bottom_offset_px: float


@dataclass(frozen=True)
class LinearComparisonEndpoint:
    """What a comparison reads of one record's placement (the ``comparison_content`` key).

    The record, its row, its X scale, the y of the anchor that faces the other
    record in the comparison's own frame (0 at the comparison's top), and the
    offset of each attached feature from that anchor. Rows that move together
    leave all of these unchanged (OV-2089); the record's own extents and axis
    do not reach the ribbons.
    """

    record_index: int
    row: int
    x: float
    px_per_bp: float
    anchor_y: float
    attachment_offsets_px: Mapping[str, float]

    def x_for_position(self, position: float) -> float:
        """Return the canvas-local X coordinate for a display-local bp position."""

        return _x_for_position(self.x, self.px_per_bp, position)


@dataclass(frozen=True)
class LinearRecordPlacement:
    """Resolved geometry for one displayed Linear record."""

    record_index: int
    record_key: RecordKey
    row: int
    column: int
    x: float
    axis_y: float
    sequence_width: float
    left_inset: float
    right_inset: float
    top_extent: float
    bottom_extent: float
    comparison_top_y: float
    comparison_bottom_y: float
    px_per_bp: float

    def x_for_position(self, position: float) -> float:
        """Return the canvas-local X coordinate for a display-local bp position."""

        return _x_for_position(self.x, self.px_per_bp, position)

    def faces_bottom(self, other: LinearRecordPlacement) -> bool:
        """Whether a comparison with ``other`` meets this record at its bottom anchor (it lies on a higher row)."""

        return self.row < other.row

    def facing_comparison_anchor_y(self, other: LinearRecordPlacement) -> float:
        """The canvas y where a comparison with ``other`` meets this record."""

        return self.comparison_bottom_y if self.faces_bottom(other) else self.comparison_top_y

    def comparison_endpoint(
        self,
        other: LinearRecordPlacement,
        *,
        anchor_y: float,
        attachments: Mapping[str, LinearFeatureAttachment],
    ) -> LinearComparisonEndpoint:
        """This record's end of a comparison with ``other``, at ``anchor_y`` in the comparison's frame."""

        faces_bottom = self.faces_bottom(other)
        return LinearComparisonEndpoint(
            self.record_index,
            self.row,
            self.x,
            self.px_per_bp,
            anchor_y,
            {
                feature_id: attachment.bottom_offset_px if faces_bottom else attachment.top_offset_px
                for feature_id, attachment in attachments.items()
            },
        )


@dataclass(frozen=True)
class LinearComparisonFrame:
    """Where one comparison draws: the canvas y of its top, and the y of its query and subject ends below it."""

    top_y: float
    query_y: float
    subject_y: float


def linear_comparison_frame(query: LinearRecordPlacement, subject: LinearRecordPlacement) -> LinearComparisonFrame:
    """The frame of a comparison between ``query`` and ``subject``: the one owner of where its ribbons end.

    The top is the higher facing anchor on the canvas, where the group is
    placed; the ends lie below it on the ``comparison_frame_px`` grid.
    """

    query_anchor = query.facing_comparison_anchor_y(subject)
    subject_anchor = subject.facing_comparison_anchor_y(query)
    top_y = min(query_anchor, subject_anchor)
    return LinearComparisonFrame(
        top_y=top_y,
        query_y=comparison_frame_px(query_anchor - top_y),
        subject_y=comparison_frame_px(subject_anchor - top_y),
    )


@dataclass(frozen=True)
class LinearLayoutPlan:
    """Immutable result of resolving Linear rows and record placements."""

    placements: tuple[LinearRecordPlacement, ...]
    row_count: int
    content_left: float
    content_right: float
    content_top: float
    content_bottom: float
    row_gap_resolutions: tuple[AxisGapResolution, ...] = ()

    def placement_for_index(self, record_index: int) -> LinearRecordPlacement:
        for placement in self.placements:
            if placement.record_index == record_index:
                return placement
        raise KeyError(record_index)

    def placement_for_key(self, record_key: str) -> LinearRecordPlacement:
        for placement in self.placements:
            if str(placement.record_key) == str(record_key):
                return placement
        raise KeyError(record_key)


def record_key_at(record: SeqRecord, index: int) -> RecordKey:
    """The key a request gave ``record``, or its position (``record-<n>``)."""

    return RecordKey(str(record.annotations.get("gbdraw_record_key") or f"record-{index + 1}"))


def stable_record_keys(records: Sequence[SeqRecord]) -> tuple[RecordKey, ...]:
    """Return deterministic opaque keys without relying on unique record IDs."""

    keys = tuple(record_key_at(record, index) for index, record in enumerate(records))
    if len(set(keys)) != len(keys):
        raise ValidationError("Linear record keys must be unique.")
    return keys


def record_pairs_between_adjacent_rows(
    rows_by_record: Sequence[int],
) -> tuple[tuple[int, int], ...]:
    """Return every cross-row pair whose normalized rows are adjacent."""

    if any(not isinstance(row, int) or isinstance(row, bool) or row < 0 for row in rows_by_record):
        raise ValidationError("Linear rows must be non-negative integers.")
    return tuple(
        (query_index, subject_index)
        for query_index in range(len(rows_by_record))
        for subject_index in range(query_index + 1, len(rows_by_record))
        if abs(rows_by_record[query_index] - rows_by_record[subject_index]) == 1
    )


@dataclass(frozen=True)
class LinearHorizontalLayout:
    """Where each record sits along its row: x, width and the shared bp scale.

    Its placements carry zero vertical fields; ``solve_linear_layout`` places them vertically.
    The labels, the definition keep-clear boxes and the collision bands read these places before
    the rows are stacked, so the drawing and the vertical solve read the same horizontal frame.
    """

    placements: tuple[LinearRecordPlacement, ...]
    rows: tuple[tuple[int, ...], ...]
    content_left: float
    content_right: float

    def placement_for_index(self, record_index: int) -> LinearRecordPlacement:
        for placement in self.placements:
            if placement.record_index == record_index:
                return placement
        raise KeyError(record_index)


def place_record_vertically(
    horizontal: LinearRecordPlacement,
    *,
    axis_y: float,
    top_extent: float,
    bottom_extent: float,
    comparison_top_extent: float,
    comparison_bottom_extent: float,
    comparison_endpoint_gap: float,
) -> LinearRecordPlacement:
    """A record's final placement: its horizontal place on the row whose axis is ``axis_y``."""

    return replace(
        horizontal,
        axis_y=axis_y,
        top_extent=top_extent,
        bottom_extent=bottom_extent,
        comparison_top_y=axis_y - comparison_top_extent - comparison_endpoint_gap,
        comparison_bottom_y=axis_y + comparison_bottom_extent + comparison_endpoint_gap,
    )


def solve_linear_horizontal_layout(
    measurements: Sequence[LinearRecordMeasurement],
    rows_by_record: Sequence[int],
    *,
    available_width: float,
    record_gap_px: float = 24.0,
    align_center: bool = True,
    normalize_length: bool = False,
    record_order: Sequence[int] | None = None,
) -> LinearHorizontalLayout:
    """Resolve each record's x, width and bp scale in its row.

    It is the one owner of a record's x and width in both row modes. A record's width is its
    share of the narrowest row's scale (``scale_px * length / scale_bp``; a record alone in its
    row: its share of the longest record), and ``align_center`` splits a row's free space
    evenly on both sides. With ``normalize_length`` every record fills its own row.
    A record's ``px_per_bp`` is its width over its length.

    It reads only the records' lengths and horizontal insets, never their vertical extents or
    collision bands, so it runs once per drawing, before the labels are placed.
    """

    items = tuple(measurements)
    if not items:
        raise ValidationError("Linear layout requires at least one record.")
    if len(rows_by_record) != len(items):
        raise ValidationError("rows_by_record must contain one row for every record.")
    width = float(available_width)
    gap = float(record_gap_px)
    if not math.isfinite(width) or width <= 0:
        raise ValidationError("available_width must be a finite positive number.")
    if not math.isfinite(gap) or gap < 0:
        raise ValidationError("record_gap_px must be a finite non-negative number.")
    if any(not isinstance(row, int) or isinstance(row, bool) or row < 0 for row in rows_by_record):
        raise ValidationError("Linear rows must be non-negative integers.")
    if len({item.record_index for item in items}) != len(items):
        raise ValidationError("Linear layout record indices must be unique.")
    if {item.record_index for item in items} != set(range(len(items))):
        raise ValidationError("Linear layout record indices must be contiguous from zero.")

    normalized_row_values = sorted(set(rows_by_record))
    row_normalization = {row: index for index, row in enumerate(normalized_row_values)}
    normalized_rows = tuple(row_normalization[row] for row in rows_by_record)
    grouped: list[list[LinearRecordMeasurement]] = [[] for _ in normalized_row_values]
    by_index = {item.record_index: item for item in items}
    ordered_indices = tuple(record_order) if record_order is not None else tuple(range(len(items)))
    if set(ordered_indices) != set(range(len(items))) or len(ordered_indices) != len(items):
        raise ValidationError("record_order must contain every record index exactly once.")
    for record_index in ordered_indices:
        row = normalized_rows[record_index]
        grouped[row].append(by_index[record_index])

    if normalize_length and any(len(row_items) > 1 for row_items in grouped):
        raise ValidationError(
            "normalize_length=True cannot be combined with multiple records in one Linear row."
        )
    row_px: list[float] = []
    row_bp: list[int] = []
    for row_items in grouped:
        gap_width = gap * max(0, len(row_items) - 1)
        if width <= gap_width:
            raise ValidationError(
                "Linear row cannot fit a positive bp scale: "
                f"record_count={len(row_items)}, record_gap_px={gap:g}, "
                f"available_width={width:g}."
            )
        row_px.append(width - gap_width)
        row_bp.append(sum(item.sequence_length for item in row_items))
    row_scales = [px / bp for px, bp in zip(row_px, row_bp)]
    narrowest_scale = min(row_scales)
    if not math.isfinite(narrowest_scale) or narrowest_scale <= 0:
        raise ValidationError("Linear layout did not produce a positive shared bp scale.")
    # The narrowest row's pixels and bp: every row draws at its scale.
    narrowest = row_scales.index(narrowest_scale)
    scale_px, scale_bp = row_px[narrowest], row_bp[narrowest]

    placements: list[LinearRecordPlacement] = []
    content_left = 0.0
    content_right = 0.0
    for row_index, row_items in enumerate(grouped):
        effective_gaps = [
            max(gap, row_items[index].right_inset + row_items[index + 1].left_inset)
            for index in range(len(row_items) - 1)
        ]
        if normalize_length:
            sequence_widths = [width]
            cursor = 0.0
        else:
            sequence_widths = [scale_px * item.sequence_length / scale_bp for item in row_items]
            # The row's free space: the scale's pixels its bp leave, plus those its gaps leave.
            free_px = (
                scale_px * ((scale_bp - row_bp[row_index]) / scale_bp)
                + (width - scale_px)
                - compensated_sum(effective_gaps)
            )
            cursor = 0.5 * free_px if align_center else 0.0
        content_left = min(content_left, cursor - row_items[0].left_inset)
        for column, (item, sequence_width_value) in enumerate(zip(row_items, sequence_widths)):
            placements.append(
                LinearRecordPlacement(
                    record_index=item.record_index,
                    record_key=item.record_key,
                    row=row_index,
                    column=column,
                    x=cursor,
                    axis_y=0.0,
                    sequence_width=sequence_width_value,
                    left_inset=item.left_inset,
                    right_inset=item.right_inset,
                    top_extent=0.0,
                    bottom_extent=0.0,
                    comparison_top_y=0.0,
                    comparison_bottom_y=0.0,
                    px_per_bp=sequence_width_value / item.sequence_length,
                )
            )
            content_right = max(
                content_right,
                cursor + sequence_width_value + item.right_inset,
            )
            if column < len(effective_gaps):
                cursor += sequence_width_value + effective_gaps[column]

    return LinearHorizontalLayout(
        placements=tuple(placements),
        rows=tuple(tuple(item.record_index for item in row_items) for row_items in grouped),
        content_left=content_left,
        content_right=content_right,
    )


def solve_linear_layout(
    measurements: Sequence[LinearRecordMeasurement],
    horizontal: LinearHorizontalLayout,
    *,
    first_axis_y: float = 0.0,
    row_gap_px: float = 0.0,
    comparison_height: float = 0.0,
    comparison_endpoint_gap_px: float = 0.0,
    definition_clear_gap: float = 0.0,
    comparison_record_indices_by_boundary: Mapping[int, Sequence[int]] | None = None,
) -> LinearLayoutPlan:
    """Stack the rows of ``horizontal`` from the records' vertical extents and collision bands.

    The records keep the horizontal places ``horizontal`` gave them; their bands are local to
    their x, which this adds back.
    """

    items = tuple(measurements)
    # One measurement per placed record, the same record and insets: a caller error, not input.
    if sorted((item.record_index, item.record_key, item.left_inset, item.right_inset) for item in items) != sorted(
        (place.record_index, place.record_key, place.left_inset, place.right_inset)
        for place in horizontal.placements
    ):
        raise ValueError("Linear layout measurements must match the records' horizontal places.")
    by_index = {item.record_index: item for item in items}
    ordinary_gap = float(row_gap_px)
    comparison_gap = float(comparison_height)
    comparison_endpoint_gap = float(comparison_endpoint_gap_px)
    definition_gap = float(definition_clear_gap)
    for name, value in (
        ("row_gap_px", ordinary_gap),
        ("comparison_height", comparison_gap),
        ("comparison_endpoint_gap_px", comparison_endpoint_gap),
        ("definition_clear_gap", definition_gap),
    ):
        if not math.isfinite(value) or value < 0:
            raise ValidationError(f"{name} must be a finite non-negative number.")
    comparison_clearance = comparison_gap + (2.0 * comparison_endpoint_gap)
    rows = horizontal.rows

    boundary_records: dict[int, frozenset[int]] = {}
    if comparison_record_indices_by_boundary is None:
        if comparison_clearance > 0.0:
            for boundary in range(max(0, len(rows) - 1)):
                boundary_records[boundary] = frozenset((*rows[boundary], *rows[boundary + 1]))
    else:
        for raw_boundary, raw_indices in comparison_record_indices_by_boundary.items():
            boundary = int(raw_boundary)
            if boundary < 0 or boundary >= len(rows) - 1:
                raise ValidationError(
                    "comparison boundary indices must identify adjacent resolved rows."
                )
            indices = frozenset(int(index) for index in raw_indices)
            if not indices <= {*rows[boundary], *rows[boundary + 1]}:
                raise ValidationError(
                    "comparison boundary record indices must belong to its adjacent rows."
                )
            if indices:
                boundary_records[boundary] = indices

    def comparison_extents(item: LinearRecordMeasurement) -> tuple[float, float]:
        return (
            item.top_extent if item.comparison_top_extent is None else item.comparison_top_extent,
            item.bottom_extent
            if item.comparison_bottom_extent is None
            else item.comparison_bottom_extent,
        )

    collision_bands_by_record: dict[int, tuple[CollisionBand, ...]] = {}
    for placement in horizontal.placements:
        item = by_index[placement.record_index]
        local_bands = item.collision_bands
        if not local_bands:
            comparison_top_extent, comparison_bottom_extent = comparison_extents(item)
            local_bands = (
                CollisionBand(
                    "body",
                    0.0,
                    placement.sequence_width,
                    -float(item.top_extent),
                    float(item.bottom_extent),
                ),
                CollisionBand(
                    "comparison",
                    0.0,
                    placement.sequence_width,
                    -float(comparison_top_extent),
                    float(comparison_bottom_extent),
                ),
            )
        collision_bands_by_record[item.record_index] = tuple(
            band.translate(x=placement.x) for band in local_bands
        )

    row_gap_resolutions: list[AxisGapResolution] = []
    row_axis_y: list[float] = [float(first_axis_y)]
    for boundary in range(max(0, len(rows) - 1)):
        active_records = boundary_records.get(boundary, frozenset())

        def boundary_bands(row_index: int) -> tuple[CollisionBand, ...]:
            return tuple(
                band
                for record_index in rows[row_index]
                for band in collision_bands_by_record[record_index]
                if band.kind != "comparison" or record_index in active_records
            )

        resolution = resolve_axis_gap(
            boundary_bands(boundary),
            boundary_bands(boundary + 1),
            ordinary_row_gap=ordinary_gap,
            comparison_height=comparison_clearance,
            definition_clear_gap=definition_gap,
            boundary_has_comparison=bool(active_records),
        )
        row_gap_resolutions.append(resolution)
        row_axis_y.append(row_axis_y[-1] + resolution.axis_gap)

    placements: list[LinearRecordPlacement] = []
    for placement in horizontal.placements:
        item = by_index[placement.record_index]
        comparison_top_extent, comparison_bottom_extent = comparison_extents(item)
        placements.append(
            place_record_vertically(
                placement,
                axis_y=row_axis_y[placement.row],
                top_extent=item.top_extent,
                bottom_extent=item.bottom_extent,
                comparison_top_extent=comparison_top_extent,
                comparison_bottom_extent=comparison_bottom_extent,
                comparison_endpoint_gap=comparison_endpoint_gap,
            )
        )

    placements.sort(key=lambda item: item.record_index)
    content_top = min(item.axis_y - item.top_extent for item in placements)
    content_bottom = max(item.axis_y + item.bottom_extent for item in placements)
    return LinearLayoutPlan(
        placements=tuple(placements),
        row_count=len(rows),
        content_left=horizontal.content_left,
        content_right=horizontal.content_right,
        content_top=content_top,
        content_bottom=content_bottom,
        row_gap_resolutions=tuple(row_gap_resolutions),
    )


__all__ = [
    "LinearComparisonEndpoint",
    "LinearComparisonFrame",
    "LinearFeatureAttachment",
    "LinearHorizontalLayout",
    "LinearLayoutPlan",
    "LinearRecordMeasurement",
    "LinearRecordPlacement",
    "RecordKey",
    "parse_record_row_position",
    "COMPARISON_FRAME_DECIMALS",
    "comparison_frame_px",
    "linear_comparison_frame",
    "place_record_vertically",
    "resolve_record_row_positions",
    "solve_linear_horizontal_layout",
    "solve_linear_layout",
    "record_key_at",
    "stable_record_keys",
]
