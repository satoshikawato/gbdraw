"""Pure placement contracts for Linear diagrams with multiple records per row."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import NewType, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError


RecordKey = NewType("RecordKey", str)


@dataclass(frozen=True)
class LinearRecordMeasurement:
    """Intrinsic record extents used by the Linear layout solver."""

    record_index: int
    record_key: RecordKey
    sequence_length: int
    left_inset: float = 0.0
    right_inset: float = 0.0
    top_extent: float = 0.0
    bottom_extent: float = 0.0

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

        return self.x + (float(position) * self.px_per_bp)


@dataclass(frozen=True)
class LinearLayoutPlan:
    """Immutable result of resolving Linear rows and record placements."""

    placements: tuple[LinearRecordPlacement, ...]
    px_per_bp: float
    row_count: int
    content_left: float
    content_right: float
    content_top: float
    content_bottom: float

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


def stable_record_keys(records: Sequence[SeqRecord]) -> tuple[RecordKey, ...]:
    """Return deterministic opaque keys without relying on unique record IDs."""

    keys = tuple(
        RecordKey(
            str(record.annotations.get("gbdraw_record_key") or f"record-{index + 1}")
        )
        for index, record in enumerate(records)
    )
    if len(set(keys)) != len(keys):
        raise ValidationError("Linear record keys must be unique.")
    return keys


def parse_record_row_position(value: str) -> tuple[str, int]:
    """Parse the shared ``<selector>@<row>`` placement syntax."""

    raw = str(value).strip()
    if not raw:
        raise ValidationError("multi_record_position does not allow empty entries.")
    if raw.count("@") != 1:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must be in '<selector>@<row>' format."
        )
    selector, row_text = (part.strip() for part in raw.rsplit("@", 1))
    if not selector:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must include a selector before '@'."
        )
    try:
        row = int(row_text)
    except ValueError as exc:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        ) from exc
    if row <= 0:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        )
    return selector, row


def resolve_record_row_positions(
    records: Sequence[SeqRecord],
    positions: Sequence[str] | None,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Resolve selectors to input indices and contiguous zero-based rows.

    The first returned tuple is the record order after sorting by row and token
    order.  The second tuple maps each input record index to its normalized row.
    """

    record_count = len(records)
    if not positions:
        indices = tuple(range(record_count))
        return indices, indices
    if len(positions) != record_count:
        raise ValidationError(
            f"multi_record_position must provide exactly {record_count} entry(ies)."
        )

    resolved: list[tuple[int, int, int]] = []
    seen: set[int] = set()
    ids: dict[str, list[int]] = {}
    for index, record in enumerate(records):
        ids.setdefault(str(record.id), []).append(index)

    for token_index, raw in enumerate(positions):
        selector, row = parse_record_row_position(str(raw))
        if selector.startswith("#"):
            try:
                record_index = int(selector[1:]) - 1
            except ValueError as exc:
                raise ValidationError(
                    f"multi_record_position selector '{selector}' is invalid."
                ) from exc
            if record_index < 0 or record_index >= record_count:
                raise ValidationError(
                    f"multi_record_position selector '{selector}' is out of range for "
                    f"{record_count} loaded record(s)."
                )
        else:
            matches = ids.get(selector, [])
            if not matches:
                raise ValidationError(
                    f"multi_record_position selector '{selector}' did not match any record ID."
                )
            if len(matches) > 1:
                raise ValidationError(
                    f"multi_record_position selector '{selector}' matched multiple records; "
                    "use a #index selector."
                )
            record_index = matches[0]
        if record_index in seen:
            raise ValidationError(
                f"multi_record_position selector '{selector}' was specified more than once."
            )
        seen.add(record_index)
        resolved.append((row, token_index, record_index))

    if len(seen) != record_count:
        raise ValidationError(
            "multi_record_position must include each loaded record exactly once."
        )

    row_map = {value: index for index, value in enumerate(sorted({item[0] for item in resolved}))}
    ordered = sorted(resolved, key=lambda item: (row_map[item[0]], item[1]))
    rows_by_index = [0] * record_count
    for input_row, _token_index, record_index in resolved:
        rows_by_index[record_index] = row_map[input_row]
    return tuple(item[2] for item in ordered), tuple(rows_by_index)


def solve_linear_layout(
    measurements: Sequence[LinearRecordMeasurement],
    rows_by_record: Sequence[int],
    *,
    available_width: float,
    record_gap_px: float = 24.0,
    align_center: bool = True,
    first_axis_y: float = 0.0,
    row_gap_px: float = 0.0,
    comparison_height: float = 0.0,
    record_order: Sequence[int] | None = None,
) -> LinearLayoutPlan:
    """Resolve one shared bp scale and record placements without creating SVG."""

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

    row_scales: list[float] = []
    for row_items in grouped:
        gap_width = gap * max(0, len(row_items) - 1)
        if width <= gap_width:
            raise ValidationError(
                "Linear row cannot fit a positive bp scale: "
                f"record_count={len(row_items)}, record_gap_px={gap:g}, "
                f"available_width={width:g}."
            )
        row_scales.append(
            (width - gap_width) / sum(item.sequence_length for item in row_items)
        )
    px_per_bp = min(row_scales)
    if not math.isfinite(px_per_bp) or px_per_bp <= 0:
        raise ValidationError("Linear layout did not produce a positive shared bp scale.")

    row_axis_y: list[float] = [float(first_axis_y)]
    for row_index in range(1, len(grouped)):
        previous_bottom = max(item.bottom_extent for item in grouped[row_index - 1])
        current_top = max(item.top_extent for item in grouped[row_index])
        row_axis_y.append(
            row_axis_y[-1]
            + previous_bottom
            + max(float(row_gap_px), float(comparison_height))
            + current_top
        )

    placements: list[LinearRecordPlacement] = []
    content_left = 0.0
    content_right = 0.0
    for row_index, row_items in enumerate(grouped):
        sequence_widths = [item.sequence_length * px_per_bp for item in row_items]
        effective_gaps = [
            max(gap, row_items[index].right_inset + row_items[index + 1].left_inset)
            for index in range(len(row_items) - 1)
        ]
        row_width = sum(sequence_widths) + sum(effective_gaps)
        cursor = 0.5 * (width - row_width) if align_center else 0.0
        content_left = min(content_left, cursor - row_items[0].left_inset)
        column = 0
        for item, sequence_width_value in zip(row_items, sequence_widths):
            placement = LinearRecordPlacement(
                record_index=item.record_index,
                record_key=item.record_key,
                row=row_index,
                column=column,
                x=cursor,
                axis_y=row_axis_y[row_index],
                sequence_width=sequence_width_value,
                left_inset=item.left_inset,
                right_inset=item.right_inset,
                top_extent=item.top_extent,
                bottom_extent=item.bottom_extent,
                comparison_top_y=row_axis_y[row_index] - item.top_extent,
                comparison_bottom_y=row_axis_y[row_index] + item.bottom_extent,
                px_per_bp=px_per_bp,
            )
            placements.append(placement)
            content_right = max(
                content_right,
                cursor + sequence_width_value + item.right_inset,
            )
            if column < len(effective_gaps):
                cursor += sequence_width_value + effective_gaps[column]
            column += 1

    placements.sort(key=lambda item: item.record_index)
    content_top = min(item.comparison_top_y for item in placements)
    content_bottom = max(item.comparison_bottom_y for item in placements)
    return LinearLayoutPlan(
        placements=tuple(placements),
        px_per_bp=px_per_bp,
        row_count=len(grouped),
        content_left=content_left,
        content_right=content_right,
        content_top=content_top,
        content_bottom=content_bottom,
    )


__all__ = [
    "LinearLayoutPlan",
    "LinearRecordMeasurement",
    "LinearRecordPlacement",
    "RecordKey",
    "parse_record_row_position",
    "resolve_record_row_positions",
    "solve_linear_layout",
    "stable_record_keys",
]
