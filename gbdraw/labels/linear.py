#!/usr/bin/env python
# coding: utf-8

"""Linear track-based label placement utilities."""

import math
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, NamedTuple, cast
from .coordinates import display_label_segment

from .filtering import get_label_text, label_filter_maps
from .label_style import ROW_LABEL_GAP_EM
from ..config.models import LinearRenderProfile
from ..features.coordinates import get_strand
from ..features.ids import compute_feature_object_hash
from ..core.record_metadata import _source_feature_index
from ..core.text import calculate_bbox_dimensions, calculate_svg_bbox_dimensions
from ..auto_sizes import determine_length_parameter
from ..layout.linear_coords import normalize_position_to_linear_track
from ..layout.linear import (
    LinearFeatureLaneGeometry,
    calculate_feature_position_factors_linear,
)
from ..layout.spatial import Aabb, Interval, IntervalIndex
from ..layout.text_geometry import (
    GEOMETRY_EPSILON,
    Point,
    aabb_from_points,
    anchor_x_bounds,
    convex_polygons_intersect,
    text_box_corner_offsets,
    translate_points,
)
from .fit import (
    FITS,
    SEPARATE_FEATURES_VALUES,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitOutcome,
    LabelFitVerification,
    LabelFitValue,
    LabelShape,
    crossing_pairs,
    label_fit_outcome,
    overlapping_pairs,
    raise_or_warn,
    switch_row,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

def _row_label_margin_px(label: Mapping[str, Any]) -> float:
    """The room a label keeps free at both ends of its text in its label row: half of
    ``ROW_LABEL_GAP_EM`` of its own font size, so neighbours on one line keep the gap between them
    (LP-D104), horizontal or turned. The rows are the auto rows and the rows ``above_feature``
    labels stack in; a label drawn inside its feature is in no row."""

    if label.get("is_embedded") and "above_feature_place_above" not in label:
        return 0.0
    return 0.5 * ROW_LABEL_GAP_EM * float(label["font_size"])


def check_label_overlap(label1, label2):
    """Check if two labels overlap horizontally, with the room each keeps beside it in its row"""
    margin1, margin2 = _row_label_margin_px(label1), _row_label_margin_px(label2)
    return not (label1["end"] + margin1 < label2["start"] - margin2 or label2["end"] + margin2 < label1["start"] - margin1)


def find_lowest_available_track(track_dict, label):
    """Find the lowest track (closest to track_1) where the label can be placed without overlap"""
    track_num = 1
    while True:
        track_id = f"track_{track_num}"
        # Check if this track has overlaps
        has_overlap = False
        if track_id in track_dict:
            for existing_label in track_dict[track_id]:
                if check_label_overlap(label, existing_label):
                    has_overlap = True
                    break

        if not has_overlap:
            return track_num
        track_num += 1


def _linear_label_interval(label: dict) -> Interval:
    margin = _row_label_margin_px(label)
    return Interval(float(label["start"]) - margin, float(label["end"]) + margin)


def _external_label_bucket_size(alignment_width: float) -> float:
    return max(16.0, float(alignment_width) / 256.0)


def _above_feature_label_bucket_size(labels: list[dict], min_gap_px: float) -> float:
    """Bucket width of the above-feature broad phase on x: a few label heights.

    A bucket sized by the widest label puts most of a dense record in a few
    buckets, which made the collision pass quadratic (OV-695).
    """
    max_height = max((float(label.get("height_px", 0.0)) for label in labels), default=0.0)
    return max(16.0, 4.0 * max_height + float(min_gap_px))


def _find_lowest_available_track_indexed(
    track_dict: dict,
    track_indexes: dict[str, IntervalIndex],
    label_by_id: dict[int, dict],
    label: dict,
    blocked: Callable[[int], bool] | None = None,
) -> int:
    """Find the lowest non-overlapping external label track using index candidates, skipping the
    tracks ``blocked`` refuses (a keep-clear box in the label's place there)."""
    label_interval = _linear_label_interval(label)
    track_num = 1
    while True:
        if blocked is not None and blocked(track_num):
            track_num += 1
            continue
        track_id = f"track_{track_num}"
        if track_id not in track_dict or not track_dict[track_id]:
            return track_num

        index = track_indexes.get(track_id)
        if index is None:
            candidates = track_dict[track_id]
        else:
            # Only int ids are inserted into these indexes.
            candidates = [label_by_id[cast(int, label_id)] for label_id in index.query(label_interval)]

        if not any(check_label_overlap(label, existing_label) for existing_label in candidates):
            return track_num
        track_num += 1


def _insert_external_label_index(
    track_indexes: dict[str, IntervalIndex],
    track_id: str,
    label_id: int,
    label: dict,
    bucket_size: float,
) -> None:
    index = track_indexes.setdefault(track_id, IntervalIndex(bucket_size=bucket_size))
    index.insert(label_id, _linear_label_interval(label))


@lru_cache(maxsize=4096)
def _rotated_corner_offsets(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
    along_text_px: float = 0.0,
) -> tuple[tuple[float, float], ...]:
    return text_box_corner_offsets(width_px, height_px, rotation_deg, text_anchor, along_text_px=along_text_px)


def _rotated_y_bounds_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> tuple[float, float]:
    """Return min/max y-offset after rotation around the text anchor point."""
    _, _, y_min, y_max = _rotated_bounds_from_anchor(width_px, height_px, rotation_deg, text_anchor)
    return y_min, y_max


def _rotated_bounds_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> tuple[float, float, float, float]:
    """Return min/max x/y offsets after rotation around the text anchor point."""
    offsets = _rotated_corner_offsets(width_px, height_px, rotation_deg, text_anchor)
    min_x, min_y, max_x, max_y = aabb_from_points(offsets)
    return min_x, max_x, min_y, max_y


def _rotated_extreme_y_point_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
    *,
    choose_max: bool,
) -> tuple[float, float]:
    """Return the rotated corner offset at the top/bottom edge used as a leader contact."""
    points = _rotated_corner_offsets(width_px, height_px, rotation_deg, text_anchor)
    return max(points, key=lambda point: point[1]) if choose_max else min(points, key=lambda point: point[1])


def calculate_label_bounds(label: dict) -> tuple[float, float, float, float]:
    """Return absolute left/right/top/bottom coordinates for a label after rotation."""
    width_px = float(label["width_px"])
    height_px = float(label["height_px"])
    rotation_deg = float(label.get("rotation_deg", 0.0))
    text_anchor = str(label.get("text_anchor", "middle"))
    middle_x = float(label["middle_x"])
    middle_y = float(label["middle_y"])
    if rotation_deg == 0.0:
        x_min_offset, x_max_offset = anchor_x_bounds(width_px, text_anchor)
        half_height = 0.5 * max(0.0, height_px)
        return (
            middle_x + x_min_offset,
            middle_x + x_max_offset,
            middle_y - half_height,
            middle_y + half_height,
        )
    x_min_offset, x_max_offset, y_min_offset, y_max_offset = _rotated_bounds_from_anchor(
        width_px,
        height_px,
        rotation_deg,
        text_anchor,
    )
    return (
        middle_x + x_min_offset,
        middle_x + x_max_offset,
        middle_y + y_min_offset,
        middle_y + y_max_offset,
    )


def calculate_label_y_bounds(label: dict) -> tuple[float, float]:
    """Return absolute top/bottom y coordinates for a label after rotation."""
    _, _, top_y, bottom_y = calculate_label_bounds(label)
    return top_y, bottom_y


def _label_covers_keep_clear(label: Mapping[str, object], boxes: Sequence[Aabb]) -> bool:
    """Whether a label's box intersects one of ``boxes`` (in the label's coordinates): the one test of a
    label over the definition text, for the packers and the fit count (OV-1395, OV-1398)."""

    if not boxes:
        return False
    corners = _label_corners(dict(label))
    return any(
        convex_polygons_intersect(
            corners,
            ((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y)),
            touching=False,
        )
        for box in boxes
    )


def _label_corners(label: dict, shift_y: float = 0.0, along_text_px: float = 0.0) -> list[tuple[float, float]]:
    """A label's box moved ``shift_y`` down, with ``along_text_px`` added at both ends of its text."""

    middle_x = float(label["middle_x"])
    middle_y = float(label["middle_y"]) + shift_y
    return list(
        translate_points(
            _rotated_corner_offsets(
                float(label["width_px"]),
                float(label["height_px"]),
                float(label.get("rotation_deg", 0.0)),
                str(label.get("text_anchor", "middle")),
                along_text_px,
            ),
            middle_x,
            middle_y,
        )
    )


def _project_polygon(points: Sequence[tuple[float, float]], axis_x: float, axis_y: float) -> tuple[float, float]:
    projections = [(x * axis_x) + (y * axis_y) for x, y in points]
    return min(projections), max(projections)


def _rotated_label_rects_intersect(label1: dict, label2: dict, label1_shift_y: float = 0.0) -> bool:
    return _shifted_label_overlaps(label1, _label_corners(label2), label1_shift_y, None, 0.0)


def _update_linear_label_leader_end(label: dict) -> None:
    contact_x_offset = float(label.get("label_contact_x_offset", 0.0))
    contact_y_offset = float(label.get("label_contact_y_offset", 0.0))
    label["leader_end_x"] = float(label["middle_x"]) + contact_x_offset
    label["leader_end_y"] = float(label["middle_y"]) + contact_y_offset


# Shifts this close to an edge of the overlap range are decided by the polygon test itself.
_SHIFT_DECISION_MARGIN_PX = 1e-6


class _PlacedLinearLabel(NamedTuple):
    """A label that no longer moves, or a box the labels keep clear of, with the shapes the collision
    pass reads: a label's box with its row gap (``_row_label_corners``), a keep-clear box as it is."""

    aabb: Aabb
    corners: list[tuple[float, float]]
    # A label: the moving label keeps its own row gap from it too. A keep-clear box: its bare box.
    keeps_row_gap: bool


def _vertical_overlap_shift_range(
    first: Sequence[tuple[float, float]],
    second: Sequence[tuple[float, float]],
) -> tuple[float, float] | None:
    """Return the open range of vertical shifts of ``first`` at which it overlaps ``second``.

    The range is the separating-axis test of ``convex_polygons_intersect`` solved for
    the shift. Text boxes are rectangles, so the first two edges of each give every
    axis. ``None`` means an axis is too close to horizontal for a safe bound.
    """
    low, high = -math.inf, math.inf
    for polygon in (first, second):
        for index in (0, 1):
            x_coord, y_coord = polygon[index]
            next_x, next_y = polygon[index + 1]
            edge_x = next_x - x_coord
            edge_y = next_y - y_coord
            axis_len = math.hypot(edge_x, edge_y)
            if axis_len <= _SHIFT_DECISION_MARGIN_PX:
                return None
            axis_x = -edge_y / axis_len
            axis_y = edge_x / axis_len
            min_first, max_first = _project_polygon(first, axis_x, axis_y)
            min_second, max_second = _project_polygon(second, axis_x, axis_y)
            # Overlap on this axis at shift s: shift * axis_y lies strictly inside (below, above).
            below = min_second + GEOMETRY_EPSILON - max_first
            above = max_second - GEOMETRY_EPSILON - min_first
            if abs(axis_y) <= 1e-12:
                if below > _SHIFT_DECISION_MARGIN_PX or above < -_SHIFT_DECISION_MARGIN_PX:
                    return math.inf, -math.inf
                if below >= -_SHIFT_DECISION_MARGIN_PX or above <= _SHIFT_DECISION_MARGIN_PX:
                    return None
                continue
            if abs(axis_y) < 1e-3:
                return None
            bounds = sorted((below / axis_y, above / axis_y))
            low = max(low, bounds[0])
            high = min(high, bounds[1])
    return low, high


def _shifted_label_overlaps(
    label: dict,
    other_corners: Sequence[tuple[float, float]],
    shift: float,
    overlap_range: tuple[float, float] | None,
    along_text_px: float,
) -> bool:
    """Polygon test of ``label`` moved by ``shift``, with ``along_text_px`` at both ends of its text,
    against fixed corners.

    A shift clearly inside or outside ``overlap_range`` is answered from the range;
    one within the margin of its edges, or any shift without a range, runs the test.
    """
    if overlap_range is not None:
        low, high = overlap_range
        if low + _SHIFT_DECISION_MARGIN_PX < shift < high - _SHIFT_DECISION_MARGIN_PX:
            return True
        if shift < low - _SHIFT_DECISION_MARGIN_PX or shift > high + _SHIFT_DECISION_MARGIN_PX:
            return False
    return convex_polygons_intersect(
        _label_corners(label, shift_y=shift, along_text_px=along_text_px), other_corners, touching=False
    )


def _vertical_label_shift_to_clear(
    label: dict,
    other_corners: Sequence[tuple[float, float]],
    overlap_range: tuple[float, float] | None,
    direction: float,
    fallback_shift: float,
    min_gap_px: float,
    along_text_px: float,
) -> float:
    shift_sign = -1.0 if direction < 0.0 else 1.0
    high = max(abs(float(fallback_shift)), float(min_gap_px), 1.0)
    while _shifted_label_overlaps(label, other_corners, shift_sign * high, overlap_range, along_text_px):
        high *= 2.0
        if high > max(abs(float(fallback_shift)) * 8.0, 1024.0):
            return fallback_shift

    low = 0.0
    for _ in range(20):
        mid = (low + high) / 2.0
        if _shifted_label_overlaps(label, other_corners, shift_sign * mid, overlap_range, along_text_px):
            low = mid
        else:
            high = mid
    return shift_sign * (high + float(min_gap_px))


def _place_linear_label_without_overlap(
    label: dict,
    placed_x_index: IntervalIndex,
    placed_by_id: dict[int, _PlacedLinearLabel],
    min_gap_px: float,
    direction: float,
    row_pitch: float | None = None,
) -> bool:
    """Move one linear label away from the axis until its bbox clears placed labels.

    With ``row_pitch`` (a horizontal ``above_feature`` label, which starts in the first row off its own
    feature), the label moves by whole rows: row k off its feature edge at ``linear_row_y``, the rows
    of the auto placement, instead of by the clearing shift plus ``min_gap_px``.

    The label moves only vertically (OV-841), so the work that does not depend on
    its height is done once per label instead of once per round:

    - the placed labels whose bbox x range meets the label's form its column; a
      round looks only at the column members whose y range its bbox has reached
      and not yet left, which the climb visits in order;
    - each neighbour's overlap range is solved for the start position and
      translated by the shift accumulated since. The rounding of the translation is
      far below ``_SHIFT_DECISION_MARGIN_PX``, inside which the polygon test decides.
    """
    original_y = float(label["middle_y"])
    row = 1
    row_gap = _row_label_margin_px(label)
    # The label's box at its start, bare (against keep-clear boxes) and with its row gap (against labels).
    original_corners = {0.0: _label_corners(label), row_gap: _label_corners(label, along_text_px=row_gap)}
    clearance = _SHIFT_DECISION_MARGIN_PX
    moves_up = direction < 0.0
    start_aabb = _row_label_aabb(label)
    column = [
        (placed_id, placed_by_id[placed_id])
        for placed_id in (
            cast(int, label_id)  # only int ids are inserted
            for label_id in placed_x_index.query(Interval(start_aabb.min_x, start_aabb.max_x), padding=min_gap_px)
        )
        if start_aabb.min_x < placed_by_id[placed_id].aabb.max_x + clearance
        and start_aabb.max_x > placed_by_id[placed_id].aabb.min_x - clearance
    ]
    # The order in which the climbing bbox reaches each column member's y range.
    column.sort(key=(lambda item: -item[1].aabb.max_y) if moves_up else (lambda item: item[1].aabb.min_y))
    next_member = 0
    reached: list[tuple[int, _PlacedLinearLabel]] = []
    start_ranges: dict[int, tuple[float, float] | None] = {}
    for _ in range(len(placed_by_id) + 1):
        required_shift = 0.0
        label_aabb = _row_label_aabb(label)
        accumulated_shift = float(label["middle_y"]) - original_y
        # The y clauses of ``Aabb.intersects(other, clearance)``: as the label moves
        # away from the axis, the first becomes true and stays true, the second
        # becomes false and stays false. The column already meets the x clauses.
        if moves_up:
            while next_member < len(column) and label_aabb.min_y < column[next_member][1].aabb.max_y + clearance:
                reached.append(column[next_member])
                next_member += 1
            reached = [item for item in reached if label_aabb.max_y > item[1].aabb.min_y - clearance]
        else:
            while next_member < len(column) and label_aabb.max_y > column[next_member][1].aabb.min_y - clearance:
                reached.append(column[next_member])
                next_member += 1
            reached = [item for item in reached if label_aabb.min_y < item[1].aabb.max_y + clearance]

        # ``reached`` now holds exactly the placed labels whose bbox the label's meets.
        for placed_id, other in reached:
            other_aabb = other.aabb
            along_text_px = row_gap if other.keeps_row_gap else 0.0
            if placed_id not in start_ranges:
                start_ranges[placed_id] = _vertical_overlap_shift_range(original_corners[along_text_px], other.corners)
            start_range = start_ranges[placed_id]
            overlap_range = (
                None
                if start_range is None
                else (start_range[0] - accumulated_shift, start_range[1] - accumulated_shift)
            )
            if not _shifted_label_overlaps(label, other.corners, 0.0, overlap_range, along_text_px):
                continue

            if moves_up:
                fallback_shift = (other_aabb.min_y - float(min_gap_px)) - label_aabb.max_y
            else:
                fallback_shift = (other_aabb.max_y + float(min_gap_px)) - label_aabb.min_y
            shift = _vertical_label_shift_to_clear(
                label, other.corners, overlap_range, direction, fallback_shift, min_gap_px, along_text_px
            )
            required_shift = min(required_shift, shift) if moves_up else max(required_shift, shift)

        if required_shift == 0.0:
            return abs(float(label["middle_y"]) - original_y) > 0.001

        if row_pitch is None:
            label["middle_y"] = float(label["middle_y"]) + required_shift
        else:
            # The clearing shift without the stacking gap, rounded up to whole rows.
            clear = abs(required_shift) - float(min_gap_px)
            row += max(1, math.ceil(clear / row_pitch - _SHIFT_DECISION_MARGIN_PX))
            label["middle_y"] = linear_row_y(float(label["leader_start_y"]), direction * row_pitch, row)

    return abs(float(label["middle_y"]) - original_y) > 0.001


def _resolve_above_feature_label_overlaps(
    labels: list[dict], min_gap_px: float, keep_clear: Sequence[Aabb] = (), row_pitch: float | None = None
) -> None:
    """Stack above-feature labels away from the axis when their boxes, each with its row gap along its
    text, collide, and off the ``keep_clear`` boxes (the definition text), which are placed first. With
    ``row_pitch``, a horizontal label moves by whole rows of that pitch off its feature."""
    if not labels:
        return

    for place_above in (True, False):
        side_labels = [label for label in labels if bool(label.get("above_feature_place_above", True)) is place_above]
        if not side_labels:
            continue
        placed_by_id: dict[int, _PlacedLinearLabel] = {}
        placed_x_index = IntervalIndex(bucket_size=_above_feature_label_bucket_size(side_labels, min_gap_px))
        for box in keep_clear:
            corners = [(box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y)]
            placed_by_id[len(placed_by_id)] = _PlacedLinearLabel(box, corners, keeps_row_gap=False)
            placed_x_index.insert(len(placed_by_id) - 1, Interval(box.min_x, box.max_x))
        direction = -1.0 if place_above else 1.0
        for label in sorted(side_labels, key=lambda item: float(item.get("feature_anchor_x", item["middle_x"]))):
            pitch = row_pitch if float(label.get("rotation_deg", 0.0)) == 0.0 else None
            if _place_linear_label_without_overlap(label, placed_x_index, placed_by_id, min_gap_px, direction, pitch):
                label["leader_line"] = True
                _update_linear_label_leader_end(label)
            corners = _row_label_corners(label)
            placed = _PlacedLinearLabel(Aabb(*aabb_from_points(corners)), corners, keeps_row_gap=True)
            placed_id = len(placed_by_id)
            placed_by_id[placed_id] = placed
            placed_x_index.insert(placed_id, Interval(placed.aabb.min_x, placed.aabb.max_x))


def _passes_orthogroup_label_eligibility(
    feature_object,
    member_ids: set[str | int] | None,
    top_member_ids: set[str | int] | None,
) -> bool:
    if member_ids is None:
        return True
    source_feature_index = _source_feature_index(feature_object)
    if source_feature_index is None:
        raw_source_feature_index: Any = getattr(
            feature_object,
            "source_feature_index",
            None,
        )
        try:
            source_feature_index = int(raw_source_feature_index)
        except (TypeError, ValueError):
            source_feature_index = None
    feature_svg_id = compute_feature_object_hash(feature_object)
    identities: set[str | int | None] = {feature_svg_id}
    if source_feature_index is not None:
        identities.add(source_feature_index)
    if not identities.intersection(member_ids):
        return True
    return top_member_ids is not None and bool(identities.intersection(top_member_ids))


# Room a label needs inside its feature beyond its advance, in em (both ends together), so letters keep off
# the feature's ends and arrowhead (Owner D-20).
_EMBEDDED_LABEL_MARGIN_EM = 0.25


class LinearLabelFrame(NamedTuple):
    """Where a record's feature labels are stacked and counted: its Linear row and the x of the record's
    start in that row. A record translation moves the record with its labels (OV-1188)."""

    row: int
    x: float


def linear_label_rows(frames: Sequence[LinearLabelFrame] | None, record_count: int) -> list[list[tuple[int, float]]]:
    """The records of each Linear row in row order, each with the x of its start (``frames``; each of the
    ``record_count`` records alone in its row without them)."""

    resolved = list(frames) if frames is not None else [LinearLabelFrame(index, 0.0) for index in range(record_count)]
    if len(resolved) != record_count:
        raise ValueError("label frames must give one frame for every record")
    rows: dict[int, list[tuple[int, float]]] = {}
    for index, frame in sorted(enumerate(resolved), key=lambda item: (item[1].row, item[1].x, item[0])):
        rows.setdefault(frame.row, []).append((index, frame.x))
    return list(rows.values())


class LinearRecordLabels(NamedTuple):
    """One record's feature labels before ``place_linear_label_row`` stacks them: ``labels`` in placement
    order (an external one has no row yet), and what the row needs of the record. ``keep_clear`` are
    boxes in the labels' record coordinates that no label may cover: the record's definition text
    beside its features (OV-1395)."""

    labels: list[dict]
    top_feature_y: float
    bottom_feature_y: float
    max_label_height: float
    alignment_width: float
    keep_clear: tuple[Aabb, ...] = ()


def linear_record_labels(
    feature_dict,
    genome_length,
    alignment_width,
    cds_height,
    strandedness,
    track_layout,
    track_axis_gap,
    profile: LinearRenderProfile,
    label_font_size: float | None = None,
    orthogroup_label_member_ids: set[str | int] | None = None,
    orthogroup_label_top_member_ids: set[str | int] | None = None,
    labels_overrides_only: bool = False,
    feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
    record_transform: RecordDisplayTransform | None = None,
    label_filtering: Mapping[str, Any] | None = None,
) -> LinearRecordLabels:
    """One record's feature labels, measured and put on or beside their features.

    A feature outside the label display scope (``labels_overrides_only``, or
    not selected by orthogroup eligibility) is labeled only by its per-feature
    override. ``label_filtering``: the compiled maps when the caller holds them.
    A label is as wide as its text's advance in SVG px, the box the renderer
    anchors the text on (OV-1294).
    """
    cfg = profile.config
    labels: list[dict] = []
    feature_track_positions = {}  # Store feature track positions
    axis_gap_factor = (
        (float(track_axis_gap) / float(cds_height))
        if (track_axis_gap is not None and float(cds_height) > 0.0)
        else None
    )

    # Get configuration values
    length_threshold = cfg.labels.length_threshold.linear
    length_param = determine_length_parameter(genome_length, length_threshold)
    font_family = cfg.objects.text.font_family
    font_size = (
        float(label_font_size)
        if label_font_size is not None
        else cfg.labels.font_size.linear.for_length_param(length_param)
    )
    label_rendering = profile.label_rendering
    force_above_feature = profile.label_placement == "above_feature"
    if force_above_feature and label_rendering != "auto":
        raise ValueError(
            "label_rendering embedded_only|external_only cannot be used with label_placement above_feature"
        )
    base_rotation_deg = profile.label_rotation
    interval = cfg.canvas.dpi
    if label_filtering is None:
        label_filtering = label_filter_maps(cfg.labels.filtering.as_dict())
    # First pass: Calculate feature track positions
    for feature_id, feature_object in feature_dict.items():
        if len(feature_object.coordinates) == 0:
            continue

        # Get feature track info
        coordinate = feature_object.coordinates[0]  # Use first coordinate for strand
        strand = get_strand(coordinate.strand)
        feature_track_id = feature_object.feature_track_id

        if feature_lane_geometry is not None:
            lane = feature_lane_geometry.lane_for(
                strand=strand,
                track_id=feature_track_id,
                separate_strands=strandedness,
            )
            track_top_y, track_y_position, track_bottom_y = lane.positions
        else:
            # Compatibility fallback for direct callers that have not measured lanes.
            factors = calculate_feature_position_factors_linear(
                strand,
                feature_track_id,
                strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            track_top_y, track_y_position, track_bottom_y = (
                cds_height * factors[0],
                cds_height * factors[1],
                cds_height * factors[2],
            )

        # Store track position for this feature
        feature_track_positions[feature_id] = {
            "middle_y": track_y_position,
            "top_y": track_top_y,
            "bottom_y": track_bottom_y,
        }

    if feature_track_positions:
        top_feature_y_limit = min(pos["top_y"] for pos in feature_track_positions.values())
        bottom_feature_y_limit = max(pos["bottom_y"] for pos in feature_track_positions.values())
    else:
        top_feature_y_limit = 0.0
        bottom_feature_y_limit = 0.0

    # Second pass: Process labels
    max_bbox_height = 0
    for feature_id, feature_object in reversed(list(feature_dict.items())):
        overrides_only = labels_overrides_only or not _passes_orthogroup_label_eligibility(
            feature_object,
            orthogroup_label_member_ids,
            orthogroup_label_top_member_ids,
        )
        if overrides_only and not profile.forced_labels:
            continue
        feature_label_text = get_label_text(
            feature_object, label_filtering, overrides_only=overrides_only,
            feature_override=getattr(feature_object, "feature_override", None),
        )
        if not feature_label_text:
            continue
        if feature_id not in feature_track_positions:
            continue

        # Calculate label dimensions and positions
        bbox_width_px, _ = calculate_svg_bbox_dimensions(feature_label_text, font_family, font_size, interval)
        _, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
        if bbox_height_px > max_bbox_height:
            max_bbox_height = bbox_height_px
        # Find the longest segment and its middle point
        longest_segment_length = 0
        coordinate_strand: str | None = None
        longest_segment_start: float = 0
        longest_segment_end: float = 0

        for coordinate in feature_object.coordinates:
            coordinate_strand = get_strand(coordinate.strand)
            start = int(coordinate.start)
            end = int(coordinate.end)
            segment_length = abs(end - start + 1)
            if segment_length > longest_segment_length:
                longest_segment_start = start
                longest_segment_end = end
                longest_segment_length = segment_length

        if getattr(feature_object, "display_parts", None) is not None:
            projected_segment = display_label_segment(feature_object, genome_length, record_transform=record_transform)
            if projected_segment is None:
                continue
            longest_segment_start = projected_segment.start
            longest_segment_end = projected_segment.end
            coordinate_strand = projected_segment.strand

        # Calculate normalized positions
        # `position` is only used in arithmetic but is declared int in layout/linear_coords.py,
        # while display fragments are floats (F-s5-1).
        normalized_start = normalize_position_to_linear_track(
            cast(int, longest_segment_start), genome_length, alignment_width
        )
        normalized_end = normalize_position_to_linear_track(
            cast(int, longest_segment_end), genome_length, alignment_width
        )
        longest_segment_length_in_pixels = abs(normalized_end - normalized_start) + 1

        assignment = getattr(feature_object, "placement", None)
        secondary_below = (
            assignment is not None and assignment.side == "below"
            and assignment.level > 0 and assignment.strand_pool == "combined"
        )
        label_below_feature = secondary_below or (
            strandedness and coordinate_strand != "positive"
        )
        if force_above_feature:
            # Above-feature mode tilts labels in the opposite direction of the user-provided angle.
            # In separate-strands mode, negative strand labels are mirrored and placed below features.
            is_negative_separate = label_below_feature
            label_rotation_deg = base_rotation_deg if is_negative_separate else -base_rotation_deg
            label_text_anchor = "start" if label_rotation_deg != 0.0 else "middle"
        else:
            # Keep auto mode horizontal regardless of configured rotation.
            label_rotation_deg = 0.0
            label_text_anchor = "middle"

        feature_anchor_x = (normalized_start + normalized_end) / 2.0
        label_anchor_x = feature_anchor_x
        label_contact_x_offset = 0.0
        label_contact_y_offset = 0.0
        label_leader_start_y = 0.0
        if force_above_feature:
            place_above_feature = not label_below_feature
            if label_rotation_deg == 0.0:
                y_min_offset, y_max_offset = _rotated_y_bounds_from_anchor(
                    bbox_width_px, bbox_height_px, label_rotation_deg, label_text_anchor
                )
                label_contact_y_offset = y_max_offset if place_above_feature else y_min_offset
            else:
                label_contact_x_offset, label_contact_y_offset = _rotated_extreme_y_point_from_anchor(
                    bbox_width_px,
                    bbox_height_px,
                    label_rotation_deg,
                    label_text_anchor,
                    choose_max=place_above_feature,
                )
            label_anchor_x = feature_anchor_x - label_contact_x_offset
        bbox_start = label_anchor_x - (bbox_width_px / 2)
        bbox_end = label_anchor_x + (bbox_width_px / 2)

        # Get actual feature track position
        feature_position = feature_track_positions[feature_id]
        feature_y = feature_position["middle_y"]
        feature_top_y = feature_position["top_y"]
        feature_bottom_y = feature_position["bottom_y"]

        # Create base label entry
        label_entry = {
            "label_text": feature_label_text,
            "source_feature_index": feature_object.source_feature_index,
            "middle": label_anchor_x,
            "start": bbox_start,
            "end": bbox_end,
            "middle_x": label_anchor_x,
            "width_px": bbox_width_px,
            "height_px": bbox_height_px,
            "strand": coordinate_strand,
            "secondary_below": secondary_below,
            "feature_middle_y": feature_y,  # Use actual feature position
            "feature_top_y": feature_top_y,
            "feature_bottom_y": feature_bottom_y,
            "feature_start_x": normalized_start,
            "feature_end_x": normalized_end,
            "feature_anchor_x": feature_anchor_x,
            "font_size": font_size,
            "font_family": font_family,
            "rotation_deg": label_rotation_deg,
            "text_anchor": label_text_anchor,
        }

        # Determine if label should be embedded
        if force_above_feature:
            y_min_offset, y_max_offset = _rotated_y_bounds_from_anchor(
                bbox_width_px, bbox_height_px, label_rotation_deg, label_text_anchor
            )
            label_vertical_gap = max(1.0, bbox_height_px * 0.05)
            is_negative_separate = label_below_feature
            if is_negative_separate:
                # Keep rotated label top edge below the feature bottom.
                label_y = feature_bottom_y + label_vertical_gap - y_min_offset
                label_leader_start_y = feature_bottom_y
            else:
                # Keep rotated label bottom edge above the feature top.
                label_y = feature_top_y - label_vertical_gap - y_max_offset
                label_leader_start_y = feature_top_y
            # A horizontal label moves to the first row off its feature when its row is placed
            # (``place_linear_label_row``, LP-D120); the row's tallest label sets the row pitch.
            label_entry.update(
                {
                    "middle_y": label_y,
                    "is_embedded": True,
                    "track_id": "track_0",
                    "above_feature_place_above": not is_negative_separate,
                    "label_contact_x_offset": label_contact_x_offset,
                    "label_contact_y_offset": label_contact_y_offset,
                    "leader_line": False,
                    "leader_start_x": feature_anchor_x,
                    "leader_start_y": label_leader_start_y,
                    "leader_end_x": feature_anchor_x,
                    "leader_end_y": label_y + label_contact_y_offset,
                }
            )
            labels.append(label_entry)
        elif (
            label_rendering != "external_only"
            and bbox_width_px + _EMBEDDED_LABEL_MARGIN_EM * font_size < longest_segment_length_in_pixels
        ):
            label_entry.update({"middle_y": feature_y, "is_embedded": True, "track_id": "track_0"})
            labels.append(label_entry)
        elif label_rendering != "embedded_only":
            label_entry.update({"middle_y": 0, "is_embedded": False})
            labels.append(label_entry)

    return LinearRecordLabels(
        labels, top_feature_y_limit, bottom_feature_y_limit, max_bbox_height, float(alignment_width)
    )


def place_linear_label_row(
    records: Sequence[tuple[LinearRecordLabels, float]],
    *,
    profile: LinearRenderProfile,
    track_layout: str,
) -> list[list[dict]]:
    """Stack the feature labels of the records in one Linear row, each with the x of its start in the
    row (records in row order), so that no two labels of the row overlap (OV-1188).

    The external labels take the lowest row free where the row draws them, in record order, and the
    rows of a side are one row pitch apart from the outermost feature of the row. ``above_feature``
    labels start at their own feature: a horizontal one in the first row off it, with a leader, a turned
    one on it (LP-D120); then they move away from the axis off each other. A record's labels keep
    record coordinates, and the first record's external labels take the rows they take alone. Returns
    each record's labels: those on their features first, then the external ones by row.
    """

    if not records:
        return []
    spacing_px = float(profile.config.labels.spacing.linear)
    below = str(track_layout).strip().lower() == "below"
    labelled = [record for record, _x in records if record.labels]
    max_height = max((record.max_label_height for record in labelled), default=0.0)
    top_feature_y = min((record.top_feature_y for record in labelled), default=0.0)
    bottom_feature_y = max((record.bottom_feature_y for record in labelled), default=0.0)
    bucket_size = _external_label_bucket_size(max(record.alignment_width for record, _x in records))
    track_height = max_height + spacing_px
    # The records' keep-clear boxes (definition text) in row coordinates.
    keep_clear = [box.translated(x, 0.0) for record, x in records for box in record.keep_clear]

    def side(label: Mapping[str, object]) -> tuple[float, float]:
        """The origin and pitch of the rows of a label's side."""

        if below or label["secondary_below"]:
            return bottom_feature_y, track_height
        return top_feature_y, -track_height

    # The external labels in row coordinates: the packer reads their x span, and their box in a row
    # where a keep-clear box lies (OV-1398).
    track_dict: defaultdict[str, list[dict]] = defaultdict(list)
    track_indexes: dict[str, IntervalIndex] = {}
    framed_by_id: dict[int, dict] = {}
    for record, x in records:
        for label in record.labels:
            if label["is_embedded"]:
                continue
            framed = {**label, "start": float(label["start"]) + x, "end": float(label["end"]) + x}
            blocked = _keep_clear_rows(label, x, side(label), keep_clear) if keep_clear else None
            track_id = f"track_{_find_lowest_available_track_indexed(track_dict, track_indexes, framed_by_id, framed, blocked)}"
            label["track_id"] = framed["track_id"] = track_id
            track_dict[track_id].append(framed)
            framed_id = len(framed_by_id)
            framed_by_id[framed_id] = framed
            _insert_external_label_index(track_indexes, track_id, framed_id, framed, bucket_size)

    if profile.label_placement == "above_feature":
        for record, _x in records:
            for label in record.labels:
                if float(label["rotation_deg"]) == 0.0:
                    _above_feature_label_in_first_row(label, track_height)
        # Every label is at its feature here; the resolver reads their anchor and box x.
        on_features = [(label, x) for record, x in records for label in record.labels]
        framed_labels = [
            {
                **label,
                "middle_x": float(label["middle_x"]) + x,
                "feature_anchor_x": float(label["feature_anchor_x"]) + x,
            }
            for label, x in on_features
        ]
        _resolve_above_feature_label_overlaps(framed_labels, max(1.0, max_height * 0.05), keep_clear, track_height)
        for (label, _x), framed in zip(on_features, framed_labels):
            label["middle_y"] = framed["middle_y"]
            if framed["leader_line"]:
                label["leader_line"] = True
                _update_linear_label_leader_end(label)

    placed: list[list[dict]] = []
    for record, _x in records:
        external = sorted(
            (label for label in record.labels if not label["is_embedded"]), key=lambda label: label["track_id"]
        )
        for label in external:
            # Compact vertical positioning for external labels: row N of its side
            # (``linear_label_in_row`` moves a label to another row of it).
            label["row_origin_y"], label["row_pitch"] = side(label)
            label["middle_y"] = linear_row_y(label["row_origin_y"], label["row_pitch"], linear_label_row(label))
        placed.append([label for label in record.labels if label["is_embedded"]] + external)
    return placed


def _keep_clear_rows(
    label: Mapping[str, object], x: float, side: tuple[float, float], keep_clear: Sequence[Aabb]
) -> Callable[[int], bool]:
    """Whether ``label`` of the record at ``x``, in a row of its ``side`` (origin, pitch), would lie over
    one of the row's ``keep_clear`` boxes (row coordinates)."""

    in_row = {**label, "middle_x": float(str(label["middle_x"])) + x}
    origin, pitch = side
    return lambda row: _label_covers_keep_clear({**in_row, "middle_y": linear_row_y(origin, pitch, row)}, keep_clear)


def linear_row_y(origin_y: float, pitch: float, row: int) -> float:
    """The y of a label in row ``row`` of a side that starts at ``origin_y`` (the features' edge) and steps
    ``pitch`` (the tallest label plus ``labels.spacing.linear``, away from the axis): the one owner of
    a Linear label row's position, for the auto rows and the rows of horizontal ``above_feature``
    labels (LP-D120)."""

    return origin_y + pitch * row


def _above_feature_label_in_first_row(label: dict, track_height: float) -> None:
    """Put a horizontal ``above_feature`` label in the first row off its own feature, where the auto
    placement puts its first row, with a leader from the feature edge, as an auto label in its rows has
    (LP-D120). A turned label stays on its feature."""

    above = bool(label["above_feature_place_above"])
    edge_y = float(label["leader_start_y"])
    label["middle_y"] = linear_row_y(edge_y, -track_height if above else track_height, 1)
    label["leader_line"] = True
    _update_linear_label_leader_end(label)


def linear_label_row(label: Mapping[str, object]) -> int:
    """The row of a label (0: on its feature, or a label that is not in a row)."""

    return int(str(label.get("track_id", "track_0")).split("_")[1])


def linear_label_side(label: Mapping[str, object]) -> tuple[float, float] | None:
    """The rows an external, horizontal label stacks in (``None`` for a label that has no rows to
    move to: embedded, above-feature mode, or without a source feature)."""

    if (
        label.get("is_embedded")
        or label.get("leader_line")
        or label.get("source_feature_index") is None
        or "row_pitch" not in label
    ):
        return None
    return float(cast(float, label["row_origin_y"])), float(cast(float, label["row_pitch"]))


def linear_label_in_row(label: dict, side: tuple[float, float], row: int) -> dict:
    """A copy of an external label in row ``row`` of its ``side`` (``linear_label_side``), at the same x
    (LP-D82 ``RGN-YIELD-ROW``)."""

    origin, pitch = side
    return {**label, "middle_y": linear_row_y(origin, pitch, row), "track_id": f"track_{row}"}


class FeatureLabelRow(NamedTuple):
    """A feature label moved to another row of its side for a region label (LP-D82 ``RGN-YIELD-ROW``),
    by its record and source feature."""

    record_index: int
    source_feature_index: int
    row: int


def with_feature_label_rows(all_labels: Sequence[list[dict]], rows: Sequence[FeatureLabelRow]) -> list[list[dict]]:
    """Each record's labels with the listed feature labels moved to their rows."""

    moved = {(row.record_index, row.source_feature_index): row.row for row in rows}

    def placed(record_index: int, label: dict) -> dict:
        side = linear_label_side(label)
        row = moved.get((record_index, cast(int, label.get("source_feature_index"))))
        return label if side is None or row is None else linear_label_in_row(label, side, row)

    return [[placed(record_index, label) for label in labels] for record_index, labels in enumerate(all_labels)]


def linear_leader_polyline(label: dict) -> tuple[Point, ...]:
    """The leader the record drawer draws for a Linear label (empty for an embedded label)."""

    if label.get("leader_line"):
        return (
            (label["leader_start_x"], label["leader_start_y"]),
            (label["leader_end_x"], label["leader_end_y"]),
        )
    if label["is_embedded"]:
        return ()
    label_middle_y = float(label["middle_y"])
    label_height = float(label["height_px"])
    if label_middle_y >= float(label["feature_middle_y"]):
        # Label is below the feature track, so connect to its upper edge.
        label_edge_y = label_middle_y - (0.45 * label_height)
    else:
        # Label is above the feature track, so connect to its lower edge.
        label_edge_y = label_middle_y + (0.45 * label_height)
    return ((label["middle"], label["feature_middle_y"]), (label["middle"], label_edge_y))


def _row_label_corners(label: dict) -> list[tuple[float, float]]:
    """The corners of a label's box with the room it keeps at both ends of its text in its row: the one
    shape of a row label for the packer, the above-feature resolver and the fit count."""

    return _label_corners(label, along_text_px=_row_label_margin_px(label))


def _row_label_aabb(label: dict) -> Aabb:
    return Aabb(*aabb_from_points(_row_label_corners(label)))


def linear_label_shapes(labels: Sequence[dict], x: float = 0.0) -> list[LabelShape]:
    """One record's drawn Linear labels as the fit count sees them, moved ``x`` right (the record's
    start in its row): each row label with the room it keeps beside it, so a pair the rows separate
    for their gap counts as an overlap."""

    return [
        LabelShape(
            "linear",
            (translate_points(_row_label_corners(label), x, 0.0),),
            translate_points(linear_leader_polyline(label), x, 0.0),
        )
        for label in labels
    ]


def _linear_row_label_fits(
    row: Sequence[tuple[Sequence[dict], float]],
    definition_boxes: Sequence[Sequence[Aabb]] | None = None,
) -> tuple[list[LabelFit], bool]:
    """The fit of each record's drawn feature labels in one Linear row (records in row order, each with
    the x of its start), counted together with the shared rule (their layout boxes and leaders): a pair
    of labels of two records counts for the later record (OV-1188). A label whose letters cover one of
    its record's ``definition_boxes`` (record coordinates, as its labels) counts once for its record
    (OV-1395). Also whether an overlapping pair has a label drawn on its feature."""

    owners: list[int] = []
    labels: list[dict] = []
    shapes: list[LabelShape] = []
    for position, (record_labels, x) in enumerate(row):
        owners.extend([position] * len(record_labels))
        labels.extend(record_labels)
        shapes.extend(linear_label_shapes(record_labels, x))
    counts = [[0, 0, 0] for _record in row]
    overlaps = overlapping_pairs([shape.parts for shape in shapes], touching=False)
    for left, right in overlaps:
        counts[max(owners[left], owners[right])][0] += 1
    for left, right in crossing_pairs([shape.leader for shape in shapes]):
        counts[max(owners[left], owners[right])][1] += 1
    for position, ((record_labels, _x), boxes) in enumerate(zip(row, definition_boxes or [()] * len(row), strict=True)):
        counts[position][2] += sum(_label_covers_keep_clear(label, boxes) for label in record_labels)
    embedded_overlap = any(labels[left]["is_embedded"] or labels[right]["is_embedded"] for left, right in overlaps)
    return [LabelFit(overlap, crossing, covered) for overlap, crossing, covered in counts], embedded_overlap


def linear_row_label_fit(row: Sequence[tuple[Sequence[dict], float]]) -> LabelFit:
    """The drawn feature labels of one Linear row (each record's labels with the x of its start, its
    ``LinearLabelFrame``) counted together, as ``linear_label_fit`` counts them."""

    return sum(_linear_row_label_fits(row)[0], FITS)


# The drawing attribute that carries a Linear figure's label fit to the check.
LINEAR_LABEL_FIT_ATTR = "_gbdraw_linear_label_fit"


@dataclass(frozen=True)
class RecordRegionLabels:
    """One record's drawn region (annotation) labels: their font sizes, how many pairs of
    them and other labels or lines still overlap (LP-D68), and which labels (by placement
    order) those are."""

    font_sizes: tuple[float, ...] = ()
    overlaps: int = 0
    conflicted: tuple[int, ...] = ()


@dataclass(frozen=True)
class FeatureLabelRows:
    """Feature label rows a figure needs room for (LP-D82 ``RGN-YIELD-ROW``): the figure is
    assembled again with ``rows`` and kept when that clears the ``overlaps`` of the region
    labels that asked for them (``region_labels``: record index, placement order) and makes
    no other count larger (``keeps_feature_label_rows``)."""

    rows: tuple[FeatureLabelRow, ...]
    overlaps: int
    region_labels: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class LinearLabelFit:
    """The drawn Linear labels that do not fit, by record.

    ``embedded_overlap`` marks an overlapping pair with a label drawn on its
    feature: the features are drawn on top of each other. ``new_rows`` are the
    feature label rows this figure has no room for yet; ``region_conflicts`` the
    region labels (record index, placement order) still over other labels or lines.
    """

    records: tuple[LabelFitRecord, ...] = ()
    embedded_overlap: bool = False
    new_rows: FeatureLabelRows | None = None
    region_conflicts: frozenset[tuple[int, int]] = frozenset()

    @property
    def fit(self) -> LabelFit:
        return sum((record.fit for record in self.records), FITS)

    @property
    def outcome(self) -> LabelFitOutcome:
        return label_fit_outcome([record.fit for record in self.records])


def linear_label_fit(
    record_labels: Sequence[tuple[str, Sequence[dict]]],
    font_size: float,
    region_labels: Sequence[RecordRegionLabels] | None = None,
    new_rows: FeatureLabelRows | None = None,
    *,
    frames: Sequence[LinearLabelFrame] | None = None,
    definition_boxes: Sequence[Sequence[Aabb]] | None = None,
) -> LinearLabelFit:
    """Count one figure's drawn Linear labels with the shared rule (LP-D48 SCOPE-ALL-LABELS).

    The feature labels of the records in one row (``frames``, where their
    labels were stacked; each record alone in its row without them) are
    counted together (OV-1188). ``definition_boxes`` are each record's definition text boxes in its
    labels' record coordinates; a label over one warns (OV-1395). ``region_labels`` are each record's region
    labels as their placer counted them (LP-D68); a record whose region labels
    overlap other labels or lines gets a ``linear_region`` record, which only warns.
    """

    regions = list(region_labels) if region_labels is not None else [RecordRegionLabels()] * len(record_labels)
    fits: list[LabelFit] = [FITS] * len(record_labels)
    embedded_overlap = False
    for row in linear_label_rows(frames, len(record_labels)):
        row_fits, row_embedded_overlap = _linear_row_label_fits(
            [(record_labels[index][1], x) for index, x in row],
            [definition_boxes[index] for index, _x in row] if definition_boxes is not None else None,
        )
        embedded_overlap = embedded_overlap or row_embedded_overlap
        for (index, _x), fit in zip(row, row_fits, strict=True):
            fits[index] = fit
    failing: list[LabelFitRecord] = []
    for (record_id, labels), region, fit in zip(record_labels, regions, fits, strict=True):
        if not fit.fits:
            failing.append(LabelFitRecord(record_id, "linear", len(labels), font_size, fit))
        if region.overlaps:
            failing.append(
                LabelFitRecord(
                    record_id,
                    "linear_region",
                    len(region.font_sizes),
                    max(region.font_sizes),
                    LabelFit(0, 0, region_overlaps=region.overlaps),
                )
            )
    conflicts = frozenset((index, label) for index, region in enumerate(regions) for label in region.conflicted)
    return LinearLabelFit(tuple(failing), embedded_overlap, new_rows, conflicts)


def keeps_feature_label_rows(before: LinearLabelFit, after: LinearLabelFit) -> bool:
    """Whether the figure assembled with ``before.new_rows`` (``after``) is kept: it clears the
    overlaps of the region labels that asked for the rows (LP-D77) and makes no other count larger."""

    if before.new_rows is None:
        return False
    old, new = before.fit, after.fit
    return (
        after.region_conflicts.isdisjoint(before.new_rows.region_labels)
        and new.region_overlaps <= old.region_overlaps - before.new_rows.overlaps
        and new.overlap_pairs <= old.overlap_pairs
        and new.crossing_pairs <= old.crossing_pairs
        and new.definition_overlaps <= old.definition_overlaps
    )


def check_linear_label_fit(
    label_fit: LinearLabelFit,
    verify: Callable[[Mapping[str, LabelFitValue]], LabelFitVerification],
    *,
    features_separated: bool,
) -> LabelFitReport | None:
    """Raise or warn for a drawn Linear figure whose labels do not fit.

    Placed Linear labels keep apart, so overlaps come from labels drawn on
    features drawn on top of each other; their one row puts overlapping
    features on separate tracks (LP-D52 ``E-ROW``), verified with ``verify``.
    A figure drawn with crossing leaders or with region labels over other
    labels or lines warns (LP-D68) and returns its report, which has no rows.
    """

    outcome = label_fit.outcome
    if outcome == "fits":
        return None
    rows: tuple[LabelFitRow, ...] = ()
    if outcome == "stop" and label_fit.embedded_overlap and not features_separated:
        values = SEPARATE_FEATURES_VALUES
        rows = (switch_row("separate_features", values, verify(values), offer_notice=True),)
    report = LabelFitReport(label_fit.records, rows=rows, drawn=outcome == "warn")
    raise_or_warn(report)
    return report


__all__ = [
    "LINEAR_LABEL_FIT_ATTR",
    "FeatureLabelRow",
    "FeatureLabelRows",
    "LinearLabelFit",
    "LinearLabelFrame",
    "LinearRecordLabels",
    "RecordRegionLabels",
    "calculate_label_bounds",
    "check_linear_label_fit",
    "keeps_feature_label_rows",
    "linear_label_fit",
    "linear_label_in_row",
    "linear_label_row",
    "linear_label_rows",
    "linear_label_side",
    "linear_row_label_fit",
    "linear_record_labels",
    "place_linear_label_row",
    "with_feature_label_rows",
    "calculate_label_y_bounds",
    "check_label_overlap",
    "find_lowest_available_track",
    "linear_label_shapes",
    "linear_leader_polyline",
]
