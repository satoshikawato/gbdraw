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
from ..config.models import LinearRenderProfile
from ..features.coordinates import get_strand
from ..features.ids import compute_feature_object_hash
from ..core.record_metadata import _source_feature_index
from ..core.text import calculate_bbox_dimensions
from .label_ink import (
    label_underline_counted,
    label_underline_in_box,
    line_half_width_em,
    linear_label_ink,
    underline_reach_em,
)
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
    label_fit_outcome,
    label_fits_by_kind,
    overlapping_pairs,
    raise_or_warn,
    switch_row,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

def check_label_overlap(label1, label2):
    """Check if two labels overlap horizontally"""
    return not (label1["end"] < label2["start"] or label2["end"] < label1["start"])


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
    return Interval(float(label["start"]), float(label["end"]))


def _linear_label_aabb(label: dict) -> Aabb:
    left, right, top, bottom = calculate_label_bounds(label)
    return Aabb(left, top, right, bottom)


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
) -> int:
    """Find the lowest non-overlapping external label track using index candidates."""
    label_interval = _linear_label_interval(label)
    track_num = 1
    while True:
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
) -> tuple[tuple[float, float], ...]:
    return text_box_corner_offsets(width_px, height_px, rotation_deg, text_anchor)


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


def _label_corners(label: dict, shift_y: float = 0.0) -> list[tuple[float, float]]:
    middle_x = float(label["middle_x"])
    middle_y = float(label["middle_y"]) + shift_y
    return list(
        translate_points(
            _rotated_corner_offsets(
            float(label["width_px"]),
            float(label["height_px"]),
            float(label.get("rotation_deg", 0.0)),
            str(label.get("text_anchor", "middle")),
            ),
            middle_x,
            middle_y,
        )
    )


def _project_polygon(points: Sequence[tuple[float, float]], axis_x: float, axis_y: float) -> tuple[float, float]:
    projections = [(x * axis_x) + (y * axis_y) for x, y in points]
    return min(projections), max(projections)


def _rotated_label_rects_intersect(label1: dict, label2: dict, label1_shift_y: float = 0.0) -> bool:
    return _shifted_label_overlaps(label1, _label_corners(label2), label1_shift_y, None)


def _update_linear_label_leader_end(label: dict) -> None:
    contact_x_offset = float(label.get("label_contact_x_offset", 0.0))
    contact_y_offset = float(label.get("label_contact_y_offset", 0.0))
    label["leader_end_x"] = float(label["middle_x"]) + contact_x_offset
    label["leader_end_y"] = float(label["middle_y"]) + contact_y_offset


# Shifts this close to an edge of the overlap range are decided by the polygon test itself.
_SHIFT_DECISION_MARGIN_PX = 1e-6


class _PlacedLinearLabel(NamedTuple):
    """A label that no longer moves, with the shapes the collision pass reads."""

    label: dict
    aabb: Aabb
    corners: list[tuple[float, float]]


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
) -> bool:
    """Polygon test of ``label`` moved by ``shift`` against fixed corners.

    A shift clearly inside or outside ``overlap_range`` is answered from the range;
    one within the margin of its edges, or any shift without a range, runs the test.
    """
    if overlap_range is not None:
        low, high = overlap_range
        if low + _SHIFT_DECISION_MARGIN_PX < shift < high - _SHIFT_DECISION_MARGIN_PX:
            return True
        if shift < low - _SHIFT_DECISION_MARGIN_PX or shift > high + _SHIFT_DECISION_MARGIN_PX:
            return False
    return convex_polygons_intersect(_label_corners(label, shift_y=shift), other_corners, touching=False)


def _vertical_label_shift_to_clear(
    label: dict,
    other_corners: Sequence[tuple[float, float]],
    overlap_range: tuple[float, float] | None,
    direction: float,
    fallback_shift: float,
    min_gap_px: float,
) -> float:
    shift_sign = -1.0 if direction < 0.0 else 1.0
    high = max(abs(float(fallback_shift)), float(min_gap_px), 1.0)
    while _shifted_label_overlaps(label, other_corners, shift_sign * high, overlap_range):
        high *= 2.0
        if high > max(abs(float(fallback_shift)) * 8.0, 1024.0):
            return fallback_shift

    low = 0.0
    for _ in range(20):
        mid = (low + high) / 2.0
        if _shifted_label_overlaps(label, other_corners, shift_sign * mid, overlap_range):
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
) -> bool:
    """Move one linear label away from the axis until its bbox clears placed labels.

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
    original_corners = _label_corners(label)
    clearance = _SHIFT_DECISION_MARGIN_PX
    moves_up = direction < 0.0
    start_aabb = _linear_label_aabb(label)
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
        label_aabb = _linear_label_aabb(label)
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
            if placed_id not in start_ranges:
                start_ranges[placed_id] = _vertical_overlap_shift_range(original_corners, other.corners)
            start_range = start_ranges[placed_id]
            overlap_range = (
                None
                if start_range is None
                else (start_range[0] - accumulated_shift, start_range[1] - accumulated_shift)
            )
            if not _shifted_label_overlaps(label, other.corners, 0.0, overlap_range):
                continue

            if moves_up:
                fallback_shift = (other_aabb.min_y - float(min_gap_px)) - label_aabb.max_y
            else:
                fallback_shift = (other_aabb.max_y + float(min_gap_px)) - label_aabb.min_y
            shift = _vertical_label_shift_to_clear(
                label, other.corners, overlap_range, direction, fallback_shift, min_gap_px
            )
            required_shift = min(required_shift, shift) if moves_up else max(required_shift, shift)

        if required_shift == 0.0:
            return abs(float(label["middle_y"]) - original_y) > 0.001

        label["middle_y"] = float(label["middle_y"]) + required_shift

    return abs(float(label["middle_y"]) - original_y) > 0.001


def _resolve_above_feature_label_overlaps(labels: list[dict], min_gap_px: float) -> None:
    """Stack rotated above-feature labels away from the axis when their bboxes collide."""
    if not labels:
        return

    for place_above in (True, False):
        side_labels = [label for label in labels if bool(label.get("above_feature_place_above", True)) is place_above]
        if not side_labels:
            continue
        placed_by_id: dict[int, _PlacedLinearLabel] = {}
        placed_x_index = IntervalIndex(bucket_size=_above_feature_label_bucket_size(side_labels, min_gap_px))
        direction = -1.0 if place_above else 1.0
        for label in sorted(side_labels, key=lambda item: float(item.get("feature_anchor_x", item["middle_x"]))):
            if _place_linear_label_without_overlap(label, placed_x_index, placed_by_id, min_gap_px, direction):
                label["leader_line"] = True
                _update_linear_label_leader_end(label)
            placed = _PlacedLinearLabel(label, _linear_label_aabb(label), _label_corners(label))
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


def prepare_label_list_linear(
    feature_dict,
    genome_length,
    alignment_width,
    genome_size_normalization_factor,
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
):
    """
    Prepares a list of labels for linear genome visualization with proper track organization.

    A feature outside the label display scope (``labels_overrides_only``, or
    not selected by orthogroup eligibility) is labeled only by its per-feature
    override. ``label_filtering``: the compiled maps when the caller holds them.
    """
    cfg = profile.config
    embedded_labels = []
    external_labels = []
    track_dict = defaultdict(list)
    external_track_indexes: dict[str, IntervalIndex] = {}
    external_label_by_id: dict[int, dict] = {}
    feature_track_positions = {}  # Store feature track positions
    track_layout_normalized = str(track_layout).strip().lower()
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
    label_spacing_px = float(cfg.labels.spacing.linear)
    external_bucket_size = _external_label_bucket_size(alignment_width)
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
    underline_box_em = underline_reach_em(font_family, "central") if label_underline_in_box("linear") else None
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
        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
        if underline_box_em is not None:
            # The underline belongs to the label's box: the placer leaves it room (LP-D93 study).
            reach_em = underline_box_em + line_half_width_em(float(font_size), "linear")
            bbox_height_px = max(bbox_height_px, 2.0 * reach_em * float(font_size))
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
            cast(int, longest_segment_start), genome_length, alignment_width, genome_size_normalization_factor
        )
        normalized_end = normalize_position_to_linear_track(
            cast(int, longest_segment_end), genome_length, alignment_width, genome_size_normalization_factor
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
            track_dict["track_0"].append(label_entry)
        elif label_rendering != "external_only" and bbox_width_px < longest_segment_length_in_pixels:
            label_entry.update({"middle_y": feature_y, "is_embedded": True, "track_id": "track_0"})
            track_dict["track_0"].append(label_entry)
        elif label_rendering != "embedded_only":
            # For external labels: find lowest possible track
            label_entry.update({"middle_y": 0, "is_embedded": False})

            # Find lowest track where label can be placed without overlaps
            best_track = _find_lowest_available_track_indexed(
                track_dict,
                external_track_indexes,
                external_label_by_id,
                label_entry,
            )
            label_entry["track_id"] = f"track_{best_track}"
            track_dict[f"track_{best_track}"].append(label_entry)
            external_label_id = len(external_label_by_id)
            external_label_by_id[external_label_id] = label_entry
            _insert_external_label_index(
                external_track_indexes,
                f"track_{best_track}",
                external_label_id,
                label_entry,
                external_bucket_size,
            )

    # Process embedded labels
    if "track_0" in track_dict:
        if force_above_feature:
            above_feature_collision_gap = max(1.0, max_bbox_height * 0.05)
            _resolve_above_feature_label_overlaps(track_dict["track_0"], above_feature_collision_gap)
        for label in track_dict["track_0"]:
            embedded_labels.append(label)

    track_height = max_bbox_height + label_spacing_px
    for track_id in sorted(track_dict.keys()):
        if track_id == "track_0":
            continue

        track_labels = track_dict[track_id]
        track_num = int(track_id.split("_")[1])

        for label in track_labels:
            # Compact vertical positioning for external labels: row N of its side
            # (``linear_label_in_row`` moves a label to another row of it).
            if track_layout_normalized == "below" or label["secondary_below"]:
                label["row_origin_y"], label["row_pitch"] = bottom_feature_y_limit, track_height
            else:
                label["row_origin_y"], label["row_pitch"] = top_feature_y_limit, -track_height
            label["middle_y"] = label["row_origin_y"] + label["row_pitch"] * track_num
            external_labels.append(label)

    return embedded_labels + external_labels


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
    return {**label, "middle_y": origin + pitch * row, "track_id": f"track_{row}"}


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


def linear_label_shapes(labels: list[dict]) -> list[LabelShape]:
    """One record's drawn Linear labels as the fit count sees them."""

    underline = label_underline_counted("linear")
    return [
        LabelShape(
            "linear",
            (tuple(_label_corners(label)),),
            linear_label_ink(label, leader, underline=True).drawn if underline and leader else leader,
        )
        for label, leader in ((label, linear_leader_polyline(label)) for label in labels)
    ]


def linear_record_label_fit(labels: Sequence[dict]) -> LabelFit:
    """One record's drawn feature labels counted with the shared rule (their layout boxes and leaders)."""

    return label_fits_by_kind(linear_label_shapes(list(labels)), touching=False).get("linear", FITS)


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
) -> LinearLabelFit:
    """Count one figure's drawn Linear labels with the shared rule (LP-D48 SCOPE-ALL-LABELS).

    ``region_labels`` are each record's region labels as their placer counted
    them (LP-D68); a record whose region labels overlap other labels or lines
    gets a ``linear_region`` record, which only warns.
    """

    regions = list(region_labels) if region_labels is not None else [RecordRegionLabels()] * len(record_labels)
    failing: list[LabelFitRecord] = []
    embedded_overlap = False
    for (record_id, labels), region in zip(record_labels, regions, strict=True):
        fit = linear_record_label_fit(labels)
        if not fit.fits:
            failing.append(LabelFitRecord(record_id, "linear", len(labels), font_size, fit))
            if fit.overlap_pairs and not embedded_overlap:
                embedded_overlap = any(
                    labels[left]["is_embedded"] or labels[right]["is_embedded"]
                    for left, right in overlapping_pairs(
                        [shape.parts for shape in linear_label_shapes(list(labels))], touching=False
                    )
                )
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
    "RecordRegionLabels",
    "calculate_label_bounds",
    "check_linear_label_fit",
    "keeps_feature_label_rows",
    "linear_label_fit",
    "linear_label_in_row",
    "linear_label_row",
    "linear_label_side",
    "linear_record_label_fit",
    "with_feature_label_rows",
    "calculate_label_y_bounds",
    "check_label_overlap",
    "find_lowest_available_track",
    "linear_label_shapes",
    "linear_leader_polyline",
    "prepare_label_list_linear",
]
