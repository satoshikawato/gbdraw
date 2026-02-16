#!/usr/bin/env python
# coding: utf-8

"""
Circular label placement utilities.

This module contains the circular/arc-based label placement logic that used to live in
`gbdraw.labels.placement`.
"""

import math

from .filtering import get_label_text, preprocess_label_filtering
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..features.coordinates import get_strand
from ..core.text import calculate_bbox_dimensions
from ..core.sequence import determine_length_parameter
from ..layout.common import calculate_cds_ratio
from ..layout.circular import calculate_feature_position_factors_circular

# Keep dense large-font labels from being pushed excessively far from features.
MIN_BBOX_GAP_RATIO = 0.01
MIN_BBOX_GAP_FLOOR_PX = 1.2
HEAVY_CLUSTER_RELAX_MAX_LABELS = 70
FULL_SCAN_LABEL_LIMIT = 70
MAX_EXPANDED_SHIFT_DEG = 90.0
RELAX_CENTER_DELTA_CAP_DEG = 35.0
LEGACY_PLACEMENT_LABEL_THRESHOLD = 50
HEMISPHERE_AXIS_NEUTRAL_DEG = 8.0
HEMISPHERE_AXIS_EPS_PX = 1.0
HEMISPHERE_REFINE_MAX_SHIFT_DEG = 90.0
HEMISPHERE_REFINE_STEP_DEG = 0.25
DENSE_DISTANCE_COMPACTION_MIN_LABELS = 30
DENSE_DISTANCE_COMPACTION_MIN_GAP_RELAX = 1
LEADER_START_PROXIMITY_EPS_PX = 0.75
LEADER_START_ORDER_GAP_CANDIDATES = (1.0, 0.75, 0.5, 0.25, 0.0)
MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX = 6.0
MIN_OUTER_LABEL_TEXT_CLEARANCE_PX = 13.0
OUTER_LABEL_FEATURE_CLEARANCE_SAMPLE_CAP_DEG = 18.0
PLACE_LABELS_BASE_MARGIN_PX = 2.0
IMPROVED_LABELS_BASE_MARGIN_PX = 1.5
RESOLVE_OVERLAP_MARGIN_EXTRA_PX = 4.2


def minimum_bbox_gap_px(label1: dict, label2: dict, base_margin_px: float = 0.0) -> float:
    """Return a minimum spacing margin in pixels derived from label bbox size."""
    width_scale = max(float(label1.get("width_px", 0.0)), float(label2.get("width_px", 0.0)))
    height_scale = max(float(label1.get("height_px", 0.0)), float(label2.get("height_px", 0.0)))
    ratio_gap = MIN_BBOX_GAP_RATIO * max(width_scale, height_scale)
    resolve_extra = 0.0
    if ("min_outer_start_radius_px" in label1) or ("min_outer_start_radius_px" in label2):
        resolve_extra = RESOLVE_OVERLAP_MARGIN_EXTRA_PX
    return max(float(base_margin_px) + resolve_extra, MIN_BBOX_GAP_FLOOR_PX, ratio_gap)


def _count_label_overlaps(labels: list[dict], total_length: int, use_min_gap: bool) -> int:
    """Count overlapping label pairs using the existing overlap predicates."""
    overlap_count = 0
    for idx in range(len(labels)):
        for jdx in range(idx + 1, len(labels)):
            label1 = labels[idx]
            label2 = labels[jdx]
            if use_min_gap:
                margin = minimum_bbox_gap_px(label1, label2, base_margin_px=0.0)
                y_margin = margin
                x_margin = margin
            else:
                y_margin = 0.1
                x_margin = 1.0
            if y_overlap(label1, label2, total_length, y_margin) and x_overlap(
                label1, label2, minimum_margin=x_margin
            ):
                overlap_count += 1
    return overlap_count


def _label_y_bounds(label: dict, total_len: int, minimum_margin: float) -> tuple[float, float]:
    """Return (min_y, max_y) for label bbox using the same baseline model as rendering."""
    label_start_y = float(label["start_y"])
    label_angle = (360.0 * (float(label["middle"]) / total_len)) % 360.0
    label_height = float(label["height_px"])
    margin_half = 0.5 * float(minimum_margin)

    if label.get("is_inner", False) is False:
        if 0 <= label_angle < 10:  # baseline = text-after-edge
            max_y = label_start_y
            min_y = label_start_y - 1.0 * label_height - margin_half
        elif 10 <= label_angle < 170:  # baseline = middle
            max_y = label_start_y + 0.5 * label_height + margin_half
            min_y = label_start_y - 0.5 * label_height - margin_half
        elif 170 <= label_angle < 190:  # baseline = hanging
            max_y = label_start_y + 1.0 * label_height + margin_half
            min_y = label_start_y - margin_half
        elif 190 <= label_angle < 350:  # baseline = middle
            max_y = label_start_y + 0.5 * label_height + margin_half
            min_y = label_start_y - 0.5 * label_height - margin_half
        else:  # baseline = text-after-edge
            max_y = label_start_y
            min_y = label_start_y - 1.0 * label_height - margin_half
    else:
        if 0 <= label_angle < 10:  # baseline = hanging
            max_y = label_start_y + 1.0 * label_height + margin_half
            min_y = label_start_y
        elif 10 <= label_angle < 170:  # baseline = middle
            max_y = label_start_y + 0.5 * label_height + margin_half
            min_y = label_start_y - 0.5 * label_height - margin_half
        elif 170 <= label_angle < 190:  # baseline = middle
            max_y = label_start_y + 0.5 * label_height + margin_half
            min_y = label_start_y - 0.5 * label_height - margin_half
        elif 190 <= label_angle < 350:  # baseline = middle
            max_y = label_start_y + 0.5 * label_height + margin_half
            min_y = label_start_y - 0.5 * label_height - margin_half
        else:  # baseline = hanging
            max_y = label_start_y + 1.0 * label_height + margin_half
            min_y = label_start_y
    return min_y, max_y


def _label_x_bounds(label: dict, minimum_margin: float) -> tuple[float, float]:
    """Return (min_x, max_x) for label bbox using current text-anchor semantics."""
    label_start_x = float(label["start_x"])
    label_width = float(label["width_px"])
    margin_half = 0.5 * float(minimum_margin)

    if label.get("is_inner", False) is False:
        if label_start_x > 0:
            max_x = label_start_x + label_width + margin_half
            min_x = label_start_x - margin_half
        else:
            max_x = label_start_x + margin_half
            min_x = label_start_x - label_width - margin_half
    else:
        if label_start_x > 0:
            max_x = label_start_x + margin_half
            min_x = label_start_x - label_width - margin_half
        else:
            max_x = label_start_x + label_width + margin_half
            min_x = label_start_x - margin_half
    return min_x, max_x


def _leader_anchor_candidates(label: dict, total_length: int) -> list[tuple[float, float]]:
    """
    Return discrete candidate anchor points on the label bbox perimeter:
    4 corners + 4 edge midpoints.
    """
    if "start_x" not in label or "start_y" not in label:
        return []
    if "width_px" not in label or "height_px" not in label:
        return []

    min_x, max_x = _label_x_bounds(label, minimum_margin=0.0)
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=0.0)
    mid_x = 0.5 * (min_x + max_x)
    mid_y = 0.5 * (min_y + max_y)
    return [
        (min_x, min_y),
        (min_x, max_y),
        (max_x, min_y),
        (max_x, max_y),
        (mid_x, min_y),
        (mid_x, max_y),
        (min_x, mid_y),
        (max_x, mid_y),
    ]


def _leader_start_meta(label: dict, total_length: int | None) -> tuple[str, float, float, float, float] | None:
    """
    Return leader anchor metadata as (side, fixed_coord, lower, upper, preferred).

    - side: one of "left", "right", "bottom", "top"
    - fixed_coord: x for left/right, y for bottom/top
    - lower/upper: movable coordinate bounds on the chosen edge
    - preferred: preferred movable coordinate before clamping
    """
    if total_length is None:
        return None
    if "middle_x" not in label or "middle_y" not in label:
        return None
    start_x = float(label.get("start_x", 0.0))
    start_y = float(label.get("start_y", 0.0))
    min_x, max_x = _label_x_bounds(label, minimum_margin=0.0)
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=0.0)
    if (max_x - min_x) <= 1e-9 or (max_y - min_y) <= 1e-9:
        return None

    middle_x = float(label["middle_x"])
    middle_y = float(label["middle_y"])

    # Avoid snapping exactly to corners; corner hits look noisy for dense labels.
    corner_inset = 0.75
    y_min_usable = min_y + corner_inset if (max_y - min_y) > (2.0 * corner_inset) else min_y
    y_max_usable = max_y - corner_inset if (max_y - min_y) > (2.0 * corner_inset) else max_y
    x_min_usable = min_x + corner_inset if (max_x - min_x) > (2.0 * corner_inset) else min_x
    x_max_usable = max_x - corner_inset if (max_x - min_x) > (2.0 * corner_inset) else max_x

    # Horizontal labels look cleaner when leaders attach on the feature-facing side.
    # Keep the attachment coordinate biased toward text anchor (`start_x/start_y`)
    # so each leader visually maps back to its own label row.
    anchor_x = float(label.get("start_x", start_x))
    anchor_y = float(label.get("start_y", start_y))
    if middle_x >= max_x:
        return "right", max_x, y_min_usable, y_max_usable, anchor_y
    if middle_x <= min_x:
        return "left", min_x, y_min_usable, y_max_usable, anchor_y
    if middle_y >= max_y:
        return "top", max_y, x_min_usable, x_max_usable, anchor_x
    if middle_y <= min_y:
        return "bottom", min_y, x_min_usable, x_max_usable, anchor_x

    # Fallback (point inside bbox): attach to nearest edge.
    left_d = abs(middle_x - min_x)
    right_d = abs(max_x - middle_x)
    bottom_d = abs(middle_y - min_y)
    top_d = abs(max_y - middle_y)
    best = min((left_d, "left"), (right_d, "right"), (bottom_d, "bottom"), (top_d, "top"))[1]
    if best == "left":
        return "left", min_x, y_min_usable, y_max_usable, middle_y
    if best == "right":
        return "right", max_x, y_min_usable, y_max_usable, middle_y
    if best == "bottom":
        return "bottom", min_y, x_min_usable, x_max_usable, middle_x
    return "top", max_y, x_min_usable, x_max_usable, middle_x


def _leader_start_point(label: dict, total_length: int | None = None) -> tuple[float, float]:
    """Choose a visually stable leader anchor on the label bbox perimeter."""
    start_x = float(label.get("start_x", 0.0))
    start_y = float(label.get("start_y", 0.0))
    meta = _leader_start_meta(label, total_length)
    if meta is None:
        return start_x, start_y
    side, fixed_coord, lower, upper, preferred = meta
    movable = min(max(float(preferred), float(lower)), float(upper))
    if side in ("left", "right"):
        return fixed_coord, movable
    return movable, fixed_coord


def _leader_start_group_score(group: list[dict], coords: list[float], side: str) -> tuple[int, int, float, float]:
    """Score leader-start coordinates for one edge group."""
    inversion_count = 0
    for idx in range(len(coords)):
        for jdx in range(idx + 1, len(coords)):
            if coords[idx] > (coords[jdx] + 1e-9):
                inversion_count += 1

    close_pairs = 0
    for idx in range(len(coords) - 1):
        if (coords[idx + 1] - coords[idx]) < LEADER_START_PROXIMITY_EPS_PX:
            close_pairs += 1

    leader_sum = 0.0
    preference_deviation = 0.0
    for item, coord in zip(group, coords):
        label = item["label"]
        if side in ("left", "right"):
            x = float(item["fixed"])
            y = float(coord)
        else:
            x = float(coord)
            y = float(item["fixed"])
        leader_sum += math.hypot(x - float(label["middle_x"]), y - float(label["middle_y"]))
        preference_deviation += abs(float(coord) - float(item["preferred"]))

    return inversion_count, close_pairs, leader_sum, preference_deviation


def _fit_monotonic_leader_coords(group: list[dict], min_gap: float) -> list[float] | None:
    """Fit edge coordinates under bounds while enforcing monotonic order."""
    if not group:
        return []

    coords = [min(max(float(item["preferred"]), float(item["lower"])), float(item["upper"])) for item in group]

    for _ in range(max(3, len(coords) * 2)):
        changed = False
        for idx in range(1, len(coords)):
            lower_bound = max(float(group[idx]["lower"]), coords[idx - 1] + min_gap)
            if coords[idx] < lower_bound:
                coords[idx] = lower_bound
                changed = True
            if coords[idx] > (float(group[idx]["upper"]) + 1e-9):
                return None

        for idx in range(len(coords) - 2, -1, -1):
            upper_bound = min(float(group[idx]["upper"]), coords[idx + 1] - min_gap)
            if coords[idx] > upper_bound:
                coords[idx] = upper_bound
                changed = True
            if coords[idx] < (float(group[idx]["lower"]) - 1e-9):
                return None

        if not changed:
            break

    for idx, coord in enumerate(coords):
        if coord < (float(group[idx]["lower"]) - 1e-9) or coord > (float(group[idx]["upper"]) + 1e-9):
            return None
        if idx > 0 and coord < (coords[idx - 1] + min_gap - 1e-9):
            return None
    return coords


def _refine_leader_start_points(labels: list[dict], total_length: int) -> list[dict]:
    """
    Improve readability by reducing leader-end inversions and excessive crowding.

    The optimization runs per bbox edge and keeps each anchor on the same edge.
    """
    if len(labels) < 2:
        return labels

    edge_groups: dict[str, list[dict]] = {"left": [], "right": [], "bottom": [], "top": []}
    for label in labels:
        meta = _leader_start_meta(label, total_length)
        if meta is None:
            continue
        side, fixed, lower, upper, preferred = meta
        if side not in edge_groups:
            continue
        if side in ("left", "right"):
            current_coord = float(label.get("leader_start_y", label.get("start_y", 0.0)))
            order_coord = float(label.get("middle_y", 0.0))
        else:
            current_coord = float(label.get("leader_start_x", label.get("start_x", 0.0)))
            order_coord = float(label.get("middle_x", 0.0))

        edge_groups[side].append(
            {
                "label": label,
                "fixed": fixed,
                "lower": lower,
                "upper": upper,
                "preferred": preferred,
                "current": current_coord,
                "order": order_coord,
            }
        )

    for side, group in edge_groups.items():
        if len(group) < 2:
            continue
        group.sort(key=lambda item: float(item["order"]))
        current_coords = [float(item["current"]) for item in group]
        current_score = _leader_start_group_score(group, current_coords, side)
        best_coords = current_coords
        best_score = current_score

        for min_gap in LEADER_START_ORDER_GAP_CANDIDATES:
            candidate_coords = _fit_monotonic_leader_coords(group, float(min_gap))
            if candidate_coords is None:
                continue
            candidate_score = _leader_start_group_score(group, candidate_coords, side)
            if candidate_score < best_score:
                best_coords = candidate_coords
                best_score = candidate_score

        if best_score >= current_score:
            continue

        for item, coord in zip(group, best_coords):
            label = item["label"]
            if side in ("left", "right"):
                label["leader_start_x"] = float(item["fixed"])
                label["leader_start_y"] = float(coord)
            else:
                label["leader_start_x"] = float(coord)
                label["leader_start_y"] = float(item["fixed"])
    return labels


def _assign_leader_start_points(labels: list[dict], total_length: int) -> list[dict]:
    """Populate `leader_start_x/y` used for drawing the first leader segment."""
    for label in labels:
        leader_x, leader_y = _leader_start_point(label, total_length)
        label["leader_start_x"] = leader_x
        label["leader_start_y"] = leader_y
    return _refine_leader_start_points(labels, total_length)


def assign_leader_start_points(labels: list[dict], total_length: int) -> list[dict]:
    """Public wrapper to (re)compute leader start anchors after label motion."""
    return _assign_leader_start_points(labels, total_length)


def _leader_length_px(label: dict, total_length: int | None = None) -> float:
    """Return leader length from label middle anchor to chosen label-side leader anchor."""
    if "middle_x" not in label or "middle_y" not in label:
        return 0.0
    leader_start_x, leader_start_y = _leader_start_point(label, total_length)
    return math.hypot(
        leader_start_x - float(label["middle_x"]),
        leader_start_y - float(label["middle_y"]),
    )


def _leader_distance_score(labels: list[dict], total_length: int | None = None) -> tuple[float, float]:
    """Return (sum_length, max_length) for leader segments."""
    if not labels:
        return 0.0, 0.0
    leader_lengths = [_leader_length_px(label, total_length) for label in labels]
    return sum(leader_lengths), max(leader_lengths)


def _label_bbox_min_radius(label: dict, total_length: int, minimum_margin: float = 0.0) -> float:
    """Return minimum distance from origin to the label bbox."""
    min_x, max_x = _label_x_bounds(label, minimum_margin=minimum_margin)
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=minimum_margin)
    closest_x = min(max(0.0, min_x), max_x)
    closest_y = min(max(0.0, min_y), max_y)
    return math.hypot(closest_x, closest_y)


def _build_outer_feature_radius_intervals(
    feature_dict: dict,
    total_length: int,
    radius: float,
    track_ratio: float,
    cfg: GbdrawConfig,
) -> list[tuple[float, float, float]]:
    """Build local genome intervals as (start, end, outer_radius_px)."""
    if total_length <= 0:
        return []

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = cfg.canvas.circular.track_ratio_factors[length_param][0]
    cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
    track_type = cfg.canvas.circular.track_type
    strandedness = cfg.canvas.strandedness

    intervals: list[tuple[float, float, float]] = []
    total_length_float = float(total_length)

    for feature_object in feature_dict.values():
        track_id = int(getattr(feature_object, "feature_track_id", 0))
        feature_location_list = getattr(feature_object, "location", [])
        list_of_coordinates = getattr(feature_object, "coordinates", [])

        for coordinate_idx, coordinate in enumerate(list_of_coordinates):
            if (
                coordinate_idx < len(feature_location_list)
                and getattr(feature_location_list[coordinate_idx], "kind", None) == "line"
            ):
                continue

            coordinate_start_raw = int(coordinate.start)
            coordinate_end_raw = int(coordinate.end)
            if coordinate_start_raw == coordinate_end_raw:
                continue

            coordinate_start = float(coordinate_start_raw % total_length)
            coordinate_end = float(coordinate_end_raw % total_length)
            if coordinate_end_raw > coordinate_start_raw and coordinate_end_raw <= total_length:
                coordinate_end = float(coordinate_end_raw)

            coordinate_strand = get_strand(coordinate.strand)
            factors = calculate_feature_position_factors_circular(
                total_length,
                coordinate_strand,
                track_ratio,
                cds_ratio,
                offset,
                track_type,
                strandedness,
                track_id,
            )
            feature_outer_radius = float(radius) * float(factors[2])

            if coordinate_start < coordinate_end:
                intervals.append((coordinate_start, coordinate_end, feature_outer_radius))
            else:
                intervals.append((coordinate_start, total_length_float, feature_outer_radius))
                intervals.append((0.0, coordinate_end, feature_outer_radius))

    return intervals


def _max_outer_feature_radius_at_position(
    intervals: list[tuple[float, float, float]],
    position: float,
    total_length: int,
) -> float:
    """Return the max feature outer radius that covers a genome position."""
    if not intervals or total_length <= 0:
        return 0.0
    position_wrapped = float(position % float(total_length))
    max_outer_radius = 0.0
    for start_pos, end_pos, outer_radius in intervals:
        if start_pos <= position_wrapped < end_pos:
            if outer_radius > max_outer_radius:
                max_outer_radius = outer_radius
    return max_outer_radius


def _sample_label_genome_positions(label: dict, total_length: int) -> list[float]:
    """Sample genome positions covered by a label for local feature clearance checks."""
    if total_length <= 0:
        return []

    start_x = float(label.get("start_x", 0.0))
    start_y = float(label.get("start_y", 0.0))
    radius = math.hypot(start_x, start_y)
    if radius <= 1e-9:
        middle = float(label.get("middle", 0.0))
        return [middle % float(total_length)]

    base_angle_deg = math.degrees(math.atan2(start_y, start_x))
    width_px = max(0.0, float(label.get("width_px", 0.0)))
    half_span_deg = 0.5 * math.degrees(width_px / radius) if width_px > 0.0 else 0.0
    half_span_deg = min(OUTER_LABEL_FEATURE_CLEARANCE_SAMPLE_CAP_DEG, half_span_deg)

    sampled_angles = [base_angle_deg]
    if half_span_deg >= 0.2:
        sampled_angles.extend([base_angle_deg - half_span_deg, base_angle_deg + half_span_deg])
    if half_span_deg >= 4.0:
        sampled_angles.extend(
            [
                base_angle_deg - (0.5 * half_span_deg),
                base_angle_deg + (0.5 * half_span_deg),
            ]
        )

    sampled_positions: list[float] = []
    for sampled_angle in sampled_angles:
        angle_mod = (sampled_angle + 90.0) % 360.0
        sampled_positions.append((angle_mod / 360.0) * float(total_length))
    return sampled_positions


def _update_outer_label_minimum_radii_against_features(
    labels: list[dict],
    total_length: int,
    feature_radius_intervals: list[tuple[float, float, float]],
) -> list[dict]:
    """Raise minimum outer-label radii so text clears local feature tracks."""
    if not labels or not feature_radius_intervals:
        return labels

    for label in labels:
        sampled_positions = _sample_label_genome_positions(label, total_length)
        if not sampled_positions:
            continue

        local_outer_radius = 0.0
        for sampled_pos in sampled_positions:
            local_outer_radius = max(
                local_outer_radius,
                _max_outer_feature_radius_at_position(feature_radius_intervals, sampled_pos, total_length),
            )
        if local_outer_radius <= 0.0:
            continue

        required_radius = local_outer_radius + MIN_OUTER_LABEL_TEXT_CLEARANCE_PX
        current_min_radius = float(label.get("min_outer_start_radius_px", 0.0))
        if required_radius > current_min_radius:
            label["min_outer_start_radius_px"] = required_radius

    return labels


def _enforce_outer_label_minimum_radius(labels: list[dict], total_length: int) -> list[dict]:
    """Push outer labels outward until their text bboxes clear required radii."""
    for label in labels:
        required_outer_radius = float(label.get("min_outer_start_radius_px", 0.0))
        if required_outer_radius <= 0.0:
            continue

        start_x = float(label.get("start_x", 0.0))
        start_y = float(label.get("start_y", 0.0))
        current_radius = math.hypot(start_x, start_y)
        if _label_bbox_min_radius(label, total_length, minimum_margin=0.0) >= (required_outer_radius - 1e-6):
            continue

        if current_radius > 1e-9:
            angle_rad = math.atan2(start_y, start_x)
        else:
            angle_rad = math.radians(angle_from_middle(float(label["middle"]), total_length))
        current_half = _current_half_from_x(start_x, axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
        candidate_radius = max(current_radius, required_outer_radius)

        for _ in range(8):
            new_x = candidate_radius * math.cos(angle_rad)
            new_y = candidate_radius * math.sin(angle_rad)
            if current_half > 0 and new_x < 0:
                new_x = abs(new_x)
            elif current_half < 0 and new_x > 0:
                new_x = -abs(new_x)
            label["start_x"] = float(new_x)
            label["start_y"] = float(new_y)

            bbox_min_radius = _label_bbox_min_radius(label, total_length, minimum_margin=0.0)
            deficit = required_outer_radius - bbox_min_radius
            if deficit <= 1e-3:
                break
            candidate_radius += deficit + 0.5

    return labels


def _resolve_outer_label_overlaps_with_fixed_radii(
    labels: list[dict],
    total_length: int,
    *,
    max_angle_shift_deg: float = 40.0,
    step_deg: float = 0.2,
) -> list[dict]:
    """Resolve overlaps by shifting label angles while keeping each label radius fixed."""
    if len(labels) < 2:
        return labels

    labels = sort_labels(labels)
    if len(labels) > 2:
        # Move the optimization seam into the sparsest angular gap so labels that
        # are adjacent around 0/360 are optimized as neighbors.
        raw_targets = [angle_from_middle_unwrapped(float(label["middle"]), total_length) for label in labels]
        gap_values = [raw_targets[idx + 1] - raw_targets[idx] for idx in range(len(raw_targets) - 1)]
        gap_values.append(raw_targets[0] + 360.0 - raw_targets[-1])
        largest_gap_idx = max(range(len(gap_values)), key=gap_values.__getitem__)
        if largest_gap_idx < len(labels) - 1:
            rotate_by = largest_gap_idx + 1
            labels = labels[rotate_by:] + labels[:rotate_by]

    min_order_gap_deg = 0.05
    max_pair_steps = max(1, int(max_angle_shift_deg / max(step_deg, 1e-6)))

    prev_target_angle: float | None = None
    for label in labels:
        current_x = float(label.get("start_x", 0.0))
        current_y = float(label.get("start_y", 0.0))
        raw_target_angle = angle_from_middle_unwrapped(float(label["middle"]), total_length)
        if prev_target_angle is None:
            target_angle = raw_target_angle
        else:
            target_angle = normalize_angle_near_reference(raw_target_angle, prev_target_angle)
            while target_angle <= prev_target_angle:
                target_angle += 360.0
        current_angle = math.degrees(math.atan2(current_y, current_x))
        current_unwrapped = normalize_angle_near_reference(current_angle, target_angle)
        label["_fixed_radius_px"] = math.hypot(current_x, current_y)
        label["_target_angle_unwrapped_fixed"] = target_angle
        label["_angle_unwrapped_fixed"] = current_unwrapped
        prev_target_angle = target_angle

    prev_angle: float | None = None
    for label in labels:
        target_angle = float(label["_target_angle_unwrapped_fixed"])
        low_bound = target_angle - max_angle_shift_deg
        high_bound = target_angle + max_angle_shift_deg
        candidate = float(label["_angle_unwrapped_fixed"])
        if prev_angle is not None:
            candidate = max(candidate, prev_angle + min_order_gap_deg)
        candidate = min(max(candidate, low_bound), high_bound)
        if prev_angle is not None and candidate <= prev_angle + min_order_gap_deg:
            candidate = prev_angle + min_order_gap_deg
        label["_angle_unwrapped_fixed"] = candidate
        prev_angle = candidate

    def _apply_angle(label: dict) -> None:
        radius = float(label["_fixed_radius_px"])
        angle_unwrapped = float(label["_angle_unwrapped_fixed"])
        label["start_x"] = radius * math.cos(math.radians(angle_unwrapped))
        label["start_y"] = radius * math.sin(math.radians(angle_unwrapped))

    def _pair_overlaps(left: dict, right: dict) -> bool:
        min_gap_px = minimum_bbox_gap_px(left, right, base_margin_px=0.0)
        return y_overlap(left, right, total_length, min_gap_px) and x_overlap(
            left,
            right,
            minimum_margin=min_gap_px,
        )

    for label in labels:
        _apply_angle(label)

    for _ in range(120):
        changed = False
        for idx in range(len(labels) - 1):
            left = labels[idx]
            right = labels[idx + 1]
            if not _pair_overlaps(left, right):
                continue

            for _ in range(max_pair_steps):
                left_target = float(left["_target_angle_unwrapped_fixed"])
                right_target = float(right["_target_angle_unwrapped_fixed"])
                left_low = left_target - max_angle_shift_deg
                right_high = right_target + max_angle_shift_deg

                new_left = max(left_low, float(left["_angle_unwrapped_fixed"]) - (0.5 * step_deg))
                new_right = min(right_high, float(right["_angle_unwrapped_fixed"]) + (0.5 * step_deg))

                if idx > 0:
                    new_left = max(new_left, float(labels[idx - 1]["_angle_unwrapped_fixed"]) + min_order_gap_deg)
                if idx + 2 < len(labels):
                    new_right = min(new_right, float(labels[idx + 2]["_angle_unwrapped_fixed"]) - min_order_gap_deg)

                if new_right <= new_left + min_order_gap_deg:
                    break

                if (
                    abs(new_left - float(left["_angle_unwrapped_fixed"])) <= 1e-9
                    and abs(new_right - float(right["_angle_unwrapped_fixed"])) <= 1e-9
                ):
                    break

                left["_angle_unwrapped_fixed"] = new_left
                right["_angle_unwrapped_fixed"] = new_right
                _apply_angle(left)
                _apply_angle(right)
                changed = True

                if not _pair_overlaps(left, right):
                    break

        if not changed:
            break

    for label in labels:
        label.pop("_fixed_radius_px", None)
        label.pop("_target_angle_unwrapped_fixed", None)
        label.pop("_angle_unwrapped_fixed", None)

    return sort_labels(labels)


def _angle_deviation_score(labels: list[dict], total_length: int) -> tuple[float, float]:
    """
    Return (sum_delta, max_delta) where deltas are angular deviations from each
    label's target feature angle.
    """
    sum_delta = 0.0
    max_delta = 0.0
    for label in labels:
        if "angle_unwrapped" in label:
            target_angle = float(
                label.get(
                    "target_angle_unwrapped",
                    angle_from_middle_unwrapped(float(label["middle"]), total_length),
                )
            )
            delta = abs(float(label["angle_unwrapped"]) - target_angle)
        elif "start_x" in label and "start_y" in label:
            actual_angle = math.degrees(math.atan2(float(label["start_y"]), float(label["start_x"]))) % 360.0
            target_angle = angle_from_middle(float(label["middle"]), total_length)
            delta = abs(angle_delta_deg(actual_angle, target_angle))
        else:
            delta = 0.0
        sum_delta += delta
        max_delta = max(max_delta, delta)
    return sum_delta, max_delta


def _placement_score(labels: list[dict], total_length: int) -> tuple[int, int, int, float, float, float, float, float]:
    """Return tuple score used to compare placement candidates."""
    overlap_plain = _count_label_overlaps(labels, total_length, use_min_gap=False)
    overlap_min_gap = _count_label_overlaps(labels, total_length, use_min_gap=True)
    hemisphere_mismatch_count, hemisphere_mismatch_weight = _hemisphere_mismatch_metrics(labels, total_length)
    leader_sum, leader_max = _leader_distance_score(labels, total_length)
    sum_delta, max_delta = _angle_deviation_score(labels, total_length)
    return (
        overlap_plain,
        overlap_min_gap,
        hemisphere_mismatch_count,
        hemisphere_mismatch_weight,
        leader_sum,
        leader_max,
        sum_delta,
        max_delta,
    )


def y_overlap(label1, label2, total_len, minimum_margin):
    min_y1, max_y1 = _label_y_bounds(label1, total_len, minimum_margin)
    min_y2, max_y2 = _label_y_bounds(label2, total_len, minimum_margin)
    if min_y1 < min_y2:
        return max_y1 >= min_y2
    else:
        return max_y2 >= min_y1


def x_overlap(label1, label2, minimum_margin=1.0):
    min_x1, max_x1 = _label_x_bounds(label1, minimum_margin)
    min_x2, max_x2 = _label_x_bounds(label2, minimum_margin)
    if min_x1 < min_x2:
        return max_x1 >= min_x2
    else:
        return max_x2 >= min_x1


def calculate_angle_degrees(
    center_x, center_y, x, y, middle, start_angle, end_angle, total_length, x_radius, y_radius, normalize
):
    """Calculate the angle in degrees from a point (x, y) relative to the center of an ellipse."""
    # Normalize coordinates to a unit circle
    x_normalized = (x - center_x) / x_radius
    y_normalized = (y - center_y) / y_radius

    angle_radians = math.atan2(y_normalized, x_normalized)
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees


def calculate_coordinates(center_x, center_y, x_radius, y_radius, angle_degrees, middle, total_length):
    angle_radians = math.radians(angle_degrees)
    y = center_y + y_radius * math.sin(angle_radians)
    x = center_x + x_radius * math.cos(angle_radians)

    return x, y


def calculate_angle_for_y(center_y, y_radius, y):
    """Calculate the angle in degrees for a given y-coordinate on an ellipse."""
    if center_y - y_radius <= y <= center_y + y_radius:
        # Calculate the arcsine of the normalized y-coordinate
        angle_radians = math.asin((y - center_y) / y_radius)
        angle_degrees = math.degrees(angle_radians)
        return angle_degrees
    else:
        return None  # Indicates the y-coordinate is outside the ellipse's bounds


def angle_from_middle(middle: float, total_length: int) -> float:
    """Convert genome midpoint to circular angle (degrees)."""
    return (360.0 * (middle / total_length) - 90.0) % 360.0


def angle_from_middle_unwrapped(middle: float, total_length: int) -> float:
    """Convert genome midpoint to monotonic circular angle (degrees, unwrapped)."""
    return 360.0 * (middle / total_length) - 90.0


def angle_delta_deg(a: float, b: float) -> float:
    """Return the signed smallest angle delta from b to a in [-180, 180)."""
    return ((a - b + 180.0) % 360.0) - 180.0


def normalize_angle_near_reference(angle_deg: float, reference_unwrapped: float) -> float:
    """Project an angle to the nearest 360-degree branch around a reference angle."""
    normalized = angle_deg % 360.0
    return normalized + 360.0 * round((reference_unwrapped - normalized) / 360.0)


def _vertical_axis_distance_deg(angle_mod: float) -> float:
    """Return angular distance to the nearest vertical axis (90 or 270 degrees)."""
    return min(abs(angle_delta_deg(angle_mod, 90.0)), abs(angle_delta_deg(angle_mod, 270.0)))


def _current_half_from_x(x: float, axis_eps_px: float = HEMISPHERE_AXIS_EPS_PX) -> int:
    """Return +1 (right), -1 (left), or 0 (near vertical axis) from x-coordinate."""
    if x > axis_eps_px:
        return 1
    if x < -axis_eps_px:
        return -1
    return 0


def _preferred_half_for_label(
    label: dict,
    total_length: int,
    axis_neutral_deg: float = HEMISPHERE_AXIS_NEUTRAL_DEG,
) -> int:
    """
    Return preferred label half: +1 right, -1 left, or 0 (neutral near top/bottom).

    Prefer feature midpoint x when available so preference follows actual feature side.
    """
    feature_middle_x = label.get("feature_middle_x")
    if feature_middle_x is not None:
        feature_half = _current_half_from_x(float(feature_middle_x), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
        if feature_half != 0:
            return feature_half

    target_angle_unwrapped = float(
        label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
    )
    target_angle_mod = target_angle_unwrapped % 360.0
    if _vertical_axis_distance_deg(target_angle_mod) <= axis_neutral_deg:
        return 0
    return 1 if math.cos(math.radians(target_angle_mod)) >= 0 else -1


def _hemisphere_mismatch_metrics(
    labels: list[dict],
    total_length: int,
    axis_neutral_deg: float = HEMISPHERE_AXIS_NEUTRAL_DEG,
) -> tuple[int, float]:
    """
    Return (mismatch_count, weighted_mismatch_score) for hemisphere preference.

    Labels close to top/bottom are treated as neutral and excluded.
    """
    mismatch_count = 0
    weighted_score = 0.0
    for label in labels:
        preferred_half = _preferred_half_for_label(label, total_length, axis_neutral_deg)
        if preferred_half == 0:
            continue

        current_half = _current_half_from_x(float(label.get("start_x", 0.0)), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
        if current_half in (0, preferred_half):
            continue

        mismatch_count += 1
        target_angle_unwrapped = float(
            label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
        )
        target_angle_mod = target_angle_unwrapped % 360.0
        distance_from_axis = _vertical_axis_distance_deg(target_angle_mod)
        weighted_score += max(0.0, distance_from_axis - axis_neutral_deg)

    return mismatch_count, weighted_score


def place_labels_on_arc_fc(
    labels: list[dict],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    start_angle: float,
    end_angle: float,
    total_length: int,
) -> list[dict]:
    def check_overlap(label1, label2, total_length, margin):
        effective_margin = minimum_bbox_gap_px(label1, label2, base_margin_px=max(margin, PLACE_LABELS_BASE_MARGIN_PX))
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
        )

    if not labels:
        return []

    rearranged_labels: list[dict] = []
    labels = sort_labels(labels)
    angle_step = 0.25
    max_shift = 25.0
    max_steps = int(max_shift / angle_step)
    min_order_gap = 0.05
    prev_angle_unwrapped = -float("inf")

    for label in labels:
        target_angle_unwrapped = angle_from_middle_unwrapped(label["middle"], total_length)
        label["target_angle_unwrapped"] = target_angle_unwrapped
        label["target_angle"] = target_angle_unwrapped % 360.0
        lower_bound = prev_angle_unwrapped + min_order_gap
        best_x = None
        best_y = None
        best_overlaps = None
        best_angle_unwrapped = None

        for step in range(max_steps + 1):
            if step == 0:
                candidate_offsets = [0.0]
            else:
                offset = step * angle_step
                candidate_offsets = [offset, -offset]

            for offset in candidate_offsets:
                candidate_angle_unwrapped = target_angle_unwrapped + offset
                if candidate_angle_unwrapped < lower_bound:
                    continue
                candidate_x, candidate_y = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, candidate_angle_unwrapped, label["middle"], total_length
                )
                candidate = label.copy()
                candidate["start_x"] = candidate_x
                candidate["start_y"] = candidate_y

                overlap_count = 0
                for existing in rearranged_labels:
                    if check_overlap(candidate, existing, total_length, 0.1):
                        overlap_count += 1

                if overlap_count == 0:
                    label["start_x"] = candidate_x
                    label["start_y"] = candidate_y
                    best_overlaps = 0
                    best_angle_unwrapped = candidate_angle_unwrapped
                    break

                if best_overlaps is None or overlap_count < best_overlaps:
                    best_x = candidate_x
                    best_y = candidate_y
                    best_overlaps = overlap_count
                    best_angle_unwrapped = candidate_angle_unwrapped

            if best_overlaps == 0:
                break

        if best_overlaps != 0:
            if best_angle_unwrapped is None:
                # Fallback: preserve ordering even when local placement is crowded.
                best_angle_unwrapped = lower_bound
                best_x, best_y = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, best_angle_unwrapped, label["middle"], total_length
                )
            label["start_x"] = best_x if best_x is not None else 0.0
            label["start_y"] = best_y if best_y is not None else 0.0

        label["angle_unwrapped"] = best_angle_unwrapped if best_angle_unwrapped is not None else lower_bound
        prev_angle_unwrapped = label["angle_unwrapped"]
        rearranged_labels.append(label)

    return rearranged_labels


def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def sort_labels(labels):
    return sorted(labels, key=lambda x: x["middle"])


def _legacy_place_labels_on_arc_fc(
    labels: list[dict],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    start_angle: float,
    end_angle: float,
    total_length: int,
) -> list[dict]:
    """
    Legacy placement strategy (from pre-93dc098 implementation).

    This remains as a fallback for dense cases where the new optimizer leaves
    unresolved overlaps.
    """

    def move_label(label: dict, angle: float) -> tuple[float, float]:
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def check_overlap(label1: dict, label2: dict, total_length: int, margin: float) -> bool:
        return y_overlap(label1, label2, total_length, margin) and x_overlap(
            label1, label2, minimum_margin=2
        )

    if not labels:
        return []

    rearranged_labels = []
    labels = sort_labels(labels)
    current_angle = -75.0
    increment = 0.1

    for idx, label in enumerate(labels):
        if idx == 0:
            label["start_x"], label["start_y"] = calculate_coordinates(
                center_x, center_y, x_radius, y_radius, current_angle, label["middle"], total_length
            )
            rearranged_labels.append(label)
            continue

        new_angle = current_angle + increment
        if new_angle < -75:
            new_angle = -75.0
        elif -75 <= new_angle < 85:
            if label["middle"] > (total_length / 2) or idx >= len(labels) * (2 / 3):
                new_angle = 85.0
        elif 85 < new_angle < 90:
            new_angle = 90.0
        if new_angle > 269:
            new_angle = 269.0

        label["start_x"], label["start_y"] = calculate_coordinates(
            center_x, center_y, x_radius, y_radius, new_angle, label["middle"], total_length
        )
        while rearranged_labels and check_overlap(label, rearranged_labels[-1], total_length, 0.1):
            new_angle += 0.01
            label["start_x"], label["start_y"] = calculate_coordinates(
                center_x, center_y, x_radius, y_radius, new_angle, label["middle"], total_length
            )

        rearranged_labels.append(label)
        current_angle = new_angle

    return rearranged_labels


def _legacy_improved_label_placement_fc(
    labels: list[dict],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    feature_radius: float,
    total_length: int,
    start_angle: float,
    end_angle: float,
    y_margin: float = 0.1,
    max_iterations: int = 10000,
) -> list[dict]:
    """
    Legacy local optimizer (from pre-93dc098 implementation).

    This solver is intentionally conservative and is used as a fallback to
    improve dense-label robustness.
    """

    def move_label(label: dict, angle: float) -> tuple[float, float]:
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def calculate_angle_of_three_points(
        x1: float, y1: float, x2: float, y2: float, x3: float, y3: float
    ) -> float:
        v1 = (x1 - x2, y1 - y2)
        v2 = (x3 - x2, y3 - y2)
        dot_product = v1[0] * v2[0] + v1[1] * v2[1]
        mag1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
        mag2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)
        if mag1 == 0 or mag2 == 0:
            return 0.0
        cos_angle = dot_product / (mag1 * mag2)
        angle = math.acos(max(-1.0, min(1.0, cos_angle)))
        return math.degrees(angle)

    def check_overlap(label1: dict, label2: dict, total_length: int) -> bool:
        return y_overlap(label1, label2, total_length, y_margin) and x_overlap(
            label1, label2, minimum_margin=1
        )

    labels = sort_labels(labels)
    if len(labels) < 2:
        return labels

    for _ in range(max_iterations):
        changes_made = False
        for idx, label in enumerate(reversed(labels)):
            reverse_idx = len(labels) - 1 - idx
            current_angle = calculate_angle_degrees(
                center_x,
                center_y,
                label["start_x"],
                label["start_y"],
                label["middle"],
                start_angle,
                end_angle,
                total_length,
                x_radius,
                y_radius,
                normalize=False,
            )

            current_score = calculate_angle_of_three_points(
                label["feature_middle_x"], label["feature_middle_y"], 0.0, 0.0, label["start_x"], label["start_y"]
            )

            if idx == 0:
                overlaps_prev = check_overlap(label, labels[reverse_idx - 1], total_length)
                overlaps_next = check_overlap(label, labels[0], total_length)
            elif 0 < idx < len(labels) - 1:
                overlaps_prev = check_overlap(labels[reverse_idx - 1], label, total_length)
                overlaps_next = check_overlap(label, labels[reverse_idx + 1], total_length)
            else:
                overlaps_prev = check_overlap(label, labels[-1], total_length)
                overlaps_next = check_overlap(label, labels[reverse_idx + 1], total_length)

            if overlaps_prev and overlaps_next:
                continue

            if overlaps_prev:
                direction = 1
            elif overlaps_next:
                direction = -1
            else:
                test_angle_plus = current_angle + 1.0
                test_x_plus, test_y_plus = move_label(label, test_angle_plus)
                score_plus = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0.0, 0.0, test_x_plus, test_y_plus
                )

                test_angle_minus = current_angle - 1.0
                test_x_minus, test_y_minus = move_label(label, test_angle_minus)
                score_minus = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0.0, 0.0, test_x_minus, test_y_minus
                )
                direction = 1 if abs(score_plus) < abs(score_minus) else -1

            while True:
                new_angle = current_angle + direction * 0.05
                new_x, new_y = move_label(label, new_angle)
                new_score = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0.0, 0.0, new_x, new_y
                )

                label_copy = label.copy()
                label_copy["start_x"], label_copy["start_y"] = new_x, new_y
                if idx == 0:
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, labels[0], total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_idx - 1], total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_idx - 1], total_length) or check_overlap(
                            label_copy, labels[0], total_length
                        )
                elif 0 < idx < len(labels) - 1:
                    prev_label = labels[reverse_idx - 1]
                    next_label = labels[reverse_idx + 1]
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, next_label, total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, prev_label, total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, prev_label, total_length) or check_overlap(
                            label_copy, next_label, total_length
                        )
                else:
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_idx + 1], total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, labels[-1], total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, labels[-1], total_length) or check_overlap(
                            label_copy, labels[reverse_idx + 1], total_length
                        )

                if (abs(new_score) <= abs(current_score)) and not creates_new_overlap:
                    label["start_x"], label["start_y"] = new_x, new_y
                    current_angle = new_angle
                    current_score = new_score
                    changes_made = True
                else:
                    break
        if not changes_made:
            break

    return labels


def improved_label_placement_fc(
    labels,
    center_x,
    center_y,
    x_radius,
    y_radius,
    feature_radius,
    total_length,
    start_angle,
    end_angle,
    y_margin=0.3,
    max_angle_shift_deg=25.0,
    max_iterations=60,
):
    labels = sort_labels(labels)
    if not labels:
        return []

    # Preserve original sorted order; optimization may rotate labels to avoid
    # placing a dense cluster across the 0/360-degree seam.
    for idx, label in enumerate(labels):
        label["_opt_original_idx"] = idx

    if len(labels) > 2:
        raw_targets = [angle_from_middle_unwrapped(label["middle"], total_length) for label in labels]
        gap_values = [raw_targets[i + 1] - raw_targets[i] for i in range(len(labels) - 1)]
        gap_values.append(raw_targets[0] + 360.0 - raw_targets[-1])
        largest_gap_idx = max(range(len(gap_values)), key=gap_values.__getitem__)
        if largest_gap_idx < len(labels) - 1:
            rotate_by = largest_gap_idx + 1
            labels = labels[rotate_by:] + labels[:rotate_by]

    # Reassign targets as a monotonic unwrapped sequence in optimization order.
    prev_target_unwrapped: float | None = None
    for label in labels:
        raw_target = angle_from_middle_unwrapped(label["middle"], total_length)
        if prev_target_unwrapped is None:
            target_angle_unwrapped = raw_target
        else:
            target_angle_unwrapped = normalize_angle_near_reference(raw_target, prev_target_unwrapped)
            while target_angle_unwrapped <= prev_target_unwrapped:
                target_angle_unwrapped += 360.0
        label["target_angle_unwrapped"] = target_angle_unwrapped
        label["target_angle"] = target_angle_unwrapped % 360.0
        prev_target_unwrapped = target_angle_unwrapped

    pair_margins: list[list[float]] = []

    def move_label(label, angle):
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def check_overlap(label1, label2, total_length, effective_margin: float):
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
        )

    def pair_margin(idx1: int, idx2: int) -> float:
        if idx1 < idx2:
            return pair_margins[idx1][idx2]
        return pair_margins[idx2][idx1]

    def overlap_count(label_idx, candidate_label=None, *, full_scan: bool = False):
        current = candidate_label if candidate_label is not None else labels[label_idx]
        n_labels = len(labels)
        if n_labels <= 1:
            return 0
        checked_indices: set[int] = set()
        if full_scan:
            checked_indices = set(range(n_labels))
            checked_indices.discard(label_idx)
        else:
            # In angular-sorted order, most overlaps are local. Start with nearby labels.
            neighbor_span = 6
            for delta in range(1, neighbor_span + 1):
                prev_idx = (label_idx - delta) % n_labels
                next_idx = (label_idx + delta) % n_labels
                checked_indices.add(prev_idx)
                checked_indices.add(next_idx)
        overlaps = 0
        for other_idx in checked_indices:
            if other_idx == label_idx:
                continue
            if check_overlap(current, labels[other_idx], total_length, pair_margin(label_idx, other_idx)):
                overlaps += 1
        return overlaps

    def update_label_coordinates(label: dict, angle_unwrapped: float) -> None:
        label["angle_unwrapped"] = angle_unwrapped
        label["start_x"], label["start_y"] = move_label(label, angle_unwrapped)

    def local_density_weight(label_idx: int) -> float:
        dense_threshold_deg = 8.0
        n_labels = len(labels)
        target = labels[label_idx]["target_angle_unwrapped"]
        close_neighbors = 0
        if label_idx > 0:
            prev_target = labels[label_idx - 1]["target_angle_unwrapped"]
            if (target - prev_target) <= dense_threshold_deg:
                close_neighbors += 1
        if label_idx < n_labels - 1:
            next_target = labels[label_idx + 1]["target_angle_unwrapped"]
            if (next_target - target) <= dense_threshold_deg:
                close_neighbors += 1
        return 1.0 + close_neighbors

    def is_dense_center(label_idx: int) -> bool:
        if label_idx <= 0 or label_idx >= len(labels) - 1:
            return False
        dense_threshold_deg = 8.0
        left_gap = labels[label_idx]["target_angle_unwrapped"] - labels[label_idx - 1]["target_angle_unwrapped"]
        right_gap = labels[label_idx + 1]["target_angle_unwrapped"] - labels[label_idx]["target_angle_unwrapped"]
        return left_gap <= dense_threshold_deg and right_gap <= dense_threshold_deg

    def middle_proximity_penalty(label_idx: int, candidate_angle: float) -> float:
        """
        Penalize configurations where edge labels in a dense triplet are closer to
        their features than the center label.
        """
        target = labels[label_idx]["target_angle_unwrapped"]
        edge_delta = abs(candidate_angle - target)
        penalty = 0.0

        center_right = label_idx + 1
        if center_right < len(labels) - 1 and is_dense_center(center_right):
            center_delta = abs(
                labels[center_right]["angle_unwrapped"] - labels[center_right]["target_angle_unwrapped"]
            )
            if edge_delta < center_delta:
                penalty += center_delta - edge_delta

        center_left = label_idx - 1
        if center_left > 0 and is_dense_center(center_left):
            center_delta = abs(labels[center_left]["angle_unwrapped"] - labels[center_left]["target_angle_unwrapped"])
            if edge_delta < center_delta:
                penalty += center_delta - edge_delta

        return penalty

    def enforce_order(min_gap: float) -> bool:
        changed = False
        for i in range(1, len(labels)):
            minimum_allowed = labels[i - 1]["angle_unwrapped"] + min_gap
            if labels[i]["angle_unwrapped"] < minimum_allowed:
                update_label_coordinates(labels[i], minimum_allowed)
                changed = True
        return changed

    def count_remaining_overlaps() -> int:
        overlap_count_all = 0
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if check_overlap(labels[i], labels[j], total_length, pair_margin(i, j)):
                    overlap_count_all += 1
        return overlap_count_all

    def count_overlaps_for_label_list(label_list: list[dict], overlap_cap: int | None = None) -> int:
        overlap_count_all = 0
        for i in range(len(label_list)):
            for j in range(i + 1, len(label_list)):
                if check_overlap(label_list[i], label_list[j], total_length, pair_margin(i, j)):
                    overlap_count_all += 1
                    if overlap_cap is not None and overlap_count_all > overlap_cap:
                        return overlap_count_all
        return overlap_count_all

    def find_overlap_components() -> list[list[int]]:
        n_labels = len(labels)
        adjacency: dict[int, set[int]] = {i: set() for i in range(n_labels)}
        for i in range(n_labels):
            for j in range(i + 1, n_labels):
                if check_overlap(labels[i], labels[j], total_length, pair_margin(i, j)):
                    adjacency[i].add(j)
                    adjacency[j].add(i)

        components: list[list[int]] = []
        visited: set[int] = set()
        for idx in range(n_labels):
            if idx in visited or not adjacency[idx]:
                continue
            stack = [idx]
            component: list[int] = []
            visited.add(idx)
            while stack:
                current = stack.pop()
                component.append(current)
                for neigh in adjacency[current]:
                    if neigh not in visited:
                        visited.add(neigh)
                        stack.append(neigh)
            if component:
                components.append(sorted(component))
        return components

    def merge_close_components(components: list[list[int]], max_gap: int = 2) -> list[list[int]]:
        if not components:
            return []
        sorted_components = sorted((sorted(comp) for comp in components), key=lambda comp: comp[0])
        merged: list[list[int]] = [sorted_components[0]]
        for comp in sorted_components[1:]:
            if comp[0] - merged[-1][-1] <= max_gap:
                merged[-1] = sorted(set(merged[-1] + comp))
            else:
                merged.append(comp)
        return merged

    def relax_overlapping_clusters(shift_limit_deg: float) -> bool:
        """
        When local single-label optimization stalls, nudge clusters as a group.
        This allows neighboring labels to move together while keeping the
        cluster center label relatively close to its feature.
        """
        if len(labels) > HEAVY_CLUSTER_RELAX_MAX_LABELS:
            return False

        components = merge_close_components(find_overlap_components(), max_gap=2)
        if not components:
            return False

        cluster_step_limit = min(int(shift_limit_deg / angle_step), 40)
        n_labels = len(labels)
        changed = False

        for component in components:
            component_center = component[len(component) // 2]
            cluster_padding = min(6, max(2, len(component) + 1))
            cluster_start = max(0, component[0] - cluster_padding)
            cluster_end = min(n_labels - 1, component[-1] + cluster_padding)
            cluster_indices = list(range(cluster_start, cluster_end + 1))

            current_score = (
                count_remaining_overlaps(),
                abs(labels[component_center]["angle_unwrapped"] - labels[component_center]["target_angle_unwrapped"]),
                sum(abs(labels[idx]["angle_unwrapped"] - labels[idx]["target_angle_unwrapped"]) for idx in cluster_indices),
            )
            best_score = current_score
            best_ordered: list[float] | None = None

            center_candidates = [component_center]
            if len(component) >= 2:
                center_candidates.append(component[0])
                center_candidates.append(component[-1])
            center_candidates = sorted(set(center_candidates))
            center_shift_limit_steps = min(int(shift_limit_deg / angle_step), 48)
            center_shift_steps = range(-center_shift_limit_steps, center_shift_limit_steps + 1, 4)
            center_delta_cap = min(float(shift_limit_deg), RELAX_CENTER_DELTA_CAP_DEG)

            for center_idx in center_candidates:
                for step in range(0, cluster_step_limit + 1):
                    spread = step * angle_step
                    for shift_step in center_shift_steps:
                        center_shift = shift_step * angle_step
                        ordered = [float(labels[i]["angle_unwrapped"]) for i in range(n_labels)]

                        for idx in cluster_indices:
                            rel = idx - center_idx
                            distance_weight = 1.0 + 0.25 * max(0, abs(rel) - 1)
                            target = float(labels[idx]["target_angle_unwrapped"])
                            lower_bound = target - shift_limit_deg
                            upper_bound = target + shift_limit_deg
                            candidate = target + center_shift + (rel * spread * distance_weight)
                            if candidate < lower_bound:
                                candidate = lower_bound
                            elif candidate > upper_bound:
                                candidate = upper_bound
                            ordered[idx] = candidate

                        # Keep global ordering constraints valid.
                        for i in range(1, n_labels):
                            minimum_allowed = ordered[i - 1] + min_order_gap
                            if ordered[i] < minimum_allowed:
                                ordered[i] = minimum_allowed
                        for i in range(n_labels - 2, -1, -1):
                            maximum_allowed = ordered[i + 1] - min_order_gap
                            if ordered[i] > maximum_allowed:
                                ordered[i] = maximum_allowed

                        # Respect per-label shift budget.
                        violates_budget = False
                        for idx in range(n_labels):
                            target = float(labels[idx]["target_angle_unwrapped"])
                            if abs(ordered[idx] - target) > shift_limit_deg + 1e-6:
                                violates_budget = True
                                break
                        if violates_budget:
                            continue

                        temp_labels = [label.copy() for label in labels]
                        for idx in range(n_labels):
                            temp_labels[idx]["angle_unwrapped"] = ordered[idx]
                            temp_labels[idx]["start_x"], temp_labels[idx]["start_y"] = move_label(
                                temp_labels[idx], ordered[idx]
                            )

                        center_delta = abs(
                            temp_labels[component_center]["angle_unwrapped"]
                            - temp_labels[component_center]["target_angle_unwrapped"]
                        )
                        if center_delta > center_delta_cap:
                            continue

                        overlap_total = count_overlaps_for_label_list(
                            temp_labels,
                            overlap_cap=best_score[0],
                        )
                        if overlap_total > best_score[0]:
                            continue

                        candidate_score = (
                            overlap_total,
                            center_delta,
                            sum(
                                abs(
                                    temp_labels[idx]["angle_unwrapped"]
                                    - temp_labels[idx]["target_angle_unwrapped"]
                                )
                                for idx in cluster_indices
                            ),
                        )
                        if candidate_score < best_score:
                            best_score = candidate_score
                            best_ordered = ordered
                        if best_score[0] == 0:
                            break
                    if best_score[0] == 0:
                        break
                if best_score[0] == 0:
                    break

            if best_ordered is not None and best_score < current_score:
                for idx in range(n_labels):
                    update_label_coordinates(labels[idx], best_ordered[idx])
                changed = True

        return changed

    def optimize_with_shift_limit(shift_limit_deg: float, *, full_scan_checks: bool) -> None:
        max_steps = int(shift_limit_deg / angle_step)
        for _ in range(max_iterations):
            changes_made = False
            if enforce_order(min_order_gap):
                changes_made = True

            iteration_order = sorted(range(len(labels)), key=local_density_weight)
            for label_idx in iteration_order:
                label = labels[label_idx]
                current_angle = label["angle_unwrapped"]
                target_angle = label["target_angle_unwrapped"]
                current_overlaps = overlap_count(label_idx, full_scan=full_scan_checks)
                current_target_delta = abs(current_angle - target_angle)
                current_middle_penalty = middle_proximity_penalty(label_idx, current_angle)
                current_leader_length = _leader_length_px(label, total_length)

                if current_overlaps == 0 and current_middle_penalty == 0 and current_target_delta == 0:
                    continue

                lower_bound = target_angle - shift_limit_deg
                upper_bound = target_angle + shift_limit_deg
                if label_idx > 0:
                    lower_bound = max(lower_bound, labels[label_idx - 1]["angle_unwrapped"] + min_order_gap)
                if label_idx < len(labels) - 1:
                    upper_bound = min(upper_bound, labels[label_idx + 1]["angle_unwrapped"] - min_order_gap)
                if lower_bound > upper_bound:
                    lower_bound = labels[label_idx - 1]["angle_unwrapped"] + min_order_gap if label_idx > 0 else lower_bound
                    upper_bound = labels[label_idx + 1]["angle_unwrapped"] - min_order_gap if label_idx < len(labels) - 1 else upper_bound
                    if lower_bound > upper_bound:
                        continue

                anchor_weight = local_density_weight(label_idx)
                best_angle = current_angle
                best_x = label["start_x"]
                best_y = label["start_y"]
                best_overlaps = current_overlaps
                best_target_delta = current_target_delta
                best_move = 0.0
                best_middle_penalty = current_middle_penalty
                best_leader_length = current_leader_length
                proximity_penalty_factor = 4.0
                initial_score = (
                    current_overlaps,
                    anchor_weight * current_target_delta + proximity_penalty_factor * current_middle_penalty,
                    current_leader_length,
                    0.0,
                )
                best_score = initial_score

                search_complete = False
                for step in range(max_steps + 1):
                    if step == 0:
                        offsets = [0.0]
                    else:
                        offset = step * angle_step
                        offsets = [offset, -offset]

                    for offset in offsets:
                        candidate_angle = target_angle + offset
                        if not (lower_bound <= candidate_angle <= upper_bound):
                            continue
                        candidate_x, candidate_y = move_label(label, candidate_angle)
                        label_copy = label.copy()
                        label_copy["start_x"] = candidate_x
                        label_copy["start_y"] = candidate_y

                        candidate_overlaps = overlap_count(label_idx, label_copy, full_scan=full_scan_checks)
                        candidate_target_delta = abs(candidate_angle - target_angle)
                        candidate_move = abs(candidate_angle - current_angle)
                        candidate_middle_penalty = middle_proximity_penalty(label_idx, candidate_angle)
                        candidate_leader_length = _leader_length_px(label_copy, total_length)
                        candidate_score = (
                            candidate_overlaps,
                            anchor_weight * candidate_target_delta + proximity_penalty_factor * candidate_middle_penalty,
                            candidate_leader_length,
                            candidate_move,
                        )

                        if candidate_score < best_score:
                            best_angle = candidate_angle
                            best_x = candidate_x
                            best_y = candidate_y
                            best_overlaps = candidate_overlaps
                            best_target_delta = candidate_target_delta
                            best_move = candidate_move
                            best_middle_penalty = candidate_middle_penalty
                            best_leader_length = candidate_leader_length
                            best_score = candidate_score

                        if (
                            best_overlaps == 0
                            and best_target_delta == 0
                            and best_middle_penalty == 0
                            and best_leader_length == 0
                        ):
                            search_complete = True
                            break
                    if search_complete:
                        break

                if best_score < initial_score:
                    label["start_x"] = best_x
                    label["start_y"] = best_y
                    label["angle_unwrapped"] = best_angle
                    changes_made = True
            if not changes_made:
                break
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            for _ in range(3):
                if not relax_overlapping_clusters(shift_limit_deg):
                    break
                if count_remaining_overlaps() == 0:
                    break

    def separate_remaining_adjacent_overlaps(shift_limit_deg: float) -> bool:
        changed = False
        n_labels = len(labels)
        if n_labels < 2:
            return False

        for _ in range(160):
            pass_changed = False
            for i in range(n_labels - 1):
                if not check_overlap(labels[i], labels[i + 1], total_length, pair_margin(i, i + 1)):
                    continue

                pair_changed = False
                for _ in range(80):
                    can_shift_left = True
                    for k in range(0, i + 1):
                        target = float(labels[k]["target_angle_unwrapped"])
                        if labels[k]["angle_unwrapped"] - angle_step < target - shift_limit_deg - 1e-6:
                            can_shift_left = False
                            break

                    can_shift_right = True
                    for k in range(i + 1, n_labels):
                        target = float(labels[k]["target_angle_unwrapped"])
                        if labels[k]["angle_unwrapped"] + angle_step > target + shift_limit_deg + 1e-6:
                            can_shift_right = False
                            break

                    if not can_shift_left and not can_shift_right:
                        break

                    if can_shift_left:
                        for k in range(0, i + 1):
                            update_label_coordinates(labels[k], labels[k]["angle_unwrapped"] - angle_step)
                    if can_shift_right:
                        for k in range(i + 1, n_labels):
                            update_label_coordinates(labels[k], labels[k]["angle_unwrapped"] + angle_step)

                    enforce_order(min_order_gap)
                    if not check_overlap(labels[i], labels[i + 1], total_length, pair_margin(i, i + 1)):
                        pair_changed = True
                        break

                if pair_changed:
                    pass_changed = True
                    changed = True

            if not pass_changed:
                break

        return changed

    n_labels = len(labels)
    pair_margins = [[0.0 for _ in range(n_labels)] for _ in range(n_labels)]
    base_margin = max(y_margin, IMPROVED_LABELS_BASE_MARGIN_PX)
    for i in range(n_labels):
        for j in range(i + 1, n_labels):
            pair_margins[i][j] = minimum_bbox_gap_px(labels[i], labels[j], base_margin_px=base_margin)

    angle_step = 0.25
    min_order_gap = 0.05

    def snapshot_angles() -> list[float]:
        return [float(label["angle_unwrapped"]) for label in labels]

    def restore_angles(angle_snapshot: list[float]) -> None:
        for idx, angle in enumerate(angle_snapshot):
            update_label_coordinates(labels[idx], angle)

    def reset_angles_to_targets() -> None:
        for label in labels:
            update_label_coordinates(label, float(label["target_angle_unwrapped"]))
        enforce_order(min_order_gap)

    def placement_score() -> tuple[int, int, float]:
        mismatch_count, _ = _hemisphere_mismatch_metrics(labels, total_length)
        return (
            count_remaining_overlaps(),
            mismatch_count,
            sum(abs(float(label["angle_unwrapped"]) - float(label["target_angle_unwrapped"])) for label in labels),
        )

    for label in labels:
        target_angle_unwrapped = float(label["target_angle_unwrapped"])

        current_angle_mod = calculate_angle_degrees(
            center_x,
            center_y,
            label["start_x"],
            label["start_y"],
            label["middle"],
            start_angle,
            end_angle,
            total_length,
            x_radius,
            y_radius,
            normalize=False,
        )
        angle_unwrapped = normalize_angle_near_reference(current_angle_mod, target_angle_unwrapped)
        update_label_coordinates(label, angle_unwrapped)

    enforce_order(min_order_gap)

    base_shift_limit = min(float(max_angle_shift_deg), MAX_EXPANDED_SHIFT_DEG)
    shift_schedule = [base_shift_limit]
    if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS:
        for extra_shift in (35.0, 45.0, 60.0):
            if extra_shift > shift_schedule[-1]:
                shift_schedule.append(extra_shift)

    for phase_idx, shift_limit in enumerate(shift_schedule):
        run_full_scan = (
            len(labels) <= FULL_SCAN_LABEL_LIMIT
            and len(shift_schedule) > 1
            and phase_idx == (len(shift_schedule) - 1)
        )
        optimize_with_shift_limit(shift_limit, full_scan_checks=run_full_scan)
        if count_remaining_overlaps() == 0:
            break

    if count_remaining_overlaps() > 0:
        for extra_shift in (75.0, MAX_EXPANDED_SHIFT_DEG):
            if extra_shift <= shift_schedule[-1]:
                continue
            optimize_with_shift_limit(
                extra_shift,
                full_scan_checks=(len(labels) <= FULL_SCAN_LABEL_LIMIT),
            )
            shift_schedule.append(extra_shift)
            if count_remaining_overlaps() == 0:
                break

    if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
        separate_remaining_adjacent_overlaps(max(shift_schedule))

    if count_remaining_overlaps() > 0:
        baseline_angles = snapshot_angles()
        baseline_score = placement_score()
        optimize_with_shift_limit(
            max(shift_schedule),
            full_scan_checks=(len(labels) <= FULL_SCAN_LABEL_LIMIT),
        )
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            separate_remaining_adjacent_overlaps(max(shift_schedule))
        continued_score = placement_score()
        if continued_score > baseline_score:
            restore_angles(baseline_angles)
            continued_score = baseline_score

    if count_remaining_overlaps() > 0:
        baseline_angles = snapshot_angles()
        baseline_score = placement_score()
        reset_angles_to_targets()
        optimize_with_shift_limit(
            max(shift_schedule),
            full_scan_checks=(len(labels) <= FULL_SCAN_LABEL_LIMIT),
        )
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            separate_remaining_adjacent_overlaps(max(shift_schedule))
        retry_score = placement_score()
        if retry_score > baseline_score:
            restore_angles(baseline_angles)

    enforce_order(min_order_gap)
    labels = _refine_labels_to_preferred_hemisphere(
        labels,
        total_length,
        center_x,
        center_y,
        x_radius,
        y_radius,
    )
    labels = _polish_total_leader_distance(
        labels,
        total_length,
        center_x,
        center_y,
        x_radius,
        y_radius,
        min_order_gap=min_order_gap,
        angle_step_deg=angle_step,
        max_passes=8,
    )
    labels.sort(key=lambda x: x["_opt_original_idx"])
    for label in labels:
        label.pop("_opt_original_idx", None)
    return labels


def _derive_monotonic_unwrapped_angles(labels: list[dict], total_length: int) -> list[float]:
    """Derive monotonically increasing unwrapped angles from current label coordinates."""
    unwrapped_angles: list[float] = []
    prev_unwrapped: float | None = None

    for label in labels:
        if "angle_unwrapped" in label:
            angle_raw = float(label["angle_unwrapped"])
        else:
            angle_raw = math.degrees(math.atan2(float(label["start_y"]), float(label["start_x"])))

        target_unwrapped = float(
            label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
        )
        reference = target_unwrapped if prev_unwrapped is None else prev_unwrapped
        angle_unwrapped = normalize_angle_near_reference(angle_raw, reference)
        if prev_unwrapped is not None:
            while angle_unwrapped <= prev_unwrapped:
                angle_unwrapped += 360.0
        unwrapped_angles.append(angle_unwrapped)
        prev_unwrapped = angle_unwrapped

    return unwrapped_angles


def _refine_labels_to_preferred_hemisphere(
    labels: list[dict],
    total_length: int,
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    *,
    axis_neutral_deg: float = HEMISPHERE_AXIS_NEUTRAL_DEG,
    axis_eps_px: float = HEMISPHERE_AXIS_EPS_PX,
    max_shift_deg: float = HEMISPHERE_REFINE_MAX_SHIFT_DEG,
    angle_step_deg: float = HEMISPHERE_REFINE_STEP_DEG,
) -> list[dict]:
    """
    Refine mismatched labels toward their preferred hemisphere without increasing plain overlaps.

    This is a lightweight post-pass and only touches labels that are clearly on the
    opposite side (outside the neutral top/bottom band).
    """
    if not labels:
        return labels

    min_order_gap = 0.05
    max_steps = int(max_shift_deg / angle_step_deg)
    block_max_steps = int(min(max_shift_deg, 6.0) / angle_step_deg)
    base_plain_overlaps = _count_label_overlaps(labels, total_length, use_min_gap=False)
    current_unwrapped_angles = _derive_monotonic_unwrapped_angles(labels, total_length)

    def plain_overlap_count_for_candidate(label_idx: int, candidate_label: dict) -> int:
        overlap_count = 0
        for other_idx, other in enumerate(labels):
            if other_idx == label_idx:
                continue
            if y_overlap(candidate_label, other, total_length, 0.1) and x_overlap(
                candidate_label, other, minimum_margin=1.0
            ):
                overlap_count += 1
        return overlap_count

    pass_orders = [range(len(labels) - 1, -1, -1), range(len(labels))]
    for pass_order in pass_orders:
        for label_idx in pass_order:
            label = labels[label_idx]
            preferred_half = _preferred_half_for_label(label, total_length, axis_neutral_deg=axis_neutral_deg)
            if preferred_half == 0:
                continue

            current_half = _current_half_from_x(float(label["start_x"]), axis_eps_px=axis_eps_px)
            if current_half in (0, preferred_half):
                continue

            lower_bound = current_unwrapped_angles[label_idx - 1] + min_order_gap if label_idx > 0 else -float("inf")
            upper_bound = (
                current_unwrapped_angles[label_idx + 1] - min_order_gap
                if label_idx < len(labels) - 1
                else float("inf")
            )
            if lower_bound >= upper_bound:
                continue

            target_angle_unwrapped = float(
                label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
            )
            target_angle_mod = target_angle_unwrapped % 360.0
            current_overlap_count = plain_overlap_count_for_candidate(label_idx, label)
            best_candidate: tuple[dict, float, int] | None = None
            best_score: tuple[int, float] | None = None

            for step in range(max_steps + 1):
                if step == 0:
                    offsets = [0.0]
                else:
                    step_offset = step * angle_step_deg
                    offsets = [step_offset, -step_offset]

                for offset in offsets:
                    candidate_angle_mod = (target_angle_mod + offset) % 360.0
                    candidate_angle_unwrapped = normalize_angle_near_reference(
                        candidate_angle_mod,
                        current_unwrapped_angles[label_idx],
                    )
                    while candidate_angle_unwrapped <= lower_bound:
                        candidate_angle_unwrapped += 360.0
                    while candidate_angle_unwrapped >= upper_bound and (candidate_angle_unwrapped - 360.0) > lower_bound:
                        candidate_angle_unwrapped -= 360.0
                    if not (lower_bound < candidate_angle_unwrapped < upper_bound):
                        continue

                    candidate_x, candidate_y = calculate_coordinates(
                        center_x,
                        center_y,
                        x_radius,
                        y_radius,
                        candidate_angle_unwrapped,
                        label["middle"],
                        total_length,
                    )
                    if _current_half_from_x(candidate_x, axis_eps_px=axis_eps_px) != preferred_half:
                        continue

                    candidate_label = label.copy()
                    candidate_label["start_x"] = candidate_x
                    candidate_label["start_y"] = candidate_y
                    candidate_overlap_count = plain_overlap_count_for_candidate(label_idx, candidate_label)
                    candidate_plain_total = base_plain_overlaps - current_overlap_count + candidate_overlap_count
                    if candidate_plain_total > base_plain_overlaps:
                        continue

                    candidate_score = (candidate_plain_total, abs(offset))
                    if best_score is None or candidate_score < best_score:
                        best_score = candidate_score
                        best_candidate = (candidate_label, candidate_angle_unwrapped, candidate_plain_total)

                if best_score is not None and best_score[0] <= base_plain_overlaps:
                    break

            if best_candidate is None:
                # Try a small coordinated block shift when single-label movement
                # cannot satisfy side preference without increasing overlaps.
                if len(labels) <= FULL_SCAN_LABEL_LIMIT and block_max_steps > 0:
                    # Keep this block search intentionally small: it is only a
                    # fallback for near-axis labels where tiny coordinated moves
                    # can fix side preference without overlap regressions.
                    target_distance_from_axis = _vertical_axis_distance_deg(target_angle_mod)
                    if target_distance_from_axis > (axis_neutral_deg + 4.0):
                        continue

                    raw_block_candidates = [
                        (label_idx - 1, label_idx),
                        (label_idx, label_idx + 1),
                        (label_idx - 2, label_idx),
                        (label_idx, label_idx + 2),
                        (label_idx - 1, label_idx + 1),
                    ]
                    block_candidates: list[tuple[int, int]] = []
                    for block_start_raw, block_end_raw in raw_block_candidates:
                        block_start = max(0, block_start_raw)
                        block_end = min(len(labels) - 1, block_end_raw)
                        if block_start >= block_end:
                            continue
                        if (block_start, block_end) not in block_candidates:
                            block_candidates.append((block_start, block_end))

                    best_block_score: tuple[int, float, int] | None = None
                    best_block_result: tuple[list[dict], list[float], int] | None = None
                    for block_start, block_end in block_candidates:
                        for step in range(1, block_max_steps + 1):
                            step_offset = step * angle_step_deg
                            for delta in (step_offset, -step_offset):
                                candidate_angles = current_unwrapped_angles.copy()
                                for move_idx in range(block_start, block_end + 1):
                                    candidate_angles[move_idx] += delta

                                order_valid = True
                                for check_idx in range(1, len(candidate_angles)):
                                    if candidate_angles[check_idx] <= candidate_angles[check_idx - 1] + min_order_gap:
                                        order_valid = False
                                        break
                                if not order_valid:
                                    continue

                                candidate_labels = [item.copy() for item in labels]
                                for move_idx in range(block_start, block_end + 1):
                                    moved_label = candidate_labels[move_idx]
                                    moved_x, moved_y = calculate_coordinates(
                                        center_x,
                                        center_y,
                                        x_radius,
                                        y_radius,
                                        candidate_angles[move_idx],
                                        moved_label["middle"],
                                        total_length,
                                    )
                                    moved_label["start_x"] = moved_x
                                    moved_label["start_y"] = moved_y
                                    if "angle_unwrapped" in moved_label:
                                        moved_label["angle_unwrapped"] = candidate_angles[move_idx]

                                target_half = _current_half_from_x(
                                    float(candidate_labels[label_idx]["start_x"]),
                                    axis_eps_px=axis_eps_px,
                                )
                                if target_half != preferred_half:
                                    continue

                                candidate_plain_total = _count_label_overlaps(
                                    candidate_labels,
                                    total_length,
                                    use_min_gap=False,
                                )
                                if candidate_plain_total > base_plain_overlaps:
                                    continue

                                block_width = block_end - block_start
                                candidate_score = (candidate_plain_total, abs(delta), block_width)
                                if best_block_score is None or candidate_score < best_block_score:
                                    best_block_score = candidate_score
                                    best_block_result = (
                                        candidate_labels,
                                        candidate_angles,
                                        candidate_plain_total,
                                    )
                            if best_block_score is not None and best_block_score[0] <= base_plain_overlaps:
                                break
                        if best_block_score is not None and best_block_score[0] <= base_plain_overlaps:
                            break

                    if best_block_result is None:
                        continue

                    candidate_labels, candidate_angles, chosen_plain_total = best_block_result
                    labels[:] = candidate_labels
                    current_unwrapped_angles = candidate_angles
                    base_plain_overlaps = chosen_plain_total
                    continue

                continue

            candidate_label, chosen_angle_unwrapped, chosen_plain_total = best_candidate
            label["start_x"] = candidate_label["start_x"]
            label["start_y"] = candidate_label["start_y"]
            if "angle_unwrapped" in label:
                label["angle_unwrapped"] = chosen_angle_unwrapped
            current_unwrapped_angles[label_idx] = chosen_angle_unwrapped
            base_plain_overlaps = chosen_plain_total

    return labels


def _polish_total_leader_distance(
    labels: list[dict],
    total_length: int,
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    *,
    min_order_gap: float = 0.05,
    angle_step_deg: float = 0.25,
    max_passes: int = 8,
) -> list[dict]:
    """
    Reduce total leader length while preserving overlap/order/hemisphere quality.

    This is a conservative final pass that only accepts strictly better global
    scores and never regresses plain/min-gap overlaps or hemisphere mismatches.
    """
    if not labels:
        return labels

    def _target_angle_unwrapped(label: dict) -> float:
        return float(
            label.get(
                "target_angle_unwrapped",
                angle_from_middle_unwrapped(float(label["middle"]), total_length),
            )
        )

    def _update_label_coordinates(label: dict, angle_unwrapped: float) -> None:
        label["angle_unwrapped"] = angle_unwrapped
        label["start_x"], label["start_y"] = calculate_coordinates(
            center_x,
            center_y,
            x_radius,
            y_radius,
            angle_unwrapped,
            label["middle"],
            total_length,
        )

    def _total_target_delta() -> float:
        delta_sum = 0.0
        for label in labels:
            target = _target_angle_unwrapped(label)
            current = float(label.get("angle_unwrapped", target))
            delta_sum += abs(current - target)
        return delta_sum

    def _labels_overlap(label_a: dict, label_b: dict, *, use_min_gap: bool) -> bool:
        if use_min_gap:
            margin = minimum_bbox_gap_px(label_a, label_b, base_margin_px=0.0)
            y_margin = margin
            x_margin = margin
        else:
            y_margin = 0.1
            x_margin = 1.0
        return y_overlap(label_a, label_b, total_length, y_margin) and x_overlap(
            label_a,
            label_b,
            minimum_margin=x_margin,
        )

    def _label_overlap_count(label_idx: int, candidate_label: dict, *, use_min_gap: bool) -> int:
        overlap_count = 0
        for other_idx, other in enumerate(labels):
            if other_idx == label_idx:
                continue
            if _labels_overlap(candidate_label, other, use_min_gap=use_min_gap):
                overlap_count += 1
        return overlap_count

    unwrapped = _derive_monotonic_unwrapped_angles(labels, total_length)
    for idx, angle_unwrapped in enumerate(unwrapped):
        labels[idx]["angle_unwrapped"] = angle_unwrapped

    plain_overlaps = _count_label_overlaps(labels, total_length, use_min_gap=False)
    min_gap_overlaps = _count_label_overlaps(labels, total_length, use_min_gap=True)
    baseline_min_gap_overlaps = min_gap_overlaps
    allow_min_gap_relaxation = len(labels) >= DENSE_DISTANCE_COMPACTION_MIN_LABELS
    max_allowed_min_gap_overlaps = baseline_min_gap_overlaps + (
        DENSE_DISTANCE_COMPACTION_MIN_GAP_RELAX if allow_min_gap_relaxation else 0
    )
    hemisphere_mismatch_count, _ = _hemisphere_mismatch_metrics(labels, total_length)
    leader_sum, _ = _leader_distance_score(labels, total_length)
    target_delta_sum = _total_target_delta()
    best_score = (
        plain_overlaps,
        min_gap_overlaps,
        hemisphere_mismatch_count,
        leader_sum,
        target_delta_sum,
    )

    for _ in range(max_passes):
        pass_changed = False
        label_order = sorted(
            range(len(labels)),
            key=lambda label_idx: _leader_length_px(labels[label_idx], total_length),
            reverse=True,
        )

        for label_idx in label_order:
            label = labels[label_idx]
            current_angle = float(label["angle_unwrapped"])
            target_angle = _target_angle_unwrapped(label)
            if math.isclose(current_angle, target_angle, abs_tol=1e-9):
                continue

            lower_bound = (
                float(labels[label_idx - 1]["angle_unwrapped"]) + min_order_gap
                if label_idx > 0
                else -float("inf")
            )
            upper_bound = (
                float(labels[label_idx + 1]["angle_unwrapped"]) - min_order_gap
                if label_idx < (len(labels) - 1)
                else float("inf")
            )
            if lower_bound >= upper_bound:
                continue

            direction = 1.0 if target_angle > current_angle else -1.0
            max_steps = int(abs(target_angle - current_angle) / angle_step_deg)
            candidate_angles = [current_angle + direction * angle_step_deg * step for step in range(1, max_steps + 1)]
            if not candidate_angles or not math.isclose(candidate_angles[-1], target_angle, abs_tol=1e-9):
                candidate_angles.append(target_angle)

            current_plain_for_label = _label_overlap_count(label_idx, label, use_min_gap=False)
            current_min_gap_for_label = _label_overlap_count(label_idx, label, use_min_gap=True)
            current_leader = _leader_length_px(label, total_length)
            current_target_delta = abs(current_angle - target_angle)
            original_start_x = float(label["start_x"])
            original_start_y = float(label["start_y"])
            original_angle = current_angle

            best_strict: tuple[float, tuple[int, int, int, float, float]] | None = None
            best_relaxed: tuple[float, tuple[int, int, int, float, float], tuple[float, float, int]] | None = None
            for candidate_angle in candidate_angles:
                if candidate_angle <= lower_bound + 1e-9 or candidate_angle >= upper_bound - 1e-9:
                    continue

                candidate_label = label.copy()
                candidate_label["angle_unwrapped"] = candidate_angle
                candidate_label["start_x"], candidate_label["start_y"] = calculate_coordinates(
                    center_x,
                    center_y,
                    x_radius,
                    y_radius,
                    candidate_angle,
                    label["middle"],
                    total_length,
                )

                candidate_plain_for_label = _label_overlap_count(label_idx, candidate_label, use_min_gap=False)
                candidate_plain_total = plain_overlaps - current_plain_for_label + candidate_plain_for_label
                if candidate_plain_total > plain_overlaps:
                    continue

                candidate_min_gap_for_label = _label_overlap_count(label_idx, candidate_label, use_min_gap=True)
                candidate_min_gap_total = min_gap_overlaps - current_min_gap_for_label + candidate_min_gap_for_label
                if candidate_min_gap_total > max_allowed_min_gap_overlaps:
                    continue

                candidate_leader = _leader_length_px(candidate_label, total_length)
                candidate_leader_sum = leader_sum - current_leader + candidate_leader
                candidate_target_delta = abs(candidate_angle - target_angle)
                candidate_target_delta_sum = target_delta_sum - current_target_delta + candidate_target_delta

                label["start_x"] = candidate_label["start_x"]
                label["start_y"] = candidate_label["start_y"]
                label["angle_unwrapped"] = candidate_angle
                candidate_mismatch_count, _ = _hemisphere_mismatch_metrics(labels, total_length)
                label["start_x"] = original_start_x
                label["start_y"] = original_start_y
                label["angle_unwrapped"] = original_angle
                if candidate_mismatch_count > hemisphere_mismatch_count:
                    continue

                candidate_score = (
                    candidate_plain_total,
                    candidate_min_gap_total,
                    candidate_mismatch_count,
                    candidate_leader_sum,
                    candidate_target_delta_sum,
                )
                if candidate_min_gap_total <= min_gap_overlaps and candidate_score < best_score:
                    if best_strict is None or candidate_score < best_strict[1]:
                        best_strict = (candidate_angle, candidate_score)
                    continue
                if (
                    allow_min_gap_relaxation
                    and candidate_plain_total == plain_overlaps
                    and candidate_mismatch_count <= hemisphere_mismatch_count
                    and candidate_min_gap_total > min_gap_overlaps
                ):
                    relaxed_cmp = (
                        candidate_leader_sum,
                        candidate_target_delta_sum,
                        candidate_min_gap_total,
                    )
                    if best_relaxed is None or relaxed_cmp < best_relaxed[2]:
                        best_relaxed = (candidate_angle, candidate_score, relaxed_cmp)

            accepted = best_strict
            if accepted is None and best_relaxed is not None:
                relaxed_angle, relaxed_score, _ = best_relaxed
                current_cmp = (leader_sum, target_delta_sum, min_gap_overlaps)
                if (relaxed_score[3], relaxed_score[4], relaxed_score[1]) < current_cmp:
                    accepted = (relaxed_angle, relaxed_score)

            if accepted is None:
                continue

            accepted_angle, accepted_score = accepted
            _update_label_coordinates(label, accepted_angle)
            plain_overlaps = accepted_score[0]
            min_gap_overlaps = accepted_score[1]
            hemisphere_mismatch_count = accepted_score[2]
            leader_sum = accepted_score[3]
            target_delta_sum = accepted_score[4]
            best_score = accepted_score
            pass_changed = True

        if not pass_changed:
            break

    return labels


def rearrange_labels_fc(
    labels,
    feature_radius,
    total_length,
    genome_len,
    config_dict,
    strands,
    is_outer,
    arena_outer_radius: float | None = None,
    cfg: GbdrawConfig | None = None,
):
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    track_type = cfg.canvas.circular.track_type

    if is_outer:
        offset_config = cfg.labels.unified_adjustment.outer_labels
        x_radius_factor = cfg.labels.arc_x_radius_factor[track_type][strands][genome_len]
        y_radius_factor = cfg.labels.arc_y_radius_factor[track_type][strands][genome_len]
        default_center_x = cfg.labels.arc_center_x[track_type][genome_len]
    else:
        offset_config = cfg.labels.unified_adjustment.inner_labels
        x_radius_factor = cfg.labels.inner_arc_x_radius_factor[track_type][strands][genome_len]
        y_radius_factor = cfg.labels.inner_arc_y_radius_factor[track_type][strands][genome_len]
        default_center_x = cfg.labels.inner_arc_center_x[track_type][genome_len]

    x_radius_offset = offset_config.x_radius_offset
    y_radius_offset = offset_config.y_radius_offset

    if is_outer:
        x_radius = feature_radius * x_radius_factor * x_radius_offset
        y_radius = feature_radius * y_radius_factor * y_radius_offset
    else:
        x_radius = feature_radius * x_radius_factor * (2 - x_radius_offset)
        y_radius = feature_radius * y_radius_factor * (2 - y_radius_offset)

    # Optional "arena" override: scale the ellipse so it fits the desired outer radius (px).
    if arena_outer_radius is not None:
        arena_outer_radius = float(arena_outer_radius)
        if arena_outer_radius > 0:
            current_max = max(abs(x_radius), abs(y_radius))
            if current_max == 0:
                x_radius = arena_outer_radius
                y_radius = arena_outer_radius
            else:
                scale = arena_outer_radius / current_max
                x_radius *= scale
                y_radius *= scale

    center_y = 0
    center_x = default_center_x
    start_angle = 0
    end_angle = 360

    sorted_labels = sorted(labels, key=lambda x: x["middle"])
    if not sorted_labels:
        return []

    # Very dense sets are better handled by the legacy strategy: it keeps labels
    # closer to feature angles and avoids expensive combinatorial relaxation.
    if len(sorted_labels) >= LEGACY_PLACEMENT_LABEL_THRESHOLD:
        candidate_legacy = [label.copy() for label in sorted_labels]
        candidate_legacy = _legacy_place_labels_on_arc_fc(
            candidate_legacy,
            center_x,
            center_y,
            x_radius,
            y_radius,
            start_angle,
            end_angle,
            total_length,
        )
        candidate_legacy = _legacy_improved_label_placement_fc(
            candidate_legacy,
            center_x,
            center_y,
            x_radius,
            y_radius,
            feature_radius,
            total_length,
            start_angle,
            end_angle,
        )
        candidate_legacy = sort_labels(candidate_legacy)
        candidate_legacy = _refine_labels_to_preferred_hemisphere(
            candidate_legacy,
            total_length,
            center_x,
            center_y,
            x_radius,
            y_radius,
        )
        legacy_score = _placement_score(candidate_legacy, total_length)
        if legacy_score[0] == 0 and legacy_score[2] == 0:
            return _assign_leader_start_points(candidate_legacy, total_length)

    # Candidate A: current placement pipeline.
    candidate_current = [label.copy() for label in sorted_labels]
    candidate_current = place_labels_on_arc_fc(
        candidate_current, center_x, center_y, x_radius, y_radius, start_angle, end_angle, total_length
    )
    candidate_current = improved_label_placement_fc(
        candidate_current,
        center_x,
        center_y,
        x_radius,
        y_radius,
        feature_radius,
        total_length,
        start_angle,
        end_angle,
    )
    candidate_current = _refine_labels_to_preferred_hemisphere(
        candidate_current,
        total_length,
        center_x,
        center_y,
        x_radius,
        y_radius,
    )
    best_labels = candidate_current
    best_score = _placement_score(candidate_current, total_length)

    # Candidate B: legacy fallback for dense unresolved overlaps.
    if len(sorted_labels) > FULL_SCAN_LABEL_LIMIT and (best_score[0] > 0 or best_score[2] > 0):
        candidate_legacy = [label.copy() for label in sorted_labels]
        candidate_legacy = _legacy_place_labels_on_arc_fc(
            candidate_legacy,
            center_x,
            center_y,
            x_radius,
            y_radius,
            start_angle,
            end_angle,
            total_length,
        )
        candidate_legacy = _legacy_improved_label_placement_fc(
            candidate_legacy,
            center_x,
            center_y,
            x_radius,
            y_radius,
            feature_radius,
            total_length,
            start_angle,
            end_angle,
        )
        candidate_legacy = sort_labels(candidate_legacy)
        candidate_legacy = _refine_labels_to_preferred_hemisphere(
            candidate_legacy,
            total_length,
            center_x,
            center_y,
            x_radius,
            y_radius,
        )
        legacy_score = _placement_score(candidate_legacy, total_length)
        if legacy_score < best_score:
            best_labels = candidate_legacy

    return _assign_leader_start_points(best_labels, total_length)


def prepare_label_list(
    feature_dict,
    total_length,
    radius,
    track_ratio,
    config_dict,
    cfg: GbdrawConfig | None = None,
    *,
    outer_arena: tuple[float, float] | None = None,
):
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    embedded_labels = []
    outer_labels = []
    inner_labels = []
    label_filtering = cfg.labels.filtering.as_dict()
    label_filtering = preprocess_label_filtering(label_filtering)
    length_threshold = cfg.labels.length_threshold.circular
    length_param = determine_length_parameter(total_length, length_threshold)
    track_type = cfg.canvas.circular.track_type
    strandedness = cfg.canvas.strandedness

    strands = "separate" if strandedness else "single"
    allow_inner_labels = cfg.canvas.circular.allow_inner_labels
    radius_factor = cfg.labels.radius_factor[track_type][strands][length_param]
    inner_radius_factor = cfg.labels.inner_radius_factor[track_type][strands][length_param]

    font_family = cfg.objects.text.font_family
    font_size: float = cfg.labels.font_size.for_length_param(length_param)
    interval = cfg.canvas.dpi

    track_ratio_factor = cfg.canvas.circular.track_ratio_factors[length_param][0]
    cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)

    for feature_object in feature_dict.values():
        feature_label_text = get_label_text(feature_object, label_filtering)
        if feature_label_text == "":
            continue
        else:
            label_entry = dict()
            longest_segment_length = 0
            is_embedded = False
            label_middle = 0
            coordinate_strand: str = "undefined"
            feature_location_list = feature_object.location
            list_of_coordinates = feature_object.coordinates
            # `FeatureObject.strand` is derived at creation time; keep placement pure.

            feature_location_count = 0
            for coordinate in list_of_coordinates:
                if feature_location_list[feature_location_count].kind == "line":
                    feature_location_count += 1
                    continue
                else:
                    coordinate_start = int(coordinate.start)
                    coordinate_end = int(coordinate.end)
                    coordinate_strand = get_strand(coordinate.strand)
                    interval_length = abs(int(coordinate_end - coordinate_start) + 1)
                    interval_middle = int(coordinate_end + coordinate_start) / 2
                    feature_location_count += 1
                    if interval_length > longest_segment_length:
                        longest_segment_start = coordinate_start
                        longest_segment_end = coordinate_end
                        longeset_segment_middle = interval_middle
                        longest_segment_length = interval_length

            # Get track_id for overlap resolution
            track_id = getattr(feature_object, 'feature_track_id', 0)
            factors: list[float] = calculate_feature_position_factors_circular(
                total_length, coordinate_strand, track_ratio, cds_ratio, offset, track_type, strandedness, track_id
            )
            # Store track_id in label entry for embedded label drawing
            label_entry["track_id"] = track_id
            bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
            label_middle = longeset_segment_middle
            label_as_feature_length = total_length * (1.1 * bbox_width_px) / (2 * math.pi * radius)
            label_start = label_middle - (label_as_feature_length / 2)
            label_end = label_middle + (label_as_feature_length / 2)
            feature_middle_x: float = (radius * factors[1]) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
            feature_middle_y: float = (radius * factors[1]) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            feature_anchor_x: float = feature_middle_x
            feature_anchor_y: float = feature_middle_y
            is_outer_label = (feature_object.strand == "positive") or (allow_inner_labels is False)
            feature_outer_radius = radius * float(factors[2])
            if is_outer_label:
                if cfg.canvas.resolve_overlaps and not strandedness:
                    # Raised feature tracks look disconnected when leaders target the track center.
                    # Anchor to the outer edge so feature-to-label distance reads more naturally.
                    feature_anchor_x = feature_outer_radius * math.cos(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                    feature_anchor_y = feature_outer_radius * math.sin(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                if outer_arena is not None:
                    # Use the inner edge of the arena as the "anchor" radius for leader lines.
                    anchor_radius = float(outer_arena[0])
                    middle_radius = anchor_radius
                else:
                    middle_radius = radius_factor * radius
                    if cfg.canvas.resolve_overlaps and not strandedness:
                        middle_radius = max(
                            float(middle_radius),
                            float(feature_outer_radius + MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX),
                        )
                        label_entry["min_outer_start_radius_px"] = float(
                            feature_outer_radius + MIN_OUTER_LABEL_TEXT_CLEARANCE_PX
                        )
                middle_x = middle_radius * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y = middle_radius * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            else:
                middle_x = (inner_radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y = (inner_radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            if label_start > longest_segment_start and label_end < longest_segment_end:
                is_embedded = True
            else:
                is_embedded = False
            label_entry["label_text"] = feature_label_text
            label_entry["middle"] = label_middle
            label_entry["start"] = label_start
            label_entry["end"] = label_end
            label_entry["middle_x"] = middle_x
            label_entry["middle_y"] = middle_y
            label_entry["feature_middle_x"] = feature_middle_x
            label_entry["feature_middle_y"] = feature_middle_y
            label_entry["feature_anchor_x"] = feature_anchor_x
            label_entry["feature_anchor_y"] = feature_anchor_y
            label_entry["width_px"] = bbox_width_px
            label_entry["height_px"] = bbox_height_px
            label_entry["strand"] = coordinate_strand
            label_entry["is_embedded"] = is_embedded
            label_entry["font_size"] = font_size
            label_entry["font_family"] = font_family

            if is_embedded:
                embedded_labels.append(label_entry)
            else:
                if label_entry["middle"] > (total_length / 2):
                    if feature_object.strand == "positive" or allow_inner_labels is False:
                        label_entry["is_inner"] = False
                        outer_labels.append(label_entry)
                    else:
                        label_entry["is_inner"] = True
                        inner_labels.append(label_entry)

                else:
                    if feature_object.strand == "positive" or allow_inner_labels is False:
                        label_entry["is_inner"] = False
                        outer_labels.append(label_entry)
                    else:
                        label_entry["is_inner"] = True
                        inner_labels.append(label_entry)
    outer_labels_rearranged = rearrange_labels_fc(
        outer_labels,
        radius,
        total_length,
        length_param,
        config_dict,
        strands,
        is_outer=True,
        arena_outer_radius=(float(outer_arena[1]) if outer_arena is not None else None),
        cfg=cfg,
    )
    if cfg.canvas.resolve_overlaps and not strandedness:
        feature_outer_radius_intervals = _build_outer_feature_radius_intervals(
            feature_dict,
            total_length,
            radius,
            track_ratio,
            cfg,
        )
        outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
            outer_labels_rearranged,
            total_length,
            feature_outer_radius_intervals,
        )
        outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
        for _ in range(2):
            outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
                outer_labels_rearranged,
                total_length,
            )
            outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
        outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
            outer_labels_rearranged,
            total_length,
        )
        outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
            outer_labels_rearranged,
            total_length,
            feature_outer_radius_intervals,
        )
        outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
        # Final radial enforcement can reintroduce tight angular collisions.
        outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
            outer_labels_rearranged,
            total_length,
        )
    outer_labels_rearranged = assign_leader_start_points(outer_labels_rearranged, total_length)

    inner_labels_rearranged = []
    if allow_inner_labels:
        inner_labels_rearranged = rearrange_labels_fc(
            inner_labels,
            radius,
            total_length,
            length_param,
            config_dict,
            strands,
            is_outer=False,
            cfg=cfg,
        )

    label_list_fc = embedded_labels + outer_labels_rearranged + inner_labels_rearranged
    return label_list_fc


__all__ = [
    "minimum_bbox_gap_px",
    "calculate_angle_degrees",
    "calculate_angle_for_y",
    "calculate_coordinates",
    "angle_from_middle",
    "angle_from_middle_unwrapped",
    "angle_delta_deg",
    "normalize_angle_near_reference",
    "euclidean_distance",
    "improved_label_placement_fc",
    "place_labels_on_arc_fc",
    "assign_leader_start_points",
    "prepare_label_list",
    "rearrange_labels_fc",
    "sort_labels",
    "x_overlap",
    "y_overlap",
]


