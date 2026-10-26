#!/usr/bin/env python
# coding: utf-8

"""Circular and arc-based label placement utilities."""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from functools import lru_cache
from typing import Any, Literal, NamedTuple, TypeAlias, cast, get_args

from ..features.coordinates import get_strand
from .filtering import label_filter_maps
from .circular_candidates import build_circular_label_candidates
from .circular_leaders import chosen_contact_end, leader_polyline
from .circular_radial import place_radial_labels
from .circular_types import CircularLabel, CircularLabelCandidateCache
from .fit import FITS, LabelFit, LabelKind, LabelShape, Polygon, label_fits_by_kind
from .label_style import CIRCULAR_LABEL_GAP_EM
from ..config.models import (
    CircularRenderProfile,
    GbdrawConfig,
)
from ..core.text import drawn_text_box
from ..diagrams.circular.presets import CircularTrackPreset, circular_feature_lane_direction_for_preset
from ..auto_sizes import determine_length_parameter
from ..layout.common import calculate_cds_ratio
from ..layout.circular import (
    CircularFeatureLayout,
    CircularRadialLayout,
    calculate_feature_position_factors_circular,
    feature_radius_intervals,
)
from ..layout.spatial import (
    Aabb,
    AabbIndex,
    bounding_aabb,
    candidate_aabb_pairs,
    touching_aabb_indices,
    touching_aabb_pairs,
    union_aabbs,
)
from ..svg.arrows import (
    calculate_circular_arrow_length,
    resolve_circular_arrow_head_length_bp,
)
from ..layout.text_geometry import annular_sector_parts
from ..svg.ids import instance_svg_id

# Keep dense large-font labels from being pushed excessively far from features.
from gbdraw.layout.record_coordinates import RecordDisplayTransform

MIN_BBOX_GAP_RATIO = 0.01
MIN_BBOX_GAP_FLOOR_PX = 1.2
HEAVY_CLUSTER_RELAX_MAX_LABELS = 70
FULL_SCAN_LABEL_LIMIT = 70
MAX_EXPANDED_SHIFT_DEG = 90.0
RELAX_CENTER_DELTA_CAP_DEG = 35.0
LEGACY_PLACEMENT_LABEL_THRESHOLD = 50
# Smallest row gap of the dense legacy packer: the y margin of the overlap count.
LEGACY_MIN_ROW_GAP_PX = 0.1
DenseOuterRows: TypeAlias = Literal["spaced", "legacy"]


class CircularLabelProbeStopped(Exception):
    """The dense outer candidate overlaps, so this placement attempt cannot fit."""

    def __init__(self, overlap_pairs: int) -> None:
        super().__init__(overlap_pairs)
        self.overlap_pairs = int(overlap_pairs)
DENSE_INNER_RELAX_MIN_LABELS = 40
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
OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX = 3.0
RESOLVE_OVERLAP_MARGIN_EXTRA_PX = 4.2
LEADER_LABEL_COLLISION_MARGIN_PX = 0.8
LEADER_LABEL_SHIFT_STEP_DEG = 0.25
LEADER_LABEL_MAX_SHIFT_DEG = 45.0
LEADER_LABEL_MAX_PASSES = 6
LEADER_LABEL_MIN_ORDER_GAP_DEG = 0.05
# OV-1823: an outer label that no angular shift clears of a leader may step outward, 1 px at a time,
# by up to twice its own height (enough to clear the elbow of a neighbour's leader one row up).
LEADER_LABEL_RADIAL_STEP_PX = 1.0
LEADER_LABEL_MAX_RADIAL_SHIFT_HEIGHTS = 2.0


def _normalize_circular_track_preset(raw: object) -> str:
    preset = str(raw or "tuckin").strip().lower()
    return preset if preset in get_args(CircularTrackPreset) else "tuckin"


def _effective_outer_middle_anchor_clearance_px(
    *,
    resolve_overlaps: bool,
    strandedness: bool,
    lane_direction: str,
    feature_band_width_px: float,
) -> float:
    """Return the middle-anchor clearance used for outer label elbows."""
    base_clearance = float(MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX)
    if (not resolve_overlaps) or strandedness or str(lane_direction) != "split":
        return base_clearance

    feature_band_width_px = max(0.0, float(feature_band_width_px))
    adaptive_extra = 2.0 + min(2.0, 0.1 * feature_band_width_px)
    return base_clearance + adaptive_extra


def _effective_outer_text_clearance_px(
    *,
    resolve_overlaps: bool,
    strandedness: bool,
    lane_direction: str,
    feature_band_width_px: float,
) -> float:
    """Return text-to-feature clearance used for outer label bboxes."""
    base_clearance = float(MIN_OUTER_LABEL_TEXT_CLEARANCE_PX)
    if (not resolve_overlaps) or strandedness or str(lane_direction) != "split":
        return base_clearance

    feature_band_width_px = max(0.0, float(feature_band_width_px))
    adaptive_extra = 1.0 + min(1.0, 0.05 * feature_band_width_px)
    return base_clearance + adaptive_extra


def _resolve_label_spacing_px(
    label1: CircularLabel | None = None,
    label2: CircularLabel | None = None,
    *,
    base_margin_px: float | None = None,
) -> float:
    """Return the configured circular label spacing for a label pair."""
    if base_margin_px is not None:
        return max(0.0, float(base_margin_px))
    for label in (label1, label2):
        if label is None:
            continue
        if not bool(label.get("_enforce_spacing_px", False)):
            continue
        spacing_px = label.get("_spacing_px")
        if spacing_px is None:
            continue
        return max(0.0, float(spacing_px))
    return 0.0


def _effective_group_spacing_px(
    requested_spacing_px: float,
    *,
    legacy_spacing_px: float,
) -> float:
    """Return the requested circular label spacing without legacy dampening."""
    del legacy_spacing_px
    return max(0.0, float(requested_spacing_px))


def minimum_bbox_gap_px(label1: CircularLabel, label2: CircularLabel, base_margin_px: float | None = None) -> float:
    """Return a minimum spacing margin in pixels derived from label bbox size."""
    resolved_base_margin_px = _resolve_label_spacing_px(label1, label2, base_margin_px=base_margin_px)
    width_scale = max(float(label1.get("width_px", 0.0)), float(label2.get("width_px", 0.0)))
    height_scale = max(float(label1.get("height_px", 0.0)), float(label2.get("height_px", 0.0)))
    ratio_gap = MIN_BBOX_GAP_RATIO * max(width_scale, height_scale)
    resolve_extra = 0.0
    if ("min_outer_start_radius_px" in label1) or ("min_outer_start_radius_px" in label2):
        resolve_extra = RESOLVE_OVERLAP_MARGIN_EXTRA_PX
    return max(resolved_base_margin_px + resolve_extra, MIN_BBOX_GAP_FLOOR_PX, ratio_gap)


def _all_label_pairs(label_count: int) -> set[tuple[int, int]]:
    return {(idx, jdx) for idx in range(label_count) for jdx in range(idx + 1, label_count)}


def _label_pair_index_padding(labels: list[CircularLabel], *, use_min_gap: bool) -> float:
    if not use_min_gap:
        return 1.0
    return max((minimum_bbox_gap_px(label, label) for label in labels), default=0.0)


def _label_aabb_bucket_size(bboxes: list[Aabb], padding: float) -> float:
    max_extent = max(
        (
            max(bbox.max_x - bbox.min_x, bbox.max_y - bbox.min_y)
            for bbox in bboxes
        ),
        default=16.0,
    )
    return max(16.0, min(256.0, max_extent + float(padding)))


def _candidate_label_pairs(
    labels: list[CircularLabel],
    total_length: int,
    padding: float,
) -> set[tuple[int, int]]:
    label_count = len(labels)
    if label_count < 2:
        return set()

    bboxes: list[Aabb] = []
    try:
        for label in labels:
            bboxes.append(Aabb(*_label_bbox_rect(label, total_length, margin_px=0.0)))
    except (KeyError, TypeError, ValueError):
        return _all_label_pairs(label_count)

    pairs = candidate_aabb_pairs(
        bboxes,
        padding=padding,
        bucket_size=_label_aabb_bucket_size(bboxes, padding),
    )

    # Preserve explicit angular-neighbor checks even when broad-phase buckets prune aggressively.
    for idx in range(label_count - 1):
        pairs.add((idx, idx + 1))
    if label_count > 2:
        pairs.add((0, label_count - 1))
    return pairs


def _count_label_overlaps(labels: list[CircularLabel], total_length: int, use_min_gap: bool) -> int:
    """Count overlapping label pairs using the existing overlap predicates."""
    overlap_count = 0
    candidate_pairs = _candidate_label_pairs(
        labels,
        total_length,
        padding=_label_pair_index_padding(labels, use_min_gap=use_min_gap),
    )
    for idx, jdx in sorted(candidate_pairs):
        if _labels_overlap(labels[idx], labels[jdx], total_length, use_min_gap=use_min_gap):
            overlap_count += 1
    return overlap_count


def circular_label_font_size(cfg: GbdrawConfig, total_length: int) -> float:
    """The label font size of a record of this length (its size class)."""

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    return float(cfg.labels.font_size.for_length_param(length_param))


def _label_fit_box(label: CircularLabel, total_length: int) -> Polygon:
    """A horizontal label's box for the fit count, with the overlap margins (x 1.0 px, y 0.1 px)."""

    min_x, max_x = _label_x_bounds(label, 1.0)
    min_y, max_y = _label_y_bounds(label, total_length, 0.1)
    return ((min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y))


def _embedded_label_parts(label: CircularLabel, total_length: int) -> tuple[Polygon, ...]:
    """An embedded label as the annular sector its text path covers, cut into nearly flat parts."""

    radius = float(label["feature_center_radius_px"])
    half_height = 0.5 * float(label["height_px"])
    if radius <= half_height:
        return ()
    center_deg = 360.0 * (float(label["middle"]) / total_length) - 90.0
    half_span_deg = math.degrees(0.5 * float(label["width_px"]) / radius)
    return annular_sector_parts(radius - half_height, radius + half_height, center_deg, half_span_deg)


def circular_label_shapes(labels: Sequence[CircularLabel], total_length: int) -> list[LabelShape]:
    """The drawn labels of one record as the fit count sees them.

    Horizontal labels are their boxes, radial labels their rotated text boxes,
    and embedded labels the sector of their text path; outer and inner labels
    carry the leader the label drawer draws.
    """

    return [shape for shape in (circular_label_shape(label, total_length) for label in labels) if shape.parts]


def circular_label_shape(label: CircularLabel, total_length: int) -> LabelShape:
    """One drawn label as the fit count sees it (no parts: an embedded label too small to draw)."""

    if label.get("is_embedded"):
        return LabelShape("circular_embedded", _embedded_label_parts(label, total_length))
    corners = label.get("oriented_corners")
    parts = (tuple(corners),) if corners else (_label_fit_box(label, total_length),)
    kind: LabelKind = "circular_inner" if label.get("is_inner") else "circular_outer"
    return LabelShape(kind, parts, leader_polyline(label))


def _circular_label_box(text: str, font_family: str, font_size: float) -> tuple[float, float]:
    """The box a Circular feature label is laid out in: its drawn box, ``CIRCULAR_LABEL_GAP_EM`` wider."""

    width, height = drawn_text_box(text, font_family, font_size)
    return width + CIRCULAR_LABEL_GAP_EM * float(font_size), height


def circular_label_fits(
    labels: Sequence[CircularLabel],
    total_length: int,
    *,
    definition_text: Sequence[Aabb] = (),
) -> dict[LabelKind, LabelFit]:
    """The fit of one record's drawn labels by kind (horizontal labels touch with their margins)."""

    radial = any(label.get("oriented_corners") for label in labels)
    return label_fits_by_kind(
        circular_label_shapes(labels, total_length),
        touching=not radial,
        definition_text=definition_text,
    )


def _label_y_bounds(label: CircularLabel, total_len: int, minimum_margin: float) -> tuple[float, float]:
    """Return (min_y, max_y) for label bbox using the same baseline model as rendering."""
    low, low_margin, high, high_margin = _label_y_offsets(
        float(label["middle"]),
        int(total_len),
        float(label["height_px"]),
        float(minimum_margin),
        label.get("is_inner", False) is False,
    )
    start_y = float(label["start_y"])
    return (start_y + low) + low_margin, (start_y + high) + high_margin


@lru_cache(maxsize=65536)
def _label_y_offsets(
    label_middle: float,
    total_len: int,
    label_height: float,
    minimum_margin: float,
    is_outer: bool,
) -> tuple[float, float, float, float]:
    """The baseline model of a label's box: (low, low_margin, high, high_margin).

    The box is ``min_y = (start_y + low) + low_margin`` and ``max_y = (start_y + high) + high_margin``,
    in that order of additions; ``-0.0`` adds nothing. A label's offsets do not depend on its
    position, so the placement loops look them up once per label and step the position (OV-870):
    the step loops of `_legacy_place_labels_on_arc_fc` and `_legacy_improved_label_placement_fc`
    apply them inline; change them with this rule.
    """
    label_angle = (360.0 * (label_middle / total_len)) % 360.0
    margin_half = 0.5 * minimum_margin
    half = 0.5 * label_height
    full = 1.0 * label_height
    hanging = (-0.0, -0.0, full, margin_half)  # max_y = start_y + height + margin
    middle = (-half, -margin_half, half, margin_half)

    if is_outer:
        if 0 <= label_angle < 10 or not label_angle < 350:  # baseline = text-after-edge
            return (-full, -margin_half, -0.0, -0.0)
        if 170 <= label_angle < 190:  # baseline = hanging
            return (-0.0, -margin_half, full, margin_half)
        return middle  # baseline = middle
    if 0 <= label_angle < 10 or not label_angle < 350:  # baseline = hanging
        return hanging
    return middle  # baseline = middle


def _label_x_bounds(label: CircularLabel, minimum_margin: float) -> tuple[float, float]:
    """Return (min_x, max_x) for label bbox using current text-anchor semantics."""
    start_x = float(label["start_x"])
    right, left = _label_x_offsets(
        float(label["width_px"]), float(minimum_margin), label.get("is_inner", False) is False
    )
    low, low_margin, high, high_margin = right if start_x > 0 else left
    return (start_x + low) + low_margin, (start_x + high) + high_margin


@lru_cache(maxsize=65536)
def _label_x_offsets(
    label_width: float,
    minimum_margin: float,
    is_outer: bool,
) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float]]:
    """The text-anchor model of a label's box: offsets for ``start_x > 0`` and for the other labels.

    Each is (low, low_margin, high, high_margin) as in `_label_y_offsets`: the text runs away
    from the circle's vertical axis (outer labels) or toward it (inner labels).

    Inlined copies of how the offsets are applied (the ``start_x > 0`` side choice, the
    two additions, the `_bounds_overlap` test) live in the step loops of
    `_legacy_place_labels_on_arc_fc` and `_legacy_improved_label_placement_fc`; change
    them with this rule. `_labels_overlap` inlines only the interval test.
    """
    margin_half = 0.5 * minimum_margin
    away = (-0.0, -margin_half, label_width, margin_half)  # text from start_x to start_x + width
    toward = (-label_width, -margin_half, -0.0, margin_half)
    return (away, toward) if is_outer else (toward, away)


def _leader_start_meta(label: CircularLabel, total_length: int | None) -> tuple[str, float, float, float, float] | None:
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
    anchor_y = float(label.get("start_y", start_y))
    if bool(label.get("is_inner", False)) and abs(start_x) > HEMISPHERE_AXIS_EPS_PX:
        # For inner labels, keep elbow anchors on the radial outer edge of the text bbox.
        # This avoids leaders exiting from the center-facing edge on lower-left / lower-right arcs.
        if start_x < 0:
            return "left", min_x, y_min_usable, y_max_usable, anchor_y
        return "right", max_x, y_min_usable, y_max_usable, anchor_y

    if middle_x >= max_x:
        return "right", max_x, y_min_usable, y_max_usable, anchor_y
    if middle_x <= min_x:
        return "left", min_x, y_min_usable, y_max_usable, anchor_y
    if middle_y >= max_y:
        return "top", max_y, x_min_usable, x_max_usable, middle_x
    if middle_y <= min_y:
        return "bottom", min_y, x_min_usable, x_max_usable, middle_x

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


def _leader_start_point(label: CircularLabel, total_length: int | None = None) -> tuple[float, float]:
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


class _LeaderEdgeItem(NamedTuple):
    """One label on a box edge while its leader start is refined (``_leader_start_meta``)."""

    label: CircularLabel
    fixed: float
    lower: float
    upper: float
    preferred: float
    current: float
    order: float


def _leader_start_group_score(group: list[_LeaderEdgeItem], coords: list[float], side: str) -> tuple[int, int, float, float]:
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
        label = item.label
        if side in ("left", "right"):
            x = float(item.fixed)
            y = float(coord)
        else:
            x = float(coord)
            y = float(item.fixed)
        leader_sum += math.hypot(x - float(label["middle_x"]), y - float(label["middle_y"]))
        preference_deviation += abs(float(coord) - float(item.preferred))

    return inversion_count, close_pairs, leader_sum, preference_deviation


def _fit_monotonic_leader_coords(group: list[_LeaderEdgeItem], min_gap: float) -> list[float] | None:
    """Fit edge coordinates under bounds while enforcing monotonic order."""
    if not group:
        return []

    coords = [min(max(float(item.preferred), float(item.lower)), float(item.upper)) for item in group]

    for _ in range(max(3, len(coords) * 2)):
        changed = False
        for idx in range(1, len(coords)):
            lower_bound = max(float(group[idx].lower), coords[idx - 1] + min_gap)
            if coords[idx] < lower_bound:
                coords[idx] = lower_bound
                changed = True
            if coords[idx] > (float(group[idx].upper) + 1e-9):
                return None

        for idx in range(len(coords) - 2, -1, -1):
            upper_bound = min(float(group[idx].upper), coords[idx + 1] - min_gap)
            if coords[idx] > upper_bound:
                coords[idx] = upper_bound
                changed = True
            if coords[idx] < (float(group[idx].lower) - 1e-9):
                return None

        if not changed:
            break

    for idx, coord in enumerate(coords):
        if coord < (float(group[idx].lower) - 1e-9) or coord > (float(group[idx].upper) + 1e-9):
            return None
        if idx > 0 and coord < (coords[idx - 1] + min_gap - 1e-9):
            return None
    return coords


def _refine_leader_start_points(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """
    Improve readability by reducing leader-end inversions and excessive crowding.

    The optimization runs per bbox edge and keeps each anchor on the same edge.
    """
    if len(labels) < 2:
        return labels

    edge_groups: dict[str, list[_LeaderEdgeItem]] = {"left": [], "right": [], "bottom": [], "top": []}
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
            _LeaderEdgeItem(label, fixed, lower, upper, preferred, current_coord, order_coord)
        )

    for side, group in edge_groups.items():
        if len(group) < 2:
            continue
        group.sort(key=lambda item: float(item.order))
        current_coords = [float(item.current) for item in group]
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
            label = item.label
            if side in ("left", "right"):
                label["leader_start_x"] = float(item.fixed)
                label["leader_start_y"] = float(coord)
            else:
                label["leader_start_x"] = float(coord)
                label["leader_start_y"] = float(item.fixed)
    return labels


def _optimize_top_bottom_leader_start_points(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """
    Prefer long-edge midpoint anchors for top/bottom labels when safe.

    Candidate order is:
    1) long-edge midpoint
    2) middle_x projection on the edge (feature-facing)
    3) current anchor (fallback)

    Scoring priority:
    1) local label-leader collisions (minimum)
    2) distance to long-edge midpoint (minimum)
    3) leader length from text middle anchor (minimum)
    """
    if not labels:
        return labels

    # Built on first use and kept current: a label's update below only moves its own leader start.
    rects: list[tuple[float, float, float, float]] | None = None
    segments: list[list[tuple[tuple[float, float], tuple[float, float]]]] | None = None
    for idx, label in enumerate(labels):
        meta = _leader_start_meta(label, total_length)
        if meta is None:
            continue

        side, fixed_coord, lower, upper, _ = meta
        if side not in ("top", "bottom"):
            continue

        lower_bound = float(lower)
        upper_bound = float(upper)
        if upper_bound < lower_bound:
            lower_bound, upper_bound = upper_bound, lower_bound

        midpoint_coord = 0.5 * (lower_bound + upper_bound)
        feature_projection = min(
            max(float(label.get("middle_x", midpoint_coord)), lower_bound),
            upper_bound,
        )
        current_coord = min(
            max(float(label.get("leader_start_x", label.get("start_x", midpoint_coord))), lower_bound),
            upper_bound,
        )

        fixed = float(fixed_coord)
        candidate_coords = [midpoint_coord, feature_projection, current_coord]
        unique_coords: list[float] = []
        for coord in candidate_coords:
            if not any(math.isclose(coord, seen, abs_tol=1e-9) for seen in unique_coords):
                unique_coords.append(coord)

        best_score: tuple[int, float, float] | None = None
        best_coord = current_coord
        if rects is None or segments is None:
            rects, segments = _build_label_leader_collision_geometry(
                labels, total_length, LEADER_LABEL_COLLISION_MARGIN_PX
            )
        for coord in unique_coords:
            candidate_label = label.copy()
            candidate_label["leader_start_x"] = float(coord)
            candidate_label["leader_start_y"] = fixed
            local_collision_count = _count_local_leader_line_collisions(
                labels,
                total_length,
                idx,
                margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
                candidate=candidate_label,
                rects=rects,
                segments=segments,
            )
            midpoint_distance = abs(float(coord) - midpoint_coord)
            leader_length = math.hypot(
                float(coord) - float(label.get("middle_x", 0.0)),
                fixed - float(label.get("middle_y", 0.0)),
            )
            score = (local_collision_count, midpoint_distance, leader_length)
            if best_score is None or score < best_score:
                best_score = score
                best_coord = float(coord)

        label["leader_start_x"] = float(best_coord)
        label["leader_start_y"] = fixed
        segments[idx] = _leader_segments(label)

    return labels


def _assign_leader_start_points(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """Populate `leader_start_x/y` used for drawing the first leader segment.

    A contact the local repair chose stays (LP-D52), and its label leaves the edge-group refinement.
    """
    engine_labels: list[CircularLabel] = []
    for label in labels:
        chosen = chosen_contact_end(label, total_length)
        leader_x, leader_y = chosen if chosen is not None else _leader_start_point(label, total_length)
        label["leader_start_x"] = leader_x
        label["leader_start_y"] = leader_y
        if chosen is None:
            engine_labels.append(label)
    _refine_leader_start_points(engine_labels, total_length)
    _optimize_top_bottom_leader_start_points(engine_labels, total_length)
    return labels


def assign_leader_start_points(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """Public wrapper to (re)compute leader start anchors after label motion."""
    return _assign_leader_start_points(labels, total_length)


def _leader_length_px(label: CircularLabel, total_length: int | None = None) -> float:
    """Return leader length from label middle anchor to chosen label-side leader anchor."""
    if "middle_x" not in label or "middle_y" not in label:
        return 0.0
    leader_start_x, leader_start_y = _leader_start_point(label, total_length)
    return math.hypot(
        leader_start_x - float(label["middle_x"]),
        leader_start_y - float(label["middle_y"]),
    )


def _leader_distance_score(labels: list[CircularLabel], total_length: int | None = None) -> tuple[float, float]:
    """Return (sum_length, max_length) for leader segments."""
    if not labels:
        return 0.0, 0.0
    leader_lengths = [_leader_length_px(label, total_length) for label in labels]
    return sum(leader_lengths), max(leader_lengths)


def _label_bbox_rect(label: CircularLabel, total_length: int, margin_px: float) -> tuple[float, float, float, float]:
    """Return label bbox rectangle as (min_x, min_y, max_x, max_y)."""
    min_x, max_x = _label_x_bounds(label, minimum_margin=float(margin_px))
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=float(margin_px))
    return min_x, min_y, max_x, max_y


def _leader_segments(label: CircularLabel) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Return two leader segments drawn for non-embedded circular labels."""
    middle_x = float(label.get("middle_x", 0.0))
    middle_y = float(label.get("middle_y", 0.0))
    leader_start_x = float(label.get("leader_start_x", label.get("start_x", 0.0)))
    leader_start_y = float(label.get("leader_start_y", label.get("start_y", 0.0)))
    feature_anchor_x = float(label.get("feature_anchor_x", label.get("feature_middle_x", 0.0)))
    feature_anchor_y = float(label.get("feature_anchor_y", label.get("feature_middle_y", 0.0)))
    middle_point = (middle_x, middle_y)
    return [
        (middle_point, (leader_start_x, leader_start_y)),
        (middle_point, (feature_anchor_x, feature_anchor_y)),
    ]


def _primary_leader_segment(label: CircularLabel) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return the primary leader segment (middle -> leader-start) for line-cross checks."""
    middle_x = float(label.get("middle_x", 0.0))
    middle_y = float(label.get("middle_y", 0.0))
    leader_start_x = float(label.get("leader_start_x", label.get("start_x", 0.0)))
    leader_start_y = float(label.get("leader_start_y", label.get("start_y", 0.0)))
    return ((middle_x, middle_y), (leader_start_x, leader_start_y))


def _build_label_leader_collision_geometry(
    labels: list[CircularLabel],
    total_length: int,
    margin_px: float,
) -> tuple[
    list[tuple[float, float, float, float]],
    list[list[tuple[tuple[float, float], tuple[float, float]]]],
]:
    """Precompute per-label rectangles and leader segments used by collision checks."""
    rects = [_label_bbox_rect(label, total_length, margin_px=margin_px) for label in labels]
    segments = [_leader_segments(label) for label in labels]
    return rects, segments


def _build_primary_leader_segments(
    labels: list[CircularLabel],
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Precompute primary leader segments used by line-cross checks."""
    return [_primary_leader_segment(label) for label in labels]


def _segments_properly_intersect(
    seg_a: tuple[tuple[float, float], tuple[float, float]],
    seg_b: tuple[tuple[float, float], tuple[float, float]],
) -> bool:
    """
    Return True when two segments intersect, including endpoint and boundary touch.
    """
    (ax1, ay1), (ax2, ay2) = seg_a
    (bx1, by1), (bx2, by2) = seg_b
    eps = 1e-9
    # Segments whose bounding boxes lie more than 1e-6 px apart neither cross nor touch within eps.
    if (
        (ax1 if ax1 > ax2 else ax2) < (bx1 if bx1 < bx2 else bx2) - 1e-6
        or (bx1 if bx1 > bx2 else bx2) < (ax1 if ax1 < ax2 else ax2) - 1e-6
        or (ay1 if ay1 > ay2 else ay2) < (by1 if by1 < by2 else by2) - 1e-6
        or (by1 if by1 > by2 else by2) < (ay1 if ay1 < ay2 else ay2) - 1e-6
    ):
        return False

    def _on_segment(px: float, py: float, qx: float, qy: float, rx: float, ry: float) -> bool:
        return (
            (min(px, rx) - eps) <= qx <= (max(px, rx) + eps)
            and (min(py, ry) - eps) <= qy <= (max(py, ry) + eps)
        )

    # Orientations of each end against the other segment, inlined (OV-870).
    o1 = ((ax2 - ax1) * (by1 - ay1)) - ((ay2 - ay1) * (bx1 - ax1))
    o2 = ((ax2 - ax1) * (by2 - ay1)) - ((ay2 - ay1) * (bx2 - ax1))
    o3 = ((bx2 - bx1) * (ay1 - by1)) - ((by2 - by1) * (ax1 - bx1))
    o4 = ((bx2 - bx1) * (ay2 - by1)) - ((by2 - by1) * (ax2 - bx1))

    if ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and (
        (o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps)
    ):
        return True
    if abs(o1) <= eps and _on_segment(ax1, ay1, bx1, by1, ax2, ay2):
        return True
    if abs(o2) <= eps and _on_segment(ax1, ay1, bx2, by2, ax2, ay2):
        return True
    if abs(o3) <= eps and _on_segment(bx1, by1, ax1, ay1, bx2, by2):
        return True
    if abs(o4) <= eps and _on_segment(bx1, by1, ax2, ay2, bx2, by2):
        return True
    return False


def _count_leader_line_intersections(
    labels: list[CircularLabel],
    *,
    primary_segments: list[tuple[tuple[float, float], tuple[float, float]]] | None = None,
) -> int:
    """Count intersections between primary leader segments."""
    if primary_segments is None:
        primary_segments = _build_primary_leader_segments(labels)
    intersection_count = 0
    for idx, jdx in touching_aabb_pairs([bounding_aabb(segment, 1.0) for segment in primary_segments]):
        if _segments_properly_intersect(primary_segments[idx], primary_segments[jdx]):
            intersection_count += 1
    return intersection_count


def _count_local_leader_line_intersections(
    labels: list[CircularLabel],
    idx: int,
    *,
    candidate: CircularLabel | None = None,
    primary_segments: list[tuple[tuple[float, float], tuple[float, float]]] | None = None,
    candidate_segment: tuple[tuple[float, float], tuple[float, float]] | None = None,
    peers: list[int] | None = None,
    cap: int | None = None,
) -> int:
    """Count primary leader-line intersections involving one label.

    `peers` limits the scan to labels that can touch it; `cap` stops counting once the count exceeds it.
    """
    if idx < 0 or idx >= len(labels):
        return 0
    if primary_segments is None:
        primary_segments = _build_primary_leader_segments(labels)
    if candidate_segment is not None:
        target_segment = candidate_segment
    elif candidate is not None:
        target_segment = _primary_leader_segment(candidate)
    else:
        target_segment = primary_segments[idx]
    intersection_count = 0
    for peer_idx in range(len(primary_segments)) if peers is None else peers:
        if peer_idx == idx:
            continue
        if _segments_properly_intersect(target_segment, primary_segments[peer_idx]):
            intersection_count += 1
            if cap is not None and intersection_count > cap:
                break
    return intersection_count


def _segment_intersects_rect_interior(
    p1: tuple[float, float],
    p2: tuple[float, float],
    rect: tuple[float, float, float, float],
) -> bool:
    """
    Return True when a segment intersects the rectangle interior.

    Boundary-only touches are treated as non-overlaps.
    """
    return _segments_through_rect_interior(((p1, p2),), rect) > 0


def _segments_through_rect_interior(
    segments: Iterable[tuple[tuple[float, float], tuple[float, float]]],
    rect: tuple[float, float, float, float],
) -> int:
    """How many of `segments` intersect the rectangle interior; boundary-only touches do not count.

    One call per segment list: the collision counts test about 0.6 M segments per search (OV-870).
    """
    min_x, min_y, max_x, max_y = map(float, rect)
    if max_x <= min_x or max_y <= min_y:
        return 0
    eps = 1e-9
    count = 0
    for p1, p2 in segments:
        x1, y1 = float(p1[0]), float(p1[1])
        x2, y2 = float(p2[0]), float(p2[1])
        # A segment on one side of the rectangle has no interior point.
        if (
            (x1 <= min_x and x2 <= min_x)
            or (x1 >= max_x and x2 >= max_x)
            or (y1 <= min_y and y2 <= min_y)
            or (y1 >= max_y and y2 >= max_y)
        ):
            continue

        dx = x2 - x1
        dy = y2 - y1
        p_vals = (-dx, dx, -dy, dy)
        q_vals = (x1 - min_x, max_x - x1, y1 - min_y, max_y - y1)

        t_enter = 0.0
        t_exit = 1.0
        missed = False
        for p_val, q_val in zip(p_vals, q_vals):
            if abs(p_val) <= eps:
                if q_val < 0.0:
                    missed = True
                    break
                continue
            t_hit = q_val / p_val
            if p_val < 0.0:
                t_enter = max(t_enter, t_hit)
            else:
                t_exit = min(t_exit, t_hit)
            if t_enter > t_exit:
                missed = True
                break
        if missed:
            continue

        clipped_start = max(0.0, t_enter)
        clipped_end = min(1.0, t_exit)
        if clipped_end <= (clipped_start + eps):
            continue

        t_mid = 0.5 * (clipped_start + clipped_end)
        x_mid = x1 + (dx * t_mid)
        y_mid = y1 + (dy * t_mid)
        if (min_x + eps) < x_mid < (max_x - eps) and (min_y + eps) < y_mid < (max_y - eps):
            count += 1
    return count


def _count_label_leader_line_collisions(
    labels: list[CircularLabel],
    total_length: int,
    margin_px: float,
    *,
    rects: list[tuple[float, float, float, float]] | None = None,
    segments: list[list[tuple[tuple[float, float], tuple[float, float]]]] | None = None,
) -> int:
    """Count segment-vs-label collisions for all label pairs."""
    if rects is None or segments is None:
        rects, segments = _build_label_leader_collision_geometry(labels, total_length, margin_px)
    collision_count = 0
    boxes = [
        bounding_aabb(
            [(rect[0], rect[1]), (rect[2], rect[3]), *(point for segment in label_segments for point in segment)],
            1.0,
        )
        for rect, label_segments in zip(rects, segments)
    ]
    for idx, jdx in touching_aabb_pairs(boxes):
        collision_count += _segments_through_rect_interior(segments[idx], rects[jdx])
        collision_count += _segments_through_rect_interior(segments[jdx], rects[idx])
    return collision_count


def _count_local_leader_line_collisions(
    labels: list[CircularLabel],
    total_length: int,
    idx: int,
    margin_px: float,
    *,
    candidate: CircularLabel | None = None,
    rects: list[tuple[float, float, float, float]] | None = None,
    segments: list[list[tuple[tuple[float, float], tuple[float, float]]]] | None = None,
    candidate_rect: tuple[float, float, float, float] | None = None,
    candidate_segments: list[tuple[tuple[float, float], tuple[float, float]]] | None = None,
    peers: list[int] | None = None,
    cap: int | None = None,
) -> int:
    """Count collisions involving one label (optionally with a candidate geometry).

    `peers` limits the scan to labels that can touch it; `cap` stops counting once the count exceeds it.
    """
    if idx < 0 or idx >= len(labels):
        return 0
    if rects is None or segments is None:
        rects, segments = _build_label_leader_collision_geometry(labels, total_length, margin_px)
    if candidate_rect is not None:
        target_rect = candidate_rect
    elif candidate is not None:
        target_rect = _label_bbox_rect(candidate, total_length, margin_px=margin_px)
    else:
        target_rect = rects[idx]

    if candidate_segments is not None:
        target_segments = candidate_segments
    elif candidate is not None:
        target_segments = _leader_segments(candidate)
    else:
        target_segments = segments[idx]
    collision_count = 0

    for peer_idx in range(len(segments)) if peers is None else peers:
        if peer_idx == idx:
            continue
        collision_count += _segments_through_rect_interior(target_segments, rects[peer_idx])
        collision_count += _segments_through_rect_interior(segments[peer_idx], target_rect)
        if cap is not None and collision_count > cap:
            break
    return collision_count


def _reach_pad(margin_bound_px: float) -> float:
    """Pad of a reach box for predicates whose margins are at most `margin_bound_px`.

    A bbox at margin m grows by m/2 on each side; the pad covers half the bound (at least
    1 px, the plain predicates' margin) plus 1 px that keeps disjoint reach boxes strictly apart.
    """
    return 0.5 * max(margin_bound_px, 1.0) + 1.0


def _label_reach_pad(labels: list[CircularLabel], *, extra_margin_px: float = 0.0) -> float:
    """Reach pad covering the bbox-overlap and leader predicates of these labels.

    The min-gap margin of a pair never exceeds the largest single-label min gap plus
    RESOLVE_OVERLAP_MARGIN_EXTRA_PX; `extra_margin_px` covers a caller's collision rect margin.
    """
    gap_bound = (
        max((minimum_bbox_gap_px(label, label) for label in labels), default=0.0)
        + RESOLVE_OVERLAP_MARGIN_EXTRA_PX
    )
    return _reach_pad(max(gap_bound, float(extra_margin_px)))


def _label_reach_box(
    label: CircularLabel,
    pad: float,
    *,
    anchors: Sequence[tuple[float, float]] | None = None,
    points: Sequence[tuple[float, float]] = (),
) -> Aabb | None:
    """Return a box holding the label's bbox anchored at each of `anchors`, and `points`.

    A label's bbox at margin m lies within width_px and height_px (plus m/2) of its anchor
    (start_x, start_y) whatever its text anchor and baseline, so full sizes padded by `pad`
    (see `_reach_pad`) hold it. `anchors` defaults to where the label is. ``None``: unknown
    extent (a non-finite coordinate), which touches every box.
    """
    if anchors is None:
        anchors = [(float(label["start_x"]), float(label["start_y"]))]
    box = bounding_aabb(anchors, float(label["width_px"]) + pad, float(label["height_px"]) + pad)
    if box is None or not points:
        return box
    point_box = bounding_aabb(points, pad)
    return None if point_box is None else union_aabbs((box, point_box))


def _labels_overlap(label_a: CircularLabel, label_b: CircularLabel, total_length: int, *, use_min_gap: bool) -> bool:
    """Return whether two label bboxes overlap, at their min gap or at the plain margins."""
    if use_min_gap:
        margin = minimum_bbox_gap_px(label_a, label_b)
        y_margin = margin
        x_margin = margin
    else:
        y_margin = 0.1
        x_margin = 1.0
    # `x_overlap` and then `y_overlap`, with the interval test inlined (OV-870).
    min_a, max_a = _label_x_bounds(label_a, x_margin)
    min_b, max_b = _label_x_bounds(label_b, x_margin)
    if not (max_a >= min_b if min_a < min_b else max_b >= min_a):
        return False
    min_a, max_a = _label_y_bounds(label_a, total_length, y_margin)
    min_b, max_b = _label_y_bounds(label_b, total_length, y_margin)
    return max_a >= min_b if min_a < min_b else max_b >= min_a


def _count_local_label_overlaps(
    labels: list[CircularLabel],
    total_length: int,
    label_idx: int,
    candidate_label: CircularLabel,
    *,
    use_min_gap: bool,
    peers: Iterable[int] | None = None,
    cap: int | None = None,
) -> int:
    """Count the labels (except `label_idx`) that `candidate_label` overlaps.

    `peers` limits the scan to labels that can touch it; `cap` stops counting once the count exceeds it.
    """
    overlap_count = 0
    for other_idx in range(len(labels)) if peers is None else peers:
        if other_idx == label_idx:
            continue
        if _labels_overlap(candidate_label, labels[other_idx], total_length, use_min_gap=use_min_gap):
            overlap_count += 1
            if cap is not None and overlap_count > cap:
                break
    return overlap_count


class _MovedLabelOverlaps:
    """Count the overlapping pairs of a label set after some labels move, from the pairs before.

    Only pairs with a moved label can change. `overlaps(label_i, label_j, i, j)` (i < j) is the pair
    predicate; it must be false for labels whose reach boxes at `reach_pad` (see `_reach_pad`) do
    not touch, so the peers of a moved label come from the reach boxes.
    """

    def __init__(
        self,
        labels: Sequence[CircularLabel],
        overlapping_pairs: Iterable[tuple[int, int]],
        overlaps: Callable[[CircularLabel, CircularLabel, int, int], bool],
        reach_pad: float,
    ) -> None:
        self._labels = labels
        self._overlaps = overlaps
        self._reach_pad = reach_pad
        self._partners: list[list[int]] = [[] for _ in labels]
        self.count = 0
        for idx, jdx in overlapping_pairs:
            self._partners[idx].append(jdx)
            self._partners[jdx].append(idx)
            self.count += 1
        self._boxes = [_label_reach_box(label, reach_pad) for label in labels]
        self._unknown = [idx for idx, box in enumerate(self._boxes) if box is None]
        finite = [box for box in self._boxes if box is not None]
        self._index = AabbIndex(bucket_size=_label_aabb_bucket_size(finite, 0.0))
        for idx, box in enumerate(self._boxes):
            if box is not None:
                self._index.insert(idx, box)

    def count_after_move(self, moved: Mapping[int, CircularLabel], cap: int | None = None) -> int:
        """Overlapping pairs once each label `idx` in `moved` is replaced by `moved[idx]`.

        With `cap` the count may stop once it exceeds `cap`: exact up to `cap`, above it otherwise.
        """
        total = self.count
        for idx in moved:
            for partner in self._partners[idx]:
                if partner not in moved or partner > idx:
                    total -= 1
        if cap is not None and total > cap:
            return total
        labels = self._labels
        overlaps = self._overlaps
        moved_ids = sorted(moved)
        moved_boxes = [_label_reach_box(moved[idx], self._reach_pad) for idx in moved_ids]
        for position, idx in enumerate(moved_ids):
            label = moved[idx]
            box = moved_boxes[position]
            nearby = None if box is None else cast("list[int]", self._index.query(box)) + self._unknown
            for peer in touching_aabb_indices(self._boxes, box, candidates=nearby):
                if peer in moved:
                    continue
                if overlaps(labels[peer], label, peer, idx) if peer < idx else overlaps(label, labels[peer], idx, peer):
                    total += 1
                    if cap is not None and total > cap:
                        return total
            for later in range(position + 1, len(moved_ids)):
                peer = moved_ids[later]
                if touching_aabb_indices(moved_boxes, box, candidates=(later,)) and overlaps(
                    label, moved[peer], idx, peer
                ):
                    total += 1
                    if cap is not None and total > cap:
                        return total
        return total


def _resolve_label_leader_line_collisions(
    labels: list[CircularLabel],
    total_length: int,
    *,
    margin_px: float,
    step_deg: float,
    max_shift_deg: float,
    max_passes: int,
    min_order_gap_deg: float,
) -> list[CircularLabel]:
    """Reduce label-vs-leader collisions by angular shifts without overlap regressions."""
    if len(labels) < 2:
        return labels

    step_deg = max(float(step_deg), 1e-6)
    max_shift_deg = max(float(max_shift_deg), step_deg)
    max_steps = max(1, int(max_shift_deg / step_deg))
    max_passes = max(1, int(max_passes))
    min_order_gap_deg = max(float(min_order_gap_deg), 1e-6)
    rects, segments = _build_label_leader_collision_geometry(labels, total_length, margin_px)
    primary_segments = _build_primary_leader_segments(labels)

    if (
        _count_label_leader_line_collisions(
            labels,
            total_length,
            margin_px=margin_px,
            rects=rects,
            segments=segments,
        )
        == 0
        and _count_leader_line_intersections(labels, primary_segments=primary_segments) == 0
    ):
        return labels

    # Peer pruning (output-preserving): every predicate below is false for two labels whose reach
    # boxes are disjoint, so the counts only visit labels whose reach box meets the label's sweep.
    reach_pad = _label_reach_pad(labels, extra_margin_px=margin_px)

    def _label_reach_boxes(idx: int) -> tuple[Aabb | None, Aabb | None, Aabb | None]:
        """Reach boxes of the text (bbox overlaps), the primary leader (crossings), and all (collisions)."""
        label = labels[idx]
        rect = rects[idx]
        return (
            _label_reach_box(label, reach_pad),
            bounding_aabb(primary_segments[idx], reach_pad),
            _label_reach_box(
                label,
                reach_pad,
                points=[
                    (rect[0], rect[1]),
                    (rect[2], rect[3]),
                    *(point for segment in segments[idx] for point in segment),
                    *primary_segments[idx],
                ],
            ),
        )

    def _set_reach_boxes(idx: int) -> None:
        text_boxes[idx], primary_boxes[idx], reach_boxes[idx] = _label_reach_boxes(idx)

    for _ in range(max_passes):
        current_total_bbox_collisions = _count_label_leader_line_collisions(
            labels,
            total_length,
            margin_px=margin_px,
            rects=rects,
            segments=segments,
        )
        current_total_line_intersections = _count_leader_line_intersections(
            labels,
            primary_segments=primary_segments,
        )
        if current_total_bbox_collisions == 0 and current_total_line_intersections == 0:
            break

        current_plain_overlaps = _count_label_overlaps(labels, total_length, use_min_gap=False)
        current_min_gap_overlaps = _count_label_overlaps(labels, total_length, use_min_gap=True)
        pass_changed = False
        unwrapped_angles = _derive_monotonic_unwrapped_angles(labels, total_length)
        text_boxes: list[Aabb | None] = [None for _ in labels]
        primary_boxes: list[Aabb | None] = [None for _ in labels]
        reach_boxes: list[Aabb | None] = [None for _ in labels]
        for idx in range(len(labels)):
            _set_reach_boxes(idx)
        local_scores = []
        for idx in range(len(labels)):
            crossing_peers = touching_aabb_indices(primary_boxes, primary_boxes[idx], skip=idx)
            near_peers = touching_aabb_indices(reach_boxes, reach_boxes[idx], skip=idx)
            local_scores.append(
                (
                    _count_local_leader_line_intersections(
                        labels,
                        idx,
                        primary_segments=primary_segments,
                        peers=crossing_peers,
                    ),
                    _count_local_leader_line_collisions(
                        labels,
                        total_length,
                        idx,
                        margin_px=margin_px,
                        rects=rects,
                        segments=segments,
                        peers=near_peers,
                    ),
                    idx,
                )
            )
        local_scores.sort(key=lambda item: (item[0], item[1]), reverse=True)

        for local_line_intersection_count, local_collision_count, label_idx in local_scores:
            if local_line_intersection_count <= 0 and local_collision_count <= 0:
                break

            label = labels[label_idx]
            start_x = float(label.get("start_x", 0.0))
            start_y = float(label.get("start_y", 0.0))
            radius = math.hypot(start_x, start_y)
            if radius <= 1e-9:
                continue

            current_angle = float(unwrapped_angles[label_idx])
            lower_bound, upper_bound = _cyclic_order_bounds(unwrapped_angles, label_idx, min_order_gap_deg)
            if lower_bound >= upper_bound:
                continue

            current_half = _current_half_from_x(start_x, axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
            current_plain_for_label = _count_local_label_overlaps(
                labels, total_length, label_idx, label, use_min_gap=False
            )
            current_min_gap_for_label = _count_local_label_overlaps(
                labels, total_length, label_idx, label, use_min_gap=True
            )
            current_score: tuple[int, int, int, tuple[int, float], float] = (
                local_line_intersection_count,
                local_collision_count,
                current_min_gap_for_label,
                (0, 0.0),
                _leader_length_px(label, total_length),
            )
            best_score = current_score
            best_candidate: CircularLabel | None = None
            best_candidate_rect: tuple[float, float, float, float] | None = None
            best_candidate_segments: list[tuple[tuple[float, float], tuple[float, float]]] | None = None
            best_candidate_primary_segment: tuple[tuple[float, float], tuple[float, float]] | None = None
            best_angle: float | None = None
            best_plain_total = current_plain_overlaps
            best_min_gap_total = current_min_gap_overlaps

            # Candidate anchors per step, in the search order; their sweep box selects the peers. Each
            # carries its move (0, angle offset) or (1, outward shift px): any angular shift ranks first.
            step_positions: list[list[tuple[tuple[int, float], float, float, float]]] = []
            for step_idx in range(1, max_steps + 1):
                angle_offset = step_idx * step_deg
                positions: list[tuple[tuple[int, float], float, float, float]] = []
                for direction in (1.0, -1.0):
                    candidate_angle = current_angle + (direction * angle_offset)
                    if candidate_angle <= (lower_bound + 1e-9) or candidate_angle >= (upper_bound - 1e-9):
                        continue

                    candidate_x = radius * math.cos(math.radians(candidate_angle))
                    candidate_y = radius * math.sin(math.radians(candidate_angle))
                    if current_half != 0 and _current_half_from_x(
                        candidate_x, axis_eps_px=HEMISPHERE_AXIS_EPS_PX
                    ) != current_half:
                        continue
                    positions.append(((0, angle_offset), candidate_angle, candidate_x, candidate_y))
                step_positions.append(positions)
            # OV-1823: the same angle, further out (an outer label only). An angular shift cannot clear
            # a label that sits across a neighbour's radial leader between the neighbour's feature and text.
            radial_positions: list[tuple[tuple[int, float], float, float, float]] = []
            if label.get("is_inner", False) is False:
                max_shift_px = LEADER_LABEL_MAX_RADIAL_SHIFT_HEIGHTS * float(label["height_px"])
                radial_steps = int(max_shift_px / LEADER_LABEL_RADIAL_STEP_PX)
                for step_idx in range(1, radial_steps + 1):
                    shift_px = step_idx * LEADER_LABEL_RADIAL_STEP_PX
                    scale = (radius + shift_px) / radius
                    radial_positions.append(((1, shift_px), current_angle, start_x * scale, start_y * scale))
            sweep_points = [
                (item[2], item[3]) for positions in (*step_positions, radial_positions) for item in positions
            ]
            if not sweep_points:
                continue
            middle_point, feature_anchor_point = _leader_segments(label)[1]
            text_peers = touching_aabb_indices(
                text_boxes, _label_reach_box(label, reach_pad, anchors=sweep_points), skip=label_idx
            )
            crossing_peers = touching_aabb_indices(
                primary_boxes,
                _label_reach_box(label, reach_pad, anchors=sweep_points, points=[middle_point]),
                skip=label_idx,
            )
            collision_peers = touching_aabb_indices(
                reach_boxes,
                _label_reach_box(
                    label, reach_pad, anchors=sweep_points, points=[middle_point, feature_anchor_point]
                ),
                skip=label_idx,
            )

            for positions in (*step_positions, radial_positions):
                if positions is radial_positions and best_score[1] == 0:
                    break  # an angular shift (or none) clears every leader collision of this label
                for move, candidate_angle, candidate_x, candidate_y in positions:
                    candidate = label.copy()
                    candidate["start_x"] = float(candidate_x)
                    candidate["start_y"] = float(candidate_y)
                    candidate["angle_unwrapped"] = float(candidate_angle)
                    leader_start_x, leader_start_y = _leader_start_point(candidate, total_length)
                    candidate["leader_start_x"] = float(leader_start_x)
                    candidate["leader_start_y"] = float(leader_start_y)
                    candidate_rect = _label_bbox_rect(candidate, total_length, margin_px=margin_px)
                    candidate_segments = _leader_segments(candidate)
                    candidate_primary_segment = _primary_leader_segment(candidate)

                    # Same rejections as counting everything first: a count only has to be exact
                    # while it stays within the bound that would reject the candidate.
                    line_cap = min(local_line_intersection_count, best_score[0])
                    candidate_local_line_intersection_count = _count_local_leader_line_intersections(
                        labels,
                        label_idx,
                        candidate=candidate,
                        primary_segments=primary_segments,
                        candidate_segment=candidate_primary_segment,
                        peers=crossing_peers,
                        cap=line_cap,
                    )
                    if candidate_local_line_intersection_count > line_cap:
                        continue
                    collision_cap = (
                        min(local_collision_count, best_score[1])
                        if candidate_local_line_intersection_count == best_score[0]
                        else local_collision_count
                    )
                    candidate_local_collision_count = _count_local_leader_line_collisions(
                        labels,
                        total_length,
                        label_idx,
                        margin_px=margin_px,
                        candidate=candidate,
                        rects=rects,
                        segments=segments,
                        candidate_rect=candidate_rect,
                        candidate_segments=candidate_segments,
                        peers=collision_peers,
                        cap=collision_cap,
                    )
                    if candidate_local_collision_count > collision_cap:
                        continue

                    candidate_plain_for_label = _count_local_label_overlaps(
                        labels,
                        total_length,
                        label_idx,
                        candidate,
                        use_min_gap=False,
                        peers=text_peers,
                        cap=current_plain_for_label,
                    )
                    candidate_plain_total = (
                        current_plain_overlaps - current_plain_for_label + candidate_plain_for_label
                    )
                    if candidate_plain_total > current_plain_overlaps:
                        continue

                    candidate_min_gap_for_label = _count_local_label_overlaps(
                        labels,
                        total_length,
                        label_idx,
                        candidate,
                        use_min_gap=True,
                        peers=text_peers,
                        cap=current_min_gap_for_label,
                    )
                    candidate_min_gap_total = (
                        current_min_gap_overlaps - current_min_gap_for_label + candidate_min_gap_for_label
                    )
                    if candidate_min_gap_total > current_min_gap_overlaps:
                        continue

                    candidate_score = (
                        candidate_local_line_intersection_count,
                        candidate_local_collision_count,
                        candidate_min_gap_for_label,
                        move,
                        _leader_length_px(candidate, total_length),
                    )
                    if candidate_score < best_score:
                        best_score = candidate_score
                        best_candidate = candidate
                        best_candidate_rect = candidate_rect
                        best_candidate_segments = candidate_segments
                        best_candidate_primary_segment = candidate_primary_segment
                        best_angle = candidate_angle
                        best_plain_total = candidate_plain_total
                        best_min_gap_total = candidate_min_gap_total

                if best_score[0] == 0 and best_score[1] == 0 and best_score[2] == 0:
                    break

            if best_candidate is None:
                continue
            if (
                best_score[0] >= local_line_intersection_count
                and best_score[1] >= local_collision_count
            ):
                continue

            label["start_x"] = float(best_candidate["start_x"])
            label["start_y"] = float(best_candidate["start_y"])
            label["angle_unwrapped"] = float(best_angle if best_angle is not None else current_angle)
            label["leader_start_x"] = float(best_candidate["leader_start_x"])
            label["leader_start_y"] = float(best_candidate["leader_start_y"])
            unwrapped_angles[label_idx] = float(best_angle if best_angle is not None else current_angle)
            if best_candidate_rect is not None:
                rects[label_idx] = best_candidate_rect
            if best_candidate_segments is not None:
                segments[label_idx] = best_candidate_segments
            if best_candidate_primary_segment is not None:
                primary_segments[label_idx] = best_candidate_primary_segment
            _set_reach_boxes(label_idx)
            current_plain_overlaps = best_plain_total
            current_min_gap_overlaps = best_min_gap_total
            pass_changed = True

        labels = assign_leader_start_points(labels, total_length)
        rects, segments = _build_label_leader_collision_geometry(labels, total_length, margin_px)
        primary_segments = _build_primary_leader_segments(labels)
        if (
            _count_label_leader_line_collisions(
                labels,
                total_length,
                margin_px=margin_px,
                rects=rects,
                segments=segments,
            )
            == 0
            and _count_leader_line_intersections(labels, primary_segments=primary_segments) == 0
        ):
            break
        if not pass_changed:
            break

    for label in labels:
        start_x = float(label.get("start_x", 0.0))
        start_y = float(label.get("start_y", 0.0))
        angle_raw = math.degrees(math.atan2(start_y, start_x))
        target_angle = float(
            label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
        )
        label["angle_unwrapped"] = normalize_angle_near_reference(angle_raw, target_angle)

    return labels


def _label_bbox_min_radius(label: CircularLabel, total_length: int, minimum_margin: float = 0.0) -> float:
    """Return minimum distance from origin to the label bbox."""
    min_x, max_x = _label_x_bounds(label, minimum_margin=minimum_margin)
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=minimum_margin)
    closest_x = min(max(0.0, min_x), max_x)
    closest_y = min(max(0.0, min_y), max_y)
    return math.hypot(closest_x, closest_y)


def _label_bbox_max_radius(label: CircularLabel, total_length: int, minimum_margin: float = 0.0) -> float:
    """Return maximum distance from origin to the label bbox."""
    min_x, max_x = _label_x_bounds(label, minimum_margin=minimum_margin)
    min_y, max_y = _label_y_bounds(label, total_length, minimum_margin=minimum_margin)
    corner_radii = (
        math.hypot(min_x, min_y),
        math.hypot(min_x, max_y),
        math.hypot(max_x, min_y),
        math.hypot(max_x, max_y),
    )
    return max(corner_radii)


def _build_feature_radius_intervals(
    feature_dict: dict,
    total_length: int,
    radius: float,
    track_ratio: float,
    profile: CircularRenderProfile,
    *,
    feature_track_ratio_factor_override: float | None = None,
    feature_layout: CircularFeatureLayout | None = None,
    track_preset: str | None = None,
) -> list[tuple[float, float, float, float]]:
    """Build local genome intervals as (start, end, inner_radius_px, outer_radius_px)."""
    cfg = profile.config
    if total_length <= 0:
        return []
    if feature_layout is not None:
        return feature_radius_intervals(feature_dict, total_length, feature_layout)

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = (
        float(feature_track_ratio_factor_override)
        if feature_track_ratio_factor_override is not None
        else float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    )
    cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
    track_type = _normalize_circular_track_preset(track_preset or cfg.canvas.circular.track_type)
    strandedness = profile.strandedness

    intervals: list[tuple[float, float, float, float]] = []
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
            factor_inner = min(float(factors[0]), float(factors[2]))
            factor_outer = max(float(factors[0]), float(factors[2]))
            feature_inner_radius = float(radius) * factor_inner
            feature_outer_radius = float(radius) * factor_outer

            if coordinate_start < coordinate_end:
                intervals.append((coordinate_start, coordinate_end, feature_inner_radius, feature_outer_radius))
            else:
                intervals.append((coordinate_start, total_length_float, feature_inner_radius, feature_outer_radius))
                intervals.append((0.0, coordinate_end, feature_inner_radius, feature_outer_radius))

    return intervals


def _build_outer_feature_radius_intervals(
    feature_dict: dict,
    total_length: int,
    radius: float,
    track_ratio: float,
    profile: CircularRenderProfile,
    *,
    feature_track_ratio_factor_override: float | None = None,
    feature_layout: CircularFeatureLayout | None = None,
    track_preset: str | None = None,
) -> list[tuple[float, float, float]]:
    """Build local genome intervals as (start, end, outer_radius_px)."""
    detailed_intervals = _build_feature_radius_intervals(
        feature_dict,
        total_length,
        radius,
        track_ratio,
        profile,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        feature_layout=feature_layout,
        track_preset=track_preset,
    )
    return [(start, end, outer_radius) for start, end, _inner_radius, outer_radius in detailed_intervals]


def _build_inner_feature_radius_intervals(
    feature_dict: dict,
    total_length: int,
    radius: float,
    track_ratio: float,
    profile: CircularRenderProfile,
    *,
    feature_track_ratio_factor_override: float | None = None,
    feature_layout: CircularFeatureLayout | None = None,
    track_preset: str | None = None,
) -> list[tuple[float, float, float]]:
    """Build local genome intervals as (start, end, inner_radius_px)."""
    detailed_intervals = _build_feature_radius_intervals(
        feature_dict,
        total_length,
        radius,
        track_ratio,
        profile,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        feature_layout=feature_layout,
        track_preset=track_preset,
    )
    return [(start, end, inner_radius) for start, end, inner_radius, _outer_radius in detailed_intervals]


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


def _min_inner_feature_radius_at_position(
    intervals: list[tuple[float, float, float]],
    position: float,
    total_length: int,
) -> float:
    """Return the min feature inner radius that covers a genome position."""
    if not intervals or total_length <= 0:
        return 0.0
    position_wrapped = float(position % float(total_length))
    min_inner_radius: float | None = None
    for start_pos, end_pos, inner_radius in intervals:
        if start_pos <= position_wrapped < end_pos:
            if min_inner_radius is None or inner_radius < min_inner_radius:
                min_inner_radius = inner_radius
    return 0.0 if min_inner_radius is None else float(min_inner_radius)


def _label_genome_intervals_for_clearance(label: CircularLabel, total_length: int) -> list[tuple[float, float]]:
    """Return genome intervals covered by the current label width estimate."""
    if total_length <= 0:
        return []

    total_length_f = float(total_length)
    start_x = float(label.get("start_x", 0.0))
    start_y = float(label.get("start_y", 0.0))
    label_radius = math.hypot(start_x, start_y)

    if label_radius <= 1e-9:
        middle = float(label.get("middle", 0.0)) % total_length_f
        return [(middle, middle)]

    width_px = max(0.0, float(label.get("width_px", 0.0)))
    if width_px <= 1e-9:
        anchor_angle = (math.degrees(math.atan2(start_y, start_x)) + 90.0) % 360.0
        anchor_position = (anchor_angle / 360.0) * total_length_f
        return [(anchor_position, anchor_position)]

    span_bp = total_length_f * (width_px / (2.0 * math.pi * label_radius))
    if span_bp >= total_length_f:
        return [(0.0, total_length_f)]

    anchor_angle = (math.degrees(math.atan2(start_y, start_x)) + 90.0) % 360.0
    center_position = (anchor_angle / 360.0) * total_length_f
    half_span = 0.5 * span_bp
    interval_start = center_position - half_span
    interval_end = center_position + half_span

    if interval_start < 0.0:
        return [(interval_start + total_length_f, total_length_f), (0.0, interval_end)]
    if interval_end > total_length_f:
        return [(interval_start, total_length_f), (0.0, interval_end - total_length_f)]
    return [(interval_start, interval_end)]


def _max_outer_feature_radius_for_intervals(
    feature_radius_intervals: list[tuple[float, float, float]],
    label_intervals: list[tuple[float, float]],
    total_length: int,
) -> float:
    """Return max feature outer radius intersecting any label genome interval."""
    if total_length <= 0 or not feature_radius_intervals or not label_intervals:
        return 0.0

    max_outer_radius = 0.0
    for label_start, label_end in label_intervals:
        label_start_f = float(label_start)
        label_end_f = float(label_end)

        # Degenerate interval: treat as a point lookup.
        if label_end_f <= label_start_f + 1e-9:
            max_outer_radius = max(
                max_outer_radius,
                _max_outer_feature_radius_at_position(feature_radius_intervals, label_start_f, total_length),
            )
            continue

        for feature_start, feature_end, feature_outer in feature_radius_intervals:
            if float(feature_start) < label_end_f and float(feature_end) > label_start_f:
                max_outer_radius = max(max_outer_radius, float(feature_outer))
    return max_outer_radius


def _min_inner_feature_radius_for_intervals(
    feature_radius_intervals: list[tuple[float, float, float]],
    label_intervals: list[tuple[float, float]],
    total_length: int,
) -> float:
    """Return min feature inner radius intersecting any label genome interval."""
    if total_length <= 0 or not feature_radius_intervals or not label_intervals:
        return 0.0

    min_inner_radius: float | None = None
    for label_start, label_end in label_intervals:
        label_start_f = float(label_start)
        label_end_f = float(label_end)

        if label_end_f <= label_start_f + 1e-9:
            point_inner = _min_inner_feature_radius_at_position(
                feature_radius_intervals,
                label_start_f,
                total_length,
            )
            if point_inner > 0.0:
                min_inner_radius = (
                    point_inner
                    if min_inner_radius is None
                    else min(min_inner_radius, point_inner)
                )
            continue

        for feature_start, feature_end, feature_inner in feature_radius_intervals:
            if float(feature_start) < label_end_f and float(feature_end) > label_start_f:
                min_inner_radius = (
                    float(feature_inner)
                    if min_inner_radius is None
                    else min(min_inner_radius, float(feature_inner))
                )
    return 0.0 if min_inner_radius is None else min_inner_radius


def _update_outer_label_minimum_radii_against_features(
    labels: list[CircularLabel],
    total_length: int,
    feature_radius_intervals: list[tuple[float, float, float]],
    *,
    text_clearance_px: float,
) -> list[CircularLabel]:
    """Raise minimum outer-label radii so text clears local feature tracks."""
    if not labels or not feature_radius_intervals:
        return labels

    for label in labels:
        label_intervals = _label_genome_intervals_for_clearance(label, total_length)
        if not label_intervals:
            continue

        local_outer_radius = _max_outer_feature_radius_for_intervals(
            feature_radius_intervals,
            label_intervals,
            total_length,
        )
        if local_outer_radius <= 0.0:
            continue

        required_radius = (
            local_outer_radius
            + float(text_clearance_px)
            + OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX
        )
        current_min_radius = float(label.get("min_outer_start_radius_px", 0.0))
        if required_radius > current_min_radius:
            label["min_outer_start_radius_px"] = required_radius

    return labels


def _update_outer_label_minimum_middle_radii_against_features(
    labels: list[CircularLabel],
    total_length: int,
    feature_radius_intervals: list[tuple[float, float, float]],
    *,
    anchor_clearance_px: float,
) -> list[CircularLabel]:
    """Raise leader middle-point radii so elbows stay outside local feature tracks."""
    if not labels or not feature_radius_intervals:
        return labels

    for label in labels:
        label_intervals = _label_genome_intervals_for_clearance(label, total_length)
        if not label_intervals:
            continue

        local_outer_radius = _max_outer_feature_radius_for_intervals(
            feature_radius_intervals,
            label_intervals,
            total_length,
        )
        if local_outer_radius <= 0.0:
            continue

        required_radius = (
            local_outer_radius
            + float(anchor_clearance_px)
            + OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX
        )
        current_min_radius = float(label.get("min_outer_middle_radius_px", 0.0))
        if required_radius > current_min_radius:
            label["min_outer_middle_radius_px"] = required_radius

    return labels


def _update_inner_label_maximum_radii_against_features(
    labels: list[CircularLabel],
    total_length: int,
    feature_radius_intervals: list[tuple[float, float, float]],
    *,
    text_clearance_px: float,
) -> list[CircularLabel]:
    """Lower maximum inner-label radii so text stays inside local feature tracks."""
    if not labels or not feature_radius_intervals:
        return labels

    for label in labels:
        label_intervals = _label_genome_intervals_for_clearance(label, total_length)
        if not label_intervals:
            continue

        local_inner_radius = _min_inner_feature_radius_for_intervals(
            feature_radius_intervals,
            label_intervals,
            total_length,
        )
        if local_inner_radius <= 0.0:
            continue

        allowed_radius = max(
            0.0,
            local_inner_radius
            - float(text_clearance_px)
            - OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX,
        )
        current_max_radius = float(label.get("max_inner_start_radius_px", float("inf")))
        if allowed_radius < current_max_radius:
            label["max_inner_start_radius_px"] = allowed_radius

    return labels


def _update_inner_label_maximum_middle_radii_against_features(
    labels: list[CircularLabel],
    total_length: int,
    feature_radius_intervals: list[tuple[float, float, float]],
    *,
    anchor_clearance_px: float,
) -> list[CircularLabel]:
    """Lower maximum inner-label middle radii so elbows stay inside local feature tracks."""
    if not labels or not feature_radius_intervals:
        return labels

    for label in labels:
        label_intervals = _label_genome_intervals_for_clearance(label, total_length)
        if not label_intervals:
            continue

        local_inner_radius = _min_inner_feature_radius_for_intervals(
            feature_radius_intervals,
            label_intervals,
            total_length,
        )
        if local_inner_radius <= 0.0:
            continue

        allowed_radius = max(
            0.0,
            local_inner_radius
            - float(anchor_clearance_px)
            - OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX,
        )
        current_max_radius = float(label.get("max_inner_middle_radius_px", float("inf")))
        if allowed_radius < current_max_radius:
            label["max_inner_middle_radius_px"] = allowed_radius

    return labels


def _enforce_outer_label_minimum_radius(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
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


def _enforce_outer_label_minimum_middle_radius(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """Push leader middle points outward until they clear required radii."""
    for label in labels:
        required_middle_radius = float(label.get("min_outer_middle_radius_px", 0.0))
        if required_middle_radius <= 0.0:
            continue

        middle_x = float(label.get("middle_x", 0.0))
        middle_y = float(label.get("middle_y", 0.0))
        current_middle_radius = math.hypot(middle_x, middle_y)
        if current_middle_radius >= (required_middle_radius - 1e-6):
            continue

        if current_middle_radius > 1e-9:
            middle_angle_rad = math.atan2(middle_y, middle_x)
        else:
            middle_angle_rad = math.radians(angle_from_middle(float(label["middle"]), total_length))

        label["middle_x"] = required_middle_radius * math.cos(middle_angle_rad)
        label["middle_y"] = required_middle_radius * math.sin(middle_angle_rad)

    return labels


def _enforce_inner_label_maximum_radius(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """Pull inner labels inward until their text bboxes stay inside required radii."""
    for label in labels:
        allowed_outer_radius = float(label.get("max_inner_start_radius_px", float("inf")))
        if not math.isfinite(allowed_outer_radius):
            continue

        start_x = float(label.get("start_x", 0.0))
        start_y = float(label.get("start_y", 0.0))
        current_radius = math.hypot(start_x, start_y)
        if _label_bbox_max_radius(label, total_length, minimum_margin=0.0) <= (allowed_outer_radius + 1e-6):
            continue

        if current_radius > 1e-9:
            angle_rad = math.atan2(start_y, start_x)
        else:
            angle_rad = math.radians(angle_from_middle(float(label["middle"]), total_length))
        current_half = _current_half_from_x(start_x, axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
        candidate_radius = min(current_radius, allowed_outer_radius)

        for _ in range(8):
            new_x = candidate_radius * math.cos(angle_rad)
            new_y = candidate_radius * math.sin(angle_rad)
            if current_half > 0 and new_x < 0:
                new_x = abs(new_x)
            elif current_half < 0 and new_x > 0:
                new_x = -abs(new_x)
            label["start_x"] = float(new_x)
            label["start_y"] = float(new_y)

            bbox_max_radius = _label_bbox_max_radius(label, total_length, minimum_margin=0.0)
            excess = bbox_max_radius - allowed_outer_radius
            if excess <= 1e-3:
                break
            candidate_radius = max(0.0, candidate_radius - (excess + 0.5))

    return labels


def _enforce_inner_label_maximum_middle_radius(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """Pull inner label middle points inward until they stay inside required radii."""
    for label in labels:
        allowed_middle_radius = float(label.get("max_inner_middle_radius_px", float("inf")))
        if not math.isfinite(allowed_middle_radius):
            continue

        middle_x = float(label.get("middle_x", 0.0))
        middle_y = float(label.get("middle_y", 0.0))
        current_middle_radius = math.hypot(middle_x, middle_y)
        if current_middle_radius <= (allowed_middle_radius + 1e-6):
            continue

        if current_middle_radius > 1e-9:
            middle_angle_rad = math.atan2(middle_y, middle_x)
        else:
            middle_angle_rad = math.radians(angle_from_middle(float(label["middle"]), total_length))

        bounded_middle_radius = max(0.0, allowed_middle_radius)
        label["middle_x"] = bounded_middle_radius * math.cos(middle_angle_rad)
        label["middle_y"] = bounded_middle_radius * math.sin(middle_angle_rad)

    return labels


def _stabilize_inner_label_clearance(
    labels: list[CircularLabel],
    total_length: int,
    feature_inner_radius_intervals: list[tuple[float, float, float]],
    *,
    text_clearance_px: float,
    anchor_clearance_px: float,
) -> list[CircularLabel]:
    """Apply iterative inner-label clearance passes against local feature tracks."""
    if not labels or not feature_inner_radius_intervals:
        return labels

    labels = _update_inner_label_maximum_radii_against_features(
        labels,
        total_length,
        feature_inner_radius_intervals,
        text_clearance_px=text_clearance_px,
    )
    labels = _enforce_inner_label_maximum_radius(labels, total_length)
    for _ in range(2):
        labels = _resolve_outer_label_overlaps_with_fixed_radii(labels, total_length)
        labels = _update_inner_label_maximum_radii_against_features(
            labels,
            total_length,
            feature_inner_radius_intervals,
            text_clearance_px=text_clearance_px,
        )
        labels = _enforce_inner_label_maximum_radius(labels, total_length)

    labels = _resolve_outer_label_overlaps_with_fixed_radii(labels, total_length)
    labels = _update_inner_label_maximum_radii_against_features(
        labels,
        total_length,
        feature_inner_radius_intervals,
        text_clearance_px=text_clearance_px,
    )
    labels = _enforce_inner_label_maximum_radius(labels, total_length)
    labels = _update_inner_label_maximum_middle_radii_against_features(
        labels,
        total_length,
        feature_inner_radius_intervals,
        anchor_clearance_px=anchor_clearance_px,
    )
    labels = _enforce_inner_label_maximum_middle_radius(labels, total_length)

    labels = _resolve_outer_label_overlaps_with_fixed_radii(labels, total_length)
    labels = _update_inner_label_maximum_radii_against_features(
        labels,
        total_length,
        feature_inner_radius_intervals,
        text_clearance_px=text_clearance_px,
    )
    labels = _enforce_inner_label_maximum_radius(labels, total_length)
    labels = _update_inner_label_maximum_middle_radii_against_features(
        labels,
        total_length,
        feature_inner_radius_intervals,
        anchor_clearance_px=anchor_clearance_px,
    )
    labels = _enforce_inner_label_maximum_middle_radius(labels, total_length)
    return labels


def _resolve_outer_label_overlaps_with_fixed_radii(
    labels: list[CircularLabel],
    total_length: int,
    *,
    max_angle_shift_deg: float = 40.0,
    step_deg: float = 0.2,
) -> list[CircularLabel]:
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
    label_count = len(labels)
    pair_margins = [[0.0 for _ in range(label_count)] for _ in range(label_count)]
    for idx in range(label_count):
        for jdx in range(idx + 1, label_count):
            pair_margins[idx][jdx] = minimum_bbox_gap_px(labels[idx], labels[jdx])
    max_pair_margin = max((max(row) for row in pair_margins), default=0.0)

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

    def _apply_angle(label: CircularLabel) -> None:
        radius = float(label["_fixed_radius_px"])
        angle_unwrapped = float(label["_angle_unwrapped_fixed"])
        label["start_x"] = radius * math.cos(math.radians(angle_unwrapped))
        label["start_y"] = radius * math.sin(math.radians(angle_unwrapped))

    def _pair_margin(idx: int, jdx: int) -> float:
        if idx < jdx:
            return pair_margins[idx][jdx]
        return pair_margins[jdx][idx]

    def _pair_overlaps(left: CircularLabel, right: CircularLabel, idx: int, jdx: int) -> bool:
        min_gap_px = _pair_margin(idx, jdx)
        return y_overlap(left, right, total_length, min_gap_px) and x_overlap(
            left,
            right,
            minimum_margin=min_gap_px,
        )

    def _pair_angle_bounds(label_idx: int) -> tuple[float, float]:
        target_angle = float(labels[label_idx]["_target_angle_unwrapped_fixed"])
        low_bound = target_angle - max_angle_shift_deg
        high_bound = target_angle + max_angle_shift_deg
        if label_idx > 0:
            low_bound = max(low_bound, float(labels[label_idx - 1]["_angle_unwrapped_fixed"]) + min_order_gap_deg)
        if label_idx + 1 < len(labels):
            high_bound = min(high_bound, float(labels[label_idx + 1]["_angle_unwrapped_fixed"]) - min_order_gap_deg)
        return low_bound, high_bound

    def _separate_pair(idx: int, jdx: int) -> bool:
        left = labels[idx]
        right = labels[jdx]
        pair_changed = False

        for _ in range(max_pair_steps):
            if not _pair_overlaps(left, right, idx, jdx):
                break

            left_low, left_high = _pair_angle_bounds(idx)
            right_low, right_high = _pair_angle_bounds(jdx)
            new_left = max(left_low, float(left["_angle_unwrapped_fixed"]) - (0.5 * step_deg))
            new_right = min(right_high, float(right["_angle_unwrapped_fixed"]) + (0.5 * step_deg))

            # Preserve strict ordering with labels between idx and jdx.
            if idx + 1 < jdx:
                new_left = min(new_left, float(labels[idx + 1]["_angle_unwrapped_fixed"]) - min_order_gap_deg)
                new_right = max(new_right, float(labels[jdx - 1]["_angle_unwrapped_fixed"]) + min_order_gap_deg)

            if (
                new_left < (left_low - 1e-9)
                or new_left > (left_high + 1e-9)
                or new_right < (right_low - 1e-9)
                or new_right > (right_high + 1e-9)
            ):
                break
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
            pair_changed = True

        return pair_changed

    for label in labels:
        _apply_angle(label)

    for _ in range(120):
        changed = False
        for idx in range(len(labels) - 1):
            left = labels[idx]
            right = labels[idx + 1]
            if not _pair_overlaps(left, right, idx, idx + 1):
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

                if not _pair_overlaps(left, right, idx, idx + 1):
                    break

        # Adjacent passes can still miss sparse non-adjacent bbox collisions.
        for idx, jdx in sorted(_candidate_label_pairs(labels, total_length, padding=max_pair_margin)):
            if jdx <= idx + 1:
                continue
            if not _pair_overlaps(labels[idx], labels[jdx], idx, jdx):
                continue
            if _separate_pair(idx, jdx):
                changed = True

        if not changed:
            break

    for label in labels:
        label.pop("_fixed_radius_px", None)
        label.pop("_target_angle_unwrapped_fixed", None)
        label.pop("_angle_unwrapped_fixed", None)

    return sort_labels(labels)


def _angle_deviation_score(labels: list[CircularLabel], total_length: int) -> tuple[float, float]:
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


def _placement_score(labels: list[CircularLabel], total_length: int) -> tuple[int, int, int, float, float, float, float, float]:
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


def _dense_inner_candidate_score(
    placement_score: tuple[int, int, int, float, float, float, float, float],
    leader_collision_count: int,
) -> tuple[int, int, int, int, float, float, float, float, float]:
    """
    Return comparison score for dense inner candidates.

    Dense inner clusters are readability-sensitive around the lower axis, so we
    prioritize hemisphere consistency and leader-line collisions ahead of
    minimum-gap strictness.
    """
    (
        overlap_plain,
        overlap_min_gap,
        hemisphere_mismatch_count,
        hemisphere_mismatch_weight,
        leader_sum,
        leader_max,
        sum_delta,
        max_delta,
    ) = placement_score
    return (
        overlap_plain,
        hemisphere_mismatch_count,
        int(leader_collision_count),
        overlap_min_gap,
        hemisphere_mismatch_weight,
        leader_sum,
        leader_max,
        sum_delta,
        max_delta,
    )


def _bounds_overlap(bounds1: tuple[float, float], bounds2: tuple[float, float]) -> bool:
    """Return whether two (min, max) label bounds overlap; touching counts."""
    if bounds1[0] < bounds2[0]:
        return bounds1[1] >= bounds2[0]
    return bounds2[1] >= bounds1[0]


def y_overlap(label1: CircularLabel, label2: CircularLabel, total_len: int, minimum_margin: float) -> bool:
    return _bounds_overlap(
        _label_y_bounds(label1, total_len, minimum_margin),
        _label_y_bounds(label2, total_len, minimum_margin),
    )


def x_overlap(label1: CircularLabel, label2: CircularLabel, minimum_margin: float = 1.0) -> bool:
    return _bounds_overlap(_label_x_bounds(label1, minimum_margin), _label_x_bounds(label2, minimum_margin))


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
    label: CircularLabel,
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
    labels: list[CircularLabel],
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
        weight = _label_hemisphere_mismatch_weight(label, total_length, axis_neutral_deg)
        if weight is None:
            continue
        mismatch_count += 1
        weighted_score += weight

    return mismatch_count, weighted_score


def _label_hemisphere_mismatch_weight(
    label: CircularLabel,
    total_length: int,
    axis_neutral_deg: float = HEMISPHERE_AXIS_NEUTRAL_DEG,
) -> float | None:
    """Return the mismatch weight of one label on the wrong half, or None when it is not."""
    preferred_half = _preferred_half_for_label(label, total_length, axis_neutral_deg)
    if preferred_half == 0:
        return None

    current_half = _current_half_from_x(float(label.get("start_x", 0.0)), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
    if current_half in (0, preferred_half):
        return None

    target_angle_unwrapped = float(
        label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
    )
    target_angle_mod = target_angle_unwrapped % 360.0
    distance_from_axis = _vertical_axis_distance_deg(target_angle_mod)
    return max(0.0, distance_from_axis - axis_neutral_deg)


def place_labels_on_arc_fc(
    labels: list[CircularLabel],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    start_angle: float,
    end_angle: float,
    total_length: int,
    spacing_px: float | None = None,
) -> list[CircularLabel]:
    requested_spacing_px = 2.0 if spacing_px is None else max(0.0, float(spacing_px))
    resolved_spacing_px = _effective_group_spacing_px(
        requested_spacing_px,
        legacy_spacing_px=2.0,
    )

    def check_overlap(label1: CircularLabel, label2: CircularLabel, total_length: int, margin: float) -> bool:
        effective_margin = minimum_bbox_gap_px(label1, label2, base_margin_px=max(margin, resolved_spacing_px))
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
        )

    if not labels:
        return []

    for label in labels:
        label["_spacing_px"] = resolved_spacing_px
        label["_enforce_spacing_px"] = True

    rearranged_labels: list[CircularLabel] = []
    labels = sort_labels(labels)
    angle_step = 0.25
    max_shift = 25.0
    max_steps = int(max_shift / angle_step)
    min_order_gap = 0.05
    prev_angle_unwrapped = -float("inf")

    for label in labels:
        target_angle_unwrapped = angle_from_middle_unwrapped(label["middle"], total_length)
        label["target_angle_unwrapped"] = target_angle_unwrapped
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


def sort_labels(labels: Iterable[CircularLabel]) -> list[CircularLabel]:
    return sorted(labels, key=lambda x: x["middle"])


def _legacy_place_labels_on_arc_fc(
    labels: list[CircularLabel],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    start_angle: float,
    end_angle: float,
    total_length: int,
    *,
    spacing_px: float,
) -> list[CircularLabel]:
    """
    Legacy placement strategy (from pre-93dc098 implementation).

    This remains as a fallback for dense cases where the new optimizer leaves
    unresolved overlaps. ``spacing_px`` is the circular label spacing: the
    vertical gap kept between neighbouring labels.
    """

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

        x, y = calculate_coordinates(center_x, center_y, x_radius, y_radius, new_angle, label["middle"], total_length)
        # The previous label stays put while this one steps away from it. The step loop makes
        # no call per step (about 0.4 M steps per search, OV-870): the label's box offsets are
        # looked up once, and the coordinates and the interval test (`calculate_coordinates`,
        # `_bounds_overlap`) are inlined.
        prev_min_y, prev_max_y = _label_y_bounds(rearranged_labels[-1], total_length, spacing_px)
        prev_min_x, prev_max_x = _label_x_bounds(rearranged_labels[-1], 2)
        is_outer = label.get("is_inner", False) is False
        low_y, low_y_margin, high_y, high_y_margin = _label_y_offsets(
            float(label["middle"]), int(total_length), float(label["height_px"]), float(spacing_px), is_outer
        )
        right, left = _label_x_offsets(float(label["width_px"]), 2.0, is_outer)
        while True:
            min_y, max_y = (y + low_y) + low_y_margin, (y + high_y) + high_y_margin
            if not (max_y >= prev_min_y if min_y < prev_min_y else prev_max_y >= min_y):
                break
            low_x, low_x_margin, high_x, high_x_margin = right if x > 0 else left
            min_x, max_x = (x + low_x) + low_x_margin, (x + high_x) + high_x_margin
            if not (max_x >= prev_min_x if min_x < prev_min_x else prev_max_x >= min_x):
                break
            new_angle += 0.01
            radians = math.radians(new_angle)
            y = center_y + y_radius * math.sin(radians)
            x = center_x + x_radius * math.cos(radians)
        label["start_x"], label["start_y"] = x, y

        rearranged_labels.append(label)
        current_angle = new_angle

    return rearranged_labels


def _legacy_improved_label_placement_fc(
    labels: list[CircularLabel],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    feature_radius: float,
    total_length: int,
    start_angle: float,
    end_angle: float,
    *,
    y_margin: float,
    max_iterations: int = 10000,
) -> list[CircularLabel]:
    """
    Legacy local optimizer (from pre-93dc098 implementation).

    This solver is intentionally conservative and is used as a fallback to
    improve dense-label robustness.
    """

    def feature_angle_scorer(label: CircularLabel) -> Callable[[float, float], float]:
        """Return the angle in degrees at the origin between the label's feature middle and a point."""
        v1_x = float(label["feature_middle_x"])
        v1_y = float(label["feature_middle_y"])
        mag1 = math.sqrt(v1_x**2 + v1_y**2)

        def angle_to(x: float, y: float) -> float:
            v2_x = float(x)
            v2_y = float(y)
            dot_product = v1_x * v2_x + v1_y * v2_y
            mag2 = math.sqrt(v2_x**2 + v2_y**2)
            if mag1 == 0 or mag2 == 0:
                return 0.0
            cos_angle = dot_product / (mag1 * mag2)
            return math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))

        return angle_to

    def check_overlap(label1: CircularLabel, label2: CircularLabel, total_length: int) -> bool:
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

            feature_angle = feature_angle_scorer(label)
            current_score = feature_angle(label["start_x"], label["start_y"])

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
                test_x_plus, test_y_plus = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, test_angle_plus, label["middle"], total_length
                )
                score_plus = feature_angle(test_x_plus, test_y_plus)

                test_angle_minus = current_angle - 1.0
                test_x_minus, test_y_minus = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, test_angle_minus, label["middle"], total_length
                )
                score_minus = feature_angle(test_x_minus, test_y_minus)
                direction = 1 if abs(score_plus) < abs(score_minus) else -1

            # The labels this move must not run into; they stay put while this label steps.
            if idx == 0:
                prev_label, next_label = labels[reverse_idx - 1], labels[0]
            elif 0 < idx < len(labels) - 1:
                prev_label, next_label = labels[reverse_idx - 1], labels[reverse_idx + 1]
            else:
                prev_label, next_label = labels[-1], labels[reverse_idx + 1]
            if overlaps_prev:
                fixed_labels = [next_label]
            elif overlaps_next:
                fixed_labels = [prev_label]
            else:
                fixed_labels = [prev_label, next_label]
            fixed_bounds = [
                (*_label_y_bounds(other, total_length, y_margin), *_label_x_bounds(other, 1)) for other in fixed_labels
            ]
            # The step loop makes no call per step (about 0.5 M steps per search, OV-870): the
            # label's box offsets are looked up once, and the coordinates, the feature angle and
            # the interval test (`calculate_coordinates`, `feature_angle`, `_bounds_overlap`) are inlined.
            is_outer = label.get("is_inner", False) is False
            low_y, low_y_margin, high_y, high_y_margin = _label_y_offsets(
                float(label["middle"]), int(total_length), float(label["height_px"]), float(y_margin), is_outer
            )
            right, left = _label_x_offsets(float(label["width_px"]), 1.0, is_outer)
            feature_x, feature_y = float(label["feature_middle_x"]), float(label["feature_middle_y"])
            feature_norm = math.sqrt(feature_x**2 + feature_y**2)
            while True:
                new_angle = current_angle + direction * 0.05
                radians = math.radians(new_angle)
                new_y = center_y + y_radius * math.sin(radians)
                new_x = center_x + x_radius * math.cos(radians)
                norm = math.sqrt(new_x**2 + new_y**2)
                if feature_norm == 0 or norm == 0:
                    new_score = 0.0
                else:
                    cos_angle = (feature_x * new_x + feature_y * new_y) / (feature_norm * norm)
                    new_score = math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))
                if not abs(new_score) <= abs(current_score):
                    break

                min_y, max_y = (new_y + low_y) + low_y_margin, (new_y + high_y) + high_y_margin
                creates_new_overlap = False
                for other_min_y, other_max_y, other_min_x, other_max_x in fixed_bounds:
                    if not (max_y >= other_min_y if min_y < other_min_y else other_max_y >= min_y):
                        continue
                    low_x, low_x_margin, high_x, high_x_margin = right if new_x > 0 else left
                    min_x, max_x = (new_x + low_x) + low_x_margin, (new_x + high_x) + high_x_margin
                    if max_x >= other_min_x if min_x < other_min_x else other_max_x >= min_x:
                        creates_new_overlap = True
                        break
                if creates_new_overlap:
                    break
                label["start_x"], label["start_y"] = new_x, new_y
                current_angle = new_angle
                current_score = new_score
                changes_made = True
        if not changes_made:
            break

    return labels


def improved_label_placement_fc(
    labels: list[CircularLabel],
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
    spacing_px: float | None = None,
):
    labels = sort_labels(labels)
    if not labels:
        return []

    requested_spacing_px = 1.5 if spacing_px is None else max(0.0, float(spacing_px))
    resolved_spacing_px = _effective_group_spacing_px(
        requested_spacing_px,
        legacy_spacing_px=1.5,
    )
    for label in labels:
        label["_spacing_px"] = resolved_spacing_px
        label["_enforce_spacing_px"] = True

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
        prev_target_unwrapped = target_angle_unwrapped

    pair_margins: list[list[float]] = []
    max_pair_margin = 0.0

    def move_label(label: CircularLabel, angle: float) -> tuple[float, float]:
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def check_overlap(
        label1: CircularLabel, label2: CircularLabel, total_length: int, effective_margin: float
    ) -> bool:
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
        )

    def pair_margin(idx1: int, idx2: int) -> float:
        if idx1 < idx2:
            return pair_margins[idx1][idx2]
        return pair_margins[idx2][idx1]

    def overlap_indices(label_idx: int, *, full_scan: bool) -> set[int]:
        n_labels = len(labels)
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
        return checked_indices

    def overlap_count(
        label_idx,
        candidate_label=None,
        *,
        full_scan: bool = False,
        indices: Iterable[int] | None = None,
        cap: int | None = None,
    ):
        """Count overlaps with the checked labels (`indices`: a subset holding every one that can overlap)."""
        current = candidate_label if candidate_label is not None else labels[label_idx]
        n_labels = len(labels)
        if n_labels <= 1:
            return 0
        overlaps = 0
        for other_idx in overlap_indices(label_idx, full_scan=full_scan) if indices is None else indices:
            if other_idx == label_idx:
                continue
            if check_overlap(current, labels[other_idx], total_length, pair_margin(label_idx, other_idx)):
                overlaps += 1
                if cap is not None and overlaps > cap:
                    break
        return overlaps

    def update_label_coordinates(label: CircularLabel, angle_unwrapped: float) -> None:
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
        angles = [float(label["angle_unwrapped"]) for label in labels]
        if not _keep_cyclic_order(angles, min_gap):
            return False
        for label, angle in zip(labels, angles):
            if angle != label["angle_unwrapped"]:
                update_label_coordinates(label, angle)
        return True

    def overlapping_pairs() -> list[tuple[int, int]]:
        return [
            (i, j)
            for i, j in sorted(_candidate_label_pairs(labels, total_length, padding=max_pair_margin))
            if check_overlap(labels[i], labels[j], total_length, pair_margin(i, j))
        ]

    def count_remaining_overlaps() -> int:
        return len(overlapping_pairs())

    def find_overlap_components() -> list[list[int]]:
        n_labels = len(labels)
        adjacency: dict[int, set[int]] = {i: set() for i in range(n_labels)}
        for i, j in overlapping_pairs():
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
        dense_inner_relax = (
            n_labels >= DENSE_INNER_RELAX_MIN_LABELS
            and all(bool(label.get("is_inner", False)) for label in labels)
        )
        if dense_inner_relax:
            # Inner dense sets are the primary hot path; a coarser search grid keeps
            # near-identical placements while avoiding runaway overlap evaluations.
            cluster_step_limit = min(cluster_step_limit, 10)
        changed = False

        for component in components:
            component_center = component[len(component) // 2]
            cluster_padding = min(6, max(2, len(component) + 1))
            cluster_start = max(0, component[0] - cluster_padding)
            cluster_end = min(n_labels - 1, component[-1] + cluster_padding)
            cluster_indices = list(range(cluster_start, cluster_end + 1))

            # Candidates are counted from the current pairs: only pairs with a label that moves can change.
            current_overlaps = _MovedLabelOverlaps(
                labels,
                overlapping_pairs(),
                lambda label_i, label_j, i, j: check_overlap(label_i, label_j, total_length, pair_margin(i, j)),
                _reach_pad(max_pair_margin),
            )
            current_angles = [float(label["angle_unwrapped"]) for label in labels]
            # A label whose stored position is not the one its angle gives moves in every candidate.
            at_angle = [
                move_label(label, angle) == (label["start_x"], label["start_y"])
                for label, angle in zip(labels, current_angles)
            ]
            current_score = (
                current_overlaps.count,
                abs(labels[component_center]["angle_unwrapped"] - labels[component_center]["target_angle_unwrapped"]),
                sum(abs(labels[idx]["angle_unwrapped"] - labels[idx]["target_angle_unwrapped"]) for idx in cluster_indices),
            )
            best_score = current_score
            best_ordered: list[float] | None = None

            center_candidates = [component_center]
            if not dense_inner_relax and len(component) >= 2:
                center_candidates.append(component[0])
                center_candidates.append(component[-1])
            center_candidates = sorted(set(center_candidates))
            center_shift_limit_steps = min(int(shift_limit_deg / angle_step), 48)
            center_shift_step = 8 if dense_inner_relax else 4
            if dense_inner_relax:
                center_shift_limit_steps = min(center_shift_limit_steps, 16)
            center_shift_steps = range(
                -center_shift_limit_steps,
                center_shift_limit_steps + 1,
                center_shift_step,
            )
            center_delta_cap = min(float(shift_limit_deg), RELAX_CENTER_DELTA_CAP_DEG)

            for center_idx in center_candidates:
                for step in range(0, cluster_step_limit + 1):
                    spread = step * angle_step
                    for shift_step in center_shift_steps:
                        center_shift = shift_step * angle_step
                        ordered = list(current_angles)

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
                        _keep_cyclic_order(ordered, min_order_gap)

                        # Respect per-label shift budget.
                        violates_budget = False
                        for idx in range(n_labels):
                            target = float(labels[idx]["target_angle_unwrapped"])
                            if abs(ordered[idx] - target) > shift_limit_deg + 1e-6:
                                violates_budget = True
                                break
                        if violates_budget:
                            continue

                        center_delta = abs(
                            ordered[component_center] - labels[component_center]["target_angle_unwrapped"]
                        )
                        if center_delta > center_delta_cap:
                            continue

                        moved: dict[int, CircularLabel] = {}
                        for idx in range(n_labels):
                            if ordered[idx] != current_angles[idx] or not at_angle[idx]:
                                moved_label = labels[idx].copy()
                                moved_label["angle_unwrapped"] = ordered[idx]
                                moved_label["start_x"], moved_label["start_y"] = move_label(moved_label, ordered[idx])
                                moved[idx] = moved_label
                        overlap_total = current_overlaps.count_after_move(moved, cap=best_score[0])
                        if overlap_total > best_score[0]:
                            continue

                        candidate_score = (
                            overlap_total,
                            center_delta,
                            sum(
                                abs(ordered[idx] - labels[idx]["target_angle_unwrapped"])
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
        dense_inner_set = (
            len(labels) >= DENSE_INNER_RELAX_MIN_LABELS
            and all(bool(label.get("is_inner", False)) for label in labels)
        )
        for _ in range(max_iterations):
            changes_made = False
            if enforce_order(min_order_gap):
                changes_made = True
            angles = [float(label["angle_unwrapped"]) for label in labels]

            iteration_order = sorted(range(len(labels)), key=local_density_weight)
            # Peer pruning (output-preserving): check_overlap is false for labels whose reach boxes are disjoint.
            reach_pad = _reach_pad(max_pair_margin)
            reach_boxes = [_label_reach_box(label, reach_pad) for label in labels]
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

                order_lower, order_upper = _cyclic_order_bounds(angles, label_idx, min_order_gap)
                lower_bound = max(target_angle - shift_limit_deg, order_lower)
                upper_bound = min(target_angle + shift_limit_deg, order_upper)
                if lower_bound > upper_bound:
                    lower_bound, upper_bound = order_lower, order_upper
                    if lower_bound > upper_bound:
                        continue

                anchor_weight = local_density_weight(label_idx)
                best_angle = current_angle
                best_x = label["start_x"]
                best_y = label["start_y"]
                best_overlaps = current_overlaps
                best_target_delta = current_target_delta
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

                step_positions: list[list[tuple[float, float, float]]] = []
                for step in range(max_steps + 1):
                    if step == 0:
                        offsets = [0.0]
                    else:
                        offset = step * angle_step
                        offsets = [offset, -offset]
                    positions: list[tuple[float, float, float]] = []
                    for offset in offsets:
                        candidate_angle = target_angle + offset
                        if not (lower_bound <= candidate_angle <= upper_bound):
                            continue
                        positions.append((candidate_angle, *move_label(label, candidate_angle)))
                    step_positions.append(positions)
                sweep_anchors = [(item[1], item[2]) for positions in step_positions for item in positions]
                if not sweep_anchors:
                    continue
                peer_indices = touching_aabb_indices(
                    reach_boxes,
                    _label_reach_box(label, reach_pad, anchors=sweep_anchors),
                    skip=label_idx,
                    candidates=overlap_indices(label_idx, full_scan=full_scan_checks),
                )

                search_complete = False
                for positions in step_positions:
                    for candidate_angle, candidate_x, candidate_y in positions:
                        label_copy = label.copy()
                        label_copy["start_x"] = candidate_x
                        label_copy["start_y"] = candidate_y

                        # A candidate with more overlaps than the best, or a worse (overlaps, weighted
                        # delta) prefix, cannot win; the rest of its score is then not needed.
                        candidate_overlaps = overlap_count(
                            label_idx, label_copy, indices=peer_indices, cap=best_score[0]
                        )
                        if candidate_overlaps > best_score[0]:
                            continue
                        candidate_target_delta = abs(candidate_angle - target_angle)
                        candidate_move = abs(candidate_angle - current_angle)
                        candidate_middle_penalty = middle_proximity_penalty(label_idx, candidate_angle)
                        candidate_weighted_delta = (
                            anchor_weight * candidate_target_delta + proximity_penalty_factor * candidate_middle_penalty
                        )
                        if (candidate_overlaps, candidate_weighted_delta) > best_score[:2]:
                            continue
                        candidate_leader_length = _leader_length_px(label_copy, total_length)
                        candidate_score = (
                            candidate_overlaps,
                            candidate_weighted_delta,
                            candidate_leader_length,
                            candidate_move,
                        )

                        if candidate_score < best_score:
                            best_angle = candidate_angle
                            best_x = candidate_x
                            best_y = candidate_y
                            best_overlaps = candidate_overlaps
                            best_target_delta = candidate_target_delta
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
                    angles[label_idx] = best_angle
                    reach_boxes[label_idx] = _label_reach_box(label, reach_pad)
                    changes_made = True
            if not changes_made:
                break
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            relax_passes = 1 if dense_inner_set else 3
            for _ in range(relax_passes):
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
    base_margin = max(float(y_margin), resolved_spacing_px)
    for i in range(n_labels):
        for j in range(i + 1, n_labels):
            pair_margins[i][j] = minimum_bbox_gap_px(labels[i], labels[j], base_margin_px=base_margin)
    max_pair_margin = max((max(row) for row in pair_margins), default=0.0)

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
    dense_inner_set = (
        len(labels) >= DENSE_INNER_RELAX_MIN_LABELS
        and all(bool(label.get("is_inner", False)) for label in labels)
    )
    allow_full_scan_checks = (len(labels) <= FULL_SCAN_LABEL_LIMIT) and (not dense_inner_set)
    if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS:
        for extra_shift in (35.0, 45.0, 60.0):
            if extra_shift > shift_schedule[-1]:
                shift_schedule.append(extra_shift)

    for phase_idx, shift_limit in enumerate(shift_schedule):
        run_full_scan = (
            allow_full_scan_checks
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
                full_scan_checks=allow_full_scan_checks,
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
            full_scan_checks=allow_full_scan_checks,
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
            full_scan_checks=allow_full_scan_checks,
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


def _cyclic_order_bounds(angles: Sequence[float], idx: int, min_gap: float) -> tuple[float, float]:
    """The open angle interval that keeps label ``idx`` between its neighbours in feature order.

    ``angles`` are the labels' unwrapped angles in feature order. The first and last labels are
    neighbours across the origin, so no label may pass another by wrapping around the circle.
    """
    last = len(angles) - 1
    lower = float(angles[idx - 1]) + min_gap if idx > 0 else float(angles[last]) - 360.0 + min_gap
    upper = float(angles[idx + 1]) - min_gap if idx < last else float(angles[0]) + 360.0 - min_gap
    return lower, upper


def _keep_cyclic_order(angles: list[float], min_gap: float) -> bool:
    """Move ``angles`` in place into their `_cyclic_order_bounds`, keeping the first; return whether any moved.

    A forward pass lifts each label above the one before it; a backward pass then lowers the labels
    that passed the first label's place one turn later, so the chain spans less than a full circle.
    The backward pass ignores overshoots below 1e-9 degrees: ``(a + gap) - gap`` can round below ``a``,
    and lowering by that rounding would undo the forward pass on every call.
    """
    changed = False
    last = len(angles) - 1
    for idx in range(1, last + 1):
        lower, _ = _cyclic_order_bounds(angles, idx, min_gap)
        if angles[idx] < lower:
            angles[idx] = lower
            changed = True
    for idx in range(last, 0, -1):
        _, upper = _cyclic_order_bounds(angles, idx, min_gap)
        if angles[idx] > upper + 1e-9:
            angles[idx] = upper
            changed = True
    return changed


def _derive_monotonic_unwrapped_angles(labels: list[CircularLabel], total_length: int) -> list[float]:
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


def _neighbor_min_gap_overlap_count(
    labels: list[CircularLabel],
    label_idx: int,
    candidate_label: CircularLabel,
    total_length: int,
) -> int:
    """Count min-gap overlaps for immediate neighbors only."""
    overlap_count = 0
    for neighbor_idx in (label_idx - 1, label_idx + 1):
        if neighbor_idx < 0 or neighbor_idx >= len(labels):
            continue
        neighbor_label = labels[neighbor_idx]
        min_gap_px = minimum_bbox_gap_px(candidate_label, neighbor_label)
        if y_overlap(candidate_label, neighbor_label, total_length, min_gap_px) and x_overlap(
            candidate_label,
            neighbor_label,
            minimum_margin=min_gap_px,
        ):
            overlap_count += 1
    return overlap_count


def _inner_rebalance_metrics(
    labels: list[CircularLabel],
    total_length: int,
    plain_overlap: int | None = None,
) -> tuple[int, int, int, float, int, int]:
    """Return baseline/candidate quality metrics for strict-order inner rebalance.

    `plain_overlap`, when given, is the plain overlap count of these labels.
    """
    if plain_overlap is None:
        plain_overlap = _count_label_overlaps(labels, total_length, use_min_gap=False)
    min_gap_overlap = _count_label_overlaps(labels, total_length, use_min_gap=True)
    mismatch_count, mismatch_weighted = _hemisphere_mismatch_metrics(labels, total_length)
    labels_with_leaders = _assign_leader_start_points([label.copy() for label in labels], total_length)
    leader_bbox_collisions = _count_label_leader_line_collisions(
        labels_with_leaders,
        total_length,
        margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
    )
    leader_line_intersections = _count_leader_line_intersections(labels_with_leaders)
    return (
        plain_overlap,
        min_gap_overlap,
        mismatch_count,
        mismatch_weighted,
        leader_bbox_collisions,
        leader_line_intersections,
    )


def _rebalance_inner_labels_strict_order(labels: list[CircularLabel], total_length: int) -> list[CircularLabel]:
    """
    Rebalance inner labels without overtaking, using suffix block pushes only.

    This pass preserves strict middle-order monotonicity and only accepts
    candidates that do not regress plain overlaps, leader/bbox collisions, or
    leader-line intersections.
    """
    if len(labels) < 2:
        return labels

    ordered_labels = sort_labels([label.copy() for label in labels])
    current_angles = _derive_monotonic_unwrapped_angles(ordered_labels, total_length)
    radii: list[float] = []
    for idx, label in enumerate(ordered_labels):
        current_angles[idx] = float(current_angles[idx])
        label["angle_unwrapped"] = current_angles[idx]
        radii.append(math.hypot(float(label["start_x"]), float(label["start_y"])))

    min_order_gap = 0.05
    max_passes = 6
    metrics_cache: dict[
        tuple[tuple[float, float, float, float, float, float, float, float, float, bool, bool], ...],
        tuple[int, int, int, float, int, int],
    ] = {}
    plain_overlap_cache: dict[
        tuple[tuple[float, float, float, float, float, float, float, float, float, bool, bool], ...],
        int,
    ] = {}
    mismatch_cache: dict[
        tuple[tuple[float, float, float, float, float, float, float, float, float, bool, bool], ...],
        tuple[int, float],
    ] = {}
    delta_schedule = (
        0.25,
        0.5,
        0.75,
        1.0,
        1.5,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        8.0,
        10.0,
        12.0,
        14.0,
        16.0,
        18.0,
        20.0,
        22.0,
        24.0,
        26.0,
        28.0,
        30.0,
    )

    def _metrics_cache_key(
        label_list: list[CircularLabel],
    ) -> tuple[tuple[float, float, float, float, float, float, float, float, float, bool, bool], ...]:
        return tuple(
            (
                float(label.get("middle", 0.0)),
                float(label.get("start_x", 0.0)),
                float(label.get("start_y", 0.0)),
                float(label.get("middle_x", 0.0)),
                float(label.get("middle_y", 0.0)),
                float(label.get("width_px", 0.0)),
                float(label.get("height_px", 0.0)),
                float(label.get("feature_anchor_x", label.get("feature_middle_x", 0.0))),
                float(label.get("feature_anchor_y", label.get("feature_middle_y", 0.0))),
                bool(label.get("is_inner", False)),
                ("min_outer_start_radius_px" in label),
            )
            for label in label_list
        )

    def _inner_rebalance_metrics_cached(label_list: list[CircularLabel]) -> tuple[int, int, int, float, int, int]:
        cache_key = _metrics_cache_key(label_list)
        cached = metrics_cache.get(cache_key)
        if cached is None:
            cached = _inner_rebalance_metrics(label_list, total_length, plain_overlap_cache.get(cache_key))
            metrics_cache[cache_key] = cached
        return cached

    def _plain_overlap_cached(label_list: list[CircularLabel]) -> int:
        cache_key = _metrics_cache_key(label_list)
        cached = plain_overlap_cache.get(cache_key)
        if cached is None:
            cached = _count_label_overlaps(label_list, total_length, use_min_gap=False)
            plain_overlap_cache[cache_key] = cached
        return cached

    def _mismatch_metrics_cached(label_list: list[CircularLabel]) -> tuple[int, float]:
        cache_key = _metrics_cache_key(label_list)
        cached = mismatch_cache.get(cache_key)
        if cached is None:
            cached = _hemisphere_mismatch_metrics(label_list, total_length)
            mismatch_cache[cache_key] = cached
        return cached

    def _mismatch_indices(label_list: list[CircularLabel]) -> list[int]:
        mismatch_indices: list[int] = []
        for idx, label in enumerate(label_list):
            preferred_half = _preferred_half_for_label(label, total_length)
            if preferred_half == 0:
                continue
            current_half = _current_half_from_x(float(label["start_x"]), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
            if current_half not in (0, preferred_half):
                mismatch_indices.append(idx)
        return mismatch_indices

    def _mismatch_priority(label: CircularLabel) -> tuple[float, float]:
        target_angle_unwrapped = float(
            label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
        )
        axis_distance = _vertical_axis_distance_deg(target_angle_unwrapped % 360.0)
        return axis_distance, abs(float(label["start_x"]))

    def _apply_block_shift(
        block_start: int,
        block_end: int,
        shift_deg: float,
    ) -> tuple[list[CircularLabel], list[float]]:
        candidate_labels = [item.copy() for item in ordered_labels]
        candidate_angles = current_angles.copy()
        for move_idx in range(block_start, block_end + 1):
            moved_angle = candidate_angles[move_idx] + shift_deg
            candidate_angles[move_idx] = moved_angle
            candidate_labels[move_idx]["angle_unwrapped"] = moved_angle
            radius = radii[move_idx]
            candidate_labels[move_idx]["start_x"] = radius * math.cos(math.radians(moved_angle))
            candidate_labels[move_idx]["start_y"] = radius * math.sin(math.radians(moved_angle))
        return candidate_labels, candidate_angles

    def _half_progress_non_worsening(
        before_labels: list[CircularLabel],
        after_labels: list[CircularLabel],
        block_start: int,
        block_end: int,
    ) -> bool:
        for move_idx in range(block_start, block_end + 1):
            preferred_half = _preferred_half_for_label(before_labels[move_idx], total_length)
            if preferred_half == 0:
                continue
            before_x = float(before_labels[move_idx]["start_x"])
            after_x = float(after_labels[move_idx]["start_x"])
            if (preferred_half * after_x) < ((preferred_half * before_x) - 1e-6):
                return False
        return True

    def _is_strict_no_overtake(label_list: list[CircularLabel]) -> bool:
        ordered = sort_labels(label_list)
        prev_unwrapped: float | None = None
        for label in ordered:
            angle = math.atan2(float(label["start_y"]), float(label["start_x"]))
            target = (2.0 * math.pi * (float(label["middle"]) / float(total_length))) - (0.5 * math.pi)
            angle_unwrapped = angle + (2.0 * math.pi) * round((target - angle) / (2.0 * math.pi))
            if prev_unwrapped is not None and angle_unwrapped <= prev_unwrapped:
                return False
            prev_unwrapped = angle_unwrapped
        return True

    for _ in range(max_passes):
        (
            baseline_plain_overlap,
            _baseline_min_gap_overlap,
            baseline_mismatch_count,
            _baseline_mismatch_weighted,
            baseline_leader_bbox_collisions,
            baseline_line_intersections,
        ) = _inner_rebalance_metrics_cached(ordered_labels)

        if baseline_mismatch_count <= 3 and baseline_line_intersections == 0:
            break

        mismatch_indices = _mismatch_indices(ordered_labels)
        if mismatch_indices:
            label_idx = max(mismatch_indices, key=lambda idx: _mismatch_priority(ordered_labels[idx]))
        else:
            labels_with_leaders = _assign_leader_start_points([label.copy() for label in ordered_labels], total_length)
            local_line_scores = [
                (_count_local_leader_line_intersections(labels_with_leaders, idx), idx)
                for idx in range(len(labels_with_leaders))
            ]
            local_line_scores.sort(reverse=True)
            if not local_line_scores or local_line_scores[0][0] <= 0:
                break
            label_idx = local_line_scores[0][1]

        focus_label = ordered_labels[label_idx]
        preferred_half = _preferred_half_for_label(focus_label, total_length)
        current_half = _current_half_from_x(float(focus_label["start_x"]), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
        current_half_bad = int(preferred_half != 0 and current_half not in (0, preferred_half))
        target_angle_unwrapped = float(
            focus_label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(focus_label["middle"]), total_length))
        )
        target_angle_unwrapped = normalize_angle_near_reference(target_angle_unwrapped, current_angles[label_idx])
        delta_to_target = target_angle_unwrapped - current_angles[label_idx]
        if abs(delta_to_target) <= 1e-9:
            if preferred_half < 0:
                primary_direction = -1.0
            elif preferred_half > 0:
                primary_direction = 1.0
            else:
                primary_direction = 1.0
        else:
            primary_direction = 1.0 if delta_to_target > 0.0 else -1.0
        directions = [primary_direction]
        if -primary_direction != primary_direction:
            directions.append(-primary_direction)

        best_candidate_labels: list[CircularLabel] | None = None
        best_candidate_angles: list[float] | None = None
        best_candidate_score: tuple[int, int, float, int, int, int, int, float] | None = None

        for block_end in range(label_idx, len(ordered_labels)):
            block_width = block_end - label_idx
            for direction in directions:
                for shift_abs in delta_schedule:
                    shift_deg = direction * shift_abs
                    moved_start = current_angles[label_idx] + shift_deg
                    moved_end = current_angles[block_end] + shift_deg
                    if label_idx > 0 and moved_start <= (current_angles[label_idx - 1] + min_order_gap):
                        continue
                    if block_end < (len(ordered_labels) - 1) and moved_end >= (
                        current_angles[block_end + 1] - min_order_gap
                    ):
                        continue

                    candidate_labels, candidate_angles = _apply_block_shift(label_idx, block_end, shift_deg)
                    if not _half_progress_non_worsening(ordered_labels, candidate_labels, label_idx, block_end):
                        continue

                    candidate_focus_preferred_half = _preferred_half_for_label(candidate_labels[label_idx], total_length)
                    candidate_focus_half = _current_half_from_x(
                        float(candidate_labels[label_idx]["start_x"]),
                        axis_eps_px=HEMISPHERE_AXIS_EPS_PX,
                    )
                    candidate_half_bad = int(
                        candidate_focus_preferred_half != 0 and candidate_focus_half not in (0, candidate_focus_preferred_half)
                    )
                    if candidate_half_bad > current_half_bad:
                        continue
                    candidate_plain_overlap = _plain_overlap_cached(candidate_labels)
                    if candidate_plain_overlap > baseline_plain_overlap:
                        continue
                    if candidate_half_bad >= current_half_bad and baseline_line_intersections == 0:
                        candidate_mismatch_count_cheap, _ = _mismatch_metrics_cached(candidate_labels)
                        if candidate_mismatch_count_cheap >= baseline_mismatch_count:
                            continue

                    (
                        _candidate_plain_overlap,
                        candidate_min_gap_overlap,
                        candidate_mismatch_count,
                        candidate_mismatch_weighted,
                        candidate_leader_bbox_collisions,
                        candidate_line_intersections,
                    ) = _inner_rebalance_metrics_cached(candidate_labels)

                    if candidate_plain_overlap > baseline_plain_overlap:
                        continue
                    if candidate_leader_bbox_collisions > baseline_leader_bbox_collisions:
                        continue
                    if candidate_line_intersections > baseline_line_intersections:
                        continue

                    if not (
                        candidate_mismatch_count < baseline_mismatch_count
                        or candidate_line_intersections < baseline_line_intersections
                        or candidate_half_bad < current_half_bad
                    ):
                        continue

                    candidate_score = (
                        candidate_mismatch_count,
                        candidate_half_bad,
                        candidate_mismatch_weighted,
                        candidate_line_intersections,
                        candidate_min_gap_overlap,
                        candidate_leader_bbox_collisions,
                        block_width,
                        shift_abs,
                    )
                    if best_candidate_score is None or candidate_score < best_candidate_score:
                        best_candidate_score = candidate_score
                        best_candidate_labels = candidate_labels
                        best_candidate_angles = candidate_angles

                if (
                    best_candidate_score is not None
                    and best_candidate_score[0] == 0
                    and best_candidate_score[1] == 0
                    and best_candidate_score[3] == 0
                ):
                    break
            if (
                best_candidate_score is not None
                and best_candidate_score[0] == 0
                and best_candidate_score[1] == 0
                and best_candidate_score[3] == 0
            ):
                break

        if best_candidate_labels is None or best_candidate_angles is None:
            break
        ordered_labels = best_candidate_labels
        current_angles = best_candidate_angles

    (
        baseline_plain_overlap,
        baseline_min_gap_overlap,
        baseline_mismatch_count,
        baseline_mismatch_weighted,
        baseline_leader_bbox_collisions,
        baseline_line_intersections,
    ) = _inner_rebalance_metrics_cached(ordered_labels)
    if baseline_plain_overlap == 0 and baseline_mismatch_count <= 3 and baseline_line_intersections == 0:
        return ordered_labels

    fallback_best_labels: list[CircularLabel] | None = None
    fallback_best_score: tuple[int, float, int, int, int, float] | None = None
    for max_shift_cap in (22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 30.0):
        candidate_labels = _resolve_outer_label_overlaps_with_fixed_radii(
            [label.copy() for label in ordered_labels],
            total_length,
            max_angle_shift_deg=max_shift_cap,
        )
        if not _is_strict_no_overtake(candidate_labels):
            continue
        candidate_plain_overlap = _plain_overlap_cached(candidate_labels)
        if candidate_plain_overlap > baseline_plain_overlap:
            continue

        (
            _candidate_plain_overlap,
            candidate_min_gap_overlap,
            candidate_mismatch_count,
            candidate_mismatch_weighted,
            candidate_leader_bbox_collisions,
            candidate_line_intersections,
        ) = _inner_rebalance_metrics_cached(candidate_labels)

        if candidate_plain_overlap > baseline_plain_overlap:
            continue
        if candidate_leader_bbox_collisions > baseline_leader_bbox_collisions:
            continue
        if candidate_line_intersections > baseline_line_intersections:
            continue
        if not (
            candidate_mismatch_count < baseline_mismatch_count
            or candidate_line_intersections < baseline_line_intersections
            or candidate_min_gap_overlap < baseline_min_gap_overlap
        ):
            continue

        fallback_candidate_score = (
            candidate_mismatch_count,
            candidate_mismatch_weighted,
            candidate_line_intersections,
            candidate_min_gap_overlap,
            candidate_leader_bbox_collisions,
            max_shift_cap,
        )
        if fallback_best_score is None or fallback_candidate_score < fallback_best_score:
            fallback_best_score = fallback_candidate_score
            fallback_best_labels = candidate_labels

    if fallback_best_labels is not None:
        return fallback_best_labels

    return ordered_labels


def _refine_labels_to_preferred_hemisphere(
    labels: list[CircularLabel],
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
) -> list[CircularLabel]:
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

    def plain_overlap_count_for_candidate(label_idx: int, candidate_label: CircularLabel) -> int:
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
            is_inner_label = bool(label.get("is_inner", False))
            preferred_half = _preferred_half_for_label(label, total_length, axis_neutral_deg=axis_neutral_deg)
            if preferred_half == 0:
                continue

            current_half = _current_half_from_x(float(label["start_x"]), axis_eps_px=axis_eps_px)
            if current_half in (0, preferred_half):
                continue

            lower_bound, upper_bound = _cyclic_order_bounds(current_unwrapped_angles, label_idx, min_order_gap)
            if lower_bound >= upper_bound:
                continue

            target_angle_unwrapped = float(
                label.get("target_angle_unwrapped", angle_from_middle_unwrapped(float(label["middle"]), total_length))
            )
            target_angle_mod = target_angle_unwrapped % 360.0
            current_overlap_count = plain_overlap_count_for_candidate(label_idx, label)
            current_neighbor_overlap_count = _neighbor_min_gap_overlap_count(
                labels,
                label_idx,
                label,
                total_length,
            )
            best_candidate: tuple[CircularLabel, float, int, int] | None = None
            best_score: tuple[int, int, float] | None = None

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

                    candidate_neighbor_overlap_count = _neighbor_min_gap_overlap_count(
                        labels,
                        label_idx,
                        candidate_label,
                        total_length,
                    )
                    candidate_score = (
                        candidate_plain_total,
                        candidate_neighbor_overlap_count,
                        abs(offset),
                    )
                    if best_score is None or candidate_score < best_score:
                        best_score = candidate_score
                        best_candidate = (
                            candidate_label,
                            candidate_angle_unwrapped,
                            candidate_plain_total,
                            candidate_neighbor_overlap_count,
                        )

                if (
                    best_score is not None
                    and best_score[0] <= base_plain_overlaps
                    and best_score[1] <= current_neighbor_overlap_count
                ):
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

                    best_block_score: tuple[int, int, float, int] | None = None
                    best_block_result: tuple[list[CircularLabel], list[float], int, int] | None = None
                    current_block_neighbor_overlap = _neighbor_min_gap_overlap_count(
                        labels,
                        label_idx,
                        label,
                        total_length,
                    )
                    for block_start, block_end in block_candidates:
                        for step in range(1, block_max_steps + 1):
                            step_offset = step * angle_step_deg
                            for delta in (step_offset, -step_offset):
                                candidate_angles = current_unwrapped_angles.copy()
                                for move_idx in range(block_start, block_end + 1):
                                    candidate_angles[move_idx] += delta

                                if any(
                                    candidate_angles[check_idx]
                                    <= _cyclic_order_bounds(candidate_angles, check_idx, min_order_gap)[0]
                                    for check_idx in range(len(candidate_angles))
                                ):
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
                                candidate_neighbor_overlap_count = _neighbor_min_gap_overlap_count(
                                    candidate_labels,
                                    label_idx,
                                    candidate_labels[label_idx],
                                    total_length,
                                )
                                block_candidate_score = (
                                    candidate_plain_total,
                                    candidate_neighbor_overlap_count,
                                    abs(delta),
                                    block_width,
                                )
                                if best_block_score is None or block_candidate_score < best_block_score:
                                    best_block_score = block_candidate_score
                                    best_block_result = (
                                        candidate_labels,
                                        candidate_angles,
                                        candidate_plain_total,
                                        candidate_neighbor_overlap_count,
                                    )
                            if (
                                best_block_score is not None
                                and best_block_score[0] <= base_plain_overlaps
                                and best_block_score[1] <= current_block_neighbor_overlap
                            ):
                                break
                        if (
                            best_block_score is not None
                            and best_block_score[0] <= base_plain_overlaps
                            and best_block_score[1] <= current_block_neighbor_overlap
                        ):
                            break

                    if best_block_result is None:
                        continue

                    candidate_labels, candidate_angles, chosen_plain_total, chosen_neighbor_overlap = best_block_result
                    if is_inner_label and (
                        chosen_plain_total >= base_plain_overlaps
                        and chosen_neighbor_overlap >= current_block_neighbor_overlap
                    ):
                        continue
                    labels[:] = candidate_labels
                    current_unwrapped_angles = candidate_angles
                    base_plain_overlaps = chosen_plain_total
                    continue

                continue

            candidate_label, chosen_angle_unwrapped, chosen_plain_total, chosen_neighbor_overlap = best_candidate
            if is_inner_label and (
                chosen_plain_total >= base_plain_overlaps
                and chosen_neighbor_overlap >= current_neighbor_overlap_count
            ):
                continue
            label["start_x"] = candidate_label["start_x"]
            label["start_y"] = candidate_label["start_y"]
            if "angle_unwrapped" in label:
                label["angle_unwrapped"] = chosen_angle_unwrapped
            current_unwrapped_angles[label_idx] = chosen_angle_unwrapped
            base_plain_overlaps = chosen_plain_total

    return labels


def _polish_total_leader_distance(
    labels: list[CircularLabel],
    total_length: int,
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    *,
    min_order_gap: float = 0.05,
    angle_step_deg: float = 0.25,
    max_passes: int = 8,
) -> list[CircularLabel]:
    """
    Reduce total leader length while preserving overlap/order/hemisphere quality.

    This is a conservative final pass that only accepts strictly better global
    scores and never regresses plain/min-gap overlaps or hemisphere mismatches.
    """
    if not labels:
        return labels

    def _target_angle_unwrapped(label: CircularLabel) -> float:
        return float(
            label.get(
                "target_angle_unwrapped",
                angle_from_middle_unwrapped(float(label["middle"]), total_length),
            )
        )

    def _move_label(label_idx: int, angle_unwrapped: float) -> None:
        unwrapped[label_idx] = angle_unwrapped
        label = labels[label_idx]
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

    # Peer pruning (output-preserving): both overlap predicates are false for labels whose reach boxes are disjoint.
    reach_pad = _label_reach_pad(labels)
    reach_boxes = [_label_reach_box(label, reach_pad) for label in labels]

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
            current_angle = unwrapped[label_idx]
            target_angle = _target_angle_unwrapped(label)
            if math.isclose(current_angle, target_angle, abs_tol=1e-9):
                continue

            lower_bound, upper_bound = _cyclic_order_bounds(unwrapped, label_idx, min_order_gap)
            if lower_bound >= upper_bound:
                continue

            direction = 1.0 if target_angle > current_angle else -1.0
            max_steps = int(abs(target_angle - current_angle) / angle_step_deg)
            candidate_angles = [current_angle + direction * angle_step_deg * step for step in range(1, max_steps + 1)]
            if not candidate_angles or not math.isclose(candidate_angles[-1], target_angle, abs_tol=1e-9):
                candidate_angles.append(target_angle)

            current_plain_for_label = _count_local_label_overlaps(
                labels, total_length, label_idx, label, use_min_gap=False
            )
            current_min_gap_for_label = _count_local_label_overlaps(
                labels, total_length, label_idx, label, use_min_gap=True
            )
            current_leader = _leader_length_px(label, total_length)
            current_target_delta = abs(current_angle - target_angle)
            # hemisphere_mismatch_count always counts the current labels, so a move changes it by
            # this label's own flag only.
            others_mismatch_count = hemisphere_mismatch_count - (
                _label_hemisphere_mismatch_weight(label, total_length) is not None
            )

            candidate_positions = []
            for candidate_angle in candidate_angles:
                if candidate_angle <= lower_bound + 1e-9 or candidate_angle >= upper_bound - 1e-9:
                    continue
                candidate_positions.append(
                    (
                        candidate_angle,
                        *calculate_coordinates(
                            center_x,
                            center_y,
                            x_radius,
                            y_radius,
                            candidate_angle,
                            label["middle"],
                            total_length,
                        ),
                    )
                )
            if not candidate_positions:
                continue
            peers = touching_aabb_indices(
                reach_boxes,
                _label_reach_box(label, reach_pad, anchors=[(item[1], item[2]) for item in candidate_positions]),
                skip=label_idx,
            )

            best_strict: tuple[float, tuple[int, int, int, float, float]] | None = None
            best_relaxed: tuple[float, tuple[int, int, int, float, float], tuple[float, float, int]] | None = None
            for candidate_angle, candidate_x, candidate_y in candidate_positions:
                candidate_label = label.copy()
                candidate_label["angle_unwrapped"] = candidate_angle
                candidate_label["start_x"], candidate_label["start_y"] = candidate_x, candidate_y

                candidate_plain_for_label = _count_local_label_overlaps(
                    labels,
                    total_length,
                    label_idx,
                    candidate_label,
                    use_min_gap=False,
                    peers=peers,
                    cap=current_plain_for_label,
                )
                candidate_plain_total = plain_overlaps - current_plain_for_label + candidate_plain_for_label
                if candidate_plain_total > plain_overlaps:
                    continue

                candidate_min_gap_for_label = _count_local_label_overlaps(
                    labels,
                    total_length,
                    label_idx,
                    candidate_label,
                    use_min_gap=True,
                    peers=peers,
                    cap=max_allowed_min_gap_overlaps - min_gap_overlaps + current_min_gap_for_label,
                )
                candidate_min_gap_total = min_gap_overlaps - current_min_gap_for_label + candidate_min_gap_for_label
                if candidate_min_gap_total > max_allowed_min_gap_overlaps:
                    continue

                candidate_leader = _leader_length_px(candidate_label, total_length)
                candidate_leader_sum = leader_sum - current_leader + candidate_leader
                candidate_target_delta = abs(candidate_angle - target_angle)
                candidate_target_delta_sum = target_delta_sum - current_target_delta + candidate_target_delta

                candidate_mismatch_count = others_mismatch_count + (
                    _label_hemisphere_mismatch_weight(candidate_label, total_length) is not None
                )
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
            _move_label(label_idx, accepted_angle)
            reach_boxes[label_idx] = _label_reach_box(label, reach_pad)
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
    labels: list[CircularLabel],
    feature_radius,
    total_length,
    genome_len,
    cfg: GbdrawConfig,
    strands,
    is_outer,
    inner_labels_enabled: bool,
    arena_outer_radius: float | None = None,
    spacing_px: float | None = None,
    track_preset: str | None = None,
    *,
    dense_outer_rows: DenseOuterRows = "legacy",
    probe: bool = False,
):
    """Place one side's external labels on their arc.

    ``dense_outer_rows`` chooses the row gap of the dense legacy packer for
    outer labels: ``"spaced"`` keeps the circular label spacing (OV-479),
    ``"legacy"`` keeps the 0.1 px gap of earlier releases. When the dense outer
    candidate overlaps and the caller can use that answer (a ``"spaced"``
    attempt, which falls back to ``"legacy"``, or a fit ``probe``), placement
    stops with :class:`CircularLabelProbeStopped` instead of trying the slower
    candidates.
    """
    track_type = _normalize_circular_track_preset(track_preset or cfg.canvas.circular.track_type)

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
    label_spacing_px = _resolve_label_spacing_px(sorted_labels[0], base_margin_px=spacing_px)
    # OV-479: a spaced attempt keeps the circular label spacing between dense outer rows.
    legacy_row_gap_px = (
        max(label_spacing_px, LEGACY_MIN_ROW_GAP_PX)
        if is_outer and dense_outer_rows == "spaced"
        else LEGACY_MIN_ROW_GAP_PX
    )

    def _build_legacy_candidate() -> list[CircularLabel]:
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
            spacing_px=legacy_row_gap_px,
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
            y_margin=legacy_row_gap_px,
        )
        candidate_legacy = sort_labels(candidate_legacy)
        return _refine_labels_to_preferred_hemisphere(
            candidate_legacy,
            total_length,
            center_x,
            center_y,
            x_radius,
            y_radius,
        )

    candidate_legacy = None
    legacy_score = None
    legacy_with_leaders: list[CircularLabel] | None = None
    legacy_leader_collisions: int | None = None
    is_dense = len(sorted_labels) >= LEGACY_PLACEMENT_LABEL_THRESHOLD

    # Very dense sets are better handled by the legacy strategy: it keeps labels
    # closer to feature angles and avoids expensive combinatorial relaxation.
    if is_dense:
        candidate_legacy = _build_legacy_candidate()
        legacy_score = _placement_score(candidate_legacy, total_length)
        if is_outer and legacy_score[0] > 0 and (probe or dense_outer_rows == "spaced"):
            raise CircularLabelProbeStopped(legacy_score[0])
        if (not is_outer) and legacy_score[0] == 0:
            # Dense inner labels prioritize leader-collision-free readability and runtime over min-gap polishing.
            legacy_with_leaders = _assign_leader_start_points(
                [label.copy() for label in candidate_legacy],
                total_length,
            )
            legacy_leader_collisions = _count_label_leader_line_collisions(
                legacy_with_leaders,
                total_length,
                margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
            )
            if legacy_leader_collisions == 0:
                return legacy_with_leaders
        if is_outer and (not inner_labels_enabled) and legacy_score[0] == 0:
            return _assign_leader_start_points(candidate_legacy, total_length)
        if is_outer and legacy_score[0] == 0 and legacy_score[2] == 0:
            return _assign_leader_start_points(candidate_legacy, total_length)

    # Candidate A: current placement pipeline.
    candidate_current = [label.copy() for label in sorted_labels]
    candidate_current = place_labels_on_arc_fc(
        candidate_current,
        center_x,
        center_y,
        x_radius,
        y_radius,
        start_angle,
        end_angle,
        total_length,
        spacing_px=label_spacing_px,
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
        spacing_px=label_spacing_px,
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

    # Always compare dense candidates when both are available.
    if candidate_legacy is not None and legacy_score is not None:
        if is_dense and (not is_outer) and legacy_score[0] == 0 and best_score[0] == 0:
            if legacy_with_leaders is None:
                legacy_with_leaders = _assign_leader_start_points(
                    [label.copy() for label in candidate_legacy],
                    total_length,
                )
            if legacy_leader_collisions is None:
                legacy_leader_collisions = _count_label_leader_line_collisions(
                    legacy_with_leaders,
                    total_length,
                    margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
                )
            current_with_leaders = _assign_leader_start_points(
                [label.copy() for label in candidate_current],
                total_length,
            )
            current_leader_collisions = _count_label_leader_line_collisions(
                current_with_leaders,
                total_length,
                margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
            )
            legacy_dense_score = _dense_inner_candidate_score(legacy_score, legacy_leader_collisions)
            current_dense_score = _dense_inner_candidate_score(best_score, current_leader_collisions)
            if legacy_dense_score < current_dense_score:
                best_labels = candidate_legacy
                best_score = legacy_score
        elif legacy_score < best_score:
            best_labels = candidate_legacy
            best_score = legacy_score

    return _assign_leader_start_points(best_labels, total_length)


def _segment_span_bp(
    segment_start: int,
    segment_end: int,
    total_length: int,
) -> float:
    """Return segment span in bp, supporting origin-spanning segments."""
    if segment_start <= segment_end or total_length <= 0:
        return float(abs(int(segment_end - segment_start)))
    return float((float(total_length) - float(segment_start)) + float(segment_end))


def _is_label_embedded_in_segment(
    label_start: float,
    label_end: float,
    segment_start: int,
    segment_end: int,
    total_length: int,
) -> bool:
    """Return True when label interval is strictly contained in the target segment."""
    label_span = float(label_end) - float(label_start)
    if label_span <= 0.0:
        return False

    if segment_start <= segment_end:
        return (label_start > float(segment_start)) and (label_end < float(segment_end))

    if total_length <= 0:
        return False

    feature_span = (float(total_length) - float(segment_start)) + float(segment_end)
    if feature_span <= 0.0:
        return False

    start_offset = (float(label_start) - float(segment_start)) % float(total_length)
    end_offset = start_offset + label_span
    eps = 1e-9
    return (start_offset > eps) and (end_offset < (feature_span - eps))


def prepare_label_list(
    feature_dict,
    total_length,
    radius,
    track_ratio,
    profile: CircularRenderProfile,
    *,
    outer_arena: tuple[float, float] | None = None,
    feature_track_ratio_factor_override: float | None = None,
    feature_layout: CircularFeatureLayout | None = None,
    radial_layout: CircularRadialLayout | None = None,
    track_preset: str | None = None,
    feature_lane_direction: str | None = None,
    label_font_size: float | None = None,
    _candidate_cache: CircularLabelCandidateCache | None = None,
    record_transform: RecordDisplayTransform | None = None,
    probe: bool = False,
    definition_text: Sequence[Aabb] = (),
    label_filtering: Mapping[str, Any] | None = None,
) -> list[CircularLabel]:
    """Return the record's drawn label list: embedded labels, then external labels.

    With ``probe`` (a fit search), horizontal placement raises
    :class:`CircularLabelProbeStopped` as soon as the dense outer candidate
    overlaps, because such a layout does not fit. Otherwise slight defects of
    the final horizontal layout are repaired locally (LP-D57); ``definition_text``
    is the centre definition's lines, in the labels' coordinates.
    ``label_filtering``: the compiled maps when the caller holds them.
    """
    cfg = profile.config
    embedded_labels: list[CircularLabel] = []
    outer_labels: list[CircularLabel] = []
    inner_labels: list[CircularLabel] = []
    if label_filtering is None:
        label_filtering = label_filter_maps(cfg.labels.filtering.as_dict())
    label_rendering = profile.label_rendering
    length_threshold = cfg.labels.length_threshold.circular
    length_param = determine_length_parameter(total_length, length_threshold)
    track_type = _normalize_circular_track_preset(track_preset or cfg.canvas.circular.track_type)
    lane_direction = (
        str(feature_lane_direction).strip().lower()
        if feature_lane_direction is not None
        else circular_feature_lane_direction_for_preset(track_type)
    )
    if lane_direction not in {"inside", "outside", "split"}:
        lane_direction = circular_feature_lane_direction_for_preset(track_type)
    strandedness = profile.strandedness
    resolve_overlaps = profile.resolve_overlaps
    displaced_lanes = any(
        getattr(feature, "placement", None) is not None
        and feature.placement.level > 0 and feature.placement.requested_target is not None
        for feature in feature_dict.values()
    )
    lane_clearance = (resolve_overlaps and not strandedness) or displaced_lanes

    strands = "separate" if strandedness else "single"
    inner_labels_enabled = profile.inner_labels_enabled
    radius_factor = cfg.labels.radius_factor[track_type][strands][length_param]
    inner_radius_factor = cfg.labels.inner_radius_factor[track_type][strands][length_param]

    font_family = cfg.objects.text.font_family
    font_size: float = (
        float(label_font_size)
        if label_font_size is not None
        else circular_label_font_size(cfg, total_length)
    )
    label_spacing_px = float(cfg.labels.spacing.circular)

    track_ratio_factor = (
        float(feature_track_ratio_factor_override)
        if feature_track_ratio_factor_override is not None
        else float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    )
    cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
    feature_band_width_px = float(radius) * float(cds_ratio)
    effective_anchor_clearance_px = _effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=lane_clearance,
        strandedness=bool(strandedness),
        lane_direction=lane_direction,
        feature_band_width_px=feature_band_width_px,
    )
    effective_text_clearance_px = _effective_outer_text_clearance_px(
        resolve_overlaps=lane_clearance,
        strandedness=bool(strandedness),
        lane_direction=lane_direction,
        feature_band_width_px=feature_band_width_px,
    )
    automatic_circular_arrow_length_bp = calculate_circular_arrow_length(total_length)
    arrow_head_length_ratio = cfg.objects.features.arrow_geometry.head_length_ratio
    arrow_shaft_width_ratio = cfg.objects.features.arrow_geometry.shaft_width_ratio

    cached_candidates = (
        _candidate_cache.get("candidates")
        if _candidate_cache is not None
        else None
    )
    if isinstance(cached_candidates, tuple):
        candidates = cached_candidates
    else:
        candidates = build_circular_label_candidates(
            feature_dict,
            total_length,
            label_filtering,
            font_family=font_family,
            font_size=font_size,
            measure_text=_circular_label_box,
            record_transform=record_transform,
            overrides_only=profile.labels_overrides_only,
        )
        if _candidate_cache is not None:
            _candidate_cache["candidates"] = candidates
    candidates_by_id = {candidate.stable_id: candidate for candidate in candidates}
    feature_ids = {
        key: feature.feature_hash for key, feature in feature_dict.items()
    }

    for input_order, (stable_id, feature_object) in enumerate(feature_dict.items()):
        candidate = candidates_by_id.get(str(stable_id))
        if candidate is None:
            continue
        else:
            feature_label_text = candidate.text
            longest_segment_start = int(candidate.segment_start_bp)
            longest_segment_end = int(candidate.segment_end_bp)
            longeset_segment_middle = float(candidate.segment_middle_bp)
            is_embedded = False
            label_middle: float = 0
            coordinate_strand = candidate.strand

            # Get track_id for overlap resolution
            track_id = getattr(feature_object, 'feature_track_id', 0)
            if feature_layout is None:
                factors: list[float] = calculate_feature_position_factors_circular(
                    total_length, coordinate_strand, track_ratio, cds_ratio, offset, track_type, strandedness, track_id
                )
                feature_inner_radius = radius * min(float(factors[0]), float(factors[2]))
                feature_center_radius = radius * float(factors[1])
                feature_outer_radius = radius * max(float(factors[0]), float(factors[2]))
            else:
                lane = feature_layout.lane_for_track_id(int(track_id))
                feature_inner_radius = float(lane.inner_px)
                feature_center_radius = float(lane.center_px)
                feature_outer_radius = float(lane.outer_px)
            # The laid-out box, ``CIRCULAR_LABEL_GAP_EM`` included: the 1.1x embed test reads the same width as
            # the placer and the fit count, so the gap is part of deciding which labels go inside their feature.
            bbox_width_px = candidate.width_px
            bbox_height_px = candidate.height_px
            label_middle = longeset_segment_middle
            label_as_feature_length = total_length * (1.1 * bbox_width_px) / (2 * math.pi * radius)
            label_start = label_middle - (label_as_feature_length / 2)
            label_end = label_middle + (label_as_feature_length / 2)
            feature_middle_x: float = feature_center_radius * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
            feature_middle_y: float = feature_center_radius * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            feature_anchor_x: float = feature_middle_x
            feature_anchor_y: float = feature_middle_y
            assignment = getattr(feature_object, "placement", None)
            is_outer_label = (feature_object.strand == "positive") or (inner_labels_enabled is False)
            if (
                assignment is not None and assignment.level > 0
                and lane_direction == "split" and inner_labels_enabled
            ):
                is_outer_label = assignment.side == "outward"
            longest_segment_bp = _segment_span_bp(
                longest_segment_start,
                longest_segment_end,
                total_length,
            )
            circular_arrow_length_bp = resolve_circular_arrow_head_length_bp(
                arrow_head_length_ratio,
                abs(feature_outer_radius - feature_inner_radius),
                automatic_circular_arrow_length_bp,
                total_length,
                feature_center_radius,
                arrow_shaft_width_ratio,
            )
            is_short_directional_feature = (
                bool(getattr(feature_object, "is_directional", False))
                and longest_segment_bp < circular_arrow_length_bp
            )
            if is_outer_label:
                if lane_clearance and not is_short_directional_feature:
                    # Raised feature tracks look disconnected when leaders target the track center.
                    # Anchor to the outer edge so feature-to-label distance reads more naturally.
                    feature_anchor_x = feature_outer_radius * math.cos(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                    feature_anchor_y = feature_outer_radius * math.sin(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                # With an arena, its inner edge is the "anchor" radius for leader lines.
                middle_radius = (
                    float(outer_arena[0]) if outer_arena is not None else radius_factor * radius
                )
                if lane_clearance:
                    middle_radius = max(
                        float(middle_radius),
                        float(
                            feature_outer_radius
                            + effective_anchor_clearance_px
                            + OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX
                        ),
                    )
                middle_x = middle_radius * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y = middle_radius * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            else:
                if lane_clearance and not is_short_directional_feature:
                    # For displaced inner features, anchor leaders on the inner edge so lines
                    # do not cross over raised feature bodies.
                    feature_anchor_x = feature_inner_radius * math.cos(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                    feature_anchor_y = feature_inner_radius * math.sin(
                        math.radians(360.0 * (label_middle / total_length) - 90)
                    )
                middle_x = (inner_radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y = (inner_radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            is_embedded = _is_label_embedded_in_segment(
                label_start,
                label_end,
                longest_segment_start,
                longest_segment_end,
                total_length,
            )
            if label_rendering == "external_only":
                is_embedded = False
            elif label_rendering == "embedded_only" and not is_embedded:
                continue
            feature_id = feature_ids[stable_id]
            label_entry: CircularLabel = {
                "track_id": track_id,
                "label_text": feature_label_text,
                "feature_id": (
                    instance_svg_id(feature_id, feature_object.rendered_instance_index)
                    if feature_id and feature_object.rendered_instance_index is not None
                    else feature_id
                ),
                "stable_id": str(stable_id),
                "input_order": int(input_order),
                "preferred_middle": float(label_middle),
                "preferred_angle_deg": (360.0 * (float(label_middle) / float(total_length))) % 360.0,
                "middle": label_middle,
                "start": label_start,
                "end": label_end,
                "middle_x": middle_x,
                "middle_y": middle_y,
                "feature_middle_x": feature_middle_x,
                "feature_middle_y": feature_middle_y,
                "feature_anchor_x": feature_anchor_x,
                "feature_anchor_y": feature_anchor_y,
                "feature_inner_radius_px": feature_inner_radius,
                "feature_center_radius_px": feature_center_radius,
                "feature_outer_radius_px": feature_outer_radius,
                "track_preset": track_type,
                "width_px": bbox_width_px,
                "height_px": bbox_height_px,
                "strand": coordinate_strand,
                "is_embedded": is_embedded,
                "font_size": font_size,
                "font_family": font_family,
            }
            # Feature clearance of a displaced lane, for the external label stages.
            if lane_clearance and is_outer_label:
                label_entry["min_outer_start_radius_px"] = float(
                    feature_outer_radius
                    + effective_text_clearance_px
                    + OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX
                )
            elif lane_clearance:
                label_entry["max_inner_start_radius_px"] = max(
                    0.0,
                    feature_inner_radius
                    - effective_text_clearance_px
                    - OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX,
                )
                label_entry["max_inner_middle_radius_px"] = max(
                    0.0,
                    feature_inner_radius
                    - effective_anchor_clearance_px
                    - OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX,
                )

            if is_embedded:
                embedded_labels.append(label_entry)
            else:
                if is_outer_label:
                    label_entry["is_inner"] = False
                    outer_labels.append(label_entry)
                else:
                    label_entry["is_inner"] = True
                    inner_labels.append(label_entry)

    if profile.label_placement == "radial":
        outer_reserved_radius_px = 0.0
        if outer_arena is not None:
            outer_reserved_radius_px = max(float(outer_arena[0]), float(outer_arena[1]))
        if radial_layout is not None:
            outer_reserved_radius_px = max(
                outer_reserved_radius_px,
                float(radial_layout.outer_content_radius_px),
            )

        inner_reserved_outer_radius_px = 0.0
        if radial_layout is not None:
            definition_band = radial_layout.definition_reserved_band_px
            if definition_band is not None:
                inner_reserved_outer_radius_px = max(
                    inner_reserved_outer_radius_px,
                    float(definition_band.outer_px),
                )
            feature_inner_radius = (
                float(radial_layout.features.all_band_px.inner_px)
                if radial_layout.features is not None
                else float("inf")
            )
            for resolved_slot in radial_layout.slots:
                if resolved_slot.renderer in {"features", "feature_labels"}:
                    continue
                reserved_band = resolved_slot.reserved_band_px
                if reserved_band is None:
                    continue
                if float(reserved_band.outer_px) <= feature_inner_radius + 1e-6:
                    inner_reserved_outer_radius_px = max(
                        inner_reserved_outer_radius_px,
                        float(reserved_band.outer_px),
                    )

        radial_result = place_radial_labels(
            outer_labels,
            inner_labels if inner_labels_enabled else (),
            total_length=total_length,
            spacing_px=label_spacing_px,
            outer_reserved_radius_px=outer_reserved_radius_px,
            inner_reserved_outer_radius_px=inner_reserved_outer_radius_px,
        )
        external_labels = [label.copy() for label in radial_result.labels]
        for label in external_labels:
            label["required_radius_growth_px"] = radial_result.required_radius_growth_px
            label["layout_collision_counts"] = (
                radial_result.text_collision_count,
                radial_result.leader_text_collision_count,
                radial_result.leader_crossing_count,
                radial_result.order_violation_count,
            )
        return embedded_labels + external_labels

    def place_external_labels(dense_outer_rows: DenseOuterRows) -> list[CircularLabel]:
        outer_labels_rearranged = rearrange_labels_fc(
            [label.copy() for label in outer_labels],
            radius,
            total_length,
            length_param,
            cfg,
            strands,
            is_outer=True,
            inner_labels_enabled=inner_labels_enabled,
            arena_outer_radius=(float(outer_arena[1]) if outer_arena is not None else None),
            spacing_px=label_spacing_px,
            track_preset=track_type,
            dense_outer_rows=dense_outer_rows,
            probe=probe,
        )
        feature_outer_radius_intervals: list[tuple[float, float, float]] | None = None
        feature_inner_radius_intervals: list[tuple[float, float, float]] | None = None
        if lane_clearance:
            feature_outer_radius_intervals = _build_outer_feature_radius_intervals(
                feature_dict,
                total_length,
                radius,
                track_ratio,
                profile,
                feature_track_ratio_factor_override=feature_track_ratio_factor_override,
                feature_layout=feature_layout,
                track_preset=track_type,
            )
            feature_inner_radius_intervals = _build_inner_feature_radius_intervals(
                feature_dict,
                total_length,
                radius,
                track_ratio,
                profile,
                feature_track_ratio_factor_override=feature_track_ratio_factor_override,
                feature_layout=feature_layout,
                track_preset=track_type,
            )
            outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
                text_clearance_px=effective_text_clearance_px,
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
                    text_clearance_px=effective_text_clearance_px,
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
                text_clearance_px=effective_text_clearance_px,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
            # Final overlap pass, then enforce outer clearances as the final priority.
            outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
                outer_labels_rearranged,
                total_length,
            )
            outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
                text_clearance_px=effective_text_clearance_px,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
            outer_labels_rearranged = _update_outer_label_minimum_middle_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
                anchor_clearance_px=effective_anchor_clearance_px,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_middle_radius(
                outer_labels_rearranged,
                total_length,
            )
            # Reconcile angular overlaps once more after stricter radial clearances.
            outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
                outer_labels_rearranged,
                total_length,
            )
            outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
                text_clearance_px=effective_text_clearance_px,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
            outer_labels_rearranged = _update_outer_label_minimum_middle_radii_against_features(
                outer_labels_rearranged,
                total_length,
                feature_outer_radius_intervals,
                anchor_clearance_px=effective_anchor_clearance_px,
            )
            outer_labels_rearranged = _enforce_outer_label_minimum_middle_radius(
                outer_labels_rearranged,
                total_length,
            )
        inner_labels_rearranged = []
        if inner_labels_enabled:
            inner_labels_rearranged = rearrange_labels_fc(
                [label.copy() for label in inner_labels],
                radius,
                total_length,
                length_param,
                cfg,
                strands,
                is_outer=False,
                inner_labels_enabled=inner_labels_enabled,
                spacing_px=label_spacing_px,
                track_preset=track_type,
            )
            if (not resolve_overlaps) and len(inner_labels_rearranged) >= 2:
                inner_labels_rearranged = _rebalance_inner_labels_strict_order(
                    inner_labels_rearranged,
                    total_length,
                )
            if feature_inner_radius_intervals is not None and inner_labels_rearranged:
                inner_labels_rearranged = _stabilize_inner_label_clearance(
                    inner_labels_rearranged,
                    total_length,
                    feature_inner_radius_intervals,
                    text_clearance_px=effective_text_clearance_px,
                    anchor_clearance_px=effective_anchor_clearance_px,
                )

        external_labels = outer_labels_rearranged + inner_labels_rearranged
        if external_labels:
            # Collision optimization assumes angularly ordered neighbors.
            external_labels_sorted = sorted(external_labels, key=lambda label: float(label["middle"]))
            assign_leader_start_points(external_labels_sorted, total_length)
            _resolve_label_leader_line_collisions(
                external_labels_sorted,
                total_length,
                margin_px=LEADER_LABEL_COLLISION_MARGIN_PX,
                step_deg=LEADER_LABEL_SHIFT_STEP_DEG,
                max_shift_deg=LEADER_LABEL_MAX_SHIFT_DEG,
                max_passes=LEADER_LABEL_MAX_PASSES,
                min_order_gap_deg=LEADER_LABEL_MIN_ORDER_GAP_DEG,
            )
            reanchor_required = False
            if feature_outer_radius_intervals is not None and outer_labels_rearranged:
                outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                    outer_labels_rearranged,
                    total_length,
                    feature_outer_radius_intervals,
                    text_clearance_px=effective_text_clearance_px,
                )
                outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
                outer_labels_rearranged = _update_outer_label_minimum_middle_radii_against_features(
                    outer_labels_rearranged,
                    total_length,
                    feature_outer_radius_intervals,
                    anchor_clearance_px=effective_anchor_clearance_px,
                )
                outer_labels_rearranged = _enforce_outer_label_minimum_middle_radius(
                    outer_labels_rearranged,
                    total_length,
                )
                outer_labels_rearranged = _resolve_outer_label_overlaps_with_fixed_radii(
                    outer_labels_rearranged,
                    total_length,
                )
                outer_labels_rearranged = _update_outer_label_minimum_radii_against_features(
                    outer_labels_rearranged,
                    total_length,
                    feature_outer_radius_intervals,
                    text_clearance_px=effective_text_clearance_px,
                )
                outer_labels_rearranged = _enforce_outer_label_minimum_radius(outer_labels_rearranged, total_length)
                outer_labels_rearranged = _update_outer_label_minimum_middle_radii_against_features(
                    outer_labels_rearranged,
                    total_length,
                    feature_outer_radius_intervals,
                    anchor_clearance_px=effective_anchor_clearance_px,
                )
                outer_labels_rearranged = _enforce_outer_label_minimum_middle_radius(
                    outer_labels_rearranged,
                    total_length,
                )
                reanchor_required = True
            if feature_inner_radius_intervals is not None and inner_labels_rearranged:
                inner_labels_rearranged = _stabilize_inner_label_clearance(
                    inner_labels_rearranged,
                    total_length,
                    feature_inner_radius_intervals,
                    text_clearance_px=effective_text_clearance_px,
                    anchor_clearance_px=effective_anchor_clearance_px,
                )
                reanchor_required = True
            if reanchor_required:
                assign_leader_start_points(external_labels_sorted, total_length)

        external_labels_ordered = sorted(
            outer_labels_rearranged + inner_labels_rearranged,
            key=lambda label: float(label["middle"]),
        )
        # The drawer reads leader starts from this final assignment.
        assign_leader_start_points(external_labels_ordered, total_length)
        return external_labels_ordered

    # OV-479 / LP-D45: dense outer rows keep the circular label spacing when the
    # drawn layout then fits; otherwise the record is placed as in earlier releases.
    if not (
        len(outer_labels) >= LEGACY_PLACEMENT_LABEL_THRESHOLD
        and max(label_spacing_px, LEGACY_MIN_ROW_GAP_PX) > LEGACY_MIN_ROW_GAP_PX
    ):
        final_labels = embedded_labels + place_external_labels("legacy")
    else:
        # A probe asks only whether either layout fits. It tries today's layout first
        # and stops when its dense rows overlap: wider spaced rows overlap as well.
        attempts: tuple[DenseOuterRows, ...] = ("legacy", "spaced") if probe else ("spaced", "legacy")
        final_labels = []
        for dense_outer_rows in attempts:
            try:
                attempt_labels = place_external_labels(dense_outer_rows)
            except CircularLabelProbeStopped:
                if dense_outer_rows == "legacy":
                    raise
                continue
            if circular_label_fits(attempt_labels, total_length).get("circular_outer", FITS).fits:
                final_labels = embedded_labels + attempt_labels
                break
            if dense_outer_rows == "legacy":
                final_labels = embedded_labels + attempt_labels
    if probe:
        return final_labels
    # LP-D57: repair slight defects of the final layout (Q7: not in probes; Q2: not the spaced attempt).
    from .circular_repair import repair_circular_labels

    repaired = repair_circular_labels(final_labels, total_length, definition_text=definition_text)
    return final_labels if repaired is None else repaired


__all__ = [
    "minimum_bbox_gap_px",
    "CircularLabelProbeStopped",
    "calculate_angle_degrees",
    "circular_label_fits",
    "circular_label_shape",
    "circular_label_shapes",
    "circular_label_font_size",
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
