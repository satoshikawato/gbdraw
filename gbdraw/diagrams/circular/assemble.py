#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import logging
import math
import copy
from typing import Any, Optional, Callable, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...analysis.conservation import (  # type: ignore[reportMissingImports]
    ConservationTrack,
    conservation_track_gradient_colors,
)
from ...analysis.depth_tracks import (  # type: ignore[reportMissingImports]
    DepthTrackData,
    sync_depth_track_legend_entries,
)
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...core.text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...labels.circular import (  # type: ignore[reportMissingImports]
    assign_leader_start_points,
    minimum_bbox_gap_px,
    prepare_label_list,
    x_overlap,
    y_overlap,
)
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...layout.circular import calculate_feature_position_factors_circular  # type: ignore[reportMissingImports]
from ...layout.circular_depth_axis import depth_axis_tick_font_size_px  # type: ignore[reportMissingImports]
from ...layout.common import calculate_cds_ratio  # type: ignore[reportMissingImports]
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]
from ...tracks import (  # type: ignore[reportMissingImports]
    CircularTrackSlot,
)
from ...tracks.circular import tick_sides_for_tick_label_layout  # type: ignore[reportMissingImports]

from .builders import (
    add_axis_group_on_canvas,
    add_conservation_group_on_canvas,
    add_depth_group_on_canvas,
    add_gc_content_group_on_canvas,
    add_gc_skew_group_on_canvas,
    add_labels_group_on_canvas,
    add_legend_group_on_canvas,
    add_record_definition_group_on_canvas,
    add_record_group_on_canvas,
    add_tick_group_on_canvas,
)
from ...render.groups.circular.definition import DefinitionGroup  # type: ignore[reportMissingImports]
from .radial_layout import (  # type: ignore[reportMissingImports]
    CircularRadialLayout,
    CircularResolvedSlot,
    CircularTickLayout,
    build_circular_feature_layout,
    resolve_circular_radial_layout,
)
from .presets import (  # type: ignore[reportMissingImports]
    CircularPresetContext,
    circular_feature_lane_direction_for_preset,
    circular_radial_plan_for_preset,
    circular_track_slots_from_preset_order,
    normalize_circular_track_preset,
)


LEGEND_LABEL_MARGIN_PX = 4.0
LABEL_NUDGE_STEP_PX = 6.0
MAX_LABEL_NUDGE_PX = 180.0
LEGEND_SHIFT_STEP_PX = 16.0
MAX_LEGEND_SHIFT_STEPS = 60
CANVAS_EXPAND_STEP_PX = 32.0
MAX_CANVAS_EXPAND_STEPS = 24
MIN_LABEL_ORDER_GAP_RAD = 1e-4
LABEL_CANVAS_PADDING_PX = 8.0
FEATURE_BAND_EPSILON = 1e-6
SINGLE_LEGEND_EDGE_MIN_PX = 16.0
SINGLE_LEGEND_CONTENT_GAP_MIN_PX = 12.0


logger = logging.getLogger(__name__)


def _circular_preset_for_layout(
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> str:
    raw = getattr(canvas_config, "circular_track_preset", None)
    if raw is None:
        raw = cfg.canvas.circular.track_type
    return normalize_circular_track_preset(str(raw))


def _lane_direction_for_feature_slot(
    feature_slot: CircularTrackSlot | None,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> str:
    if feature_slot is not None:
        raw = feature_slot.params.get("lane_direction", feature_slot.params.get("lanes"))
        if raw is not None:
            direction = str(raw).strip().lower()
            if direction in {"inside", "outside", "split"}:
                return direction
    return circular_feature_lane_direction_for_preset(
        _circular_preset_for_layout(canvas_config, cfg)
    )


def _sync_canvas_viewbox(canvas: Drawing, canvas_config: CircularCanvasConfigurator) -> None:
    """Sync drawing viewport attrs with mutable canvas config values."""
    canvas.attribs["width"] = f"{canvas_config.total_width}px"
    canvas.attribs["height"] = f"{canvas_config.total_height}px"
    canvas.attribs["viewBox"] = f"0 0 {canvas_config.total_width} {canvas_config.total_height}"


def _resolved_slots_from_radial_layout(
    radial_layout: CircularRadialLayout,
    layout_slots: Sequence[CircularTrackSlot],
) -> list[CircularResolvedSlot]:
    del layout_slots
    return list(radial_layout.slots)


def _legend_bbox(canvas_config: CircularCanvasConfigurator, legend_config: LegendDrawingConfigurator) -> tuple[float, float, float, float]:
    """Return legend bbox on canvas as (min_x, min_y, max_x, max_y)."""
    min_x = float(canvas_config.legend_offset_x)
    min_y = float(canvas_config.legend_offset_y) - 0.5 * float(legend_config.color_rect_size)
    max_x = min_x + float(legend_config.legend_width)
    max_y = min_y + float(legend_config.legend_height)
    return min_x, min_y, max_x, max_y


def _bbox_overlaps(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    """Return True when two axis-aligned bboxes overlap."""
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def _unwrap_angle_near_reference(angle_rad: float, reference_rad: float) -> float:
    """Project angle to the nearest 2*pi branch around reference."""
    return angle_rad + (2.0 * math.pi) * round((reference_rad - angle_rad) / (2.0 * math.pi))


def _label_target_unwrapped_angle(label: dict[str, Any], total_length: int) -> float:
    """Target label angle derived from genomic midpoint, in unwrapped radians."""
    return (2.0 * math.pi * (float(label["middle"]) / float(total_length))) - (0.5 * math.pi)


def _label_unwrapped_angle(label: dict[str, Any], total_length: int) -> float:
    """Current label angle projected near its target branch."""
    current_angle = math.atan2(float(label["start_y"]), float(label["start_x"]))
    target_angle = _label_target_unwrapped_angle(label, total_length)
    return _unwrap_angle_near_reference(current_angle, target_angle)


def _label_bbox_local(label: dict[str, Any], total_length: int, margin_px: float) -> tuple[float, float, float, float]:
    """Approximate label bbox in local (centered) coordinates."""
    half_margin = 0.5 * margin_px
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    width_px = float(label["width_px"])
    height_px = float(label["height_px"])
    is_inner = bool(label.get("is_inner", False))

    if not is_inner:
        if start_x > 0:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin
        else:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
    else:
        if start_x > 0:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
        else:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin

    angle_deg = (360.0 * (float(label["middle"]) / float(total_length))) % 360.0
    if not is_inner:
        if 0.0 <= angle_deg < 10.0 or angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        elif 170.0 <= angle_deg < 190.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    else:
        if 0.0 <= angle_deg < 10.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y
        elif angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    return min_x, min_y, max_x, max_y


def _label_bbox_on_canvas(
    label: dict[str, Any],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    margin_px: float,
) -> tuple[float, float, float, float]:
    """Return label bbox on canvas as (min_x, min_y, max_x, max_y)."""
    min_x, min_y, max_x, max_y = _label_bbox_local(label, total_length, margin_px)
    return (
        min_x + float(canvas_config.offset_x),
        min_y + float(canvas_config.offset_y),
        max_x + float(canvas_config.offset_x),
        max_y + float(canvas_config.offset_y),
    )


def _external_label_bounds_on_canvas(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    *,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> tuple[float, float, float, float] | None:
    """Return bounds of all non-embedded labels in canvas coordinates."""
    min_x: float | None = None
    min_y: float | None = None
    max_x: float | None = None
    max_y: float | None = None

    for label in labels:
        if label.get("is_embedded"):
            continue
        box_min_x, box_min_y, box_max_x, box_max_y = _label_bbox_on_canvas(
            label,
            total_length,
            canvas_config,
            margin_px=margin_px,
        )
        min_x = box_min_x if min_x is None else min(min_x, box_min_x)
        min_y = box_min_y if min_y is None else min(min_y, box_min_y)
        max_x = box_max_x if max_x is None else max(max_x, box_max_x)
        max_y = box_max_y if max_y is None else max(max_y, box_max_y)

    if min_x is None or min_y is None or max_x is None or max_y is None:
        return None
    return min_x, min_y, max_x, max_y


def _expand_canvas_to_fit_external_labels(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    *,
    padding_px: float = LABEL_CANVAS_PADDING_PX,
) -> bool:
    """Expand canvas so all non-embedded labels are inside with padding."""
    bounds = _external_label_bounds_on_canvas(labels, total_length, canvas_config)
    if bounds is None:
        return False

    total_width = float(canvas_config.total_width)
    total_height = float(canvas_config.total_height)

    grow_left = max(0.0, float(padding_px) - float(bounds[0]))
    grow_top = max(0.0, float(padding_px) - float(bounds[1]))
    grow_right = max(0.0, float(bounds[2]) - (total_width - float(padding_px)))
    grow_bottom = max(0.0, float(bounds[3]) - (total_height - float(padding_px)))

    if grow_left <= 1e-6 and grow_top <= 1e-6 and grow_right <= 1e-6 and grow_bottom <= 1e-6:
        return False

    canvas_config.total_width = total_width + grow_left + grow_right
    canvas_config.total_height = total_height + grow_top + grow_bottom
    canvas_config.offset_x = float(canvas_config.offset_x) + grow_left
    canvas_config.offset_y = float(canvas_config.offset_y) + grow_top

    if hasattr(canvas_config, "legend_offset_x"):
        canvas_config.legend_offset_x = float(canvas_config.legend_offset_x) + grow_left
    if hasattr(canvas_config, "legend_offset_y"):
        canvas_config.legend_offset_y = float(canvas_config.legend_offset_y) + grow_top

    return True


def _expand_canvas_to_fit_radius(
    canvas_config: CircularCanvasConfigurator,
    required_radius_px: float,
    *,
    padding_px: float = LABEL_CANVAS_PADDING_PX,
) -> bool:
    """Expand canvas so a centered circle with required radius fully fits."""
    required_extent = float(required_radius_px) + float(padding_px)
    if required_extent <= 0:
        return False

    total_width = float(canvas_config.total_width)
    total_height = float(canvas_config.total_height)
    offset_x = float(canvas_config.offset_x)
    offset_y = float(canvas_config.offset_y)

    grow_left = max(0.0, required_extent - offset_x)
    grow_right = max(0.0, required_extent - (total_width - offset_x))
    grow_top = max(0.0, required_extent - offset_y)
    grow_bottom = max(0.0, required_extent - (total_height - offset_y))

    if grow_left <= 1e-6 and grow_top <= 1e-6 and grow_right <= 1e-6 and grow_bottom <= 1e-6:
        return False

    canvas_config.total_width = total_width + grow_left + grow_right
    canvas_config.total_height = total_height + grow_top + grow_bottom
    canvas_config.offset_x = offset_x + grow_left
    canvas_config.offset_y = offset_y + grow_top

    if hasattr(canvas_config, "legend_offset_x"):
        canvas_config.legend_offset_x = float(canvas_config.legend_offset_x) + grow_left
    if hasattr(canvas_config, "legend_offset_y"):
        canvas_config.legend_offset_y = float(canvas_config.legend_offset_y) + grow_top

    return True


def _legend_collision_indices(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> list[int]:
    """Return indices of labels that currently collide with the legend."""
    legend_box = _legend_bbox(canvas_config, legend_config)
    collided: list[int] = []
    for idx, label in enumerate(labels):
        if label.get("is_embedded"):
            continue
        label_box = _label_bbox_on_canvas(label, total_length, canvas_config, margin_px=margin_px)
        if _bbox_overlaps(label_box, legend_box):
            collided.append(idx)
    return collided


def _labels_collide_with_legend(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> bool:
    """Whether any external label overlaps with the legend bbox."""
    return bool(_legend_collision_indices(labels, total_length, canvas_config, legend_config))


def _label_overlaps_other_labels(candidate: dict[str, Any], labels: list[dict[str, Any]], idx: int, total_length: int) -> bool:
    """Check candidate against peer labels with the existing overlap predicates."""
    for peer_idx, peer in enumerate(labels):
        if peer_idx == idx or peer.get("is_embedded"):
            continue
        min_gap_px = minimum_bbox_gap_px(candidate, peer)
        if y_overlap(candidate, peer, total_length, min_gap_px) and x_overlap(candidate, peer, minimum_margin=min_gap_px):
            return True
    return False


def _legend_center_local(
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> tuple[float, float]:
    """Return legend center in local (circle-centered) coordinates."""
    legend_min_x, legend_min_y, legend_max_x, legend_max_y = _legend_bbox(canvas_config, legend_config)
    legend_center_x = 0.5 * (legend_min_x + legend_max_x)
    legend_center_y = 0.5 * (legend_min_y + legend_max_y)
    return legend_center_x - float(canvas_config.offset_x), legend_center_y - float(canvas_config.offset_y)


def _preferred_angular_shift_sign(
    label: dict[str, Any],
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> int:
    """Return preferred angular direction (+1 ccw / -1 cw) to move away from legend."""
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    radius = math.hypot(start_x, start_y)
    if radius <= 1e-6:
        return 1
    legend_local_x, legend_local_y = _legend_center_local(canvas_config, legend_config)
    away_x = start_x - legend_local_x
    away_y = start_y - legend_local_y
    # Unit tangents at current point (CCW and CW).
    ccw_tx = -start_y / radius
    ccw_ty = start_x / radius
    cw_tx = -ccw_tx
    cw_ty = -ccw_ty
    ccw_score = ccw_tx * away_x + ccw_ty * away_y
    cw_score = cw_tx * away_x + cw_ty * away_y
    return 1 if ccw_score >= cw_score else -1


def _expand_shift_block_to_preserve_order(
    labels: list[dict[str, Any]],
    total_length: int,
    side_indices: list[int],
    center_pos: int,
    delta_angle: float,
) -> tuple[int, int]:
    """Expand move block so angular order cannot be inverted by this shift."""
    left = center_pos
    right = center_pos
    while True:
        expanded = False
        moved_left_idx = side_indices[left]
        moved_right_idx = side_indices[right]
        moved_left_angle = _label_unwrapped_angle(labels[moved_left_idx], total_length) + delta_angle
        moved_right_angle = _label_unwrapped_angle(labels[moved_right_idx], total_length) + delta_angle

        if left > 0:
            prev_idx = side_indices[left - 1]
            prev_angle = _label_unwrapped_angle(labels[prev_idx], total_length)
            if moved_left_angle <= prev_angle + MIN_LABEL_ORDER_GAP_RAD:
                left -= 1
                expanded = True

        if right < len(side_indices) - 1:
            next_idx = side_indices[right + 1]
            next_angle = _label_unwrapped_angle(labels[next_idx], total_length)
            if moved_right_angle >= next_angle - MIN_LABEL_ORDER_GAP_RAD:
                right += 1
                expanded = True

        if not expanded:
            break
    return left, right


def _try_shift_labels_away_from_legend(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> bool:
    """Resolve collisions by moving labels first (preferred strategy)."""
    if not labels:
        return True

    max_passes = 4
    for _ in range(max_passes):
        collided_indices = _legend_collision_indices(labels, total_length, canvas_config, legend_config)
        if not collided_indices:
            return True

        changed = False
        legend_box = _legend_bbox(canvas_config, legend_config)
        for idx in collided_indices:
            label = labels[idx]
            if label.get("is_embedded"):
                continue
            start_x = float(label["start_x"])
            start_y = float(label["start_y"])
            radius = math.hypot(start_x, start_y)
            if radius <= 1e-6:
                continue
            moving_side = 1 if start_x >= 0 else -1
            preferred_sign = _preferred_angular_shift_sign(label, canvas_config, legend_config)
            side_indices = [
                side_idx
                for side_idx, side_label in enumerate(labels)
                if not side_label.get("is_embedded") and ((float(side_label["start_x"]) >= 0) == (moving_side > 0))
            ]
            side_indices.sort(key=lambda side_idx: float(labels[side_idx]["middle"]))
            if idx not in side_indices:
                continue
            center_pos = side_indices.index(idx)

            placed = False
            shift_px = LABEL_NUDGE_STEP_PX
            while shift_px <= MAX_LABEL_NUDGE_PX and not placed:
                angle_delta = shift_px / radius
                for direction_sign in (preferred_sign, -preferred_sign):
                    signed_delta = direction_sign * angle_delta
                    block_left, block_right = _expand_shift_block_to_preserve_order(
                        labels, total_length, side_indices, center_pos, signed_delta
                    )
                    block_indices = side_indices[block_left : block_right + 1]

                    candidate_positions: dict[int, tuple[float, float]] = {}
                    candidate_labels: dict[int, dict[str, Any]] = {}
                    candidate_valid = True

                    for block_idx in block_indices:
                        block_label = labels[block_idx]
                        block_radius = math.hypot(float(block_label["start_x"]), float(block_label["start_y"]))
                        if block_radius <= 1e-6:
                            candidate_valid = False
                            break
                        block_angle = math.atan2(float(block_label["start_y"]), float(block_label["start_x"]))
                        candidate_angle = block_angle + signed_delta
                        candidate_x = block_radius * math.cos(candidate_angle)
                        candidate_y = block_radius * math.sin(candidate_angle)

                        block_sign = 1 if float(block_label["start_x"]) >= 0 else -1
                        if block_sign > 0 and candidate_x <= 0:
                            candidate_valid = False
                            break
                        if block_sign < 0 and candidate_x >= 0:
                            candidate_valid = False
                            break

                        candidate = block_label.copy()
                        candidate["start_x"] = float(candidate_x)
                        candidate["start_y"] = float(candidate_y)
                        candidate_positions[block_idx] = (float(candidate_x), float(candidate_y))
                        candidate_labels[block_idx] = candidate

                    if not candidate_valid:
                        continue

                    for candidate in candidate_labels.values():
                        candidate_box = _label_bbox_on_canvas(
                            candidate, total_length, canvas_config, margin_px=LEGEND_LABEL_MARGIN_PX
                        )
                        if _bbox_overlaps(candidate_box, legend_box):
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, candidate in candidate_labels.items():
                        for peer_idx, peer_label in enumerate(labels):
                            if peer_idx == block_idx or peer_label.get("is_embedded"):
                                continue
                            peer = candidate_labels.get(peer_idx, peer_label)
                            min_gap_px = minimum_bbox_gap_px(candidate, peer)
                            if y_overlap(candidate, peer, total_length, min_gap_px) and x_overlap(
                                candidate, peer, minimum_margin=min_gap_px
                            ):
                                candidate_valid = False
                                break
                        if not candidate_valid:
                            break
                    if not candidate_valid:
                        continue

                    def _unwrapped_for_order(check_idx: int) -> float:
                        candidate = candidate_labels.get(check_idx)
                        if candidate is not None:
                            return _label_unwrapped_angle(candidate, total_length)
                        return _label_unwrapped_angle(labels[check_idx], total_length)

                    for order_pos in range(1, len(side_indices)):
                        prev_idx = side_indices[order_pos - 1]
                        curr_idx = side_indices[order_pos]
                        if _unwrapped_for_order(curr_idx) <= _unwrapped_for_order(prev_idx) + MIN_LABEL_ORDER_GAP_RAD:
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, (candidate_x, candidate_y) in candidate_positions.items():
                        labels[block_idx]["start_x"] = candidate_x
                        labels[block_idx]["start_y"] = candidate_y
                    changed = True
                    placed = True
                    break
                shift_px += LABEL_NUDGE_STEP_PX

        if not changed:
            break

    return not _labels_collide_with_legend(labels, total_length, canvas_config, legend_config)


def _legend_offset_bounds(
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> tuple[float, float, float, float]:
    """Return min/max bounds for legend offsets (x_min, x_max, y_min, y_max)."""
    x_min = 0.0
    x_max = max(0.0, float(canvas_config.total_width) - float(legend_config.legend_width))
    y_min = 0.5 * float(legend_config.color_rect_size)
    y_max = max(y_min, float(canvas_config.total_height) - float(legend_config.legend_height) + 0.5 * float(legend_config.color_rect_size))
    return x_min, x_max, y_min, y_max


def _clamp_legend_offsets(canvas_config: CircularCanvasConfigurator, legend_config: LegendDrawingConfigurator) -> None:
    """Keep legend offsets inside the current canvas."""
    x_min, x_max, y_min, y_max = _legend_offset_bounds(canvas_config, legend_config)
    canvas_config.legend_offset_x = min(max(float(canvas_config.legend_offset_x), x_min), x_max)
    canvas_config.legend_offset_y = min(max(float(canvas_config.legend_offset_y), y_min), y_max)


def _legend_push_direction(canvas_config: CircularCanvasConfigurator, legend_config: LegendDrawingConfigurator) -> tuple[float, float]:
    """Direction that moves legend farther away from the circular map center."""
    legend_min_x, legend_min_y, legend_max_x, legend_max_y = _legend_bbox(canvas_config, legend_config)
    legend_center_x = 0.5 * (legend_min_x + legend_max_x)
    legend_center_y = 0.5 * (legend_min_y + legend_max_y)
    dx = legend_center_x - float(canvas_config.offset_x)
    dy = legend_center_y - float(canvas_config.offset_y)
    norm = math.hypot(dx, dy)
    if norm > 1e-6:
        return dx / norm, dy / norm
    position = str(getattr(canvas_config, "legend_position", "right"))
    if "left" in position:
        return -1.0, 0.0
    if "upper" in position:
        return 0.0, -1.0
    if "lower" in position:
        return 0.0, 1.0
    return 1.0, 0.0


def _try_move_legend_away_from_labels(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> bool:
    """Resolve collisions by moving legend within the current canvas bounds."""
    if not _labels_collide_with_legend(labels, total_length, canvas_config, legend_config):
        return True
    dir_x, dir_y = _legend_push_direction(canvas_config, legend_config)
    for _ in range(MAX_LEGEND_SHIFT_STEPS):
        prev_x = float(canvas_config.legend_offset_x)
        prev_y = float(canvas_config.legend_offset_y)
        canvas_config.legend_offset_x = prev_x + dir_x * LEGEND_SHIFT_STEP_PX
        canvas_config.legend_offset_y = prev_y + dir_y * LEGEND_SHIFT_STEP_PX
        _clamp_legend_offsets(canvas_config, legend_config)
        if (
            abs(float(canvas_config.legend_offset_x) - prev_x) < 1e-6
            and abs(float(canvas_config.legend_offset_y) - prev_y) < 1e-6
        ):
            break
        if not _labels_collide_with_legend(labels, total_length, canvas_config, legend_config):
            return True
    return False


def _expand_canvas_for_legend(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> bool:
    """Expand canvas in legend direction until legend-label collisions disappear."""
    if not _labels_collide_with_legend(labels, total_length, canvas_config, legend_config):
        return True

    dir_x, dir_y = _legend_push_direction(canvas_config, legend_config)
    sign_x = 1 if dir_x > 0.25 else (-1 if dir_x < -0.25 else 0)
    sign_y = 1 if dir_y > 0.25 else (-1 if dir_y < -0.25 else 0)
    if sign_x == 0 and sign_y == 0:
        sign_x = 1

    for _ in range(MAX_CANVAS_EXPAND_STEPS):
        if sign_x > 0:
            canvas_config.total_width = float(canvas_config.total_width) + CANVAS_EXPAND_STEP_PX
            canvas_config.legend_offset_x = float(canvas_config.legend_offset_x) + CANVAS_EXPAND_STEP_PX
        elif sign_x < 0:
            canvas_config.total_width = float(canvas_config.total_width) + CANVAS_EXPAND_STEP_PX
            canvas_config.offset_x = float(canvas_config.offset_x) + CANVAS_EXPAND_STEP_PX

        if sign_y > 0:
            canvas_config.total_height = float(canvas_config.total_height) + CANVAS_EXPAND_STEP_PX
            canvas_config.legend_offset_y = float(canvas_config.legend_offset_y) + CANVAS_EXPAND_STEP_PX
        elif sign_y < 0:
            canvas_config.total_height = float(canvas_config.total_height) + CANVAS_EXPAND_STEP_PX
            canvas_config.offset_y = float(canvas_config.offset_y) + CANVAS_EXPAND_STEP_PX

        _clamp_legend_offsets(canvas_config, legend_config)
        if not _labels_collide_with_legend(labels, total_length, canvas_config, legend_config):
            return True
    return False


def _translate_canvas_top_level_groups(canvas: Drawing, *, dy: float) -> None:
    """Translate top-level canvas elements in Y (except defs)."""
    if abs(float(dy)) <= 1e-6:
        return
    for element in getattr(canvas, "elements", []):
        class_name = element.__class__.__name__.lower()
        element_name = str(getattr(element, "elementname", "")).lower()
        if class_name == "defs" or element_name == "defs":
            continue
        translate = getattr(element, "translate", None)
        if callable(translate):
            translate(0, float(dy))


def _default_feature_outer_radius_px(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    feature_track_ratio_factor_override: float | None,
) -> float:
    """Return a conservative outer feature radius when no precomputed feature band is available."""
    length_param = str(canvas_config.length_param)
    track_ratio_factor = (
        float(feature_track_ratio_factor_override)
        if feature_track_ratio_factor_override is not None
        else float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    )
    cds_ratio, offset = calculate_cds_ratio(
        float(canvas_config.track_ratio), length_param, track_ratio_factor
    )
    factors_positive = calculate_feature_position_factors_circular(
        total_length,
        "positive",
        float(canvas_config.track_ratio),
        cds_ratio,
        offset,
        _circular_preset_for_layout(canvas_config, cfg),
        bool(cfg.canvas.strandedness),
        0,
    )
    factors_negative = calculate_feature_position_factors_circular(
        total_length,
        "negative",
        float(canvas_config.track_ratio),
        cds_ratio,
        offset,
        _circular_preset_for_layout(canvas_config, cfg),
        bool(cfg.canvas.strandedness),
        0,
    )
    max_factor = max(
        [abs(float(factor)) for factor in [*factors_positive, *factors_negative]],
        default=1.0,
    )
    return float(canvas_config.radius) * float(max_factor)


def _resolve_content_vertical_bounds_for_single_legend(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    show_features: bool,
    precomputed_feature_dict: dict | None,
    rendered_feature_band_all_tracks: tuple[float, float] | None,
    feature_track_ratio_factor_override: float | None,
    tick_label_annulus: tuple[float, float] | None,
    external_label_bounds: tuple[float, float, float, float] | None,
) -> tuple[float, float]:
    """Return content top/bottom bounds in canvas coordinates for top/bottom legend placement."""
    content_outer_radius = float(canvas_config.radius)

    feature_band = rendered_feature_band_all_tracks
    if show_features and feature_band is None and precomputed_feature_dict is not None:
        length_param = str(canvas_config.length_param)
        track_ratio_factor = (
            float(feature_track_ratio_factor_override)
            if feature_track_ratio_factor_override is not None
            else float(cfg.canvas.circular.track_ratio_factors[length_param][0])
        )
        feature_band = _compute_feature_band_bounds_px(
            precomputed_feature_dict,
            total_length,
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=length_param,
            track_ratio_factor=track_ratio_factor,
            cfg=cfg,
        )
    if show_features and feature_band is None:
        default_feature_outer = _default_feature_outer_radius_px(
            total_length=total_length,
            canvas_config=canvas_config,
            cfg=cfg,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        )
        content_outer_radius = max(content_outer_radius, default_feature_outer)
    elif feature_band is not None:
        content_outer_radius = max(
            content_outer_radius,
            abs(float(feature_band[0])),
            abs(float(feature_band[1])),
        )

    if tick_label_annulus is not None:
        content_outer_radius = max(
            content_outer_radius,
            abs(float(tick_label_annulus[0])),
            abs(float(tick_label_annulus[1])),
        )

    content_top = float(canvas_config.offset_y) - float(content_outer_radius)
    content_bottom = float(canvas_config.offset_y) + float(content_outer_radius)

    if external_label_bounds is not None:
        content_top = min(content_top, float(external_label_bounds[1]))
        content_bottom = max(content_bottom, float(external_label_bounds[3]))

    return content_top, content_bottom


def _position_single_top_bottom_legend_between_edge_and_content(
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
    *,
    position: str,
    content_top: float,
    content_bottom: float,
    edge_min_px: float = SINGLE_LEGEND_EDGE_MIN_PX,
    content_gap_px: float = SINGLE_LEGEND_CONTENT_GAP_MIN_PX,
) -> None:
    """Place top/bottom legend at the midpoint between content edge and canvas edge."""
    if position not in {"top", "bottom"}:
        return

    legend_height = float(legend_config.legend_height)
    legend_local_top = -0.5 * float(legend_config.color_rect_size)
    canvas_config.legend_offset_x = (
        float(canvas_config.total_width) - float(legend_config.legend_width)
    ) / 2.0

    if position == "top":
        lane_top = float(edge_min_px)
        lane_bottom = float(content_top) - float(content_gap_px)
        free_height = lane_bottom - lane_top
        if free_height < legend_height:
            missing = legend_height - free_height
            _translate_canvas_top_level_groups(canvas, dy=missing)
            canvas_config.offset_y = float(canvas_config.offset_y) + missing
            canvas_config.total_height = float(canvas_config.total_height) + missing
            content_top = float(content_top) + missing
            _sync_canvas_viewbox(canvas, canvas_config)
            lane_bottom = float(content_top) - float(content_gap_px)
        legend_top = lane_top + max(0.0, 0.5 * (lane_bottom - lane_top - legend_height))
        canvas_config.legend_offset_y = legend_top - legend_local_top
        legend_bottom = legend_top + legend_height
        required_canvas_bottom = legend_bottom + float(edge_min_px)
        if required_canvas_bottom > float(canvas_config.total_height):
            canvas_config.total_height = required_canvas_bottom
            _sync_canvas_viewbox(canvas, canvas_config)
        return

    lane_top = float(content_bottom) + float(content_gap_px)
    lane_bottom = float(canvas_config.total_height) - float(edge_min_px)
    free_height = lane_bottom - lane_top
    if free_height < legend_height:
        missing = legend_height - free_height
        canvas_config.total_height = float(canvas_config.total_height) + missing
        _sync_canvas_viewbox(canvas, canvas_config)
        lane_bottom = float(canvas_config.total_height) - float(edge_min_px)
    legend_top = lane_top + max(0.0, 0.5 * (lane_bottom - lane_top - legend_height))
    canvas_config.legend_offset_y = legend_top - legend_local_top


def _resolve_label_legend_collisions(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> None:
    """Resolve label-vs-legend collisions with ordered fallbacks."""
    if canvas_config.legend_position == "none":
        return
    if float(legend_config.legend_width) <= 0 or float(legend_config.legend_height) <= 0:
        return

    external_labels = [label for label in labels if not label.get("is_embedded")]
    if not external_labels:
        return
    if not _labels_collide_with_legend(external_labels, total_length, canvas_config, legend_config):
        return

    if _try_shift_labels_away_from_legend(external_labels, total_length, canvas_config, legend_config):
        return
    if _try_move_legend_away_from_labels(external_labels, total_length, canvas_config, legend_config):
        return
    _expand_canvas_for_legend(external_labels, total_length, canvas_config, legend_config)


def _slot_width_ratio_factor(
    slot: CircularTrackSlot | None,
    *,
    base_radius_px: float,
    base_track_ratio: float,
) -> float | None:
    if slot is None or slot.width is None:
        return None
    denominator = float(base_radius_px) * float(base_track_ratio)
    if denominator <= 0:
        return None
    return float(slot.width.resolve(float(base_radius_px))) / denominator


def _feature_track_ratio_factor_from_draw_width(
    draw_width_px: float,
    *,
    precomputed_feature_dict: dict | None,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float | None:
    """Convert a measured feature draw band width back to renderer track_ratio_factor."""
    target_width = max(0.0, float(draw_width_px))
    if target_width <= FEATURE_BAND_EPSILON:
        return None

    if precomputed_feature_dict is not None:
        zero_band = _compute_feature_band_bounds_px(
            precomputed_feature_dict,
            total_length,
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=str(canvas_config.length_param),
            track_ratio_factor=0.0,
            cfg=cfg,
        )
        unit_band = _compute_feature_band_bounds_px(
            precomputed_feature_dict,
            total_length,
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=str(canvas_config.length_param),
            track_ratio_factor=1.0,
            cfg=cfg,
        )
        if zero_band is not None and unit_band is not None:
            zero_width = abs(float(zero_band[1]) - float(zero_band[0]))
            unit_width = abs(float(unit_band[1]) - float(unit_band[0]))
            scalable_width = unit_width - zero_width
            if scalable_width > FEATURE_BAND_EPSILON:
                return max(0.0, (target_width - zero_width) / scalable_width)

    denominator = float(canvas_config.radius) * float(canvas_config.track_ratio)
    if denominator <= FEATURE_BAND_EPSILON:
        return None
    return target_width / denominator


def _slot_dinucleotide(slot_or_resolved: CircularTrackSlot | CircularResolvedSlot, default: str) -> str:
    params = getattr(slot_or_resolved, "params", {}) or {}
    raw = params.get("nt", params.get("dinucleotide", default))
    nt = str(raw or default).upper()
    return nt if len(nt) >= 2 else str(default or "GC").upper()


def _svg_number(value: object, *, default: float = 0.0) -> float:
    if value is None:
        return float(default)
    raw = str(value).strip()
    if raw.endswith("px"):
        raw = raw[:-2]
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default)


def _text_element_plain_text(element: Any) -> str:
    parts: list[str] = []
    text = getattr(element, "text", None)
    if text:
        parts.append(str(text))
    for child in getattr(element, "elements", []) or []:
        child_text = getattr(child, "text", None)
        if child_text:
            parts.append(str(child_text))
    return "".join(parts)


def _definition_reserved_radius_px(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    species: str | None,
    strain: str | None,
    config_dict: dict,
    *,
    cfg: GbdrawConfig,
    plot_title: str | None,
    definition_profile: str,
) -> float:
    definition_group = DefinitionGroup(
        gb_record,
        canvas_config,
        config_dict,
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_profile=definition_profile,
        cfg=cfg,
    ).get_group()

    max_extent = 0.0
    for element in getattr(definition_group, "elements", []) or []:
        attribs = getattr(element, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        x_value = _svg_number(attribs.get("x"), default=0.0)
        y_value = _svg_number(attribs.get("y"), default=0.0)
        font_size = _svg_number(attribs.get("font-size"), default=cfg.objects.definition.circular.font_size)
        text = _text_element_plain_text(element)
        if not text.strip():
            continue
        width_px, height_px = calculate_bbox_dimensions(
            text,
            str(attribs.get("font-family", cfg.objects.text.font_family)),
            font_size,
            int(canvas_config.dpi),
        )
        half_width = 0.5 * float(width_px)
        half_height = 0.5 * float(height_px if height_px else font_size)
        max_extent = max(
            max_extent,
            abs(x_value - half_width),
            abs(x_value + half_width),
            abs(y_value - half_height),
            abs(y_value + half_height),
        )

    if max_extent <= FEATURE_BAND_EPSILON:
        return 0.0
    padding_px = max(8.0, 0.02 * float(canvas_config.radius))
    return float(max_extent + padding_px)


def _slot_config_with_dinucleotide(config: Any, nt: str) -> Any:
    if str(getattr(config, "dinucleotide", "")).upper() == str(nt).upper():
        return config
    cloned = copy.copy(config)
    cloned.dinucleotide = str(nt).upper()
    return cloned


def _slot_gc_config_with_axis_font_size(
    config: GcContentConfigurator,
    nt: str,
    tick_font_size: float | None,
) -> GcContentConfigurator:
    slot_config = _slot_config_with_dinucleotide(config, nt)
    if tick_font_size is None or getattr(slot_config, "tick_font_size", None) is not None:
        return slot_config
    cloned = copy.copy(slot_config)
    cloned.tick_font_size = float(tick_font_size)
    return cloned


def _gc_content_matches_depth_axis_font_size(
    *,
    gc_config: GcContentConfigurator,
    depth_config: DepthConfigurator | None,
    resolved_track_slots: Sequence[CircularResolvedSlot],
) -> float | None:
    if depth_config is None:
        return None
    if str(getattr(gc_config, "mode", "deviation")).strip().lower() != "percent":
        return None
    if getattr(gc_config, "tick_font_size", None) is not None:
        return None
    if not bool(getattr(depth_config, "show_axis", True)) or not bool(getattr(depth_config, "show_ticks", True)):
        return None
    depth_slot = next(
        (
            slot
            for slot in resolved_track_slots
            if str(slot.renderer) == "depth" and float(slot.draw_width_px) > 0
        ),
        None,
    )
    if depth_slot is None:
        return None
    return depth_axis_tick_font_size_px(depth_config, float(depth_slot.draw_width_px))


def _slot_dataframe_for_nt(
    *,
    nt: str,
    default_df: DataFrame,
    default_nt: str,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
) -> DataFrame:
    normalized_nt = str(nt).upper()
    if dinucleotide_dataframes and normalized_nt in dinucleotide_dataframes:
        return dinucleotide_dataframes[normalized_nt]
    if normalized_nt == str(default_nt).upper():
        return default_df
    return DataFrame()


def _unique_legend_key(legend_table: dict, preferred: str) -> str:
    if preferred not in legend_table:
        return preferred
    suffix = 2
    while f"{preferred} ({suffix})" in legend_table:
        suffix += 1
    return f"{preferred} ({suffix})"


def _slot_legend_label(slot: CircularTrackSlot, fallback: str) -> str:
    params = slot.params or {}
    raw_label = params.get("legend_label", params.get("label"))
    label = str(raw_label).strip() if raw_label is not None else ""
    return label or fallback


def _sync_legend_table_for_circular_slots(
    legend_table: dict,
    *,
    circular_track_slots: list[CircularTrackSlot] | None,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
    depth_df: DataFrame | None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    cfg: GbdrawConfig,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float | None = None,
) -> dict:
    """Replace singleton numeric legend entries with slot-aware entries."""
    if circular_track_slots is None:
        return legend_table

    out = dict(legend_table)
    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    removable = {
        "Depth",
        f"{default_nt} content",
        f"{default_nt} content (+)",
        f"{default_nt} content (-)",
        f"{default_nt} skew",
        f"{default_nt} skew (+)",
        f"{default_nt} skew (-)",
    }
    for key in removable:
        out.pop(key, None)

    for slot in circular_track_slots:
        if not slot.enabled:
            continue
        renderer = str(slot.renderer)
        if renderer == "depth":
            if depth_tracks:
                try:
                    track_index = int((slot.params or {}).get("track_index", 0) or 0)
                except (TypeError, ValueError):
                    track_index = 0
                if track_index < 0 or track_index >= len(depth_tracks):
                    continue
                track = depth_tracks[track_index]
                label = _slot_legend_label(slot, track.label)
                out[_unique_legend_key(out, label)] = {
                    "type": "solid",
                    "fill": track.config.fill_color,
                    "stroke": track.config.stroke_color,
                    "width": track.config.stroke_width,
                }
                continue
            if depth_config is None or depth_df is None:
                continue
            label = _slot_legend_label(slot, "Depth")
            out[_unique_legend_key(out, label)] = {
                "type": "solid",
                "fill": depth_config.fill_color,
                "stroke": depth_config.stroke_color,
                "width": depth_config.stroke_width,
            }
        elif renderer == "dinucleotide_content":
            nt = _slot_dinucleotide(slot, default_nt)
            label = _slot_legend_label(slot, f"{nt} content")
            if gc_config.high_fill_color == gc_config.low_fill_color:
                out[_unique_legend_key(out, label)] = {
                    "type": "solid",
                    "fill": gc_config.high_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
            else:
                out[_unique_legend_key(out, f"{label} (+)")] = {
                    "type": "solid",
                    "fill": gc_config.high_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
                out[_unique_legend_key(out, f"{label} (-)")] = {
                    "type": "solid",
                    "fill": gc_config.low_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
        elif renderer == "dinucleotide_skew":
            nt = _slot_dinucleotide(slot, default_nt)
            label = _slot_legend_label(slot, f"{nt} skew")
            if skew_config.high_fill_color == skew_config.low_fill_color:
                out[_unique_legend_key(out, label)] = {
                    "type": "solid",
                    "fill": skew_config.high_fill_color,
                    "stroke": skew_config.stroke_color,
                    "width": skew_config.stroke_width,
                }
            else:
                out[_unique_legend_key(out, f"{label} (+)")] = {
                    "type": "solid",
                    "fill": skew_config.high_fill_color,
                    "stroke": skew_config.stroke_color,
                    "width": skew_config.stroke_width,
                }
                out[_unique_legend_key(out, f"{label} (-)")] = {
                    "type": "solid",
                    "fill": skew_config.low_fill_color,
                    "stroke": skew_config.stroke_color,
                    "width": skew_config.stroke_width,
                }
    if conservation_tracks:
        if any(track.track_color for track in conservation_tracks):
            for track in conservation_tracks:
                min_color, max_color = conservation_track_gradient_colors(
                    track.track_color,
                    default_min_color=cfg.objects.conservation.min_color,
                    default_max_color=cfg.objects.conservation.max_color,
                )
                out[_unique_legend_key(out, track.track_label)] = {
                    "type": "gradient",
                    "min_color": min_color,
                    "max_color": max_color,
                    "stroke": "none",
                    "width": 0,
                    "min_value": float(conservation_min_identity or 0.0),
                }
        else:
            out[_unique_legend_key(out, "Conservation identity")] = {
                "type": "gradient",
                "min_color": cfg.objects.conservation.min_color,
                "max_color": cfg.objects.conservation.max_color,
                "stroke": "none",
                "width": 0,
                "min_value": float(conservation_min_identity or 0.0),
            }
    return out


def _draw_resolved_circular_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    config_dict: dict,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks: Sequence[DepthTrackData] | None,
    conservation_tracks: Sequence[ConservationTrack] | None,
    conservation_min_identity: float,
    cfg: GbdrawConfig,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
    gc_content_tick_font_size_override: float | None = None,
    precomputed_feature_dict: dict | None = None,
    precalculated_labels: list[dict] | None = None,
    _tick_track_channel_override: str | None = None,
    use_slot_group_id: bool = True,
    use_slot_tick_options: bool = True,
    use_feature_anchor_override: bool = True,
) -> Drawing:
    """Draw one resolved circular slot."""
    renderer = str(resolved_slot.renderer)
    norm_factor_override = float(resolved_slot.anchor_radius_px) / float(canvas_config.radius)
    if renderer == "spacer":
        return canvas

    if renderer == "features":
        base_width = (
            float(canvas_config.radius)
            * float(canvas_config.track_ratio)
            * float(cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)][0])
        )
        resolved_feature_width = float(resolved_slot.draw_width_px)
        radial_layout = getattr(canvas_config, "circular_radial_layout", None)
        if radial_layout is not None and getattr(radial_layout, "features", None) is not None:
            resolved_feature_width = float(radial_layout.features.width_px)
        ratio_override = None
        if base_width > 0 and resolved_feature_width > 0:
            ratio_override = resolved_feature_width / max(
                FEATURE_BAND_EPSILON,
                float(canvas_config.radius) * float(canvas_config.track_ratio),
            )
        feature_kwargs: dict[str, Any] = {
            "cfg": cfg,
            "precomputed_feature_dict": precomputed_feature_dict,
            "precalculated_labels": precalculated_labels,
            "feature_track_ratio_factor_override": ratio_override,
        }
        if use_feature_anchor_override or not math.isclose(
            float(resolved_slot.anchor_radius_px),
            float(canvas_config.radius),
            rel_tol=1e-9,
            abs_tol=1e-9,
        ):
            feature_kwargs["feature_anchor_radius_px"] = float(resolved_slot.anchor_radius_px)
        return add_record_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            **feature_kwargs,
        )

    if renderer == "ticks":
        tick_layout = (
            resolved_slot.payload
            if isinstance(resolved_slot.payload, CircularTickLayout)
            else None
        )
        label_side = (
            str(tick_layout.label_side)
            if tick_layout is not None
            else tick_sides_for_tick_label_layout(
                resolved_slot.params.get("tick_label_layout"),
                side=resolved_slot.side,
            )[0]
        ).strip().lower()
        tick_side = (
            str(tick_layout.tick_side)
            if tick_layout is not None
            else tick_sides_for_tick_label_layout(
                resolved_slot.params.get("tick_label_layout"),
                side=resolved_slot.side,
            )[1]
        ).strip().lower()
        if label_side != "none" or tick_side != "none":
            tick_group_kwargs: dict[str, Any] = {
                "radius_override": float(resolved_slot.anchor_radius_px),
                "cfg": cfg,
            }
            if use_slot_tick_options or tick_layout is not None:
                tick_group_kwargs["label_side"] = label_side
                tick_group_kwargs["tick_side"] = tick_side
                tick_group_kwargs["tick_length_px"] = (
                    tick_layout.tick_length_px
                    if tick_layout is not None
                    else (
                        float(resolved_slot.draw_width_px)
                        if resolved_slot.explicit_width and resolved_slot.draw_width_px > 0
                        else None
                    )
                )
                tick_group_kwargs["track_preset"] = (
                    tick_layout.track_preset
                    if tick_layout is not None
                    else normalize_circular_track_preset(
                        str(resolved_slot.params.get("preset", resolved_slot.params.get("track_preset", "tuckin")))
                    )
                )
            if _tick_track_channel_override is not None:
                tick_group_kwargs["tick_track_channel_override"] = _tick_track_channel_override
            canvas = add_tick_group_on_canvas(
                canvas,
                gb_record,
                canvas_config,
                config_dict,
                **tick_group_kwargs,
            )
        return canvas

    if renderer == "depth":
        depth_group_id = str(resolved_slot.id)
        selected_depth_df = depth_df
        selected_depth_config = depth_config
        if depth_tracks:
            try:
                track_index = int(resolved_slot.params.get("track_index", 0) or 0)
            except (TypeError, ValueError):
                track_index = 0
            if track_index < 0 or track_index >= len(depth_tracks):
                logger.warning(
                    "Skipping circular depth slot '%s' because track_index=%s is unavailable.",
                    resolved_slot.id,
                    track_index,
                )
                return canvas
            depth_track = depth_tracks[track_index]
            selected_depth_df = depth_track.df
            selected_depth_config = depth_track.config
            if not use_slot_group_id:
                depth_group_id = depth_track.id
        if selected_depth_config is None or selected_depth_df is None:
            logger.warning("Skipping circular depth slot '%s' because depth data are unavailable.", resolved_slot.id)
            return canvas
        depth_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
            "cfg": cfg,
        }
        if use_slot_group_id:
            depth_kwargs["group_id"] = depth_group_id
        return add_depth_group_on_canvas(
            canvas,
            gb_record,
            selected_depth_df,
            canvas_config,
            selected_depth_config,
            config_dict,
            **depth_kwargs,
        )

    if renderer == "sequence_conservation":
        track_index = int(resolved_slot.params.get("track_index", 0) or 0)
        track = next(
            (
                conservation_track
                for conservation_track in (conservation_tracks or ())
                if int(conservation_track.track_index) == track_index
            ),
            None,
        )
        if track is None:
            logger.warning(
                "Skipping circular conservation slot '%s' because conservation data are unavailable.",
                resolved_slot.id,
            )
            return canvas
        return add_conservation_group_on_canvas(
            canvas,
            gb_record,
            track,
            canvas_config,
            inner_radius_px=float(resolved_slot.draw_inner_radius_px),
            outer_radius_px=float(resolved_slot.draw_outer_radius_px),
            min_identity=float(conservation_min_identity),
            cfg=cfg,
        )

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = _slot_dinucleotide(resolved_slot, default_nt)
    slot_df = _slot_dataframe_for_nt(
        nt=nt,
        default_df=gc_df,
        default_nt=default_nt,
        dinucleotide_dataframes=dinucleotide_dataframes,
    )

    if renderer == "dinucleotide_content":
        if f"{nt} content" not in slot_df.columns:
            logger.warning(
                "Skipping circular content slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        gc_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
            "cfg": cfg,
        }
        if use_slot_group_id:
            gc_kwargs["group_id"] = str(resolved_slot.id)
        slot_gc_config = _slot_gc_config_with_axis_font_size(
            gc_config,
            nt,
            gc_content_tick_font_size_override,
        )
        return add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_gc_config,
            config_dict,
            **gc_kwargs,
        )

    if renderer == "dinucleotide_skew":
        if f"{nt} skew" not in slot_df.columns:
            logger.warning(
                "Skipping circular skew slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        skew_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
            "cfg": cfg,
        }
        if use_slot_group_id:
            skew_kwargs["group_id"] = str(resolved_slot.id)
        return add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            _slot_config_with_dinucleotide(skew_config, nt),
            config_dict,
            **skew_kwargs,
        )

    logger.warning("Skipping unsupported circular track slot renderer '%s'.", renderer)
    return canvas


def _default_gc_skew_track_ids_without_depth(*, show_gc: bool, show_skew: bool) -> dict[str, int]:
    """Return built-in GC/skew track ids for the same visibility with depth disabled."""
    track_ids: dict[str, int] = {}
    if show_gc:
        track_ids["gc_content"] = 2
    if show_skew:
        track_ids["gc_skew"] = 3 if show_gc else 2
    return track_ids


def _default_gc_skew_track_width_px(
    kind: str,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    """Return the default width for a built-in circular GC or skew track."""
    length_param = str(canvas_config.length_param)
    factor_index = 1 if kind == "gc_content" else 2
    return (
        float(canvas_config.radius)
        * float(canvas_config.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors[length_param][factor_index])
    )


def _default_depth_track_width_px(
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    length_param = str(canvas_config.length_param)
    return (
        float(canvas_config.radius)
        * float(canvas_config.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors[length_param][1])
        * 0.5
    )


def _default_numeric_slot_width_px(
    slot: CircularTrackSlot,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    renderer = str(slot.renderer)
    if renderer == "depth":
        return _default_depth_track_width_px(canvas_config=canvas_config, cfg=cfg)
    if renderer == "dinucleotide_skew":
        return _default_gc_skew_track_width_px(
            "gc_skew",
            canvas_config=canvas_config,
            cfg=cfg,
        )
    return _default_gc_skew_track_width_px(
        "gc_content",
        canvas_config=canvas_config,
        cfg=cfg,
    )


def _default_gc_skew_layout_without_depth(
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    show_gc: bool,
    show_skew: bool,
    radius_mapper: Callable[[float], float] | None = None,
) -> dict[str, tuple[float, float]]:
    """Return default no-depth GC/skew layout as {kind: (center_px, width_px)}."""
    track_ids = _default_gc_skew_track_ids_without_depth(
        show_gc=show_gc,
        show_skew=show_skew,
    )
    if not track_ids:
        return {}

    length_param = str(canvas_config.length_param)
    track_type = _circular_preset_for_layout(canvas_config, cfg)
    track_dict = cfg.canvas.circular.track_dict[length_param][track_type]
    layout: dict[str, tuple[float, float]] = {}
    for kind, track_id in track_ids.items():
        center_px = float(canvas_config.radius) * float(track_dict[str(track_id)])
        if radius_mapper is not None:
            center_px = float(radius_mapper(center_px))
        width_px = _default_gc_skew_track_width_px(
            kind,
            canvas_config=canvas_config,
            cfg=cfg,
        )
        layout[kind] = (float(center_px), float(width_px))
    return layout


def _default_gc_skew_gap_px(
    layout: dict[str, tuple[float, float]],
    *,
    canvas_config: CircularCanvasConfigurator,
) -> float:
    """Return the default radial gap between adjacent GC/skew tracks."""
    fallback_gap = max(1.0, 0.01 * float(canvas_config.radius))
    if len(layout) < 2:
        return fallback_gap

    annuli = sorted(
        (
            (
                float(center_px) - (0.5 * float(width_px)),
                float(center_px) + (0.5 * float(width_px)),
            )
            for center_px, width_px in layout.values()
        ),
        key=lambda annulus: annulus[0],
    )
    gaps = [
        max(0.0, float(annuli[idx + 1][0]) - float(annuli[idx][1]))
        for idx in range(len(annuli) - 1)
    ]
    positive_gaps = [gap for gap in gaps if gap > FEATURE_BAND_EPSILON]
    if not positive_gaps:
        return 0.0
    return min(positive_gaps)


def _distributed_lane_specs_between_bounds(
    *,
    inner_radius_px: float,
    outer_radius_px: float,
    default_widths_outer_to_inner: Sequence[float],
    desired_gap_px: float,
    include_outer_boundary_gap: bool = False,
    fill_available: bool = False,
) -> list[tuple[float, float]]:
    """Return outer-to-inner lane specs as (center_px, width_scale)."""
    lane_count = len(default_widths_outer_to_inner)
    if lane_count <= 0:
        return []

    inner = float(inner_radius_px)
    outer = float(outer_radius_px)
    if outer < inner:
        inner, outer = outer, inner
    available_width = max(0.0, outer - inner)

    default_widths_inner_to_outer = [
        max(0.0, float(width_px))
        for width_px in reversed(default_widths_outer_to_inner)
    ]
    total_default_width = sum(default_widths_inner_to_outer)
    if total_default_width <= FEATURE_BAND_EPSILON:
        return []

    gap_count = lane_count if include_outer_boundary_gap else max(0, lane_count - 1)
    desired_gap = max(0.0, float(desired_gap_px))
    if gap_count > 0:
        compressed_gap_px = min(
            desired_gap,
            max(0.0, available_width / (3.0 * float(gap_count))),
        )
    else:
        compressed_gap_px = 0.0

    compressed_width_budget = max(0.0, available_width - (float(gap_count) * compressed_gap_px))
    if compressed_width_budget < total_default_width - FEATURE_BAND_EPSILON:
        gap_px = compressed_gap_px
        width_scale = max(0.0, compressed_width_budget / total_default_width)
    else:
        width_scale = 1.0
        if fill_available and gap_count > 0:
            gap_px = max(0.0, (available_width - total_default_width) / float(gap_count))
        else:
            gap_px = compressed_gap_px

    specs_inner_to_outer: list[tuple[float, float]] = []
    cursor = inner
    for width_px in default_widths_inner_to_outer:
        scaled_width_px = width_px * width_scale
        center_px = cursor + (0.5 * scaled_width_px)
        specs_inner_to_outer.append((float(center_px), float(width_scale)))
        cursor += scaled_width_px + gap_px

    return list(reversed(specs_inner_to_outer))


def _compute_feature_band_bounds_px(
    feature_dict: dict | None,
    total_length: int,
    *,
    base_radius_px: float,
    track_ratio: float,
    length_param: str,
    track_ratio_factor: float,
    cfg: GbdrawConfig,
    track_id_whitelist: set[int] | None = None,
    preset: str | None = None,
) -> tuple[float, float] | None:
    """Compute (inner_radius_px, outer_radius_px) over all drawn feature blocks."""
    if not feature_dict or total_length <= 0:
        return None

    filtered_feature_dict = feature_dict
    if track_id_whitelist is not None:
        filtered_feature_dict = {
            feature_id: feature_object
            for feature_id, feature_object in feature_dict.items()
            if int(getattr(feature_object, "feature_track_id", 0)) in track_id_whitelist
        }
        if not filtered_feature_dict:
            return None

    feature_width_px = float(base_radius_px) * float(track_ratio) * float(track_ratio_factor)
    feature_layout = build_circular_feature_layout(
        filtered_feature_dict,
        axis_radius_px=float(base_radius_px),
        width_px=feature_width_px,
        lane_direction=circular_feature_lane_direction_for_preset(
            normalize_circular_track_preset(preset or cfg.canvas.circular.track_type)
        ),
        strandedness=bool(cfg.canvas.strandedness),
    )
    if feature_layout is None:
        return None
    return float(feature_layout.all_band_px.inner_px), float(feature_layout.all_band_px.outer_px)


def _map_radius_across_feature_bands(
    radius_px: float,
    old_band: tuple[float, float],
    new_band: tuple[float, float],
) -> float:
    """Map radius while preserving distance outside/inside feature band edges."""
    old_inner, old_outer = old_band
    new_inner, new_outer = new_band

    if old_outer <= old_inner + FEATURE_BAND_EPSILON:
        return float(radius_px) + (new_outer - old_outer)

    if radius_px >= old_outer:
        return new_outer + (float(radius_px) - old_outer)
    if radius_px <= old_inner:
        return new_inner - (old_inner - float(radius_px))

    rel = (float(radius_px) - old_inner) / (old_outer - old_inner)
    return new_inner + rel * (new_outer - new_inner)


def _build_feature_radius_mapper(
    feature_dict: dict | None,
    total_length: int,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    feature_track_ratio_factor_override: float,
) -> tuple[Callable[[float], float] | None, tuple[float, float] | None, tuple[float, float] | None]:
    """Create radius mapper from default feature band -> overridden feature band.

    Use primary feature track (track_id=0) as baseline to avoid sparse
    overlap-resolved outer tracks from disproportionately shifting unrelated
    auto tracks (labels/ticks/gc/skew).
    """
    if feature_dict is None:
        return None, None, None

    length_param = str(canvas_config.length_param)
    default_ratio_factor = float(cfg.canvas.circular.track_ratio_factors[length_param][0])

    old_band = _compute_feature_band_bounds_px(
        feature_dict,
        total_length,
        base_radius_px=float(canvas_config.radius),
        track_ratio=float(canvas_config.track_ratio),
        length_param=str(canvas_config.length_param),
        track_ratio_factor=default_ratio_factor,
        cfg=cfg,
        track_id_whitelist={0},
    )
    new_band = _compute_feature_band_bounds_px(
        feature_dict,
        total_length,
        base_radius_px=float(canvas_config.radius),
        track_ratio=float(canvas_config.track_ratio),
        length_param=str(canvas_config.length_param),
        track_ratio_factor=float(feature_track_ratio_factor_override),
        cfg=cfg,
        track_id_whitelist={0},
    )

    if old_band is None or new_band is None:
        # Fallback for datasets where all features were assigned to non-zero tracks.
        old_band = _compute_feature_band_bounds_px(
            feature_dict,
            total_length,
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=str(canvas_config.length_param),
            track_ratio_factor=default_ratio_factor,
            cfg=cfg,
        )
        new_band = _compute_feature_band_bounds_px(
            feature_dict,
            total_length,
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=str(canvas_config.length_param),
            track_ratio_factor=float(feature_track_ratio_factor_override),
            cfg=cfg,
        )

    if old_band is None or new_band is None:
        return None, old_band, new_band

    def mapper(radius_px: float) -> float:
        return _map_radius_across_feature_bands(float(radius_px), old_band, new_band)

    return mapper, old_band, new_band


def _default_track_center_radius_px(
    kind: str,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float | None:
    """Return default center radius for a built-in circular track kind."""
    if kind in {"axis", "ticks"}:
        return float(canvas_config.radius)

    if kind == "depth":
        track_id = canvas_config.track_ids.get("depth_track")
    elif kind == "gc_content":
        track_id = canvas_config.track_ids.get("gc_track")
    elif kind == "gc_skew":
        track_id = canvas_config.track_ids.get("skew_track")
    else:
        track_id = None

    if track_id is None:
        return None

    length_param = str(canvas_config.length_param)
    track_type = _circular_preset_for_layout(canvas_config, cfg)
    norm_factor = float(cfg.canvas.circular.track_dict[length_param][track_type][str(track_id)])
    return float(canvas_config.radius) * norm_factor


def _default_outer_label_arena(
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[float, float]:
    """Return default (anchor_radius_px, arc_outer_radius_px) for outer labels."""
    length_param = str(canvas_config.length_param)
    track_type = _circular_preset_for_layout(canvas_config, cfg)
    strands = "separate" if cfg.canvas.strandedness else "single"
    base_radius = float(canvas_config.radius)

    anchor_radius = float(cfg.labels.radius_factor[track_type][strands][length_param]) * base_radius
    offset_cfg = cfg.labels.unified_adjustment.outer_labels
    arc_x = (
        base_radius
        * float(cfg.labels.arc_x_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.x_radius_offset)
    )
    arc_y = (
        base_radius
        * float(cfg.labels.arc_y_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.y_radius_offset)
    )
    arc_outer = max(abs(arc_x), abs(arc_y))

    if arc_outer < anchor_radius:
        arc_outer = anchor_radius
    return anchor_radius, arc_outer


def _arena_from_center_and_width(center_px: float | None, width_px: float | None) -> tuple[float, float] | None:
    """Convert center+width to (inner, outer) arena bounds."""
    if center_px is None or width_px is None:
        return None
    inner_px = float(center_px) - (float(width_px) / 2.0)
    outer_px = float(center_px) + (float(width_px) / 2.0)
    if outer_px < inner_px:
        inner_px, outer_px = outer_px, inner_px
    return inner_px, outer_px


def _annulus_from_center_and_width(center_px: float, width_px: float) -> tuple[float, float]:
    """Return annulus bounds from center radius and width."""
    half_width = max(0.0, 0.5 * float(width_px))
    inner = float(center_px) - half_width
    outer = float(center_px) + half_width
    if outer < inner:
        inner, outer = outer, inner
    return inner, outer


def _annulus_overlaps_band(
    annulus: tuple[float, float],
    band: tuple[float, float],
) -> bool:
    """Whether annulus intersects a forbidden band."""
    annulus_inner, annulus_outer = sorted((float(annulus[0]), float(annulus[1])))
    band_inner, band_outer = sorted((float(band[0]), float(band[1])))
    return (
        annulus_inner < (band_outer - FEATURE_BAND_EPSILON)
        and annulus_outer > (band_inner + FEATURE_BAND_EPSILON)
    )


def _annulus_overlaps_any_band(
    annulus: tuple[float, float],
    forbidden_bands: list[tuple[float, float]],
) -> bool:
    """Whether annulus intersects any forbidden band."""
    return any(_annulus_overlaps_band(annulus, band) for band in forbidden_bands)


def _resolve_ratio_annulus_center_avoiding_forbidden_bands(
    *,
    current_center_px: float,
    ratio_bounds: tuple[float, float],
    forbidden_bands: list[tuple[float, float]],
    prefer_inside: bool = True,
) -> float | None:
    """Resolve center for ratio-scaled annulus (ticks), prioritizing inward moves."""
    ratio_min, ratio_max = sorted((float(ratio_bounds[0]), float(ratio_bounds[1])))
    if ratio_min <= 0.0 or ratio_max <= 0.0:
        return None

    if not forbidden_bands:
        return float(current_center_px)

    def annulus_for_center(center_px: float) -> tuple[float, float]:
        return _tick_annulus_for_center(center_px, (ratio_min, ratio_max))

    current = float(current_center_px)
    if current > FEATURE_BAND_EPSILON and not _annulus_overlaps_any_band(
        annulus_for_center(current), forbidden_bands
    ):
        return current

    inside_candidates: list[float] = []
    outside_candidates: list[float] = []
    for band_inner, band_outer in forbidden_bands:
        low, high = sorted((float(band_inner), float(band_outer)))
        inside_center = (low - FEATURE_BAND_EPSILON) / ratio_max
        outside_center = (high + FEATURE_BAND_EPSILON) / ratio_min
        if inside_center > FEATURE_BAND_EPSILON:
            inside_candidates.append(float(inside_center))
        if outside_center > FEATURE_BAND_EPSILON:
            outside_candidates.append(float(outside_center))

    def pick_valid(candidates: list[float]) -> float | None:
        for candidate in sorted(candidates, key=lambda value: (abs(value - current), value)):
            if _annulus_overlaps_any_band(annulus_for_center(candidate), forbidden_bands):
                continue
            return float(candidate)
        return None

    if prefer_inside:
        resolved = pick_valid(inside_candidates)
        if resolved is not None:
            return resolved
        return pick_valid(outside_candidates)

    combined_candidates = inside_candidates + outside_candidates
    return pick_valid(combined_candidates)


def _resolve_fixed_width_annulus_center_avoiding_forbidden_bands(
    *,
    current_center_px: float,
    width_px: float,
    forbidden_bands: list[tuple[float, float]],
    prefer_inside: bool = True,
) -> float | None:
    """Resolve center for fixed-width annulus (GC tracks), prioritizing inward moves."""
    if width_px < 0:
        return None
    if not forbidden_bands:
        return float(current_center_px)

    current = float(current_center_px)
    annulus = _annulus_from_center_and_width(current, float(width_px))
    if current > FEATURE_BAND_EPSILON and not _annulus_overlaps_any_band(annulus, forbidden_bands):
        return current

    half_width = max(0.0, 0.5 * float(width_px))
    inside_candidates: list[float] = []
    outside_candidates: list[float] = []
    for band_inner, band_outer in forbidden_bands:
        low, high = sorted((float(band_inner), float(band_outer)))
        inside_center = low - half_width - FEATURE_BAND_EPSILON
        outside_center = high + half_width + FEATURE_BAND_EPSILON
        if inside_center > FEATURE_BAND_EPSILON:
            inside_candidates.append(float(inside_center))
        if outside_center > FEATURE_BAND_EPSILON:
            outside_candidates.append(float(outside_center))

    def pick_valid(candidates: list[float]) -> float | None:
        for candidate in sorted(candidates, key=lambda value: (abs(value - current), value)):
            candidate_annulus = _annulus_from_center_and_width(candidate, float(width_px))
            if _annulus_overlaps_any_band(candidate_annulus, forbidden_bands):
                continue
            return float(candidate)
        return None

    if prefer_inside:
        resolved = pick_valid(inside_candidates)
        if resolved is not None:
            return resolved
        return pick_valid(outside_candidates)

    return pick_valid(inside_candidates + outside_candidates)


def _shrink_fixed_width_annulus_to_avoid_forbidden_bands(
    *,
    center_px: float,
    current_width_px: float,
    forbidden_bands: list[tuple[float, float]],
) -> float | None:
    """Shrink fixed-width annulus around center so it avoids all forbidden bands."""
    if current_width_px < 0:
        return None
    if not forbidden_bands:
        return float(current_width_px)

    center = float(center_px)
    if center <= FEATURE_BAND_EPSILON:
        return None

    max_half_width = float("inf")
    for band_inner, band_outer in forbidden_bands:
        low, high = sorted((float(band_inner), float(band_outer)))
        if center <= low:
            max_half_width = min(max_half_width, low - FEATURE_BAND_EPSILON - center)
        elif center >= high:
            max_half_width = min(max_half_width, center - (high + FEATURE_BAND_EPSILON))
        else:
            # Center inside forbidden band: shrinking alone cannot resolve overlap.
            return None

    if not math.isfinite(max_half_width):
        return float(current_width_px)
    if max_half_width < 0.0:
        return None

    shrunk_width = min(float(current_width_px), 2.0 * max_half_width)
    if shrunk_width < 0.0:
        return None
    return float(shrunk_width)


def _tick_annulus_for_center(
    center_radius_px: float,
    tick_ratio_bounds: tuple[float, float],
) -> tuple[float, float]:
    """Return tick annulus (inner, outer) for a center radius."""
    min_ratio, max_ratio = sorted((float(tick_ratio_bounds[0]), float(tick_ratio_bounds[1])))
    inner_radius = float(center_radius_px) * min_ratio
    outer_radius = float(center_radius_px) * max_ratio
    if outer_radius < inner_radius:
        inner_radius, outer_radius = outer_radius, inner_radius
    return inner_radius, outer_radius


def _tick_annulus_overlaps_feature_band(
    center_radius_px: float,
    tick_ratio_bounds: tuple[float, float],
    feature_band_px: tuple[float, float],
) -> bool:
    """Whether the tick annulus intersects the feature band."""
    tick_inner, tick_outer = _tick_annulus_for_center(center_radius_px, tick_ratio_bounds)
    feature_inner, feature_outer = sorted((float(feature_band_px[0]), float(feature_band_px[1])))
    return (
        tick_inner < (feature_outer - FEATURE_BAND_EPSILON)
        and tick_outer > (feature_inner + FEATURE_BAND_EPSILON)
    )


def _resolve_ticks_center_radius_avoiding_feature_band(
    *,
    current_ticks_radius_px: float,
    default_ticks_radius_px: float,
    feature_band_px: tuple[float, float],
    default_feature_band_px: tuple[float, float],
    tick_ratio_bounds: tuple[float, float],
) -> float | None:
    """Resolve a non-overlapping tick center while preserving default band gaps when possible."""
    tick_min_ratio, tick_max_ratio = sorted((float(tick_ratio_bounds[0]), float(tick_ratio_bounds[1])))
    if tick_min_ratio <= 0.0 or tick_max_ratio <= 0.0:
        return None

    feature_inner, feature_outer = sorted((float(feature_band_px[0]), float(feature_band_px[1])))
    default_feature_inner, default_feature_outer = sorted(
        (float(default_feature_band_px[0]), float(default_feature_band_px[1]))
    )

    default_tick_inner, default_tick_outer = _tick_annulus_for_center(default_ticks_radius_px, tick_ratio_bounds)
    base_inside_gap = default_feature_inner - default_tick_outer
    base_outside_gap = default_tick_inner - default_feature_outer
    inside_gap = max(0.0, base_inside_gap)
    outside_gap = max(0.0, base_outside_gap)

    candidates: list[float] = []

    inside_center = (feature_inner - inside_gap - FEATURE_BAND_EPSILON) / tick_max_ratio
    if inside_center > FEATURE_BAND_EPSILON:
        if not _tick_annulus_overlaps_feature_band(inside_center, tick_ratio_bounds, feature_band_px):
            candidates.append(float(inside_center))

    outside_center = (feature_outer + outside_gap + FEATURE_BAND_EPSILON) / tick_min_ratio
    if outside_center > FEATURE_BAND_EPSILON:
        if not _tick_annulus_overlaps_feature_band(outside_center, tick_ratio_bounds, feature_band_px):
            candidates.append(float(outside_center))

    if not candidates:
        return None
    return min(candidates, key=lambda candidate: abs(float(candidate) - float(current_ticks_radius_px)))


def add_record_on_circular_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    gc_df: DataFrame,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    config_dict: dict,
    legend_config,
    legend_table,
    *,
    depth_config: DepthConfigurator | None = None,
    depth_df: DataFrame | None = None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    cfg: GbdrawConfig | None = None,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    _tick_track_channel_override: str | None = None,
) -> Drawing:
    """
    Adds various record-related groups to a circular canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    species (str): Species name.
    strain (str): Strain name.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with all record-related groups added.
    """
    cfg = cfg or canvas_config._cfg
    effective_circular_track_slots = circular_track_slots
    user_slot_mode = effective_circular_track_slots is not None
    user_active_slot_renderers = {
        str(slot.renderer)
        for slot in (effective_circular_track_slots or [])
        if slot.enabled
    }

    raw_show_labels = cfg.canvas.show_labels
    show_labels_base = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)
    depth_enabled = bool(
        canvas_config.show_depth
        and (
            bool(depth_tracks)
            or (depth_config is not None and depth_df is not None)
        )
    )
    if user_slot_mode:
        show_depth_track = bool("depth" in user_active_slot_renderers and depth_enabled)
        show_gc_track = bool("dinucleotide_content" in user_active_slot_renderers)
        show_skew_track = bool("dinucleotide_skew" in user_active_slot_renderers)
        show_ticks_track = bool("ticks" in user_active_slot_renderers)
    else:
        show_depth_track = bool(depth_enabled)
        show_gc_track = bool(canvas_config.show_gc)
        show_skew_track = bool(canvas_config.show_skew)
        show_ticks_track = True
    circular_preset = normalize_circular_track_preset(cfg.canvas.circular.track_type)
    setattr(canvas_config, "circular_track_preset", circular_preset)

    if user_slot_mode:
        show_features = "features" in user_active_slot_renderers
        radial_plan = circular_track_slots_from_preset_order(
            list(effective_circular_track_slots or []),
            circular_preset,
            CircularPresetContext(
                cfg=cfg,
                canvas_config=canvas_config,
                total_length=len(gb_record.seq),
                strandedness=bool(cfg.canvas.strandedness),
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
                tick_track_channel_override=_tick_track_channel_override,
            ),
            axis_index=circular_track_axis_index,
        )
        layout_slots = list(radial_plan.slots)
        preferred_anchor_slot_ids = radial_plan.preferred_anchor_slot_ids
    else:
        show_features = True
        radial_plan = circular_radial_plan_for_preset(
            circular_preset,
            CircularPresetContext(
                cfg=cfg,
                canvas_config=canvas_config,
                total_length=len(gb_record.seq),
                strandedness=bool(cfg.canvas.strandedness),
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
                tick_track_channel_override=_tick_track_channel_override,
            ),
        )
        layout_slots = list(radial_plan.slots)
        preferred_anchor_slot_ids = radial_plan.preferred_anchor_slot_ids
    feature_slot = next(
        (
            slot
            for slot in layout_slots
            if slot.enabled and str(slot.renderer) == "features"
        ),
        None,
    )
    feature_lane_direction = _lane_direction_for_feature_slot(
        feature_slot,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    setattr(canvas_config, "circular_feature_lane_direction", feature_lane_direction)

    show_external_labels = show_labels_base and show_features
    core_track_overlap_relayout_enabled = (
        show_features
        and bool(cfg.canvas.resolve_overlaps)
        and (not bool(cfg.canvas.strandedness))
    )
    split_overlaps_by_strand = (
        bool(cfg.canvas.resolve_overlaps)
        and (not bool(cfg.canvas.strandedness))
        and feature_lane_direction == "split"
    )

    feature_track_ratio_factor_override: float | None = None
    if show_features:
        feature_track_ratio_factor_override = _slot_width_ratio_factor(
            feature_slot,
            base_radius_px=float(canvas_config.radius),
            base_track_ratio=float(canvas_config.track_ratio),
        )
    feature_width_override_requested = feature_track_ratio_factor_override is not None

    precomputed_feature_dict: dict | None = None
    should_precompute_feature_dict = show_features and (
        show_external_labels
        or feature_track_ratio_factor_override is not None
        or core_track_overlap_relayout_enabled
        or bool(feature_slot is not None)
    )
    if should_precompute_feature_dict:
        compute_label_text = show_external_labels
        label_filtering = (
            preprocess_label_filtering(cfg.labels.filtering.as_dict())
            if compute_label_text
            else {}
        )
        color_table, default_colors = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
        precomputed_feature_dict, _ = create_feature_dict(
            gb_record,
            color_table,
            feature_config.selected_features_set,
            default_colors,
            cfg.canvas.strandedness,
            cfg.canvas.resolve_overlaps,
            label_filtering,
            split_overlaps_by_strand=split_overlaps_by_strand,
            directional_feature_types=feature_config.directional_feature_types,
            feature_visibility_rules=feature_config.feature_visibility_rules,
            compute_label_text=compute_label_text,
        )

    tick_label_annulus_for_legend_bounds: tuple[float, float] | None = None
    resolved_track_slots: list[CircularResolvedSlot] = []
    resolved_feature_anchor_radius_px: float | None = None
    definition_reserved_radius_px: float | None = None
    if center_reserved_radius is not None:
        definition_reserved_radius_px = max(0.0, float(center_reserved_radius))
    elif str(definition_position).strip().lower() == "center":
        definition_reserved_radius_px = _definition_reserved_radius_px(
            gb_record,
            canvas_config,
            species,
            strain,
            config_dict,
            cfg=cfg,
            plot_title=plot_title,
            definition_profile=definition_profile,
        )
    radial_layout = resolve_circular_radial_layout(
        total_length=len(gb_record.seq),
        canvas_config=canvas_config,
        cfg=cfg,
        slots=layout_slots,
        feature_dict=precomputed_feature_dict,
        show_features=show_features,
        show_ticks=show_ticks_track,
        definition_reserved_radius_px=definition_reserved_radius_px,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        tick_track_channel_override=_tick_track_channel_override,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
        depth_config=depth_config if show_depth_track else None,
    )
    setattr(canvas_config, "circular_radial_layout", radial_layout)
    setattr(canvas_config, "circular_feature_layout", radial_layout.features)
    resolved_track_slots = _resolved_slots_from_radial_layout(radial_layout, layout_slots)
    gc_content_tick_font_size_override = _gc_content_matches_depth_axis_font_size(
        gc_config=gc_config,
        depth_config=depth_config if show_depth_track else None,
        resolved_track_slots=resolved_track_slots,
    )
    if radial_layout.outer_content_radius_px > float(canvas_config.radius):
        if _expand_canvas_to_fit_radius(canvas_config, radial_layout.outer_content_radius_px):
            _sync_canvas_viewbox(canvas, canvas_config)

    rendered_feature_band_all_tracks: tuple[float, float] | None = None
    if radial_layout.features is not None:
        rendered_feature_band_all_tracks = (
            float(radial_layout.features.all_band_px.inner_px),
            float(radial_layout.features.all_band_px.outer_px),
        )
        resolved_feature_anchor_radius_px = float(radial_layout.features.anchor_radius_px)
        feature_track_ratio_factor_override = (
            float(radial_layout.features.width_px)
            / max(FEATURE_BAND_EPSILON, float(canvas_config.radius) * float(canvas_config.track_ratio))
        )
    for resolved_slot in resolved_track_slots:
        if resolved_slot.renderer != "ticks":
            continue
        tick_label_annulus = (
            float(resolved_slot.reserved_inner_radius_px),
            float(resolved_slot.reserved_outer_radius_px),
        )
        if tick_label_annulus_for_legend_bounds is None:
            tick_label_annulus_for_legend_bounds = tick_label_annulus
        else:
            tick_label_annulus_for_legend_bounds = (
                min(float(tick_label_annulus_for_legend_bounds[0]), float(tick_label_annulus[0])),
                max(float(tick_label_annulus_for_legend_bounds[1]), float(tick_label_annulus[1])),
            )

    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features.
    outer_arena: tuple[float, float] | None = None
    if show_external_labels:
        default_anchor_px, default_arc_outer_px = _default_outer_label_arena(
            canvas_config=canvas_config,
            cfg=cfg,
        )
        if default_arc_outer_px < default_anchor_px:
            default_anchor_px, default_arc_outer_px = default_arc_outer_px, default_anchor_px

        if feature_width_override_requested:
            outer_arena = (default_anchor_px, default_arc_outer_px)

    precalculated_labels: list[dict[str, Any]] | None = None
    if show_external_labels and precomputed_feature_dict is not None:
        precalculated_labels = prepare_label_list(
            precomputed_feature_dict,
            len(gb_record.seq),
            (
                float(resolved_feature_anchor_radius_px)
                if resolved_feature_anchor_radius_px is not None
                else canvas_config.radius
            ),
            canvas_config.track_ratio,
            config_dict,
            cfg=cfg,
            outer_arena=outer_arena,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            feature_layout=radial_layout.features,
            track_preset=_circular_preset_for_layout(canvas_config, cfg),
            feature_lane_direction=feature_lane_direction,
        )
        _resolve_label_legend_collisions(
            precalculated_labels,
            len(gb_record.seq),
            canvas_config,
            legend_config,
        )
        assign_leader_start_points(
            [label for label in precalculated_labels if not label.get("is_embedded")],
            len(gb_record.seq),
        )
        _expand_canvas_to_fit_external_labels(
            precalculated_labels,
            len(gb_record.seq),
            canvas_config,
        )
        _sync_canvas_viewbox(canvas, canvas_config)

    canvas = add_axis_group_on_canvas(
        canvas,
        canvas_config,
        config_dict,
        radius_override=None,
        cfg=cfg,
    )

    if show_external_labels:
        labels_group_kwargs: dict[str, Any] = {
            "outer_arena": outer_arena,
            "cfg": cfg,
            "precomputed_feature_dict": precomputed_feature_dict,
            "precalculated_labels": precalculated_labels,
            "feature_track_ratio_factor_override": feature_track_ratio_factor_override,
        }
        if resolved_feature_anchor_radius_px is not None and (
            precalculated_labels is None
            and (
                user_slot_mode
                or not math.isclose(
                    float(resolved_feature_anchor_radius_px),
                    float(canvas_config.radius),
                    rel_tol=1e-9,
                    abs_tol=1e-9,
                )
            )
        ):
            labels_group_kwargs["feature_anchor_radius_px"] = resolved_feature_anchor_radius_px
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            phase="leaders",
            **labels_group_kwargs,
        )

    drawable_slots = [
        slot for slot in resolved_track_slots
        if str(slot.renderer) != "feature_labels"
    ]
    for resolved_slot in sorted(drawable_slots, key=lambda item: (int(item.z), int(item.slot_index))):
        canvas = _draw_resolved_circular_slot(
            canvas,
            resolved_slot,
            gb_record=gb_record,
            canvas_config=canvas_config,
            feature_config=feature_config,
            config_dict=config_dict,
            gc_df=gc_df,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_df=depth_df,
            depth_config=depth_config,
            depth_tracks=depth_tracks,
            conservation_tracks=conservation_tracks,
            conservation_min_identity=conservation_min_identity,
            cfg=cfg,
            dinucleotide_dataframes=dinucleotide_dataframes,
            gc_content_tick_font_size_override=gc_content_tick_font_size_override,
            precomputed_feature_dict=precomputed_feature_dict,
            precalculated_labels=precalculated_labels,
            _tick_track_channel_override=_tick_track_channel_override,
            use_slot_group_id=user_slot_mode,
            use_slot_tick_options=user_slot_mode,
            use_feature_anchor_override=user_slot_mode,
        )

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            phase="text",
            **labels_group_kwargs,
        )

    definition_kwargs: dict[str, Any] = {"cfg": cfg}
    if plot_title is not None:
        definition_kwargs["plot_title"] = plot_title
    if str(definition_profile) != "full":
        definition_kwargs["definition_profile"] = definition_profile
    if str(definition_position) != "center":
        definition_kwargs["definition_position"] = definition_position
    if definition_group_id is not None:
        definition_kwargs["definition_group_id"] = definition_group_id
    canvas = add_record_definition_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        species,
        strain,
        config_dict,
        **definition_kwargs,
    )

    if canvas_config.legend_position != "none":
        if canvas_config.legend_position in {"top", "bottom"}:
            external_label_bounds: tuple[float, float, float, float] | None = None
            if precalculated_labels is not None:
                external_label_bounds = _external_label_bounds_on_canvas(
                    precalculated_labels,
                    len(gb_record.seq),
                    canvas_config,
                    margin_px=0.0,
                )
            content_top, content_bottom = _resolve_content_vertical_bounds_for_single_legend(
                total_length=len(gb_record.seq),
                canvas_config=canvas_config,
                cfg=cfg,
                show_features=show_features,
                precomputed_feature_dict=precomputed_feature_dict,
                rendered_feature_band_all_tracks=rendered_feature_band_all_tracks,
                feature_track_ratio_factor_override=feature_track_ratio_factor_override,
                tick_label_annulus=tick_label_annulus_for_legend_bounds,
                external_label_bounds=external_label_bounds,
            )
            _position_single_top_bottom_legend_between_edge_and_content(
                canvas,
                canvas_config,
                legend_config,
                position=canvas_config.legend_position,
                content_top=content_top,
                content_bottom=content_bottom,
            )
        canvas = add_legend_group_on_canvas(canvas, canvas_config, legend_config, legend_table)
    return canvas


def assemble_circular_diagram(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    feature_config: FeatureDrawingConfigurator,
    species: Optional[str],
    strain: Optional[str],
    plot_title: Optional[str],
    config_dict: dict,
    legend_config: LegendDrawingConfigurator,
    depth_df: DataFrame | None = None,
    depth_config: DepthConfigurator | None = None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    cfg: GbdrawConfig | None = None,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    _tick_track_channel_override: str | None = None,
) -> Drawing:
    """
    Assembles a circular diagram for a GenBank record and returns the SVG canvas.

    Parameters:
    gb_record (SeqRecord): The GenBank record to plot.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    species (str): Species name.
    strain (str): Strain name.
    config_dict (dict): Configuration dictionary for drawing parameters.
    out_formats (list): List of formats to save the output (e.g., ['png', 'svg']).

    Returns:
    Drawing: The assembled SVG canvas (not saved).
    """
    # Configure and create canvas

    # Prefer a pre-parsed config model when available to avoid repeated from_dict() calls.
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    effective_circular_track_slots = circular_track_slots

    legend_table: dict = {}
    if canvas_config.legend_position != "none":
        features_present = check_feature_presence(
            gb_record,
            feature_config.selected_features_set,
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )

        color_map, default_color_map = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            gb_record,
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        legend_table = prepare_legend_table(
            gc_config,
            skew_config,
            feature_config,
            features_present,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=depth_config if (
                depth_config is not None
                and (depth_df is not None or len(depth_tracks or ()) == 1)
            ) else None,
        )
        if depth_tracks and effective_circular_track_slots is None:
            legend_table = sync_depth_track_legend_entries(legend_table, depth_tracks)
        legend_table = _sync_legend_table_for_circular_slots(
            legend_table,
            circular_track_slots=effective_circular_track_slots,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_config=depth_config,
            depth_df=depth_df,
            depth_tracks=depth_tracks,
            cfg=cfg,
            conservation_tracks=conservation_tracks,
            conservation_min_identity=conservation_min_identity,
        )
        legend_config = legend_config.recalculate_legend_dimensions(legend_table, canvas_config)
        canvas_config.recalculate_canvas_dimensions(legend_config)
    canvas: Drawing = canvas_config.create_svg_canvas()
    canvas = add_record_on_circular_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        gc_config,
        skew_config,
        gc_df,
        species,
        strain,
        plot_title,
        config_dict,
        legend_config,
        legend_table,
        depth_config=depth_config,
        depth_df=depth_df,
        depth_tracks=depth_tracks,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=conservation_min_identity,
        cfg=cfg,
        circular_track_slots=effective_circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_dataframes=dinucleotide_dataframes,
        definition_position=definition_position,
        definition_profile=definition_profile,
        definition_group_id=definition_group_id,
        center_reserved_radius=center_reserved_radius,
        _tick_track_channel_override=_tick_track_channel_override,
    )
    return canvas


def plot_circular_diagram(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    feature_config: FeatureDrawingConfigurator,
    species: Optional[str],
    strain: Optional[str],
    plot_title: Optional[str],
    config_dict: dict,
    out_formats: list,
    legend_config: LegendDrawingConfigurator,
    depth_df: DataFrame | None = None,
    depth_config: DepthConfigurator | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    cfg: GbdrawConfig | None = None,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
) -> Drawing:
    """
    Backwards-compatible wrapper that assembles and saves a circular diagram.
    """
    canvas = assemble_circular_diagram(
        gb_record=gb_record,
        canvas_config=canvas_config,
        gc_df=gc_df,
        depth_df=depth_df,
        gc_config=gc_config,
        depth_config=depth_config,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=conservation_min_identity,
        skew_config=skew_config,
        feature_config=feature_config,
        species=species,
        strain=strain,
        plot_title=plot_title,
        config_dict=config_dict,
        legend_config=legend_config,
        cfg=cfg,
        circular_track_slots=circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_dataframes=dinucleotide_dataframes,
        definition_position=definition_position,
        definition_profile=definition_profile,
        definition_group_id=definition_group_id,
        center_reserved_radius=center_reserved_radius,
    )
    save_figure(canvas, out_formats)
    return canvas


__all__ = [
    "assemble_circular_diagram",
    "plot_circular_diagram",
    "add_record_on_circular_canvas",
]
