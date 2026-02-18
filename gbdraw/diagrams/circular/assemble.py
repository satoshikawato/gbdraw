#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Optional, Callable

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...features.coordinates import get_strand  # type: ignore[reportMissingImports]
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
from ...layout.common import calculate_cds_ratio  # type: ignore[reportMissingImports]
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]
from ...svg.circular_ticks import get_circular_tick_path_ratio_bounds  # type: ignore[reportMissingImports]
from ...tracks import TrackSpec  # type: ignore[reportMissingImports]

from .builders import (
    add_axis_group_on_canvas,
    add_gc_content_group_on_canvas,
    add_gc_skew_group_on_canvas,
    add_labels_group_on_canvas,
    add_legend_group_on_canvas,
    add_record_definition_group_on_canvas,
    add_record_group_on_canvas,
    add_tick_group_on_canvas,
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


logger = logging.getLogger(__name__)


def _sync_canvas_viewbox(canvas: Drawing, canvas_config: CircularCanvasConfigurator) -> None:
    """Sync drawing viewport attrs with mutable canvas config values."""
    canvas.attribs["width"] = f"{canvas_config.total_width}px"
    canvas.attribs["height"] = f"{canvas_config.total_height}px"
    canvas.attribs["viewBox"] = f"0 0 {canvas_config.total_width} {canvas_config.total_height}"


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
        min_gap_px = minimum_bbox_gap_px(candidate, peer, base_margin_px=0.0)
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
                            min_gap_px = minimum_bbox_gap_px(candidate, peer, base_margin_px=0.0)
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


def _track_specs_by_kind(track_specs: list[TrackSpec] | None) -> dict[str, TrackSpec]:
    """Index TrackSpec list by kind (first occurrence wins)."""
    if not track_specs:
        return {}
    out: dict[str, TrackSpec] = {}
    for ts in track_specs:
        out.setdefault(str(ts.kind), ts)
    return out


def _resolve_circular_track_center_and_width_px(
    ts: TrackSpec | None, *, base_radius_px: float
) -> tuple[float | None, float | None]:
    """Resolve a circular placement into (center_radius_px, width_px).

    Supports:
    - radius + width
    - inner_radius + outer_radius
    """
    if ts is None or ts.placement is None:
        return None, None

    placement = ts.placement
    # Only handle circular placements (best-effort; keep permissive).
    if not hasattr(placement, "radius"):
        return None, None

    # ScalarSpec has .resolve(reference) method.
    inner_spec = getattr(placement, "inner_radius", None)
    outer_spec = getattr(placement, "outer_radius", None)
    radius_spec = getattr(placement, "radius", None)
    width_spec = getattr(placement, "width", None)

    inner_px = inner_spec.resolve(base_radius_px) if inner_spec is not None else None
    outer_px = outer_spec.resolve(base_radius_px) if outer_spec is not None else None
    radius_px = radius_spec.resolve(base_radius_px) if radius_spec is not None else None
    width_px = width_spec.resolve(base_radius_px) if width_spec is not None else None

    if inner_px is not None and outer_px is not None:
        if outer_px < inner_px:
            inner_px, outer_px = outer_px, inner_px
        return (inner_px + outer_px) / 2.0, (outer_px - inner_px)

    if radius_px is not None and width_px is not None:
        return radius_px, width_px

    # Partial specs: width-only keeps width and lets caller decide center.
    if width_px is not None:
        return None, width_px

    # Partial specs: treat "radius" alone as center; width unknown.
    if radius_px is not None:
        return radius_px, None

    return None, None


def _track_spec_has_explicit_center(ts: TrackSpec | None) -> bool:
    """Whether TrackSpec placement explicitly sets center radius."""
    if ts is None or ts.placement is None:
        return False
    placement = ts.placement
    if not hasattr(placement, "radius"):
        return False
    radius_spec = getattr(placement, "radius", None)
    inner_spec = getattr(placement, "inner_radius", None)
    outer_spec = getattr(placement, "outer_radius", None)
    return (radius_spec is not None) or (inner_spec is not None and outer_spec is not None)


def _resolve_feature_track_ratio_factor_override(
    ts: TrackSpec | None,
    *,
    base_radius_px: float,
    base_track_ratio: float,
) -> float | None:
    """Resolve feature width override into track_ratio_factor (width-only support)."""
    if ts is None or ts.placement is None:
        return None
    placement = ts.placement
    if not hasattr(placement, "radius"):
        return None

    radius_spec = getattr(placement, "radius", None)
    inner_spec = getattr(placement, "inner_radius", None)
    outer_spec = getattr(placement, "outer_radius", None)
    width_spec = getattr(placement, "width", None)

    if radius_spec is not None or inner_spec is not None or outer_spec is not None:
        logger.warning(
            "TrackSpec kind='features' supports width only. Center placement (r/ri/ro) is ignored."
        )

    if width_spec is None:
        return None

    width_px = float(width_spec.resolve(base_radius_px))
    if width_px <= 0:
        logger.warning("Ignoring non-positive features track width override: %.6f", width_px)
        return None
    if base_radius_px <= 0 or base_track_ratio <= 0:
        logger.warning(
            "Ignoring features track width override because base radius/track_ratio is non-positive."
        )
        return None
    return width_px / (base_radius_px * base_track_ratio)


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
) -> tuple[float, float] | None:
    """Compute (inner_radius_px, outer_radius_px) over all drawn feature blocks."""
    if not feature_dict or total_length <= 0:
        return None

    cds_ratio, offset = calculate_cds_ratio(track_ratio, str(length_param), track_ratio_factor)
    track_type = cfg.canvas.circular.track_type
    strandedness = cfg.canvas.strandedness

    min_inner: float | None = None
    max_outer: float | None = None

    for feature_object in feature_dict.values():
        track_id = int(getattr(feature_object, "feature_track_id", 0))
        if track_id_whitelist is not None and track_id not in track_id_whitelist:
            continue
        feature_location_list = list(getattr(feature_object, "location", []))
        list_of_coordinates = list(getattr(feature_object, "coordinates", []))

        for coordinate_idx, coordinate in enumerate(list_of_coordinates):
            if (
                coordinate_idx < len(feature_location_list)
                and getattr(feature_location_list[coordinate_idx], "kind", None) == "line"
            ):
                continue

            coord_start = int(coordinate.start)
            coord_end = int(coordinate.end)
            if coord_start == coord_end:
                continue

            coord_strand = get_strand(coordinate.strand)
            factors = calculate_feature_position_factors_circular(
                total_length,
                coord_strand,
                track_ratio,
                cds_ratio,
                offset,
                track_type,
                strandedness,
                track_id,
            )
            inner_px = base_radius_px * min(float(factors[0]), float(factors[2]))
            outer_px = base_radius_px * max(float(factors[0]), float(factors[2]))

            min_inner = inner_px if min_inner is None else min(min_inner, inner_px)
            max_outer = outer_px if max_outer is None else max(max_outer, outer_px)

    if min_inner is None or max_outer is None:
        return None
    return float(min_inner), float(max_outer)


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

    if kind == "gc_content":
        track_id = canvas_config.track_ids.get("gc_track")
    elif kind == "gc_skew":
        track_id = canvas_config.track_ids.get("skew_track")
    else:
        track_id = None

    if track_id is None:
        return None

    length_param = str(canvas_config.length_param)
    track_type = str(cfg.canvas.circular.track_type)
    norm_factor = float(cfg.canvas.circular.track_dict[length_param][track_type][str(track_id)])
    return float(canvas_config.radius) * norm_factor


def _default_outer_label_arena(
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[float, float]:
    """Return default (anchor_radius_px, arc_outer_radius_px) for outer labels."""
    length_param = str(canvas_config.length_param)
    track_type = str(cfg.canvas.circular.track_type)
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
    species: str,
    strain: str,
    config_dict: dict,
    legend_config,
    legend_table,
    *,
    cfg: GbdrawConfig | None = None,
    track_specs: list[TrackSpec] | None = None,
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
    ts_by_kind = _track_specs_by_kind(track_specs)

    raw_show_labels = cfg.canvas.show_labels
    show_labels_base = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)
    features_ts = ts_by_kind.get("features")
    labels_ts = ts_by_kind.get("labels")

    show_features = features_ts is None or features_ts.show
    show_external_labels = show_labels_base and (labels_ts is None or labels_ts.show) and show_features

    feature_track_ratio_factor_override: float | None = None
    if show_features:
        feature_track_ratio_factor_override = _resolve_feature_track_ratio_factor_override(
            features_ts,
            base_radius_px=float(canvas_config.radius),
            base_track_ratio=float(canvas_config.track_ratio),
        )

    precomputed_feature_dict: dict | None = None
    should_precompute_feature_dict = show_features and (
        show_external_labels or feature_track_ratio_factor_override is not None
    )
    if should_precompute_feature_dict:
        label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
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
        )

    feature_radius_mapper: Callable[[float], float] | None = None
    default_primary_feature_band: tuple[float, float] | None = None
    overridden_primary_feature_band: tuple[float, float] | None = None
    auto_relayout_active = False
    if (
        show_features
        and feature_track_ratio_factor_override is not None
        and precomputed_feature_dict is not None
    ):
        feature_radius_mapper, default_primary_feature_band, overridden_primary_feature_band = _build_feature_radius_mapper(
            precomputed_feature_dict,
            len(gb_record.seq),
            canvas_config=canvas_config,
            cfg=cfg,
            feature_track_ratio_factor_override=float(feature_track_ratio_factor_override),
        )
        auto_relayout_active = feature_radius_mapper is not None
        all_tracks_feature_band = _compute_feature_band_bounds_px(
            precomputed_feature_dict,
            len(gb_record.seq),
            base_radius_px=float(canvas_config.radius),
            track_ratio=float(canvas_config.track_ratio),
            length_param=str(canvas_config.length_param),
            track_ratio_factor=float(feature_track_ratio_factor_override),
            cfg=cfg,
        )
        if all_tracks_feature_band is not None:
            max_feature_radius = max(abs(float(all_tracks_feature_band[0])), abs(float(all_tracks_feature_band[1])))
            if _expand_canvas_to_fit_radius(canvas_config, max_feature_radius):
                _sync_canvas_viewbox(canvas, canvas_config)

    def _resolve_track_center_and_width_with_autorelayout(
        kind: str, ts: TrackSpec | None
    ) -> tuple[float | None, float | None]:
        center_px, width_px = _resolve_circular_track_center_and_width_px(
            ts,
            base_radius_px=float(canvas_config.radius),
        )
        has_explicit_center = _track_spec_has_explicit_center(ts)
        should_resolve_default_center = (
            center_px is None
            and (auto_relayout_active or (width_px is not None and not has_explicit_center))
        )
        if should_resolve_default_center:
            default_center_px = _default_track_center_radius_px(
                kind,
                canvas_config=canvas_config,
                cfg=cfg,
            )
            if default_center_px is not None:
                if auto_relayout_active and (not has_explicit_center) and feature_radius_mapper is not None:
                    center_px = float(feature_radius_mapper(default_center_px))
                else:
                    center_px = float(default_center_px)
        return center_px, width_px

    axis_ts = ts_by_kind.get("axis")
    if axis_ts is None or axis_ts.show:
        axis_radius_px, _ = _resolve_track_center_and_width_with_autorelayout("axis", axis_ts)
        canvas = add_axis_group_on_canvas(
            canvas,
            canvas_config,
            config_dict,
            radius_override=axis_radius_px,
            cfg=cfg,
        )

    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features.
    outer_arena: tuple[float, float] | None = None
    if show_external_labels:
        labels_center_px, labels_width_px = _resolve_circular_track_center_and_width_px(
            labels_ts,
            base_radius_px=float(canvas_config.radius),
        )
        labels_has_explicit_center = _track_spec_has_explicit_center(labels_ts)

        default_anchor_px, default_arc_outer_px = _default_outer_label_arena(
            canvas_config=canvas_config,
            cfg=cfg,
        )
        if auto_relayout_active and (not labels_has_explicit_center) and feature_radius_mapper is not None:
            default_anchor_px = float(feature_radius_mapper(default_anchor_px))
            default_arc_outer_px = float(feature_radius_mapper(default_arc_outer_px))
        if default_arc_outer_px < default_anchor_px:
            default_anchor_px, default_arc_outer_px = default_arc_outer_px, default_anchor_px

        if labels_center_px is None and labels_width_px is None:
            if auto_relayout_active and (not labels_has_explicit_center):
                outer_arena = (default_anchor_px, default_arc_outer_px)
        else:
            if labels_center_px is None:
                labels_center_px = float(default_anchor_px)
            if labels_width_px is None:
                labels_width_px = max(0.0, float(default_arc_outer_px - default_anchor_px))
            outer_arena = _arena_from_center_and_width(labels_center_px, labels_width_px)

    precalculated_labels: list[dict[str, Any]] | None = None
    if show_external_labels and precomputed_feature_dict is not None:
        precalculated_labels = prepare_label_list(
            precomputed_feature_dict,
            len(gb_record.seq),
            canvas_config.radius,
            canvas_config.track_ratio,
            config_dict,
            cfg=cfg,
            outer_arena=outer_arena,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
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

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            outer_arena=outer_arena,
            cfg=cfg,
            precomputed_feature_dict=precomputed_feature_dict,
            precalculated_labels=precalculated_labels,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        )

    if show_features:
        canvas = add_record_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            cfg=cfg,
            precomputed_feature_dict=precomputed_feature_dict,
            precalculated_labels=precalculated_labels,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        )

    definition_ts = ts_by_kind.get("definition")
    if definition_ts is None or definition_ts.show:
        canvas = add_record_definition_group_on_canvas(
            canvas, gb_record, canvas_config, species, strain, config_dict, cfg=cfg
        )

    ticks_ts = ts_by_kind.get("ticks")
    if ticks_ts is None or ticks_ts.show:
        ticks_radius_px, _ = _resolve_track_center_and_width_with_autorelayout("ticks", ticks_ts)
        ticks_has_explicit_center = _track_spec_has_explicit_center(ticks_ts)
        if (
            auto_relayout_active
            and show_features
            and (feature_track_ratio_factor_override is not None)
            and (default_primary_feature_band is not None)
            and (overridden_primary_feature_band is not None)
            and (not ticks_has_explicit_center)
        ):
            tick_ratio_bounds = get_circular_tick_path_ratio_bounds(
                len(gb_record.seq),
                str(cfg.canvas.circular.track_type),
                bool(cfg.canvas.strandedness),
            )
            effective_ticks_radius = (
                float(ticks_radius_px) if ticks_radius_px is not None else float(canvas_config.radius)
            )
            if _tick_annulus_overlaps_feature_band(
                effective_ticks_radius,
                tick_ratio_bounds,
                overridden_primary_feature_band,
            ):
                adjusted_ticks_radius = _resolve_ticks_center_radius_avoiding_feature_band(
                    current_ticks_radius_px=effective_ticks_radius,
                    default_ticks_radius_px=float(canvas_config.radius),
                    feature_band_px=overridden_primary_feature_band,
                    default_feature_band_px=default_primary_feature_band,
                    tick_ratio_bounds=tick_ratio_bounds,
                )
                if adjusted_ticks_radius is not None:
                    ticks_radius_px = float(adjusted_ticks_radius)
        canvas = add_tick_group_on_canvas(
            canvas, gb_record, canvas_config, config_dict, radius_override=ticks_radius_px, cfg=cfg
        )

    legend_ts = ts_by_kind.get("legend")
    if canvas_config.legend_position != "none" and (legend_ts is None or legend_ts.show):
        canvas = add_legend_group_on_canvas(canvas, canvas_config, legend_config, legend_table)

    # Add GC content group if configured to show.
    gc_ts = ts_by_kind.get("gc_content")
    if canvas_config.show_gc and (gc_ts is None or gc_ts.show):
        gc_center_px, gc_width_px = _resolve_track_center_and_width_with_autorelayout("gc_content", gc_ts)
        norm_factor_override = (gc_center_px / canvas_config.radius) if gc_center_px is not None else None
        canvas = add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            gc_df,
            canvas_config,
            gc_config,
            config_dict,
            track_width_override=gc_width_px,
            norm_factor_override=norm_factor_override,
            cfg=cfg,
        )

    skew_ts = ts_by_kind.get("gc_skew")
    if canvas_config.show_skew and (skew_ts is None or skew_ts.show):
        skew_center_px, skew_width_px = _resolve_track_center_and_width_with_autorelayout("gc_skew", skew_ts)
        norm_factor_override = (skew_center_px / canvas_config.radius) if skew_center_px is not None else None
        canvas = add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            gc_df,
            canvas_config,
            skew_config,
            config_dict,
            track_width_override=skew_width_px,
            norm_factor_override=norm_factor_override,
            cfg=cfg,
        )
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
    config_dict: dict,
    legend_config: LegendDrawingConfigurator,
    cfg: GbdrawConfig | None = None,
    track_specs: list[TrackSpec] | None = None,
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

    features_present = check_feature_presence(gb_record, feature_config.selected_features_set)

    # Pre-compute which color rules are actually used for accurate legend
    color_map, default_color_map = preprocess_color_tables(
        feature_config.color_table, feature_config.default_colors
    )
    used_color_rules, default_used_features = precompute_used_color_rules(
        gb_record, color_map, default_color_map, set(feature_config.selected_features_set)
    )
    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        features_present,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
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
        config_dict,
        legend_config,
        legend_table,
        cfg=cfg,
        track_specs=track_specs,
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
    config_dict: dict,
    out_formats: list,
    legend_config: LegendDrawingConfigurator,
    cfg: GbdrawConfig | None = None,
    track_specs: list[TrackSpec] | None = None,
) -> Drawing:
    """
    Backwards-compatible wrapper that assembles and saves a circular diagram.
    """
    canvas = assemble_circular_diagram(
        gb_record=gb_record,
        canvas_config=canvas_config,
        gc_df=gc_df,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        species=species,
        strain=strain,
        config_dict=config_dict,
        legend_config=legend_config,
        cfg=cfg,
        track_specs=track_specs,
    )
    save_figure(canvas, out_formats)
    return canvas


__all__ = [
    "assemble_circular_diagram",
    "plot_circular_diagram",
    "add_record_on_circular_canvas",
]
