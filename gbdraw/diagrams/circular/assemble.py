#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import math
from typing import Any, Optional

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
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...labels.circular import minimum_bbox_gap_px, prepare_label_list, x_overlap, y_overlap  # type: ignore[reportMissingImports]
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]
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


def _candidate_shift_directions(
    label: dict[str, Any],
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
) -> list[tuple[float, float]]:
    """Return normalized shift directions to test for one colliding label."""
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    directions: list[tuple[float, float]] = []

    radial_norm = math.hypot(start_x, start_y)
    if radial_norm > 1e-6:
        # First preference: move label inward toward the circular map.
        directions.append((-start_x / radial_norm, -start_y / radial_norm))

    legend_min_x, legend_min_y, legend_max_x, legend_max_y = _legend_bbox(canvas_config, legend_config)
    legend_center_x = 0.5 * (legend_min_x + legend_max_x)
    legend_center_y = 0.5 * (legend_min_y + legend_max_y)
    start_abs_x = float(canvas_config.offset_x) + start_x
    start_abs_y = float(canvas_config.offset_y) + start_y
    away_x = start_abs_x - legend_center_x
    away_y = start_abs_y - legend_center_y
    away_norm = math.hypot(away_x, away_y)
    if away_norm > 1e-6:
        directions.append((away_x / away_norm, away_y / away_norm))

    return directions


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
        for idx in collided_indices:
            label = labels[idx]
            if label.get("is_embedded"):
                continue
            start_x = float(label["start_x"])
            start_y = float(label["start_y"])
            x_sign = 1 if start_x >= 0 else -1

            placed = False
            shift_px = LABEL_NUDGE_STEP_PX
            directions = _candidate_shift_directions(label, canvas_config, legend_config)
            while shift_px <= MAX_LABEL_NUDGE_PX and not placed:
                for dir_x, dir_y in directions:
                    candidate = label.copy()
                    candidate["start_x"] = start_x + dir_x * shift_px
                    candidate["start_y"] = start_y + dir_y * shift_px
                    # Avoid anchor-side flips that can make labels jump visually.
                    if x_sign > 0 and float(candidate["start_x"]) <= 0:
                        continue
                    if x_sign < 0 and float(candidate["start_x"]) >= 0:
                        continue
                    candidate_box = _label_bbox_on_canvas(candidate, total_length, canvas_config, margin_px=LEGEND_LABEL_MARGIN_PX)
                    if _bbox_overlaps(candidate_box, _legend_bbox(canvas_config, legend_config)):
                        continue
                    if _label_overlaps_other_labels(candidate, labels, idx, total_length):
                        continue
                    label["start_x"] = float(candidate["start_x"])
                    label["start_y"] = float(candidate["start_y"])
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

    # Partial specs: treat "radius" alone as center; width unknown.
    if radius_px is not None:
        return radius_px, None

    return None, None


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

    axis_ts = ts_by_kind.get("axis")
    if axis_ts is None or axis_ts.show:
        axis_radius_px, _ = _resolve_circular_track_center_and_width_px(axis_ts, base_radius_px=canvas_config.radius)
        canvas = add_axis_group_on_canvas(canvas, canvas_config, config_dict, radius_override=axis_radius_px, cfg=cfg)

    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features
    labels_ts = ts_by_kind.get("labels")
    raw_show_labels = cfg.canvas.show_labels
    show_labels_base = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)
    features_ts = ts_by_kind.get("features")
    show_features = features_ts is None or features_ts.show
    show_external_labels = show_labels_base and (labels_ts is None or labels_ts.show) and show_features

    outer_arena = None
    if show_external_labels and labels_ts is not None:
        center_px, width_px = _resolve_circular_track_center_and_width_px(labels_ts, base_radius_px=canvas_config.radius)
        if center_px is not None and width_px is not None:
            inner_px = center_px - (width_px / 2.0)
            outer_px = center_px + (width_px / 2.0)
            if outer_px < inner_px:
                inner_px, outer_px = outer_px, inner_px
            outer_arena = (inner_px, outer_px)

    precomputed_feature_dict = None
    precalculated_labels = None
    if show_external_labels:
        label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
        color_table, default_colors = preprocess_color_tables(feature_config.color_table, feature_config.default_colors)
        precomputed_feature_dict, _ = create_feature_dict(
            gb_record,
            color_table,
            feature_config.selected_features_set,
            default_colors,
            cfg.canvas.strandedness,
            cfg.canvas.resolve_overlaps,
            label_filtering,
        )
        precalculated_labels = prepare_label_list(
            precomputed_feature_dict,
            len(gb_record.seq),
            canvas_config.radius,
            canvas_config.track_ratio,
            config_dict,
            cfg=cfg,
            outer_arena=outer_arena,
        )
        _resolve_label_legend_collisions(
            precalculated_labels,
            len(gb_record.seq),
            canvas_config,
            legend_config,
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
        )

    definition_ts = ts_by_kind.get("definition")
    if definition_ts is None or definition_ts.show:
        canvas = add_record_definition_group_on_canvas(
            canvas, gb_record, canvas_config, species, strain, config_dict, cfg=cfg
        )

    ticks_ts = ts_by_kind.get("ticks")
    if ticks_ts is None or ticks_ts.show:
        ticks_radius_px, _ = _resolve_circular_track_center_and_width_px(ticks_ts, base_radius_px=canvas_config.radius)
        canvas = add_tick_group_on_canvas(
            canvas, gb_record, canvas_config, config_dict, radius_override=ticks_radius_px, cfg=cfg
        )

    legend_ts = ts_by_kind.get("legend")
    if canvas_config.legend_position != "none" and (legend_ts is None or legend_ts.show):
        canvas = add_legend_group_on_canvas(canvas, canvas_config, legend_config, legend_table)

    # Add GC content group if configured to show
    gc_ts = ts_by_kind.get("gc_content")
    if canvas_config.show_gc and (gc_ts is None or gc_ts.show):
        gc_center_px, gc_width_px = _resolve_circular_track_center_and_width_px(gc_ts, base_radius_px=canvas_config.radius)
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
        skew_center_px, skew_width_px = _resolve_circular_track_center_and_width_px(skew_ts, base_radius_px=canvas_config.radius)
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
