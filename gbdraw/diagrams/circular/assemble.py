#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from typing import Optional

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
from ...labels.circular import prepare_label_list  # type: ignore[reportMissingImports]
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
