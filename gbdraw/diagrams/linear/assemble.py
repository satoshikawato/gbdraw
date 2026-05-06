#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from dataclasses import replace
import math

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.text import Text  # type: ignore[reportMissingImports]

from ...analysis.skew import skew_df  # type: ignore[reportMissingImports]
from ...analysis.depth import depth_df as build_depth_df  # type: ignore[reportMissingImports]
from ...analysis.protein_colinearity import OrthogroupResult  # type: ignore[reportMissingImports]
from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    LegendDrawingConfigurator,
)
from ...core.text import calculate_bbox_dimensions
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...render.groups.linear import LengthBarGroup, LegendGroup, PlotTitleGroup  # type: ignore[reportMissingImports]
from ...render.groups.linear.length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
)
from ...io.comparisons import filter_comparison_dataframe, load_comparisons
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]
from ...layout.linear import calculate_feature_position_factors_linear  # type: ignore[reportMissingImports]
from ...labels.linear import calculate_label_y_bounds  # type: ignore[reportMissingImports]
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...labels.linear import prepare_label_list_linear  # type: ignore[reportMissingImports]

from .builders import (
    add_comparison_on_linear_canvas,
    add_depth_group,
    add_gc_content_group,
    add_gc_skew_group,
    add_legends_on_linear_canvas,
    add_length_bar_on_linear_canvas,
    add_record_definition_group,
    add_record_group,
    add_row_comparison_on_linear_canvas,
)
from .orthogroup_alignment import (
    calculate_orthogroup_alignment_canvas_adjustment,
    calculate_orthogroup_alignment_offsets,
)
from .precalc import (
    FeatureDict,
    _precalculate_definition_metrics,
    _precalculate_feature_dicts,
    _precalculate_label_dimensions,
)
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...exceptions import ValidationError  # type: ignore[reportMissingImports]
from ...layout.linear_rows import (
    LinearInputRow,
    LinearLayoutIndex,
    LinearRecordLayout,
    build_linear_layout_index,
    clone_rows_for_rendering,
    flatten_linear_rows,
    segment_canvas_config,
    segment_cfg,
)


def _is_axis_ruler_enabled(canvas_config: LinearCanvasConfigurator, cfg: GbdrawConfig) -> bool:
    track_layout = str(canvas_config.track_layout).strip().lower()
    scale_style = str(cfg.objects.scale.style).strip().lower()
    return (
        bool(canvas_config.ruler_on_axis)
        and scale_style == "ruler"
        and track_layout in {"above", "below"}
    )


def _axis_ruler_extents(canvas_config: LinearCanvasConfigurator, cfg: GbdrawConfig) -> tuple[float, float]:
    """
    Return (height_above_axis, height_below_axis) required by axis-based ruler labels.
    """
    if not _is_axis_ruler_enabled(canvas_config, cfg):
        return 0.0, 0.0

    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param(canvas_config.length_param)
    label_height = calculate_bbox_dimensions(
        "0",
        cfg.objects.text.font_family,
        font_size,
        cfg.canvas.dpi,
    )[1]
    protrusion = max(
        0.5 * float(cfg.objects.scale.stroke_width),
        float(RULER_TICK_LENGTH),
        float(RULER_LABEL_OFFSET) + float(label_height),
    )
    track_layout = str(canvas_config.track_layout).strip().lower()
    if track_layout == "above":
        return 0.0, protrusion
    if track_layout == "below":
        return protrusion, 0.0
    return 0.0, 0.0


def _precalculate_feature_track_heights(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
    precomputed_feature_dicts: list[FeatureDict] | None = None,
) -> tuple[dict[str, float], dict[str, float], dict[str, float], dict[str, float], dict[str, float]]:
    """
    Pre-calculates the height required for feature tracks for each record.
    This is needed when resolve_overlaps is enabled as features may span multiple tracks.
    
    Returns:
        - dict mapping record_id -> height below the axis line (for lower tracks)
        - dict mapping record_id -> height above the axis line (for upper tracks)
        - dict mapping record_id -> minimum above-axis extent required to keep top visible
        - dict mapping record_id -> above-axis extent with non-displaced track positioning
        - dict mapping record_id -> middle-layout above-axis extent with non-displaced tracks
    """
    record_heights_below: dict[str, float] = {}
    record_heights_above: dict[str, float] = {}
    record_top_guard_above: dict[str, float] = {}
    record_top_guard_undisplaced: dict[str, float] = {}
    record_top_guard_middle_undisplaced: dict[str, float] = {}
    track_layout = str(canvas_config.track_layout).strip().lower()
    axis_gap_factor = (
        (float(canvas_config.track_axis_gap) / float(canvas_config.cds_height))
        if (canvas_config.track_axis_gap is not None and float(canvas_config.cds_height) > 0.0)
        else None
    )
    axis_ruler_above, axis_ruler_below = _axis_ruler_extents(canvas_config, cfg)
    
    if precomputed_feature_dicts is None:
        color_table, default_colors = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
    
    for index, record in enumerate(records):
        if precomputed_feature_dicts is not None:
            feature_dict = precomputed_feature_dicts[index]
        else:
            feature_dict, _ = create_feature_dict(
                record,
                color_table,
                feature_config.selected_features_set,
                default_colors,
                canvas_config.strandedness,
                canvas_config.resolve_overlaps,
                {},
                directional_feature_types=feature_config.directional_feature_types,
                feature_visibility_rules=feature_config.feature_visibility_rules,
                compute_label_text=False,
            )
        min_top_y = 0.0
        max_bottom_y = 0.0
        min_top_y_undisplaced = 0.0
        min_top_y_middle_undisplaced = 0.0
        for feature_obj in feature_dict.values():
            track_id = int(getattr(feature_obj, "feature_track_id", 0))
            strand = str(getattr(feature_obj, "strand", "undefined"))
            factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=track_id,
                separate_strands=canvas_config.strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            top_y = canvas_config.cds_height * float(factors[0])
            bottom_y = canvas_config.cds_height * float(factors[2])
            if top_y < min_top_y:
                min_top_y = top_y
            if bottom_y > max_bottom_y:
                max_bottom_y = bottom_y

            undisplaced_track_id = -1 if (canvas_config.strandedness and strand == "negative") else 0
            undisplaced_factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=undisplaced_track_id,
                separate_strands=canvas_config.strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            undisplaced_top_y = canvas_config.cds_height * float(undisplaced_factors[0])
            if undisplaced_top_y < min_top_y_undisplaced:
                min_top_y_undisplaced = undisplaced_top_y
            middle_undisplaced_factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=undisplaced_track_id,
                separate_strands=canvas_config.strandedness,
                track_layout="middle",
                axis_gap_factor=axis_gap_factor,
            )
            middle_undisplaced_top_y = canvas_config.cds_height * float(middle_undisplaced_factors[0])
            if middle_undisplaced_top_y < min_top_y_middle_undisplaced:
                min_top_y_middle_undisplaced = middle_undisplaced_top_y

        precise_height_above = max(0.0, -min_top_y)
        precise_height_below = max(0.0, max_bottom_y)
        undisplaced_height_above = max(0.0, -min_top_y_undisplaced)
        middle_undisplaced_height_above = max(0.0, -min_top_y_middle_undisplaced)
        
        if track_layout == "middle":
            # Keep existing middle-mode sizing behavior for backward compatibility.
            max_positive_track = 0
            min_negative_track = 0
            
            for feature_obj in feature_dict.values():
                track_id = feature_obj.feature_track_id
                if track_id > max_positive_track:
                    max_positive_track = track_id
                if track_id < min_negative_track:
                    min_negative_track = track_id
            
            if canvas_config.strandedness:
                # Stranded mode: positive tracks above axis, negative tracks below
                num_tracks_above = max_positive_track + 1
                num_tracks_below = abs(min_negative_track) + 1 if min_negative_track < 0 else 1
                height_above = num_tracks_above * canvas_config.cds_height * 1.1
                height_below = num_tracks_below * canvas_config.cds_height * 1.1
            else:
                # Non-stranded mode: track 0 is at axis, higher track IDs go below
                # Each track needs full cds_height of space plus some padding
                # Track 0 extends both above and below the axis (0.6 each direction)
                # Additional tracks (1, 2, ...) add cds_height * 1.1 below for breathing room
                height_above = canvas_config.cds_height * 0.6
                height_below = canvas_config.cds_height * (0.6 + max_positive_track * 1.1)
        else:
            # For above/below layouts, use actual positioned factors so downstream spacing,
            # GC/skew placement, and comparison ribbons align with feature extents.
            height_above = precise_height_above
            height_below = precise_height_below

        if axis_ruler_above > 0.0:
            height_above = max(height_above, axis_ruler_above)
            precise_height_above = max(precise_height_above, axis_ruler_above)
            undisplaced_height_above = max(undisplaced_height_above, axis_ruler_above)
            middle_undisplaced_height_above = max(middle_undisplaced_height_above, axis_ruler_above)
        if axis_ruler_below > 0.0:
            height_below = max(height_below, axis_ruler_below)
        
        record_heights_above[record.id] = height_above
        record_heights_below[record.id] = height_below
        record_top_guard_above[record.id] = precise_height_above
        record_top_guard_undisplaced[record.id] = undisplaced_height_above
        record_top_guard_middle_undisplaced[record.id] = middle_undisplaced_height_above

    return (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_undisplaced,
        record_top_guard_middle_undisplaced,
    )


def _precalculate_gc_dataframes(
    records: list[SeqRecord],
    *,
    window: int,
    step: int,
    dinucleotide: str,
    enabled: bool,
) -> list[DataFrame | None]:
    """Build GC/skew data once per record and share it across linear groups."""
    if not enabled:
        return [None for _ in records]
    return [skew_df(record, window, step, dinucleotide) for record in records]


def _precalculate_depth_dataframes(
    records: list[SeqRecord],
    *,
    depth_tables: list[DataFrame | None] | None,
    depth_config: DepthConfigurator | None,
    enabled: bool,
) -> list[DataFrame | None]:
    """Build depth data once per record and share it across linear groups."""
    if not enabled or depth_config is None or not depth_tables:
        return [None for _ in records]
    out: list[DataFrame | None] = []
    for index, record in enumerate(records):
        table = depth_tables[index] if index < len(depth_tables) else None
        if table is None:
            out.append(None)
            continue
        out.append(
            build_depth_df(
                record,
                table,
                int(depth_config.window),
                int(depth_config.step),
                normalize=bool(depth_config.normalize),
                min_depth=depth_config.min_depth,
                max_depth=depth_config.max_depth,
            )
        )
    return out


def _apply_shared_depth_axis(
    record_depth_dfs: list[DataFrame | None],
    depth_config: DepthConfigurator | None,
) -> None:
    """Use one automatically determined depth axis max across all records."""

    if depth_config is None or not bool(getattr(depth_config, "share_axis", False)):
        return
    if depth_config.max_depth is not None:
        return

    max_depth_values = [
        float(df["depth"].max())
        for df in record_depth_dfs
        if df is not None and not df.empty and "depth" in df.columns
    ]
    if max_depth_values:
        depth_config.max_depth = max(max_depth_values)


def _precalculate_label_heights_below(all_labels_by_record: dict[str, list[dict]]) -> dict[str, float]:
    """Return per-record label extents that protrude below the axis."""
    record_label_heights_below: dict[str, float] = {}
    for record_id, labels in all_labels_by_record.items():
        max_bottom_y = 0.0
        for label in labels:
            _, label_bottom_y = calculate_label_y_bounds(label)
            if bool(label.get("is_embedded")):
                feature_bottom_y = float(label.get("feature_bottom_y", 0.0))
                # Ignore labels that do not protrude outside their feature track.
                if label_bottom_y <= feature_bottom_y:
                    continue
            if label_bottom_y > max_bottom_y:
                max_bottom_y = label_bottom_y
        record_label_heights_below[record_id] = max_bottom_y
    return record_label_heights_below


def _row_label_metrics(
    rows: list[LinearInputRow],
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[float, dict[int, float], dict[int, float]]:
    font_size = cfg.objects.definition.linear.font_size.for_length_param(canvas_config.length_param)
    font_family = cfg.objects.text.font_family
    max_width = 0.0
    heights: dict[int, float] = {}
    half_heights: dict[int, float] = {}
    for row in rows:
        width, height = calculate_bbox_dimensions(
            str(row.label),
            font_family,
            font_size,
            cfg.canvas.dpi,
        )
        max_width = max(max_width, float(width))
        heights[row.row_index] = float(height)
        half_heights[row.row_index] = 0.5 * float(height)
    return max_width, heights, half_heights


def _add_row_definition_group(
    canvas: Drawing,
    row: LinearInputRow,
    offset_y: float,
    canvas_config: LinearCanvasConfigurator,
    max_def_width: float,
    cfg: GbdrawConfig,
    definition_column_shift_x: float = 0.0,
) -> None:
    def_cfg = cfg.objects.definition.linear
    font_size = def_cfg.font_size.for_length_param(canvas_config.length_param)
    definition_gap = min(canvas_config.canvas_padding, 20)
    definition_offset_x = (max_def_width / 2.0) + definition_gap + max(0.0, float(definition_column_shift_x))
    group = Group(id=f"{row.source_id}_definition")
    group.add(
        Text(
            str(row.label),
            insert=(0, 0),
            stroke=def_cfg.stroke,
            fill=def_cfg.fill,
            font_size=font_size,
            font_weight=def_cfg.font_weight,
            font_family=cfg.objects.text.font_family,
            text_anchor=def_cfg.text_anchor,
            dominant_baseline=def_cfg.dominant_baseline,
        )
    )
    group.translate(canvas_config.horizontal_offset - definition_offset_x, offset_y)
    canvas.add(group)


def _add_segment_label(
    canvas: Drawing,
    segment: LinearRecordLayout,
    offset_y: float,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> None:
    label = str(segment.ref.label or segment.ref.record_id).strip()
    if not label:
        return
    font_size = max(6.0, 0.75 * cfg.objects.definition.linear.font_size.for_length_param(canvas_config.length_param))
    label_gap = float(getattr(canvas_config, "record_label_gap", cfg.canvas.linear.record_label_gap))
    label_y = -label_gap
    dominant_baseline = "text-after-edge"
    if str(canvas_config.track_layout).strip().lower() == "above":
        label_y = label_gap
        dominant_baseline = "hanging"
    group = Group(id=f"{segment.ref.unique_key}_segment_label")
    group.add(
        Text(
            label,
            insert=(segment.x + (0.5 * segment.width), label_y),
            stroke="none",
            fill=cfg.objects.definition.linear.fill,
            font_size=font_size,
            font_weight=cfg.objects.definition.linear.font_weight,
            font_family=cfg.objects.text.font_family,
            text_anchor="middle",
            dominant_baseline=dominant_baseline,
        )
    )
    group.translate(canvas_config.horizontal_offset, offset_y)
    canvas.add(group)


def _segment_label_extents(
    rows: list[LinearInputRow],
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[dict[str, float], dict[str, float]]:
    font_size = max(6.0, 0.75 * cfg.objects.definition.linear.font_size.for_length_param(canvas_config.length_param))
    label_gap = float(getattr(canvas_config, "record_label_gap", cfg.canvas.linear.record_label_gap))
    track_layout = str(canvas_config.track_layout).strip().lower()
    above: dict[str, float] = {}
    below: dict[str, float] = {}
    for record in flatten_linear_rows(rows):
        label = str((getattr(record, "annotations", None) or {}).get("gbdraw_segment_label", record.id))
        _, label_height = calculate_bbox_dimensions(
            label,
            cfg.objects.text.font_family,
            font_size,
            cfg.canvas.dpi,
        )
        if track_layout == "above":
            above[record.id] = 0.0
            below[record.id] = label_gap + float(label_height)
        else:
            above[record.id] = label_gap + float(label_height)
            below[record.id] = 0.0
    return above, below


def _precalculate_row_segment_data(
    rows: list[LinearInputRow],
    layout_index: LinearLayoutIndex,
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig,
) -> tuple[list[FeatureDict], dict[str, list[dict]], dict[str, float], dict[str, float]]:
    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")
    color_table, default_colors = preprocess_color_tables(
        feature_config.color_table,
        feature_config.default_colors,
    )
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    segment_render_cfg = segment_cfg(cfg)

    feature_dicts: list[FeatureDict] = []
    all_labels_by_record: dict[str, list[dict]] = {}
    record_label_heights_above: dict[str, float] = {}
    record_label_heights_below: dict[str, float] = {}

    for row in rows:
        row_layout = layout_index.resolve_row(row.row_index)
        for record, segment in zip(row.records, row_layout.records):
            labels_enabled = (
                show_labels_mode == "all"
                or (show_labels_mode == "first" and row.row_index == 0)
            )
            feature_dict, _ = create_feature_dict(
                record,
                color_table,
                feature_config.selected_features_set,
                default_colors,
                canvas_config.strandedness,
                canvas_config.resolve_overlaps,
                label_filtering if labels_enabled else {},
                directional_feature_types=feature_config.directional_feature_types,
                feature_visibility_rules=feature_config.feature_visibility_rules,
                compute_label_text=labels_enabled,
            )
            feature_dicts.append(feature_dict)
            label_list: list[dict] = []
            if labels_enabled:
                segment_config = segment_canvas_config(canvas_config, segment)
                label_list = prepare_label_list_linear(
                    feature_dict,
                    len(record.seq),
                    segment_config.alignment_width,
                    1.0,
                    segment_config.cds_height,
                    segment_config.strandedness,
                    segment_config.track_layout,
                    segment_config.track_axis_gap,
                    config_dict,
                    cfg=segment_render_cfg,
                )
            all_labels_by_record[record.id] = label_list

            min_y_coord_record = 0.0
            max_bottom_y = 0.0
            for label in label_list:
                label_top_y, label_bottom_y = calculate_label_y_bounds(label)
                is_above_feature_embedded = (
                    bool(label.get("is_embedded"))
                    and label_top_y < float(label.get("feature_top_y", 0.0))
                )
                if bool(label.get("is_embedded")) and not is_above_feature_embedded:
                    pass
                elif label_top_y < min_y_coord_record:
                    min_y_coord_record = label_top_y

                if bool(label.get("is_embedded")):
                    feature_bottom_y = float(label.get("feature_bottom_y", 0.0))
                    if label_bottom_y <= feature_bottom_y:
                        continue
                if label_bottom_y > max_bottom_y:
                    max_bottom_y = label_bottom_y
            record_label_heights_above[record.id] = abs(min_y_coord_record)
            record_label_heights_below[record.id] = max_bottom_y
    segment_label_above, segment_label_below = _segment_label_extents(rows, canvas_config, cfg)
    for record_id, height in segment_label_above.items():
        record_label_heights_above[record_id] = max(record_label_heights_above.get(record_id, 0.0), height)
    for record_id, height in segment_label_below.items():
        record_label_heights_below[record_id] = max(record_label_heights_below.get(record_id, 0.0), height)
    return feature_dicts, all_labels_by_record, record_label_heights_above, record_label_heights_below


def assemble_linear_diagram_from_rows(
    rows: list[LinearInputRow],
    blast_files,
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    config_dict: dict,
    legend_config: LegendDrawingConfigurator,
    skew_config,
    depth_config: DepthConfigurator | None = None,
    depth_tables: list[DataFrame | None] | None = None,
    plot_title: str | None = None,
    plot_title_position: str = "bottom",
    plot_title_font_size: float = 32.0,
    comparison_dataframes: list[DataFrame] | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Assemble a row-aware linear diagram where one input source is one row."""

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    if not rows:
        raise ValidationError("rows is empty")
    render_rows, _id_map = clone_rows_for_rendering(rows)
    records = flatten_linear_rows(render_rows)
    if not records:
        raise ValidationError("rows contain no records")

    track_layout = str(canvas_config.track_layout).strip().lower()
    non_middle_layout = track_layout in {"above", "below"}
    axis_ruler_enabled = _is_axis_ruler_enabled(canvas_config, cfg)
    normalized_plot_title = str(plot_title or "").strip()
    normalized_plot_title_position = str(plot_title_position or "bottom").strip().lower()
    if normalized_plot_title_position not in {"center", "top", "bottom"}:
        raise ValueError("plot_title_position must be one of: center, top, bottom")

    plot_title_obj: PlotTitleGroup | None = None
    plot_title_edge_margin = 24.0
    plot_title_vertical_gap = float(canvas_config.vertical_padding)
    plot_title_top_reserve = 0.0
    plot_title_bottom_reserve = 0.0
    if normalized_plot_title:
        plot_title_obj = PlotTitleGroup(
            normalized_plot_title,
            config_dict,
            font_size=float(plot_title_font_size),
            cfg=cfg,
        )
        plot_title_text_height = max(float(plot_title_obj.text_bbox_height), float(plot_title_font_size))
        reserve = plot_title_edge_margin + plot_title_text_height + plot_title_vertical_gap
        if normalized_plot_title_position == "top":
            plot_title_top_reserve = reserve
        elif normalized_plot_title_position == "bottom":
            plot_title_bottom_reserve = reserve

    horizontal_layout = build_linear_layout_index(
        render_rows,
        alignment_width=canvas_config.alignment_width,
        record_gap=float(getattr(canvas_config, "record_gap", cfg.canvas.linear.record_gap)),
        normalize_length=bool(cfg.canvas.linear.normalize_length),
        align_center=bool(canvas_config.align_center),
    )
    segment_render_cfg = segment_cfg(cfg)
    record_feature_dicts, all_labels, record_label_heights_above, record_label_heights_below = _precalculate_row_segment_data(
        render_rows,
        horizontal_layout,
        feature_config,
        canvas_config,
        config_dict,
        cfg,
    )
    max_def_width, _definition_heights, definition_half_heights = _row_label_metrics(
        render_rows,
        canvas_config,
        cfg,
    )
    (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_undisplaced,
        record_top_guard_middle_undisplaced,
    ) = _precalculate_feature_track_heights(
        records,
        feature_config,
        canvas_config,
        cfg,
        precomputed_feature_dicts=record_feature_dicts,
    )

    row_record_ids: dict[int, list[str]] = {
        row.row_index: [record.id for record in row.records] for row in render_rows
    }

    def _row_max(row_index: int, values: dict[str, float], default: float = 0.0) -> float:
        return max((values.get(record_id, default) for record_id in row_record_ids[row_index]), default=default)

    row_feature_below = {row.row_index: _row_max(row.row_index, record_heights_below, canvas_config.cds_padding) for row in render_rows}
    row_feature_above = {row.row_index: _row_max(row.row_index, record_heights_above, canvas_config.cds_padding) for row in render_rows}
    row_label_above = {row.row_index: _row_max(row.row_index, record_label_heights_above, 0.0) for row in render_rows}
    row_label_below = {row.row_index: _row_max(row.row_index, record_label_heights_below, 0.0) for row in render_rows}
    row_top_guard_above = {row.row_index: _row_max(row.row_index, record_top_guard_above, canvas_config.cds_padding) for row in render_rows}
    row_top_guard_undisplaced = {row.row_index: _row_max(row.row_index, record_top_guard_undisplaced, canvas_config.cds_padding) for row in render_rows}
    row_top_guard_middle_undisplaced = {
        row.row_index: _row_max(row.row_index, record_top_guard_middle_undisplaced, canvas_config.cds_padding)
        for row in render_rows
    }

    required_label_height = max(row_label_above.values(), default=0.0)
    if required_label_height > 0:
        if canvas_config.vertical_offset < required_label_height:
            canvas_config.vertical_offset = (
                required_label_height + canvas_config.original_vertical_offset + canvas_config.cds_padding
            )
    else:
        canvas_config.vertical_offset = canvas_config.original_vertical_offset + canvas_config.cds_padding

    first_row_id = render_rows[0].row_index
    base_axis_y = canvas_config.vertical_offset
    first_actual_extent = row_top_guard_above.get(first_row_id, canvas_config.cds_padding)
    first_normal_extent = row_top_guard_undisplaced.get(first_row_id, first_actual_extent)
    normal_top_margin = max(canvas_config.vertical_padding, base_axis_y - first_normal_extent)
    if track_layout == "above":
        middle_floor_extent = row_top_guard_middle_undisplaced.get(first_row_id, first_actual_extent)
        middle_floor_margin = max(canvas_config.vertical_padding, base_axis_y - middle_floor_extent)
        normal_top_margin = max(normal_top_margin, middle_floor_margin)
    if track_layout == "above" or (canvas_config.resolve_overlaps and track_layout == "middle"):
        canvas_config.vertical_offset = max(base_axis_y, first_actual_extent + normal_top_margin)
    else:
        canvas_config.vertical_offset = max(
            base_axis_y,
            first_actual_extent + canvas_config.vertical_padding,
        )
    canvas_config.vertical_offset = max(
        canvas_config.vertical_offset,
        definition_half_heights.get(first_row_id, 0.0) + canvas_config.vertical_padding,
    )

    record_gc_dfs = _precalculate_gc_dataframes(
        records,
        window=int(gc_config.window),
        step=int(gc_config.step),
        dinucleotide=str(gc_config.dinucleotide),
        enabled=bool(canvas_config.show_gc or canvas_config.show_skew),
    )
    depth_enabled = bool(canvas_config.show_depth and depth_config is not None and depth_tables)
    record_depth_dfs = _precalculate_depth_dataframes(
        records,
        depth_tables=depth_tables,
        depth_config=depth_config,
        enabled=depth_enabled,
    )
    _apply_shared_depth_axis(record_depth_dfs, depth_config)
    if not depth_enabled:
        canvas_config.show_depth = False
        canvas_config.set_gc_height_and_gc_padding()

    has_blast = bool(blast_files or comparison_dataframes)
    comparisons: list[DataFrame] = []
    legend_table: dict = {}
    legend_group: LegendGroup | None = None
    required_legend_height = 0.0
    if canvas_config.legend_position != "none":
        features_present = check_feature_presence(
            records,
            feature_config.selected_features_set,
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        color_map, default_color_map = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            records,
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        legend_table = prepare_legend_table(
            gc_config, skew_config, feature_config, features_present, blast_config, has_blast,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=depth_config if depth_enabled else None,
        )
        legend_config = legend_config.recalculate_legend_dimensions(legend_table, canvas_config)
        legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table, cfg=cfg)
        required_legend_height = float(legend_group.legend_height)
        canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)
    else:
        canvas_config.alignment_width = canvas_config.fig_width
        canvas_config.horizontal_offset = 2 * canvas_config.canvas_padding + max_def_width
        canvas_config.total_width = (
            canvas_config.horizontal_offset
            + canvas_config.alignment_width
            + 2 * canvas_config.canvas_padding
        )
        canvas_config.legend_offset_x = 0
        canvas_config.legend_offset_y = 0

    vertical_shift = 0
    if canvas_config.legend_position not in ["top", "bottom"]:
        if required_legend_height > canvas_config.total_height:
            height_difference = required_legend_height - canvas_config.total_height
            canvas_config.total_height = int(required_legend_height)
            vertical_shift = height_difference / 2

    row_offsets: list[float] = []
    if canvas_config.legend_position == "top":
        current_y = canvas_config.original_vertical_offset + required_legend_height + canvas_config.vertical_offset
    else:
        current_y = canvas_config.vertical_offset if canvas_config.vertical_offset > vertical_shift else canvas_config.original_vertical_offset + vertical_shift
    current_y += plot_title_top_reserve

    for i, row in enumerate(render_rows):
        row_offsets.append(current_y)
        if i < len(render_rows) - 1:
            current_row_id = row.row_index
            next_row_id = render_rows[i + 1].row_index
            height_below_axis = (
                row_feature_below.get(current_row_id, canvas_config.cds_padding)
                + canvas_config.plot_tracks_height
                + row_label_below.get(current_row_id, 0.0)
            )
            height_above_next_axis = max(
                row_label_above.get(next_row_id, 0.0),
                row_feature_above.get(next_row_id, canvas_config.cds_padding),
            )
            min_gap = canvas_config.comparison_height if has_blast else canvas_config.cds_padding * 1.5
            inter_record_space = height_below_axis + min_gap + height_above_next_axis
            inter_record_space = max(
                inter_record_space,
                definition_half_heights.get(current_row_id, 0.0)
                + definition_half_heights.get(next_row_id, 0.0),
            )
            if definition_half_heights.get(current_row_id, 0.0) > 0.0 and definition_half_heights.get(next_row_id, 0.0) > 0.0:
                inter_record_space = max(
                    inter_record_space,
                    definition_half_heights.get(current_row_id, 0.0)
                    + definition_half_heights.get(next_row_id, 0.0)
                    + max(1.0, 0.5 * float(canvas_config.vertical_padding)),
                )
            current_y += inter_record_space

    layout_index = horizontal_layout.with_y_axes(row_offsets)

    length_bar_group: LengthBarGroup | None = None
    if not canvas_config.normalize_length and not axis_ruler_enabled:
        length_bar_group = LengthBarGroup(
            canvas_config.fig_width,
            canvas_config.alignment_width,
            canvas_config.longest_genome,
            config_dict,
            canvas_config,
            cfg=cfg,
        )
    length_bar_height = float(length_bar_group.scale_group_height) if length_bar_group is not None else 0.0

    final_row_id = render_rows[-1].row_index
    final_feature_height_below = (
        row_feature_below.get(final_row_id, canvas_config.cds_padding)
        if non_middle_layout
        else canvas_config.cds_padding
    )
    final_record_height_below = max(
        final_feature_height_below
        + row_label_below.get(final_row_id, 0.0)
        + canvas_config.plot_tracks_height,
        definition_half_heights.get(final_row_id, 0.0),
    )
    final_height = (
        current_y
        + final_record_height_below
        + length_bar_height
        + 4 * canvas_config.vertical_padding
        + canvas_config.original_vertical_offset
        + plot_title_bottom_reserve
    )
    canvas_config.height_below_final_record = (
        current_y
        + final_record_height_below
        + 4 * canvas_config.vertical_padding
    )
    if canvas_config.legend_position in ["top", "bottom"]:
        final_height += int(required_legend_height)
    canvas_config.total_height = max(final_height, canvas_config.total_height)

    if legend_group is not None:
        canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)

    if has_blast:
        comparison_sources: list[list[DataFrame]] = []
        if blast_files:
            comparison_sources.append(load_comparisons(blast_files, blast_config))
        if comparison_dataframes:
            comparison_sources.append(
                [
                    filter_comparison_dataframe(comparison, blast_config)
                    for comparison in comparison_dataframes
                ]
            )
        max_source_len = max((len(source) for source in comparison_sources), default=0)
        if max_source_len > max(0, len(render_rows) - 1):
            raise ValidationError(
                f"Too many comparison sources ({max_source_len}) for {len(render_rows)} source rows."
            )
        for index in range(max_source_len):
            frames = [
                source[index]
                for source in comparison_sources
                if index < len(source)
            ]
            if len(frames) == 1:
                comparisons.append(frames[0])
            elif frames:
                comparisons.append(pd.concat(frames, ignore_index=True))

    canvas: Drawing = canvas_config.create_svg_canvas()

    vertical_vb_width = canvas_config.total_width
    vertical_vb_height = canvas_config.total_height
    horizontal_vb_width = canvas_config.total_width
    horizontal_vb_height = canvas_config.total_height
    if legend_group is not None:
        h_legend_width, h_legend_height = legend_group.get_horizontal_dimensions()
        v_legend_width, v_legend_height = legend_group.get_vertical_dimensions()
        if canvas_config.legend_position in ["top", "bottom"]:
            vertical_vb_width = canvas_config.total_width - h_legend_width + v_legend_width
            vertical_vb_height = canvas_config.total_height - h_legend_height
        if canvas_config.legend_position in ["left", "right"]:
            horizontal_vb_width = canvas_config.total_width - v_legend_width
            horizontal_vb_height = canvas_config.total_height + h_legend_height
    canvas.attribs["data-vertical-viewbox"] = f"0 0 {vertical_vb_width} {vertical_vb_height}"
    canvas.attribs["data-horizontal-viewbox"] = f"0 0 {horizontal_vb_width} {horizontal_vb_height}"

    if canvas_config.legend_position == "top" and plot_title_top_reserve > 0:
        canvas_config.legend_offset_y += plot_title_top_reserve
    if canvas_config.legend_position != "none":
        canvas = add_legends_on_linear_canvas(canvas, config_dict, canvas_config, legend_group, legend_table)
    if length_bar_group is not None:
        canvas = add_length_bar_on_linear_canvas(canvas, canvas_config, config_dict, length_bar_group, legend_group)

    if has_blast:
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(render_rows) - 1):
            if non_middle_layout:
                ribbon_start_y = row_offsets[i]
                ribbon_end_y = row_offsets[i + 1]
            else:
                height_below_axis = canvas_config.cds_padding + canvas_config.plot_tracks_height
                ribbon_start_y = row_offsets[i] + height_below_axis
                ribbon_end_y = row_offsets[i + 1] - canvas_config.cds_padding
            comparison_offsets.append(ribbon_start_y)
            actual_comparison_heights.append(max(0.0, ribbon_end_y - ribbon_start_y))
        canvas = add_row_comparison_on_linear_canvas(
            canvas,
            comparisons,
            canvas_config,
            blast_config,
            layout_index,
            comparison_offsets,
            actual_comparison_heights,
        )

    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")
    cfg_labels_on = segment_cfg(replace(cfg, canvas=replace(cfg.canvas, show_labels=True)))
    cfg_labels_off = segment_cfg(replace(cfg, canvas=replace(cfg.canvas, show_labels=False)))

    flat_records = flatten_linear_rows(render_rows)
    flat_index_by_id = {record.id: index for index, record in enumerate(flat_records)}
    row_offset_by_index = {row.row_index: row_offsets[index] for index, row in enumerate(render_rows)}
    for row in render_rows:
        offset_y = row_offset_by_index[row.row_index]
        _add_row_definition_group(
            canvas,
            row,
            offset_y,
            canvas_config,
            max_def_width,
            cfg,
        )
        row_layout = layout_index.resolve_row(row.row_index)
        for record, segment in zip(row.records, row_layout.records):
            flat_index = flat_index_by_id[record.id]
            should_show_labels = False
            if show_labels_mode == "all":
                should_show_labels = True
            elif show_labels_mode == "first" and row.row_index == 0:
                should_show_labels = True
            record_cfg = cfg_labels_on if should_show_labels else cfg_labels_off
            segment_config = segment_canvas_config(canvas_config, segment)
            labels_for_record = all_labels.get(record.id)
            add_record_group(
                canvas,
                record,
                offset_y,
                segment.x,
                segment_config,
                feature_config,
                config_dict,
                precalculated_labels=labels_for_record,
                cfg=record_cfg,
                precomputed_feature_dict=record_feature_dicts[flat_index],
            )
            _add_segment_label(canvas, segment, offset_y, canvas_config, cfg)

            gc_offset_y = offset_y
            if non_middle_layout:
                current_feature_height_below = record_heights_below.get(record.id, canvas_config.cds_padding)
                gc_offset_y = offset_y + (current_feature_height_below - canvas_config.cds_padding)
            gc_offset_y += record_label_heights_below.get(record.id, 0.0)
            shared_gc_df = record_gc_dfs[flat_index] if flat_index < len(record_gc_dfs) else None
            shared_depth_df = record_depth_dfs[flat_index] if flat_index < len(record_depth_dfs) else None
            if depth_enabled and shared_depth_df is not None and depth_config is not None:
                add_depth_group(
                    canvas,
                    record,
                    gc_offset_y,
                    segment.x,
                    segment_config,
                    depth_config,
                    config_dict,
                    cfg=record_cfg,
                    depth_df=shared_depth_df,
                )
            if canvas_config.show_gc:
                add_gc_content_group(
                    canvas,
                    record,
                    gc_offset_y,
                    segment.x,
                    segment_config,
                    gc_config,
                    config_dict,
                    cfg=record_cfg,
                    gc_df=shared_gc_df,
                )
            if canvas_config.show_skew:
                add_gc_skew_group(
                    canvas,
                    record,
                    gc_offset_y,
                    segment.x,
                    segment_config,
                    skew_config,
                    config_dict,
                    cfg=record_cfg,
                    gc_df=shared_gc_df,
                )

    if plot_title_obj is not None:
        title_group = plot_title_obj.get_group()
        title_height = float(plot_title_obj.text_bbox_height)
        title_y = 0.5 * float(canvas_config.total_height)
        if normalized_plot_title_position == "top":
            title_y = plot_title_edge_margin + (0.5 * title_height)
        elif normalized_plot_title_position == "bottom":
            title_y = float(canvas_config.total_height) - plot_title_edge_margin - (0.5 * title_height)
        title_group.translate(0.5 * float(canvas_config.total_width), title_y)
        canvas.add(title_group)

    return canvas


def assemble_linear_diagram(
    records: list[SeqRecord],
    blast_files,
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    config_dict: dict,
    legend_config: LegendDrawingConfigurator,
    skew_config,
    depth_config: DepthConfigurator | None = None,
    depth_tables: list[DataFrame | None] | None = None,
    plot_title: str | None = None,
    plot_title_position: str = "bottom",
    plot_title_font_size: float = 32.0,
    comparison_dataframes: list[DataFrame] | None = None,
    orthogroups: OrthogroupResult | None = None,
    align_orthogroup_feature: str | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """
    Assembles a linear diagram of genomic records with optional BLAST comparison data,
    and returns the SVG canvas (not saved).
    """
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    track_layout = str(canvas_config.track_layout).strip().lower()
    non_middle_layout = track_layout in {"above", "below"}
    axis_ruler_enabled = _is_axis_ruler_enabled(canvas_config, cfg)
    normalized_plot_title = str(plot_title or "").strip()
    normalized_plot_title_position = str(plot_title_position or "bottom").strip().lower()
    if normalized_plot_title_position not in {"center", "top", "bottom"}:
        raise ValueError("plot_title_position must be one of: center, top, bottom")

    plot_title_obj: PlotTitleGroup | None = None
    plot_title_edge_margin = 24.0
    plot_title_vertical_gap = float(canvas_config.vertical_padding)
    plot_title_top_reserve = 0.0
    plot_title_bottom_reserve = 0.0
    if normalized_plot_title:
        plot_title_obj = PlotTitleGroup(
            normalized_plot_title,
            config_dict,
            font_size=float(plot_title_font_size),
            cfg=cfg,
        )
        plot_title_text_height = max(float(plot_title_obj.text_bbox_height), float(plot_title_font_size))
        reserve = plot_title_edge_margin + plot_title_text_height + plot_title_vertical_gap
        if normalized_plot_title_position == "top":
            plot_title_top_reserve = reserve
        elif normalized_plot_title_position == "bottom":
            plot_title_bottom_reserve = reserve

    record_feature_dicts = _precalculate_feature_dicts(
        records,
        feature_config,
        canvas_config,
        config_dict,
        cfg=cfg,
    )
    required_label_height, all_labels, record_label_heights_above = _precalculate_label_dimensions(
        records,
        feature_config,
        canvas_config,
        config_dict,
        cfg=cfg,
        precomputed_feature_dicts=record_feature_dicts,
    )
    record_label_heights_below = _precalculate_label_heights_below(all_labels)
    max_def_width, _definition_heights, definition_half_heights = _precalculate_definition_metrics(
        records,
        config_dict,
        canvas_config,
        cfg=cfg,
    )
    
    # Pre-calculate feature track heights for each record (needed for resolve_overlaps)
    (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_undisplaced,
        record_top_guard_middle_undisplaced,
    ) = _precalculate_feature_track_heights(
        records,
        feature_config,
        canvas_config,
        cfg,
        precomputed_feature_dicts=record_feature_dicts,
    )
    
    if required_label_height > 0:
        if canvas_config.vertical_offset < required_label_height:
            canvas_config.vertical_offset = (
                required_label_height + canvas_config.original_vertical_offset + canvas_config.cds_padding
            )
    else:
        canvas_config.vertical_offset = canvas_config.original_vertical_offset + canvas_config.cds_padding

    if records:
        first_record_id = records[0].id
        base_axis_y = canvas_config.vertical_offset
        first_actual_extent = record_top_guard_above.get(first_record_id, canvas_config.cds_padding)
        first_normal_extent = record_top_guard_undisplaced.get(first_record_id, first_actual_extent)
        normal_top_margin = max(canvas_config.vertical_padding, base_axis_y - first_normal_extent)
        if track_layout == "above":
            middle_floor_extent = record_top_guard_middle_undisplaced.get(first_record_id, first_actual_extent)
            middle_floor_margin = max(canvas_config.vertical_padding, base_axis_y - middle_floor_extent)
            normal_top_margin = max(normal_top_margin, middle_floor_margin)
        if track_layout == "above" or (canvas_config.resolve_overlaps and track_layout == "middle"):
            required_axis_y = first_actual_extent + normal_top_margin
            # In above layout, keep at least the middle-layout top-margin floor so
            # non-stranded tracks do not end up visually too close to the top edge.
            canvas_config.vertical_offset = max(base_axis_y, required_axis_y)
        else:
            minimum_axis_y_for_features = first_actual_extent + canvas_config.vertical_padding
            canvas_config.vertical_offset = max(base_axis_y, minimum_axis_y_for_features)
        canvas_config.vertical_offset = max(
            canvas_config.vertical_offset,
            definition_half_heights.get(first_record_id, 0.0) + canvas_config.vertical_padding,
        )

    normalize_length = cfg.canvas.linear.normalize_length
    record_gc_dfs = _precalculate_gc_dataframes(
        records,
        window=int(gc_config.window),
        step=int(gc_config.step),
        dinucleotide=str(gc_config.dinucleotide),
        enabled=bool(canvas_config.show_gc or canvas_config.show_skew),
    )
    depth_enabled = bool(canvas_config.show_depth and depth_config is not None and depth_tables)
    record_depth_dfs = _precalculate_depth_dataframes(
        records,
        depth_tables=depth_tables,
        depth_config=depth_config,
        enabled=depth_enabled,
    )
    _apply_shared_depth_axis(record_depth_dfs, depth_config)
    if not depth_enabled:
        canvas_config.show_depth = False
        canvas_config.set_gc_height_and_gc_padding()

    # Prepare legend group
    has_blast = bool(blast_files or comparison_dataframes)
    comparisons: list[DataFrame] = []
    legend_table: dict = {}
    legend_group: LegendGroup | None = None
    required_legend_height = 0.0
    if canvas_config.legend_position != "none":
        features_present = check_feature_presence(
            records,
            feature_config.selected_features_set,
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        color_map, default_color_map = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            records,
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        legend_table = prepare_legend_table(
            gc_config, skew_config, feature_config, features_present, blast_config, has_blast,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=depth_config if depth_enabled else None,
        )
        legend_config = legend_config.recalculate_legend_dimensions(legend_table, canvas_config)
        legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table, cfg=cfg)
        required_legend_height = float(legend_group.legend_height)
        canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)
    else:
        canvas_config.alignment_width = canvas_config.fig_width
        canvas_config.horizontal_offset = 2 * canvas_config.canvas_padding + max_def_width
        canvas_config.total_width = (
            canvas_config.horizontal_offset
            + canvas_config.alignment_width
            + 2 * canvas_config.canvas_padding
        )
        canvas_config.legend_offset_x = 0
        canvas_config.legend_offset_y = 0
    # Vertical shift: how much the records should be moved downward in order to place the records in the middle of the canvas
    vertical_shift = 0
    if canvas_config.legend_position in ["top", "bottom"]:
        pass  # If the legend is placed at the top or bottom of the canvas, no need to care about this
    else:
        # the height of the legend might be larger than that of the canvas if the legend is stacked vertically.
        if required_legend_height > canvas_config.total_height:
            height_difference = required_legend_height - canvas_config.total_height
            canvas_config.total_height = int(required_legend_height)
            vertical_shift = height_difference / 2

    record_ids = [r.id for r in records]
    record_offsets = []

    if canvas_config.legend_position == "top":
        current_y = canvas_config.original_vertical_offset + required_legend_height + canvas_config.vertical_offset
    else:
        if canvas_config.vertical_offset > vertical_shift:
            current_y = canvas_config.vertical_offset
        else:
            current_y = canvas_config.original_vertical_offset + vertical_shift
    current_y += plot_title_top_reserve

    for i, _ in enumerate(record_ids):
        record_offsets.append(current_y)

        if i < len(record_ids) - 1:
            current_record_id = record_ids[i]
            next_record_id = record_ids[i + 1]
            
            # Get the height below axis for the current record (feature tracks + GC/skew)
            current_feature_height_below = record_heights_below.get(current_record_id, canvas_config.cds_padding)
            height_below_axis = (
                current_feature_height_below
                + canvas_config.plot_tracks_height
            )
            current_label_height_below = record_label_heights_below.get(current_record_id, 0.0)
            height_below_axis += current_label_height_below
            
            # Get the height above axis for the next record (labels or feature tracks)
            next_label_height = record_label_heights_above.get(next_record_id, 0)
            next_feature_height_above = record_heights_above.get(next_record_id, canvas_config.cds_padding)
            # For above-axis height, we need to consider both labels and the upper part of features
            height_above_next_axis = max(next_label_height, next_feature_height_above)
            # For BLAST comparisons, use comparison_height as minimum space between records
            if has_blast:
                min_gap = canvas_config.comparison_height
            else:
                # Use cds_padding * 1.5 as minimum gap to ensure clean separation with some breathing room
                min_gap = canvas_config.cds_padding * 1.5
            
            # Total inter-record space: below current + gap + above next
            inter_record_space = height_below_axis + min_gap + height_above_next_axis
            inter_record_space = max(
                inter_record_space,
                definition_half_heights.get(current_record_id, 0.0)
                + definition_half_heights.get(next_record_id, 0.0),
            )
            if definition_half_heights.get(current_record_id, 0.0) > 0.0 and definition_half_heights.get(next_record_id, 0.0) > 0.0:
                inter_record_space = max(
                    inter_record_space,
                    definition_half_heights.get(current_record_id, 0.0)
                    + definition_half_heights.get(next_record_id, 0.0)
                    + max(1.0, 0.5 * float(canvas_config.vertical_padding)),
                )
            current_y += inter_record_space

    length_bar_group: LengthBarGroup | None = None
    if not canvas_config.normalize_length and not axis_ruler_enabled:
        length_bar_group = LengthBarGroup(
            canvas_config.fig_width,
            canvas_config.alignment_width,
            canvas_config.longest_genome,
            config_dict,
            canvas_config,
            cfg=cfg,
        )
    length_bar_height = (
        float(length_bar_group.scale_group_height)
        if length_bar_group is not None
        else 0.0
    )

    final_record_id = record_ids[-1] if record_ids else ""
    final_feature_height_below = (
        record_heights_below.get(final_record_id, canvas_config.cds_padding)
        if non_middle_layout
        else canvas_config.cds_padding
    )
    final_label_height_below = (
        record_label_heights_below.get(final_record_id, 0.0)
    )
    final_definition_height_below = definition_half_heights.get(final_record_id, 0.0)
    final_record_height_below = max(
        final_feature_height_below
        + final_label_height_below
        + canvas_config.plot_tracks_height,
        final_definition_height_below,
    )

    final_height = (
        current_y
        + final_record_height_below
        + length_bar_height
        + 4 * canvas_config.vertical_padding
        + canvas_config.original_vertical_offset
        + plot_title_bottom_reserve
    )
    canvas_config.height_below_final_record = (
        current_y
        + final_record_height_below
        + 4 * canvas_config.vertical_padding
    )
    if canvas_config.legend_position in ["top", "bottom"]:
        final_height += int(required_legend_height)
    canvas_config.total_height = max(final_height, canvas_config.total_height)

    if legend_group is not None:
        canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)

    if has_blast:
        comparison_sources: list[list[DataFrame]] = []
        if blast_files:
            comparison_sources.append(load_comparisons(blast_files, blast_config))
        if comparison_dataframes:
            comparison_sources.append(
                [
                    filter_comparison_dataframe(comparison, blast_config)
                    for comparison in comparison_dataframes
                ]
            )
        max_source_len = max((len(source) for source in comparison_sources), default=0)
        for index in range(max_source_len):
            frames = [
                source[index]
                for source in comparison_sources
                if index < len(source)
            ]
            if len(frames) == 1:
                comparisons.append(frames[0])
            elif frames:
                comparisons.append(pd.concat(frames, ignore_index=True))

    orthogroup_alignment_offsets = calculate_orthogroup_alignment_offsets(
        records,
        comparisons,
        canvas_config,
        align_orthogroup_feature,
        orthogroups=orthogroups,
    )
    alignment_shift_x, alignment_width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        orthogroup_alignment_offsets,
    )
    definition_column_shift_x = alignment_shift_x
    if alignment_shift_x or alignment_width_extension:
        width_extension_px = math.ceil(alignment_width_extension)
        canvas_config.horizontal_offset += alignment_shift_x
        canvas_config.total_width += width_extension_px
        if canvas_config.legend_position == "right":
            canvas_config.legend_offset_x += width_extension_px
        elif canvas_config.legend_position in {"top", "bottom"}:
            canvas_config.legend_offset_x += 0.5 * width_extension_px
    canvas: Drawing = canvas_config.create_svg_canvas()

    # Embed both viewBox configurations as data attributes for JavaScript repositioning
    # This allows switching between horizontal and vertical legend layouts without accumulation errors
    vertical_vb_width = canvas_config.total_width
    vertical_vb_height = canvas_config.total_height
    horizontal_vb_width = canvas_config.total_width
    horizontal_vb_height = canvas_config.total_height
    if legend_group is not None:
        h_legend_width, h_legend_height = legend_group.get_horizontal_dimensions()
        v_legend_width, v_legend_height = legend_group.get_vertical_dimensions()

        if canvas_config.legend_position in ["top", "bottom"]:
            vertical_vb_width = canvas_config.total_width - h_legend_width + v_legend_width
            vertical_vb_height = canvas_config.total_height - h_legend_height

        if canvas_config.legend_position in ["left", "right"]:
            horizontal_vb_width = canvas_config.total_width - v_legend_width
            horizontal_vb_height = canvas_config.total_height + h_legend_height

    canvas.attribs["data-vertical-viewbox"] = f"0 0 {vertical_vb_width} {vertical_vb_height}"
    canvas.attribs["data-horizontal-viewbox"] = f"0 0 {horizontal_vb_width} {horizontal_vb_height}"

    if canvas_config.legend_position == "top" and plot_title_top_reserve > 0:
        canvas_config.legend_offset_y += plot_title_top_reserve

    if canvas_config.legend_position != "none":
        canvas = add_legends_on_linear_canvas(canvas, config_dict, canvas_config, legend_group, legend_table)
    if length_bar_group is not None:
        canvas = add_length_bar_on_linear_canvas(canvas, canvas_config, config_dict, length_bar_group, legend_group)

    if has_blast:
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(records) - 1):
            if non_middle_layout:
                # In above/below layouts, anchor ribbons to consecutive record axes.
                ribbon_start_y = record_offsets[i]
                ribbon_end_y = record_offsets[i + 1]
            else:
                height_below_axis = (
                    canvas_config.cds_padding
                    + canvas_config.plot_tracks_height
                )
                ribbon_start_y = record_offsets[i] + height_below_axis
                ribbon_end_y = record_offsets[i + 1] - canvas_config.cds_padding
            comparison_offsets.append(ribbon_start_y)
            height = max(0.0, ribbon_end_y - ribbon_start_y)
            actual_comparison_heights.append(height)

        canvas = add_comparison_on_linear_canvas(
            canvas,
            comparisons,
            canvas_config,
            blast_config,
            config_dict,
            records,
            comparison_offsets,
            actual_comparison_heights,
            orthogroup_alignment_offsets,
        )

    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")
    cfg_labels_on = replace(cfg, canvas=replace(cfg.canvas, show_labels=True))
    cfg_labels_off = replace(cfg, canvas=replace(cfg.canvas, show_labels=False))

    for count, record in enumerate(records, start=1):
        offset_y = record_offsets[count - 1]

        if normalize_length:
            offset_x = 0
        else:
            offset_x = (
                (canvas_config.alignment_width * ((canvas_config.longest_genome - len(record.seq)) / canvas_config.longest_genome) / 2)
                if canvas_config.align_center
                else 0
            )
        offset_x += orthogroup_alignment_offsets.get(count - 1, 0.0)

        labels_for_record = all_labels.get(record.id)
        should_show_labels = False
        if show_labels_mode == "all":
            should_show_labels = True
        elif show_labels_mode == "first" and count == 1:
            should_show_labels = True
        record_cfg = cfg_labels_on if should_show_labels else cfg_labels_off

        add_record_group(
            canvas,
            record,
            offset_y,
            offset_x,
            canvas_config,
            feature_config,
            config_dict,
            precalculated_labels=labels_for_record,
            cfg=record_cfg,
            precomputed_feature_dict=record_feature_dicts[count - 1],
        )
        add_record_definition_group(
            canvas,
            record,
            offset_y,
            offset_x,
            canvas_config,
            config_dict,
            max_def_width,
            cfg=record_cfg,
            definition_column_shift_x=definition_column_shift_x,
        )
        gc_offset_y = offset_y
        if non_middle_layout:
            current_feature_height_below = record_heights_below.get(record.id, canvas_config.cds_padding)
            gc_offset_y = offset_y + (current_feature_height_below - canvas_config.cds_padding)
        gc_offset_y += record_label_heights_below.get(record.id, 0.0)
        shared_gc_df = record_gc_dfs[count - 1] if (count - 1) < len(record_gc_dfs) else None
        shared_depth_df = record_depth_dfs[count - 1] if (count - 1) < len(record_depth_dfs) else None

        if depth_enabled and shared_depth_df is not None and depth_config is not None:
            add_depth_group(
                canvas,
                record,
                gc_offset_y,
                offset_x,
                canvas_config,
                depth_config,
                config_dict,
                cfg=record_cfg,
                depth_df=shared_depth_df,
            )
        if canvas_config.show_gc:
            add_gc_content_group(
                canvas,
                record,
                gc_offset_y,
                offset_x,
                canvas_config,
                gc_config,
                config_dict,
                cfg=record_cfg,
                gc_df=shared_gc_df,
            )
        if canvas_config.show_skew:
            add_gc_skew_group(
                canvas,
                record,
                gc_offset_y,
                offset_x,
                canvas_config,
                skew_config,
                config_dict,
                cfg=record_cfg,
                gc_df=shared_gc_df,
            )

    if plot_title_obj is not None:
        title_group = plot_title_obj.get_group()
        title_height = float(plot_title_obj.text_bbox_height)
        title_y = 0.5 * float(canvas_config.total_height)
        if normalized_plot_title_position == "top":
            title_y = plot_title_edge_margin + (0.5 * title_height)
        elif normalized_plot_title_position == "bottom":
            title_y = float(canvas_config.total_height) - plot_title_edge_margin - (0.5 * title_height)
        title_group.translate(0.5 * float(canvas_config.total_width), title_y)
        canvas.add(title_group)

    return canvas


def plot_linear_diagram(
    records: list[SeqRecord],
    blast_files,
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    config_dict: dict,
    out_formats,
    legend_config,
    skew_config,
    depth_config: DepthConfigurator | None = None,
    depth_tables: list[DataFrame | None] | None = None,
    cfg: GbdrawConfig | None = None,
    comparison_dataframes: list[DataFrame] | None = None,
    orthogroups: OrthogroupResult | None = None,
    align_orthogroup_feature: str | None = None,
) -> Drawing:
    """Backwards-compatible wrapper that assembles and saves a linear diagram."""
    canvas = assemble_linear_diagram(
        records=records,
        blast_files=blast_files,
        canvas_config=canvas_config,
        blast_config=blast_config,
        feature_config=feature_config,
        gc_config=gc_config,
        config_dict=config_dict,
        legend_config=legend_config,
        skew_config=skew_config,
        depth_config=depth_config,
        depth_tables=depth_tables,
        comparison_dataframes=comparison_dataframes,
        orthogroups=orthogroups,
        align_orthogroup_feature=align_orthogroup_feature,
        cfg=cfg,
    )
    save_figure(canvas, out_formats)
    return canvas


__all__ = [
    "assemble_linear_diagram",
    "assemble_linear_diagram_from_rows",
    "plot_linear_diagram",
]
