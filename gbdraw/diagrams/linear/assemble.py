#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from dataclasses import replace

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    LegendDrawingConfigurator,
)
from ...core.text import calculate_bbox_dimensions
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...render.groups.linear import LengthBarGroup, LegendGroup  # type: ignore[reportMissingImports]
from ...render.groups.linear.length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
)
from ...io.comparisons import load_comparisons
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]
from ...layout.linear import calculate_feature_position_factors_linear  # type: ignore[reportMissingImports]
from ...labels.linear import calculate_label_y_bounds  # type: ignore[reportMissingImports]

from .builders import (
    add_comparison_on_linear_canvas,
    add_gc_content_group,
    add_gc_skew_group,
    add_legends_on_linear_canvas,
    add_length_bar_on_linear_canvas,
    add_record_definition_group,
    add_record_group,
)
from .precalc import _precalculate_definition_widths, _precalculate_label_dimensions
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]


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
) -> tuple[dict[str, float], dict[str, float], dict[str, float], dict[str, float]]:
    """
    Pre-calculates the height required for feature tracks for each record.
    This is needed when resolve_overlaps is enabled as features may span multiple tracks.
    
    Returns:
        - dict mapping record_id -> height below the axis line (for lower tracks)
        - dict mapping record_id -> height above the axis line (for upper tracks)
        - dict mapping record_id -> minimum above-axis extent required to keep top visible
        - dict mapping record_id -> above-axis extent with middle-layout positioning
    """
    record_heights_below: dict[str, float] = {}
    record_heights_above: dict[str, float] = {}
    record_top_guard_above: dict[str, float] = {}
    record_top_guard_middle_ref: dict[str, float] = {}
    track_layout = str(canvas_config.track_layout).strip().lower()
    axis_gap_factor = (
        (float(canvas_config.track_axis_gap) / float(canvas_config.cds_height))
        if (canvas_config.track_axis_gap is not None and float(canvas_config.cds_height) > 0.0)
        else None
    )
    axis_ruler_above, axis_ruler_below = _axis_ruler_extents(canvas_config, cfg)
    
    color_table, default_colors = preprocess_color_tables(
        feature_config.color_table, feature_config.default_colors
    )
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    
    for record in records:
        feature_dict, _ = create_feature_dict(
            record,
            color_table,
            feature_config.selected_features_set,
            default_colors,
            canvas_config.strandedness,
            canvas_config.resolve_overlaps,
            label_filtering,
            directional_feature_types=feature_config.directional_feature_types,
        )
        min_top_y = 0.0
        max_bottom_y = 0.0
        min_top_y_middle_ref = 0.0
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
            if track_layout == "above":
                middle_ref_factors = calculate_feature_position_factors_linear(
                    strand=strand,
                    track_id=track_id,
                    separate_strands=canvas_config.strandedness,
                    track_layout="middle",
                    axis_gap_factor=axis_gap_factor,
                )
                middle_ref_top_y = canvas_config.cds_height * float(middle_ref_factors[0])
                if middle_ref_top_y < min_top_y_middle_ref:
                    min_top_y_middle_ref = middle_ref_top_y

        precise_height_above = max(0.0, -min_top_y)
        precise_height_below = max(0.0, max_bottom_y)
        middle_ref_height_above = (
            max(0.0, -min_top_y_middle_ref) if track_layout == "above" else precise_height_above
        )
        
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
            middle_ref_height_above = max(middle_ref_height_above, axis_ruler_above)
        if axis_ruler_below > 0.0:
            height_below = max(height_below, axis_ruler_below)
        
        record_heights_above[record.id] = height_above
        record_heights_below[record.id] = height_below
        record_top_guard_above[record.id] = precise_height_above
        record_top_guard_middle_ref[record.id] = middle_ref_height_above

    return record_heights_below, record_heights_above, record_top_guard_above, record_top_guard_middle_ref


def _precalculate_label_heights_below(all_labels_by_record: dict[str, list[dict]]) -> dict[str, float]:
    """Return per-record label extents below the axis for spacing in below layout."""
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

    required_label_height, all_labels, record_label_heights_above = _precalculate_label_dimensions(
        records, feature_config, canvas_config, config_dict, cfg=cfg
    )
    record_label_heights_below = (
        _precalculate_label_heights_below(all_labels)
        if track_layout == "below"
        else {}
    )
    
    # Pre-calculate feature track heights for each record (needed for resolve_overlaps)
    (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_middle_ref,
    ) = _precalculate_feature_track_heights(records, feature_config, canvas_config, cfg)
    
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
        first_above_extent = record_top_guard_above.get(first_record_id, canvas_config.cds_padding)
        if track_layout == "above":
            first_middle_extent = record_top_guard_middle_ref.get(first_record_id, first_above_extent)
            middle_axis_y = max(base_axis_y, first_middle_extent + canvas_config.vertical_padding)
            target_top_margin = middle_axis_y - first_middle_extent
            required_axis_y = first_above_extent + target_top_margin
            # Match above-layout top margin to middle-layout feature-top margin for the first record.
            canvas_config.vertical_offset = max(base_axis_y, required_axis_y)
        else:
            minimum_axis_y_for_features = first_above_extent + canvas_config.vertical_padding
            canvas_config.vertical_offset = max(base_axis_y, minimum_axis_y_for_features)

    normalize_length = cfg.canvas.linear.normalize_length

    # Prepare legend group
    has_blast = bool(blast_files)
    # Determine which features should be displayed in the legend
    features_present = check_feature_presence(records, feature_config.selected_features_set)
    # Pre-compute which color rules are actually used for accurate legend
    color_map, default_color_map = preprocess_color_tables(
        feature_config.color_table, feature_config.default_colors
    )
    used_color_rules, default_used_features = precompute_used_color_rules(
        records, color_map, default_color_map, set(feature_config.selected_features_set)
    )
    # Prepare legend table
    legend_table = prepare_legend_table(
        gc_config, skew_config, feature_config, features_present, blast_config, has_blast,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
    )
    # Predetermine legend dimensions (number of columns etc.)
    legend_config = legend_config.recalculate_legend_dimensions(legend_table, canvas_config)
    # Draw legend group to determine the actual dimensions
    legend_group: Group = LegendGroup(config_dict, canvas_config, legend_config, legend_table, cfg=cfg)
    # Get the legend height
    required_legend_height = legend_group.legend_height
    max_def_width = _precalculate_definition_widths(records, config_dict, canvas_config, cfg=cfg)

    canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)
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
        current_y = canvas_config.original_vertical_offset + legend_group.legend_height + canvas_config.vertical_offset
    else:
        if canvas_config.vertical_offset > vertical_shift:
            current_y = canvas_config.vertical_offset
        else:
            current_y = canvas_config.original_vertical_offset + vertical_shift

    for i, _ in enumerate(record_ids):
        record_offsets.append(current_y)

        if i < len(record_ids) - 1:
            current_record_id = record_ids[i]
            next_record_id = record_ids[i + 1]
            
            # Get the height below axis for the current record (feature tracks + GC/skew)
            current_feature_height_below = record_heights_below.get(current_record_id, canvas_config.cds_padding)
            height_below_axis = current_feature_height_below + canvas_config.gc_padding + canvas_config.skew_padding
            
            # Get the height above axis for the next record (labels or feature tracks)
            next_label_height = record_label_heights_above.get(next_record_id, 0)
            next_feature_height_above = record_heights_above.get(next_record_id, canvas_config.cds_padding)
            # For above-axis height, we need to consider both labels and the upper part of features
            height_above_next_axis = max(next_label_height, next_feature_height_above)
            if track_layout == "below":
                current_label_height_below = record_label_heights_below.get(current_record_id, 0.0)
                height_below_axis += current_label_height_below
            
            # For BLAST comparisons, use comparison_height as minimum space between records
            if has_blast:
                min_gap = canvas_config.comparison_height
            else:
                # Use cds_padding * 1.5 as minimum gap to ensure clean separation with some breathing room
                min_gap = canvas_config.cds_padding * 1.5
            
            # Total inter-record space: below current + gap + above next
            inter_record_space = height_below_axis + min_gap + height_above_next_axis
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
        if track_layout == "below"
        else 0.0
    )

    final_height = (
        current_y
        + final_feature_height_below
        + final_label_height_below
        + canvas_config.gc_padding
        + canvas_config.skew_padding
        + length_bar_height
        + 4 * canvas_config.vertical_padding
        + canvas_config.original_vertical_offset
    )
    canvas_config.height_below_final_record = (
        current_y
        + final_feature_height_below
        + final_label_height_below
        + canvas_config.gc_padding
        + canvas_config.skew_padding
        + 4 * canvas_config.vertical_padding
    )
    if canvas_config.legend_position in ["top", "bottom"]:
        final_height += int(required_legend_height)
    canvas_config.total_height = max(final_height, canvas_config.total_height)

    canvas_config.recalculate_canvas_dimensions(legend_group, max_def_width)
    canvas: Drawing = canvas_config.create_svg_canvas()

    # Embed both viewBox configurations as data attributes for JavaScript repositioning
    # This allows switching between horizontal and vertical legend layouts without accumulation errors
    h_legend_width, h_legend_height = legend_group.get_horizontal_dimensions()
    v_legend_width, v_legend_height = legend_group.get_vertical_dimensions()

    # Calculate viewBox for vertical layout (left/right legend positions)
    vertical_vb_width = canvas_config.total_width
    vertical_vb_height = canvas_config.total_height
    if canvas_config.legend_position in ["top", "bottom"]:
        # Current dimensions are for horizontal layout, calculate vertical
        vertical_vb_width = canvas_config.total_width - h_legend_width + v_legend_width
        vertical_vb_height = canvas_config.total_height - h_legend_height

    # Calculate viewBox for horizontal layout (top/bottom legend positions)
    horizontal_vb_width = canvas_config.total_width
    horizontal_vb_height = canvas_config.total_height
    if canvas_config.legend_position in ["left", "right"]:
        # Current dimensions are for vertical layout, calculate horizontal
        horizontal_vb_width = canvas_config.total_width - v_legend_width
        horizontal_vb_height = canvas_config.total_height + h_legend_height

    canvas.attribs["data-vertical-viewbox"] = f"0 0 {vertical_vb_width} {vertical_vb_height}"
    canvas.attribs["data-horizontal-viewbox"] = f"0 0 {horizontal_vb_width} {horizontal_vb_height}"

    if canvas_config.legend_position != "none":
        canvas = add_legends_on_linear_canvas(canvas, config_dict, canvas_config, legend_group, legend_table)
    if length_bar_group is not None:
        canvas = add_length_bar_on_linear_canvas(canvas, canvas_config, config_dict, length_bar_group, legend_group)

    if blast_files:
        comparisons = load_comparisons(blast_files, blast_config)
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(records) - 1):
            if non_middle_layout:
                # In above/below layouts, anchor ribbons to consecutive record axes.
                ribbon_start_y = record_offsets[i]
                ribbon_end_y = record_offsets[i + 1]
            else:
                height_below_axis = canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding
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
        )
        gc_offset_y = offset_y
        if non_middle_layout:
            current_feature_height_below = record_heights_below.get(record.id, canvas_config.cds_padding)
            gc_offset_y = offset_y + (current_feature_height_below - canvas_config.cds_padding)
            if track_layout == "below":
                gc_offset_y += record_label_heights_below.get(record.id, 0.0)

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
            )

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
    cfg: GbdrawConfig | None = None,
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
        cfg=cfg,
    )
    save_figure(canvas, out_formats)
    return canvas


__all__ = [
    "assemble_linear_diagram",
    "plot_linear_diagram",
]


