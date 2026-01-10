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
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...render.groups.linear import LengthBarGroup, LegendGroup  # type: ignore[reportMissingImports]
from ...io.comparisons import load_comparisons
from ...legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from ...render.export import save_figure  # type: ignore[reportMissingImports]

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


def _precalculate_feature_track_heights(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Pre-calculates the height required for feature tracks for each record.
    This is needed when resolve_overlaps is enabled as features may span multiple tracks.
    
    Returns:
        - dict mapping record_id -> height below the axis line (for lower tracks)
        - dict mapping record_id -> height above the axis line (for upper tracks)
    """
    record_heights_below: dict[str, float] = {}
    record_heights_above: dict[str, float] = {}
    
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
        )
        
        # Find the maximum and minimum track IDs
        max_positive_track = 0
        min_negative_track = 0
        
        for feature_obj in feature_dict.values():
            track_id = feature_obj.feature_track_id
            if track_id > max_positive_track:
                max_positive_track = track_id
            if track_id < min_negative_track:
                min_negative_track = track_id
        
        # Calculate heights based on number of tracks
        # Use cds_height as the base unit for track spacing
        # Add generous padding to prevent any overlap
        
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
        
        record_heights_above[record.id] = height_above
        record_heights_below[record.id] = height_below
    
    return record_heights_below, record_heights_above


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

    required_label_height, all_labels, record_label_heights = _precalculate_label_dimensions(
        records, feature_config, canvas_config, config_dict, cfg=cfg
    )
    
    # Pre-calculate feature track heights for each record (needed for resolve_overlaps)
    record_heights_below, record_heights_above = _precalculate_feature_track_heights(
        records, feature_config, canvas_config, cfg
    )
    
    if required_label_height > 0:
        if canvas_config.vertical_offset < required_label_height:
            canvas_config.vertical_offset = (
                required_label_height + canvas_config.original_vertical_offset + canvas_config.cds_padding
            )
    else:
        canvas_config.vertical_offset = canvas_config.original_vertical_offset + canvas_config.cds_padding

    normalize_length = cfg.canvas.linear.normalize_length

    # Prepare legend group
    has_blast = bool(blast_files)
    # Determine which features should be displayed in the legend
    features_present = check_feature_presence(records, feature_config.selected_features_set)
    # Pre-compute which color rules are actually used for accurate legend
    color_map, default_color_map = preprocess_color_tables(
        feature_config.color_table, feature_config.default_colors
    )
    used_color_rules = precompute_used_color_rules(
        records, color_map, default_color_map, set(feature_config.selected_features_set)
    )
    # Prepare legend table
    legend_table = prepare_legend_table(
        gc_config, skew_config, feature_config, features_present, blast_config, has_blast,
        used_color_rules=used_color_rules
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
            next_label_height = record_label_heights.get(next_record_id, 0)
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
            current_y += inter_record_space

    length_bar_group: Group = LengthBarGroup(
        canvas_config.fig_width,
        canvas_config.alignment_width,
        canvas_config.longest_genome,
        config_dict,
        canvas_config,
        cfg=cfg,
    )

    final_height = (
        current_y
        + canvas_config.cds_padding
        + canvas_config.gc_padding
        + canvas_config.skew_padding
        + length_bar_group.scale_group_height
        + 4 * canvas_config.vertical_padding
        + canvas_config.original_vertical_offset
    )
    canvas_config.height_below_final_record = (
        current_y
        + canvas_config.cds_padding
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
    if not canvas_config.normalize_length:
        canvas = add_length_bar_on_linear_canvas(canvas, canvas_config, config_dict, length_bar_group, legend_group)

    if blast_files:
        comparisons = load_comparisons(blast_files, blast_config)
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(records) - 1):
            height_below_axis = canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding
            ribbon_start_y = record_offsets[i] + height_below_axis
            comparison_offsets.append(ribbon_start_y)
            next_record_id = records[i + 1].id
            next_label_height = record_label_heights.get(next_record_id, 0)
            ribbon_end_y = record_offsets[i + 1] - canvas_config.cds_padding
            height = ribbon_end_y - ribbon_start_y
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
        if canvas_config.show_gc:
            add_gc_content_group(canvas, record, offset_y, offset_x, canvas_config, gc_config, config_dict, cfg=record_cfg)
        if canvas_config.show_skew:
            add_gc_skew_group(canvas, record, offset_y, offset_x, canvas_config, skew_config, config_dict, cfg=record_cfg)

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


