#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import logging
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
from ...groups.linear import LengthBarGroup, LegendGroup  # type: ignore[reportMissingImports]
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

logger = logging.getLogger(__name__)


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
    # Prepare legend table
    legend_table = prepare_legend_table(gc_config, skew_config, feature_config, features_present, blast_config, has_blast)
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
            height_below_axis = canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding
            next_record_id = record_ids[i + 1]
            next_label_height = record_label_heights.get(next_record_id, 0)
            if next_label_height > canvas_config.comparison_height:
                inter_record_space = height_below_axis + next_label_height + canvas_config.cds_padding
            else:
                inter_record_space = height_below_axis + canvas_config.comparison_height + canvas_config.cds_padding
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


