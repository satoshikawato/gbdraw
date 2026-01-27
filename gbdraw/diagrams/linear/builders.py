#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for linear diagram assembly."""

from __future__ import annotations

from typing import Optional

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
)
from ...render.groups.linear import (  # type: ignore[reportMissingImports]
    DefinitionGroup,
    GcContentGroup,
    GcSkewGroup,
    PairWiseMatchGroup,
    SeqRecordGroup,
)
from .positioning import (
    position_gc_content_group,
    position_gc_skew_group,
    position_record_definition_group,
    position_record_group,
)


def add_record_group(
    canvas: Drawing,
    record: SeqRecord,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    config_dict: dict,
    precalculated_labels: Optional[list],
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Adds a record group to the linear canvas."""
    record_group: Group = SeqRecordGroup(
        gb_record=record,
        canvas_config=canvas_config,
        feature_config=feature_config,
        config_dict=config_dict,
        precalculated_labels=precalculated_labels,
        cfg=cfg,
    ).get_group()
    position_record_group(record_group, offset_y, offset_x, canvas_config)
    canvas.add(record_group)
    return canvas


def add_gc_content_group(
    canvas: Drawing,
    record: SeqRecord,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    gc_config: GcContentConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Adds a GC content group to the linear canvas."""
    gc_content_group: Group = GcContentGroup(
        gb_record=record,
        alignment_width=canvas_config.alignment_width,
        longest_record_len=canvas_config.longest_genome,
        track_height=canvas_config.gc_height,
        gc_config=gc_config,
        config_dict=config_dict,
        cfg=cfg,
    ).get_group()
    position_gc_content_group(gc_content_group, offset_y, offset_x, canvas_config)
    canvas.add(gc_content_group)
    return canvas


def add_gc_skew_group(
    canvas: Drawing,
    record: SeqRecord,
    offset: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    skew_config: GcSkewConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Adds a GC skew group to the linear canvas."""
    gc_skew_group: Group = GcSkewGroup(
        gb_record=record,
        alignment_width=canvas_config.alignment_width,
        longest_record_len=canvas_config.longest_genome,
        track_height=canvas_config.skew_height,
        skew_config=skew_config,
        config_dict=config_dict,
        cfg=cfg,
    ).get_group()
    position_gc_skew_group(gc_skew_group, offset, offset_x, canvas_config)
    canvas.add(gc_skew_group)
    return canvas


def add_record_definition_group(
    canvas: Drawing,
    record: SeqRecord,
    record_offset_y: float,
    record_offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    max_def_width,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Adds a record definition group to the linear canvas."""
    definition_group_obj = DefinitionGroup(record, config_dict, canvas_config, cfg=cfg)
    definition_gap = min(canvas_config.canvas_padding, 20)
    definition_offset_x = (definition_group_obj.definition_bounding_box_width / 2) + definition_gap
    record_definition_group: Group = definition_group_obj.get_group()

    position_record_definition_group(
        record_definition_group,
        record_offset_y,
        (definition_offset_x - record_offset_x),
        canvas_config,
    )

    canvas.add(record_definition_group)
    return canvas


def add_comparison_on_linear_canvas(
    canvas: Drawing,
    comparisons,
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    config_dict: dict,
    records: list,
    comparison_offsets: list,
    actual_comparison_heights: list,
) -> Drawing:
    """Adds comparison groups at specified y-offsets with dynamic height."""
    for comparison_count, comparison in enumerate(comparisons, start=1):
        if comparison_count > len(comparison_offsets):
            break

        height = actual_comparison_heights[comparison_count - 1]
        offset = comparison_offsets[comparison_count - 1]

        match_group: Group = PairWiseMatchGroup(
            canvas_config,
            blast_config.sequence_length_dict,
            comparison,
            height,
            comparison_count,
            blast_config,
            records,
        ).get_group()

        match_group.translate(canvas_config.horizontal_offset, offset)
        canvas.add(match_group)
    return canvas


def add_length_bar_on_linear_canvas(
    canvas: Drawing, canvas_config: LinearCanvasConfigurator, config_dict: dict, scale_group, legend_group
) -> Drawing:
    """Adds a length bar to the linear canvas."""
    if canvas_config.legend_position == "bottom" or canvas_config.legend_position == "top":
        offset_for_length_bar = canvas_config.height_below_final_record
    else:
        offset_for_length_bar = canvas_config.height_below_final_record
    scale_group = scale_group.get_group()
    scale_group.translate(canvas_config.horizontal_offset, offset_for_length_bar)
    canvas.add(scale_group)
    return canvas


def add_legends_on_linear_canvas(canvas: Drawing, config_dict, canvas_config: LinearCanvasConfigurator, legend_group, legend_table):
    legend_group = legend_group.get_group()
    offset_x = canvas_config.legend_offset_x
    offset_y = canvas_config.legend_offset_y
    legend_group.translate(offset_x, offset_y)
    canvas.add(legend_group)
    return canvas


__all__ = [
    "add_record_group",
    "add_gc_content_group",
    "add_gc_skew_group",
    "add_record_definition_group",
    "add_comparison_on_linear_canvas",
    "add_length_bar_on_linear_canvas",
    "add_legends_on_linear_canvas",
]


