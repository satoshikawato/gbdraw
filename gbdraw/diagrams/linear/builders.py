#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for linear diagram assembly."""

from __future__ import annotations

from typing import Optional

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...layout.linear_multi_record import LinearRecordPlacement
from ...linear_comparison import LinearComparison
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
)
from ...render.groups.linear import (  # type: ignore[reportMissingImports]
    DefinitionGroup,
    DepthGroup,
    GcContentGroup,
    GcSkewGroup,
    PairWiseMatchGroup,
    SeqRecordGroup,
)
from .positioning import (
    position_depth_group,
    position_gc_content_group,
    position_gc_skew_group,
    position_linear_track_group,
    position_record_definition_group,
    position_record_group,
)
from .precalc import FeatureDict


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
    precomputed_feature_dict: FeatureDict | None = None,
    feature_track_layout: str | None = None,
    draw_features: bool = True,
    label_font_size: float | None = None,
    orthogroup_label_member_ids: set[str] | None = None,
    orthogroup_label_top_member_ids: set[str] | None = None,
    record_index: int = 0,
    record_count: int = 1,
    group_id: str | None = None,
    sequence_width: float | None = None,
    record_local_ruler: bool = False,
    placement: LinearRecordPlacement | None = None,
) -> Drawing:
    """Adds a record group to the linear canvas."""
    if placement is not None:
        sequence_width = placement.sequence_width
        record_local_ruler = True
    record_group: Group = SeqRecordGroup(
        gb_record=record,
        canvas_config=canvas_config,
        feature_config=feature_config,
        config_dict=config_dict,
        precalculated_labels=precalculated_labels,
        cfg=cfg,
        precomputed_feature_dict=precomputed_feature_dict,
        feature_track_layout=feature_track_layout,
        draw_features=draw_features,
        label_font_size=label_font_size,
        orthogroup_label_member_ids=orthogroup_label_member_ids,
        orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
        record_index=record_index,
        record_count=record_count,
        group_id=group_id,
        sequence_width=sequence_width,
        record_local_ruler=record_local_ruler,
    ).get_group()
    if placement is not None:
        record_group.attribs["data-record-index"] = placement.record_index
        record_group.attribs["data-record-row"] = placement.row
        record_group.attribs["data-record-column"] = placement.column
        record_group.attribs["data-record-key"] = str(placement.record_key)
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
    gc_df: DataFrame | None = None,
    track_height: float | None = None,
    track_offset_y: float | None = None,
    group_id: str = "gc_content",
    sequence_width: float | None = None,
) -> Drawing:
    """Adds a GC content group to the linear canvas."""
    gc_content_group: Group = GcContentGroup(
        gb_record=record,
        alignment_width=canvas_config.alignment_width,
        longest_record_len=canvas_config.longest_genome,
        track_height=canvas_config.gc_height if track_height is None else float(track_height),
        gc_config=gc_config,
        config_dict=config_dict,
        cfg=cfg,
        gc_df=gc_df,
        group_id=group_id,
        sequence_width=sequence_width,
    ).get_group()
    if track_offset_y is None:
        position_gc_content_group(gc_content_group, offset_y, offset_x, canvas_config)
    else:
        position_linear_track_group(gc_content_group, offset_y, offset_x, canvas_config, track_offset_y)
    canvas.add(gc_content_group)
    return canvas


def add_depth_group(
    canvas: Drawing,
    record: SeqRecord,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    depth_config: DepthConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
    depth_df: DataFrame | None = None,
    depth_track_index: int = 0,
    group_id: str = "depth",
    axis_group_id: str = "depth_axis",
    track_height: float | None = None,
    track_offset_y: float | None = None,
    sequence_width: float | None = None,
) -> Drawing:
    """Adds a depth coverage group to the linear canvas."""
    depth_group: Group = DepthGroup(
        gb_record=record,
        alignment_width=canvas_config.alignment_width,
        longest_record_len=canvas_config.longest_genome,
        track_height=canvas_config.depth_height if track_height is None else float(track_height),
        depth_config=depth_config,
        config_dict=config_dict,
        cfg=cfg,
        depth_df=depth_df,
        group_id=group_id,
        axis_group_id=axis_group_id,
        sequence_width=sequence_width,
    ).get_group()
    if track_offset_y is None:
        position_depth_group(depth_group, offset_y, offset_x, canvas_config, depth_track_index)
    else:
        position_linear_track_group(depth_group, offset_y, offset_x, canvas_config, track_offset_y)
    canvas.add(depth_group)
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
    gc_df: DataFrame | None = None,
    track_height: float | None = None,
    track_offset_y: float | None = None,
    group_id: str = "gc_skew",
    sequence_width: float | None = None,
) -> Drawing:
    """Adds a GC skew group to the linear canvas."""
    gc_skew_group: Group = GcSkewGroup(
        gb_record=record,
        alignment_width=canvas_config.alignment_width,
        longest_record_len=canvas_config.longest_genome,
        track_height=canvas_config.skew_height if track_height is None else float(track_height),
        skew_config=skew_config,
        config_dict=config_dict,
        cfg=cfg,
        gc_df=gc_df,
        group_id=group_id,
        sequence_width=sequence_width,
    ).get_group()
    if track_offset_y is None:
        position_gc_skew_group(gc_skew_group, offset, offset_x, canvas_config)
    else:
        position_linear_track_group(gc_skew_group, offset, offset_x, canvas_config, track_offset_y)
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
    group_id: str | None = None,
    placement: LinearRecordPlacement | None = None,
    row_definition_width: float | None = None,
    definition_center_y: float | None = None,
) -> Drawing:
    """Adds a record definition group to the linear canvas."""
    keep_definition_left_aligned = bool(getattr(canvas_config, "keep_definition_left_aligned", False))
    try:
        definition_gap = max(0.0, float(getattr(canvas_config, "definition_gap", 20.0)))
    except (TypeError, ValueError):
        definition_gap = 20.0

    if placement is not None:
        split_row_definition = (
            keep_definition_left_aligned and placement.column == 0
        )
        local_line_kinds = (
            {"replicon", "accession", "length"}
            if split_row_definition
            else None
        )
        definition_group_obj = DefinitionGroup(
            record,
            config_dict,
            canvas_config,
            cfg=cfg,
            text_anchor="middle",
            text_x=0.0,
            group_id=group_id,
            line_kinds=local_line_kinds,
        )
        record_definition_group = definition_group_obj.get_group()
        header_y = (
            placement.axis_y
            - placement.top_extent
            + 0.5 * definition_group_obj.definition_bounding_box_height
        )
        record_definition_group.translate(
            canvas_config.horizontal_offset
            + placement.x
            + 0.5 * placement.sequence_width,
            header_y,
        )
        canvas.add(record_definition_group)
        if split_row_definition:
            row_group_obj = DefinitionGroup(
                record,
                config_dict,
                canvas_config,
                cfg=cfg,
                text_anchor="start",
                text_x=0.0,
                group_id=f"{group_id or str(record.id)}_row",
                line_kinds={"name", "subtitle"},
            )
            reserved_width = (
                max(0.0, float(row_definition_width))
                if row_definition_width is not None
                else row_group_obj.definition_bounding_box_width
            )
            row_group = row_group_obj.get_group()
            row_group.translate(
                canvas_config.horizontal_offset - definition_gap - reserved_width,
                (
                    placement.axis_y
                    if definition_center_y is None
                    else float(definition_center_y)
                ),
            )
            canvas.add(row_group)
        return canvas

    if keep_definition_left_aligned:
        try:
            definition_column_width = max(0.0, float(max_def_width))
        except (TypeError, ValueError):
            definition_column_width = 0.0

        if definition_column_width == 0.0:
            provisional_group_obj = DefinitionGroup(
                record,
                config_dict,
                canvas_config,
                cfg=cfg,
                group_id=group_id,
            )
            definition_column_width = provisional_group_obj.definition_bounding_box_width

        definition_group_obj = DefinitionGroup(
            record,
            config_dict,
            canvas_config,
            cfg=cfg,
            text_anchor="start",
            text_x=0.0,
            group_id=group_id,
        )
        positioned_definition_offset_x = definition_column_width + definition_gap
    else:
        definition_group_obj = DefinitionGroup(
            record,
            config_dict,
            canvas_config,
            cfg=cfg,
            group_id=group_id,
        )
        definition_offset_x = (definition_group_obj.definition_bounding_box_width / 2) + definition_gap
        positioned_definition_offset_x = definition_offset_x - record_offset_x

    record_definition_group: Group = definition_group_obj.get_group()

    position_record_definition_group(
        record_definition_group,
        (
            record_offset_y
            if definition_center_y is None
            else float(definition_center_y)
        ),
        positioned_definition_offset_x,
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
    record_offsets_x: dict[int, float] | None = None,
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
            record_offsets_x=record_offsets_x,
        ).get_group()

        match_group.translate(canvas_config.horizontal_offset, offset)
        canvas.add(match_group)
    return canvas


def add_explicit_comparisons_on_linear_canvas(
    canvas: Drawing,
    comparisons: list[LinearComparison],
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    records: list[SeqRecord],
    placements: dict[int, LinearRecordPlacement],
) -> Drawing:
    """Draw comparisons using explicit endpoint placements."""

    for comparison_count, comparison in enumerate(comparisons, start=1):
        query = placements[comparison.query_record_index]
        subject = placements[comparison.subject_record_index]
        query_anchor = (
            query.comparison_bottom_y
            if query.row < subject.row
            else query.comparison_top_y
        )
        subject_anchor = (
            subject.comparison_bottom_y
            if subject.row < query.row
            else subject.comparison_top_y
        )
        top_y = min(query_anchor, subject_anchor)
        height = abs(subject_anchor - query_anchor)
        match_group = PairWiseMatchGroup(
            canvas_config,
            blast_config.sequence_length_dict,
            comparison.matches,
            height,
            comparison_count,
            blast_config,
            records,
            query_record_index=comparison.query_record_index,
            subject_record_index=comparison.subject_record_index,
            query_placement=query,
            subject_placement=subject,
            query_y=query_anchor - top_y,
            subject_y=subject_anchor - top_y,
        ).get_group()
        match_group.translate(canvas_config.horizontal_offset, top_y)
        canvas.add(match_group)
    return canvas


def add_length_bar_on_linear_canvas(
    canvas: Drawing,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    scale_group,
    legend_group,
    offset_x: float = 0.0,
) -> Drawing:
    """Adds a length bar to the linear canvas."""
    if canvas_config.legend_position == "bottom" or canvas_config.legend_position == "top":
        offset_for_length_bar = canvas_config.height_below_final_record
    else:
        offset_for_length_bar = canvas_config.height_below_final_record
    scale_group = scale_group.get_group()
    scale_group.translate(canvas_config.horizontal_offset + offset_x, offset_for_length_bar)
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
    "add_depth_group",
    "add_gc_content_group",
    "add_gc_skew_group",
    "add_record_definition_group",
    "add_comparison_on_linear_canvas",
    "add_explicit_comparisons_on_linear_canvas",
    "add_length_bar_on_linear_canvas",
    "add_legends_on_linear_canvas",
]


