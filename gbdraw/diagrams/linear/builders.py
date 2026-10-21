#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for linear diagram assembly."""

from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING, Optional, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

from ...canvas import LinearCanvasConfigurator
from ...features.factory import FeatureBuildResult
from ...layout.linear_multi_record import LinearRecordPlacement
from ...layout.linear import LinearFeatureLaneGeometry, LinearRecordRenderContext, place_linear_definition
from ...linear_comparison import LinearComparison
from ...configurators import (
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
)
from ...render.groups.linear import (
    DefinitionGroup,
    DepthGroup,
    GcContentGroup,
    GcSkewGroup,
    LinearFeatureDomIndex,
    PairWiseMatchGroup,
    SeqRecordGroup,
)
from ...render.groups.linear.definition import LinearDefinitionStyle
from ...render.group_content import (
    CANVAS_SCOPE,
    DefinitionMeasurement,
    GroupContent,
    GroupKey,
    GroupKind,
    build_group_content,
    place_content,
)
from ...render_stages import render_stage
from .positioning import (
    position_linear_track_group,
    position_record_group,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation

def add_record_group(
    canvas: Drawing,
    record: SeqRecord,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    feature_layers: FeatureBuildResult,
    render_context: LinearRecordRenderContext,
    feature_dom_index: LinearFeatureDomIndex,
    precalculated_labels: Optional[list],
    draw_features: bool = True,
    label_font_size: float | None = None,
    orthogroup_label_member_ids: set[str | int] | None = None,
    orthogroup_label_top_member_ids: set[str | int] | None = None,
    record_index: int = 0,
    record_count: int = 1,
    group_id: str | None = None,
    dom_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_local_ruler: bool = False,
    multi_record_layout: bool = False,
    feature_offset_y: float = 0.0,
    feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
    record_transform: RecordDisplayTransform | None = None,
    record_translation_x: float = 0.0,
    record_translation_y: float = 0.0,
    *,
    placement: LinearRecordPlacement,
    group_key: GroupKey,
) -> Drawing:
    """Adds a record group to the linear canvas."""
    content = render_stage(
        "record_content",
        build_group_content,
        SeqRecordGroup,
        gb_record=record,
        geometry=canvas_config.geometry(),
        feature_config=feature_config,
        feature_layers=feature_layers,
        profile=render_context.profile,
        feature_track_layout=render_context.feature_track_layout,
        axis_ruler_enabled=render_context.axis_ruler_enabled,
        feature_dom_index=feature_dom_index,
        precalculated_labels=precalculated_labels,
        draw_features=draw_features,
        label_font_size=label_font_size,
        orthogroup_label_member_ids=orthogroup_label_member_ids,
        orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
        record_index=record_index,
        record_count=record_count,
        group_id=dom_group_id or group_id,
        sequence_width=placement.sequence_width,
        record_local_ruler=record_local_ruler,
        feature_offset_y=feature_offset_y,
        feature_lane_geometry=feature_lane_geometry,
        record_transform=record_transform,
    )
    record_group = place_content(canvas, content, group_key)
    if slot_id:
        record_group.attribs["data-gbdraw-slot-id"] = str(slot_id)
    if slot_renderer:
        record_group.attribs["data-gbdraw-slot-renderer"] = str(slot_renderer)
    if record_count > 1:
        record_group.attribs["data-record-key"] = str(placement.record_key)
        record_group.attribs["data-record-translation-x"] = str(
            float(record_translation_x)
        )
        record_group.attribs["data-record-translation-y"] = str(
            float(record_translation_y)
        )
    if multi_record_layout:
        record_group.attribs["data-record-index"] = placement.record_index
        record_group.attribs["data-record-row"] = placement.row
        record_group.attribs["data-record-column"] = placement.column
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
    track_offset_y: float,
    gc_df: DataFrame | None = None,
    track_height: float | None = None,
    group_id: str = "gc_content",
    dom_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
    *,
    sequence_width: float,
    group_key: GroupKey,
) -> Drawing:
    """Adds a GC content group to the linear canvas."""
    content = render_stage(
        "gc_content_content",
        build_group_content,
        GcContentGroup,
        gb_record=record,
        track_height=canvas_config.gc_height if track_height is None else float(track_height),
        gc_config=gc_config,
        gc_df=gc_df,
        group_id=dom_group_id or group_id,
        sequence_width=sequence_width,
        record_transform=record_transform,
    )
    gc_content_group = place_content(canvas, content, group_key)
    if slot_id:
        gc_content_group.attribs["data-gbdraw-slot-id"] = str(slot_id)
    if slot_renderer:
        gc_content_group.attribs["data-gbdraw-slot-renderer"] = str(slot_renderer)
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
    track_offset_y: float,
    depth_df: DataFrame | None = None,
    group_id: str = "depth",
    axis_group_id: str = "depth_axis",
    dom_group_id: str | None = None,
    dom_axis_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    track_height: float | None = None,
    record_transform: RecordDisplayTransform | None = None,
    *,
    sequence_width: float,
    group_key: GroupKey,
) -> Drawing:
    """Adds a depth coverage group to the linear canvas."""
    content = render_stage(
        "depth_content",
        build_group_content,
        DepthGroup,
        gb_record=record,
        track_height=canvas_config.depth_height if track_height is None else float(track_height),
        depth_config=depth_config,
        depth_df=depth_df,
        group_id=dom_group_id or group_id,
        axis_group_id=dom_axis_group_id or axis_group_id,
        sequence_width=sequence_width,
        record_transform=record_transform,
    )
    depth_group = place_content(canvas, content, group_key)
    if slot_id:
        depth_group.attribs["data-gbdraw-slot-id"] = str(slot_id)
    if slot_renderer:
        depth_group.attribs["data-gbdraw-slot-renderer"] = str(slot_renderer)
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
    track_offset_y: float,
    gc_df: DataFrame | None = None,
    track_height: float | None = None,
    group_id: str = "gc_skew",
    dom_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
    *,
    sequence_width: float,
    group_key: GroupKey,
) -> Drawing:
    """Adds a GC skew group to the linear canvas."""
    content = render_stage(
        "gc_skew_content",
        build_group_content,
        GcSkewGroup,
        gb_record=record,
        track_height=canvas_config.skew_height if track_height is None else float(track_height),
        skew_config=skew_config,
        gc_df=gc_df,
        group_id=dom_group_id or group_id,
        sequence_width=sequence_width,
        record_transform=record_transform,
    )
    gc_skew_group = place_content(canvas, content, group_key)
    if slot_id:
        gc_skew_group.attribs["data-gbdraw-slot-id"] = str(slot_id)
    if slot_renderer:
        gc_skew_group.attribs["data-gbdraw-slot-renderer"] = str(slot_renderer)
    position_linear_track_group(gc_skew_group, offset, offset_x, canvas_config, track_offset_y)
    canvas.add(gc_skew_group)
    return canvas


def _measurement(content: GroupContent) -> DefinitionMeasurement:
    if content.definition is None:  # pragma: no cover - DefinitionGroup measures itself
        raise RuntimeError("a record definition content has no measurement")
    return content.definition


def add_record_definition_group(
    canvas: Drawing,
    record: SeqRecord,
    record_offset_y: float,
    record_offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    max_def_width,
    group_id: str | None = None,
    placement: LinearRecordPlacement | None = None,
    row_definition_width: float | None = None,
    definition_center_y: float | None = None,
    definition_header_center_y: float | None = None,
    multi_record_layout: bool = False,
    local_line_kinds: Collection[str] | None = None,
    row_line_kinds: Collection[str] | None = None,
    record_index: int = 0,
    record_count: int = 1,
    record_transform: RecordDisplayTransform | None = None,
    record_presentation: RecordPresentation | None = None,
    definition_column_left: float = 0.0,
    *,
    record_key: str,
) -> Drawing:
    """Adds a record definition group to the linear canvas.

    In multi-record layouts the caller owns the split between the row-level and
    record-local definition parts, so that measurement and drawing agree on the
    same line kinds.
    """
    style = LinearDefinitionStyle.of(canvas_config.profile.config)
    length_param = str(canvas_config.length_param)
    keep_definition_left_aligned = bool(getattr(canvas_config, "keep_definition_left_aligned", False))
    definition_gap = max(0.0, float(canvas_config.definition_column_gap))

    if multi_record_layout:
        if placement is None:
            raise ValueError("multi_record_layout requires a resolved record placement")
        # A row heading is drawn from the leading record, and only when
        # something describes the row as a whole.
        split_row_definition = placement.column == 0 and bool(row_line_kinds)
        local_content = render_stage(
            "definition_content",
            build_group_content,
            DefinitionGroup,
            record,
            length_param,
            style=style,
            text_anchor="middle",
            text_x=0.0,
            group_id=group_id,
            line_kinds=local_line_kinds,
            record_index=record_index,
            record_count=record_count,
            record_transform=record_transform,
            record_presentation=record_presentation,
        )
        local_measurement = _measurement(local_content)
        header_y = (
            float(definition_header_center_y)
            if definition_header_center_y is not None
            else (
                placement.axis_y
                - placement.top_extent
                + 0.5 * local_measurement.height
            )
        )
        local_position = place_linear_definition(
            width=local_measurement.width,
            column_width=0.0,
            record_x=placement.x, gap=definition_gap,
            keep_left=keep_definition_left_aligned, sequence_width=placement.sequence_width,
            column_left=definition_column_left,
        )
        record_definition_group = place_content(
            canvas, local_content, GroupKey(record_key, GroupKind.DEFINITION, "main")
        )
        record_definition_group.translate(
            canvas_config.horizontal_offset + local_position.x,
            header_y,
        )
        canvas.add(record_definition_group)
        if split_row_definition:
            row_content = render_stage(
                "definition_content",
                build_group_content,
                DefinitionGroup,
                record,
                length_param,
                style=style,
                text_anchor="start" if keep_definition_left_aligned else "middle",
                text_x=0.0,
                group_id=f"{group_id or str(record.id)}_row",
                line_kinds=row_line_kinds,
                record_index=record_index,
                record_count=record_count,
                definition_part="row",
                record_transform=record_transform,
                record_presentation=record_presentation,
            )
            row_measurement = _measurement(row_content)
            if not row_measurement.line_count:
                # Nothing describes this row as a whole, so it gets no heading.
                return canvas
            reserved_width = (
                max(0.0, float(row_definition_width))
                if row_definition_width is not None
                else row_measurement.width
            )
            row_position = place_linear_definition(
                width=row_measurement.width,
                column_width=reserved_width,
                record_x=placement.x, gap=definition_gap,
                keep_left=keep_definition_left_aligned,
                column_left=definition_column_left,
            )
            row_group = place_content(
                canvas, row_content, GroupKey(record_key, GroupKind.DEFINITION, "row")
            )
            row_group.translate(
                canvas_config.horizontal_offset + row_position.x,
                (
                    placement.axis_y
                    if definition_center_y is None
                    else float(definition_center_y)
                ),
            )
            canvas.add(row_group)
        return canvas

    content = render_stage(
        "definition_content",
        build_group_content,
        DefinitionGroup,
        record,
        length_param,
        style=style,
        text_anchor="start" if keep_definition_left_aligned else None,
        group_id=group_id,
        record_index=record_index,
        record_count=record_count,
        record_transform=record_transform,
        record_presentation=record_presentation,
    )
    measurement = _measurement(content)
    definition_position = place_linear_definition(
        width=measurement.width,
        column_width=max(0.0, float(max_def_width)),
        record_x=record_offset_x,
        gap=definition_gap, keep_left=keep_definition_left_aligned,
        text_anchor=measurement.text_anchor,
        column_left=definition_column_left,
    )
    record_definition_group = place_content(
        canvas, content, GroupKey(record_key, GroupKind.DEFINITION, "main")
    )
    record_definition_group.translate(
        canvas_config.horizontal_offset + definition_position.x,
        record_offset_y if definition_center_y is None else float(definition_center_y),
    )

    canvas.add(record_definition_group)
    return canvas


def linear_comparison_content(
    comparisons: tuple[LinearComparison, ...],
    index: int,
    blast_config: BlastMatchConfigurator,
    records: tuple[SeqRecord, ...],
    *,
    height: float,
    query_placement: LinearRecordPlacement,
    subject_placement: LinearRecordPlacement,
    query_y: float,
    subject_y: float,
    feature_dom_index: LinearFeatureDomIndex | None,
    record_transforms: tuple[RecordDisplayTransform, ...] | None,
) -> GroupContent:
    """The ribbons of comparison ``index`` (render stage ``comparison_content``).

    It takes the whole ``comparison_rows`` result and an index, because that
    result, not each comparison in it, has the stage identity.
    """

    comparison = comparisons[index]
    return build_group_content(
        PairWiseMatchGroup,
        blast_config.sequence_length_dict,
        comparison.matches,
        height,
        index + 1,
        blast_config,
        records,
        query_placement=query_placement,
        subject_placement=subject_placement,
        query_y=query_y,
        subject_y=subject_y,
        feature_dom_index=feature_dom_index,
        record_transforms=record_transforms,
    )


def add_explicit_comparisons_on_linear_canvas(
    canvas: Drawing,
    comparisons: tuple[LinearComparison, ...],
    canvas_config: LinearCanvasConfigurator,
    blast_config: BlastMatchConfigurator,
    records: list[SeqRecord],
    placements: dict[int, LinearRecordPlacement],
    feature_dom_index: LinearFeatureDomIndex | None = None,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    *,
    record_keys: Sequence[str],
) -> Drawing:
    """Draw comparisons using explicit endpoint placements."""

    for index, comparison in enumerate(comparisons):
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
        content = render_stage(
            "comparison_content",
            linear_comparison_content,
            comparisons,
            index,
            blast_config,
            tuple(records),
            height=abs(subject_anchor - query_anchor),
            query_placement=query,
            subject_placement=subject,
            query_y=query_anchor - top_y,
            subject_y=subject_anchor - top_y,
            feature_dom_index=feature_dom_index,
            record_transforms=tuple(record_transforms) if record_transforms is not None else None,
        )
        match_group = place_content(
            canvas,
            content,
            GroupKey(
                (
                    str(record_keys[comparison.query_record_index]),
                    str(record_keys[comparison.subject_record_index]),
                ),
                GroupKind.COMPARISON,
            ),
        )
        match_group.translate(canvas_config.horizontal_offset, top_y)
        canvas.add(match_group)
    return canvas


def add_length_bar_on_linear_canvas(
    canvas: Drawing,
    canvas_config: LinearCanvasConfigurator,
    content: GroupContent,
    offset_x: float = 0.0,
) -> Drawing:
    """Adds a length bar to the linear canvas."""
    offset_for_length_bar = canvas_config.height_below_final_record
    scale_group = place_content(canvas, content, GroupKey(CANVAS_SCOPE, GroupKind.LENGTH_BAR))
    scale_group.translate(canvas_config.horizontal_offset + offset_x, offset_for_length_bar)
    canvas.add(scale_group)
    return canvas


__all__ = [
    "linear_comparison_content",
    "add_record_group",
    "add_depth_group",
    "add_gc_content_group",
    "add_gc_skew_group",
    "add_record_definition_group",
    "add_explicit_comparisons_on_linear_canvas",
    "add_length_bar_on_linear_canvas",
]
