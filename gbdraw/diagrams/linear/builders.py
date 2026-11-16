#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for linear diagram assembly."""

from __future__ import annotations

from collections.abc import Collection, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

from ...canvas import LinearCanvasConfigurator
from ...features.factory import FeatureBuildResult
from ...render.drawers.linear.feature_density import FeatureDensityHeatmapGroup
from ...labels.linear_types import LinearLabel
from ...layout.linear_multi_record import (
    LinearComparisonEndpoint,
    LinearFeatureAttachment,
    LinearRecordPlacement,
    linear_comparison_frame,
)
from ...layout.linear import (
    LinearFeatureLaneGeometry,
    LinearRecordRenderContext,
    LinearResolvedTrack,
    place_linear_definition,
)
from ...linear_comparison import LinearComparison
from ...configurators import (
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
)
from ...render.groups.linear import (
    DefinitionGroup,
    LinearFeatureDomIndex,
    NumericSeriesGroup,
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
from ...render.numeric_series import NumericSeriesConfig, NumericSeriesStage
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
    precalculated_labels: list[LinearLabel] | None,
    label_font_size: float | None = None,
    orthogroup_label_member_ids: set[str | int] | None = None,
    orthogroup_label_top_member_ids: set[str | int] | None = None,
    record_index: int = 0,
    record_count: int = 1,
    record_local_ruler: bool = False,
    multi_record_layout: bool = False,
    feature_offset_y: float = 0.0,
    feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
    record_transform: RecordDisplayTransform | None = None,
    record_translation_x: float = 0.0,
    record_translation_y: float = 0.0,
    axis_row_offsets: tuple[float, ...] = (),
    feature_row_shifts: tuple[float, ...] = (),
    *,
    placement: LinearRecordPlacement,
    group_key: GroupKey,
    group_id: str,
    feature_slot: LinearResolvedTrack | None,
) -> Drawing:
    """Adds a record group to the linear canvas.

    ``feature_slot`` is the features slot the group draws the record's features in, and the
    group carries its slot tags; ``None`` draws the record's axis alone, for a figure whose
    slots hold no feature track.
    """
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
        draw_features=feature_slot is not None,
        label_font_size=label_font_size,
        orthogroup_label_member_ids=orthogroup_label_member_ids,
        orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
        record_index=record_index,
        record_count=record_count,
        group_id=group_id,
        sequence_width=placement.sequence_width,
        record_local_ruler=record_local_ruler,
        feature_offset_y=feature_offset_y,
        feature_lane_geometry=feature_lane_geometry,
        record_transform=record_transform,
        axis_row_offsets=axis_row_offsets,
        feature_row_shifts=feature_row_shifts,
    )
    record_group = place_content(canvas, content, group_key)
    if feature_slot is not None:
        record_group.attribs["data-gbdraw-slot-id"] = str(feature_slot.id)
        record_group.attribs["data-gbdraw-slot-renderer"] = str(feature_slot.renderer)
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


@dataclass(frozen=True)
class LinearSeriesSlot:
    """A slot's numeric series track as its group is built and placed: the record's offsets on the canvas, the
    track's offset below the record's axis and its height, the record's drawn axis width, transform and row
    shifts, the group's id and key, and the slot it tags."""

    offset_y: float
    offset_x: float
    track_offset_y: float
    track_height: float
    sequence_width: float
    record_transform: RecordDisplayTransform | None
    row_shifts: tuple[float, ...]
    group_id: str
    group_key: GroupKey
    slot_id: str
    slot_renderer: str


def add_numeric_series_group(
    canvas: Drawing,
    stage: NumericSeriesStage,
    record: SeqRecord,
    canvas_config: LinearCanvasConfigurator,
    config: NumericSeriesConfig,
    frame: DataFrame,
    slot: LinearSeriesSlot,
) -> Drawing:
    """Add a numeric series track (GC content, GC skew, depth) to the linear canvas, built in the render stage
    ``stage``."""

    content = render_stage(
        stage,
        build_group_content,
        NumericSeriesGroup,
        gb_record=record,
        track_height=slot.track_height,
        config=config,
        frame=frame,
        group_id=slot.group_id,
        sequence_width=slot.sequence_width,
        record_transform=slot.record_transform,
        row_shifts=slot.row_shifts,
    )
    group = place_content(canvas, content, slot.group_key)
    group.attribs["data-gbdraw-slot-id"] = slot.slot_id
    group.attribs["data-gbdraw-slot-renderer"] = slot.slot_renderer
    position_linear_track_group(group, slot.offset_y, slot.offset_x, canvas_config, slot.track_offset_y)
    canvas.add(group)
    return canvas


def add_feature_density_heatmap_group(
    canvas: Drawing,
    heatmap: FeatureDensityHeatmapGroup,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    *,
    track_offset_y: float,
    group_key: GroupKey,
    slot_id: str,
    slot_renderer: str,
) -> Drawing:
    """Adds a feature density heatmap band to the linear canvas."""
    heatmap_group = place_content(canvas, GroupContent.of(heatmap.get_group()), group_key)
    heatmap_group.attribs["data-gbdraw-slot-id"] = slot_id
    heatmap_group.attribs["data-gbdraw-slot-renderer"] = slot_renderer
    position_linear_track_group(heatmap_group, offset_y, offset_x, canvas_config, track_offset_y)
    canvas.add(heatmap_group)
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
    query_placement: LinearComparisonEndpoint,
    subject_placement: LinearComparisonEndpoint,
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
        index + 1,
        blast_config,
        records,
        query_placement=query_placement,
        subject_placement=subject_placement,
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
    feature_attachments: Mapping[int, Mapping[str, LinearFeatureAttachment]],
) -> Drawing:
    """Draw comparisons using explicit endpoint placements.

    Each comparison draws in its own frame (``linear_comparison_frame``), and the group is placed at the frame's
    top, so rows that move together change only where the group is placed (OV-2089).
    """

    for index, comparison in enumerate(comparisons):
        query = placements[comparison.query_record_index]
        subject = placements[comparison.subject_record_index]
        frame = linear_comparison_frame(query, subject)
        content = render_stage(
            "comparison_content",
            linear_comparison_content,
            comparisons,
            index,
            blast_config,
            tuple(records),
            query_placement=query.comparison_endpoint(
                subject, anchor_y=frame.query_y, attachments=feature_attachments.get(query.record_index, {})
            ),
            subject_placement=subject.comparison_endpoint(
                query, anchor_y=frame.subject_y, attachments=feature_attachments.get(subject.record_index, {})
            ),
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
        match_group.translate(canvas_config.horizontal_offset, frame.top_y)
        canvas.add(match_group)
    return canvas


def add_row_window_heading_group(
    canvas: Drawing,
    content: GroupContent,
    canvas_config: LinearCanvasConfigurator,
    placement: LinearRecordPlacement,
    *,
    row_definition_width: float,
    center_y: float,
    definition_column_left: float,
    record_key: str,
) -> Drawing:
    """Adds a windowed row's heading (``RowWindowHeadingGroup``) in the definition column, its name line on
    ``center_y``, where a row heading of records stands (``place_linear_definition``)."""

    measurement = _measurement(content)
    position = place_linear_definition(
        width=measurement.width,
        column_width=max(0.0, float(row_definition_width)),
        record_x=placement.x,
        gap=max(0.0, float(canvas_config.definition_column_gap)),
        keep_left=bool(getattr(canvas_config, "keep_definition_left_aligned", False)),
        column_left=definition_column_left,
    )
    group = place_content(canvas, content, GroupKey(record_key, GroupKind.DEFINITION, "row"))
    group.translate(canvas_config.horizontal_offset + position.x, center_y)
    canvas.add(group)
    return canvas


def add_row_window_marks_group(
    canvas: Drawing, content: GroupContent, canvas_config: LinearCanvasConfigurator, *, record_key: str
) -> Drawing:
    """Adds a windowed row's marks (``RowWindowMarksGroup``), drawn in canvas coordinates, above the records."""

    group = place_content(canvas, content, GroupKey(record_key, GroupKind.DEFINITION, "marks"))
    group.translate(canvas_config.horizontal_offset, 0.0)
    canvas.add(group)
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
    "LinearSeriesSlot",
    "add_feature_density_heatmap_group",
    "add_numeric_series_group",
    "add_record_definition_group",
    "add_explicit_comparisons_on_linear_canvas",
    "add_length_bar_on_linear_canvas",
    "add_row_window_heading_group",
    "add_row_window_marks_group",
]
