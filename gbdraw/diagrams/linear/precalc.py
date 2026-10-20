#!/usr/bin/env python
# coding: utf-8

"""Pre-calculation helpers for linear diagram assembly.

These functions compute layout constraints (definition widths, label heights) before
final canvas sizing.
"""

from __future__ import annotations

from gbdraw.features.placement import FeaturePlacementSlot

from collections.abc import Collection
import math
from typing import TYPE_CHECKING, Sequence

from Bio.SeqRecord import SeqRecord
from ...canvas import LinearCanvasConfigurator
from ...config.models import (
    GbdrawConfig,
    LinearRenderProfile,
)
from ...configurators import FeatureDrawingConfigurator
from ...auto_sizes import linear_label_font_size, linear_label_size_lengths
from ...features.factory import (
    FeatureBuildResult,
    create_feature_layers,
)
from ...features.objects import FeatureObject
from ...render.groups.linear import DefinitionGroup, LinearDefinitionStyle
from ...labels.filtering import label_filter_maps
from ...render_stages import render_stage
from ...labels.linear import (
    LinearLabelFrame,
    LinearRecordLabels,
    calculate_label_y_bounds,
    linear_label_rows,
    linear_record_labels,
    place_linear_label_row,
)
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.linear import LinearFeatureLaneGeometry, LinearRecordRenderContext
from ...layout.spatial import Aabb
from .orthogroup_alignment import OrthogroupLabelEligibility, orthogroup_label_sets_for_record

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation


FeatureDict = dict[str, FeatureObject]


def _resolve_linear_diagram_label_font_size(
    records: list[SeqRecord],
    *,
    canvas_config: LinearCanvasConfigurator,
    profile: LinearRenderProfile,
) -> float:
    """Resolve one auto label size for all labels drawn in a linear diagram.

    The longest record in the label display scope sets the size. With no record
    in scope, labels shown by per-feature overrides use the size of all records.
    """

    cfg = profile.config
    sized_lengths = linear_label_size_lengths([len(record.seq) for record in records], profile.label_scope)
    if sized_lengths:
        return linear_label_font_size(sized_lengths, cfg)
    return cfg.labels.font_size.linear.for_length_param(canvas_config.length_param)


def _precalculate_definition_metrics(
    records: list[SeqRecord],
    canvas_config,
    cfg: GbdrawConfig,
    line_kinds_by_record: Sequence[Collection[str] | None] | None = None,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    record_presentations: Sequence[RecordPresentation | None] | None = None,
) -> tuple[float, list[float], list[float]]:
    """
    Pre-calculate definition widths and heights for all records.
    """
    max_definition_width: float = 0
    definition_heights: list[float] = []
    definition_widths: list[float] = []
    if not records:
        return 0.0, definition_widths, definition_heights

    style = LinearDefinitionStyle.of(cfg)
    if line_kinds_by_record is not None and len(line_kinds_by_record) != len(records):
        raise ValueError("line_kinds_by_record must match the number of records")
    for index, record in enumerate(records):
        line_kinds = (
            line_kinds_by_record[index]
            if line_kinds_by_record is not None
            else None
        )
        def_group = DefinitionGroup(
            record,
            str(canvas_config.length_param),
            style=style,
            line_kinds=line_kinds,
            record_transform=(record_transforms[index] if record_transforms is not None else None),
            record_presentation=(record_presentations[index] if record_presentations is not None else None),
        )
        if def_group.definition_bounding_box_width > max_definition_width:
            max_definition_width = def_group.definition_bounding_box_width
        definition_height = float(def_group.definition_bounding_box_height)
        definition_heights.append(definition_height)
        definition_widths.append(float(def_group.definition_bounding_box_width))
    return math.ceil(max_definition_width), definition_widths, definition_heights


def _precalculate_feature_layers(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    profile: LinearRenderProfile,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    placement_slot: FeaturePlacementSlot | None = None,
) -> list[FeatureBuildResult]:
    """Build feature objects once per record for the linear assembly pipeline."""

    cfg = profile.config
    label_filtering = label_filter_maps(cfg.labels.filtering.as_dict())

    feature_layers: list[FeatureBuildResult] = []
    for i, record in enumerate(records):
        compute_label_text = profile.labels_in_scope(i)
        result = render_stage(
            "feature_layers",
            create_feature_layers,
            record,
            feature_config.specific_color_rules,
            feature_config.selected_features_set,
            feature_config.default_color_map,
            profile.strandedness,
            profile.resolve_overlaps,
            label_filtering if compute_label_text else {},
            feature_shapes=feature_config.feature_shapes,
            feature_visibility_rules=feature_config.feature_visibility_rules,
            compute_label_text=compute_label_text,
            record_features=feature_config.record_features[i] if feature_config.record_features else None,
            placement_slot=placement_slot,
            feature_overlap_tolerance_bp=profile.feature_overlap_tolerance_bp,
            record_transform=(record_transforms[i] if record_transforms is not None else None),
        )
        feature_layers.append(result)
    return feature_layers


def _precalculate_label_dimensions(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    render_context: LinearRecordRenderContext,
    precomputed_feature_dicts: list[FeatureDict] | None = None,
    orthogroup_label_eligibility: OrthogroupLabelEligibility | None = None,
    sequence_widths: Sequence[float] | None = None,
    feature_lane_geometries: Sequence[LinearFeatureLaneGeometry] | None = None,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    label_frames: Sequence[LinearLabelFrame] | None = None,
    label_keep_clear: Sequence[tuple[Aabb, ...]] | None = None,
) -> tuple[float, list[list[dict]], list[float]]:
    """Pre-calculates label placements for all records to determine the required canvas height.

    The labels of the records in one row (``label_frames``; each record alone in its row without them)
    are stacked together (OV-1188). ``label_keep_clear`` are each record's boxes its labels keep
    clear of, in their coordinates (``LinearRecordLabels.keep_clear``). ``sequence_widths`` are the
    records' drawn axis widths (their placements' ``sequence_width``); without them each record is
    as wide as it is alone in its row (``LinearCanvasConfigurator.record_sequence_width``).
    """

    profile = render_context.profile
    cfg = profile.config
    label_scope = profile.label_scope

    if not any(profile.record_has_labels(i) for i in range(len(records))):
        return 0, [[] for _ in records], [0.0 for _ in records]

    label_font_size = _resolve_linear_diagram_label_font_size(
        records,
        canvas_config=canvas_config,
        profile=profile,
    )
    max_required_height: float = 0
    record_labels: list[LinearRecordLabels | None] = []
    all_labels_by_record: list[list[dict]] = []
    record_label_heights: list[float] = []
    label_filtering = label_filter_maps(cfg.labels.filtering.as_dict())

    for i, record in enumerate(records):
        if not profile.record_has_labels(i):
            record_labels.append(None)
            continue
        member_ids, top_member_ids = orthogroup_label_sets_for_record(
            orthogroup_label_eligibility if label_scope == "orthogroup_top" else None,
            i,
        )

        if precomputed_feature_dicts is not None:
            feature_dict = precomputed_feature_dicts[i]
        else:
            feature_result = create_feature_layers(
                record,
                feature_config.specific_color_rules,
                feature_config.selected_features_set,
                feature_config.default_color_map,
                profile.strandedness,
                profile.resolve_overlaps,
                label_filtering,
                feature_shapes=feature_config.feature_shapes,
                feature_overlap_tolerance_bp=profile.feature_overlap_tolerance_bp,
                feature_visibility_rules=feature_config.feature_visibility_rules,
            )
            feature_dict = feature_result.foreground_features

        record_length = len(record.seq)
        # Label x-positions use the record's drawn axis width (its placement's ``sequence_width``).
        labels = linear_record_labels(
            feature_dict,
            record_length,
            (
                float(sequence_widths[i])
                if sequence_widths is not None
                else canvas_config.record_sequence_width(record_length)
            ),
            1.0,
            canvas_config.cds_height,
            profile.strandedness,
            render_context.feature_track_layout,
            canvas_config.track_axis_gap,
            profile=profile,
            label_font_size=label_font_size,
            orthogroup_label_member_ids=member_ids,
            orthogroup_label_top_member_ids=top_member_ids,
            labels_overrides_only=not profile.labels_in_scope(i),
            feature_lane_geometry=(
                feature_lane_geometries[i]
                if feature_lane_geometries is not None
                else None
            ),
            record_transform=(record_transforms[i] if record_transforms is not None else None),
            label_filtering=label_filtering,
        )
        if label_keep_clear is not None:
            labels = labels._replace(keep_clear=label_keep_clear[i])
        record_labels.append(labels)

    placed_by_record: dict[int, list[dict]] = {}
    for row in linear_label_rows(label_frames, len(records)):
        members = [(i, x, unplaced) for i, x in row if (unplaced := record_labels[i]) is not None]
        placed = place_linear_label_row(
            [(unplaced, x) for _i, x, unplaced in members],
            profile=profile,
            track_layout=render_context.feature_track_layout,
        )
        placed_by_record.update(zip((i for i, _x, _unplaced in members), placed, strict=True))

    for i in range(len(records)):
        label_list = placed_by_record.get(i, [])
        all_labels_by_record.append(label_list)
        if i not in placed_by_record:
            record_label_heights.append(0.0)
            continue

        min_y_coord_record = 0.0
        for label in label_list:
            label_top_y, _ = calculate_label_y_bounds(label)
            is_above_feature_embedded = (
                bool(label.get("is_embedded"))
                and label_top_y < float(label.get("feature_top_y", 0.0))
            )
            if bool(label.get("is_embedded")) and not is_above_feature_embedded:
                continue
            if label_top_y < min_y_coord_record:
                min_y_coord_record = label_top_y

        record_height = abs(min_y_coord_record)
        record_label_heights.append(record_height)

        if record_height > max_required_height:
            max_required_height = record_height

    return max_required_height, all_labels_by_record, record_label_heights


__all__ = [
    "FeatureDict",
    "_precalculate_definition_metrics",
    "_precalculate_feature_layers",
    "_precalculate_label_dimensions",
    "_resolve_linear_diagram_label_font_size",
]


