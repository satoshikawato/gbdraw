#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from gbdraw.features.placement import FeaturePlacementSlot
from gbdraw.render_stages import render_stage
from ...layout.linear_wrap import drawn_axis_length, refuse_feature_track_on_the_axis, refuse_wrapped_anchor_moves
from ...render.track_slot_metadata import (
    RecordAxisInput,
    attach_feature_placement_targets,
    attach_track_slot_geometry,
    linear_axis_geometry,
)
from ...layout.record_labels import coordinate_map_ticks, source_ruler_ticks

import copy
import math
from dataclasses import dataclass, replace
from functools import partial
import logging
from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Literal, Mapping, NamedTuple, Sequence, cast

from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pandas import DataFrame
from svgwrite import Drawing
from svgwrite.container import Group

from ...analysis.skew import skew_df
from ...analysis.depth import depth_df as build_depth_df
from ...analysis.depth_tracks import (
    DepthTrackData,
    DepthTrackSpec,
    build_depth_track_dataframes,
    depth_track_count,
    depth_track_heights,
    index_depth_track_row,
    representative_depth_tracks,
    sync_depth_track_legend_entries,
)
from ...auto_sizes import linear_tick_interval
from ...canvas import LinearCanvasConfigurator
from ...config.models import (
    GbdrawConfig,
    LinearRenderProfile,
)
from ...exceptions import ValidationError
from ...configurators import (
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
)
from ...configurators.bundle import DiagramConfigurators
from ...configurators.gc import _slot_skew_config
from ...configurators.legend import measure_legend
from ...core.text import calculate_svg_bbox_dimensions, drawn_text_box
from ...render.group_content import (
    GroupContent,
    GroupKey,
    GroupKind,
    build_group_content,
    place_content,
)
from ...render.groups.linear import (
    LengthBarGroup,
    linear_legend_content,
)
from ...render.groups.linear.definition import LinearDefinitionStyle
from ...render.groups.linear.row_windows import RowWindowHeadingGroup, RowWindowMarksGroup, RowWindowPlace
from ...render.groups.linear.length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
    format_linear_tick_label,
)
from ...io.comparisons import comparison_row_tally, read_comparison_resource
from ...legend.row_facts import LegendRowFacts, attach_legend_row_facts
from ...legend.sources import (
    ComparisonLegendSource,
    DepthLegendSource,
    LegendSlotSource,
    LegendSources,
    NumericLegendSource,
    depth_legend_sources,
    feature_legend_presence,
    legend_slot_sources,
)
from ...legend.table import (
    LegendRowEdits,
    LegendTable,
    PairwiseIdentityLegend,
    apply_legend_row_edits,
    pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
    sync_feature_density_legend_entries,
)
from ...layout.linear import (
    place_linear_definition,
    AxisGapResolution,
    CollisionBand,
    LinearFeatureLaneGeometry,
    LinearRecordRenderContext,
    LinearResolvedTrack,
    LinearTrackLayout,
    VerticalBand,
    measure_linear_feature_lanes,
    measure_linear_label_band,
    resolve_axis_gap,
)
from ...layout.composition import (
    CompositionInputs,
    LegendPlacement,
)
from ...layout.spatial import Aabb, union_aabbs
from ...layout.linear import union_vertical_bands
from ...layout.linear_multi_record import (
    LinearFeatureAttachment,
    LinearRecordMeasurement,
    LinearRecordPlacement,
    RecordKey,
    place_record_vertically,
    solve_linear_horizontal_layout,
    solve_linear_layout,
    comparison_frame_px,
    stable_record_keys,
)
from ...layout.linear_row_windows import (
    LinearRowWindows,
    RowWindowStyle,
    plan_row_windows,
    row_window_insets,
)
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.record_transform import bind_record_transforms
from ...layout.record_units import RecordUnit, record_units
from ...layout.record_placement import FeatureAnchorPlacement, resolve_record_row_positions
from ...layout.similarity_alignment import (
    AlignmentDecisionStatus,
    SimilarityAlignmentPlan,
)
from ...mode_profiles import ComparisonThresholds
from ...linear_comparison import (
    ComparisonRecordIdWarning,
    ComparisonTableTally,
    LinearComparison,
    ProteinComparisonFrames,
    project_search_frame_comparisons,
    search_frame_linear_comparisons,
    validate_linear_comparison_record_ids,
)
from ...layout.circular_depth_axis import ScalarAxisSource
from ...layout.scalar_axis import linear_scalar_axis_tick_font_size_px
from ...layout.scalar_axis import (
    LINEAR_SCALAR_AXIS_LABEL_GAP_PX,
    LINEAR_SCALAR_AXIS_STROKE_PX,
    LINEAR_SCALAR_AXIS_TICK_PX,
    ScalarAxisLabels,
    _depth_axis_bounds,
    depth_axis_labels,
    percent_axis_labels,
)
from ...labels.label_style import FeatureLabelStyle, label_leader_width_px
from ...labels.linear_types import LinearLabel
from ...labels.linear import (
    LINEAR_LABEL_FIT_ATTR,
    FeatureLabelRows,
    LinearLabelFrame,
    calculate_label_bounds,
    linear_label_fit,
    linear_label_on_its_wrap_row,
    linear_label_rows,
    linear_paint_order,
    with_feature_label_rows,
)
from ...tracks import (
    LinearTrackSlot,
    ScalarSpec,
    default_linear_track_slots,
    normalize_linear_track_slots_with_axis,
)
from ...tracks.parsing import (
    depth_slot_track_index,
    slot_dinucleotide,
    slot_legend_label,
    slots_draw_features,
    check_feature_density_slots,
    feature_density_slot_id,
    with_feature_density_slots,
)
from ...analysis.feature_density import FeatureDensityLayer
from ..feature_density import (
    feature_density_record_layers,
    feature_density_color,
    feature_density_histogram_source,
    feature_density_legend,
    slot_feature_density_layer,
)
from ...annotations import (
    AnnotationOptions,
    LabelSizing,
    ResolvedAnnotationBundle,
    ResolvedAnnotationTrack,
    annotation_track_params_from_mapping,
    effective_annotation_style,
    feature_underlay_anchor_slot_id,
    feature_underlay_slot_id,
    LinearAnnotationLabelSide,
    layout_annotation_track,
    linear_annotation_label_side,
    merge_feature_underlays,
    sync_annotation_legend_entries,
)
from ...annotations.planning import prepare_annotation_track_slots
from ...annotations.models import slots_in_paint_order
from ...render.drawers.linear.feature_density import FeatureDensityHeatmapGroup
from ...render.drawers.linear.annotations import (
    linear_annotation_label_extent,
    linear_annotation_label_room,
    linear_annotation_mark_past_lane,
    linear_annotation_track_content,
    linear_region_label_front_content,
)
from .annotation_labels import (
    AnnotationLabelsBySlot,
    LabelObstacle,
    PlacedAnnotationLabels,
    linear_annotation_label_inks,
    linear_annotation_slot_frame,
    linear_feature_glyph_boxes,
    linear_feature_label_obstacles,
    place_linear_annotation_labels,
)
from ...render.groups.linear import build_linear_feature_dom_index
from ...render.composition import CompositionPrimary, compose_drawing, place_legend, place_title

from .builders import (
    add_explicit_comparisons_on_linear_canvas,
    add_depth_group,
    add_feature_density_heatmap_group,
    add_gc_content_group,
    add_gc_skew_group,
    add_length_bar_on_linear_canvas,
    add_record_definition_group,
    add_record_group,
    add_row_window_heading_group,
    add_row_window_marks_group,
)
from .parts import (
    LinearAlignmentInputs,
    LinearComparisonInputs,
    LinearRecordInputs,
    LinearTrackInputs,
)
from .orthogroup_alignment import (
    build_orthogroup_label_eligibility,
    orthogroup_label_sets_for_record,
)
from .precalc import (
    FeatureDict,
    _precalculate_definition_metrics,
    _precalculate_feature_layers,
    _precalculate_label_dimensions,
    _resolve_linear_diagram_label_font_size,
)
from ...svg.ids import (
    definition_group_svg_id,
    instance_svg_id,
    record_group_svg_id,
    track_slot_svg_id,
)
from .track_slots import (
    LinearRecordVerticalPlan,
    LinearResolvedSlot,
    LinearSlotFootprint,
    resolve_linear_record_vertical_plan,
    resolve_linear_track_layout,
    stack_linear_wrap_rows,
)


logger = logging.getLogger(__name__)
_COMPARISON_ENDPOINT_GAP_PX = 4.0

if TYPE_CHECKING:
    from ...api.options import LinearMultiRecordOptions
    from ...api.requests import RecordPresentation
    from ...features.factory import FeatureBuildResult
    from ...layout.composition import CompositionPlan
    from ...layout.linear_multi_record import LinearHorizontalLayout, LinearLayoutPlan
    from ...labels.linear import LinearLabelFit
    from ...render.drawers.linear.annotations import LinearAnnotationLabel
    from ...render.groups.linear.feature_identity import LinearFeatureDomIndex
    from ...tracks import NormalizedLinearTrackSlot
    from .orthogroup_alignment import OrthogroupLabelEligibility


def _linear_feature_placement_slot(side: str, separate_strands: bool) -> FeaturePlacementSlot:
    # FeaturePlacementSlot.__post_init__ validates the side for linear mode.
    return FeaturePlacementSlot(
        "linear",
        cast(Literal["overlay", "above", "below"], side),
        separate_strands,
    )


def _drawn_axis_lengths(
    records: Sequence[SeqRecord], record_transforms: Sequence[RecordDisplayTransform] | None
) -> list[int]:
    """The bases on each record's drawn axis, its row length when it wraps."""

    return [
        drawn_axis_length(record, record_transforms[index] if record_transforms is not None else None)
        for index, record in enumerate(records)
    ]


def _prepare_linear_annotation_tracks(
    records: list[SeqRecord],
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    slots: list[LinearTrackSlot] | None,
    *,
    profile: LinearRenderProfile,
    record_depth_tracks: list[list[DepthTrackSpec]] | None,
    label_font_size: float,
    sequence_widths: Sequence[float],
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    axis_index: int | None = None,
) -> tuple[
    list[LinearTrackSlot] | None,
    ResolvedAnnotationBundle,
    dict[str, ResolvedAnnotationTrack],
    frozenset[str],
]:
    slots, bundle, auto_slot_ids = prepare_annotation_track_slots(
        annotations,
        records,
        slots,
        mode="linear",
        record_transforms=record_transforms,
        default_slots=lambda: default_linear_track_slots(
            show_features=True,
            show_depth=bool(record_depth_tracks),
            depth_track_count=max(1, depth_track_count(record_depth_tracks)),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            track_layout=profile.track_layout,
        ),
        slot_factory=LinearTrackSlot,
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, {}, frozenset()
    # The planner generates the slot list whenever annotations are present.
    assert slots is not None
    updated_slots, layouts = _layout_linear_annotation_tracks(
        records,
        slots,
        bundle,
        auto_slot_ids=auto_slot_ids,
        label_font_size=label_font_size,
        font_family=str(profile.config.objects.text.font_family),
        label_sides=_linear_annotation_label_sides(normalize_linear_track_slots_with_axis(slots, axis_index)),
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
    )
    return updated_slots, bundle, layouts, auto_slot_ids


def _linear_annotation_label_sides(slots: Sequence[NormalizedLinearTrackSlot]) -> dict[str, LinearAnnotationLabelSide]:
    """Where each drawn annotation slot puts its labels (``linear_annotation_label_side``), found by its anchor's slot
    id among the drawn features slots."""

    feature_sides = {slot.id: slot.side for slot in slots if slot.renderer == "features"}
    return {
        slot.id: linear_annotation_label_side(slot.side, annotation_track_params_from_mapping(slot.params), feature_sides)
        for slot in slots
        if slot.renderer == "annotations"
    }


def _layout_linear_annotation_tracks(
    records: list[SeqRecord],
    slots: list[LinearTrackSlot],
    bundle: ResolvedAnnotationBundle,
    *,
    auto_slot_ids: frozenset[str],
    label_font_size: float,
    font_family: str,
    label_sides: Mapping[str, LinearAnnotationLabelSide],
    sequence_widths: Sequence[float],
    record_transforms: Sequence[RecordDisplayTransform] | None,
) -> tuple[list[LinearTrackSlot], dict[str, ResolvedAnnotationTrack]]:
    """Resolve annotation lanes for the records' drawn axis widths (their placements'
    ``sequence_width``).

    ``label_font_size`` is the feature label font size of the diagram, the size of a region label whose style
    sets none; ``font_family`` measures the labels' letters. ``label_sides`` is where each drawn slot puts its
    labels (``linear_annotation_label_side``): the track keeps it (``ResolvedAnnotationTrack.labels_above``,
    ``labels_beyond_slot``) for the slot's footprint, its drawer and its label placement."""

    if len(sequence_widths) != len(records):
        raise ValueError("sequence_widths must contain one width for every record")
    record_lengths = dict(enumerate(_drawn_axis_lengths(records, record_transforms)))
    bp_per_px = {index: length / max(1.0, float(sequence_widths[index])) for index, length in record_lengths.items()}
    updated_slots: list[LinearTrackSlot] = []
    layouts: dict[str, ResolvedAnnotationTrack] = {}
    for slot in slots:
        if str(slot.renderer).strip().lower() != "annotations":
            updated_slots.append(slot)
            continue
        params = annotation_track_params_from_mapping(slot.params)
        available = (
            None
            if slot.id in auto_slot_ids
            else (slot.height.resolve(1.0) if slot.height is not None else None)
        )
        # A disabled slot is not drawn.
        label_side = label_sides.get(slot.id, LinearAnnotationLabelSide(above=False))
        layout = layout_annotation_track(
            slot.id,
            params.set_id,
            bundle.annotations,
            record_lengths=record_lengths,
            params=params,
            available_extent_px=available,
            labels=LabelSizing(
                font_size=label_font_size,
                extent=partial(linear_annotation_label_extent, font_family=font_family, bp_per_px=bp_per_px),
                reach=linear_annotation_mark_past_lane,
                bp_per_px=bp_per_px,
                labels_above=label_side.above,
            ),
            join_origin=False,
        )
        layouts[slot.id] = replace(layout, labels_beyond_slot=label_side.beyond_slot)
        updated_slots.append(
            replace(slot, height=ScalarSpec(layout.required_extent_px, "px"))
            if slot.id in auto_slot_ids
            else slot
        )
    return updated_slots, layouts


def _add_linear_feature_underlays(
    records: list[SeqRecord],
    feature_layers: list,
    slots: list[LinearTrackSlot],
    bundle: ResolvedAnnotationBundle,
    layouts: dict[str, ResolvedAnnotationTrack],
    auto_slot_ids: frozenset[str],
    *,
    label_font_size: float,
    font_family: str,
    sequence_widths: Sequence[float],
    record_transforms: Sequence[RecordDisplayTransform] | None,
) -> tuple[
    list[LinearTrackSlot],
    ResolvedAnnotationBundle,
    dict[str, ResolvedAnnotationTrack],
    frozenset[str],
    int,
]:
    merged, set_id = merge_feature_underlays(
        bundle,
        [result.underlay_features for result in feature_layers],
        records,
        mode="linear",
    )
    anchor_id = (
        feature_underlay_anchor_slot_id(slots, (f for result in feature_layers for f in result.underlay_features))
        if set_id is not None
        else None
    )
    if anchor_id is None:
        return slots, bundle, layouts, auto_slot_ids, 0
    anchor = next(slot for slot in slots if str(slot.id) == anchor_id)
    slot_id = feature_underlay_slot_id(slots)
    underlay_slot = LinearTrackSlot(
        id=slot_id,
        renderer="annotations",
        side="overlay",
        z=int(anchor.z) - 1,
        params={
            "set_id": set_id,
            "marks": ("highlight",),
            "anchor_slot": anchor_id,
            "layer": "underlay",
            "cover_anchor": True,
            "padding_px": 0.0,
            "show_labels": False,
        },
    )
    # The underlay set is new (collision-safe id), so the slots laid out before keep their layouts.
    sized_underlay, underlay_layouts = _layout_linear_annotation_tracks(
        records,
        [underlay_slot],
        merged,
        auto_slot_ids=frozenset((slot_id,)),
        label_font_size=label_font_size,
        font_family=font_family,
        label_sides={
            slot_id: linear_annotation_label_side(
                str(underlay_slot.side), annotation_track_params_from_mapping(underlay_slot.params), {}
            )
        },
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
    )
    return (
        [*sized_underlay, *slots],
        merged,
        {**underlay_layouts, **layouts},
        frozenset((*auto_slot_ids, slot_id)),
        len(sized_underlay),
    )


def _apply_depth_track_heights_to_linear_slots(
    slots: list,
    record_depth_tracks: list[list[DepthTrackSpec]] | None,
) -> list:
    depth_heights = depth_track_heights(record_depth_tracks)
    if not depth_heights:
        return slots
    out = []
    for slot in slots:
        if slot.renderer != "depth" or slot.height is not None:
            out.append(slot)
            continue
        track_index = _depth_slot_track_index(slot)
        track_height = depth_heights[track_index] if 0 <= track_index < len(depth_heights) else None
        if track_height is not None:
            out.append(replace(slot, height=ScalarSpec(float(track_height), "px")))
        else:
            out.append(slot)
    return out


def _depth_slot_track_index(slot) -> int:
    return depth_slot_track_index(str(getattr(slot, "id", "")), getattr(slot, "params", {}) or {})


def _clone_gc_config_with_dinucleotide(gc_config: GcContentConfigurator, dinucleotide: str) -> GcContentConfigurator:
    cloned = copy.copy(gc_config)
    cloned.dinucleotide = str(dinucleotide).upper()
    return cloned










def _sync_legend_table_for_linear_slots(
    legend_table: LegendTable,
    *,
    slots: Sequence[LegendSlotSource] | None,
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    depth_tracks: Sequence[DepthLegendSource] | None,
) -> LegendTable:
    """Replace singleton numeric legends with slot-authoritative entries.

    Each content and skew slot draws its own rows, named by its
    ``legend_label`` (alias ``label``), as in Circular (OV-846, Owner LP-D67).
    The rows keep the Linear order: content rows, then Depth rows, then skew rows.
    """
    if slots is None:
        return legend_table

    default_nt = str(gc.dinucleotide).upper()
    for kind in ("content", "skew"):
        for key in (f"{default_nt} {kind}", f"{default_nt} {kind} (+)", f"{default_nt} {kind} (-)"):
            legend_table.discard_default_track(key)

    def add_slot_rows(table: LegendTable, renderer: str) -> None:
        for slot in slots:
            if slot.renderer != renderer:
                continue
            nt = slot_dinucleotide(slot.params, default_nt)
            if renderer == "dinucleotide_content":
                table.add_numeric_source(slot_legend_label(slot, f"{nt} content"), gc)
            else:
                table.add_numeric_source(slot_legend_label(slot, f"{nt} skew"), skew.for_skew_slot(slot, nt))

    add_slot_rows(legend_table, "dinucleotide_content")
    out = sync_depth_track_legend_entries(legend_table, depth_tracks, slots=slots)
    add_slot_rows(out, "dinucleotide_skew")
    return out


def _serialize_linear_track_slot_geometry(
    *,
    record_axes: Sequence[RecordAxisInput],
    units: Sequence[RecordUnit],
    layout: LinearTrackLayout,
    record_plans: list[LinearRecordVerticalPlan],
    record_offsets: list[float],
    record_placements: dict[int, LinearRecordPlacement],
    composition_offset_x: float,
    horizontal_offset_x: float,
    record_collision_bands: list[tuple[CollisionBand, ...]],
    boundary_gap_resolutions: list[AxisGapResolution],
) -> dict[str, Any]:
    def band_payload(band: VerticalBand | None, *, axis_y: float) -> dict[str, float] | None:
        if band is None:
            return None
        return {
            "topPx": float(band.top_y),
            "bottomPx": float(band.bottom_y),
            "absoluteTopPx": axis_y + float(band.top_y),
            "absoluteBottomPx": axis_y + float(band.bottom_y),
        }

    def collision_payload(
        band: CollisionBand,
        *,
        axis_y: float,
    ) -> dict[str, float | str]:
        return {
            "kind": str(band.kind),
            "xStartPx": float(band.x_start),
            "xEndPx": float(band.x_end),
            "topPx": float(band.top_y),
            "bottomPx": float(band.bottom_y),
            "absoluteTopPx": axis_y + float(band.top_y),
            "absoluteBottomPx": axis_y + float(band.bottom_y),
        }

    base_by_id = {slot.id: slot for slot in layout.slots}
    records_payload: list[dict[str, Any]] = []
    for record_index, (drawn, unit) in enumerate(zip(record_axes, units, strict=True)):
        axis_y = float(record_offsets[record_index]) if record_index < len(record_offsets) else 0.0
        placement = record_placements[record_index]
        plan = record_plans[record_index]
        slots_payload: list[dict[str, Any]] = []
        for slot in plan.slots:
            base_slot = base_by_id[slot.id]
            slots_payload.append(
                {
                    "slotIndex": int(slot.slot_index),
                    "slotId": str(slot.id),
                    "renderer": str(slot.renderer),
                    "side": str(slot.side),
                    "heightPx": float(slot.height),
                    "spacingAfterPx": float(slot.spacing_after_px),
                    "baseYOffsetPx": float(base_slot.y_offset),
                    "resolvedOriginPx": float(slot.origin_y),
                    "finalYOffsetPx": axis_y + float(slot.origin_y),
                    "dataAvailable": bool(slot.data_available),
                    "paintBand": band_payload(slot.paint_band, axis_y=axis_y),
                    "reserveBand": band_payload(slot.reserve_band, axis_y=axis_y),
                    "source": "resolved",
                }
            )
        record_id = str(getattr(drawn.record, "id", "") or "")
        row_bands = [(row.axis_offset_y, row.plan.feature_lanes_band.translate(row.axis_offset_y)) for row in plan.rows]
        # A window of a record with several names its place among the record's windows.
        window_entry: dict[str, int] = {} if unit.window_index is None else {"windowIndex": unit.window_index}
        records_payload.append(
            {
                "recordIndex": int(record_index),
                "recordId": record_id,
                "recordLabel": record_id,
                # A window of a record names its record (D-01); any other names itself.
                "recordKey": unit.record_key,
                **window_entry,
                "axis": linear_axis_geometry(
                    drawn,
                    rows=[(axis, band.top_y, band.bottom_y) for axis, band in row_bands],
                    row_length_bp=drawn_axis_length(drawn.record, drawn.transform),
                    px_per_bp=float(placement.px_per_bp),
                ),
                "axisYpx": axis_y,
                "axisXpx": (
                    float(placement.x)
                    + float(horizontal_offset_x)
                    + float(composition_offset_x)
                ),
                "sequenceWidthPx": float(placement.sequence_width),
                "recordBodyBand": band_payload(plan.record_body_band, axis_y=axis_y),
                "comparisonExclusionBand": band_payload(
                    plan.comparison_exclusion_band,
                    axis_y=axis_y,
                ),
                "canvasBand": band_payload(plan.canvas_band, axis_y=axis_y),
                "collisionBands": [
                    collision_payload(band, axis_y=axis_y)
                    for band in (
                        record_collision_bands[record_index]
                        if record_index < len(record_collision_bands)
                        else ()
                    )
                ],
                "slots": slots_payload,
                # One record: the rows of a wrapped one are only where its slots repeat below its axis.
                **(
                    {
                        "wrap": {
                            "rowOffsets": list(plan.row_offsets),
                            "rowCanvasBands": [
                                band_payload(row.plan.canvas_band.translate(row.axis_offset_y), axis_y=axis_y)
                                for row in plan.wrap_rows
                            ],
                        }
                    }
                    if plan.wrap_rows
                    else {}
                ),
            }
        )
    return {
        "schema": 2,
        "mode": "linear",
        "source": "resolved",
        "records": records_payload,
        "axisGapConstraints": [
            {
                "boundaryRow": int(boundary),
                "nextRow": int(boundary + 1),
                "axisGapPx": float(resolution.axis_gap),
                "clearGapPx": float(resolution.clear_gap),
                "currentKind": (
                    str(resolution.current_band.kind)
                    if resolution.current_band is not None
                    else None
                ),
                "nextKind": (
                    str(resolution.next_band.kind)
                    if resolution.next_band is not None
                    else None
                ),
            }
            for boundary, resolution in enumerate(boundary_gap_resolutions)
        ],
    }


def _content_bounds(content: GroupContent) -> Aabb:
    """The measured bounds of a content whose group class measures them."""

    if content.local_bounds is None:  # pragma: no cover - the group classes measure them
        raise ValueError("this group content has no measured bounds")
    return content.local_bounds


def _read_comparison_files(
    blast_files: Sequence[str] | None,
    thresholds: ComparisonThresholds,
) -> tuple[DataFrame, ...]:
    """Read each comparison file of the request and log the thresholds that will filter it.

    A file the request staged as a resource is parsed once per resource and
    gives its table that resource's identity, so the ``comparison_rows`` stage
    that filters it can be reused. This runs outside that stage, so every
    request logs what it reads.
    """

    if not blast_files:
        return ()
    logger.info(
        "BLAST output visualization settings: e-value threshold: {}; bitscore threshold: {}; "
        "identity threshold: {}; alignment length threshold: {}".format(
            thresholds.evalue, thresholds.bitscore, thresholds.identity, thresholds.alignment_length
        )
    )
    logger.info("Loading comparison file(s)...")
    tables: list[DataFrame] = []
    for comparison_file in blast_files:
        logger.info("Loading {}".format(comparison_file))
        tables.append(read_comparison_resource(comparison_file))
    logger.info("... finished loading comparison file(s)")
    return tuple(tables)


@dataclass(frozen=True)
class LinearComparisonRows:
    """The comparisons a Linear diagram draws and the record-ID warnings they raise."""

    comparisons: tuple[LinearComparison, ...]
    record_id_warnings: tuple[ComparisonRecordIdWarning, ...]
    # One per table the thresholds filtered, in the order the tables were given.
    tallies: tuple[ComparisonTableTally, ...] = ()


def linear_comparison_rows(
    blast_tables: tuple[DataFrame, ...],
    protein_comparisons: ProteinComparisonFrames | None,
    linear_comparisons: tuple[LinearComparison, ...],
    thresholds: ComparisonThresholds,
    records: tuple[SeqRecord, ...],
    rows_by_record: tuple[int, ...],
    multi_record_enabled: bool,
    blast_names: tuple[str, ...] = (),
) -> LinearComparisonRows:
    """Filter, merge, validate and project the comparisons (render stage ``comparison_rows``).

    :func:`search_frame_linear_comparisons` filters, merges and validates them;
    this stage projects them onto the records. The ordered tables (comparison
    files, then protein comparisons) pair record ``i`` with ``i + 1``;
    explicit comparisons (given, then protein pairs) name their records.
    ``blast_names`` names the comparison files.
    """

    proteins = protein_comparisons or ProteinComparisonFrames()
    tallies = (
        *(
            ComparisonTableTally(
                blast_names[index] if index < len(blast_names) else None,
                index, index + 1, comparison_row_tally(table, thresholds),
            )
            for index, table in enumerate(blast_tables)
        ),
        *(
            ComparisonTableTally(None, index, index + 1, comparison_row_tally(table, thresholds), protein_matches=True)
            for index, table in enumerate(proteins.ordered or ())
        ),
        *(
            ComparisonTableTally(
                item.source_name, item.query_record_index, item.subject_record_index,
                comparison_row_tally(item.matches, thresholds), protein_matches=item.protein_matches,
            )
            for item in (*linear_comparisons, *proteins.pairs)
        ),
    )
    normalized = search_frame_linear_comparisons(
        blast_tables, protein_comparisons, linear_comparisons, thresholds,
        rows_by_record, multi_record_enabled,
    )
    warnings = validate_linear_comparison_record_ids(normalized, records)
    return LinearComparisonRows(
        comparisons=tuple(project_search_frame_comparisons(normalized, records)),
        record_id_warnings=tuple(warnings),
        tallies=tallies,
    )


def _gc_config_matching_linear_depth_axis_font_size(
    *,
    gc_config: GcContentConfigurator,
    depth_config: DepthConfigurator | None,
    depth_height: float,
    depth_enabled: bool,
) -> GcContentConfigurator:
    if (
        not depth_enabled
        or depth_config is None
        or str(getattr(gc_config, "mode", "deviation")).strip().lower() != "percent"
        or getattr(gc_config, "tick_font_size", None) is not None
        or not bool(getattr(depth_config, "show_axis", True))
        or not bool(getattr(depth_config, "show_ticks", True))
    ):
        return gc_config

    cloned = copy.copy(gc_config)
    cloned.tick_font_size = linear_scalar_axis_tick_font_size_px(depth_config, depth_height)
    return cloned


# The GC content configuration a record's slot (by its params) is drawn with, given the record's depth
# tracks: the drawing, the slot's vertical footprint and the record's left axis room read the same one.
LinearGcConfigOf = Callable[[Mapping[str, object], Mapping[int, DepthTrackData]], GcContentConfigurator]


def _axis_ruler_extents(
    canvas_config: LinearCanvasConfigurator,
    render_context: LinearRecordRenderContext,
) -> tuple[float, float]:
    """
    Return (height_above_axis, height_below_axis) required by axis-based ruler labels.
    """
    if not render_context.axis_ruler_enabled:
        return 0.0, 0.0

    cfg = render_context.profile.config
    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param(canvas_config.length_param)
    label_height = drawn_text_box("0", cfg.objects.text.font_family, font_size, cfg.objects.scale.font_weight)[1]
    protrusion = max(
        0.5 * float(cfg.objects.scale.stroke_width),
        float(RULER_TICK_LENGTH),
        float(RULER_LABEL_OFFSET) + float(label_height),
    )
    if render_context.feature_track_layout == "above":
        return 0.0, protrusion
    if render_context.feature_track_layout == "below":
        return protrusion, 0.0
    return 0.0, 0.0


def _linear_axis_band(
    canvas_config: LinearCanvasConfigurator,
    render_context: LinearRecordRenderContext,
) -> VerticalBand:
    cfg = render_context.profile.config
    axis_stroke_width = cfg.objects.axis.linear.stroke_width.for_length_param(
        canvas_config.length_param
    )
    half_stroke = 0.5 * max(0.0, float(axis_stroke_width))
    ruler_above, ruler_below = _axis_ruler_extents(
        canvas_config,
        render_context,
    )
    return VerticalBand(
        -max(half_stroke, float(ruler_above)),
        max(half_stroke, float(ruler_below)),
    )


def _feature_reserve_band(
    paint_band: VerticalBand,
    *,
    side: str,
    minimum_height: float,
) -> VerticalBand:
    minimum = max(0.0, float(minimum_height))
    if paint_band.height >= minimum:
        return paint_band
    if side == "above":
        return VerticalBand(paint_band.bottom_y - minimum, paint_band.bottom_y)
    if side == "below":
        return VerticalBand(paint_band.top_y, paint_band.top_y + minimum)
    center_y = 0.5 * (paint_band.top_y + paint_band.bottom_y)
    return VerticalBand(center_y - 0.5 * minimum, center_y + 0.5 * minimum)


def _slot_body_band(slot: LinearResolvedTrack) -> VerticalBand:
    return VerticalBand(-float(slot.top_extent), float(slot.bottom_extent))


def _axis_text_overhang(axis_config: ScalarAxisSource, track_height: float) -> float:
    """How far a Linear scalar axis paints above and below its track: half its stroke, or half a label line when
    numbered tick labels sit centred on the top and bottom ticks (the range label sits centred on the track)."""

    if not bool(axis_config.show_axis):
        return 0.0
    axis_stroke = 0.5 * LINEAR_SCALAR_AXIS_STROKE_PX
    if not bool(axis_config.show_ticks) or axis_config.axis_labels != "numbered":
        return axis_stroke
    return max(
        axis_stroke,
        0.5 * linear_scalar_axis_tick_font_size_px(axis_config, track_height),
    )


def _linear_slot_footprints_for_record(
    *,
    record_index: int,
    layout: LinearTrackLayout,
    feature_dict: FeatureDict,
    feature_geometry: LinearFeatureLaneGeometry,
    labels: list[LinearLabel],
    canvas_config: LinearCanvasConfigurator,
    gc_config_of: LinearGcConfigOf,
    skew_config: GcSkewConfigurator,
    depth_config_by_index: dict[int, DepthConfigurator],
    depth_by_index: Mapping[int, DepthTrackData],
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    wrap_row: int | None = None,
) -> tuple[dict[str, LinearSlotFootprint], float]:
    """Each slot's footprint on the record, or on one row of a wrapped record (``wrap_row``, with that row's
    ``labels``)."""

    footprints: dict[str, LinearSlotFootprint] = {}
    font_family = str(canvas_config._cfg.objects.text.font_family)
    label_band = measure_linear_label_band(labels)

    for slot in layout.slots:
        body_band = _slot_body_band(slot)
        if slot.renderer == "features":
            paint_band = feature_geometry.occupied_band
            if label_band is not None:
                paint_band = paint_band.union(label_band)
            reserve_band = _feature_reserve_band(
                paint_band,
                side=slot.side,
                minimum_height=float(slot.height),
            )
            lanes = feature_geometry.occupied_band
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=paint_band if feature_dict or labels else None,
                reserve_band=reserve_band,
                data_available=bool(feature_dict or labels),
                # Nominal lanes under underlay features are a body without lane features.
                body_band=lanes if lanes.height > 0.0 else None,
                strand_bands=feature_geometry.strand_lane_bands(),
            )
            continue

        if slot.renderer == "spacer":
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=None,
                reserve_band=body_band,
                data_available=False,
            )
            continue

        if slot.renderer == "depth":
            track_index = _depth_slot_track_index(slot)
            axis_config = depth_config_by_index.get(track_index)
            overhang = (
                _axis_text_overhang(axis_config, float(slot.height))
                if axis_config is not None
                else 0.0
            )
            stroke_overhang = 0.5 * max(
                0.0,
                float(getattr(axis_config, "stroke_width", 0.0)),
            )
            expanded = body_band.expand(max(overhang, stroke_overhang))
            available = track_index in depth_by_index
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=expanded if available else None,
                reserve_band=expanded if available else VerticalBand(0.0, 0.0),
                data_available=available,
            )
            continue

        if slot.renderer == "dinucleotide_content":
            slot_gc_config = gc_config_of(slot.params, depth_by_index)
            overhang = max(
                0.5 * max(0.0, float(slot_gc_config.stroke_width)),
                0.5 * max(0.0, float(slot_gc_config.percent_border_width)),
            )
            if str(slot_gc_config.mode) == "percent":
                overhang = max(overhang, _axis_text_overhang(slot_gc_config, float(slot.height)))
            expanded = body_band.expand(overhang)
            footprints[slot.id] = LinearSlotFootprint(expanded, expanded)
            continue

        if slot.renderer == "dinucleotide_skew":
            overhang = 0.5 * max(0.0, float(skew_config.stroke_width))
            expanded = body_band.expand(overhang)
            footprints[slot.id] = LinearSlotFootprint(expanded, expanded)
            continue

        if slot.renderer == "annotations":
            track = annotation_track_layouts.get(slot.id)
            footprints[slot.id] = _linear_annotation_slot_footprint(
                slot,
                record_index=record_index,
                annotation_layout=track,
                font_family=font_family,
                wrap_row=wrap_row,
                beyond_band=_linear_labels_beyond_band(layout, slot, track, feature_geometry),
            )
            continue

        footprints[slot.id] = LinearSlotFootprint(body_band, body_band)
    feature_center_y = 0.5 * (
        feature_geometry.occupied_band.top_y
        + feature_geometry.occupied_band.bottom_y
    )
    return footprints, feature_center_y


def _linear_labels_beyond_band(
    layout: LinearTrackLayout,
    slot: LinearResolvedTrack,
    track: ResolvedAnnotationTrack | None,
    feature_geometry: LinearFeatureLaneGeometry,
) -> VerticalBand | None:
    """The record's feature lanes (``LinearFeatureLaneGeometry.occupied_band``) where an overlay's labels stand
    beyond them (``ResolvedAnnotationTrack.labels_beyond_slot``), in the overlay's own frame: it keeps its origin
    from that features slot. None where its labels stand beside their marks or the record has no lanes."""

    anchor_id = track.labels_beyond_slot if track is not None else None
    anchor = next((item for item in layout.slots if item.id == anchor_id), None)
    lanes = feature_geometry.occupied_band
    if anchor is None or lanes.height <= 0.0:
        return None
    return lanes.translate(float(anchor.y_offset) - float(slot.y_offset))


def _linear_annotation_slot_footprint(
    slot: LinearResolvedTrack,
    *,
    record_index: int,
    annotation_layout: ResolvedAnnotationTrack | None,
    font_family: str,
    span_height: float | None = None,
    wrap_row: int | None = None,
    beyond_band: VerticalBand | None = None,
) -> LinearSlotFootprint:
    """An annotation slot's own band with its labels' room on their side (on one row of a wrapped record:
    that row's labels, and paint only where the row holds a mark).

    ``span_height`` is the stretched height of an overlay the record's plan
    spans over slots (``spans_slots``): its labels then keep clear of every
    slot it covers, measured on that height; otherwise on the slot's own.
    ``beyond_band`` is the lanes its labels stand beyond (``_linear_labels_beyond_band``), which the footprint
    keeps for its drawer (``LinearSlotFootprint.label_beyond_band``).
    """

    body_band = _slot_body_band(slot)
    params = annotation_track_params_from_mapping(slot.params)
    label_overhang = (
        linear_annotation_label_room(
            annotation_layout,
            record_index=record_index,
            covers_slots=span_height is not None,
            height_px=float(body_band.height if span_height is None else span_height),
            font_family=font_family,
            params=params,
            wrap_row=wrap_row,
            beyond_band=beyond_band,
        )
        if annotation_layout is not None
        else 0.0
    )
    if annotation_layout is not None and annotation_layout.labels_above:
        expanded = body_band.expand(label_overhang, 0.0)
    else:
        expanded = body_band.expand(0.0, label_overhang)
    placed = (
        [item for item in annotation_layout.placements if item.annotation.record_index == record_index]
        if annotation_layout is not None
        else []
    )
    available = bool(placed)
    painted = available and (wrap_row is None or any(wrap_row in item.annotation.geometry_rows for item in placed))
    return LinearSlotFootprint(
        paint_band=expanded if painted else None,
        reserve_band=expanded,
        data_available=available,
        label_beyond_band=beyond_band,
    )


def _resolved_linear_depth_heights(layout: LinearTrackLayout) -> dict[int, float]:
    heights: dict[int, float] = {}
    for slot in layout.slots:
        if slot.renderer != "depth":
            continue
        track_index = _depth_slot_track_index(slot)
        heights.setdefault(track_index, float(slot.height))
    return heights


@dataclass(frozen=True)
class _LinearRecordDefinitionGeometry:
    """Record-axis-local bands used by definition measurement and rendering."""

    local_band: VerticalBand | None = None
    row_band: VerticalBand | None = None
    local_width: float = 0.0
    row_width: float = 0.0

    @property
    def local_center_y(self) -> float | None:
        if self.local_band is None:
            return None
        return 0.5 * (self.local_band.top_y + self.local_band.bottom_y)

    @property
    def row_center_y(self) -> float | None:
        if self.row_band is None:
            return None
        return 0.5 * (self.row_band.top_y + self.row_band.bottom_y)


def _centered_vertical_band(center_y: float, height: float) -> VerticalBand | None:
    resolved_height = max(0.0, float(height))
    if resolved_height <= 0.0:
        return None
    half_height = 0.5 * resolved_height
    return VerticalBand(float(center_y) - half_height, float(center_y) + half_height)


def _definition_geometry_beside_features(
    feature_center_y: float,
    *,
    local_height: float,
    local_width: float,
    row_height: float,
    row_width: float,
    multi_record_enabled: bool,
    split_row_definitions: bool,
    row_top: float | None = None,
) -> _LinearRecordDefinitionGeometry:
    """The definition text drawn beside a record's features, centered on them: a record alone in its
    row has its own text there, and the leading record of a shared row the row's heading. Its place does
    not depend on the feature labels, so the labels can keep clear of it (OV-1395); the record-local
    text of a shared row is a header above the record's canvas band, added by the vertical plan. A
    windowed row's heading (``row_top`` given) has its name line on the features' center and its list
    below it."""

    local_band = None if multi_record_enabled else _centered_vertical_band(feature_center_y, local_height)
    row_band = (
        (
            _centered_vertical_band(feature_center_y, row_height)
            if row_top is None
            else VerticalBand(feature_center_y + row_top, feature_center_y + row_top + max(0.0, float(row_height)))
        )
        if multi_record_enabled and split_row_definitions
        else None
    )
    return _LinearRecordDefinitionGeometry(
        local_band=local_band,
        row_band=row_band,
        local_width=float(local_width),
        row_width=float(row_width),
    )


_ROW_ELIGIBLE_LINE_KINDS: tuple[str, ...] = ("name", "subtitle")
_RECORD_LOCAL_LINE_KINDS = frozenset({"replicon", "accession", "length"})


def _record_definition_text(presentation: RecordPresentation | None, kind: str) -> str:
    """Return a record's explicit label ("name") or subtitle, or an empty string."""
    if presentation is None:
        return ""
    return str((presentation.label if kind == "name" else presentation.subtitle) or "").strip()


def _row_wide_definition_line_kinds(
    presentations: Sequence[RecordPresentation | None],
    row_record_indices: list[int],
    *,
    leading_index: int,
) -> frozenset[str]:
    """Return the definition line kinds that describe a whole row rather than a record.

    A kind qualifies when the record leading the row names it and no other record
    of the row contradicts that name. An empty value means the record has nothing
    of its own to say, which is how a records table writes a row identity once on
    its leading record; a different non-empty value makes the kind record-local, so
    a per-record replicon name is never promoted to a row heading.

    The leading record anchors the comparison because the row part is drawn from
    that record, so a name that only a following record carries could not be
    rendered beside the row at all.

    The kinds are tested in stacking order and the scan stops at the first one
    that does not qualify: a subtitle belongs under its own label, so it must not
    be left alone beside the row once the label has moved above the records.
    """

    row_wide: set[str] = set()
    for kind in _ROW_ELIGIBLE_LINE_KINDS:
        leading_text = _record_definition_text(presentations[leading_index], kind)
        if not leading_text:
            break
        if any(
            text and text != leading_text
            for text in (
                _record_definition_text(presentations[index], kind)
                for index in row_record_indices
            )
        ):
            break
        row_wide.add(kind)
    return frozenset(row_wide)


def _split_definition_line_kinds(
    presentations: Sequence[RecordPresentation | None],
    *,
    rows_by_record: tuple[int, ...],
    row_leading_indices: set[int],
) -> tuple[list[frozenset[str]], list[frozenset[str]]]:
    """Assign definition line kinds to the row part and local part of each record.

    The row part carries only what describes the whole row: a label or subtitle
    that no record of the row contradicts, which is what a Web file-level default
    applied to each record of a source, or a records table naming a row once on
    its leading record, produces. Anything that varies within the row stays above
    its own record, including for the record that leads the row.

    Both lists are explicit, so the caller measures exactly the lines it draws.
    """

    row_record_indices: dict[int, list[int]] = {}
    for index in range(len(presentations)):
        row_record_indices.setdefault(rows_by_record[index], []).append(index)

    leading_index_by_row = {
        rows_by_record[index]: index for index in row_leading_indices
    }
    row_wide_by_row = {
        row: _row_wide_definition_line_kinds(
            presentations,
            indices,
            leading_index=leading_index_by_row[row],
        )
        for row, indices in row_record_indices.items()
    }

    local_kinds: list[frozenset[str]] = []
    row_kinds: list[frozenset[str]] = []
    for index in range(len(presentations)):
        row_wide = row_wide_by_row[rows_by_record[index]]
        local_kinds.append(
            _RECORD_LOCAL_LINE_KINDS | (frozenset(_ROW_ELIGIBLE_LINE_KINDS) - row_wide)
        )
        row_kinds.append(row_wide if index in row_leading_indices else frozenset())
    return local_kinds, row_kinds


def _linear_record_vertical_offset(
    plans: list[LinearRecordVerticalPlan],
    rows_by_record: tuple[int, ...],
    canvas_config: LinearCanvasConfigurator,
) -> float:
    """Place the first row below the top edge using final measured canvas bands."""

    first_row_top_extent = max(
        (
            plan.canvas_top_extent
            for index, plan in enumerate(plans)
            if index < len(rows_by_record) and rows_by_record[index] == 0
        ),
        default=0.0,
    )
    return max(
        float(canvas_config.original_vertical_offset) + float(canvas_config.cds_padding),
        first_row_top_extent + float(canvas_config.vertical_padding),
    )


def _build_linear_record_vertical_plans(
    *,
    records: list[SeqRecord],
    layout: LinearTrackLayout,
    feature_dicts: list[FeatureDict],
    feature_geometries: list[LinearFeatureLaneGeometry],
    labels_by_record: list[list[LinearLabel]],
    canvas_config: LinearCanvasConfigurator,
    gc_config_of: LinearGcConfigOf,
    skew_config: GcSkewConfigurator,
    record_depth_by_index: list[dict[int, DepthTrackData]],
    representative_depth_configs: dict[int, DepthConfigurator],
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    definitions: _LinearDefinitionMetrics,
    multi_record_enabled: bool,
    split_row_definitions: bool,
    render_context: LinearRecordRenderContext,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
) -> tuple[list[LinearRecordVerticalPlan], list[_LinearRecordDefinitionGeometry]]:
    axis_band = _linear_axis_band(
        canvas_config,
        render_context,
    )
    plans: list[LinearRecordVerticalPlan] = []
    definition_geometries: list[_LinearRecordDefinitionGeometry] = []
    feature_slot_id = next(
        (slot.id for slot in layout.slots if slot.renderer == "features"),
        None,
    )

    def row_plan(record_index: int, wrap_row: int | None) -> tuple[LinearRecordVerticalPlan, float]:
        """The packed slots of a record, or of one row of a wrapped record, and its feature lanes' center."""

        depth_by_index = (
            record_depth_by_index[record_index]
            if record_index < len(record_depth_by_index)
            else {}
        )
        depth_configs = dict(representative_depth_configs)
        depth_configs.update(
            {index: track.config for index, track in depth_by_index.items()}
        )
        labels = labels_by_record[record_index]
        footprints, feature_center_y = _linear_slot_footprints_for_record(
            record_index=record_index,
            layout=layout,
            feature_dict=feature_dicts[record_index],
            feature_geometry=feature_geometries[record_index],
            labels=labels if wrap_row is None else [label for label in labels if label.get("wrap_row") == wrap_row],
            canvas_config=canvas_config,
            gc_config_of=gc_config_of,
            skew_config=skew_config,
            depth_config_by_index=depth_configs,
            depth_by_index=depth_by_index,
            annotation_track_layouts=annotation_track_layouts,
            wrap_row=wrap_row,
        )
        plan = resolve_linear_record_vertical_plan(
            layout,
            axis_band=axis_band,
            footprints=footprints,
        )
        # The plan decides which overlays span slots on this record (OV-1158). Their labels keep
        # clear of the slots they cover: measure their room again on the stretched height. A
        # single span's height does not depend on its room, which lies beyond the slots it covers.
        # Where one span covers a slot that holds another span's label room, this pass changes
        # that height and the room stays measured on the first one: only where the mark's paint
        # passes the span's edge does the far-side gap then differ from label_offset, by at most
        # 0.08 (band) or 0.5 (bracket, line) of the height change per lane (review 1177).
        span_heights = {slot.id: float(slot.height) for slot in plan.slots if slot.spans_slots}
        if span_heights:
            font_family = str(canvas_config._cfg.objects.text.font_family)
            footprints = {
                **footprints,
                **{
                    slot.id: _linear_annotation_slot_footprint(
                        slot,
                        record_index=record_index,
                        annotation_layout=annotation_track_layouts.get(slot.id),
                        font_family=font_family,
                        span_height=span_heights[slot.id],
                        wrap_row=wrap_row,
                    )
                    for slot in layout.slots
                    if slot.id in span_heights
                },
            }
            plan = resolve_linear_record_vertical_plan(
                layout,
                axis_band=axis_band,
                footprints=footprints,
            )
        return plan, feature_center_y

    for record_index, _record in enumerate(records):
        bound = record_transforms[record_index] if record_transforms is not None else None
        transform = bound if bound is not None and bound.wrap_length is not None else None
        plan, feature_center_y = row_plan(record_index, 0 if transform is not None else None)
        feature_origin_y = (
            plan.slot_by_id(feature_slot_id).origin_y
            if feature_slot_id is not None
            else 0.0
        )
        feature_center_y += feature_origin_y
        canvas_band = plan.canvas_band
        # A wrapped record's name and definition stand beside its first row, once.
        local_definition_height = float(definitions.definition_heights[record_index])
        beside = _definition_geometry_beside_features(
            feature_center_y,
            local_height=local_definition_height,
            local_width=float(definitions.definition_widths[record_index]),
            row_height=float(definitions.row_definition_heights[record_index]),
            row_width=float(definitions.row_definition_widths[record_index]),
            multi_record_enabled=multi_record_enabled,
            split_row_definitions=split_row_definitions,
            row_top=definitions.row_definition_tops[record_index],
        )
        if beside.row_band is not None:
            canvas_band = canvas_band.union(beside.row_band)
        local_band = beside.local_band
        if multi_record_enabled and local_definition_height > 0.0:
            header_bottom_y = (
                canvas_band.top_y - float(canvas_config.vertical_padding)
            )
            local_band = VerticalBand(
                header_bottom_y - local_definition_height,
                header_bottom_y,
            )

        if local_band is not None:
            canvas_band = canvas_band.union(local_band)
        plan = replace(plan, canvas_band=canvas_band)
        plans.append(
            stack_linear_wrap_rows(
                [plan, *(row_plan(record_index, row)[0] for row in range(1, transform.row_count))],
                _ordinary_row_gap(canvas_config),
            )
            if transform is not None
            else plan
        )
        definition_geometries.append(replace(beside, local_band=local_band))
    return plans, definition_geometries


def _ordinary_row_gap(canvas_config: LinearCanvasConfigurator) -> float:
    """The gap between stacked records, and between the rows of a wrapped record."""

    return float(canvas_config.cds_padding) * 1.5


def _labels_in_first_row_frame(
    labels: list[LinearLabel], plan: LinearRecordVerticalPlan, feature_slot_id: str | None
) -> list[LinearLabel]:
    """A wrapped record's feature labels moved from each row's frame into its first row's."""

    if not plan.wrap_rows or feature_slot_id is None:
        return labels
    shifts = plan.slot_row_shifts(feature_slot_id)
    return [linear_label_on_its_wrap_row(label, shifts[label["wrap_row"]]) for label in labels]


def _comparison_records_by_boundary(
    comparisons: list[LinearComparison],
    rows_by_record: tuple[int, ...],
) -> dict[int, tuple[int, ...]]:
    records_by_boundary: dict[int, set[int]] = {}
    for comparison in comparisons:
        query_index = int(comparison.query_record_index)
        subject_index = int(comparison.subject_record_index)
        query_row = int(rows_by_record[query_index])
        subject_row = int(rows_by_record[subject_index])
        boundary = min(query_row, subject_row)
        records_by_boundary.setdefault(boundary, set()).update(
            (query_index, subject_index)
        )
    return {
        boundary: tuple(sorted(record_indices))
        for boundary, record_indices in records_by_boundary.items()
    }


def _record_collision_bands(
    *,
    plan: LinearRecordVerticalPlan,
    definition_geometry: _LinearRecordDefinitionGeometry,
    sequence_width: float,
    record_x: float,
    multi_record_enabled: bool,
    keep_definition_left_aligned: bool,
    definition_column_width: float,
    row_definition_width: float,
    definition_gap: float,
    text_anchor: str,
    definition_column_left: float,
) -> tuple[CollisionBand, ...]:
    """Build alignment-local collision domains for one placed record."""

    width = max(0.0, float(sequence_width))
    x = float(record_x)
    bands: list[CollisionBand] = [
        CollisionBand(
            "body",
            x,
            x + width,
            plan.record_body_band.top_y,
            plan.record_body_band.bottom_y,
        ),
        CollisionBand(
            "comparison",
            x,
            x + width,
            plan.comparison_exclusion_band.top_y,
            plan.comparison_exclusion_band.bottom_y,
        ),
    ]

    bands.extend(
        _record_definition_bands(
            definition_geometry,
            sequence_width=width,
            record_x=x,
            multi_record_enabled=multi_record_enabled,
            keep_definition_left_aligned=keep_definition_left_aligned,
            definition_column_width=definition_column_width,
            row_definition_width=row_definition_width,
            definition_gap=definition_gap,
            text_anchor=text_anchor,
            definition_column_left=definition_column_left,
        )
    )
    return tuple(bands)


def _record_definition_bands(
    definition_geometry: _LinearRecordDefinitionGeometry,
    *,
    sequence_width: float,
    record_x: float,
    multi_record_enabled: bool,
    keep_definition_left_aligned: bool,
    definition_column_width: float,
    row_definition_width: float,
    definition_gap: float,
    text_anchor: str,
    definition_column_left: float,
) -> list[CollisionBand]:
    """The boxes of a placed record's definition text (``place_linear_definition`` on x, the
    definition geometry's bands on y), in the coordinates of ``record_x``."""

    width = max(0.0, float(sequence_width))
    x = float(record_x)
    bands: list[CollisionBand] = []
    local_band = definition_geometry.local_band
    local_width = max(0.0, float(definition_geometry.local_width))
    if local_band is not None and local_width > 0.0:
        placement = place_linear_definition(
            width=local_width, column_width=definition_column_width,
            record_x=x,
            gap=max(0.0, float(definition_gap)), keep_left=keep_definition_left_aligned,
            text_anchor=text_anchor, sequence_width=width if multi_record_enabled else None,
            column_left=definition_column_left,
        )
        bands.append(
            CollisionBand(
                "definition",
                placement.left,
                placement.right,
                local_band.top_y,
                local_band.bottom_y,
            )
        )

    row_band = definition_geometry.row_band
    actual_row_width = max(0.0, float(definition_geometry.row_width))
    if row_band is not None and actual_row_width > 0.0:
        placement = place_linear_definition(
            width=actual_row_width, column_width=row_definition_width,
            record_x=x,
            gap=max(0.0, float(definition_gap)), keep_left=keep_definition_left_aligned,
            column_left=definition_column_left,
        )
        bands.append(
            CollisionBand(
                "definition",
                placement.left,
                placement.right,
                row_band.top_y,
                row_band.bottom_y,
            )
        )
    return bands


def _definition_label_boxes(
    bands: Sequence[CollisionBand], *, record_x: float, features_origin_y: float
) -> tuple[Aabb, ...]:
    """A record's definition boxes (``_record_definition_bands`` at ``record_x``) in the coordinates
    of its feature labels: from the record's start and the features slot's origin."""

    return tuple(
        Aabb(band.x_start - record_x, band.top_y - features_origin_y, band.x_end - record_x, band.bottom_y - features_origin_y)
        for band in bands
        if band.kind == "definition" and band.x_end > band.x_start
    )


def _linear_axis_ruler_bounds(
    record: SeqRecord,
    placement: LinearRecordPlacement,
    plan: LinearRecordVerticalPlan,
    *,
    canvas_config: LinearCanvasConfigurator,
    render_context: LinearRecordRenderContext,
    cfg: GbdrawConfig,
    record_local_ruler: bool,
    record_transform: RecordDisplayTransform | None = None,
) -> Aabb | None:
    """Return exact horizontal text extents for a ruler drawn on the record axis."""

    if not render_context.axis_ruler_enabled or placement.sequence_width <= 0.0:
        return None
    interval = cfg.objects.scale.interval
    tick_interval = (
        int(interval)
        if interval is not None and int(interval) > 0
        else linear_tick_interval(max(1, int(canvas_config.longest_genome)))
    )
    if tick_interval <= 0:
        return None
    axis_length = drawn_axis_length(record, record_transform)
    # Each tick: the base its label names, its place along the axis (0..1), and its wrap row.
    ticks: list[tuple[int, float, int]]
    if record_transform is not None and record_transform.projects:
        ticks = [
            (tick.base, tick.offset / axis_length, tick.row)
            for tick in source_ruler_ticks(record_transform, tick_interval)
        ]
    else:
        ticks = [
            (tick.base, tick.position, 0)
            for tick in coordinate_map_ticks(record, len(record.seq), tick_interval, record_local=record_local_ruler)
        ]
    if not ticks:
        return None

    bounds: list[Aabb] = []
    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param(
        canvas_config.length_param
    )
    row_offsets = plan.row_offsets
    for coordinate, offset, row in ticks:
        x = placement.sequence_width * offset
        row_y = row_offsets[row] if record_transform is not None and record_transform.projects else 0.0
        label = format_linear_tick_label(
            coordinate,
            context_length=max(1, int(canvas_config.longest_genome)),
            tick_interval=tick_interval,
        )
        width, _height = drawn_text_box(label, cfg.objects.text.font_family, font_size, cfg.objects.scale.font_weight)
        bounds.append(
            Aabb(
                placement.x + x - (0.5 * width),
                placement.axis_y + row_y + plan.axis_band.top_y,
                placement.x + x + (0.5 * width),
                placement.axis_y + row_y + plan.axis_band.bottom_y,
            )
        )
    return union_aabbs(bounds)


class LinearScalarAxisLeftRoom(NamedTuple):
    """How far left of a record's axis start its left-facing scalar axes (depth, GC%) draw: their line
    (``axis_px``, 0 without an axis) and their tick labels (``tick_labels_px``, 0 without labels)."""

    axis_px: float
    tick_labels_px: float


def _linear_scalar_axis_left_room(
    slots: Sequence[LinearResolvedTrack],
    *,
    gc_config_of: LinearGcConfigOf,
    depth_by_index: Mapping[int, DepthTrackData],
    cfg: GbdrawConfig,
) -> LinearScalarAxisLeftRoom:
    """One record's left-facing scalar axes as the drawers draw them, measured once from the resolved
    data with the drawn (SVG) text width: the figure bound and the definition column read it (OV-1399,
    OV-1403). ``gc_config_of`` gives the GC content configuration a slot (by its params) of the record
    (by its depth tracks) is drawn with."""

    axis_px = 0.0
    tick_labels_px = 0.0
    for slot in slots:
        axis_labels: ScalarAxisLabels
        if slot.renderer == "depth":
            depth_track = depth_by_index.get(_depth_slot_track_index(slot))
            if depth_track is None or not bool(depth_track.config.show_axis):
                continue
            depth_config = depth_track.config
            font_size = linear_scalar_axis_tick_font_size_px(depth_config, float(slot.height))
            axis_min, axis_max = _depth_axis_bounds(depth_track.df, depth_config.min_depth, depth_config.max_depth)
            axis_labels = depth_axis_labels(
                axis_min,
                axis_max,
                label_style=depth_config.axis_labels,
                min_depth=depth_config.min_depth,
                normalize=bool(depth_config.normalize),
                show_ticks=bool(depth_config.show_ticks),
                large_tick_interval=depth_config.large_tick_interval,
                length_px=float(slot.height),
                label_line_px=font_size,
            )
        elif slot.renderer == "dinucleotide_content":
            gc_slot_config = gc_config_of(slot.params, depth_by_index)
            if str(gc_slot_config.mode) != "percent" or not bool(gc_slot_config.show_axis):
                continue
            font_size = linear_scalar_axis_tick_font_size_px(gc_slot_config, float(slot.height))
            axis_labels = percent_axis_labels(
                float(0.0 if gc_slot_config.min_percent is None else gc_slot_config.min_percent),
                float(100.0 if gc_slot_config.max_percent is None else gc_slot_config.max_percent),
                label_style=gc_slot_config.axis_labels,
                show_ticks=bool(gc_slot_config.show_ticks),
                large_tick_interval=gc_slot_config.large_tick_interval,
                length_px=float(slot.height),
                label_line_px=font_size,
            )
        else:
            continue
        axis_px = max(axis_px, 0.5 * LINEAR_SCALAR_AXIS_STROKE_PX)
        texts = [tick.label for tick in axis_labels.ticks] + [axis_labels.range_label]
        for label in (text for text in texts if text is not None):
            width, _height = calculate_svg_bbox_dimensions(
                label,
                cfg.objects.text.font_family,
                font_size,
            )
            tick_labels_px = max(
                tick_labels_px, LINEAR_SCALAR_AXIS_TICK_PX + LINEAR_SCALAR_AXIS_LABEL_GAP_PX + float(width)
            )
    return LinearScalarAxisLeftRoom(axis_px, tick_labels_px)


def _linear_annotation_bounds(
    record_index: int,
    record_length: int,
    placement: LinearRecordPlacement,
    plan: LinearRecordVerticalPlan,
    *,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    annotation_labels: AnnotationLabelsBySlot,
) -> list[Aabb]:
    """Return structured mark and label bounds for rendered annotation slots."""

    scale = placement.sequence_width / max(1, int(record_length))
    bounds: list[Aabb] = []
    for slot in plan.slots:
        if slot.renderer != "annotations":
            continue
        painted = [
            band.translate(row.axis_offset_y)
            for row in plan.rows
            if (band := row.plan.slot_by_id(slot.id).paint_band) is not None
        ]
        track = annotation_track_layouts.get(slot.id)
        if track is None or not painted:
            continue
        paint_band = union_vertical_bands(painted)
        params = annotation_track_params_from_mapping(slot.params)
        labels = {label.annotation_id: label for label in annotation_labels.get(slot.id, ())}
        for item in track.placements:
            annotation = item.annotation
            if annotation.record_index != record_index:
                continue
            style = effective_annotation_style(annotation, params)
            half_stroke = 0.5 * max(0.0, float(style.stroke_width))
            min_x = min(float(start) for start, _end in annotation.geometry_segments) * scale
            max_x = max(float(end) for _start, end in annotation.geometry_segments) * scale
            label = labels.get(annotation.id)
            if label is not None:
                min_x = min(min_x, label.footprint.min_x)
                max_x = max(max_x, label.footprint.max_x)
            bounds.append(
                Aabb(
                    placement.x + min_x - half_stroke,
                    placement.axis_y + paint_band.top_y,
                    placement.x + max_x + half_stroke,
                    placement.axis_y + paint_band.bottom_y,
                )
            )
    return bounds


def _place_linear_annotation_labels_by_record(
    *,
    records: Sequence[SeqRecord],
    axis_lengths: Sequence[int],
    placements: Mapping[int, LinearRecordPlacement],
    plans: Sequence[LinearRecordVerticalPlan],
    collision_bands: Sequence[tuple[CollisionBand, ...]],
    feature_labels: list[list[LinearLabel]],
    feature_dicts: Sequence[FeatureDict],
    feature_geometries: Sequence[LinearFeatureLaneGeometry],
    separate_strands: bool,
    feature_stroke_width: float,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack],
    font_family: str,
    label_style: FeatureLabelStyle,
    gap: float,
    new_rows: bool,
    label_frames: Sequence[LinearLabelFrame] | None,
    definition_boxes: Sequence[tuple[Aabb, ...]],
) -> tuple[list[PlacedAnnotationLabels], list[list[LinearLabel]]]:
    """Each record's annotation labels, placed along their spans in record order, and every
    record's feature labels after the moves made for them (LP-D82), which the label fit of their
    row (``label_frames``, as ``linear_label_fit`` reads them, with each record's ``definition_boxes``)
    counts no larger.

    Obstacles are the record's own feature labels and leaders, every record's
    name and definition text, and, from the other records of the same row,
    their bars (the ``body`` collision bands the figure bounds read), feature
    labels, leaders and the annotation labels placed before (LP-D42). Each
    record's region labels come with the pairs that still overlap (LP-D68).
    The record's feature glyphs bound the empty area inside a filled mark (LP-D75).
    """

    def labels_of(index: int) -> list[LinearLabel]:
        # Read when used: an earlier record's region labels may have moved its feature labels.
        return feature_labels[index] if index < len(feature_labels) else []

    def own_ink(index: int) -> list[LabelObstacle]:
        return linear_feature_label_obstacles(labels_of(index), plans[index], index, label_style=label_style)

    label_row_of = {index: row for row in linear_label_rows(label_frames, len(records)) for index, _x in row}
    placed: list[PlacedAnnotationLabels] = []
    for index, _record in enumerate(records):
        here = placements[index]
        obstacles = own_ink(index)
        label_row = label_row_of[index]
        for other, there in placements.items():
            same_row = other != index and there.row == here.row
            dx, dy = there.x - here.x, there.axis_y - here.axis_y
            obstacles.extend(
                LabelObstacle(
                    Aabb(band.x_start - here.x, dy + band.top_y, band.x_end - here.x, dy + band.bottom_y),
                    (other, "band", band_index),
                )
                for band_index, band in enumerate(collision_bands[other])
                if band.x_end > band.x_start and (band.kind == "definition" or (same_row and band.kind == "body"))
            )
            if same_row:
                inks = own_ink(other) + (
                    linear_annotation_label_inks(placed[other].by_slot, plans[other], other) if other < index else []
                )
                obstacles.extend(item._replace(box=item.box.translated(dx, dy)) for item in inks)
        placed.append(
            place_linear_annotation_labels(
                record_index=index,
                record_length=axis_lengths[index],
                sequence_width=here.sequence_width,
                plan=plans[index],
                annotation_track_layouts=annotation_track_layouts,
                obstacles=obstacles,
                feature_labels=labels_of(index),
                # Only a record with annotation labels to place reads its glyphs.
                glyphs=linear_feature_glyph_boxes(
                    feature_dicts[index],
                    feature_geometries[index],
                    plans[index],
                    separate_strands=separate_strands,
                    scale=here.sequence_width / max(1, axis_lengths[index]),
                    stroke_width=feature_stroke_width,
                )
                if any(slot.renderer == "annotations" and slot.id in annotation_track_layouts for slot in plans[index].slots)
                else (),
                font_family=font_family,
                label_style=label_style,
                gap=gap,
                new_rows=new_rows,
                row=[(labels_of(other), x) for other, x in label_row],
                row_position=next(position for position, (other, _x) in enumerate(label_row) if other == index),
                definition_boxes=definition_boxes[index],
            )
        )
        feature_labels = with_feature_label_rows(feature_labels, placed[-1].rows)
    return placed, feature_labels


def _collect_linear_primary_bounds(
    *,
    records: list[SeqRecord],
    record_placements: dict[int, LinearRecordPlacement],
    record_plans: list[LinearRecordVerticalPlan],
    record_collision_bands: list[tuple[CollisionBand, ...]],
    labels_by_record: list[list[LinearLabel]],
    record_scalar_axis_rooms: Sequence[LinearScalarAxisLeftRoom],
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    annotation_labels_by_record: Sequence[AnnotationLabelsBySlot],
    normalized_comparisons: list[LinearComparison],
    canvas_config: LinearCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    blast_config: object,
    render_context: LinearRecordRenderContext,
    cfg: GbdrawConfig,
    multi_record_enabled: bool,
    length_bar_bounds: Aabb | None,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
) -> Aabb:
    """Collect one authoritative plot-space box without inspecting rendered SVG."""

    axis_stroke = cfg.objects.axis.linear.stroke_width.for_length_param(
        canvas_config.length_param
    )
    half_stroke = 0.5 * max(
        0.0,
        float(axis_stroke),
        float(feature_config.block_stroke_width),
        float(feature_config.line_stroke_width),
        float(gc_config.stroke_width),
        float(skew_config.stroke_width),
    )
    horizontal_offset = float(canvas_config.horizontal_offset)
    bounds: list[Aabb] = []
    for record_index, record in enumerate(records):
        placement = record_placements[record_index]
        plan = record_plans[record_index]
        scalar_axis = record_scalar_axis_rooms[record_index]
        bounds.append(
            Aabb(
                horizontal_offset
                + placement.x
                - max(half_stroke, scalar_axis.axis_px, scalar_axis.tick_labels_px),
                placement.axis_y + plan.paint_band.top_y - half_stroke,
                horizontal_offset + placement.x + placement.sequence_width + half_stroke,
                placement.axis_y + plan.paint_band.bottom_y + half_stroke,
            )
        )

        for band in record_collision_bands[record_index]:
            if band.kind != "definition" or band.x_end <= band.x_start:
                continue
            bounds.append(
                Aabb(
                    horizontal_offset + band.x_start,
                    placement.axis_y + band.top_y,
                    horizontal_offset + band.x_end,
                    placement.axis_y + band.bottom_y,
                )
            )

        feature_origin_y = plan.features_origin_y
        for label in labels_by_record[record_index]:
            left, right, top, bottom = calculate_label_bounds(label)
            bounds.append(
                Aabb(
                    horizontal_offset + placement.x + left,
                    placement.axis_y + feature_origin_y + top,
                    horizontal_offset + placement.x + right,
                    placement.axis_y + feature_origin_y + bottom,
                )
            )
            # An ``above_feature`` label carries both ends of its leader.
            leader_points = (
                [
                    (float(label["leader_start_x"]), float(label["leader_start_y"])),
                    (float(label["leader_end_x"]), float(label["leader_end_y"])),
                ]
                if "leader_start_x" in label
                else []
            )
            if leader_points:
                label_half_stroke = 0.5 * label_leader_width_px(label)
                bounds.append(
                    Aabb(
                        horizontal_offset
                        + placement.x
                        + min(point[0] for point in leader_points)
                        - label_half_stroke,
                        placement.axis_y
                        + feature_origin_y
                        + min(point[1] for point in leader_points)
                        - label_half_stroke,
                        horizontal_offset
                        + placement.x
                        + max(point[0] for point in leader_points)
                        + label_half_stroke,
                        placement.axis_y
                        + feature_origin_y
                        + max(point[1] for point in leader_points)
                        + label_half_stroke,
                    )
                )

        ruler_bounds = _linear_axis_ruler_bounds(
            record,
            placement,
            plan,
            canvas_config=canvas_config,
            render_context=render_context,
            cfg=cfg,
            record_local_ruler=multi_record_enabled,
            record_transform=(record_transforms[record_index] if record_transforms is not None else None),
        )
        if ruler_bounds is not None:
            bounds.append(ruler_bounds.translated(horizontal_offset, 0.0))
        bounds.extend(
            item.translated(horizontal_offset, 0.0)
            for item in _linear_annotation_bounds(
                record_index,
                drawn_axis_length(record, record_transforms[record_index] if record_transforms is not None else None),
                placement,
                plan,
                annotation_track_layouts=annotation_track_layouts,
                annotation_labels=annotation_labels_by_record[record_index],
            )
        )

    comparison_half_stroke = 0.5 * max(
        0.0,
        float(getattr(blast_config, "stroke_width", 0.0)),
    )
    for comparison in normalized_comparisons:
        query = record_placements[comparison.query_record_index]
        subject = record_placements[comparison.subject_record_index]
        query_y = query.facing_comparison_anchor_y(subject)
        subject_y = subject.facing_comparison_anchor_y(query)
        bounds.append(
            Aabb(
                horizontal_offset + min(query.x, subject.x) - comparison_half_stroke,
                min(query_y, subject_y) - comparison_half_stroke,
                horizontal_offset
                + max(
                    query.x + query.sequence_width,
                    subject.x + subject.sequence_width,
                )
                + comparison_half_stroke,
                max(query_y, subject_y) + comparison_half_stroke,
            )
        )

    if length_bar_bounds is not None:
        bounds.append(length_bar_bounds)
    primary_bounds = union_aabbs(bounds)
    if primary_bounds is None or primary_bounds.width <= 0.0 or primary_bounds.height <= 0.0:
        raise RuntimeError("linear assembly produced no primary painted bounds")
    return primary_bounds


def _precalculate_gc_dataframes(
    records: list[SeqRecord],
    record_transforms: Sequence[RecordDisplayTransform] | None,
    *,
    window: int,
    step: int,
    dinucleotide: str,
    enabled: bool,
) -> list[DataFrame | None]:
    """Build GC/skew data once per record and share it across linear groups.

    A record's windows wrap round its ends when it is drawn as circular (its transform)."""
    if not enabled:
        return [None for _ in records]
    return [
        render_stage("gc_skew_data", skew_df, record, window, step, dinucleotide, circular=transform.is_circular)
        for record, transform in zip(records, bind_record_transforms(records, record_transforms), strict=True)
    ]






def _final_record_translations(
    *,
    record_keys: Sequence[str],
    placements: Mapping[int, LinearRecordPlacement],
    layout: LinearMultiRecordOptions | None,
    similarity_alignment: SimilarityAlignmentPlan | None,
    anchor_centers: Sequence[float | None] | None,
    feature_anchor: FeatureAnchorPlacement | None = None,
) -> tuple[tuple[float, float], ...]:
    """Resolve base translations and absolute anchor alignment exactly once.

    The records align at display positions: the Similarity alignment anchor
    centres (``anchor_centers``) with ``similarity_alignment``, else the
    ``anchor_feature`` starts and reference the planner decided
    (``feature_anchor``). An aligned or rotated record moves so that its
    position meets the reference's; a rotated reference is drawn at x 0 (its y
    kept), any other reference and every other record keep their base
    translation.

    ``record_keys`` name the record each drawn unit draws: the windows of one
    record share its key, its base translation and its move, which the window
    holding its anchor centre decides (D-50).
    """

    base_by_key = {
        str(item.record_key): (float(item.x), float(item.y))
        for item in (layout.record_translations if layout is not None else ())
    }
    base = tuple(base_by_key.get(str(key), (0.0, 0.0)) for key in record_keys)
    if similarity_alignment is None:
        if feature_anchor is None:
            return base
        centers = feature_anchor.starts
        reference_index = feature_anchor.reference_index
        if feature_anchor.statuses[reference_index] == "rotated":
            base = tuple(
                (0.0, base_y) if index == reference_index else (base_x, base_y)
                for index, (base_x, base_y) in enumerate(base)
            )
        # Each moved unit and the unit whose position it meets the reference's at.
        aligned = {
            index: index for index, status in enumerate(feature_anchor.statuses)
            if status in ("rotated", "aligned") and index != reference_index
        }
    else:
        centers = tuple(anchor_centers or ())
        similarity_alignment.validate_record_coverage(tuple(dict.fromkeys(str(key) for key in record_keys)))
        if len(centers) != len(record_keys):
            raise ValidationError(
                "Similarity alignment requires one projected anchor center per drawn unit."
            )
        anchor_units = {str(key): index for index, key in enumerate(record_keys) if centers[index] is not None}
        if similarity_alignment.reference.record_key not in anchor_units:
            raise ValidationError("Similarity alignment reference center is missing.")
        reference_index = anchor_units[similarity_alignment.reference.record_key]
        decisions = {decision.record_key: decision for decision in similarity_alignment.records}
        for decision in decisions.values():
            if decision.status is AlignmentDecisionStatus.ALIGNED and decision.record_key not in anchor_units:
                raise ValidationError(
                    f"Similarity alignment center is missing for record {decision.record_key!r}."
                )
        aligned = {
            index: anchor_units[str(record_key)] for index, record_key in enumerate(record_keys)
            if decisions[str(record_key)].status is AlignmentDecisionStatus.ALIGNED
        }
    reference_world_x = (
        placements[reference_index].x_for_position(cast(float, centers[reference_index]))
        + base[reference_index][0]
    )
    return tuple(
        (
            reference_world_x - placements[aligned[index]].x_for_position(cast(float, centers[aligned[index]]))
            if index in aligned
            else base_x,
            base_y,
        )
        for index, (base_x, base_y) in enumerate(base)
    )


def build_linear_legend_table(
    records: Sequence[SeqRecord],
    *,
    sources: LegendSources,
    comparison: ComparisonLegendSource | None,
    pairwise_identity_legend: PairwiseIdentityLegend,
    depth: DepthLegendSource | None,
    show_gc: bool,
    show_skew: bool,
    slots: tuple[LegendSlotSource, ...] | None,
    depth_tracks: tuple[DepthLegendSource, ...] | None,
    annotations: ResolvedAnnotationBundle,
    edits: LegendRowEdits = LegendRowEdits(),
) -> LegendTable:
    """Build the Linear legend table (render stage ``legend_table``).

    ``depth`` is given only when one depth track has a singleton row. The
    request's Legend row edits are applied last.
    """

    features_present, used_color_rules, default_used_features = feature_legend_presence(
        records, sources.features, slots
    )
    legend_table = prepare_legend_table(
        sources.gc,
        sources.skew,
        sources.features,
        features_present,
        comparison=comparison,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
        depth=depth,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=depth is not None,
        pairwise_identity_legend=pairwise_identity_legend,
    )
    legend_table = _sync_legend_table_for_linear_slots(
        legend_table,
        slots=slots,
        gc=sources.gc,
        skew=sources.skew,
        depth_tracks=depth_tracks,
    )
    legend_table = sync_feature_density_legend_entries(legend_table, sources.feature_density, slots=slots)
    return apply_legend_row_edits(sync_annotation_legend_entries(legend_table, annotations, slots), edits)


@dataclass(frozen=True, slots=True)
class _LinearRowFrame:
    """Where each record sits horizontally, solved once before its tracks and labels.

    ``rows_by_record`` is each record's row and ``ordered_record_indices`` the
    records in row and column order (the figure draws them in ``linear_paint_order``);
    ``multi_record_enabled`` is true when a row holds more than
    one record, and then each record's definition is split into its own lines
    and its row's, and its feature labels are stacked and counted in its row's
    frame (``label_frames``: its row and its start in it, OV-1188; ``None``
    when every record is alone in its row). The placements have zero vertical
    extents, which the row stacking fills in.
    """

    presentations: tuple[RecordPresentation | None, ...]
    record_keys: tuple[RecordKey, ...]
    # The drawn unit of each record: the record, or the window of it, it shows (``record_units``).
    units: tuple[RecordUnit, ...]
    ordered_record_indices: tuple[int, ...]
    rows_by_record: tuple[int, ...]
    multi_record_enabled: bool
    row_leading_indices: set[int]
    label_frames: list[LinearLabelFrame] | None
    axis_lengths: list[int]
    horizontal_layout: LinearHorizontalLayout
    horizontal_placements: dict[int, LinearRecordPlacement]
    sequence_widths: list[float]
    # The rows that hold record windows (PD-OI-122), their captions and marks, and the look of both.
    row_windows: tuple[LinearRowWindows, ...]
    row_window_style: RowWindowStyle

    def row_windows_of(self, record_index: int) -> LinearRowWindows | None:
        """The windowed row ``record_index`` is drawn in, or ``None`` when its row keeps each record's text."""

        return next((row for row in self.row_windows if record_index in row.record_indices), None)


def _plan_linear_row_windows(
    units: tuple[RecordUnit, ...],
    presentations: tuple[RecordPresentation | None, ...],
    ordered_record_indices: tuple[int, ...],
    rows_by_record: tuple[int, ...],
) -> tuple[LinearRowWindows, ...]:
    """Each windowed row's captions and marks, its units in drawing order (``plan_row_windows``)."""

    units_by_row: dict[int, list[int]] = {}
    for record_index in ordered_record_indices:
        units_by_row.setdefault(rows_by_record[record_index], []).append(record_index)
    plans = (
        plan_row_windows(
            [units[index] for index in indices],
            [
                (presentation.label or "") if (presentation := presentations[index]) is not None else ""
                for index in indices
            ],
        )
        for indices in units_by_row.values()
    )
    return tuple(plan for plan in plans if plan is not None)


def _linear_row_frame(canvas_config: LinearCanvasConfigurator, record_inputs: LinearRecordInputs) -> _LinearRowFrame:
    """Each record's row and horizontal place: its x, width and scale, in both row modes."""

    records = record_inputs.records
    linear_layout = record_inputs.linear_layout
    presentations: tuple[RecordPresentation | None, ...] = (
        tuple(record_inputs.record_presentations) if record_inputs.record_presentations is not None else (None,) * len(records)
    )
    record_keys = stable_record_keys(records)
    units = record_units(records, record_inputs.record_transforms)
    ordered_record_indices, rows_by_record = resolve_record_row_positions(
        records,
        linear_layout.multi_record_positions if linear_layout is not None else None,
    )
    row_counts: dict[int, int] = {}
    for row in rows_by_record:
        row_counts[row] = row_counts.get(row, 0) + 1
    multi_record_enabled = any(count > 1 for count in row_counts.values())
    row_leading_indices: set[int] = set()
    seen_rows: set[int] = set()
    for record_index in ordered_record_indices:
        row = rows_by_record[record_index]
        if row in seen_rows:
            continue
        seen_rows.add(row)
        row_leading_indices.add(record_index)
    axis_lengths = _drawn_axis_lengths(records, record_inputs.record_transforms)
    row_windows = _plan_linear_row_windows(units, presentations, ordered_record_indices, rows_by_record)
    row_window_style = RowWindowStyle.of(
        canvas_config.profile.config, str(canvas_config.length_param), float(canvas_config.default_cds_height)
    )
    insets: dict[int, tuple[float, float]] = {}
    for windows in row_windows:
        insets.update(row_window_insets(windows, row_window_style))
    horizontal_layout = solve_linear_horizontal_layout(
        [
            LinearRecordMeasurement(
                record_index=index,
                record_key=RecordKey(str(record_keys[index])),
                sequence_length=axis_lengths[index],
                left_inset=insets.get(index, (0.0, 0.0))[0],
                right_inset=insets.get(index, (0.0, 0.0))[1],
            )
            for index, _record in enumerate(records)
        ],
        rows_by_record,
        available_width=float(canvas_config.alignment_width),
        record_gap_px=(linear_layout.record_gap_px if linear_layout is not None else 24.0),
        align_center=bool(canvas_config.align_center),
        normalize_length=bool(canvas_config.normalize_length),
        record_order=ordered_record_indices,
    )
    horizontal_placements: dict[int, LinearRecordPlacement] = {
        index: horizontal_layout.placement_for_index(index) for index in range(len(records))
    }
    return _LinearRowFrame(
        presentations=presentations,
        record_keys=record_keys,
        units=units,
        ordered_record_indices=ordered_record_indices,
        rows_by_record=rows_by_record,
        multi_record_enabled=multi_record_enabled,
        row_leading_indices=row_leading_indices,
        label_frames=(
            [LinearLabelFrame(rows_by_record[index], horizontal_placements[index].x) for index in range(len(records))]
            if multi_record_enabled
            else None
        ),
        axis_lengths=axis_lengths,
        horizontal_layout=horizontal_layout,
        horizontal_placements=horizontal_placements,
        sequence_widths=[horizontal_placements[index].sequence_width for index in range(len(records))],
        row_windows=row_windows,
        row_window_style=row_window_style,
    )


@dataclass(frozen=True, slots=True)
class _LinearSlotPlan:
    """The requested slots after annotation planning, or the default slots, before the feature underlays."""

    linear_track_slots: list[LinearTrackSlot]
    resolved_annotations: ResolvedAnnotationBundle
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack]
    auto_annotation_slot_ids: frozenset[str]
    normalized_linear_track_slots: list[NormalizedLinearTrackSlot]


def _linear_slot_plan(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    track_inputs: LinearTrackInputs,
    frame: _LinearRowFrame,
    label_font_size: float,
) -> _LinearSlotPlan:
    """Plan the annotation slots on the requested slots (or the default ones) and normalize them once."""

    profile = canvas_config.profile
    records = record_inputs.records
    record_depth_tracks = track_inputs.record_depth_tracks
    linear_track_slots = track_inputs.linear_track_slots
    (
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
    ) = _prepare_linear_annotation_tracks(
        records,
        track_inputs.annotations,
        linear_track_slots,
        profile=profile,
        record_depth_tracks=record_depth_tracks,
        label_font_size=label_font_size,
        sequence_widths=frame.sequence_widths,
        record_transforms=record_inputs.record_transforms,
        axis_index=track_inputs.linear_track_axis_index,
    )
    if linear_track_slots is None:
        linear_track_slots = default_linear_track_slots(
            show_features=True,
            show_depth=bool(profile.show_depth and record_depth_tracks),
            depth_track_count=max(1, depth_track_count(record_depth_tracks)),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            dinucleotide=str(configurators.gc.dinucleotide),
            track_layout=profile.track_layout,
        )
    density_ids = tuple(track.id for track in track_inputs.feature_density_tracks)
    if track_inputs.linear_track_slots is None:
        linear_track_slots = with_feature_density_slots(
            linear_track_slots,
            density_ids,
            lambda track_id: LinearTrackSlot(
                id=feature_density_slot_id(track_id), renderer="feature_density", side="below", params={"track": track_id}
            ),
        )
    check_feature_density_slots(linear_track_slots, density_ids)
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        track_inputs.linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        record_depth_tracks,
    )
    return _LinearSlotPlan(
        linear_track_slots=linear_track_slots,
        resolved_annotations=resolved_annotations,
        annotation_track_layouts=annotation_track_layouts,
        auto_annotation_slot_ids=auto_annotation_slot_ids,
        normalized_linear_track_slots=normalized_linear_track_slots,
    )


@dataclass(frozen=True, slots=True)
class _LinearComparisons:
    """The comparisons the diagram draws and the orthogroup labels they allow.

    ``matches`` are the drawn comparisons' match tables, in order;
    ``orthogroup_label_eligibility`` is set only for the ``orthogroup_top``
    label scope.
    """

    rows: LinearComparisonRows
    normalized: list[LinearComparison]
    matches: list[DataFrame]
    has_blast: bool
    label_scope: str
    orthogroup_label_eligibility: OrthogroupLabelEligibility | None


def _linear_comparisons(
    canvas_config: LinearCanvasConfigurator,
    blast_config: BlastMatchConfigurator,
    record_inputs: LinearRecordInputs,
    comparison_inputs: LinearComparisonInputs,
    frame: _LinearRowFrame,
) -> _LinearComparisons:
    """Read and filter the comparisons (stage ``comparison_rows``) and find which labels orthogroups allow."""

    records = record_inputs.records
    comparison_thresholds = ComparisonThresholds(
        evalue=blast_config.evalue,
        bitscore=blast_config.bitscore,
        identity=blast_config.identity,
        alignment_length=blast_config.alignment_length,
    )
    comparison_rows = render_stage(
        "comparison_rows",
        linear_comparison_rows,
        _read_comparison_files(comparison_inputs.blast_files, comparison_thresholds),
        comparison_inputs.protein_comparisons,
        tuple(comparison_inputs.linear_comparisons or ()),
        comparison_thresholds,
        tuple(records),
        tuple(frame.rows_by_record),
        bool(frame.multi_record_enabled),
        tuple(Path(path).name for path in comparison_inputs.blast_files or ()),
    )
    normalized_comparisons = list(comparison_rows.comparisons)
    comparisons = [item.matches for item in normalized_comparisons]

    label_scope = canvas_config.profile.label_scope
    orthogroup_label_eligibility = None
    if label_scope == "orthogroup_top":
        orthogroup_label_eligibility = build_orthogroup_label_eligibility(
            orthogroups=comparison_inputs.orthogroups,
            comparisons=comparisons,
            records=records,
        )
    if (
        label_scope == "orthogroup_top"
        and (
            orthogroup_label_eligibility is None
            or not orthogroup_label_eligibility.member_ids_by_record
        )
    ):
        logger.warning(
            "--show_labels orthogroup_top requires orthogroup metadata; no orthogroup-specific label suppression was applied."
        )
    return _LinearComparisons(
        rows=comparison_rows,
        normalized=normalized_comparisons,
        matches=comparisons,
        has_blast=bool(normalized_comparisons),
        label_scope=label_scope,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
    )


@dataclass(frozen=True, slots=True)
class _LinearTrackPlan:
    """Each record's feature layers and lanes, and the track layout after the feature underlays.

    The annotations, annotation track layouts and normalized slots here
    supersede the slot plan's: the underlays add slots and annotation sets.
    ``label_font_size`` is the diagram's one feature label font size, which
    region labels without their own font take too.
    """

    label_font_size: float
    feature_placement_slot: FeaturePlacementSlot | None
    record_feature_layers: list[FeatureBuildResult]
    record_feature_dicts: list[FeatureDict]
    feature_dom_index: LinearFeatureDomIndex
    resolved_annotations: ResolvedAnnotationBundle
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack]
    normalized_linear_track_slots: list[NormalizedLinearTrackSlot]
    linear_track_layout: LinearTrackLayout
    render_context: LinearRecordRenderContext
    resolved_depth_heights_by_index: dict[int, float]
    primary_depth_track_index: int
    record_feature_lane_geometries: list[LinearFeatureLaneGeometry]
    # Each record's feature density layers by track id.
    record_feature_density: Sequence[Mapping[str, FeatureDensityLayer]]


def _linear_track_plan(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    track_inputs: LinearTrackInputs,
    frame: _LinearRowFrame,
    slot_plan: _LinearSlotPlan,
    label_font_size: float,
) -> _LinearTrackPlan:
    """Build the feature layers (stage ``feature_dom_index``), add their underlay slots, resolve the track
    layout and measure each record's feature lanes. Refuses a wrapped record whose feature track is the axis."""

    profile = canvas_config.profile
    cfg = profile.config
    records = record_inputs.records
    record_transforms = record_inputs.record_transforms
    feature_config = configurators.feature
    resolved_feature_slot = next(
        (slot for slot in slot_plan.normalized_linear_track_slots if slot.renderer == "features"),
        None,
    )
    feature_placement_slot = (
        _linear_feature_placement_slot(resolved_feature_slot.side, profile.strandedness)
        if resolved_feature_slot is not None
        else None
    )
    record_feature_layers = _precalculate_feature_layers(
        records,
        feature_config,
        profile,
        record_transforms=record_transforms,
        placement_slot=feature_placement_slot,
    )
    record_feature_dicts = [
        result.foreground_features for result in record_feature_layers
    ]
    feature_dom_index = render_stage(
        "feature_dom_index", build_linear_feature_dom_index, tuple(record_feature_dicts)
    )
    (
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        _auto_annotation_slot_ids,
        private_underlay_prefix_count,
    ) = _add_linear_feature_underlays(
        records,
        record_feature_layers,
        slot_plan.linear_track_slots,
        slot_plan.resolved_annotations,
        slot_plan.annotation_track_layouts,
        slot_plan.auto_annotation_slot_ids,
        label_font_size=label_font_size,
        font_family=str(cfg.objects.text.font_family),
        sequence_widths=frame.sequence_widths,
        record_transforms=record_transforms,
    )
    effective_linear_track_axis_index = (
        track_inputs.linear_track_axis_index + private_underlay_prefix_count
        if track_inputs.linear_track_axis_index is not None
        else None
    )
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        effective_linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        track_inputs.record_depth_tracks,
    )
    linear_track_layout = resolve_linear_track_layout(
        normalized_linear_track_slots,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    depth_available = bool(configurators.depth is not None and track_inputs.record_depth_tracks)
    render_context = LinearRecordRenderContext(
        profile=profile,
        track_layout=linear_track_layout,
        depth_available=depth_available,
    )
    if any(transform.row_count > 1 for transform in record_transforms or ()):
        refuse_feature_track_on_the_axis(render_context.feature_track_layout)
    resolved_depth_heights_by_index = _resolved_linear_depth_heights(
        linear_track_layout
    )
    record_feature_lane_geometries = [
        measure_linear_feature_lanes(
            feature_dict,
            cds_height=float(canvas_config.cds_height),
            separate_strands=profile.strandedness,
            track_layout=render_context.feature_track_layout,
            axis_gap=canvas_config.track_axis_gap,
            stroke_width=max(
                float(feature_config.block_stroke_width),
                float(feature_config.line_stroke_width),
            ),
            include_nominal_lanes=bool(
                record_feature_layers[index].underlay_features
            ),
        )
        for index, feature_dict in enumerate(record_feature_dicts)
    ]
    return _LinearTrackPlan(
        label_font_size=label_font_size,
        feature_placement_slot=feature_placement_slot,
        record_feature_layers=record_feature_layers,
        record_feature_dicts=record_feature_dicts,
        feature_dom_index=feature_dom_index,
        resolved_annotations=resolved_annotations,
        annotation_track_layouts=annotation_track_layouts,
        normalized_linear_track_slots=normalized_linear_track_slots,
        linear_track_layout=linear_track_layout,
        render_context=render_context,
        resolved_depth_heights_by_index=resolved_depth_heights_by_index,
        primary_depth_track_index=min(resolved_depth_heights_by_index, default=0),
        record_feature_lane_geometries=record_feature_lane_geometries,
        # The records share one window and one scale per track.
        record_feature_density=feature_density_record_layers(
            records,
            track_inputs.feature_density_tracks,
            windows=track_inputs.feature_density_windows,
            transforms=bind_record_transforms(records, record_transforms),
        ),
    )


@dataclass(frozen=True, slots=True)
class _LinearDefinitionMetrics:
    """Each record's definition size: its own lines and, with records sharing a row, its row's lines.

    ``local_definition_line_kinds`` is ``None`` when definitions are not split
    by row; then the row-wide sizes are zero.
    """

    local_definition_line_kinds: list[frozenset[str]] | None
    row_definition_line_kinds: list[frozenset[str]]
    max_def_width: float
    definition_widths: list[float]
    definition_heights: list[float]
    row_definition_width: float
    row_definition_widths: list[float]
    row_definition_heights: list[float]
    # Where a row heading's box starts above the line it is placed on: ``None`` for a heading centred there, a
    # windowed row's heading hangs its list below its name line (``RowWindowHeadingGroup``).
    row_definition_tops: list[float | None]


def _row_window_heading_content(
    canvas_config: LinearCanvasConfigurator,
    record: SeqRecord,
    record_count: int,
    frame: _LinearRowFrame,
    windows: LinearRowWindows,
) -> GroupContent:
    """A windowed row's heading. It carries the identity (record index and ID) of the record whose label it shows,
    so a click on the name finds that record; without a label in the row, the leading record's. ``record`` is the
    record at ``windows.heading_record_index``."""

    index = windows.heading_record_index
    group_id = definition_group_svg_id(record.id, mode="linear", record_index=index, record_count=record_count)
    return build_group_content(
        RowWindowHeadingGroup,
        record,
        str(canvas_config.length_param),
        style=LinearDefinitionStyle.of(canvas_config.profile.config),
        windows=windows,
        window_style=frame.row_window_style,
        text_anchor="start" if bool(canvas_config.keep_definition_left_aligned) else "middle",
        group_id=f"{group_id}_row",
        record_index=index,
        record_presentation=frame.presentations[index] if windows.name_record_index is not None else None,
    )


def _linear_definition_metrics(
    canvas_config: LinearCanvasConfigurator,
    record_inputs: LinearRecordInputs,
    frame: _LinearRowFrame,
) -> _LinearDefinitionMetrics:
    """Measure each record's definition, split into its own lines and its row's when rows hold several records."""

    cfg = canvas_config.profile.config
    records = record_inputs.records
    local_definition_line_kinds: list[frozenset[str]] | None = None
    row_definition_line_kinds: list[frozenset[str]] = [
        frozenset() for _record in records
    ]
    if frame.multi_record_enabled:
        (
            local_definition_line_kinds,
            row_definition_line_kinds,
        ) = _split_definition_line_kinds(
            frame.presentations,
            rows_by_record=frame.rows_by_record,
            row_leading_indices=frame.row_leading_indices,
        )
        # A windowed row draws no text above its windows and its heading from the row's windows (PD-OI-122).
        for windows in frame.row_windows:
            for index in windows.record_indices:
                local_definition_line_kinds[index] = frozenset()
                row_definition_line_kinds[index] = frozenset()
    max_def_width, definition_widths, definition_heights = _precalculate_definition_metrics(
        records,
        canvas_config,
        cfg=cfg,
        line_kinds_by_record=local_definition_line_kinds,
        record_presentations=frame.presentations,
    )
    row_definition_width = 0.0
    row_definition_widths = [0.0 for _record in records]
    row_definition_heights = [0.0 for _record in records]
    if frame.multi_record_enabled:
        (
            row_definition_width,
            row_definition_widths,
            row_definition_heights,
        ) = _precalculate_definition_metrics(
            records,
            canvas_config,
            cfg=cfg,
            line_kinds_by_record=list(row_definition_line_kinds),
            record_presentations=frame.presentations,
        )
    row_definition_tops: list[float | None] = [None for _record in records]
    for windows in frame.row_windows:
        content = _row_window_heading_content(
            canvas_config, records[windows.heading_record_index], len(records), frame, windows
        )
        index = windows.leading_record_index
        measured = content.definition
        bounds = content.local_bounds
        if measured is None or bounds is None:  # pragma: no cover - the heading measures itself
            raise RuntimeError("a windowed row heading has no measurement")
        row_definition_widths[index] = float(measured.width)
        row_definition_heights[index] = float(measured.height)
        row_definition_tops[index] = float(bounds.min_y)
        row_definition_width = max(row_definition_width, float(math.ceil(measured.width)))
    return _LinearDefinitionMetrics(
        local_definition_line_kinds=local_definition_line_kinds,
        row_definition_line_kinds=row_definition_line_kinds,
        max_def_width=max_def_width,
        definition_widths=definition_widths,
        definition_heights=definition_heights,
        row_definition_width=row_definition_width,
        row_definition_widths=row_definition_widths,
        row_definition_heights=row_definition_heights,
        row_definition_tops=row_definition_tops,
    )


@dataclass(frozen=True, slots=True)
class _LinearDepthData:
    """Each record's depth tracks by track index, and what its left-facing scalar axes draw left of its axis.

    ``drawn_gc_config`` gives the GC content configuration a record's slot
    (by its params) is drawn with.
    """

    record_depth_by_index: list[dict[int, DepthTrackData]]
    representative_depth_data: list[DepthTrackData]
    representative_depth_configs: dict[int, DepthConfigurator]
    drawn_gc_config: LinearGcConfigOf
    record_scalar_axis_rooms: list[LinearScalarAxisLeftRoom]


def _linear_depth_data(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    track_inputs: LinearTrackInputs,
    tracks: _LinearTrackPlan,
) -> _LinearDepthData:
    """Build each record's depth frames and measure its scalar axes' left room.

    Sets the canvas's definition column gap from the widest tick labels: the
    definition column keeps its gap from them (LP-D122, OV-1399).
    """

    cfg = canvas_config.profile.config
    records = record_inputs.records
    gc_config = configurators.gc
    depth_config = configurators.depth
    render_context = tracks.render_context
    primary_depth_track_index = tracks.primary_depth_track_index
    resolved_depth_heights_by_index = tracks.resolved_depth_heights_by_index
    record_depth_data: list[list[DepthTrackData]] = (
        build_depth_track_dataframes(
            records,
            track_inputs.record_depth_tracks,
            record_inputs.record_transforms,
            base_config=depth_config,
            depth_df_builder=build_depth_df,
        )
        if render_context.depth_enabled
        else [[] for _ in records]
    )
    record_depth_by_index = [
        index_depth_track_row(row, record_index=record_index)
        for record_index, row in enumerate(record_depth_data)
    ]
    representative_depth_data = representative_depth_tracks(record_depth_data)
    representative_depth_configs = {
        track.track_index: track.config
        for track in representative_depth_data
    }

    def drawn_gc_config(
        slot_params: Mapping[str, object], depth_by_index: Mapping[int, DepthTrackData]
    ) -> GcContentConfigurator:
        """The GC content configuration a record's slot (by its params) is drawn with: its dinucleotide,
        and the depth axis font size when it shows a percent axis beside a depth axis."""

        primary_depth_track = depth_by_index.get(primary_depth_track_index)
        return _gc_config_matching_linear_depth_axis_font_size(
            gc_config=_clone_gc_config_with_dinucleotide(
                gc_config, slot_dinucleotide(slot_params, str(gc_config.dinucleotide))
            ),
            depth_config=(
                primary_depth_track.config
                if primary_depth_track is not None
                else representative_depth_configs.get(primary_depth_track_index, depth_config)
            ),
            depth_height=resolved_depth_heights_by_index.get(
                primary_depth_track_index,
                float(canvas_config.default_depth_height),
            ),
            depth_enabled=render_context.depth_enabled,
        )

    # What each record's left-facing scalar axes draw left of its axis start, measured once: the figure
    # bound reads it, and the definition column keeps its gap from the tick labels (LP-D122, OV-1399).
    record_scalar_axis_rooms = [
        _linear_scalar_axis_left_room(
            tracks.linear_track_layout.slots,
            gc_config_of=drawn_gc_config,
            depth_by_index=depth_by_index,
            cfg=cfg,
        )
        for depth_by_index in record_depth_by_index
    ]
    canvas_config.configure_definition_column_gap(
        max((room.tick_labels_px for room in record_scalar_axis_rooms), default=0.0)
    )
    return _LinearDepthData(
        record_depth_by_index=record_depth_by_index,
        representative_depth_data=representative_depth_data,
        representative_depth_configs=representative_depth_configs,
        drawn_gc_config=drawn_gc_config,
        record_scalar_axis_rooms=record_scalar_axis_rooms,
    )


@dataclass(frozen=True, slots=True)
class _LinearHorizontalFrame:
    """Each record's final translation, the definition column's left, and the definition boxes its labels keep clear of."""

    final_translations: tuple[tuple[float, float], ...]
    definition_column_left: float
    label_keep_clear: list[tuple[Aabb, ...]]


def _linear_horizontal_frame(
    canvas_config: LinearCanvasConfigurator,
    record_inputs: LinearRecordInputs,
    alignment: LinearAlignmentInputs,
    frame: _LinearRowFrame,
    tracks: _LinearTrackPlan,
    definitions: _LinearDefinitionMetrics,
) -> _LinearHorizontalFrame:
    """Each record's final translation, the definition column's left (the leftmost translated record), and the
    definition boxes its feature labels keep clear of (OV-1395) in their coordinates.

    They read only the records' horizontal places, which the labels do not
    move, so the drawing and the fit count read the same values. Reads the
    definition column gap the depth data set; refuses a feature anchor that
    would move a wrapped record.
    """

    cfg = canvas_config.profile.config
    records = record_inputs.records
    record_transforms = record_inputs.record_transforms
    feature_anchor = alignment.feature_anchor
    placements = frame.horizontal_placements
    translations = _final_record_translations(
        # The record each unit draws: the windows of a record move as one (D-50).
        record_keys=[unit.record_key for unit in frame.units],
        placements=placements,
        layout=record_inputs.linear_layout,
        similarity_alignment=alignment.similarity_alignment,
        anchor_centers=alignment.alignment_anchor_centers,
        feature_anchor=feature_anchor,
    )
    if feature_anchor is not None and record_transforms is not None:
        refuse_wrapped_anchor_moves(
            feature_anchor,
            [x for x, _y in translations],
            [transform.row_count > 1 for transform in record_transforms],
            [str(record.id) for record in records],
        )
    record_xs = [placements[index].x + translations[index][0] for index in range(len(records))]
    # The definition column keeps its gap from a window's marks too (its left inset).
    column_left = min((x - placements[index].left_inset for index, x in enumerate(record_xs)), default=0.0)
    keep_clear = [
        _definition_label_boxes(
            _record_definition_bands(
                _definition_geometry_beside_features(
                    0.5 * (geometry.occupied_band.top_y + geometry.occupied_band.bottom_y),
                    local_height=definitions.definition_heights[index],
                    local_width=definitions.definition_widths[index],
                    row_height=definitions.row_definition_heights[index],
                    row_width=definitions.row_definition_widths[index],
                    multi_record_enabled=frame.multi_record_enabled,
                    split_row_definitions=frame.multi_record_enabled,
                    row_top=definitions.row_definition_tops[index],
                ),
                sequence_width=placements[index].sequence_width,
                record_x=record_xs[index],
                multi_record_enabled=frame.multi_record_enabled,
                keep_definition_left_aligned=bool(canvas_config.keep_definition_left_aligned),
                definition_column_width=definitions.max_def_width,
                row_definition_width=definitions.row_definition_width,
                definition_gap=float(canvas_config.definition_column_gap),
                text_anchor=cfg.objects.definition.linear.text_anchor,
                definition_column_left=column_left,
            ),
            record_x=record_xs[index],
            features_origin_y=0.0,
        )
        for index, geometry in enumerate(tracks.record_feature_lane_geometries)
    ]
    return _LinearHorizontalFrame(
        final_translations=translations,
        definition_column_left=column_left,
        label_keep_clear=keep_clear,
    )


def _linear_feature_labels(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    frame: _LinearRowFrame,
    comparisons: _LinearComparisons,
    tracks: _LinearTrackPlan,
    horizontal: _LinearHorizontalFrame,
    feature_label_rows: FeatureLabelRows | None,
) -> list[list[LinearLabel]]:
    """Place each record's feature labels in its row's frame (OV-1188), in the rows an earlier assembly proposed
    when given; the row stacking reads their heights."""

    # The labels are placed once, at the records' drawn widths; the row stacking reads their heights.
    _unused_label_height, all_labels, _record_label_heights_above = _precalculate_label_dimensions(
        record_inputs.records,
        configurators.feature,
        canvas_config,
        render_context=tracks.render_context,
        precomputed_feature_dicts=tracks.record_feature_dicts,
        orthogroup_label_eligibility=comparisons.orthogroup_label_eligibility,
        feature_lane_geometries=tracks.record_feature_lane_geometries,
        sequence_widths=frame.sequence_widths,
        record_transforms=record_inputs.record_transforms,
        label_frames=frame.label_frames,
        label_keep_clear=horizontal.label_keep_clear,
    )
    if feature_label_rows is not None:
        all_labels = with_feature_label_rows(all_labels, feature_label_rows.rows)
    return all_labels


@dataclass(frozen=True, slots=True)
class _LinearVerticalPlans:
    """Each record's vertical plan and definition geometry, and its feature labels counted in its first row frame.

    ``labels_by_record`` are the labels before the annotation labels move them.
    """

    labels_by_record: list[list[LinearLabel]]
    record_vertical_plans: list[LinearRecordVerticalPlan]
    record_definition_geometries: list[_LinearRecordDefinitionGeometry]
    ribbon_clear_bands: dict[int, VerticalBand]
    # Each window of a windowed row: the feature lanes of its row (``_shared_by_row_windows``), whose top its
    # number and "//" length stand on (PD-OI-122).
    row_window_lanes: dict[int, VerticalBand]
    # Each record: where its ribbons attach (``_comparison_attachments``).
    comparison_attachments: dict[int, _ComparisonAttachment]


@dataclass(frozen=True, slots=True)
class _ComparisonAttachment:
    """Where a record's ribbons attach, from its axis (``band``: for a window of a windowed row, its row's), and the
    paint of its tracks other than features (``None`` without any)."""

    band: VerticalBand
    non_feature_band: VerticalBand | None


def _shared_by_row_windows(
    bands: Mapping[int, VerticalBand], row_windows: Sequence[LinearRowWindows]
) -> dict[int, VerticalBand]:
    """``bands`` by record index, each window of a windowed row given the union of its row's bands; other records
    keep their own. The windows of a row share one top and one bottom, so their numbers stand on one line and
    the ribbons to and from them end on one line, whichever strands each window's features take (OV-3374)."""

    shared = dict(bands)
    for windows in row_windows:
        row_band = union_vertical_bands([bands[index] for index in windows.record_indices])
        shared.update(dict.fromkeys(windows.record_indices, row_band))
    return shared


def _ribbon_clear_bands(
    frame: _LinearRowFrame,
    comparisons: _LinearComparisons,
    definition_geometries: Sequence[_LinearRecordDefinitionGeometry],
    marks_bands: Mapping[int, VerticalBand],
) -> dict[int, VerticalBand]:
    """The row text of each record of a shared row that a comparison ribbon reaches, by record index: the text
    above the record (OV-2712), and a window's number and "//" length (PD-OI-122: no row text is drawn over the
    comparison ribbons)."""

    if not frame.multi_record_enabled:
        return {}
    endpoints = {
        int(index)
        for comparison in comparisons.normalized
        for index in (comparison.query_record_index, comparison.subject_record_index)
    }
    clear: dict[int, VerticalBand] = {}
    for index in sorted(endpoints):
        bands = [
            band for band in (definition_geometries[index].local_band, marks_bands.get(index)) if band is not None
        ]
        if bands:
            clear[index] = bands[0] if len(bands) == 1 else bands[0].union(bands[1])
    return clear


def _linear_vertical_plans(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    frame: _LinearRowFrame,
    tracks: _LinearTrackPlan,
    definitions: _LinearDefinitionMetrics,
    depth: _LinearDepthData,
    feature_labels: list[list[LinearLabel]],
    comparisons: _LinearComparisons,
) -> _LinearVerticalPlans:
    """Plan each record's slots, labels and definition vertically. Sets the canvas's vertical offset."""

    records = record_inputs.records
    record_vertical_plans, record_definition_geometries = _build_linear_record_vertical_plans(
        records=records,
        layout=tracks.linear_track_layout,
        feature_dicts=tracks.record_feature_dicts,
        feature_geometries=tracks.record_feature_lane_geometries,
        labels_by_record=feature_labels,
        canvas_config=canvas_config,
        gc_config_of=depth.drawn_gc_config,
        skew_config=configurators.skew,
        record_depth_by_index=depth.record_depth_by_index,
        representative_depth_configs=depth.representative_depth_configs,
        annotation_track_layouts=tracks.annotation_track_layouts,
        definitions=definitions,
        multi_record_enabled=frame.multi_record_enabled,
        split_row_definitions=frame.multi_record_enabled,
        render_context=tracks.render_context,
        record_transforms=record_inputs.record_transforms,
    )
    # A window's marks, number and gap text belong to its record's body (PD-OI-122), on its row's lanes.
    row_window_lanes = _shared_by_row_windows(
        {
            index: record_vertical_plans[index].feature_lanes_band
            for windows in frame.row_windows
            for index in windows.record_indices
        },
        frame.row_windows,
    )
    marks_bands = {index: frame.row_window_style.marks_band(lanes.top_y) for index, lanes in row_window_lanes.items()}
    for index, marks in marks_bands.items():
        plan = record_vertical_plans[index]
        record_vertical_plans[index] = replace(
            plan,
            record_body_band=plan.record_body_band.union(marks),
            canvas_band=plan.canvas_band.union(marks),
        )
    # The row text of a record that a ribbon reaches belongs to its body: the ribbons end above it, and the rows
    # move apart by its height (OV-2712, PD-OI-122).
    ribbon_clear_bands = _ribbon_clear_bands(frame, comparisons, record_definition_geometries, marks_bands)
    for index, band in ribbon_clear_bands.items():
        plan = record_vertical_plans[index]
        record_vertical_plans[index] = replace(
            plan,
            record_body_band=plan.record_body_band.union(band),
            comparison_exclusion_band=plan.comparison_exclusion_band.union(band),
        )
    # The rows are spaced from a windowed row's shared edge, where its ribbons attach: a ribbon from its shallowest
    # window keeps the comparison height (OV-3374).
    row_comparison_bands = _shared_by_row_windows(
        {index: record_vertical_plans[index].comparison_exclusion_band for index in row_window_lanes},
        frame.row_windows,
    )
    for index, band in row_comparison_bands.items():
        plan = record_vertical_plans[index]
        record_vertical_plans[index] = replace(
            plan,
            record_body_band=plan.record_body_band.union(band),
            comparison_exclusion_band=band,
            canvas_band=plan.canvas_band.union(band),
        )
    # Every later reader counts a wrapped record's labels in one frame, where rows do not overlap.
    feature_slot_id = next((slot.id for slot in tracks.linear_track_layout.slots if slot.renderer == "features"), None)
    all_labels = [
        _labels_in_first_row_frame(labels, plan, feature_slot_id)
        for labels, plan in zip(feature_labels, record_vertical_plans, strict=True)
    ]

    canvas_config.vertical_offset = _linear_record_vertical_offset(
        record_vertical_plans,
        frame.rows_by_record,
        canvas_config,
    )
    return _LinearVerticalPlans(
        labels_by_record=all_labels,
        record_vertical_plans=record_vertical_plans,
        record_definition_geometries=record_definition_geometries,
        ribbon_clear_bands=ribbon_clear_bands,
        row_window_lanes=row_window_lanes,
        comparison_attachments=_comparison_attachments(
            record_vertical_plans, tracks, ribbon_clear_bands, frame.row_windows
        ),
    )


@dataclass(frozen=True, slots=True)
class _LinearGcFrames:
    """Each record's GC frame per dinucleotide a slot draws, and the frames of the request's dinucleotide."""

    by_dinucleotide: dict[str, list[DataFrame | None]]
    default: list[DataFrame | None]


def _linear_gc_frames(
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    tracks: _LinearTrackPlan,
) -> _LinearGcFrames:
    """Compute each record's GC frame once per dinucleotide its content and skew slots draw."""

    gc_config = configurators.gc
    records = record_inputs.records
    needed_nts = {
        slot_dinucleotide(slot.params, str(gc_config.dinucleotide))
        for slot in tracks.linear_track_layout.slots
        if slot.renderer in {"dinucleotide_content", "dinucleotide_skew"}
    }
    record_gc_dfs_by_nt = {
        nt: _precalculate_gc_dataframes(
            records,
            record_inputs.record_transforms,
            window=int(gc_config.window),
            step=int(gc_config.step),
            dinucleotide=nt,
            enabled=True,
        )
        for nt in sorted(needed_nts)
    }
    record_gc_dfs = record_gc_dfs_by_nt.get(
        str(gc_config.dinucleotide).upper(),
        [None for _ in records],
    )
    return _LinearGcFrames(by_dinucleotide=record_gc_dfs_by_nt, default=record_gc_dfs)


@dataclass(frozen=True, slots=True)
class _LinearLegend:
    """The legend's sources and its table (empty when the diagram draws no legend)."""

    sources: LegendSources
    table: LegendTable


def _linear_legend_table(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    blast_config: BlastMatchConfigurator,
    record_inputs: LinearRecordInputs,
    track_inputs: LinearTrackInputs,
    comparison_inputs: LinearComparisonInputs,
    comparisons: _LinearComparisons,
    tracks: _LinearTrackPlan,
    depth: _LinearDepthData,
) -> _LinearLegend:
    """Build the legend table (stage ``legend_table``) unless the legend is hidden."""

    profile = canvas_config.profile
    cfg = profile.config
    depth_config = configurators.depth
    render_context = tracks.render_context
    # The records share one scale, so the first record's layers name every record's rows.
    legend_sources = LegendSources.of(
        configurators.feature,
        configurators.gc,
        configurators.skew,
        feature_density_legend(tracks.record_feature_density[0] if tracks.record_feature_density else None, profile),
    )
    legend_table = LegendTable()
    if canvas_config.legend_position != "none":
        singleton_depth = render_context.depth_enabled and depth_track_count(track_inputs.record_depth_tracks) == 1
        comparison_identity_legend = pairwise_identity_legend_from_comparisons(blast_config, comparisons.matches)
        legend_table = render_stage(
            "legend_table",
            build_linear_legend_table,
            record_inputs.records,
            sources=legend_sources,
            comparison=ComparisonLegendSource.of(blast_config) if comparisons.has_blast and blast_config else None,
            pairwise_identity_legend=(
                comparison_identity_legend
                or comparison_inputs.pairwise_identity_legend
                or PairwiseIdentityLegend()
            ),
            depth=(
                DepthLegendSource.of_config(depth_config)
                if singleton_depth and depth_config is not None
                else None
            ),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            slots=legend_slot_sources(tracks.normalized_linear_track_slots),
            depth_tracks=(
                depth_legend_sources(depth.representative_depth_data) if render_context.depth_enabled else None
            ),
            annotations=tracks.resolved_annotations,
            edits=LegendRowEdits.of_config(cfg.objects.legends),
        )
    return _LinearLegend(sources=legend_sources, table=legend_table)


@dataclass(frozen=True, slots=True)
class _LinearRowStacking:
    """Each record's axis y before its final translation, and the collision bands the rows were stacked by.

    ``multi_record_plan`` is the solved layout when rows hold several
    records; ``boundary_gap_resolutions`` are the gaps between rows, in
    stacking order.
    """

    record_collision_bands: list[tuple[CollisionBand, ...]]
    axis_ys: list[float]
    multi_record_plan: LinearLayoutPlan | None
    boundary_gap_resolutions: list[AxisGapResolution]
    comparison_endpoint_gap_px: float


def _row_window_mark_bands(
    frame: _LinearRowFrame, vertical: _LinearVerticalPlans, record_index: int, record_x: float
) -> tuple[CollisionBand, ...]:
    """A window's marks beside it, in its insets (PD-OI-122), as a body band: rows keep clear of them and the
    figure bounds hold them."""

    lanes = vertical.row_window_lanes.get(record_index)
    if lanes is None:
        return ()
    placement = frame.horizontal_placements[record_index]
    marks = frame.row_window_style.marks_band(lanes.top_y)
    return (
        CollisionBand(
            "body",
            record_x - placement.left_inset,
            record_x + placement.sequence_width + placement.right_inset,
            marks.top_y,
            marks.bottom_y,
        ),
    )


def _active_collision_bands(
    bands: tuple[CollisionBand, ...], index: int, active_records: frozenset[int]
) -> tuple[CollisionBand, ...]:
    """A record's collision bands at one row boundary: its comparison band only when it draws a comparison there."""

    return tuple(band for band in bands if band.kind != "comparison" or index in active_records)


def _linear_row_stacking(
    canvas_config: LinearCanvasConfigurator,
    record_inputs: LinearRecordInputs,
    frame: _LinearRowFrame,
    comparisons: _LinearComparisons,
    definitions: _LinearDefinitionMetrics,
    horizontal: _LinearHorizontalFrame,
    vertical: _LinearVerticalPlans,
) -> _LinearRowStacking:
    """Stack the rows from the drawn collision bands, so the gap between rows is solved from the drawn boxes.

    Sets the canvas's plot width from the definition room it reserves.
    """

    cfg = canvas_config.profile.config
    records = record_inputs.records
    placements = frame.horizontal_placements
    plans = vertical.record_vertical_plans
    definition_reserve_width = (
        definitions.row_definition_width if frame.multi_record_enabled else definitions.max_def_width
    )
    canvas_config.configure_plot_width(definition_reserve_width)
    ordinary_row_gap = _ordinary_row_gap(canvas_config)
    definition_clear_gap = max(1.0, 0.5 * float(canvas_config.vertical_padding))
    comparison_endpoint_gap_px = (
        _COMPARISON_ENDPOINT_GAP_PX if comparisons.has_blast else 0.0
    )
    comparison_boundary_records = _comparison_records_by_boundary(
        comparisons.normalized,
        frame.rows_by_record,
    )
    current_y = max(canvas_config.vertical_offset, canvas_config.original_vertical_offset)
    # Each record's collision bands where the drawing puts them (OV-1400): at the record's final
    # translation, with the definition column's left and the text anchor the drawing uses.
    record_collision_bands = [
        _record_collision_bands(
            plan=plans[index],
            definition_geometry=vertical.record_definition_geometries[index],
            sequence_width=placements[index].sequence_width,
            record_x=placements[index].x + horizontal.final_translations[index][0],
            multi_record_enabled=frame.multi_record_enabled,
            text_anchor=cfg.objects.definition.linear.text_anchor,
            keep_definition_left_aligned=bool(canvas_config.keep_definition_left_aligned),
            definition_column_width=definitions.max_def_width,
            row_definition_width=definitions.row_definition_width,
            definition_gap=float(canvas_config.definition_column_gap),
            definition_column_left=horizontal.definition_column_left,
        )
        + _row_window_mark_bands(frame, vertical, index, placements[index].x + horizontal.final_translations[index][0])
        for index in range(len(records))
    ]
    record_offsets: list[float] = [0.0 for _record in records]
    boundary_gap_resolutions: list[AxisGapResolution] = []
    multi_record_plan = None
    if frame.multi_record_enabled:
        multi_record_plan = solve_linear_layout(
            [
                LinearRecordMeasurement(
                    record_index=index,
                    record_key=RecordKey(str(frame.record_keys[index])),
                    sequence_length=frame.axis_lengths[index],
                    left_inset=placements[index].left_inset,
                    right_inset=placements[index].right_inset,
                    top_extent=plans[index].canvas_top_extent,
                    bottom_extent=plans[index].canvas_bottom_extent,
                    comparison_top_extent=plans[index].comparison_top_extent,
                    comparison_bottom_extent=plans[index].comparison_bottom_extent,
                    # The drawn bands, local to the record's x, which the solve adds back.
                    collision_bands=tuple(
                        band.translate(x=-placements[index].x)
                        for band in record_collision_bands[index]
                    ),
                )
                for index, _record in enumerate(records)
            ],
            frame.horizontal_layout,
            first_axis_y=current_y,
            row_gap_px=ordinary_row_gap,
            comparison_height=float(canvas_config.configured_comparison_height),
            comparison_endpoint_gap_px=comparison_endpoint_gap_px,
            definition_clear_gap=definition_clear_gap,
            comparison_record_indices_by_boundary=comparison_boundary_records,
        )
        record_offsets = [
            multi_record_plan.placement_for_index(index).axis_y
            for index in range(len(records))
        ]
        boundary_gap_resolutions = list(multi_record_plan.row_gap_resolutions)
    else:
        ordered_rows = tuple(frame.ordered_record_indices)
        for row_position, record_index in enumerate(ordered_rows):
            record_offsets[record_index] = current_y
            if row_position >= len(ordered_rows) - 1:
                continue
            next_record_index = ordered_rows[row_position + 1]
            boundary = int(frame.rows_by_record[record_index])
            active_records = frozenset(
                comparison_boundary_records.get(boundary, ())
            )
            resolution = resolve_axis_gap(
                _active_collision_bands(record_collision_bands[record_index], record_index, active_records),
                _active_collision_bands(record_collision_bands[next_record_index], next_record_index, active_records),
                ordinary_row_gap=ordinary_row_gap,
                comparison_height=(
                    float(canvas_config.configured_comparison_height)
                    + (2.0 * comparison_endpoint_gap_px)
                ),
                definition_clear_gap=definition_clear_gap,
                boundary_has_comparison=bool(active_records),
            )
            boundary_gap_resolutions.append(resolution)
            current_y += resolution.axis_gap
    return _LinearRowStacking(
        record_collision_bands=record_collision_bands,
        axis_ys=record_offsets,
        multi_record_plan=multi_record_plan,
        boundary_gap_resolutions=boundary_gap_resolutions,
        comparison_endpoint_gap_px=comparison_endpoint_gap_px,
    )


@dataclass(frozen=True, slots=True)
class _LinearPlacements:
    """Each record's final placement (translated, with its comparison endpoints) and its feature attachments.

    A placement's ``x`` and ``axis_y`` are where the record is drawn;
    ``alignment_placements`` are the placements before the final
    translations and the composition. ``length_bar_offset_y`` is where the
    scale is drawn when ``draw_length_bar``.
    """

    record_placements: dict[int, LinearRecordPlacement]
    alignment_placements: tuple[LinearRecordPlacement, ...]
    feature_attachments: dict[int, dict[str, LinearFeatureAttachment]]
    draw_length_bar: bool
    length_bar_offset_y: float


def _comparison_attachments(
    plans: Sequence[LinearRecordVerticalPlan],
    tracks: _LinearTrackPlan,
    clear_bands: Mapping[int, VerticalBand],
    row_windows: Sequence[LinearRowWindows],
) -> dict[int, _ComparisonAttachment]:
    """Where each record's ribbons attach: its axis, its features by their occupied lanes, its other track paint
    and the row text its ribbons keep clear of (``clear_bands``). The windows of a windowed row attach on their row's
    band (``_shared_by_row_windows``), the band the rows are spaced from."""

    bands: dict[int, VerticalBand] = {}
    non_feature_bands: dict[int, VerticalBand | None] = {}
    for record_index, plan in enumerate(plans):
        band = plan.axis_band
        non_feature_band: VerticalBand | None = None
        for plan_slot in plan.slots:
            paint = plan_slot.paint_band
            if paint is None:
                continue
            if plan_slot.renderer == "features":
                paint = tracks.record_feature_lane_geometries[record_index].occupied_band.translate(plan_slot.origin_y)
            elif not plan_slot.spans_slots:
                # A span's band comes from the slots it covers, counted here or as lanes.
                non_feature_band = paint if non_feature_band is None else non_feature_band.union(paint)
            band = band.union(paint)
        clear_band = clear_bands.get(record_index)
        bands[record_index] = band if clear_band is None else band.union(clear_band)
        non_feature_bands[record_index] = non_feature_band
    return {
        record_index: _ComparisonAttachment(band, non_feature_bands[record_index])
        for record_index, band in _shared_by_row_windows(bands, row_windows).items()
    }


def _attach_comparison_endpoints(
    record_placements: dict[int, LinearRecordPlacement],
    tracks: _LinearTrackPlan,
    vertical: _LinearVerticalPlans,
    comparison_endpoint_gap_px: float,
) -> dict[int, dict[str, LinearFeatureAttachment]]:
    """Attach each record's comparison endpoints to its track paint, replacing its placement in
    ``record_placements``, and return the attachments of features that moved off the first lane.

    Row spacing reserves labels, but ribbon endpoints attach to track paint
    (``_LinearVerticalPlans.comparison_attachments``); this runs after both
    layout paths have fixed every slot and record. A feature's attachment is
    kept as its offset from the record's anchors, so moving rows leaves it as
    it is. The ribbons of a windowed row's windows end on the row's edge, so
    they follow no feature.
    """

    feature_attachments: dict[int, dict[str, LinearFeatureAttachment]] = {}
    for record_index, placement in record_placements.items():
        plan = vertical.record_vertical_plans[record_index]
        attachment_band = vertical.comparison_attachments[record_index].band
        non_feature_band = vertical.comparison_attachments[record_index].non_feature_band
        clear_band = vertical.ribbon_clear_bands.get(record_index)
        record_placements[record_index] = replace(
            placement,
            comparison_top_y=(
                placement.axis_y + attachment_band.top_y - comparison_endpoint_gap_px
            ),
            comparison_bottom_y=(
                placement.axis_y + attachment_band.bottom_y + comparison_endpoint_gap_px
            ),
        )
        # Feature-associated matches follow their endpoint features. A moved
        # feature must not shorten matches attached to another occupied lane.
        features = tracks.record_feature_dicts[record_index]
        if record_index in vertical.row_window_lanes:
            continue
        if not any(feature.placement is not None and feature.placement.requested_target is not None
                   for feature in features.values()):
            continue
        if not any(feature.placement is not None and feature.placement.level > 0
                   for feature in features.values()):
            continue
        for plan_slot in plan.slots:
            if plan_slot.renderer != "features" or plan_slot.paint_band is None:
                continue
            geometry = tracks.record_feature_lane_geometries[record_index]
            for feature in features.values():
                if feature.source_feature_index is None:
                    continue  # by_source_index only holds integer source indexes
                dom_id = tracks.feature_dom_index.by_source_index.get((record_index, feature.source_feature_index))
                if dom_id is None:
                    continue
                lane = geometry.lane_for(strand=feature.strand, track_id=feature.feature_track_id,
                                         separate_strands=tracks.render_context.profile.strandedness)
                band = lane.band.translate(plan_slot.origin_y)
                if non_feature_band is not None:
                    band = band.union(non_feature_band)
                feature_attachments.setdefault(record_index, {})[dom_id] = LinearFeatureAttachment(
                    # A ribbon that meets a record below its text ends above the text (OV-2712).
                    0.0 if clear_band is not None else comparison_frame_px(band.top_y - attachment_band.top_y),
                    comparison_frame_px(band.bottom_y - attachment_band.bottom_y),
                )
    return feature_attachments


def _linear_vertical_placement(
    canvas_config: LinearCanvasConfigurator,
    frame: _LinearRowFrame,
    tracks: _LinearTrackPlan,
    horizontal: _LinearHorizontalFrame,
    vertical: _LinearVerticalPlans,
    stacking: _LinearRowStacking,
) -> _LinearPlacements:
    """Place each record at its stacked axis y and final translation, and attach its comparison endpoints.

    Sets the canvas's height below the final record and its total height
    from the painted content.
    """

    cfg = canvas_config.profile.config
    record_count = len(frame.axis_lengths)
    record_plans = vertical.record_vertical_plans
    # Every record and window of the figure draws at one bp-per-pixel scale unless lengths are normalized
    # (``solve_linear_horizontal_layout``), so several records in a row keep the scale (OV-2714).
    draw_length_bar = (
        cfg.objects.scale.show
        and not canvas_config.normalize_length
        and not tracks.render_context.axis_ruler_enabled
    )
    record_placements: dict[int, LinearRecordPlacement] = {}
    for record_index in range(record_count):
        if stacking.multi_record_plan is not None:
            record_placements[record_index] = stacking.multi_record_plan.placement_for_index(
                record_index
            )
            continue
        record_plan = record_plans[record_index]
        record_placements[record_index] = place_record_vertically(
            frame.horizontal_placements[record_index],
            axis_y=stacking.axis_ys[record_index],
            top_extent=record_plan.canvas_top_extent,
            bottom_extent=record_plan.canvas_bottom_extent,
            comparison_top_extent=record_plan.comparison_top_extent,
            comparison_bottom_extent=record_plan.comparison_bottom_extent,
            comparison_endpoint_gap=stacking.comparison_endpoint_gap_px,
        )
    # Keep raw layout placements available to the existing alignment helper.
    # They exclude canonical translations and composition/viewBox fitting.
    alignment_placements = tuple(record_placements.values())
    for record_index, (translation_x, translation_y) in enumerate(
        horizontal.final_translations
    ):
        placement = record_placements[record_index]
        record_placements[record_index] = replace(
            placement,
            x=placement.x + translation_x,
            axis_y=placement.axis_y + translation_y,
            comparison_top_y=placement.comparison_top_y + translation_y,
            comparison_bottom_y=placement.comparison_bottom_y + translation_y,
        )
    painted_content_bottom = max(
        (
            float(record_placements[index].axis_y) + plan.canvas_band.bottom_y
            for index, plan in enumerate(record_plans)
        ),
        default=0.0,
    )
    canvas_config.height_below_final_record = (
        painted_content_bottom + 4 * canvas_config.vertical_padding
    )
    length_bar_offset_y = float(canvas_config.height_below_final_record)
    canvas_config.total_height = max(
        1.0,
        painted_content_bottom + float(canvas_config.vertical_padding),
    )
    feature_attachments = _attach_comparison_endpoints(
        record_placements,
        tracks,
        vertical,
        stacking.comparison_endpoint_gap_px,
    )
    return _LinearPlacements(
        record_placements=record_placements,
        alignment_placements=alignment_placements,
        feature_attachments=feature_attachments,
        draw_length_bar=bool(draw_length_bar),
        length_bar_offset_y=length_bar_offset_y,
    )


@dataclass(frozen=True, slots=True)
class _LinearLengthBar:
    """The scale's content and where it is drawn (``group`` is ``None`` when the figure draws no scale)."""

    group: GroupContent | None
    offset_x: float
    offset_y: float


def _linear_length_bar(
    canvas_config: LinearCanvasConfigurator, frame: _LinearRowFrame, placements: _LinearPlacements
) -> _LinearLengthBar:
    """Build the scale (stage ``length_bar_content``) under the records when the figure draws one.

    Sets the canvas's total height to include it.
    """

    if not placements.draw_length_bar:
        return _LinearLengthBar(group=None, offset_x=0.0, offset_y=placements.length_bar_offset_y)
    cfg = canvas_config.profile.config
    record_placements = placements.record_placements
    min_record_x = min(
        placement.x for placement in record_placements.values()
    )
    max_record_x = max(
        placement.x + placement.sequence_width
        for placement in record_placements.values()
    )
    widest_record = max(placement.sequence_width for placement in record_placements.values())
    length_bar_group = render_stage(
        "length_bar_content",
        build_group_content,
        LengthBarGroup,
        # The longest record's drawn width: the bar maps bp to px as that record does.
        widest_record,
        canvas_config.longest_genome,
        canvas_config.length_param,
        scale=cfg.objects.scale,
        font_family=cfg.objects.text.font_family,
        # A ruler across a row of several records would number their bases on from one record to the next;
        # it spans the longest record (PD-OI-122, OV-2714).
        ruler_width=widest_record if frame.multi_record_enabled else max_record_x - min_record_x,
    )
    length_bar_bottom = placements.length_bar_offset_y + float(
        _content_bounds(length_bar_group).max_y
    )
    canvas_config.total_height = max(
        1.0,
        length_bar_bottom + float(canvas_config.vertical_padding),
    )
    return _LinearLengthBar(
        group=length_bar_group, offset_x=min_record_x, offset_y=placements.length_bar_offset_y
    )


@dataclass(frozen=True, slots=True)
class _LinearPlan:
    """Everything the Linear figure is drawn from, solved before its canvas exists, one field per phase."""

    frame: _LinearRowFrame
    comparisons: _LinearComparisons
    tracks: _LinearTrackPlan
    definitions: _LinearDefinitionMetrics
    depth: _LinearDepthData
    horizontal: _LinearHorizontalFrame
    vertical: _LinearVerticalPlans
    gc_frames: _LinearGcFrames
    legend: _LinearLegend
    stacking: _LinearRowStacking
    placements: _LinearPlacements
    length_bar: _LinearLengthBar


def _linear_canvas(
    canvas_config: LinearCanvasConfigurator,
    blast_config: BlastMatchConfigurator,
    record_inputs: LinearRecordInputs,
    plan: _LinearPlan,
) -> tuple[Drawing, int]:
    """Create the canvas and draw the scale and the comparisons on it.

    Returns the canvas and the index of its first element that the
    composition moves with the records.
    """

    canvas: Drawing = canvas_config.create_svg_canvas()
    primary_target_start = len(getattr(canvas, "elements", []))
    if plan.length_bar.group is not None:
        canvas = add_length_bar_on_linear_canvas(
            canvas,
            canvas_config,
            plan.length_bar.group,
            offset_x=plan.length_bar.offset_x,
        )
    if plan.comparisons.has_blast:
        canvas = add_explicit_comparisons_on_linear_canvas(
            canvas,
            plan.comparisons.rows.comparisons,
            canvas_config,
            blast_config,
            record_inputs.records,
            plan.placements.record_placements,
            feature_dom_index=plan.tracks.feature_dom_index,
            record_transforms=record_inputs.record_transforms,
            record_keys=plan.frame.record_keys,
            feature_attachments=plan.placements.feature_attachments,
        )
    return canvas, primary_target_start


@dataclass(frozen=True, slots=True)
class _LinearAnnotationLabels:
    """The annotation labels by record and slot, and each record's feature labels after the moves made for them.

    Later phases read ``labels_by_record`` here, not the vertical plans'.
    """

    labels_by_record: list[list[LinearLabel]]
    annotation_labels_by_record: list[dict[str, tuple[LinearAnnotationLabel, ...]]]
    label_fit: LinearLabelFit


def _linear_annotation_labels(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    plan: _LinearPlan,
    feature_label_rows: FeatureLabelRows | None,
) -> _LinearAnnotationLabels:
    """Place the annotation labels along their spans, off the feature labels and leaders (OV-845), and count
    the drawn labels' fit, which the API checks and offers rows from (LP-D48, LP-D52)."""

    cfg = canvas_config.profile.config
    records = record_inputs.records
    tracks = plan.tracks
    record_placements = plan.placements.record_placements
    record_collision_bands = plan.stacking.record_collision_bands
    feature_config = configurators.feature
    # Each record's drawn definition text, in its feature labels' coordinates (OV-1395).
    definition_label_boxes = [
        _definition_label_boxes(
            record_collision_bands[index],
            record_x=record_placements[index].x,
            features_origin_y=plan.vertical.record_vertical_plans[index].features_origin_y,
        )
        for index in range(len(records))
    ]
    placed_annotation_labels, all_labels = _place_linear_annotation_labels_by_record(
        records=records,
        axis_lengths=plan.frame.axis_lengths,
        placements=record_placements,
        plans=plan.vertical.record_vertical_plans,
        collision_bands=record_collision_bands,
        feature_labels=plan.vertical.labels_by_record,
        feature_dicts=tracks.record_feature_dicts,
        feature_geometries=tracks.record_feature_lane_geometries,
        separate_strands=tracks.render_context.profile.strandedness,
        feature_stroke_width=max(float(feature_config.block_stroke_width), float(feature_config.line_stroke_width)),
        annotation_track_layouts=tracks.annotation_track_layouts,
        font_family=str(cfg.objects.text.font_family),
        label_style=tracks.render_context.profile.label_style,
        gap=float(cfg.labels.spacing.linear),
        # One assembly more at most: one that places proposed rows proposes none.
        new_rows=feature_label_rows is None,
        label_frames=plan.frame.label_frames,
        definition_boxes=definition_label_boxes,
    )
    label_fit = linear_label_fit(
        [(str(record.id), labels) for record, labels in zip(records, all_labels)],
        tracks.label_font_size,
        [placed.region_labels for placed in placed_annotation_labels],
        FeatureLabelRows(
            tuple(row for placed in placed_annotation_labels for row in (*placed.rows, *placed.new_rows)),
            sum(placed.new_row_overlaps for placed in placed_annotation_labels),
            tuple(
                (record_index, label)
                for record_index, placed in enumerate(placed_annotation_labels)
                for label in placed.proposers
            ),
        )
        if any(placed.new_rows for placed in placed_annotation_labels)
        else None,
        frames=plan.frame.label_frames,
        definition_boxes=definition_label_boxes,
    )
    return _LinearAnnotationLabels(
        labels_by_record=all_labels,
        annotation_labels_by_record=[placed.by_slot for placed in placed_annotation_labels],
        label_fit=label_fit,
    )


@dataclass(frozen=True, slots=True)
class _LinearRecordDrawing:
    """What drawing one record's slots, axis and definition reads.

    Values of the whole figure come from ``plan``. Each of this record's
    values is one field here, taken from the plan by the record's index
    once, and the drawers read it only here. ``placement`` is where the
    record is drawn (its ``x`` and ``axis_y``); ``translation`` is its final
    translation; ``axis_row_offsets`` are its wrap rows' offsets (empty when
    it does not wrap); ``gc_df_by_dinucleotide`` holds its GC frame per
    dinucleotide a slot draws, and ``default_gc_df`` the frame of the
    request's dinucleotide.
    """

    canvas_config: LinearCanvasConfigurator
    configurators: DiagramConfigurators
    plan: _LinearPlan
    record: SeqRecord
    record_index: int
    record_count: int
    record_key: str
    record_transform: RecordDisplayTransform | None
    presentation: RecordPresentation | None
    axis_length: int
    placement: LinearRecordPlacement
    translation: tuple[float, float]
    vertical_plan: LinearRecordVerticalPlan
    definition_geometry: _LinearRecordDefinitionGeometry
    local_line_kinds: frozenset[str] | None
    row_line_kinds: frozenset[str]
    feature_layers: FeatureBuildResult
    feature_lane_geometry: LinearFeatureLaneGeometry
    feature_labels: list[LinearLabel]
    annotation_labels: Mapping[str, tuple[LinearAnnotationLabel, ...]]
    depth_by_index: Mapping[int, DepthTrackData]
    feature_density: Mapping[str, FeatureDensityLayer]
    gc_df_by_dinucleotide: Mapping[str, DataFrame | None]
    default_gc_df: DataFrame | None
    orthogroup_label_member_ids: set[str | int] | None
    orthogroup_label_top_member_ids: set[str | int] | None
    record_group_id: str
    axis_row_offsets: tuple[float, ...]
    # The record whose identity the heading of this record's windowed row carries, when this record leads the row.
    row_heading_record: SeqRecord | None


def _linear_record_drawing(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    record_inputs: LinearRecordInputs,
    plan: _LinearPlan,
    labels: _LinearAnnotationLabels,
    record_index: int,
) -> _LinearRecordDrawing:
    """Collect what drawing one record reads, by its index."""

    records = record_inputs.records
    record = records[record_index]
    definitions = plan.definitions
    vertical_plan = plan.vertical.record_vertical_plans[record_index]
    windows = plan.frame.row_windows_of(record_index)
    member_ids, top_member_ids = orthogroup_label_sets_for_record(
        plan.comparisons.orthogroup_label_eligibility if plan.comparisons.label_scope == "orthogroup_top" else None,
        record_index,
    )
    return _LinearRecordDrawing(
        canvas_config=canvas_config,
        configurators=configurators,
        plan=plan,
        record=record,
        record_index=record_index,
        record_count=len(records),
        record_key=str(plan.frame.record_keys[record_index]),
        record_transform=(
            record_inputs.record_transforms[record_index] if record_inputs.record_transforms is not None else None
        ),
        presentation=plan.frame.presentations[record_index],
        axis_length=plan.frame.axis_lengths[record_index],
        placement=plan.placements.record_placements[record_index],
        translation=plan.horizontal.final_translations[record_index],
        vertical_plan=vertical_plan,
        definition_geometry=plan.vertical.record_definition_geometries[record_index],
        local_line_kinds=(
            definitions.local_definition_line_kinds[record_index]
            if definitions.local_definition_line_kinds is not None
            else None
        ),
        row_line_kinds=definitions.row_definition_line_kinds[record_index],
        feature_layers=plan.tracks.record_feature_layers[record_index],
        feature_lane_geometry=plan.tracks.record_feature_lane_geometries[record_index],
        feature_labels=labels.labels_by_record[record_index],
        annotation_labels=labels.annotation_labels_by_record[record_index],
        depth_by_index=plan.depth.record_depth_by_index[record_index],
        feature_density=plan.tracks.record_feature_density[record_index],
        gc_df_by_dinucleotide={nt: frames[record_index] for nt, frames in plan.gc_frames.by_dinucleotide.items()},
        default_gc_df=plan.gc_frames.default[record_index],
        orthogroup_label_member_ids=member_ids,
        orthogroup_label_top_member_ids=top_member_ids,
        record_group_id=record_group_svg_id(
            record.id,
            mode="linear",
            record_index=record_index,
            record_count=len(records),
        ),
        axis_row_offsets=vertical_plan.row_offsets if vertical_plan.wrap_rows else (),
        row_heading_record=(
            records[windows.heading_record_index]
            if windows is not None and windows.leading_record_index == record_index
            else None
        ),
    )


class _LinearSlotDrawing(NamedTuple):
    """One slot of one record as the drawing reads it: the slot, where the record's plan puts it, how far it
    lies on each wrap row below the first, and its DOM id."""

    slot: LinearResolvedTrack
    resolved: LinearResolvedSlot
    row_shifts: tuple[float, ...]
    group_id: str


def _linear_slot_group_id(drawing: _LinearRecordDrawing, slot: LinearResolvedTrack) -> str:
    """The DOM id of one record's slot group, the same whether the request listed the slots or took the
    default ones, so a figure drawn without a stack and the same stack written out carry the same ids
    (OV-3012). It is the Circular slot id (``track_slot_svg_id``), namespaced by record when the figure has
    more than one."""

    group_id = track_slot_svg_id(str(slot.id), renderer=str(slot.renderer))
    if int(drawing.record_count) <= 1:
        return group_id
    return instance_svg_id(group_id, ("linear-record", int(drawing.record_index), int(drawing.record_count)))


def _draw_linear_feature_slot(canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing) -> Drawing:
    """Draw one record's features, labels and axis in its feature slot."""

    plan = drawing.plan
    translation_x, translation_y = drawing.translation
    return add_record_group(
        canvas,
        drawing.record,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        drawing.configurators.feature,
        feature_layers=drawing.feature_layers,
        render_context=plan.tracks.render_context,
        feature_dom_index=plan.tracks.feature_dom_index,
        precalculated_labels=drawing.feature_labels,
        label_font_size=plan.tracks.label_font_size,
        orthogroup_label_member_ids=drawing.orthogroup_label_member_ids,
        orthogroup_label_top_member_ids=drawing.orthogroup_label_top_member_ids,
        record_index=drawing.record_index,
        record_count=drawing.record_count,
        group_id=slot.group_id,
        feature_slot=slot.slot,
        placement=drawing.placement,
        multi_record_layout=plan.frame.multi_record_enabled,
        record_local_ruler=plan.frame.multi_record_enabled,
        feature_offset_y=slot.resolved.origin_y,
        feature_lane_geometry=drawing.feature_lane_geometry,
        record_transform=drawing.record_transform,
        record_translation_x=translation_x,
        record_translation_y=translation_y,
        axis_row_offsets=drawing.axis_row_offsets,
        feature_row_shifts=slot.row_shifts,
        group_key=GroupKey(drawing.record_key, GroupKind.FEATURES, str(slot.slot.id)),
    )


def _draw_linear_annotation_slot(canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing) -> Drawing:
    """Draw one record's region annotations in an annotation slot (stage ``annotation_content``)."""

    slot_id = slot.slot.id
    annotation_layout = drawing.plan.tracks.annotation_track_layouts.get(slot_id)
    if annotation_layout is None:
        return canvas
    annotation_content = render_stage(
        "annotation_content",
        linear_annotation_track_content,
        annotation_layout,
        record_id=str(drawing.record.id),
        record_index=drawing.record_index,
        record_length=drawing.axis_length,
        bar_length_px=float(drawing.placement.sequence_width),
        y_offset_px=0.0,
        frame=linear_annotation_slot_frame(drawing.vertical_plan, slot_id),
        font_family=str(drawing.canvas_config.profile.config.objects.text.font_family),
        params=annotation_track_params_from_mapping(slot.slot.params),
        labels=drawing.annotation_labels.get(slot_id, ()),
        group_id=slot.group_id,
        slot_id=str(slot_id),
    )
    annotation_group = place_content(
        canvas,
        annotation_content,
        GroupKey(drawing.record_key, GroupKind.ANNOTATIONS, str(slot_id)),
    )
    annotation_group.translate(
        drawing.placement.x + drawing.canvas_config.horizontal_offset,
        drawing.placement.axis_y + float(slot.resolved.origin_y),
    )
    canvas.add(annotation_group)
    return canvas


def _draw_linear_depth_slot(canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing) -> Drawing:
    """Draw one record's depth track in a depth slot, when the record has that track."""

    depth_track = drawing.depth_by_index.get(_depth_slot_track_index(slot.slot))
    if depth_track is None:
        return canvas
    return add_depth_group(
        canvas,
        drawing.record,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        depth_track.config,
        depth_df=depth_track.df,
        group_id=slot.group_id,
        slot_id=str(slot.slot.id),
        slot_renderer=str(slot.slot.renderer),
        track_height=slot.slot.height,
        track_offset_y=slot.resolved.origin_y,
        sequence_width=drawing.placement.sequence_width,
        record_transform=drawing.record_transform,
        row_shifts=slot.row_shifts,
        group_key=GroupKey(drawing.record_key, GroupKind.DEPTH, str(slot.slot.id)),
    )


def _draw_linear_feature_density_slot(
    canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing
) -> Drawing:
    """Draw one record's feature density track in its slot: a histogram as a depth area, a heatmap as a band."""

    layer = slot_feature_density_layer(str(slot.slot.id), slot.slot, drawing.feature_density)
    if layer is None:
        return canvas
    profile = drawing.canvas_config.profile
    group_key = GroupKey(drawing.record_key, GroupKind.DEPTH, str(slot.slot.id))
    slot_id = str(slot.slot.id)
    slot_renderer = str(slot.slot.renderer)
    if layer.track.style == "histogram":
        frame, histogram_config = feature_density_histogram_source(layer, profile)
        return add_depth_group(
            canvas,
            drawing.record,
            drawing.placement.axis_y,
            drawing.placement.x,
            drawing.canvas_config,
            histogram_config,
            depth_df=frame,
            group_id=slot.group_id,
            slot_id=slot_id,
            slot_renderer=slot_renderer,
            track_height=slot.slot.height,
            track_offset_y=slot.resolved.origin_y,
            sequence_width=drawing.placement.sequence_width,
            record_transform=drawing.record_transform,
            row_shifts=slot.row_shifts,
            group_key=group_key,
        )
    heatmap = FeatureDensityHeatmapGroup(
        layer,
        len(drawing.record.seq),
        float(slot.slot.height),
        feature_density_color(layer, profile),
        sequence_width=drawing.placement.sequence_width,
        group_id=slot.group_id,
        record_transform=drawing.record_transform,
        row_shifts=slot.row_shifts,
    )
    return add_feature_density_heatmap_group(
        canvas,
        heatmap,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        track_offset_y=slot.resolved.origin_y,
        group_key=group_key,
        slot_id=slot_id,
        slot_renderer=slot_renderer,
    )


def _slot_gc_df(drawing: _LinearRecordDrawing, nt: str) -> DataFrame | None:
    """The record's GC frame of one dinucleotide, or of the request's dinucleotide when no slot asked for it."""

    return drawing.gc_df_by_dinucleotide.get(nt, drawing.default_gc_df)


def _draw_linear_gc_content_slot(canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing) -> Drawing:
    """Draw one record's dinucleotide content in its slot."""

    nt = slot_dinucleotide(slot.slot.params, str(drawing.configurators.gc.dinucleotide))
    return add_gc_content_group(
        canvas,
        drawing.record,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        drawing.plan.depth.drawn_gc_config(slot.slot.params, drawing.depth_by_index),
        gc_df=_slot_gc_df(drawing, nt),
        track_height=slot.slot.height,
        track_offset_y=slot.resolved.origin_y,
        group_id=slot.group_id,
        slot_id=str(slot.slot.id),
        slot_renderer=str(slot.slot.renderer),
        sequence_width=drawing.placement.sequence_width,
        record_transform=drawing.record_transform,
        row_shifts=slot.row_shifts,
        group_key=GroupKey(drawing.record_key, GroupKind.GC_CONTENT, str(slot.slot.id)),
    )


def _draw_linear_gc_skew_slot(canvas: Drawing, drawing: _LinearRecordDrawing, slot: _LinearSlotDrawing) -> Drawing:
    """Draw one record's dinucleotide skew in its slot."""

    skew_config = drawing.configurators.skew
    nt = slot_dinucleotide(slot.slot.params, str(skew_config.dinucleotide))
    return add_gc_skew_group(
        canvas,
        drawing.record,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        _slot_skew_config(skew_config, slot.slot, nt),
        gc_df=_slot_gc_df(drawing, nt),
        track_height=slot.slot.height,
        track_offset_y=slot.resolved.origin_y,
        group_id=slot.group_id,
        slot_id=str(slot.slot.id),
        slot_renderer=str(slot.slot.renderer),
        sequence_width=drawing.placement.sequence_width,
        record_transform=drawing.record_transform,
        row_shifts=slot.row_shifts,
        group_key=GroupKey(drawing.record_key, GroupKind.GC_SKEW, str(slot.slot.id)),
    )


_LinearSlotDrawer = Callable[[Drawing, _LinearRecordDrawing, _LinearSlotDrawing], Drawing]
# Each drawn slot renderer's drawer; a spacer and an unknown renderer draw nothing.
_LINEAR_SLOT_DRAWERS: dict[str, _LinearSlotDrawer] = {
    "features": _draw_linear_feature_slot,
    "annotations": _draw_linear_annotation_slot,
    "depth": _draw_linear_depth_slot,
    "feature_density": _draw_linear_feature_density_slot,
    "dinucleotide_content": _draw_linear_gc_content_slot,
    "dinucleotide_skew": _draw_linear_gc_skew_slot,
}


def _draw_linear_region_label_fronts(canvas: Drawing, drawing: _LinearRecordDrawing) -> None:
    """Draw one record's backed region labels in front of its leader lines (LP-D82 RGN-FRONT)."""

    tracks = drawing.plan.tracks
    for slot in tracks.linear_track_layout.slots:
        annotation_layout = tracks.annotation_track_layouts.get(slot.id)
        front = (
            linear_region_label_front_content(
                annotation_layout,
                record_index=drawing.record_index,
                font_family=str(drawing.canvas_config.profile.config.objects.text.font_family),
                labels=drawing.annotation_labels.get(slot.id, ()),
                group_id=_linear_slot_group_id(drawing, slot),
            )
            if slot.renderer == "annotations" and annotation_layout is not None
            else None
        )
        if front is None:
            continue
        front_group = place_content(
            canvas, front, GroupKey(drawing.record_key, GroupKind.LABEL_TEXT, str(slot.id))
        )
        front_group.translate(
            drawing.placement.x + drawing.canvas_config.horizontal_offset,
            drawing.placement.axis_y + float(drawing.vertical_plan.slot_by_id(slot.id).origin_y),
        )
        canvas.add(front_group)


def _draw_linear_axis_only(canvas: Drawing, drawing: _LinearRecordDrawing) -> None:
    """Draw one record's axis without features, for a figure whose slots hold no feature track."""

    plan = drawing.plan
    translation_x, translation_y = drawing.translation
    add_record_group(
        canvas,
        drawing.record,
        drawing.placement.axis_y,
        drawing.placement.x,
        drawing.canvas_config,
        drawing.configurators.feature,
        feature_layers=drawing.feature_layers,
        render_context=plan.tracks.render_context,
        feature_dom_index=plan.tracks.feature_dom_index,
        precalculated_labels=None,
        label_font_size=plan.tracks.label_font_size,
        orthogroup_label_member_ids=drawing.orthogroup_label_member_ids,
        orthogroup_label_top_member_ids=drawing.orthogroup_label_top_member_ids,
        record_index=drawing.record_index,
        record_count=drawing.record_count,
        group_id=drawing.record_group_id,
        feature_slot=None,
        placement=drawing.placement,
        multi_record_layout=plan.frame.multi_record_enabled,
        record_local_ruler=plan.frame.multi_record_enabled,
        feature_lane_geometry=drawing.feature_lane_geometry,
        record_transform=drawing.record_transform,
        record_translation_x=translation_x,
        record_translation_y=translation_y,
        axis_row_offsets=drawing.axis_row_offsets,
        group_key=GroupKey(drawing.record_key, GroupKind.AXIS),
    )


def _draw_linear_record_definition(canvas: Drawing, drawing: _LinearRecordDrawing) -> None:
    """Draw one record's definition beside its tracks, or its own and its row's lines when rows hold several records."""

    plan = drawing.plan
    definitions = plan.definitions
    multi_record_enabled = plan.frame.multi_record_enabled
    geometry = drawing.definition_geometry
    offset_y = drawing.placement.axis_y
    record_index = drawing.record_index
    add_record_definition_group(
        canvas,
        drawing.record,
        offset_y,
        drawing.placement.x,
        drawing.canvas_config,
        definitions.max_def_width,
        group_id=definition_group_svg_id(
            drawing.record.id,
            mode="linear",
            record_index=record_index,
            record_count=drawing.record_count,
        ),
        placement=drawing.placement,
        row_definition_width=definitions.row_definition_width,
        definition_center_y=(
            offset_y + geometry.row_center_y
            if multi_record_enabled
            and geometry.row_center_y is not None
            else (
                offset_y + geometry.local_center_y
                if not multi_record_enabled
                and geometry.local_center_y is not None
                else None
            )
        ),
        definition_header_center_y=(
            offset_y + geometry.local_center_y
            if multi_record_enabled
            and geometry.local_center_y is not None
            else None
        ),
        multi_record_layout=multi_record_enabled,
        local_line_kinds=drawing.local_line_kinds,
        row_line_kinds=drawing.row_line_kinds,
        record_index=record_index,
        record_count=drawing.record_count,
        record_presentation=drawing.presentation,
        definition_column_left=plan.horizontal.definition_column_left,
        record_key=drawing.record_key,
    )
    windows = plan.frame.row_windows_of(record_index)
    if (
        windows is not None
        and drawing.row_heading_record is not None
        and geometry.row_center_y is not None
    ):
        add_row_window_heading_group(
            canvas,
            _row_window_heading_content(
                drawing.canvas_config, drawing.row_heading_record, drawing.record_count, plan.frame, windows
            ),
            drawing.canvas_config,
            drawing.placement,
            row_definition_width=definitions.row_definition_width,
            # The heading's name line stands on the features' center, where its band starts ``row_top`` above.
            center_y=offset_y + geometry.row_band.top_y - float(definitions.row_definition_tops[record_index] or 0.0)
            if geometry.row_band is not None
            else offset_y,
            definition_column_left=plan.horizontal.definition_column_left,
            record_key=drawing.record_key,
        )


def _draw_linear_row_window_marks(canvas: Drawing, canvas_config: LinearCanvasConfigurator, plan: _LinearPlan) -> None:
    """Draw each windowed row's marks above its records (PD-OI-122)."""

    for windows in plan.frame.row_windows:
        places = tuple(
            RowWindowPlace(
                record_index=index,
                x=float(plan.placements.record_placements[index].x),
                width=float(plan.placements.record_placements[index].sequence_width),
                axis_y=float(plan.placements.record_placements[index].axis_y),
                feature_top_y=float(plan.vertical.row_window_lanes[index].top_y),
            )
            for index in windows.record_indices
        )
        leading = windows.leading_record_index
        add_row_window_marks_group(
            canvas,
            build_group_content(
                RowWindowMarksGroup,
                windows,
                places,
                plan.frame.row_window_style,
                f"linear_row_window_marks_{plan.frame.rows_by_record[leading] + 1}",
            ),
            canvas_config,
            record_key=str(plan.frame.record_keys[leading]),
        )


def _draw_linear_record(canvas: Drawing, drawing: _LinearRecordDrawing) -> Drawing:
    """Draw one record: its slots by depth order, its region label fronts, its axis when no slot draws
    features, and its definition."""

    vertical_plan = drawing.vertical_plan
    for slot in slots_in_paint_order(drawing.plan.tracks.linear_track_layout.slots):
        slot_drawing = _LinearSlotDrawing(
            slot,
            vertical_plan.slot_by_id(slot.id),
            vertical_plan.slot_row_shifts(slot.id),
            _linear_slot_group_id(drawing, slot),
        )
        drawer = _LINEAR_SLOT_DRAWERS.get(str(slot.renderer))
        if drawer is None:
            continue
        canvas = drawer(canvas, drawing, slot_drawing)
    _draw_linear_region_label_fronts(canvas, drawing)
    if not slots_draw_features(drawing.plan.tracks.normalized_linear_track_slots):
        _draw_linear_axis_only(canvas, drawing)
    _draw_linear_record_definition(canvas, drawing)
    return canvas


@dataclass(frozen=True, slots=True)
class _LinearPrimaryBounds:
    """The drawn records' groups the composition moves, the scale's painted bounds, and the records' bounds.

    ``length_bar_bounds`` is where ``add_length_bar_on_linear_canvas`` drew
    the scale: one value for the primary bounds and the composition metadata.
    """

    targets: tuple[Group, ...]
    length_bar_bounds: Aabb | None
    source_bounds: Aabb


def _linear_primary_bounds(
    canvas: Drawing,
    primary_target_start: int,
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    blast_config: BlastMatchConfigurator,
    record_inputs: LinearRecordInputs,
    plan: _LinearPlan,
    labels: _LinearAnnotationLabels,
) -> _LinearPrimaryBounds:
    """Collect the groups drawn since the canvas was created and the bounds of everything they paint."""

    primary_targets = tuple(
        element
        for element in getattr(canvas, "elements", [])[primary_target_start:]
        if isinstance(element, Group)
    )
    length_bar_group = plan.length_bar.group
    length_bar_bounds = (
        _content_bounds(length_bar_group).translated(
            float(canvas_config.horizontal_offset) + float(plan.length_bar.offset_x),
            float(plan.length_bar.offset_y),
        )
        if length_bar_group is not None and _content_bounds(length_bar_group).width > 0.0
        else None
    )
    source_primary_bounds = _collect_linear_primary_bounds(
        records=record_inputs.records,
        record_placements=plan.placements.record_placements,
        record_plans=plan.vertical.record_vertical_plans,
        record_collision_bands=plan.stacking.record_collision_bands,
        labels_by_record=labels.labels_by_record,
        record_scalar_axis_rooms=plan.depth.record_scalar_axis_rooms,
        annotation_track_layouts=plan.tracks.annotation_track_layouts,
        annotation_labels_by_record=labels.annotation_labels_by_record,
        normalized_comparisons=plan.comparisons.normalized,
        canvas_config=canvas_config,
        feature_config=configurators.feature,
        gc_config=configurators.gc,
        skew_config=configurators.skew,
        blast_config=blast_config,
        render_context=plan.tracks.render_context,
        cfg=canvas_config.profile.config,
        multi_record_enabled=plan.frame.multi_record_enabled,
        length_bar_bounds=length_bar_bounds,
        record_transforms=record_inputs.record_transforms,
    )
    return _LinearPrimaryBounds(
        targets=primary_targets,
        length_bar_bounds=length_bar_bounds,
        source_bounds=source_primary_bounds,
    )


def _compose_linear_canvas(
    canvas: Drawing,
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    composition: CompositionInputs,
    legend: _LinearLegend,
    primary: _LinearPrimaryBounds,
) -> CompositionPlan:
    """Draw the legend (stage ``legend_content``) and the title, plan where they and the records go (stage
    ``composition_plan``) and move them there."""

    legend_placement = LegendPlacement(str(canvas_config.legend_position))
    legend_measurement = measure_legend(
        configurators.legend,
        legend.table,
        placement=canvas_config.legend_position,
        wrap_width=primary.source_bounds.width,
    )
    placed_legend = place_legend(canvas, linear_legend_content, legend_placement, legend_measurement, legend.table)
    title_target = place_title(canvas, composition)
    return compose_drawing(
        canvas,
        composition,
        CompositionPrimary(primary.source_bounds, primary.targets, length_bar_bounds=primary.length_bar_bounds),
        legend_placement,
        placed_legend,
        title_target,
    )


def _attach_linear_facts(
    canvas: Drawing,
    canvas_config: LinearCanvasConfigurator,
    record_inputs: LinearRecordInputs,
    plan: _LinearPlan,
    composition_plan: CompositionPlan,
    primary: _LinearPrimaryBounds,
    label_fit: LinearLabelFit,
) -> None:
    """Attach what the API and the Web read from the drawn figure to the canvas: the track slot geometry,
    the legend rows, the placements and translations, the bounds, the comparison facts and the label fit."""

    primary_placement = composition_plan.placement_for("primary")
    if primary_placement is None:  # pragma: no cover - primary is required
        raise RuntimeError("linear composition has no primary placement")
    placements = plan.placements
    track_slot_geometry = _serialize_linear_track_slot_geometry(
        record_axes=RecordAxisInput.of_records(
            record_inputs.records, plan.frame.record_keys,
            bind_record_transforms(record_inputs.records, record_inputs.record_transforms),
        ),
        units=plan.frame.units,
        layout=plan.tracks.linear_track_layout,
        record_plans=plan.vertical.record_vertical_plans,
        record_offsets=[
            float(placements.record_placements[index].axis_y) + primary_placement.dy
            for index in range(len(record_inputs.records))
        ],
        record_placements=placements.record_placements,
        composition_offset_x=primary_placement.dx,
        horizontal_offset_x=float(canvas_config.horizontal_offset),
        record_collision_bands=plan.stacking.record_collision_bands,
        boundary_gap_resolutions=plan.stacking.boundary_gap_resolutions,
    )
    attach_feature_placement_targets(track_slot_geometry["records"], plan.tracks.feature_placement_slot)
    attach_track_slot_geometry(canvas, track_slot_geometry)
    attach_legend_row_facts(
        canvas,
        LegendRowFacts(
            records=record_inputs.records,
            sources=plan.legend.sources,
            legend_table=plan.legend.table,
        ),
    )
    setattr(canvas, "_gbdraw_alignment_placements", tuple(
        replace(placement, x=placement.x + float(canvas_config.horizontal_offset))
        for placement in placements.alignment_placements
    ))
    # The translation each record was drawn with, by record key: what a Session
    # stores for an anchor_feature figure (LinearDiagramMetadata.record_translations).
    setattr(canvas, "_gbdraw_record_translations", tuple(
        (str(key), x, y) for key, (x, y) in zip(plan.frame.record_keys, plan.horizontal.final_translations, strict=True)
    ))
    setattr(canvas, "_gbdraw_comparison_record_id_warnings", plan.comparisons.rows.record_id_warnings)
    setattr(canvas, "_gbdraw_comparison_tallies", plan.comparisons.rows.tallies)
    setattr(canvas, LINEAR_LABEL_FIT_ATTR, label_fit)


def assemble_linear_diagram(
    canvas_config: LinearCanvasConfigurator,
    configurators: DiagramConfigurators,
    blast_config: BlastMatchConfigurator,
    record_inputs: LinearRecordInputs,
    track_inputs: LinearTrackInputs,
    comparison_inputs: LinearComparisonInputs,
    alignment: LinearAlignmentInputs,
    composition: CompositionInputs,
    feature_label_rows: FeatureLabelRows | None,
) -> Drawing:
    """
    Assembles a linear diagram of genomic records with optional BLAST comparison data,
    and returns the SVG canvas (not saved).

    ``record_inputs.record_presentations`` hold each record's label and subtitle.
    ``feature_label_rows`` are feature label rows an earlier assembly of the
    same figure proposed for its region labels (``LinearLabelFit.new_rows``,
    LP-D82): this assembly places those labels there and plans room for them.

    The phases run in this order; the canvas configurator's definition
    column gap, vertical offset, plot width and heights are set by the
    phases that name them, and the later phases read them.
    """
    frame = _linear_row_frame(canvas_config, record_inputs)
    # One feature label font size for the diagram; region labels without their own font take it too.
    label_font_size = _resolve_linear_diagram_label_font_size(
        record_inputs.records, canvas_config=canvas_config, profile=canvas_config.profile
    )
    slot_plan = _linear_slot_plan(canvas_config, configurators, record_inputs, track_inputs, frame, label_font_size)
    comparisons = _linear_comparisons(canvas_config, blast_config, record_inputs, comparison_inputs, frame)
    tracks = _linear_track_plan(
        canvas_config, configurators, record_inputs, track_inputs, frame, slot_plan, label_font_size
    )
    definitions = _linear_definition_metrics(canvas_config, record_inputs, frame)
    depth = _linear_depth_data(canvas_config, configurators, record_inputs, track_inputs, tracks)
    horizontal = _linear_horizontal_frame(canvas_config, record_inputs, alignment, frame, tracks, definitions)
    feature_labels = _linear_feature_labels(
        canvas_config, configurators, record_inputs, frame, comparisons, tracks, horizontal, feature_label_rows
    )
    vertical = _linear_vertical_plans(
        canvas_config, configurators, record_inputs, frame, tracks, definitions, depth, feature_labels, comparisons
    )
    gc_frames = _linear_gc_frames(configurators, record_inputs, tracks)
    legend = _linear_legend_table(
        canvas_config, configurators, blast_config, record_inputs, track_inputs, comparison_inputs, comparisons,
        tracks, depth,
    )
    stacking = _linear_row_stacking(canvas_config, record_inputs, frame, comparisons, definitions, horizontal, vertical)
    placements = _linear_vertical_placement(canvas_config, frame, tracks, horizontal, vertical, stacking)
    plan = _LinearPlan(
        frame=frame,
        comparisons=comparisons,
        tracks=tracks,
        definitions=definitions,
        depth=depth,
        horizontal=horizontal,
        vertical=vertical,
        gc_frames=gc_frames,
        legend=legend,
        stacking=stacking,
        placements=placements,
        length_bar=_linear_length_bar(canvas_config, frame, placements),
    )

    canvas, primary_target_start = _linear_canvas(canvas_config, blast_config, record_inputs, plan)
    labels = _linear_annotation_labels(canvas_config, configurators, record_inputs, plan, feature_label_rows)
    # The label rows were placed for this order (``linear_label_rows``).
    for record_index in linear_paint_order(len(record_inputs.records)):
        canvas = _draw_linear_record(
            canvas, _linear_record_drawing(canvas_config, configurators, record_inputs, plan, labels, record_index)
        )
    _draw_linear_row_window_marks(canvas, canvas_config, plan)
    primary = _linear_primary_bounds(
        canvas, primary_target_start, canvas_config, configurators, blast_config, record_inputs, plan, labels
    )
    composition_plan = _compose_linear_canvas(canvas, canvas_config, configurators, composition, plan.legend, primary)
    _attach_linear_facts(canvas, canvas_config, record_inputs, plan, composition_plan, primary, labels.label_fit)
    return canvas


__all__ = ["assemble_linear_diagram"]
