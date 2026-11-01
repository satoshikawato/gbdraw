#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from gbdraw.features.placement import FeaturePlacementSlot
from gbdraw.render_stages import render_stage
from ...layout.linear_wrap import drawn_axis_length, refuse_feature_track_on_the_axis
from ...layout.record_labels import source_ruler_ticks

import copy
from dataclasses import dataclass, replace
import logging
import math
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
    normalize_depth_tracks,
    representative_depth_tracks,
    sync_depth_track_legend_entries,
)
from ...analysis.protein_colinearity import OrthogroupResult, OrthogroupGraphResult
from ...auto_sizes import linear_tick_interval
from ...canvas import LinearCanvasConfigurator
from ...config.models import (
    GbdrawConfig,
    LinearRenderProfile,
)
from ...exceptions import ValidationError
from ...configurators import (
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
    LegendMeasurement,
)
from ...configurators.gc import _slot_skew_config
from ...configurators.legend import measure_legend
from ...core.text import calculate_svg_bbox_dimensions, drawn_text_box
from ...render.group_content import (
    CANVAS_SCOPE,
    GroupContent,
    GroupKey,
    GroupKind,
    build_group_content,
    place_content,
)
from ...render.groups.linear import (
    LengthBarGroup,
    linear_legend_content,
    plot_title_content,
)
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
    LegendTable,
    PairwiseIdentityLegend,
    pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
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
    CompositionItem,
    CompositionRequest,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
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
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.record_placement import resolve_record_row_positions
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
from ...layout.scalar_axis import linear_scalar_axis_tick_font_size_px
from ...layout.scalar_axis import (
    LINEAR_SCALAR_AXIS_LABEL_GAP_PX,
    LINEAR_SCALAR_AXIS_STROKE_PX,
    LINEAR_SCALAR_AXIS_TICK_PX,
    _depth_axis_bounds,
    _format_depth_tick,
    format_percent_tick,
    scalar_axis_tick_values,
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
    with_feature_label_rows,
)
from ...tracks import (
    LinearTrackSlot,
    ScalarSpec,
    default_linear_track_slots,
    normalize_linear_track_slots_with_axis,
    parse_nonnegative_integer,
)
from ...tracks.parsing import slot_dinucleotide, slot_legend_label
from ...annotations import (
    AnnotationOptions,
    ResolvedAnnotationBundle,
    ResolvedAnnotationTrack,
    annotation_track_params_from_mapping,
    effective_annotation_style,
    feature_underlay_anchor_slot_id,
    feature_underlay_slot_id,
    layout_annotation_track,
    linear_annotation_labels_above,
    merge_feature_underlays,
    sync_annotation_legend_entries,
)
from ...annotations.planning import prepare_annotation_track_slots
from ...render.drawers.linear.annotations import (
    linear_annotation_label_room,
    linear_annotation_track_content,
    linear_region_label_front_content,
)
from .annotation_labels import (
    AnnotationLabelsBySlot,
    LabelObstacle,
    PlacedAnnotationLabels,
    linear_annotation_label_inks,
    linear_annotation_row_frames,
    linear_feature_glyph_boxes,
    linear_feature_label_obstacles,
    place_linear_annotation_labels,
)
from ...render.groups.linear import build_linear_feature_dom_index
from ...render.composition import apply_composition_plan, painted_title_content

from .builders import (
    add_explicit_comparisons_on_linear_canvas,
    add_depth_group,
    add_gc_content_group,
    add_gc_skew_group,
    add_length_bar_on_linear_canvas,
    add_record_definition_group,
    add_record_group,
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
from ...features.ids import make_linear_dom_id
from ...svg.ids import (
    definition_group_svg_id,
    instance_svg_id,
    record_group_svg_id,
    track_slot_svg_id,
)
from .track_slots import (
    LinearRecordVerticalPlan,
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


def _linear_feature_placement_slot(side: str, separate_strands: bool) -> FeaturePlacementSlot:
    # FeaturePlacementSlot.__post_init__ validates the side for linear mode.
    return FeaturePlacementSlot(
        "linear",
        cast(Literal["overlay", "above", "below"], side),
        separate_strands,
    )


def _linear_track_slot_dom_id(
    *,
    slot_id: str,
    renderer: str,
    slot_index: int,
    record_index: int,
    record_count: int,
) -> str:
    """Return the deterministic DOM ID for one rendered Linear user slot."""

    group_id = track_slot_svg_id(
        slot_id,
        renderer=renderer,
        slot_index=slot_index,
    )
    if int(record_count) <= 1:
        return group_id
    return instance_svg_id(
        group_id,
        ("linear-record", int(record_index), int(record_count)),
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
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
    )
    return updated_slots, bundle, layouts, auto_slot_ids


def _layout_linear_annotation_tracks(
    records: list[SeqRecord],
    slots: list[LinearTrackSlot],
    bundle: ResolvedAnnotationBundle,
    *,
    auto_slot_ids: frozenset[str],
    label_font_size: float,
    sequence_widths: Sequence[float],
    record_transforms: Sequence[RecordDisplayTransform] | None,
) -> tuple[list[LinearTrackSlot], dict[str, ResolvedAnnotationTrack]]:
    """Resolve annotation lanes for the records' drawn axis widths (their placements'
    ``sequence_width``).

    ``label_font_size`` is the feature label font size of the diagram, the size of a region label whose style
    sets none."""

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
        layout = layout_annotation_track(
            slot.id,
            params.set_id,
            bundle.annotations,
            record_lengths=record_lengths,
            params=params,
            available_extent_px=available,
            bp_per_px=bp_per_px,
            label_font_size=label_font_size,
        )
        layouts[slot.id] = layout
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
    if set_id is None:
        return slots, bundle, layouts, auto_slot_ids, 0
    anchor_id = feature_underlay_anchor_slot_id(slots)
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
    params = getattr(slot, "params", {}) or {}
    return parse_nonnegative_integer(
        params.get("track_index", 0),
        field_name=f"depth slot '{getattr(slot, 'id', '')}' track_index",
    )


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
    records: list[SeqRecord],
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
    for record_index, record in enumerate(records):
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
        record_id = str(getattr(record, "id", "") or "")
        records_payload.append(
            {
                "recordIndex": int(record_index),
                "recordId": record_id,
                "recordLabel": record_id,
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
        "INFO: BLAST output visualization settings: e-value threshold: {}; bitscore threshold: {}; "
        "identity threshold: {}; alignment length threshold: {}".format(
            thresholds.evalue, thresholds.bitscore, thresholds.identity, thresholds.alignment_length
        )
    )
    logger.info("INFO: Loading comparison file(s)...")
    tables: list[DataFrame] = []
    for comparison_file in blast_files:
        logger.info("INFO: Loading {}".format(comparison_file))
        tables.append(read_comparison_resource(comparison_file))
    logger.info("INFO:             ... finished loading comparison file(s)")
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
                comparison_row_tally(item.matches, thresholds), protein_matches=protein,
            )
            for items, protein in ((linear_comparisons, False), (proteins.pairs, True))
            for item in items
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


def _axis_text_overhang(axis_config: object, track_height: float) -> float:
    if not bool(getattr(axis_config, "show_axis", False)):
        return 0.0
    axis_stroke = 0.5 * LINEAR_SCALAR_AXIS_STROKE_PX
    if not bool(getattr(axis_config, "show_ticks", False)):
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
            footprints[slot.id] = _linear_annotation_slot_footprint(
                slot,
                record_index=record_index,
                annotation_layout=annotation_track_layouts.get(slot.id),
                font_family=font_family,
                wrap_row=wrap_row,
            )
            continue

        footprints[slot.id] = LinearSlotFootprint(body_band, body_band)
    feature_center_y = 0.5 * (
        feature_geometry.occupied_band.top_y
        + feature_geometry.occupied_band.bottom_y
    )
    return footprints, feature_center_y


def _linear_annotation_slot_footprint(
    slot: LinearResolvedTrack,
    *,
    record_index: int,
    annotation_layout: ResolvedAnnotationTrack | None,
    font_family: str,
    span_height: float | None = None,
    wrap_row: int | None = None,
) -> LinearSlotFootprint:
    """An annotation slot's own band with its labels' room on their side (on one row of a wrapped record:
    that row's labels, and paint only where the row holds a mark).

    ``span_height`` is the stretched height of an overlay the record's plan
    spans over slots (``spans_slots``): its labels then keep clear of every
    slot it covers, measured on that height; otherwise on the slot's own.
    """

    body_band = _slot_body_band(slot)
    params = annotation_track_params_from_mapping(slot.params)
    label_overhang = (
        linear_annotation_label_room(
            annotation_layout,
            record_index=record_index,
            side=slot.side,
            covers_slots=span_height is not None,
            height_px=float(body_band.height if span_height is None else span_height),
            font_family=font_family,
            params=params,
            wrap_row=wrap_row,
        )
        if annotation_layout is not None
        else 0.0
    )
    if linear_annotation_labels_above(slot.side, params):
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
) -> _LinearRecordDefinitionGeometry:
    """The definition text drawn beside a record's features, centered on them: a record alone in its
    row has its own text there, and the leading record of a shared row the row's heading. Its place does
    not depend on the feature labels, so the labels can keep clear of it (OV-1395); the record-local
    text of a shared row is a header above the record's canvas band, added by the vertical plan."""

    local_band = None if multi_record_enabled else _centered_vertical_band(feature_center_y, local_height)
    row_band = (
        _centered_vertical_band(feature_center_y, row_height)
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
    definition_heights: list[float],
    definition_widths: list[float],
    row_definition_heights: list[float],
    row_definition_widths: list[float],
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
        local_definition_height = (
            float(definition_heights[record_index])
            if record_index < len(definition_heights)
            else 0.0
        )
        beside = _definition_geometry_beside_features(
            feature_center_y,
            local_height=local_definition_height,
            local_width=(
                float(definition_widths[record_index])
                if record_index < len(definition_widths)
                else 0.0
            ),
            row_height=(
                float(row_definition_heights[record_index])
                if record_index < len(row_definition_heights)
                else 0.0
            ),
            row_width=(
                float(row_definition_widths[record_index])
                if record_index < len(row_definition_widths)
                else 0.0
            ),
            multi_record_enabled=multi_record_enabled,
            split_row_definitions=split_row_definitions,
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
    if record_local_ruler:
        start_coord, end_coord = 0, len(record.seq)
    else:
        annotations = getattr(record, "annotations", None) or {}
        try:
            start_coord = int(annotations.get("gbdraw_coord_base", 1))
        except (TypeError, ValueError):
            start_coord = 1
        try:
            step = int(annotations.get("gbdraw_coord_step", 1))
        except (TypeError, ValueError):
            step = 1
        step = 1 if step >= 0 else -1
        end_coord = start_coord + (step * max(0, len(record.seq) - 1))

    span = abs(end_coord - start_coord)
    interval = cfg.objects.scale.interval
    tick_interval = (
        int(interval)
        if interval is not None and int(interval) > 0
        else linear_tick_interval(max(1, int(canvas_config.longest_genome)))
    )
    if span <= 0 or tick_interval <= 0:
        return None
    minimum = min(start_coord, end_coord)
    maximum = max(start_coord, end_coord)
    tick = (minimum // tick_interval) * tick_interval
    if tick < minimum:
        tick += tick_interval
    ticks: list[int] = []
    while tick < maximum:
        if tick > minimum:
            ticks.append(int(tick))
        tick += tick_interval
    projected_ticks = None
    if record_transform is not None and record_transform.projects:
        projected_ticks = {tick.base: tick for tick in source_ruler_ticks(record_transform, tick_interval)}
        ticks = list(projected_ticks)
    if not ticks:
        return None

    bounds: list[Aabb] = []
    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param(
        canvas_config.length_param
    )
    axis_length = drawn_axis_length(record, record_transform)
    row_offsets = plan.row_offsets
    for coordinate in ticks:
        offset = (projected_ticks[coordinate].offset / axis_length if projected_ticks is not None
                  else abs(coordinate - start_coord) / float(span))
        x = placement.sequence_width * offset
        row_y = row_offsets[projected_ticks[coordinate].row] if projected_ticks is not None else 0.0
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
        labels: list[str] = []
        axis_config: object
        if slot.renderer == "depth":
            depth_track = depth_by_index.get(_depth_slot_track_index(slot))
            if depth_track is None or not bool(depth_track.config.show_axis):
                continue
            axis_config = depth_track.config
            min_depth = depth_track.config.min_depth
            axis_min, axis_max = _depth_axis_bounds(depth_track.df, min_depth, depth_track.config.max_depth)
            labels = [
                _format_depth_tick(value)
                for value in scalar_axis_tick_values(
                    axis_min,
                    axis_max,
                    show_ticks=bool(depth_track.config.show_ticks),
                    large_tick_interval=depth_track.config.large_tick_interval,
                )
                if not (min_depth is None and math.isclose(value, axis_min, rel_tol=1e-9, abs_tol=1e-9))
            ]
        elif slot.renderer == "dinucleotide_content":
            gc_slot_config = gc_config_of(slot.params, depth_by_index)
            if str(gc_slot_config.mode) != "percent" or not bool(gc_slot_config.show_axis):
                continue
            axis_config = gc_slot_config
            labels = [
                format_percent_tick(value)
                for value in scalar_axis_tick_values(
                    float(0.0 if gc_slot_config.min_percent is None else gc_slot_config.min_percent),
                    float(100.0 if gc_slot_config.max_percent is None else gc_slot_config.max_percent),
                    show_ticks=bool(gc_slot_config.show_ticks),
                    large_tick_interval=gc_slot_config.large_tick_interval,
                )
            ]
        else:
            continue
        axis_px = max(axis_px, 0.5 * LINEAR_SCALAR_AXIS_STROKE_PX)
        font_size = linear_scalar_axis_tick_font_size_px(axis_config, float(slot.height))
        for label in labels:
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
                label_x=next(x for other, x in label_row if other == index),
                row_labels=[(labels_of(other), x) for other, x in label_row if other != index],
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
    *,
    window: int,
    step: int,
    dinucleotide: str,
    enabled: bool,
) -> list[DataFrame | None]:
    """Build GC/skew data once per record and share it across linear groups."""
    if not enabled:
        return [None for _ in records]
    return [render_stage("gc_skew_data", skew_df, record, window, step, dinucleotide) for record in records]






def _linear_depth_group_id(
    base_id: str,
    *,
    record_index: int,
    record_count: int,
) -> str:
    """Return a stable linear depth group id."""

    return make_linear_dom_id(
        base_id or "depth",
        record_index=record_index,
        record_count=record_count,
    )


def _final_record_translations(
    *,
    record_keys: Sequence[str],
    placements: Mapping[int, LinearRecordPlacement],
    layout: LinearMultiRecordOptions | None,
    similarity_alignment: SimilarityAlignmentPlan | None,
    anchor_centers: Sequence[float | None] | None,
) -> tuple[tuple[float, float], ...]:
    """Resolve base translations and absolute anchor alignment exactly once.

    ``anchor_centers`` are the display positions the records align at: the
    Similarity alignment anchor centres with ``similarity_alignment``, else the
    ``anchor_feature`` starts, whose reference is the first record that has one.
    An aligned record moves so that its position meets the reference's; any
    other record keeps its base translation.
    """

    base_by_key = {
        str(item.record_key): (float(item.x), float(item.y))
        for item in (layout.record_translations if layout is not None else ())
    }
    base = tuple(base_by_key.get(str(key), (0.0, 0.0)) for key in record_keys)
    centers = tuple(anchor_centers or ())
    if similarity_alignment is None:
        reference_index = next(
            (index for index, center in enumerate(centers) if center is not None), None
        )
        if reference_index is None:
            return base
        aligned = {
            index for index, center in enumerate(centers)
            if center is not None and index != reference_index
        }
    else:
        similarity_alignment.validate_record_coverage(record_keys)
        if len(centers) != len(record_keys):
            raise ValidationError(
                "Similarity alignment requires one projected anchor center per record."
            )
        index_by_key = {str(key): index for index, key in enumerate(record_keys)}
        reference_index = index_by_key[similarity_alignment.reference.record_key]
        if centers[reference_index] is None:
            raise ValidationError("Similarity alignment reference center is missing.")
        decisions = {decision.record_key: decision for decision in similarity_alignment.records}
        aligned = {
            index for index, record_key in enumerate(record_keys)
            if decisions[str(record_key)].status is AlignmentDecisionStatus.ALIGNED
        }
        for index in aligned:
            if centers[index] is None:
                raise ValidationError(
                    f"Similarity alignment center is missing for record {record_keys[index]!r}."
                )
    reference_world_x = (
        placements[reference_index].x_for_position(cast(float, centers[reference_index]))
        + base[reference_index][0]
    )
    return tuple(
        (
            reference_world_x - placements[index].x_for_position(cast(float, centers[index]))
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
) -> LegendTable:
    """Build the Linear legend table (render stage ``legend_table``).

    ``depth`` is given only when one depth track has a singleton row.
    """

    features_present, used_color_rules, default_used_features = feature_legend_presence(
        records, sources.features
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
    return sync_annotation_legend_entries(legend_table, annotations, slots)


def assemble_linear_diagram(
    records: list[SeqRecord],
    blast_files,
    canvas_config: LinearCanvasConfigurator,
    blast_config,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    legend_config: LegendDrawingConfigurator,
    skew_config,
    depth_config: DepthConfigurator | None = None,
    depth_tables: list[DataFrame | None] | None = None,
    record_depth_tracks: list[list[DepthTrackSpec]] | None = None,
    linear_track_slots: list[LinearTrackSlot] | None = None,
    linear_track_axis_index: int | None = None,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None,
    plot_title: str | None = None,
    plot_title_position: str = "bottom",
    plot_title_font_size: float = 32.0,
    protein_comparisons: ProteinComparisonFrames | None = None,
    pairwise_identity_legend: PairwiseIdentityLegend | None = None,
    linear_comparisons: list[LinearComparison] | None = None,
    linear_layout: LinearMultiRecordOptions | None = None,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None,
    similarity_alignment: SimilarityAlignmentPlan | None = None,
    alignment_anchor_centers: Sequence[float | None] | None = None,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    record_presentations: Sequence[RecordPresentation] | None = None,
    feature_label_rows: FeatureLabelRows | None = None,
) -> Drawing:
    """
    Assembles a linear diagram of genomic records with optional BLAST comparison data,
    and returns the SVG canvas (not saved).

    ``record_presentations`` hold each record's label and subtitle.
    ``feature_label_rows`` are feature label rows an earlier assembly of the
    same figure proposed for its region labels (``LinearLabelFit.new_rows``,
    LP-D82): this assembly places those labels there and plans room for them.
    """
    profile = canvas_config.profile
    cfg = profile.config
    presentations: tuple[RecordPresentation | None, ...] = (
        tuple(record_presentations) if record_presentations is not None else (None,) * len(records)
    )
    record_keys = stable_record_keys(records)
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
    split_row_definitions = bool(multi_record_enabled)
    # Every record's horizontal place, solved once before its tracks and labels, in both row
    # modes: its x, width and scale, with zero vertical extents that the row stacking fills in.
    axis_lengths = _drawn_axis_lengths(records, record_transforms)
    horizontal_layout = solve_linear_horizontal_layout(
        [
            LinearRecordMeasurement(
                record_index=index,
                record_key=RecordKey(str(record_keys[index])),
                sequence_length=axis_lengths[index],
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
    sequence_widths = [horizontal_placements[index].sequence_width for index in range(len(records))]
    if record_depth_tracks is None and depth_tables:
        record_depth_tracks = normalize_depth_tracks(
            records,
            depth_track_tables=[[table] for table in depth_tables],
        )
    user_slot_mode = linear_track_slots is not None
    # One feature label font size for the diagram; region labels without their own font take it too.
    label_font_size = _resolve_linear_diagram_label_font_size(
        records,
        canvas_config=canvas_config,
        profile=profile,
    )
    (
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
    ) = _prepare_linear_annotation_tracks(
        records,
        annotations,
        linear_track_slots,
        profile=profile,
        record_depth_tracks=record_depth_tracks,
        label_font_size=label_font_size,
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
    )
    if linear_track_slots is None:
        linear_track_slots = default_linear_track_slots(
            show_features=True,
            show_depth=bool(profile.show_depth and record_depth_tracks),
            depth_track_count=max(1, depth_track_count(record_depth_tracks)),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            dinucleotide=str(gc_config.dinucleotide),
            track_layout=profile.track_layout,
        )
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        record_depth_tracks,
    )
    normalized_plot_title = str(plot_title or "").strip()
    normalized_plot_title_position = str(plot_title_position or "bottom").strip().lower()
    if normalized_plot_title_position not in {"center", "top", "bottom"}:
        raise ValueError("plot_title_position must be one of: center, top, bottom")

    # A typed title that paints nothing ("<i></i>") gets no title (OV-1845).
    title_content: GroupContent | None = None
    if normalized_plot_title:
        title_content = painted_title_content(render_stage(
            "plot_title_content",
            plot_title_content,
            normalized_plot_title,
            font_family=cfg.objects.text.font_family,
            font_size=float(plot_title_font_size),
        ))

    comparison_thresholds = ComparisonThresholds(
        evalue=blast_config.evalue,
        bitscore=blast_config.bitscore,
        identity=blast_config.identity,
        alignment_length=blast_config.alignment_length,
    )
    comparison_rows = render_stage(
        "comparison_rows",
        linear_comparison_rows,
        _read_comparison_files(blast_files, comparison_thresholds),
        protein_comparisons,
        tuple(linear_comparisons or ()),
        comparison_thresholds,
        tuple(records),
        tuple(rows_by_record),
        bool(multi_record_enabled),
        tuple(Path(path).name for path in blast_files or ()),
    )
    normalized_comparisons = list(comparison_rows.comparisons)
    comparison_record_id_warnings = comparison_rows.record_id_warnings
    comparisons = [item.matches for item in normalized_comparisons]
    has_blast = bool(normalized_comparisons)

    label_scope = profile.label_scope
    orthogroup_label_eligibility = None
    if label_scope == "orthogroup_top":
        orthogroup_label_eligibility = build_orthogroup_label_eligibility(
            orthogroups=orthogroups,
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
            "WARNING: --show_labels orthogroup_top requires orthogroup metadata; no orthogroup-specific label suppression was applied."
        )

    resolved_feature_slot = next(
        (slot for slot in normalized_linear_track_slots if slot.renderer == "features"),
        None,
    )
    record_feature_layers = _precalculate_feature_layers(
        records,
        feature_config,
        profile,
        record_transforms=record_transforms,
        placement_slot=(
            _linear_feature_placement_slot(resolved_feature_slot.side, profile.strandedness)
            if resolved_feature_slot is not None
            else None
        ),
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
        auto_annotation_slot_ids,
        private_underlay_prefix_count,
    ) = _add_linear_feature_underlays(
        records,
        record_feature_layers,
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
        label_font_size=label_font_size,
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
    )
    effective_linear_track_axis_index = (
        linear_track_axis_index + private_underlay_prefix_count
        if linear_track_axis_index is not None
        else None
    )
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        effective_linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        record_depth_tracks,
    )
    linear_track_layout = resolve_linear_track_layout(
        normalized_linear_track_slots,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    depth_available = bool(depth_config is not None and record_depth_tracks)
    render_context = LinearRecordRenderContext(
        profile=profile,
        track_layout=linear_track_layout,
        depth_available=depth_available,
    )
    axis_ruler_enabled = render_context.axis_ruler_enabled
    if any(transform.row_count > 1 for transform in record_transforms or ()):
        refuse_feature_track_on_the_axis(render_context.feature_track_layout)
    resolved_depth_heights_by_index = _resolved_linear_depth_heights(
        linear_track_layout
    )
    primary_depth_track_index = min(resolved_depth_heights_by_index, default=0)
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
    local_definition_line_kinds: list[frozenset[str]] | None = None
    row_definition_line_kinds: list[frozenset[str]] = [
        frozenset() for _record in records
    ]
    if split_row_definitions:
        (
            local_definition_line_kinds,
            row_definition_line_kinds,
        ) = _split_definition_line_kinds(
            presentations,
            rows_by_record=rows_by_record,
            row_leading_indices=row_leading_indices,
        )
    max_def_width, definition_widths, definition_heights = _precalculate_definition_metrics(
        records,
        canvas_config,
        cfg=cfg,
        line_kinds_by_record=local_definition_line_kinds,
        record_presentations=presentations,
    )
    row_definition_width = 0.0
    row_definition_widths = [0.0 for _record in records]
    row_definition_heights = [0.0 for _record in records]
    if split_row_definitions:
        (
            row_definition_width,
            row_definition_widths,
            row_definition_heights,
        ) = _precalculate_definition_metrics(
            records,
            canvas_config,
            cfg=cfg,
            line_kinds_by_record=list(row_definition_line_kinds),
            record_presentations=presentations,
        )
    record_depth_data: list[list[DepthTrackData]] = (
        build_depth_track_dataframes(
            records,
            record_depth_tracks,
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
            linear_track_layout.slots,
            gc_config_of=drawn_gc_config,
            depth_by_index=depth_by_index,
            cfg=cfg,
        )
        for depth_by_index in record_depth_by_index
    ]
    canvas_config.configure_definition_column_gap(
        max((room.tick_labels_px for room in record_scalar_axis_rooms), default=0.0)
    )

    def resolve_horizontal_frame(
        placements: Mapping[int, LinearRecordPlacement],
    ) -> tuple[tuple[tuple[float, float], ...], float, list[tuple[Aabb, ...]]]:
        """Each record's final translation, the definition column's left (the leftmost translated
        record), and the definition boxes its feature labels keep clear of (OV-1395) in their
        coordinates. They read only the records' horizontal places, which the labels do not move,
        so the drawing and the fit count read the same values."""

        translations = _final_record_translations(
            record_keys=record_keys,
            placements=placements,
            layout=linear_layout,
            similarity_alignment=similarity_alignment,
            anchor_centers=alignment_anchor_centers,
        )
        record_xs = [placements[index].x + translations[index][0] for index in range(len(records))]
        column_left = min(record_xs, default=0.0)
        keep_clear = [
            _definition_label_boxes(
                _record_definition_bands(
                    _definition_geometry_beside_features(
                        0.5 * (geometry.occupied_band.top_y + geometry.occupied_band.bottom_y),
                        local_height=definition_heights[index],
                        local_width=definition_widths[index],
                        row_height=row_definition_heights[index],
                        row_width=row_definition_widths[index],
                        multi_record_enabled=multi_record_enabled,
                        split_row_definitions=split_row_definitions,
                    ),
                    sequence_width=placements[index].sequence_width,
                    record_x=record_xs[index],
                    multi_record_enabled=multi_record_enabled,
                    keep_definition_left_aligned=bool(canvas_config.keep_definition_left_aligned),
                    definition_column_width=max_def_width,
                    row_definition_width=row_definition_width,
                    definition_gap=float(canvas_config.definition_column_gap),
                    text_anchor=cfg.objects.definition.linear.text_anchor,
                    definition_column_left=column_left,
                ),
                record_x=record_xs[index],
                features_origin_y=0.0,
            )
            for index, geometry in enumerate(record_feature_lane_geometries)
        ]
        return translations, column_left, keep_clear

    final_translations, definition_column_left, label_keep_clear = resolve_horizontal_frame(
        horizontal_placements
    )
    # Where each record's feature labels are stacked and counted: its row and its start in it
    # (OV-1188); a record alone in its row needs no frame.
    label_frames: list[LinearLabelFrame] | None = (
        [
            LinearLabelFrame(rows_by_record[index], horizontal_placements[index].x)
            for index in range(len(records))
        ]
        if multi_record_enabled
        else None
    )
    # The labels are placed once, at the records' drawn widths; the row stacking reads their heights.
    _unused_label_height, all_labels, _record_label_heights_above = _precalculate_label_dimensions(
        records,
        feature_config,
        canvas_config,
        render_context=render_context,
        precomputed_feature_dicts=record_feature_dicts,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
        feature_lane_geometries=record_feature_lane_geometries,
        sequence_widths=sequence_widths,
        record_transforms=record_transforms,
        label_frames=label_frames,
        label_keep_clear=label_keep_clear,
    )
    if feature_label_rows is not None:
        all_labels = with_feature_label_rows(all_labels, feature_label_rows.rows)

    record_vertical_plans, record_definition_geometries = _build_linear_record_vertical_plans(
        records=records,
        layout=linear_track_layout,
        feature_dicts=record_feature_dicts,
        feature_geometries=record_feature_lane_geometries,
        labels_by_record=all_labels,
        canvas_config=canvas_config,
        gc_config_of=drawn_gc_config,
        skew_config=skew_config,
        record_depth_by_index=record_depth_by_index,
        representative_depth_configs=representative_depth_configs,
        annotation_track_layouts=annotation_track_layouts,
        definition_heights=definition_heights,
        definition_widths=definition_widths,
        row_definition_heights=row_definition_heights,
        row_definition_widths=row_definition_widths,
        multi_record_enabled=multi_record_enabled,
        split_row_definitions=split_row_definitions,
        render_context=render_context,
        record_transforms=record_transforms,
    )
    # Every later reader counts a wrapped record's labels in one frame, where rows do not overlap.
    feature_slot_id = next((slot.id for slot in linear_track_layout.slots if slot.renderer == "features"), None)
    all_labels = [
        _labels_in_first_row_frame(labels, plan, feature_slot_id)
        for labels, plan in zip(all_labels, record_vertical_plans, strict=True)
    ]

    canvas_config.vertical_offset = _linear_record_vertical_offset(
        record_vertical_plans,
        rows_by_record,
        canvas_config,
    )

    needed_nts = {
        slot_dinucleotide(slot.params, str(gc_config.dinucleotide))
        for slot in linear_track_layout.slots
        if slot.renderer in {"dinucleotide_content", "dinucleotide_skew"}
    }
    record_gc_dfs_by_nt = {
        nt: _precalculate_gc_dataframes(
            records,
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
    # Prepare legend group
    legend_sources = LegendSources.of(feature_config, gc_config, skew_config)
    legend_table = LegendTable()
    if canvas_config.legend_position != "none":
        singleton_depth = render_context.depth_enabled and depth_track_count(record_depth_tracks) == 1
        comparison_identity_legend = pairwise_identity_legend_from_comparisons(blast_config, comparisons)
        legend_table = render_stage(
            "legend_table",
            build_linear_legend_table,
            records,
            sources=legend_sources,
            comparison=ComparisonLegendSource.of(blast_config) if has_blast and blast_config else None,
            pairwise_identity_legend=(
                comparison_identity_legend
                or pairwise_identity_legend
                or PairwiseIdentityLegend()
            ),
            depth=(
                DepthLegendSource.of_config(depth_config)
                if singleton_depth and depth_config is not None
                else None
            ),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            slots=legend_slot_sources(normalized_linear_track_slots),
            depth_tracks=(
                depth_legend_sources(representative_depth_data) if render_context.depth_enabled else None
            ),
            annotations=resolved_annotations,
        )
    definition_reserve_width = (
        row_definition_width
        if split_row_definitions
        else (0.0 if multi_record_enabled else max_def_width)
    )
    canvas_config.configure_plot_width(definition_reserve_width)
    ordinary_row_gap = _ordinary_row_gap(canvas_config)
    definition_clear_gap = max(1.0, 0.5 * float(canvas_config.vertical_padding))
    comparison_endpoint_gap_px = (
        _COMPARISON_ENDPOINT_GAP_PX if has_blast else 0.0
    )
    comparison_boundary_records = _comparison_records_by_boundary(
        normalized_comparisons,
        rows_by_record,
    )


    definition_column_width = max_def_width

    def drawn_collision_bands() -> list[tuple[CollisionBand, ...]]:
        """Each record's collision bands where the drawing puts them (OV-1400): at the record's
        final translation, with the definition column's left and the text anchor the drawing uses.
        The row stacking reads these, so the gap between rows is solved from the drawn boxes."""

        return [
            _record_collision_bands(
                plan=record_vertical_plans[index],
                definition_geometry=record_definition_geometries[index],
                sequence_width=horizontal_placements[index].sequence_width,
                record_x=horizontal_placements[index].x + final_translations[index][0],
                multi_record_enabled=multi_record_enabled,
                text_anchor=cfg.objects.definition.linear.text_anchor,
                keep_definition_left_aligned=bool(canvas_config.keep_definition_left_aligned),
                definition_column_width=definition_column_width,
                row_definition_width=row_definition_width,
                definition_gap=float(canvas_config.definition_column_gap),
                definition_column_left=definition_column_left,
            )
            for index in range(len(records))
        ]

    record_offsets: list[float] = [0.0 for _record in records]
    record_offsets_x: list[float] = [0.0 for _record in records]
    record_collision_bands: list[tuple[CollisionBand, ...]]
    boundary_gap_resolutions: list[AxisGapResolution] = []
    multi_record_plan = None

    def resolve_first_axis_y() -> float:
        return max(
            canvas_config.vertical_offset,
            canvas_config.original_vertical_offset,
        )

    current_y = resolve_first_axis_y()

    if multi_record_enabled:
        record_collision_bands = drawn_collision_bands()
        multi_record_plan = solve_linear_layout(
            [
                LinearRecordMeasurement(
                    record_index=index,
                    record_key=RecordKey(str(record_keys[index])),
                    sequence_length=axis_lengths[index],
                    top_extent=record_vertical_plans[index].canvas_top_extent,
                    bottom_extent=record_vertical_plans[index].canvas_bottom_extent,
                    comparison_top_extent=record_vertical_plans[index].comparison_top_extent,
                    comparison_bottom_extent=record_vertical_plans[index].comparison_bottom_extent,
                    # The drawn bands, local to the record's x, which the solve adds back.
                    collision_bands=tuple(
                        band.translate(x=-horizontal_placements[index].x)
                        for band in record_collision_bands[index]
                    ),
                )
                for index, _record in enumerate(records)
            ],
            horizontal_layout,
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
        current_y = max(record_offsets)
    else:
        record_collision_bands = drawn_collision_bands()

        ordered_rows = tuple(ordered_record_indices)
        for row_position, record_index in enumerate(ordered_rows):
            record_offsets[record_index] = current_y
            if row_position >= len(ordered_rows) - 1:
                continue
            next_record_index = ordered_rows[row_position + 1]
            boundary = int(rows_by_record[record_index])
            active_records = frozenset(
                comparison_boundary_records.get(boundary, ())
            )

            def active_bands(index: int) -> tuple[CollisionBand, ...]:
                return tuple(
                    band
                    for band in record_collision_bands[index]
                    if band.kind != "comparison" or index in active_records
                )

            resolution = resolve_axis_gap(
                active_bands(record_index),
                active_bands(next_record_index),
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
        current_y = max(record_offsets, default=current_y)

    length_bar_group: GroupContent | None = None
    length_bar_offset_x = 0.0
    draw_length_bar = (
        cfg.objects.scale.show
        and not multi_record_enabled
        and not canvas_config.normalize_length
        and not axis_ruler_enabled
    )
    record_placements: dict[int, LinearRecordPlacement] = {}
    for record_index in range(len(records)):
        if multi_record_plan is not None:
            record_placements[record_index] = multi_record_plan.placement_for_index(
                record_index
            )
            continue
        record_plan = record_vertical_plans[record_index]
        record_placements[record_index] = place_record_vertically(
            horizontal_placements[record_index],
            axis_y=record_offsets[record_index],
            top_extent=record_plan.canvas_top_extent,
            bottom_extent=record_plan.canvas_bottom_extent,
            comparison_top_extent=record_plan.comparison_top_extent,
            comparison_bottom_extent=record_plan.comparison_bottom_extent,
            comparison_endpoint_gap=comparison_endpoint_gap_px,
        )
    # Keep raw layout placements available to the existing alignment helper.
    # They exclude canonical translations and composition/viewBox fitting.
    alignment_placements = tuple(record_placements.values())
    for record_index, (translation_x, translation_y) in enumerate(
        final_translations
    ):
        placement = record_placements[record_index]
        record_placements[record_index] = replace(
            placement,
            x=placement.x + translation_x,
            axis_y=placement.axis_y + translation_y,
            comparison_top_y=placement.comparison_top_y + translation_y,
            comparison_bottom_y=placement.comparison_bottom_y + translation_y,
        )
        record_offsets_x[record_index] = record_placements[record_index].x
        record_offsets[record_index] = record_placements[record_index].axis_y
    painted_content_bottom = max(
        (
            float(record_offsets[index]) + plan.canvas_band.bottom_y
            for index, plan in enumerate(record_vertical_plans)
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
    # Row spacing reserves labels, but ribbon endpoints attach to track paint.
    # Resolve this after both layout paths have fixed every slot and record.
    # A feature's attachment is kept as its offset from the record's anchors, so moving rows leaves it as it is.
    feature_attachments: dict[int, dict[str, LinearFeatureAttachment]] = {}
    for record_index, placement in record_placements.items():
        plan = record_vertical_plans[record_index]
        attachment_band = plan.axis_band
        non_feature_band: VerticalBand | None = None
        for plan_slot in plan.slots:
            band = plan_slot.paint_band
            if band is None:
                continue
            if plan_slot.renderer == "features":
                band = record_feature_lane_geometries[
                    record_index
                ].occupied_band.translate(plan_slot.origin_y)
            elif not plan_slot.spans_slots:
                # A span's band comes from the slots it covers, counted here or as lanes.
                non_feature_band = band if non_feature_band is None else non_feature_band.union(band)
            attachment_band = attachment_band.union(band)
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
        features = record_feature_dicts[record_index]
        if not any(feature.placement is not None and feature.placement.requested_target is not None
                   for feature in features.values()):
            continue
        if not any(feature.placement is not None and feature.placement.level > 0
                   for feature in features.values()):
            continue
        for plan_slot in plan.slots:
            if plan_slot.renderer != "features" or plan_slot.paint_band is None:
                continue
            geometry = record_feature_lane_geometries[record_index]
            for feature in features.values():
                if feature.source_feature_index is None:
                    continue  # by_source_index only holds integer source indexes
                dom_id = feature_dom_index.by_source_index.get((record_index, feature.source_feature_index))
                if dom_id is None:
                    continue
                lane = geometry.lane_for(strand=feature.strand, track_id=feature.feature_track_id,
                                         separate_strands=render_context.profile.strandedness)
                band = lane.band.translate(plan_slot.origin_y)
                if non_feature_band is not None:
                    band = band.union(non_feature_band)
                feature_attachments.setdefault(record_index, {})[dom_id] = LinearFeatureAttachment(
                    comparison_frame_px(band.top_y - attachment_band.top_y),
                    comparison_frame_px(band.bottom_y - attachment_band.bottom_y),
                )

    if draw_length_bar:
        min_record_x = min(
            placement.x for placement in record_placements.values()
        )
        max_record_x = max(
            placement.x + placement.sequence_width
            for placement in record_placements.values()
        )
        length_bar_offset_x = min_record_x
        length_bar_group = render_stage(
            "length_bar_content",
            build_group_content,
            LengthBarGroup,
            # The longest record's drawn width: the bar maps bp to px as that record does.
            max(placement.sequence_width for placement in record_placements.values()),
            canvas_config.longest_genome,
            canvas_config.length_param,
            scale=cfg.objects.scale,
            font_family=cfg.objects.text.font_family,
            ruler_width=max_record_x - min_record_x,
        )
        length_bar_bottom = length_bar_offset_y + float(
            _content_bounds(length_bar_group).max_y
        )
        canvas_config.total_height = max(
            1.0,
            length_bar_bottom + float(canvas_config.vertical_padding),
        )

    canvas: Drawing = canvas_config.create_svg_canvas()
    primary_target_start = len(getattr(canvas, "elements", []))
    if length_bar_group is not None:
        canvas = add_length_bar_on_linear_canvas(
            canvas,
            canvas_config,
            length_bar_group,
            offset_x=length_bar_offset_x,
        )

    comparison_placements = record_placements if has_blast else None

    if has_blast and comparison_placements is not None:
        canvas = add_explicit_comparisons_on_linear_canvas(
            canvas,
            comparison_rows.comparisons,
            canvas_config,
            blast_config,
            records,
            comparison_placements,
            feature_dom_index=feature_dom_index,
            record_transforms=record_transforms,
            record_keys=record_keys,
            feature_attachments=feature_attachments,
        )

    # Each record's drawn definition text, in its feature labels' coordinates (OV-1395).
    definition_label_boxes = [
        _definition_label_boxes(
            record_collision_bands[index],
            record_x=record_placements[index].x,
            features_origin_y=record_vertical_plans[index].features_origin_y,
        )
        for index in range(len(records))
    ]
    # The API checks the drawn labels and offers the rows (LP-D48, LP-D52).
    # Annotation labels move along their spans off the feature labels and leaders (OV-845).
    placed_annotation_labels, all_labels = _place_linear_annotation_labels_by_record(
        records=records,
        axis_lengths=axis_lengths,
        placements=record_placements,
        plans=record_vertical_plans,
        collision_bands=record_collision_bands,
        feature_labels=all_labels,
        feature_dicts=record_feature_dicts,
        feature_geometries=record_feature_lane_geometries,
        separate_strands=render_context.profile.strandedness,
        feature_stroke_width=max(float(feature_config.block_stroke_width), float(feature_config.line_stroke_width)),
        annotation_track_layouts=annotation_track_layouts,
        font_family=str(cfg.objects.text.font_family),
        label_style=render_context.profile.label_style,
        gap=float(cfg.labels.spacing.linear),
        # One assembly more at most: one that places proposed rows proposes none.
        new_rows=feature_label_rows is None,
        label_frames=label_frames,
        definition_boxes=definition_label_boxes,
    )
    annotation_labels_by_record = [placed.by_slot for placed in placed_annotation_labels]
    label_fit = linear_label_fit(
        [(str(record.id), labels) for record, labels in zip(records, all_labels)],
        label_font_size,
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
        frames=label_frames,
        definition_boxes=definition_label_boxes,
    )

    total_records = len(records)
    for count, record in enumerate(records, start=1):
        record_index = count - 1
        depth_by_index = (
            record_depth_by_index[record_index]
            if record_index < len(record_depth_by_index)
            else {}
        )
        offset_y = record_offsets[record_index]
        record_placement = record_placements[record_index]
        sequence_width = record_placement.sequence_width
        record_vertical_plan = record_vertical_plans[record_index]
        record_definition_geometry = record_definition_geometries[record_index]

        offset_x = record_offsets_x[record_index] if record_index < len(record_offsets_x) else 0.0

        labels_for_record = all_labels[record_index] if record_index < len(all_labels) else []
        orthogroup_label_member_ids, orthogroup_label_top_member_ids = orthogroup_label_sets_for_record(
            orthogroup_label_eligibility if label_scope == "orthogroup_top" else None,
            record_index,
        )
        record_group_id = record_group_svg_id(
            record.id,
            mode="linear",
            record_index=record_index,
            record_count=total_records,
        )
        definition_group_id = definition_group_svg_id(
            record.id,
            mode="linear",
            record_index=record_index,
            record_count=total_records,
        )

        if linear_track_layout is not None:
            feature_rendered = False
            axis_row_offsets = record_vertical_plan.row_offsets if record_vertical_plan.wrap_rows else ()

            for slot in sorted(linear_track_layout.slots, key=lambda item: (item.z, item.slot_index)):
                resolved_slot = record_vertical_plan.slot_by_id(slot.id)
                row_shifts = record_vertical_plan.slot_row_shifts(slot.id)
                slot_group_id = (
                    _linear_track_slot_dom_id(
                        slot_id=str(slot.id),
                        renderer=str(slot.renderer),
                        slot_index=int(slot.slot_index),
                        record_index=record_index,
                        record_count=total_records,
                    )
                    if user_slot_mode
                    else None
                )
                if slot.renderer == "spacer":
                    continue
                if slot.renderer == "features":
                    canvas = add_record_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        feature_config,
                        feature_layers=record_feature_layers[record_index],
                        render_context=render_context,
                        feature_dom_index=feature_dom_index,
                        precalculated_labels=labels_for_record,
                        label_font_size=label_font_size,
                        orthogroup_label_member_ids=orthogroup_label_member_ids,
                        orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
                        record_index=record_index,
                        record_count=total_records,
                        group_id=record_group_id,
                        dom_group_id=slot_group_id,
                        slot_id=str(slot.id) if slot_group_id is not None else None,
                        slot_renderer=(
                            str(slot.renderer) if slot_group_id is not None else None
                        ),
                        placement=record_placement,
                        multi_record_layout=multi_record_enabled,
                        record_local_ruler=multi_record_enabled,
                        feature_offset_y=resolved_slot.origin_y,
                        feature_lane_geometry=record_feature_lane_geometries[record_index],
                        record_transform=(record_transforms[record_index] if record_transforms is not None else None),
                        record_translation_x=final_translations[record_index][0],
                        record_translation_y=final_translations[record_index][1],
                        axis_row_offsets=axis_row_offsets,
                        feature_row_shifts=row_shifts,
                        group_key=GroupKey(str(record_keys[record_index]), GroupKind.FEATURES, str(slot.id)),
                    )
                    feature_rendered = True
                    continue

                track_offset_y = resolved_slot.origin_y
                if slot.renderer == "annotations":
                    annotation_layout = annotation_track_layouts.get(slot.id)
                    if annotation_layout is None:
                        continue
                    annotation_content = render_stage(
                        "annotation_content",
                        linear_annotation_track_content,
                        annotation_layout,
                        record_id=str(record.id),
                        record_index=record_index,
                        record_length=axis_lengths[record_index],
                        bar_length_px=float(sequence_width),
                        y_offset_px=0.0,
                        height_px=float(resolved_slot.height),
                        font_family=str(cfg.objects.text.font_family),
                        params=annotation_track_params_from_mapping(slot.params),
                        labels=annotation_labels_by_record[record_index].get(slot.id, ()),
                        dom_group_id=slot_group_id,
                        semantic_slot_id=(
                            str(slot.id) if slot_group_id is not None else None
                        ),
                        row_frames=linear_annotation_row_frames(record_vertical_plan, slot.id),
                    )
                    annotation_group = place_content(
                        canvas,
                        annotation_content,
                        GroupKey(str(record_keys[record_index]), GroupKind.ANNOTATIONS, str(slot.id)),
                    )
                    annotation_group.translate(
                        offset_x + canvas_config.horizontal_offset,
                        offset_y + float(track_offset_y),
                    )
                    canvas.add(annotation_group)
                    continue
                if slot.renderer == "depth":
                    track_index = _depth_slot_track_index(slot)
                    depth_track = depth_by_index.get(track_index)
                    if depth_track is None:
                        continue
                    group_id = _linear_depth_group_id(
                        slot.id or depth_track.id,
                        record_index=record_index,
                        record_count=total_records,
                    )
                    canvas = add_depth_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        depth_track.config,
                        depth_df=depth_track.df,
                        group_id=group_id,
                        axis_group_id=f"{group_id}_axis",
                        dom_group_id=slot_group_id,
                        dom_axis_group_id=(
                            f"{slot_group_id}_axis"
                            if slot_group_id is not None
                            else None
                        ),
                        slot_id=str(slot.id) if slot_group_id is not None else None,
                        slot_renderer=(
                            str(slot.renderer) if slot_group_id is not None else None
                        ),
                        track_height=slot.height,
                        track_offset_y=track_offset_y,
                        sequence_width=sequence_width,
                        record_transform=(record_transforms[record_index] if record_transforms is not None else None),
                        row_shifts=row_shifts,
                        group_key=GroupKey(str(record_keys[record_index]), GroupKind.DEPTH, str(slot.id)),
                    )
                    continue

                if slot.renderer == "dinucleotide_content":
                    nt = slot_dinucleotide(slot.params, str(gc_config.dinucleotide))
                    per_nt_gc_dfs = record_gc_dfs_by_nt.get(nt, record_gc_dfs)
                    shared_gc_df = per_nt_gc_dfs[record_index] if record_index < len(per_nt_gc_dfs) else None
                    slot_gc_config = drawn_gc_config(slot.params, depth_by_index)
                    canvas = add_gc_content_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        slot_gc_config,
                        gc_df=shared_gc_df,
                        track_height=slot.height,
                        track_offset_y=track_offset_y,
                        group_id=make_linear_dom_id(
                            slot.id,
                            record_index=record_index,
                            record_count=total_records,
                        ),
                        dom_group_id=slot_group_id,
                        slot_id=str(slot.id),
                        slot_renderer=str(slot.renderer),
                        sequence_width=sequence_width,
                        record_transform=(record_transforms[record_index] if record_transforms is not None else None),
                        row_shifts=row_shifts,
                        group_key=GroupKey(str(record_keys[record_index]), GroupKind.GC_CONTENT, str(slot.id)),
                    )
                    continue

                if slot.renderer == "dinucleotide_skew":
                    nt = slot_dinucleotide(slot.params, str(skew_config.dinucleotide))
                    per_nt_gc_dfs = record_gc_dfs_by_nt.get(nt, record_gc_dfs)
                    shared_gc_df = per_nt_gc_dfs[record_index] if record_index < len(per_nt_gc_dfs) else None
                    canvas = add_gc_skew_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        _slot_skew_config(skew_config, slot, nt),
                        gc_df=shared_gc_df,
                        track_height=slot.height,
                        track_offset_y=track_offset_y,
                        group_id=make_linear_dom_id(
                            slot.id,
                            record_index=record_index,
                            record_count=total_records,
                        ),
                        dom_group_id=slot_group_id,
                        slot_id=str(slot.id),
                        slot_renderer=str(slot.renderer),
                        sequence_width=sequence_width,
                        record_transform=(record_transforms[record_index] if record_transforms is not None else None),
                        row_shifts=row_shifts,
                        group_key=GroupKey(str(record_keys[record_index]), GroupKind.GC_SKEW, str(slot.id)),
                    )

            # Backed region labels, in front of the record's leader lines (LP-D82 RGN-FRONT).
            for slot in linear_track_layout.slots:
                annotation_layout = annotation_track_layouts.get(slot.id)
                front = (
                    linear_region_label_front_content(
                        annotation_layout,
                        record_index=record_index,
                        font_family=str(cfg.objects.text.font_family),
                        labels=annotation_labels_by_record[record_index].get(slot.id, ()),
                        dom_group_id=(
                            _linear_track_slot_dom_id(
                                slot_id=str(slot.id),
                                renderer=str(slot.renderer),
                                slot_index=int(slot.slot_index),
                                record_index=record_index,
                                record_count=total_records,
                            )
                            if user_slot_mode
                            else None
                        ),
                    )
                    if slot.renderer == "annotations" and annotation_layout is not None
                    else None
                )
                if front is None:
                    continue
                front_group = place_content(
                    canvas, front, GroupKey(str(record_keys[record_index]), GroupKind.LABEL_TEXT, str(slot.id))
                )
                front_group.translate(
                    offset_x + canvas_config.horizontal_offset,
                    offset_y + float(record_vertical_plan.slot_by_id(slot.id).origin_y),
                )
                canvas.add(front_group)

            if not feature_rendered:
                add_record_group(
                    canvas,
                    record,
                    offset_y,
                    offset_x,
                    canvas_config,
                    feature_config,
                    feature_layers=record_feature_layers[record_index],
                    render_context=render_context,
                    feature_dom_index=feature_dom_index,
                    precalculated_labels=None,
                    draw_features=False,
                    label_font_size=label_font_size,
                    orthogroup_label_member_ids=orthogroup_label_member_ids,
                    orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
                    record_index=record_index,
                    record_count=total_records,
                    group_id=record_group_id,
                    placement=record_placement,
                    multi_record_layout=multi_record_enabled,
                    record_local_ruler=multi_record_enabled,
                    feature_lane_geometry=record_feature_lane_geometries[record_index],
                    record_transform=(record_transforms[record_index] if record_transforms is not None else None),
                    record_translation_x=final_translations[record_index][0],
                    record_translation_y=final_translations[record_index][1],
                    axis_row_offsets=axis_row_offsets,
                    group_key=GroupKey(str(record_keys[record_index]), GroupKind.AXIS),
                )
            add_record_definition_group(
                canvas,
                record,
                offset_y,
                offset_x,
                canvas_config,
                definition_column_width,
                group_id=definition_group_id,
                placement=record_placement,
                row_definition_width=row_definition_width,
                definition_center_y=(
                    offset_y + record_definition_geometry.row_center_y
                    if multi_record_enabled
                    and record_definition_geometry.row_center_y is not None
                    else (
                        offset_y + record_definition_geometry.local_center_y
                        if not multi_record_enabled
                        and record_definition_geometry.local_center_y is not None
                        else None
                    )
                ),
                definition_header_center_y=(
                    offset_y + record_definition_geometry.local_center_y
                    if multi_record_enabled
                    and record_definition_geometry.local_center_y is not None
                    else None
                ),
                multi_record_layout=multi_record_enabled,
                local_line_kinds=(
                    local_definition_line_kinds[record_index]
                    if local_definition_line_kinds is not None
                    else None
                ),
                row_line_kinds=row_definition_line_kinds[record_index],
                record_index=record_index,
                record_count=total_records,
                record_presentation=presentations[record_index],
                definition_column_left=definition_column_left,
                record_key=str(record_keys[record_index]),
            )
            continue

    primary_targets = tuple(
        element
        for element in getattr(canvas, "elements", [])[primary_target_start:]
        if isinstance(element, Group)
    )
    # The scale's painted bounds where `add_length_bar_on_linear_canvas` drew
    # it: one value for the primary bounds and the composition metadata.
    length_bar_bounds = (
        _content_bounds(length_bar_group).translated(
            float(canvas_config.horizontal_offset) + float(length_bar_offset_x),
            float(length_bar_offset_y),
        )
        if length_bar_group is not None and _content_bounds(length_bar_group).width > 0.0
        else None
    )
    source_primary_bounds = _collect_linear_primary_bounds(
        records=records,
        record_placements=record_placements,
        record_plans=record_vertical_plans,
        record_collision_bands=record_collision_bands,
        labels_by_record=all_labels,
        record_scalar_axis_rooms=record_scalar_axis_rooms,
        annotation_track_layouts=annotation_track_layouts,
        annotation_labels_by_record=annotation_labels_by_record,
        normalized_comparisons=normalized_comparisons,
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        blast_config=blast_config,
        render_context=render_context,
        cfg=cfg,
        multi_record_enabled=multi_record_enabled,
        length_bar_bounds=length_bar_bounds,
        record_transforms=record_transforms,
    )

    legend_measurement: LegendMeasurement | None = None
    legend_target: Group | None = None
    if canvas_config.legend_position != "none":
        legend_measurement = measure_legend(
            legend_config,
            legend_table,
            placement=canvas_config.legend_position,
            wrap_width=source_primary_bounds.width,
        )
        if (
            legend_measurement.local_bounds.width > 0.0
            and legend_measurement.local_bounds.height > 0.0
        ):
            legend_target = place_content(
                canvas,
                render_stage(
                    "legend_content", linear_legend_content, legend_measurement, legend_table
                ),
                GroupKey(CANVAS_SCOPE, GroupKind.LEGEND),
            )
            canvas.add(legend_target)

    title_target = (
        place_content(canvas, title_content, GroupKey(CANVAS_SCOPE, GroupKind.PLOT_TITLE))
        if title_content is not None
        else None
    )
    if title_target is not None:
        canvas.add(title_target)

    legend_placement = LegendPlacement(str(canvas_config.legend_position))
    title_placement = (
        TitlePlacement(normalized_plot_title_position)
        if title_target is not None
        else TitlePlacement.NONE
    )
    composition_plan = render_stage(
        "composition_plan",
        plan_composition,
        CompositionRequest(
            primary=CompositionItem("primary", source_primary_bounds),
            legend=(
                CompositionItem("legend", legend_measurement.local_bounds)
                if legend_target is not None and legend_measurement is not None
                else None
            ),
            title=(
                CompositionItem("title", title_content.local_bounds)
                if title_target is not None
                and title_content is not None
                and title_content.local_bounds is not None
                else None
            ),
            legend_placement=legend_placement,
            title_placement=title_placement,
        ),
    )
    apply_composition_plan(
        canvas,
        composition_plan,
        primary_targets=primary_targets,
        legend_target=legend_target,
        legend_side=legend_placement,
        legend_reflow_metrics=(
            legend_measurement.reflow_metrics
            if legend_target is not None and legend_measurement is not None
            else None
        ),
        title_target=title_target,
        title_side=title_placement,
        length_bar_local_bounds=length_bar_bounds,
    )

    primary_placement = composition_plan.placement_for("primary")
    if primary_placement is None:  # pragma: no cover - primary is required
        raise RuntimeError("linear composition has no primary placement")

    if linear_track_layout is not None:
        setattr(
            canvas,
            "_gbdraw_track_slot_geometry",
            _serialize_linear_track_slot_geometry(
                records=records,
                layout=linear_track_layout,
                record_plans=record_vertical_plans,
                record_offsets=[
                    float(offset) + primary_placement.dy
                    for offset in record_offsets
                ],
                record_placements=record_placements,
                composition_offset_x=primary_placement.dx,
                horizontal_offset_x=float(canvas_config.horizontal_offset),
                record_collision_bands=record_collision_bands,
                boundary_gap_resolutions=boundary_gap_resolutions,
            ),
        )
    if resolved_feature_slot is not None:
        targets = _linear_feature_placement_slot(
            resolved_feature_slot.side, profile.strandedness
        ).supported_targets()
        for record_geometry in getattr(canvas, "_gbdraw_track_slot_geometry", {}).get("records", []):
            record_geometry["featurePlacementTargets"] = targets
    attach_legend_row_facts(
        canvas,
        LegendRowFacts(
            records=records,
            sources=legend_sources,
            legend_table=legend_table,
        ),
    )
    setattr(canvas, "_gbdraw_alignment_placements", tuple(
        replace(placement, x=placement.x + float(canvas_config.horizontal_offset))
        for placement in alignment_placements
    ))
    setattr(canvas, "_gbdraw_linear_source_content_bounds", source_primary_bounds)
    setattr(canvas, "_gbdraw_linear_composition_plan", composition_plan)
    setattr(canvas, "_gbdraw_comparison_record_id_warnings", comparison_record_id_warnings)
    setattr(canvas, "_gbdraw_comparison_tallies", comparison_rows.tallies)
    setattr(canvas, LINEAR_LABEL_FIT_ATTR, label_fit)

    return canvas


__all__ = ["assemble_linear_diagram"]
