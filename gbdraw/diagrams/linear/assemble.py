#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, replace
import logging
import math
from typing import TYPE_CHECKING, Any

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...analysis.skew import skew_df  # type: ignore[reportMissingImports]
from ...analysis.depth import depth_df as build_depth_df  # type: ignore[reportMissingImports]
from ...analysis.depth_tracks import (  # type: ignore[reportMissingImports]
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
from ...analysis.protein_colinearity import OrthogroupResult  # type: ignore[reportMissingImports]
from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import (  # type: ignore[reportMissingImports]
    GbdrawConfig,
    LinearRenderProfile,
)
from ...exceptions import ValidationError
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
    LegendMeasurement,
)
from ...configurators.gc import _slot_skew_config
from ...core.text import calculate_bbox_dimensions
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...render.groups.linear import LengthBarGroup, LegendGroup, PlotTitleGroup  # type: ignore[reportMissingImports]
from ...render.groups.linear.length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
    auto_linear_tick_interval,
    format_linear_tick_label,
)
from ...io.comparisons import filter_comparison_dataframe, load_comparisons
from ...legend.table import (  # type: ignore[reportMissingImports]
    _unique_legend_key,
    configure_pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
)
from ...layout.linear import (  # type: ignore[reportMissingImports]
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
from ...layout.linear_multi_record import (
    LinearRecordMeasurement,
    LinearRecordPlacement,
    RecordKey,
    solve_linear_layout,
    stable_record_keys,
)
from ...layout.record_placement import resolve_record_row_positions
from ...linear_comparison import (
    LinearComparison,
    merge_linear_comparisons,
    validate_linear_comparison_topology,
)
from ...layout.scalar_axis import linear_scalar_axis_tick_font_size_px  # type: ignore[reportMissingImports]
from ...layout.scalar_axis import (
    _depth_axis_bounds,
    _format_depth_tick,
    format_percent_tick,
    scalar_axis_tick_values,
)
from ...labels.linear import calculate_label_bounds
from ...tracks import (
    LinearTrackSlot,
    ScalarSpec,
    default_linear_track_slots,
    normalize_linear_track_slots_with_axis,
    parse_nonnegative_integer,
)
from ...annotations import (
    AnnotationOptions,
    ResolvedAnnotationBundle,
    ResolvedAnnotationTrack,
    annotation_track_params_from_mapping,
    effective_annotation_style,
    feature_underlay_anchor_slot_id,
    feature_underlay_slot_id,
    layout_annotation_track,
    merge_feature_underlays,
    sync_annotation_legend_entries,
)
from ...annotations.planning import prepare_annotation_track_slots
from ...render.drawers.linear.annotations import draw_linear_annotation_track
from ...render.groups.linear import build_linear_feature_dom_index
from ...render.composition import apply_composition_plan

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
    calculate_orthogroup_alignment_canvas_extents,
    calculate_orthogroup_alignment_offsets,
    orthogroup_label_sets_for_record,
)
from .precalc import (
    FeatureDict,
    _precalculate_definition_metrics,
    _precalculate_feature_layers,
    _precalculate_label_dimensions,
    _resolve_linear_diagram_label_font_size,
)
from ...features.colors import precompute_used_color_rules  # type: ignore[reportMissingImports]
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
)


logger = logging.getLogger(__name__)
_COMPARISON_ENDPOINT_GAP_PX = 4.0

if TYPE_CHECKING:
    from ...api.options import LinearMultiRecordOptions


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


def _prepare_linear_annotation_tracks(
    records: list[SeqRecord],
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    slots: list[LinearTrackSlot] | None,
    *,
    canvas_config: LinearCanvasConfigurator,
    profile: LinearRenderProfile,
    record_depth_tracks: list[list[DepthTrackSpec]] | None,
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
    updated_slots, layouts = _layout_linear_annotation_tracks(
        records,
        slots,
        bundle,
        canvas_config=canvas_config,
        auto_slot_ids=auto_slot_ids,
    )
    return updated_slots, bundle, layouts, auto_slot_ids


def _layout_linear_annotation_tracks(
    records: list[SeqRecord],
    slots: list[LinearTrackSlot],
    bundle: ResolvedAnnotationBundle,
    *,
    canvas_config: LinearCanvasConfigurator,
    auto_slot_ids: frozenset[str],
    sequence_widths: list[float] | None = None,
) -> tuple[list[LinearTrackSlot], dict[str, ResolvedAnnotationTrack]]:
    """Resolve annotation lanes for the axis widths used by final rendering."""

    if sequence_widths is not None and len(sequence_widths) != len(records):
        raise ValueError("sequence_widths must contain one width for every record")
    record_lengths = {index: len(record.seq) for index, record in enumerate(records)}
    if sequence_widths is None:
        bp_per_px = {
            index: len(record.seq)
            / max(
                1.0,
                float(canvas_config.alignment_width)
                * (
                    1.0
                    if canvas_config.normalize_length
                    else len(record.seq) / max(1, canvas_config.longest_genome)
                ),
            )
            for index, record in enumerate(records)
        }
    else:
        bp_per_px = {
            index: len(record.seq) / max(1.0, float(sequence_widths[index]))
            for index, record in enumerate(records)
        }
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
    canvas_config: LinearCanvasConfigurator,
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
    updated_slots = [underlay_slot, *slots]
    private_prefix_count = len(updated_slots) - len(slots)
    updated_auto_ids = frozenset((*auto_slot_ids, slot_id))
    updated_slots, updated_layouts = _layout_linear_annotation_tracks(
        records,
        updated_slots,
        merged,
        canvas_config=canvas_config,
        auto_slot_ids=updated_auto_ids,
    )
    return (
        updated_slots,
        merged,
        updated_layouts,
        updated_auto_ids,
        private_prefix_count,
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
        if 0 <= track_index < len(depth_heights) and depth_heights[track_index] is not None:
            out.append(replace(slot, height=ScalarSpec(float(depth_heights[track_index]), "px")))
        else:
            out.append(slot)
    return out


def _slot_nt(slot: LinearResolvedTrack, default_nt: str) -> str:
    params = slot.params or {}
    return str(params.get("nt", params.get("dinucleotide", default_nt)) or default_nt).upper()


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










def _slot_legend_label(slot, fallback: str) -> str:
    params = getattr(slot, "params", {}) or {}
    raw_label = params.get("legend_label", params.get("label"))
    label = str(raw_label).strip() if raw_label is not None else ""
    return label or fallback


def _sync_legend_table_for_linear_slots(
    legend_table: dict,
    *,
    linear_track_slots: list | None,
    skew_config: GcSkewConfigurator,
    depth_tracks: list[DepthTrackData] | None,
) -> dict:
    """Replace singleton numeric legends with slot-authoritative entries."""
    if linear_track_slots is None:
        return legend_table

    out = sync_depth_track_legend_entries(
        legend_table,
        depth_tracks,
        slots=linear_track_slots,
    )
    default_nt = str(getattr(skew_config, "dinucleotide", "GC")).upper()
    for key in (
        f"{default_nt} skew",
        f"{default_nt} skew (+)",
        f"{default_nt} skew (-)",
    ):
        out.pop(key, None)

    for slot in linear_track_slots:
        if not getattr(slot, "enabled", True):
            continue
        if str(getattr(slot, "renderer", "")) != "dinucleotide_skew":
            continue
        nt = _slot_nt(slot, default_nt)
        label = _slot_legend_label(slot, f"{nt} skew")
        slot_skew_config = _slot_skew_config(skew_config, slot, nt)
        if slot_skew_config.high_fill_color == slot_skew_config.low_fill_color:
            out[_unique_legend_key(out, label)] = {
                "type": "solid",
                "fill": slot_skew_config.high_fill_color,
                "stroke": slot_skew_config.stroke_color,
                "width": slot_skew_config.stroke_width,
            }
        else:
            out[_unique_legend_key(out, f"{label} (+)")] = {
                "type": "solid",
                "fill": slot_skew_config.high_fill_color,
                "stroke": slot_skew_config.stroke_color,
                "width": slot_skew_config.stroke_width,
            }
            out[_unique_legend_key(out, f"{label} (-)")] = {
                "type": "solid",
                "fill": slot_skew_config.low_fill_color,
                "stroke": slot_skew_config.stroke_color,
                "width": slot_skew_config.stroke_width,
            }
    return out


def _serialize_linear_track_slot_geometry(
    *,
    records: list[SeqRecord],
    layout: LinearTrackLayout,
    record_plans: list[LinearRecordVerticalPlan],
    record_offsets: list[float],
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


def _load_linear_comparison_dataframes(
    blast_files,
    comparison_dataframes: list[DataFrame] | None,
    blast_config,
) -> list[DataFrame]:
    if not (blast_files or comparison_dataframes):
        return []

    comparison_sources: list[list[DataFrame]] = []
    if blast_files:
        comparison_sources.append(load_comparisons(blast_files, blast_config))
    if comparison_dataframes:
        comparison_sources.append(
            [
                filter_comparison_dataframe(comparison, blast_config)
                for comparison in comparison_dataframes
            ]
        )

    comparisons: list[DataFrame] = []
    max_source_len = max((len(source) for source in comparison_sources), default=0)
    for index in range(max_source_len):
        frames = [
            source[index]
            for source in comparison_sources
            if index < len(source)
        ]
        if len(frames) == 1:
            comparisons.append(frames[0])
        elif frames:
            comparisons.append(pd.concat(frames, ignore_index=True))
    return comparisons


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
    axis_stroke = 0.4
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
    labels: list[dict],
    canvas_config: LinearCanvasConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config_by_index: dict[int, DepthConfigurator],
    available_depth_indices: set[int],
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
) -> tuple[dict[str, LinearSlotFootprint], float]:
    footprints: dict[str, LinearSlotFootprint] = {}
    label_band = measure_linear_label_band(
        labels,
        leader_stroke_width=float(canvas_config._cfg.labels.stroke_width.long),
    )

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
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=paint_band if feature_dict or labels else None,
                reserve_band=reserve_band,
                data_available=bool(feature_dict or labels),
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
            available = track_index in available_depth_indices
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=expanded if available else None,
                reserve_band=expanded if available else VerticalBand(0.0, 0.0),
                data_available=available,
            )
            continue

        if slot.renderer == "dinucleotide_content":
            overhang = max(
                0.5 * max(0.0, float(gc_config.stroke_width)),
                0.5 * max(0.0, float(gc_config.percent_border_width)),
            )
            if str(getattr(gc_config, "mode", "deviation")).strip().lower() == "percent":
                overhang = max(overhang, _axis_text_overhang(gc_config, float(slot.height)))
            expanded = body_band.expand(overhang)
            footprints[slot.id] = LinearSlotFootprint(expanded, expanded)
            continue

        if slot.renderer == "dinucleotide_skew":
            overhang = 0.5 * max(0.0, float(skew_config.stroke_width))
            expanded = body_band.expand(overhang)
            footprints[slot.id] = LinearSlotFootprint(expanded, expanded)
            continue

        if slot.renderer == "annotations":
            annotation_layout = annotation_track_layouts.get(slot.id)
            placements = (
                [
                    item
                    for item in annotation_layout.placements
                    if item.annotation.record_index == record_index
                ]
                if annotation_layout is not None
                else []
            )
            params = annotation_track_params_from_mapping(slot.params)
            top_overhang = 0.0
            bottom_overhang = 0.0
            if params.show_labels:
                for placement in placements:
                    style = placement.annotation.style
                    label_overhang = float(style.label_offset) + 0.5 * float(
                        style.label_font_size or 10.0
                    )
                    labels_above = slot.side == "above" or (
                        slot.side == "overlay" and params.layer == "underlay"
                    )
                    if labels_above:
                        top_overhang = max(top_overhang, label_overhang)
                    else:
                        bottom_overhang = max(bottom_overhang, label_overhang)
            expanded = body_band.expand(top_overhang, bottom_overhang)
            available = bool(placements)
            footprints[slot.id] = LinearSlotFootprint(
                paint_band=expanded if available else None,
                reserve_band=expanded,
                data_available=available,
            )
            continue

        footprints[slot.id] = LinearSlotFootprint(body_band, body_band)
    feature_center_y = 0.5 * (
        feature_geometry.occupied_band.top_y
        + feature_geometry.occupied_band.bottom_y
    )
    return footprints, feature_center_y


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


def _definition_metrics_by_record(
    records: list[SeqRecord],
    canvas_config: LinearCanvasConfigurator,
    *,
    cfg: GbdrawConfig,
    line_kinds_by_record: list[frozenset[str] | None] | None = None,
) -> tuple[float, list[float], list[float]]:
    """Return the maximum width plus record-local definition widths and heights."""

    widths: list[float] = []
    heights: list[float] = []
    for index, record in enumerate(records):
        line_kinds = (
            line_kinds_by_record[index]
            if line_kinds_by_record is not None
            else None
        )
        width, record_heights, _half_heights = _precalculate_definition_metrics(
            [record],
            canvas_config,
            cfg=cfg,
            line_kinds_by_record=[line_kinds],
        )
        widths.append(float(width))
        heights.append(float(record_heights[0]))
    return max(widths, default=0.0), widths, heights


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
    labels_by_record: list[list[dict]],
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
    gc_config: GcContentConfigurator,
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
    for record_index, _record in enumerate(records):
        depth_by_index = (
            record_depth_by_index[record_index]
            if record_index < len(record_depth_by_index)
            else {}
        )
        depth_configs = dict(representative_depth_configs)
        depth_configs.update(
            {index: track.config for index, track in depth_by_index.items()}
        )
        footprints, feature_center_y = _linear_slot_footprints_for_record(
            record_index=record_index,
            layout=layout,
            feature_dict=feature_dicts[record_index],
            feature_geometry=feature_geometries[record_index],
            labels=labels_by_record[record_index],
            canvas_config=canvas_config,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_config_by_index=depth_configs,
            available_depth_indices=set(depth_by_index),
            annotation_track_layouts=annotation_track_layouts,
        )
        plan = resolve_linear_record_vertical_plan(
            layout,
            axis_band=axis_band,
            footprints=footprints,
        )
        feature_origin_y = (
            plan.slot_by_id(feature_slot_id).origin_y
            if feature_slot_id is not None
            else 0.0
        )
        feature_center_y += feature_origin_y
        canvas_band = plan.canvas_band
        local_definition_height = (
            float(definition_heights[record_index])
            if record_index < len(definition_heights)
            else 0.0
        )
        row_definition_height = (
            float(row_definition_heights[record_index])
            if record_index < len(row_definition_heights)
            else 0.0
        )
        if multi_record_enabled:
            row_band = (
                _centered_vertical_band(feature_center_y, row_definition_height)
                if split_row_definitions
                else None
            )
            if row_band is not None:
                canvas_band = canvas_band.union(row_band)
            local_band = None
            if local_definition_height > 0.0:
                header_bottom_y = (
                    canvas_band.top_y - float(canvas_config.vertical_padding)
                )
                local_band = VerticalBand(
                    header_bottom_y - local_definition_height,
                    header_bottom_y,
                )
        else:
            local_band = _centered_vertical_band(
                feature_center_y,
                local_definition_height,
            )
            row_band = None

        for definition_band in (local_band,):
            if definition_band is not None:
                canvas_band = canvas_band.union(definition_band)
        plans.append(replace(plan, canvas_band=canvas_band))
        definition_geometries.append(
            _LinearRecordDefinitionGeometry(
                local_band=local_band,
                row_band=row_band,
                local_width=(
                    float(definition_widths[record_index])
                    if record_index < len(definition_widths)
                    else 0.0
                ),
                row_width=(
                    float(row_definition_widths[record_index])
                    if record_index < len(row_definition_widths)
                    else 0.0
                ),
            )
        )
    return plans, definition_geometries


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

    local_band = definition_geometry.local_band
    local_width = max(0.0, float(definition_geometry.local_width))
    if local_band is not None and local_width > 0.0:
        if multi_record_enabled:
            center_x = x + (0.5 * width)
            definition_start = center_x - (0.5 * local_width)
            definition_end = center_x + (0.5 * local_width)
        elif keep_definition_left_aligned:
            definition_start = -(
                max(0.0, float(definition_column_width))
                + max(0.0, float(definition_gap))
            )
            definition_end = definition_start + local_width
        else:
            definition_end = x - max(0.0, float(definition_gap))
            definition_start = definition_end - local_width
        bands.append(
            CollisionBand(
                "definition",
                definition_start,
                definition_end,
                local_band.top_y,
                local_band.bottom_y,
            )
        )

    row_band = definition_geometry.row_band
    actual_row_width = max(0.0, float(definition_geometry.row_width))
    if row_band is not None and actual_row_width > 0.0:
        row_start = -(
            max(0.0, float(definition_gap))
            + max(0.0, float(row_definition_width))
        )
        bands.append(
            CollisionBand(
                "definition",
                row_start,
                row_start + actual_row_width,
                row_band.top_y,
                row_band.bottom_y,
            )
        )
    return tuple(bands)


def _linear_axis_ruler_bounds(
    record: SeqRecord,
    placement: LinearRecordPlacement,
    plan: LinearRecordVerticalPlan,
    *,
    canvas_config: LinearCanvasConfigurator,
    render_context: LinearRecordRenderContext,
    cfg: GbdrawConfig,
    record_local_ruler: bool,
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
        else auto_linear_tick_interval(max(1, int(canvas_config.longest_genome)))
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
    if not ticks:
        return None

    bounds: list[Aabb] = []
    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param(
        canvas_config.length_param
    )
    for coordinate in ticks:
        x = placement.sequence_width * (
            abs(coordinate - start_coord) / float(span)
        )
        label = format_linear_tick_label(
            coordinate,
            context_length=max(1, int(canvas_config.longest_genome)),
            tick_interval=tick_interval,
        )
        width, _height = calculate_bbox_dimensions(
            label,
            cfg.objects.text.font_family,
            font_size,
            cfg.canvas.dpi,
        )
        bounds.append(
            Aabb(
                placement.x + x - (0.5 * width),
                placement.axis_y + plan.axis_band.top_y,
                placement.x + x + (0.5 * width),
                placement.axis_y + plan.axis_band.bottom_y,
            )
        )
    return union_aabbs(bounds)


def _linear_scalar_axis_left_overhang(
    plan: LinearRecordVerticalPlan,
    *,
    gc_config: GcContentConfigurator,
    depth_by_index: dict[int, DepthTrackData],
    cfg: GbdrawConfig,
) -> float:
    """Measure left-facing scalar-axis ticks and labels from resolved data."""

    overhang = 0.0
    for slot in plan.slots:
        if slot.paint_band is None:
            continue
        labels: list[str] = []
        axis_config: object | None = None
        if slot.renderer == "depth":
            depth_track = depth_by_index.get(_depth_slot_track_index(slot))
            if depth_track is None:
                continue
            axis_config = depth_track.config
            if not bool(getattr(axis_config, "show_axis", False)):
                continue
            axis_min, axis_max = _depth_axis_bounds(
                depth_track.df,
                getattr(axis_config, "min_depth", None),
                getattr(axis_config, "max_depth", None),
            )
            labels = [
                _format_depth_tick(value)
                for value in scalar_axis_tick_values(
                    axis_min,
                    axis_max,
                    show_ticks=bool(getattr(axis_config, "show_ticks", False)),
                    large_tick_interval=getattr(
                        axis_config,
                        "large_tick_interval",
                        None,
                    ),
                )
                if not (
                    getattr(axis_config, "min_depth", None) is None
                    and math.isclose(value, axis_min, rel_tol=1e-9, abs_tol=1e-9)
                )
            ]
        elif (
            slot.renderer == "dinucleotide_content"
            and str(getattr(gc_config, "mode", "deviation")).strip().lower()
            == "percent"
        ):
            axis_config = gc_config
            if not bool(getattr(axis_config, "show_axis", False)):
                continue
            labels = [
                format_percent_tick(value)
                for value in scalar_axis_tick_values(
                    float(getattr(gc_config, "min_percent", 0.0) or 0.0),
                    float(getattr(gc_config, "max_percent", 100.0) or 100.0),
                    show_ticks=bool(getattr(gc_config, "show_ticks", False)),
                    large_tick_interval=getattr(
                        gc_config,
                        "large_tick_interval",
                        None,
                    ),
                )
            ]
        if axis_config is None:
            continue
        overhang = max(overhang, 0.4)
        font_size = linear_scalar_axis_tick_font_size_px(
            axis_config,
            float(slot.height),
        )
        for label in labels:
            width, _height = calculate_bbox_dimensions(
                label,
                cfg.objects.text.font_family,
                font_size,
                cfg.canvas.dpi,
            )
            overhang = max(overhang, 4.5 + float(width))
    return overhang


def _linear_annotation_bounds(
    record_index: int,
    record_length: int,
    placement: LinearRecordPlacement,
    plan: LinearRecordVerticalPlan,
    *,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    cfg: GbdrawConfig,
) -> list[Aabb]:
    """Return structured mark and label bounds for rendered annotation slots."""

    scale = placement.sequence_width / max(1, int(record_length))
    bounds: list[Aabb] = []
    for slot in plan.slots:
        if slot.renderer != "annotations" or slot.paint_band is None:
            continue
        track = annotation_track_layouts.get(slot.id)
        if track is None:
            continue
        params = annotation_track_params_from_mapping(slot.params)
        for item in track.placements:
            annotation = item.annotation
            if annotation.record_index != record_index:
                continue
            style = effective_annotation_style(annotation, params)
            half_stroke = 0.5 * max(0.0, float(style.stroke_width))
            min_x = min(float(start) for start, _end in annotation.segments) * scale
            max_x = max(float(end) for _start, end in annotation.segments) * scale
            if params.show_labels and annotation.label:
                font_size = float(style.label_font_size or 10.0)
                width, _height = calculate_bbox_dimensions(
                    annotation.label,
                    cfg.objects.text.font_family,
                    font_size,
                    cfg.canvas.dpi,
                )
                if style.label_position == "start":
                    label_x = float(annotation.segments[0][0]) * scale
                    label_min_x, label_max_x = label_x, label_x + width
                elif style.label_position == "end":
                    label_x = float(annotation.segments[-1][1]) * scale
                    label_min_x, label_max_x = label_x - width, label_x
                else:
                    label_x = float(annotation.midpoint_bp) * scale
                    label_min_x = label_x - (0.5 * width)
                    label_max_x = label_x + (0.5 * width)
                min_x = min(min_x, label_min_x)
                max_x = max(max_x, label_max_x)
            bounds.append(
                Aabb(
                    placement.x + min_x - half_stroke,
                    placement.axis_y + slot.paint_band.top_y,
                    placement.x + max_x + half_stroke,
                    placement.axis_y + slot.paint_band.bottom_y,
                )
            )
    return bounds


def _collect_linear_primary_bounds(
    *,
    records: list[SeqRecord],
    record_placements: dict[int, LinearRecordPlacement],
    record_plans: list[LinearRecordVerticalPlan],
    record_collision_bands: list[tuple[CollisionBand, ...]],
    labels_by_record: list[list[dict]],
    record_depth_by_index: list[dict[int, DepthTrackData]],
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    normalized_comparisons: list[LinearComparison],
    canvas_config: LinearCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    blast_config: object,
    render_context: LinearRecordRenderContext,
    cfg: GbdrawConfig,
    multi_record_enabled: bool,
    length_bar_group: LengthBarGroup | None,
    length_bar_offset_x: float,
    length_bar_offset_y: float,
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
    label_half_stroke = 0.5 * max(
        0.0,
        float(cfg.labels.stroke_width.long),
    )
    bounds: list[Aabb] = []
    for record_index, record in enumerate(records):
        placement = record_placements[record_index]
        plan = record_plans[record_index]
        scalar_overhang = _linear_scalar_axis_left_overhang(
            plan,
            gc_config=gc_config,
            depth_by_index=record_depth_by_index[record_index],
            cfg=cfg,
        )
        bounds.append(
            Aabb(
                horizontal_offset + placement.x - max(half_stroke, scalar_overhang),
                placement.axis_y + plan.paint_band.top_y - half_stroke,
                horizontal_offset + placement.x + placement.sequence_width + half_stroke,
                placement.axis_y + plan.paint_band.bottom_y + half_stroke,
            )
        )

        for band in record_collision_bands[record_index]:
            if band.kind != "definition" or band.width <= 0.0:
                continue
            bounds.append(
                Aabb(
                    horizontal_offset + band.x_start,
                    placement.axis_y + band.top_y,
                    horizontal_offset + band.x_end,
                    placement.axis_y + band.bottom_y,
                )
            )

        feature_origin_y = next(
            (
                float(slot.origin_y)
                for slot in plan.slots
                if slot.renderer == "features"
            ),
            0.0,
        )
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
            leader_points = [
                (
                    float(label[name_x]),
                    float(label[name_y]),
                )
                for name_x, name_y in (
                    ("leader_start_x", "leader_start_y"),
                    ("leader_end_x", "leader_end_y"),
                )
                if label.get(name_x) is not None and label.get(name_y) is not None
            ]
            if leader_points:
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
        )
        if ruler_bounds is not None:
            bounds.append(ruler_bounds.translated(horizontal_offset, 0.0))
        bounds.extend(
            item.translated(horizontal_offset, 0.0)
            for item in _linear_annotation_bounds(
                record_index,
                len(record.seq),
                placement,
                plan,
                annotation_track_layouts=annotation_track_layouts,
                cfg=cfg,
            )
        )

    comparison_half_stroke = 0.5 * max(
        0.0,
        float(getattr(blast_config, "stroke_width", 0.0)),
    )
    for comparison in normalized_comparisons:
        query = record_placements[comparison.query_record_index]
        subject = record_placements[comparison.subject_record_index]
        query_y = (
            query.comparison_bottom_y
            if query.row < subject.row
            else query.comparison_top_y
        )
        subject_y = (
            subject.comparison_bottom_y
            if subject.row < query.row
            else subject.comparison_top_y
        )
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

    if length_bar_group is not None and length_bar_group.local_bounds.width > 0.0:
        bounds.append(
            length_bar_group.local_bounds.translated(
                horizontal_offset + float(length_bar_offset_x),
                float(length_bar_offset_y),
            )
        )
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
    return [skew_df(record, window, step, dinucleotide) for record in records]






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
    comparison_dataframes: list[DataFrame] | None = None,
    linear_comparisons: list[LinearComparison] | None = None,
    linear_layout: LinearMultiRecordOptions | None = None,
    orthogroups: OrthogroupResult | None = None,
    align_orthogroup_feature: str | None = None,
) -> Drawing:
    """
    Assembles a linear diagram of genomic records with optional BLAST comparison data,
    and returns the SVG canvas (not saved).
    """
    profile = canvas_config.profile
    cfg = profile.config
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
    split_row_definitions = (
        multi_record_enabled and bool(canvas_config.keep_definition_left_aligned)
    )
    if multi_record_enabled and bool(cfg.canvas.linear.normalize_length):
        raise ValidationError(
            "normalize_length=True cannot be combined with multiple records in one Linear row."
        )
    if record_depth_tracks is None and depth_tables:
        record_depth_tracks = normalize_depth_tracks(
            records,
            depth_track_tables=[[table] for table in depth_tables],
        )
    user_slot_mode = linear_track_slots is not None
    (
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
    ) = _prepare_linear_annotation_tracks(
        records,
        annotations,
        linear_track_slots,
        canvas_config=canvas_config,
        profile=profile,
        record_depth_tracks=record_depth_tracks,
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

    plot_title_obj: PlotTitleGroup | None = None
    if normalized_plot_title:
        plot_title_obj = PlotTitleGroup(
            normalized_plot_title,
            font_size=float(plot_title_font_size),
            cfg=cfg,
        )

    legacy_comparison_frames = _load_linear_comparison_dataframes(
        blast_files,
        comparison_dataframes,
        blast_config,
    )
    normalized_comparisons: list[LinearComparison] = []
    if linear_comparisons:
        normalized_comparisons.extend(
            LinearComparison(
                item.query_record_index,
                item.subject_record_index,
                filter_comparison_dataframe(item.matches, blast_config),
            )
            for item in linear_comparisons
        )
    if legacy_comparison_frames:
        if multi_record_enabled:
            raise ValidationError(
                "Ordered BLAST/protein comparison inputs are ambiguous when a Linear row "
                "contains multiple records; provide explicit LinearComparison endpoints."
            )
        normalized_comparisons.extend(
            LinearComparison(index, index + 1, frame)
            for index, frame in enumerate(legacy_comparison_frames)
            if index + 1 < len(records)
        )
    normalized_comparisons = list(merge_linear_comparisons(normalized_comparisons))
    validate_linear_comparison_topology(normalized_comparisons, rows_by_record)
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

    record_feature_layers = _precalculate_feature_layers(
        records,
        feature_config,
        canvas_config,
        profile,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
    )
    record_feature_dicts = [
        result.foreground_features for result in record_feature_layers
    ]
    feature_dom_index = build_linear_feature_dom_index(record_feature_dicts)
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
        canvas_config=canvas_config,
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
    _unused_label_height, all_labels, _record_label_heights_above = _precalculate_label_dimensions(
        records,
        feature_config,
        canvas_config,
        render_context=render_context,
        precomputed_feature_dicts=record_feature_dicts,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
        feature_lane_geometries=record_feature_lane_geometries,
    )
    local_definition_line_kinds = (
        [
            (
                frozenset({"replicon", "accession", "length"})
                if index in row_leading_indices
                else None
            )
            for index in range(len(records))
        ]
        if split_row_definitions
        else None
    )
    max_def_width, definition_widths, definition_heights = _definition_metrics_by_record(
        records,
        canvas_config,
        cfg=cfg,
        line_kinds_by_record=local_definition_line_kinds,
    )
    row_definition_width = 0.0
    row_definition_widths = [0.0 for _record in records]
    row_definition_heights = [0.0 for _record in records]
    if split_row_definitions:
        (
            row_definition_width,
            row_definition_widths,
            row_definition_heights,
        ) = _definition_metrics_by_record(
            records,
            canvas_config,
            cfg=cfg,
            line_kinds_by_record=[
                frozenset({"name", "subtitle"})
                if index in row_leading_indices
                else frozenset()
                for index in range(len(records))
            ],
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
    record_vertical_plans, record_definition_geometries = _build_linear_record_vertical_plans(
        records=records,
        layout=linear_track_layout,
        feature_dicts=record_feature_dicts,
        feature_geometries=record_feature_lane_geometries,
        labels_by_record=all_labels,
        canvas_config=canvas_config,
        cfg=cfg,
        gc_config=gc_config,
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
    )

    canvas_config.vertical_offset = _linear_record_vertical_offset(
        record_vertical_plans,
        rows_by_record,
        canvas_config,
    )

    normalize_length = cfg.canvas.linear.normalize_length
    needed_nts = {
        _slot_nt(slot, str(gc_config.dinucleotide))
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
    configure_pairwise_identity_legend_from_comparisons(blast_config, comparisons)
    legend_table: dict = {}
    if canvas_config.legend_position != "none":
        color_map = feature_config.specific_color_rules
        default_color_map = feature_config.default_color_map
        features_present = check_feature_presence(
            records,
            feature_config.selected_features_set,
            feature_visibility_rules=feature_config.feature_visibility_rules,
            specific_color_rules=color_map,
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            records,
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        legend_table = prepare_legend_table(
            gc_config, skew_config, feature_config, features_present, blast_config, has_blast,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=(
                depth_config
                if render_context.depth_enabled
                and depth_track_count(record_depth_tracks) == 1
                else None
            ),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            show_depth=(
                render_context.depth_enabled
                and depth_track_count(record_depth_tracks) == 1
            ),
        )
        legend_table = _sync_legend_table_for_linear_slots(
            legend_table,
            linear_track_slots=normalized_linear_track_slots,
            skew_config=skew_config,
            depth_tracks=(
                representative_depth_data
                if render_context.depth_enabled
                else None
            ),
        )
        legend_table = sync_annotation_legend_entries(
            legend_table,
            resolved_annotations,
            normalized_linear_track_slots,
        )
    definition_reserve_width = (
        row_definition_width
        if split_row_definitions
        else (0.0 if multi_record_enabled else max_def_width)
    )
    canvas_config.configure_plot_width(definition_reserve_width)
    ordinary_row_gap = float(canvas_config.cds_padding) * 1.5
    definition_clear_gap = max(1.0, 0.5 * float(canvas_config.vertical_padding))
    comparison_endpoint_gap_px = (
        _COMPARISON_ENDPOINT_GAP_PX if has_blast else 0.0
    )
    comparison_boundary_records = _comparison_records_by_boundary(
        normalized_comparisons,
        rows_by_record,
    )

    if multi_record_enabled:
        if align_orthogroup_feature:
            raise ValidationError(
                "align_orthogroup_feature is not supported with multiple records in one Linear row."
            )
        orthogroup_alignment_offsets = {index: 0.0 for index in range(len(records))}
    else:
        orthogroup_alignment_offsets = calculate_orthogroup_alignment_offsets(
            records,
            comparisons,
            canvas_config,
            align_orthogroup_feature,
            orthogroups=orthogroups,
        )
    alignment_extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        orthogroup_alignment_offsets,
    )

    record_offsets_x: list[float] = []
    if not multi_record_enabled:
        for record_index, record in enumerate(records):
            if normalize_length:
                record_offset_x = 0.0
            elif canvas_config.align_center:
                record_offset_x = (
                    canvas_config.alignment_width
                    * (
                        (canvas_config.longest_genome - len(record.seq))
                        / canvas_config.longest_genome
                    )
                    / 2
                )
            else:
                record_offset_x = 0.0
            record_offset_x += orthogroup_alignment_offsets.get(record_index, 0.0)
            record_offsets_x.append(record_offset_x)

    definition_column_width = max_def_width
    if canvas_config.keep_definition_left_aligned and record_offsets_x:
        definition_column_width = max(
            0.0,
            float(max_def_width) - min(record_offsets_x),
        )

    record_offsets: list[float] = [0.0 for _record in records]
    record_collision_bands: list[tuple[CollisionBand, ...]] = [
        () for _record in records
    ]
    boundary_gap_resolutions: list[AxisGapResolution] = []
    multi_record_plan = None

    def resolve_first_axis_y() -> float:
        return max(
            canvas_config.vertical_offset,
            canvas_config.original_vertical_offset,
        )

    current_y = resolve_first_axis_y()

    if multi_record_enabled:
        def build_measurements(
            horizontal_plan=None,
        ) -> list[LinearRecordMeasurement]:
            measurements: list[LinearRecordMeasurement] = []
            for index, record in enumerate(records):
                record_plan = record_vertical_plans[index]
                placement = (
                    horizontal_plan.placement_for_index(index)
                    if horizontal_plan is not None
                    else None
                )
                collision_bands: tuple[CollisionBand, ...] = ()
                if placement is not None:
                    absolute_bands = _record_collision_bands(
                        plan=record_plan,
                        definition_geometry=record_definition_geometries[index],
                        sequence_width=placement.sequence_width,
                        record_x=placement.x,
                        multi_record_enabled=True,
                        keep_definition_left_aligned=bool(
                            canvas_config.keep_definition_left_aligned
                        ),
                        definition_column_width=definition_column_width,
                        row_definition_width=row_definition_width,
                        definition_gap=float(canvas_config.definition_gap),
                    )
                    collision_bands = tuple(
                        band.translate(x=-placement.x)
                        for band in absolute_bands
                    )
                measurements.append(
                    LinearRecordMeasurement(
                        record_index=index,
                        record_key=RecordKey(str(record_keys[index])),
                        sequence_length=len(record.seq),
                        left_inset=0.0,
                        right_inset=0.0,
                        top_extent=record_plan.canvas_top_extent,
                        bottom_extent=record_plan.canvas_bottom_extent,
                        comparison_top_extent=record_plan.comparison_top_extent,
                        comparison_bottom_extent=record_plan.comparison_bottom_extent,
                        collision_bands=collision_bands,
                    )
                )
            return measurements

        def solve_measurements(
            measurements: list[LinearRecordMeasurement],
            *,
            first_axis_y: float,
        ):
            return solve_linear_layout(
                measurements,
                rows_by_record,
                available_width=float(canvas_config.alignment_width),
                record_gap_px=(
                    linear_layout.record_gap_px if linear_layout is not None else 24.0
                ),
                align_center=bool(canvas_config.align_center),
                first_axis_y=first_axis_y,
                row_gap_px=float(canvas_config.cds_padding) * 1.5,
                comparison_height=float(canvas_config.configured_comparison_height),
                comparison_endpoint_gap_px=comparison_endpoint_gap_px,
                definition_clear_gap=definition_clear_gap,
                comparison_record_indices_by_boundary=comparison_boundary_records,
                record_order=ordered_record_indices,
            )

        multi_record_plan = solve_measurements(
            build_measurements(),
            first_axis_y=current_y,
        )
        # Horizontal widths are independent of vertical extents. Re-resolve all
        # width-sensitive lanes before the authoritative vertical solve.
        final_sequence_widths = [
            multi_record_plan.placement_for_index(index).sequence_width
            for index in range(len(records))
        ]
        if annotation_track_layouts:
            linear_track_slots, annotation_track_layouts = (
                _layout_linear_annotation_tracks(
                    records,
                    linear_track_slots,
                    resolved_annotations,
                    canvas_config=canvas_config,
                    auto_slot_ids=auto_annotation_slot_ids,
                    sequence_widths=final_sequence_widths,
                )
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
            render_context = LinearRecordRenderContext(
                profile=profile,
                track_layout=linear_track_layout,
                depth_available=depth_available,
            )
            axis_ruler_enabled = render_context.axis_ruler_enabled
            resolved_depth_heights_by_index = _resolved_linear_depth_heights(
                linear_track_layout
            )
            primary_depth_track_index = min(
                resolved_depth_heights_by_index,
                default=0,
            )
        _unused_height, all_labels, _record_label_heights_above = _precalculate_label_dimensions(
            records,
            feature_config,
            canvas_config,
            render_context=render_context,
            precomputed_feature_dicts=record_feature_dicts,
            orthogroup_label_eligibility=orthogroup_label_eligibility,
            feature_lane_geometries=record_feature_lane_geometries,
            sequence_widths=final_sequence_widths,
        )
        record_vertical_plans, record_definition_geometries = _build_linear_record_vertical_plans(
            records=records,
            layout=linear_track_layout,
            feature_dicts=record_feature_dicts,
            feature_geometries=record_feature_lane_geometries,
            labels_by_record=all_labels,
            canvas_config=canvas_config,
            cfg=cfg,
            gc_config=gc_config,
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
        )
        canvas_config.vertical_offset = _linear_record_vertical_offset(
            record_vertical_plans,
            rows_by_record,
            canvas_config,
        )
        current_y = resolve_first_axis_y()
        multi_record_plan = solve_measurements(
            build_measurements(multi_record_plan),
            first_axis_y=current_y,
        )
        record_offsets = [
            multi_record_plan.placement_for_index(index).axis_y
            for index in range(len(records))
        ]
        record_offsets_x = [
            multi_record_plan.placement_for_index(index).x
            for index in range(len(records))
        ]
        record_collision_bands = [
            _record_collision_bands(
                plan=record_vertical_plans[index],
                definition_geometry=record_definition_geometries[index],
                sequence_width=multi_record_plan.placement_for_index(index).sequence_width,
                record_x=multi_record_plan.placement_for_index(index).x,
                multi_record_enabled=True,
                keep_definition_left_aligned=bool(
                    canvas_config.keep_definition_left_aligned
                ),
                definition_column_width=definition_column_width,
                row_definition_width=row_definition_width,
                definition_gap=float(canvas_config.definition_gap),
            )
            for index in range(len(records))
        ]
        boundary_gap_resolutions = list(multi_record_plan.row_gap_resolutions)
        current_y = max(record_offsets)
    else:
        for index, record in enumerate(records):
            sequence_width = (
                float(canvas_config.alignment_width)
                if normalize_length
                else float(canvas_config.alignment_width)
                * len(record.seq)
                / max(1, canvas_config.longest_genome)
            )
            record_collision_bands[index] = _record_collision_bands(
                plan=record_vertical_plans[index],
                definition_geometry=record_definition_geometries[index],
                sequence_width=sequence_width,
                record_x=record_offsets_x[index],
                multi_record_enabled=False,
                keep_definition_left_aligned=bool(
                    canvas_config.keep_definition_left_aligned
                ),
                definition_column_width=definition_column_width,
                row_definition_width=row_definition_width,
                definition_gap=float(canvas_config.definition_gap),
            )

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

    length_bar_group: LengthBarGroup | None = None
    length_bar_offset_x = 0.0
    if (
        cfg.objects.scale.show
        and not multi_record_enabled
        and not canvas_config.normalize_length
        and not axis_ruler_enabled
    ):
        length_bar_offset_x = alignment_extents.ruler_offset_x
        length_bar_group = LengthBarGroup(
            canvas_config.fig_width,
            canvas_config.alignment_width,
            canvas_config.longest_genome,
            canvas_config,
            cfg=cfg,
            ruler_width=alignment_extents.ruler_width,
        )
    painted_content_bottom = max(
        (
            float(record_offsets[index]) + plan.canvas_band.bottom_y
            for index, plan in enumerate(record_vertical_plans)
        ),
        default=0.0,
    )
    canvas_config.height_below_final_record = (
        painted_content_bottom
        + 4 * canvas_config.vertical_padding
    )
    length_bar_offset_y = float(canvas_config.height_below_final_record)
    length_bar_bottom = (
        length_bar_offset_y + float(length_bar_group.local_bounds.max_y)
        if length_bar_group is not None
        else painted_content_bottom
    )
    canvas_config.total_height = max(
        1.0,
        length_bar_bottom + float(canvas_config.vertical_padding),
    )
    record_placements: dict[int, LinearRecordPlacement] = {}
    for record_index, record in enumerate(records):
        if multi_record_plan is not None:
            record_placements[record_index] = multi_record_plan.placement_for_index(
                record_index
            )
            continue
        axis_y = record_offsets[record_index]
        record_plan = record_vertical_plans[record_index]
        sequence_width = (
            float(canvas_config.alignment_width)
            if canvas_config.normalize_length
            else float(canvas_config.alignment_width)
            * len(record.seq)
            / max(1, canvas_config.longest_genome)
        )
        record_placements[record_index] = LinearRecordPlacement(
            record_index=record_index,
            record_key=RecordKey(str(record_keys[record_index])),
            row=rows_by_record[record_index],
            column=0,
            x=record_offsets_x[record_index],
            axis_y=axis_y,
            sequence_width=sequence_width,
            left_inset=0.0,
            right_inset=0.0,
            top_extent=record_plan.canvas_top_extent,
            bottom_extent=record_plan.canvas_bottom_extent,
            comparison_top_y=(
                axis_y
                - record_plan.comparison_top_extent
                - comparison_endpoint_gap_px
            ),
            comparison_bottom_y=(
                axis_y
                + record_plan.comparison_bottom_extent
                + comparison_endpoint_gap_px
            ),
            px_per_bp=sequence_width / max(1, len(record.seq)),
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
            normalized_comparisons,
            canvas_config,
            blast_config,
            records,
            comparison_placements,
            feature_dom_index=feature_dom_index,
        )

    label_font_size = _resolve_linear_diagram_label_font_size(
        records,
        canvas_config=canvas_config,
        profile=profile,
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

            for slot in sorted(linear_track_layout.slots, key=lambda item: (item.z, item.slot_index)):
                resolved_slot = record_vertical_plan.slot_by_id(slot.id)
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
                    )
                    feature_rendered = True
                    continue

                track_offset_y = resolved_slot.origin_y
                if slot.renderer == "annotations":
                    annotation_layout = annotation_track_layouts.get(slot.id)
                    if annotation_layout is None:
                        continue
                    params = annotation_track_params_from_mapping(slot.params)
                    annotation_height = float(slot.height)
                    if params.cover_anchor and any(
                        anchor.renderer == "features"
                        and anchor.id == params.anchor_slot
                        for anchor in linear_track_layout.slots
                    ):
                        feature_band = record_feature_lane_geometries[
                            record_index
                        ].occupied_band
                        if feature_band.height > 0.0:
                            track_offset_y = float(feature_band.top_y)
                            annotation_height = float(feature_band.height)
                    bar_length = (
                        float(sequence_width)
                        if sequence_width is not None
                        else float(canvas_config.alignment_width)
                        * (
                            1.0
                            if canvas_config.normalize_length
                            else len(record.seq) / max(1, canvas_config.longest_genome)
                        )
                    )
                    annotation_group = draw_linear_annotation_track(
                        canvas,
                        annotation_layout,
                        record_id=str(record.id),
                        record_index=record_index,
                        record_length=len(record.seq),
                        bar_length_px=bar_length,
                        y_offset_px=0.0,
                        side=slot.side,
                        height_px=annotation_height,
                        font_family=str(cfg.objects.text.font_family),
                        params=params,
                        dom_group_id=slot_group_id,
                        semantic_slot_id=(
                            str(slot.id) if slot_group_id is not None else None
                        ),
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
                    )
                    continue

                if slot.renderer == "dinucleotide_content":
                    nt = _slot_nt(slot, str(gc_config.dinucleotide))
                    per_nt_gc_dfs = record_gc_dfs_by_nt.get(nt, record_gc_dfs)
                    shared_gc_df = per_nt_gc_dfs[record_index] if record_index < len(per_nt_gc_dfs) else None
                    primary_depth_track = depth_by_index.get(primary_depth_track_index)
                    slot_gc_config = _gc_config_matching_linear_depth_axis_font_size(
                        gc_config=_clone_gc_config_with_dinucleotide(gc_config, nt),
                        depth_config=(
                            primary_depth_track.config
                            if primary_depth_track is not None
                            else representative_depth_configs.get(
                                primary_depth_track_index,
                                depth_config,
                            )
                        ),
                        depth_height=resolved_depth_heights_by_index.get(
                            primary_depth_track_index,
                            float(canvas_config.default_depth_height),
                        ),
                        depth_enabled=render_context.depth_enabled,
                    )
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
                    )
                    continue

                if slot.renderer == "dinucleotide_skew":
                    nt = _slot_nt(slot, str(skew_config.dinucleotide))
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
                    )

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
                record_index=record_index,
                record_count=total_records,
            )
            continue

    primary_targets = tuple(
        element
        for element in getattr(canvas, "elements", [])[primary_target_start:]
        if isinstance(element, Group)
    )
    source_primary_bounds = _collect_linear_primary_bounds(
        records=records,
        record_placements=record_placements,
        record_plans=record_vertical_plans,
        record_collision_bands=record_collision_bands,
        labels_by_record=all_labels,
        record_depth_by_index=record_depth_by_index,
        annotation_track_layouts=annotation_track_layouts,
        normalized_comparisons=normalized_comparisons,
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        blast_config=blast_config,
        render_context=render_context,
        cfg=cfg,
        multi_record_enabled=multi_record_enabled,
        length_bar_group=length_bar_group,
        length_bar_offset_x=length_bar_offset_x,
        length_bar_offset_y=length_bar_offset_y,
    )

    legend_measurement: LegendMeasurement | None = None
    legend_target: Group | None = None
    if canvas_config.legend_position != "none":
        legend_measurement = legend_config.measure_legend(
            legend_table,
            placement=canvas_config.legend_position,
            wrap_width=source_primary_bounds.width,
        )
        if (
            legend_measurement.local_bounds.width > 0.0
            and legend_measurement.local_bounds.height > 0.0
        ):
            legend_target = LegendGroup(
                canvas_config,
                legend_measurement,
                legend_table,
                cfg=cfg,
            ).get_group()
            canvas.add(legend_target)

    title_target = (
        plot_title_obj.get_group()
        if plot_title_obj is not None
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
    composition_plan = plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", source_primary_bounds),
            legend=(
                CompositionItem("legend", legend_measurement.local_bounds)
                if legend_target is not None and legend_measurement is not None
                else None
            ),
            title=(
                CompositionItem("title", plot_title_obj.local_bounds)
                if title_target is not None and plot_title_obj is not None
                else None
            ),
            legend_placement=legend_placement,
            title_placement=title_placement,
        )
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
                record_collision_bands=record_collision_bands,
                boundary_gap_resolutions=boundary_gap_resolutions,
            ),
        )
    setattr(canvas, "_gbdraw_linear_source_content_bounds", source_primary_bounds)
    setattr(canvas, "_gbdraw_linear_composition_plan", composition_plan)

    return canvas


__all__ = ["assemble_linear_diagram"]
