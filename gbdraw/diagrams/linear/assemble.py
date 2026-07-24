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
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...exceptions import ValidationError
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from ...configurators.gc import _slot_skew_config
from ...core.text import calculate_bbox_dimensions
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...render.groups.linear import LengthBarGroup, LegendGroup, PlotTitleGroup  # type: ignore[reportMissingImports]
from ...render.groups.linear.length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
)
from ...io.comparisons import filter_comparison_dataframe, load_comparisons
from ...legend.table import (  # type: ignore[reportMissingImports]
    _unique_legend_key,
    configure_pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
)
from ...render.export import save_figure  # type: ignore[reportMissingImports]
from ...layout.linear import (  # type: ignore[reportMissingImports]
    AxisGapResolution,
    CollisionBand,
    LinearFeatureLaneGeometry,
    VerticalBand,
    measure_linear_feature_lanes,
    measure_linear_label_band,
    resolve_axis_gap,
)
from ...layout.linear_multi_record import (
    LinearRecordMeasurement,
    LinearRecordPlacement,
    RecordKey,
    resolve_record_row_positions,
    solve_linear_layout,
    stable_record_keys,
)
from ...linear_comparison import (
    LinearComparison,
    merge_linear_comparisons,
    validate_linear_comparison_topology,
)
from ...layout.scalar_axis import linear_scalar_axis_tick_font_size_px  # type: ignore[reportMissingImports]
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
    feature_underlay_anchor_slot_id,
    feature_underlay_slot_id,
    layout_annotation_track,
    merge_feature_underlays,
    resolve_annotations,
    sync_annotation_legend_entries,
)
from ...render.drawers.linear.annotations import draw_linear_annotation_track

from .builders import (
    add_explicit_comparisons_on_linear_canvas,
    add_depth_group,
    add_gc_content_group,
    add_gc_skew_group,
    add_legends_on_linear_canvas,
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
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.ids import make_linear_dom_id
from .track_slots import (
    LinearRecordVerticalPlan,
    LinearResolvedTrack,
    LinearSlotFootprint,
    LinearTrackLayout,
    resolve_linear_record_vertical_plan,
    resolve_linear_track_layout,
)


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from ...api.options import LinearMultiRecordOptions


def _annotation_marks_for_set(
    bundle: ResolvedAnnotationBundle,
    set_id: str,
) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(item.mark for item in bundle.annotations if item.set_id == set_id)
    )


def _prepare_linear_annotation_tracks(
    records: list[SeqRecord],
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    slots: list[LinearTrackSlot] | None,
    *,
    canvas_config: LinearCanvasConfigurator,
    record_depth_tracks: list[list[DepthTrackSpec]] | None,
) -> tuple[
    list[LinearTrackSlot] | None,
    ResolvedAnnotationBundle,
    dict[str, ResolvedAnnotationTrack],
    frozenset[str],
]:
    bundle = (
        annotations
        if isinstance(annotations, ResolvedAnnotationBundle)
        else resolve_annotations(annotations, records, mode="linear")
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, {}, frozenset()
    set_ids = bundle.set_ids or tuple(dict.fromkeys(item.set_id for item in bundle.annotations))
    if slots is None:
        slots = default_linear_track_slots(
            show_features=True,
            show_depth=bool(record_depth_tracks),
            depth_track_count=max(1, depth_track_count(record_depth_tracks)),
            show_gc=bool(canvas_config.show_gc),
            show_skew=bool(canvas_config.show_skew),
            track_layout=str(canvas_config.track_layout),
        )
        auto_annotation_slots = []
        for index, set_id in enumerate(set_ids):
            marks = _annotation_marks_for_set(bundle, set_id)
            lane_marks = tuple(mark for mark in marks if mark != "highlight")
            if lane_marks:
                auto_annotation_slots.append(
                    LinearTrackSlot(
                        id=f"annotations_{index + 1}",
                        renderer="annotations",
                        side="above",
                        params={"set_id": set_id, "marks": lane_marks},
                    )
                )
            if "highlight" in marks:
                highlight_slot_id = (
                    f"annotations_{index + 1}_highlight"
                    if lane_marks
                    else f"annotations_{index + 1}"
                )
                auto_annotation_slots.append(
                    LinearTrackSlot(
                        id=highlight_slot_id,
                        renderer="annotations",
                        side="overlay",
                        z=-1,
                        params={
                            "set_id": set_id,
                            "marks": ("highlight",),
                            "anchor_slot": "features",
                            "layer": "underlay",
                            "cover_anchor": True,
                            "padding_px": 0.0,
                        },
                    )
                )
        slots = auto_annotation_slots + slots

    requested_set_ids = {
        str(slot.params.get("set_id", "")).strip()
        for slot in slots
        if str(slot.renderer).strip().lower() == "annotations"
    }
    unknown = requested_set_ids - set(set_ids)
    if unknown:
        raise ValueError(f"Annotation track references unknown set_id(s): {', '.join(sorted(unknown))}")
    for set_id in set_ids:
        if set_id not in requested_set_ids:
            logger.warning("Annotation set %s is not referenced by a linear track slot.", set_id)

    auto_slot_ids = frozenset(
        slot.id
        for slot in slots
        if str(slot.renderer).strip().lower() == "annotations" and slot.height is None
    )
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
]:
    merged, set_id = merge_feature_underlays(
        bundle,
        [result.underlay_features for result in feature_layers],
        records,
        mode="linear",
    )
    if set_id is None:
        return slots, bundle, layouts, auto_slot_ids
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
    updated_auto_ids = frozenset((*auto_slot_ids, slot_id))
    updated_slots, updated_layouts = _layout_linear_annotation_tracks(
        records,
        updated_slots,
        merged,
        canvas_config=canvas_config,
        auto_slot_ids=updated_auto_ids,
    )
    return updated_slots, merged, updated_layouts, updated_auto_ids


def _is_axis_ruler_enabled(canvas_config: LinearCanvasConfigurator, cfg: GbdrawConfig) -> bool:
    track_layout = str(canvas_config.track_layout).strip().lower()
    scale_style = str(cfg.objects.scale.style).strip().lower()
    return (
        bool(canvas_config.ruler_on_axis)
        and scale_style == "ruler"
        and track_layout in {"above", "below"}
    )


def _feature_track_layout_for_linear_slots(
    slots: list,
    fallback: str,
) -> str:
    for slot in slots:
        if slot.renderer != "features":
            continue
        if slot.side == "above":
            return "above"
        if slot.side == "below":
            return "below"
        return "middle"
    return str(fallback or "middle").strip().lower()


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


def _axis_ruler_extents(canvas_config: LinearCanvasConfigurator, cfg: GbdrawConfig) -> tuple[float, float]:
    """
    Return (height_above_axis, height_below_axis) required by axis-based ruler labels.
    """
    if not _is_axis_ruler_enabled(canvas_config, cfg):
        return 0.0, 0.0

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
    track_layout = str(canvas_config.track_layout).strip().lower()
    if track_layout == "above":
        return 0.0, protrusion
    if track_layout == "below":
        return protrusion, 0.0
    return 0.0, 0.0


def _linear_axis_band(
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> VerticalBand:
    axis_stroke_width = cfg.objects.axis.linear.stroke_width.for_length_param(
        canvas_config.length_param
    )
    half_stroke = 0.5 * max(0.0, float(axis_stroke_width))
    ruler_above, ruler_below = _axis_ruler_extents(canvas_config, cfg)
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


def _layout_with_preferred_origins(
    layout: LinearTrackLayout,
    preferred_origins: dict[str, float],
) -> LinearTrackLayout:
    if not preferred_origins:
        return layout
    return replace(
        layout,
        slots=tuple(
            replace(slot, y_offset=float(preferred_origins.get(slot.id, slot.y_offset)))
            for slot in layout.slots
        ),
    )


def _default_preferred_origins(
    layout: LinearTrackLayout,
    canvas_config: LinearCanvasConfigurator,
) -> dict[str, float]:
    base = float(canvas_config.cds_padding) + float(canvas_config.vertical_padding)
    axis_clearance = float(canvas_config.vertical_padding)
    return {
        slot.id: base + float(slot.y_offset) - axis_clearance
        for slot in layout.slots
        if slot.renderer in {"depth", "dinucleotide_content", "dinucleotide_skew"}
    }


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
    config_dict: dict,
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
            config_dict,
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
) -> tuple[list[LinearRecordVerticalPlan], list[_LinearRecordDefinitionGeometry]]:
    axis_band = _linear_axis_band(canvas_config, cfg)
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
    config_dict: dict,
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
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """
    Assembles a linear diagram of genomic records with optional BLAST comparison data,
    and returns the SVG canvas (not saved).
    """
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
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
    custom_linear_track_slots_requested = linear_track_slots is not None
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
        record_depth_tracks=record_depth_tracks,
    )
    if linear_track_slots is None:
        linear_track_slots = default_linear_track_slots(
            show_features=True,
            show_depth=bool(canvas_config.show_depth and record_depth_tracks),
            depth_track_count=max(1, depth_track_count(record_depth_tracks)),
            show_gc=bool(canvas_config.show_gc),
            show_skew=bool(canvas_config.show_skew),
            dinucleotide=str(gc_config.dinucleotide),
            track_layout=str(canvas_config.track_layout),
        )
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        record_depth_tracks,
    )
    canvas_config.track_layout = _feature_track_layout_for_linear_slots(
        normalized_linear_track_slots,
        canvas_config.track_layout,
    )
    axis_ruler_enabled = _is_axis_ruler_enabled(canvas_config, cfg)
    normalized_plot_title = str(plot_title or "").strip()
    normalized_plot_title_position = str(plot_title_position or "bottom").strip().lower()
    if normalized_plot_title_position not in {"center", "top", "bottom"}:
        raise ValueError("plot_title_position must be one of: center, top, bottom")

    plot_title_obj: PlotTitleGroup | None = None
    plot_title_edge_margin = 24.0
    plot_title_vertical_gap = float(canvas_config.vertical_padding)
    plot_title_top_reserve = 0.0
    plot_title_bottom_reserve = 0.0
    if normalized_plot_title:
        plot_title_obj = PlotTitleGroup(
            normalized_plot_title,
            config_dict,
            font_size=float(plot_title_font_size),
            cfg=cfg,
        )
        plot_title_text_height = max(float(plot_title_obj.text_bbox_height), float(plot_title_font_size))
        reserve = plot_title_edge_margin + plot_title_text_height + plot_title_vertical_gap
        if normalized_plot_title_position == "top":
            plot_title_top_reserve = reserve
        elif normalized_plot_title_position == "bottom":
            plot_title_bottom_reserve = reserve

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

    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")
    orthogroup_label_eligibility = None
    if show_labels_mode == "orthogroup_top":
        orthogroup_label_eligibility = build_orthogroup_label_eligibility(
            orthogroups=orthogroups,
            comparisons=comparisons,
        )
    if (
        show_labels_mode == "orthogroup_top"
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
        config_dict,
        cfg=cfg,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
    )
    record_feature_dicts = [
        result.foreground_features for result in record_feature_layers
    ]
    (
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
    ) = _add_linear_feature_underlays(
        records,
        record_feature_layers,
        linear_track_slots,
        resolved_annotations,
        annotation_track_layouts,
        auto_annotation_slot_ids,
        canvas_config=canvas_config,
    )
    normalized_linear_track_slots = normalize_linear_track_slots_with_axis(
        linear_track_slots,
        linear_track_axis_index,
    )
    normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
        normalized_linear_track_slots,
        record_depth_tracks,
    )
    record_feature_lane_geometries = [
        measure_linear_feature_lanes(
            feature_dict,
            cds_height=float(canvas_config.cds_height),
            separate_strands=bool(canvas_config.strandedness),
            track_layout=str(canvas_config.track_layout),
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
        config_dict,
        cfg=cfg,
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
        config_dict,
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
            config_dict,
            canvas_config,
            cfg=cfg,
            line_kinds_by_record=[
                frozenset({"name", "subtitle"})
                if index in row_leading_indices
                else frozenset()
                for index in range(len(records))
            ],
        )

    linear_track_layout = resolve_linear_track_layout(
        normalized_linear_track_slots,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    if not custom_linear_track_slots_requested:
        linear_track_layout = _layout_with_preferred_origins(
            linear_track_layout,
            _default_preferred_origins(linear_track_layout, canvas_config),
        )
    resolved_depth_heights_by_index = _resolved_linear_depth_heights(
        linear_track_layout
    )
    primary_depth_track_index = min(resolved_depth_heights_by_index, default=0)
    depth_enabled = bool(canvas_config.show_depth and depth_config is not None and record_depth_tracks)
    record_depth_data: list[list[DepthTrackData]] = (
        build_depth_track_dataframes(
            records,
            record_depth_tracks,
            base_config=depth_config,
            depth_df_builder=build_depth_df,
        )
        if depth_enabled
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
    if not depth_enabled:
        canvas_config.show_depth = False

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
    legend_group: LegendGroup | None = None
    required_legend_height = 0.0
    if canvas_config.legend_position != "none":
        color_map, default_color_map = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
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
            depth_config=depth_config if depth_enabled and depth_track_count(record_depth_tracks) == 1 else None,
        )
        legend_table = _sync_legend_table_for_linear_slots(
            legend_table,
            linear_track_slots=normalized_linear_track_slots,
            skew_config=skew_config,
            depth_tracks=representative_depth_data if depth_enabled else None,
        )
        legend_table = sync_annotation_legend_entries(
            legend_table,
            resolved_annotations,
            normalized_linear_track_slots,
        )
        legend_config = legend_config.recalculate_legend_dimensions(legend_table, canvas_config)
        legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table, cfg=cfg)
        required_legend_height = float(legend_group.legend_height)
        definition_reserve_width = (
            row_definition_width
            if split_row_definitions
            else (0.0 if multi_record_enabled else max_def_width)
        )
        canvas_config.recalculate_canvas_dimensions(
            legend_group,
            definition_reserve_width,
        )
    else:
        definition_reserve_width = (
            row_definition_width
            if split_row_definitions
            else (0.0 if multi_record_enabled else max_def_width)
        )
        canvas_config.alignment_width = canvas_config.fig_width
        canvas_config.horizontal_offset = (
            2 * canvas_config.canvas_padding
            + definition_reserve_width
            + canvas_config.definition_gap
        )
        canvas_config.total_width = (
            canvas_config.horizontal_offset
            + canvas_config.alignment_width
            + 2 * canvas_config.canvas_padding
        )
        canvas_config.legend_offset_x = 0
        canvas_config.legend_offset_y = 0
    ordinary_row_gap = float(canvas_config.cds_padding) * 1.5
    definition_clear_gap = max(1.0, 0.5 * float(canvas_config.vertical_padding))
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
        if canvas_config.legend_position == "top":
            first_axis_y = (
                canvas_config.original_vertical_offset
                + required_legend_height
                + canvas_config.vertical_offset
            )
        else:
            first_axis_y = max(
                canvas_config.vertical_offset,
                canvas_config.original_vertical_offset,
            )
        return first_axis_y + plot_title_top_reserve

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
                linear_track_axis_index,
            )
            normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
                normalized_linear_track_slots,
                record_depth_tracks,
            )
            canvas_config.track_layout = _feature_track_layout_for_linear_slots(
                normalized_linear_track_slots,
                canvas_config.track_layout,
            )
            linear_track_layout = resolve_linear_track_layout(
                normalized_linear_track_slots,
                canvas_config=canvas_config,
                cfg=cfg,
            )
            if not custom_linear_track_slots_requested:
                linear_track_layout = _layout_with_preferred_origins(
                    linear_track_layout,
                    _default_preferred_origins(linear_track_layout, canvas_config),
                )
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
            config_dict,
            cfg=cfg,
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
                comparison_height=float(
                    canvas_config.configured_comparison_height
                ),
                definition_clear_gap=definition_clear_gap,
                boundary_has_comparison=bool(active_records),
            )
            boundary_gap_resolutions.append(resolution)
            current_y += resolution.axis_gap
        current_y = max(record_offsets, default=current_y)

    length_bar_group: LengthBarGroup | None = None
    length_bar_offset_x = 0.0
    if not multi_record_enabled and not canvas_config.normalize_length and not axis_ruler_enabled:
        length_bar_offset_x = alignment_extents.ruler_offset_x
        length_bar_group = LengthBarGroup(
            canvas_config.fig_width,
            canvas_config.alignment_width,
            canvas_config.longest_genome,
            config_dict,
            canvas_config,
            cfg=cfg,
            ruler_width=alignment_extents.ruler_width,
        )
    length_bar_height = (
        float(length_bar_group.scale_group_height)
        if length_bar_group is not None
        else 0.0
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
    bottom_title_stack = (
        plot_title_obj is not None
        and normalized_plot_title_position == "bottom"
        and canvas_config.legend_position == "bottom"
        and legend_group is not None
    )
    bottom_stack_gap = float(canvas_config.vertical_padding)
    if bottom_title_stack:
        final_height = (
            canvas_config.height_below_final_record
            + length_bar_height
            + bottom_stack_gap
            + required_legend_height
            + bottom_stack_gap
            + float(plot_title_obj.text_bbox_height)
            + plot_title_edge_margin
        )
    else:
        final_height = (
            canvas_config.height_below_final_record
            + length_bar_height
            + canvas_config.original_vertical_offset
            + plot_title_bottom_reserve
        )
        if canvas_config.legend_position in ["top", "bottom"]:
            final_height += int(required_legend_height)
    canvas_config.total_height = max(final_height, required_legend_height)

    if legend_group is not None:
        canvas_config.recalculate_canvas_dimensions(
            legend_group, definition_reserve_width
        )

    alignment_shift_x = alignment_extents.horizontal_shift
    alignment_width_extension = alignment_extents.width_extension
    if alignment_shift_x or alignment_width_extension:
        width_extension_px = math.ceil(alignment_width_extension)
        canvas_config.horizontal_offset += alignment_shift_x
        canvas_config.total_width += width_extension_px
        if canvas_config.legend_position == "right":
            canvas_config.legend_offset_x += width_extension_px
        elif canvas_config.legend_position in {"top", "bottom"}:
            canvas_config.legend_offset_x += 0.5 * width_extension_px
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
            comparison_top_y=axis_y - record_plan.comparison_top_extent,
            comparison_bottom_y=axis_y + record_plan.comparison_bottom_extent,
            px_per_bp=sequence_width / max(1, len(record.seq)),
        )
    canvas: Drawing = canvas_config.create_svg_canvas()

    # Embed both viewBox configurations as data attributes for JavaScript repositioning
    # This allows switching between horizontal and vertical legend layouts without accumulation errors
    vertical_vb_width = canvas_config.total_width
    vertical_vb_height = canvas_config.total_height
    horizontal_vb_width = canvas_config.total_width
    horizontal_vb_height = canvas_config.total_height
    if legend_group is not None:
        h_legend_width, h_legend_height = legend_group.get_horizontal_dimensions()
        v_legend_width, v_legend_height = legend_group.get_vertical_dimensions()

        if canvas_config.legend_position in ["top", "bottom"]:
            vertical_vb_width = canvas_config.total_width - h_legend_width + v_legend_width
            vertical_vb_height = canvas_config.total_height - h_legend_height

        if canvas_config.legend_position in ["left", "right"]:
            horizontal_vb_width = canvas_config.total_width - v_legend_width
            horizontal_vb_height = canvas_config.total_height + h_legend_height

    canvas.attribs["data-vertical-viewbox"] = f"0 0 {vertical_vb_width} {vertical_vb_height}"
    canvas.attribs["data-horizontal-viewbox"] = f"0 0 {horizontal_vb_width} {horizontal_vb_height}"

    if canvas_config.legend_position == "top" and plot_title_top_reserve > 0:
        canvas_config.legend_offset_y += plot_title_top_reserve
    if bottom_title_stack:
        canvas_config.legend_offset_y = (
            float(canvas_config.total_height)
            - plot_title_edge_margin
            - float(plot_title_obj.text_bbox_height)
            - bottom_stack_gap
            - required_legend_height
        )

    if canvas_config.legend_position != "none":
        canvas = add_legends_on_linear_canvas(canvas, config_dict, canvas_config, legend_group, legend_table)
    if length_bar_group is not None:
        canvas = add_length_bar_on_linear_canvas(
            canvas,
            canvas_config,
            config_dict,
            length_bar_group,
            legend_group,
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
        )

    label_font_size = _resolve_linear_diagram_label_font_size(
        records,
        show_labels_mode=show_labels_mode,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    cfg_labels_on = replace(cfg, canvas=replace(cfg.canvas, show_labels=True))
    cfg_labels_off = replace(cfg, canvas=replace(cfg.canvas, show_labels=False))

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
        should_show_labels = False
        if show_labels_mode == "all":
            should_show_labels = True
        elif show_labels_mode == "first" and count == 1:
            should_show_labels = True
        elif show_labels_mode == "orthogroup_top":
            should_show_labels = True
        record_cfg = cfg_labels_on if should_show_labels else cfg_labels_off
        orthogroup_label_member_ids, orthogroup_label_top_member_ids = orthogroup_label_sets_for_record(
            orthogroup_label_eligibility if show_labels_mode == "orthogroup_top" else None,
            record_index,
        )
        record_group_id = make_linear_dom_id(
            record.id,
            record_index=record_index,
            record_count=total_records,
        )
        definition_group_id = make_linear_dom_id(
            record.id,
            record_index=record_index,
            record_count=total_records,
            suffix="definition",
        )

        if linear_track_layout is not None:
            feature_rendered = False

            for slot in sorted(linear_track_layout.slots, key=lambda item: (item.z, item.slot_index)):
                resolved_slot = record_vertical_plan.slot_by_id(slot.id)
                if slot.renderer == "spacer":
                    continue
                if slot.renderer == "features":
                    slot_feature_layout = "middle" if slot.side == "overlay" else slot.side
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
                        precomputed_feature_dict=record_feature_dicts[record_index],
                        feature_track_layout=slot_feature_layout,
                        label_font_size=label_font_size,
                        orthogroup_label_member_ids=orthogroup_label_member_ids,
                        orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
                        record_index=record_index,
                        record_count=total_records,
                        group_id=record_group_id,
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
                    add_depth_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        depth_track.config,
                        config_dict,
                        cfg=record_cfg,
                        depth_df=depth_track.df,
                        group_id=group_id,
                        axis_group_id=f"{group_id}_axis",
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
                        depth_enabled=depth_enabled,
                    )
                    add_gc_content_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        slot_gc_config,
                        config_dict,
                        cfg=record_cfg,
                        gc_df=shared_gc_df,
                        track_height=slot.height,
                        track_offset_y=track_offset_y,
                        group_id=make_linear_dom_id(
                            slot.id,
                            record_index=record_index,
                            record_count=total_records,
                        ),
                        sequence_width=sequence_width,
                    )
                    continue

                if slot.renderer == "dinucleotide_skew":
                    nt = _slot_nt(slot, str(skew_config.dinucleotide))
                    per_nt_gc_dfs = record_gc_dfs_by_nt.get(nt, record_gc_dfs)
                    shared_gc_df = per_nt_gc_dfs[record_index] if record_index < len(per_nt_gc_dfs) else None
                    add_gc_skew_group(
                        canvas,
                        record,
                        offset_y,
                        offset_x,
                        canvas_config,
                        _slot_skew_config(skew_config, slot, nt),
                        config_dict,
                        cfg=record_cfg,
                        gc_df=shared_gc_df,
                        track_height=slot.height,
                        track_offset_y=track_offset_y,
                        group_id=make_linear_dom_id(
                            slot.id,
                            record_index=record_index,
                            record_count=total_records,
                        ),
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
                    config_dict,
                    precalculated_labels=None,
                    cfg=record_cfg,
                    precomputed_feature_dict=record_feature_dicts[record_index],
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
                config_dict,
                definition_column_width,
                cfg=record_cfg,
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
            )
            continue

    if plot_title_obj is not None:
        title_group = plot_title_obj.get_group()
        title_height = float(plot_title_obj.text_bbox_height)
        title_y = 0.5 * float(canvas_config.total_height)
        if normalized_plot_title_position == "top":
            title_y = plot_title_edge_margin + (0.5 * title_height)
        elif normalized_plot_title_position == "bottom":
            if bottom_title_stack:
                title_y = (
                    canvas_config.legend_offset_y
                    + required_legend_height
                    + bottom_stack_gap
                    + (0.5 * title_height)
                )
            else:
                title_y = float(canvas_config.total_height) - plot_title_edge_margin - (0.5 * title_height)
        title_group.translate(0.5 * float(canvas_config.total_width), title_y)
        canvas.add(title_group)

    if linear_track_layout is not None:
        setattr(
            canvas,
            "_gbdraw_track_slot_geometry",
            _serialize_linear_track_slot_geometry(
                records=records,
                layout=linear_track_layout,
                record_plans=record_vertical_plans,
                record_offsets=record_offsets,
                record_collision_bands=record_collision_bands,
                boundary_gap_resolutions=boundary_gap_resolutions,
            ),
        )

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
    depth_config: DepthConfigurator | None = None,
    depth_tables: list[DataFrame | None] | None = None,
    cfg: GbdrawConfig | None = None,
    comparison_dataframes: list[DataFrame] | None = None,
    orthogroups: OrthogroupResult | None = None,
    align_orthogroup_feature: str | None = None,
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
        depth_config=depth_config,
        depth_tables=depth_tables,
        comparison_dataframes=comparison_dataframes,
        orthogroups=orthogroups,
        align_orthogroup_feature=align_orthogroup_feature,
        cfg=cfg,
    )
    save_figure(canvas, out_formats)
    return canvas


__all__ = [
    "assemble_linear_diagram",
    "plot_linear_diagram",
]
