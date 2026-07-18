#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly (implementation).

This module was extracted from `gbdraw.linear_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import copy
from dataclasses import replace
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
    normalize_depth_tracks,
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
    calculate_feature_position_factors_linear,
    resolve_feature_axis_gap_linear,
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
from ...labels.linear import calculate_label_y_bounds  # type: ignore[reportMissingImports]
from ...tracks import (
    LinearTrackSlot,
    ScalarSpec,
    default_linear_track_slots,
    normalize_linear_track_slots_with_axis,
)
from ...annotations import (
    AnnotationOptions,
    ResolvedAnnotationBundle,
    ResolvedAnnotationTrack,
    annotation_track_params_from_mapping,
    layout_annotation_track,
    resolve_annotations,
    sync_annotation_legend_entries,
)
from ...render.drawers.linear.annotations import draw_linear_annotation_track

from .builders import (
    add_comparison_on_linear_canvas,
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
    _precalculate_feature_dicts,
    _precalculate_label_dimensions,
    _resolve_linear_diagram_label_font_size,
)
from ...features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.ids import make_linear_dom_id
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from .track_slots import LinearResolvedTrack, LinearTrackLayout, resolve_linear_track_layout


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from ...api.options import LinearMultiRecordOptions


def _prepare_linear_annotation_tracks(
    records: list[SeqRecord],
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    slots: list[LinearTrackSlot] | None,
    *,
    canvas_config: LinearCanvasConfigurator,
    record_depth_tracks: list[list[DepthTrackSpec]] | None,
) -> tuple[list[LinearTrackSlot] | None, ResolvedAnnotationBundle, dict[str, ResolvedAnnotationTrack]]:
    bundle = (
        annotations
        if isinstance(annotations, ResolvedAnnotationBundle)
        else resolve_annotations(annotations, records, mode="linear")
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, {}
    set_ids = bundle.set_ids or tuple(dict.fromkeys(item.set_id for item in bundle.annotations))
    if slots is None:
        slots = default_linear_track_slots(
            show_features=True,
            show_depth=bool(record_depth_tracks),
            depth_track_count=max((len(items) for items in (record_depth_tracks or ())), default=1),
            show_gc=bool(canvas_config.show_gc),
            show_skew=bool(canvas_config.show_skew),
            track_layout=str(canvas_config.track_layout),
        )
        slots = [
            LinearTrackSlot(
                id=f"annotations_{index + 1}",
                renderer="annotations",
                side="above",
                params={"set_id": set_id},
            )
            for index, set_id in enumerate(set_ids)
        ] + slots

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

    record_lengths = {index: len(record.seq) for index, record in enumerate(records)}
    bp_per_px = {
        index: len(record.seq)
        / max(
            1.0,
            float(canvas_config.alignment_width)
            * (1.0 if canvas_config.normalize_length else len(record.seq) / max(1, canvas_config.longest_genome)),
        )
        for index, record in enumerate(records)
    }
    updated_slots: list[LinearTrackSlot] = []
    layouts: dict[str, ResolvedAnnotationTrack] = {}
    for slot in slots:
        if str(slot.renderer).strip().lower() != "annotations":
            updated_slots.append(slot)
            continue
        params = annotation_track_params_from_mapping(slot.params)
        available = slot.height.resolve(1.0) if slot.height is not None else None
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
            slot
            if slot.height is not None
            else replace(slot, height=ScalarSpec(layout.required_extent_px, "px"))
        )
    return updated_slots, bundle, layouts


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


def _feature_slot_for_linear_slots(slots: list) -> object | None:
    return next((slot for slot in slots if slot.renderer == "features"), None)


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
        try:
            track_index = int(slot.params.get("track_index", 0))
        except (AttributeError, TypeError, ValueError):
            track_index = 0
        if 0 <= track_index < len(depth_heights) and depth_heights[track_index] is not None:
            out.append(replace(slot, height=ScalarSpec(float(depth_heights[track_index]), "px")))
        else:
            out.append(slot)
    return out


def _slot_nt(slot: LinearResolvedTrack, default_nt: str) -> str:
    params = slot.params or {}
    return str(params.get("nt", params.get("dinucleotide", default_nt)) or default_nt).upper()


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
) -> dict:
    """Replace singleton skew legend entries with slot-aware entries."""
    if linear_track_slots is None:
        return legend_table

    out = dict(legend_table)
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


def _custom_record_bottom_extent(
    *,
    record_index: int,
    layout: LinearTrackLayout,
    feature_slot,
    canvas_config: LinearCanvasConfigurator,
    record_heights_below: list[float],
    record_label_heights_below: list[float],
) -> float:
    extent = float(layout.bottom_extent)
    feature_side = getattr(feature_slot, "side", None) if feature_slot is not None else None
    feature_extent = (
        record_heights_below[record_index]
        + record_label_heights_below[record_index]
    )
    if feature_side == "below":
        return extent + feature_extent + canvas_config.vertical_padding
    if feature_side is None or feature_side == "overlay":
        return max(extent, feature_extent)
    return extent


def _custom_record_top_extent(
    *,
    record_index: int,
    layout: LinearTrackLayout,
    feature_slot,
    canvas_config: LinearCanvasConfigurator,
    record_heights_above: list[float],
    record_label_heights_above: list[float],
) -> float:
    extent = float(layout.top_extent)
    feature_side = getattr(feature_slot, "side", None) if feature_slot is not None else None
    feature_extent = max(
        record_heights_above[record_index],
        record_label_heights_above[record_index],
    )
    if feature_side == "above":
        return extent + feature_extent + canvas_config.vertical_padding
    if feature_side is None or feature_side == "overlay":
        return max(extent, feature_extent)
    return extent


def _custom_slot_track_offset_y(
    slot: LinearResolvedTrack,
    *,
    feature_slot,
    record_index: int,
    canvas_config: LinearCanvasConfigurator,
    record_heights_below: list[float],
    record_heights_above: list[float],
    record_label_heights_below: list[float],
    record_label_heights_above: list[float],
) -> float:
    offset = float(slot.y_offset)
    if feature_slot is None or slot.renderer in {"features", "spacer"}:
        return offset
    feature_side = getattr(feature_slot, "side", None)
    feature_index = int(getattr(feature_slot, "slot_index", -1))
    if slot.side == "below" and feature_side == "below" and feature_index < int(slot.slot_index):
        return (
            offset
            + record_heights_below[record_index]
            + record_label_heights_below[record_index]
            + canvas_config.vertical_padding
        )
    if slot.side == "above" and feature_side == "above" and feature_index > int(slot.slot_index):
        return (
            offset
            - record_heights_above[record_index]
            - record_label_heights_above[record_index]
            - canvas_config.vertical_padding
        )
    return offset


def _serialize_linear_track_slot_geometry(
    *,
    records: list[SeqRecord],
    layout: LinearTrackLayout,
    record_offsets: list[float],
    feature_slot,
    canvas_config: LinearCanvasConfigurator,
    record_heights_below: list[float],
    record_heights_above: list[float],
    record_label_heights_below: list[float],
    record_label_heights_above: list[float],
) -> dict[str, Any]:
    records_payload: list[dict[str, Any]] = []
    for record_index, record in enumerate(records):
        axis_y = float(record_offsets[record_index]) if record_index < len(record_offsets) else 0.0
        slots_payload: list[dict[str, Any]] = []
        for slot in layout.slots:
            final_offset = _custom_slot_track_offset_y(
                slot,
                feature_slot=feature_slot,
                record_index=record_index,
                canvas_config=canvas_config,
                record_heights_below=record_heights_below,
                record_heights_above=record_heights_above,
                record_label_heights_below=record_label_heights_below,
                record_label_heights_above=record_label_heights_above,
            )
            slots_payload.append(
                {
                    "slotIndex": int(slot.slot_index),
                    "slotId": str(slot.id),
                    "renderer": str(slot.renderer),
                    "side": str(slot.side),
                    "heightPx": float(slot.height),
                    "spacingAfterPx": float(slot.spacing_after_px),
                    "baseYOffsetPx": float(slot.y_offset),
                    "finalYOffsetPx": axis_y + float(final_offset),
                    "source": "resolved",
                }
            )
        record_id = str(getattr(record, "id", "") or "")
        records_payload.append(
            {
                "recordIndex": int(record_index),
                "recordId": record_id,
                "recordLabel": record_id,
                "slots": slots_payload,
            }
        )
    return {
        "schema": 1,
        "mode": "linear",
        "source": "resolved",
        "records": records_payload,
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
    canvas_config: LinearCanvasConfigurator,
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
    cloned.tick_font_size = linear_scalar_axis_tick_font_size_px(depth_config, canvas_config.depth_height)
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


def _record_plot_track_stack_bottom_y(
    *,
    axis_y: float,
    record_index: int,
    canvas_config: LinearCanvasConfigurator,
    record_heights_below: list[float],
    record_label_heights_below: list[float],
    non_middle_layout: bool,
) -> float | None:
    """Return the absolute y-coordinate for the drawn plot-track stack bottom."""
    visual_bottom = float(getattr(canvas_config, "plot_tracks_visual_bottom", 0.0))
    if visual_bottom <= 0.0:
        return None

    plot_offset_y = float(axis_y)
    if non_middle_layout:
        current_feature_height_below = record_heights_below[record_index]
        plot_offset_y += current_feature_height_below - canvas_config.cds_padding
    plot_offset_y += record_label_heights_below[record_index]
    return (
        plot_offset_y
        + canvas_config.cds_padding
        + canvas_config.vertical_padding
        + visual_bottom
    )


def _precalculate_feature_track_heights(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
    precomputed_feature_dicts: list[FeatureDict] | None = None,
) -> tuple[list[float], list[float], list[float], list[float], list[float], list[float]]:
    """
    Pre-calculates the height required for feature tracks for each record.
    This is needed when resolve_overlaps is enabled as features may span multiple tracks.

    Returns:
        - list mapping record index -> height below the axis line (for lower tracks)
        - list mapping record index -> height above the axis line (for upper tracks)
        - list mapping record index -> minimum above-axis extent required to keep top visible
        - list mapping record index -> above-axis extent with non-displaced track positioning
        - list mapping record index -> middle-layout above-axis extent with non-displaced tracks
        - list mapping record index -> feature-only visual center relative to the axis
    """
    record_heights_below: list[float] = []
    record_heights_above: list[float] = []
    record_top_guard_above: list[float] = []
    record_top_guard_undisplaced: list[float] = []
    record_top_guard_middle_undisplaced: list[float] = []
    record_feature_center_offsets: list[float] = []
    track_layout = str(canvas_config.track_layout).strip().lower()
    axis_gap_factor = (
        (float(canvas_config.track_axis_gap) / float(canvas_config.cds_height))
        if (canvas_config.track_axis_gap is not None and float(canvas_config.cds_height) > 0.0)
        else None
    )
    axis_ruler_above, axis_ruler_below = _axis_ruler_extents(canvas_config, cfg)

    if precomputed_feature_dicts is None:
        color_table, default_colors = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )

    for index, record in enumerate(records):
        if precomputed_feature_dicts is not None:
            feature_dict = precomputed_feature_dicts[index]
        else:
            feature_dict, _ = create_feature_dict(
                record,
                color_table,
                feature_config.selected_features_set,
                default_colors,
                canvas_config.strandedness,
                canvas_config.resolve_overlaps,
                {},
                directional_feature_types=feature_config.directional_feature_types,
                feature_visibility_rules=feature_config.feature_visibility_rules,
                compute_label_text=False,
            )
        min_top_y = 0.0
        max_bottom_y = 0.0
        feature_top_y: float | None = None
        feature_bottom_y: float | None = None
        min_top_y_undisplaced = 0.0
        min_top_y_middle_undisplaced = 0.0
        for feature_obj in feature_dict.values():
            track_id = int(getattr(feature_obj, "feature_track_id", 0))
            strand = str(getattr(feature_obj, "strand", "undefined"))
            factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=track_id,
                separate_strands=canvas_config.strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            top_y = canvas_config.cds_height * float(factors[0])
            bottom_y = canvas_config.cds_height * float(factors[2])
            feature_top_y = top_y if feature_top_y is None else min(feature_top_y, top_y)
            feature_bottom_y = (
                bottom_y if feature_bottom_y is None else max(feature_bottom_y, bottom_y)
            )
            if top_y < min_top_y:
                min_top_y = top_y
            if bottom_y > max_bottom_y:
                max_bottom_y = bottom_y

            undisplaced_track_id = -1 if (canvas_config.strandedness and strand == "negative") else 0
            undisplaced_factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=undisplaced_track_id,
                separate_strands=canvas_config.strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            undisplaced_top_y = canvas_config.cds_height * float(undisplaced_factors[0])
            if undisplaced_top_y < min_top_y_undisplaced:
                min_top_y_undisplaced = undisplaced_top_y
            middle_undisplaced_factors = calculate_feature_position_factors_linear(
                strand=strand,
                track_id=undisplaced_track_id,
                separate_strands=canvas_config.strandedness,
                track_layout="middle",
                axis_gap_factor=axis_gap_factor,
            )
            middle_undisplaced_top_y = canvas_config.cds_height * float(middle_undisplaced_factors[0])
            if middle_undisplaced_top_y < min_top_y_middle_undisplaced:
                min_top_y_middle_undisplaced = middle_undisplaced_top_y

        precise_height_above = max(0.0, -min_top_y)
        precise_height_below = max(0.0, max_bottom_y)
        undisplaced_height_above = max(0.0, -min_top_y_undisplaced)
        middle_undisplaced_height_above = max(0.0, -min_top_y_middle_undisplaced)

        if track_layout == "middle":
            # Keep existing middle-mode sizing behavior for backward compatibility.
            max_positive_track = 0
            min_negative_track = 0

            for feature_obj in feature_dict.values():
                track_id = feature_obj.feature_track_id
                if track_id > max_positive_track:
                    max_positive_track = track_id
                if track_id < min_negative_track:
                    min_negative_track = track_id

            if canvas_config.strandedness:
                # Stranded mode: positive tracks above axis, negative tracks below
                num_tracks_above = max_positive_track + 1
                num_tracks_below = abs(min_negative_track) + 1 if min_negative_track < 0 else 1
                height_above = num_tracks_above * canvas_config.cds_height * 1.1
                height_below = num_tracks_below * canvas_config.cds_height * 1.1
            else:
                # Non-stranded mode: track 0 is at axis, higher track IDs go below
                # Each track needs full cds_height of space plus some padding
                # Track 0 extends both above and below the axis (0.6 each direction)
                # Additional tracks (1, 2, ...) add cds_height * 1.1 below for breathing room
                height_above = canvas_config.cds_height * 0.6
                height_below = canvas_config.cds_height * (0.6 + max_positive_track * 1.1)
        else:
            # For above/below layouts, use actual positioned factors so downstream spacing,
            # GC/skew placement, and comparison ribbons align with feature extents.
            height_above = precise_height_above
            height_below = precise_height_below

        if axis_ruler_above > 0.0:
            height_above = max(height_above, axis_ruler_above)
            precise_height_above = max(precise_height_above, axis_ruler_above)
            undisplaced_height_above = max(undisplaced_height_above, axis_ruler_above)
            middle_undisplaced_height_above = max(middle_undisplaced_height_above, axis_ruler_above)
        if axis_ruler_below > 0.0:
            height_below = max(height_below, axis_ruler_below)

        record_heights_above.append(height_above)
        record_heights_below.append(height_below)
        record_top_guard_above.append(precise_height_above)
        record_top_guard_undisplaced.append(undisplaced_height_above)
        record_top_guard_middle_undisplaced.append(middle_undisplaced_height_above)
        record_feature_center_offsets.append(
            0.0
            if feature_top_y is None or feature_bottom_y is None
            else 0.5 * (feature_top_y + feature_bottom_y)
        )

    return (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_undisplaced,
        record_top_guard_middle_undisplaced,
        record_feature_center_offsets,
    )


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


def _precalculate_label_heights_below(all_labels_by_record: list[list[dict]]) -> list[float]:
    """Return per-record label extents that protrude below the axis."""
    record_label_heights_below: list[float] = []
    for labels in all_labels_by_record:
        max_bottom_y = 0.0
        for label in labels:
            _, label_bottom_y = calculate_label_y_bounds(label)
            if bool(label.get("is_embedded")):
                feature_bottom_y = float(label.get("feature_bottom_y", 0.0))
                # Ignore labels that do not protrude outside their feature track.
                if label_bottom_y <= feature_bottom_y:
                    continue
            if label_bottom_y > max_bottom_y:
                max_bottom_y = label_bottom_y
        record_label_heights_below.append(max_bottom_y)
    return record_label_heights_below


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
    linear_track_slots, resolved_annotations, annotation_track_layouts = _prepare_linear_annotation_tracks(
        records,
        annotations,
        linear_track_slots,
        canvas_config=canvas_config,
        record_depth_tracks=record_depth_tracks,
    )
    normalized_linear_track_slots = (
        normalize_linear_track_slots_with_axis(
            linear_track_slots,
            linear_track_axis_index,
        )
        if linear_track_slots is not None
        else None
    )
    if normalized_linear_track_slots is not None:
        normalized_linear_track_slots = _apply_depth_track_heights_to_linear_slots(
            normalized_linear_track_slots,
            record_depth_tracks,
        )
    if normalized_linear_track_slots is not None:
        canvas_config.track_layout = _feature_track_layout_for_linear_slots(
            normalized_linear_track_slots,
            canvas_config.track_layout,
        )
    track_layout = str(canvas_config.track_layout).strip().lower()
    non_middle_layout = track_layout in {"above", "below"}
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

    record_feature_dicts = _precalculate_feature_dicts(
        records,
        feature_config,
        canvas_config,
        config_dict,
        cfg=cfg,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
    )
    required_label_height, all_labels, record_label_heights_above = _precalculate_label_dimensions(
        records,
        feature_config,
        canvas_config,
        config_dict,
        cfg=cfg,
        precomputed_feature_dicts=record_feature_dicts,
        orthogroup_label_eligibility=orthogroup_label_eligibility,
    )
    record_label_heights_below = _precalculate_label_heights_below(all_labels)
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
    max_def_width, _definition_heights, definition_half_heights = _precalculate_definition_metrics(
        records,
        config_dict,
        canvas_config,
        cfg=cfg,
        line_kinds_by_record=local_definition_line_kinds,
    )
    row_definition_width = 0.0
    if split_row_definitions:
        row_definition_width, _row_definition_heights, _row_definition_half_heights = (
            _precalculate_definition_metrics(
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
        )

    # Pre-calculate feature track heights for each record (needed for resolve_overlaps)
    (
        record_heights_below,
        record_heights_above,
        record_top_guard_above,
        record_top_guard_undisplaced,
        record_top_guard_middle_undisplaced,
        record_feature_center_offsets,
    ) = _precalculate_feature_track_heights(
        records,
        feature_config,
        canvas_config,
        cfg,
        precomputed_feature_dicts=record_feature_dicts,
    )

    linear_track_layout: LinearTrackLayout | None = None
    feature_slot = None
    if normalized_linear_track_slots is not None:
        linear_track_layout = resolve_linear_track_layout(
            normalized_linear_track_slots,
            canvas_config=canvas_config,
            cfg=cfg,
        )
        canvas_config.set_linear_track_layout(linear_track_layout)
        feature_slot = _feature_slot_for_linear_slots(normalized_linear_track_slots)

    if required_label_height > 0:
        if canvas_config.vertical_offset < required_label_height:
            canvas_config.vertical_offset = (
                required_label_height + canvas_config.original_vertical_offset + canvas_config.cds_padding
            )
    else:
        canvas_config.vertical_offset = canvas_config.original_vertical_offset + canvas_config.cds_padding

    if records:
        base_axis_y = canvas_config.vertical_offset
        first_actual_extent = record_top_guard_above[0]
        first_normal_extent = record_top_guard_undisplaced[0]
        normal_top_margin = max(canvas_config.vertical_padding, base_axis_y - first_normal_extent)
        if track_layout == "above":
            middle_floor_extent = record_top_guard_middle_undisplaced[0]
            middle_floor_margin = max(canvas_config.vertical_padding, base_axis_y - middle_floor_extent)
            normal_top_margin = max(normal_top_margin, middle_floor_margin)
        if track_layout == "above" or (canvas_config.resolve_overlaps and track_layout == "middle"):
            required_axis_y = first_actual_extent + normal_top_margin
            # In above layout, keep at least the middle-layout top-margin floor so
            # non-stranded tracks do not end up visually too close to the top edge.
            canvas_config.vertical_offset = max(base_axis_y, required_axis_y)
        else:
            minimum_axis_y_for_features = first_actual_extent + canvas_config.vertical_padding
            canvas_config.vertical_offset = max(base_axis_y, minimum_axis_y_for_features)
        canvas_config.vertical_offset = max(
            canvas_config.vertical_offset,
            definition_half_heights[0] + canvas_config.vertical_padding,
        )
        if linear_track_layout is not None:
            custom_top_extent = float(linear_track_layout.top_extent)
            if feature_slot is not None and getattr(feature_slot, "side", "") == "above":
                custom_top_extent += record_heights_above[0]
                custom_top_extent += record_label_heights_above[0]
            canvas_config.vertical_offset = max(
                canvas_config.vertical_offset,
                custom_top_extent + canvas_config.vertical_padding,
            )

    normalize_length = cfg.canvas.linear.normalize_length
    if linear_track_layout is not None:
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
        record_gc_dfs = record_gc_dfs_by_nt.get(str(gc_config.dinucleotide).upper(), [None for _ in records])
    else:
        record_gc_dfs_by_nt = {}
        record_gc_dfs = _precalculate_gc_dataframes(
            records,
            window=int(gc_config.window),
            step=int(gc_config.step),
            dinucleotide=str(gc_config.dinucleotide),
            enabled=bool(canvas_config.show_gc or canvas_config.show_skew),
        )
    depth_enabled = bool(canvas_config.show_depth and depth_config is not None and record_depth_tracks)
    record_depth_data: list[list[DepthTrackData]] = build_depth_track_dataframes(
        records,
        record_depth_tracks,
        base_config=depth_config,
        depth_df_builder=build_depth_df,
    ) if depth_enabled else [[] for _ in records]
    if not depth_enabled:
        canvas_config.show_depth = False
        if linear_track_layout is None:
            canvas_config.set_gc_height_and_gc_padding()
    gc_axis_config = _gc_config_matching_linear_depth_axis_font_size(
        gc_config=gc_config,
        depth_config=depth_config,
        canvas_config=canvas_config,
        depth_enabled=depth_enabled,
    )

    # Prepare legend group
    has_blast = bool(normalized_comparisons)
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
        if depth_enabled:
            first_depth_row = next((row for row in record_depth_data if row), [])
            legend_table = sync_depth_track_legend_entries(legend_table, first_depth_row)
        legend_table = _sync_legend_table_for_linear_slots(
            legend_table,
            linear_track_slots=normalized_linear_track_slots,
            skew_config=skew_config,
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

    record_offsets: list[float] = []
    multi_record_plan = None

    if canvas_config.legend_position == "top":
        current_y = canvas_config.original_vertical_offset + required_legend_height + canvas_config.vertical_offset
    else:
        if canvas_config.vertical_offset > vertical_shift:
            current_y = canvas_config.vertical_offset
        else:
            current_y = canvas_config.original_vertical_offset + vertical_shift
    current_y += plot_title_top_reserve

    if multi_record_enabled:
        measurements: list[LinearRecordMeasurement] = []
        for i, record in enumerate(records):
            if linear_track_layout is not None:
                bottom_extent = _custom_record_bottom_extent(
                    record_index=i,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                )
                top_extent = _custom_record_top_extent(
                    record_index=i,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_above=record_heights_above,
                    record_label_heights_above=record_label_heights_above,
                )
            else:
                top_extent = max(record_label_heights_above[i], record_heights_above[i])
                bottom_extent = (
                    record_heights_below[i]
                    + record_label_heights_below[i]
                    + canvas_config.plot_tracks_height
                )
            comparison_endpoint_gap = 0.0
            if has_blast:
                comparison_endpoint_gap = float(canvas_config.vertical_padding)
                if non_middle_layout:
                    comparison_endpoint_gap = resolve_feature_axis_gap_linear(
                        cds_height=float(canvas_config.cds_height),
                        separate_strands=bool(canvas_config.strandedness),
                        axis_gap=canvas_config.track_axis_gap,
                    )
            comparison_top_extent = top_extent + comparison_endpoint_gap
            comparison_bottom_extent = bottom_extent + comparison_endpoint_gap
            top_extent = comparison_top_extent
            bottom_extent = comparison_bottom_extent
            # Multi-record definitions are record-local headers stacked above
            # every other record-local element. They reserve layout height but
            # may overlay comparison ribbons to avoid an empty header band.
            if float(_definition_heights[i]) > 0.0:
                top_extent += (
                    float(_definition_heights[i])
                    + float(canvas_config.vertical_padding)
                )
            measurements.append(
                LinearRecordMeasurement(
                    record_index=i,
                    record_key=RecordKey(str(record_keys[i])),
                    sequence_length=len(record.seq),
                    left_inset=0.0,
                    right_inset=0.0,
                    top_extent=top_extent,
                    bottom_extent=bottom_extent,
                    comparison_top_extent=comparison_top_extent,
                    comparison_bottom_extent=comparison_bottom_extent,
                )
            )
        multi_record_plan = solve_linear_layout(
            measurements,
            rows_by_record,
            available_width=float(canvas_config.alignment_width),
            record_gap_px=(
                linear_layout.record_gap_px if linear_layout is not None else 24.0
            ),
            align_center=bool(canvas_config.align_center),
            first_axis_y=current_y,
            row_gap_px=float(canvas_config.cds_padding) * 1.5,
            comparison_height=(
                float(canvas_config.comparison_height) if has_blast else 0.0
            ),
            record_order=ordered_record_indices,
        )
        record_offsets = [
            multi_record_plan.placement_for_index(index).axis_y
            for index in range(len(records))
        ]
        current_y = max(record_offsets)
        # Re-run label X placement with the final record-local sequence widths.
        _unused_height, all_labels, _unused_record_heights = _precalculate_label_dimensions(
            records,
            feature_config,
            canvas_config,
            config_dict,
            cfg=cfg,
            precomputed_feature_dicts=record_feature_dicts,
            orthogroup_label_eligibility=orthogroup_label_eligibility,
            sequence_widths=[
                multi_record_plan.placement_for_index(index).sequence_width
                for index in range(len(records))
            ],
        )
    else:
        for i, _record in enumerate(records):
            record_offsets.append(current_y)

            if i >= len(records) - 1:
                continue
            # Get the height below axis for the current record (feature tracks + GC/skew)
            current_feature_height_below = record_heights_below[i]
            current_label_height_below = record_label_heights_below[i]
            if linear_track_layout is not None:
                height_below_axis = float(linear_track_layout.bottom_extent)
                if feature_slot is not None and getattr(feature_slot, "side", "") == "below":
                    height_below_axis += (
                        current_feature_height_below
                        + current_label_height_below
                        + canvas_config.vertical_padding
                    )
                elif feature_slot is None or getattr(feature_slot, "side", "") == "overlay":
                    height_below_axis = max(
                        height_below_axis,
                        current_feature_height_below + current_label_height_below,
                    )
            else:
                height_below_axis = (
                    current_feature_height_below
                    + current_label_height_below
                    + canvas_config.plot_tracks_height
                )
                current_plot_stack_bottom_y = _record_plot_track_stack_bottom_y(
                    axis_y=current_y,
                    record_index=i,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                    non_middle_layout=non_middle_layout,
                )
                if current_plot_stack_bottom_y is not None:
                    height_below_axis = max(
                        height_below_axis,
                        current_plot_stack_bottom_y - current_y,
                    )

            # Get the height above axis for the next record (labels or feature tracks)
            next_label_height = record_label_heights_above[i + 1]
            next_feature_height_above = record_heights_above[i + 1]
            # For above-axis height, we need to consider both labels and the upper part of features
            height_above_next_axis = max(next_label_height, next_feature_height_above)
            if linear_track_layout is not None:
                custom_height_above = float(linear_track_layout.top_extent)
                if feature_slot is not None and getattr(feature_slot, "side", "") == "above":
                    custom_height_above += (
                        next_feature_height_above
                        + next_label_height
                        + canvas_config.vertical_padding
                    )
                elif feature_slot is None or getattr(feature_slot, "side", "") == "overlay":
                    custom_height_above = max(custom_height_above, next_label_height, next_feature_height_above)
                height_above_next_axis = max(height_above_next_axis, custom_height_above)
            # For BLAST comparisons, use comparison_height as minimum space between records
            if has_blast:
                min_gap = canvas_config.comparison_height
            else:
                # Use cds_padding * 1.5 as minimum gap to ensure clean separation with some breathing room
                min_gap = canvas_config.cds_padding * 1.5

            # Total inter-record space: below current + gap + above next
            inter_record_space = height_below_axis + min_gap + height_above_next_axis
            inter_record_space = max(
                inter_record_space,
                definition_half_heights[i]
                + definition_half_heights[i + 1],
            )
            if definition_half_heights[i] > 0.0 and definition_half_heights[i + 1] > 0.0:
                inter_record_space = max(
                    inter_record_space,
                    definition_half_heights[i]
                    + definition_half_heights[i + 1]
                    + max(1.0, 0.5 * float(canvas_config.vertical_padding)),
                )
            current_y += inter_record_space

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

    final_record_index = len(records) - 1 if records else -1
    final_feature_height_below = (
        record_heights_below[final_record_index]
        if non_middle_layout and final_record_index >= 0
        else canvas_config.cds_padding
    )
    final_label_height_below = (
        record_label_heights_below[final_record_index] if final_record_index >= 0 else 0.0
    )
    final_definition_height_below = (
        definition_half_heights[final_record_index] if final_record_index >= 0 else 0.0
    )
    if multi_record_plan is not None:
        final_record_height_below = max(
            placement.bottom_extent
            for placement in multi_record_plan.placements
            if placement.row == multi_record_plan.row_count - 1
        )
    elif linear_track_layout is not None:
        custom_final_height_below = float(linear_track_layout.bottom_extent)
        if feature_slot is not None and getattr(feature_slot, "side", "") == "below":
            custom_final_height_below += (
                final_feature_height_below
                + final_label_height_below
                + canvas_config.vertical_padding
            )
        elif feature_slot is None or getattr(feature_slot, "side", "") == "overlay":
            custom_final_height_below = max(
                custom_final_height_below,
                final_feature_height_below + final_label_height_below,
            )
        final_record_height_below = max(custom_final_height_below, final_definition_height_below)
    else:
        final_record_height_below = max(
            final_feature_height_below
            + final_label_height_below
            + canvas_config.plot_tracks_height,
            final_definition_height_below,
        )

    canvas_config.height_below_final_record = (
        current_y
        + final_record_height_below
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
    if multi_record_plan is not None:
        # The configurator's initial height assumes one record per row. Once a
        # multi-record plan exists, retaining that estimate leaves an empty
        # band for every record that shares a row. Size the canvas from the
        # resolved rows while still allowing a side legend to set the minimum.
        canvas_config.total_height = max(final_height, required_legend_height)
    else:
        canvas_config.total_height = max(final_height, canvas_config.total_height)

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
    record_offsets_x: list[float] = []
    for record_index, record in enumerate(records):
        if multi_record_plan is not None:
            record_offsets_x.append(
                multi_record_plan.placement_for_index(record_index).x
            )
            continue
        if normalize_length:
            record_offset_x = 0.0
        elif canvas_config.align_center:
            record_offset_x = (
                canvas_config.alignment_width
                * ((canvas_config.longest_genome - len(record.seq)) / canvas_config.longest_genome)
                / 2
            )
        else:
            record_offset_x = 0.0
        record_offset_x += orthogroup_alignment_offsets.get(len(record_offsets_x), 0.0)
        record_offsets_x.append(record_offset_x)
    definition_column_width = max_def_width
    if canvas_config.keep_definition_left_aligned and record_offsets_x:
        definition_column_width = max(0.0, float(max_def_width) - min(record_offsets_x))
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

    explicit_placements: dict[int, LinearRecordPlacement] | None = None
    if has_blast and linear_comparisons and multi_record_plan is None:
        explicit_placements = {}
        for record_index, record in enumerate(records):
            axis_y = record_offsets[record_index]
            if linear_track_layout is not None:
                top_y = axis_y - _custom_record_top_extent(
                    record_index=record_index,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_above=record_heights_above,
                    record_label_heights_above=record_label_heights_above,
                )
                bottom_y = axis_y + _custom_record_bottom_extent(
                    record_index=record_index,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                )
            elif non_middle_layout:
                top_y = bottom_y = axis_y
            else:
                top_y = axis_y - canvas_config.cds_padding
                bottom_y = _record_plot_track_stack_bottom_y(
                    axis_y=axis_y,
                    record_index=record_index,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                    non_middle_layout=non_middle_layout,
                ) or (axis_y + canvas_config.cds_padding)
            sequence_width = (
                float(canvas_config.alignment_width)
                if canvas_config.normalize_length
                else float(canvas_config.alignment_width)
                * len(record.seq)
                / max(1, canvas_config.longest_genome)
            )
            explicit_placements[record_index] = LinearRecordPlacement(
                record_index=record_index,
                record_key=RecordKey(str(record_keys[record_index])),
                row=rows_by_record[record_index],
                column=0,
                x=record_offsets_x[record_index],
                axis_y=axis_y,
                sequence_width=sequence_width,
                left_inset=0.0,
                right_inset=0.0,
                top_extent=max(0.0, axis_y - top_y),
                bottom_extent=max(0.0, bottom_y - axis_y),
                comparison_top_y=top_y,
                comparison_bottom_y=bottom_y,
                px_per_bp=sequence_width / len(record.seq),
            )

    if has_blast and (multi_record_plan is not None or explicit_placements is not None):
        canvas = add_explicit_comparisons_on_linear_canvas(
            canvas,
            normalized_comparisons,
            canvas_config,
            blast_config,
            records,
            explicit_placements or {
                placement.record_index: placement
                for placement in multi_record_plan.placements
            },
        )
    elif has_blast:
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(records) - 1):
            if linear_track_layout is not None:
                ribbon_start_y = record_offsets[i] + _custom_record_bottom_extent(
                    record_index=i,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                )
                ribbon_end_y = record_offsets[i + 1] - _custom_record_top_extent(
                    record_index=i + 1,
                    layout=linear_track_layout,
                    feature_slot=feature_slot,
                    canvas_config=canvas_config,
                    record_heights_above=record_heights_above,
                    record_label_heights_above=record_label_heights_above,
                )
            elif non_middle_layout:
                # In above/below layouts, anchor ribbons to consecutive record axes.
                ribbon_start_y = record_offsets[i]
                ribbon_end_y = record_offsets[i + 1]
            else:
                plot_stack_bottom_y = _record_plot_track_stack_bottom_y(
                    axis_y=record_offsets[i],
                    record_index=i,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_label_heights_below=record_label_heights_below,
                    non_middle_layout=non_middle_layout,
                )
                if plot_stack_bottom_y is None:
                    ribbon_start_y = record_offsets[i] + canvas_config.cds_padding
                else:
                    ribbon_start_y = plot_stack_bottom_y
                ribbon_end_y = record_offsets[i + 1] - canvas_config.cds_padding
            comparison_offsets.append(ribbon_start_y)
            height = max(0.0, ribbon_end_y - ribbon_start_y)
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
            orthogroup_alignment_offsets,
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
        offset_y = record_offsets[record_index]
        record_placement = (
            multi_record_plan.placement_for_index(record_index)
            if multi_record_plan is not None
            else None
        )
        sequence_width = (
            record_placement.sequence_width if record_placement is not None else None
        )

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
            shared_depth_tracks = record_depth_data[record_index] if record_index < len(record_depth_data) else []
            feature_rendered = False

            for slot in sorted(linear_track_layout.slots, key=lambda item: (item.z, item.slot_index)):
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
                    )
                    feature_rendered = True
                    continue

                track_offset_y = _custom_slot_track_offset_y(
                    slot,
                    feature_slot=feature_slot,
                    record_index=record_index,
                    canvas_config=canvas_config,
                    record_heights_below=record_heights_below,
                    record_heights_above=record_heights_above,
                    record_label_heights_below=record_label_heights_below,
                    record_label_heights_above=record_label_heights_above,
                )
                if slot.renderer == "annotations":
                    annotation_layout = annotation_track_layouts.get(slot.id)
                    if annotation_layout is None:
                        continue
                    params = annotation_track_params_from_mapping(slot.params)
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
                        height_px=slot.height,
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
                    track_index = int(slot.params.get("track_index", 0))
                    if track_index < 0 or track_index >= len(shared_depth_tracks):
                        raise ValueError(
                            f"Depth slot '{slot.id}' track_index={track_index} is outside the available depth tracks."
                        )
                    depth_track = shared_depth_tracks[track_index]
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
                        depth_track_index=track_index,
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
                    slot_gc_config = _gc_config_matching_linear_depth_axis_font_size(
                        gc_config=_clone_gc_config_with_dinucleotide(gc_config, nt),
                        depth_config=depth_config,
                        canvas_config=canvas_config,
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
                    offset_y + record_feature_center_offsets[record_index]
                ),
            )
            continue

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
            label_font_size=label_font_size,
            orthogroup_label_member_ids=orthogroup_label_member_ids,
            orthogroup_label_top_member_ids=orthogroup_label_top_member_ids,
            record_index=record_index,
            record_count=total_records,
            group_id=record_group_id,
            placement=record_placement,
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
                offset_y + record_feature_center_offsets[record_index]
            ),
        )
        gc_offset_y = offset_y
        if non_middle_layout:
            current_feature_height_below = record_heights_below[record_index]
            gc_offset_y = offset_y + (current_feature_height_below - canvas_config.cds_padding)
        gc_offset_y += record_label_heights_below[record_index]
        shared_gc_df = record_gc_dfs[record_index] if record_index < len(record_gc_dfs) else None
        shared_depth_tracks = record_depth_data[record_index] if record_index < len(record_depth_data) else []

        if depth_enabled:
            for depth_track_index, depth_track in enumerate(shared_depth_tracks):
                group_id = _linear_depth_group_id(
                    depth_track.id,
                    record_index=record_index,
                    record_count=total_records,
                )
                add_depth_group(
                    canvas,
                    record,
                    gc_offset_y,
                    offset_x,
                    canvas_config,
                    depth_track.config,
                    config_dict,
                    cfg=record_cfg,
                    depth_df=depth_track.df,
                    depth_track_index=depth_track_index,
                    group_id=group_id,
                    axis_group_id=f"{group_id}_axis",
                    track_height=depth_track.height,
                    sequence_width=sequence_width,
                )
        if canvas_config.show_gc:
            gc_group_id = make_linear_dom_id(
                "gc_content",
                record_index=record_index,
                record_count=total_records,
            )
            add_gc_content_group(
                canvas,
                record,
                gc_offset_y,
                offset_x,
                canvas_config,
                gc_axis_config,
                config_dict,
                cfg=record_cfg,
                gc_df=shared_gc_df,
                group_id=gc_group_id,
                sequence_width=sequence_width,
            )
        if canvas_config.show_skew:
            skew_group_id = make_linear_dom_id(
                "gc_skew",
                record_index=record_index,
                record_count=total_records,
            )
            add_gc_skew_group(
                canvas,
                record,
                gc_offset_y,
                offset_x,
                canvas_config,
                skew_config,
                config_dict,
                cfg=record_cfg,
                gc_df=shared_gc_df,
                group_id=skew_group_id,
                sequence_width=sequence_width,
            )

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
                record_offsets=record_offsets,
                feature_slot=feature_slot,
                canvas_config=canvas_config,
                record_heights_below=record_heights_below,
                record_heights_above=record_heights_above,
                record_label_heights_below=record_label_heights_below,
                record_label_heights_above=record_label_heights_above,
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
