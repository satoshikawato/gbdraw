#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

import logging
import math
import copy
from dataclasses import dataclass, replace
from typing import Any, Optional, Mapping, Sequence, cast

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...analysis.conservation import (  # type: ignore[reportMissingImports]
    ConservationTrack,
    conservation_track_gradient_colors,
)
from ...analysis.depth_tracks import (  # type: ignore[reportMissingImports]
    DepthTrackData,
    depth_track_data_count,
    index_depth_track_row,
    sync_depth_track_legend_entries,
)
from ...config.models import (  # type: ignore[reportMissingImports]
    GbdrawConfig,
)
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
    LegendMeasurement,
)
from ...configurators.gc import _slot_skew_config
from ...core.sequence import check_feature_presence  # type: ignore[reportMissingImports]
from ...core.text import calculate_bbox_dimensions, calculate_svg_bbox_dimensions
from ...exceptions import ValidationError
from ...features.colors import precompute_used_color_rules  # type: ignore[reportMissingImports]
from ...features.factory import FeatureBuildResult, create_feature_layers  # type: ignore[reportMissingImports]
from ...labels.circular import (  # type: ignore[reportMissingImports]
    assign_leader_start_points,
    minimum_bbox_gap_px,
    prepare_label_list,
    x_overlap,
    y_overlap,
)
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...layout.circular import (  # type: ignore[reportMissingImports]
    CircularFeatureLaneDirection,
    CircularRadialLayout,
    CircularRecordRenderContext,
    CircularResolvedSlot,
    CircularTickLayout,
)
from ...layout.circular_depth_axis import depth_axis_tick_font_size_px  # type: ignore[reportMissingImports]
from ...layout.composition import (  # type: ignore[reportMissingImports]
    CompositionItem,
    CompositionPlan,
    CompositionRequest,
    DEFAULT_COMPOSITION_SPACING,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from ...layout.spatial import Aabb, union_aabbs  # type: ignore[reportMissingImports]
from ...legend.table import _unique_legend_key, prepare_legend_table  # type: ignore[reportMissingImports]
from ...tracks import (  # type: ignore[reportMissingImports]
    CircularTrackSlot,
    ScalarSpec,
    default_circular_track_slots,
    normalize_circular_track_slots,
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
    sync_annotation_legend_entries,
)
from ...annotations.planning import prepare_annotation_track_slots
from ...render.drawers.circular.annotations import draw_circular_annotation_track
from ...render.composition import apply_composition_plan
from ...svg.ids import stable_svg_id, track_slot_svg_id
from .positioning import _parse_svg_number as _svg_number, center_group_on_canvas
from ...tracks.circular import tick_sides_for_tick_label_layout  # type: ignore[reportMissingImports]

from .builders import (
    add_axis_group_on_canvas,
    add_conservation_group_on_canvas,
    add_depth_group_on_canvas,
    add_gc_content_group_on_canvas,
    add_gc_skew_group_on_canvas,
    add_labels_group_on_canvas,
    add_legend_group_on_canvas,
    add_record_definition_group_on_canvas,
    add_record_group_on_canvas,
    add_tick_group_on_canvas,
)
from ...render.groups.circular.definition import DefinitionGroup  # type: ignore[reportMissingImports]
from .radial_layout import resolve_circular_radial_layout  # type: ignore[reportMissingImports]
from .presets import (  # type: ignore[reportMissingImports]
    CircularPresetContext,
    circular_feature_lane_direction_for_preset,
    circular_radial_plan_for_preset,
    circular_track_slots_from_preset_order,
    normalize_circular_track_preset,
)


LEGEND_LABEL_MARGIN_PX = 2.0 * DEFAULT_COMPOSITION_SPACING.overlay_clearance_px
LABEL_NUDGE_STEP_PX = 6.0
MAX_LABEL_NUDGE_PX = 180.0
MIN_LABEL_ORDER_GAP_RAD = 1e-4
FEATURE_BAND_EPSILON = 1e-6
logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class _CircularPlotAssembly:
    """Uncomposed Circular plot geometry and its existing SVG targets."""

    drawing: Drawing
    canvas_config: CircularCanvasConfigurator
    legend_measurement: LegendMeasurement
    legend_table: Mapping[object, Mapping[str, object]]
    source_content_bounds: Aabb
    source_overlay_obstacles: tuple[Aabb, ...]
    primary_targets: tuple[Group, ...]
    legend_target: Group | None
    track_slot_geometry: Mapping[str, Any]


@dataclass(frozen=True, slots=True)
class CircularAssemblyResult:
    """Final single-record Circular assembly with authoritative geometry."""

    drawing: Drawing
    legend_measurement: LegendMeasurement
    legend_table: Mapping[object, Mapping[str, object]]
    source_content_bounds: Aabb
    content_bounds: Aabb
    overlay_obstacles: tuple[Aabb, ...]
    primary_targets: tuple[Group, ...]
    legend_target: Group | None
    title_target: Group | None
    track_slot_geometry: Mapping[str, Any]
    composition_plan: CompositionPlan


def _mark_circular_track_slot_group(
    canvas: Drawing,
    *,
    source_group_id: str,
    group_id: str,
    slot_id: str,
    renderer: str,
) -> Drawing:
    """Namespace a newly added slot group and its renderer-owned child IDs."""

    for element in reversed(getattr(canvas, "elements", [])):
        attribs = getattr(element, "attribs", None)
        if isinstance(attribs, dict) and str(attribs.get("id", "")) == source_group_id:
            subtree: list[Any] = []

            def collect(node: Any) -> None:
                subtree.append(node)
                for child in getattr(node, "elements", []) or []:
                    collect(child)

            collect(element)
            id_map: dict[str, str] = {}
            for child_index, child in enumerate(subtree):
                child_attribs = getattr(child, "attribs", None)
                if not isinstance(child_attribs, dict):
                    continue
                old_id = str(child_attribs.get("id", ""))
                if not old_id:
                    continue
                new_id = (
                    group_id
                    if child is element
                    else stable_svg_id(
                        "track_slot_child",
                        group_id,
                        old_id,
                        child_index,
                    )
                )
                child_attribs["id"] = new_id
                id_map[old_id] = new_id

            for child in subtree:
                child_attribs = getattr(child, "attribs", None)
                if not isinstance(child_attribs, dict):
                    continue
                for name, value in tuple(child_attribs.items()):
                    if not isinstance(value, str):
                        continue
                    updated = value
                    for old_id, new_id in id_map.items():
                        updated = updated.replace(
                            f"url(#{old_id})",
                            f"url(#{new_id})",
                        )
                        if updated == f"#{old_id}":
                            updated = f"#{new_id}"
                    if updated != value:
                        child_attribs[name] = updated

            attribs["data-gbdraw-slot-id"] = slot_id
            attribs["data-gbdraw-slot-renderer"] = renderer
            break
    return canvas


def _tag_circular_track_slot_group(
    canvas: Drawing,
    *,
    source_group_id: str,
    slot_id: str,
    renderer: str,
) -> Drawing:
    """Attach semantic slot metadata without changing legacy default IDs."""

    for element in reversed(getattr(canvas, "elements", [])):
        attribs = getattr(element, "attribs", None)
        if isinstance(attribs, dict) and str(attribs.get("id", "")) == source_group_id:
            attribs["data-gbdraw-slot-id"] = str(slot_id)
            attribs["data-gbdraw-slot-renderer"] = str(renderer)
            break
    return canvas


def _prepare_circular_annotation_tracks(
    gb_record: SeqRecord,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    slots: list[CircularTrackSlot] | None,
    *,
    canvas_config: CircularCanvasConfigurator,
    record_index: int,
    show_ticks: bool,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
    depth_track_count: int,
) -> tuple[list[CircularTrackSlot] | None, ResolvedAnnotationBundle, dict[str, ResolvedAnnotationTrack]]:
    slots, bundle, _auto_slot_ids = prepare_annotation_track_slots(
        annotations,
        [gb_record],
        slots,
        mode="circular",
        default_slots=lambda: default_circular_track_slots(
            show_features=True,
            show_ticks=show_ticks,
            show_depth=show_depth,
            depth_track_count=depth_track_count,
            show_gc=show_gc,
            show_skew=show_skew,
        ),
        slot_factory=CircularTrackSlot,
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, {}
    relevant = tuple(item for item in bundle.annotations if item.record_index == record_index)

    record_lengths = {record_index: len(gb_record.seq)}
    bp_per_px = {
        record_index: len(gb_record.seq) / max(1.0, 2.0 * math.pi * float(canvas_config.radius))
    }
    updated: list[CircularTrackSlot] = []
    layouts: dict[str, ResolvedAnnotationTrack] = {}
    for slot in slots:
        if str(slot.renderer).strip().lower() != "annotations":
            updated.append(slot)
            continue
        params = annotation_track_params_from_mapping(slot.params)
        available = slot.width.resolve(float(canvas_config.radius)) if slot.width is not None else None
        layout = layout_annotation_track(
            slot.id,
            params.set_id,
            relevant,
            record_lengths=record_lengths,
            params=params,
            available_extent_px=available,
            bp_per_px=bp_per_px,
        )
        layouts[slot.id] = layout
        updated.append(
            slot
            if slot.width is not None
            else replace(slot, width=ScalarSpec(layout.required_extent_px, "px"))
        )
    return updated, bundle, layouts


def _add_circular_feature_underlays(
    gb_record: SeqRecord,
    underlay_features: Sequence,
    slots: list[CircularTrackSlot],
    bundle: ResolvedAnnotationBundle,
    *,
    canvas_config: CircularCanvasConfigurator,
    record_index: int,
    show_ticks: bool,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
    depth_track_count: int,
) -> tuple[
    list[CircularTrackSlot],
    ResolvedAnnotationBundle,
    dict[str, ResolvedAnnotationTrack],
]:
    merged, set_id = merge_feature_underlays(
        bundle,
        [underlay_features],
        [gb_record],
        mode="circular",
        record_indices=[record_index],
    )
    if set_id is None:
        return slots, bundle, {}
    anchor_id = feature_underlay_anchor_slot_id(slots)
    anchor = next(slot for slot in slots if str(slot.id) == anchor_id)
    slot_id = feature_underlay_slot_id(slots)
    underlay_slot = CircularTrackSlot(
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
    updated, merged, layouts = _prepare_circular_annotation_tracks(
        gb_record,
        merged,
        [underlay_slot, *slots],
        canvas_config=canvas_config,
        record_index=record_index,
        show_ticks=show_ticks,
        show_depth=show_depth,
        show_gc=show_gc,
        show_skew=show_skew,
        depth_track_count=depth_track_count,
    )
    return list(updated or []), merged, layouts


def _circular_preset_for_layout(
    cfg: GbdrawConfig,
) -> str:
    return normalize_circular_track_preset(str(cfg.canvas.circular.track_type))


def _lane_direction_for_feature_slot(
    feature_slot: CircularTrackSlot | None,
    *,
    track_preset: str,
) -> CircularFeatureLaneDirection:
    if feature_slot is not None:
        raw = feature_slot.params.get("lane_direction", feature_slot.params.get("lanes"))
        if raw is not None:
            direction = str(raw).strip().lower()
            if direction in {"inside", "outside", "split"}:
                return cast(CircularFeatureLaneDirection, direction)
    return circular_feature_lane_direction_for_preset(
        track_preset
    )


def _serialize_circular_track_slot_geometry(
    *,
    gb_record: SeqRecord,
    radial_layout: CircularRadialLayout,
    layout_slots: Sequence[CircularTrackSlot],
    base_radius_px: float,
    record_index: int = 0,
) -> dict[str, Any]:
    source_by_index = {
        int(getattr(slot, "slot_index", index)): slot
        for index, slot in enumerate(layout_slots)
    }
    safe_base_radius = max(FEATURE_BAND_EPSILON, float(base_radius_px))
    slots_payload: list[dict[str, Any]] = []
    for resolved in radial_layout.slots:
        source_slot = source_by_index.get(int(resolved.slot_index))
        width_spec = getattr(source_slot, "width", None)
        width_factor = None
        if width_spec is not None and getattr(width_spec, "unit", None) == "factor":
            width_factor = float(width_spec.value)
        resolved_width = (
            float(resolved.resolved_width_px)
            if resolved.resolved_width_px is not None
            else float(resolved.draw_width_px)
        )
        slots_payload.append(
            {
                "slotIndex": int(resolved.slot_index),
                "slotId": str(resolved.id),
                "renderer": str(resolved.renderer),
                "side": str(resolved.side),
                "widthPx": resolved_width,
                "widthFactor": width_factor,
                "radiusFactor": float(resolved.center_radius_px) / safe_base_radius,
                "innerGapPx": float(resolved.inner_gap_px),
                "outerGapPx": float(resolved.outer_gap_px),
                "compressed": bool(resolved.compressed),
                "source": "resolved",
            }
        )
    record_id = str(getattr(gb_record, "id", "") or "")
    return {
        "schema": 1,
        "mode": "circular",
        "source": "resolved",
        "records": [
            {
                "recordIndex": int(record_index),
                "recordId": record_id,
                "recordLabel": record_id,
                "axisRadiusPx": float(base_radius_px),
                "slots": slots_payload,
            }
        ],
    }


def _legend_bbox(legend_bounds: Aabb) -> tuple[float, float, float, float]:
    """Return legend bounds as a collision tuple."""
    return (
        legend_bounds.min_x,
        legend_bounds.min_y,
        legend_bounds.max_x,
        legend_bounds.max_y,
    )


def _bbox_overlaps(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    """Return True when two axis-aligned bboxes overlap."""
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def _unwrap_angle_near_reference(angle_rad: float, reference_rad: float) -> float:
    """Project angle to the nearest 2*pi branch around reference."""
    return angle_rad + (2.0 * math.pi) * round((reference_rad - angle_rad) / (2.0 * math.pi))


def _label_target_unwrapped_angle(label: dict[str, Any], total_length: int) -> float:
    """Target label angle derived from genomic midpoint, in unwrapped radians."""
    return (2.0 * math.pi * (float(label["middle"]) / float(total_length))) - (0.5 * math.pi)


def _label_unwrapped_angle(label: dict[str, Any], total_length: int) -> float:
    """Current label angle projected near its target branch."""
    current_angle = math.atan2(float(label["start_y"]), float(label["start_x"]))
    target_angle = _label_target_unwrapped_angle(label, total_length)
    return _unwrap_angle_near_reference(current_angle, target_angle)


def _label_bbox_local(label: dict[str, Any], total_length: int, margin_px: float) -> tuple[float, float, float, float]:
    """Approximate label bbox in local (centered) coordinates."""
    authoritative_aabb = label.get("aabb_local")
    if authoritative_aabb is not None and len(authoritative_aabb) == 4:
        half_margin = 0.5 * float(margin_px)
        return (
            float(authoritative_aabb[0]) - half_margin,
            float(authoritative_aabb[1]) - half_margin,
            float(authoritative_aabb[2]) + half_margin,
            float(authoritative_aabb[3]) + half_margin,
        )
    half_margin = 0.5 * margin_px
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    width_px = float(label["width_px"])
    height_px = float(label["height_px"])
    is_inner = bool(label.get("is_inner", False))

    if not is_inner:
        if start_x > 0:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin
        else:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
    else:
        if start_x > 0:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
        else:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin

    angle_deg = (360.0 * (float(label["middle"]) / float(total_length))) % 360.0
    if not is_inner:
        if 0.0 <= angle_deg < 10.0 or angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        elif 170.0 <= angle_deg < 190.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    else:
        if 0.0 <= angle_deg < 10.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y
        elif angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    return min_x, min_y, max_x, max_y


def _label_bbox_on_canvas(
    label: dict[str, Any],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    margin_px: float,
) -> tuple[float, float, float, float]:
    """Return label bbox on canvas as (min_x, min_y, max_x, max_y)."""
    min_x, min_y, max_x, max_y = _label_bbox_local(label, total_length, margin_px)
    return (
        min_x + float(canvas_config.offset_x),
        min_y + float(canvas_config.offset_y),
        max_x + float(canvas_config.offset_x),
        max_y + float(canvas_config.offset_y),
    )


def _move_label_text_anchor(
    label: dict[str, Any],
    *,
    start_x: float,
    start_y: float,
) -> None:
    """Move label text geometry while keeping authoritative bounds coherent."""
    old_x = float(label.get("start_x", 0.0))
    old_y = float(label.get("start_y", 0.0))
    dx = float(start_x) - old_x
    dy = float(start_y) - old_y
    label["start_x"] = float(start_x)
    label["start_y"] = float(start_y)

    for x_key, y_key in (("text_x", "text_y"), ("leader_start_x", "leader_start_y")):
        if x_key in label:
            label[x_key] = float(label[x_key]) + dx
        if y_key in label:
            label[y_key] = float(label[y_key]) + dy

    authoritative_aabb = label.get("aabb_local")
    if authoritative_aabb is not None and len(authoritative_aabb) == 4:
        label["aabb_local"] = (
            float(authoritative_aabb[0]) + dx,
            float(authoritative_aabb[1]) + dy,
            float(authoritative_aabb[2]) + dx,
            float(authoritative_aabb[3]) + dy,
        )
    oriented_corners = label.get("oriented_corners")
    if oriented_corners is not None:
        label["oriented_corners"] = tuple(
            (float(point[0]) + dx, float(point[1]) + dy)
            for point in oriented_corners
        )


def _external_label_bounds_on_canvas(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    *,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> tuple[float, float, float, float] | None:
    """Return bounds of all non-embedded labels in canvas coordinates."""
    min_x: float | None = None
    min_y: float | None = None
    max_x: float | None = None
    max_y: float | None = None

    for label in labels:
        if label.get("is_embedded"):
            continue
        box_min_x, box_min_y, box_max_x, box_max_y = _label_bbox_on_canvas(
            label,
            total_length,
            canvas_config,
            margin_px=margin_px,
        )
        min_x = box_min_x if min_x is None else min(min_x, box_min_x)
        min_y = box_min_y if min_y is None else min(min_y, box_min_y)
        max_x = box_max_x if max_x is None else max(max_x, box_max_x)
        max_y = box_max_y if max_y is None else max(max_y, box_max_y)

    if min_x is None or min_y is None or max_x is None or max_y is None:
        return None
    return min_x, min_y, max_x, max_y


def _external_label_paint_bounds_local(
    label: Mapping[str, Any],
    total_length: int,
    *,
    leader_half_stroke_px: float,
) -> Aabb:
    """Return exact text plus leader paint bounds in circle-local coordinates."""
    measured_label = dict(label)
    if measured_label.get("aabb_local") is None:
        text_width, text_height = calculate_svg_bbox_dimensions(
            str(measured_label.get("label_text", "")),
            str(measured_label.get("font_family", "sans-serif")),
            float(measured_label.get("font_size", 0.0)),
            96,
            font_weight="normal",
            font_style="normal",
        )
        if text_width > 0.0:
            measured_label["width_px"] = float(text_width)
        if text_height > 0.0:
            measured_label["height_px"] = float(text_height)
    text_bounds = Aabb(
        *_label_bbox_local(measured_label, total_length, margin_px=0.0)
    )
    leader_points = (
        (
            float(label.get("middle_x", label.get("start_x", 0.0))),
            float(label.get("middle_y", label.get("start_y", 0.0))),
        ),
        (
            float(label.get("leader_start_x", label.get("start_x", 0.0))),
            float(label.get("leader_start_y", label.get("start_y", 0.0))),
        ),
        (
            float(label.get("feature_anchor_x", label.get("feature_middle_x", 0.0))),
            float(label.get("feature_anchor_y", label.get("feature_middle_y", 0.0))),
        ),
    )
    leader_bounds = Aabb(
        min(point[0] for point in leader_points),
        min(point[1] for point in leader_points),
        max(point[0] for point in leader_points),
        max(point[1] for point in leader_points),
    ).expanded(float(leader_half_stroke_px))
    combined = union_aabbs((text_bounds, leader_bounds))
    if combined is None:  # pragma: no cover - both inputs are present
        raise RuntimeError("external label paint bounds unexpectedly empty")
    return combined


def _external_label_obstacles_on_canvas(
    labels: Sequence[Mapping[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
) -> tuple[Aabb, ...]:
    """Return authoritative per-label overlay obstacles in plot coordinates."""
    half_stroke = 0.5 * float(
        canvas_config.profile.config.labels.stroke_width.for_length_param(
            canvas_config.length_param
        )
    )
    return tuple(
        _external_label_paint_bounds_local(
            label,
            total_length,
            leader_half_stroke_px=half_stroke,
        ).translated(float(canvas_config.offset_x), float(canvas_config.offset_y))
        for label in labels
        if not label.get("is_embedded")
    )


def _primary_radial_half_stroke_px(
    *,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
) -> float:
    """Return conservative stroke safety for radial record/track geometry."""
    cfg = canvas_config.profile.config
    widths = [
        float(
            cfg.objects.axis.circular.stroke_width.for_length_param(
                canvas_config.length_param
            )
        ),
        float(feature_config.block_stroke_width),
        float(feature_config.line_stroke_width),
        float(gc_config.stroke_width),
        float(skew_config.stroke_width),
    ]
    if depth_config is not None:
        widths.append(float(depth_config.stroke_width))
    conservation = getattr(cfg.objects, "conservation", None)
    if conservation is not None:
        widths.append(float(getattr(conservation, "stroke_width", 0.0)))
    return 0.5 * max(widths, default=0.0)


def _collect_circular_primary_bounds(
    *,
    radial_layout: CircularRadialLayout,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
    overlay_obstacles: Sequence[Aabb],
    definition_target: Group | None = None,
) -> Aabb:
    """Collect authoritative plot-space bounds before outer composition."""
    stroke_safety = _primary_radial_half_stroke_px(
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
    )
    radius = float(radial_layout.outer_content_radius_px) + stroke_safety
    center_x = float(canvas_config.offset_x)
    center_y = float(canvas_config.offset_y)
    bounds: list[Aabb] = [
        Aabb(
            center_x - radius,
            center_y - radius,
            center_x + radius,
            center_y + radius,
        )
    ]
    bounds.extend(overlay_obstacles)
    if definition_target is not None:
        local_bounds = getattr(definition_target, "_gbdraw_local_bounds", None)
        translation = getattr(definition_target, "_gbdraw_plot_translation", None)
        if isinstance(local_bounds, Aabb) and translation is not None:
            bounds.append(
                local_bounds.translated(float(translation[0]), float(translation[1]))
            )
    combined = union_aabbs(bounds)
    if combined is None:  # pragma: no cover - radial bounds are always present
        raise RuntimeError("circular primary bounds unexpectedly empty")
    return combined


def _overlay_legend_bounds_for_collision(
    *,
    legend_measurement: LegendMeasurement,
    primary_bounds: Aabb,
    placement: LegendPlacement,
) -> Aabb:
    """Return the composer's requested corner bounds before label adjustment."""
    if not placement.is_overlay:
        raise ValueError("collision bounds require an overlay legend placement")
    plan = plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", primary_bounds),
            legend=CompositionItem("legend", legend_measurement.local_bounds),
            legend_placement=placement,
        )
    )
    primary = plan.placement_for("primary")
    legend = plan.placement_for("legend")
    if primary is None or legend is None:  # pragma: no cover - request has both items
        raise RuntimeError("overlay collision plan is incomplete")
    return legend.final_bounds.translated(-primary.dx, -primary.dy)


def _legend_collision_indices(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> list[int]:
    """Return indices of labels that currently collide with the legend."""
    legend_box = _legend_bbox(legend_bounds)
    collided: list[int] = []
    for idx, label in enumerate(labels):
        if label.get("is_embedded"):
            continue
        label_box = _label_bbox_on_canvas(label, total_length, canvas_config, margin_px=margin_px)
        if _bbox_overlaps(label_box, legend_box):
            collided.append(idx)
    return collided


def _labels_collide_with_legend(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> bool:
    """Whether any external label overlaps with the legend bbox."""
    return bool(
        _legend_collision_indices(
            labels,
            total_length,
            canvas_config,
            legend_bounds,
        )
    )




def _legend_center_local(
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> tuple[float, float]:
    """Return legend center in local (circle-centered) coordinates."""
    legend_min_x, legend_min_y, legend_max_x, legend_max_y = _legend_bbox(legend_bounds)
    legend_center_x = 0.5 * (legend_min_x + legend_max_x)
    legend_center_y = 0.5 * (legend_min_y + legend_max_y)
    return legend_center_x - float(canvas_config.offset_x), legend_center_y - float(canvas_config.offset_y)


def _preferred_angular_shift_sign(
    label: dict[str, Any],
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> int:
    """Return preferred angular direction (+1 ccw / -1 cw) to move away from legend."""
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    radius = math.hypot(start_x, start_y)
    if radius <= 1e-6:
        return 1
    legend_local_x, legend_local_y = _legend_center_local(
        canvas_config,
        legend_bounds,
    )
    away_x = start_x - legend_local_x
    away_y = start_y - legend_local_y
    # Unit tangents at current point (CCW and CW).
    ccw_tx = -start_y / radius
    ccw_ty = start_x / radius
    cw_tx = -ccw_tx
    cw_ty = -ccw_ty
    ccw_score = ccw_tx * away_x + ccw_ty * away_y
    cw_score = cw_tx * away_x + cw_ty * away_y
    return 1 if ccw_score >= cw_score else -1


def _expand_shift_block_to_preserve_order(
    labels: list[dict[str, Any]],
    total_length: int,
    side_indices: list[int],
    center_pos: int,
    delta_angle: float,
) -> tuple[int, int]:
    """Expand move block so angular order cannot be inverted by this shift."""
    left = center_pos
    right = center_pos
    while True:
        expanded = False
        moved_left_idx = side_indices[left]
        moved_right_idx = side_indices[right]
        moved_left_angle = _label_unwrapped_angle(labels[moved_left_idx], total_length) + delta_angle
        moved_right_angle = _label_unwrapped_angle(labels[moved_right_idx], total_length) + delta_angle

        if left > 0:
            prev_idx = side_indices[left - 1]
            prev_angle = _label_unwrapped_angle(labels[prev_idx], total_length)
            if moved_left_angle <= prev_angle + MIN_LABEL_ORDER_GAP_RAD:
                left -= 1
                expanded = True

        if right < len(side_indices) - 1:
            next_idx = side_indices[right + 1]
            next_angle = _label_unwrapped_angle(labels[next_idx], total_length)
            if moved_right_angle >= next_angle - MIN_LABEL_ORDER_GAP_RAD:
                right += 1
                expanded = True

        if not expanded:
            break
    return left, right


def _try_shift_labels_away_from_legend(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> bool:
    """Resolve collisions by moving labels first (preferred strategy)."""
    if not labels:
        return True

    max_passes = 4
    for _ in range(max_passes):
        collided_indices = _legend_collision_indices(
            labels,
            total_length,
            canvas_config,
            legend_bounds,
        )
        if not collided_indices:
            return True

        changed = False
        legend_box = _legend_bbox(legend_bounds)
        for idx in collided_indices:
            label = labels[idx]
            if label.get("is_embedded"):
                continue
            start_x = float(label["start_x"])
            start_y = float(label["start_y"])
            radius = math.hypot(start_x, start_y)
            if radius <= 1e-6:
                continue
            moving_side = 1 if start_x >= 0 else -1
            preferred_sign = _preferred_angular_shift_sign(
                label,
                canvas_config,
                legend_bounds,
            )
            side_indices = [
                side_idx
                for side_idx, side_label in enumerate(labels)
                if not side_label.get("is_embedded") and ((float(side_label["start_x"]) >= 0) == (moving_side > 0))
            ]
            side_indices.sort(key=lambda side_idx: float(labels[side_idx]["middle"]))
            if idx not in side_indices:
                continue
            center_pos = side_indices.index(idx)

            placed = False
            shift_px = LABEL_NUDGE_STEP_PX
            while shift_px <= MAX_LABEL_NUDGE_PX and not placed:
                angle_delta = shift_px / radius
                for direction_sign in (preferred_sign, -preferred_sign):
                    signed_delta = direction_sign * angle_delta
                    block_left, block_right = _expand_shift_block_to_preserve_order(
                        labels, total_length, side_indices, center_pos, signed_delta
                    )
                    block_indices = side_indices[block_left : block_right + 1]

                    candidate_positions: dict[int, tuple[float, float]] = {}
                    candidate_labels: dict[int, dict[str, Any]] = {}
                    candidate_valid = True

                    for block_idx in block_indices:
                        block_label = labels[block_idx]
                        block_radius = math.hypot(float(block_label["start_x"]), float(block_label["start_y"]))
                        if block_radius <= 1e-6:
                            candidate_valid = False
                            break
                        block_angle = math.atan2(float(block_label["start_y"]), float(block_label["start_x"]))
                        candidate_angle = block_angle + signed_delta
                        candidate_x = block_radius * math.cos(candidate_angle)
                        candidate_y = block_radius * math.sin(candidate_angle)

                        block_sign = 1 if float(block_label["start_x"]) >= 0 else -1
                        if block_sign > 0 and candidate_x <= 0:
                            candidate_valid = False
                            break
                        if block_sign < 0 and candidate_x >= 0:
                            candidate_valid = False
                            break

                        candidate = block_label.copy()
                        _move_label_text_anchor(
                            candidate,
                            start_x=float(candidate_x),
                            start_y=float(candidate_y),
                        )
                        candidate_positions[block_idx] = (float(candidate_x), float(candidate_y))
                        candidate_labels[block_idx] = candidate

                    if not candidate_valid:
                        continue

                    for candidate in candidate_labels.values():
                        candidate_box = _label_bbox_on_canvas(
                            candidate, total_length, canvas_config, margin_px=LEGEND_LABEL_MARGIN_PX
                        )
                        if _bbox_overlaps(candidate_box, legend_box):
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, candidate in candidate_labels.items():
                        for peer_idx, peer_label in enumerate(labels):
                            if peer_idx == block_idx or peer_label.get("is_embedded"):
                                continue
                            peer = candidate_labels.get(peer_idx, peer_label)
                            min_gap_px = minimum_bbox_gap_px(candidate, peer)
                            if y_overlap(candidate, peer, total_length, min_gap_px) and x_overlap(
                                candidate, peer, minimum_margin=min_gap_px
                            ):
                                candidate_valid = False
                                break
                        if not candidate_valid:
                            break
                    if not candidate_valid:
                        continue

                    def _unwrapped_for_order(check_idx: int) -> float:
                        candidate = candidate_labels.get(check_idx)
                        if candidate is not None:
                            return _label_unwrapped_angle(candidate, total_length)
                        return _label_unwrapped_angle(labels[check_idx], total_length)

                    for order_pos in range(1, len(side_indices)):
                        prev_idx = side_indices[order_pos - 1]
                        curr_idx = side_indices[order_pos]
                        if _unwrapped_for_order(curr_idx) <= _unwrapped_for_order(prev_idx) + MIN_LABEL_ORDER_GAP_RAD:
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, (candidate_x, candidate_y) in candidate_positions.items():
                        _move_label_text_anchor(
                            labels[block_idx],
                            start_x=candidate_x,
                            start_y=candidate_y,
                        )
                    changed = True
                    placed = True
                    break
                shift_px += LABEL_NUDGE_STEP_PX

        if not changed:
            break

    return not _labels_collide_with_legend(
        labels,
        total_length,
        canvas_config,
        legend_bounds,
    )


def _resolve_label_legend_collisions(
    labels: list[dict[str, Any]],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> None:
    """Shift horizontal labels only for the four corner-overlay placements."""
    try:
        placement = LegendPlacement(str(canvas_config.legend_position))
    except ValueError:
        return
    if not placement.is_overlay:
        return
    if (
        legend_bounds.width <= 0
        or legend_bounds.height <= 0
    ):
        return

    external_labels = [label for label in labels if not label.get("is_embedded")]
    if not external_labels:
        return
    if not _labels_collide_with_legend(
        external_labels,
        total_length,
        canvas_config,
        legend_bounds,
    ):
        return

    if any(label.get("placement") == "radial" for label in external_labels):
        return
    _try_shift_labels_away_from_legend(
        external_labels,
        total_length,
        canvas_config,
        legend_bounds,
    )


def _slot_width_ratio_factor(
    slot: CircularTrackSlot | None,
    *,
    base_radius_px: float,
    base_track_ratio: float,
) -> float | None:
    if slot is None or slot.width is None:
        return None
    denominator = float(base_radius_px) * float(base_track_ratio)
    if denominator <= 0:
        return None
    return float(slot.width.resolve(float(base_radius_px))) / denominator




def _slot_dinucleotide(slot_or_resolved: CircularTrackSlot | CircularResolvedSlot, default: str) -> str:
    params = getattr(slot_or_resolved, "params", {}) or {}
    raw = params.get("nt", params.get("dinucleotide", default))
    nt = str(raw or default).upper()
    return nt if len(nt) >= 2 else str(default or "GC").upper()




def _text_element_plain_text(element: Any) -> str:
    parts: list[str] = []
    text = getattr(element, "text", None)
    if text:
        parts.append(str(text))
    for child in getattr(element, "elements", ()):
        child_text = getattr(child, "text", None)
        if child_text:
            parts.append(str(child_text))
    return "".join(parts)


def _definition_reserved_radius_px(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    species: str | None,
    strain: str | None,
    *,
    plot_title: str | None,
    definition_profile: str,
) -> float:
    cfg = canvas_config.profile.config
    definition_group = DefinitionGroup(
        gb_record,
        canvas_config,
        cfg=cfg,
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_profile=definition_profile,
    ).get_group()

    # Track-slot reservation is a record-local compatibility contract. Keep its
    # historical unstyled text metric independent from the exact painted bounds
    # used by outer composition, so composition metadata cannot reflow tracks.
    max_extent = 0.0
    for element in getattr(definition_group, "elements", ()):
        attribs = getattr(element, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        x_value = _svg_number(attribs.get("x"), default=0.0)
        y_value = _svg_number(attribs.get("y"), default=0.0)
        font_size = _svg_number(
            attribs.get("font-size"),
            default=cfg.objects.definition.circular.font_size,
        )
        text = _text_element_plain_text(element)
        if not text.strip():
            continue
        width_px, height_px = calculate_bbox_dimensions(
            text,
            str(attribs.get("font-family", cfg.objects.text.font_family)),
            font_size,
            int(canvas_config.dpi),
        )
        half_width = 0.5 * float(width_px)
        half_height = 0.5 * float(height_px if height_px else font_size)
        max_extent = max(
            max_extent,
            abs(x_value - half_width),
            abs(x_value + half_width),
            abs(y_value - half_height),
            abs(y_value + half_height),
        )

    if max_extent <= FEATURE_BAND_EPSILON:
        return 0.0
    padding_px = max(8.0, 0.02 * float(canvas_config.radius))
    return float(max_extent + padding_px)


def _slot_config_with_dinucleotide(config: Any, nt: str) -> Any:
    if str(getattr(config, "dinucleotide", "")).upper() == str(nt).upper():
        return config
    cloned = copy.copy(config)
    cloned.dinucleotide = str(nt).upper()
    return cloned






def _slot_gc_config_with_axis_font_size(
    config: GcContentConfigurator,
    nt: str,
    tick_font_size: float | None,
) -> GcContentConfigurator:
    slot_config = _slot_config_with_dinucleotide(config, nt)
    if tick_font_size is None or getattr(slot_config, "tick_font_size", None) is not None:
        return slot_config
    cloned = copy.copy(slot_config)
    cloned.tick_font_size = float(tick_font_size)
    return cloned


def _gc_content_matches_depth_axis_font_size(
    *,
    gc_config: GcContentConfigurator,
    depth_config: DepthConfigurator | None,
    resolved_track_slots: Sequence[CircularResolvedSlot],
) -> float | None:
    if depth_config is None:
        return None
    if str(getattr(gc_config, "mode", "deviation")).strip().lower() != "percent":
        return None
    if getattr(gc_config, "tick_font_size", None) is not None:
        return None
    if not bool(getattr(depth_config, "show_axis", True)) or not bool(getattr(depth_config, "show_ticks", True)):
        return None
    depth_slot = next(
        (
            slot
            for slot in resolved_track_slots
            if str(slot.renderer) == "depth" and float(slot.draw_width_px) > 0
        ),
        None,
    )
    if depth_slot is None:
        return None
    return depth_axis_tick_font_size_px(depth_config, float(depth_slot.draw_width_px))


def _slot_dataframe_for_nt(
    *,
    nt: str,
    default_df: DataFrame,
    default_nt: str,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
) -> DataFrame:
    normalized_nt = str(nt).upper()
    if dinucleotide_dataframes and normalized_nt in dinucleotide_dataframes:
        return dinucleotide_dataframes[normalized_nt]
    if normalized_nt == str(default_nt).upper():
        return default_df
    return DataFrame()




def _slot_legend_label(slot: CircularTrackSlot, fallback: str) -> str:
    params = slot.params or {}
    raw_label = params.get("legend_label", params.get("label"))
    label = str(raw_label).strip() if raw_label is not None else ""
    return label or fallback


def _sync_legend_table_for_circular_slots(
    legend_table: dict,
    *,
    circular_track_slots: list[CircularTrackSlot] | None,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
    depth_df: DataFrame | None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    cfg: GbdrawConfig,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float | None = None,
) -> dict:
    """Replace singleton numeric legend entries with slot-aware entries."""
    if circular_track_slots is None:
        return legend_table

    out = dict(legend_table)
    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    removable = {
        "Depth",
        f"{default_nt} content",
        f"{default_nt} content (+)",
        f"{default_nt} content (-)",
        f"{default_nt} skew",
        f"{default_nt} skew (+)",
        f"{default_nt} skew (-)",
    }
    for key in removable:
        out.pop(key, None)

    depth_by_index = index_depth_track_row(depth_tracks or ())
    for slot in circular_track_slots:
        if not slot.enabled:
            continue
        renderer = str(slot.renderer)
        if renderer == "depth":
            if depth_by_index:
                track_index = parse_nonnegative_integer(
                    (slot.params or {}).get("track_index", 0),
                    field_name=f"depth slot '{slot.id}' track_index",
                )
                track = depth_by_index.get(track_index)
                if track is None:
                    continue
                label = _slot_legend_label(slot, track.label)
                out[_unique_legend_key(out, label)] = {
                    "type": "solid",
                    "fill": track.config.fill_color,
                    "stroke": track.config.stroke_color,
                    "width": track.config.stroke_width,
                }
                continue
            if depth_config is None or depth_df is None:
                continue
            label = _slot_legend_label(slot, "Depth")
            out[_unique_legend_key(out, label)] = {
                "type": "solid",
                "fill": depth_config.fill_color,
                "stroke": depth_config.stroke_color,
                "width": depth_config.stroke_width,
            }
        elif renderer == "dinucleotide_content":
            nt = _slot_dinucleotide(slot, default_nt)
            label = _slot_legend_label(slot, f"{nt} content")
            if gc_config.high_fill_color == gc_config.low_fill_color:
                out[_unique_legend_key(out, label)] = {
                    "type": "solid",
                    "fill": gc_config.high_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
            else:
                out[_unique_legend_key(out, f"{label} (+)")] = {
                    "type": "solid",
                    "fill": gc_config.high_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
                out[_unique_legend_key(out, f"{label} (-)")] = {
                    "type": "solid",
                    "fill": gc_config.low_fill_color,
                    "stroke": gc_config.stroke_color,
                    "width": gc_config.stroke_width,
                }
        elif renderer == "dinucleotide_skew":
            nt = _slot_dinucleotide(slot, default_nt)
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
    if conservation_tracks:
        if any(track.track_color for track in conservation_tracks):
            for track in conservation_tracks:
                min_color, max_color = conservation_track_gradient_colors(
                    track.track_color,
                    default_min_color=cfg.objects.conservation.min_color,
                    default_max_color=cfg.objects.conservation.max_color,
                )
                out[_unique_legend_key(out, track.track_label)] = {
                    "type": "gradient",
                    "min_color": min_color,
                    "max_color": max_color,
                    "stroke": "none",
                    "width": 0,
                    "min_value": float(conservation_min_identity or 0.0),
                }
        else:
            out[_unique_legend_key(out, "Conservation identity")] = {
                "type": "gradient",
                "min_color": cfg.objects.conservation.min_color,
                "max_color": cfg.objects.conservation.max_color,
                "stroke": "none",
                "width": 0,
                "min_value": float(conservation_min_identity or 0.0),
            }
    return out


def _draw_resolved_circular_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
    conservation_tracks: Sequence[ConservationTrack] | None,
    conservation_min_identity: float,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None = None,
    gc_content_tick_font_size_override: float | None = None,
    feature_layers: FeatureBuildResult,
    render_context: CircularRecordRenderContext,
    precalculated_labels: list[dict] | None = None,
    _tick_track_channel_override: str | None = None,
    use_slot_group_id: bool = True,
    use_slot_tick_options: bool = True,
    use_feature_anchor_override: bool = True,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack] | None = None,
    annotation_record_index: int = 0,
) -> Drawing:
    """Draw one resolved circular slot."""
    cfg = render_context.profile.config
    renderer = str(resolved_slot.renderer)
    norm_factor_override = float(resolved_slot.anchor_radius_px) / float(canvas_config.radius)
    if renderer == "spacer":
        return canvas

    if renderer == "annotations":
        annotation_layout = (annotation_track_layouts or {}).get(resolved_slot.id)
        if annotation_layout is None:
            return canvas
        params = annotation_track_params_from_mapping(resolved_slot.params)
        group = draw_circular_annotation_track(
            canvas,
            annotation_layout,
            record_id=str(gb_record.id),
            record_index=annotation_record_index,
            record_length=len(gb_record.seq),
            inner_radius_px=float(resolved_slot.draw_inner_radius_px),
            outer_radius_px=float(resolved_slot.draw_outer_radius_px),
            side=str(resolved_slot.side),
            font_family=str(cfg.objects.text.font_family),
            params=params,
            dom_group_id=(
                track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                )
                if use_slot_group_id
                else None
            ),
            semantic_slot_id=(
                str(resolved_slot.id) if use_slot_group_id else None
            ),
        )
        group = center_group_on_canvas(group, canvas_config)
        canvas.add(group)
        return canvas

    if renderer == "features":
        base_width = (
            float(canvas_config.radius)
            * float(canvas_config.track_ratio)
            * float(cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)][0])
        )
        resolved_feature_width = float(resolved_slot.draw_width_px)
        if render_context.feature_layout is not None:
            resolved_feature_width = float(render_context.feature_layout.width_px)
        ratio_override = None
        if base_width > 0 and resolved_feature_width > 0:
            ratio_override = resolved_feature_width / max(
                FEATURE_BAND_EPSILON,
                float(canvas_config.radius) * float(canvas_config.track_ratio),
            )
        feature_kwargs: dict[str, Any] = {
            "feature_layers": feature_layers,
            "render_context": render_context,
            "precalculated_labels": precalculated_labels,
            "feature_track_ratio_factor_override": ratio_override,
        }
        if annotation_record_index:
            feature_kwargs["record_index"] = annotation_record_index
        if use_slot_group_id:
            feature_kwargs["group_id"] = track_slot_svg_id(
                resolved_slot.id,
                renderer=renderer,
                slot_index=int(resolved_slot.slot_index),
            )
            feature_kwargs["slot_id"] = str(resolved_slot.id)
            feature_kwargs["feature_dom_namespace"] = (
                f"{int(resolved_slot.slot_index) + 1}_{resolved_slot.id}"
            )
        if use_feature_anchor_override or not math.isclose(
            float(resolved_slot.anchor_radius_px),
            float(canvas_config.radius),
            rel_tol=1e-9,
            abs_tol=1e-9,
        ):
            feature_kwargs["feature_anchor_radius_px"] = float(resolved_slot.anchor_radius_px)
        return add_record_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            **feature_kwargs,
        )

    if renderer == "ticks":
        tick_layout = (
            resolved_slot.payload
            if isinstance(resolved_slot.payload, CircularTickLayout)
            else None
        )
        label_side = (
            str(tick_layout.label_side)
            if tick_layout is not None
            else tick_sides_for_tick_label_layout(
                resolved_slot.params.get("tick_label_layout"),
                side=resolved_slot.side,
            )[0]
        ).strip().lower()
        tick_side = (
            str(tick_layout.tick_side)
            if tick_layout is not None
            else tick_sides_for_tick_label_layout(
                resolved_slot.params.get("tick_label_layout"),
                side=resolved_slot.side,
            )[1]
        ).strip().lower()
        if label_side != "none" or tick_side != "none":
            tick_group_kwargs: dict[str, Any] = {
                "radius_override": float(resolved_slot.anchor_radius_px),
            }
            if use_slot_tick_options or tick_layout is not None:
                tick_group_kwargs["label_side"] = label_side
                tick_group_kwargs["tick_side"] = tick_side
                tick_group_kwargs["tick_length_px"] = (
                    tick_layout.tick_length_px
                    if tick_layout is not None
                    else (
                        float(resolved_slot.draw_width_px)
                        if resolved_slot.explicit_width and resolved_slot.draw_width_px > 0
                        else None
                    )
                )
                tick_group_kwargs["track_preset"] = (
                    tick_layout.track_preset
                    if tick_layout is not None
                    else normalize_circular_track_preset(
                        str(resolved_slot.params.get("preset", resolved_slot.params.get("track_preset", "tuckin")))
                    )
                )
            if _tick_track_channel_override is not None:
                tick_group_kwargs["tick_track_channel_override"] = _tick_track_channel_override
            tick_group_kwargs["slot_id"] = str(resolved_slot.id)
            if use_slot_group_id:
                tick_group_kwargs["group_id"] = track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                )
            canvas = add_tick_group_on_canvas(
                canvas,
                gb_record,
                canvas_config,
                **tick_group_kwargs,
            )
        return canvas

    if renderer == "depth":
        depth_group_id = str(resolved_slot.id)
        selected_depth_df = depth_df
        selected_depth_config = depth_config
        if depth_tracks_by_index is not None:
            track_index = parse_nonnegative_integer(
                resolved_slot.params.get("track_index", 0),
                field_name=f"depth slot '{resolved_slot.id}' track_index",
            )
            depth_track = depth_tracks_by_index.get(track_index)
            if depth_track is None:
                return canvas
            selected_depth_df = depth_track.df
            selected_depth_config = depth_track.config
            if not use_slot_group_id:
                depth_group_id = depth_track.id
        if selected_depth_config is None or selected_depth_df is None:
            logger.warning("Skipping circular depth slot '%s' because depth data are unavailable.", resolved_slot.id)
            return canvas
        depth_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
        }
        if use_slot_group_id:
            depth_kwargs["group_id"] = depth_group_id
        canvas = add_depth_group_on_canvas(
            canvas,
            gb_record,
            selected_depth_df,
            canvas_config,
            selected_depth_config,
            **depth_kwargs,
        )
        if use_slot_group_id:
            canvas = _mark_circular_track_slot_group(
                canvas,
                source_group_id=depth_group_id,
                group_id=track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                ),
                slot_id=str(resolved_slot.id),
                renderer=renderer,
            )
        return canvas

    if renderer == "sequence_conservation":
        track_index = int(resolved_slot.params.get("track_index", 0) or 0)
        track = next(
            (
                conservation_track
                for conservation_track in (conservation_tracks or ())
                if int(conservation_track.track_index) == track_index
            ),
            None,
        )
        if track is None:
            logger.warning(
                "Skipping circular conservation slot '%s' because conservation data are unavailable.",
                resolved_slot.id,
            )
            return canvas
        conservation_kwargs: dict[str, Any] = {}
        if use_slot_group_id:
            conservation_kwargs["group_id"] = track_slot_svg_id(
                resolved_slot.id,
                renderer=renderer,
                slot_index=int(resolved_slot.slot_index),
            )
            conservation_kwargs["slot_id"] = str(resolved_slot.id)
        return add_conservation_group_on_canvas(
            canvas,
            gb_record,
            track,
            canvas_config,
            inner_radius_px=float(resolved_slot.draw_inner_radius_px),
            outer_radius_px=float(resolved_slot.draw_outer_radius_px),
            min_identity=float(conservation_min_identity),
            **conservation_kwargs,
        )

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = _slot_dinucleotide(resolved_slot, default_nt)
    if renderer == "dinucleotide_content":
        slot_df = _slot_dataframe_for_nt(
            nt=nt,
            default_df=gc_df,
            default_nt=default_nt,
            dinucleotide_dataframes=dinucleotide_content_dataframes or dinucleotide_dataframes,
        )
        if f"{nt} content" not in slot_df.columns:
            logger.warning(
                "Skipping circular content slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        gc_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
        }
        if use_slot_group_id:
            gc_kwargs["group_id"] = str(resolved_slot.id)
        slot_gc_config = _slot_gc_config_with_axis_font_size(
            gc_config,
            nt,
            gc_content_tick_font_size_override,
        )
        canvas = add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_gc_config,
            **gc_kwargs,
        )
        if use_slot_group_id:
            canvas = _mark_circular_track_slot_group(
                canvas,
                source_group_id=str(gc_kwargs["group_id"]),
                group_id=track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                ),
                slot_id=str(resolved_slot.id),
                renderer=renderer,
            )
        else:
            canvas = _tag_circular_track_slot_group(
                canvas,
                source_group_id="gc_content",
                slot_id=str(resolved_slot.id),
                renderer=renderer,
            )
        return canvas

    if renderer == "dinucleotide_skew":
        slot_df = _slot_dataframe_for_nt(
            nt=nt,
            default_df=gc_df,
            default_nt=default_nt,
            dinucleotide_dataframes=dinucleotide_skew_dataframes or dinucleotide_dataframes,
        )
        if f"{nt} skew" not in slot_df.columns:
            logger.warning(
                "Skipping circular skew slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        skew_kwargs: dict[str, Any] = {
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
        }
        if use_slot_group_id:
            skew_kwargs["group_id"] = str(resolved_slot.id)
        canvas = add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            _slot_skew_config(skew_config, resolved_slot, nt),
            **skew_kwargs,
        )
        if use_slot_group_id:
            canvas = _mark_circular_track_slot_group(
                canvas,
                source_group_id=str(skew_kwargs["group_id"]),
                group_id=track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                ),
                slot_id=str(resolved_slot.id),
                renderer=renderer,
            )
        else:
            canvas = _tag_circular_track_slot_group(
                canvas,
                source_group_id="skew",
                slot_id=str(resolved_slot.id),
                renderer=renderer,
            )
        return canvas

    logger.warning("Skipping unsupported circular track slot renderer '%s'.", renderer)
    return canvas


def _default_outer_label_arena(
    *,
    canvas_config: CircularCanvasConfigurator,
) -> tuple[float, float]:
    """Return default (anchor_radius_px, arc_outer_radius_px) for outer labels."""
    profile = canvas_config.profile
    cfg = profile.config
    length_param = str(canvas_config.length_param)
    track_type = _circular_preset_for_layout(cfg)
    strands = "separate" if profile.strandedness else "single"
    base_radius = float(canvas_config.radius)

    anchor_radius = float(cfg.labels.radius_factor[track_type][strands][length_param]) * base_radius
    offset_cfg = cfg.labels.unified_adjustment.outer_labels
    arc_x = (
        base_radius
        * float(cfg.labels.arc_x_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.x_radius_offset)
    )
    arc_y = (
        base_radius
        * float(cfg.labels.arc_y_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.y_radius_offset)
    )
    arc_outer = max(abs(arc_x), abs(arc_y))

    if arc_outer < anchor_radius:
        arc_outer = anchor_radius
    return anchor_radius, arc_outer




def add_record_on_circular_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    gc_df: DataFrame,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    legend_measurement: LegendMeasurement,
    legend_table,
    *,
    legend_config: LegendDrawingConfigurator | None = None,
    depth_config: DepthConfigurator | None = None,
    depth_df: DataFrame | None = None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    depth_track_count_value: int | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    _tick_track_channel_override: str | None = None,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None,
    annotation_record_index: int = 0,
    definition_record_count: int = 1,
) -> _CircularPlotAssembly:
    """
    Adds various record-related groups to a circular canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    species (str): Species name.
    strain (str): Strain name.
    Returns:
    Drawing: The updated SVG drawing with all record-related groups added.
    """
    primary_target_start = len(getattr(canvas, "elements", []))
    profile = canvas_config.profile
    cfg = profile.config
    depth_tracks_by_index = (
        index_depth_track_row(depth_tracks)
        if depth_tracks is not None
        else None
    )
    resolved_depth_track_count = (
        int(depth_track_count_value)
        if depth_track_count_value is not None
        else depth_track_data_count([depth_tracks or ()])
    )

    show_labels_base = profile.labels_enabled
    depth_enabled = bool(profile.show_depth and depth_config is not None)
    simple_show_ticks = bool(cfg.objects.scale.show)
    effective_circular_track_slots, resolved_annotations, annotation_track_layouts = _prepare_circular_annotation_tracks(
        gb_record,
        annotations,
        circular_track_slots,
        canvas_config=canvas_config,
        record_index=annotation_record_index,
        show_ticks=simple_show_ticks,
        show_depth=depth_enabled,
        show_gc=profile.show_gc,
        show_skew=profile.show_skew,
        depth_track_count=max(1, resolved_depth_track_count),
    )
    user_slot_mode = effective_circular_track_slots is not None
    user_active_slot_renderers = {
        str(slot.renderer)
        for slot in (effective_circular_track_slots or [])
        if slot.enabled
    }
    if user_slot_mode:
        show_depth_track = bool("depth" in user_active_slot_renderers and depth_enabled)
        show_gc_track = bool("dinucleotide_content" in user_active_slot_renderers)
        show_skew_track = bool("dinucleotide_skew" in user_active_slot_renderers)
        show_ticks_track = bool("ticks" in user_active_slot_renderers)
    else:
        show_depth_track = bool(depth_enabled)
        show_gc_track = profile.show_gc
        show_skew_track = profile.show_skew
        show_ticks_track = simple_show_ticks
    circular_preset = normalize_circular_track_preset(cfg.canvas.circular.track_type)

    if user_slot_mode:
        show_features = "features" in user_active_slot_renderers
        radial_plan = circular_track_slots_from_preset_order(
            list(effective_circular_track_slots or []),
            circular_preset,
            CircularPresetContext(
                cfg=cfg,
                canvas_config=canvas_config,
                total_length=len(gb_record.seq),
                strandedness=profile.strandedness,
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
                tick_track_channel_override=_tick_track_channel_override,
            ),
            axis_index=circular_track_axis_index,
        )
        layout_slots = list(radial_plan.slots)
        preferred_anchor_slot_ids = radial_plan.preferred_anchor_slot_ids
    else:
        show_features = True
        radial_plan = circular_radial_plan_for_preset(
            circular_preset,
            CircularPresetContext(
                cfg=cfg,
                canvas_config=canvas_config,
                total_length=len(gb_record.seq),
                strandedness=profile.strandedness,
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
                tick_track_channel_override=_tick_track_channel_override,
            ),
        )
        layout_slots = list(radial_plan.slots)
        preferred_anchor_slot_ids = radial_plan.preferred_anchor_slot_ids
    feature_slot = next(
        (
            slot
            for slot in layout_slots
            if slot.enabled and str(slot.renderer) == "features"
        ),
        None,
    )
    feature_lane_direction = _lane_direction_for_feature_slot(
        feature_slot,
        track_preset=circular_preset,
    )

    show_external_labels = show_labels_base and show_features
    split_overlaps_by_strand = (
        profile.resolve_overlaps
        and (not profile.strandedness)
        and feature_lane_direction == "split"
    )

    feature_track_ratio_factor_override: float | None = None
    if show_features:
        feature_track_ratio_factor_override = _slot_width_ratio_factor(
            feature_slot,
            base_radius_px=float(canvas_config.radius),
            base_track_ratio=float(canvas_config.track_ratio),
        )
    feature_width_override_requested = feature_track_ratio_factor_override is not None

    compute_label_text = show_external_labels
    label_filtering = (
        preprocess_label_filtering(cfg.labels.filtering.as_dict())
        if compute_label_text
        else {}
    )
    feature_layers = create_feature_layers(
        gb_record,
        feature_config.specific_color_rules,
        feature_config.selected_features_set,
        feature_config.default_color_map,
        profile.strandedness,
        profile.resolve_overlaps,
        label_filtering,
        split_overlaps_by_strand=split_overlaps_by_strand,
        feature_shapes=feature_config.feature_shapes,
        feature_visibility_rules=feature_config.feature_visibility_rules,
        compute_label_text=compute_label_text,
        feature_instance_identity_plan=(
            feature_config.feature_instance_identity_plan
        ),
        feature_semantic_selector_context=(
            feature_config.feature_semantic_selector_context
        ),
    )
    precomputed_feature_dict: dict = feature_layers.foreground_features
    layout_feature_dict = {
        f"foreground:{feature_id}": feature
        for feature_id, feature in feature_layers.foreground_features.items()
    }
    layout_feature_dict.update(
        {
            f"underlay:{feature_index}": feature
            for feature_index, feature in enumerate(feature_layers.underlay_features)
        }
    )
    if feature_layers.underlay_features:
        (
            layout_slots,
            resolved_annotations,
            annotation_track_layouts,
        ) = _add_circular_feature_underlays(
            gb_record,
            feature_layers.underlay_features,
            layout_slots,
            resolved_annotations,
            canvas_config=canvas_config,
            record_index=annotation_record_index,
            show_ticks=show_ticks_track,
            show_depth=show_depth_track,
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            depth_track_count=max(1, resolved_depth_track_count),
        )
    resolved_track_slots: list[CircularResolvedSlot] = []
    resolved_feature_anchor_radius_px: float | None = None
    definition_reserved_radius_px: float | None = None
    if center_reserved_radius is not None:
        definition_reserved_radius_px = max(0.0, float(center_reserved_radius))
    elif str(definition_position).strip().lower() == "center":
        definition_reserved_radius_px = _definition_reserved_radius_px(
            gb_record,
            canvas_config,
            species,
            strain,
            plot_title=plot_title,
            definition_profile=definition_profile,
        )
    radial_layout = resolve_circular_radial_layout(
        total_length=len(gb_record.seq),
        canvas_config=canvas_config,
        slots=layout_slots,
        feature_dict=layout_feature_dict,
        definition_reserved_radius_px=definition_reserved_radius_px,
        tick_track_channel_override=_tick_track_channel_override,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
        depth_config=depth_config if show_depth_track else None,
    )
    resolved_track_slots = list(radial_layout.slots)
    setattr(
        canvas,
        "_gbdraw_track_slot_geometry",
        _serialize_circular_track_slot_geometry(
            gb_record=gb_record,
            radial_layout=radial_layout,
            layout_slots=layout_slots,
            base_radius_px=float(canvas_config.radius),
        ),
    )
    gc_content_tick_font_size_override = _gc_content_matches_depth_axis_font_size(
        gc_config=gc_config,
        depth_config=depth_config if show_depth_track else None,
        resolved_track_slots=resolved_track_slots,
    )
    if radial_layout.features is not None:
        resolved_feature_anchor_radius_px = float(radial_layout.features.anchor_radius_px)
        feature_track_ratio_factor_override = (
            float(radial_layout.features.width_px)
            / max(FEATURE_BAND_EPSILON, float(canvas_config.radius) * float(canvas_config.track_ratio))
        )
    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features.
    outer_arena: tuple[float, float] | None = None
    if show_external_labels:
        default_anchor_px, default_arc_outer_px = _default_outer_label_arena(
            canvas_config=canvas_config,
        )
        if default_arc_outer_px < default_anchor_px:
            default_anchor_px, default_arc_outer_px = default_arc_outer_px, default_anchor_px

        if feature_width_override_requested:
            outer_arena = (default_anchor_px, default_arc_outer_px)

    precalculated_labels: list[dict[str, Any]] | None = None
    if show_external_labels and precomputed_feature_dict is not None:
        label_candidate_cache: dict[str, object] = {}
        precalculated_labels = prepare_label_list(
            precomputed_feature_dict,
            len(gb_record.seq),
            (
                float(resolved_feature_anchor_radius_px)
                if resolved_feature_anchor_radius_px is not None
                else canvas_config.radius
            ),
            canvas_config.track_ratio,
            profile,
            outer_arena=outer_arena,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            feature_layout=radial_layout.features,
            radial_layout=radial_layout,
            track_preset=_circular_preset_for_layout(cfg),
            feature_lane_direction=feature_lane_direction,
            _candidate_cache=label_candidate_cache,
        )
        if profile.label_placement == "radial":
            previous_growth: float | None = None
            stagnant_passes = 0
            preflight_tracks_frozen = False
            while True:
                required_growth = max(
                    (
                        float(label.get("required_radius_growth_px", 0.0))
                        for label in precalculated_labels
                        if not label.get("is_embedded")
                    ),
                    default=0.0,
                )
                if required_growth <= 1e-3:
                    break
                if previous_growth is not None and required_growth >= previous_growth - 1e-3:
                    stagnant_passes += 1
                else:
                    stagnant_passes = 0
                if stagnant_passes >= 2:
                    fixed_inside_slots = [
                        slot.id
                        for slot in radial_layout.slots
                        if slot.side == "inside" and (slot.explicit_anchor or slot.explicit_width)
                    ]
                    slot_context = ", ".join(fixed_inside_slots) or "center/feature interval"
                    raise ValidationError(
                        "radial inner labels cannot fit the fixed circular geometry: "
                        f"slot={slot_context}, required_additional_px={required_growth:.3f}"
                    )
                previous_growth = required_growth

                if not preflight_tracks_frozen:
                    frozen_slots: list[CircularTrackSlot] = []
                    for slot_index, slot in enumerate(layout_slots):
                        resolved_slot = next(
                            (
                                candidate
                                for candidate in radial_layout.slots
                                if int(candidate.slot_index) == int(slot_index)
                            ),
                            None,
                        )
                        if (
                            str(slot.renderer) != "features"
                            and resolved_slot is not None
                            and resolved_slot.anchor_radius_px is not None
                        ):
                            frozen_slots.append(
                                replace(
                                    slot,
                                    radius=ScalarSpec(
                                        float(resolved_slot.anchor_radius_px),
                                        unit="px",
                                    ),
                                )
                            )
                        elif str(slot.renderer) == "features" and slot.radius is None and slot.side == "inside":
                            # Radial inner labels need a real arena between features and
                            # the frozen inner tracks. Move the automatic tuck-in feature
                            # slot to the outside of the enlarged axis while preserving
                            # its lane-direction payload.
                            feature_params = dict(slot.params)
                            feature_params["lane_direction"] = "outside"
                            frozen_slots.append(
                                replace(slot, side="outside", params=feature_params)
                            )
                        else:
                            frozen_slots.append(slot)
                    layout_slots = frozen_slots
                    feature_slot = next(
                        (
                            slot
                            for slot in layout_slots
                            if slot.enabled and str(slot.renderer) == "features"
                        ),
                        None,
                    )
                    if feature_slot is not None:
                        feature_lane_direction = _lane_direction_for_feature_slot(
                            feature_slot,
                            track_preset=circular_preset,
                        )
                    preflight_tracks_frozen = True

                canvas_config.radius = float(canvas_config.radius) + required_growth
                feature_track_ratio_factor_override = _slot_width_ratio_factor(
                    feature_slot,
                    base_radius_px=float(canvas_config.radius),
                    base_track_ratio=float(canvas_config.track_ratio),
                )
                radial_layout = resolve_circular_radial_layout(
                    total_length=len(gb_record.seq),
                    canvas_config=canvas_config,
                    slots=layout_slots,
                    feature_dict=layout_feature_dict,
                    definition_reserved_radius_px=definition_reserved_radius_px,
                    tick_track_channel_override=_tick_track_channel_override,
                    preferred_anchor_slot_ids=preferred_anchor_slot_ids,
                    depth_config=depth_config if show_depth_track else None,
                )
                resolved_track_slots = list(radial_layout.slots)
                setattr(
                    canvas,
                    "_gbdraw_track_slot_geometry",
                    _serialize_circular_track_slot_geometry(
                        gb_record=gb_record,
                        radial_layout=radial_layout,
                        layout_slots=layout_slots,
                        base_radius_px=float(canvas_config.radius),
                    ),
                )
                if radial_layout.features is not None:
                    resolved_feature_anchor_radius_px = float(radial_layout.features.anchor_radius_px)
                    feature_track_ratio_factor_override = (
                        float(radial_layout.features.width_px)
                        / max(
                            FEATURE_BAND_EPSILON,
                            float(canvas_config.radius) * float(canvas_config.track_ratio),
                        )
                    )
                default_anchor_px, default_arc_outer_px = _default_outer_label_arena(
                    canvas_config=canvas_config,
                )
                if default_arc_outer_px < default_anchor_px:
                    default_anchor_px, default_arc_outer_px = default_arc_outer_px, default_anchor_px
                outer_arena = (
                    (default_anchor_px, default_arc_outer_px)
                    if feature_width_override_requested
                    else None
                )
                precalculated_labels = prepare_label_list(
                    precomputed_feature_dict,
                    len(gb_record.seq),
                    (
                        float(resolved_feature_anchor_radius_px)
                        if resolved_feature_anchor_radius_px is not None
                        else canvas_config.radius
                    ),
                    canvas_config.track_ratio,
                    profile,
                    outer_arena=outer_arena,
                    feature_track_ratio_factor_override=feature_track_ratio_factor_override,
                    feature_layout=radial_layout.features,
                    radial_layout=radial_layout,
                    track_preset=_circular_preset_for_layout(cfg),
                    feature_lane_direction=feature_lane_direction,
                    _candidate_cache=label_candidate_cache,
                )

            gc_content_tick_font_size_override = _gc_content_matches_depth_axis_font_size(
                gc_config=gc_config,
                depth_config=depth_config if show_depth_track else None,
                resolved_track_slots=resolved_track_slots,
            )
        if profile.label_placement == "horizontal":
            assign_leader_start_points(
                [label for label in precalculated_labels if not label.get("is_embedded")],
                len(gb_record.seq),
            )

    overlay_obstacles = _external_label_obstacles_on_canvas(
        precalculated_labels or (),
        len(gb_record.seq),
        canvas_config,
    )
    source_content_bounds = _collect_circular_primary_bounds(
        radial_layout=radial_layout,
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
        overlay_obstacles=overlay_obstacles,
    )
    if legend_config is not None:
        legend_measurement = legend_config.measure_legend(
            legend_table,
            placement=canvas_config.legend_position,
            wrap_width=source_content_bounds.width,
        )
    legend_placement = LegendPlacement(str(canvas_config.legend_position))
    if (
        legend_placement.is_overlay
        and legend_measurement.local_bounds.width > 0
        and legend_measurement.local_bounds.height > 0
    ):
        legend_collision_bounds = _overlay_legend_bounds_for_collision(
            legend_measurement=legend_measurement,
            primary_bounds=source_content_bounds,
            placement=legend_placement,
        )
        if precalculated_labels:
            _resolve_label_legend_collisions(
                precalculated_labels,
                len(gb_record.seq),
                canvas_config,
                legend_collision_bounds,
            )
            if profile.label_placement == "horizontal":
                assign_leader_start_points(
                    [
                        label
                        for label in precalculated_labels
                        if not label.get("is_embedded")
                    ],
                    len(gb_record.seq),
                )
            overlay_obstacles = _external_label_obstacles_on_canvas(
                precalculated_labels,
                len(gb_record.seq),
                canvas_config,
            )
            source_content_bounds = _collect_circular_primary_bounds(
                radial_layout=radial_layout,
                canvas_config=canvas_config,
                feature_config=feature_config,
                gc_config=gc_config,
                skew_config=skew_config,
                depth_config=depth_config,
                overlay_obstacles=overlay_obstacles,
            )

    render_context = CircularRecordRenderContext(
        profile=profile,
        track_preset=circular_preset,
        feature_lane_direction=feature_lane_direction,
        radial_layout=radial_layout,
    )

    canvas = add_axis_group_on_canvas(
        canvas,
        canvas_config,
        radius_override=None,
    )

    if show_external_labels:
        labels_group_kwargs: dict[str, Any] = {
            "outer_arena": outer_arena,
            "feature_layers": feature_layers,
            "render_context": render_context,
            "precalculated_labels": precalculated_labels,
            "feature_track_ratio_factor_override": feature_track_ratio_factor_override,
        }
        if resolved_feature_anchor_radius_px is not None and (
            precalculated_labels is None
            and (
                user_slot_mode
                or not math.isclose(
                    float(resolved_feature_anchor_radius_px),
                    float(canvas_config.radius),
                    rel_tol=1e-9,
                    abs_tol=1e-9,
                )
            )
        ):
            labels_group_kwargs["feature_anchor_radius_px"] = resolved_feature_anchor_radius_px
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            phase="leaders",
            **labels_group_kwargs,
        )

    drawable_slots = [
        slot for slot in resolved_track_slots
        if str(slot.renderer) != "feature_labels"
    ]
    for resolved_slot in sorted(drawable_slots, key=lambda item: (int(item.z), int(item.slot_index))):
        canvas = _draw_resolved_circular_slot(
            canvas,
            resolved_slot,
            gb_record=gb_record,
            canvas_config=canvas_config,
            feature_config=feature_config,
            gc_df=gc_df,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_df=depth_df,
            depth_config=depth_config,
            depth_tracks_by_index=depth_tracks_by_index,
            conservation_tracks=conservation_tracks,
            conservation_min_identity=conservation_min_identity,
            dinucleotide_dataframes=dinucleotide_dataframes,
            dinucleotide_content_dataframes=dinucleotide_content_dataframes,
            dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
            gc_content_tick_font_size_override=gc_content_tick_font_size_override,
            feature_layers=feature_layers,
            render_context=render_context,
            precalculated_labels=precalculated_labels,
            _tick_track_channel_override=_tick_track_channel_override,
            use_slot_group_id=user_slot_mode,
            use_slot_tick_options=user_slot_mode,
            use_feature_anchor_override=user_slot_mode,
            annotation_track_layouts=annotation_track_layouts,
            annotation_record_index=annotation_record_index,
        )

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            phase="text",
            **labels_group_kwargs,
        )

    definition_kwargs: dict[str, Any] = {}
    definition_kwargs["record_index"] = int(annotation_record_index)
    definition_kwargs["record_count"] = int(definition_record_count)
    if plot_title is not None:
        definition_kwargs["plot_title"] = plot_title
    if str(definition_profile) != "full":
        definition_kwargs["definition_profile"] = definition_profile
    if str(definition_position) != "center":
        definition_kwargs["definition_position"] = definition_position
    if definition_group_id is not None:
        definition_kwargs["definition_group_id"] = definition_group_id
    definition_target_start = len(getattr(canvas, "elements", []))
    canvas = add_record_definition_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        species,
        strain,
        **definition_kwargs,
    )
    definition_target = next(
        (
            element
            for element in getattr(canvas, "elements", [])[definition_target_start:]
            if isinstance(element, Group)
        ),
        None,
    )
    source_content_bounds = _collect_circular_primary_bounds(
        radial_layout=radial_layout,
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
        overlay_obstacles=overlay_obstacles,
        definition_target=definition_target,
    )
    primary_targets = tuple(
        element
        for element in getattr(canvas, "elements", [])[primary_target_start:]
        if isinstance(element, Group)
    )

    legend_target: Group | None = None
    if (
        canvas_config.legend_position != "none"
        and legend_measurement.local_bounds.width > 0.0
        and legend_measurement.local_bounds.height > 0.0
    ):
        legend_target_start = len(getattr(canvas, "elements", []))
        canvas = add_legend_group_on_canvas(
            canvas,
            canvas_config,
            legend_measurement,
            legend_table,
        )
        legend_target = next(
            (
                element
                for element in getattr(canvas, "elements", [])[legend_target_start:]
                if isinstance(element, Group)
            ),
            None,
        )
    track_slot_geometry = getattr(canvas, "_gbdraw_track_slot_geometry", {})
    return _CircularPlotAssembly(
        drawing=canvas,
        canvas_config=canvas_config,
        legend_measurement=legend_measurement,
        legend_table=legend_table,
        source_content_bounds=source_content_bounds,
        source_overlay_obstacles=overlay_obstacles,
        primary_targets=primary_targets,
        legend_target=legend_target,
        track_slot_geometry=track_slot_geometry,
    )


def _compose_circular_plot(
    plot: _CircularPlotAssembly,
    *,
    title_target: Group | None = None,
    title_bounds: Aabb | None = None,
    title_placement: TitlePlacement | str = TitlePlacement.NONE,
) -> CircularAssemblyResult:
    """Compose one fully measured Circular plot exactly once."""
    resolved_title_placement = (
        title_placement
        if isinstance(title_placement, TitlePlacement)
        else TitlePlacement(str(title_placement))
    )
    if (title_target is None) != (title_bounds is None):
        raise ValueError("title_target and title_bounds must either both exist or be absent")
    if title_target is None and resolved_title_placement is not TitlePlacement.NONE:
        raise ValueError("title placement requires a title target and bounds")
    if title_target is not None and resolved_title_placement is TitlePlacement.NONE:
        raise ValueError("title target requires a non-none title placement")
    if title_target is not None and title_target not in getattr(plot.drawing, "elements", []):
        plot.drawing.add(title_target)

    legend_placement = LegendPlacement(str(plot.canvas_config.legend_position))
    legend_item = (
        CompositionItem("legend", plot.legend_measurement.local_bounds)
        if plot.legend_target is not None
        else None
    )
    title_item = (
        CompositionItem("title", title_bounds)
        if title_target is not None and title_bounds is not None
        else None
    )
    plan = plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", plot.source_content_bounds),
            legend=legend_item,
            title=title_item,
            legend_placement=legend_placement,
            title_placement=resolved_title_placement,
            overlay_obstacles=plot.source_overlay_obstacles,
        )
    )

    apply_composition_plan(
        plot.drawing,
        plan,
        primary_targets=plot.primary_targets,
        legend_target=plot.legend_target,
        legend_side=legend_placement,
        legend_reflow_metrics=(
            plot.legend_measurement.reflow_metrics
            if plot.legend_target is not None
            else None
        ),
        title_target=title_target,
        title_side=resolved_title_placement,
    )

    primary_placement = plan.placement_for("primary")
    if primary_placement is None:  # pragma: no cover - request always has primary
        raise RuntimeError("circular composition has no primary placement")
    final_obstacles = tuple(
        obstacle.translated(primary_placement.dx, primary_placement.dy)
        for obstacle in plot.source_overlay_obstacles
    )
    plot.canvas_config.total_width = plan.width
    plot.canvas_config.total_height = plan.height
    plot.canvas_config.offset_x = float(plot.canvas_config.offset_x) + primary_placement.dx
    plot.canvas_config.offset_y = float(plot.canvas_config.offset_y) + primary_placement.dy
    legend_plan_placement = plan.placement_for("legend")
    if legend_plan_placement is not None:
        plot.canvas_config.legend_offset_x = legend_plan_placement.dx
        plot.canvas_config.legend_offset_y = legend_plan_placement.dy

    result = CircularAssemblyResult(
        drawing=plot.drawing,
        legend_measurement=plot.legend_measurement,
        legend_table=plot.legend_table,
        source_content_bounds=plot.source_content_bounds,
        content_bounds=plan.primary_bounds,
        overlay_obstacles=final_obstacles,
        primary_targets=plot.primary_targets,
        legend_target=plot.legend_target,
        title_target=title_target,
        track_slot_geometry=plot.track_slot_geometry,
        composition_plan=plan,
    )
    setattr(plot.drawing, "_gbdraw_circular_assembly_result", result)
    return result


def _assemble_circular_diagram_result(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    feature_config: FeatureDrawingConfigurator,
    species: Optional[str],
    strain: Optional[str],
    plot_title: Optional[str],
    legend_config: LegendDrawingConfigurator,
    depth_df: DataFrame | None = None,
    depth_config: DepthConfigurator | None = None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    depth_track_count_value: int | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    _tick_track_channel_override: str | None = None,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None,
    annotation_record_index: int = 0,
    definition_record_count: int = 1,
    title_target: Group | None = None,
    title_bounds: Aabb | None = None,
    title_placement: TitlePlacement | str = TitlePlacement.NONE,
) -> CircularAssemblyResult:
    """
    Assembles a circular diagram for a GenBank record and returns the SVG canvas.

    Parameters:
    gb_record (SeqRecord): The GenBank record to plot.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    species (str): Species name.
    strain (str): Strain name.
    out_formats (list): List of formats to save the output (e.g., ['png', 'svg']).

    Returns:
    CircularAssemblyResult: The composed canvas and authoritative geometry.
    """
    # Configure and create canvas

    profile = canvas_config.profile
    cfg = profile.config
    resolved_depth_track_count = (
        int(depth_track_count_value)
        if depth_track_count_value is not None
        else depth_track_data_count([depth_tracks or ()])
    )
    effective_circular_track_slots, resolved_annotations, _annotation_layouts = _prepare_circular_annotation_tracks(
        gb_record,
        annotations,
        circular_track_slots,
        canvas_config=canvas_config,
        record_index=annotation_record_index,
        show_ticks=bool(cfg.objects.scale.show),
        show_depth=bool(profile.show_depth and depth_config is not None),
        show_gc=profile.show_gc,
        show_skew=profile.show_skew,
        depth_track_count=max(1, resolved_depth_track_count),
    )

    legend_table: dict = {}
    legend_measurement = legend_config.measure_legend(
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )
    if canvas_config.legend_position != "none":
        color_map = feature_config.specific_color_rules
        default_color_map = feature_config.default_color_map
        features_present = check_feature_presence(
            gb_record,
            feature_config.selected_features_set,
            feature_visibility_rules=feature_config.feature_visibility_rules,
            specific_color_rules=color_map,
            feature_instance_identity_plan=(
                feature_config.feature_instance_identity_plan
            ),
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            gb_record,
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
            feature_instance_identity_plan=(
                feature_config.feature_instance_identity_plan
            ),
        )
        legend_table = prepare_legend_table(
            gc_config,
            skew_config,
            feature_config,
            features_present,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=depth_config if (
                depth_config is not None
                and (depth_df is not None or resolved_depth_track_count == 1)
            ) else None,
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            show_depth=bool(
                profile.show_depth
                and depth_config is not None
                and (depth_df is not None or resolved_depth_track_count == 1)
            ),
        )
        if depth_tracks and effective_circular_track_slots is None:
            legend_table = sync_depth_track_legend_entries(legend_table, depth_tracks)
        legend_table = _sync_legend_table_for_circular_slots(
            legend_table,
            circular_track_slots=effective_circular_track_slots,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_config=depth_config,
            depth_df=depth_df,
            depth_tracks=depth_tracks,
            cfg=cfg,
            conservation_tracks=conservation_tracks,
            conservation_min_identity=conservation_min_identity,
        )
        normalized_annotation_slots = (
            normalize_circular_track_slots(effective_circular_track_slots)
            if effective_circular_track_slots is not None
            else None
        )
        legend_table = sync_annotation_legend_entries(
            legend_table,
            resolved_annotations,
            normalized_annotation_slots,
        )
    canvas: Drawing = canvas_config.create_svg_canvas()
    plot = add_record_on_circular_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        gc_config,
        skew_config,
        gc_df,
        species,
        strain,
        plot_title,
        legend_measurement,
        legend_table,
        legend_config=legend_config,
        depth_config=depth_config,
        depth_df=depth_df,
        depth_tracks=depth_tracks,
        depth_track_count_value=resolved_depth_track_count,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=conservation_min_identity,
        circular_track_slots=effective_circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_dataframes=dinucleotide_dataframes,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position=definition_position,
        definition_profile=definition_profile,
        definition_group_id=definition_group_id,
        center_reserved_radius=center_reserved_radius,
        _tick_track_channel_override=_tick_track_channel_override,
        annotations=resolved_annotations,
        annotation_record_index=annotation_record_index,
        definition_record_count=definition_record_count,
    )
    return _compose_circular_plot(
        plot,
        title_target=title_target,
        title_bounds=title_bounds,
        title_placement=title_placement,
    )


def assemble_circular_diagram(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    feature_config: FeatureDrawingConfigurator,
    species: Optional[str],
    strain: Optional[str],
    plot_title: Optional[str],
    legend_config: LegendDrawingConfigurator,
    depth_df: DataFrame | None = None,
    depth_config: DepthConfigurator | None = None,
    depth_tracks: Sequence[DepthTrackData] | None = None,
    depth_track_count_value: int | None = None,
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: str = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    _tick_track_channel_override: str | None = None,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None,
    annotation_record_index: int = 0,
    definition_record_count: int = 1,
) -> tuple[Drawing, LegendMeasurement]:
    """Assemble and compose a no-title Circular diagram.

    This release-backed wrapper intentionally preserves its historical return
    tuple.  Internal callers that need authoritative geometry or a plot title
    use :func:`_assemble_circular_diagram_result`.
    """
    result = _assemble_circular_diagram_result(
        gb_record=gb_record,
        canvas_config=canvas_config,
        gc_df=gc_df,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        species=species,
        strain=strain,
        plot_title=plot_title,
        legend_config=legend_config,
        depth_df=depth_df,
        depth_config=depth_config,
        depth_tracks=depth_tracks,
        depth_track_count_value=depth_track_count_value,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=conservation_min_identity,
        circular_track_slots=circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_dataframes=dinucleotide_dataframes,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position=definition_position,
        definition_profile=definition_profile,
        definition_group_id=definition_group_id,
        center_reserved_radius=center_reserved_radius,
        _tick_track_channel_override=_tick_track_channel_override,
        annotations=annotations,
        annotation_record_index=annotation_record_index,
        definition_record_count=definition_record_count,
    )
    return result.drawing, result.legend_measurement


__all__ = [
    "CircularAssemblyResult",
    "assemble_circular_diagram",
    "add_record_on_circular_canvas",
]
