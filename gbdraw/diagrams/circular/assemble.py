#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from gbdraw.features.placement import FeaturePlacementSlot, ResolvedRecordFeatureInputs

import logging
import math
import copy
from dataclasses import dataclass, replace
from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Collection, Optional, Mapping, Sequence, TypedDict, cast

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing
from svgwrite.container import Group

from ...auto_sizes import circular_track_width_px
from ...canvas import CircularCanvasConfigurator, CircularCanvasGeometry
from ...analysis.conservation import (
    ConservationLegendSource,
    ConservationTrack,
    conservation_legend_sources,
    conservation_track_gradient_colors,
)
from ...analysis.depth_tracks import (
    DepthTrackData,
    depth_track_data_count,
    index_depth_track_row,
    sync_depth_track_legend_entries,
)
from ...config.models import (
    CircularRenderProfile,
    GbdrawConfig,
)
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
from ...render_stages import _pre_keyed, render_stage
from ...core.text import calculate_bbox_dimensions, calculate_svg_bbox_dimensions
from ...config.modify import with_config_overrides, with_fit_row_overrides
from ...exceptions import ValidationError
from ...features.factory import FeatureBuildResult, create_feature_layers
from ...features.objects import FeatureObject
from ...labels.label_ink import circular_label_ink
from ...labels.label_style import label_leader_width_px
from ...render.drawers.label_backing import feature_label_drawn_parts
from ...labels.circular import (
    CircularLabelProbeStopped,
    assign_leader_start_points,
    circular_label_fits,
    circular_label_font_size,
    circular_label_shape,
    minimum_bbox_gap_px,
    prepare_label_list,
    x_overlap,
    y_overlap,
)
from ...labels.circular_fit import CircularDefinitionFitProbe, CircularLabelFitProbe, CircularLabelLayout
from ...labels.circular_radial import RADIAL_GROWTH_TOLERANCE_PX
from ...labels.circular_types import CircularLabel, CircularLabelCandidate, CircularLabelCandidateCache
from ...labels.fit import DEFINITION_FONT_PATH, FITS, LabelFit, LabelFitValue, LabelKind
from ...labels.filtering import label_filter_maps
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.circular import (
    CircularAxisLayout,
    CircularFeatureLaneDirection,
    CircularFeatureLayout,
    CircularRadialLayout,
    CircularRecordRenderContext,
    CircularResolvedSlot,
    CircularTickLayout,
    RadialBand,
    feature_radius_intervals,
)
from ...layout.circular_depth_axis import DepthAxisStyle, depth_axis_tick_font_size_px
from ...layout.circular_ink import CircularInk, CircularInkRing, box_meets_ink, stroke_pad_px
from ...layout.composition import (
    CompositionItem,
    CompositionPlan,
    CompositionRequest,
    DEFAULT_COMPOSITION_SPACING,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from ...layout.spatial import Aabb, union_aabbs
from ...legend.row_facts import LegendRowFacts, attach_legend_row_facts
from ...legend.sources import (
    DepthLegendSource,
    LegendSlotSource,
    LegendSources,
    NumericLegendSource,
    depth_legend_sources,
    feature_legend_presence,
    legend_slot_sources,
)
from ...legend.table import LegendTable, prepare_legend_table
from ...tracks import (
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
from ...render.drawers.circular.depth import DepthDrawer
from ...render.drawers.circular.gc_content import GcContentDrawer
from ...render.groups.circular.conservation import conservation_ring_ink
from ...render.groups.circular.ticks import TickGroup
from ...render.drawers.circular.gc_skew import SkewDrawer
from ...render.drawers.circular.annotations import (
    CircularAnnotationLabel,
    CircularLabelObstacles,
    LeaderSegment,
    circular_annotation_label_room,
    circular_annotation_labels,
    circular_annotation_track_content,
    circular_annotation_track_ink,
)
from ...render.group_content import CANVAS_SCOPE, GroupContent, GroupKey, GroupKind, place_content
from ...layout.linear_multi_record import record_key_at
from ...render.composition import apply_composition_plan
from ...svg.ids import track_slot_svg_id
from .positioning import _parse_svg_number as _svg_number, center_group_on_canvas
from ...tracks.circular import tick_sides_for_tick_label_layout
from ...tracks.parsing import slot_dinucleotide, slot_legend_label

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
from ...render.groups.circular.definition import CircularDefinitionProfile, CircularDefinitionStyle, DefinitionGroup
from .radial_layout import AnnotationLabelRoom, resolve_circular_radial_layout
from .presets import (
    CircularPresetContext,
    circular_feature_lane_direction_for_preset,
    circular_radial_plan_for_preset,
    circular_track_slots_from_preset_order,
    normalize_circular_track_preset,
)

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation


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
    label_layout: CircularLabelLayout | None


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
    # The drawn external labels: their fit by kind and a probe for the fit search (None without them).
    label_layout: CircularLabelLayout | None = None


def plan_circular_annotation_slots(
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    records: Sequence[SeqRecord],
    slots: list[CircularTrackSlot] | None,
    *,
    show_ticks: bool,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
    depth_track_count: int,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    record_indices: Sequence[int] | None = None,
) -> tuple[list[CircularTrackSlot] | None, ResolvedAnnotationBundle]:
    """Bind annotation sets to the Circular slot list that the records use."""

    slots, bundle, _auto_slot_ids = prepare_annotation_track_slots(
        annotations,
        records,
        slots,
        mode="circular",
        record_transforms=record_transforms,
        record_indices=record_indices,
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
    return slots, bundle


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
    record_transform: RecordDisplayTransform | None = None,
) -> tuple[list[CircularTrackSlot] | None, ResolvedAnnotationBundle, dict[str, ResolvedAnnotationTrack]]:
    slots, bundle = plan_circular_annotation_slots(
        annotations,
        [gb_record],
        slots,
        show_ticks=show_ticks,
        show_depth=show_depth,
        show_gc=show_gc,
        show_skew=show_skew,
        depth_track_count=depth_track_count,
        record_transforms=(record_transform,) if record_transform is not None else None,
        record_indices=(record_index,),
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, {}
    # The planner generates the slot list whenever annotations are present.
    assert slots is not None
    updated, layouts = _layout_circular_annotation_tracks(
        gb_record, slots, bundle, canvas_config=canvas_config, record_index=record_index
    )
    return updated, bundle, layouts


def _layout_circular_annotation_tracks(
    gb_record: SeqRecord,
    slots: list[CircularTrackSlot],
    bundle: ResolvedAnnotationBundle,
    *,
    canvas_config: CircularCanvasConfigurator,
    record_index: int,
) -> tuple[list[CircularTrackSlot], dict[str, ResolvedAnnotationTrack]]:
    """Resolve the lanes of each annotation slot; a slot without a width takes its required extent."""

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
            label_font_size=circular_label_font_size(canvas_config.profile.config, len(gb_record.seq)),
        )
        layouts[slot.id] = layout
        updated.append(
            slot
            if slot.width is not None
            else replace(slot, width=ScalarSpec(layout.required_extent_px, "px"))
        )
    return updated, layouts


def _add_circular_feature_underlays(
    gb_record: SeqRecord,
    underlay_features: Sequence,
    slots: list[CircularTrackSlot],
    bundle: ResolvedAnnotationBundle,
    layouts: dict[str, ResolvedAnnotationTrack],
    *,
    canvas_config: CircularCanvasConfigurator,
    record_index: int,
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
        return slots, bundle, layouts
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
    # The underlay set is new (collision-safe id), so the slots laid out before keep their layouts.
    sized_underlay, underlay_layouts = _layout_circular_annotation_tracks(
        gb_record, [underlay_slot], merged, canvas_config=canvas_config, record_index=record_index
    )
    return [*sized_underlay, *slots], merged, {**underlay_layouts, **layouts}


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


def _label_target_unwrapped_angle(label: CircularLabel, total_length: int) -> float:
    """Target label angle derived from genomic midpoint, in unwrapped radians."""
    return (2.0 * math.pi * (float(label["middle"]) / float(total_length))) - (0.5 * math.pi)


def _label_unwrapped_angle(label: CircularLabel, total_length: int) -> float:
    """Current label angle projected near its target branch."""
    current_angle = math.atan2(float(label["start_y"]), float(label["start_x"]))
    target_angle = _label_target_unwrapped_angle(label, total_length)
    return _unwrap_angle_near_reference(current_angle, target_angle)


def _label_bbox_local(label: CircularLabel, total_length: int, margin_px: float) -> tuple[float, float, float, float]:
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
    label: CircularLabel,
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
    label: CircularLabel,
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

    if "text_x" in label:
        label["text_x"] = float(label["text_x"]) + dx
    if "text_y" in label:
        label["text_y"] = float(label["text_y"]) + dy
    if "leader_start_x" in label:
        label["leader_start_x"] = float(label["leader_start_x"]) + dx
    if "leader_start_y" in label:
        label["leader_start_y"] = float(label["leader_start_y"]) + dy

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
    labels: list[CircularLabel],
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


def _external_label_paint_bounds_local(label: CircularLabel, total_length: int) -> Aabb:
    """Return exact text plus leader paint bounds (the leader at its label's own width) in circle-local
    coordinates."""
    measured_label = label.copy()
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
    ).expanded(0.5 * label_leader_width_px(label))
    combined = union_aabbs((text_bounds, leader_bounds))
    if combined is None:  # pragma: no cover - both inputs are present
        raise RuntimeError("external label paint bounds unexpectedly empty")
    return combined


def _external_label_obstacles_on_canvas(
    labels: Sequence[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
) -> tuple[Aabb, ...]:
    """Return authoritative per-label overlay obstacles in plot coordinates."""
    return tuple(
        _external_label_paint_bounds_local(label, total_length).translated(
            float(canvas_config.offset_x), float(canvas_config.offset_y)
        )
        for label in labels
        if not label.get("is_embedded")
    )


def _circular_annotation_labels_by_slot(
    radial_layout: CircularRadialLayout,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None,
    *,
    record_index: int,
    record_length: int,
    font_family: str,
    obstacles: CircularLabelObstacles,
) -> dict[str, tuple[CircularAnnotationLabel, ...]]:
    """The labels of each annotation slot of a record, measured once for the drawer and the canvas bounds.

    The letters clear the axis circle's paint too where they would meet it, and
    move along their regions off ``obstacles`` and the labels of earlier slots.
    """

    axis_paint = radial_layout.axis.paint_band_px
    labels: dict[str, tuple[CircularAnnotationLabel, ...]] = {}
    for slot in radial_layout.slots:
        layout = (annotation_track_layouts or {}).get(slot.id)
        if str(slot.renderer) != "annotations" or layout is None:
            continue
        labels[slot.id] = circular_annotation_labels(
            layout,
            record_index=record_index,
            record_length=record_length,
            inner_radius_px=float(slot.draw_inner_radius_px),
            outer_radius_px=float(slot.draw_outer_radius_px),
            side=str(slot.side),
            font_family=font_family,
            params=annotation_track_params_from_mapping(slot.params),
            rings=(axis_paint,),
            obstacles=obstacles,
            covers_slots=slot.spans_slots,
        )
        obstacles = obstacles.with_labels(labels[slot.id])
    return labels


# Letters lie on a track when they reach more than this into its band, as a label covers the definition text.
REGION_TRACK_MIN_OVERLAP_PX = 0.5


def _region_labels_on_other_tracks(
    radial_layout: CircularRadialLayout, labels: Mapping[str, Sequence[CircularAnnotationLabel]]
) -> int:
    """How many annotation labels' letters lie on the band another track takes.

    The layout keeps a room for these letters beyond an inside or outside
    annotation slot and beyond a span overlay (OV-1192), so they lie on a
    track only where set radii leave no room (OV-1207); the figure is then
    drawn and warns (LP-D86). A track takes the band the layout keeps for it
    (``own_reserved_band_px``): its drawing with its tick or axis labels and
    its own labels' room. Annotation overlays draw over other tracks and take
    none. An overlay the layout keeps no room for (``span_keeps_label_room``:
    a span whose outer edge is a plot overlay) has its labels drawn over the
    tracks around it whatever the radii, so they are not counted.
    """

    tracks = _track_bands(radial_layout)
    count = 0
    for slot in radial_layout.slots:
        own = slot.draw_band_px
        if own is None or (slot.side == "overlay" and not slot.span_keeps_label_room):
            continue
        others = [
            band for track, band in tracks if track is not slot and not _reaches_band(own.inner_px, own.outer_px, band)
        ]
        count += sum(
            1 for label in labels.get(slot.id, ()) if any(_reaches_band(*label.ink_radii, band) for band in others)
        )
    return count


def _track_bands(radial_layout: CircularRadialLayout) -> list[tuple[CircularResolvedSlot, RadialBand]]:
    """The band each drawn track takes (``own_reserved_band_px``: its drawing with its tick or axis labels and its
    own labels' room); annotation overlays draw over other tracks and take none."""

    return [
        (slot, band)
        for slot in radial_layout.slots
        if slot.draw_band_px is not None
        and not (slot.renderer == "annotations" and slot.side == "overlay")
        and (band := slot.own_reserved_band_px) is not None
    ]


def _reaches_band(inner: float, outer: float, band: RadialBand) -> bool:
    """Whether ink between radii ``inner`` and ``outer`` reaches more than the threshold into ``band``."""

    return min(outer - band.inner_px, band.outer_px - inner) > REGION_TRACK_MIN_OVERLAP_PX


def _record_track_ink(
    radial_layout: CircularRadialLayout,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    record_transform: RecordDisplayTransform | None,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
    conservation_tracks: Sequence[ConservationTrack] | None,
    conservation_min_identity: float,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None,
    gc_content_tick_font_size_override: float | None,
    feature_dict: Mapping[str, FeatureObject],
    feature_config: FeatureDrawingConfigurator,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None,
    annotation_labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    annotation_record_index: int,
    use_slot_tick_options: bool,
) -> tuple[CircularInk, ...]:
    """What the record's tracks paint, from the owners that draw them (D-18), read with the same slot sources
    the drawing reads (``_draw_resolved_circular_slot``).

    The GC content, GC skew and depth areas and their axes come from their drawers, the feature blocks from the
    intervals the labels clear (``feature_radius_intervals``), the conservation hits and background from the
    ring's owner, the annotation marks and labels from the annotation track's owner, the tick marks and labels
    from the tick group, and the axis circle from its paint band. A spacer paints nothing.
    """

    canvas_radius = float(canvas_config.radius)
    record_length = len(gb_record.seq)
    font_family = str(canvas_config.profile.config.objects.text.font_family)
    inks: list[CircularInk] = [CircularInkRing(radial_layout.axis.paint_band_px)]
    for slot in radial_layout.slots:
        renderer = str(slot.renderer)
        if slot.draw_band_px is None:
            continue
        norm_factor = _slot_norm_factor(slot, canvas_radius)
        if renderer == "dinucleotide_content":
            _, content_df, content_config = _slot_content_source(
                slot,
                gc_df=gc_df,
                gc_config=gc_config,
                dataframes=dinucleotide_content_dataframes or dinucleotide_dataframes,
                tick_font_size_override=gc_content_tick_font_size_override,
            )
            if content_df is not None:
                inks.extend(
                    GcContentDrawer(content_config).ink(
                        canvas_radius, content_df, record_length, slot.draw_width_px, norm_factor,
                        content_config.dinucleotide, record_transform, font_family=font_family,
                    )
                )
        elif renderer == "dinucleotide_skew":
            _, skew_df, slot_skew_config = _slot_skew_source(
                slot,
                gc_df=gc_df,
                gc_config=gc_config,
                skew_config=skew_config,
                dataframes=dinucleotide_skew_dataframes or dinucleotide_dataframes,
            )
            if skew_df is not None:
                inks.extend(
                    SkewDrawer(slot_skew_config).ink(
                        canvas_radius, skew_df, record_length, slot.draw_width_px, norm_factor,
                        slot_skew_config.dinucleotide, record_transform,
                    )
                )
        elif renderer == "depth":
            depth_source = _slot_depth_source(slot, depth_df, depth_config, depth_tracks_by_index)
            if depth_source is not None and depth_source[0] is not None and depth_source[1] is not None:
                inks.extend(
                    DepthDrawer(depth_source[1]).ink(
                        canvas_radius, depth_source[0], record_length, float(slot.draw_width_px), norm_factor,
                        font_family=font_family, record_transform=record_transform,
                    )
                )
        elif renderer == "sequence_conservation":
            track = _slot_conservation_track(slot, conservation_tracks)
            if track is not None:
                inks.extend(
                    conservation_ring_ink(
                        track,
                        total_length=record_length,
                        inner_radius_px=float(slot.draw_inner_radius_px),
                        outer_radius_px=float(slot.draw_outer_radius_px),
                        min_identity=float(conservation_min_identity),
                        cfg=canvas_config.profile.config,
                        record_transform=record_transform,
                    )
                )
        elif renderer == "annotations":
            annotation_layout = (annotation_track_layouts or {}).get(slot.id)
            if annotation_layout is not None:
                inks.extend(
                    circular_annotation_track_ink(
                        annotation_layout,
                        record_index=annotation_record_index,
                        record_length=record_length,
                        inner_radius_px=float(slot.draw_inner_radius_px),
                        outer_radius_px=float(slot.draw_outer_radius_px),
                        side=str(slot.side),
                        params=annotation_track_params_from_mapping(slot.params),
                        labels=annotation_labels.get(slot.id, ()),
                    )
                )
        elif renderer == "ticks":
            tick_options = _slot_tick_options(slot, use_slot_tick_options)
            if tick_options is not None:
                inks.extend(
                    TickGroup(
                        gb_record,
                        canvas_config.geometry(),
                        profile=canvas_config.profile,
                        radius=float(slot.anchor_radius_px or 0.0),
                        record_transform=record_transform,
                        **tick_options,
                    ).ink()
                )
        elif renderer == "features" and isinstance(slot.payload, CircularFeatureLayout):
            pad = stroke_pad_px(feature_config.block_stroke_color, feature_config.block_stroke_width)
            inks.extend(
                CircularInkRing(RadialBand(inner - pad, outer + pad), start / record_length, end / record_length)
                for start, end, inner, outer in feature_radius_intervals(feature_dict, record_length, slot.payload)
            )
    return tuple(inks)


def _definition_lines_on_tracks(inks: Sequence[CircularInk], definition_text: Sequence[Aabb]) -> int:
    """How many lines of the centre definition text meet what a track paints (D-18).

    ``definition_text`` holds the painted line boxes centred on the record
    (``_centre_definition_text``); a line counts when its box reaches more than
    the threshold that region labels use (LP-D86) into a track's paint
    (``_record_track_ink``). The figure is drawn as before and warns.
    """

    return sum(1 for box in definition_text if box_meets_ink(box, inks, threshold_px=REGION_TRACK_MIN_OVERLAP_PX))


def _feature_label_obstacles(
    labels: Sequence[CircularLabel], total_length: int, canvas_config: CircularCanvasConfigurator
) -> CircularLabelObstacles:
    """The record's drawn feature labels and their leaders, as the fit count sees them, for annotation labels to
    keep the label spacing from (OV-1173). A label with a leader also covers what its style draws with its letters,
    its underline and its white box (OV-1358, the Linear rule of OV-1357)."""

    style = canvas_config.profile.label_style
    shapes = [
        (
            label,
            shape,
            ()
            if label.get("is_embedded")
            else feature_label_drawn_parts(circular_label_ink(label, total_length, style), style),
        )
        for label in labels
        if (shape := circular_label_shape(label, total_length)).parts
    ]
    return CircularLabelObstacles(
        shapes=tuple((*shape.parts, *drawn) for _label, shape, drawn in shapes),
        segments=tuple(
            LeaderSegment(start, end, 0.5 * label_leader_width_px(label))
            for label, shape, _drawn in shapes
            for start, end in zip(shape.leader, shape.leader[1:])
        ),
        clearance=float(canvas_config.profile.config.labels.spacing.circular),
    )


@dataclass(frozen=True)
class _AnnotationLabelRoom:
    """The room one annotation slot's labels need beyond a band (``AnnotationLabelRoom``), measured by the labels'
    own owner; a frozen value, so the label stage that replans the radial layout stays cacheable."""

    track: ResolvedAnnotationTrack
    record_index: int
    record_length: int
    font_family: str

    def __call__(
        self, band: RadialBand, axis: CircularAxisLayout, side: str, params: Mapping[str, Any], covers_slots: bool
    ) -> float:
        return circular_annotation_label_room(
            self.track,
            record_index=self.record_index,
            record_length=self.record_length,
            band=band,
            side=side,
            font_family=self.font_family,
            params=annotation_track_params_from_mapping(params),
            rings=(axis.paint_band_px,),
            covers_slots=covers_slots,
        )


def _circular_annotation_label_rooms(
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None,
    *,
    record_index: int,
    record_length: int,
    font_family: str,
) -> dict[str, AnnotationLabelRoom]:
    """The room each annotation slot's labels need beyond its band, measured by the labels' own owner."""

    return {
        slot_id: _AnnotationLabelRoom(track, record_index, record_length, font_family)
        for slot_id, track in (annotation_track_layouts or {}).items()
    }


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
    record_obstacles: Sequence[Aabb],
) -> Aabb:
    """Return where the composer puts a corner legend before label adjustment.

    ``record_obstacles`` are the obstacles of the final plan that labels cannot
    clear (the record and its tracks), so the predicted corner is one the final
    plan can take; labels are left out because they move off the prediction
    (OV-1609).
    """
    if not placement.is_overlay:
        raise ValueError("collision bounds require an overlay legend placement")
    plan = render_stage(
        "composition_plan",
        plan_composition,
        CompositionRequest(
            primary=CompositionItem("primary", primary_bounds),
            legend=CompositionItem("legend", legend_measurement.local_bounds),
            legend_placement=placement,
            overlay_obstacles=tuple(record_obstacles),
        ),
    )
    primary = plan.placement_for("primary")
    legend = plan.placement_for("legend")
    if primary is None or legend is None:  # pragma: no cover - request has both items
        raise RuntimeError("overlay collision plan is incomplete")
    return legend.final_bounds.translated(-primary.dx, -primary.dy)


def _legend_collision_indices(
    labels: list[CircularLabel],
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
    labels: list[CircularLabel],
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
    label: CircularLabel,
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
    labels: list[CircularLabel],
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
    labels: list[CircularLabel],
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
                    candidate_labels: dict[int, CircularLabel] = {}
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
    labels: list[CircularLabel],
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
    geometry: CircularCanvasGeometry,
    species: str | None,
    strain: str | None,
    *,
    cfg: GbdrawConfig,
    plot_title: str | None,
    definition_profile: CircularDefinitionProfile,
    record_presentation: RecordPresentation | None = None,
) -> float:
    definition_group = DefinitionGroup(
        gb_record,
        geometry,
        style=CircularDefinitionStyle.of(cfg, definition_profile),
        species=species,
        strain=strain,
        plot_title=plot_title,
        record_presentation=record_presentation,
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
            int(geometry.dpi),
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
    padding_px = max(8.0, 0.02 * float(geometry.radius))
    return float(max_extent + padding_px)


def _definition_clearance_px(layout: CircularRadialLayout) -> float:
    """How far the centre definition's band may reach before it enters a drawn track's band (``_track_bands``).

    This is the narrow row's reading of "the definition over a track": the
    definition's band (``_definition_reserved_radius_px``: the largest
    ``max(|x +/- w/2|, |y +/- h/2|)`` of its lines, plus padding) against each
    track's band (its drawing, its tick or axis labels, and its labels' room).
    The notice reads another one: the painted line boxes against what the
    tracks paint (``_definition_lines_on_tracks`` over ``_record_track_ink``).
    The band is not the boxes' corner distance, so a figure narrowed to it may
    still have letters on a track and warn (MERS-CoV at 231.99 px).
    """

    return min((float(band.inner_px) for _slot, band in _track_bands(layout)), default=float("inf"))


@dataclass(frozen=True)
class _CentreDefinitionFit:
    """The tracks and the centre definition drawn among them (PD-OI-078).

    ``reserved_radius_px`` is the definition's band (None without a centre definition); radial labels keep out of
    it. ``narrow_radius_px`` is that band when it enters a track's band (``_definition_clearance_px``), the
    ``center_reserved_radius`` that narrows the inside tracks around it (the "narrow the inside tracks" row); 0 when
    the band clears the tracks' bands or the radius is explicit already.
    """

    layout: CircularRadialLayout
    reserved_radius_px: float | None
    narrow_radius_px: float = 0.0


def _resolve_radial_layout_fitting_definition(
    resolve_layout: Callable[[float | None], CircularRadialLayout],
    *,
    gb_record: SeqRecord,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    definition_profile: CircularDefinitionProfile,
    layout_definition_profile: CircularDefinitionProfile,
    definition_position: str,
    center_reserved_radius: float | None,
    record_presentation: RecordPresentation | None,
) -> _CentreDefinitionFit:
    """Place the tracks around the centre definition, then draw it among them (PD-OI-078, rule L).

    The definition is never wrapped (DF-D25). The radial layout places the
    tracks around its band (``resolve_circular_radial_layout``): it moves the
    inside tracks outward, then narrows their Auto widths down to
    ``canvas.circular.track_narrowing_min_px``, and when that does not clear
    the band it keeps the widths and the definition is drawn over the tracks;
    the label fit check warns when its painted lines meet a track's ink
    (``_definition_lines_on_tracks``, PD-OI-097). An explicit
    ``center_reserved_radius`` keeps the tracks out of the centre (DF-D05).

    The tracks are placed around the band of ``layout_definition_profile``,
    the definition the record draws without a plot title or, for rows that
    pack on a Multi-Record Canvas, alone (OV-1710, OV-1711), so a title or the
    canvas changes only the centre text. The band and the narrow row read the
    drawn definition (``definition_profile``).
    """

    if center_reserved_radius is not None:
        reserved = max(0.0, float(center_reserved_radius))
        return _CentreDefinitionFit(resolve_layout(reserved), reserved)
    if str(definition_position).strip().lower() != "center":
        return _CentreDefinitionFit(resolve_layout(None), None)

    def band_px(profile: CircularDefinitionProfile) -> float:
        return _definition_reserved_radius_px(
            gb_record,
            geometry,
            species,
            strain,
            cfg=cfg,
            plot_title=plot_title,
            definition_profile=profile,
            record_presentation=record_presentation,
        )

    reserved = band_px(definition_profile)
    layout_band = reserved if layout_definition_profile == definition_profile else band_px(layout_definition_profile)
    if max(reserved, layout_band) <= FEATURE_BAND_EPSILON:
        return _CentreDefinitionFit(resolve_layout(None), None)
    layout = resolve_layout(float(layout_band) if layout_band > FEATURE_BAND_EPSILON else None)
    narrow_radius = 0.0 if reserved <= _definition_clearance_px(layout) + FEATURE_BAND_EPSILON else reserved
    return _CentreDefinitionFit(layout, reserved, narrow_radius)


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


def _slot_content_source(
    resolved_slot: CircularResolvedSlot,
    *,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    dataframes: dict[str, DataFrame] | None,
    tick_font_size_override: float | None,
) -> tuple[str, DataFrame | None, GcContentConfigurator]:
    """The dinucleotide, data (``None`` when there are none) and settings a content slot draws with."""

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = slot_dinucleotide(resolved_slot.params, default_nt)
    slot_df = _slot_dataframe_for_nt(nt=nt, default_df=gc_df, default_nt=default_nt, dinucleotide_dataframes=dataframes)
    config = _slot_gc_config_with_axis_font_size(gc_config, nt, tick_font_size_override)
    return nt, slot_df if f"{nt} content" in slot_df.columns else None, config


def _slot_skew_source(
    resolved_slot: CircularResolvedSlot,
    *,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    dataframes: dict[str, DataFrame] | None,
) -> tuple[str, DataFrame | None, GcSkewConfigurator]:
    """The dinucleotide, data (``None`` when there are none) and settings a skew slot draws with."""

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = slot_dinucleotide(resolved_slot.params, default_nt)
    slot_df = _slot_dataframe_for_nt(nt=nt, default_df=gc_df, default_nt=default_nt, dinucleotide_dataframes=dataframes)
    config = _slot_skew_config(skew_config, resolved_slot, nt)
    return nt, slot_df if f"{nt} skew" in slot_df.columns else None, config


def _slot_norm_factor(resolved_slot: CircularResolvedSlot, canvas_radius: float) -> float:
    """The slot's anchor radius as a share of the canvas radius, as its drawer takes it."""

    # Radial layout always resolves a radius for a drawn slot.
    anchor_radius_px = resolved_slot.anchor_radius_px
    assert anchor_radius_px is not None
    return float(anchor_radius_px) / float(canvas_radius)


def _slot_depth_source(
    resolved_slot: CircularResolvedSlot,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
) -> tuple[DataFrame | None, DepthConfigurator | None, DepthTrackData | None] | None:
    """The data and settings a depth slot draws with (``None`` where they are missing) and its track when there
    are several; ``None`` when the slot names a track there is not."""

    if depth_tracks_by_index is None:
        return depth_df, depth_config, None
    track_index = parse_nonnegative_integer(
        resolved_slot.params.get("track_index", 0),
        field_name=f"depth slot '{resolved_slot.id}' track_index",
    )
    depth_track = depth_tracks_by_index.get(track_index)
    if depth_track is None:
        return None
    return depth_track.df, depth_track.config, depth_track


def _slot_conservation_track(
    resolved_slot: CircularResolvedSlot, conservation_tracks: Sequence[ConservationTrack] | None
) -> ConservationTrack | None:
    """The conservation track a conservation slot draws."""

    track_index = int(resolved_slot.params.get("track_index", 0) or 0)
    return next(
        (track for track in (conservation_tracks or ()) if int(track.track_index) == track_index),
        None,
    )


class _TickOptions(TypedDict):
    label_side: str
    tick_side: str
    tick_length_px: float | None
    track_preset: str | None


def _slot_tick_options(resolved_slot: CircularResolvedSlot, use_slot_tick_options: bool) -> _TickOptions | None:
    """The sides, tick length and preset a tick slot draws with; ``None`` when it draws neither ticks nor
    labels."""

    tick_layout = resolved_slot.payload if isinstance(resolved_slot.payload, CircularTickLayout) else None
    if tick_layout is not None:
        label_side, tick_side = str(tick_layout.label_side), str(tick_layout.tick_side)
    else:
        label_side, tick_side = tick_sides_for_tick_label_layout(
            resolved_slot.params.get("tick_label_layout"),
            side=resolved_slot.side,
        )[:2]
    label_side, tick_side = label_side.strip().lower(), tick_side.strip().lower()
    if label_side == "none" and tick_side == "none":
        return None
    if not (use_slot_tick_options or tick_layout is not None):
        return {"label_side": "legacy", "tick_side": "legacy", "tick_length_px": None, "track_preset": None}
    return {
        "label_side": label_side,
        "tick_side": tick_side,
        "tick_length_px": (
            tick_layout.tick_length_px
            if tick_layout is not None
            else (
                float(resolved_slot.draw_width_px)
                if resolved_slot.explicit_width and resolved_slot.draw_width_px > 0
                else None
            )
        ),
        "track_preset": (
            tick_layout.track_preset
            if tick_layout is not None
            else normalize_circular_track_preset(
                str(resolved_slot.params.get("preset", resolved_slot.params.get("track_preset", "tuckin")))
            )
        ),
    }




def _sync_legend_table_for_circular_slots(
    out: LegendTable,
    *,
    slots: Sequence[LegendSlotSource] | None,
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    depth: DepthLegendSource | None,
    depth_tracks: Sequence[DepthLegendSource] | None = None,
    conservation_legend: Sequence[ConservationLegendSource] = (),
    conservation_colors: tuple[str, str] = ("", ""),
    conservation_min_identity: float | None = None,
) -> LegendTable:
    """Replace singleton numeric legend entries with slot-aware entries.

    ``depth`` is the lone depth track's swatch when its frame was built;
    ``conservation_colors`` the default gradient ends (min, max).
    """
    if slots is None:
        return out

    default_nt = str(gc.dinucleotide).upper()
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
        out.discard_default_track(key)

    depth_by_index = index_depth_track_row(depth_tracks or ())
    for slot in slots:
        renderer = slot.renderer
        if renderer == "depth":
            if depth_by_index:
                track_index = parse_nonnegative_integer(
                    slot.params.get("track_index", 0),
                    field_name=f"depth slot '{slot.id}' track_index",
                )
                track = depth_by_index.get(track_index)
                if track is None:
                    continue
                swatch = track
                label = slot_legend_label(slot, track.label)
            elif depth is None:
                continue
            else:
                swatch = depth
                label = slot_legend_label(slot, "Depth")
            out.add_source(label, {
                "type": "solid",
                "fill": swatch.fill_color,
                "stroke": swatch.stroke_color,
                "width": swatch.stroke_width,
            })
        elif renderer == "dinucleotide_content":
            nt = slot_dinucleotide(slot.params, default_nt)
            out.add_numeric_source(slot_legend_label(slot, f"{nt} content"), gc)
        elif renderer == "dinucleotide_skew":
            nt = slot_dinucleotide(slot.params, default_nt)
            out.add_numeric_source(slot_legend_label(slot, f"{nt} skew"), skew.for_skew_slot(slot, nt))
    if conservation_legend:
        default_min_color, default_max_color = conservation_colors
        if any(track.track_color for track in conservation_legend):
            for conservation_track in conservation_legend:
                min_color, max_color = conservation_track_gradient_colors(
                    conservation_track.track_color,
                    default_min_color=default_min_color,
                    default_max_color=default_max_color,
                )
                out.add_source(conservation_track.track_label, {
                    "type": "gradient",
                    "min_color": min_color,
                    "max_color": max_color,
                    "stroke": "none",
                    "width": 0,
                    "min_value": float(conservation_min_identity or 0.0),
                })
        else:
            out.add_source("Conservation identity", {
                "type": "gradient",
                "min_color": default_min_color,
                "max_color": default_max_color,
                "stroke": "none",
                "width": 0,
                "min_value": float(conservation_min_identity or 0.0),
            })
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
    precalculated_labels: list[CircularLabel] | None = None,
    use_slot_group_id: bool = True,
    use_slot_tick_options: bool = True,
    use_feature_anchor_override: bool = True,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack] | None = None,
    annotation_labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    annotation_record_index: int = 0,
    record_transform: RecordDisplayTransform | None = None,
    record_key: str,
) -> Drawing:
    """Draw one resolved circular slot."""
    cfg = render_context.profile.config
    renderer = str(resolved_slot.renderer)
    # Radial layout always resolves a radius; float(None) already failed here before.
    anchor_radius_px = resolved_slot.anchor_radius_px
    assert anchor_radius_px is not None
    norm_factor_override = _slot_norm_factor(resolved_slot, float(canvas_config.radius))
    if renderer == "spacer":
        return canvas

    if renderer == "annotations":
        annotation_layout = (annotation_track_layouts or {}).get(resolved_slot.id)
        if annotation_layout is None:
            return canvas
        params = annotation_track_params_from_mapping(resolved_slot.params)
        content = render_stage(
            "annotation_content",
            circular_annotation_track_content,
            annotation_layout,
            record_id=str(gb_record.id),
            record_index=annotation_record_index,
            record_length=len(gb_record.seq),
            inner_radius_px=float(resolved_slot.draw_inner_radius_px),
            outer_radius_px=float(resolved_slot.draw_outer_radius_px),
            side=str(resolved_slot.side),
            font_family=str(cfg.objects.text.font_family),
            params=params,
            labels=annotation_labels.get(resolved_slot.id, ()),
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
        group = place_content(canvas, content, GroupKey(record_key, GroupKind.ANNOTATIONS, str(resolved_slot.id)))
        canvas.add(center_group_on_canvas(group, canvas_config))
        return canvas

    if renderer == "features":
        base_width = circular_track_width_px(
            "features",
            radius=canvas_config.radius,
            track_ratio=canvas_config.track_ratio,
            factors=cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)],
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
            "group_key": GroupKey(record_key, GroupKind.FEATURES, str(resolved_slot.id)),
            "feature_layers": feature_layers,
            "render_context": render_context,
            # The record group draws the embedded labels; the legend shift moves only external ones.
            "precalculated_labels": tuple(
                label for label in precalculated_labels or () if label.get("is_embedded")
            ),
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
            float(anchor_radius_px),
            float(canvas_config.radius),
            rel_tol=1e-9,
            abs_tol=1e-9,
        ):
            feature_kwargs["feature_anchor_radius_px"] = float(anchor_radius_px)
        return add_record_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            **feature_kwargs,
        )

    if renderer == "ticks":
        tick_options = _slot_tick_options(resolved_slot, use_slot_tick_options)
        if tick_options is not None:
            tick_group_kwargs: dict[str, Any] = {
                "group_key": GroupKey(record_key, GroupKind.TICKS, str(resolved_slot.id)),
                "radius_override": float(anchor_radius_px),
                **tick_options,
            }
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
                record_transform=record_transform,
            )
        return canvas

    if renderer == "depth":
        depth_source = _slot_depth_source(resolved_slot, depth_df, depth_config, depth_tracks_by_index)
        if depth_source is None:
            return canvas
        selected_depth_df, selected_depth_config, depth_track = depth_source
        depth_group_id = depth_track.id if depth_track is not None and not use_slot_group_id else str(resolved_slot.id)
        if selected_depth_config is None or selected_depth_df is None:
            logger.warning("Skipping circular depth slot '%s' because depth data are unavailable.", resolved_slot.id)
            return canvas
        depth_kwargs: dict[str, Any] = {
            "group_key": GroupKey(record_key, GroupKind.DEPTH, str(resolved_slot.id)),
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
        }
        if use_slot_group_id:
            depth_kwargs.update(
                group_id=depth_group_id,
                slot_group_id=track_slot_svg_id(
                    resolved_slot.id,
                    renderer=renderer,
                    slot_index=int(resolved_slot.slot_index),
                ),
                slot_id=str(resolved_slot.id),
                slot_renderer=renderer,
            )
        return add_depth_group_on_canvas(
            canvas,
            gb_record,
            selected_depth_df,
            canvas_config,
            selected_depth_config,
            **depth_kwargs,
            record_transform=record_transform,
        )

    if renderer == "sequence_conservation":
        track = _slot_conservation_track(resolved_slot, conservation_tracks)
        if track is None:
            logger.warning(
                "Skipping circular conservation slot '%s' because conservation data are unavailable.",
                resolved_slot.id,
            )
            return canvas
        conservation_kwargs: dict[str, Any] = {
            "group_key": GroupKey(record_key, GroupKind.CONSERVATION, str(resolved_slot.id)),
        }
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
            record_index=annotation_record_index,
            record_transform=record_transform,
            **conservation_kwargs,
        )

    if renderer == "dinucleotide_content":
        nt, slot_df, slot_gc_config = _slot_content_source(
            resolved_slot,
            gc_df=gc_df,
            gc_config=gc_config,
            dataframes=dinucleotide_content_dataframes or dinucleotide_dataframes,
            tick_font_size_override=gc_content_tick_font_size_override,
        )
        if slot_df is None:
            logger.warning(
                "Skipping circular content slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        gc_kwargs: dict[str, Any] = {
            "group_key": GroupKey(record_key, GroupKind.GC_CONTENT, str(resolved_slot.id)),
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
            "slot_id": str(resolved_slot.id),
            "slot_renderer": renderer,
        }
        if use_slot_group_id:
            gc_kwargs["group_id"] = str(resolved_slot.id)
            gc_kwargs["slot_group_id"] = track_slot_svg_id(
                resolved_slot.id,
                renderer=renderer,
                slot_index=int(resolved_slot.slot_index),
            )
        return add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_gc_config,
            **gc_kwargs,
            record_transform=record_transform,
        )

    if renderer == "dinucleotide_skew":
        nt, slot_df, slot_skew_config = _slot_skew_source(
            resolved_slot,
            gc_df=gc_df,
            gc_config=gc_config,
            skew_config=skew_config,
            dataframes=dinucleotide_skew_dataframes or dinucleotide_dataframes,
        )
        if slot_df is None:
            logger.warning(
                "Skipping circular skew slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        skew_kwargs: dict[str, Any] = {
            "group_key": GroupKey(record_key, GroupKind.GC_SKEW, str(resolved_slot.id)),
            "track_width_override": float(resolved_slot.draw_width_px),
            "norm_factor_override": norm_factor_override,
            "slot_id": str(resolved_slot.id),
            "slot_renderer": renderer,
        }
        if use_slot_group_id:
            skew_kwargs["group_id"] = str(resolved_slot.id)
            skew_kwargs["slot_group_id"] = track_slot_svg_id(
                resolved_slot.id,
                renderer=renderer,
                slot_index=int(resolved_slot.slot_index),
            )
        return add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_skew_config,
            **skew_kwargs,
            record_transform=record_transform,
        )

    logger.warning("Skipping unsupported circular track slot renderer '%s'.", renderer)
    return canvas


def _default_outer_label_arena(
    *,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    strandedness: bool,
) -> tuple[float, float]:
    """Return default (anchor_radius_px, arc_outer_radius_px) for outer labels.

    ``cfg`` gives the labels settings (a fit probe passes its own).
    """
    length_param = str(geometry.length_param)
    track_type = _circular_preset_for_layout(cfg)
    strands = "separate" if strandedness else "single"
    base_radius = float(geometry.radius)

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




def _record_label_layout(
    gb_record: SeqRecord,
    profile: CircularRenderProfile,
    drawn_labels: Sequence[CircularLabel],
    place_labels: Callable[..., list[CircularLabel]] | None,
    *,
    record_index: int,
    definition_text: Sequence[Aabb],
    radial_refused: bool,
    region_labels: Sequence[CircularAnnotationLabel] = (),
    region_track_overlaps: int = 0,
    definition_track_overlaps: int = 0,
    definition_probe: CircularDefinitionFitProbe | None = None,
    definition_narrow_radius_px: float = 0.0,
) -> CircularLabelLayout | None:
    """The drawn labels' fit by kind, and a probe that places horizontal outer labels with other settings.

    Region labels cut to fit their circle (LP-D87) or lying on another track
    (LP-D86) are a ``circular_region`` kind that only warns, and so are the
    lines of the centre definition text that reach into a track (D-12,
    ``circular_definition``).
    """

    total_length = len(gb_record.seq)
    label_counts: dict[LabelKind, int] = {}
    for label in drawn_labels:
        kind: LabelKind = (
            "circular_embedded"
            if label.get("is_embedded")
            else "circular_inner"
            if label.get("is_inner")
            else "circular_outer"
        )
        label_counts[kind] = label_counts.get(kind, 0) + 1
    placed = bool(label_counts) and place_labels is not None
    cuts = sum(label.cut for label in region_labels)
    if not placed and not cuts and not region_track_overlaps and not definition_track_overlaps:
        return None
    cfg = profile.config
    horizontal = profile.label_placement == "horizontal"

    probe: CircularLabelFitProbe | None = None
    if place_labels is not None and horizontal and label_counts.get("circular_outer"):
        placer = place_labels

        def probe_outer(overrides: Mapping[str, LabelFitValue]) -> LabelFit:
            probe_profile = replace(profile, config=with_config_overrides(cfg, overrides))
            try:
                labels = placer(probe_profile, probe=True)
            except CircularLabelProbeStopped as stopped:
                # The dense candidate already overlaps; crossings are not counted.
                return LabelFit(overlap_pairs=stopped.overlap_pairs, crossing_pairs=0)
            return circular_label_fits(labels, total_length).get("circular_outer", FITS)

        probe = probe_outer

    counts: dict[LabelKind, int] = dict(label_counts) if placed else {}
    kind_fits = circular_label_fits(drawn_labels, total_length, definition_text=definition_text) if placed else {}
    if cuts or region_track_overlaps:
        counts["circular_region"] = len(region_labels)
        kind_fits["circular_region"] = LabelFit(0, 0, region_cuts=cuts, region_track_overlaps=region_track_overlaps)
    if definition_track_overlaps:
        counts["circular_definition"] = definition_track_overlaps
        kind_fits["circular_definition"] = LabelFit(0, 0, definition_track_overlaps=definition_track_overlaps)

    outer_offsets = cfg.labels.unified_adjustment.outer_labels
    return CircularLabelLayout(
        record_id=str(gb_record.id),
        record_index=record_index,
        placement="horizontal" if horizontal else "radial",
        font_size=circular_label_font_size(cfg, total_length),
        label_counts=counts,
        kind_fits=kind_fits,
        x_radius_offset=float(outer_offsets.x_radius_offset),
        y_radius_offset=float(outer_offsets.y_radius_offset),
        label_spacing=float(cfg.labels.spacing.circular),
        refused=radial_refused,
        features_separated=bool(profile.resolve_overlaps),
        probe=probe,
        region_font_size=max((label.font_size for label in region_labels), default=0.0),
        definition_font_size=float(cfg.objects.definition.circular.font_size) if definition_track_overlaps else 0.0,
        definition_probe=definition_probe if definition_track_overlaps else None,
        definition_narrow_radius_px=definition_narrow_radius_px if definition_track_overlaps else 0.0,
    )


def _centre_definition_text(definition: DefinitionGroup) -> tuple[Aabb, ...]:
    """The painted lines of a centred definition in the coordinates of the record's labels.

    ``place_definition_group_on_canvas`` centres the lines' union on the
    record's centre, the origin of the label coordinates.
    """

    bounds = definition.local_bounds
    dx = -0.5 * (bounds.min_x + bounds.max_x)
    dy = -0.5 * (bounds.min_y + bounds.max_y)
    return tuple(box.translated(dx, dy) for box in definition.line_bounds)


def _radial_layout_feature_dict(feature_layers: FeatureBuildResult) -> dict[str, FeatureObject]:
    """The features the radial layout measures: foreground, then underlays."""

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
    return layout_feature_dict


@dataclass(frozen=True)
class _CircularLabelState:
    """What a placement of one record's labels reads besides the label settings.

    The geometry is the one the radial growth left; the fit probes of the
    label-fit check place the labels again with it.
    """

    feature_dict: Mapping[str, FeatureObject]
    total_length: int
    geometry: CircularCanvasGeometry
    strandedness: bool
    track_preset: str
    resolved_feature_anchor_radius_px: float | None
    outer_arena_enabled: bool
    feature_track_ratio_factor_override: float | None
    radial_layout: CircularRadialLayout
    feature_lane_direction: CircularFeatureLaneDirection
    record_transform: RecordDisplayTransform | None
    definition_text: tuple[Aabb, ...]
    label_filtering: Mapping[str, Any]


def _place_circular_labels(
    state: _CircularLabelState,
    label_profile: CircularRenderProfile,
    *,
    probe: bool = False,
    candidate_caches: dict[float, CircularLabelCandidateCache],
) -> list[CircularLabel]:
    """Place one record's labels with ``label_profile`` and the geometry of ``state``."""

    font_size = circular_label_font_size(label_profile.config, state.total_length)
    return prepare_label_list(
        state.feature_dict,
        state.total_length,
        (
            float(state.resolved_feature_anchor_radius_px)
            if state.resolved_feature_anchor_radius_px is not None
            else state.geometry.radius
        ),
        state.geometry.track_ratio,
        label_profile,
        # The arena follows the outer label offsets, which a probe may change.
        outer_arena=(
            _default_outer_label_arena(
                geometry=state.geometry,
                cfg=label_profile.config,
                strandedness=state.strandedness,
            )
            if state.outer_arena_enabled
            else None
        ),
        feature_track_ratio_factor_override=state.feature_track_ratio_factor_override,
        feature_layout=state.radial_layout.features,
        radial_layout=state.radial_layout,
        track_preset=state.track_preset,
        feature_lane_direction=state.feature_lane_direction,
        _candidate_cache=candidate_caches.setdefault(font_size, {}),
        record_transform=state.record_transform,
        probe=probe,
        definition_text=state.definition_text,
        label_filtering=state.label_filtering,
    )


@dataclass(frozen=True)
class CircularLabelPlacement:
    """The drawn labels of one record and the geometry their radial growth left.

    Result of the ``circular_labels`` stage. The caller writes ``state.geometry.radius``
    into its canvas configurator and draws copies of ``labels``.
    """

    labels: tuple[CircularLabel, ...]
    state: _CircularLabelState
    feature_layers: FeatureBuildResult
    layout_slots: tuple[CircularTrackSlot, ...]
    feature_slot: CircularTrackSlot | None
    radial_refused: bool
    # Text measurements per label font size, for the fit probes.
    label_candidates: tuple[tuple[float, tuple[CircularLabelCandidate, ...]], ...]


def circular_labels(
    feature_layers: FeatureBuildResult,
    profile: CircularRenderProfile,
    *,
    total_length: int,
    geometry: CircularCanvasGeometry,
    radial_layout: CircularRadialLayout,
    layout_slots: tuple[CircularTrackSlot, ...],
    feature_slot: CircularTrackSlot | None,
    feature_lane_direction: CircularFeatureLaneDirection,
    circular_preset: str,
    track_preset: str,
    resolved_feature_anchor_radius_px: float | None,
    feature_track_ratio_factor_override: float | None,
    feature_width_override_requested: bool,
    record_features: ResolvedRecordFeatureInputs | None,
    definition_reserved_radius_px: float | None,
    preferred_anchor_slot_ids: Collection[str],
    depth_axis: DepthAxisStyle | None,
    center_reserved_radius_explicit: bool,
    record_transform: RecordDisplayTransform | None,
    definition_text: tuple[Aabb, ...],
    label_filtering: Mapping[str, Any],
    annotation_label_rooms: Mapping[str, AnnotationLabelRoom],
) -> CircularLabelPlacement:
    """Place one record's external and embedded labels; radial labels grow the radius (stage ``circular_labels``).

    The growth loop replans the radial layout with a larger radius until the
    radial labels fit or growing stops making room. It works on copies: the
    feature layers it is given and the canvas configurator are not changed.
    """

    cfg = profile.config
    candidate_caches: dict[float, CircularLabelCandidateCache] = {}

    def state() -> _CircularLabelState:
        return _CircularLabelState(
            feature_dict=feature_layers.foreground_features,
            total_length=total_length,
            geometry=geometry,
            strandedness=profile.strandedness,
            track_preset=track_preset,
            resolved_feature_anchor_radius_px=resolved_feature_anchor_radius_px,
            outer_arena_enabled=feature_width_override_requested,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            radial_layout=radial_layout,
            feature_lane_direction=feature_lane_direction,
            record_transform=record_transform,
            definition_text=definition_text,
            label_filtering=label_filtering,
        )

    slots = list(layout_slots)
    labels = _place_circular_labels(state(), profile, candidate_caches=candidate_caches)
    # Radial labels that collide once the canvas stops growing are refused (OV-694, RADIAL-JOIN).
    radial_refused = False
    if profile.label_placement == "radial":
        previous_growth: float | None = None
        stagnant_passes = 0
        preflight_tracks_frozen = False
        while True:
            required_growth = max(
                (
                    float(label.get("required_radius_growth_px", 0.0))
                    for label in labels
                    if not label.get("is_embedded")
                ),
                default=0.0,
            )
            if required_growth <= RADIAL_GROWTH_TOLERANCE_PX:
                break
            if (
                previous_growth is not None
                and required_growth >= previous_growth - RADIAL_GROWTH_TOLERANCE_PX
            ):
                stagnant_passes += 1
            else:
                stagnant_passes = 0
            if stagnant_passes >= 2:
                # Growing the canvas no longer makes room (LP-D49 RADIAL-JOIN).
                radial_refused = True
                break
            previous_growth = required_growth

            if not preflight_tracks_frozen:
                frozen_slots: list[CircularTrackSlot] = []
                for slot_index, slot in enumerate(slots):
                    resolved_slot = next(
                        (
                            candidate
                            for candidate in radial_layout.slots
                            if int(candidate.slot_index) == int(slot_index)
                        ),
                        None,
                    )
                    # Outside rows are not frozen: they stay beyond the enlarged
                    # axis and the feature slot moved outside, in stack order.
                    if (
                        str(slot.renderer) != "features"
                        and resolved_slot is not None
                        and resolved_slot.anchor_radius_px is not None
                        and resolved_slot.side != "outside"
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
                        # the remaining slot parameters.
                        feature_params = dict(slot.params)
                        feature_params["lane_direction"] = "outside"
                        frozen_slots.append(
                            replace(slot, side="outside", params=feature_params)
                        )
                    else:
                        frozen_slots.append(slot)
                slots = frozen_slots
                feature_slot = next(
                    (
                        slot
                        for slot in slots
                        if slot.enabled and str(slot.renderer) == "features"
                    ),
                    None,
                )
                if feature_slot is not None:
                    feature_lane_direction = _lane_direction_for_feature_slot(
                        feature_slot,
                        track_preset=circular_preset,
                    )
                # Reflow changed the physical slot. Rebuild occupancy through
                # the same planner before measuring any final geometry.
                from ...features.placement import plan_feature_placements

                if feature_slot is not None:
                    # The layers are a stage result: replan on copies of the features.
                    replanned = {
                        feature_id: copy.copy(feature)
                        for feature_id, feature in feature_layers.foreground_features.items()
                    }
                    feature_layers = replace(feature_layers, foreground_features=replanned)
                    plan_feature_placements(
                        replanned,
                        slot=FeaturePlacementSlot(
                            "circular", feature_lane_direction, profile.strandedness,
                        ),
                        record_features=record_features,
                        resolve_overlaps=profile.resolve_overlaps,
                        tolerance_bp=profile.feature_overlap_tolerance_bp,
                        genome_length=total_length,
                    )
                preflight_tracks_frozen = True

            geometry = replace(geometry, radius=float(geometry.radius) + required_growth)
            feature_track_ratio_factor_override = _slot_width_ratio_factor(
                feature_slot,
                base_radius_px=float(geometry.radius),
                base_track_ratio=float(geometry.track_ratio),
            )
            radial_layout = resolve_circular_radial_layout(
                total_length=total_length,
                geometry=geometry,
                cfg=cfg,
                slots=slots,
                feature_dict=_radial_layout_feature_dict(feature_layers),
                definition_reserved_radius_px=definition_reserved_radius_px,
                preferred_anchor_slot_ids=preferred_anchor_slot_ids,
                depth_config=depth_axis,
                center_reserved_radius_explicit=center_reserved_radius_explicit,
                annotation_label_rooms=annotation_label_rooms,
            )
            if radial_layout.features is not None:
                resolved_feature_anchor_radius_px = float(radial_layout.features.anchor_radius_px)
                feature_track_ratio_factor_override = (
                    float(radial_layout.features.width_px)
                    / max(
                        FEATURE_BAND_EPSILON,
                        float(geometry.radius) * float(geometry.track_ratio),
                    )
                )
            labels = _place_circular_labels(state(), profile, candidate_caches=candidate_caches)
        radial_refused = radial_refused or any(
            any(label.get("layout_collision_counts", ()))
            for label in labels
            if not label.get("is_embedded")
        )

    return CircularLabelPlacement(
        labels=tuple(labels),
        state=state(),
        feature_layers=feature_layers,
        layout_slots=tuple(slots),
        feature_slot=feature_slot,
        radial_refused=radial_refused,
        label_candidates=tuple(
            (font_size, cache["candidates"])
            for font_size, cache in candidate_caches.items()
            if "candidates" in cache
        ),
    )


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
    conservation_tracks: Sequence[ConservationTrack] | None = None,
    conservation_min_identity: float = 0.0,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    dinucleotide_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_content_dataframes: dict[str, DataFrame] | None = None,
    dinucleotide_skew_dataframes: dict[str, DataFrame] | None = None,
    definition_position: str = "center",
    definition_profile: CircularDefinitionProfile = "full",
    layout_definition_profile: CircularDefinitionProfile | None = None,
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    annotations: ResolvedAnnotationBundle,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack],
    annotation_record_index: int = 0,
    definition_record_count: int = 1,
    record_transform: RecordDisplayTransform | None = None,
    record_presentation: RecordPresentation | None = None,
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
    record_key = str(record_key_at(gb_record, annotation_record_index))
    profile = canvas_config.profile
    cfg = profile.config
    depth_tracks_by_index = (
        index_depth_track_row(depth_tracks)
        if depth_tracks is not None
        else None
    )

    show_labels_base = profile.labels_enabled
    depth_enabled = bool(profile.show_depth and depth_config is not None)
    simple_show_ticks = bool(cfg.objects.scale.show)
    # The caller planned and laid out the annotation slots (_prepare_circular_annotation_tracks).
    effective_circular_track_slots = circular_track_slots
    resolved_annotations = annotations
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
                geometry=canvas_config.geometry(),
                total_length=len(gb_record.seq),
                strandedness=profile.strandedness,
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                track_ids=tuple(canvas_config.track_ids.items()),
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
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
                geometry=canvas_config.geometry(),
                total_length=len(gb_record.seq),
                strandedness=profile.strandedness,
                show_features=show_features,
                show_ticks=show_ticks_track,
                show_depth=show_depth_track,
                show_gc=show_gc_track,
                show_skew=show_skew_track,
                track_ids=tuple(canvas_config.track_ids.items()),
                dinucleotide=str(getattr(gc_config, "dinucleotide", "GC")),
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
        label_filter_maps(cfg.labels.filtering.as_dict())
        if compute_label_text
        else {}
    )
    feature_layers = render_stage(
        "feature_layers",
        create_feature_layers,
        gb_record,
        feature_config.specific_color_rules,
        feature_config.selected_features_set,
        feature_config.default_color_map,
        profile.strandedness,
        profile.resolve_overlaps,
        label_filtering,
        split_overlaps_by_strand=split_overlaps_by_strand,
        record_features=feature_config.record_features[0] if feature_config.record_features else None,
        placement_slot=(FeaturePlacementSlot(
            "circular", feature_lane_direction, profile.strandedness,
        ) if feature_slot is not None else None),
        feature_overlap_tolerance_bp=profile.feature_overlap_tolerance_bp,
        feature_shapes=feature_config.feature_shapes,
        feature_visibility_rules=feature_config.feature_visibility_rules,
        compute_label_text=compute_label_text,
        record_transform=record_transform,
    )
    layout_feature_dict = _radial_layout_feature_dict(feature_layers)
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
            annotation_track_layouts,
            canvas_config=canvas_config,
            record_index=annotation_record_index,
        )
    resolved_track_slots: list[CircularResolvedSlot] = []
    resolved_feature_anchor_radius_px: float | None = None
    annotation_label_rooms = _circular_annotation_label_rooms(
        annotation_track_layouts,
        record_index=annotation_record_index,
        record_length=len(gb_record.seq),
        font_family=str(cfg.objects.text.font_family),
    )
    # Rows that pack from the feature track are placed around the full definition the record draws alone, so a
    # record summary on a Multi-Record Canvas changes only the center text (OV-1711). The preset lanes keep the
    # drawn definition there (D-08 covers the plot title).
    reserved_definition_profile: CircularDefinitionProfile = (
        "full"
        if definition_profile == "record_summary" and not preferred_anchor_slot_ids
        else definition_profile if layout_definition_profile is None else layout_definition_profile
    )
    # The tracks and the definition text drawn among them, for a config: the figure's own, or a smaller definition
    # font that the definition row probes (D-12).
    def place_around_definition(layout_cfg: GbdrawConfig) -> _CentreDefinitionFit:
        return _resolve_radial_layout_fitting_definition(
            lambda reserved_radius_px: resolve_circular_radial_layout(
                total_length=len(gb_record.seq),
                geometry=canvas_config.geometry(),
                cfg=layout_cfg,
                slots=layout_slots,
                feature_dict=layout_feature_dict,
                definition_reserved_radius_px=reserved_radius_px,
                preferred_anchor_slot_ids=preferred_anchor_slot_ids,
                depth_config=(
                    DepthAxisStyle.of(depth_config) if show_depth_track and depth_config is not None else None
                ),
                center_reserved_radius_explicit=center_reserved_radius is not None,
                annotation_label_rooms=annotation_label_rooms,
            ),
            gb_record=gb_record,
            geometry=canvas_config.geometry(),
            cfg=layout_cfg,
            species=species,
            strain=strain,
            plot_title=plot_title,
            definition_profile=definition_profile,
            layout_definition_profile=reserved_definition_profile,
            definition_position=definition_position,
            center_reserved_radius=center_reserved_radius,
            record_presentation=record_presentation,
        )

    def centre_definition_text(text_cfg: GbdrawConfig) -> tuple[Aabb, ...]:
        if str(definition_position) != "center" or str(definition_profile) == "shared_common":
            return ()
        return _centre_definition_text(
            DefinitionGroup(
                gb_record,
                canvas_config.geometry(),
                style=CircularDefinitionStyle.of(text_cfg, definition_profile),
                species=species,
                strain=strain,
                plot_title=plot_title,
                record_index=int(annotation_record_index),
                record_count=int(definition_record_count),
                record_presentation=record_presentation,
            )
        )

    def definition_lines_on_tracks(
        layout: CircularRadialLayout,
        definition_lines: Sequence[Aabb],
        labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    ) -> int:
        inks = _record_track_ink(
            layout,
            gb_record=gb_record,
            canvas_config=canvas_config,
            record_transform=record_transform,
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
            feature_dict=layout_feature_dict,
            feature_config=feature_config,
            annotation_track_layouts=annotation_track_layouts,
            annotation_labels=labels,
            annotation_record_index=annotation_record_index,
            use_slot_tick_options=user_slot_mode,
        )
        return _definition_lines_on_tracks(inks, definition_lines)

    def definition_clears_tracks(font_size: float) -> bool:
        """A hint for the definition row: the tracks and the text at ``font_size`` before the labels are placed."""

        probe_cfg = with_fit_row_overrides(cfg, {DEFINITION_FONT_PATH: font_size})
        try:
            probe = place_around_definition(probe_cfg)
        except ValidationError:
            return False
        return definition_lines_on_tracks(probe.layout, centre_definition_text(probe_cfg), {}) == 0

    centre_fit = place_around_definition(cfg)
    radial_layout = centre_fit.layout
    definition_reserved_radius_px = centre_fit.reserved_radius_px
    if radial_layout.features is not None:
        resolved_feature_anchor_radius_px = float(radial_layout.features.anchor_radius_px)
        feature_track_ratio_factor_override = (
            float(radial_layout.features.width_px)
            / max(FEATURE_BAND_EPSILON, float(canvas_config.radius) * float(canvas_config.track_ratio))
        )
    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features.
    precalculated_labels: list[CircularLabel] | None = None
    # LP-D49 DEF-WARN and the LP-D57 repair: labels over the centre definition text (not a plot title); D-12: the
    # text over a track.
    definition_text = centre_definition_text(cfg)
    radial_refused = False
    label_state: _CircularLabelState | None = None
    # Text measurements per label font size, shared by the drawn placement and fit probes.
    label_candidate_caches: dict[float, CircularLabelCandidateCache] = {}
    if show_external_labels:
        placement = render_stage(
            "circular_radial_labels" if profile.label_placement == "radial" else "circular_labels",
            circular_labels,
            feature_layers,
            profile,
            total_length=len(gb_record.seq),
            geometry=canvas_config.geometry(),
            radial_layout=radial_layout,
            layout_slots=tuple(layout_slots),
            feature_slot=feature_slot,
            feature_lane_direction=feature_lane_direction,
            circular_preset=circular_preset,
            track_preset=_circular_preset_for_layout(cfg),
            resolved_feature_anchor_radius_px=resolved_feature_anchor_radius_px,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            feature_width_override_requested=feature_width_override_requested,
            record_features=feature_config.record_features[0] if feature_config.record_features else None,
            definition_reserved_radius_px=definition_reserved_radius_px,
            preferred_anchor_slot_ids=preferred_anchor_slot_ids,
            # The radial replan reads the depth axis and ticks only.
            depth_axis=(
                DepthAxisStyle.of(depth_config) if show_depth_track and depth_config is not None else None
            ),
            center_reserved_radius_explicit=center_reserved_radius is not None,
            record_transform=record_transform,
            definition_text=definition_text,
            label_filtering=label_filtering,
            annotation_label_rooms=annotation_label_rooms,
        )
        # The one write of the radius the labels grew (a stage result is never changed).
        canvas_config.radius = placement.state.geometry.radius
        label_state = placement.state
        radial_layout = label_state.radial_layout
        feature_layers = placement.feature_layers
        layout_slots = list(placement.layout_slots)
        feature_slot = placement.feature_slot
        feature_lane_direction = label_state.feature_lane_direction
        resolved_feature_anchor_radius_px = label_state.resolved_feature_anchor_radius_px
        feature_track_ratio_factor_override = label_state.feature_track_ratio_factor_override
        radial_refused = placement.radial_refused
        label_candidate_caches.update(
            (font_size, {"candidates": candidates}) for font_size, candidates in placement.label_candidates
        )
        # Copies: the legend collision shift and the leader pass move them.
        precalculated_labels = [label.copy() for label in placement.labels]

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

    # Annotation labels lie outside the radial bands too: the canvas holds them and a legend or title avoids them.
    # They move along their regions off the feature labels and leaders, so they follow the feature labels.
    def place_annotation_labels() -> tuple[dict[str, tuple[CircularAnnotationLabel, ...]], tuple[Aabb, ...]]:
        labels_by_slot = _circular_annotation_labels_by_slot(
            radial_layout,
            annotation_track_layouts,
            record_index=annotation_record_index,
            record_length=len(gb_record.seq),
            font_family=str(cfg.objects.text.font_family),
            obstacles=_feature_label_obstacles(precalculated_labels or (), len(gb_record.seq), canvas_config),
        )
        return labels_by_slot, tuple(
            label.ink_bounds.translated(float(canvas_config.offset_x), float(canvas_config.offset_y))
            for labels in labels_by_slot.values()
            for label in labels
        )

    # A corner legend must also clear the record and its tracks, even when
    # every label has moved outside and the interior has no label obstacles.
    # The record is an obstacle whatever this record's own legend position:
    # a grid places its legend over records drawn without one, and the Web
    # moves a drawn legend to a corner with these obstacles. The label pass
    # predicts the corner with it before the centre definition is drawn.
    def record_obstacle(definition: Group | None = None) -> Aabb:
        return _collect_circular_primary_bounds(
            radial_layout=radial_layout,
            canvas_config=canvas_config,
            feature_config=feature_config,
            gc_config=gc_config,
            skew_config=skew_config,
            depth_config=depth_config,
            overlay_obstacles=(),
            definition_target=definition,
        )

    annotation_labels, annotation_label_obstacles = place_annotation_labels()
    overlay_obstacles = (
        *_external_label_obstacles_on_canvas(precalculated_labels or (), len(gb_record.seq), canvas_config),
        *annotation_label_obstacles,
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
        legend_measurement = measure_legend(
            legend_config,
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
            record_obstacles=(record_obstacle(),),
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
            annotation_labels, annotation_label_obstacles = place_annotation_labels()
            overlay_obstacles = (
                *_external_label_obstacles_on_canvas(precalculated_labels, len(gb_record.seq), canvas_config),
                *annotation_label_obstacles,
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

    # The labels do not change after this: both label phases take one keyed copy.
    placed_labels: Sequence[CircularLabel] = (
        _pre_keyed(label.copy() for label in precalculated_labels or ()) if show_external_labels else ()
    )
    render_context = CircularRecordRenderContext(
        profile=profile,
        track_preset=circular_preset,
        feature_lane_direction=feature_lane_direction,
        radial_layout=radial_layout,
    )

    canvas = add_axis_group_on_canvas(canvas, canvas_config, group_key=GroupKey(record_key, GroupKind.AXIS_CIRCLE))

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            group_key=GroupKey(record_key, GroupKind.LABEL_LEADERS),
            render_context=render_context,
            precalculated_labels=placed_labels,
            phase="leaders",
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
            use_slot_group_id=user_slot_mode,
            use_slot_tick_options=user_slot_mode,
            use_feature_anchor_override=user_slot_mode,
            annotation_track_layouts=annotation_track_layouts,
            annotation_labels=annotation_labels,
            annotation_record_index=annotation_record_index,
            record_transform=record_transform,
            record_key=record_key,
        )

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            group_key=GroupKey(record_key, GroupKind.LABEL_TEXT),
            render_context=render_context,
            precalculated_labels=placed_labels,
            phase="text",
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
        group_key=GroupKey(record_key, GroupKind.DEFINITION),
        **definition_kwargs,
        record_presentation=record_presentation,
    )
    definition_target = next(
        (
            element
            for element in getattr(canvas, "elements", [])[definition_target_start:]
            if isinstance(element, Group)
        ),
        None,
    )
    # The check counts the labels as drawn: after the legend shift and the last leader pass.
    label_layout = _record_label_layout(
        gb_record,
        profile,
        precalculated_labels or [],
        (
            partial(_place_circular_labels, label_state, candidate_caches=label_candidate_caches)
            if label_state is not None
            else None
        ),
        # The index the record's definition group carries (``definition_kwargs``).
        record_index=int(annotation_record_index),
        definition_text=definition_text,
        radial_refused=radial_refused,
        region_labels=tuple(label for labels in annotation_labels.values() for label in labels),
        region_track_overlaps=_region_labels_on_other_tracks(radial_layout, annotation_labels),
        definition_track_overlaps=(
            definition_lines_on_tracks(radial_layout, definition_text, annotation_labels) if definition_text else 0
        ),
        definition_probe=definition_clears_tracks,
        definition_narrow_radius_px=centre_fit.narrow_radius_px,
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
    overlay_obstacles = (*overlay_obstacles, record_obstacle(definition_target))
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
        canvas = add_legend_group_on_canvas(canvas, legend_measurement, legend_table)
        legend_target = next(
            (
                element
                for element in getattr(canvas, "elements", [])[legend_target_start:]
                if isinstance(element, Group)
            ),
            None,
        )
    track_slot_geometry = getattr(canvas, "_gbdraw_track_slot_geometry", {})
    if feature_slot is not None:
        targets = FeaturePlacementSlot("circular", feature_lane_direction, profile.strandedness).supported_targets()
        for record_geometry in track_slot_geometry.get("records", []):
            record_geometry["featurePlacementTargets"] = targets

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
        label_layout=label_layout,
    )


def _compose_circular_plot(
    plot: _CircularPlotAssembly,
    *,
    title_content: GroupContent | None = None,
    title_placement: TitlePlacement | str = TitlePlacement.NONE,
) -> CircularAssemblyResult:
    """Compose one fully measured Circular plot exactly once."""
    resolved_title_placement = (
        title_placement
        if isinstance(title_placement, TitlePlacement)
        else TitlePlacement(str(title_placement))
    )
    title_bounds = title_content.local_bounds if title_content is not None else None
    if title_content is not None and title_bounds is None:
        raise ValueError("a title content needs its local bounds")
    if title_content is None and resolved_title_placement is not TitlePlacement.NONE:
        raise ValueError("title placement requires a title content")
    if title_content is not None and resolved_title_placement is TitlePlacement.NONE:
        raise ValueError("title content requires a non-none title placement")
    title_target: Group | None = None
    if title_content is not None:
        title_target = place_content(
            plot.drawing, title_content, GroupKey(CANVAS_SCOPE, GroupKind.PLOT_TITLE)
        )
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
    plan = render_stage(
        "composition_plan",
        plan_composition,
        CompositionRequest(
            primary=CompositionItem("primary", plot.source_content_bounds),
            legend=legend_item,
            title=title_item,
            legend_placement=legend_placement,
            title_placement=resolved_title_placement,
            overlay_obstacles=plot.source_overlay_obstacles,
        ),
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
        label_layout=plot.label_layout,
    )
    setattr(plot.drawing, "_gbdraw_circular_assembly_result", result)
    return result


def build_circular_legend_table(
    records: Sequence[SeqRecord],
    *,
    sources: LegendSources,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
    depth: DepthLegendSource | None = None,
    has_depth_frame: bool = False,
    depth_tracks: tuple[DepthLegendSource, ...] | None = None,
    depth_track_count_value: int = 0,
    slots: tuple[LegendSlotSource, ...] | None = None,
    annotation_slots: tuple[LegendSlotSource, ...] | None = None,
    annotations: ResolvedAnnotationBundle | None = None,
    conservation_legend: Sequence[ConservationLegendSource] = (),
    conservation_colors: tuple[str, str] = ("", ""),
    conservation_min_identity: float = 0.0,
) -> LegendTable:
    """Build the one Circular legend table for the records drawn on a canvas.

    Single-record and Multi-Record Canvas output both call this owner, so slot
    renderer parameters, region annotation and depth ``legend_label`` values,
    and conservation rows reach every Circular legend in the same order.
    ``circular_legend_table`` translates the drawing settings into these
    legend-only arguments.
    """

    features_present, used_color_rules, default_used_features = feature_legend_presence(
        records, sources.features
    )
    singleton_depth = bool(
        show_depth
        and depth is not None
        and (has_depth_frame or int(depth_track_count_value) == 1)
    )
    legend_table = prepare_legend_table(
        sources.gc,
        sources.skew,
        sources.features,
        features_present,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
        depth=depth if singleton_depth else None,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=singleton_depth,
    )
    if depth_tracks and slots is None:
        legend_table = sync_depth_track_legend_entries(legend_table, depth_tracks)
    legend_table = _sync_legend_table_for_circular_slots(
        legend_table,
        slots=slots,
        gc=sources.gc,
        skew=sources.skew,
        depth=depth if has_depth_frame else None,
        depth_tracks=depth_tracks,
        conservation_legend=conservation_legend,
        conservation_colors=conservation_colors,
        conservation_min_identity=conservation_min_identity,
    )
    if annotations is not None:
        legend_table = sync_annotation_legend_entries(legend_table, annotations, annotation_slots)
    return legend_table


def circular_legend_table(
    records: Sequence[SeqRecord],
    sources: LegendSources,
    *,
    profile: CircularRenderProfile,
    depth_config: DepthConfigurator | None,
    depth_df: DataFrame | None,
    depth_tracks: Sequence[DepthTrackData] | None,
    depth_track_count_value: int,
    circular_track_slots: list[CircularTrackSlot] | None,
    annotations: ResolvedAnnotationBundle | None,
    conservation_legend: Sequence[ConservationLegendSource],
    conservation_min_identity: float,
) -> LegendTable:
    """Call the ``legend_table`` stage with the legend-only parts of the drawing settings.

    ``circular_track_slots`` is the slot list after annotation planning;
    annotation rows read it normalized.
    """

    conservation = profile.config.objects.conservation
    return render_stage(
        "legend_table",
        build_circular_legend_table,
        list(records),
        sources=sources,
        show_gc=profile.show_gc,
        show_skew=profile.show_skew,
        show_depth=profile.show_depth,
        depth=DepthLegendSource.of_config(depth_config) if depth_config is not None else None,
        has_depth_frame=depth_df is not None,
        depth_tracks=depth_legend_sources(depth_tracks) if depth_tracks is not None else None,
        depth_track_count_value=depth_track_count_value,
        slots=legend_slot_sources(circular_track_slots),
        annotation_slots=(
            legend_slot_sources(normalize_circular_track_slots(circular_track_slots))
            if circular_track_slots is not None and annotations is not None
            else None
        ),
        annotations=annotations,
        conservation_legend=tuple(conservation_legend),
        conservation_colors=(conservation.min_color, conservation.max_color),
        conservation_min_identity=conservation_min_identity,
    )


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
    definition_profile: CircularDefinitionProfile = "full",
    layout_definition_profile: CircularDefinitionProfile | None = None,
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None,
    annotation_record_index: int = 0,
    definition_record_count: int = 1,
    title_content: GroupContent | None = None,
    title_placement: TitlePlacement | str = TitlePlacement.NONE,
    record_transform: RecordDisplayTransform | None = None,
    record_presentation: RecordPresentation | None = None,
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
    effective_circular_track_slots, resolved_annotations, annotation_track_layouts = _prepare_circular_annotation_tracks(
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
        record_transform=record_transform,
    )

    legend_sources = LegendSources.of(feature_config, gc_config, skew_config)
    legend_table = LegendTable()
    legend_measurement = measure_legend(
        legend_config,
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )
    if canvas_config.legend_position != "none":
        legend_table = circular_legend_table(
            [gb_record],
            legend_sources,
            profile=profile,
            depth_config=depth_config,
            depth_df=depth_df,
            depth_tracks=depth_tracks,
            depth_track_count_value=resolved_depth_track_count,
            circular_track_slots=effective_circular_track_slots,
            annotations=resolved_annotations,
            conservation_legend=conservation_legend_sources(conservation_tracks),
            conservation_min_identity=conservation_min_identity,
        )
    canvas: Drawing = canvas_config.create_svg_canvas()
    attach_legend_row_facts(
        canvas,
        LegendRowFacts(
            records=[gb_record],
            sources=legend_sources,
            legend_table=legend_table,
        ),
    )
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
        conservation_tracks=conservation_tracks,
        conservation_min_identity=conservation_min_identity,
        circular_track_slots=effective_circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_dataframes=dinucleotide_dataframes,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position=definition_position,
        definition_profile=definition_profile,
        layout_definition_profile=layout_definition_profile,
        definition_group_id=definition_group_id,
        center_reserved_radius=center_reserved_radius,
        annotations=resolved_annotations,
        annotation_track_layouts=annotation_track_layouts,
        annotation_record_index=annotation_record_index,
        definition_record_count=definition_record_count,
        record_transform=record_transform,
        record_presentation=record_presentation,
    )
    return _compose_circular_plot(
        plot,
        title_content=title_content,
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
    definition_profile: CircularDefinitionProfile = "full",
    definition_group_id: str | None = None,
    center_reserved_radius: float | None = None,
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
        annotations=annotations,
        annotation_record_index=annotation_record_index,
        definition_record_count=definition_record_count,
    )
    return result.drawing, result.legend_measurement


__all__ = [
    "CircularAssemblyResult",
    "assemble_circular_diagram",
    "add_record_on_circular_canvas",
    "build_circular_legend_table",
    "plan_circular_annotation_slots",
]

