"""Unified radial slot resolver for circular diagrams."""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import Any, Collection, Literal, Mapping, Sequence, TypedDict

from ...canvas import CircularCanvasConfigurator
from ...config.models import GbdrawConfig
from ...configurators import DepthConfigurator
from ...exceptions import ValidationError
from ...layout.circular import (
    CircularAxisLayout,
    CircularFeatureLane,
    CircularFeatureLayout,
    CircularFeatureStackMetrics,
    CircularRadialLayout,
    CircularResolvedSlot,
    CircularTickLayout,
    LAYOUT_EPSILON,
    RadialBand,
    band_union,
    bands_overlap,
)
from ...layout.circular_depth_axis import (
    DepthAxisFootprint,
    resolve_depth_axis_footprint,
)
from ...tracks.circular import (
    NUMERIC_CIRCULAR_TRACK_RENDERERS,
    CircularTrackSlot,
    NormalizedCircularTrackSlot,
    normalize_circular_track_slots,
    tick_sides_for_tick_label_layout,
)
from ...svg.circular_ticks import (
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_radius_bounds,
)
from .presets import normalize_circular_track_preset


logger = logging.getLogger(__name__)

TRACK_SPACING_MULTIPLIER = 1.2
MIN_NUMERIC_WIDTH_PX = 10.0
MIN_NUMERIC_WIDTH_FRACTION = 0.55
MIN_SKEW_WIDTH_PX = 12.0
MIN_SKEW_WIDTH_FRACTION = 0.65
MIN_CONSERVATION_WIDTH_PX = 6.0
MIN_CONSERVATION_WIDTH_FRACTION = 0.35
MIN_DENSE_CONSERVATION_WIDTH_PX = 4.0
MIN_AUTO_STACK_SPACING_PX = 1.0
PREFERRED_MIN_NUMERIC_WIDTH_FRACTION = 0.4
PREFERRED_MIN_SKEW_WIDTH_FRACTION = 0.4

PlacementPolicy = Literal["hard", "preferred", "auto", "overlay"]


@dataclass(frozen=True)
class PlacementWindow:
    inner_px: float
    outer_px: float
    # The reserved band (tick labels) may reach below inner_px to this edge of
    # a pinned row, as it may touch the row below it in top-down packing.
    reserved_inner_px: float | None = None


@dataclass(frozen=True)
class _RadialSlotIntent:
    slot: NormalizedCircularTrackSlot
    slot_index: int
    slot_id: str
    renderer: str
    side: str
    anchor_offset_px: float | None
    width_px: float
    explicit_anchor: bool
    explicit_width: bool
    explicit_spacing: bool
    spacing_px: float
    explicit_inner_gap: bool
    explicit_outer_gap: bool
    inner_gap_px: float
    outer_gap_px: float
    z: int
    compress: bool
    reserve: bool
    placement_policy: PlacementPolicy
    params: Mapping[str, Any]


def _cannot_fit_diagnostic(intent: "_RadialSlotIntent | None", window: Any = None) -> dict[str, object]:
    """Track row and usable band of a fit failure; the slot ID stays private.

    A pinned slot has no band to report, and an unbounded band is omitted.
    """

    diagnostic: dict[str, object] = {"code": "TRACK_LAYOUT", "reason": "CANNOT_FIT"}
    if window is not None and float(window.outer_px) < float("inf"):
        diagnostic["innerPx"] = max(0, round(float(window.inner_px)))
        diagnostic["outerPx"] = max(0, round(float(window.outer_px)))
    if intent is not None:
        diagnostic["slotIndex"] = int(intent.slot_index)
    return diagnostic


def _center_reserved_diagnostic(
    diagnostic: Mapping[str, object], *, explicit_radius: bool
) -> dict[str, object]:
    """The same fit failure, limited by the center reservation.

    ``DEFINITION_RESERVED``: the band that the center definition text needs.
    ``CENTER_RESERVED``: an explicit ``center_reserved_radius``.
    """

    if explicit_radius:
        return {**diagnostic, "code": "TRACK_LAYOUT", "reason": "CENTER_RESERVED"}
    return {**diagnostic, "code": "TRACK_LAYOUT", "reason": "DEFINITION_RESERVED"}


def _band_from_center_width(center_px: float, width_px: float) -> RadialBand:
    half = max(0.0, float(width_px) / 2.0)
    return RadialBand(float(center_px) - half, float(center_px) + half)


def _depth_reserved_band_for_draw_band(
    draw_band_px: RadialBand,
    footprint: DepthAxisFootprint | None,
) -> RadialBand:
    if footprint is None:
        return draw_band_px
    return RadialBand(
        float(draw_band_px.inner_px) - float(footprint.radial_inner_extra_px),
        float(draw_band_px.outer_px) + float(footprint.radial_outer_extra_px),
    )


def _lane_direction_from_legacy_track_type(track_type: str | None) -> str:
    preset = normalize_circular_track_preset(track_type)
    if preset == "middle":
        return "split"
    if preset == "spreadout":
        return "outside"
    return "inside"


def _preset_from_lane_direction(lane_direction: str | None) -> str:
    direction = str(lane_direction or "inside").strip().lower()
    if direction == "split":
        return "middle"
    if direction == "outside":
        return "spreadout"
    return "tuckin"


def _feature_track_ids(
    feature_dict: Mapping[str, Any] | None,
    *,
    strandedness: bool = False,
    include_nominal_lanes: bool = False,
) -> tuple[int, ...]:
    if not feature_dict:
        return (
            (0, -1)
            if include_nominal_lanes and strandedness
            else (0,)
        )
    track_ids = {
        int(getattr(feature_object, "feature_track_id", 0))
        for feature_object in feature_dict.values()
    }
    return tuple(sorted(track_ids or {0}))


def _feature_lane_strand_group(track_id: int, strandedness: bool) -> str:
    track = int(track_id)
    if track < 0:
        return "negative"
    if bool(strandedness):
        return "positive"
    if track > 0:
        return "positive"
    return "combined"


def _inside_lane_order(track_ids: Sequence[int], *, strandedness: bool) -> tuple[int, ...]:
    ids = [int(track_id) for track_id in track_ids]
    has_negative = any(track_id < 0 for track_id in ids)
    if bool(strandedness) or has_negative:
        negative = sorted((track_id for track_id in ids if track_id < 0), key=lambda value: abs(value), reverse=True)
        positive = sorted((track_id for track_id in ids if track_id >= 0), reverse=True)
        return tuple(negative + positive)
    return tuple(sorted(ids, reverse=True))


def _outside_lane_order(track_ids: Sequence[int], *, strandedness: bool) -> tuple[int, ...]:
    ids = [int(track_id) for track_id in track_ids]
    has_negative = any(track_id < 0 for track_id in ids)
    if bool(strandedness) or has_negative:
        negative = sorted((track_id for track_id in ids if track_id < 0), key=lambda value: abs(value), reverse=True)
        positive = sorted((track_id for track_id in ids if track_id >= 0), reverse=True)
        return tuple(negative + positive)
    return tuple(sorted(ids))


def _feature_stack_band_width(lane_count: int, lane_width_px: float, lane_spacing_px: float) -> float:
    count = max(1, int(lane_count))
    width = max(0.0, float(lane_width_px))
    spacing = max(0.0, float(lane_spacing_px))
    return (count * width) + ((count - 1) * spacing)


def _single_band_lane_centers(
    *,
    ordered_track_ids_inner_to_outer: Sequence[int],
    center_radius_px: float,
    lane_width_px: float,
    lane_spacing_px: float,
) -> dict[int, float]:
    ordered = tuple(int(track_id) for track_id in ordered_track_ids_inner_to_outer) or (0,)
    band_width = _feature_stack_band_width(len(ordered), lane_width_px, lane_spacing_px)
    step = max(0.0, float(lane_width_px)) + max(0.0, float(lane_spacing_px))
    first_center = float(center_radius_px) - (0.5 * band_width) + (0.5 * max(0.0, float(lane_width_px)))
    return {
        int(track_id): max(0.0, first_center + (idx * step))
        for idx, track_id in enumerate(ordered)
    }


def measure_circular_feature_stack(
    *,
    axis_radius_px: float,
    lane_width_px: float,
    lane_spacing_px: float | None = None,
    preset: str | None = None,
    lane_direction: str | None = None,
    strandedness: bool,
    track_ids: Sequence[int] = (0,),
    center_radius_px: float | None = None,
    placements_by_track_id: Mapping[int, Any] | None = None,
) -> CircularFeatureStackMetrics:
    """Measure feature-lane centers from the resolved feature-stack preset."""

    width = max(0.0, float(lane_width_px))
    spacing = (
        _default_spacing_px(axis_radius_px)
        if lane_spacing_px is None
        else max(0.0, float(lane_spacing_px))
    )
    axis = float(axis_radius_px)
    ids = tuple(sorted({int(track_id) for track_id in track_ids})) or (0,)
    lane_count = len(ids)
    band_width = _feature_stack_band_width(lane_count, width, spacing)
    normalized_preset = (
        normalize_circular_track_preset(preset)
        if preset is not None
        else normalize_circular_track_preset(_preset_from_lane_direction(lane_direction))
    )
    center = float(center_radius_px) if center_radius_px is not None else axis

    if center_radius_px is None:
        if normalized_preset == "tuckin":
            center = axis - spacing - (0.5 * band_width)
        elif normalized_preset == "spreadout":
            center = axis + spacing + (0.5 * band_width)

    lane_centers: dict[int, float]
    if normalized_preset == "middle" and placements_by_track_id and not strandedness:
        # Semantic levels reserve the empty Main band, even for a lone lane 1.
        # Track IDs only associate the existing geometry consumers with this result.
        step = width + spacing
        main_center = center
        lane_centers = {
            track_id: max(0.0, (
                main_center - (assignment.level * step)
                if assignment.side == "inward"
                else main_center + (assignment.level * step)
            ))
            for track_id, assignment in placements_by_track_id.items()
        }
        lane_centers.setdefault(0, main_center)
    elif normalized_preset == "middle":
        step = width + spacing
        ids_set = set(ids)
        has_negative = any(track_id < 0 for track_id in ids)
        if bool(strandedness) or has_negative:
            lane_centers = {}
            split_seed = (0.5 * width) + (0.5 * spacing)
            for track_id in (tid for tid in ids if tid < 0):
                lane_centers[int(track_id)] = max(0.0, center - split_seed - ((abs(track_id) - 1) * step))
            for track_id in (tid for tid in ids if tid >= 0):
                lane_centers[int(track_id)] = max(0.0, center + split_seed + (track_id * step))
            if not lane_centers:
                lane_centers[0] = max(0.0, center)
        else:
            lane_centers = {}
            anchor_id = 0 if 0 in ids_set else ids[0]
            lane_centers[int(anchor_id)] = max(0.0, center)
            for track_id in ids:
                if track_id == anchor_id:
                    continue
                if track_id < 0:
                    lane_centers[int(track_id)] = max(0.0, center - (abs(track_id) * step))
                else:
                    lane_centers[int(track_id)] = max(0.0, center + (abs(track_id) * step))
    else:
        ordered = (
            _outside_lane_order(ids, strandedness=strandedness)
            if normalized_preset == "spreadout"
            else _inside_lane_order(ids, strandedness=strandedness)
        )
        lane_centers = _single_band_lane_centers(
            ordered_track_ids_inner_to_outer=ordered,
            center_radius_px=center,
            lane_width_px=width,
            lane_spacing_px=spacing,
        )

    lane_bands = [
        RadialBand(center_px - (0.5 * width), center_px + (0.5 * width))
        for center_px in lane_centers.values()
    ]
    if normalized_preset == "middle" and placements_by_track_id and not strandedness:
        band_width = (band_union(lane_bands) or RadialBand(center, center)).width_px
    all_band = band_union(lane_bands) or RadialBand(center, center)
    return CircularFeatureStackMetrics(
        lane_count=lane_count,
        lane_width_px=width,
        lane_spacing_px=spacing,
        band_width_px=band_width,
        center_radius_px=float(center),
        inner_radius_px=float(all_band.inner_px),
        outer_radius_px=float(all_band.outer_px),
        lane_centers_by_track_id=lane_centers,
    )


def build_circular_feature_layout(
    feature_dict: Mapping[str, Any] | None,
    *,
    axis_radius_px: float,
    width_px: float,
    track_type: str | None = None,
    lane_direction: str | None = None,
    strandedness: bool,
    anchor_radius_px: float | None = None,
    lane_spacing_px: float | None = None,
    include_nominal_lanes: bool = False,
) -> CircularFeatureLayout | None:
    width = max(0.0, float(width_px))
    direction = (
        str(lane_direction).strip().lower()
        if lane_direction is not None
        else _lane_direction_from_legacy_track_type(track_type)
    )
    track_ids = _feature_track_ids(
        feature_dict,
        strandedness=bool(strandedness),
        include_nominal_lanes=include_nominal_lanes,
    )
    placements = {
        int(feature.feature_track_id): feature.placement
        for feature in (feature_dict or {}).values()
        if getattr(feature, "placement", None) is not None
    }
    metrics = measure_circular_feature_stack(
        axis_radius_px=float(axis_radius_px),
        lane_width_px=width,
        lane_spacing_px=lane_spacing_px,
        preset=track_type if lane_direction is None else None,
        lane_direction=direction,
        strandedness=bool(strandedness),
        track_ids=track_ids,
        center_radius_px=anchor_radius_px,
        placements_by_track_id=placements,
    )
    anchor = float(metrics.center_radius_px)

    lanes: dict[int, CircularFeatureLane] = {}
    for track_id in sorted(track_ids):
        center = float(metrics.lane_centers_by_track_id[int(track_id)])
        assignment = placements.get(int(track_id))
        strand_group = assignment.strand_pool if assignment is not None else _feature_lane_strand_group(int(track_id), bool(strandedness))
        half_width = width / 2.0
        lanes[int(track_id)] = CircularFeatureLane(
            track_id=int(track_id),
            strand_group=strand_group,
            inner_px=max(0.0, float(center) - half_width),
            center_px=max(0.0, float(center)),
            outer_px=max(0.0, float(center) + half_width),
        )

    all_band = RadialBand(metrics.inner_radius_px, metrics.outer_radius_px)
    primary_lanes = [lane.band_px for tid, lane in lanes.items() if tid in {0, -1}]
    primary_band = (
        RadialBand(
            metrics.lane_centers_by_track_id[0] - 0.5 * width,
            metrics.lane_centers_by_track_id[0] + 0.5 * width,
        )
        if placements and direction == "split" and not strandedness
        else band_union(primary_lanes) or all_band
    )
    return CircularFeatureLayout(
        anchor_radius_px=anchor,
        width_px=width,
        lanes_by_track_id=lanes,
        primary_band_px=primary_band,
        all_band_px=all_band,
    )


def _tick_layout_from_params(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    params: Mapping[str, Any],
    anchor_radius_px: float,
    width_px: float,
    explicit_width: bool,
    tick_track_channel_override: str | None = None,
) -> CircularTickLayout:
    base_radius = float(canvas_config.radius)
    tick_preset = normalize_circular_track_preset(str(params.get("preset", params.get("track_preset", "tuckin"))))
    tick_label_layout = str(params.get("tick_label_layout", "label_out_tick_in")).strip().lower()
    label_side, tick_side = tick_sides_for_tick_label_layout(
        tick_label_layout,
        side=params.get("_slot_side"),
    )
    tick_length_px = float(width_px) if explicit_width and width_px > 0 else None

    if tick_side in {"none", ""}:
        tick_band = RadialBand(anchor_radius_px, anchor_radius_px)
    else:
        tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
            center_radius_px=float(anchor_radius_px),
            total_len=int(total_length),
            size="large",
            track_type=tick_preset,
            strandedness=canvas_config.profile.strandedness,
            tick_track_channel_override=tick_track_channel_override,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            length_reference_radius_px=base_radius,
        )
        tick_band = RadialBand(tick_inner, tick_outer)

    label_band: RadialBand | None = None
    label_bounds = get_circular_tick_label_radius_bounds(
        center_radius_px=float(anchor_radius_px),
        total_len=int(total_length),
        track_type=tick_preset,
        strandedness=canvas_config.profile.strandedness,
        font_size=float(cfg.objects.ticks.tick_labels.font_size),
        font_family=str(cfg.objects.text.font_family),
        dpi=int(canvas_config.dpi),
        manual_interval=cfg.objects.scale.interval,
        tick_track_channel_override=tick_track_channel_override,
        label_side=label_side,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
        tick_width=float(cfg.objects.ticks.tick_width),
        length_reference_radius_px=base_radius,
    )
    if label_bounds is not None:
        label_band = RadialBand(*label_bounds)
    reserved = band_union([band for band in (tick_band, label_band) if band is not None]) or tick_band
    return CircularTickLayout(
        anchor_radius_px=float(anchor_radius_px),
        tick_band_px=tick_band,
        label_band_px=label_band,
        reserved_band_px=reserved,
        track_preset=tick_preset,
        label_side=label_side,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
    )


def _default_width_px(
    renderer: str,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    length_param = str(canvas_config.length_param)
    base = float(canvas_config.radius) * float(canvas_config.track_ratio)
    if renderer == "features":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    if renderer == "sequence_conservation":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    if renderer == "depth":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][1]) * 0.5
    if renderer == "dinucleotide_skew":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][2])
    if renderer == "ticks":
        return 0.0
    return base * float(cfg.canvas.circular.track_ratio_factors[length_param][1])


def _default_spacing_px(axis_radius_px: float) -> float:
    return max(1.0, 0.01 * float(axis_radius_px))


def _slot_intents(
    slots: Sequence[CircularTrackSlot],
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    preferred_anchor_slot_ids: Collection[str] = (),
) -> list[_RadialSlotIntent]:
    axis_radius_px = float(canvas_config.radius)
    preferred_ids = {str(slot_id) for slot_id in preferred_anchor_slot_ids}
    intents: list[_RadialSlotIntent] = []
    for slot in normalize_circular_track_slots(slots):
        radius_px = slot.radius.resolve(axis_radius_px) if slot.radius is not None else None
        preferred_anchor_px = (
            slot.preferred_anchor_radius.resolve(axis_radius_px)
            if (
                radius_px is None
                and slot.preferred_anchor_radius is not None
                and slot.id in preferred_ids
                and slot.side == "inside"
                and slot.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
            )
            else None
        )
        width_px = slot.width.resolve(axis_radius_px) if slot.width is not None else None
        if width_px is None:
            width_px = _default_width_px(slot.renderer, canvas_config=canvas_config, cfg=cfg)
        default_gap_px = _default_spacing_px(axis_radius_px)
        spacing_px = (
            float(slot.legacy_spacing.resolve(axis_radius_px))
            if slot.legacy_spacing is not None
            else default_gap_px
        )
        inner_gap_px = (
            float(slot.inner_gap_px)
            if slot.inner_gap_px is not None
            else spacing_px
        )
        outer_gap_px = (
            float(slot.outer_gap_px)
            if slot.outer_gap_px is not None
            else spacing_px
        )
        legacy_spacing_px = max(float(inner_gap_px), float(outer_gap_px))
        explicit_anchor = radius_px is not None
        preferred_anchor_available = radius_px is not None or preferred_anchor_px is not None
        placement_policy: PlacementPolicy
        if slot.side == "overlay":
            placement_policy = "overlay"
        elif (
            slot.id in preferred_ids
            and slot.side == "inside"
            and slot.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
            and preferred_anchor_available
        ):
            placement_policy = "preferred"
        elif explicit_anchor:
            placement_policy = "hard"
        else:
            if slot.id in preferred_ids:
                logger.debug(
                    "Ignoring preferred-anchor intent for circular slot '%s' without an inside numeric/depth radius.",
                    slot.id,
                )
            placement_policy = "auto"
        anchor_px = radius_px if radius_px is not None else preferred_anchor_px
        intents.append(
            _RadialSlotIntent(
                slot=slot,
                slot_index=slot.slot_index,
                slot_id=slot.id,
                renderer=slot.renderer,
                side=slot.side,
                anchor_offset_px=(
                    float(anchor_px) - axis_radius_px
                    if anchor_px is not None
                    else None
                ),
                width_px=max(0.0, float(width_px)),
                explicit_anchor=explicit_anchor,
                explicit_width=slot.width is not None,
                explicit_spacing=slot.legacy_spacing is not None,
                spacing_px=max(0.0, float(legacy_spacing_px)),
                explicit_inner_gap=slot.legacy_spacing is not None or slot.inner_gap_px is not None,
                explicit_outer_gap=slot.legacy_spacing is not None or slot.outer_gap_px is not None,
                inner_gap_px=max(0.0, float(inner_gap_px)),
                outer_gap_px=max(0.0, float(outer_gap_px)),
                z=slot.z,
                compress=slot.compress,
                reserve=slot.reserve,
                placement_policy=placement_policy,
                params=slot.params,
            )
        )
    return intents


def _measure_radial_slot(
    intent: _RadialSlotIntent,
    *,
    anchor_offset_px: float,
    width_px: float | None = None,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
    compressed: bool = False,
) -> CircularResolvedSlot:
    anchor_radius_px = float(axis_radius_px) + float(anchor_offset_px)
    resolved_width = max(0.0, float(intent.width_px if width_px is None else width_px))
    renderer = intent.renderer

    if renderer == "features":
        feature_preset = intent.params.get("stack_preset", intent.params.get("preset"))
        feature_layout = build_circular_feature_layout(
            feature_dict,
            axis_radius_px=float(axis_radius_px),
            width_px=resolved_width,
            track_type=str(feature_preset) if feature_preset is not None else None,
            lane_direction=str(intent.params.get("lane_direction", "inside")),
            strandedness=canvas_config.profile.strandedness,
            anchor_radius_px=anchor_radius_px,
            lane_spacing_px=float(intent.spacing_px),
        )
        if feature_layout is not None:
            anchor_radius_px = float(feature_layout.anchor_radius_px)
            anchor_offset_px = float(anchor_radius_px) - float(axis_radius_px)
        band = feature_layout.all_band_px if feature_layout is not None else RadialBand(anchor_radius_px, anchor_radius_px)
        draw_band: RadialBand | None = (
            feature_layout.primary_band_px if feature_layout is not None else band
        )
        return CircularResolvedSlot(
            slot_index=int(intent.slot_index),
            id=intent.slot_id,
            renderer=renderer,
            side=intent.side,
            z=intent.z,
            anchor_radius_px=anchor_radius_px,
            anchor_offset_px=float(anchor_offset_px),
            requested_width_px=float(intent.width_px),
            resolved_width_px=resolved_width,
            packing_band_px=band,
            draw_band_px=draw_band,
            reserved_band_px=band,
            inner_gap_px=float(intent.inner_gap_px),
            outer_gap_px=float(intent.outer_gap_px),
            params=dict(intent.params),
            payload=feature_layout,
            explicit_anchor=bool(intent.explicit_anchor),
            explicit_width=bool(intent.explicit_width),
            compressed=bool(compressed),
        )

    if renderer == "ticks":
        tick_params = dict(intent.params)
        tick_params["_slot_side"] = intent.side
        tick_layout = _tick_layout_from_params(
            total_length=int(total_length),
            canvas_config=canvas_config,
            cfg=cfg,
            params=tick_params,
            anchor_radius_px=anchor_radius_px,
            width_px=resolved_width,
            explicit_width=bool(intent.explicit_width),
            tick_track_channel_override=tick_track_channel_override,
        )
        return CircularResolvedSlot(
            slot_index=int(intent.slot_index),
            id=intent.slot_id,
            renderer=renderer,
            side=intent.side,
            z=intent.z,
            anchor_radius_px=anchor_radius_px,
            anchor_offset_px=float(anchor_offset_px),
            requested_width_px=float(intent.width_px),
            resolved_width_px=resolved_width,
            packing_band_px=tick_layout.tick_band_px,
            draw_band_px=tick_layout.tick_band_px,
            reserved_band_px=tick_layout.reserved_band_px,
            inner_gap_px=float(intent.inner_gap_px),
            outer_gap_px=float(intent.outer_gap_px),
            params=dict(intent.params),
            payload=tick_layout,
            explicit_anchor=bool(intent.explicit_anchor),
            explicit_width=bool(intent.explicit_width),
            compressed=bool(compressed),
        )

    band = _band_from_center_width(anchor_radius_px, resolved_width)
    draw_band = None if renderer == "spacer" else band
    reserved_band = band
    if renderer == "depth" and depth_config is not None:
        reserved_band = _depth_reserved_band_for_draw_band(
            band,
            resolve_depth_axis_footprint(depth_config, resolved_width),
        )
    return CircularResolvedSlot(
        slot_index=int(intent.slot_index),
        id=intent.slot_id,
        renderer=renderer,
        side=intent.side,
        z=intent.z,
        anchor_radius_px=anchor_radius_px,
        anchor_offset_px=float(anchor_offset_px),
        requested_width_px=float(intent.width_px),
        resolved_width_px=resolved_width,
        packing_band_px=band,
        draw_band_px=draw_band,
        reserved_band_px=reserved_band,
        inner_gap_px=float(intent.inner_gap_px),
        outer_gap_px=float(intent.outer_gap_px),
        params=dict(intent.params),
        payload=None,
        explicit_anchor=bool(intent.explicit_anchor),
        explicit_width=bool(intent.explicit_width),
        compressed=bool(compressed),
    )


def _reserved_overlap_any(band: RadialBand, occupied: Sequence[tuple[str, RadialBand]]) -> tuple[str, RadialBand] | None:
    for owner, other in occupied:
        if bands_overlap(band, other):
            return owner, other
    return None


def _min_readable_numeric_width_px(renderer: str, default_width_px: float) -> float:
    width = max(0.0, float(default_width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if renderer == "sequence_conservation":
        return min(width, max(MIN_CONSERVATION_WIDTH_PX, MIN_CONSERVATION_WIDTH_FRACTION * width))
    if renderer == "dinucleotide_skew":
        return min(width, max(MIN_SKEW_WIDTH_PX, MIN_SKEW_WIDTH_FRACTION * width))
    return min(width, max(MIN_NUMERIC_WIDTH_PX, MIN_NUMERIC_WIDTH_FRACTION * width))


def _min_readable_preferred_width_px(renderer: str, default_width_px: float) -> float:
    width = max(0.0, float(default_width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if renderer == "sequence_conservation":
        return min(width, max(MIN_CONSERVATION_WIDTH_PX, PREFERRED_MIN_NUMERIC_WIDTH_FRACTION * width))
    if renderer == "dinucleotide_skew":
        return min(width, max(MIN_SKEW_WIDTH_PX, PREFERRED_MIN_SKEW_WIDTH_FRACTION * width))
    return min(width, max(MIN_NUMERIC_WIDTH_PX, PREFERRED_MIN_NUMERIC_WIDTH_FRACTION * width))


def _min_dense_stack_numeric_width_px(renderer: str, default_width_px: float) -> float:
    width = max(0.0, float(default_width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if renderer == "sequence_conservation":
        return min(width, MIN_DENSE_CONSERVATION_WIDTH_PX)
    if renderer == "dinucleotide_skew":
        return min(width, MIN_SKEW_WIDTH_PX)
    return min(width, MIN_NUMERIC_WIDTH_PX)


def _candidate_widths(intent: _RadialSlotIntent) -> list[tuple[float, bool]]:
    width = max(0.0, float(intent.width_px))
    if not intent.compress or intent.renderer not in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        return [(width, False)]
    min_width = _min_readable_numeric_width_px(intent.renderer, width)
    if min_width >= width - LAYOUT_EPSILON:
        return [(width, False)]
    values = [width]
    steps = 8
    for idx in range(1, steps + 1):
        fraction = idx / float(steps)
        values.append(width - ((width - min_width) * fraction))
    return [(max(0.0, value), idx > 0) for idx, value in enumerate(values)]


def _slot_reserves(intent: _RadialSlotIntent) -> bool:
    if intent.renderer == "annotations" and intent.side == "overlay":
        return False
    return True


def _place_outside_auto(
    intent: _RadialSlotIntent,
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> CircularResolvedSlot:
    for width_px, compressed in _candidate_widths(intent):
        anchor_offset = max(0.0, float(intent.anchor_offset_px or 0.0))
        for _ in range(256):
            resolved = _measure_radial_slot(
                intent,
                anchor_offset_px=anchor_offset,
                width_px=width_px,
                axis_radius_px=axis_radius_px,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=total_length,
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
                compressed=compressed,
            )
            if resolved.packing_band_px is not None and resolved.packing_band_px.inner_px < placement_window.inner_px - LAYOUT_EPSILON:
                anchor_offset += float(placement_window.inner_px) - float(resolved.packing_band_px.inner_px)
                continue
            if (
                resolved.packing_band_px is not None
                and resolved.packing_band_px.outer_px > placement_window.outer_px + LAYOUT_EPSILON
            ):
                break
            if resolved.reserved_band_px is not None:
                conflict = _reserved_overlap_any(resolved.reserved_band_px, occupied)
                if conflict is not None:
                    _owner, band = conflict
                    anchor_offset += (
                        float(band.outer_px)
                        - float(resolved.reserved_band_px.inner_px)
                        + max(0.0, float(intent.inner_gap_px))
                    )
                    continue
            return resolved
    raise ValidationError(
        f"Circular track slot '{intent.slot_id}' cannot be placed outside without overlap.",
        diagnostic=_cannot_fit_diagnostic(intent, placement_window),
    )


def _free_intervals(
    occupied: Sequence[tuple[str, RadialBand]],
    *,
    inner_limit_px: float,
    outer_limit_px: float,
) -> list[tuple[float, float]]:
    lower_limit = max(0.0, float(inner_limit_px))
    upper_limit = max(lower_limit, float(outer_limit_px))
    intervals: list[tuple[float, float]] = []
    cursor = lower_limit
    for _owner, band in sorted(occupied, key=lambda item: item[1].inner_px):
        if band.outer_px <= lower_limit + LAYOUT_EPSILON:
            continue
        if band.inner_px >= upper_limit - LAYOUT_EPSILON:
            break
        band_inner = max(lower_limit, float(band.inner_px))
        band_outer = min(upper_limit, float(band.outer_px))
        if band_inner > cursor + LAYOUT_EPSILON:
            intervals.append((cursor, band_inner))
        cursor = max(cursor, band_outer)
    if cursor < upper_limit - LAYOUT_EPSILON:
        intervals.append((cursor, upper_limit))
    return intervals


def _try_measure_inside_interval(
    intent: _RadialSlotIntent,
    *,
    width_px: float,
    interval: tuple[float, float],
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
    compressed: bool,
    packing_inner_px: float = 0.0,
) -> CircularResolvedSlot | None:
    lower, upper = float(interval[0]), float(interval[1])
    if intent.renderer == "features":
        seeds = [
            upper,
            upper - (0.5 * float(width_px)),
            (lower + upper) / 2.0,
            lower + (0.5 * float(width_px)),
        ]
    elif intent.renderer == "ticks" and upper >= float(axis_radius_px) - LAYOUT_EPSILON:
        seeds = [
            upper,
            upper - (0.5 * float(width_px)),
            (lower + upper) / 2.0,
            lower + (0.5 * float(width_px)),
        ]
    else:
        seeds = [
            upper - (0.5 * float(width_px)),
            (lower + upper) / 2.0,
            upper,
            lower + (0.5 * float(width_px)),
        ]
    for seed_radius in seeds:
        anchor_offset = float(seed_radius) - float(axis_radius_px)
        for _ in range(32):
            resolved = _measure_radial_slot(
                intent,
                anchor_offset_px=anchor_offset,
                width_px=width_px,
                axis_radius_px=axis_radius_px,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=total_length,
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
                compressed=compressed,
            )
            bands = [
                band
                for band in (resolved.packing_band_px, resolved.reserved_band_px)
                if band is not None
            ]
            if not bands:
                return resolved
            min_inner = min(float(band.inner_px) for band in bands)
            if min_inner < lower - LAYOUT_EPSILON:
                anchor_offset += lower - min_inner
                continue
            packing = resolved.packing_band_px
            if packing is not None and float(packing.inner_px) < packing_inner_px - LAYOUT_EPSILON:
                anchor_offset += packing_inner_px - float(packing.inner_px)
                continue
            max_outer = max(float(band.outer_px) for band in bands)
            if max_outer > upper + LAYOUT_EPSILON:
                anchor_offset -= max_outer - upper
                continue
            if resolved.reserved_band_px is not None and _reserved_overlap_any(resolved.reserved_band_px, occupied) is not None:
                break
            return resolved
    return None


def _place_inside_auto(
    intent: _RadialSlotIntent,
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> CircularResolvedSlot:
    for width_px, compressed in _candidate_widths(intent):
        resolved = _place_inside_auto_fixed_width(
            intent,
            width_px=width_px,
            compressed=compressed,
            occupied=occupied,
            axis_radius_px=axis_radius_px,
            placement_window=placement_window,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
        )
        if resolved is not None:
            return resolved
    raise _InsideFitError(
        f"Circular track slot '{intent.slot_id}'",
        intent=intent,
        placement_window=placement_window,
    )


def _place_inside_auto_fixed_width(
    intent: _RadialSlotIntent,
    *,
    width_px: float,
    compressed: bool,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> CircularResolvedSlot | None:
    outer_limit = max(0.0, float(placement_window.outer_px))
    reserved_inner = placement_window.reserved_inner_px
    intervals = _free_intervals(
        occupied,
        inner_limit_px=float(placement_window.inner_px if reserved_inner is None else reserved_inner),
        outer_limit_px=outer_limit,
    )
    for interval in sorted(intervals, key=lambda item: item[1], reverse=True):
        resolved = _try_measure_inside_interval(
            intent,
            width_px=width_px,
            interval=interval,
            occupied=occupied,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
            compressed=compressed,
            packing_inner_px=float(placement_window.inner_px),
        )
        if resolved is None or resolved.reserved_band_px is None:
            continue
        if _reserved_overlap_any(resolved.reserved_band_px, occupied) is None:
            return resolved
    return None


def _measure_anchored_inside(
    intent: _RadialSlotIntent,
    *,
    width_px: float,
    compressed: bool,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> CircularResolvedSlot | None:
    """A pinned row at ``width_px`` on its radius, or ``None`` when it leaves the window or overlaps."""

    resolved = _measure_radial_slot(
        intent,
        anchor_offset_px=float(intent.anchor_offset_px or 0.0),
        width_px=width_px,
        axis_radius_px=axis_radius_px,
        feature_dict=feature_dict,
        canvas_config=canvas_config,
        cfg=cfg,
        total_length=total_length,
        tick_track_channel_override=tick_track_channel_override,
        depth_config=depth_config,
        compressed=compressed,
    )
    bands = [band for band in (resolved.packing_band_px, resolved.reserved_band_px) if band is not None]
    reserved_inner = placement_window.inner_px if placement_window.reserved_inner_px is None else placement_window.reserved_inner_px
    if any(float(band.outer_px) > float(placement_window.outer_px) + LAYOUT_EPSILON for band in bands):
        return None
    if any(float(band.inner_px) < float(reserved_inner) - LAYOUT_EPSILON for band in bands):
        return None
    packing = resolved.packing_band_px
    if packing is not None and float(packing.inner_px) < float(placement_window.inner_px) - LAYOUT_EPSILON:
        return None
    if resolved.reserved_band_px is not None and _reserved_overlap_any(resolved.reserved_band_px, occupied) is not None:
        return None
    return resolved


def _shrinkable_inside_numeric(intent: _RadialSlotIntent) -> bool:
    """An inside numeric row with an Auto width, pinned or not."""

    return (
        intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
        and intent.side == "inside"
        and intent.compress
    )


def _anchored_inside(intent: _RadialSlotIntent) -> bool:
    """A pinned row that compresses like Auto, centred on its radius (GX-17).

    It is placed with the unpinned rows around it, at the same width scale. A
    row pinned at or beyond the Axis lies outside Auto's inside window and keeps
    its width.
    """

    return (
        intent.placement_policy == "hard"
        and _shrinkable_inside_numeric(intent)
        and float(intent.anchor_offset_px or 0.0) < 0.0
    )


def _inside_auto_stack_group_from(
    ordered_intents: Sequence[_RadialSlotIntent],
    start_pos: int,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> list[_RadialSlotIntent]:
    group: list[_RadialSlotIntent] = []
    for future in ordered_intents[start_pos:]:
        if future.slot_index in resolved_by_index:
            break
        if (
            future.side != "inside"
            or future.placement_policy != "auto"
            or future.explicit_anchor
        ):
            break
        group.append(future)
    return group


def _outside_auto_stack_group_from(
    ordered_intents: Sequence[_RadialSlotIntent],
    start_pos: int,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> list[_RadialSlotIntent]:
    group: list[_RadialSlotIntent] = []
    for future in ordered_intents[start_pos:]:
        if future.slot_index in resolved_by_index:
            break
        if (
            future.side != "outside"
            or future.placement_policy != "auto"
            or future.explicit_anchor
        ):
            break
        group.append(future)
    return group


def _inside_movable_stack_group_from(
    ordered_intents: Sequence[_RadialSlotIntent],
    start_pos: int,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> list[_RadialSlotIntent]:
    group: list[_RadialSlotIntent] = []
    for future in ordered_intents[start_pos:]:
        if future.slot_index in resolved_by_index:
            break
        if future.side != "inside" or (future.explicit_anchor and not _anchored_inside(future)):
            break
        if future.placement_policy == "auto" or _anchored_inside(future):
            group.append(future)
            continue
        if future.placement_policy == "preferred" and future.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
            group.append(future)
            continue
        break
    return group


def _inside_auto_stack_width_scales(intents: Sequence[_RadialSlotIntent]) -> list[float]:
    min_scales = [
        _min_dense_stack_numeric_width_px(intent.renderer, intent.width_px) / intent.width_px
        for intent in intents
        if _shrinkable_inside_numeric(intent) and intent.width_px > LAYOUT_EPSILON
    ]
    if not min_scales:
        return [1.0]
    return _linear_scales_from_1_to_min(min(min_scales), steps=16)


def _scaled_inside_auto_width(intent: _RadialSlotIntent, scale: float) -> tuple[float, bool]:
    width = max(0.0, float(intent.width_px))
    if not _shrinkable_inside_numeric(intent) or width <= LAYOUT_EPSILON:
        return width, False
    min_width = _min_dense_stack_numeric_width_px(intent.renderer, width)
    scaled_width = max(min_width, width * float(scale))
    return scaled_width, scaled_width < width - LAYOUT_EPSILON


def _scaled_inside_auto_gap_px(
    intent: _RadialSlotIntent,
    *,
    gap_px: float,
    explicit_gap: bool,
    scale: float,
) -> float:
    spacing = max(0.0, float(gap_px))
    if (
        spacing <= LAYOUT_EPSILON
        or bool(explicit_gap)
        or not _shrinkable_inside_numeric(intent)
    ):
        return spacing
    return max(min(spacing, MIN_AUTO_STACK_SPACING_PX), spacing * float(scale))


def _scaled_inside_auto_inner_gap(intent: _RadialSlotIntent, scale: float) -> float:
    return _scaled_inside_auto_gap_px(
        intent,
        gap_px=float(intent.inner_gap_px),
        explicit_gap=bool(intent.explicit_inner_gap),
        scale=scale,
    )


def _scaled_inside_auto_outer_gap(intent: _RadialSlotIntent, scale: float) -> float:
    return _scaled_inside_auto_gap_px(
        intent,
        gap_px=float(intent.outer_gap_px),
        explicit_gap=bool(intent.explicit_outer_gap),
        scale=scale,
    )




def _gap_between_inner_outer_tracks(
    inner_intent: _RadialSlotIntent,
    outer_intent: _RadialSlotIntent,
    *,
    inner_scale: float = 1.0,
    outer_scale: float = 1.0,
) -> float:
    """Facing gap of two adjacent rows: the inner row's outer gap or the outer row's inner gap, whichever is larger.

    Every placement and the order check use this rule, also after a pinned row
    and across a group boundary.
    """

    return max(
        _scaled_inside_auto_outer_gap(inner_intent, inner_scale),
        _scaled_inside_auto_inner_gap(outer_intent, outer_scale),
    )


def _outer_limit_below_rows(
    limit_px: float,
    rows_above: Sequence[tuple[_RadialSlotIntent, float]],
    intent: _RadialSlotIntent,
) -> float:
    """Outer limit of ``intent`` below placed rows, each with the edge facing it."""

    return min(
        [float(limit_px)]
        + [float(edge_px) - _gap_between_inner_outer_tracks(intent, row) for row, edge_px in rows_above]
    )


def _inner_limit_above_rows(
    limit_px: float,
    rows_below: Sequence[tuple[_RadialSlotIntent, float]],
    intent: _RadialSlotIntent,
) -> float:
    """Inner limit of ``intent`` above placed rows, each with the edge facing it."""

    return max(
        [float(limit_px)]
        + [float(edge_px) + _gap_between_inner_outer_tracks(row, intent) for row, edge_px in rows_below]
    )


class _InsideFitError(ValidationError):
    """An inside slot or group that cannot be placed in its window."""

    def __init__(
        self,
        subject: str,
        *,
        intent: _RadialSlotIntent | None,
        placement_window: PlacementWindow,
        hint: str = "",
    ) -> None:
        self.subject = subject
        self.placement_window = placement_window
        self.hint = hint
        super().__init__(
            f"{self._between()}. Move the slot, reduce widths, disable conflicting labels, "
            f"or use side=outside.{hint}",
            diagnostic=_cannot_fit_diagnostic(intent, placement_window),
        )

    def _between(self) -> str:
        return (
            f"{self.subject} cannot fit inside between "
            f"{self.placement_window.inner_px:.1f}px and {self.placement_window.outer_px:.1f}px"
        )

    def limited_by_center(self, reserved_radius_px: float, *, explicit_radius: bool) -> ValidationError:
        """Name the center reservation (definition text or explicit radius) as the cause."""

        reserved = f"{float(reserved_radius_px):.1f}px"
        if explicit_radius:
            message = (
                f"{self._between()} because center_reserved_radius reserves {reserved}. "
                f"Set a smaller center_reserved_radius or place tracks outside.{self.hint}"
            )
        else:
            message = (
                f"{self._between()} because the center definition text reserves {reserved}. "
                "Shorten the species or strain text, reduce the definition font size, "
                f"set a smaller center_reserved_radius, or place tracks outside.{self.hint}"
            )
        return ValidationError(
            message,
            diagnostic=_center_reserved_diagnostic(self.diagnostic or {}, explicit_radius=explicit_radius),
        )


def _inside_stack_failure_hint(intents: Sequence[_RadialSlotIntent]) -> str:
    if any(intent.renderer == "sequence_conservation" for intent in intents):
        return (
            " For many similarity rings, reduce Ring Width/Ring Gap, disable GC/skew/depth tracks, "
            "move some tracks outside, or reduce the number of comparison files."
        )
    return ""


def _place_inside_auto_stack_group(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> tuple[CircularResolvedSlot, ...]:
    failure: _RadialSlotIntent | None = None
    for scale in _inside_auto_stack_width_scales(intents):
        working_occupied = list(occupied)
        working_outer = float(placement_window.outer_px)
        resolved_group: list[CircularResolvedSlot] = []
        failed = False
        for intent_index, intent in enumerate(intents):
            width_px, compressed = _scaled_inside_auto_width(intent, scale)
            place = _measure_anchored_inside if _anchored_inside(intent) else _place_inside_auto_fixed_width
            resolved = place(
                intent,
                width_px=width_px,
                compressed=compressed,
                occupied=working_occupied,
                axis_radius_px=axis_radius_px,
                placement_window=replace(placement_window, outer_px=working_outer),
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=total_length,
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
            )
            if resolved is None:
                failure = intent
                failed = True
                break
            resolved_group.append(resolved)
            if _slot_reserves(intent) and resolved.reserved_band_px is not None:
                working_occupied.append((intent.slot_id, resolved.reserved_band_px))
            next_inner_intent = intents[intent_index + 1] if intent_index + 1 < len(intents) else None
            if resolved.packing_band_px is not None and next_inner_intent is not None:
                gap_to_next = _gap_between_inner_outer_tracks(
                    next_inner_intent,
                    intent,
                    inner_scale=scale,
                    outer_scale=scale,
                )
                working_outer = min(
                    working_outer,
                    float(resolved.packing_band_px.inner_px) - gap_to_next,
                )
        if not failed:
            return tuple(resolved_group)

    # Name the slot that failed at the smallest scale; the band is the group's.
    failed_intent = failure
    failed_slot_id = failed_intent.slot_id if failed_intent is not None else "<empty>"
    raise _InsideFitError(
        f"Circular track slot '{failed_slot_id}'",
        intent=failed_intent,
        placement_window=placement_window,
        hint=_inside_stack_failure_hint(intents),
    )


def _place_outside_auto_stack_group(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    placement_window: PlacementWindow,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> tuple[CircularResolvedSlot, ...]:
    """Place ``intents``, given from the axis outward, from the window's inner edge."""

    working_occupied = list(occupied)
    working_inner = float(placement_window.inner_px)
    resolved_group: list[CircularResolvedSlot] = []
    for intent_index, intent in enumerate(intents):
        resolved = _place_outside_auto(
            intent,
            occupied=working_occupied,
            axis_radius_px=axis_radius_px,
            placement_window=PlacementWindow(working_inner, float(placement_window.outer_px)),
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
        )
        resolved_group.append(resolved)
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            working_occupied.append((intent.slot_id, resolved.reserved_band_px))
        next_outer_intent = intents[intent_index + 1] if intent_index + 1 < len(intents) else None
        if resolved.packing_band_px is not None and next_outer_intent is not None:
            gap_to_next = _gap_between_inner_outer_tracks(intent, next_outer_intent)
            working_inner = max(working_inner, float(resolved.packing_band_px.outer_px) + gap_to_next)

    return tuple(resolved_group)


def _linear_scales_from_1_to_min(min_scale: float, *, steps: int = 8) -> list[float]:
    lower = max(0.0, min(1.0, float(min_scale)))
    if lower >= 1.0 - LAYOUT_EPSILON:
        return [1.0]
    return [1.0 - ((1.0 - lower) * (idx / float(steps))) for idx in range(0, steps + 1)]


def _preferred_group_width_scales(intents: Sequence[_RadialSlotIntent]) -> list[float]:
    min_scales = [
        _min_readable_preferred_width_px(intent.renderer, intent.width_px) / intent.width_px
        for intent in intents
        if intent.width_px > LAYOUT_EPSILON
    ]
    if not min_scales:
        return [1.0]
    return _linear_scales_from_1_to_min(min(min_scales), steps=8)


def _scaled_preferred_width(intent: _RadialSlotIntent, scale: float) -> float:
    width = max(0.0, float(intent.width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    return max(_min_readable_preferred_width_px(intent.renderer, width), width * float(scale))


def _preferred_numeric_group_from(
    ordered_intents: Sequence[_RadialSlotIntent],
    start_pos: int,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> list[_RadialSlotIntent]:
    group: list[_RadialSlotIntent] = []
    for future in ordered_intents[start_pos:]:
        if future.slot_index in resolved_by_index:
            break
        if (
            future.placement_policy != "preferred"
            or future.side != "inside"
            or future.renderer not in NUMERIC_CIRCULAR_TRACK_RENDERERS
        ):
            break
        group.append(future)
    return group


def _group_packing_span(resolved_group: Sequence[CircularResolvedSlot]) -> RadialBand | None:
    return band_union(
        [slot.packing_band_px for slot in resolved_group if slot.packing_band_px is not None]
    )


def _group_fits_window_and_order(
    intents: Sequence[_RadialSlotIntent],
    resolved_group: Sequence[CircularResolvedSlot],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    placement_window: PlacementWindow,
) -> bool:
    if len(intents) != len(resolved_group):
        return False
    span = _group_packing_span(resolved_group)
    if span is not None:
        if span.inner_px < placement_window.inner_px - LAYOUT_EPSILON:
            return False
        if span.outer_px > placement_window.outer_px + LAYOUT_EPSILON:
            return False

    for previous_intent, current_intent, previous, current in zip(
        intents,
        intents[1:],
        resolved_group,
        resolved_group[1:],
    ):
        if previous.packing_band_px is None or current.packing_band_px is None:
            continue
        gap = _gap_between_inner_outer_tracks(current_intent, previous_intent)
        if current.packing_band_px.outer_px > previous.packing_band_px.inner_px - gap + LAYOUT_EPSILON:
            return False

    working_occupied = list(occupied)
    for intent, resolved in zip(intents, resolved_group):
        if not _slot_reserves(intent) or resolved.reserved_band_px is None:
            continue
        if _reserved_overlap_any(resolved.reserved_band_px, working_occupied) is not None:
            return False
        working_occupied.append((intent.slot_id, resolved.reserved_band_px))
    return True


def _try_place_preferred_numeric_group_at_anchors(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    placement_window: PlacementWindow,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> tuple[CircularResolvedSlot, ...] | None:
    resolved_group = tuple(
        _measure_radial_slot(
            intent,
            anchor_offset_px=float(intent.anchor_offset_px or 0.0),
            width_px=float(intent.width_px),
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
            compressed=False,
        )
        for intent in intents
    )
    if _group_fits_window_and_order(
        intents,
        resolved_group,
        occupied=occupied,
        placement_window=placement_window,
    ):
        return resolved_group
    return None


def _place_inside_auto_group_with_width_scale(
    intents: Sequence[_RadialSlotIntent],
    *,
    width_scale: float,
    occupied: Sequence[tuple[str, RadialBand]],
    placement_window: PlacementWindow,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> tuple[CircularResolvedSlot, ...] | None:
    working_occupied = list(occupied)
    working_outer = float(placement_window.outer_px)
    resolved_group: list[CircularResolvedSlot] = []
    for intent_index, intent in enumerate(intents):
        width_px = _scaled_preferred_width(intent, width_scale)
        resolved = _place_inside_auto_fixed_width(
            intent,
            width_px=width_px,
            compressed=width_px < float(intent.width_px) - LAYOUT_EPSILON,
            occupied=working_occupied,
            axis_radius_px=axis_radius_px,
            placement_window=replace(placement_window, outer_px=working_outer),
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
        )
        if resolved is None:
            return None
        resolved_group.append(resolved)
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            working_occupied.append((intent.slot_id, resolved.reserved_band_px))
        next_inner_intent = intents[intent_index + 1] if intent_index + 1 < len(intents) else None
        if resolved.packing_band_px is not None and next_inner_intent is not None:
            gap_to_next = _gap_between_inner_outer_tracks(next_inner_intent, intent)
            working_outer = min(working_outer, float(resolved.packing_band_px.inner_px) - gap_to_next)

    if _group_fits_window_and_order(
        intents,
        resolved_group,
        occupied=occupied,
        placement_window=placement_window,
    ):
        return tuple(resolved_group)
    return None


def _place_preferred_numeric_group(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    placement_window: PlacementWindow,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    depth_config: DepthConfigurator | None,
) -> tuple[CircularResolvedSlot, ...]:
    anchored = _try_place_preferred_numeric_group_at_anchors(
        intents,
        occupied=occupied,
        placement_window=placement_window,
        axis_radius_px=axis_radius_px,
        feature_dict=feature_dict,
        canvas_config=canvas_config,
        cfg=cfg,
        total_length=total_length,
        tick_track_channel_override=tick_track_channel_override,
        depth_config=depth_config,
    )
    if anchored is not None:
        return anchored

    for scale in _preferred_group_width_scales(intents):
        placed = _place_inside_auto_group_with_width_scale(
            intents,
            width_scale=scale,
            occupied=occupied,
            placement_window=placement_window,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
        )
        if placed is not None:
            return placed

    group_name = ",".join(intent.slot_id for intent in intents) or "<empty>"
    raise _InsideFitError(
        f"Preferred numeric group '{group_name}'",
        intent=intents[0] if intents else None,
        placement_window=placement_window,
    )


def _validate_same_side_order(
    slots: Sequence[CircularResolvedSlot],
    intent_by_index: Mapping[int, _RadialSlotIntent],
    movable_by_index: Mapping[int, bool],
) -> None:
    for side in ("outside", "inside"):
        side_slots = [
            slot for slot in sorted(slots, key=lambda item: item.slot_index)
            if slot.side == side and slot.renderer != "features" and slot.packing_band_px is not None
        ]
        for previous, current in zip(side_slots, side_slots[1:]):
            if not bool(movable_by_index.get(previous.slot_index)) and not bool(movable_by_index.get(current.slot_index)):
                continue
            previous_intent = intent_by_index.get(previous.slot_index)
            current_intent = intent_by_index.get(current.slot_index)
            if previous_intent is None or current_intent is None:
                continue
            spacing = (
                _gap_between_inner_outer_tracks(
                    current_intent,
                    previous_intent,
                    inner_scale=0.0,
                    outer_scale=0.0,
                )
                if side == "inside"
                else _gap_between_inner_outer_tracks(current_intent, previous_intent)
            )
            # side_slots above keeps only slots with a packing band.
            assert previous.packing_band_px is not None
            assert current.packing_band_px is not None
            # On both sides a later row lies inside the previous row.
            if current.packing_band_px.outer_px > previous.packing_band_px.inner_px - spacing + LAYOUT_EPSILON:
                raise ValidationError(
                    "Circular track slot order cannot be honored with the supplied pinned geometry: "
                    f"'{current.id}' would overlap or move outside '{previous.id}'.",
                    diagnostic=_cannot_fit_diagnostic(current_intent),
                )


def _next_future_hard_slot(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    side: str,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
    below_px: float | None = None,
) -> tuple[_RadialSlotIntent, CircularResolvedSlot] | None:
    """The next pinned row on ``side``; with ``below_px``, the next one below that edge.

    A row pinned wholly above ``below_px`` contradicts the stack order and
    cannot bound the rows before it from below; the order check judges it.
    """

    for future in ordered_intents[start_pos + 1:]:
        if future.side != side or future.placement_policy != "hard":
            continue
        resolved = resolved_by_index.get(future.slot_index)
        if resolved is None or resolved.packing_band_px is None:
            continue
        if below_px is not None and float(resolved.packing_band_px.inner_px) >= float(below_px) - LAYOUT_EPSILON:
            continue
        return future, resolved
    return None


def _minimum_future_inside_width_px(
    intent: _RadialSlotIntent,
    *,
    axis_radius_px: float,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthConfigurator | None,
) -> float:
    width = max(0.0, float(intent.width_px))
    if intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        if intent.placement_policy == "preferred":
            width = _min_readable_preferred_width_px(intent.renderer, width)
        elif intent.compress:
            width = _min_readable_numeric_width_px(intent.renderer, width)

    if intent.renderer in {"features", "ticks", "depth"}:
        resolved = _measure_radial_slot(
            intent,
            anchor_offset_px=0.0,
            width_px=width,
            axis_radius_px=float(axis_radius_px),
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=int(total_length),
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
            compressed=width < float(intent.width_px) - LAYOUT_EPSILON,
        )
        footprint_widths = [width]
        if resolved.packing_band_px is not None:
            footprint_widths.append(float(resolved.packing_band_px.width_px))
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            footprint_widths.append(float(resolved.reserved_band_px.width_px))
        return max(footprint_widths)

    return width


def _future_unresolved_inside_span_px(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    current_outer_intent: _RadialSlotIntent,
    axis_radius_px: float,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthConfigurator | None,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> float:
    span = 0.0
    outer_neighbor = current_outer_intent
    for future in ordered_intents[start_pos + 1:]:
        if future.slot_index in resolved_by_index:
            break
        if future.side != "inside" or future.placement_policy == "hard":
            break
        span += _gap_between_inner_outer_tracks(future, outer_neighbor) + _minimum_future_inside_width_px(
            future,
            axis_radius_px=axis_radius_px,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            feature_dict=feature_dict,
            depth_config=depth_config,
        )
        outer_neighbor = future
    return span


def _inner_limit_with_reserved_future_span(
    occupied: Sequence[tuple[str, RadialBand]],
    *,
    inner_limit_px: float,
    outer_limit_px: float,
    required_span_px: float,
) -> float:
    required = max(0.0, float(required_span_px))
    inner_limit = max(0.0, float(inner_limit_px))
    if required <= LAYOUT_EPSILON:
        return inner_limit

    intervals = _free_intervals(
        occupied,
        inner_limit_px=inner_limit,
        outer_limit_px=float(outer_limit_px),
    )
    for lower, upper in sorted(intervals, key=lambda item: item[1], reverse=True):
        if float(upper) - float(lower) >= required - LAYOUT_EPSILON:
            return max(inner_limit, float(lower) + required)
    return max(inner_limit, float(outer_limit_px))


def _inside_placement_window(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    current_outer_intent: _RadialSlotIntent,
    inside_max_outer: float,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthConfigurator | None,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
    group: Sequence[_RadialSlotIntent] = (),
) -> PlacementWindow:
    inner_limit = 0.0
    reserved_inner_limit: float | None = None
    # A pinned row listed after the group but lying above its window, or above
    # a pinned row of the group, contradicts the stack order: no lower bound.
    future_hard = _next_future_hard_slot(
        ordered_intents,
        start_pos=start_pos,
        side="inside",
        resolved_by_index=resolved_by_index,
        below_px=min(
            [float(inside_max_outer)]
            + [float(axis_radius_px) + float(member.anchor_offset_px or 0.0) for member in group if _anchored_inside(member)]
        ),
    )
    if future_hard is not None:
        future_hard_intent, future_hard_slot = future_hard
        # _next_future_hard_slot only returns slots that have a packing band.
        assert future_hard_slot.packing_band_px is not None
        gap_to_future_hard = _gap_between_inner_outer_tracks(future_hard_intent, current_outer_intent)
        # Mirror of top-down packing: the gap separates the packing bands, while
        # a reserved band may touch the pinned row.
        reserved_inner_limit = max(
            float(band.outer_px)
            for band in (future_hard_slot.packing_band_px, future_hard_slot.reserved_band_px)
            if band is not None
        )
        inner_limit = reserved_inner_limit + gap_to_future_hard
    future_span = _future_unresolved_inside_span_px(
        ordered_intents,
        start_pos=start_pos,
        current_outer_intent=current_outer_intent,
        axis_radius_px=axis_radius_px,
        canvas_config=canvas_config,
        cfg=cfg,
        total_length=total_length,
        tick_track_channel_override=tick_track_channel_override,
        feature_dict=feature_dict,
        depth_config=depth_config,
        resolved_by_index=resolved_by_index,
    )
    if future_span > LAYOUT_EPSILON:
        # The span kept free for the rows after this group bounds both bands.
        reserved_inner_limit = None
        inner_limit = max(
            inner_limit,
            _inner_limit_with_reserved_future_span(
                occupied,
                inner_limit_px=inner_limit,
                outer_limit_px=float(inside_max_outer),
                required_span_px=future_span,
            ),
        )
    return PlacementWindow(inner_limit, float(inside_max_outer), reserved_inner_limit)


def _outside_placement_window(
    outside_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    current_inner_intent: _RadialSlotIntent,
    outside_min_inner: float,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> PlacementWindow:
    """Window of an outside group; ``outside_intents`` run from the axis outward."""

    outer_limit = float("inf")
    future_hard = _next_future_hard_slot(
        outside_intents,
        start_pos=start_pos,
        side="outside",
        resolved_by_index=resolved_by_index,
    )
    if future_hard is not None:
        future_hard_intent, future_hard_slot = future_hard
        # _next_future_hard_slot only returns slots that have a packing band.
        assert future_hard_slot.packing_band_px is not None
        outer_limit = float(future_hard_slot.packing_band_px.inner_px) - _gap_between_inner_outer_tracks(
            current_inner_intent, future_hard_intent
        )
    return PlacementWindow(float(outside_min_inner), outer_limit)


class _RadialLayoutInputs(TypedDict):
    """Keyword arguments shared by both _resolve_circular_radial_layout attempts."""

    total_length: int
    canvas_config: CircularCanvasConfigurator
    slots: Sequence[CircularTrackSlot]
    feature_dict: Mapping[str, Any] | None
    tick_track_channel_override: str | None
    preferred_anchor_slot_ids: Collection[str]
    depth_config: DepthConfigurator | None
    center_reserved_radius_explicit: bool


def resolve_circular_radial_layout(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    slots: Sequence[CircularTrackSlot],
    feature_dict: Mapping[str, Any] | None = None,
    definition_reserved_radius_px: float | None = None,
    tick_track_channel_override: str | None = None,
    preferred_anchor_slot_ids: Collection[str] = (),
    depth_config: DepthConfigurator | None = None,
    center_reserved_radius_explicit: bool = False,
) -> CircularRadialLayout:
    """Resolve every slot band; a failed inside placement names its cause.

    When an inside slot cannot be placed but the same slots fit without the
    center reservation, the ``TRACK_LAYOUT`` diagnostic reason is
    ``DEFINITION_RESERVED`` (the definition text band) or ``CENTER_RESERVED``
    (``center_reserved_radius_explicit``); otherwise it is ``CANNOT_FIT``.
    """

    layout_inputs = _RadialLayoutInputs(
        total_length=total_length,
        canvas_config=canvas_config,
        slots=slots,
        feature_dict=feature_dict,
        tick_track_channel_override=tick_track_channel_override,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
        depth_config=depth_config,
        center_reserved_radius_explicit=center_reserved_radius_explicit,
    )
    try:
        return _resolve_circular_radial_layout(
            definition_reserved_radius_px=definition_reserved_radius_px,
            **layout_inputs,
        )
    except _InsideFitError as error:
        if definition_reserved_radius_px is None or definition_reserved_radius_px <= LAYOUT_EPSILON:
            raise
        try:
            _resolve_circular_radial_layout(definition_reserved_radius_px=None, **layout_inputs)
        except ValidationError:
            raise error from None
        raise error.limited_by_center(
            float(definition_reserved_radius_px),
            explicit_radius=center_reserved_radius_explicit,
        ) from None


def _resolve_circular_radial_layout(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    slots: Sequence[CircularTrackSlot],
    feature_dict: Mapping[str, Any] | None = None,
    definition_reserved_radius_px: float | None = None,
    tick_track_channel_override: str | None = None,
    preferred_anchor_slot_ids: Collection[str] = (),
    depth_config: DepthConfigurator | None = None,
    center_reserved_radius_explicit: bool = False,
) -> CircularRadialLayout:
    cfg = canvas_config.profile.config
    axis_radius_px = float(canvas_config.radius)
    axis = CircularAxisLayout(
        radius_px=axis_radius_px,
        stroke_width_px=float(cfg.objects.axis.circular.stroke_width.for_length_param(str(canvas_config.length_param))),
    )
    definition_band = (
        RadialBand(0.0, float(definition_reserved_radius_px))
        if definition_reserved_radius_px is not None and definition_reserved_radius_px > LAYOUT_EPSILON
        else None
    )
    occupied: list[tuple[str, RadialBand]] = []
    if definition_band is not None and definition_band.width_px > LAYOUT_EPSILON:
        occupied.append(("definition", definition_band))

    intents = _slot_intents(
        slots,
        canvas_config=canvas_config,
        cfg=cfg,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
    )
    resolved_by_index: dict[int, CircularResolvedSlot] = {}
    intent_by_index = {intent.slot_index: intent for intent in intents}
    movable_by_index = {
        intent.slot_index: intent.placement_policy in {"auto", "preferred"}
        for intent in intents
    }

    # Manual anchors and overlays become blockers before movable placement. A
    # pinned row with an Auto width is placed with the inside rows around it; at
    # full width it bounds only the outside rows.
    anchored_blockers: list[tuple[str, RadialBand]] = []
    for intent in intents:
        if intent.slot_index in resolved_by_index:
            continue
        if intent.placement_policy not in {"hard", "overlay"}:
            continue
        anchor_offset = float(intent.anchor_offset_px or 0.0)
        resolved = _measure_radial_slot(
            intent,
            anchor_offset_px=anchor_offset,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=int(total_length),
            tick_track_channel_override=tick_track_channel_override,
            depth_config=depth_config,
        )
        if _anchored_inside(intent):
            if resolved.reserved_band_px is not None:
                anchored_blockers.append((intent.slot_id, resolved.reserved_band_px))
            continue
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            conflict = _reserved_overlap_any(resolved.reserved_band_px, occupied)
            if conflict is not None:
                diagnostic = _cannot_fit_diagnostic(intent)
                if conflict[0] == "definition":
                    diagnostic = _center_reserved_diagnostic(diagnostic, explicit_radius=center_reserved_radius_explicit)
                raise ValidationError(
                    f"Pinned circular track slot '{intent.slot_id}' overlaps reserved circular slot '{conflict[0]}'.",
                    diagnostic=diagnostic,
                )
        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            band = (
                resolved.reserved_band_px.expanded_sides(intent.inner_gap_px, intent.outer_gap_px)
                if intent.side == "overlay"
                else resolved.reserved_band_px
            )
            occupied.append((intent.slot_id, band))

    default_axis_gap_px = _default_spacing_px(axis_radius_px)
    outside_axis_gap_px = max(
        [float(intent.inner_gap_px) for intent in intents if intent.side == "outside"] or [default_axis_gap_px]
    )
    inside_axis_gap_px = max(
        [float(intent.outer_gap_px) for intent in intents if intent.side == "inside"] or [default_axis_gap_px]
    )
    outside_min_inner = axis_radius_px + outside_axis_gap_px
    inside_max_outer = axis_radius_px - inside_axis_gap_px
    # Placed rows that bound the next row of a side, each with its edge that
    # faces that row; the next row keeps the facing gap from every one of them.
    rows_below_outside: list[tuple[_RadialSlotIntent, float]] = []
    rows_above_inside: list[tuple[_RadialSlotIntent, float]] = []
    # Stack order runs from the outermost row to the innermost. A pinned row
    # bounds the rows after it on its own side (below); a pinned feature row
    # also bounds every row on the other side of the axis.
    for slot_index, resolved in resolved_by_index.items():
        slot_intent = intent_by_index.get(slot_index)
        if slot_intent is None or resolved.packing_band_px is None or resolved.renderer != "features":
            continue
        if resolved.side != "outside":
            rows_below_outside.append((slot_intent, max(axis_radius_px, float(resolved.packing_band_px.outer_px))))
        if resolved.side != "inside":
            rows_above_inside.append((slot_intent, min(axis_radius_px, float(resolved.packing_band_px.inner_px))))

    ordered_intents = sorted(intents, key=lambda item: item.slot_index)

    # Outside rows are placed from the axis outward: each row's window starts at
    # the row below it and ends at the next pinned row above it.
    outside_intents = [intent for intent in reversed(ordered_intents) if intent.side == "outside"]
    for outside_pos, intent in enumerate(outside_intents):
        if intent.slot_index not in resolved_by_index:
            outside_group = _outside_auto_stack_group_from(outside_intents, outside_pos, resolved_by_index)
            resolved_group = _place_outside_auto_stack_group(
                outside_group,
                occupied=[*occupied, *anchored_blockers],
                axis_radius_px=axis_radius_px,
                placement_window=_outside_placement_window(
                    outside_intents,
                    start_pos=outside_pos + len(outside_group) - 1,
                    current_inner_intent=outside_group[-1],
                    outside_min_inner=_inner_limit_above_rows(outside_min_inner, rows_below_outside, outside_group[0]),
                    resolved_by_index=resolved_by_index,
                ),
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
            )
            for group_intent, group_resolved in zip(outside_group, resolved_group):
                resolved_by_index[group_intent.slot_index] = group_resolved
                if _slot_reserves(group_intent) and group_resolved.reserved_band_px is not None:
                    occupied.append((group_intent.slot_id, group_resolved.reserved_band_px))
        outside_band = resolved_by_index[intent.slot_index].packing_band_px
        if outside_band is not None:
            outer_px = float(outside_band.outer_px)
            if intent.renderer == "features":
                outer_px = max(axis_radius_px, outer_px)
            rows_below_outside.append((intent, outer_px))

    for intent_pos, intent in enumerate(ordered_intents):
        already_resolved = resolved_by_index.get(intent.slot_index)
        if already_resolved is not None:
            if (
                intent.placement_policy == "hard"
                and intent.side == "inside"
                and already_resolved.packing_band_px is not None
            ):
                inner_px = float(already_resolved.packing_band_px.inner_px)
                if intent.renderer == "features":
                    inner_px = min(axis_radius_px, inner_px)
                rows_above_inside.append((intent, inner_px))
            continue

        if intent.side == "overlay":
            resolved = _measure_radial_slot(
                intent,
                anchor_offset_px=0.0,
                axis_radius_px=axis_radius_px,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
            )
        else:
            # Every group starts with this row, which faces the rows above it.
            outer_limit = _outer_limit_below_rows(inside_max_outer, rows_above_inside, intent)
            movable_group = _inside_movable_stack_group_from(
                ordered_intents,
                intent_pos,
                resolved_by_index,
            )
            if len(movable_group) > 1 or _anchored_inside(intent):
                placement_window = _inside_placement_window(
                    ordered_intents,
                    start_pos=intent_pos + len(movable_group) - 1,
                    current_outer_intent=movable_group[-1],
                    inside_max_outer=outer_limit,
                    occupied=occupied,
                    axis_radius_px=axis_radius_px,
                    canvas_config=canvas_config,
                    cfg=cfg,
                    total_length=int(total_length),
                    tick_track_channel_override=tick_track_channel_override,
                    feature_dict=feature_dict,
                    depth_config=depth_config,
                    resolved_by_index=resolved_by_index,
                    group=movable_group,
                )
                resolved_group = _place_inside_auto_stack_group(
                    movable_group,
                    occupied=occupied,
                    axis_radius_px=axis_radius_px,
                    placement_window=placement_window,
                    feature_dict=feature_dict,
                    canvas_config=canvas_config,
                    cfg=cfg,
                    total_length=int(total_length),
                    tick_track_channel_override=tick_track_channel_override,
                    depth_config=depth_config,
                )
                for group_intent, group_resolved in zip(movable_group, resolved_group):
                    resolved_by_index[group_intent.slot_index] = group_resolved
                    if _slot_reserves(group_intent) and group_resolved.reserved_band_px is not None:
                        occupied.append((group_intent.slot_id, group_resolved.reserved_band_px))
                    if group_resolved.packing_band_px is not None:
                        rows_above_inside.append((group_intent, float(group_resolved.packing_band_px.inner_px)))
                continue

            preferred_group = _preferred_numeric_group_from(
                ordered_intents,
                intent_pos,
                resolved_by_index,
            )
            if preferred_group:
                placement_window = _inside_placement_window(
                    ordered_intents,
                    start_pos=intent_pos + len(preferred_group) - 1,
                    current_outer_intent=preferred_group[-1],
                    inside_max_outer=outer_limit,
                    occupied=occupied,
                    axis_radius_px=axis_radius_px,
                    canvas_config=canvas_config,
                    cfg=cfg,
                    total_length=int(total_length),
                    tick_track_channel_override=tick_track_channel_override,
                    feature_dict=feature_dict,
                    depth_config=depth_config,
                    resolved_by_index=resolved_by_index,
                )
                resolved_group = _place_preferred_numeric_group(
                    preferred_group,
                    occupied=occupied,
                    placement_window=placement_window,
                    axis_radius_px=axis_radius_px,
                    feature_dict=feature_dict,
                    canvas_config=canvas_config,
                    cfg=cfg,
                    total_length=int(total_length),
                    tick_track_channel_override=tick_track_channel_override,
                    depth_config=depth_config,
                )
                for group_intent, group_resolved in zip(preferred_group, resolved_group):
                    resolved_by_index[group_intent.slot_index] = group_resolved
                    if _slot_reserves(group_intent) and group_resolved.reserved_band_px is not None:
                        occupied.append((group_intent.slot_id, group_resolved.reserved_band_px))
                    if group_resolved.packing_band_px is not None:
                        rows_above_inside.append((group_intent, float(group_resolved.packing_band_px.inner_px)))
                continue

            inside_group = _inside_auto_stack_group_from(
                ordered_intents,
                intent_pos,
                resolved_by_index,
            )
            if not inside_group:
                inside_group = [intent]
            placement_window = _inside_placement_window(
                ordered_intents,
                start_pos=intent_pos + len(inside_group) - 1,
                current_outer_intent=inside_group[-1],
                inside_max_outer=outer_limit,
                occupied=occupied,
                axis_radius_px=axis_radius_px,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
                feature_dict=feature_dict,
                depth_config=depth_config,
                resolved_by_index=resolved_by_index,
            )
            if len(inside_group) > 1:
                try:
                    resolved_group = _place_inside_auto_stack_group(
                        inside_group,
                        occupied=occupied,
                        axis_radius_px=axis_radius_px,
                        placement_window=placement_window,
                        feature_dict=feature_dict,
                        canvas_config=canvas_config,
                        cfg=cfg,
                        total_length=int(total_length),
                        tick_track_channel_override=tick_track_channel_override,
                        depth_config=depth_config,
                    )
                except ValidationError:
                    fallback_group = _inside_movable_stack_group_from(
                        ordered_intents,
                        intent_pos,
                        resolved_by_index,
                    )
                    if len(fallback_group) <= len(inside_group):
                        raise
                    fallback_window = _inside_placement_window(
                        ordered_intents,
                        start_pos=intent_pos + len(fallback_group) - 1,
                        current_outer_intent=fallback_group[-1],
                        inside_max_outer=outer_limit,
                        occupied=occupied,
                        axis_radius_px=axis_radius_px,
                        canvas_config=canvas_config,
                        cfg=cfg,
                        total_length=int(total_length),
                        tick_track_channel_override=tick_track_channel_override,
                        feature_dict=feature_dict,
                        depth_config=depth_config,
                        resolved_by_index=resolved_by_index,
                    )
                    resolved_group = _place_inside_auto_stack_group(
                        fallback_group,
                        occupied=occupied,
                        axis_radius_px=axis_radius_px,
                        placement_window=fallback_window,
                        feature_dict=feature_dict,
                        canvas_config=canvas_config,
                        cfg=cfg,
                        total_length=int(total_length),
                        tick_track_channel_override=tick_track_channel_override,
                        depth_config=depth_config,
                    )
                    inside_group = fallback_group
                for group_intent, group_resolved in zip(inside_group, resolved_group):
                    resolved_by_index[group_intent.slot_index] = group_resolved
                    if _slot_reserves(group_intent) and group_resolved.reserved_band_px is not None:
                        occupied.append((group_intent.slot_id, group_resolved.reserved_band_px))
                    if group_resolved.packing_band_px is not None:
                        rows_above_inside.append((group_intent, float(group_resolved.packing_band_px.inner_px)))
                continue
            resolved = _place_inside_auto(
                intent,
                occupied=occupied,
                axis_radius_px=axis_radius_px,
                placement_window=placement_window,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
            )

        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            band = (
                resolved.reserved_band_px.expanded_sides(intent.inner_gap_px, intent.outer_gap_px)
                if intent.side == "overlay"
                else resolved.reserved_band_px
            )
            occupied.append((intent.slot_id, band))
        if resolved.packing_band_px is not None and resolved.side == "inside":
            rows_above_inside.append((intent, float(resolved.packing_band_px.inner_px)))

    by_id = {slot.id: slot for slot in resolved_by_index.values()}
    for intent in intents:
        if intent.renderer != "annotations" or intent.side != "overlay":
            continue
        overlay_resolved = resolved_by_index.get(intent.slot_index)
        anchor_id = str(intent.params.get("anchor_slot", "")).strip()
        anchor = by_id.get(anchor_id)
        if overlay_resolved is None or anchor is None:
            continue
        anchor_band = anchor.draw_band_px or anchor.packing_band_px
        anchor_radius = (
            float(anchor_band.center_px)
            if anchor_band is not None
            else float(anchor.anchor_radius_px or axis_radius_px)
        )
        cover_anchor = str(intent.params.get("cover_anchor", "false")).strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        width = (
            float(anchor_band.width_px)
            if cover_anchor and anchor_band is not None
            else max(0.0, float(overlay_resolved.resolved_width_px or intent.width_px))
        )
        draw_band = _band_from_center_width(anchor_radius, width)
        anchored = replace(
            overlay_resolved,
            anchor_radius_px=anchor_radius,
            anchor_offset_px=anchor_radius - axis_radius_px,
            resolved_width_px=width,
            packing_band_px=None,
            draw_band_px=draw_band,
            reserved_band_px=None,
        )
        resolved_by_index[intent.slot_index] = anchored
        by_id[intent.slot_id] = anchored

    resolved_slots = tuple(resolved_by_index[index] for index in sorted(resolved_by_index))
    _validate_same_side_order(resolved_slots, intent_by_index, movable_by_index)
    outer_content_radius = max(
        [axis_radius_px]
        + ([definition_band.outer_px] if definition_band is not None else [])
        + [
            float(slot.reserved_band_px.outer_px)
            for slot in resolved_slots
            if slot.reserved_band_px is not None
        ]
    )
    return CircularRadialLayout(
        axis=axis,
        slots=resolved_slots,
        definition_reserved_band_px=definition_band,
        outer_content_radius_px=float(outer_content_radius),
    )


__all__ = [
    "build_circular_feature_layout",
    "resolve_circular_radial_layout",
]
