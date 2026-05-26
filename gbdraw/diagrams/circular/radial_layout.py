"""Unified radial slot resolver for circular diagrams."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Collection, Literal, Mapping, Sequence, TypeAlias

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import DepthConfigurator  # type: ignore[reportMissingImports]
from ...exceptions import ValidationError  # type: ignore[reportMissingImports]
from ...layout.circular_depth_axis import (  # type: ignore[reportMissingImports]
    DepthAxisFootprint,
    resolve_depth_axis_footprint,
)
from ...tracks.circular import (  # type: ignore[reportMissingImports]
    NUMERIC_CIRCULAR_TRACK_RENDERERS,
    CircularTrackSlot,
    NormalizedCircularTrackSlot,
    normalize_circular_track_slots,
    tick_sides_for_tick_label_layout,
)
from ...svg.circular_ticks import (  # type: ignore[reportMissingImports]
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_radius_bounds,
)
from .presets import normalize_circular_track_preset  # type: ignore[reportMissingImports]


logger = logging.getLogger(__name__)

LAYOUT_EPSILON = 1e-6
TRACK_SPACING_MULTIPLIER = 1.2
MIN_NUMERIC_WIDTH_PX = 10.0
MIN_NUMERIC_WIDTH_FRACTION = 0.55
MIN_SKEW_WIDTH_PX = 12.0
MIN_SKEW_WIDTH_FRACTION = 0.65
PREFERRED_MIN_NUMERIC_WIDTH_FRACTION = 0.4
PREFERRED_MIN_SKEW_WIDTH_FRACTION = 0.4

PlacementPolicy = Literal["hard", "preferred", "auto", "overlay"]


@dataclass(frozen=True)
class RadialBand:
    inner_px: float
    outer_px: float

    def __post_init__(self) -> None:
        inner, outer = sorted((max(0.0, float(self.inner_px)), max(0.0, float(self.outer_px))))
        object.__setattr__(self, "inner_px", inner)
        object.__setattr__(self, "outer_px", outer)

    @property
    def width_px(self) -> float:
        return max(0.0, float(self.outer_px) - float(self.inner_px))

    @property
    def center_px(self) -> float:
        return (float(self.inner_px) + float(self.outer_px)) / 2.0

    def expanded(self, spacing_px: float) -> "RadialBand":
        spacing = max(0.0, float(spacing_px))
        return RadialBand(max(0.0, self.inner_px - spacing), self.outer_px + spacing)


@dataclass(frozen=True)
class CircularFeatureLane:
    track_id: int
    strand_group: str
    inner_px: float
    center_px: float
    outer_px: float

    @property
    def band_px(self) -> RadialBand:
        return RadialBand(self.inner_px, self.outer_px)


@dataclass(frozen=True)
class CircularFeatureStackMetrics:
    lane_count: int
    lane_width_px: float
    lane_spacing_px: float
    band_width_px: float
    center_radius_px: float
    inner_radius_px: float
    outer_radius_px: float
    lane_centers_by_track_id: Mapping[int, float]


@dataclass(frozen=True)
class CircularFeatureLayout:
    anchor_radius_px: float
    width_px: float
    lanes_by_track_id: Mapping[int, CircularFeatureLane]
    primary_band_px: RadialBand
    all_band_px: RadialBand
    stack_metrics: CircularFeatureStackMetrics | None = None

    def lane_for_track_id(self, track_id: int) -> CircularFeatureLane:
        if int(track_id) in self.lanes_by_track_id:
            return self.lanes_by_track_id[int(track_id)]
        if 0 in self.lanes_by_track_id:
            return self.lanes_by_track_id[0]
        return next(iter(self.lanes_by_track_id.values()))


@dataclass(frozen=True)
class CircularAxisLayout:
    radius_px: float
    stroke_width_px: float


@dataclass(frozen=True)
class CircularTickLayout:
    anchor_radius_px: float
    tick_band_px: RadialBand
    label_band_px: RadialBand | None
    reserved_band_px: RadialBand
    track_preset: str
    label_side: str
    tick_side: str
    tick_length_px: float | None


@dataclass(frozen=True)
class CircularFeatureLabelLayout:
    side: str
    reserved_band_px: RadialBand
    labels: tuple[Mapping[str, Any], ...] = ()


CircularSlotPayload: TypeAlias = CircularFeatureLayout | CircularTickLayout | CircularFeatureLabelLayout


@dataclass(frozen=True)
class CircularResolvedSlot:
    slot_index: int
    id: str
    renderer: str
    side: str
    z: int
    anchor_radius_px: float | None
    anchor_offset_px: float | None
    requested_width_px: float | None
    resolved_width_px: float | None
    packing_band_px: RadialBand | None
    draw_band_px: RadialBand | None
    reserved_band_px: RadialBand | None
    params: Mapping[str, Any]
    payload: CircularSlotPayload | None = None
    explicit_anchor: bool = False
    explicit_width: bool = False
    compressed: bool = False

    @property
    def center_radius_px(self) -> float:
        if self.draw_band_px is not None:
            return float(self.draw_band_px.center_px)
        return float(self.anchor_radius_px or 0.0)

    @property
    def draw_inner_radius_px(self) -> float:
        return float(self.draw_band_px.inner_px) if self.draw_band_px is not None else 0.0

    @property
    def draw_outer_radius_px(self) -> float:
        return float(self.draw_band_px.outer_px) if self.draw_band_px is not None else 0.0

    @property
    def reserved_inner_radius_px(self) -> float:
        return float(self.reserved_band_px.inner_px) if self.reserved_band_px is not None else 0.0

    @property
    def reserved_outer_radius_px(self) -> float:
        return float(self.reserved_band_px.outer_px) if self.reserved_band_px is not None else 0.0

    @property
    def draw_width_px(self) -> float:
        return float(self.draw_band_px.width_px) if self.draw_band_px is not None else 0.0

    @property
    def reserved_width_px(self) -> float:
        return float(self.reserved_band_px.width_px) if self.reserved_band_px is not None else 0.0


@dataclass(frozen=True)
class CircularRadialLayout:
    axis: CircularAxisLayout
    slots: tuple[CircularResolvedSlot, ...]
    definition_reserved_band_px: RadialBand | None
    outer_content_radius_px: float

    @property
    def features(self) -> CircularFeatureLayout | None:
        for slot in self.slots:
            if slot.renderer == "features" and isinstance(slot.payload, CircularFeatureLayout):
                return slot.payload
        return None

    @property
    def ticks(self) -> CircularTickLayout | None:
        for slot in self.slots:
            if slot.renderer == "ticks" and isinstance(slot.payload, CircularTickLayout):
                return slot.payload
        return None

    @property
    def tracks(self) -> tuple[CircularResolvedSlot, ...]:
        return tuple(slot for slot in self.slots if slot.renderer not in {"features", "ticks", "feature_labels"})


@dataclass(frozen=True)
class PlacementWindow:
    inner_px: float
    outer_px: float


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
    spacing_px: float
    z: int
    compress: bool
    reserve: bool
    placement_policy: PlacementPolicy
    params: Mapping[str, Any]


def band_union(bands: Sequence[RadialBand]) -> RadialBand | None:
    non_empty = [band for band in bands if band.width_px > LAYOUT_EPSILON]
    if not non_empty:
        return None
    return RadialBand(
        min(float(band.inner_px) for band in non_empty),
        max(float(band.outer_px) for band in non_empty),
    )


def bands_overlap(a: RadialBand, b: RadialBand) -> bool:
    return (
        float(a.inner_px) < (float(b.outer_px) - LAYOUT_EPSILON)
        and float(a.outer_px) > (float(b.inner_px) + LAYOUT_EPSILON)
    )


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


def _feature_track_ids(feature_dict: Mapping[str, Any] | None) -> tuple[int, ...]:
    if not feature_dict:
        return (0,)
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
        positive = sorted((track_id for track_id in ids if track_id >= 0), key=lambda value: (value != 0, value))
        negative = sorted((track_id for track_id in ids if track_id < 0), key=lambda value: abs(value))
        return tuple(positive + negative)
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
    if normalized_preset == "middle":
        step = width + spacing
        ids_set = set(ids)
        has_negative = any(track_id < 0 for track_id in ids)
        if bool(strandedness) or has_negative:
            lane_centers = {}
            split_seed = (0.5 * width) + (0.5 * spacing)
            for idx, track_id in enumerate(sorted((tid for tid in ids if tid < 0), key=lambda value: abs(value))):
                lane_centers[int(track_id)] = max(0.0, center - split_seed - (idx * step))
            for idx, track_id in enumerate(sorted((tid for tid in ids if tid >= 0), key=lambda value: (value != 0, value))):
                lane_centers[int(track_id)] = max(0.0, center + split_seed + (idx * step))
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
) -> CircularFeatureLayout | None:
    width = max(0.0, float(width_px))
    direction = (
        str(lane_direction).strip().lower()
        if lane_direction is not None
        else _lane_direction_from_legacy_track_type(track_type)
    )
    track_ids = _feature_track_ids(feature_dict)
    metrics = measure_circular_feature_stack(
        axis_radius_px=float(axis_radius_px),
        lane_width_px=width,
        lane_spacing_px=lane_spacing_px,
        preset=track_type,
        lane_direction=direction,
        strandedness=bool(strandedness),
        track_ids=track_ids,
        center_radius_px=anchor_radius_px,
    )
    anchor = float(metrics.center_radius_px)

    lanes: dict[int, CircularFeatureLane] = {}
    for track_id in sorted(track_ids):
        center = float(metrics.lane_centers_by_track_id[int(track_id)])
        strand_group = _feature_lane_strand_group(int(track_id), bool(strandedness))
        half_width = width / 2.0
        lanes[int(track_id)] = CircularFeatureLane(
            track_id=int(track_id),
            strand_group=strand_group,
            inner_px=max(0.0, float(center) - half_width),
            center_px=max(0.0, float(center)),
            outer_px=max(0.0, float(center) + half_width),
        )

    all_band = band_union([lane.band_px for lane in lanes.values()]) or RadialBand(anchor, anchor)
    primary_lanes = [lane.band_px for tid, lane in lanes.items() if tid in {0, -1}]
    primary_band = band_union(primary_lanes) or all_band
    return CircularFeatureLayout(
        anchor_radius_px=anchor,
        width_px=width,
        lanes_by_track_id=lanes,
        primary_band_px=primary_band,
        all_band_px=all_band,
        stack_metrics=metrics,
    )


def feature_radii_for_object(
    feature_object: Any,
    feature_layout: CircularFeatureLayout,
) -> tuple[float, float, float]:
    lane = feature_layout.lane_for_track_id(int(getattr(feature_object, "feature_track_id", 0)))
    return float(lane.inner_px), float(lane.center_px), float(lane.outer_px)


def feature_radii_for_coordinate(
    feature_object: Any,
    coordinate: Any,
    feature_layout: CircularFeatureLayout,
) -> tuple[float, float, float]:
    del coordinate
    return feature_radii_for_object(feature_object, feature_layout)


def _tick_layout_from_params(
    *,
    axis_radius_px: float,
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
            strandedness=bool(cfg.canvas.strandedness),
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
        strandedness=bool(cfg.canvas.strandedness),
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


def _preferred_anchor_radius_px(params: Mapping[str, Any], axis_radius_px: float) -> float | None:
    raw = params.get("_preferred_anchor_radius")
    if raw is None:
        return None
    if hasattr(raw, "resolve"):
        try:
            return float(raw.resolve(float(axis_radius_px)))
        except (TypeError, ValueError):
            return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


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
            _preferred_anchor_radius_px(slot.params, axis_radius_px)
            if (
                radius_px is None
                and slot.id in preferred_ids
                and slot.side == "inside"
                and slot.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
            )
            else None
        )
        width_px = slot.width.resolve(axis_radius_px) if slot.width is not None else None
        if width_px is None:
            width_px = _default_width_px(slot.renderer, canvas_config=canvas_config, cfg=cfg)
        spacing_px = (
            float(slot.spacing.resolve(axis_radius_px))
            if slot.spacing is not None
            else _default_spacing_px(axis_radius_px)
        )
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
        intents.append(
            _RadialSlotIntent(
                slot=slot,
                slot_index=slot.slot_index,
                slot_id=slot.id,
                renderer=slot.renderer,
                side=slot.side,
                anchor_offset_px=(
                    float(radius_px if radius_px is not None else preferred_anchor_px) - axis_radius_px
                    if (radius_px is not None or preferred_anchor_px is not None)
                    else None
                ),
                width_px=max(0.0, float(width_px)),
                explicit_anchor=explicit_anchor,
                explicit_width=slot.width is not None,
                spacing_px=max(0.0, float(spacing_px)),
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
            strandedness=bool(cfg.canvas.strandedness),
            anchor_radius_px=anchor_radius_px,
            lane_spacing_px=float(intent.spacing_px),
        )
        if feature_layout is not None:
            anchor_radius_px = float(feature_layout.anchor_radius_px)
            anchor_offset_px = float(anchor_radius_px) - float(axis_radius_px)
        band = feature_layout.all_band_px if feature_layout is not None else RadialBand(anchor_radius_px, anchor_radius_px)
        draw_band = feature_layout.primary_band_px if feature_layout is not None else band
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
            axis_radius_px=float(axis_radius_px),
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
    if renderer == "dinucleotide_skew":
        return min(width, max(MIN_SKEW_WIDTH_PX, MIN_SKEW_WIDTH_FRACTION * width))
    return min(width, max(MIN_NUMERIC_WIDTH_PX, MIN_NUMERIC_WIDTH_FRACTION * width))


def _min_readable_preferred_width_px(renderer: str, default_width_px: float) -> float:
    width = max(0.0, float(default_width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if renderer == "dinucleotide_skew":
        return min(width, max(MIN_SKEW_WIDTH_PX, PREFERRED_MIN_SKEW_WIDTH_FRACTION * width))
    return min(width, max(MIN_NUMERIC_WIDTH_PX, PREFERRED_MIN_NUMERIC_WIDTH_FRACTION * width))


def _min_dense_stack_numeric_width_px(renderer: str, default_width_px: float) -> float:
    width = max(0.0, float(default_width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
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
    del intent
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
                    anchor_offset += float(band.outer_px) - float(resolved.reserved_band_px.inner_px) + intent.spacing_px
                    continue
            return resolved
    raise ValidationError(f"Circular track slot '{intent.slot_id}' cannot be placed outside without overlap.")


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
    raise ValidationError(
        f"Circular track slot '{intent.slot_id}' cannot fit inside between "
        f"{placement_window.inner_px:.1f}px and {placement_window.outer_px:.1f}px. "
        "Move the slot, reduce widths, disable conflicting labels, or use side=outside."
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
    intervals = _free_intervals(
        occupied,
        inner_limit_px=float(placement_window.inner_px),
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
        )
        if resolved is None or resolved.reserved_band_px is None:
            continue
        if _reserved_overlap_any(resolved.reserved_band_px, occupied) is None:
            return resolved
    return None


def _shrinkable_inside_numeric(intent: _RadialSlotIntent) -> bool:
    return (
        intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
        and intent.side == "inside"
        and intent.compress
        and not intent.explicit_anchor
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
        if future.side != "inside" or future.explicit_anchor:
            break
        if future.placement_policy == "auto":
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
    for scale in _inside_auto_stack_width_scales(intents):
        working_occupied = list(occupied)
        working_outer = float(placement_window.outer_px)
        resolved_group: list[CircularResolvedSlot] = []
        failed = False
        for intent in intents:
            width_px, compressed = _scaled_inside_auto_width(intent, scale)
            resolved = _place_inside_auto_fixed_width(
                intent,
                width_px=width_px,
                compressed=compressed,
                occupied=working_occupied,
                axis_radius_px=axis_radius_px,
                placement_window=PlacementWindow(float(placement_window.inner_px), working_outer),
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=total_length,
                tick_track_channel_override=tick_track_channel_override,
                depth_config=depth_config,
            )
            if resolved is None:
                failed = True
                break
            resolved_group.append(resolved)
            if _slot_reserves(intent) and resolved.reserved_band_px is not None:
                working_occupied.append((intent.slot_id, resolved.reserved_band_px))
            if resolved.packing_band_px is not None:
                working_outer = min(working_outer, float(resolved.packing_band_px.inner_px) - intent.spacing_px)
        if not failed:
            return tuple(resolved_group)

    first_unplaced = intents[-1].slot_id if intents else "<empty>"
    raise ValidationError(
        f"Circular track slot '{first_unplaced}' cannot fit inside between "
        f"{placement_window.inner_px:.1f}px and {placement_window.outer_px:.1f}px. "
        "Move the slot, reduce widths, disable conflicting labels, or use side=outside."
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
    working_occupied = list(occupied)
    working_inner = float(placement_window.inner_px)
    resolved_by_slot_index: dict[int, CircularResolvedSlot] = {}

    # Outside rows are displayed from outermost to innermost. The packer fills
    # from the axis outward, so place that visual stack in reverse.
    for intent in reversed(tuple(intents)):
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
        resolved_by_slot_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            working_occupied.append((intent.slot_id, resolved.reserved_band_px))
        if resolved.packing_band_px is not None:
            working_inner = max(working_inner, float(resolved.packing_band_px.outer_px) + intent.spacing_px)

    return tuple(resolved_by_slot_index[intent.slot_index] for intent in intents)


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

    for previous_intent, previous, current in zip(intents, resolved_group, resolved_group[1:]):
        if previous.packing_band_px is None or current.packing_band_px is None:
            continue
        if current.packing_band_px.outer_px > previous.packing_band_px.inner_px - previous_intent.spacing_px + LAYOUT_EPSILON:
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
    for intent in intents:
        width_px = _scaled_preferred_width(intent, width_scale)
        resolved = _place_inside_auto_fixed_width(
            intent,
            width_px=width_px,
            compressed=width_px < float(intent.width_px) - LAYOUT_EPSILON,
            occupied=working_occupied,
            axis_radius_px=axis_radius_px,
            placement_window=PlacementWindow(float(placement_window.inner_px), working_outer),
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
        if resolved.packing_band_px is not None:
            working_outer = min(working_outer, float(resolved.packing_band_px.inner_px) - intent.spacing_px)

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
    raise ValidationError(
        f"Preferred numeric group '{group_name}' cannot fit inside between "
        f"{placement_window.inner_px:.1f}px and {placement_window.outer_px:.1f}px."
    )


def _validate_same_side_order(
    slots: Sequence[CircularResolvedSlot],
    spacing_by_index: Mapping[int, float],
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
            spacing = max(0.0, float(spacing_by_index.get(previous.slot_index, 0.0)))
            if side == "outside":
                if previous.packing_band_px.inner_px < current.packing_band_px.outer_px + spacing - LAYOUT_EPSILON:
                    raise ValidationError(
                        "Circular track slot order cannot be honored with the supplied pinned geometry: "
                        f"'{current.id}' would overlap or move outside '{previous.id}'."
                    )
            else:
                if current.packing_band_px.outer_px > previous.packing_band_px.inner_px - spacing + LAYOUT_EPSILON:
                    raise ValidationError(
                        "Circular track slot order cannot be honored with the supplied pinned geometry: "
                        f"'{current.id}' would overlap or move outside '{previous.id}'."
                    )


def _next_future_hard_slot(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    side: str,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> CircularResolvedSlot | None:
    for future in ordered_intents[start_pos + 1:]:
        if future.side != side or future.placement_policy != "hard":
            continue
        resolved = resolved_by_index.get(future.slot_index)
        if resolved is not None and resolved.packing_band_px is not None:
            return resolved
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
        if resolved.reserved_band_px is not None:
            footprint_widths.append(float(resolved.reserved_band_px.width_px))
        return max(footprint_widths)

    return width


def _future_unresolved_inside_span_px(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    current_spacing_px: float,
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
    spacing_before = max(0.0, float(current_spacing_px))
    for future in ordered_intents[start_pos + 1:]:
        if future.slot_index in resolved_by_index:
            break
        if future.side != "inside" or future.placement_policy == "hard":
            break
        span += spacing_before + _minimum_future_inside_width_px(
            future,
            axis_radius_px=axis_radius_px,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            feature_dict=feature_dict,
            depth_config=depth_config,
        )
        spacing_before = max(0.0, float(future.spacing_px))
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
    current_spacing_px: float,
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
) -> PlacementWindow:
    inner_limit = 0.0
    future_hard = _next_future_hard_slot(
        ordered_intents,
        start_pos=start_pos,
        side="inside",
        resolved_by_index=resolved_by_index,
    )
    if future_hard is not None and future_hard.packing_band_px is not None:
        inner_limit = max(
            inner_limit,
            float(future_hard.packing_band_px.outer_px) + max(0.0, float(current_spacing_px)),
        )
    future_span = _future_unresolved_inside_span_px(
        ordered_intents,
        start_pos=start_pos,
        current_spacing_px=current_spacing_px,
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
        inner_limit = max(
            inner_limit,
            _inner_limit_with_reserved_future_span(
                occupied,
                inner_limit_px=inner_limit,
                outer_limit_px=float(inside_max_outer),
                required_span_px=future_span,
            ),
        )
    return PlacementWindow(inner_limit, float(inside_max_outer))


def _outside_placement_window(
    ordered_intents: Sequence[_RadialSlotIntent],
    *,
    start_pos: int,
    current_spacing_px: float,
    outside_min_inner: float,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> PlacementWindow:
    outer_limit = float("inf")
    future_hard = _next_future_hard_slot(
        ordered_intents,
        start_pos=start_pos,
        side="outside",
        resolved_by_index=resolved_by_index,
    )
    if future_hard is not None and future_hard.packing_band_px is not None:
        outer_limit = float(future_hard.packing_band_px.inner_px) - max(0.0, float(current_spacing_px))
    return PlacementWindow(float(outside_min_inner), outer_limit)


def resolve_circular_radial_layout(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    slots: Sequence[CircularTrackSlot],
    feature_dict: Mapping[str, Any] | None = None,
    show_features: bool = True,
    show_ticks: bool = True,
    definition_reserved_radius_px: float | None = None,
    feature_track_ratio_factor_override: float | None = None,
    tick_track_channel_override: str | None = None,
    preferred_anchor_slot_ids: Collection[str] = (),
    depth_config: DepthConfigurator | None = None,
) -> CircularRadialLayout:
    del show_features, show_ticks, feature_track_ratio_factor_override
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
    spacing_by_index = {intent.slot_index: intent.spacing_px for intent in intents}
    movable_by_index = {
        intent.slot_index: intent.placement_policy in {"auto", "preferred"}
        for intent in intents
    }

    # Manual anchors and overlays become blockers before movable placement.
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
        if resolved.reserved_band_px is not None:
            conflict = _reserved_overlap_any(resolved.reserved_band_px, occupied)
            if conflict is not None:
                message = f"Pinned circular track slot '{intent.slot_id}' overlaps reserved circular slot '{conflict[0]}'."
                raise ValidationError(message)
        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            band = resolved.reserved_band_px.expanded(intent.spacing_px) if intent.side == "overlay" else resolved.reserved_band_px
            occupied.append((intent.slot_id, band))

    outside_min_inner = axis_radius_px + _default_spacing_px(axis_radius_px)
    inside_max_outer = axis_radius_px - _default_spacing_px(axis_radius_px)
    for slot_index, resolved in resolved_by_index.items():
        intent = intent_by_index.get(slot_index)
        if intent is None or resolved.packing_band_px is None:
            continue
        if resolved.renderer == "features":
            outside_min_inner = max(
                outside_min_inner,
                max(axis_radius_px, float(resolved.packing_band_px.outer_px)) + intent.spacing_px,
            )
            inside_max_outer = min(
                inside_max_outer,
                min(axis_radius_px, float(resolved.packing_band_px.inner_px)) - intent.spacing_px,
            )
        elif resolved.side == "outside":
            outside_min_inner = max(outside_min_inner, float(resolved.packing_band_px.outer_px) + intent.spacing_px)
        elif resolved.side == "inside":
            inside_max_outer = min(inside_max_outer, float(resolved.packing_band_px.inner_px) - intent.spacing_px)

    ordered_intents = sorted(intents, key=lambda item: item.slot_index)

    for intent_pos, intent in enumerate(ordered_intents):
        resolved = resolved_by_index.get(intent.slot_index)
        if resolved is not None:
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
        elif intent.side == "outside":
            outside_group = _outside_auto_stack_group_from(
                ordered_intents,
                intent_pos,
                resolved_by_index,
            )
            if len(outside_group) > 1:
                placement_window = _outside_placement_window(
                    ordered_intents,
                    start_pos=intent_pos + len(outside_group) - 1,
                    current_spacing_px=outside_group[-1].spacing_px,
                    outside_min_inner=outside_min_inner,
                    resolved_by_index=resolved_by_index,
                )
                resolved_group = _place_outside_auto_stack_group(
                    outside_group,
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
                for group_intent, group_resolved in zip(outside_group, resolved_group):
                    resolved_by_index[group_intent.slot_index] = group_resolved
                    if _slot_reserves(group_intent) and group_resolved.reserved_band_px is not None:
                        occupied.append((group_intent.slot_id, group_resolved.reserved_band_px))
                    if group_resolved.packing_band_px is not None:
                        outside_min_inner = max(
                            outside_min_inner,
                            float(group_resolved.packing_band_px.outer_px) + group_intent.spacing_px,
                        )
                continue
            placement_window = _outside_placement_window(
                ordered_intents,
                start_pos=intent_pos,
                current_spacing_px=intent.spacing_px,
                outside_min_inner=outside_min_inner,
                resolved_by_index=resolved_by_index,
            )
            resolved = _place_outside_auto(
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
        else:
            movable_group = _inside_movable_stack_group_from(
                ordered_intents,
                intent_pos,
                resolved_by_index,
            )
            if len(movable_group) > 1:
                placement_window = _inside_placement_window(
                    ordered_intents,
                    start_pos=intent_pos + len(movable_group) - 1,
                    current_spacing_px=movable_group[-1].spacing_px,
                    inside_max_outer=inside_max_outer,
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
                        inside_max_outer = min(
                            inside_max_outer,
                            float(group_resolved.packing_band_px.inner_px) - group_intent.spacing_px,
                        )
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
                    current_spacing_px=preferred_group[-1].spacing_px,
                    inside_max_outer=inside_max_outer,
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
                        inside_max_outer = min(
                            inside_max_outer,
                            float(group_resolved.packing_band_px.inner_px) - group_intent.spacing_px,
                        )
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
                current_spacing_px=inside_group[-1].spacing_px,
                inside_max_outer=inside_max_outer,
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
                        current_spacing_px=fallback_group[-1].spacing_px,
                        inside_max_outer=inside_max_outer,
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
                        inside_max_outer = min(
                            inside_max_outer,
                            float(group_resolved.packing_band_px.inner_px) - group_intent.spacing_px,
                        )
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
            band = resolved.reserved_band_px.expanded(intent.spacing_px) if intent.side == "overlay" else resolved.reserved_band_px
            occupied.append((intent.slot_id, band))
        if resolved.packing_band_px is not None:
            if resolved.side == "outside":
                outside_min_inner = max(outside_min_inner, float(resolved.packing_band_px.outer_px) + intent.spacing_px)
            elif resolved.side == "inside":
                inside_max_outer = min(inside_max_outer, float(resolved.packing_band_px.inner_px) - intent.spacing_px)

    resolved_slots = tuple(resolved_by_index[index] for index in sorted(resolved_by_index))
    _validate_same_side_order(resolved_slots, spacing_by_index, movable_by_index)
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


def feature_band_bounds_px(feature_layout: CircularFeatureLayout | None) -> tuple[float, float] | None:
    if feature_layout is None:
        return None
    return float(feature_layout.all_band_px.inner_px), float(feature_layout.all_band_px.outer_px)


def feature_radius_intervals(
    feature_dict: Mapping[str, Any],
    total_length: int,
    feature_layout: CircularFeatureLayout,
) -> list[tuple[float, float, float, float]]:
    if total_length <= 0:
        return []

    intervals: list[tuple[float, float, float, float]] = []
    total_length_float = float(total_length)
    for feature_object in feature_dict.values():
        lane = feature_layout.lane_for_track_id(int(getattr(feature_object, "feature_track_id", 0)))
        feature_location_list = getattr(feature_object, "location", [])
        list_of_coordinates = getattr(feature_object, "coordinates", [])

        for coordinate_idx, coordinate in enumerate(list_of_coordinates):
            if (
                coordinate_idx < len(feature_location_list)
                and getattr(feature_location_list[coordinate_idx], "kind", None) == "line"
            ):
                continue

            coordinate_start_raw = int(coordinate.start)
            coordinate_end_raw = int(coordinate.end)
            if coordinate_start_raw == coordinate_end_raw:
                continue

            coordinate_start = float(coordinate_start_raw % total_length)
            coordinate_end = float(coordinate_end_raw % total_length)
            if coordinate_end_raw > coordinate_start_raw and coordinate_end_raw <= total_length:
                coordinate_end = float(coordinate_end_raw)

            if coordinate_start < coordinate_end:
                intervals.append((coordinate_start, coordinate_end, lane.inner_px, lane.outer_px))
            else:
                intervals.append((coordinate_start, total_length_float, lane.inner_px, lane.outer_px))
                intervals.append((0.0, coordinate_end, lane.inner_px, lane.outer_px))
    return intervals


__all__ = [
    "RadialBand",
    "CircularFeatureLane",
    "CircularFeatureLayout",
    "CircularAxisLayout",
    "CircularTickLayout",
    "CircularFeatureLabelLayout",
    "CircularResolvedSlot",
    "CircularSlotPayload",
    "CircularRadialLayout",
    "band_union",
    "bands_overlap",
    "build_circular_feature_layout",
    "feature_band_bounds_px",
    "feature_radii_for_coordinate",
    "feature_radii_for_object",
    "feature_radius_intervals",
    "resolve_circular_radial_layout",
]
