"""Unified radial slot resolver for circular diagrams."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Mapping, Sequence, TypeAlias

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...exceptions import ValidationError  # type: ignore[reportMissingImports]
from ...tracks.circular import (  # type: ignore[reportMissingImports]
    NUMERIC_CIRCULAR_TRACK_RENDERERS,
    CircularTrackSlot,
    NormalizedCircularTrackSlot,
    normalize_circular_track_slots,
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
class CircularFeatureLayout:
    anchor_radius_px: float
    width_px: float
    lanes_by_track_id: Mapping[int, CircularFeatureLane]
    primary_band_px: RadialBand
    all_band_px: RadialBand

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
    strict: bool
    compress: bool
    reserve: bool
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


def _lane_direction_from_legacy_track_type(track_type: str | None) -> str:
    preset = normalize_circular_track_preset(track_type)
    if preset == "middle":
        return "split"
    if preset == "spreadout":
        return "outside"
    return "inside"


def _feature_lane_center(
    *,
    axis_radius_px: float,
    width_px: float,
    lane_direction: str,
    strandedness: bool,
    track_id: int,
) -> tuple[float, str]:
    width = max(0.0, float(width_px))
    step = width * TRACK_SPACING_MULTIPLIER
    axis = float(axis_radius_px)
    layout = str(lane_direction or "inside").strip().lower()
    track = int(track_id)

    if strandedness:
        if layout == "outside":
            if track < 0:
                return axis + (0.9 * width) + ((abs(track) - 1) * step), "negative"
            return axis + (1.9 * width) + (max(0, track) * step), "positive"
        if layout == "inside":
            if track < 0:
                return axis - (2.0 * width) - ((abs(track) - 1) * step), "negative"
            return axis - (1.0 * width) - (max(0, track) * step), "positive"
        if track < 0:
            return axis - (0.5 * width) - ((abs(track) - 1) * step), "negative"
        return axis + (0.5 * width) + (max(0, track) * step), "positive"

    if layout == "outside":
        return axis + (0.9 * width) + (abs(track) * step), "combined"
    if layout == "inside":
        return axis - (0.75 * width) - (abs(track) * step), "combined"
    if track < 0:
        return axis - (abs(track) * step), "negative"
    if track > 0:
        return axis + (track * step), "positive"
    return axis, "combined"


def build_circular_feature_layout(
    feature_dict: Mapping[str, Any] | None,
    *,
    axis_radius_px: float,
    width_px: float,
    track_type: str | None = None,
    lane_direction: str | None = None,
    strandedness: bool,
    anchor_radius_px: float | None = None,
) -> CircularFeatureLayout | None:
    width = max(0.0, float(width_px))
    anchor = float(axis_radius_px if anchor_radius_px is None else anchor_radius_px)
    direction = (
        str(lane_direction).strip().lower()
        if lane_direction is not None
        else _lane_direction_from_legacy_track_type(track_type)
    )
    track_ids: set[int] = {0}
    if feature_dict:
        track_ids = {
            int(getattr(feature_object, "feature_track_id", 0))
            for feature_object in feature_dict.values()
        } or {0}

    lanes: dict[int, CircularFeatureLane] = {}
    for track_id in sorted(track_ids):
        center, strand_group = _feature_lane_center(
            axis_radius_px=anchor,
            width_px=width,
            lane_direction=direction,
            strandedness=strandedness,
            track_id=int(track_id),
        )
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
    label_side = str(params.get("label_side", "inside")).strip().lower()
    tick_side = str(params.get("tick_side", "inside")).strip().lower()
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
) -> list[_RadialSlotIntent]:
    axis_radius_px = float(canvas_config.radius)
    intents: list[_RadialSlotIntent] = []
    for slot in normalize_circular_track_slots(slots):
        radius_px = slot.radius.resolve(axis_radius_px) if slot.radius is not None else None
        width_px = slot.width.resolve(axis_radius_px) if slot.width is not None else None
        if width_px is None:
            width_px = _default_width_px(slot.renderer, canvas_config=canvas_config, cfg=cfg)
        spacing_px = (
            float(slot.spacing.resolve(axis_radius_px))
            if slot.spacing is not None
            else _default_spacing_px(axis_radius_px)
        )
        intents.append(
            _RadialSlotIntent(
                slot=slot,
                slot_index=slot.slot_index,
                slot_id=slot.id,
                renderer=slot.renderer,
                side=slot.side,
                anchor_offset_px=(float(radius_px) - axis_radius_px) if radius_px is not None else None,
                width_px=max(0.0, float(width_px)),
                explicit_anchor=radius_px is not None,
                explicit_width=slot.width is not None,
                spacing_px=max(0.0, float(spacing_px)),
                z=slot.z,
                strict=slot.strict,
                compress=slot.compress,
                reserve=slot.reserve,
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
    compressed: bool = False,
) -> CircularResolvedSlot:
    anchor_radius_px = float(axis_radius_px) + float(anchor_offset_px)
    resolved_width = max(0.0, float(intent.width_px if width_px is None else width_px))
    renderer = intent.renderer

    if renderer == "features":
        feature_layout = build_circular_feature_layout(
            feature_dict,
            axis_radius_px=float(axis_radius_px),
            width_px=resolved_width,
            lane_direction=str(intent.params.get("lane_direction", "inside")),
            strandedness=bool(cfg.canvas.strandedness),
            anchor_radius_px=anchor_radius_px,
        )
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
        tick_layout = _tick_layout_from_params(
            axis_radius_px=float(axis_radius_px),
            total_length=int(total_length),
            canvas_config=canvas_config,
            cfg=cfg,
            params=intent.params,
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
    return intent.side != "overlay" or bool(intent.reserve)


def _place_outside_auto(
    intent: _RadialSlotIntent,
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    min_packing_inner_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
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
                compressed=compressed,
            )
            if resolved.packing_band_px is not None and resolved.packing_band_px.inner_px < min_packing_inner_px - LAYOUT_EPSILON:
                anchor_offset += min_packing_inner_px - float(resolved.packing_band_px.inner_px)
                continue
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
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
    compressed: bool,
) -> CircularResolvedSlot | None:
    lower, upper = float(interval[0]), float(interval[1])
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
            max_outer = max(float(band.outer_px) for band in bands)
            if max_outer > upper + LAYOUT_EPSILON:
                anchor_offset -= max_outer - upper
                continue
            if min_inner < lower - LAYOUT_EPSILON:
                anchor_offset += lower - min_inner
                continue
            return resolved
    return None


def _place_inside_auto(
    intent: _RadialSlotIntent,
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    max_packing_outer_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
) -> CircularResolvedSlot:
    for width_px, compressed in _candidate_widths(intent):
        resolved = _place_inside_auto_fixed_width(
            intent,
            width_px=width_px,
            compressed=compressed,
            occupied=occupied,
            axis_radius_px=axis_radius_px,
            max_packing_outer_px=max_packing_outer_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
        )
        if resolved is not None:
            return resolved
    raise ValidationError(
        f"Circular track slot '{intent.slot_id}' cannot fit inside. "
        "Move the slot, reduce widths, disable conflicting labels, or use side=outside."
    )


def _place_inside_auto_fixed_width(
    intent: _RadialSlotIntent,
    *,
    width_px: float,
    compressed: bool,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    max_packing_outer_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
) -> CircularResolvedSlot | None:
    outer_limit = max(0.0, float(max_packing_outer_px))
    intervals = _free_intervals(occupied, inner_limit_px=0.0, outer_limit_px=outer_limit)
    for interval in sorted(intervals, key=lambda item: item[1], reverse=True):
        resolved = _try_measure_inside_interval(
            intent,
            width_px=width_px,
            interval=interval,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            canvas_config=canvas_config,
            cfg=cfg,
            total_length=total_length,
            tick_track_channel_override=tick_track_channel_override,
            compressed=compressed,
        )
        if resolved is None or resolved.reserved_band_px is None:
            continue
        if _reserved_overlap_any(resolved.reserved_band_px, occupied) is None:
            return resolved
    return None


def _place_inside_auto_group(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    axis_radius_px: float,
    max_packing_outer_px: float,
    feature_dict: Mapping[str, Any] | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    total_length: int,
    tick_track_channel_override: str | None,
) -> tuple[CircularResolvedSlot, ...]:
    candidate_lists = [_candidate_widths(intent) for intent in intents]
    max_candidates = max(len(candidates) for candidates in candidate_lists)
    for candidate_index in range(max_candidates):
        working_occupied = list(occupied)
        working_outer = float(max_packing_outer_px)
        resolved_group: list[CircularResolvedSlot] = []
        failed = False
        for intent, candidates in zip(intents, candidate_lists):
            width_px, compressed = candidates[min(candidate_index, len(candidates) - 1)]
            resolved = _place_inside_auto_fixed_width(
                intent,
                width_px=width_px,
                compressed=compressed,
                occupied=working_occupied,
                axis_radius_px=axis_radius_px,
                max_packing_outer_px=working_outer,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=total_length,
                tick_track_channel_override=tick_track_channel_override,
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
        f"Circular track slot '{first_unplaced}' cannot fit inside. "
        "Move the slot, reduce widths, disable conflicting labels, or use side=outside."
    )


def _validate_same_side_order(slots: Sequence[CircularResolvedSlot], spacing_by_index: Mapping[int, float]) -> None:
    def is_preset_generated(slot: CircularResolvedSlot) -> bool:
        return bool(slot.params.get("_preset_generated"))

    for side in ("outside", "inside"):
        side_slots = [
            slot for slot in sorted(slots, key=lambda item: item.slot_index)
            if slot.side == side and slot.packing_band_px is not None
        ]
        for previous, current in zip(side_slots, side_slots[1:]):
            if is_preset_generated(previous) and is_preset_generated(current):
                continue
            spacing = max(0.0, float(spacing_by_index.get(previous.slot_index, 0.0)))
            if side == "outside":
                if current.packing_band_px.inner_px < previous.packing_band_px.outer_px + spacing - LAYOUT_EPSILON:
                    raise ValidationError(
                        "Circular track slot order cannot be honored with the supplied pinned geometry: "
                        f"'{current.id}' would overlap or move inside '{previous.id}'."
                    )
            else:
                if current.packing_band_px.outer_px > previous.packing_band_px.inner_px - spacing + LAYOUT_EPSILON:
                    raise ValidationError(
                        "Circular track slot order cannot be honored with the supplied pinned geometry: "
                        f"'{current.id}' would overlap or move outside '{previous.id}'."
                    )


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

    intents = _slot_intents(slots, canvas_config=canvas_config, cfg=cfg)
    resolved_by_index: dict[int, CircularResolvedSlot] = {}
    spacing_by_index = {intent.slot_index: intent.spacing_px for intent in intents}

    # Fixed anchors and reserving overlays become blockers before auto placement.
    for intent in intents:
        if not intent.explicit_anchor and not (intent.side == "overlay" and intent.reserve):
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
        )
        if resolved.reserved_band_px is not None:
            conflict = _reserved_overlap_any(resolved.reserved_band_px, occupied)
            if conflict is not None:
                message = f"Pinned circular track slot '{intent.slot_id}' overlaps reserved circular slot '{conflict[0]}'."
                if intent.strict:
                    raise ValidationError(message)
                if not (bool(intent.params.get("_preset_generated")) and conflict[0] == "definition"):
                    logger.warning(message)
        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and resolved.reserved_band_px is not None:
            band = resolved.reserved_band_px.expanded(intent.spacing_px) if intent.side == "overlay" else resolved.reserved_band_px
            occupied.append((intent.slot_id, band))

    outside_min_inner = axis_radius_px
    inside_max_outer = axis_radius_px
    ordered_intents = sorted(intents, key=lambda item: item.slot_index)
    for intent_pos, intent in enumerate(ordered_intents):
        resolved = resolved_by_index.get(intent.slot_index)
        if resolved is not None and resolved.packing_band_px is not None:
            if resolved.side == "outside":
                outside_min_inner = max(outside_min_inner, float(resolved.packing_band_px.outer_px) + intent.spacing_px)
            elif resolved.side == "inside":
                inside_max_outer = min(inside_max_outer, float(resolved.packing_band_px.inner_px) - intent.spacing_px)
            continue
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
            )
        elif intent.side == "outside":
            resolved = _place_outside_auto(
                intent,
                occupied=occupied,
                axis_radius_px=axis_radius_px,
                min_packing_inner_px=outside_min_inner,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
            )
        else:
            inside_group = [intent]
            if intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS and intent.compress and not intent.explicit_anchor:
                for future in ordered_intents[intent_pos + 1:]:
                    if future.slot_index in resolved_by_index:
                        break
                    if (
                        future.side != "inside"
                        or future.explicit_anchor
                        or future.renderer not in NUMERIC_CIRCULAR_TRACK_RENDERERS
                        or not future.compress
                    ):
                        break
                    inside_group.append(future)
            if len(inside_group) > 1:
                resolved_group = _place_inside_auto_group(
                    inside_group,
                    occupied=occupied,
                    axis_radius_px=axis_radius_px,
                    max_packing_outer_px=inside_max_outer,
                    feature_dict=feature_dict,
                    canvas_config=canvas_config,
                    cfg=cfg,
                    total_length=int(total_length),
                    tick_track_channel_override=tick_track_channel_override,
                )
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
                max_packing_outer_px=inside_max_outer,
                feature_dict=feature_dict,
                canvas_config=canvas_config,
                cfg=cfg,
                total_length=int(total_length),
                tick_track_channel_override=tick_track_channel_override,
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
    _validate_same_side_order(resolved_slots, spacing_by_index)
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
