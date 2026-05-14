"""Resolved radial layout for circular diagrams.

This module is the single place where circular diagram radii are converted into
absolute pixel geometry. Renderers should consume the resolved values and avoid
re-deriving placement from track type or legacy track dictionaries.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...tracks.circular import CircularTrackSlot  # type: ignore[reportMissingImports]
from ...tracks.spec import CircularTrackPlacement, TrackSpec  # type: ignore[reportMissingImports]
from ...svg.circular_ticks import (  # type: ignore[reportMissingImports]
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_radius_bounds,
)


logger = logging.getLogger(__name__)

LAYOUT_EPSILON = 1e-6
TRACK_SPACING_MULTIPLIER = 1.2
NUMERIC_RENDERERS = frozenset({"dinucleotide_content", "dinucleotide_skew", "depth"})


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

    def expanded(self, gap_px: float) -> RadialBand:
        gap = max(0.0, float(gap_px))
        return RadialBand(max(0.0, self.inner_px - gap), self.outer_px + gap)


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
class CircularResolvedTrack:
    id: str
    renderer: str
    channel: str
    anchor_radius_px: float | None
    draw_band_px: RadialBand
    reserved_band_px: RadialBand
    z: int
    params: Mapping[str, Any]
    explicit_width: bool = False


@dataclass(frozen=True)
class CircularRadialLayout:
    axis: CircularAxisLayout
    features: CircularFeatureLayout | None
    ticks: CircularTickLayout | None
    tracks: tuple[CircularResolvedTrack, ...]
    definition_reserved_band_px: RadialBand | None
    outer_content_radius_px: float


@dataclass(frozen=True)
class _SlotIntent:
    slot: CircularTrackSlot
    center_px: float | None
    width_px: float
    explicit_anchor: bool
    explicit_annulus: bool
    explicit_width: bool
    channel: str
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


def _normalize_side(raw: object, default: str = "inside") -> str:
    side = str(raw or default).strip().lower()
    if side in {"inside", "outside", "overlay"}:
        return side
    logger.warning("Unknown circular track side '%s'; using '%s'.", raw, default)
    return default


def _parse_bool_param(value: object, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    return default


def _track_specs_by_match(track_specs: Sequence[TrackSpec] | None) -> tuple[dict[str, TrackSpec], dict[str, TrackSpec]]:
    by_id: dict[str, TrackSpec] = {}
    by_kind: dict[str, TrackSpec] = {}
    for ts in track_specs or []:
        by_id.setdefault(str(ts.id), ts)
        by_kind.setdefault(str(ts.kind), ts)
    return by_id, by_kind


def _legacy_kind_for_slot(slot: CircularTrackSlot) -> str | None:
    if slot.id == "gc_content":
        return "gc_content"
    if slot.id == "gc_skew":
        return "gc_skew"
    if slot.id in {"features", "ticks", "depth"}:
        return str(slot.id)
    renderer = str(slot.renderer)
    if renderer == "dinucleotide_content":
        return "gc_content"
    if renderer == "dinucleotide_skew":
        return "gc_skew"
    if renderer in {"features", "ticks", "depth"}:
        return renderer
    return None


def _track_spec_for_slot(
    slot: CircularTrackSlot,
    track_specs: Sequence[TrackSpec] | None,
) -> TrackSpec | None:
    by_id, by_kind = _track_specs_by_match(track_specs)
    if str(slot.id) in by_id:
        return by_id[str(slot.id)]
    legacy_kind = _legacy_kind_for_slot(slot)
    if legacy_kind is not None:
        return by_kind.get(legacy_kind)
    return None


def _resolve_placement_center_and_width(
    placement: CircularTrackPlacement | None,
    *,
    base_radius_px: float,
) -> tuple[float | None, float | None, bool, bool, bool]:
    if placement is None:
        return None, None, False, False, False

    inner_px = placement.inner_radius.resolve(base_radius_px) if placement.inner_radius is not None else None
    outer_px = placement.outer_radius.resolve(base_radius_px) if placement.outer_radius is not None else None
    radius_px = placement.radius.resolve(base_radius_px) if placement.radius is not None else None
    width_px = placement.width.resolve(base_radius_px) if placement.width is not None else None
    explicit_annulus = inner_px is not None and outer_px is not None

    if explicit_annulus:
        if outer_px < inner_px:
            inner_px, outer_px = outer_px, inner_px
        return (inner_px + outer_px) / 2.0, outer_px - inner_px, False, True, True
    return radius_px, width_px, radius_px is not None, False, width_px is not None


def _resolve_track_spec_center_and_width(
    ts: TrackSpec | None,
    *,
    base_radius_px: float,
) -> tuple[float | None, float | None, bool, bool, bool]:
    if ts is None or not isinstance(ts.placement, CircularTrackPlacement):
        return None, None, False, False, False
    return _resolve_placement_center_and_width(ts.placement, base_radius_px=base_radius_px)


def _feature_width_from_specs(
    *,
    feature_slot: CircularTrackSlot | None,
    feature_ts: TrackSpec | None,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    feature_track_ratio_factor_override: float | None = None,
) -> tuple[float, float]:
    base_radius = float(canvas_config.radius)
    base_track_ratio = float(canvas_config.track_ratio)
    length_param = str(canvas_config.length_param)
    default_ratio_factor = float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    ratio_factor = (
        float(feature_track_ratio_factor_override)
        if feature_track_ratio_factor_override is not None
        else default_ratio_factor
    )

    width_px: float | None = None
    if feature_ts is not None and isinstance(feature_ts.placement, CircularTrackPlacement):
        _center, spec_width, _explicit_anchor, _explicit_annulus, explicit_width = _resolve_placement_center_and_width(
            feature_ts.placement,
            base_radius_px=base_radius,
        )
        if explicit_width and spec_width is not None:
            width_px = float(spec_width)

    if feature_slot is not None:
        if feature_slot.width is not None:
            width_px = float(feature_slot.width.resolve(base_radius))
        elif feature_slot.placement is not None:
            _center, slot_width, _explicit_anchor, _explicit_annulus, explicit_width = _resolve_placement_center_and_width(
                feature_slot.placement,
                base_radius_px=base_radius,
            )
            if explicit_width and slot_width is not None:
                width_px = float(slot_width)

    if width_px is not None and width_px > 0.0:
        ratio_factor = width_px / max(LAYOUT_EPSILON, base_radius * base_track_ratio)
        return float(width_px), float(ratio_factor)

    width_px = base_radius * base_track_ratio * ratio_factor
    return float(width_px), float(ratio_factor)


def _feature_anchor_from_slot(
    feature_slot: CircularTrackSlot | None,
    *,
    base_radius_px: float,
) -> float:
    if feature_slot is None:
        return float(base_radius_px)
    center_px, _width_px, explicit_anchor, explicit_annulus, _explicit_width = _resolve_placement_center_and_width(
        feature_slot.placement,
        base_radius_px=base_radius_px,
    )
    if (explicit_anchor or explicit_annulus) and center_px is not None:
        return float(center_px)
    return float(base_radius_px)


def _feature_lane_center(
    *,
    axis_radius_px: float,
    width_px: float,
    track_type: str,
    strandedness: bool,
    track_id: int,
) -> tuple[float, str]:
    width = max(0.0, float(width_px))
    step = width * TRACK_SPACING_MULTIPLIER
    axis = float(axis_radius_px)
    layout = str(track_type or "middle").strip().lower()
    track = int(track_id)

    if strandedness:
        if layout == "spreadout":
            if track < 0:
                return axis + (0.9 * width) + ((abs(track) - 1) * step), "negative"
            return axis + (1.9 * width) + (max(0, track) * step), "positive"
        if layout == "tuckin":
            if track < 0:
                return axis - (2.0 * width) - ((abs(track) - 1) * step), "negative"
            return axis - (1.0 * width) - (max(0, track) * step), "positive"
        if track < 0:
            return axis - (0.5 * width) - ((abs(track) - 1) * step), "negative"
        return axis + (0.5 * width) + (max(0, track) * step), "positive"

    if layout == "spreadout":
        return axis + (0.9 * width) + (abs(track) * step), "combined"
    if layout == "tuckin":
        return axis - (1.2 * width) - (abs(track) * step), "combined"
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
    track_type: str,
    strandedness: bool,
    anchor_radius_px: float | None = None,
) -> CircularFeatureLayout | None:
    width = max(0.0, float(width_px))
    anchor = float(axis_radius_px if anchor_radius_px is None else anchor_radius_px)
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
            track_type=track_type,
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


def _tick_layout_from_slot(
    *,
    axis_radius_px: float,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    tick_slot: CircularTrackSlot | None,
    tick_ts: TrackSpec | None,
    tick_track_channel_override: str | None = None,
) -> CircularTickLayout | None:
    if tick_ts is not None and not tick_ts.show:
        return None

    base_radius = float(canvas_config.radius)
    params: dict[str, Any] = {}
    if tick_slot is not None:
        params.update(tick_slot.params)
    if tick_ts is not None and tick_ts.params:
        params.update(tick_ts.params)

    center_px, width_px, explicit_anchor, _explicit_annulus, explicit_width = _resolve_track_spec_center_and_width(
        tick_ts,
        base_radius_px=base_radius,
    )
    if tick_slot is not None:
        slot_center, slot_width, slot_explicit_anchor, _slot_annulus, slot_explicit_width = _resolve_placement_center_and_width(
            tick_slot.placement,
            base_radius_px=base_radius,
        )
        if slot_center is not None:
            center_px = slot_center
            explicit_anchor = slot_explicit_anchor
        if tick_slot.width is not None:
            width_px = float(tick_slot.width.resolve(base_radius))
            explicit_width = True
        elif slot_width is not None:
            width_px = slot_width
            explicit_width = slot_explicit_width

    anchor = float(center_px) if explicit_anchor and center_px is not None else float(axis_radius_px)
    label_side = str(params.get("label_side", "legacy")).strip().lower()
    tick_side = str(params.get("tick_side", "legacy")).strip().lower()
    tick_length_px = float(width_px) if explicit_width and width_px is not None and width_px > 0 else None

    if tick_side in {"none", ""}:
        tick_band = RadialBand(anchor, anchor)
    else:
        tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
            center_radius_px=anchor,
            total_len=total_length,
            size="large",
            track_type=str(cfg.canvas.circular.track_type),
            strandedness=bool(cfg.canvas.strandedness),
            tick_track_channel_override=tick_track_channel_override,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            length_reference_radius_px=base_radius,
        )
        tick_band = RadialBand(tick_inner, tick_outer)

    label_band: RadialBand | None = None
    label_bounds = get_circular_tick_label_radius_bounds(
        center_radius_px=anchor,
        total_len=total_length,
        track_type=str(cfg.canvas.circular.track_type),
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
        anchor_radius_px=anchor,
        tick_band_px=tick_band,
        label_band_px=label_band,
        reserved_band_px=reserved,
    )


def _default_numeric_width_px(
    renderer: str,
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    length_param = str(canvas_config.length_param)
    base = float(canvas_config.radius) * float(canvas_config.track_ratio)
    if renderer == "depth":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][1]) * 0.5
    if renderer == "dinucleotide_skew":
        return base * float(cfg.canvas.circular.track_ratio_factors[length_param][2])
    return base * float(cfg.canvas.circular.track_ratio_factors[length_param][1])


def _slot_intents(
    slots: Sequence[CircularTrackSlot],
    *,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    track_specs: Sequence[TrackSpec] | None,
) -> list[_SlotIntent]:
    base_radius = float(canvas_config.radius)
    intents: list[_SlotIntent] = []
    for slot in slots:
        renderer = str(slot.renderer)
        if not slot.enabled or renderer in {"features", "ticks"}:
            continue
        ts = _track_spec_for_slot(slot, track_specs)
        if ts is not None and not ts.show:
            continue

        params = dict(slot.params)
        if ts is not None and ts.params:
            params.update(ts.params)

        center_px, width_px, explicit_anchor, explicit_annulus, explicit_width = _resolve_track_spec_center_and_width(
            ts,
            base_radius_px=base_radius,
        )
        slot_center, slot_width, slot_explicit_anchor, slot_explicit_annulus, slot_explicit_width = _resolve_placement_center_and_width(
            slot.placement,
            base_radius_px=base_radius,
        )
        if slot_center is not None:
            center_px = slot_center
            explicit_anchor = slot_explicit_anchor
            explicit_annulus = slot_explicit_annulus
        if slot_width is not None:
            width_px = slot_width
            explicit_width = slot_explicit_width
        if slot.width is not None:
            width_px = float(slot.width.resolve(base_radius))
            explicit_width = True
        if width_px is None:
            width_px = _default_numeric_width_px(renderer, canvas_config=canvas_config, cfg=cfg)

        channel = _normalize_side(params.get("side"), default="inside")
        intents.append(
            _SlotIntent(
                slot=slot,
                center_px=center_px,
                width_px=max(0.0, float(width_px)),
                explicit_anchor=bool(explicit_anchor),
                explicit_annulus=bool(explicit_annulus),
                explicit_width=bool(explicit_width),
                channel=channel,
                params=params,
            )
        )
    return intents


def _slot_gap_px(slot: CircularTrackSlot, *, base_radius_px: float) -> float:
    if slot.gap_after is None:
        return max(1.0, 0.01 * float(base_radius_px))
    return max(0.0, float(slot.gap_after.resolve(float(base_radius_px))))


def _inside_start_outer(occupied: Sequence[RadialBand], *, guard_px: float, axis_radius_px: float) -> float:
    candidates = [
        float(band.inner_px)
        for band in occupied
        if band.width_px > LAYOUT_EPSILON and float(band.inner_px) > float(guard_px) + LAYOUT_EPSILON
    ]
    if candidates:
        return max(candidates)
    return max(float(guard_px), float(axis_radius_px))


def _free_intervals(
    occupied: Sequence[RadialBand],
    *,
    inner_limit_px: float,
    outer_limit_px: float,
) -> list[tuple[float, float]]:
    lower_limit = max(0.0, float(inner_limit_px))
    upper_limit = max(lower_limit, float(outer_limit_px))
    intervals: list[tuple[float, float]] = []
    cursor = lower_limit
    for band in sorted(occupied, key=lambda item: item.inner_px):
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


def _place_inside(
    intent: _SlotIntent,
    *,
    occupied: Sequence[RadialBand],
    cursor_outer_px: float,
    guard_px: float,
    base_radius_px: float,
) -> tuple[RadialBand, float]:
    gap_px = _slot_gap_px(intent.slot, base_radius_px=base_radius_px)
    width_px = max(0.0, float(intent.width_px))
    allow_compress = (
        _parse_bool_param(intent.params.get("compress"), default=False)
        or (
            str(intent.slot.renderer) in NUMERIC_RENDERERS
            and not intent.explicit_width
            and not intent.explicit_annulus
        )
    )
    best_compressed: tuple[RadialBand, float] | None = None
    for raw_lower, raw_upper in sorted(
        _free_intervals(occupied, inner_limit_px=guard_px, outer_limit_px=cursor_outer_px),
        key=lambda interval: interval[1],
        reverse=True,
    ):
        lower = float(raw_lower) + gap_px
        upper = float(raw_upper) - gap_px
        available = max(0.0, upper - lower)
        if available <= LAYOUT_EPSILON:
            continue
        if width_px <= available + LAYOUT_EPSILON:
            band = RadialBand(upper - width_px, upper)
            return band, float(band.inner_px)
        if allow_compress and available > LAYOUT_EPSILON:
            compressed_width = available
            candidate = RadialBand(upper - compressed_width, upper)
            if best_compressed is None or candidate.outer_px > best_compressed[0].outer_px:
                best_compressed = (candidate, compressed_width)

    if best_compressed is not None:
        band, compressed_width = best_compressed
        logger.info(
            "Auto-compressed circular track '%s' to %.3fpx to fit reserved radial space.",
            intent.slot.id,
            compressed_width,
        )
        return band, float(band.inner_px)

    raise ValueError(
        f"circular track slot '{intent.slot.id}' cannot fit without overlapping the center definition"
    )


def _resolve_tracks(
    intents: Sequence[_SlotIntent],
    *,
    reference_bands: Sequence[RadialBand],
    definition_band: RadialBand | None,
    axis_radius_px: float,
    base_radius_px: float,
) -> tuple[CircularResolvedTrack, ...]:
    occupied: list[RadialBand] = [
        band for band in reference_bands if band.width_px > LAYOUT_EPSILON
    ]
    if definition_band is not None and definition_band.width_px > LAYOUT_EPSILON:
        occupied.append(definition_band)

    resolved_by_id: dict[str, CircularResolvedTrack] = {}

    for intent in intents:
        if intent.center_px is None:
            continue
        draw_band = _band_from_center_width(float(intent.center_px), float(intent.width_px))
        reserved_band = draw_band
        track = CircularResolvedTrack(
            id=str(intent.slot.id),
            renderer=str(intent.slot.renderer),
            channel=intent.channel,
            anchor_radius_px=draw_band.center_px,
            draw_band_px=draw_band,
            reserved_band_px=reserved_band,
            z=int(intent.slot.z),
            params=dict(intent.params),
            explicit_width=bool(intent.explicit_width),
        )
        strict = _parse_bool_param(intent.params.get("strict"), default=False)
        for band in occupied:
            if bands_overlap(reserved_band, band):
                message = f"Pinned circular track slot '{intent.slot.id}' overlaps a reserved circular layout band."
                if strict:
                    raise ValueError(message)
                logger.warning(message)
                break
        resolved_by_id[str(intent.slot.id)] = track
        if intent.channel != "overlay" or _parse_bool_param(intent.params.get("reserve"), default=False):
            occupied.append(reserved_band)

    guard_px = float(definition_band.outer_px) if definition_band is not None else 0.0
    inside_cursor = _inside_start_outer(occupied, guard_px=guard_px, axis_radius_px=axis_radius_px)
    outside_cursor = max([float(band.outer_px) for band in occupied] + [float(axis_radius_px)])

    for intent in intents:
        if str(intent.slot.id) in resolved_by_id:
            continue
        if intent.channel == "outside":
            gap_px = _slot_gap_px(intent.slot, base_radius_px=base_radius_px)
            draw_band = RadialBand(outside_cursor + gap_px, outside_cursor + gap_px + float(intent.width_px))
            outside_cursor = draw_band.outer_px
        elif intent.channel == "overlay":
            center = float(axis_radius_px)
            draw_band = _band_from_center_width(center, float(intent.width_px))
        else:
            try:
                draw_band, inside_cursor = _place_inside(
                    intent,
                    occupied=occupied,
                    cursor_outer_px=inside_cursor,
                    guard_px=guard_px,
                    base_radius_px=base_radius_px,
                )
            except ValueError:
                if _parse_bool_param(intent.params.get("strict"), default=False):
                    raise
                gap_px = _slot_gap_px(intent.slot, base_radius_px=base_radius_px)
                draw_band = RadialBand(outside_cursor + gap_px, outside_cursor + gap_px + float(intent.width_px))
                outside_cursor = draw_band.outer_px
                logger.warning(
                    "Placed circular track slot '%s' outside because the inside channel has no free radial space.",
                    intent.slot.id,
                )

        reserved_band = draw_band
        track = CircularResolvedTrack(
            id=str(intent.slot.id),
            renderer=str(intent.slot.renderer),
            channel=intent.channel,
            anchor_radius_px=draw_band.center_px,
            draw_band_px=draw_band,
            reserved_band_px=reserved_band,
            z=int(intent.slot.z),
            params=dict(intent.params),
            explicit_width=bool(intent.explicit_width),
        )
        resolved_by_id[str(intent.slot.id)] = track
        if intent.channel != "overlay" or _parse_bool_param(intent.params.get("reserve"), default=False):
            occupied.append(reserved_band)

    return tuple(resolved_by_id[str(intent.slot.id)] for intent in intents if str(intent.slot.id) in resolved_by_id)


def resolve_circular_radial_layout(
    *,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    slots: Sequence[CircularTrackSlot],
    track_specs: Sequence[TrackSpec] | None = None,
    feature_dict: Mapping[str, Any] | None = None,
    show_features: bool = True,
    show_ticks: bool = True,
    definition_reserved_radius_px: float | None = None,
    feature_track_ratio_factor_override: float | None = None,
    tick_track_channel_override: str | None = None,
) -> CircularRadialLayout:
    axis_ts = _track_specs_by_match(track_specs)[1].get("axis")
    axis_radius_px = float(canvas_config.radius)
    axis_center_px, _axis_width_px, axis_explicit, _axis_annulus, _axis_explicit_width = _resolve_track_spec_center_and_width(
        axis_ts,
        base_radius_px=float(canvas_config.radius),
    )
    if axis_explicit and axis_center_px is not None:
        axis_radius_px = float(axis_center_px)

    axis = CircularAxisLayout(
        radius_px=axis_radius_px,
        stroke_width_px=float(cfg.objects.axis.circular.stroke_width.for_length_param(str(canvas_config.length_param))),
    )

    slots_by_renderer = {
        str(slot.renderer): slot
        for slot in slots
        if slot.enabled and str(slot.renderer) in {"features", "ticks"}
    }
    ts_by_kind = _track_specs_by_match(track_specs)[1]
    features_ts = ts_by_kind.get("features")
    ticks_ts = ts_by_kind.get("ticks")

    feature_layout: CircularFeatureLayout | None = None
    if show_features and (features_ts is None or features_ts.show):
        feature_width_px, _ratio_factor = _feature_width_from_specs(
            feature_slot=slots_by_renderer.get("features"),
            feature_ts=features_ts,
            canvas_config=canvas_config,
            cfg=cfg,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        )
        feature_anchor_px = _feature_anchor_from_slot(
            slots_by_renderer.get("features"),
            base_radius_px=axis_radius_px,
        )
        feature_layout = build_circular_feature_layout(
            feature_dict,
            axis_radius_px=axis_radius_px,
            width_px=feature_width_px,
            track_type=str(cfg.canvas.circular.track_type),
            strandedness=bool(cfg.canvas.strandedness),
            anchor_radius_px=feature_anchor_px,
        )

    tick_layout: CircularTickLayout | None = None
    if show_ticks and (ticks_ts is None or ticks_ts.show):
        tick_layout = _tick_layout_from_slot(
            axis_radius_px=axis_radius_px,
            total_length=int(total_length),
            canvas_config=canvas_config,
            cfg=cfg,
            tick_slot=slots_by_renderer.get("ticks"),
            tick_ts=ticks_ts,
            tick_track_channel_override=tick_track_channel_override,
        )

    definition_band = (
        RadialBand(0.0, float(definition_reserved_radius_px))
        if definition_reserved_radius_px is not None and definition_reserved_radius_px > LAYOUT_EPSILON
        else None
    )
    reference_bands: list[RadialBand] = []
    if feature_layout is not None:
        reference_bands.append(feature_layout.all_band_px)
    if tick_layout is not None:
        reference_bands.append(tick_layout.reserved_band_px)

    tracks = _resolve_tracks(
        _slot_intents(slots, canvas_config=canvas_config, cfg=cfg, track_specs=track_specs),
        reference_bands=reference_bands,
        definition_band=definition_band,
        axis_radius_px=axis_radius_px,
        base_radius_px=float(canvas_config.radius),
    )

    outer_content_radius = max(
        [axis_radius_px]
        + [band.outer_px for band in reference_bands]
        + [track.reserved_band_px.outer_px for track in tracks]
    )

    return CircularRadialLayout(
        axis=axis,
        features=feature_layout,
        ticks=tick_layout,
        tracks=tuple(tracks),
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
    "CircularResolvedTrack",
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
