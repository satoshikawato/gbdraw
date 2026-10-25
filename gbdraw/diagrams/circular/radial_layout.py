"""Unified radial slot resolver for circular diagrams."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, replace
from typing import Any, Callable, Collection, Literal, Mapping, Protocol, Sequence, TypedDict, TypeVar

from ...annotations.models import AnnotationTrackParams
from ...auto_sizes import circular_track_width_px
from ...canvas import CircularCanvasGeometry
from ...config.models import GbdrawConfig
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
    DepthAxisSource,
    DepthAxisStyle,
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
from .presets import (
    circular_feature_lane_direction_for_preset,
    circular_preset_for_feature_lane_direction,
    normalize_circular_track_preset,
)


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
# The band a span overlay's labels start from, on the slot that holds their room: ``to_slot`` stretches over the
# slots' whole bands ("whole"), ``cover_anchor`` over the anchor's draw band ("draw"), and an overlay that spans
# nothing draws its own width on the anchor's draw band ("own").
SpanExtent = Literal["whole", "draw", "own"]
# How far an annotation slot's labels reach beyond its band, gap included: for the band, the axis, the slot's side and
# params, and whether the band is a span over other slots.
AnnotationLabelRoom = Callable[[RadialBand, CircularAxisLayout, str, Mapping[str, Any], bool], float]


@dataclass(frozen=True)
class _SwitchedTrack:
    """A track row with an off option, in the words a fit failure offers it with (OV-1274, OV-1276)."""

    # The diagnostic identifier (``keptInsideTracks``, ``slotTrack``; utils/error-normalization.js words it).
    diagnostic_id: str
    words: str
    # How the CLI turns the row off, as the remedy names it.
    off_remedy: str
    # The option that narrows the row, or None.
    width_option: str | None = None


# The rows with an off option, by renderer. Other rows are named by slot ID only.
_SWITCHED_TRACKS: Mapping[str, _SwitchedTrack] = {
    "dinucleotide_content": _SwitchedTrack("gc_content", "GC content", "--no-gc", "--gc_content_width"),
    "dinucleotide_skew": _SwitchedTrack("gc_skew", "GC skew", "--no-skew", "--gc_skew_width"),
    "ticks": _SwitchedTrack("ticks", "scale", "--hide_scale; with an explicit track list, remove the ticks slot"),
}
# Their diagnostic identifiers: the vocabulary of the Web error adapter.
SWITCHED_TRACK_IDS: frozenset[str] = frozenset(track.diagnostic_id for track in _SWITCHED_TRACKS.values())


def _switched_track_remedy(tracks: Sequence[_SwitchedTrack]) -> str:
    """Turning ``tracks`` off (and narrowing those with a width option), as list items ending in ", ", or ""."""

    narrowable = [(track, track.width_option) for track in tracks if track.width_option is not None]
    remedy = ""
    if narrowable:
        words = " and ".join(track.words for track, _width in narrowable)
        options = ", ".join([*(track.off_remedy for track, _width in narrowable), *(width for _track, width in narrowable)])
        remedy = f"turn off or narrow the {words} track{'s' if len(narrowable) > 1 else ''} ({options}), "
    return remedy + "".join(f"hide the {track.words} ({track.off_remedy}), " for track in tracks if track.width_option is None)


# The owner of the center reservation band among the occupied bands: an explicit ``center_reserved_radius``
# (DF-D05), or the centre definition text's band in the "move" pass of rule L (``_resolve_around_definition``).
_CENTER_RESERVATION = "definition"

# The passes of rule L (PD-OI-078) that place the tracks around the centre definition text's band:
# "move" keeps the inside tracks out of the band, moving them outward and narrowing Auto widths down to the floor
# (``canvas.circular.track_narrowing_min_px``); "limit" packs the inside rows on preset lanes outward at their own
# widths and draws the definition over them; "none" only records the band (no definition, or the tracks clear it).
_DefinitionPass = Literal["none", "move", "limit"]

# Inside rows that are not tracks of the inside stack: the outermost stack track starts right under them.
_ROWS_ABOVE_THE_INSIDE_STACK = frozenset({"features", "ticks"})


@dataclass(frozen=True)
class _KeptInsideSlot:
    slot_id: str
    renderer: str


@dataclass(frozen=True)
class PlacementWindow:
    inner_px: float
    outer_px: float
    # The reserved band (tick labels) may reach below inner_px to this edge of
    # a pinned row, as it may touch the row below it in top-down packing.
    reserved_inner_px: float | None = None
    # The room kept below inner_px for the inside rows placed after this window's rows, and those rows.
    kept_inside_px: float = 0.0
    kept_inside_slots: tuple[_KeptInsideSlot, ...] = ()
    # The outer edge of the center reservation (an explicit center_reserved_radius, or the definition band in the
    # "move" pass). Placement reads it from
    # the occupied bands; a fit failure reports the window from it (OV-1276).
    center_reserved_px: float = 0.0

    @property
    def usable_inner_px(self) -> float:
        """The lower edge of the band the window's rows may use: inner_px, or the center reservation above it."""

        return max(float(self.inner_px), float(self.center_reserved_px))

    @property
    def kept_switched_tracks(self) -> tuple[_SwitchedTrack, ...]:
        """The narrowable rows (GC) among the rows kept inside, each once, in slot order; the scale is offered only when it fails."""

        if self.kept_inside_px <= LAYOUT_EPSILON:
            return ()
        tracks = (_SWITCHED_TRACKS.get(slot.renderer) for slot in self.kept_inside_slots)
        return tuple(dict.fromkeys(track for track in tracks if track is not None and track.width_option is not None))


@dataclass(frozen=True)
class _SpanLabelRoom:
    """The labels' room of a span overlay past the slot that holds it (``_span_room_holder``), for a band of it."""

    room: Callable[[RadialBand], float]
    extent: SpanExtent
    # The overlay's own width (its intent's, never negative), the band its labels start from when it spans nothing
    # ("own"); the final overlay loop draws the same width.
    width_px: float


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
    # The room of an annotation slot's labels on their side, for a band (``AnnotationLabelRoom`` on this axis).
    label_room: Callable[[RadialBand], float] | None = None
    # The rooms of the span overlays whose outermost covered slot this is: reserved past its band while it is packed.
    span_label_rooms: tuple[_SpanLabelRoom, ...] = ()
    # The narrowest width (px) an Auto width may take: 0 keeps the readable minimums alone; the "move" pass of rule L
    # sets ``canvas.circular.track_narrowing_min_px`` (infinite with ``track_narrowing`` off).
    narrowing_floor_px: float = 0.0
    # An explicit center_reserved_radius narrows a typed width on a preset lane, as before rule L (DF-D38).
    typed_lane_width_narrows: bool = False


def _cannot_fit_diagnostic(
    intent: "_RadialSlotIntent | None", window: PlacementWindow | None = None
) -> dict[str, object]:
    """Track row and usable band of a fit failure; the slot ID stays private.

    A pinned slot has no band to report, and an unbounded band is omitted. The
    room kept below the band and the rows with an off option it is kept for (OV-1274) are named.
    """

    diagnostic: dict[str, object] = {"code": "TRACK_LAYOUT", "reason": "CANNOT_FIT"}
    if window is not None and float(window.outer_px) < float("inf"):
        diagnostic["innerPx"] = max(0, round(window.usable_inner_px))
        diagnostic["outerPx"] = max(0, round(float(window.outer_px)))
    if window is not None and window.kept_inside_px > LAYOUT_EPSILON:
        diagnostic["keptInsidePx"] = round(float(window.kept_inside_px))
        kept_tracks = window.kept_switched_tracks
        if kept_tracks:
            diagnostic["keptInsideTracks"] = [track.diagnostic_id for track in kept_tracks]
    if intent is not None:
        diagnostic["slotIndex"] = int(intent.slot_index)
    return diagnostic


def track_layout_reason(error: ValidationError) -> str | None:
    """The ``TRACK_LAYOUT`` diagnostic reason of a placement failure, or None for any other error."""

    diagnostic = error.diagnostic or {}
    if diagnostic.get("code") != "TRACK_LAYOUT":
        return None
    reason = diagnostic.get("reason")
    return reason if isinstance(reason, str) else None


def _center_reserved_diagnostic(diagnostic: Mapping[str, object]) -> dict[str, object]:
    """The same fit failure, limited by an explicit ``center_reserved_radius`` (``CENTER_RESERVED``)."""

    return {**diagnostic, "code": "TRACK_LAYOUT", "reason": "CENTER_RESERVED"}


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
        else circular_preset_for_feature_lane_direction(lane_direction)
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
        else circular_feature_lane_direction_for_preset(track_type)
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    params: Mapping[str, Any],
    anchor_radius_px: float,
    width_px: float,
    explicit_width: bool,
) -> CircularTickLayout:
    base_radius = float(geometry.radius)
    tick_preset = normalize_circular_track_preset(str(params.get("preset", params.get("track_preset", "tuckin"))))
    tick_label_layout = str(params.get("tick_label_layout", "label_out_tick_in")).strip().lower()
    label_side, tick_side = tick_sides_for_tick_label_layout(
        tick_label_layout,
        side=params.get("_slot_side"),
    )
    tick_length_px = float(width_px) if explicit_width and width_px > 0 else None
    length_threshold = int(cfg.labels.length_threshold.circular)

    if tick_side in {"none", ""}:
        tick_band = RadialBand(anchor_radius_px, anchor_radius_px)
    else:
        tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
            center_radius_px=float(anchor_radius_px),
            total_len=int(total_length),
            size="large",
            track_type=tick_preset,
            strandedness=bool(cfg.canvas.strandedness),
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            length_reference_radius_px=base_radius,
            length_threshold=length_threshold,
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
        manual_interval=cfg.objects.scale.interval,
        label_side=label_side,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
        tick_width=float(cfg.objects.ticks.tick_width),
        length_reference_radius_px=base_radius,
        length_threshold=length_threshold,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
) -> float:
    if renderer == "ticks":
        return 0.0
    return circular_track_width_px(
        renderer,
        radius=geometry.radius,
        track_ratio=geometry.track_ratio,
        factors=geometry.track_ratio_factors,
    )


def _default_spacing_px(axis_radius_px: float) -> float:
    return max(1.0, 0.01 * float(axis_radius_px))


def _slot_intents(
    slots: Sequence[CircularTrackSlot],
    *,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    preferred_anchor_slot_ids: Collection[str] = (),
    axis: CircularAxisLayout,
    annotation_label_rooms: Mapping[str, AnnotationLabelRoom],
) -> list[_RadialSlotIntent]:
    axis_radius_px = float(geometry.radius)
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
            width_px = _default_width_px(slot.renderer, geometry=geometry, cfg=cfg)
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
                label_room=_bound_label_room(annotation_label_rooms.get(slot.id), axis, slot)
                if slot.renderer == "annotations" and slot.side != "overlay"
                else None,
            )
        )
    span_rooms = _span_label_rooms_by_slot(intents, axis, annotation_label_rooms)
    return [
        replace(intent, span_label_rooms=tuple(span_rooms[intent.slot_id])) if intent.slot_id in span_rooms else intent
        for intent in intents
    ]


def _span_label_rooms_by_slot(
    intents: Sequence[_RadialSlotIntent], axis: CircularAxisLayout, rooms: Mapping[str, AnnotationLabelRoom]
) -> dict[str, list[_SpanLabelRoom]]:
    """Each span overlay's label room, keyed by the slot that holds it (``_span_room_holder``)."""

    by_id = {intent.slot_id: intent for intent in intents}
    span_rooms: dict[str, list[_SpanLabelRoom]] = {}
    for intent in intents:
        annotation = intent.slot.annotation
        room = rooms.get(intent.slot_id)
        if intent.renderer != "annotations" or intent.side != "overlay" or annotation is None or room is None:
            continue
        holder = _span_room_holder(annotation, by_id)
        if holder is None:
            continue
        extent = _span_extent(annotation)
        span_rooms.setdefault(holder.slot_id, []).append(
            _SpanLabelRoom(_bound_span_label_room(room, axis, intent.slot, extent), extent, float(intent.width_px))
        )
    return span_rooms


def _span_extent(annotation: AnnotationTrackParams) -> SpanExtent:
    if annotation.to_slot is not None:
        return "whole"
    return "draw" if annotation.cover_anchor else "own"


class _SpanCoveredSlot(Protocol):
    @property
    def side(self) -> str: ...

    @property
    def renderer(self) -> str: ...

    @property
    def slot_index(self) -> int: ...


_SpanCoveredSlotT = TypeVar("_SpanCoveredSlotT", bound=_SpanCoveredSlot)


def _span_room_holder(
    annotation: AnnotationTrackParams, by_id: Mapping[str, _SpanCoveredSlotT]
) -> _SpanCoveredSlotT | None:
    """The slot whose reserved band holds a span overlay's labels' room: the outermost slot it covers, or its anchor
    when it spans nothing (OD-R1174 F1). An ``outside`` slot lies beyond every other, split feature lanes (a features
    overlay) lie on the axis, beyond every ``inside`` slot, and on one side the slot listed first lies outermost.
    Another covered overlay (a plot overlay) holds none."""

    covered = [
        by_id[slot_id]
        for slot_id in (annotation.anchor_slot, annotation.to_slot)
        if slot_id is not None
        and slot_id in by_id
        and (by_id[slot_id].side != "overlay" or by_id[slot_id].renderer == "features")
    ]
    rank = {"outside": 0, "overlay": 1}
    return min(covered, key=lambda item: (rank.get(item.side, 2), item.slot_index)) if covered else None


def _span_room_band(resolved: CircularResolvedSlot, extent: SpanExtent, width_px: float) -> RadialBand | None:
    """The band of a span's room holder that the room lies past: its whole band for ``to_slot``, its draw band for
    ``cover_anchor``, and the overlay's own width on its draw band for an overlay that spans nothing. Without
    ``to_slot`` the overlay draws this band on its anchor."""

    if extent == "whole":
        return resolved.packing_band_px or resolved.draw_band_px
    band = resolved.draw_band_px or resolved.packing_band_px
    if band is None:
        return None
    return _band_from_center_width(band.center_px, band.width_px if extent == "draw" else width_px)


def _bound_span_label_room(
    room: AnnotationLabelRoom, axis: CircularAxisLayout, overlay: NormalizedCircularTrackSlot, extent: SpanExtent
) -> Callable[[RadialBand], float]:
    # An overlay that spans nothing places its labels from its own mark, as the drawer does (``spans_slots`` unset).
    covers_slots = extent != "own"
    return lambda band: room(band, axis, overlay.side, overlay.params, covers_slots)


def _with_span_label_room(intent: _RadialSlotIntent, resolved: CircularResolvedSlot) -> CircularResolvedSlot:
    """``resolved`` with the labels' room of the spans it holds (``_span_room_holder``) reserved past its band.

    The room is measured on the band the span stretches over from this slot (``_span_room_band``). A ``to_slot``
    span also covers slots further in; its letters' reach past the edge does not grow with the span's width, so this
    slot's own band bounds it. A span whose labels need no room (none shown) reserves nothing, even when it draws
    wider than this slot.
    """

    reserved = resolved.reserved_band_px
    reach = 0.0
    for span in intent.span_label_rooms:
        band = _span_room_band(resolved, span.extent, span.width_px)
        if band is None:
            continue
        room = span.room(band)
        if room > LAYOUT_EPSILON:
            reach = max(reach, float(band.outer_px) + room)
    if reserved is None or reach <= float(reserved.outer_px) + LAYOUT_EPSILON:
        return resolved
    return replace(
        resolved,
        reserved_band_px=RadialBand(reserved.inner_px, reach),
        span_label_room_px=reach - float(reserved.outer_px),
    )


def _bound_label_room(
    room: AnnotationLabelRoom | None, axis: CircularAxisLayout, slot: NormalizedCircularTrackSlot
) -> Callable[[RadialBand], float] | None:
    return None if room is None else lambda band: room(band, axis, slot.side, slot.params, False)


def _measure_radial_slot(
    intent: _RadialSlotIntent,
    *,
    anchor_offset_px: float,
    width_px: float | None = None,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
        draw_band: RadialBand | None = (
            feature_layout.primary_band_px if feature_layout is not None else band
        )
        return _with_span_label_room(intent, CircularResolvedSlot(
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
        ))

    if renderer == "ticks":
        tick_params = dict(intent.params)
        tick_params["_slot_side"] = intent.side
        tick_layout = _tick_layout_from_params(
            total_length=int(total_length),
            geometry=geometry,
            cfg=cfg,
            params=tick_params,
            anchor_radius_px=anchor_radius_px,
            width_px=resolved_width,
            explicit_width=bool(intent.explicit_width),
        )
        return _with_span_label_room(intent, CircularResolvedSlot(
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
        ))

    band = _band_from_center_width(anchor_radius_px, resolved_width)
    draw_band = None if renderer == "spacer" else band
    reserved_band = band
    if intent.label_room is not None:
        # An annotation slot reserves its labels' room on their side, so other rings keep clear of the letters.
        room = intent.label_room(band)
        reserved_band = (
            RadialBand(band.inner_px - room, band.outer_px)
            if intent.side == "inside"
            else RadialBand(band.inner_px, band.outer_px + room)
        )
    if renderer == "depth" and depth_config is not None:
        reserved_band = _depth_reserved_band_for_draw_band(
            band,
            resolve_depth_axis_footprint(depth_config, resolved_width),
        )
    return _with_span_label_room(intent, CircularResolvedSlot(
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
    ))


def _reserved_overlap_any(band: RadialBand, occupied: Sequence[tuple[str, RadialBand]]) -> tuple[str, RadialBand] | None:
    for owner, other in occupied:
        if bands_overlap(band, other):
            return owner, other
    return None


def _floored_min_width_px(intent: _RadialSlotIntent, min_width_px: float) -> float:
    """A readable minimum raised to the intent's narrowing floor; a track narrower than the floor keeps its width."""

    return max(float(min_width_px), min(max(0.0, float(intent.width_px)), float(intent.narrowing_floor_px)))


def _min_readable_numeric_width_px(intent: _RadialSlotIntent) -> float:
    width = max(0.0, float(intent.width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if intent.renderer == "sequence_conservation":
        minimum = min(width, max(MIN_CONSERVATION_WIDTH_PX, MIN_CONSERVATION_WIDTH_FRACTION * width))
    elif intent.renderer == "dinucleotide_skew":
        minimum = min(width, max(MIN_SKEW_WIDTH_PX, MIN_SKEW_WIDTH_FRACTION * width))
    else:
        minimum = min(width, max(MIN_NUMERIC_WIDTH_PX, MIN_NUMERIC_WIDTH_FRACTION * width))
    return _floored_min_width_px(intent, minimum)


def _min_readable_preferred_width_px(intent: _RadialSlotIntent) -> float:
    width = max(0.0, float(intent.width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if intent.renderer == "sequence_conservation":
        minimum = min(width, max(MIN_CONSERVATION_WIDTH_PX, PREFERRED_MIN_NUMERIC_WIDTH_FRACTION * width))
    elif intent.renderer == "dinucleotide_skew":
        minimum = min(width, max(MIN_SKEW_WIDTH_PX, PREFERRED_MIN_SKEW_WIDTH_FRACTION * width))
    else:
        minimum = min(width, max(MIN_NUMERIC_WIDTH_PX, PREFERRED_MIN_NUMERIC_WIDTH_FRACTION * width))
    return _floored_min_width_px(intent, minimum)


def _min_dense_stack_numeric_width_px(intent: _RadialSlotIntent) -> float:
    width = max(0.0, float(intent.width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if intent.renderer == "sequence_conservation":
        minimum = min(width, MIN_DENSE_CONSERVATION_WIDTH_PX)
    elif intent.renderer == "dinucleotide_skew":
        minimum = min(width, MIN_SKEW_WIDTH_PX)
    else:
        minimum = min(width, MIN_NUMERIC_WIDTH_PX)
    return _floored_min_width_px(intent, minimum)


def _candidate_widths(intent: _RadialSlotIntent) -> list[tuple[float, bool]]:
    width = max(0.0, float(intent.width_px))
    if not intent.compress or intent.renderer not in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        return [(width, False)]
    min_width = _min_readable_numeric_width_px(intent)
    if min_width >= width - LAYOUT_EPSILON:
        return [(width, False)]
    values = [width]
    steps = 8
    for idx in range(1, steps + 1):
        fraction = idx / float(steps)
        values.append(width - ((width - min_width) * fraction))
    return [(max(0.0, value), idx > 0) for idx, value in enumerate(values)]


def _occupied_band(intent: _RadialSlotIntent, resolved: CircularResolvedSlot) -> RadialBand | None:
    """The band other slots keep out of: the reserved band, and for an overlay (split feature lanes, a plot overlay)
    its own band widened by its gaps, as no row is packed next to it. The span label room past it is not widened:
    it ends ``label_offset`` beyond the letters, as past any other slot."""

    reserved = resolved.reserved_band_px
    own = resolved.own_reserved_band_px
    if intent.side != "overlay" or reserved is None or own is None:
        return reserved
    widened = own.expanded_sides(intent.inner_gap_px, intent.outer_gap_px)
    return RadialBand(widened.inner_px, max(widened.outer_px, float(reserved.outer_px)))


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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
                geometry=geometry,
                cfg=cfg,
                total_length=total_length,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
                geometry=geometry,
                cfg=cfg,
                total_length=total_length,
                depth_config=depth_config,
                compressed=compressed,
            )
            bands = [band for band in (resolved.packing_band_px, resolved.own_reserved_band_px) if band is not None]
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
            reserved = resolved.reserved_band_px
            conflict = _reserved_overlap_any(reserved, occupied) if reserved is not None else None
            if reserved is not None and conflict is not None:
                room_inner = float(reserved.outer_px) - resolved.span_label_room_px
                if resolved.span_label_room_px > 0.0 and float(conflict[1].inner_px) >= room_inner - LAYOUT_EPSILON:
                    # Only the span labels' room meets the row above: the slot moves in by the overlap.
                    anchor_offset -= float(reserved.outer_px) - float(conflict[1].inner_px)
                    continue
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
) -> CircularResolvedSlot | None:
    """A pinned row at ``width_px`` on its radius, or ``None`` when it leaves the window or overlaps."""

    resolved = _measure_radial_slot(
        intent,
        anchor_offset_px=float(intent.anchor_offset_px or 0.0),
        width_px=width_px,
        axis_radius_px=axis_radius_px,
        feature_dict=feature_dict,
        geometry=geometry,
        cfg=cfg,
        total_length=total_length,
        depth_config=depth_config,
        compressed=compressed,
    )
    bands = [band for band in (resolved.packing_band_px, resolved.own_reserved_band_px) if band is not None]
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
    *,
    keep_preset_lanes: bool,
) -> list[_RadialSlotIntent]:
    """Auto rows and pinned Auto-width rows that pack as one stack.

    With ``keep_preset_lanes`` the group ends where a preset-lane group starts
    (GX-19); without it, preset-lane rows without a radius reflow in the stack.
    """

    group: list[_RadialSlotIntent] = []
    for future_pos in range(start_pos, len(ordered_intents)):
        future = ordered_intents[future_pos]
        if future.slot_index in resolved_by_index or future.side != "inside":
            break
        if keep_preset_lanes and _preferred_numeric_group_from(
            ordered_intents, future_pos, resolved_by_index, keep_preset_lanes=True
        ):
            break
        reflowing_preset_lane = (
            not keep_preset_lanes and _preferred_inside_numeric(future) and not future.explicit_anchor
        )
        if not (future.placement_policy == "auto" or _anchored_inside(future) or reflowing_preset_lane):
            break
        group.append(future)
    return group


def _inside_auto_stack_width_scales(intents: Sequence[_RadialSlotIntent]) -> list[float]:
    min_scales = [
        _min_dense_stack_numeric_width_px(intent) / intent.width_px
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
    min_width = _min_dense_stack_numeric_width_px(intent)
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
        # The failed row, when it has an off option (OV-1276: the scale ticks).
        self.slot_track = _SWITCHED_TRACKS.get(intent.renderer) if intent is not None else None
        kept = self._kept_inside()
        diagnostic = _cannot_fit_diagnostic(intent, placement_window)
        if self.slot_track is not None:
            diagnostic["slotTrack"] = self.slot_track.diagnostic_id
        super().__init__(
            f"{self._between()}{f'; below it, {kept}' if kept else ''}. Move the slot, reduce widths, "
            f"disable conflicting labels, {self._switch_remedy()}or use side=outside.{hint}",
            diagnostic=diagnostic,
        )

    def _between(self) -> str:
        return (
            f"{self.subject} cannot fit inside between "
            f"{self.placement_window.usable_inner_px:.1f}px and {self.placement_window.outer_px:.1f}px"
        )

    def _kept_inside(self) -> str:
        """The room the window keeps for the inside rows after it (OV-1268), or ""."""

        window = self.placement_window
        if window.kept_inside_px <= LAYOUT_EPSILON:
            return ""
        slot_ids = ", ".join(f"'{slot.slot_id}'" for slot in window.kept_inside_slots)
        return f"{window.kept_inside_px:.1f}px are kept for the tracks inside it ({slot_ids})"

    def _switch_remedy(self) -> str:
        """Turning off the failed row or the rows the window keeps room for (OV-1274, OV-1276), or ""."""

        failed = () if self.slot_track is None else (self.slot_track,)
        return _switched_track_remedy(tuple(dict.fromkeys((*failed, *self.placement_window.kept_switched_tracks))))

    def limited_by_center(self, reserved_radius_px: float) -> ValidationError:
        """Name an explicit ``center_reserved_radius`` as the cause."""

        kept = self._kept_inside()
        reserved = f"{float(reserved_radius_px):.1f}px{f' and {kept}' if kept else ''}"
        switch_remedy = self._switch_remedy()
        message = (
            f"{self._between()} because center_reserved_radius reserves {reserved}. "
            f"Set a smaller center_reserved_radius{f', {switch_remedy}' if switch_remedy else ' '}"
            f"or place tracks outside.{self.hint}"
        )
        return ValidationError(message, diagnostic=_center_reserved_diagnostic(self.diagnostic or {}))


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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
                geometry=geometry,
                cfg=cfg,
                total_length=total_length,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
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
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
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


def _lane_width_narrows(intent: _RadialSlotIntent) -> bool:
    """Whether a row on a preset lane may narrow: an Auto width, or any width under an explicit center_reserved_radius."""

    return intent.compress or intent.typed_lane_width_narrows


def _preferred_group_width_scales(intents: Sequence[_RadialSlotIntent]) -> list[float]:
    min_scales = [
        _min_readable_preferred_width_px(intent) / intent.width_px
        for intent in intents
        if _lane_width_narrows(intent) and intent.width_px > LAYOUT_EPSILON
    ]
    if not min_scales:
        return [1.0]
    return _linear_scales_from_1_to_min(min(min_scales), steps=8)


def _scaled_preferred_width(intent: _RadialSlotIntent, scale: float) -> float:
    """A row's width at ``scale``; rule L never narrows a typed width (``compress`` off, OV-1853)."""

    width = max(0.0, float(intent.width_px))
    if width <= LAYOUT_EPSILON:
        return 0.0
    if not _lane_width_narrows(intent):
        return width
    return max(_min_readable_preferred_width_px(intent), width * float(scale))


def _preferred_numeric_group_from(
    ordered_intents: Sequence[_RadialSlotIntent],
    start_pos: int,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
    *,
    keep_preset_lanes: bool,
) -> list[_RadialSlotIntent]:
    """Inside numeric rows on preset lanes, placed at their preset anchors or packed.

    The default stack places its numeric rows this way, so a preset stack with
    nothing typed draws the default figure (GX-19). With ``keep_preset_lanes``
    a pinned row with an Auto width joins the group it touches and keeps its
    radius (GX-17); a run without a preset-lane row is not a preferred group.
    """

    group: list[_RadialSlotIntent] = []
    for future in ordered_intents[start_pos:]:
        if future.slot_index in resolved_by_index:
            break
        if not (_preferred_inside_numeric(future) or (keep_preset_lanes and _anchored_inside(future))):
            break
        group.append(future)
    if not any(_preferred_inside_numeric(member) for member in group):
        return []
    return group


def _preferred_inside_numeric(intent: _RadialSlotIntent) -> bool:
    return (
        intent.placement_policy == "preferred"
        and intent.side == "inside"
        and intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
    )


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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
) -> tuple[CircularResolvedSlot, ...] | None:
    resolved_group = tuple(
        _measure_radial_slot(
            intent,
            anchor_offset_px=float(intent.anchor_offset_px or 0.0),
            width_px=float(intent.width_px),
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
) -> tuple[CircularResolvedSlot, ...] | None:
    working_occupied = list(occupied)
    working_outer = float(placement_window.outer_px)
    resolved_group: list[CircularResolvedSlot] = []
    for intent_index, intent in enumerate(intents):
        width_px = _scaled_preferred_width(intent, width_scale)
        place = _measure_anchored_inside if _anchored_inside(intent) else _place_inside_auto_fixed_width
        resolved = place(
            intent,
            width_px=width_px,
            compressed=width_px < float(intent.width_px) - LAYOUT_EPSILON,
            occupied=working_occupied,
            axis_radius_px=axis_radius_px,
            placement_window=replace(placement_window, outer_px=working_outer),
            feature_dict=feature_dict,
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
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


def _same_band(a: RadialBand | None, b: RadialBand | None) -> bool:
    if a is None or b is None:
        return a is b
    return abs(float(a.inner_px) - float(b.inner_px)) <= LAYOUT_EPSILON and abs(
        float(a.outer_px) - float(b.outer_px)
    ) <= LAYOUT_EPSILON


def _place_preferred_numeric_group(
    intents: Sequence[_RadialSlotIntent],
    *,
    occupied: Sequence[tuple[str, RadialBand]],
    placement_window: PlacementWindow,
    axis_radius_px: float,
    feature_dict: Mapping[str, Any] | None,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    depth_config: DepthAxisSource | None,
    under_a_track: bool,
    pack: bool,
) -> tuple[CircularResolvedSlot, ...]:
    """Place a group of inside rows on preset lanes (rule L, PD-OI-078).

    Without an inside stack track above the group (``under_a_track``), its lanes start right under the tick labels:
    they move outward together, keeping their spacing, until the outermost lane row's band meets the rows above it,
    as the Auto packing places it. A pinned row keeps its radius. Then the rows stay on their lanes when the lanes fit; otherwise they pack outward
    against each other at their own widths, and then at narrower Auto widths (``_preferred_group_width_scales``).
    With ``pack`` (the "limit" pass) the rows pack outward at their own widths without trying the lanes.
    """

    def at_anchors(group: Sequence[_RadialSlotIntent]) -> tuple[CircularResolvedSlot, ...] | None:
        return _try_place_preferred_numeric_group_at_anchors(
            group,
            occupied=occupied,
            placement_window=placement_window,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
            depth_config=depth_config,
        )

    def packed_at(group: Sequence[_RadialSlotIntent], scale: float) -> tuple[CircularResolvedSlot, ...] | None:
        return _place_inside_auto_group_with_width_scale(
            group,
            width_scale=scale,
            occupied=occupied,
            placement_window=placement_window,
            axis_radius_px=axis_radius_px,
            feature_dict=feature_dict,
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
            depth_config=depth_config,
        )

    lanes = list(intents)
    # The first row on a lane, under the pinned rows of the group above it (GX-17): its packed place against them
    # sets the shift of every lane.
    first = next((pos for pos, intent in enumerate(intents) if _preferred_inside_numeric(intent)), None)
    if first is not None and not under_a_track:
        first_lane = at_anchors(intents[: first + 1])
        first_packed = packed_at(intents[: first + 1], 1.0)
        lane_band = first_lane[first].draw_band_px if first_lane else None
        packed_band = first_packed[first].draw_band_px if first_packed else None
        shift = 0.0 if lane_band is None or packed_band is None else float(packed_band.outer_px) - float(lane_band.outer_px)
        if shift > LAYOUT_EPSILON:
            lanes = [
                replace(intent, anchor_offset_px=float(intent.anchor_offset_px or 0.0) + shift)
                if _preferred_inside_numeric(intent)
                else intent
                for intent in intents
            ]
    on_lanes = at_anchors(lanes)
    if on_lanes is not None and not pack:
        return on_lanes
    if pack:
        packed = packed_at(intents, 1.0)
        # Where the packed places are the lanes, keep the lanes' measurement: the packed radii differ from the
        # lanes only in the last bits, which would change the SVG IDs hashed from them (OV-1852).
        if packed is not None and on_lanes is not None and all(
            _same_band(a.draw_band_px, b.draw_band_px) for a, b in zip(packed, on_lanes)
        ):
            return on_lanes
        if packed is not None:
            return packed

    for scale in _preferred_group_width_scales(intents):
        placed = packed_at(intents, scale)
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthAxisSource | None,
) -> float:
    width = max(0.0, float(intent.width_px))
    if intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        if intent.placement_policy == "preferred":
            if _lane_width_narrows(intent):
                width = _min_readable_preferred_width_px(intent)
        elif intent.compress:
            width = _min_readable_numeric_width_px(intent)

    if intent.renderer in {"features", "ticks", "depth"}:
        resolved = _measure_radial_slot(
            intent,
            anchor_offset_px=0.0,
            width_px=width,
            axis_radius_px=float(axis_radius_px),
            feature_dict=feature_dict,
            geometry=geometry,
            cfg=cfg,
            total_length=int(total_length),
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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthAxisSource | None,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
) -> tuple[float, tuple[_KeptInsideSlot, ...]]:
    """The smallest span the unresolved inside rows after ``start_pos`` need, and those rows."""

    span = 0.0
    slots: list[_KeptInsideSlot] = []
    outer_neighbor = current_outer_intent
    for future in ordered_intents[start_pos + 1:]:
        if future.slot_index in resolved_by_index:
            break
        if future.side != "inside" or future.placement_policy == "hard":
            break
        slots.append(_KeptInsideSlot(future.slot_id, future.renderer))
        span += _gap_between_inner_outer_tracks(future, outer_neighbor) + _minimum_future_inside_width_px(
            future,
            axis_radius_px=axis_radius_px,
            geometry=geometry,
            cfg=cfg,
            total_length=total_length,
            feature_dict=feature_dict,
            depth_config=depth_config,
        )
        outer_neighbor = future
    return span, tuple(slots)


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
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    total_length: int,
    feature_dict: Mapping[str, Any] | None,
    depth_config: DepthAxisSource | None,
    resolved_by_index: Mapping[int, CircularResolvedSlot],
    group: Sequence[_RadialSlotIntent] = (),
) -> PlacementWindow:
    inner_limit = 0.0
    reserved_inner_limit: float | None = None
    center_reserved = max(
        (float(band.outer_px) for owner, band in occupied if owner == _CENTER_RESERVATION), default=0.0
    )
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
    future_span, future_slots = _future_unresolved_inside_span_px(
        ordered_intents,
        start_pos=start_pos,
        current_outer_intent=current_outer_intent,
        axis_radius_px=axis_radius_px,
        geometry=geometry,
        cfg=cfg,
        total_length=total_length,
        feature_dict=feature_dict,
        depth_config=depth_config,
        resolved_by_index=resolved_by_index,
    )
    if future_span <= LAYOUT_EPSILON:
        return PlacementWindow(
            inner_limit, float(inside_max_outer), reserved_inner_limit, center_reserved_px=center_reserved
        )
    # The span kept free for the rows after this group bounds both bands.
    inner_limit = max(
        inner_limit,
        _inner_limit_with_reserved_future_span(
            occupied,
            inner_limit_px=inner_limit,
            outer_limit_px=float(inside_max_outer),
            required_span_px=future_span,
        ),
    )
    return PlacementWindow(
        inner_limit,
        float(inside_max_outer),
        kept_inside_px=float(future_span),
        kept_inside_slots=future_slots,
        center_reserved_px=center_reserved,
    )


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
    geometry: CircularCanvasGeometry
    cfg: GbdrawConfig
    slots: Sequence[CircularTrackSlot]
    feature_dict: Mapping[str, Any] | None
    preferred_anchor_slot_ids: Collection[str]
    depth_config: DepthAxisStyle | None
    center_reserved_radius_explicit: bool
    keep_preset_lanes: bool
    annotation_label_rooms: Mapping[str, AnnotationLabelRoom]
    definition_pass: _DefinitionPass


def resolve_circular_radial_layout(
    *,
    total_length: int,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    slots: Sequence[CircularTrackSlot],
    feature_dict: Mapping[str, Any] | None = None,
    definition_reserved_radius_px: float | None = None,
    preferred_anchor_slot_ids: Collection[str] = (),
    depth_config: DepthAxisStyle | None = None,
    center_reserved_radius_explicit: bool = False,
    annotation_label_rooms: Mapping[str, AnnotationLabelRoom] | None = None,
) -> CircularRadialLayout:
    """Resolve every slot band; a failed inside placement names its cause.

    ``annotation_label_rooms`` measure, per annotation slot ID, the room its
    labels need beyond its band; an annotation slot reserves that room, and an
    annotation overlay's room is reserved past the slot that holds it
    (``_span_room_holder``) while that slot is packed, so rings placed around them keep clear of the
    letters. Two pinned slots keep their radii even when one reaches into the
    other's label room, and when the rings on one side of the axis cannot
    make the room at all, they keep the places they take without that side's
    rooms.

    ``definition_reserved_radius_px`` is the centre's band: with
    ``center_reserved_radius_explicit`` the tracks keep out of it and the
    inside tracks narrow as far as their readable minimums (DF-D05);
    otherwise it is the centre definition text's band, placed by rule L
    (``_resolve_around_definition``). Radial labels keep out of it either way.
    When an inside slot cannot be placed but the same slots fit without an
    explicit center reservation, the ``TRACK_LAYOUT`` diagnostic reason is
    ``CENTER_RESERVED``; otherwise it is ``CANNOT_FIT``.
    """

    layout_inputs = _RadialLayoutInputs(
        total_length=total_length,
        geometry=geometry,
        cfg=cfg,
        slots=slots,
        feature_dict=feature_dict,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
        depth_config=depth_config,
        center_reserved_radius_explicit=center_reserved_radius_explicit,
        keep_preset_lanes=True,
        annotation_label_rooms=annotation_label_rooms or {},
        definition_pass="none",
    )

    try:
        return _resolve_explaining_failure(definition_reserved_radius_px, layout_inputs.copy())
    except ValidationError as error:
        rooms = layout_inputs["annotation_label_rooms"]
        kept = _label_rooms_off_failed_side(rooms, slots, error)
        if len(kept) == len(rooms):
            raise
        # Rings that cannot make the labels' room (a pinned ring already at its smallest width) keep their places,
        # and the labels on that side may meet them: explicit radii win (OD-1192, OV-1207). The other side keeps
        # its rooms. Whether Circular should warn about the overlap is open (OV-1208).
        without_rooms = layout_inputs.copy()
        without_rooms["annotation_label_rooms"] = kept
        try:
            return _resolve_explaining_failure(definition_reserved_radius_px, without_rooms)
        except ValidationError:
            raise error from None


def _label_rooms_off_failed_side(
    rooms: Mapping[str, AnnotationLabelRoom], slots: Sequence[CircularTrackSlot], error: ValidationError
) -> dict[str, AnnotationLabelRoom]:
    """``rooms`` without those that reach the side (inside or outside of the axis) of the slot ``error`` names: an
    annotation slot's room reaches its own side; a span's room reaches the side of the slot that holds it
    (``_span_room_holder``), and the outside when that slot is split feature lanes or the inside row next to the axis
    (its room then reaches past the axis, as in the placement order). An error that names no slot drops every
    room."""

    normalized = normalize_circular_track_slots(slots)
    slot_index = (error.diagnostic or {}).get("slotIndex")
    failed = next((slot.side for slot in normalized if slot.slot_index == slot_index), None)
    by_id = {slot.id: slot for slot in normalized}
    first_inside = min(
        (slot for slot in normalized if slot.side == "inside"), key=lambda slot: slot.slot_index, default=None
    )

    def sides(slot: NormalizedCircularTrackSlot) -> set[str]:
        annotation = slot.annotation
        if slot.side != "overlay" or annotation is None:
            return {slot.side}
        holder = _span_room_holder(annotation, by_id)
        if holder is None:
            return set()
        past_axis = holder.side == "overlay" or holder is first_inside
        return {holder.side} - {"overlay"} | ({"outside"} if past_axis else set())

    return {
        slot_id: room
        for slot_id, room in rooms.items()
        if failed is not None and slot_id in by_id and failed not in sides(by_id[slot_id])
    }


def _resolve_explaining_failure(
    definition_reserved_radius_px: float | None, layout_inputs: _RadialLayoutInputs
) -> CircularRadialLayout:
    if (
        not layout_inputs["center_reserved_radius_explicit"]
        and definition_reserved_radius_px is not None
        and definition_reserved_radius_px > LAYOUT_EPSILON
    ):
        return _resolve_around_definition(float(definition_reserved_radius_px), layout_inputs)
    return _resolve_in_pass(definition_reserved_radius_px, layout_inputs, "none")


def _resolve_around_definition(definition_reserved_radius_px: float, layout_inputs: _RadialLayoutInputs) -> CircularRadialLayout:
    """Rule L (PD-OI-078): place the tracks around the centre definition text's band.

    1. "move": the band is occupied; the inside rows move outward at their own widths, then narrow, each Auto width
       to the floor at the least (``canvas.circular.track_narrowing_min_px``; no narrowing with
       ``canvas.circular.track_narrowing`` off). A typed width, a track narrower than the floor, and a pinned radius
       are kept.
    2. When that fails and the band does not reach the tracks placed without it, those tracks (they failed for
       another reason, as without a definition).
    3. Otherwise "limit": the inside rows on preset lanes pack outward at their own widths, and the definition is
       drawn over them (the label fit check warns, PD-OI-097). If even that fails, the tracks placed without the band.
    """

    try:
        return _resolve_in_pass(definition_reserved_radius_px, layout_inputs, "move")
    except ValidationError:
        pass
    free = _resolve_in_pass(definition_reserved_radius_px, layout_inputs, "none")
    band = RadialBand(0.0, definition_reserved_radius_px)
    reserved = [(slot.id, slot.reserved_band_px) for slot in free.slots if slot.reserved_band_px is not None]
    if _reserved_overlap_any(band, reserved) is None:
        return free
    try:
        return _resolve_in_pass(definition_reserved_radius_px, layout_inputs, "limit")
    except ValidationError:
        return free


def _resolve_in_pass(
    definition_reserved_radius_px: float | None, layout_inputs: _RadialLayoutInputs, definition_pass: _DefinitionPass
) -> CircularRadialLayout:
    # Each pass starts from the caller's inputs: the retry below turns the preset lanes off in its own copy.
    layout_inputs = layout_inputs.copy()
    layout_inputs["definition_pass"] = definition_pass
    try:
        return _resolve_circular_radial_layout(
            definition_reserved_radius_px=definition_reserved_radius_px,
            **layout_inputs,
        )
    except ValidationError:
        # A preset lane is a preference (GX-19): when keeping the preset-lane
        # rows at their anchors leaves the stack no room, they reflow with the
        # Auto rows around them, and a failure names that layout's cause.
        layout_inputs["keep_preset_lanes"] = False
    try:
        return _resolve_circular_radial_layout(
            definition_reserved_radius_px=definition_reserved_radius_px,
            **layout_inputs,
        )
    except _InsideFitError as error:
        if (
            not layout_inputs["center_reserved_radius_explicit"]
            or definition_reserved_radius_px is None
            or definition_reserved_radius_px <= LAYOUT_EPSILON
        ):
            raise
        try:
            _resolve_circular_radial_layout(definition_reserved_radius_px=None, **layout_inputs)
        except ValidationError:
            raise error from None
        raise error.limited_by_center(float(definition_reserved_radius_px)) from None


def _resolve_circular_radial_layout(
    *,
    total_length: int,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    slots: Sequence[CircularTrackSlot],
    feature_dict: Mapping[str, Any] | None = None,
    definition_reserved_radius_px: float | None = None,
    preferred_anchor_slot_ids: Collection[str] = (),
    depth_config: DepthAxisSource | None = None,
    center_reserved_radius_explicit: bool = False,
    keep_preset_lanes: bool = True,
    annotation_label_rooms: Mapping[str, AnnotationLabelRoom],
    definition_pass: _DefinitionPass = "none",
) -> CircularRadialLayout:
    axis_radius_px = float(geometry.radius)
    axis = CircularAxisLayout(
        radius_px=axis_radius_px,
        stroke_width_px=float(cfg.objects.axis.circular.stroke_width.for_length_param(str(geometry.length_param))),
    )
    definition_band = (
        RadialBand(0.0, float(definition_reserved_radius_px))
        if definition_reserved_radius_px is not None and definition_reserved_radius_px > LAYOUT_EPSILON
        else None
    )
    occupied: list[tuple[str, RadialBand]] = []
    # An explicit center_reserved_radius keeps the tracks out of the centre (DF-D05), and so does the centre
    # definition text's band in the "move" pass of rule L; otherwise the band is only recorded (radial labels keep
    # out of it).
    if (
        (center_reserved_radius_explicit or definition_pass == "move")
        and definition_band is not None
        and definition_band.width_px > LAYOUT_EPSILON
    ):
        occupied.append((_CENTER_RESERVATION, definition_band))

    intents = _slot_intents(
        slots,
        geometry=geometry,
        cfg=cfg,
        preferred_anchor_slot_ids=preferred_anchor_slot_ids,
        axis=axis,
        annotation_label_rooms=annotation_label_rooms,
    )
    if definition_pass == "move":
        circular = cfg.canvas.circular
        floor_px = float(circular.track_narrowing_min_px) if circular.track_narrowing else math.inf
        intents = [replace(intent, narrowing_floor_px=floor_px) for intent in intents]
    elif center_reserved_radius_explicit:
        intents = [replace(intent, typed_lane_width_narrows=True) for intent in intents]
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
    # Pinned slots are checked against each other without annotation or span label rooms: explicit radii win.
    pinned_occupied = list(occupied)
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
            geometry=geometry,
            cfg=cfg,
            total_length=int(total_length),
            depth_config=depth_config,
        )
        if _anchored_inside(intent):
            if resolved.reserved_band_px is not None:
                anchored_blockers.append((intent.slot_id, resolved.reserved_band_px))
            continue
        pinned_band = resolved.packing_band_px if intent.label_room is not None else resolved.own_reserved_band_px
        if _slot_reserves(intent) and pinned_band is not None:
            conflict = _reserved_overlap_any(pinned_band, pinned_occupied)
            if conflict is not None:
                diagnostic = _cannot_fit_diagnostic(intent)
                if conflict[0] == _CENTER_RESERVATION:
                    diagnostic = _center_reserved_diagnostic(diagnostic)
                raise ValidationError(
                    f"Pinned circular track slot '{intent.slot_id}' overlaps reserved circular slot '{conflict[0]}'.",
                    diagnostic=diagnostic,
                )
        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and (band := _occupied_band(intent, resolved)) is not None:
            occupied.append((intent.slot_id, band))
            if pinned_band is not None:
                pinned_occupied.append((intent.slot_id, pinned_band))

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

    def place_outside_rows() -> None:
        """Outside rows are placed from the axis outward: each row's window starts at the row below it and ends at
        the next pinned row above it."""

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
                    geometry=geometry,
                    cfg=cfg,
                    total_length=int(total_length),
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

    # The inside row next to the axis may carry a span's label room past the axis; then the inside rows are placed
    # first and the outside rows keep out of that room. Otherwise the outside rows go first, as their reserved
    # bands (outside tick labels) can reach inside.
    first_inside = next((intent for intent in ordered_intents if intent.side == "inside"), None)
    inside_rows_first = first_inside is not None and bool(first_inside.span_label_rooms)
    if not inside_rows_first:
        place_outside_rows()

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
        # Outside rows are placed by ``place_outside_rows``, and every overlay (placement policy "overlay") was
        # resolved with the pinned slots above, so only inside rows remain.
        if intent.side == "outside":
            continue

        # Every group starts with this row, which faces the rows above it.
        outer_limit = _outer_limit_below_rows(inside_max_outer, rows_above_inside, intent)
        movable_group = _inside_movable_stack_group_from(
            ordered_intents,
            intent_pos,
            resolved_by_index,
            keep_preset_lanes=keep_preset_lanes,
        )
        if movable_group and (len(movable_group) > 1 or _anchored_inside(intent)):
            placement_window = _inside_placement_window(
                ordered_intents,
                start_pos=intent_pos + len(movable_group) - 1,
                current_outer_intent=movable_group[-1],
                inside_max_outer=outer_limit,
                occupied=occupied,
                axis_radius_px=axis_radius_px,
                geometry=geometry,
                cfg=cfg,
                total_length=int(total_length),
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
                geometry=geometry,
                cfg=cfg,
                total_length=int(total_length),
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
            keep_preset_lanes=keep_preset_lanes,
        )
        if preferred_group:
            placement_window = _inside_placement_window(
                ordered_intents,
                start_pos=intent_pos + len(preferred_group) - 1,
                current_outer_intent=preferred_group[-1],
                inside_max_outer=outer_limit,
                occupied=occupied,
                axis_radius_px=axis_radius_px,
                geometry=geometry,
                cfg=cfg,
                total_length=int(total_length),
                feature_dict=feature_dict,
                depth_config=depth_config,
                resolved_by_index=resolved_by_index,
                group=preferred_group,
            )
            resolved_group = _place_preferred_numeric_group(
                preferred_group,
                occupied=occupied,
                placement_window=placement_window,
                axis_radius_px=axis_radius_px,
                feature_dict=feature_dict,
                geometry=geometry,
                cfg=cfg,
                total_length=int(total_length),
                depth_config=depth_config,
                under_a_track=any(
                    row.renderer not in _ROWS_ABOVE_THE_INSIDE_STACK for row, _edge in rows_above_inside
                ),
                pack=definition_pass == "limit",
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
            geometry=geometry,
            cfg=cfg,
            total_length=int(total_length),
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
                    geometry=geometry,
                    cfg=cfg,
                    total_length=int(total_length),
                    depth_config=depth_config,
                )
            except ValidationError:
                fallback_group = _inside_movable_stack_group_from(
                    ordered_intents,
                    intent_pos,
                    resolved_by_index,
                    keep_preset_lanes=keep_preset_lanes,
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
                    geometry=geometry,
                    cfg=cfg,
                    total_length=int(total_length),
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
                    geometry=geometry,
                    cfg=cfg,
                    total_length=int(total_length),
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
            geometry=geometry,
            cfg=cfg,
            total_length=int(total_length),
            depth_config=depth_config,
        )

        resolved_by_index[intent.slot_index] = resolved
        if _slot_reserves(intent) and (band := _occupied_band(intent, resolved)) is not None:
            occupied.append((intent.slot_id, band))
        if resolved.packing_band_px is not None and resolved.side == "inside":
            rows_above_inside.append((intent, float(resolved.packing_band_px.inner_px)))

    if inside_rows_first:
        place_outside_rows()

    by_id = {slot.id: slot for slot in resolved_by_index.values()}
    for intent in intents:
        if intent.renderer != "annotations" or intent.side != "overlay":
            continue
        overlay_resolved = resolved_by_index.get(intent.slot_index)
        annotation = intent.slot.annotation
        if overlay_resolved is None or annotation is None:
            continue
        referenced = [
            by_id[slot_id]
            for slot_id in (annotation.anchor_slot, annotation.to_slot)
            if slot_id in by_id
        ]
        if not referenced:
            continue
        anchor = referenced[0]
        anchor_band = anchor.draw_band_px or anchor.packing_band_px
        # to_slot: the hull of the two slots' whole bands, every feature lane included.
        # cover_anchor alone keeps the draw band (the first feature lane of each strand).
        whole_bands = [slot.packing_band_px or slot.draw_band_px for slot in referenced]
        span_band = (
            band_union([item for item in whole_bands if item is not None])
            if annotation.to_slot is not None
            else None
        )
        if span_band is not None:
            anchor_radius, width, draw_band = span_band.center_px, span_band.width_px, span_band
        else:
            anchor_radius = (
                float(anchor_band.center_px)
                if anchor_band is not None
                else float(anchor.anchor_radius_px or axis_radius_px)
            )
            # The overlay's own width is the intent's, the one its labels' room is measured from (``_SpanLabelRoom``).
            width = (
                float(anchor_band.width_px)
                if annotation.cover_anchor and anchor_band is not None
                else float(intent.width_px)
            )
            # The same band its labels' room is measured from when the anchor holds that room.
            anchor_draw_band = _span_room_band(anchor, _span_extent(annotation), width)
            draw_band = (
                anchor_draw_band if anchor_draw_band is not None else _band_from_center_width(anchor_radius, width)
            )
        # The room lies past the holder, so it keeps the labels off other tracks only when the holder sets the span's
        # outer edge.
        holder = _span_room_holder(annotation, by_id)
        holder_band = None if holder is None else _span_room_band(holder, _span_extent(annotation), width)
        anchored = replace(
            overlay_resolved,
            anchor_radius_px=anchor_radius,
            anchor_offset_px=anchor_radius - axis_radius_px,
            resolved_width_px=width,
            packing_band_px=None,
            draw_band_px=draw_band,
            reserved_band_px=None,
            spans_slots=span_band is not None or (annotation.cover_anchor and anchor_band is not None),
            span_keeps_label_room=holder_band is not None
            and float(holder_band.outer_px) >= float(draw_band.outer_px) - LAYOUT_EPSILON,
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
    "track_layout_reason",
]
