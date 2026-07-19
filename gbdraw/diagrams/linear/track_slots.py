from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from ...canvas import LinearCanvasConfigurator
from ...config.models import GbdrawConfig
from ...layout.linear import VerticalBand, union_vertical_bands
from ...tracks import NormalizedLinearTrackSlot


@dataclass(frozen=True)
class LinearResolvedTrack:
    slot_index: int
    id: str
    renderer: str
    side: str
    y_offset: float
    height: float
    spacing_after_px: float
    top_extent: float
    bottom_extent: float
    z: int
    params: Mapping[str, object]


@dataclass(frozen=True)
class LinearTrackLayout:
    slots: tuple[LinearResolvedTrack, ...]
    top_extent: float
    bottom_extent: float
    plot_tracks_height: float
    plot_tracks_visual_bottom: float
    depth_track_offsets: tuple[float, ...]
    depth_track_heights: tuple[float, ...]
    gc_content_track_offset: float
    gc_skew_track_offset: float


@dataclass(frozen=True)
class LinearSlotFootprint:
    """Record-local paint and reservation geometry for one slot."""

    paint_band: VerticalBand | None
    reserve_band: VerticalBand
    data_available: bool = True


@dataclass(frozen=True)
class LinearResolvedSlot:
    """One slot after per-record vertical packing."""

    slot_index: int
    id: str
    renderer: str
    side: str
    origin_y: float
    height: float
    spacing_after_px: float
    paint_band: VerticalBand | None
    reserve_band: VerticalBand
    data_available: bool
    z: int
    params: Mapping[str, object]

@dataclass(frozen=True)
class LinearRecordVerticalPlan:
    """Immutable vertical geometry for one Linear record."""

    axis_band: VerticalBand
    slots: tuple[LinearResolvedSlot, ...]
    record_body_band: VerticalBand
    comparison_exclusion_band: VerticalBand
    canvas_band: VerticalBand

    @property
    def top_extent(self) -> float:
        return max(0.0, -self.record_body_band.top_y)

    @property
    def bottom_extent(self) -> float:
        return max(0.0, self.record_body_band.bottom_y)

    @property
    def comparison_top_extent(self) -> float:
        return max(0.0, -self.comparison_exclusion_band.top_y)

    @property
    def comparison_bottom_extent(self) -> float:
        return max(0.0, self.comparison_exclusion_band.bottom_y)

    @property
    def canvas_top_extent(self) -> float:
        return max(0.0, -self.canvas_band.top_y)

    @property
    def canvas_bottom_extent(self) -> float:
        return max(0.0, self.canvas_band.bottom_y)

    def slot_by_id(self, slot_id: str) -> LinearResolvedSlot:
        for slot in self.slots:
            if slot.id == slot_id:
                return slot
        raise KeyError(slot_id)


def _resolve_scalar_px(value, default: float) -> float:
    if value is None:
        return float(default)
    return float(value.resolve(1.0))


def _slot_height(
    slot: NormalizedLinearTrackSlot,
    *,
    canvas_config: LinearCanvasConfigurator,
) -> float:
    if slot.renderer == "features":
        return max(0.0, _resolve_scalar_px(slot.height, 0.0))
    if slot.renderer == "spacer":
        return max(0.0, _resolve_scalar_px(slot.height, 0.0))
    if slot.renderer == "depth":
        return max(0.0, _resolve_scalar_px(slot.height, canvas_config.default_depth_height))
    return max(0.0, _resolve_scalar_px(slot.height, canvas_config.default_gc_height))


def _slot_spacing(
    slot: NormalizedLinearTrackSlot,
    *,
    canvas_config: LinearCanvasConfigurator,
) -> float:
    if slot.spacing is not None:
        return max(0.0, _resolve_scalar_px(slot.spacing, 0.0))
    if slot.renderer == "depth":
        return max(0.0, float(canvas_config.configured_depth_padding))
    if slot.renderer == "features":
        return max(0.0, float(canvas_config.vertical_padding))
    return 0.0


def _axis_clearance_px(canvas_config: LinearCanvasConfigurator) -> float:
    return max(0.0, float(canvas_config.vertical_padding))


def _slot_needs_axis_clearance(slot: NormalizedLinearTrackSlot) -> bool:
    return slot.renderer != "features"


def _track_extents_for_slot(
    slot: NormalizedLinearTrackSlot,
    *,
    height: float,
    cfg: GbdrawConfig,
) -> tuple[float, float]:
    if slot.renderer == "dinucleotide_content":
        mode = str(getattr(cfg.objects.gc_content, "mode", "deviation")).strip().lower()
        if mode == "percent":
            return 0.0, float(height)
        return 0.5 * float(height), 0.5 * float(height)
    if slot.renderer == "dinucleotide_skew":
        return 0.5 * float(height), 0.5 * float(height)
    if slot.renderer in {"depth", "annotations", "spacer"}:
        return 0.0, float(height)
    return 0.0, 0.0


def _resolve_below_tracks(
    slots: list[NormalizedLinearTrackSlot],
    *,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[list[LinearResolvedTrack], float]:
    resolved: list[LinearResolvedTrack] = []
    cursor = (
        _axis_clearance_px(canvas_config)
        if slots and _slot_needs_axis_clearance(slots[0])
        else 0.0
    )
    for slot in slots:
        height = _slot_height(slot, canvas_config=canvas_config)
        spacing = _slot_spacing(slot, canvas_config=canvas_config)
        top_extent, bottom_extent = _track_extents_for_slot(slot, height=height, cfg=cfg)
        y_offset = cursor + top_extent
        resolved.append(
            LinearResolvedTrack(
                slot_index=slot.slot_index,
                id=slot.id,
                renderer=slot.renderer,
                side=slot.side,
                y_offset=float(y_offset),
                height=float(height),
                spacing_after_px=float(spacing),
                top_extent=float(top_extent),
                bottom_extent=float(bottom_extent),
                z=slot.z,
                params=dict(slot.params),
            )
        )
        if slot.renderer != "features":
            cursor = y_offset + bottom_extent + spacing
    return resolved, float(cursor)


def _resolve_above_tracks(
    slots: list[NormalizedLinearTrackSlot],
    *,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> tuple[list[LinearResolvedTrack], float]:
    by_slot_index: dict[int, LinearResolvedTrack] = {}
    axis_adjacent_slot = slots[-1] if slots else None
    cursor = (
        _axis_clearance_px(canvas_config)
        if axis_adjacent_slot is not None and _slot_needs_axis_clearance(axis_adjacent_slot)
        else 0.0
    )
    for slot in reversed(slots):
        height = _slot_height(slot, canvas_config=canvas_config)
        spacing = _slot_spacing(slot, canvas_config=canvas_config)
        top_extent, bottom_extent = _track_extents_for_slot(slot, height=height, cfg=cfg)
        y_offset = -(cursor + bottom_extent)
        resolved = LinearResolvedTrack(
            slot_index=slot.slot_index,
            id=slot.id,
            renderer=slot.renderer,
            side=slot.side,
            y_offset=float(y_offset),
            height=float(height),
            spacing_after_px=float(spacing),
            top_extent=float(top_extent),
            bottom_extent=float(bottom_extent),
            z=slot.z,
            params=dict(slot.params),
        )
        by_slot_index[slot.slot_index] = resolved
        if slot.renderer != "features":
            cursor = cursor + top_extent + bottom_extent + spacing
    return [by_slot_index[slot.slot_index] for slot in slots], float(cursor)


def resolve_linear_track_layout(
    slots: list[NormalizedLinearTrackSlot],
    *,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> LinearTrackLayout:
    """Resolve normalized linear slots into concrete y offsets around the axis."""

    above_slots = [slot for slot in slots if slot.side == "above"]
    below_slots = [slot for slot in slots if slot.side == "below"]
    overlay_slots = [slot for slot in slots if slot.side == "overlay"]

    above_resolved, top_extent = _resolve_above_tracks(
        above_slots,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    below_resolved, bottom_extent = _resolve_below_tracks(
        below_slots,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    overlay_resolved: list[LinearResolvedTrack] = []
    anchor_by_id = {track.id: track for track in (*above_resolved, *below_resolved)}
    for slot in overlay_slots:
        height = _slot_height(slot, canvas_config=canvas_config) if slot.renderer == "annotations" else 0.0
        spacing = (
            _slot_spacing(slot, canvas_config=canvas_config)
            if slot.renderer == "features"
            else 0.0
        )
        anchor = anchor_by_id.get(str(slot.params.get("anchor_slot", "")))
        anchor_center = 0.0
        if anchor is not None:
            anchor_center = 0.5 * (
                (float(anchor.y_offset) - float(anchor.top_extent))
                + (float(anchor.y_offset) + float(anchor.bottom_extent))
            )
        resolved_overlay = LinearResolvedTrack(
            slot_index=slot.slot_index,
            id=slot.id,
            renderer=slot.renderer,
            side=slot.side,
            y_offset=anchor_center - (0.5 * height),
            height=float(height),
            spacing_after_px=float(spacing),
            top_extent=0.0,
            bottom_extent=float(height),
            z=slot.z,
            params=dict(slot.params),
        )
        overlay_resolved.append(resolved_overlay)
        anchor_by_id[slot.id] = resolved_overlay

    resolved_by_index = {
        track.slot_index: track
        for track in (*above_resolved, *below_resolved, *overlay_resolved)
    }
    resolved = tuple(
        resolved_by_index[slot.slot_index]
        for slot in sorted(slots, key=lambda item: item.slot_index)
        if slot.slot_index in resolved_by_index
    )

    depth_offsets_by_index: dict[int, float] = {}
    depth_heights_by_index: dict[int, float] = {}
    gc_content_track_offset = 0.0
    gc_skew_track_offset = 0.0
    for track in resolved:
        if track.renderer == "depth":
            track_index = int(track.params.get("track_index", 0))
            depth_offsets_by_index[track_index] = float(track.y_offset)
            depth_heights_by_index[track_index] = float(track.height)
        elif track.renderer == "dinucleotide_content":
            gc_content_track_offset = float(track.y_offset)
        elif track.renderer == "dinucleotide_skew":
            gc_skew_track_offset = float(track.y_offset)

    depth_track_count = max(
        max(depth_offsets_by_index.keys(), default=-1),
        max(depth_heights_by_index.keys(), default=-1),
    ) + 1
    depth_track_offsets = tuple(
        depth_offsets_by_index.get(index, 0.0)
        for index in range(depth_track_count)
    )
    depth_track_heights = tuple(
        depth_heights_by_index.get(index, 0.0)
        for index in range(depth_track_count)
    )

    return LinearTrackLayout(
        slots=resolved,
        top_extent=float(top_extent),
        bottom_extent=float(bottom_extent),
        plot_tracks_height=float(bottom_extent),
        plot_tracks_visual_bottom=float(bottom_extent),
        depth_track_offsets=depth_track_offsets,
        depth_track_heights=depth_track_heights,
        gc_content_track_offset=float(gc_content_track_offset),
        gc_skew_track_offset=float(gc_skew_track_offset),
    )


def _default_slot_footprint(track: LinearResolvedTrack) -> LinearSlotFootprint:
    reserve_band = VerticalBand(-float(track.top_extent), float(track.bottom_extent))
    return LinearSlotFootprint(
        paint_band=None if track.renderer == "spacer" else reserve_band,
        reserve_band=reserve_band,
        data_available=track.renderer != "spacer",
    )


def _translate_optional_band(band: VerticalBand | None, origin_y: float) -> VerticalBand | None:
    return None if band is None else band.translate(origin_y)


def _combine_slot_footprint(
    footprint: LinearSlotFootprint,
    overlays: Sequence[tuple[float, LinearSlotFootprint]],
) -> LinearSlotFootprint:
    reserve_bands = [footprint.reserve_band]
    paint_bands = [footprint.paint_band] if footprint.paint_band is not None else []
    for relative_origin, overlay in overlays:
        reserve_bands.append(overlay.reserve_band.translate(relative_origin))
        if overlay.data_available and overlay.paint_band is not None:
            paint_bands.append(overlay.paint_band.translate(relative_origin))
    return LinearSlotFootprint(
        paint_band=(union_vertical_bands(paint_bands) if paint_bands else None),
        reserve_band=union_vertical_bands(reserve_bands),
        data_available=footprint.data_available,
    )


def _resolved_slot(
    track: LinearResolvedTrack,
    footprint: LinearSlotFootprint,
    origin_y: float,
) -> LinearResolvedSlot:
    return LinearResolvedSlot(
        slot_index=track.slot_index,
        id=track.id,
        renderer=track.renderer,
        side=track.side,
        origin_y=float(origin_y),
        height=float(track.height),
        spacing_after_px=float(track.spacing_after_px),
        paint_band=(
            _translate_optional_band(footprint.paint_band, origin_y)
            if footprint.data_available
            else None
        ),
        reserve_band=footprint.reserve_band.translate(origin_y),
        data_available=bool(footprint.data_available),
        z=track.z,
        params=dict(track.params),
    )


def _pack_side_slots(
    tracks: Sequence[LinearResolvedTrack],
    *,
    direction: str,
    seed_band: VerticalBand,
    axis_band: VerticalBand,
    footprints: Mapping[str, LinearSlotFootprint],
    overlays_by_anchor: Mapping[str, Sequence[tuple[float, LinearSlotFootprint]]],
    initial_spacing: float = 0.0,
) -> list[LinearResolvedSlot]:
    if direction not in {"above", "below"}:
        raise ValueError("direction must be 'above' or 'below'")

    ordered = list(reversed(tracks)) if direction == "above" else list(tracks)
    resolved: list[LinearResolvedSlot] = []
    occupied_edge = seed_band.top_y if direction == "above" else seed_band.bottom_y
    preferred_inner_edge = axis_band.top_y if direction == "above" else axis_band.bottom_y
    spacing_from_inner = max(0.0, float(initial_spacing))

    for track in ordered:
        base_footprint = footprints.get(track.id, _default_slot_footprint(track))
        footprint = _combine_slot_footprint(
            base_footprint,
            overlays_by_anchor.get(track.id, ()),
        )
        preferred_band = footprint.reserve_band.translate(track.y_offset)
        if direction == "above":
            preferred_gap = max(0.0, preferred_inner_edge - preferred_band.bottom_y)
            gap = max(spacing_from_inner, preferred_gap)
            origin_y = min(
                float(track.y_offset),
                occupied_edge - gap - footprint.reserve_band.bottom_y,
            )
        else:
            preferred_gap = max(0.0, preferred_band.top_y - preferred_inner_edge)
            gap = max(spacing_from_inner, preferred_gap)
            origin_y = max(
                float(track.y_offset),
                occupied_edge + gap - footprint.reserve_band.top_y,
            )

        item = _resolved_slot(track, footprint, origin_y)
        resolved.append(item)
        if direction == "above":
            occupied_edge = item.reserve_band.top_y
            preferred_inner_edge = preferred_band.top_y
        else:
            occupied_edge = item.reserve_band.bottom_y
            preferred_inner_edge = preferred_band.bottom_y
        spacing_from_inner = max(0.0, float(track.spacing_after_px))

    if direction == "above":
        resolved.reverse()
    return resolved


def resolve_linear_record_vertical_plan(
    layout: LinearTrackLayout,
    *,
    axis_band: VerticalBand,
    footprints: Mapping[str, LinearSlotFootprint] | None = None,
    comparison_gap: float = 0.0,
) -> LinearRecordVerticalPlan:
    """Pack one record's slots around its axis from measured footprints."""

    footprint_map = dict(footprints or {})
    tracks_by_id = {track.id: track for track in layout.slots}
    overlay_tracks = [track for track in layout.slots if track.side == "overlay"]
    overlays_by_anchor: dict[str, list[tuple[float, LinearSlotFootprint]]] = {}
    standalone_overlays: list[LinearResolvedSlot] = []
    relative_overlay_origins: dict[str, tuple[str, float]] = {}

    for overlay in overlay_tracks:
        if overlay.renderer == "features":
            continue
        anchor_id = str(overlay.params.get("anchor_slot", "")).strip()
        anchor = tracks_by_id.get(anchor_id)
        overlay_footprint = footprint_map.get(overlay.id, _default_slot_footprint(overlay))
        if anchor is None:
            standalone_overlays.append(
                _resolved_slot(overlay, overlay_footprint, overlay.y_offset)
            )
            continue
        relative_origin = float(overlay.y_offset) - float(anchor.y_offset)
        overlays_by_anchor.setdefault(anchor.id, []).append((relative_origin, overlay_footprint))
        relative_overlay_origins[overlay.id] = (anchor.id, relative_origin)

    structural_slots: list[LinearResolvedSlot] = []
    seed_band = axis_band
    for track in overlay_tracks:
        if track.renderer != "features":
            continue
        footprint = _combine_slot_footprint(
            footprint_map.get(track.id, _default_slot_footprint(track)),
            overlays_by_anchor.get(track.id, ()),
        )
        item = _resolved_slot(track, footprint, 0.0)
        structural_slots.append(item)
        seed_band = seed_band.union(item.reserve_band)

    above_tracks = [track for track in layout.slots if track.side == "above"]
    below_tracks = [track for track in layout.slots if track.side == "below"]
    structural_spacing = max(
        (slot.spacing_after_px for slot in structural_slots),
        default=0.0,
    )
    above_slots = _pack_side_slots(
        above_tracks,
        direction="above",
        seed_band=seed_band,
        axis_band=axis_band,
        footprints=footprint_map,
        overlays_by_anchor=overlays_by_anchor,
        initial_spacing=structural_spacing,
    )
    below_slots = _pack_side_slots(
        below_tracks,
        direction="below",
        seed_band=seed_band,
        axis_band=axis_band,
        footprints=footprint_map,
        overlays_by_anchor=overlays_by_anchor,
        initial_spacing=structural_spacing,
    )

    anchors = {
        slot.id: slot
        for slot in (*structural_slots, *above_slots, *below_slots, *standalone_overlays)
    }
    resolved_overlays = list(standalone_overlays)
    for overlay in overlay_tracks:
        relation = relative_overlay_origins.get(overlay.id)
        if relation is None:
            continue
        anchor_id, relative_origin = relation
        anchor = anchors.get(anchor_id)
        if anchor is None:
            continue
        footprint = footprint_map.get(overlay.id, _default_slot_footprint(overlay))
        item = _resolved_slot(overlay, footprint, anchor.origin_y + relative_origin)
        resolved_overlays.append(item)
        anchors[overlay.id] = item

    resolved_by_index = {
        slot.slot_index: slot
        for slot in (*structural_slots, *above_slots, *below_slots, *resolved_overlays)
    }
    resolved_slots = tuple(
        resolved_by_index[track.slot_index]
        for track in sorted(layout.slots, key=lambda item: item.slot_index)
        if track.slot_index in resolved_by_index
    )
    reserve_bands = [axis_band, *(slot.reserve_band for slot in resolved_slots)]
    body_band = union_vertical_bands(reserve_bands)
    gap = max(0.0, float(comparison_gap))
    comparison_band = body_band.expand(gap)
    return LinearRecordVerticalPlan(
        axis_band=axis_band,
        slots=resolved_slots,
        record_body_band=body_band,
        comparison_exclusion_band=comparison_band,
        canvas_band=body_band,
    )


__all__ = [
    "LinearRecordVerticalPlan",
    "LinearResolvedTrack",
    "LinearResolvedSlot",
    "LinearSlotFootprint",
    "LinearTrackLayout",
    "resolve_linear_track_layout",
    "resolve_linear_record_vertical_plan",
]
