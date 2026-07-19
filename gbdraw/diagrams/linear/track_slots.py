from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from ...canvas import LinearCanvasConfigurator
from ...config.models import GbdrawConfig
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
        return 0.0
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
            spacing_after_px=0.0,
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
            try:
                track_index = int(track.params.get("track_index", 0))
            except (TypeError, ValueError):
                track_index = 0
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


__all__ = [
    "LinearResolvedTrack",
    "LinearTrackLayout",
    "resolve_linear_track_layout",
]
