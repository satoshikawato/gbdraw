from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from ...tracks.circular import (
    SUPPORTED_CIRCULAR_TRACK_RENDERERS,
    CircularTrackSlot,
)
from ...tracks.spec import CircularTrackPlacement, TrackSpec
from ...svg.circular_ticks import (
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_radius_bounds,
)


logger = logging.getLogger(__name__)

LAYOUT_EPSILON = 1e-6
COMPRESSIBLE_NUMERIC_RENDERERS = frozenset({"dinucleotide_content", "dinucleotide_skew", "depth"})


class CircularSlotAutoPackError(ValueError):
    """Raised when an auto circular slot cannot be packed in the available radius."""


@dataclass(frozen=True)
class CircularSlotFootprint:
    anchor_radius_px: float | None
    draw_inner_px: float
    draw_outer_px: float
    hard_inner_px: float
    hard_outer_px: float
    soft_inner_px: float
    soft_outer_px: float
    explicit_anchor: bool
    explicit_annulus: bool
    explicit_width: bool

    @property
    def reserved_inner_px(self) -> float:
        return self.hard_inner_px

    @property
    def reserved_outer_px(self) -> float:
        return self.hard_outer_px


@dataclass(frozen=True)
class ResolvedCircularTrackSlot:
    id: str
    renderer: str
    anchor_radius_px: float
    center_radius_px: float
    width_px: float
    draw_inner_radius_px: float
    draw_outer_radius_px: float
    hard_inner_radius_px: float
    hard_outer_radius_px: float
    soft_inner_radius_px: float
    soft_outer_radius_px: float
    reserved_inner_radius_px: float
    reserved_outer_radius_px: float
    z: int
    params: Mapping[str, Any]

    @property
    def inner_radius_px(self) -> float:
        return self.draw_inner_radius_px

    @property
    def outer_radius_px(self) -> float:
        return self.draw_outer_radius_px


@dataclass(frozen=True)
class CircularTrackLayoutContext:
    """Data needed to measure and pack circular custom slots."""

    base_radius_px: float
    legacy_centers_px: Mapping[str, float] = field(default_factory=dict)
    legacy_widths_px: Mapping[str, float] = field(default_factory=dict)
    preferred_layouts_px: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    default_gap_px: float = 0.0
    auto_start_radius_px: float | None = None
    feature_band_offsets_px: tuple[float, float] | None = None
    tick_path_ratio_bounds: tuple[float, float] | None = None
    tick_label_radius_ratio: float | None = None
    tick_label_extent_px: float = 0.0
    tick_label_offsets_px: tuple[float, float] | None = None
    tick_labels_hard: bool = False
    tick_font_size_px: float = 10.0
    tick_width_px: float = 0.0
    tick_total_len: int | None = None
    tick_track_type: str = ""
    tick_strandedness: bool = True
    tick_font_family: str = ""
    tick_dpi: int = 72
    tick_manual_interval: int | None = None
    tick_track_channel_override: str | None = None
    reserved_bands_px: tuple[tuple[float, float], ...] = ()
    min_auto_inner_radius_px: float | None = None


def _legacy_kind_for_slot(slot: CircularTrackSlot) -> str | None:
    if slot.id == "gc_content":
        return "gc_content"
    if slot.id == "gc_skew":
        return "gc_skew"
    if slot.id in {"features", "ticks", "depth"}:
        return slot.id
    renderer = str(slot.renderer)
    if renderer == "dinucleotide_content":
        return "gc_content"
    if renderer == "dinucleotide_skew":
        return "gc_skew"
    return renderer if renderer in {"features", "ticks", "depth"} else None


def _track_specs_by_match(track_specs: Sequence[TrackSpec] | None) -> tuple[dict[str, TrackSpec], dict[str, TrackSpec]]:
    by_id: dict[str, TrackSpec] = {}
    by_kind: dict[str, TrackSpec] = {}
    for ts in track_specs or []:
        by_id.setdefault(str(ts.id), ts)
        by_kind.setdefault(str(ts.kind), ts)
    return by_id, by_kind


def _track_spec_for_slot(
    slot: CircularTrackSlot,
    track_specs: Sequence[TrackSpec] | None,
) -> TrackSpec | None:
    by_id, by_kind = _track_specs_by_match(track_specs)
    if slot.id in by_id:
        return by_id[slot.id]
    legacy_kind = _legacy_kind_for_slot(slot)
    if legacy_kind is not None:
        return by_kind.get(legacy_kind)
    return None


def _track_spec_has_explicit_circular_width(ts: TrackSpec | None) -> bool:
    if ts is None or not isinstance(ts.placement, CircularTrackPlacement):
        return False
    placement = ts.placement
    return placement.width is not None or (
        placement.inner_radius is not None and placement.outer_radius is not None
    )


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


def _legacy_lookup(
    mapping: Mapping[str, float],
    slot: CircularTrackSlot,
) -> float | None:
    if slot.id in mapping:
        return float(mapping[slot.id])
    legacy_kind = _legacy_kind_for_slot(slot)
    if legacy_kind is not None and legacy_kind in mapping:
        return float(mapping[legacy_kind])
    if str(slot.renderer) in mapping:
        return float(mapping[str(slot.renderer)])
    return None


def _preferred_layout_lookup(
    mapping: Mapping[str, tuple[float, float]],
    slot: CircularTrackSlot,
) -> tuple[float, float] | None:
    if slot.id in mapping:
        center_px, width_px = mapping[slot.id]
        return float(center_px), float(width_px)
    return None


def _annulus(center_px: float, width_px: float) -> tuple[float, float]:
    half_width = max(0.0, 0.5 * float(width_px))
    inner = float(center_px) - half_width
    outer = float(center_px) + half_width
    return (inner, outer) if inner <= outer else (outer, inner)


def _slot_gap_after_px(slot: CircularTrackSlot, context: CircularTrackLayoutContext) -> float:
    if slot.gap_after is None:
        return float(context.default_gap_px)
    return max(0.0, float(slot.gap_after.resolve(float(context.base_radius_px))))


def _is_compressible_numeric_width(
    slot: CircularTrackSlot,
    *,
    explicit_annulus: bool,
    explicit_width: bool,
) -> bool:
    return (
        str(slot.renderer) in COMPRESSIBLE_NUMERIC_RENDERERS
        and not bool(explicit_annulus)
        and not bool(explicit_width)
    )


def _uses_auto_inner_guard(slot: CircularTrackSlot) -> bool:
    return str(slot.renderer) in COMPRESSIBLE_NUMERIC_RENDERERS


def _entry_width_px(entry: Mapping[str, Any], width_scale: float) -> float:
    width = max(0.0, float(entry["width_px"]))
    if entry.get("compressible_width"):
        return max(0.0, width * float(width_scale))
    return width


def _entry_gap_after_px(
    entry: Mapping[str, Any],
    *,
    context: CircularTrackLayoutContext,
    implicit_gap_scale: float,
) -> float:
    slot = entry["slot"]
    gap = _slot_gap_after_px(slot, context)
    if entry.get("explicit_gap"):
        return gap
    return max(0.0, gap * float(implicit_gap_scale))


def _overlaps(a: tuple[float, float], b: tuple[float, float]) -> bool:
    a_inner, a_outer = sorted((float(a[0]), float(a[1])))
    b_inner, b_outer = sorted((float(b[0]), float(b[1])))
    return a_inner < (b_outer - LAYOUT_EPSILON) and a_outer > (b_inner + LAYOUT_EPSILON)


def _normalize_band(band: tuple[float, float]) -> tuple[float, float]:
    inner, outer = sorted((float(band[0]), float(band[1])))
    return max(0.0, inner), max(0.0, outer)


def _reserved_auto_inner_guard(context: CircularTrackLayoutContext) -> float:
    guard = 0.0
    if context.min_auto_inner_radius_px is not None:
        guard = max(guard, float(context.min_auto_inner_radius_px))
    for band in context.reserved_bands_px:
        inner, outer = _normalize_band(band)
        if inner <= LAYOUT_EPSILON:
            guard = max(guard, outer)
    return guard


def _auto_pack_error(slot_id: str, *, inner_guard_px: float) -> ValueError:
    if inner_guard_px > LAYOUT_EPSILON:
        return CircularSlotAutoPackError(
            f"circular track slot '{slot_id}' cannot fit without overlapping the center definition; "
            "reduce width, remove inner slots, or set an explicit radius."
        )
    return CircularSlotAutoPackError(f"circular track slot '{slot_id}' could not be packed without overlap")


def _free_intervals_px(
    occupied: Sequence[tuple[str, tuple[float, float]]],
    *,
    inner_limit_px: float,
    outer_limit_px: float,
) -> list[tuple[float, float]]:
    lower_limit = max(0.0, float(inner_limit_px))
    upper_limit = max(lower_limit, float(outer_limit_px))
    intervals: list[tuple[float, float]] = []
    cursor = lower_limit
    for _slot_id, raw_band in sorted(occupied, key=lambda item: _normalize_band(item[1])[0]):
        band_inner, band_outer = _normalize_band(raw_band)
        if band_outer <= lower_limit + LAYOUT_EPSILON:
            continue
        if band_inner >= upper_limit - LAYOUT_EPSILON:
            break
        band_inner = max(lower_limit, band_inner)
        band_outer = min(upper_limit, band_outer)
        if band_inner > cursor + LAYOUT_EPSILON:
            intervals.append((cursor, band_inner))
        cursor = max(cursor, band_outer)
    if cursor < upper_limit - LAYOUT_EPSILON:
        intervals.append((cursor, upper_limit))
    return intervals


def _find_auto_gap_footprint(
    slot: CircularTrackSlot,
    *,
    preferred_anchor_px: float,
    width_px: float,
    params: Mapping[str, Any],
    context: CircularTrackLayoutContext,
    occupied: Sequence[tuple[str, tuple[float, float]]],
    inner_guard_px: float,
    explicit_width: bool,
    gap_px: float | None = None,
    outer_limit_px: float | None = None,
) -> CircularSlotFootprint | None:
    global_outer_limit_px = (
        float(context.auto_start_radius_px)
        if context.auto_start_radius_px is not None
        else float(context.base_radius_px)
    )
    ordered_outer_limit_px = (
        global_outer_limit_px
        if outer_limit_px is None
        else min(global_outer_limit_px, float(outer_limit_px))
    )
    gap_px = max(0.0, float(_slot_gap_after_px(slot, context) if gap_px is None else gap_px))
    candidates: list[tuple[float, CircularSlotFootprint]] = []

    for raw_lower, raw_upper in _free_intervals_px(
        occupied,
        inner_limit_px=inner_guard_px,
        outer_limit_px=ordered_outer_limit_px,
    ):
        lower = float(raw_lower) + gap_px
        upper = float(raw_upper) - gap_px
        if upper < lower - LAYOUT_EPSILON:
            continue

        seed_anchors = [
            float(preferred_anchor_px),
            (lower + upper) / 2.0,
            lower + (0.5 * float(width_px)),
            upper - (0.5 * float(width_px)),
        ]
        for seed_anchor in seed_anchors:
            initial = _measure_slot(
                slot,
                anchor_radius_px=float(seed_anchor),
                width_px=float(width_px),
                params=params,
                context=context,
                explicit_anchor=False,
                explicit_annulus=False,
                explicit_width=bool(explicit_width),
            )
            anchor = float(seed_anchor)
            if initial.reserved_inner_px < lower:
                anchor += lower - initial.reserved_inner_px
            if initial.reserved_outer_px > upper:
                anchor -= initial.reserved_outer_px - upper

            footprint = _measure_slot(
                slot,
                anchor_radius_px=float(anchor),
                width_px=float(width_px),
                params=params,
                context=context,
                explicit_anchor=False,
                explicit_annulus=False,
                explicit_width=bool(explicit_width),
            )
            band = (footprint.reserved_inner_px, footprint.reserved_outer_px)
            if band[0] < lower - LAYOUT_EPSILON or band[1] > upper + LAYOUT_EPSILON:
                continue
            if any(_overlaps(band, other_band) for _other_id, other_band in occupied):
                continue
            candidates.append((abs(float(anchor) - float(preferred_anchor_px)), footprint))

    if not candidates:
        return None
    return min(candidates, key=lambda item: item[0])[1]


def _make_slot_footprint(
    *,
    anchor: float,
    draw_inner: float,
    draw_outer: float,
    hard_inner: float,
    hard_outer: float,
    soft_inner: float,
    soft_outer: float,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    draw_inner, draw_outer = sorted((float(draw_inner), float(draw_outer)))
    hard_inner, hard_outer = sorted((float(hard_inner), float(hard_outer)))
    soft_inner, soft_outer = sorted((float(soft_inner), float(soft_outer)))
    return CircularSlotFootprint(
        anchor_radius_px=float(anchor),
        draw_inner_px=max(0.0, draw_inner),
        draw_outer_px=max(0.0, draw_outer),
        hard_inner_px=max(0.0, hard_inner),
        hard_outer_px=max(0.0, hard_outer),
        soft_inner_px=max(0.0, soft_inner),
        soft_outer_px=max(0.0, soft_outer),
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _legacy_tick_label_soft_band(
    anchor: float,
    context: CircularTrackLayoutContext,
) -> tuple[float, float] | None:
    if context.tick_label_radius_ratio is not None:
        label_center = float(anchor) * float(context.tick_label_radius_ratio)
        extent = max(0.0, float(context.tick_label_extent_px))
        return label_center - extent, label_center + extent
    if context.tick_label_offsets_px is not None:
        return (
            float(anchor) + float(context.tick_label_offsets_px[0]),
            float(anchor) + float(context.tick_label_offsets_px[1]),
        )
    return None


def _tick_label_soft_band(
    *,
    anchor: float,
    width: float,
    params: Mapping[str, Any],
    context: CircularTrackLayoutContext,
    explicit_width: bool,
) -> tuple[float, float] | None:
    label_side = str(params.get("label_side", "legacy")).strip().lower()
    if label_side in {"none", ""}:
        return None
    tick_side = str(params.get("tick_side", "legacy")).strip().lower()

    if context.tick_total_len is not None and context.tick_track_type and context.tick_font_family:
        tick_length_px = float(width) if explicit_width and width > 0 else None
        return get_circular_tick_label_radius_bounds(
            center_radius_px=float(anchor),
            total_len=int(context.tick_total_len),
            track_type=str(context.tick_track_type),
            strandedness=bool(context.tick_strandedness),
            font_size=float(context.tick_font_size_px),
            font_family=str(context.tick_font_family),
            dpi=int(context.tick_dpi),
            manual_interval=context.tick_manual_interval,
            tick_track_channel_override=context.tick_track_channel_override,
            label_side=label_side,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            tick_width=float(context.tick_width_px),
            length_reference_radius_px=float(context.base_radius_px),
        )

    if label_side == "legacy":
        return _legacy_tick_label_soft_band(anchor, context)
    return None


def _measure_feature_slot(
    *,
    anchor: float,
    width: float,
    context: CircularTrackLayoutContext,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    if context.feature_band_offsets_px is not None:
        offset_inner, offset_outer = context.feature_band_offsets_px
        draw_inner = float(anchor) + float(offset_inner)
        draw_outer = float(anchor) + float(offset_outer)
    else:
        draw_inner, draw_outer = _annulus(anchor, width)
    return _make_slot_footprint(
        anchor=anchor,
        draw_inner=draw_inner,
        draw_outer=draw_outer,
        hard_inner=draw_inner,
        hard_outer=draw_outer,
        soft_inner=draw_inner,
        soft_outer=draw_outer,
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _measure_tick_slot(
    *,
    anchor: float,
    width: float,
    params: Mapping[str, Any],
    context: CircularTrackLayoutContext,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    label_side = str(params.get("label_side", "legacy")).strip().lower()
    tick_side = str(params.get("tick_side", "legacy")).strip().lower()

    inner_ext = 0.0
    outer_ext = 0.0
    if context.tick_total_len is not None and context.tick_track_type:
        tick_length_px = float(width) if explicit_width and width > 0 else None
        draw_inner, draw_outer = get_circular_tick_path_radius_bounds(
            center_radius_px=float(anchor),
            total_len=int(context.tick_total_len),
            size="large",
            track_type=str(context.tick_track_type),
            strandedness=bool(context.tick_strandedness),
            tick_track_channel_override=context.tick_track_channel_override,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            length_reference_radius_px=float(context.base_radius_px),
        )
        if tick_side not in {"legacy", "none", ""}:
            inner_ext = max(0.0, float(anchor) - float(draw_inner))
            outer_ext = max(0.0, float(draw_outer) - float(anchor))
    elif tick_side == "legacy":
        if context.tick_path_ratio_bounds is not None:
            draw_inner = float(anchor) * min(context.tick_path_ratio_bounds)
            draw_outer = float(anchor) * max(context.tick_path_ratio_bounds)
        else:
            draw_inner, draw_outer = _annulus(anchor, width)
    else:
        tick_len = width if explicit_width and width > 0 else max(6.0, 0.025 * float(context.base_radius_px))
        if tick_side in {"inside", "both"}:
            inner_ext = max(inner_ext, tick_len)
        if tick_side in {"outside", "both"}:
            outer_ext = max(outer_ext, tick_len)
        if tick_side in {"none", ""}:
            inner_ext = outer_ext = 0.0
        draw_inner = float(anchor) - inner_ext
        draw_outer = float(anchor) + outer_ext

    hard_inner, hard_outer = draw_inner, draw_outer
    soft_inner, soft_outer = draw_inner, draw_outer
    label_band = _tick_label_soft_band(
        anchor=anchor,
        width=width,
        params=params,
        context=context,
        explicit_width=explicit_width,
    )
    if label_band is not None:
        soft_inner = min(soft_inner, label_band[0])
        soft_outer = max(soft_outer, label_band[1])
    else:
        label_pad = max(float(context.tick_font_size_px) * 1.8, 10.0)
        if label_side == "inside":
            soft_inner = min(soft_inner, float(anchor) - inner_ext - label_pad)
        elif label_side == "outside":
            soft_outer = max(soft_outer, float(anchor) + outer_ext + label_pad)

    if context.tick_labels_hard:
        hard_inner, hard_outer = soft_inner, soft_outer
    return _make_slot_footprint(
        anchor=anchor,
        draw_inner=draw_inner,
        draw_outer=draw_outer,
        hard_inner=hard_inner,
        hard_outer=hard_outer,
        soft_inner=soft_inner,
        soft_outer=soft_outer,
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _measure_numeric_slot(
    *,
    anchor: float,
    width: float,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    draw_inner, draw_outer = _annulus(anchor, width)
    return _make_slot_footprint(
        anchor=anchor,
        draw_inner=draw_inner,
        draw_outer=draw_outer,
        hard_inner=draw_inner,
        hard_outer=draw_outer,
        soft_inner=draw_inner,
        soft_outer=draw_outer,
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _measure_spacer_slot(
    *,
    anchor: float,
    width: float,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    return _measure_numeric_slot(
        anchor=anchor,
        width=width,
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _measure_slot(
    slot: CircularTrackSlot,
    *,
    anchor_radius_px: float,
    width_px: float,
    params: Mapping[str, Any],
    context: CircularTrackLayoutContext,
    explicit_anchor: bool,
    explicit_annulus: bool,
    explicit_width: bool,
) -> CircularSlotFootprint:
    renderer = str(slot.renderer)
    anchor = float(anchor_radius_px)
    width = max(0.0, float(width_px))

    if renderer == "features":
        return _measure_feature_slot(
            anchor=anchor,
            width=width,
            context=context,
            explicit_anchor=explicit_anchor,
            explicit_annulus=explicit_annulus,
            explicit_width=explicit_width,
        )
    if renderer == "ticks":
        return _measure_tick_slot(
            anchor=anchor,
            width=width,
            params=params,
            context=context,
            explicit_anchor=explicit_anchor,
            explicit_annulus=explicit_annulus,
            explicit_width=explicit_width,
        )
    if renderer == "spacer":
        return _measure_spacer_slot(
            anchor=anchor,
            width=width,
            explicit_anchor=explicit_anchor,
            explicit_annulus=explicit_annulus,
            explicit_width=explicit_width,
        )
    return _measure_numeric_slot(
        anchor=anchor,
        width=width,
        explicit_anchor=explicit_anchor,
        explicit_annulus=explicit_annulus,
        explicit_width=explicit_width,
    )


def _resolved_from_footprint(
    slot: CircularTrackSlot,
    footprint: CircularSlotFootprint,
    *,
    width_px: float,
    params: Mapping[str, Any],
) -> ResolvedCircularTrackSlot:
    anchor = float(footprint.anchor_radius_px if footprint.anchor_radius_px is not None else 0.0)
    center = anchor if str(slot.renderer) in {"features", "ticks"} else (
        (float(footprint.draw_inner_px) + float(footprint.draw_outer_px)) / 2.0
    )
    return ResolvedCircularTrackSlot(
        id=str(slot.id),
        renderer=str(slot.renderer),
        anchor_radius_px=anchor,
        center_radius_px=center,
        width_px=max(0.0, float(width_px)),
        draw_inner_radius_px=float(footprint.draw_inner_px),
        draw_outer_radius_px=float(footprint.draw_outer_px),
        hard_inner_radius_px=float(footprint.hard_inner_px),
        hard_outer_radius_px=float(footprint.hard_outer_px),
        soft_inner_radius_px=float(footprint.soft_inner_px),
        soft_outer_radius_px=float(footprint.soft_outer_px),
        reserved_inner_radius_px=float(footprint.hard_inner_px),
        reserved_outer_radius_px=float(footprint.hard_outer_px),
        z=int(slot.z),
        params=dict(params),
    )


def _warn_soft_annotation_overlaps(resolved_slots: Sequence[ResolvedCircularTrackSlot]) -> None:
    for slot in resolved_slots:
        soft_band = _normalize_band((slot.soft_inner_radius_px, slot.soft_outer_radius_px))
        hard_band = _normalize_band((slot.hard_inner_radius_px, slot.hard_outer_radius_px))
        if (
            abs(soft_band[0] - hard_band[0]) <= LAYOUT_EPSILON
            and abs(soft_band[1] - hard_band[1]) <= LAYOUT_EPSILON
        ):
            continue

        for other in resolved_slots:
            if other.id == slot.id:
                continue
            other_hard_band = _normalize_band((other.hard_inner_radius_px, other.hard_outer_radius_px))
            if _overlaps(soft_band, other_hard_band) and not _overlaps(hard_band, other_hard_band):
                logger.warning(
                    "Soft annotation footprint for circular track slot '%s' overlaps hard footprint of slot '%s'.",
                    slot.id,
                    other.id,
                )


def _pack_circular_slot_entries(
    entries: Sequence[Mapping[str, Any]],
    *,
    context: CircularTrackLayoutContext,
    width_scale: float = 1.0,
    implicit_gap_scale: float = 1.0,
) -> list[ResolvedCircularTrackSlot]:
    normalized_reserved_bands = tuple(_normalize_band(band) for band in context.reserved_bands_px)
    inner_guard_px = _reserved_auto_inner_guard(context)
    occupied: list[tuple[str, tuple[float, float]]] = [
        (f"reserved:{idx}", band)
        for idx, band in enumerate(normalized_reserved_bands)
        if band[1] > band[0] + LAYOUT_EPSILON and band[0] > LAYOUT_EPSILON
    ]
    resolved_by_id: dict[str, ResolvedCircularTrackSlot] = {}

    for entry in entries:
        if not entry["pinned"]:
            continue
        slot = entry["slot"]
        width_px = _entry_width_px(entry, width_scale)
        footprint = _measure_slot(
            slot,
            anchor_radius_px=float(entry["center_px"]),
            width_px=width_px,
            params=entry["params"],
            context=context,
            explicit_anchor=bool(entry["explicit_anchor"]),
            explicit_annulus=bool(entry["explicit_annulus"]),
            explicit_width=bool(entry["explicit_width"]),
        )
        band = (footprint.reserved_inner_px, footprint.reserved_outer_px)
        for other_id, other_band in occupied:
            if _overlaps(band, other_band):
                logger.warning(
                    "Pinned circular track slot '%s' overlaps pinned slot '%s'.",
                    slot.id,
                    other_id,
                )
        for reserved_band in normalized_reserved_bands:
            if _overlaps(band, reserved_band):
                logger.warning(
                    "Pinned circular track slot '%s' overlaps a reserved circular layout band.",
                    slot.id,
                )
        occupied.append((str(slot.id), band))
        resolved_by_id[str(slot.id)] = _resolved_from_footprint(
            slot,
            footprint,
            width_px=width_px,
            params=entry["params"],
        )

    cursor_outer = (
        float(context.auto_start_radius_px)
        if context.auto_start_radius_px is not None
        else float(context.base_radius_px)
    )
    for entry in entries:
        slot = entry["slot"]
        width_px = _entry_width_px(entry, width_scale)
        gap_after_px = _entry_gap_after_px(entry, context=context, implicit_gap_scale=implicit_gap_scale)
        if entry["pinned"]:
            resolved = resolved_by_id[str(slot.id)]
            cursor_outer = min(
                cursor_outer,
                float(resolved.reserved_inner_radius_px) - gap_after_px,
            )
            continue

        preferred_anchor = entry.get("preferred_center_px")
        if preferred_anchor is None and str(slot.renderer) in {"features", "ticks"}:
            preferred_anchor = _legacy_lookup(context.legacy_centers_px, slot)
        if preferred_anchor is None:
            preferred_anchor = cursor_outer - (0.5 * width_px)
        preferred_anchor = float(preferred_anchor)
        anchor = preferred_anchor
        slot_inner_guard_px = inner_guard_px if _uses_auto_inner_guard(slot) else 0.0

        for _ in range(128):
            footprint = _measure_slot(
                slot,
                anchor_radius_px=anchor,
                width_px=width_px,
                params=entry["params"],
                context=context,
                explicit_anchor=False,
                explicit_annulus=False,
                explicit_width=bool(entry["explicit_width"]),
            )
            if footprint.reserved_inner_px < (slot_inner_guard_px - LAYOUT_EPSILON):
                anchor += slot_inner_guard_px - footprint.reserved_inner_px
                continue
            if footprint.reserved_outer_px > cursor_outer:
                anchor -= footprint.reserved_outer_px - cursor_outer
                continue
            conflicts = [
                band
                for _, band in occupied
                if _overlaps((footprint.reserved_inner_px, footprint.reserved_outer_px), band)
            ]
            if not conflicts:
                break
            nearest_inner = min(float(band[0]) for band in conflicts)
            anchor -= footprint.reserved_outer_px - (nearest_inner - gap_after_px)
        else:
            gap_footprint = _find_auto_gap_footprint(
                slot,
                preferred_anchor_px=preferred_anchor,
                width_px=width_px,
                params=entry["params"],
                context=context,
                occupied=occupied,
                inner_guard_px=slot_inner_guard_px,
                explicit_width=bool(entry["explicit_width"]),
                gap_px=gap_after_px,
                outer_limit_px=cursor_outer,
            )
            if gap_footprint is None:
                raise _auto_pack_error(str(slot.id), inner_guard_px=slot_inner_guard_px)
            footprint = gap_footprint

        if footprint.reserved_inner_px < (slot_inner_guard_px - LAYOUT_EPSILON):
            gap_footprint = _find_auto_gap_footprint(
                slot,
                preferred_anchor_px=preferred_anchor,
                width_px=width_px,
                params=entry["params"],
                context=context,
                occupied=occupied,
                inner_guard_px=slot_inner_guard_px,
                explicit_width=bool(entry["explicit_width"]),
                gap_px=gap_after_px,
                outer_limit_px=cursor_outer,
            )
            if gap_footprint is None:
                raise _auto_pack_error(str(slot.id), inner_guard_px=slot_inner_guard_px)
            footprint = gap_footprint
        if footprint.reserved_inner_px < -LAYOUT_EPSILON:
            raise ValueError(
                f"circular track slot '{slot.id}' cannot fit: inner radius is {footprint.reserved_inner_px:.3f}px"
            )
        occupied.append((str(slot.id), (footprint.reserved_inner_px, footprint.reserved_outer_px)))
        resolved = _resolved_from_footprint(
            slot,
            footprint,
            width_px=width_px,
            params=entry["params"],
        )
        resolved_by_id[str(slot.id)] = resolved
        cursor_outer = min(cursor_outer, float(resolved.reserved_inner_radius_px) - gap_after_px)

    return [resolved_by_id[str(entry["slot"].id)] for entry in entries if str(entry["slot"].id) in resolved_by_id]


def _try_pack_circular_slot_entries(
    entries: Sequence[Mapping[str, Any]],
    *,
    context: CircularTrackLayoutContext,
    width_scale: float,
    implicit_gap_scale: float,
) -> list[ResolvedCircularTrackSlot] | None:
    try:
        return _pack_circular_slot_entries(
            entries,
            context=context,
            width_scale=width_scale,
            implicit_gap_scale=implicit_gap_scale,
        )
    except ValueError:
        return None


def _pack_with_auto_numeric_compression(
    entries: Sequence[Mapping[str, Any]],
    *,
    context: CircularTrackLayoutContext,
    initial_error: ValueError,
) -> list[ResolvedCircularTrackSlot]:
    if not any(bool(entry.get("compressible_width")) for entry in entries):
        raise initial_error

    lower_scale = 0.05
    lower_zero_gap_layout = _try_pack_circular_slot_entries(
        entries,
        context=context,
        width_scale=lower_scale,
        implicit_gap_scale=0.0,
    )
    if lower_zero_gap_layout is None:
        raise initial_error

    lower_scaled_gap_layout = _try_pack_circular_slot_entries(
        entries,
        context=context,
        width_scale=lower_scale,
        implicit_gap_scale=lower_scale,
    )
    scale_gaps_with_width = lower_scaled_gap_layout is not None

    low = lower_scale
    high = 1.0
    best_scale = lower_scale
    best_layout = lower_scaled_gap_layout if scale_gaps_with_width else lower_zero_gap_layout

    for _ in range(20):
        mid = (low + high) / 2.0
        gap_scale = mid if scale_gaps_with_width else 0.0
        layout = _try_pack_circular_slot_entries(
            entries,
            context=context,
            width_scale=mid,
            implicit_gap_scale=gap_scale,
        )
        if layout is None:
            high = mid
            continue
        low = mid
        best_scale = mid
        best_layout = layout

    if best_scale < 0.999:
        logger.info(
            "Auto-compressed circular numeric track widths to %.1f%% to fit reserved radial space.",
            best_scale * 100.0,
        )
    return best_layout


def resolve_circular_track_slots(
    slots: Sequence[CircularTrackSlot],
    *,
    context: CircularTrackLayoutContext,
    legacy_track_specs: Sequence[TrackSpec] | None = None,
    compatibility_mode: bool = False,
) -> list[ResolvedCircularTrackSlot]:
    """Resolve circular slots into measured, non-overlapping radial footprints."""

    active_slots: list[CircularTrackSlot] = []
    seen: set[str] = set()
    for slot in slots:
        if slot.id in seen:
            raise ValueError(f"duplicate circular track slot id: {slot.id}")
        seen.add(slot.id)
        if str(slot.renderer) not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
            raise ValueError(f"unknown circular track renderer: {slot.renderer}")
        if slot.enabled:
            active_slots.append(slot)

    base_radius_px = float(context.base_radius_px)
    entries: list[dict[str, Any]] = []
    for slot in active_slots:
        ts = _track_spec_for_slot(slot, legacy_track_specs)
        if ts is not None and not ts.show:
            continue

        track_spec_explicit_width = _track_spec_has_explicit_circular_width(ts)
        placement = slot.placement
        if ts is not None and isinstance(ts.placement, CircularTrackPlacement) and placement is None:
            placement = ts.placement

        center_px, placement_width_px, explicit_anchor, explicit_annulus, explicit_width = (
            _resolve_placement_center_and_width(placement, base_radius_px=base_radius_px)
        )
        explicit_width = explicit_width or track_spec_explicit_width
        width_px = placement_width_px
        if slot.width is not None:
            width_px = float(slot.width.resolve(base_radius_px))
            explicit_width = True

        preferred_layout = None
        if center_px is None and not explicit_annulus:
            preferred_layout = _preferred_layout_lookup(context.preferred_layouts_px, slot)
            if preferred_layout is not None and width_px is None:
                width_px = preferred_layout[1]

        if width_px is None:
            width_px = _legacy_lookup(context.legacy_widths_px, slot)
        if width_px is None:
            width_px = 0.0

        if center_px is None and (compatibility_mode or explicit_anchor):
            center_px = _legacy_lookup(context.legacy_centers_px, slot)

        params = dict(slot.params)
        if ts is not None and ts.params:
            params.update(ts.params)

        entries.append(
            {
                "slot": slot,
                "center_px": center_px,
                "preferred_center_px": preferred_layout[0] if preferred_layout is not None else None,
                "width_px": max(0.0, float(width_px)),
                "params": params,
                "explicit_anchor": bool(explicit_anchor),
                "explicit_annulus": bool(explicit_annulus),
                "explicit_width": bool(explicit_width),
                "explicit_gap": slot.gap_after is not None,
                "pinned": center_px is not None,
                "compressible_width": _is_compressible_numeric_width(
                    slot,
                    explicit_annulus=bool(explicit_annulus),
                    explicit_width=bool(explicit_width),
                ),
            }
        )

    if compatibility_mode:
        resolved_compat: list[ResolvedCircularTrackSlot] = []
        for entry in entries:
            slot = entry["slot"]
            center_px = entry["center_px"]
            if center_px is None:
                center_px = _legacy_lookup(context.legacy_centers_px, slot)
            if center_px is None:
                center_px = float(context.base_radius_px)
            footprint = _measure_slot(
                slot,
                anchor_radius_px=float(center_px),
                width_px=float(entry["width_px"]),
                params=entry["params"],
                context=context,
                explicit_anchor=bool(entry["explicit_anchor"]),
                explicit_annulus=bool(entry["explicit_annulus"]),
                explicit_width=bool(entry["explicit_width"]),
            )
            resolved_compat.append(
                _resolved_from_footprint(
                    slot,
                    footprint,
                    width_px=float(entry["width_px"]),
                    params=entry["params"],
                )
            )
        _warn_soft_annotation_overlaps(resolved_compat)
        return resolved_compat

    try:
        resolved = _pack_circular_slot_entries(entries, context=context, width_scale=1.0)
    except CircularSlotAutoPackError as exc:
        resolved = _pack_with_auto_numeric_compression(entries, context=context, initial_error=exc)
    _warn_soft_annotation_overlaps(resolved)
    return resolved


__all__ = [
    "CircularSlotFootprint",
    "CircularTrackLayoutContext",
    "ResolvedCircularTrackSlot",
    "resolve_circular_track_slots",
]
