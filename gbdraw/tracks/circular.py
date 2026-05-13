from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Literal, Mapping, Sequence

from .parser import TrackSpecParseError, _parse_bool, _split_kv_list
from .spec import CircularTrackPlacement, ScalarSpec, TrackSpec


CircularTrackRendererName = Literal[
    "features",
    "ticks",
    "dinucleotide_content",
    "dinucleotide_skew",
    "depth",
    "spacer",
]

SUPPORTED_CIRCULAR_TRACK_RENDERERS: frozenset[str] = frozenset(
    {
        "features",
        "ticks",
        "dinucleotide_content",
        "dinucleotide_skew",
        "depth",
        "spacer",
    }
)

_LEGACY_SLOT_RENDERER_BY_ID: dict[str, str] = {
    "features": "features",
    "ticks": "ticks",
    "depth": "depth",
    "gc_content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
}

_LEGACY_KIND_BY_RENDERER: dict[str, str] = {
    "features": "features",
    "ticks": "ticks",
    "depth": "depth",
    "dinucleotide_content": "gc_content",
    "dinucleotide_skew": "gc_skew",
}

_TICK_AXIS_PARAM_ERROR = "ticks slots no longer accept 'axis'; the circular axis is fixed and not a slot"


@dataclass(frozen=True)
class CircularTrackSlot:
    """One radial slot in a circular diagram.

    `id` identifies the slot instance. `renderer` selects the drawing
    implementation, so several slots may share a renderer.
    """

    id: str
    renderer: CircularTrackRendererName | str
    enabled: bool = True
    placement: CircularTrackPlacement | None = None
    width: ScalarSpec | None = None
    gap_after: ScalarSpec | None = None
    z: int = 0
    params: Mapping[str, Any] = field(default_factory=dict)


def _normalize_renderer(raw: str) -> str:
    renderer = str(raw).strip().lower()
    aliases = {
        "gc_content": "dinucleotide_content",
        "content": "dinucleotide_content",
        "gc_skew": "dinucleotide_skew",
        "skew": "dinucleotide_skew",
    }
    return aliases.get(renderer, renderer)


def _validate_circular_track_slot_params(
    *,
    renderer: str,
    params: Mapping[str, Any],
    original: str,
) -> None:
    if renderer == "ticks" and any(str(key).strip().lower() == "axis" for key in params):
        raise TrackSpecParseError(_TICK_AXIS_PARAM_ERROR, original)


def _parse_slot_head(head: str, original: str) -> tuple[str, str]:
    if ":" not in head:
        slot_id = head.strip()
        if not slot_id:
            raise TrackSpecParseError("missing circular track slot id", original)
        renderer = _LEGACY_SLOT_RENDERER_BY_ID.get(slot_id, slot_id)
        return slot_id, renderer

    slot_id_raw, renderer_raw = head.split(":", 1)
    slot_id = slot_id_raw.strip()
    renderer = renderer_raw.strip()
    if not slot_id:
        raise TrackSpecParseError("missing circular track slot id", original)
    if not renderer:
        raise TrackSpecParseError("missing circular track renderer", original)
    return slot_id, renderer


def parse_circular_track_slot(raw: str) -> CircularTrackSlot:
    """Parse `<slot_id>:<renderer>@key=value,...` into a slot object."""

    original = raw
    s = str(raw).strip()
    if not s or s.startswith("#"):
        raise TrackSpecParseError("empty/comment line", original)

    if "#" in s:
        s = s.split("#", 1)[0].strip()

    if "@" in s:
        head, opts = s.split("@", 1)
        opts = opts.strip()
    else:
        head, opts = s, ""

    slot_id, renderer_raw = _parse_slot_head(head.strip(), original)
    renderer = _normalize_renderer(renderer_raw)
    if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
        raise TrackSpecParseError(f"unknown circular track renderer '{renderer}'", original)

    enabled = True
    placement: CircularTrackPlacement | None = None
    width: ScalarSpec | None = None
    gap_after: ScalarSpec | None = None
    z = 0
    params: dict[str, Any] = {}

    if opts:
        try:
            for raw_key, raw_value in _split_kv_list(opts):
                key = raw_key.strip().lower()
                value = raw_value.strip()
                if key in {"id"}:
                    slot_id = value
                elif key in {"renderer", "type"}:
                    renderer = _normalize_renderer(value)
                    if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
                        raise ValueError(f"unknown circular track renderer '{renderer}'")
                elif key in {"enabled", "show", "visible"}:
                    enabled = _parse_bool(value)
                elif key in {"z", "z_index", "zindex"}:
                    z = int(value)
                elif key in {"gap", "gap_after"}:
                    if value.strip().lower() not in {"", "auto", "legacy", "none"}:
                        gap_after = ScalarSpec.parse(value)
                elif key in {"w", "width"}:
                    if value.strip().lower() not in {"", "auto", "legacy"}:
                        width = ScalarSpec.parse(value)
                elif key in {"r", "radius", "ri", "inner", "inner_radius", "ro", "outer", "outer_radius"}:
                    if placement is None:
                        placement = CircularTrackPlacement()
                    if key in {"r", "radius"}:
                        placement = replace(placement, radius=ScalarSpec.parse(value))
                    elif key in {"ri", "inner", "inner_radius"}:
                        placement = replace(placement, inner_radius=ScalarSpec.parse(value))
                    else:
                        placement = replace(placement, outer_radius=ScalarSpec.parse(value))
                elif key in {"placement"}:
                    params["placement"] = value
                elif key in {"nt", "dinucleotide"}:
                    params["nt"] = value.upper()
                else:
                    params[key] = value
        except Exception as exc:
            raise TrackSpecParseError(str(exc), original) from exc

    if placement is not None:
        placement = replace(placement, z=z)
    if not slot_id:
        raise TrackSpecParseError("missing circular track slot id", original)
    _validate_circular_track_slot_params(renderer=renderer, params=params, original=original)

    return CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        enabled=enabled,
        placement=placement,
        width=width,
        gap_after=gap_after,
        z=z,
        params=params,
    )


def parse_circular_track_slots(specs: Sequence[str | CircularTrackSlot]) -> list[CircularTrackSlot]:
    """Parse and validate circular slot specs."""

    out: list[CircularTrackSlot] = []
    seen: set[str] = set()
    for item in specs:
        if isinstance(item, CircularTrackSlot):
            renderer = _normalize_renderer(str(item.renderer))
            slot = replace(item, renderer=renderer) if renderer != str(item.renderer) else item
        else:
            slot = parse_circular_track_slot(str(item))
        if slot.id in seen:
            raise TrackSpecParseError("duplicate circular track slot id", slot.id)
        if str(slot.renderer) not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
            raise TrackSpecParseError(
                f"unknown circular track renderer '{slot.renderer}'",
                str(slot.id),
            )
        _validate_circular_track_slot_params(
            renderer=str(slot.renderer),
            params=slot.params,
            original=str(slot.id),
        )
        seen.add(slot.id)
        out.append(slot)
    return out


def default_circular_track_slots(
    *,
    show_features: bool = True,
    show_ticks: bool = True,
    show_depth: bool = False,
    show_gc: bool = True,
    show_skew: bool = True,
    dinucleotide: str = "GC",
) -> list[CircularTrackSlot]:
    """Return the legacy-compatible default slot list."""

    slots: list[CircularTrackSlot] = []
    nt = str(dinucleotide or "GC").upper()
    if show_features:
        slots.append(CircularTrackSlot(id="features", renderer="features"))
    if show_ticks:
        slots.append(
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"placement": "legacy_axis", "label_side": "legacy", "tick_side": "legacy"},
            )
        )
    if show_depth:
        slots.append(CircularTrackSlot(id="depth", renderer="depth"))
    if show_gc:
        slots.append(CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": nt}))
    if show_skew:
        slots.append(CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": nt}))
    return slots


def circular_track_slots_from_order(order: str, *, dinucleotide: str = "GC") -> list[CircularTrackSlot]:
    """Build slots from a comma-separated legacy/default slot id order."""

    specs: list[str] = []
    nt = str(dinucleotide or "GC").upper()
    for raw_part in str(order).split(","):
        part = raw_part.strip()
        if not part:
            continue
        renderer = _LEGACY_SLOT_RENDERER_BY_ID.get(part, _normalize_renderer(part))
        suffix = ""
        if renderer in {"dinucleotide_content", "dinucleotide_skew"}:
            suffix = f"@nt={nt}"
        specs.append(f"{part}:{renderer}{suffix}")
    return parse_circular_track_slots(specs)


def _legacy_kind_for_slot(slot: CircularTrackSlot) -> str | None:
    if slot.id in _LEGACY_SLOT_RENDERER_BY_ID:
        if slot.id == "gc_content":
            return "gc_content"
        if slot.id == "gc_skew":
            return "gc_skew"
        return slot.id
    return _LEGACY_KIND_BY_RENDERER.get(str(slot.renderer))


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


def _resolve_placement_center_and_width(
    placement: CircularTrackPlacement | None,
    *,
    base_radius_px: float,
) -> tuple[float | None, float | None]:
    if placement is None:
        return None, None
    inner_px = placement.inner_radius.resolve(base_radius_px) if placement.inner_radius is not None else None
    outer_px = placement.outer_radius.resolve(base_radius_px) if placement.outer_radius is not None else None
    radius_px = placement.radius.resolve(base_radius_px) if placement.radius is not None else None
    width_px = placement.width.resolve(base_radius_px) if placement.width is not None else None

    if inner_px is not None and outer_px is not None:
        if outer_px < inner_px:
            inner_px, outer_px = outer_px, inner_px
        return (inner_px + outer_px) / 2.0, outer_px - inner_px
    return radius_px, width_px


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


def _annulus(center_px: float, width_px: float) -> tuple[float, float]:
    half_width = max(0.0, 0.5 * float(width_px))
    return float(center_px) - half_width, float(center_px) + half_width


def _slot_gap_after_px(slot: CircularTrackSlot, context: CircularTrackLayoutContext) -> float:
    if slot.gap_after is None:
        return float(context.default_gap_px)
    return float(slot.gap_after.resolve(float(context.base_radius_px)))


def _resolve_circular_track_slots_legacy(
    slots: Sequence[CircularTrackSlot],
    *,
    context: CircularTrackLayoutContext,
    legacy_track_specs: Sequence[TrackSpec] | None = None,
    compatibility_mode: bool = False,
) -> list[ResolvedCircularTrackSlot]:
    """Resolve circular slots into pixel annuli.

    In compatibility mode, legacy center/width values are preferred for known
    built-in slots. Outside compatibility mode, slots without an explicit center
    are packed outer-to-inner in the order provided.
    """

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

    entries: list[dict[str, Any]] = []
    base_radius_px = float(context.base_radius_px)
    for slot in active_slots:
        ts = _track_spec_for_slot(slot, legacy_track_specs)
        if ts is not None and not ts.show:
            continue

        placement = slot.placement
        if ts is not None and isinstance(ts.placement, CircularTrackPlacement):
            placement = ts.placement

        center_px, placement_width_px = _resolve_placement_center_and_width(
            placement,
            base_radius_px=base_radius_px,
        )
        width_px = placement_width_px
        if slot.width is not None:
            width_px = float(slot.width.resolve(base_radius_px))

        if width_px is None:
            width_px = _legacy_lookup(context.legacy_widths_px, slot)
        if width_px is None:
            width_px = 0.0

        if compatibility_mode and center_px is None:
            center_px = _legacy_lookup(context.legacy_centers_px, slot)

        params = dict(slot.params)
        if ts is not None and ts.params:
            params.update(ts.params)

        entries.append(
            {
                "slot": slot,
                "center_px": center_px,
                "width_px": max(0.0, float(width_px)),
                "params": params,
            }
        )

    cursor_outer = (
        float(context.auto_start_radius_px)
        if context.auto_start_radius_px is not None
        else float(context.base_radius_px)
    )
    resolved: list[ResolvedCircularTrackSlot] = []
    for entry in entries:
        slot = entry["slot"]
        width_px = float(entry["width_px"])
        center_px = entry["center_px"]

        if center_px is None:
            center_px = cursor_outer - (0.5 * width_px)

        inner_px, outer_px = _annulus(float(center_px), width_px)
        if inner_px < -1e-6:
            raise ValueError(
                f"circular track slot '{slot.id}' cannot fit: inner radius is {inner_px:.3f}px"
            )
        inner_px = max(0.0, inner_px)

        resolved.append(
            ResolvedCircularTrackSlot(
                id=str(slot.id),
                renderer=str(slot.renderer),
                center_radius_px=float(center_px),
                width_px=width_px,
                inner_radius_px=float(inner_px),
                outer_radius_px=float(outer_px),
                z=int(slot.z),
                params=entry["params"],
            )
        )
        cursor_outer = min(cursor_outer, float(inner_px) - _slot_gap_after_px(slot, context))

    return resolved


__all__ = [
    "CircularTrackLayoutContext",
    "CircularTrackRendererName",
    "CircularSlotFootprint",
    "CircularTrackSlot",
    "ResolvedCircularTrackSlot",
    "SUPPORTED_CIRCULAR_TRACK_RENDERERS",
    "circular_track_slots_from_order",
    "default_circular_track_slots",
    "parse_circular_track_slot",
    "parse_circular_track_slots",
    "resolve_circular_track_slots",
]


from ..diagrams.circular.slot_layout import (  # noqa: E402  # type: ignore[reportMissingImports]
    CircularSlotFootprint,
    CircularTrackLayoutContext,
    ResolvedCircularTrackSlot,
    resolve_circular_track_slots,
)
