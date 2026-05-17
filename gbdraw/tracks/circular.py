from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Literal, Mapping, Sequence

from .parsing import CircularTrackSlotParseError, parse_bool, split_kv_list
from .scalars import ScalarSpec


CircularTrackRendererName = Literal[
    "features",
    "ticks",
    "dinucleotide_content",
    "dinucleotide_skew",
    "depth",
    "spacer",
]

CircularTrackSide = Literal["inside", "outside", "overlay"]

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

NUMERIC_CIRCULAR_TRACK_RENDERERS: frozenset[str] = frozenset(
    {"dinucleotide_content", "dinucleotide_skew", "depth"}
)

_RENDERER_ALIASES: dict[str, str] = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
}

_ORDER_RENDERER_BY_ID: dict[str, str] = {
    "features": "features",
    "ticks": "ticks",
    "depth": "depth",
    "gc_content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
}

_OBSOLETE_GEOMETRY_KEYS = {
    "ri",
    "inner",
    "inner_radius",
    "ro",
    "outer",
    "outer_radius",
}

_OBSOLETE_SPACING_KEYS = {"gap", "gap_after"}
_OBSOLETE_CAMEL_KEYS = {
    "gapafter",
    "innerradius",
    "outerradius",
}
_GENERIC_LAYOUT_KEYS = {
    "side",
    "r",
    "radius",
    "w",
    "width",
    "spacing",
    "z",
    "z_index",
    "zindex",
    "strict",
    "compress",
    "reserve",
    "enabled",
    "show",
    "visible",
}
_TICK_SIDE_VALUES = {"inside", "outside", "both", "none"}
_FEATURE_LANE_VALUES = {"inside", "outside", "split"}
_SIDE_VALUES = {"inside", "outside", "overlay"}


@dataclass(frozen=True)
class CircularTrackSlot:
    """One public radial slot input for a circular diagram."""

    id: str
    renderer: CircularTrackRendererName | str
    enabled: bool = True
    side: str | None = None
    radius: ScalarSpec | None = None
    width: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z: int = 0
    strict: bool | None = None
    compress: bool | None = None
    reserve: bool | None = None
    params: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class NormalizedCircularTrackSlot:
    slot_index: int
    id: str
    renderer: str
    enabled: bool
    side: str
    radius: ScalarSpec | None
    width: ScalarSpec | None
    spacing: ScalarSpec | None
    z: int
    strict: bool
    compress: bool
    reserve: bool
    params: Mapping[str, Any]


def _normalize_renderer(raw: str) -> str:
    renderer = str(raw).strip().lower()
    return _RENDERER_ALIASES.get(renderer, renderer)


def _normalize_side_value(raw: object, *, field_name: str = "side") -> str:
    side = str(raw).strip().lower()
    if side not in _SIDE_VALUES:
        raise ValueError(f"{field_name} must be one of inside, outside, overlay")
    return side


def _normalize_tick_side(raw: object, *, field_name: str) -> str:
    side = str(raw).strip().lower()
    if side not in _TICK_SIDE_VALUES:
        raise ValueError(f"{field_name} must be one of inside, outside, both, none")
    return side


def _normalize_feature_lane(raw: object) -> str:
    lane = str(raw).strip().lower()
    if lane not in _FEATURE_LANE_VALUES:
        raise ValueError("lane_direction must be one of inside, outside, split")
    return lane


def _parse_slot_head(head: str, original: str) -> tuple[str, str]:
    if ":" not in head:
        raise CircularTrackSlotParseError(
            "circular track slots require '<slot_id>:<renderer>@...'; shortcut slot strings are no longer supported",
            original,
        )
    slot_id_raw, renderer_raw = head.split(":", 1)
    slot_id = slot_id_raw.strip()
    renderer = renderer_raw.strip()
    if not slot_id:
        raise CircularTrackSlotParseError("missing circular track slot id", original)
    if not renderer:
        raise CircularTrackSlotParseError("missing circular track renderer", original)
    return slot_id, renderer


def parse_circular_track_slot(raw: str) -> CircularTrackSlot:
    """Parse `<slot_id>:<renderer>@key=value,...` into a slot input object."""

    original = raw
    s = str(raw).strip()
    if not s or s.startswith("#"):
        raise CircularTrackSlotParseError("empty/comment line", original)
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
        raise CircularTrackSlotParseError(f"unknown circular track renderer '{renderer}'", original)

    enabled = True
    side: str | None = None
    radius: ScalarSpec | None = None
    width: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z = 0
    strict: bool | None = None
    compress: bool | None = None
    reserve: bool | None = None
    params: dict[str, Any] = {}

    if opts:
        try:
            for raw_key, raw_value in split_kv_list(opts):
                key = raw_key.strip().lower()
                value = raw_value.strip()
                if key in {"id"}:
                    slot_id = value
                elif key in {"renderer", "type"}:
                    renderer = _normalize_renderer(value)
                    if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
                        raise ValueError(f"unknown circular track renderer '{renderer}'")
                elif key in {"enabled", "show", "visible"}:
                    enabled = parse_bool(value)
                elif key in {"z", "z_index", "zindex"}:
                    z = int(value)
                elif key in {"r", "radius"}:
                    radius = ScalarSpec.parse(value)
                elif key in {"w", "width"}:
                    width = ScalarSpec.parse(value)
                elif key == "spacing":
                    spacing = ScalarSpec.parse(value)
                elif key in _OBSOLETE_GEOMETRY_KEYS:
                    raise ValueError(f"'{key}' is no longer supported; use r=<radius> with w=<width>")
                elif key in {"innerradius", "outerradius"}:
                    raise ValueError(f"'{raw_key}' is no longer supported; use r=<radius> with w=<width>")
                elif key in _OBSOLETE_SPACING_KEYS or key == "gapafter":
                    raise ValueError(f"'{key}' is no longer supported; use spacing=<ScalarSpec>")
                elif key == "side":
                    side = _normalize_side_value(value)
                elif key == "strict":
                    strict = parse_bool(value)
                elif key == "compress":
                    compress = parse_bool(value)
                elif key == "reserve":
                    reserve = parse_bool(value)
                elif key in {"nt", "dinucleotide"}:
                    params["nt"] = value.upper()
                else:
                    params[key] = value
        except CircularTrackSlotParseError:
            raise
        except Exception as exc:
            raise CircularTrackSlotParseError(str(exc), original) from exc

    if not slot_id:
        raise CircularTrackSlotParseError("missing circular track slot id", original)

    slot = CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        enabled=enabled,
        side=side,
        radius=radius,
        width=width,
        spacing=spacing,
        z=z,
        strict=strict,
        compress=compress,
        reserve=reserve,
        params=params,
    )
    try:
        normalize_circular_track_slots([slot])
    except Exception as exc:
        raise CircularTrackSlotParseError(str(exc), original) from exc
    return slot


def parse_circular_track_slots(specs: Sequence[str | CircularTrackSlot]) -> list[CircularTrackSlot]:
    """Parse and validate circular slot specs."""

    out: list[CircularTrackSlot] = []
    seen: set[str] = set()
    for item in specs:
        if isinstance(item, CircularTrackSlot):
            renderer = _normalize_renderer(str(item.renderer))
            slot = replace(item, renderer=renderer) if renderer != str(item.renderer) else item
            try:
                normalize_circular_track_slots([slot])
            except Exception as exc:
                raise CircularTrackSlotParseError(str(exc), str(slot.id)) from exc
        else:
            slot = parse_circular_track_slot(str(item))
        if slot.id in seen:
            raise CircularTrackSlotParseError("duplicate circular track slot id", slot.id)
        seen.add(slot.id)
        out.append(slot)
    normalize_circular_track_slots(out)
    return out


def _normalized_feature_side_and_params(slot: CircularTrackSlot, params: dict[str, Any]) -> tuple[str, bool, dict[str, Any]]:
    raw_lane = params.get("lane_direction", params.get("lanes"))
    raw_side = slot.side
    reserve = bool(slot.reserve) if slot.reserve is not None else False

    lane: str | None = _normalize_feature_lane(raw_lane) if raw_lane is not None else None
    side: str | None = _normalize_side_value(raw_side) if raw_side is not None else None

    if lane is None and side is None:
        lane = "inside"
        side = "inside"
    elif lane is None:
        if side == "inside":
            lane = "inside"
        elif side == "outside":
            lane = "outside"
        else:
            lane = "split"
            reserve = True if slot.reserve is None else reserve
    elif side is None:
        if lane == "inside":
            side = "inside"
        elif lane == "outside":
            side = "outside"
        else:
            side = "overlay"
            reserve = True if slot.reserve is None else reserve
    else:
        expected = "overlay" if lane == "split" else lane
        if side != expected:
            raise ValueError(
                f"features slot '{slot.id}' has conflicting side={side!r} and lane_direction={lane!r}"
            )
        if side == "overlay":
            reserve = True if slot.reserve is None else reserve

    params.pop("lanes", None)
    params["lane_direction"] = lane
    return str(side), reserve, params


def _normalized_tick_params(slot: CircularTrackSlot, params: dict[str, Any]) -> dict[str, Any]:
    if "axis" in {str(key).strip().lower() for key in params}:
        raise ValueError("ticks slots no longer accept 'axis'; the circular axis is fixed and not a slot")
    if "label_side" in params:
        params["label_side"] = _normalize_tick_side(params["label_side"], field_name="label_side")
    else:
        params["label_side"] = "inside"
    if "tick_side" in params:
        params["tick_side"] = _normalize_tick_side(params["tick_side"], field_name="tick_side")
    else:
        params["tick_side"] = "inside"
    return params


def normalize_circular_track_slots(slots: Sequence[CircularTrackSlot]) -> list[NormalizedCircularTrackSlot]:
    """Validate and normalize slots for the circular radial resolver."""

    normalized: list[NormalizedCircularTrackSlot] = []
    seen: set[str] = set()
    for slot_index, slot in enumerate(slots):
        if slot.id in seen:
            raise ValueError(f"duplicate circular track slot id: {slot.id}")
        seen.add(str(slot.id))
        renderer = _normalize_renderer(str(slot.renderer))
        if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
            raise ValueError(f"unknown circular track renderer: {slot.renderer}")
        if not slot.enabled:
            continue

        raw_params = {str(key): value for key, value in dict(slot.params or {}).items()}
        for raw_key in raw_params:
            key = raw_key.strip()
            normalized_key = key.lower()
            if (
                normalized_key in _OBSOLETE_CAMEL_KEYS
                or normalized_key in _GENERIC_LAYOUT_KEYS
                or normalized_key in _OBSOLETE_GEOMETRY_KEYS
                or normalized_key in _OBSOLETE_SPACING_KEYS
            ):
                raise ValueError(
                    f"circular track slot '{slot.id}' stores generic layout field '{raw_key}' in params; "
                    "use slot-level radius, width, spacing, side, strict, compress, reserve, and z fields"
                )
        params = dict(raw_params)
        side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
        reserve = bool(slot.reserve) if slot.reserve is not None else False

        if renderer == "features":
            side, reserve, params = _normalized_feature_side_and_params(slot, params)
        elif renderer == "ticks":
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
            params = _normalized_tick_params(slot, params)
        elif renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
        elif renderer == "spacer":
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"

        compress = bool(slot.compress) if slot.compress is not None else False
        explicit_anchor = slot.radius is not None
        if compress:
            if renderer not in NUMERIC_CIRCULAR_TRACK_RENDERERS:
                raise ValueError(f"compress is not supported for {renderer} circular slots")
            if explicit_anchor:
                raise ValueError("compress is not supported on pinned circular numeric/depth slots")
            if side != "inside":
                raise ValueError("compress is valid only for auto inside numeric/depth circular slots")

        normalized.append(
            NormalizedCircularTrackSlot(
                slot_index=int(slot_index),
                id=str(slot.id),
                renderer=renderer,
                enabled=True,
                side=side,
                radius=slot.radius,
                width=slot.width,
                spacing=slot.spacing,
                z=int(slot.z),
                strict=bool(slot.strict) if slot.strict is not None else False,
                compress=compress,
                reserve=reserve,
                params=params,
            )
        )
    return normalized


def default_circular_track_slots(
    *,
    show_features: bool = True,
    show_ticks: bool = True,
    show_depth: bool = False,
    show_gc: bool = True,
    show_skew: bool = True,
    dinucleotide: str = "GC",
) -> list[CircularTrackSlot]:
    """Return the default circular slot input list."""

    slots: list[CircularTrackSlot] = []
    nt = str(dinucleotide or "GC").upper()
    if show_features:
        slots.append(
            CircularTrackSlot(
                id="features",
                renderer="features",
            )
        )
    if show_ticks:
        slots.append(
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
            )
        )
    if show_depth:
        slots.append(CircularTrackSlot(id="depth", renderer="depth"))
    if show_gc:
        slots.append(CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": nt}))
    if show_skew:
        slots.append(CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": nt}))
    return slots


def circular_track_slots_from_order(
    order: str | Sequence[str],
    *,
    show_features: bool = True,
    show_ticks: bool = True,
    show_depth: bool = False,
    show_gc: bool = True,
    show_skew: bool = True,
    dinucleotide: str = "GC",
) -> list[CircularTrackSlot]:
    """Expand a comma-separated slot order into explicit slot inputs."""

    enabled = {
        "features": show_features,
        "ticks": show_ticks,
        "depth": show_depth,
        "gc_content": show_gc,
        "gc_skew": show_skew,
    }
    nt = str(dinucleotide or "GC").upper()
    items = order.split(",") if isinstance(order, str) else list(order)
    slots: list[CircularTrackSlot] = []
    seen: set[str] = set()
    for raw in items:
        slot_id = str(raw).strip()
        if not slot_id:
            continue
        if slot_id in seen:
            raise ValueError(f"duplicate circular slot id in order: {slot_id}")
        if slot_id not in _ORDER_RENDERER_BY_ID:
            raise ValueError(f"unknown circular slot id in order: {slot_id}")
        seen.add(slot_id)
        if not enabled.get(slot_id, True):
            continue
        params: dict[str, Any] = {}
        renderer = _ORDER_RENDERER_BY_ID[slot_id]
        if renderer in {"dinucleotide_content", "dinucleotide_skew"}:
            params["nt"] = nt
        slots.append(
            CircularTrackSlot(
                id=slot_id,
                renderer=renderer,
                params=params,
            )
        )
    return slots


__all__ = [
    "CircularTrackRendererName",
    "CircularTrackSide",
    "CircularTrackSlot",
    "CircularTrackSlotParseError",
    "NormalizedCircularTrackSlot",
    "SUPPORTED_CIRCULAR_TRACK_RENDERERS",
    "circular_track_slots_from_order",
    "default_circular_track_slots",
    "normalize_circular_track_slots",
    "parse_circular_track_slot",
    "parse_circular_track_slots",
]
