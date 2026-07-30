from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import isfinite
from typing import Any, Literal, Mapping, Sequence

from .parsing import (
    CircularTrackSlotParseError,
    normalize_dinucleotide_skew_color_params as _normalize_dinucleotide_skew_color_params,
    parse_bool,
    parse_nonnegative_integer,
    parse_track_slot_text,
    split_kv_list,
    validate_overlay_annotation_anchors,
)
from .scalars import ScalarSpec
from gbdraw.annotations.models import AnnotationTrackParams, annotation_track_params_from_mapping


CircularTrackRendererName = Literal[
    "features",
    "ticks",
    "dinucleotide_content",
    "dinucleotide_skew",
    "depth",
    "sequence_conservation",
    "annotations",
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
        "sequence_conservation",
        "annotations",
        "spacer",
    }
)

NUMERIC_CIRCULAR_TRACK_RENDERERS: frozenset[str] = frozenset(
    {"dinucleotide_content", "dinucleotide_skew", "depth", "sequence_conservation"}
)

CONSERVATION_CIRCULAR_TRACK_RENDERERS: frozenset[str] = frozenset(
    {"sequence_conservation"}
)

_RENDERER_ALIASES: dict[str, str] = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
    "conservation": "sequence_conservation",
    "circular_comparison": "sequence_conservation",
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
    "inner_gap_px",
    "outer_gap_px",
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
_TICK_LABEL_LAYOUT_VALUES = {
    "label_out_tick_in",
    "label_in_tick_out",
    "tick_only",
    "label_only",
}
DEFAULT_TICK_LABEL_LAYOUT = "label_out_tick_in"
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
    z: int = 0
    params: Mapping[str, Any] = field(default_factory=dict)
    inner_gap_px: float | None = None
    outer_gap_px: float | None = None


@dataclass(frozen=True)
class _InternalCircularTrackSlot(CircularTrackSlot):
    """Reader/planner-only state that must never be accepted as public params."""

    legacy_spacing: ScalarSpec | None = None
    auto_compress: bool = False
    preferred_anchor_radius: ScalarSpec | None = None


def _internal_circular_track_slot(
    slot: CircularTrackSlot,
    *,
    legacy_spacing: ScalarSpec | None = None,
    auto_compress: bool = False,
    preferred_anchor_radius: ScalarSpec | None = None,
) -> _InternalCircularTrackSlot:
    """Attach trusted reader/planner state without exposing private param keys."""

    existing = (
        slot
        if isinstance(slot, _InternalCircularTrackSlot)
        else None
    )
    return _InternalCircularTrackSlot(
        id=slot.id,
        renderer=slot.renderer,
        enabled=slot.enabled,
        side=slot.side,
        radius=slot.radius,
        width=slot.width,
        z=slot.z,
        params=slot.params,
        inner_gap_px=slot.inner_gap_px,
        outer_gap_px=slot.outer_gap_px,
        legacy_spacing=(
            legacy_spacing
            if legacy_spacing is not None
            else (existing.legacy_spacing if existing is not None else None)
        ),
        auto_compress=(
            auto_compress
            or (existing.auto_compress if existing is not None else False)
        ),
        preferred_anchor_radius=(
            preferred_anchor_radius
            if preferred_anchor_radius is not None
            else (
                existing.preferred_anchor_radius
                if existing is not None
                else None
            )
        ),
    )


@dataclass(frozen=True)
class NormalizedCircularTrackSlot:
    slot_index: int
    id: str
    renderer: str
    enabled: bool
    side: str
    radius: ScalarSpec | None
    width: ScalarSpec | None
    legacy_spacing: ScalarSpec | None
    inner_gap_px: float | None
    outer_gap_px: float | None
    z: int
    compress: bool
    reserve: bool
    params: Mapping[str, Any]
    annotation: AnnotationTrackParams | None = None
    preferred_anchor_radius: ScalarSpec | None = None


def _normalize_renderer(raw: str) -> str:
    renderer = str(raw).strip().lower()
    return _RENDERER_ALIASES.get(renderer, renderer)




def _parse_gap_px(raw: object, *, field_name: str) -> float:
    text = str(raw).strip()
    if not text:
        raise ValueError(f"{field_name} must be a numeric pixel value")
    if text.endswith("px") or text.endswith("%"):
        raise ValueError(f"{field_name} must be a numeric pixel value without a unit")
    value = float(text)
    if not isfinite(value) or value < 0:
        raise ValueError(f"{field_name} must be a nonnegative numeric pixel value")
    return value


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


def normalize_tick_label_layout(raw: object | None) -> str:
    layout = str(raw or DEFAULT_TICK_LABEL_LAYOUT).strip().lower()
    if layout not in _TICK_LABEL_LAYOUT_VALUES:
        raise ValueError(
            "tick_label_layout must be one of "
            "label_out_tick_in, label_in_tick_out, tick_only, label_only"
        )
    return layout


def tick_sides_for_tick_label_layout(layout: object | None, side: object | None = None) -> tuple[str, str]:
    normalized_layout = normalize_tick_label_layout(layout)
    normalized_side = _normalize_side_value(side) if side is not None else "inside"
    single_side = "outside" if normalized_side == "outside" else "inside"
    if normalized_layout == "label_out_tick_in":
        return "outside", "inside"
    if normalized_layout == "label_in_tick_out":
        return "inside", "outside"
    if normalized_layout == "tick_only":
        return "none", single_side
    return single_side, "none"


def tick_label_layout_from_sides(label_side: object, tick_side: object) -> str:
    label = _normalize_tick_side(label_side, field_name="label_side")
    tick = _normalize_tick_side(tick_side, field_name="tick_side")
    if label == "outside" and tick == "inside":
        return "label_out_tick_in"
    if label == "inside" and tick == "outside":
        return "label_in_tick_out"
    if label == "none" and tick != "none":
        return "tick_only"
    if label != "none" and tick == "none":
        return "label_only"
    return DEFAULT_TICK_LABEL_LAYOUT


def _normalize_feature_lane(raw: object) -> str:
    lane = str(raw).strip().lower()
    if lane not in _FEATURE_LANE_VALUES:
        raise ValueError("lane_direction must be one of inside, outside, split")
    return lane


def parse_circular_track_slot(
    raw: str,
    *,
    _allow_legacy_transport: bool = False,
) -> CircularTrackSlot:
    """Parse `<slot_id>:<renderer>@key=value,...` into a slot input object."""

    original, slot_id, renderer_raw, opts = parse_track_slot_text(
        raw,
        mode="circular",
        error_type=CircularTrackSlotParseError,
    )
    renderer = _normalize_renderer(renderer_raw)
    if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
        raise CircularTrackSlotParseError(f"unknown circular track renderer '{renderer}'", original)

    enabled = True
    side: str | None = None
    radius: ScalarSpec | None = None
    width: ScalarSpec | None = None
    inner_gap_px: float | None = None
    outer_gap_px: float | None = None
    legacy_spacing: ScalarSpec | None = None
    z = 0
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
                    if not _allow_legacy_transport:
                        raise ValueError(
                            "'spacing' is no longer supported; use inner_gap_px and outer_gap_px"
                        )
                    legacy_spacing = ScalarSpec.parse(value)
                elif key == "__gbdraw_legacy_spacing":
                    if not _allow_legacy_transport:
                        raise ValueError(
                            "private legacy spacing transport is not accepted by fresh track slots"
                        )
                    legacy_spacing = ScalarSpec.parse(value)
                elif key == "inner_gap_px":
                    inner_gap_px = _parse_gap_px(value, field_name="inner_gap_px")
                elif key == "outer_gap_px":
                    outer_gap_px = _parse_gap_px(value, field_name="outer_gap_px")
                elif key in _OBSOLETE_GEOMETRY_KEYS:
                    raise ValueError(f"'{key}' is no longer supported; use r=<radius> with w=<width>")
                elif key in {"innerradius", "outerradius"}:
                    raise ValueError(f"'{raw_key}' is no longer supported; use r=<radius> with w=<width>")
                elif key in _OBSOLETE_SPACING_KEYS or key == "gapafter":
                    raise ValueError(
                        f"'{key}' is no longer supported; use inner_gap_px and outer_gap_px"
                    )
                elif key == "side":
                    side = _normalize_side_value(value)
                elif key in {"strict", "compress", "reserve"}:
                    if not _allow_legacy_transport:
                        raise ValueError(
                            f"'{key}' is no longer supported; geometry and reservation are derived from side"
                        )
                elif key in {"nt", "dinucleotide"}:
                    params["nt"] = value.upper()
                elif key.startswith("_"):
                    raise ValueError(
                        f"private circular track parameter '{raw_key}' is not supported"
                    )
                else:
                    params[key] = value
        except CircularTrackSlotParseError:
            raise
        except Exception as exc:
            raise CircularTrackSlotParseError(str(exc), original) from exc

    if not slot_id:
        raise CircularTrackSlotParseError("missing circular track slot id", original)
    if renderer == "dinucleotide_skew":
        params = _normalize_dinucleotide_skew_color_params(params)

    slot = CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        enabled=enabled,
        side=side,
        radius=radius,
        width=width,
        z=z,
        params=params,
        inner_gap_px=inner_gap_px,
        outer_gap_px=outer_gap_px,
    )
    if legacy_spacing is not None:
        slot = _internal_circular_track_slot(
            slot,
            legacy_spacing=legacy_spacing,
        )
    try:
        _normalize_circular_track_slots([slot], validate_anchor_references=False)
    except Exception as exc:
        raise CircularTrackSlotParseError(str(exc), original) from exc
    return slot


def parse_circular_track_slots(
    specs: Sequence[str | CircularTrackSlot],
    *,
    _allow_legacy_transport: bool = False,
) -> list[CircularTrackSlot]:
    """Parse and validate circular slot specs."""

    out: list[CircularTrackSlot] = []
    seen: set[str] = set()
    for item in specs:
        if isinstance(item, CircularTrackSlot):
            renderer = _normalize_renderer(str(item.renderer))
            params = {str(key): value for key, value in dict(item.params or {}).items()}
            legacy_spacing: ScalarSpec | None = None
            if _allow_legacy_transport:
                legacy_spacing_raw = params.pop("__gbdraw_legacy_spacing", None)
                if legacy_spacing_raw is not None:
                    legacy_spacing = (
                        legacy_spacing_raw
                        if isinstance(legacy_spacing_raw, ScalarSpec)
                        else ScalarSpec.parse(legacy_spacing_raw)
                    )
            if renderer == "dinucleotide_skew":
                params = _normalize_dinucleotide_skew_color_params(params)
            slot = replace(item, renderer=renderer, params=params)
            if legacy_spacing is not None:
                slot = _internal_circular_track_slot(
                    slot,
                    legacy_spacing=legacy_spacing,
                )
            try:
                _normalize_circular_track_slots([slot], validate_anchor_references=False)
            except Exception as exc:
                raise CircularTrackSlotParseError(str(exc), str(slot.id)) from exc
        else:
            slot = parse_circular_track_slot(
                str(item),
                _allow_legacy_transport=_allow_legacy_transport,
            )
        if slot.id in seen:
            raise CircularTrackSlotParseError("duplicate circular track slot id", slot.id)
        seen.add(slot.id)
        out.append(slot)
    normalize_circular_track_slots(out)
    return out


def _normalized_feature_side_and_params(slot: CircularTrackSlot, params: dict[str, Any]) -> tuple[str, bool, dict[str, Any]]:
    raw_lane = params.get("lane_direction", params.get("lanes"))
    raw_side = slot.side
    reserve = False

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
            reserve = True
    elif side is None:
        if lane == "inside":
            side = "inside"
        elif lane == "outside":
            side = "outside"
        else:
            side = "overlay"
            reserve = True
    else:
        expected = "overlay" if lane == "split" else lane
        if side != expected:
            raise ValueError(
                f"features slot '{slot.id}' has conflicting side={side!r} and lane_direction={lane!r}"
            )
        if side == "overlay":
            reserve = True

    params.pop("lanes", None)
    params["lane_direction"] = lane
    return str(side), reserve, params


def _axis_derived_side(slot_index: int, axis_index: int) -> str:
    return "outside" if int(slot_index) < int(axis_index) else "inside"


def _axis_side_conflict_message(slot: CircularTrackSlot, derived_side: str, explicit_side: str) -> str:
    return (
        f"--circular_track_axis_index places slot '{slot.id}' {derived_side}, "
        f"but the slot specifies side={explicit_side}. Remove side= or move the Axis boundary."
    )


def _slot_with_axis_derived_side(slot: CircularTrackSlot, slot_index: int, axis_index: int) -> CircularTrackSlot:
    """Return a transient slot whose side/lane matches an explicit axis boundary."""

    derived_side = _axis_derived_side(slot_index, axis_index)
    renderer = _normalize_renderer(str(slot.renderer))
    params = {str(key): value for key, value in dict(slot.params or {}).items()}
    explicit_side = _normalize_side_value(slot.side) if slot.side is not None else None

    if renderer == "features":
        raw_lane = params.get("lane_direction", params.get("lanes"))
        lane = _normalize_feature_lane(raw_lane) if raw_lane is not None else None
        if lane == "split":
            if explicit_side is not None and explicit_side != "overlay":
                raise ValueError(
                    f"--circular_track_axis_index places split feature slot '{slot.id}' on the Axis, "
                    f"but the slot specifies side={explicit_side}. Use side=overlay or remove side=."
                )
            params.pop("lanes", None)
            params["lane_direction"] = "split"
            return replace(slot, renderer=renderer, side="overlay", params=params)
        if lane in {"inside", "outside"} and lane != derived_side:
            raise ValueError(
                f"--circular_track_axis_index places feature slot '{slot.id}' {derived_side}, "
                f"but the slot specifies lane_direction={lane}. Remove lane_direction= or move the Axis boundary."
            )
        if explicit_side is not None and explicit_side != derived_side:
            raise ValueError(_axis_side_conflict_message(slot, derived_side, explicit_side))
        params.pop("lanes", None)
        params["lane_direction"] = lane or derived_side
        return replace(slot, renderer=renderer, side=derived_side, params=params)

    if renderer in {"ticks", "annotations"} and explicit_side == "overlay":
        return replace(slot, renderer=renderer, side="overlay", params=params)

    if explicit_side is not None and explicit_side != derived_side:
        raise ValueError(_axis_side_conflict_message(slot, derived_side, explicit_side))
    return replace(slot, renderer=renderer, side=derived_side, params=params)


def circular_track_slots_with_axis_side(
    slots: Sequence[CircularTrackSlot],
    axis_index: int,
) -> list[CircularTrackSlot]:
    """Derive transient slot sides from an explicit circular Axis boundary."""

    if not isinstance(axis_index, int):
        raise ValueError("--circular_track_axis_index must be an integer")
    if axis_index < 0 or axis_index > len(slots):
        raise ValueError(
            f"--circular_track_axis_index must be between 0 and the number of circular track slots ({len(slots)})"
        )
    return [
        _slot_with_axis_derived_side(slot, slot_index, axis_index)
        for slot_index, slot in enumerate(slots)
    ]


def _normalized_tick_params(params: dict[str, Any]) -> dict[str, Any]:
    if "axis" in {str(key).strip().lower() for key in params}:
        raise ValueError("ticks slots no longer accept 'axis'; the circular axis is fixed and not a slot")
    lowered_keys = {str(key).strip().lower() for key in params}
    if "label_side" in lowered_keys or "tick_side" in lowered_keys:
        raise ValueError("ticks slots use tick_label_layout; label_side and tick_side are no longer supported")
    params["tick_label_layout"] = normalize_tick_label_layout(params.get("tick_label_layout"))
    return params


def _normalize_circular_track_slots(
    slots: Sequence[CircularTrackSlot],
    *,
    validate_anchor_references: bool,
) -> list[NormalizedCircularTrackSlot]:
    """Normalize slots, optionally deferring aggregate anchor validation."""

    normalized: list[NormalizedCircularTrackSlot] = []
    seen: set[str] = set()
    for slot_index, slot in enumerate(slots):
        if slot.id in seen:
            raise ValueError(f"duplicate circular track slot id: {slot.id}")
        seen.add(str(slot.id))
        renderer = _normalize_renderer(str(slot.renderer))
        if renderer not in SUPPORTED_CIRCULAR_TRACK_RENDERERS:
            raise ValueError(f"unknown circular track renderer: {slot.renderer}")

        raw_params = {str(key): value for key, value in dict(slot.params or {}).items()}
        internal_slot = (
            slot
            if isinstance(slot, _InternalCircularTrackSlot)
            else None
        )
        legacy_spacing = (
            internal_slot.legacy_spacing
            if internal_slot is not None
            else None
        )
        for raw_key in raw_params:
            key = raw_key.strip()
            normalized_key = key.lower()
            if normalized_key.startswith("_"):
                raise ValueError(
                    f"circular track slot '{slot.id}' contains private "
                    f"parameter '{raw_key}', which is not accepted by fresh slots"
                )
            if (
                normalized_key in _OBSOLETE_CAMEL_KEYS
                or normalized_key in _GENERIC_LAYOUT_KEYS
                or normalized_key in _OBSOLETE_GEOMETRY_KEYS
                or normalized_key in _OBSOLETE_SPACING_KEYS
            ):
                raise ValueError(
                    f"circular track slot '{slot.id}' stores generic layout field '{raw_key}' in params; "
                    "use slot-level radius, width, inner_gap_px, outer_gap_px, side, and z fields"
                )
        if not slot.enabled:
            continue
        params = dict(raw_params)
        side = _normalize_side_value(slot.side) if slot.side is not None else (
            "outside" if renderer == "annotations" else "inside"
        )
        reserve = side == "overlay"
        annotation_params: AnnotationTrackParams | None = None

        if renderer == "features":
            side, reserve, params = _normalized_feature_side_and_params(slot, params)
        elif renderer == "ticks":
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
            reserve = side == "overlay"
            params = _normalized_tick_params(params)
        elif renderer == "sequence_conservation":
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
            if side == "overlay":
                raise ValueError("sequence_conservation slots cannot use side=overlay")
        elif renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
            if renderer == "dinucleotide_skew":
                params = _normalize_dinucleotide_skew_color_params(params)
            elif renderer == "depth" and "track_index" in params:
                raw_track_index = params["track_index"]
                try:
                    params["track_index"] = parse_nonnegative_integer(
                        raw_track_index,
                        field_name=f"depth slot '{slot.id}' track_index",
                    )
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"depth slot '{slot.id}' has invalid track_index={raw_track_index!r}"
                    ) from exc
        elif renderer == "spacer":
            side = _normalize_side_value(slot.side) if slot.side is not None else "inside"
        elif renderer == "annotations":
            side = _normalize_side_value(slot.side) if slot.side is not None else "outside"
            annotation_params = annotation_track_params_from_mapping(params)
            reserve = side != "overlay"
            if side == "overlay" and annotation_params.anchor_slot is None:
                raise ValueError(f"annotation slot '{slot.id}' with side=overlay requires anchor_slot")
            if side != "overlay" and annotation_params.anchor_slot is not None:
                raise ValueError(f"annotation slot '{slot.id}' uses anchor_slot without side=overlay")

        auto_compress = (
            internal_slot.auto_compress
            if internal_slot is not None
            else False
        )
        inner_gap_px = (
            _parse_gap_px(slot.inner_gap_px, field_name="inner_gap_px")
            if slot.inner_gap_px is not None
            else None
        )
        outer_gap_px = (
            _parse_gap_px(slot.outer_gap_px, field_name="outer_gap_px")
            if slot.outer_gap_px is not None
            else None
        )
        compress = (
            renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
            and side == "inside"
            and slot.radius is None
            and (slot.width is None or auto_compress)
        )

        normalized.append(
            NormalizedCircularTrackSlot(
                slot_index=int(slot_index),
                id=str(slot.id),
                renderer=renderer,
                enabled=True,
                side=side,
                radius=slot.radius,
                width=slot.width,
                legacy_spacing=legacy_spacing,
                inner_gap_px=inner_gap_px,
                outer_gap_px=outer_gap_px,
                z=int(slot.z),
                compress=compress,
                reserve=reserve,
                params=params,
                annotation=annotation_params,
                preferred_anchor_radius=(
                    internal_slot.preferred_anchor_radius
                    if internal_slot is not None
                    else None
                ),
            )
        )
    if validate_anchor_references:
        validate_overlay_annotation_anchors(
            normalized,
            anchorless_renderers={"ticks", "spacer"},
        )
    return normalized


def normalize_circular_track_slots(
    slots: Sequence[CircularTrackSlot],
) -> list[NormalizedCircularTrackSlot]:
    """Validate and normalize a complete circular track-slot list."""

    return _normalize_circular_track_slots(slots, validate_anchor_references=True)


def normalize_circular_track_slots_with_axis(
    slots: Sequence[CircularTrackSlot],
    axis_index: int | None = None,
) -> list[NormalizedCircularTrackSlot]:
    """Normalize circular slots, optionally deriving side from an explicit Axis index."""

    if axis_index is None:
        return normalize_circular_track_slots(slots)
    return normalize_circular_track_slots(circular_track_slots_with_axis_side(slots, axis_index))


def default_circular_track_slots(
    *,
    show_features: bool = True,
    show_ticks: bool = True,
    show_depth: bool = False,
    depth_track_count: int = 1,
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
        count = max(1, int(depth_track_count))
        if count == 1:
            slots.append(CircularTrackSlot(id="depth", renderer="depth", params={"track_index": 0}))
        else:
            for index in range(count):
                slots.append(
                    CircularTrackSlot(
                        id=f"depth_{index + 1}",
                        renderer="depth",
                        params={"track_index": index},
                    )
                )
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
    depth_track_count: int = 1,
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
        if renderer == "depth" and int(depth_track_count) > 1:
            for track_index in range(int(depth_track_count)):
                slots.append(
                    CircularTrackSlot(
                        id=f"depth_{track_index + 1}",
                        renderer=renderer,
                        params={"track_index": track_index},
                    )
                )
        else:
            if renderer == "depth":
                params["track_index"] = 0
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
    "CONSERVATION_CIRCULAR_TRACK_RENDERERS",
    "SUPPORTED_CIRCULAR_TRACK_RENDERERS",
    "circular_track_slots_with_axis_side",
    "circular_track_slots_from_order",
    "default_circular_track_slots",
    "normalize_circular_track_slots",
    "normalize_circular_track_slots_with_axis",
    "parse_circular_track_slot",
    "parse_circular_track_slots",
]
