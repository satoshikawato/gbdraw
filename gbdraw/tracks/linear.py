from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Literal, Mapping, Sequence

from .parsing import parse_bool, split_kv_list
from .scalars import ScalarSpec


LinearTrackRendererName = Literal[
    "features",
    "dinucleotide_content",
    "dinucleotide_skew",
    "depth",
    "spacer",
]

LinearTrackSide = Literal["above", "below", "overlay"]

SUPPORTED_LINEAR_TRACK_RENDERERS: frozenset[str] = frozenset(
    {
        "features",
        "dinucleotide_content",
        "dinucleotide_skew",
        "depth",
        "spacer",
    }
)

NUMERIC_LINEAR_TRACK_RENDERERS: frozenset[str] = frozenset(
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
    "depth": "depth",
    "gc_content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
}

_GENERIC_LAYOUT_KEYS = {
    "id",
    "renderer",
    "type",
    "side",
    "h",
    "height",
    "spacing",
    "z",
    "z_index",
    "zindex",
    "enabled",
    "show",
    "visible",
}

_SIDE_VALUES = {"above", "below", "overlay"}


@dataclass
class LinearTrackSlotParseError(ValueError):
    message: str
    raw: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.message}: {self.raw!r}"


@dataclass(frozen=True)
class LinearTrackSlot:
    """One public vertical slot input for a linear diagram."""

    id: str
    renderer: LinearTrackRendererName | str
    enabled: bool = True
    side: str | None = None
    height: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z: int = 0
    params: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class NormalizedLinearTrackSlot:
    slot_index: int
    id: str
    renderer: str
    enabled: bool
    side: str
    height: ScalarSpec | None
    spacing: ScalarSpec | None
    z: int
    reserve: bool
    params: Mapping[str, Any]


def _normalize_renderer(raw: str) -> str:
    renderer = str(raw).strip().lower()
    return _RENDERER_ALIASES.get(renderer, renderer)


def _normalize_side_value(raw: object, *, field_name: str = "side") -> str:
    side = str(raw).strip().lower()
    if side not in _SIDE_VALUES:
        raise ValueError(f"{field_name} must be one of above, below, overlay")
    return side


def _parse_px_scalar(raw: object, *, field_name: str) -> ScalarSpec:
    text = str(raw).strip()
    if not text:
        raise ValueError(f"{field_name} cannot be empty")
    if text.endswith("%"):
        raise ValueError(f"{field_name} only accepts px or unitless px values")
    try:
        if text.endswith("px"):
            value = float(text[:-2])
        else:
            value = float(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be a number of pixels") from exc
    return ScalarSpec(value=value, unit="px")


def _validate_px_scalar(spec: ScalarSpec | None, *, field_name: str, allow_zero: bool) -> None:
    if spec is None:
        return
    if spec.unit != "px":
        raise ValueError(f"{field_name} only accepts px or unitless px values")
    value = float(spec.value)
    if value < 0 or (value == 0 and not allow_zero):
        relation = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{field_name} must be {relation}")


def _parse_slot_head(head: str, original: str) -> tuple[str, str]:
    if ":" not in head:
        raise LinearTrackSlotParseError(
            "linear track slots require '<slot_id>:<renderer>@...'",
            original,
        )
    slot_id_raw, renderer_raw = head.split(":", 1)
    slot_id = slot_id_raw.strip()
    renderer = renderer_raw.strip()
    if not slot_id:
        raise LinearTrackSlotParseError("missing linear track slot id", original)
    if not renderer:
        raise LinearTrackSlotParseError("missing linear track renderer", original)
    return slot_id, renderer


def parse_linear_track_slot(raw: str) -> LinearTrackSlot:
    """Parse `<slot_id>:<renderer>@key=value,...` into a slot input object."""

    original = raw
    s = str(raw).strip()
    if not s or s.startswith("#"):
        raise LinearTrackSlotParseError("empty/comment line", original)
    if "#" in s:
        s = s.split("#", 1)[0].strip()

    if "@" in s:
        head, opts = s.split("@", 1)
        opts = opts.strip()
    else:
        head, opts = s, ""

    slot_id, renderer_raw = _parse_slot_head(head.strip(), original)
    renderer = _normalize_renderer(renderer_raw)
    if renderer not in SUPPORTED_LINEAR_TRACK_RENDERERS:
        raise LinearTrackSlotParseError(f"unknown linear track renderer '{renderer}'", original)

    enabled = True
    side: str | None = None
    height: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z = 0
    params: dict[str, Any] = {}

    if opts:
        try:
            for raw_key, raw_value in split_kv_list(opts):
                key = raw_key.strip().lower()
                value = raw_value.strip()
                if key == "id":
                    slot_id = value
                elif key in {"renderer", "type"}:
                    renderer = _normalize_renderer(value)
                    if renderer not in SUPPORTED_LINEAR_TRACK_RENDERERS:
                        raise ValueError(f"unknown linear track renderer '{renderer}'")
                elif key in {"enabled", "show", "visible"}:
                    enabled = parse_bool(value)
                elif key in {"z", "z_index", "zindex"}:
                    z = int(value)
                elif key == "side":
                    side = _normalize_side_value(value)
                elif key in {"h", "height"}:
                    height = _parse_px_scalar(value, field_name="height")
                elif key == "spacing":
                    spacing = _parse_px_scalar(value, field_name="spacing")
                elif key in {"nt", "dinucleotide"}:
                    params["nt"] = value.upper()
                elif key == "track_index":
                    params["track_index"] = int(value)
                else:
                    params[key] = value
        except LinearTrackSlotParseError:
            raise
        except Exception as exc:
            raise LinearTrackSlotParseError(str(exc), original) from exc

    if not slot_id:
        raise LinearTrackSlotParseError("missing linear track slot id", original)

    slot = LinearTrackSlot(
        id=slot_id,
        renderer=renderer,
        enabled=enabled,
        side=side,
        height=height,
        spacing=spacing,
        z=z,
        params=params,
    )
    try:
        normalize_linear_track_slots([slot])
    except Exception as exc:
        raise LinearTrackSlotParseError(str(exc), original) from exc
    return slot


def parse_linear_track_slots(specs: Sequence[str | LinearTrackSlot]) -> list[LinearTrackSlot]:
    """Parse and validate linear slot specs."""

    if len(specs) == 0:
        raise LinearTrackSlotParseError("linear track slot list cannot be empty", "")
    out: list[LinearTrackSlot] = []
    seen: set[str] = set()
    for item in specs:
        if isinstance(item, LinearTrackSlot):
            renderer = _normalize_renderer(str(item.renderer))
            slot = replace(item, renderer=renderer) if renderer != str(item.renderer) else item
            try:
                normalize_linear_track_slots([slot])
            except Exception as exc:
                raise LinearTrackSlotParseError(str(exc), str(slot.id)) from exc
        else:
            slot = parse_linear_track_slot(str(item))
        if slot.id in seen:
            raise LinearTrackSlotParseError("duplicate linear track slot id", slot.id)
        seen.add(slot.id)
        out.append(slot)
    normalize_linear_track_slots(out)
    return out


def _axis_derived_side(slot_index: int, axis_index: int) -> str:
    return "above" if int(slot_index) < int(axis_index) else "below"


def _axis_side_conflict_message(
    *,
    slot_id: str,
    slot_index: int,
    explicit_side: str,
    derived_side: str,
    axis_index: int,
) -> str:
    return (
        f"linear track slot '{slot_id}' at index {slot_index} has side={explicit_side!r}, "
        f"which conflicts with --linear_track_axis_index {axis_index} "
        f"(axis-derived side is {derived_side!r})"
    )


def linear_track_slots_with_axis_side(
    slots: Sequence[LinearTrackSlot],
    axis_index: int,
) -> list[LinearTrackSlot]:
    """Derive transient slot sides from an explicit linear Axis boundary."""

    if not isinstance(axis_index, int):
        raise ValueError("--linear_track_axis_index must be an integer")
    if axis_index < 0 or axis_index > len(slots):
        raise ValueError(
            f"--linear_track_axis_index must be between 0 and the number of linear track slots ({len(slots)})"
        )
    out: list[LinearTrackSlot] = []
    for slot_index, slot in enumerate(slots):
        renderer = _normalize_renderer(str(slot.renderer))
        params = {str(key): value for key, value in dict(slot.params or {}).items()}
        derived_side = _axis_derived_side(slot_index, axis_index)
        explicit_side = _normalize_side_value(slot.side) if slot.side is not None else None
        if explicit_side == "overlay":
            if renderer != "features":
                raise ValueError("side=overlay is only supported for features slots")
            if slot_index != axis_index:
                raise ValueError(
                    _axis_side_conflict_message(
                        slot_id=str(slot.id),
                        slot_index=slot_index,
                        explicit_side=explicit_side,
                        derived_side=derived_side,
                        axis_index=axis_index,
                    )
                )
        elif explicit_side is not None and explicit_side != derived_side:
            raise ValueError(
                _axis_side_conflict_message(
                    slot_id=str(slot.id),
                    slot_index=slot_index,
                    explicit_side=explicit_side,
                    derived_side=derived_side,
                    axis_index=axis_index,
                )
            )
        out.append(
            replace(
                slot,
                renderer=renderer,
                side=explicit_side or derived_side,
                params=params,
            )
        )
    return out


def normalize_linear_track_slots(slots: Sequence[LinearTrackSlot]) -> list[NormalizedLinearTrackSlot]:
    """Validate and normalize slots for the linear vertical resolver."""

    if len(slots) == 0:
        raise ValueError("linear track slot list cannot be empty")

    normalized: list[NormalizedLinearTrackSlot] = []
    seen: set[str] = set()
    feature_slot_count = 0
    for slot_index, slot in enumerate(slots):
        if slot.id in seen:
            raise ValueError(f"duplicate linear track slot id: {slot.id}")
        seen.add(str(slot.id))
        renderer = _normalize_renderer(str(slot.renderer))
        if renderer not in SUPPORTED_LINEAR_TRACK_RENDERERS:
            raise ValueError(f"unknown linear track renderer: {slot.renderer}")

        raw_params = {str(key): value for key, value in dict(slot.params or {}).items()}
        for raw_key in raw_params:
            normalized_key = raw_key.strip().lower()
            if normalized_key in _GENERIC_LAYOUT_KEYS:
                raise ValueError(
                    f"linear track slot '{slot.id}' stores generic layout field '{raw_key}' in params; "
                    "use slot-level height, spacing, side, and z fields"
                )
        params = dict(raw_params)

        _validate_px_scalar(slot.height, field_name="height", allow_zero=False)
        _validate_px_scalar(slot.spacing, field_name="spacing", allow_zero=True)

        if not slot.enabled:
            continue

        side = _normalize_side_value(slot.side) if slot.side is not None else "below"
        if side == "overlay" and renderer != "features":
            raise ValueError("side=overlay is only supported for features slots")
        if renderer == "features":
            feature_slot_count += 1
            if feature_slot_count > 1:
                raise ValueError("linear custom track slots support only one features slot")

        if renderer in {"dinucleotide_content", "dinucleotide_skew"}:
            if "dinucleotide" in params and "nt" not in params:
                params["nt"] = str(params.pop("dinucleotide")).upper()
            elif "nt" in params:
                params["nt"] = str(params["nt"]).upper()
        elif renderer == "depth":
            raw_track_index = params.get("track_index", 0)
            try:
                track_index = int(raw_track_index)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"depth slot '{slot.id}' has invalid track_index={raw_track_index!r}"
                ) from exc
            if track_index < 0:
                raise ValueError(f"depth slot '{slot.id}' track_index must be >= 0")
            params["track_index"] = track_index
        elif renderer in {"features", "spacer"}:
            pass

        normalized.append(
            NormalizedLinearTrackSlot(
                slot_index=int(slot_index),
                id=str(slot.id),
                renderer=renderer,
                enabled=True,
                side=side,
                height=slot.height,
                spacing=slot.spacing,
                z=int(slot.z),
                reserve=(renderer != "spacer"),
                params=params,
            )
        )
    return normalized


def normalize_linear_track_slots_with_axis(
    slots: Sequence[LinearTrackSlot],
    axis_index: int | None = None,
) -> list[NormalizedLinearTrackSlot]:
    """Normalize linear slots, optionally deriving side from an explicit Axis index."""

    if axis_index is None:
        return normalize_linear_track_slots(slots)
    return normalize_linear_track_slots(linear_track_slots_with_axis_side(slots, axis_index))


def _feature_side_for_track_layout(track_layout: str) -> str:
    normalized = str(track_layout or "middle").strip().lower()
    if normalized == "above":
        return "above"
    if normalized == "below":
        return "below"
    return "overlay"


def default_linear_track_slots(
    *,
    show_features: bool = True,
    show_depth: bool = False,
    depth_track_count: int = 1,
    show_gc: bool = True,
    show_skew: bool = True,
    dinucleotide: str = "GC",
    track_layout: str = "middle",
) -> list[LinearTrackSlot]:
    """Return the default linear slot input list."""

    slots: list[LinearTrackSlot] = []
    nt = str(dinucleotide or "GC").upper()
    if show_features:
        slots.append(
            LinearTrackSlot(
                id="features",
                renderer="features",
                side=_feature_side_for_track_layout(track_layout),
            )
        )
    if show_depth:
        count = max(1, int(depth_track_count))
        if count == 1:
            slots.append(LinearTrackSlot(id="depth", renderer="depth", side="below", params={"track_index": 0}))
        else:
            for index in range(count):
                slots.append(
                    LinearTrackSlot(
                        id=f"depth_{index + 1}",
                        renderer="depth",
                        side="below",
                        params={"track_index": index},
                    )
                )
    if show_gc:
        slots.append(
            LinearTrackSlot(
                id="gc_content",
                renderer="dinucleotide_content",
                side="below",
                params={"nt": nt},
            )
        )
    if show_skew:
        slots.append(
            LinearTrackSlot(
                id="gc_skew",
                renderer="dinucleotide_skew",
                side="below",
                params={"nt": nt},
            )
        )
    return slots


def linear_track_slots_from_order(
    order: str | Sequence[str],
    *,
    show_features: bool = True,
    show_depth: bool = False,
    depth_track_count: int = 1,
    show_gc: bool = True,
    show_skew: bool = True,
    dinucleotide: str = "GC",
    track_layout: str = "middle",
) -> list[LinearTrackSlot]:
    """Expand a comma-separated slot order into explicit slot inputs."""

    enabled = {
        "features": show_features,
        "depth": show_depth,
        "gc_content": show_gc,
        "gc_skew": show_skew,
    }
    nt = str(dinucleotide or "GC").upper()
    items = order.split(",") if isinstance(order, str) else list(order)
    slots: list[LinearTrackSlot] = []
    seen: set[str] = set()
    for raw in items:
        slot_id = str(raw).strip()
        if not slot_id:
            continue
        if slot_id in seen:
            raise ValueError(f"duplicate linear slot id in order: {slot_id}")
        if slot_id not in _ORDER_RENDERER_BY_ID:
            raise ValueError(f"unknown linear slot id in order: {slot_id}")
        seen.add(slot_id)
        if not enabled.get(slot_id, True):
            continue
        renderer = _ORDER_RENDERER_BY_ID[slot_id]
        params: dict[str, Any] = {}
        side = "below"
        if renderer == "features":
            side = _feature_side_for_track_layout(track_layout)
        if renderer in {"dinucleotide_content", "dinucleotide_skew"}:
            params["nt"] = nt
        if renderer == "depth" and int(depth_track_count) > 1:
            for track_index in range(int(depth_track_count)):
                slots.append(
                    LinearTrackSlot(
                        id=f"depth_{track_index + 1}",
                        renderer=renderer,
                        side=side,
                        params={"track_index": track_index},
                    )
                )
        else:
            if renderer == "depth":
                params["track_index"] = 0
            slots.append(
                LinearTrackSlot(
                    id=slot_id,
                    renderer=renderer,
                    side=side,
                    params=params,
                )
            )
    if not slots:
        raise ValueError("linear track order did not enable any slots")
    return slots


__all__ = [
    "LinearTrackRendererName",
    "LinearTrackSide",
    "LinearTrackSlot",
    "LinearTrackSlotParseError",
    "NormalizedLinearTrackSlot",
    "NUMERIC_LINEAR_TRACK_RENDERERS",
    "SUPPORTED_LINEAR_TRACK_RENDERERS",
    "default_linear_track_slots",
    "linear_track_slots_from_order",
    "linear_track_slots_with_axis_side",
    "normalize_linear_track_slots",
    "normalize_linear_track_slots_with_axis",
    "parse_linear_track_slot",
    "parse_linear_track_slots",
]
