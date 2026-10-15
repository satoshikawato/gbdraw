from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import TYPE_CHECKING, Any, Collection, Literal, Mapping, Protocol, Sequence

from gbdraw.core.numeric import parse_ascii_integer
from gbdraw.io.colors import resolve_color_to_hex
from gbdraw.mode_profiles import validate_dinucleotide

if TYPE_CHECKING:
    from gbdraw.annotations.models import AnnotationTrackParams


@dataclass
class CircularTrackSlotParseError(ValueError):
    message: str
    raw: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.message}: {self.raw!r}"


def parse_track_slot_text(
    raw: str,
    *,
    mode: Literal["circular", "linear"],
    error_type: type[ValueError],
) -> tuple[str, str, str, str]:
    """Parse the mode-neutral raw/head portion of one track-slot token."""

    original = raw
    text = str(raw).strip()
    if not text or text.startswith("#"):
        raise error_type("empty/comment line", original)
    text = strip_inline_comment(text)
    if "@" in text:
        head, options = text.split("@", 1)
        options = options.strip()
    else:
        head, options = text, ""
    if ":" not in head:
        suffix = (
            "; shortcut slot strings are no longer supported"
            if mode == "circular"
            else ""
        )
        raise error_type(
            f"{mode} track slots require '<slot_id>:<renderer>@...'{suffix}",
            original,
        )
    slot_id_raw, renderer_raw = head.split(":", 1)
    slot_id = slot_id_raw.strip()
    renderer = renderer_raw.strip()
    if not slot_id:
        raise error_type(f"missing {mode} track slot id", original)
    if not renderer:
        raise error_type(f"missing {mode} track renderer", original)
    return original, slot_id, renderer, options


def parse_bool(raw: object) -> bool:
    s = str(raw).strip().lower()
    if s in {"1", "true", "yes", "on"}:
        return True
    if s in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"invalid bool: {raw!r}")


def parse_nonnegative_integer(raw: object, *, field_name: str) -> int:
    """Parse an integer identity without silently truncating numeric values."""

    if isinstance(raw, bool):
        raise ValueError(f"{field_name} must be a nonnegative integer")
    if isinstance(raw, Integral):
        value = int(raw)
    else:
        parsed = parse_ascii_integer(raw.strip()) if isinstance(raw, str) else None
        if parsed is None:
            raise ValueError(f"{field_name} must be a nonnegative integer")
        value = parsed
    if value < 0:
        raise ValueError(f"{field_name} must be >= 0")
    return value


# An overlay annotation slot refers to other slots by id: parameter name, and its role in messages.
_OVERLAY_SLOT_REFERENCES = (("anchor_slot", "anchor"), ("to_slot", "to_slot"))


def validate_annotation_slot_side(slot_id: str, side: str, annotation: AnnotationTrackParams) -> None:
    """Require an anchor on an overlay slot, and no slot reference on any other side."""

    if side == "overlay":
        if annotation.anchor_slot is None:
            raise ValueError(f"annotation slot '{slot_id}' with side=overlay requires anchor_slot")
        return
    for name, _role in _OVERLAY_SLOT_REFERENCES:
        if getattr(annotation, name) is not None:
            raise ValueError(f"annotation slot '{slot_id}' uses {name} without side=overlay")


def validate_overlay_annotation_anchors(
    slots: Sequence[Any],
    *,
    anchorless_renderers: Collection[str],
) -> None:
    """Validate aggregate overlay references after every slot is normalized."""

    by_id = {str(slot.id): slot for slot in slots}
    for slot in slots:
        annotation = getattr(slot, "annotation", None)
        if slot.renderer != "annotations" or slot.side != "overlay" or annotation is None:
            continue
        if annotation.to_slot is not None and annotation.to_slot == annotation.anchor_slot:
            raise ValueError(
                f"annotation slot '{slot.id}' to_slot must name a slot other than "
                f"its anchor_slot '{annotation.anchor_slot}'"
            )
        for name, role in _OVERLAY_SLOT_REFERENCES:
            reference = getattr(annotation, name)
            if reference is None:
                continue
            target = by_id.get(str(reference))
            if target is None:
                raise ValueError(
                    f"annotation slot '{slot.id}' references unknown {name}={reference!r}"
                )
            if target.renderer in anchorless_renderers:
                raise ValueError(
                    f"annotation slot '{slot.id}' {role} '{target.id}' has no drawable band"
                )
            if target.renderer == "annotations":
                raise ValueError(
                    f"annotation slot '{slot.id}' cannot use annotation slot "
                    f"'{target.id}' as {role}"
                )
            if annotation.layer == "underlay" and slot.z >= target.z:
                raise ValueError(
                    f"annotation underlay slot '{slot.id}' must have z less than {role} '{target.id}'"
                )
            if annotation.layer == "foreground" and slot.z <= target.z:
                raise ValueError(
                    f"annotation foreground slot '{slot.id}' must have z greater than {role} '{target.id}'"
                )


def split_kv_list(raw: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"expected key=value, got: {part!r}")
        key, value = part.split("=", 1)
        out.append((key.strip(), value.strip()))
    return out


def normalize_dinucleotide_skew_color_params(params: dict) -> dict:
    if "positive_color" not in params and "high_color" in params:
        params["positive_color"] = params["high_color"]
    if "negative_color" not in params and "low_color" in params:
        params["negative_color"] = params["low_color"]
    params.pop("high_color", None)
    params.pop("low_color", None)
    return params


class SlotWithParams(Protocol):
    """A track slot as its params readers see it (any slot type of either mode)."""

    @property
    def params(self) -> Mapping[str, object]: ...


class SlotView(SlotWithParams, Protocol):
    """A track slot as the legend reads it: identity, renderer, params, and whether it draws."""

    @property
    def id(self) -> str: ...

    @property
    def renderer(self) -> str: ...

    @property
    def enabled(self) -> bool: ...


def _slot_param_text(slot: SlotWithParams, *names: str) -> str | None:
    params = slot.params
    for name in names:
        value = params.get(name)
        text = str(value).strip() if value is not None else ""
        if text:
            return text
    return None


def slot_skew_colors(slot: SlotWithParams) -> tuple[str | None, str | None]:
    """The high and low fill colors a skew slot sets in its params, as hex, or ``None``."""

    positive_color = _slot_param_text(slot, "positive_color", "high_color")
    negative_color = _slot_param_text(slot, "negative_color", "low_color")
    return (
        resolve_color_to_hex(positive_color) if positive_color else None,
        resolve_color_to_hex(negative_color) if negative_color else None,
    )


def slot_dinucleotide(params: Mapping[str, Any] | None, default: str) -> str:
    """Return a numeric slot's validated pair; an unset ``nt`` uses ``default``."""

    raw = (params or {}).get("nt", (params or {}).get("dinucleotide"))
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        raw = default
    return validate_dinucleotide(raw, field_name="nt")


def normalize_dinucleotide_params(params: dict) -> dict:
    """Store an explicit slot pair once, as validated ``nt``."""

    if "dinucleotide" in params and "nt" not in params:
        params["nt"] = params.pop("dinucleotide")
    if "nt" in params:
        params["nt"] = validate_dinucleotide(params["nt"], field_name="nt")
    return params


def strip_inline_comment(raw: str) -> str:
    for marker in (" #", "\t#"):
        index = raw.find(marker)
        if index >= 0:
            return raw[:index].strip()
    return raw


__all__ = [
    "CircularTrackSlotParseError",
    "SlotView",
    "SlotWithParams",
    "normalize_dinucleotide_params",
    "parse_bool",
    "parse_nonnegative_integer",
    "slot_skew_colors",
    "parse_track_slot_text",
    "slot_dinucleotide",
    "split_kv_list",
    "validate_annotation_slot_side",
    "validate_overlay_annotation_anchors",
]
