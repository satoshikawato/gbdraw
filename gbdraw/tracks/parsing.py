from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from numbers import Integral
from typing import TYPE_CHECKING, Any, Collection, Iterable, Literal, Mapping, Protocol, Sequence

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


def depth_slot_track_index(slot_id: str, params: Mapping[str, object]) -> int:
    """The logical depth track a depth slot draws: its ``track_index`` param, 0 when it has none."""

    return parse_nonnegative_integer(params.get("track_index", 0), field_name=f"depth slot '{slot_id}' track_index")


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
    listed_ids: Collection[str] | None = None,
) -> None:
    """Validate aggregate overlay references after every slot is normalized.

    ``slots`` is the drawing's list unless ``listed_ids`` is given: then it is a
    list the drawing may still add slots to (Circular comparison rings), and a
    reference to an id not in ``listed_ids`` (every listed id, enabled or not)
    is left to the check of the drawing's complete list.
    """

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
            if target is None and listed_ids is not None and str(reference) not in listed_ids:
                continue
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


# The params keys each numeric or spacer renderer reads (parser aliases
# included), the same in both modes; each mode adds its own renderers
# (`_*_SLOT_PARAM_KEYS` in tracks/circular.py and tracks/linear.py). Annotation
# slots check their own keys (`annotation_track_params_from_mapping`).
LEGEND_LABEL_PARAM_KEYS = ("legend_label", "label")
_DINUCLEOTIDE_PARAM_KEYS = frozenset({"nt", "dinucleotide"})
SHARED_SLOT_PARAM_KEYS: Mapping[str, frozenset[str]] = MappingProxyType({
    "dinucleotide_content": _DINUCLEOTIDE_PARAM_KEYS | set(LEGEND_LABEL_PARAM_KEYS),
    "dinucleotide_skew": _DINUCLEOTIDE_PARAM_KEYS | set(LEGEND_LABEL_PARAM_KEYS) | {
        "positive_color", "negative_color", "high_color", "low_color"
    },
    "depth": frozenset({"track_index", *LEGEND_LABEL_PARAM_KEYS}),
    "spacer": frozenset(),
})
# Keys every non-annotation slot takes and ignores where its renderer does not
# read them (Owner LP-D70): a legend label left on a row whose renderer changed
# in a released Web, and a comparison ring's label, color and fileName.
CARRIED_SLOT_PARAM_KEYS = frozenset({*LEGEND_LABEL_PARAM_KEYS, "color", "fileName"})


def unread_slot_param_keys(
    read: Mapping[str, frozenset[str]],
    carried: Mapping[str, frozenset[str]] | None = None,
) -> Mapping[str, tuple[str, ...]]:
    """The params keys each renderer of ``read`` takes without reading them: the keys every slot carries
    (``CARRIED_SLOT_PARAM_KEYS``) and the renderer's own ``carried`` ones, less the keys it reads. They do not change
    the drawing, so a Source recipe may leave them out (OV-3009)."""

    return MappingProxyType({
        renderer: tuple(sorted((CARRIED_SLOT_PARAM_KEYS | (carried or {}).get(renderer, frozenset())) - keys))
        for renderer, keys in read.items()
    })


def reject_unread_slot_params(
    params: Mapping[str, object],
    *,
    mode: Literal["circular", "linear"],
    slot_id: str,
    renderer: str,
    accepted: Mapping[str, frozenset[str]],
) -> None:
    """Raise for a params key the slot's renderer neither reads nor carries (OV-477, OV-1186)."""

    read = accepted[renderer]
    unknown = sorted(str(key) for key in params if str(key) not in read and str(key) not in CARRIED_SLOT_PARAM_KEYS)
    if not unknown:
        return
    references = [name for name, _role in _OVERLAY_SLOT_REFERENCES if name in unknown]
    if references:
        raise ValueError(
            f"{mode} track slot '{slot_id}' ({renderer}) cannot use {', '.join(sorted(references))}: "
            "only annotation slots take anchor_slot and to_slot"
        )
    raise ValueError(
        f"{mode} track slot '{slot_id}' ({renderer}) has unknown parameter(s): {', '.join(unknown)}; "
        f"it takes {', '.join(sorted(read)) or 'none'}"
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


def slots_draw_features(slots: Iterable[SlotView] | None) -> bool:
    """Whether a drawing on ``slots`` draws its features, in either mode.

    The mode's default track stack (``None``) draws them; a custom stack draws
    them only with an enabled ``features`` slot. The drawing, its Legend and its
    input summary read this one answer.
    """

    return slots is None or any(slot.enabled and str(slot.renderer) == "features" for slot in slots)


def _slot_param_text(slot: SlotWithParams, *names: str) -> str | None:
    params = slot.params
    for name in names:
        value = params.get(name)
        text = str(value).strip() if value is not None else ""
        if text:
            return text
    return None


def slot_legend_label(slot: SlotWithParams, fallback: str) -> str:
    """The legend row name a slot sets (``legend_label``, else ``label``), or ``fallback``."""

    params = slot.params
    raw_label = params.get(LEGEND_LABEL_PARAM_KEYS[0], params.get(LEGEND_LABEL_PARAM_KEYS[1]))
    label = str(raw_label).strip() if raw_label is not None else ""
    return label or fallback


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
    "CARRIED_SLOT_PARAM_KEYS",
    "CircularTrackSlotParseError",
    "LEGEND_LABEL_PARAM_KEYS",
    "SHARED_SLOT_PARAM_KEYS",
    "SlotView",
    "SlotWithParams",
    "depth_slot_track_index",
    "normalize_dinucleotide_params",
    "parse_bool",
    "parse_nonnegative_integer",
    "reject_unread_slot_params",
    "slot_legend_label",
    "slot_skew_colors",
    "slots_draw_features",
    "parse_track_slot_text",
    "unread_slot_param_keys",
    "slot_dinucleotide",
    "split_kv_list",
    "validate_annotation_slot_side",
    "validate_overlay_annotation_anchors",
]
