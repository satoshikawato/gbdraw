from __future__ import annotations

import re
from dataclasses import dataclass
from numbers import Integral
from typing import Any, Collection, Sequence


@dataclass
class CircularTrackSlotParseError(ValueError):
    message: str
    raw: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.message}: {self.raw!r}"


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
    elif isinstance(raw, str) and re.fullmatch(r"[+-]?\d+", raw.strip()):
        value = int(raw.strip())
    else:
        raise ValueError(f"{field_name} must be a nonnegative integer")
    if value < 0:
        raise ValueError(f"{field_name} must be >= 0")
    return value


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
        anchor = by_id.get(str(annotation.anchor_slot))
        if anchor is None:
            raise ValueError(
                f"annotation slot '{slot.id}' references unknown "
                f"anchor_slot={annotation.anchor_slot!r}"
            )
        if anchor.renderer in anchorless_renderers:
            raise ValueError(
                f"annotation slot '{slot.id}' anchor '{anchor.id}' has no drawable band"
            )
        if anchor.renderer == "annotations":
            raise ValueError(
                f"annotation slot '{slot.id}' cannot use annotation slot "
                f"'{anchor.id}' as anchor"
            )
        if annotation.layer == "underlay" and slot.z >= anchor.z:
            raise ValueError(
                f"annotation underlay slot '{slot.id}' must have z less than anchor '{anchor.id}'"
            )
        if annotation.layer == "foreground" and slot.z <= anchor.z:
            raise ValueError(
                f"annotation foreground slot '{slot.id}' must have z greater than anchor '{anchor.id}'"
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


def strip_inline_comment(raw: str) -> str:
    for marker in (" #", "\t#"):
        index = raw.find(marker)
        if index >= 0:
            return raw[:index].strip()
    return raw


__all__ = [
    "CircularTrackSlotParseError",
    "parse_bool",
    "parse_nonnegative_integer",
    "split_kv_list",
    "validate_overlay_annotation_anchors",
]
