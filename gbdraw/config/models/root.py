from __future__ import annotations

import math
from dataclasses import dataclass, fields, is_dataclass
from numbers import Real
from typing import TYPE_CHECKING, Any, Callable, Iterator, Mapping

from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import USER_COLOR_FORMS, is_user_color

from .canvas import CanvasConfig
from .labels import LabelsConfig
from .objects import ObjectsConfig

if TYPE_CHECKING:
    from _typeshed import DataclassInstance


# D-27 (PD-OI-081): only values whose meaning SVG/CSS fixes are constrained.
# A path segment ending in one of the names selects the domain; ``tick_width``
# is the SVG stroke-width of the Circular tick paths.
_STYLE_LEAF_DOMAINS = (
    (("font_size",), "POSITIVE", "font sizes must be finite numbers greater than zero"),
    (("stroke_width", "tick_width"), "NONNEGATIVE", "stroke widths must be finite numbers of zero or greater"),
)
# D-04: None keeps the automatic interval. A request (a Session, the Web) that
# carries <= 0 is read as None before it reaches here (session_request_codec.py).
_EXACT_LEAF_DOMAINS = {
    "objects.scale.interval": ("POSITIVE_INTEGER_OR_AUTO", "expected a positive integer or None"),
}
# P10c (D-13, D-14): text leaves written into SVG attributes. A color leaf is a
# ``fill``, ``stroke`` or ``*_color`` field (``is_user_color``, any letter case);
# the keyword leaves below match their SVG keywords exactly, as svgwrite did,
# because CairoSVG compares them case-sensitively.
_FONT_WEIGHT_KEYWORDS = frozenset({"normal", "bold", "bolder", "lighter", "inherit"})
_TEXT_ANCHOR_KEYWORDS = frozenset({"start", "middle", "end", "inherit"})
_DOMINANT_BASELINE_KEYWORDS = frozenset(
    "auto use-script no-change reset-size ideographic alphabetic hanging mathematical "
    "central middle text-after-edge text-before-edge text-top text-bottom inherit".split()
)


def _is_font_weight(text: str) -> bool:
    if text in _FONT_WEIGHT_KEYWORDS:
        return True
    # Digits only, the form CairoSVG reads; "700.0" and "7e2" errored under svgwrite.
    return text.isascii() and text.isdigit() and 1 <= int(text) <= 1000


_KEYWORD_LEAF_DOMAINS: dict[str, tuple[str, Callable[[str], bool]]] = {
    "font_weight": ("expected normal, bold, bolder, lighter, inherit, or a number from 1 to 1000", _is_font_weight),
    "text_anchor": ("expected start, middle, end, or inherit", _TEXT_ANCHOR_KEYWORDS.__contains__),
    "dominant_baseline": (
        "expected an SVG dominant-baseline keyword such as auto, central, middle, or hanging",
        _DOMINANT_BASELINE_KEYWORDS.__contains__,
    ),
    "font_family": ("expected a non-empty font family list", lambda text: bool(text.strip())),
}


def style_leaf_domain(path: str) -> tuple[str, str] | None:
    """Return ``(reason, requirement)`` for a constrained SVG style leaf.

    Numeric leaves: font sizes, stroke widths, the scale interval. Text leaves:
    colors (``COLOR``) and the SVG keywords font-weight, text-anchor,
    dominant-baseline, and a non-empty font-family (``KEYWORD``).
    """

    if path in _EXACT_LEAF_DOMAINS:
        return _EXACT_LEAF_DOMAINS[path]
    for suffixes, reason, requirement in _STYLE_LEAF_DOMAINS:
        if any(part.endswith(suffixes) for part in path.split(".")):
            return reason, requirement
    leaf = path.rsplit(".", 1)[-1]
    if leaf in {"fill", "stroke"} or leaf.endswith("_color"):
        return "COLOR", f"colors must be {USER_COLOR_FORMS}"
    if leaf in _KEYWORD_LEAF_DOMAINS:
        return "KEYWORD", _KEYWORD_LEAF_DOMAINS[leaf][0]
    return None


def validate_style_leaf(path: str, value: object, *, prefix: str) -> None:
    """Reject a style leaf outside its SVG/CSS domain (see ``style_leaf_domain``)."""

    domain = style_leaf_domain(path)
    if domain is None or value is None or isinstance(value, bool):
        return
    reason, requirement = domain
    if reason in {"COLOR", "KEYWORD"}:
        if isinstance(value, str) and (
            is_user_color(value)
            if reason == "COLOR"
            else _KEYWORD_LEAF_DOMAINS[path.rsplit(".", 1)[-1]][1](value)
        ):
            return
        # COLOR is in the Web's reason vocabulary; a keyword error names its path only.
        diagnostic: dict[str, object] = {"code": "INPUT_INVALID", "configPath": path}
        if reason == "COLOR":
            diagnostic["reason"] = reason
        raise ValidationError(f"{prefix} {path!r}: {value!r}; {requirement}.", diagnostic=diagnostic)
    if not isinstance(value, Real):
        return
    number = float(value)
    if math.isfinite(number) and (number >= 0 if reason == "NONNEGATIVE" else number > 0):
        return
    raise ValidationError(
        f"{prefix} {path!r}; {requirement}.",
        diagnostic={"code": "INPUT_INVALID", "reason": reason, "configPath": path},
    )


def _typed_leaves(value: DataclassInstance, prefix: str = "") -> Iterator[tuple[str, object]]:
    for config_field in fields(value):
        child = getattr(value, config_field.name)
        path = f"{prefix}{config_field.name}"
        if is_dataclass(child) and not isinstance(child, type):
            yield from _typed_leaves(child, f"{path}.")
        else:
            yield path, child


@dataclass(frozen=True)
class GbdrawConfig:
    canvas: CanvasConfig
    labels: LabelsConfig
    objects: ObjectsConfig

    def __post_init__(self) -> None:
        for path, value in _typed_leaves(self):
            validate_style_leaf(path, value, prefix="Invalid configuration value at")

    @classmethod
    def from_dict(cls, config_dict: Mapping[str, Any]) -> "GbdrawConfig":
        return cls(
            canvas=CanvasConfig.from_dict(config_dict["canvas"]),
            labels=LabelsConfig.from_dict(config_dict["labels"]),
            objects=ObjectsConfig.from_dict(config_dict["objects"]),
        )


__all__ = ["GbdrawConfig", "style_leaf_domain", "validate_style_leaf"]
