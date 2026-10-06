from __future__ import annotations

import math
from dataclasses import dataclass, fields, is_dataclass
from numbers import Real
from typing import TYPE_CHECKING, Any, Iterator, Mapping

from gbdraw.exceptions import ValidationError

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


def style_leaf_domain(path: str) -> tuple[str, str] | None:
    """Return ``(reason, requirement)`` for a font-size or stroke-width leaf."""

    for suffixes, reason, requirement in _STYLE_LEAF_DOMAINS:
        if any(part.endswith(suffixes) for part in path.split(".")):
            return reason, requirement
    return None


def validate_style_leaf(path: str, value: object, *, prefix: str) -> None:
    """Reject a numeric font size <= 0 or stroke width < 0 at ``path``."""

    domain = style_leaf_domain(path)
    if domain is None or value is None or isinstance(value, bool) or not isinstance(value, Real):
        return
    reason, requirement = domain
    number = float(value)
    if math.isfinite(number) and (number > 0 if reason == "POSITIVE" else number >= 0):
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
