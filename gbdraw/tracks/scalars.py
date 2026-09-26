from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
import re
from typing import Literal

ScalarUnit = Literal["px", "factor"]


# Match JavaScript whitespace, including BOM, at both text boundaries.
_PIXEL_WHITESPACE = (
    " \t\n\r\f\v\u00a0\u1680"
    + "".join(map(chr, range(0x2000, 0x200B)))
    + "\u2028\u2029\u202f\u205f\u3000\ufeff"
)
_PIXEL_TEXT = re.compile(
    r"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)(?:["
    + re.escape(_PIXEL_WHITESPACE)
    + r"]*[pP][xX])?"
)


def parse_optional_pixel(
    raw: object, *, field_name: str, allow_zero: bool
) -> float | None:
    """Parse pure pixel text before projecting it to a typed slot value."""
    if raw is None:
        return None
    relation = "nonnegative" if allow_zero else "positive"
    message = f"{field_name} must be {relation} finite number of pixels (px optional)"
    if isinstance(raw, str):
        text = raw.strip(_PIXEL_WHITESPACE)
        if not text:
            return None
        match = _PIXEL_TEXT.fullmatch(text)
        if match is None:
            raise ValueError(message)
        value = float(match[1])
    elif isinstance(raw, (int, float)) and not isinstance(raw, bool):
        try:
            value = float(raw)
        except OverflowError as exc:
            raise ValueError(message) from exc
    else:
        raise ValueError(message)
    if not isfinite(value) or value < 0 or (not allow_zero and value == 0):
        raise ValueError(message)
    return value


@dataclass(frozen=True)
class ScalarSpec:
    """A scalar value with an explicit unit."""

    value: float
    unit: ScalarUnit = "factor"

    @classmethod
    def parse(cls, raw: str) -> "ScalarSpec":
        s = str(raw).strip()
        if not s:
            raise ValueError("empty scalar")
        if s.endswith("px"):
            return cls(value=float(s[:-2]), unit="px")
        if s.endswith("%"):
            return cls(value=float(s[:-1]) / 100.0, unit="factor")
        return cls(value=float(s), unit="factor")

    def resolve(self, reference: float) -> float:
        return self.value if self.unit == "px" else (self.value * reference)


__all__ = ["ScalarSpec", "ScalarUnit"]
