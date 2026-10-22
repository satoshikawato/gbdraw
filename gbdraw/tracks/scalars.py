from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
import re
from typing import Literal

from gbdraw.core.numeric import JS_WHITESPACE, parse_ascii_finite_float

ScalarUnit = Literal["px", "factor"]


_PIXEL_TEXT = re.compile(
    r"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)(?:["
    + re.escape(JS_WHITESPACE)
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
        text = raw.strip(JS_WHITESPACE)
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


def _finite_scalar(text: str, raw: object) -> float:
    value = parse_ascii_finite_float(text)
    if value is None:
        raise ValueError(f"track scalar must be a finite number (px or % optional), not {raw!r}")
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
            return cls(value=_finite_scalar(s[:-2], raw), unit="px")
        if s.endswith("%"):
            return cls(value=_finite_scalar(s[:-1], raw) / 100.0, unit="factor")
        return cls(value=_finite_scalar(s, raw), unit="factor")

    def resolve(self, reference: float) -> float:
        return self.value if self.unit == "px" else (self.value * reference)


__all__ = ["ScalarSpec", "ScalarUnit"]
