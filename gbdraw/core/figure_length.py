"""A requested figure length with its unit (``180mm``, ``7in``): the one reader of ``canvas.figure_width`` and
``canvas.figure_height``.

A length is a positive decimal number and a unit, ``mm``, ``cm``, ``in`` or ``px`` (any case, an optional space).
A bare number is refused: every other gbdraw size is in px, while journals give column widths in mm, so a bare
``180`` would be misread either way. A px is a CSS px, 96 per inch.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from decimal import Decimal
from typing import Literal, TypeAlias, cast

from gbdraw.exceptions import ValidationError

FigureLengthUnit: TypeAlias = Literal["mm", "cm", "in", "px"]
CSS_PX_PER_INCH = 96.0
_PX_PER_UNIT: dict[FigureLengthUnit, float] = {
    "mm": CSS_PX_PER_INCH / 25.4,
    "cm": CSS_PX_PER_INCH / 2.54,
    "in": CSS_PX_PER_INCH,
    "px": 1.0,
}
# The decimals a size in this unit is written with in a message: a tenth of a mm or cm, a hundredth of an inch.
_MESSAGE_DECIMALS: dict[FigureLengthUnit, int] = {"mm": 1, "cm": 1, "in": 2, "px": 0}
FIGURE_LENGTH_AUTO = "auto"
_FIGURE_LENGTH_RE = re.compile(r"\s*([0-9]+(?:\.[0-9]*)?|\.[0-9]+)\s*([A-Za-z]+)\s*", re.ASCII)
_REFUSAL = "must be a length with a unit (mm, cm, in, or px), for example 180mm or 7in, or auto"
# A caller that knows the config path raises again with it (``figure_length_from_config``).
_INVALID = {"code": "INPUT_INVALID"}


def _plain_number(value: float) -> str:
    """``value`` as a plain decimal with no exponent and no trailing zeros: ``180``, ``7.5``."""

    return format(Decimal(repr(value)).normalize(), "f")


@dataclass(frozen=True)
class FigureLength:
    """A positive length in one unit; ``text`` is its canonical spelling (``180mm``)."""

    value: float
    unit: FigureLengthUnit

    def __post_init__(self) -> None:
        if self.unit not in _PX_PER_UNIT:
            raise ValidationError(f"figure length unit {self.unit!r} {_REFUSAL}.", diagnostic=_INVALID)
        if isinstance(self.value, bool) or not math.isfinite(self.value) or self.value <= 0:
            raise ValidationError(f"figure length {self.value!r} {_REFUSAL}.", diagnostic=_INVALID)

    @property
    def px(self) -> float:
        return self.value * _PX_PER_UNIT[self.unit]

    @property
    def text(self) -> str:
        return f"{_plain_number(self.value)}{self.unit}"

    def of_px(self, px: float) -> FigureLength:
        """The length ``px`` measures in this length's unit."""

        return FigureLength(px / _PX_PER_UNIT[self.unit], self.unit)

    def message_text(self, px: float, *, round_up: bool = False) -> str:
        """``px`` in this length's unit as a message writes it: ``53.9 mm``, ``180 mm``, ``2.13 in``, ``204 px``.

        With ``round_up`` the value is rounded up, as a minimum is.
        """

        decimals = _MESSAGE_DECIMALS[self.unit]
        scaled = px / _PX_PER_UNIT[self.unit] * 10**decimals
        # A float a hair over a whole step (53.900000001) is still that step.
        rounded = math.ceil(round(scaled, 6)) if round_up else round(scaled)
        return f"{_plain_number(rounded / 10**decimals)} {self.unit}"


def parse_figure_length(value: str | None) -> FigureLength | None:
    """The length ``value`` spells, or ``None`` for ``"auto"`` (any case) and ``None``.

    ``"180mm"``, ``"180 MM"``, ``"7.5in"``, ``"18cm"``, ``"680px"``. Raises ``ValidationError`` for anything else:
    a bare number, zero, a negative or non-finite number, an exponent, another unit, ``"none"`` and empty text.
    """

    if value is None:
        return None
    if not isinstance(value, str):
        raise ValidationError(f"{_REFUSAL}; got {value!r}.", diagnostic=_INVALID)
    if value.strip().lower() == FIGURE_LENGTH_AUTO:
        return None
    match = _FIGURE_LENGTH_RE.fullmatch(value)
    unit = match[2].lower() if match is not None else ""
    number = float(match[1]) if match is not None and len(match[1]) <= 32 else 0.0
    if unit not in _PX_PER_UNIT or not math.isfinite(number) or number <= 0:
        raise ValidationError(f"{_REFUSAL}; got {value!r}.", diagnostic=_INVALID)
    return FigureLength(number, cast(FigureLengthUnit, unit))


def argparse_figure_length(text: str) -> str:
    """argparse ``type=`` of ``--figure_width`` and ``--figure_height``: the canonical text, or ``"auto"``."""

    try:
        length = parse_figure_length(text)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(f"invalid figure length: {text!r} ({_REFUSAL})") from exc
    return FIGURE_LENGTH_AUTO if length is None else length.text


@dataclass(frozen=True)
class RequestedFigureSize:
    """The figure size a request asks for; a side that is ``None`` is auto."""

    width: FigureLength | None = None
    height: FigureLength | None = None

    @property
    def requested(self) -> bool:
        return self.width is not None or self.height is not None

    @property
    def width_px(self) -> float | None:
        return None if self.width is None else self.width.px

    @property
    def height_px(self) -> float | None:
        return None if self.height is None else self.height.px


__all__ = [
    "CSS_PX_PER_INCH",
    "FIGURE_LENGTH_AUTO",
    "FigureLength",
    "FigureLengthUnit",
    "RequestedFigureSize",
    "argparse_figure_length",
    "parse_figure_length",
]
