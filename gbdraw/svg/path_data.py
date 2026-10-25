"""The one writer of the numbers in SVG path data.

Every ``<path d=...>`` that gbdraw draws writes its coordinates, radii and
control points through :func:`path_number`: rounded to ``PATH_DECIMALS``
decimals in the shortest form (no trailing zeros, no ``-0``, no exponent).
At 3 decimals a point moves by at most 0.0005 SVG px.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike

PATH_DECIMALS = 3
_FORMAT = f".{PATH_DECIMALS}f"


def path_number(value: float) -> str:
    """``value`` as path data writes it: ``PATH_DECIMALS`` decimals, shortest form."""

    text = format(value, _FORMAT).rstrip("0")
    if text[-1] == ".":
        text = text[:-1]
    return "0" if text == "-0" else text


def path_numbers(values: ArrayLike) -> list[str]:
    """:func:`path_number` of each value of a numeric array."""

    return [path_number(value) for value in np.asarray(values, dtype=float).tolist()]


def path_point(x: float, y: float) -> str:
    """One ``x,y`` pair of path data."""

    return f"{path_number(x)},{path_number(y)}"


def path_line_to(xs: ArrayLike, ys: ArrayLike) -> str:
    """``L x y`` for each point, joined (the polyline of a track path)."""

    return "".join(f"L{x} {y}" for x, y in zip(path_numbers(xs), path_numbers(ys)))


def join_path_numbers(values: Iterable[float], separator: str = " ") -> str:
    """Several path numbers joined by ``separator``."""

    return separator.join(path_number(value) for value in values)


__all__ = [
    "PATH_DECIMALS",
    "join_path_numbers",
    "path_line_to",
    "path_number",
    "path_numbers",
    "path_point",
]
