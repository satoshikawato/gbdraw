"""Numeric helpers whose results must not depend on the Python version."""

from __future__ import annotations

import math
from collections.abc import Iterable


def compensated_sum(values: Iterable[float]) -> float:
    """Sum floats as CPython 3.12+ ``sum()`` does, on every supported Python.

    Python 3.12 made the built-in ``sum()`` of floats compensated (Neumaier's
    variant of Kahan summation, ``Python/bltinmodule.c``); 3.10 and 3.11 add
    naively, so the same layout came out one ulp apart. The Web port of the
    Legend layout (``pythonFloatSum`` in ``gbdraw/web/js/services/legend-layout.js``)
    computes exactly this: the first value added to ``0``, each further value
    compensated, and the compensation added at the end when it is nonzero and
    finite. Every value is taken as a float, as the port does.
    """
    iterator = iter(values)
    try:
        first = next(iterator)
    except StopIteration:
        return 0.0
    total = 0 + float(first)
    compensation = 0.0
    for value in iterator:
        item = float(value)
        partial = total + item
        if abs(total) >= abs(item):
            compensation += (total - partial) + item
        else:
            compensation += (item - partial) + total
        total = partial
    if compensation and math.isfinite(compensation):
        total += compensation
    return total


def scaled_tick_text(position: int, divisor: int, interval: int | None) -> str:
    """Return ``position / divisor`` with the decimals that ``interval`` needs.

    The decimals are the fewest (at most 6) that write ``interval / divisor``
    exactly, so a 250 bp interval in kbp gives ``0.25``, ``0.5``, ``0.75``;
    trailing zeros are dropped. Without an interval the value is rounded.
    """
    decimals = 0
    step = abs(int(interval or 0))
    while step and decimals < 6 and (step * 10**decimals) % divisor:
        decimals += 1
    value = float(position) / float(divisor)
    if decimals == 0:
        return f"{value:.0f}"
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")


__all__ = ["compensated_sum", "scaled_tick_text"]
