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


__all__ = ["compensated_sum"]
