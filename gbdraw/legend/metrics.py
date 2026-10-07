"""Shared legend-entry metrics used by renderers and composition metadata."""

from __future__ import annotations

import math
from collections.abc import Iterable


LEGEND_LINE_HEIGHT_RATIO = 24.0 / 14.0
LEGEND_TEXT_OFFSET_RATIO = 22.0 / 14.0


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


def _positive_size(value: float) -> float:
    size = float(value)
    if not math.isfinite(size) or size <= 0.0:
        raise ValueError("legend color rectangle size must be finite and positive")
    return size


def legend_line_height(color_rect_size: float) -> float:
    """Return the authoritative baseline step for legend entries."""

    return LEGEND_LINE_HEIGHT_RATIO * _positive_size(color_rect_size)


def legend_text_x_offset(color_rect_size: float) -> float:
    """Return the authoritative rectangle-to-text offset for legend entries."""

    return LEGEND_TEXT_OFFSET_RATIO * _positive_size(color_rect_size)


__all__ = ["compensated_sum", "legend_line_height", "legend_text_x_offset"]
