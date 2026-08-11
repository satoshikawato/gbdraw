"""Shared legend-entry metrics used by renderers and composition metadata."""

from __future__ import annotations

import math


LEGEND_LINE_HEIGHT_RATIO = 24.0 / 14.0
LEGEND_TEXT_OFFSET_RATIO = 22.0 / 14.0


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


__all__ = ["legend_line_height", "legend_text_x_offset"]
