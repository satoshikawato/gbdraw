"""Shared legend-entry metrics used by renderers and composition metadata."""

from __future__ import annotations

import math

from gbdraw.exceptions import ValidationError


LEGEND_LINE_HEIGHT_RATIO = 24.0 / 14.0
LEGEND_TEXT_OFFSET_RATIO = 22.0 / 14.0


LEGEND_SWATCH_SIZE_PATH = "objects.legends.color_rect_size"


def legend_swatch_size(value: float, *, length_param: str | None = None) -> float:
    """The swatch size a drawn Legend measures its rows with.

    OV-1247 (D-10): the setting accepts any number; a Legend is drawn only with
    a finite size greater than zero. ``length_param`` (``short`` or ``long``)
    names the config leaf the size came from.
    """

    size = float(value)
    if not math.isfinite(size) or size <= 0.0:
        path = f"{LEGEND_SWATCH_SIZE_PATH}.{length_param}" if length_param else LEGEND_SWATCH_SIZE_PATH
        raise ValidationError(
            f"{path} must be a finite number greater than zero to draw a Legend, not {value!r}.",
            diagnostic={"code": "INPUT_INVALID", "reason": "POSITIVE", "configPath": path},
        )
    return size


def legend_line_height(color_rect_size: float) -> float:
    """Return the authoritative baseline step for legend entries."""

    return LEGEND_LINE_HEIGHT_RATIO * legend_swatch_size(color_rect_size)


def legend_text_x_offset(color_rect_size: float) -> float:
    """Return the authoritative rectangle-to-text offset for legend entries."""

    return LEGEND_TEXT_OFFSET_RATIO * legend_swatch_size(color_rect_size)


__all__ = ["LEGEND_SWATCH_SIZE_PATH", "legend_line_height", "legend_swatch_size", "legend_text_x_offset"]
