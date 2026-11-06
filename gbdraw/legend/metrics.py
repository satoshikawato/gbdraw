"""Shared legend-entry metrics used by renderers and composition metadata."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

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


def single_gradient_block_x(bar_width: float, title_width: float) -> tuple[float, float]:
    """The width of a one-gradient Legend block and the x of its bar.

    The title is centred over the bar and may be wider than it; the block is
    then as wide as the title and the bar is centred under it, so the whole
    block starts at x=0 (OV-1261 Circular, OV-1640 Linear).
    """

    width = max(float(bar_width), float(title_width))
    return width, (width - float(bar_width)) / 2.0


def gradient_label_text(value: object) -> str:
    """The Legend text of one end of an identity scale: ``70%``, ``72.5%``."""
    # Non-numeric values were never accepted; legend properties are numeric or numeric strings.
    identity = float(cast("float | str", value or 0))
    if identity == int(identity):
        return f"{int(identity)}%"
    return f"{identity}%"


def gradient_scale_texts(properties: Mapping[str, object]) -> tuple[str, str]:
    """The two texts of a gradient row's scale, its start and its end.

    A row that names its ends (``min_label``, ``max_label``: an annotation gradient's ``start`` and ``end``) shows
    those names; an identity row shows its ``min_value`` and ``max_value`` as percentages. Every Legend layout reads
    a gradient or ring row's scale texts here.
    """

    min_label, max_label = properties.get("min_label"), properties.get("max_label")
    return (
        min_label if isinstance(min_label, str) else gradient_label_text(properties.get("min_value", 0)),
        max_label if isinstance(max_label, str) else gradient_label_text(properties["max_value"]),
    )


# The attributes of a drawn gradient row that names its scale's ends, which the Web reads (generated into its
# ``services/drawn-attributes.generated.js``).
GRADIENT_SCALE_NAME_ATTRIBUTES = ("data-legend-min-label", "data-legend-max-label")


def write_gradient_scale_names(attributes: dict[str, object], properties: Mapping[str, object]) -> None:
    """Write a gradient row's named scale ends (``min_label``, ``max_label``) on its drawn entry; an identity row,
    whose ends are its percentages, gets none."""

    for attribute, name in zip(GRADIENT_SCALE_NAME_ATTRIBUTES, ("min_label", "max_label"), strict=True):
        value = properties.get(name)
        if isinstance(value, str):
            attributes[attribute] = value


@dataclass(frozen=True, slots=True)
class GradientScaleLayout:
    """One scale line below gradient bars: its start and end texts, under the bars' two ends."""

    min_label_text: str
    max_label_text: str
    min_label_x: float
    max_label_x: float
    scale_label_y: float


@dataclass(frozen=True, slots=True)
class CompactGradientRows:
    """Gradient rows stacked by scale: ``rows`` holds (index of the input row, centre y of its bar and caption) in
    drawn order, ``scales`` one scale line per group, and ``height`` the bottom of the last scale line."""

    rows: tuple[tuple[int, float], ...]
    scales: tuple[GradientScaleLayout, ...]
    height: float


def stack_compact_gradient_rows(
    scale_texts: Sequence[tuple[str, str]],
    *,
    color_rect_size: float,
    line_height: float,
    label_height: float,
    min_label_x: float,
    max_label_x: float,
) -> CompactGradientRows:
    """Stack the rows of a compact gradient block: the rows that share scale texts, in order of first appearance,
    then their one scale line; the next group starts half a legend line below that line."""

    groups: dict[tuple[str, str], list[int]] = {}
    for index, texts in enumerate(scale_texts):
        groups.setdefault(texts, []).append(index)
    rows: list[tuple[int, float]] = []
    scales: list[GradientScaleLayout] = []
    top = 0.0
    for (min_text, max_text), indexes in groups.items():
        if scales:
            top = scales[-1].scale_label_y + label_height + (line_height / 2.0)
        rows.extend(
            (index, top + (color_rect_size / 2.0) + (position * line_height))
            for position, index in enumerate(indexes)
        )
        scale_label_y = top + color_rect_size + ((len(indexes) - 1) * line_height) + 2.0
        scales.append(GradientScaleLayout(min_text, max_text, min_label_x, max_label_x, scale_label_y))
    return CompactGradientRows(tuple(rows), tuple(scales), scales[-1].scale_label_y + label_height)


__all__ = [
    "CompactGradientRows",
    "GRADIENT_SCALE_NAME_ATTRIBUTES",
    "GradientScaleLayout",
    "LEGEND_SWATCH_SIZE_PATH",
    "gradient_label_text",
    "gradient_scale_texts",
    "legend_line_height",
    "legend_swatch_size",
    "legend_text_x_offset",
    "single_gradient_block_x",
    "stack_compact_gradient_rows",
    "write_gradient_scale_names",
]
