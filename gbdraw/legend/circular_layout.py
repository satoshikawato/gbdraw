#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, cast

from ..core.text import drawn_text_box
from ..core.numeric import compensated_sum
from .metrics import legend_line_height, legend_text_x_offset
from .table import DrawnLegendName, LegendRowKey, as_legend_table


GRADIENT_BAR_WIDTH_RATIO = 10
# The similarity-ring block (PD-OI-110): its captions, and the label of its identity band.
RING_OUTER_CAPTION = "outer"
RING_INNER_CAPTION = "inner"
RING_BAND_LABEL = "identity"
# The band sits this many lines below the ``inner`` caption.
RING_BAND_LINE_GAP = 1.5


@dataclass(frozen=True)
class CircularLegendEntryLayout:
    key: LegendRowKey
    properties: Mapping[str, object]
    text_width: float
    entry_width: float
    rect_x: float
    rect_y: float
    text_x: float
    text_y: float


@dataclass(frozen=True)
class CircularRingRowLayout:
    """One drawn ring's row: its swatch (when the block shows swatches) and its name, centred on ``y``."""

    key: LegendRowKey
    properties: Mapping[str, object]
    swatch_x: float
    text_x: float
    y: float


@dataclass(frozen=True)
class CircularRingLegendLayout:
    """The similarity-ring block (PD-OI-110): caption ``outer``, the drawn rings outermost first, caption
    ``inner``, then the identity band labelled ``identity`` with its minimum and 100% labels below it."""

    width: float
    height: float
    show_swatches: bool
    outer_caption_y: float
    inner_caption_y: float
    rows: tuple[CircularRingRowLayout, ...]
    bar_width: float
    band_y: float
    band_label_x: float
    min_label_text: str
    max_label_text: str
    min_label_x: float
    max_label_x: float
    scale_label_y: float


@dataclass(frozen=True)
class CircularLegendLayout:
    horizontal: bool
    width: float
    height: float
    feature_width: float
    feature_height: float
    pairwise_legend_width: float
    line_margin: float
    x_margin: float
    has_gradient: bool
    num_lines: int
    num_columns: int
    num_items_per_line: int
    solid_entries: tuple[CircularLegendEntryLayout, ...]
    # The similarity-ring block, beside or below the solid rows.
    gradient: CircularRingLegendLayout | None
    gradient_x: float
    gradient_y: float


@dataclass(frozen=True)
class _MeasuredLegendEntry:
    key: LegendRowKey
    text: DrawnLegendName
    properties: Mapping[str, object]
    text_width: float
    entry_width: float


def gradient_label_text(value: object) -> str:
    """The Legend text of one end of an identity scale: ``70%``, ``72.5%``."""
    # Non-numeric values were never accepted; legend properties are numeric or numeric strings.
    identity = float(cast("float | str", value or 0))
    if identity == int(identity):
        return f"{int(identity)}%"
    return f"{identity}%"


def _measure_entries(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    font_family: str,
    font_size: float,
    x_margin: float,
) -> tuple[tuple[_MeasuredLegendEntry, ...], tuple[_MeasuredLegendEntry, ...]]:
    solid_entries: list[_MeasuredLegendEntry] = []
    ring_entries: list[_MeasuredLegendEntry] = []
    table = as_legend_table(legend_table)
    drawn_names = table.names()
    for key, properties in table.items():
        text = drawn_names[key].drawn
        text_width, _ = drawn_text_box(
            text,
            font_family,
            font_size,
        )
        measured = _MeasuredLegendEntry(
            key=key,
            text=text,
            properties=properties,
            text_width=float(text_width),
            entry_width=float(text_width) + (2 * float(x_margin)),
        )
        if properties.get("type") == "ring":
            ring_entries.append(measured)
        elif properties.get("type") == "solid":
            solid_entries.append(measured)
    return tuple(solid_entries), tuple(ring_entries)


def _build_ring_layout(
    ring_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    font_family: str,
    font_size: float,
    color_rect_size: float,
) -> CircularRingLegendLayout | None:
    if not ring_entries:
        return None
    rect = float(color_rect_size)
    line = legend_line_height(rect)
    show_swatches = any(entry.properties.get("fill") for entry in ring_entries)
    text_x = legend_text_x_offset(rect) if show_swatches else 0.0
    bar_width = GRADIENT_BAR_WIDTH_RATIO * rect
    band_label_x = bar_width + (legend_text_x_offset(rect) - rect)
    outer_caption_y = rect / 2
    rows = tuple(
        CircularRingRowLayout(
            key=entry.key,
            properties=entry.properties,
            swatch_x=0.0,
            text_x=text_x,
            y=outer_caption_y + ((index + 1) * line),
        )
        for index, entry in enumerate(ring_entries)
    )
    inner_caption_y = outer_caption_y + ((len(ring_entries) + 1) * line)
    band_y = inner_caption_y + (RING_BAND_LINE_GAP * line)
    scale_label_y = band_y + (rect / 2) + 2
    _, label_height = drawn_text_box("100%", font_family, font_size)
    caption_width = max(
        float(drawn_text_box(text, font_family, font_size)[0]) for text in (RING_OUTER_CAPTION, RING_INNER_CAPTION)
    )
    band_width = band_label_x + float(drawn_text_box(RING_BAND_LABEL, font_family, font_size)[0])
    width = max(caption_width, band_width, *(text_x + entry.text_width for entry in ring_entries))
    height = scale_label_y + float(label_height)
    return CircularRingLegendLayout(
        width=width,
        height=height,
        show_swatches=show_swatches,
        outer_caption_y=outer_caption_y,
        inner_caption_y=inner_caption_y,
        rows=rows,
        bar_width=bar_width,
        band_y=band_y,
        band_label_x=band_label_x,
        min_label_text=gradient_label_text(ring_entries[0].properties.get("min_value", 0)),
        max_label_text=gradient_label_text(ring_entries[0].properties["max_value"]),
        min_label_x=0.0,
        max_label_x=bar_width,
        scale_label_y=scale_label_y,
    )


def _wrap_horizontal_rows(
    entries: tuple[_MeasuredLegendEntry, ...],
    wrap_width: float,
    *,
    x_margin: float,
) -> tuple[tuple[_MeasuredLegendEntry, ...], ...]:
    if not entries:
        return ()
    if wrap_width <= 0:
        wrap_width = float("inf")

    rows: list[list[_MeasuredLegendEntry]] = []
    current_row: list[_MeasuredLegendEntry] = []
    current_row_width = 0.0
    for entry in entries:
        exceeds_wrap = (
            math.isfinite(wrap_width)
            and current_row
            and ((current_row_width + float(x_margin) + entry.entry_width) > wrap_width)
        )
        if exceeds_wrap:
            rows.append(current_row)
            current_row = []
            current_row_width = 0.0
        current_row.append(entry)
        current_row_width += entry.entry_width
    if current_row:
        rows.append(current_row)
    return tuple(tuple(row) for row in rows)


def _horizontal_legend_width(
    solid_entries: tuple[_MeasuredLegendEntry, ...],
    gradient_layout: CircularRingLegendLayout | None,
    *,
    pairwise_legend_width: float,
    canvas_width: float,
    color_rect_size: float,
    x_margin: float,
) -> float:
    solid_single_row_width = compensated_sum(entry.entry_width for entry in solid_entries)
    solid_desired_width = solid_single_row_width + (2 * float(x_margin) if solid_entries else 0.0)
    desired_width = solid_desired_width
    if gradient_layout is not None:
        desired_width += max(float(pairwise_legend_width), float(gradient_layout.width))
    if desired_width <= 0:
        return 0.0

    min_solid_width = (
        max((entry.entry_width + (2 * float(x_margin)) for entry in solid_entries), default=0.0)
    )
    min_width = min_solid_width
    if gradient_layout is not None:
        min_width += max(float(pairwise_legend_width), float(gradient_layout.width))
        if solid_entries:
            min_width += float(x_margin)
    if min_width <= 0:
        min_width = float(color_rect_size)
    return max(min_width, min(desired_width, float(canvas_width)))


def _build_horizontal_layout(
    solid_entries: tuple[_MeasuredLegendEntry, ...],
    gradient_layout: CircularRingLegendLayout | None,
    *,
    width: float,
    pairwise_legend_width: float,
    color_rect_size: float,
    line_margin: float,
    x_margin: float,
) -> tuple[tuple[CircularLegendEntryLayout, ...], float, float, float, int, int, int]:
    wrap_width = float(width)
    if gradient_layout is not None:
        wrap_width = max(float(color_rect_size), wrap_width - float(pairwise_legend_width) - float(x_margin))
    if wrap_width <= 0:
        wrap_width = float("inf")

    rows = _wrap_horizontal_rows(solid_entries, wrap_width, x_margin=x_margin)
    feature_height = (
        float(color_rect_size) + ((len(rows) - 1) * float(line_margin))
        if rows
        else 0.0
    )
    gradient_height = float(gradient_layout.height) if gradient_layout is not None else 0.0
    feature_y_offset = max(0.0, (gradient_height - feature_height) / 2.0)

    laid_out_entries: list[CircularLegendEntryLayout] = []
    for row_index, row_entries in enumerate(rows):
        row_y = (
            feature_y_offset + (float(color_rect_size) / 2.0) + (row_index * float(line_margin))
            if gradient_layout is not None
            else row_index * float(line_margin)
        )
        row_width = compensated_sum(float(entry.entry_width) for entry in row_entries)
        if math.isfinite(wrap_width):
            row_start_x = float(x_margin) + max(0.0, (wrap_width - row_width) * 0.5)
        else:
            row_start_x = float(x_margin)
        current_x = row_start_x
        for entry in row_entries:
            laid_out_entries.append(
                CircularLegendEntryLayout(
                    key=entry.key,
                    properties=entry.properties,
                    text_width=entry.text_width,
                    entry_width=entry.entry_width,
                    rect_x=current_x - float(x_margin),
                    rect_y=row_y,
                    text_x=current_x,
                    text_y=row_y,
                )
            )
            current_x += float(entry.entry_width)

    num_lines = max(1, len(rows)) if solid_entries or gradient_layout is not None else 0
    num_items_per_line = max((len(row) for row in rows), default=0)
    num_columns = num_items_per_line + (1 if gradient_layout is not None else 0)
    return (
        tuple(laid_out_entries),
        wrap_width if math.isfinite(wrap_width) else width,
        feature_height,
        feature_y_offset,
        num_lines,
        num_columns,
        num_items_per_line,
    )


def build_circular_legend_layout(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    legend_position: str,
    canvas_width: float,
    font_family: str,
    font_size: float,
    color_rect_size: float,
    legend_width: float | None = None,
    pairwise_legend_width: float | None = None,
) -> CircularLegendLayout:
    line_margin = legend_line_height(color_rect_size)
    x_margin = legend_text_x_offset(color_rect_size)
    solid_entries, ring_entries = _measure_entries(
        legend_table,
        font_family=font_family,
        font_size=font_size,
        x_margin=x_margin,
    )
    gradient_layout = _build_ring_layout(
        ring_entries,
        font_family=font_family,
        font_size=font_size,
        color_rect_size=float(color_rect_size),
    )
    resolved_pairwise_width = (
        float(pairwise_legend_width)
        if pairwise_legend_width is not None
        else gradient_layout.width
        if gradient_layout is not None
        else 0.0
    )
    horizontal = str(legend_position) in {"top", "bottom"}

    if horizontal:
        width = (
            float(legend_width)
            if legend_width is not None
            else _horizontal_legend_width(
                solid_entries,
                gradient_layout,
                pairwise_legend_width=resolved_pairwise_width,
                canvas_width=float(canvas_width),
                color_rect_size=float(color_rect_size),
                x_margin=x_margin,
            )
        )
        (
            laid_out_entries,
            feature_width,
            feature_height,
            _feature_y_offset,
            num_lines,
            num_columns,
            num_items_per_line,
        ) = _build_horizontal_layout(
            solid_entries,
            gradient_layout,
            width=width,
            pairwise_legend_width=resolved_pairwise_width,
            color_rect_size=float(color_rect_size),
            line_margin=line_margin,
            x_margin=x_margin,
        )
        gradient_x = 0.0
        gradient_y = 0.0
        content_bottom = feature_height
        if gradient_layout is not None:
            gradient_y = max(0.0, (feature_height - float(gradient_layout.height)) / 2.0)
            # Without feature rows the gradient block is the whole Legend (OV-1612).
            gradient_x = max(
                feature_width + x_margin if solid_entries else 0.0,
                width - float(gradient_layout.width),
            )
            content_bottom = max(
                (max((entry.rect_y for entry in laid_out_entries), default=0.0) + float(color_rect_size) / 2.0)
                if laid_out_entries
                else 0.0,
                gradient_y + float(gradient_layout.height),
            )
            height = content_bottom + (float(color_rect_size) / 2.0)
        else:
            height = feature_height
        return CircularLegendLayout(
            horizontal=True,
            width=width,
            height=height,
            feature_width=feature_width,
            feature_height=feature_height,
            pairwise_legend_width=resolved_pairwise_width,
            line_margin=line_margin,
            x_margin=x_margin,
            has_gradient=gradient_layout is not None,
            num_lines=num_lines,
            num_columns=num_columns,
            num_items_per_line=num_items_per_line,
            solid_entries=laid_out_entries,
            gradient=gradient_layout,
            gradient_x=gradient_x,
            gradient_y=gradient_y,
        )

    feature_block_width = max(
        (x_margin + entry.text_width for entry in solid_entries),
        default=0.0,
    )
    width = (
        float(legend_width)
        if legend_width is not None
        else max(
            feature_block_width,
            resolved_pairwise_width,
            float(gradient_layout.width) if gradient_layout is not None else 0.0,
        )
    )
    alignment_width = max(
        width,
        feature_block_width,
        float(gradient_layout.width) if gradient_layout is not None else 0.0,
    )
    feature_x_offset = (
        max(0.0, (alignment_width - feature_block_width) / 2.0)
        if gradient_layout is not None and feature_block_width > 0
        else 0.0
    )
    laid_out_entries = tuple(
        CircularLegendEntryLayout(
            key=entry.key,
            properties=entry.properties,
            text_width=entry.text_width,
            entry_width=entry.entry_width,
            rect_x=feature_x_offset,
            rect_y=idx * line_margin,
            text_x=feature_x_offset + x_margin,
            text_y=idx * line_margin,
        )
        for idx, entry in enumerate(solid_entries)
    )
    feature_height = (
        float(color_rect_size) + ((len(solid_entries) - 1) * line_margin)
        if solid_entries
        else 0.0
    )
    gradient_x = 0.0
    gradient_y = 0.0
    if gradient_layout is not None:
        gradient_x = max(0.0, (alignment_width - float(gradient_layout.width)) / 2.0)
        gradient_y = len(solid_entries) * line_margin + (
            line_margin * 0.5 if solid_entries else 0.0
        )
    content_bottom = feature_height
    if gradient_layout is not None:
        content_bottom = max(content_bottom, gradient_y + float(gradient_layout.height))
        height = content_bottom + (float(color_rect_size) / 2.0)
    else:
        height = feature_height
    num_lines = len(solid_entries)
    if gradient_layout is not None:
        # outer, the rings, inner, the band and its scale labels
        num_lines += len(ring_entries) + 4
    return CircularLegendLayout(
        horizontal=False,
        width=width,
        height=height,
        feature_width=feature_block_width,
        feature_height=feature_height,
        pairwise_legend_width=resolved_pairwise_width,
        line_margin=line_margin,
        x_margin=x_margin,
        has_gradient=gradient_layout is not None,
        num_lines=num_lines,
        num_columns=1,
        num_items_per_line=1 if solid_entries else 0,
        solid_entries=laid_out_entries,
        gradient=gradient_layout,
        gradient_x=gradient_x,
        gradient_y=gradient_y,
    )


__all__ = [
    "CircularLegendEntryLayout",
    "CircularLegendLayout",
    "CircularRingLegendLayout",
    "CircularRingRowLayout",
    "RING_BAND_LABEL",
    "RING_INNER_CAPTION",
    "RING_OUTER_CAPTION",
    "build_circular_legend_layout",
    "gradient_label_text",
]
