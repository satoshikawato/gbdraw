#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

from ..core.text import calculate_bbox_dimensions


LINE_MARGIN_RATIO = 24 / 14
X_MARGIN_RATIO = 22 / 14
GRADIENT_BAR_WIDTH_RATIO = 10
SINGLE_GRADIENT_TRAILING_GAP_RATIO = 0.35


@dataclass(frozen=True)
class CircularLegendEntryLayout:
    key: str
    properties: Mapping[str, object]
    text_width: float
    entry_width: float
    rect_x: float
    rect_y: float
    text_x: float
    text_y: float


@dataclass(frozen=True)
class CircularCompactGradientEntryLayout:
    key: str
    properties: Mapping[str, object]
    label_y: float
    bar_y: float


@dataclass(frozen=True)
class CircularSingleGradientEntryLayout:
    key: str
    properties: Mapping[str, object]
    title_x: float
    title_y: float
    bar_x: float
    bar_y: float
    min_label_x: float
    max_label_x: float
    scale_label_y: float


@dataclass(frozen=True)
class CircularGradientLegendLayout:
    compact: bool
    width: float
    height: float
    bar_width: float
    bar_x: float
    min_label_text: str
    scale_y: float
    compact_entries: tuple[CircularCompactGradientEntryLayout, ...]
    single_entries: tuple[CircularSingleGradientEntryLayout, ...]


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
    gradient: CircularGradientLegendLayout | None
    gradient_x: float
    gradient_y: float


@dataclass(frozen=True)
class _MeasuredLegendEntry:
    key: str
    properties: Mapping[str, object]
    text_width: float
    entry_width: float


def min_gradient_label_text(value: object) -> str:
    min_identity = float(value or 0)
    if min_identity == int(min_identity):
        return f"{int(min_identity)}%"
    return f"{min_identity}%"


def _measure_entries(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    font_family: str,
    font_size: float,
    dpi: int,
    x_margin: float,
) -> tuple[tuple[_MeasuredLegendEntry, ...], tuple[_MeasuredLegendEntry, ...]]:
    solid_entries: list[_MeasuredLegendEntry] = []
    gradient_entries: list[_MeasuredLegendEntry] = []
    for key, properties in legend_table.items():
        normalized_key = str(key)
        text_width, _ = calculate_bbox_dimensions(
            normalized_key,
            font_family,
            font_size,
            dpi,
        )
        measured = _MeasuredLegendEntry(
            key=normalized_key,
            properties=properties,
            text_width=float(text_width),
            entry_width=float(text_width) + (2 * float(x_margin)),
        )
        if properties.get("type") == "gradient":
            gradient_entries.append(measured)
        elif properties.get("type") == "solid":
            solid_entries.append(measured)
    return tuple(solid_entries), tuple(gradient_entries)


def _gradient_reserved_width(
    gradient_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    color_rect_size: float,
    x_margin: float,
) -> float:
    if not gradient_entries:
        return 0.0
    width = GRADIENT_BAR_WIDTH_RATIO * float(color_rect_size)
    if len(gradient_entries) > 1:
        width += max(entry.text_width for entry in gradient_entries) + float(x_margin)
    return width


def _build_gradient_layout(
    gradient_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    font_family: str,
    font_size: float,
    dpi: int,
    color_rect_size: float,
) -> CircularGradientLegendLayout | None:
    if not gradient_entries:
        return None

    color_rect_size = float(color_rect_size)
    bar_width = GRADIENT_BAR_WIDTH_RATIO * color_rect_size
    row_height = LINE_MARGIN_RATIO * color_rect_size
    first_min_label = min_gradient_label_text(gradient_entries[0].properties.get("min_value", 0))

    if len(gradient_entries) > 1:
        label_width = max(entry.text_width for entry in gradient_entries)
        bar_x = label_width + color_rect_size
        scale_y = color_rect_size + ((len(gradient_entries) - 1) * row_height) + 2
        _, min_label_height = calculate_bbox_dimensions(
            first_min_label,
            font_family,
            font_size,
            dpi,
        )
        _, max_label_height = calculate_bbox_dimensions(
            "100%",
            font_family,
            font_size,
            dpi,
        )
        compact_entries = tuple(
            CircularCompactGradientEntryLayout(
                key=entry.key,
                properties=entry.properties,
                label_y=(color_rect_size / 2) + (idx * row_height),
                bar_y=(color_rect_size / 2) + (idx * row_height),
            )
            for idx, entry in enumerate(gradient_entries)
        )
        return CircularGradientLegendLayout(
            compact=True,
            width=bar_x + bar_width,
            height=scale_y + max(float(min_label_height), float(max_label_height)),
            bar_width=bar_width,
            bar_x=bar_x,
            min_label_text=first_min_label,
            scale_y=scale_y,
            compact_entries=compact_entries,
            single_entries=(),
        )

    entry = gradient_entries[0]
    title_width, title_height = calculate_bbox_dimensions(
        entry.key,
        font_family,
        font_size,
        dpi,
    )
    bar_y = float(title_height) + (color_rect_size / 2)
    scale_label_y = bar_y + (color_rect_size / 2) + 2
    _, label_height = calculate_bbox_dimensions(
        "100%",
        font_family,
        font_size,
        dpi,
    )
    height = (
        scale_label_y
        + float(label_height)
        + (row_height * SINGLE_GRADIENT_TRAILING_GAP_RATIO)
    )
    single_entry = CircularSingleGradientEntryLayout(
        key=entry.key,
        properties=entry.properties,
        title_x=bar_width / 2.0,
        title_y=0.0,
        bar_x=0.0,
        bar_y=bar_y,
        min_label_x=0.0,
        max_label_x=bar_width,
        scale_label_y=scale_label_y,
    )
    return CircularGradientLegendLayout(
        compact=False,
        width=max(bar_width, float(title_width)),
        height=max(0.0, height),
        bar_width=bar_width,
        bar_x=0.0,
        min_label_text=min_gradient_label_text(entry.properties.get("min_value", 0)),
        scale_y=scale_label_y,
        compact_entries=(),
        single_entries=(single_entry,),
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
    gradient_layout: CircularGradientLegendLayout | None,
    *,
    pairwise_legend_width: float,
    canvas_width: float,
    color_rect_size: float,
    x_margin: float,
) -> float:
    solid_single_row_width = sum(entry.entry_width for entry in solid_entries)
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
    gradient_layout: CircularGradientLegendLayout | None,
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
        row_width = sum(float(entry.entry_width) for entry in row_entries)
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
    dpi: int,
    color_rect_size: float,
    legend_width: float | None = None,
    pairwise_legend_width: float | None = None,
) -> CircularLegendLayout:
    line_margin = LINE_MARGIN_RATIO * float(color_rect_size)
    x_margin = X_MARGIN_RATIO * float(color_rect_size)
    solid_entries, gradient_entries = _measure_entries(
        legend_table,
        font_family=font_family,
        font_size=font_size,
        dpi=dpi,
        x_margin=x_margin,
    )
    gradient_layout = _build_gradient_layout(
        gradient_entries,
        font_family=font_family,
        font_size=font_size,
        dpi=dpi,
        color_rect_size=float(color_rect_size),
    )
    resolved_pairwise_width = (
        float(pairwise_legend_width)
        if pairwise_legend_width is not None
        else _gradient_reserved_width(
            gradient_entries,
            color_rect_size=float(color_rect_size),
            x_margin=x_margin,
        )
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
            gradient_x = max(
                feature_width + x_margin,
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
        num_lines += len(gradient_entries) + (1 if len(gradient_entries) > 1 else 2)
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
    "CircularCompactGradientEntryLayout",
    "CircularGradientLegendLayout",
    "CircularLegendEntryLayout",
    "CircularLegendLayout",
    "CircularSingleGradientEntryLayout",
    "build_circular_legend_layout",
    "min_gradient_label_text",
]
