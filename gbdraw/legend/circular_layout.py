#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Mapping

from ..core.text import drawn_text_box
from ..core.numeric import compensated_sum
from .metrics import (
    GRADIENT_TEXT_GAP,
    GradientScaleLayout,
    gradient_scale_texts,
    legend_caption_box,
    legend_line_height,
    legend_text_x_offset,
    single_gradient_block_x,
    single_gradient_row_ys,
    stack_compact_gradient_rows,
)
from .table import DrawnLegendName, LegendRowKey, as_legend_table


GRADIENT_BAR_WIDTH_RATIO = 10
# The similarity-ring block (PD-OI-110): its captions, and the label of its identity scale.
RING_OUTER_CAPTION = "outer"
RING_INNER_CAPTION = "inner"
RING_BAND_LABEL = "identity"
# One colour: the band sits this many lines below the ``inner`` caption.
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
class CircularCompactGradientEntryLayout:
    key: LegendRowKey
    properties: Mapping[str, object]
    label_x: float
    label_y: float
    bar_y: float


@dataclass(frozen=True)
class CircularSingleGradientEntryLayout:
    key: LegendRowKey
    properties: Mapping[str, object]
    title_x: float
    title_y: float
    bar_x: float
    bar_y: float


@dataclass(frozen=True)
class CircularGradientLegendLayout:
    """The gradient rows (an annotation gradient's bar): two or more rows with their captions left of the bars
    (``compact``), or one row with its title centred over its bar.

    The rows that share scale texts share one scale line below them. The ``min_label_*``, ``max_label_*`` and
    ``scale_label_y`` fields are the first group's scale line; ``more_scales`` holds the later groups' lines, each
    below its own rows. Every y is in the side block (``CircularLegendLayout.gradient_x``, ``gradient_y``): below
    a similarity-ring block it includes the ring block's height and the gap; ``width`` and ``height`` are the rows'
    own."""

    compact: bool
    width: float
    height: float
    bar_width: float
    bar_x: float
    min_label_text: str
    max_label_text: str
    min_label_x: float
    max_label_x: float
    scale_label_y: float
    compact_entries: tuple[CircularCompactGradientEntryLayout, ...]
    single_entries: tuple[CircularSingleGradientEntryLayout, ...]
    more_scales: tuple[GradientScaleLayout, ...] = ()

    @property
    def scales(self) -> tuple[GradientScaleLayout, ...]:
        """Every scale line of the rows, the first group's first."""

        first = GradientScaleLayout(
            self.min_label_text, self.max_label_text, self.min_label_x, self.max_label_x, self.scale_label_y
        )
        return (first, *self.more_scales)

    def shifted(self, dx: float) -> CircularGradientLegendLayout:
        """The same rows ``dx`` to the right in the side block."""

        if dx == 0:
            return self
        return replace(
            self,
            bar_x=self.bar_x + dx,
            min_label_x=self.min_label_x + dx,
            max_label_x=self.max_label_x + dx,
            more_scales=tuple(
                replace(scale, min_label_x=scale.min_label_x + dx, max_label_x=scale.max_label_x + dx)
                for scale in self.more_scales
            ),
            compact_entries=tuple(replace(entry, label_x=entry.label_x + dx) for entry in self.compact_entries),
            single_entries=tuple(
                replace(entry, title_x=entry.title_x + dx, bar_x=entry.bar_x + dx) for entry in self.single_entries
            ),
        )


@dataclass(frozen=True)
class CircularRingRowLayout:
    """One drawn ring's row: its name at the block's ``names_x`` and, when the rings differ in colour, its own
    identity ramp at the block's ``bar_x``, both centred on ``y``."""

    key: LegendRowKey
    properties: Mapping[str, object]
    y: float


@dataclass(frozen=True)
class CircularRingLegendLayout:
    """The similarity-ring block (PD-OI-110, scenario revision 2), in one of two forms.

    - One colour: caption ``outer``, the drawn rings' names outermost first, caption
      ``inner``, then one identity band at ``band_y``. The names block (captions and names, left-aligned at
      ``names_x``) and the band (at ``bar_x``) are two blocks; a vertical Legend centres each on the block's width,
      a horizontal one starts both at the block's left edge.
    - Colours that differ (``per_ring_ramps``): each row is the ring's name at ``names_x`` and its own identity
      ramp at ``bar_x``, as one block, and no band (``band_y`` is None). The captions stand above and below the
      names only when two or more drawn rings share a colour (``show_captions``).

    Either form has one scale line under its bars: the minimum at the bars' left end, ``identity`` centred on
    them, 100% at their right end."""

    width: float
    height: float
    show_captions: bool
    outer_caption_y: float
    inner_caption_y: float
    rows: tuple[CircularRingRowLayout, ...]
    names_x: float
    bar_x: float
    bar_width: float
    band_y: float | None
    min_label_text: str
    max_label_text: str
    scale_label_y: float

    @property
    def min_label_x(self) -> float:
        return self.bar_x

    @property
    def max_label_x(self) -> float:
        return self.bar_x + self.bar_width

    @property
    def identity_label_x(self) -> float:
        """``identity`` stands centred under the bars, on the scale line."""

        return self.bar_x + (self.bar_width / 2)

    @property
    def identity_label_y(self) -> float:
        return self.scale_label_y

    def shifted(self, dx: float) -> CircularRingLegendLayout:
        """The same block ``dx`` to the right in the side block."""

        return self if dx == 0 else replace(self, names_x=self.names_x + dx, bar_x=self.bar_x + dx)

    @property
    def per_ring_ramps(self) -> bool:
        """Whether each ring draws its own ramp (the rings differ in colour) instead of the one band."""

        return self.band_y is None

    @property
    def bar_ys(self) -> tuple[float, ...]:
        """The centre of each bar the block draws at ``bar_x``: every ring's own ramp, or the one band."""

        return (self.band_y,) if self.band_y is not None else tuple(row.y for row in self.rows)

    @property
    def lines(self) -> int:
        """The block's lines in a vertical Legend: the captions, the rings, and the band and its scale labels (one
        colour) or the one scale line (colours that differ)."""

        return (2 if self.show_captions else 0) + len(self.rows) + (1 if self.per_ring_ramps else 2)


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
    # The similarity-ring block, beside or below the solid rows, at the side block's origin.
    gradient: CircularRingLegendLayout | None
    # Where the side block sits: the ring block, then the gradient rows below it.
    gradient_x: float
    gradient_y: float
    # The gradient rows (annotation gradients), below the ring block in the side block.
    gradient_rows: CircularGradientLegendLayout | None = None


@dataclass(frozen=True)
class _SideBlock:
    """The block beside or below the solid rows: the ring block, then the gradient rows half a legend line below
    it; its size, the width it reserves when no pairwise Legend width is given, and its lines in a vertical Legend."""

    ring: CircularRingLegendLayout | None
    rows: CircularGradientLegendLayout | None
    width: float
    height: float
    reserved_width: float
    lines: int


@dataclass(frozen=True)
class _MeasuredLegendEntry:
    key: LegendRowKey
    text: DrawnLegendName
    properties: Mapping[str, object]
    text_width: float
    entry_width: float


def _measure_entries(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    font_family: str,
    font_size: float,
    x_margin: float,
) -> tuple[tuple[_MeasuredLegendEntry, ...], tuple[_MeasuredLegendEntry, ...], tuple[_MeasuredLegendEntry, ...]]:
    solid_entries: list[_MeasuredLegendEntry] = []
    gradient_entries: list[_MeasuredLegendEntry] = []
    ring_entries: list[_MeasuredLegendEntry] = []
    table = as_legend_table(legend_table)
    drawn_names = table.names()
    for key, properties in table.items():
        text = drawn_names[key].drawn
        text_width, _ = legend_caption_box(text, font_family, font_size)
        measured = _MeasuredLegendEntry(
            key=key,
            text=text,
            properties=properties,
            text_width=float(text_width),
            entry_width=float(text_width) + (2 * float(x_margin)),
        )
        if properties.get("type") == "gradient":
            gradient_entries.append(measured)
        elif properties.get("type") == "ring":
            ring_entries.append(measured)
        elif properties.get("type") == "solid":
            solid_entries.append(measured)
    return tuple(solid_entries), tuple(gradient_entries), tuple(ring_entries)


def _build_ring_layout(
    ring_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    font_family: str,
    font_size: float,
    color_rect_size: float,
    centred: bool,
) -> CircularRingLegendLayout | None:
    """The similarity-ring block. A row's ``fill`` is its ring's colour when the drawn rings differ in colour
    (``_add_conservation_ring_rows``), else None: rows with a colour draw their own ramps, and rows that share
    a colour (one colour included) keep the captions that give the order outer to inner. ``centred`` (a vertical
    Legend) centres the one-colour names block and band on the block's width."""

    if not ring_entries:
        return None
    rect = float(color_rect_size)
    line = legend_line_height(rect)
    fills = [entry.properties.get("fill") for entry in ring_entries]
    per_ring_ramps = any(fills)
    show_captions = not per_ring_ramps or len(set(fills)) < len(fills)
    bar_width = GRADIENT_BAR_WIDTH_RATIO * rect
    outer_caption_y = rect / 2
    # The rows start a line below ``outer`` when it is drawn; ``inner`` stands a line below the last row.
    lead = 1 if show_captions else 0
    rows = tuple(
        CircularRingRowLayout(key=entry.key, properties=entry.properties, y=outer_caption_y + ((index + lead) * line))
        for index, entry in enumerate(ring_entries)
    )
    inner_caption_y = outer_caption_y + ((len(ring_entries) + lead) * line)
    _, label_height = drawn_text_box("100%", font_family, font_size)
    caption_width = (
        max(float(drawn_text_box(text, font_family, font_size)[0]) for text in (RING_OUTER_CAPTION, RING_INNER_CAPTION))
        if show_captions
        else 0.0
    )
    names_width = max(caption_width, *(entry.text_width for entry in ring_entries))
    band_y: float | None
    names_x = 0.0
    if per_ring_ramps:
        bar_x = names_width + rect
        band_y = None
        scale_label_y = rows[-1].y + (rect / 2) + GRADIENT_TEXT_GAP
        width = bar_x + bar_width
        height = max(scale_label_y + float(label_height), inner_caption_y + (rect / 2) if show_captions else 0.0)
    else:
        width = max(names_width, bar_width)
        bar_x = 0.0
        if centred:
            names_x, bar_x = (width - names_width) / 2, (width - bar_width) / 2
        band_y = inner_caption_y + (RING_BAND_LINE_GAP * line)
        scale_label_y = band_y + (rect / 2) + GRADIENT_TEXT_GAP
        height = scale_label_y + float(label_height)
    min_label_text, max_label_text = gradient_scale_texts(ring_entries[0].properties)
    return CircularRingLegendLayout(
        width=width,
        height=height,
        show_captions=show_captions,
        outer_caption_y=outer_caption_y,
        inner_caption_y=inner_caption_y,
        rows=rows,
        names_x=names_x,
        bar_x=bar_x,
        bar_width=bar_width,
        band_y=band_y,
        min_label_text=min_label_text,
        max_label_text=max_label_text,
        scale_label_y=scale_label_y,
    )


def _gradient_reserved_width(
    gradient_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    color_rect_size: float,
    x_margin: float,
) -> float:
    width = GRADIENT_BAR_WIDTH_RATIO * float(color_rect_size)
    if len(gradient_entries) > 1:
        return width + (max(entry.text_width for entry in gradient_entries) + float(x_margin))
    # One gradient reserves its block, title included (OV-1611).
    return single_gradient_block_x(width, gradient_entries[0].text_width)[0]


def _build_gradient_layout(
    gradient_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    top: float,
    font_family: str,
    font_size: float,
    color_rect_size: float,
) -> CircularGradientLegendLayout | None:
    """The gradient rows, their first row ``top`` below the side block's origin."""

    if not gradient_entries:
        return None

    color_rect_size = float(color_rect_size)
    bar_width = GRADIENT_BAR_WIDTH_RATIO * color_rect_size
    row_height = legend_line_height(color_rect_size)
    # The scale labels share one line box, whatever their digits.
    _, label_height = drawn_text_box("100%", font_family, font_size)

    if len(gradient_entries) > 1:
        label_width = max(entry.text_width for entry in gradient_entries)
        bar_x = label_width + color_rect_size
        stack = stack_compact_gradient_rows(
            [gradient_scale_texts(entry.properties) for entry in gradient_entries],
            color_rect_size=color_rect_size,
            line_height=row_height,
            label_height=float(label_height),
            min_label_x=bar_x,
            max_label_x=bar_x + bar_width,
        )
        scales = tuple(
            GradientScaleLayout(
                scale.min_label_text, scale.max_label_text, scale.min_label_x, scale.max_label_x,
                top + scale.scale_label_y,
            )
            for scale in stack.scales
        )
        compact_entries = tuple(
            CircularCompactGradientEntryLayout(
                key=gradient_entries[index].key,
                properties=gradient_entries[index].properties,
                # The block's left edge; an int keeps the drawn translate(0,y).
                label_x=0,
                label_y=top + y,
                bar_y=top + y,
            )
            for index, y in stack.rows
        )
        return CircularGradientLegendLayout(
            compact=True,
            width=bar_x + bar_width,
            height=stack.height,
            bar_width=bar_width,
            bar_x=bar_x,
            min_label_text=scales[0].min_label_text,
            max_label_text=scales[0].max_label_text,
            min_label_x=scales[0].min_label_x,
            max_label_x=scales[0].max_label_x,
            scale_label_y=scales[0].scale_label_y,
            compact_entries=compact_entries,
            single_entries=(),
            more_scales=scales[1:],
        )

    entry = gradient_entries[0]
    title_width, title_height = legend_caption_box(entry.text, font_family, font_size)
    ys = single_gradient_row_ys(float(title_height), color_rect_size, float(label_height))
    width, bar_x = single_gradient_block_x(bar_width, float(title_width))
    min_label_text, max_label_text = gradient_scale_texts(entry.properties)
    return CircularGradientLegendLayout(
        compact=False,
        width=width,
        height=ys.height,
        bar_width=bar_width,
        bar_x=bar_x,
        min_label_text=min_label_text,
        max_label_text=max_label_text,
        min_label_x=bar_x,
        max_label_x=bar_x + bar_width,
        scale_label_y=top + ys.scale_label_y,
        compact_entries=(),
        single_entries=(
            CircularSingleGradientEntryLayout(
                key=entry.key,
                properties=entry.properties,
                title_x=width / 2.0,
                title_y=top,
                bar_x=bar_x,
                bar_y=top + ys.bar_y,
            ),
        ),
    )


def _build_side_block(
    ring_entries: tuple[_MeasuredLegendEntry, ...],
    gradient_entries: tuple[_MeasuredLegendEntry, ...],
    *,
    font_family: str,
    font_size: float,
    color_rect_size: float,
    line_margin: float,
    x_margin: float,
    centred: bool,
) -> _SideBlock | None:
    """The side block of the ring rows and the gradient rows, or ``None`` without either. ``centred`` (a vertical
    Legend) centres the ring block and the gradient rows on the side block's width."""

    ring = _build_ring_layout(
        ring_entries, font_family=font_family, font_size=font_size, color_rect_size=color_rect_size, centred=centred
    )
    rows_top = ring.height + (line_margin / 2.0) if ring is not None else 0.0
    rows = _build_gradient_layout(
        gradient_entries, top=rows_top, font_family=font_family, font_size=font_size, color_rect_size=color_rect_size
    )
    if ring is None and rows is None:
        return None
    lines = ring.lines if ring is not None else 0
    if rows is not None:
        # the rows and one scale line per group, or a title, its bar and its scale line
        lines += len(gradient_entries) + (len(rows.scales) if rows.compact else 2)
    width = max(ring.width if ring is not None else 0.0, rows.width if rows is not None else 0.0)
    if centred:
        ring = ring.shifted((width - ring.width) / 2) if ring is not None else None
        rows = rows.shifted((width - rows.width) / 2) if rows is not None else None
    return _SideBlock(
        ring=ring,
        rows=rows,
        width=width,
        height=rows_top + rows.height if rows is not None else ring.height if ring is not None else 0.0,
        reserved_width=max(
            ring.width if ring is not None else 0.0,
            _gradient_reserved_width(gradient_entries, color_rect_size=color_rect_size, x_margin=x_margin)
            if gradient_entries
            else 0.0,
        ),
        lines=lines,
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
    side_block: _SideBlock | None,
    *,
    pairwise_legend_width: float,
    canvas_width: float,
    color_rect_size: float,
    x_margin: float,
) -> float:
    solid_single_row_width = compensated_sum(entry.entry_width for entry in solid_entries)
    solid_desired_width = solid_single_row_width + (2 * float(x_margin) if solid_entries else 0.0)
    desired_width = solid_desired_width
    if side_block is not None:
        desired_width += max(float(pairwise_legend_width), float(side_block.width))
    if desired_width <= 0:
        return 0.0

    min_solid_width = (
        max((entry.entry_width + (2 * float(x_margin)) for entry in solid_entries), default=0.0)
    )
    min_width = min_solid_width
    if side_block is not None:
        min_width += max(float(pairwise_legend_width), float(side_block.width))
        if solid_entries:
            min_width += float(x_margin)
    if min_width <= 0:
        min_width = float(color_rect_size)
    return max(min_width, min(desired_width, float(canvas_width)))


def _build_horizontal_layout(
    solid_entries: tuple[_MeasuredLegendEntry, ...],
    side_block: _SideBlock | None,
    *,
    width: float,
    pairwise_legend_width: float,
    color_rect_size: float,
    line_margin: float,
    x_margin: float,
) -> tuple[tuple[CircularLegendEntryLayout, ...], float, float, float, int, int, int]:
    wrap_width = float(width)
    if side_block is not None:
        wrap_width = max(float(color_rect_size), wrap_width - float(pairwise_legend_width) - float(x_margin))
    if wrap_width <= 0:
        wrap_width = float("inf")

    rows = _wrap_horizontal_rows(solid_entries, wrap_width, x_margin=x_margin)
    feature_height = (
        float(color_rect_size) + ((len(rows) - 1) * float(line_margin))
        if rows
        else 0.0
    )
    gradient_height = float(side_block.height) if side_block is not None else 0.0
    feature_y_offset = max(0.0, (gradient_height - feature_height) / 2.0)

    laid_out_entries: list[CircularLegendEntryLayout] = []
    for row_index, row_entries in enumerate(rows):
        row_y = (
            feature_y_offset + (float(color_rect_size) / 2.0) + (row_index * float(line_margin))
            if side_block is not None
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

    num_lines = max(1, len(rows)) if solid_entries or side_block is not None else 0
    num_items_per_line = max((len(row) for row in rows), default=0)
    num_columns = num_items_per_line + (1 if side_block is not None else 0)
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
    solid_entries, gradient_entries, ring_entries = _measure_entries(
        legend_table,
        font_family=font_family,
        font_size=font_size,
        x_margin=x_margin,
    )
    horizontal = str(legend_position) in {"top", "bottom"}
    side_block = _build_side_block(
        ring_entries, gradient_entries, font_family=font_family, font_size=font_size,
        color_rect_size=float(color_rect_size), line_margin=line_margin, x_margin=x_margin, centred=not horizontal,
    )
    ring_layout, gradient_rows = (side_block.ring, side_block.rows) if side_block is not None else (None, None)
    resolved_pairwise_width = (
        float(pairwise_legend_width)
        if pairwise_legend_width is not None
        else side_block.reserved_width if side_block is not None else 0.0
    )

    if horizontal:
        width = (
            float(legend_width)
            if legend_width is not None
            else _horizontal_legend_width(
                solid_entries,
                side_block,
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
            side_block,
            width=width,
            pairwise_legend_width=resolved_pairwise_width,
            color_rect_size=float(color_rect_size),
            line_margin=line_margin,
            x_margin=x_margin,
        )
        gradient_x = 0.0
        gradient_y = 0.0
        content_bottom = feature_height
        if side_block is not None:
            gradient_y = max(0.0, (feature_height - float(side_block.height)) / 2.0)
            # Without feature rows the gradient block is the whole Legend (OV-1612).
            gradient_x = max(
                feature_width + x_margin if solid_entries else 0.0,
                width - float(side_block.width),
            )
            content_bottom = max(
                (max((entry.rect_y for entry in laid_out_entries), default=0.0) + float(color_rect_size) / 2.0)
                if laid_out_entries
                else 0.0,
                gradient_y + float(side_block.height),
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
            has_gradient=side_block is not None,
            num_lines=num_lines,
            num_columns=num_columns,
            num_items_per_line=num_items_per_line,
            solid_entries=laid_out_entries,
            gradient=ring_layout,
            gradient_x=gradient_x,
            gradient_y=gradient_y,
            gradient_rows=gradient_rows,
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
            float(side_block.width) if side_block is not None else 0.0,
        )
    )
    # With a side block, every block of a vertical Legend is centred on one axis (OV-2835): the feature rows here,
    # the side block's blocks on its own width, and the side block on the Legend's.
    alignment_width = max(width, feature_block_width, float(side_block.width) if side_block is not None else 0.0)
    feature_x = (alignment_width - feature_block_width) / 2.0 if side_block is not None and solid_entries else 0.0
    laid_out_entries = tuple(
        CircularLegendEntryLayout(
            key=entry.key,
            properties=entry.properties,
            text_width=entry.text_width,
            entry_width=entry.entry_width,
            rect_x=feature_x,
            rect_y=idx * line_margin,
            text_x=feature_x + x_margin,
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
    if side_block is not None:
        gradient_x = (alignment_width - float(side_block.width)) / 2.0
        gradient_y = len(solid_entries) * line_margin + (
            line_margin * 0.5 if solid_entries else 0.0
        )
    content_bottom = feature_height
    if side_block is not None:
        content_bottom = max(content_bottom, gradient_y + float(side_block.height))
        height = content_bottom + (float(color_rect_size) / 2.0)
    else:
        height = feature_height
    num_lines = len(solid_entries) + (side_block.lines if side_block is not None else 0)
    return CircularLegendLayout(
        horizontal=False,
        width=width,
        height=height,
        feature_width=feature_block_width,
        feature_height=feature_height,
        pairwise_legend_width=resolved_pairwise_width,
        line_margin=line_margin,
        x_margin=x_margin,
        has_gradient=side_block is not None,
        num_lines=num_lines,
        num_columns=1,
        num_items_per_line=1 if solid_entries else 0,
        solid_entries=laid_out_entries,
        gradient=ring_layout,
        gradient_x=gradient_x,
        gradient_y=gradient_y,
        gradient_rows=gradient_rows,
    )


__all__ = [
    "CircularCompactGradientEntryLayout",
    "CircularGradientLegendLayout",
    "CircularLegendEntryLayout",
    "CircularLegendLayout",
    "CircularRingLegendLayout",
    "CircularRingRowLayout",
    "CircularSingleGradientEntryLayout",
    "RING_BAND_LABEL",
    "RING_INNER_CAPTION",
    "RING_OUTER_CAPTION",
    "build_circular_legend_layout",
]
