from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

from ..core.text import calculate_bbox_dimensions
from .circular_layout import min_gradient_label_text


LINE_MARGIN_RATIO = 24 / 14
X_MARGIN_RATIO = 22 / 14
GRADIENT_BAR_WIDTH_RATIO = 10
GRADIENT_LABEL_GAP_RATIO = 0.2


@dataclass(frozen=True, slots=True)
class LinearLegendEntryLayout:
    key: str
    properties: Mapping[str, object]
    rect_x: float
    rect_y: float
    text_x: float
    text_y: float


@dataclass(frozen=True, slots=True)
class LinearFeatureLegendLayout:
    entries: tuple[LinearLegendEntryLayout, ...]
    width: float
    height: float
    num_lines: int


@dataclass(frozen=True, slots=True)
class LinearGradientEntryLayout:
    key: str
    properties: Mapping[str, object]
    title_x: float
    title_y: float
    bar_x: float
    bar_y: float


@dataclass(frozen=True, slots=True)
class LinearGradientLegendLayout:
    compact: bool
    entries: tuple[LinearGradientEntryLayout, ...]
    width: float
    height: float
    bar_width: float
    min_label_text: str
    min_label_x: float
    max_label_x: float
    scale_label_y: float


@dataclass(frozen=True, slots=True)
class LinearOrientationLegendLayout:
    feature: LinearFeatureLegendLayout
    gradient: LinearGradientLegendLayout | None
    feature_x: float
    feature_y: float
    gradient_x: float
    gradient_y: float
    width: float
    height: float


@dataclass(frozen=True, slots=True)
class LinearLegendLayout:
    horizontal: LinearOrientationLegendLayout
    vertical: LinearOrientationLegendLayout
    active_orientation: Literal["horizontal", "vertical"]

    @property
    def active(self) -> LinearOrientationLegendLayout:
        return (
            self.horizontal
            if self.active_orientation == "horizontal"
            else self.vertical
        )

    @property
    def has_gradient(self) -> bool:
        return self.horizontal.gradient is not None


@dataclass(frozen=True, slots=True)
class _MeasuredSolidEntry:
    key: str
    properties: Mapping[str, object]
    text_width: float


def _measure_solid_entries(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    font_family: str,
    font_size: float,
    dpi: int,
) -> tuple[_MeasuredSolidEntry, ...]:
    return tuple(
        _MeasuredSolidEntry(
            key=str(key),
            properties=properties,
            text_width=calculate_bbox_dimensions(
                str(key),
                font_family,
                font_size,
                dpi,
            )[0],
        )
        for key, properties in legend_table.items()
        if properties.get("type") == "solid"
    )


def _build_gradient_layout(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    font_family: str,
    font_size: float,
    dpi: int,
    color_rect_size: float,
    line_height: float,
) -> LinearGradientLegendLayout | None:
    gradient_entries = tuple(
        (str(key), properties)
        for key, properties in legend_table.items()
        if properties.get("type") == "gradient"
    )
    if not gradient_entries:
        return None

    bar_width = GRADIENT_BAR_WIDTH_RATIO * color_rect_size
    min_label_text = min_gradient_label_text(
        gradient_entries[0][1].get("min_value", 0)
    )
    if len(gradient_entries) > 1:
        label_width = max(
            calculate_bbox_dimensions(
                key,
                font_family,
                font_size,
                dpi,
            )[0]
            for key, _ in gradient_entries
        )
        bar_x = label_width + (GRADIENT_LABEL_GAP_RATIO * color_rect_size)
        entries = tuple(
            LinearGradientEntryLayout(
                key=key,
                properties=properties,
                title_x=0.0,
                title_y=(color_rect_size / 2.0) + (index * line_height),
                bar_x=bar_x,
                bar_y=(color_rect_size / 2.0) + (index * line_height),
            )
            for index, (key, properties) in enumerate(gradient_entries)
        )
        scale_label_y = (
            color_rect_size
            + ((len(gradient_entries) - 1) * line_height)
            + 2.0
        )
        scale_label_height = max(
            calculate_bbox_dimensions(
                min_label_text,
                font_family,
                font_size,
                dpi,
            )[1],
            calculate_bbox_dimensions(
                "100%",
                font_family,
                font_size,
                dpi,
            )[1],
        )
        return LinearGradientLegendLayout(
            compact=True,
            entries=entries,
            width=bar_x + bar_width,
            height=scale_label_y + scale_label_height,
            bar_width=bar_width,
            min_label_text=min_label_text,
            min_label_x=bar_x,
            max_label_x=bar_x + bar_width,
            scale_label_y=scale_label_y,
        )

    key, properties = gradient_entries[0]
    title_height = calculate_bbox_dimensions(
        key,
        font_family,
        font_size,
        dpi,
    )[1]
    bar_y = title_height + (color_rect_size / 2.0)
    scale_label_y = bar_y + (color_rect_size / 2.0) + 2.0
    scale_label_height = max(
        calculate_bbox_dimensions(
            min_label_text,
            font_family,
            font_size,
            dpi,
        )[1],
        calculate_bbox_dimensions(
            "100%",
            font_family,
            font_size,
            dpi,
        )[1],
    )
    entry = LinearGradientEntryLayout(
        key=key,
        properties=properties,
        title_x=bar_width / 2.0,
        title_y=0.0,
        bar_x=0.0,
        bar_y=bar_y,
    )
    return LinearGradientLegendLayout(
        compact=False,
        entries=(entry,),
        width=bar_width,
        height=scale_label_y + scale_label_height,
        bar_width=bar_width,
        min_label_text=min_label_text,
        min_label_x=0.0,
        max_label_x=bar_width,
        scale_label_y=scale_label_y,
    )


def _build_horizontal_feature_layout(
    entries: tuple[_MeasuredSolidEntry, ...],
    *,
    canvas_width: float,
    gradient_width: float,
    color_rect_size: float,
    line_height: float,
    text_x_offset: float,
) -> LinearFeatureLegendLayout:
    wrap_width = float(canvas_width)
    if gradient_width > 0.0 and wrap_width > 0.0:
        min_width = color_rect_size + (2.0 * text_x_offset)
        wrap_width = max(
            wrap_width - gradient_width - text_x_offset,
            min_width,
        )

    y_offset = color_rect_size / 2.0
    current_x = 0.0
    current_row_width = 0.0
    max_row_width = 0.0
    height = line_height
    num_lines = 1 if entries else 0
    laid_out: list[LinearLegendEntryLayout] = []
    for entry in entries:
        entry_width = (
            color_rect_size
            + text_x_offset
            + entry.text_width
            + text_x_offset
        )
        if wrap_width > 0.0 and current_x + entry_width > wrap_width:
            max_row_width = max(max_row_width, current_row_width)
            current_x = 0.0
            current_row_width = 0.0
            y_offset += line_height
            height += line_height
            num_lines += 1
        laid_out.append(
            LinearLegendEntryLayout(
                key=entry.key,
                properties=entry.properties,
                rect_x=current_x,
                rect_y=y_offset,
                text_x=current_x + text_x_offset,
                text_y=y_offset,
            )
        )
        current_x += entry_width
        current_row_width += entry_width

    return LinearFeatureLegendLayout(
        entries=tuple(laid_out),
        width=max(max_row_width, current_row_width),
        height=height,
        num_lines=num_lines,
    )


def _build_vertical_feature_layout(
    entries: tuple[_MeasuredSolidEntry, ...],
    *,
    color_rect_size: float,
    line_height: float,
    text_x_offset: float,
) -> LinearFeatureLegendLayout:
    y_offset = color_rect_size / 2.0
    max_width = 0.0
    laid_out: list[LinearLegendEntryLayout] = []
    for entry in entries:
        laid_out.append(
            LinearLegendEntryLayout(
                key=entry.key,
                properties=entry.properties,
                rect_x=0.0,
                rect_y=y_offset,
                text_x=text_x_offset,
                text_y=y_offset,
            )
        )
        max_width = max(max_width, text_x_offset + entry.text_width)
        y_offset += line_height
    return LinearFeatureLegendLayout(
        entries=tuple(laid_out),
        width=max_width,
        height=y_offset,
        num_lines=len(entries),
    )


def build_linear_legend_layout(
    legend_table: Mapping[object, Mapping[str, object]],
    *,
    legend_position: str,
    canvas_width: float,
    font_family: str,
    font_size: float,
    dpi: int,
    color_rect_size: float,
) -> LinearLegendLayout:
    active_orientation: Literal["horizontal", "vertical"] = (
        "horizontal"
        if legend_position in {"top", "bottom"}
        else "vertical"
    )
    if not legend_table:
        empty_feature = LinearFeatureLegendLayout(
            entries=(),
            width=0.0,
            height=0.0,
            num_lines=0,
        )
        empty_orientation = LinearOrientationLegendLayout(
            feature=empty_feature,
            gradient=None,
            feature_x=0.0,
            feature_y=0.0,
            gradient_x=0.0,
            gradient_y=0.0,
            width=0.0,
            height=0.0,
        )
        return LinearLegendLayout(
            horizontal=empty_orientation,
            vertical=empty_orientation,
            active_orientation=active_orientation,
        )

    line_height = LINE_MARGIN_RATIO * color_rect_size
    text_x_offset = X_MARGIN_RATIO * color_rect_size
    solid_entries = _measure_solid_entries(
        legend_table,
        font_family=font_family,
        font_size=font_size,
        dpi=dpi,
    )
    gradient = _build_gradient_layout(
        legend_table,
        font_family=font_family,
        font_size=font_size,
        dpi=dpi,
        color_rect_size=color_rect_size,
        line_height=line_height,
    )
    gradient_width = gradient.width if gradient is not None else 0.0
    horizontal_feature = _build_horizontal_feature_layout(
        solid_entries,
        canvas_width=canvas_width,
        gradient_width=gradient_width,
        color_rect_size=color_rect_size,
        line_height=line_height,
        text_x_offset=text_x_offset,
    )
    vertical_feature = _build_vertical_feature_layout(
        solid_entries,
        color_rect_size=color_rect_size,
        line_height=line_height,
        text_x_offset=text_x_offset,
    )

    if gradient is None:
        horizontal = LinearOrientationLegendLayout(
            feature=horizontal_feature,
            gradient=None,
            feature_x=0.0,
            feature_y=0.0,
            gradient_x=0.0,
            gradient_y=0.0,
            width=horizontal_feature.width,
            height=horizontal_feature.height,
        )
        vertical = LinearOrientationLegendLayout(
            feature=vertical_feature,
            gradient=None,
            feature_x=0.0,
            feature_y=0.0,
            gradient_x=0.0,
            gradient_y=0.0,
            width=vertical_feature.width,
            height=vertical_feature.height,
        )
    else:
        horizontal = LinearOrientationLegendLayout(
            feature=horizontal_feature,
            gradient=gradient,
            feature_x=0.0,
            feature_y=max(
                0.0,
                (gradient.height - horizontal_feature.height) / 2.0,
            ),
            gradient_x=horizontal_feature.width + text_x_offset,
            gradient_y=max(
                0.0,
                (horizontal_feature.height - gradient.height) / 2.0,
            ),
            width=horizontal_feature.width + text_x_offset + gradient.width,
            height=max(horizontal_feature.height, gradient.height),
        )
        vertical = LinearOrientationLegendLayout(
            feature=vertical_feature,
            gradient=gradient,
            feature_x=max(
                0.0,
                (gradient.width - vertical_feature.width) / 2.0,
            ),
            feature_y=0.0,
            gradient_x=0.0,
            gradient_y=vertical_feature.height + (line_height / 2.0),
            width=max(vertical_feature.width, gradient.width),
            height=(
                vertical_feature.height
                + (line_height / 2.0)
                + gradient.height
            ),
        )
    return LinearLegendLayout(
        horizontal=horizontal,
        vertical=vertical,
        active_orientation=active_orientation,
    )


__all__ = [
    "LinearFeatureLegendLayout",
    "LinearGradientEntryLayout",
    "LinearGradientLegendLayout",
    "LinearLegendEntryLayout",
    "LinearLegendLayout",
    "LinearOrientationLegendLayout",
    "build_linear_legend_layout",
]
