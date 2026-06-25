#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.shapes import Circle, Line
from svgwrite.text import Text

from ....analysis.gc import gc_content_percent_df
from ....layout.circular_depth_axis import (
    DEPTH_AXIS_SMALL_TICK_SIZE_PX,
    DEPTH_AXIS_STROKE_WIDTH_PX,
    DEPTH_AXIS_TICK_SIZE_PX,
    depth_axis_tick_font_size_px,
)
from ....layout.scalar_axis import (
    format_percent_tick,
    scalar_axis_small_tick_values,
    scalar_axis_tick_values,
    scaled_scalar_fraction,
)
from ....svg.circular_tracks import (  # type: ignore[reportMissingImports]
    generate_circular_gc_content_path_desc,
    generate_circular_scalar_area_path_desc,
)


class GcContentDrawer:
    """
    Draws GC content track on a circular canvas.
    """

    def __init__(self, gc_config) -> None:
        self._gc_config = gc_config
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity
        self.mode: str = str(getattr(gc_config, "mode", "deviation"))
        min_percent = getattr(gc_config, "min_percent", None)
        max_percent = getattr(gc_config, "max_percent", None)
        self.min_percent: float = float(0.0 if min_percent is None else min_percent)
        self.max_percent: float = float(100.0 if max_percent is None else max_percent)
        self.show_axis: bool = bool(getattr(gc_config, "show_axis", True))
        self.show_ticks: bool = bool(getattr(gc_config, "show_ticks", True))
        self.large_tick_interval: float | None = getattr(gc_config, "large_tick_interval", 20.0)
        self.small_tick_interval: float | None = getattr(gc_config, "small_tick_interval", None)
        self.tick_font_size: float | None = getattr(gc_config, "tick_font_size", None)
        self.percent_background_color: str = str(getattr(gc_config, "percent_background_color", "#737373"))
        self.percent_background_opacity: float = float(getattr(gc_config, "percent_background_opacity", 1.0))
        self.percent_border_color: str = str(getattr(gc_config, "percent_border_color", "#4b5563"))
        self.percent_border_width: float = float(getattr(gc_config, "percent_border_width", 0.8))
        self.axis_stroke: str = "#4b5563"
        self.axis_stroke_width: float = DEPTH_AXIS_STROKE_WIDTH_PX
        self.axis_tick_size: float = DEPTH_AXIS_TICK_SIZE_PX
        self.axis_small_tick_size: float = DEPTH_AXIS_SMALL_TICK_SIZE_PX

    def _tick_font_size(self, track_width: float) -> float:
        return depth_axis_tick_font_size_px(self._gc_config, track_width)

    @staticmethod
    def _percent_radii(radius: float, track_width: float, norm_factor: float) -> tuple[float, float]:
        baseline_radius = max(
            0.0,
            (float(radius) * float(norm_factor)) - (0.5 * float(track_width)),
        )
        return baseline_radius, baseline_radius + float(track_width)

    def _add_percent_container(
        self,
        group: Group,
        radius: float,
        track_width: float,
        norm_factor: float,
        group_id: str,
    ) -> Group:
        baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
        middle_radius = baseline_radius + (0.5 * float(track_width))
        group.add(
            Circle(
                id=f"{group_id}_percent_background",
                center=(0, 0),
                r=middle_radius,
                fill="none",
                stroke=self.percent_background_color,
                stroke_width=float(track_width),
                stroke_opacity=self.percent_background_opacity,
            )
        )
        return group

    def _add_percent_border(
        self,
        group: Group,
        baseline_radius: float,
        outer_radius: float,
        group_id: str,
    ) -> Group:
        if self.percent_border_width <= 0:
            return group
        border_group = Group(id=f"{group_id}_percent_border")
        for border_radius in (baseline_radius, outer_radius):
            border_group.add(
                Circle(
                    center=(0, 0),
                    r=border_radius,
                    fill="none",
                    stroke=self.percent_border_color,
                    stroke_width=self.percent_border_width,
                )
            )
        group.add(border_group)
        return group

    def _add_percent_axis(
        self,
        group: Group,
        radius: float,
        track_width: float,
        norm_factor: float,
        axis_group_id: str = "gc_content_axis",
    ) -> Group:
        if not self.show_axis:
            return group

        baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
        axis_group = Group(id=axis_group_id)
        axis_group.add(
            Circle(
                center=(0, 0),
                r=baseline_radius,
                fill="none",
                stroke=self.axis_stroke,
                stroke_width=self.axis_stroke_width,
            )
        )
        axis_group.add(
            Line(
                start=(0, -baseline_radius),
                end=(0, -outer_radius),
                stroke=self.axis_stroke,
                stroke_width=self.axis_stroke_width,
            )
        )

        font_size = self._tick_font_size(track_width)
        large_tick_values = scalar_axis_tick_values(
            self.min_percent,
            self.max_percent,
            show_ticks=self.show_ticks,
            large_tick_interval=self.large_tick_interval,
        )
        for tick_value in scalar_axis_small_tick_values(
            self.min_percent,
            self.max_percent,
            show_ticks=self.show_ticks,
            small_tick_interval=self.small_tick_interval,
            large_tick_values=large_tick_values,
        ):
            normalized = scaled_scalar_fraction(tick_value, self.min_percent, self.max_percent)
            tick_radius = baseline_radius + (float(track_width) * normalized)
            axis_group.add(
                Line(
                    start=(0, -tick_radius),
                    end=(self.axis_small_tick_size, -tick_radius),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
        for tick_value in large_tick_values:
            normalized = scaled_scalar_fraction(tick_value, self.min_percent, self.max_percent)
            tick_radius = baseline_radius + (float(track_width) * normalized)
            axis_group.add(
                Line(
                    start=(0, -tick_radius),
                    end=(self.axis_tick_size, -tick_radius),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
            axis_group.add(
                Text(
                    format_percent_tick(tick_value),
                    insert=(self.axis_tick_size + 1.5, -tick_radius),
                    fill=self.axis_stroke,
                    font_size=font_size,
                    text_anchor="start",
                    dominant_baseline="central",
                )
            )

        group.add(axis_group)
        return group

    def draw(
        self,
        radius: float,
        group: Group,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
        group_id: str = "gc_content",
    ) -> Group:
        if self.mode == "percent":
            plot_df = gc_content_percent_df(
                gc_df,
                dinucleotide=dinucleotide,
                min_percent=self.min_percent,
                max_percent=self.max_percent,
            )
            gc_path_desc: str = generate_circular_scalar_area_path_desc(
                radius,
                record_len,
                plot_df,
                track_width,
                norm_factor,
                close_at_record_len=True,
            )
        else:
            gc_path_desc = generate_circular_gc_content_path_desc(
                radius, record_len, gc_df, track_width, norm_factor, dinucleotide
            )
        if not gc_path_desc:
            return group
        if self.mode == "percent":
            self._add_percent_container(group, radius, track_width, norm_factor, group_id)
        gc_path: Path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd",
        )
        group.add(gc_path)
        if self.mode == "percent":
            baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
            self._add_percent_border(group, baseline_radius, outer_radius, group_id)
            self._add_percent_axis(group, radius, track_width, norm_factor, f"{group_id}_axis")
        return group


__all__ = ["GcContentDrawer"]


