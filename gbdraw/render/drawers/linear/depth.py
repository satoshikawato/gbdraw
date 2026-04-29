#!/usr/bin/env python
# coding: utf-8

import math

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.shapes import Line
from svgwrite.text import Text

from ....configurators import DepthConfigurator
from ....layout.linear_coords import normalize_position_to_linear_track
from ....svg.linear_tracks import calculate_depth_path_desc


class DepthDrawer:
    """Draws a depth coverage track on a linear canvas."""

    def __init__(self, depth_config: DepthConfigurator) -> None:
        self.fill_color: str = depth_config.fill_color
        self.stroke_color: str = depth_config.stroke_color
        self.stroke_width: float = depth_config.stroke_width
        self.fill_opacity: float = depth_config.fill_opacity
        self.min_depth: float | None = depth_config.min_depth
        self.max_depth: float | None = depth_config.max_depth
        self.normalize: bool = depth_config.normalize
        self.show_axis: bool = depth_config.show_axis
        self.show_ticks: bool = depth_config.show_ticks
        self.large_tick_interval: float | None = depth_config.large_tick_interval
        self.tick_interval: float | None = self.large_tick_interval
        self.small_tick_interval: float | None = depth_config.small_tick_interval
        self.tick_font_size: float | None = depth_config.tick_font_size
        self.axis_stroke: str = "#4b5563"
        self.axis_stroke_width: float = 0.8
        self.axis_tick_size: float = 3.0
        self.axis_small_tick_size: float = 2.0

    @staticmethod
    def _format_depth_tick(value: float) -> str:
        value = 0.0 if not math.isfinite(value) else float(value)
        if abs(value - round(value)) < 0.05 or abs(value) >= 100.0:
            return f"{value:.0f}x"
        if abs(value) >= 10.0:
            return f"{value:.1f}x"
        return f"{value:.2f}".rstrip("0").rstrip(".") + "x"

    def _format_depth_tick_label(self, value: float, axis_min: float) -> str | None:
        if self.min_depth is None and math.isclose(value, axis_min, rel_tol=1e-9, abs_tol=1e-9):
            return None
        return self._format_depth_tick(value)

    def _axis_bounds(self, depth_df: DataFrame) -> tuple[float, float]:
        axis_min = float(self.min_depth) if self.min_depth is not None else 0.0
        if self.max_depth is not None:
            axis_max = float(self.max_depth)
        elif "depth" in depth_df.columns and not depth_df.empty:
            axis_max = float(depth_df["depth"].max())
        else:
            axis_max = axis_min
        if not math.isfinite(axis_max):
            axis_max = axis_min
        return axis_min, max(axis_min, axis_max)

    def _scale_depth_value(self, value: float) -> float:
        value = 0.0 if not math.isfinite(value) else float(value)
        if not self.normalize:
            return value
        return math.log10(value) + 1.0 if value > 0 else 0.0

    def _scaled_fraction(self, value: float, axis_min: float, axis_max: float) -> float:
        scaled_min = self._scale_depth_value(axis_min)
        scaled_max = self._scale_depth_value(axis_max)
        if scaled_max <= scaled_min:
            return 0.0
        scaled_value = self._scale_depth_value(value)
        return max(0.0, min(1.0, (scaled_value - scaled_min) / (scaled_max - scaled_min)))

    @staticmethod
    def _has_tick_value(values: list[float], candidate: float) -> bool:
        return any(math.isclose(candidate, value, rel_tol=1e-9, abs_tol=1e-9) for value in values)

    def _tick_values(self, axis_min: float, axis_max: float) -> list[float]:
        if not self.show_ticks:
            return []
        if axis_max <= axis_min:
            return [axis_min]
        if self.large_tick_interval is None:
            return [axis_min, axis_max]

        interval = float(self.large_tick_interval)
        if interval <= 0:
            return [axis_min, axis_max]
        tick_values: list[float] = []
        start = math.ceil((axis_min - 1e-9) / interval) * interval
        if start > axis_min + 1e-9:
            tick_values.append(axis_min)
        tick = start
        while tick <= axis_max + 1e-9 and len(tick_values) < 100:
            if tick >= axis_min - 1e-9:
                tick_values.append(0.0 if abs(tick) < 1e-9 else float(tick))
            tick += interval
        if not tick_values or not math.isclose(tick_values[-1], axis_max, rel_tol=1e-9, abs_tol=1e-9):
            tick_values.append(axis_max)

        deduplicated: list[float] = []
        for value in tick_values:
            if not deduplicated or not math.isclose(value, deduplicated[-1], rel_tol=1e-9, abs_tol=1e-9):
                deduplicated.append(value)
        return deduplicated

    def _small_tick_values(
        self,
        axis_min: float,
        axis_max: float,
        large_tick_values: list[float],
    ) -> list[float]:
        if not self.show_ticks or self.small_tick_interval is None or axis_max <= axis_min:
            return []
        interval = float(self.small_tick_interval)
        if interval <= 0:
            return []

        tick_values: list[float] = []
        tick = math.ceil((axis_min - 1e-9) / interval) * interval
        while tick <= axis_max + 1e-9 and len(tick_values) < 500:
            normalized_tick = 0.0 if abs(tick) < 1e-9 else float(tick)
            if (
                normalized_tick >= axis_min - 1e-9
                and not self._has_tick_value(large_tick_values, normalized_tick)
                and not self._has_tick_value(tick_values, normalized_tick)
            ):
                tick_values.append(normalized_tick)
            tick += interval
        return tick_values

    def _plot_depth_df(self, depth_df: DataFrame) -> DataFrame:
        plot_df = depth_df.copy()
        axis_min, axis_max = self._axis_bounds(plot_df)
        if axis_max <= axis_min:
            plot_df["depth_normalized"] = 0.0
            return plot_df
        plot_df["depth_normalized"] = plot_df["depth"].map(
            lambda value: self._scaled_fraction(float(value), axis_min, axis_max)
        )
        return plot_df

    def _add_axes(
        self,
        group: Group,
        depth_df: DataFrame,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
    ) -> Group:
        if not self.show_axis:
            return group
        baseline_y = float(start_y) + float(track_height)
        final_x = normalize_position_to_linear_track(
            record_len,
            record_len,
            alignment_width,
            genome_size_normalization_factor,
        )
        axis_group = Group(id="depth_axis")
        axis_group.add(
            Line(
                start=(start_x, baseline_y),
                end=(final_x, baseline_y),
                stroke=self.axis_stroke,
                stroke_width=self.axis_stroke_width,
            )
        )
        axis_group.add(
            Line(
                start=(start_x, start_y),
                end=(start_x, baseline_y),
                stroke=self.axis_stroke,
                stroke_width=self.axis_stroke_width,
            )
        )

        axis_min, axis_max = self._axis_bounds(depth_df)
        font_size = float(self.tick_font_size) if self.tick_font_size is not None else max(5.0, min(8.0, float(track_height) * 0.7))
        large_tick_values = self._tick_values(axis_min, axis_max)
        for tick_value in self._small_tick_values(axis_min, axis_max, large_tick_values):
            normalized = self._scaled_fraction(tick_value, axis_min, axis_max)
            tick_y = baseline_y - (float(track_height) * normalized)
            axis_group.add(
                Line(
                    start=(start_x - self.axis_small_tick_size, tick_y),
                    end=(start_x, tick_y),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
        for tick_value in large_tick_values:
            normalized = self._scaled_fraction(tick_value, axis_min, axis_max)
            tick_y = baseline_y - (float(track_height) * normalized)
            axis_group.add(
                Line(
                    start=(start_x - self.axis_tick_size, tick_y),
                    end=(start_x, tick_y),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
            tick_label = self._format_depth_tick_label(tick_value, axis_min)
            if tick_label is not None:
                axis_group.add(
                    Text(
                        tick_label,
                        insert=(start_x - self.axis_tick_size - 1.5, tick_y),
                        fill=self.axis_stroke,
                        font_size=font_size,
                        text_anchor="end",
                        dominant_baseline="central",
                    )
                )

        group.add(axis_group)
        return group

    def draw(
        self,
        group: Group,
        depth_df: DataFrame,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
    ) -> Group:
        plot_df = self._plot_depth_df(depth_df)
        depth_path_desc = calculate_depth_path_desc(
            start_x,
            start_y,
            plot_df,
            record_len,
            alignment_width,
            genome_size_normalization_factor,
            track_height,
        )
        if not depth_path_desc:
            return group
        depth_path = Path(
            d=depth_path_desc,
            fill=self.fill_color,
            stroke=self.stroke_color,
            stroke_width=self.stroke_width,
            fill_opacity=self.fill_opacity,
            fill_rule="evenodd",
        )
        group.add(depth_path)
        self._add_axes(
            group,
            depth_df,
            record_len,
            alignment_width,
            genome_size_normalization_factor,
            track_height,
            start_x,
            start_y,
        )
        return group


__all__ = ["DepthDrawer"]
