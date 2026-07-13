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
from ....layout.scalar_axis import (
    _depth_axis_bounds,
    _format_depth_tick,
    _prepare_depth_plot_dataframe,
    _scaled_depth_fraction,
    linear_scalar_axis_tick_font_size_px,
    scalar_axis_small_tick_values,
    scalar_axis_tick_values,
)
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
        axis_group_id: str = "depth_axis",
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
        axis_group = Group(id=axis_group_id)
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

        axis_min, axis_max = _depth_axis_bounds(depth_df, self.min_depth, self.max_depth)
        font_size = linear_scalar_axis_tick_font_size_px(self, track_height)
        large_tick_values = scalar_axis_tick_values(
            axis_min,
            axis_max,
            show_ticks=self.show_ticks,
            large_tick_interval=self.large_tick_interval,
        )
        small_tick_values = scalar_axis_small_tick_values(
            axis_min,
            axis_max,
            show_ticks=self.show_ticks,
            small_tick_interval=self.small_tick_interval,
            large_tick_values=large_tick_values,
        )
        for tick_value in small_tick_values:
            normalized = _scaled_depth_fraction(
                tick_value, axis_min, axis_max, normalize=self.normalize
            )
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
            normalized = _scaled_depth_fraction(
                tick_value, axis_min, axis_max, normalize=self.normalize
            )
            tick_y = baseline_y - (float(track_height) * normalized)
            axis_group.add(
                Line(
                    start=(start_x - self.axis_tick_size, tick_y),
                    end=(start_x, tick_y),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
            tick_label = (
                None
                if self.min_depth is None and math.isclose(tick_value, axis_min, rel_tol=1e-9, abs_tol=1e-9)
                else _format_depth_tick(tick_value)
            )
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
        axis_group_id: str = "depth_axis",
    ) -> Group:
        plot_df = _prepare_depth_plot_dataframe(
            depth_df,
            self.min_depth,
            self.max_depth,
            normalize=self.normalize,
        )
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
            axis_group_id,
        )
        return group


__all__ = ["DepthDrawer"]
