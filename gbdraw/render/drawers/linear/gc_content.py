#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.shapes import Line, Rect
from svgwrite.text import Text

from ....analysis.gc import gc_content_percent_df
from ....layout.linear_coords import normalize_position_to_linear_track
from ....layout.scalar_axis import (
    format_percent_tick,
    linear_scalar_axis_tick_font_size_px,
    scalar_axis_small_tick_values,
    scalar_axis_tick_values,
    scaled_scalar_fraction,
)
from ....svg.linear_tracks import calculate_gc_content_path_desc, calculate_linear_scalar_area_path_desc


class GcContentDrawer:
    """
    Draws GC content track on a linear canvas.
    """

    def __init__(self, gc_config):
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
        self.axis_stroke_width: float = 0.8
        self.axis_tick_size: float = 3.0
        self.axis_small_tick_size: float = 2.0

    def _tick_font_size(self, track_height: float) -> float:
        return linear_scalar_axis_tick_font_size_px(self._gc_config, track_height)

    def _add_percent_container(
        self,
        group: Group,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
        group_id: str,
    ) -> Group:
        final_x = normalize_position_to_linear_track(
            record_len,
            record_len,
            alignment_width,
            genome_size_normalization_factor,
        )
        width = max(0.0, float(final_x) - float(start_x))
        group.add(
            Rect(
                id=f"{group_id}_percent_background",
                insert=(start_x, start_y),
                size=(width, track_height),
                fill=self.percent_background_color,
                fill_opacity=self.percent_background_opacity,
                stroke="none",
            )
        )
        return group

    def _add_percent_border(
        self,
        group: Group,
        start_x: float,
        start_y: float,
        width: float,
        track_height: float,
        group_id: str,
    ) -> Group:
        if self.percent_border_width <= 0:
            return group
        group.add(
            Rect(
                id=f"{group_id}_percent_border",
                insert=(start_x, start_y),
                size=(width, track_height),
                fill="none",
                stroke=self.percent_border_color,
                stroke_width=self.percent_border_width,
            )
        )
        return group

    def _add_percent_axis(
        self,
        group: Group,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
        axis_group_id: str = "gc_content_axis",
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

        font_size = self._tick_font_size(track_height)
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
            normalized = scaled_scalar_fraction(tick_value, self.min_percent, self.max_percent)
            tick_y = baseline_y - (float(track_height) * normalized)
            axis_group.add(
                Line(
                    start=(start_x - self.axis_tick_size, tick_y),
                    end=(start_x, tick_y),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
            axis_group.add(
                Text(
                    format_percent_tick(tick_value),
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
        gc_df: DataFrame,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
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
            gc_path_desc: str = calculate_linear_scalar_area_path_desc(
                start_x,
                start_y,
                plot_df,
                record_len,
                alignment_width,
                genome_size_normalization_factor,
                track_height,
            )
        else:
            gc_path_desc = calculate_gc_content_path_desc(
                start_x,
                start_y,
                gc_df,
                record_len,
                alignment_width,
                genome_size_normalization_factor,
                track_height,
                dinucleotide,
            )
        if not gc_path_desc:
            return group
        if self.mode == "percent":
            self._add_percent_container(
                group,
                record_len,
                alignment_width,
                genome_size_normalization_factor,
                track_height,
                start_x,
                start_y,
                group_id,
            )
        gc_path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd",
        )
        group.add(gc_path)
        if self.mode == "percent":
            final_x = normalize_position_to_linear_track(
                record_len,
                record_len,
                alignment_width,
                genome_size_normalization_factor,
            )
            self._add_percent_border(
                group,
                start_x,
                start_y,
                max(0.0, float(final_x) - float(start_x)),
                track_height,
                group_id,
            )
            self._add_percent_axis(
                group,
                record_len,
                alignment_width,
                genome_size_normalization_factor,
                track_height,
                start_x,
                start_y,
                f"{group_id}_axis",
            )
        return group


__all__ = ["GcContentDrawer"]


