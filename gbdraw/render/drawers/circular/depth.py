#!/usr/bin/env python
# coding: utf-8

import math

from pandas import DataFrame
from gbdraw.svg.elements import Circle, Group, Line, Path, Text

from ....configurators import DepthConfigurator
from ....layout.circular_depth_axis import (
    DEPTH_AXIS_SMALL_TICK_SIZE_PX,
    DEPTH_AXIS_STROKE_WIDTH_PX,
    DEPTH_AXIS_TICK_SIZE_PX,
    SCALAR_AXIS_LABEL_GAP_PX,
    ScalarAxisTick,
    depth_axis_tick_font_size_px,
    scalar_axis_ink,
)
from ....layout.circular_ink import CircularInk, CircularInkProfile, stroke_pad_px
from ....layout.scalar_axis import (
    _depth_axis_bounds,
    _format_depth_tick,
    _prepare_depth_plot_dataframe,
    _scaled_depth_fraction,
    scalar_axis_small_tick_values,
    scalar_axis_tick_values,
)
from ....svg.circular_tracks import circular_scalar_area_line, generate_circular_depth_path_desc


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class DepthDrawer:
    """Draws a depth coverage track on a circular canvas."""

    def __init__(self, depth_config: DepthConfigurator) -> None:
        self._depth_config = depth_config
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
        self.small_tick_interval: float | None = depth_config.small_tick_interval
        self.tick_font_size: float | None = depth_config.tick_font_size
        self.axis_stroke: str = "#4b5563"
        self.axis_stroke_width: float = DEPTH_AXIS_STROKE_WIDTH_PX
        self.axis_tick_size: float = DEPTH_AXIS_TICK_SIZE_PX
        self.axis_small_tick_size: float = DEPTH_AXIS_SMALL_TICK_SIZE_PX

    def _axis_ticks(
        self, depth_df: DataFrame, baseline_radius: float, track_width: float
    ) -> list[ScalarAxisTick]:
        """The axis ticks as drawn: small ticks first, then large ticks with their labels."""

        axis_min, axis_max = _depth_axis_bounds(depth_df, self.min_depth, self.max_depth)
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

        def tick_radius(tick_value: float) -> float:
            normalized = _scaled_depth_fraction(tick_value, axis_min, axis_max, normalize=self.normalize)
            return baseline_radius + (float(track_width) * normalized)

        ticks = [ScalarAxisTick(tick_radius(value), self.axis_small_tick_size, None) for value in small_tick_values]
        ticks.extend(
            ScalarAxisTick(
                tick_radius(value),
                self.axis_tick_size,
                None
                if self.min_depth is None and math.isclose(value, axis_min, rel_tol=1e-9, abs_tol=1e-9)
                else _format_depth_tick(value),
            )
            for value in large_tick_values
        )
        return ticks

    @staticmethod
    def _radii(radius: float, track_width: float, norm_factor: float) -> tuple[float, float]:
        baseline_radius = max(
            0.0,
            (float(radius) * float(norm_factor)) - (0.5 * float(track_width)),
        )
        return baseline_radius, baseline_radius + float(track_width)

    def _add_axes(
        self,
        group: Group,
        depth_df: DataFrame,
        radius: float,
        track_width: float,
        norm_factor: float,
        axis_group_id: str = "depth_axis",
    ) -> Group:
        if not self.show_axis:
            return group
        baseline_radius, outer_radius = self._radii(radius, track_width, norm_factor)
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
        font_size = depth_axis_tick_font_size_px(self._depth_config, track_width)
        for tick in self._axis_ticks(depth_df, baseline_radius, track_width):
            axis_group.add(
                Line(
                    start=(0, -tick.radius_px),
                    end=(tick.length_px, -tick.radius_px),
                    stroke=self.axis_stroke,
                    stroke_width=self.axis_stroke_width,
                )
            )
            if tick.label is not None:
                axis_group.add(
                    Text(
                        tick.label,
                        insert=(tick.length_px + SCALAR_AXIS_LABEL_GAP_PX, -tick.radius_px),
                        fill=self.axis_stroke,
                        font_size=font_size,
                        text_anchor="start",
                        dominant_baseline="central",
                    )
                )
        group.add(axis_group)
        return group

    def ink(
        self,
        radius: float,
        depth_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        *,
        font_family: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> tuple[CircularInk, ...]:
        """What ``draw`` paints with the same arguments: the area down to the baseline with half its stroke, and
        the axis (``font_family`` is the family the axis labels inherit)."""

        plot_df = _prepare_depth_plot_dataframe(depth_df, self.min_depth, self.max_depth, normalize=self.normalize)
        line = circular_scalar_area_line(
            radius, record_len, plot_df, track_width, norm_factor,
            value_column="depth_normalized", record_transform=record_transform,
        )
        if line is None:
            return ()
        positions, outer_radii, baseline_radius = line
        pad = stroke_pad_px(self.stroke_color, self.stroke_width)
        inks: list[CircularInk] = [
            CircularInkProfile.between(positions, outer_radii, record_length=record_len, base_px=baseline_radius,
                                       pad_px=pad)
        ]
        if self.show_axis:
            inks.extend(
                scalar_axis_ink(
                    *self._radii(radius, track_width, norm_factor),
                    self._axis_ticks(depth_df, baseline_radius, track_width),
                    font_size=depth_axis_tick_font_size_px(self._depth_config, track_width),
                    font_family=font_family,
                )
            )
        return tuple(inks)

    def draw(
        self,
        radius: float,
        group: Group,
        depth_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        axis_group_id: str = "depth_axis",
        record_transform: RecordDisplayTransform | None = None,
    ) -> Group:
        plot_df = _prepare_depth_plot_dataframe(
            depth_df,
            self.min_depth,
            self.max_depth,
            normalize=self.normalize,
        )
        depth_path_desc = generate_circular_depth_path_desc(
            radius, record_len, plot_df, track_width, norm_factor,
            record_transform=record_transform,
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
        self._add_axes(group, depth_df, radius, track_width, norm_factor, axis_group_id)
        return group


__all__ = ["DepthDrawer"]
