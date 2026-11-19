#!/usr/bin/env python
# coding: utf-8

"""Depth track configurator."""

from __future__ import annotations

from gbdraw.config.models import RenderProfile
from gbdraw.config.models.objects import ScalarAxisLabelStyle
from gbdraw.layout.series_style import (
    DEFAULT_SERIES_BASELINE,
    DEFAULT_SERIES_STYLE,
    DEPTH_UNIT,
    ScalarAxisStyle,
    SeriesBaseline,
    SeriesPaint,
    SeriesStyle,
)

# The stroke of a ``line`` series when the Depth stroke width is 0, which would draw nothing (D-25).
LINE_SERIES_MIN_STROKE_PX = 1.0


class DepthConfigurator:
    """Collect drawing and normalization options for depth coverage tracks.

    One series' own values (``depth_track_config``) replace the colour, the value range (``min_depth``,
    ``max_depth``), the ticks, and set ``series_style``, ``negative_color``, ``baseline`` and ``unit``; the defaults draw as
    every depth track drew before they existed.
    """

    def __init__(
        self,
        window: int,
        step: int,
        profile: RenderProfile,
        *,
        fill_color: str | None = None,
        min_depth: float | None = None,
        max_depth: float | None = None,
        normalize: bool | None = None,
        show_axis: bool | None = None,
        show_ticks: bool | None = None,
        large_tick_interval: float | None = None,
        small_tick_interval: float | None = None,
        tick_font_size: float | None = None,
        share_axis: bool | None = None,
    ) -> None:
        cfg = profile.config
        self.window = int(window)
        self.step = int(step)
        self.fill_color: str = str(fill_color or cfg.objects.depth.fill_color)
        self.stroke_color: str = cfg.objects.depth.stroke_color
        self.stroke_width: float = cfg.objects.depth.stroke_width
        self.fill_opacity: float = cfg.objects.depth.fill_opacity
        self.normalize: bool = bool(cfg.objects.depth.normalize if normalize is None else normalize)
        self.min_depth: float | None = cfg.objects.depth.min_depth if min_depth is None else float(min_depth)
        self.max_depth: float | None = cfg.objects.depth.max_depth if max_depth is None else float(max_depth)
        self.show_axis: bool = bool(cfg.objects.depth.show_axis if show_axis is None else show_axis)
        self.show_ticks: bool = bool(cfg.objects.depth.show_ticks if show_ticks is None else show_ticks)
        self.large_tick_interval: float | None = (
            cfg.objects.depth.large_tick_interval
            if large_tick_interval is None
            else float(large_tick_interval)
        )
        self.small_tick_interval: float | None = (
            cfg.objects.depth.small_tick_interval
            if small_tick_interval is None
            else float(small_tick_interval)
        )
        self.tick_font_size: float | None = (
            cfg.objects.depth.tick_font_size if tick_font_size is None else float(tick_font_size)
        )
        self.share_axis: bool = bool(cfg.objects.depth.share_axis if share_axis is None else share_axis)
        # The axis labels are measured and drawn in the figure's font (OV-2757).
        self.font_family: str = str(cfg.objects.text.font_family)
        self.axis_labels: ScalarAxisLabelStyle = cfg.objects.depth.axis_labels
        self.series_style: SeriesStyle = DEFAULT_SERIES_STYLE
        # The paint of the part below the baseline; None: the series colour.
        self.negative_color: str | None = None
        self.baseline: SeriesBaseline = DEFAULT_SERIES_BASELINE
        # The text after each axis number and in the range label; empty for none.
        self.unit: str = DEPTH_UNIT

    @property
    def series_paint(self) -> SeriesPaint:
        """How the track's path is painted: an ``area`` filled in the series colour, a ``line`` stroked in it at the
        Depth stroke width (1 px when that is 0); the part below the baseline in the negative colour when it is
        another colour."""

        below = None if self.negative_color in (None, self.fill_color) else self.negative_color
        if self.series_style == "line":
            width = self.stroke_width if self.stroke_width > 0 else LINE_SERIES_MIN_STROKE_PX
            return SeriesPaint("none", self.fill_color, width, self.fill_opacity, None, below)
        return SeriesPaint(self.fill_color, self.stroke_color, self.stroke_width, self.fill_opacity, below, None)

    @property
    def legend_paint(self) -> SeriesPaint:
        """The paint of the track's one Legend row: its colour as an area, whatever its style and negative colour
        (one row per series)."""

        return SeriesPaint(self.fill_color, self.stroke_color, self.stroke_width, self.fill_opacity, None, None)

    @property
    def axis_style(self) -> ScalarAxisStyle:
        """The axis the track draws and keeps room for."""

        return ScalarAxisStyle.of(self)

__all__ = ["LINE_SERIES_MIN_STROKE_PX", "DepthConfigurator"]
