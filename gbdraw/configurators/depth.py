#!/usr/bin/env python
# coding: utf-8

"""Depth track configurator."""

from __future__ import annotations

from gbdraw.config.models import RenderProfile  # type: ignore[reportMissingImports]


class DepthConfigurator:
    """Collect drawing and normalization options for depth coverage tracks."""

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

__all__ = ["DepthConfigurator"]
