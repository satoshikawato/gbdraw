"""Shared circular depth-axis radial sizing helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from gbdraw.auto_sizes import scalar_axis_tick_font_size


DEPTH_AXIS_STROKE_WIDTH_PX = 0.8
DEPTH_AXIS_TICK_SIZE_PX = 3.0
DEPTH_AXIS_SMALL_TICK_SIZE_PX = 2.0


@dataclass(frozen=True)
class DepthAxisFootprint:
    radial_inner_extra_px: float
    radial_outer_extra_px: float


def depth_axis_tick_font_size_px(
    depth_config: Any,
    track_width_px: float,
) -> float:
    tick_font_size = getattr(depth_config, "tick_font_size", None)
    if tick_font_size is not None:
        return max(0.0, float(tick_font_size))
    return scalar_axis_tick_font_size("circular", track_width_px)


def resolve_depth_axis_footprint(
    depth_config: Any,
    track_width_px: float,
) -> DepthAxisFootprint:
    if not bool(depth_config.show_axis):
        return DepthAxisFootprint(0.0, 0.0)

    extra = 0.5 * DEPTH_AXIS_STROKE_WIDTH_PX
    if bool(depth_config.show_ticks):
        extra = max(extra, 0.5 * depth_axis_tick_font_size_px(depth_config, track_width_px))
    return DepthAxisFootprint(extra, extra)


__all__ = [
    "DEPTH_AXIS_SMALL_TICK_SIZE_PX",
    "DEPTH_AXIS_STROKE_WIDTH_PX",
    "DEPTH_AXIS_TICK_SIZE_PX",
    "DepthAxisFootprint",
    "depth_axis_tick_font_size_px",
    "resolve_depth_axis_footprint",
]
