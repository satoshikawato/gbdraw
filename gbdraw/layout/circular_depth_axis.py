"""Shared circular depth-axis radial sizing helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol

from gbdraw.auto_sizes import scalar_axis_tick_font_size


DEPTH_AXIS_STROKE_WIDTH_PX = 0.8
DEPTH_AXIS_TICK_SIZE_PX = 3.0
DEPTH_AXIS_SMALL_TICK_SIZE_PX = 2.0


class DepthAxisSource(Protocol):
    """What the radial layout reads of a depth track: its axis, its ticks and their font size."""

    @property
    def show_axis(self) -> bool: ...

    @property
    def show_ticks(self) -> bool: ...

    @property
    def tick_font_size(self) -> float | None: ...


@dataclass(frozen=True)
class DepthAxisStyle:
    """A depth track's axis and ticks as a value: a render stage that replans the
    radial layout takes this instead of the depth configurator, so a depth
    window or scale change does not rebuild it."""

    show_axis: bool
    show_ticks: bool
    tick_font_size: float | None

    @classmethod
    def of(cls, source: DepthAxisSource) -> DepthAxisStyle:
        return cls(
            show_axis=bool(source.show_axis),
            show_ticks=bool(source.show_ticks),
            tick_font_size=source.tick_font_size,
        )


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
    depth_config: DepthAxisSource,
    track_width_px: float,
) -> DepthAxisFootprint:
    if not bool(depth_config.show_axis):
        return DepthAxisFootprint(0.0, 0.0)

    extra = 0.5 * DEPTH_AXIS_STROKE_WIDTH_PX
    if bool(depth_config.show_ticks):
        extra = max(extra, 0.5 * depth_axis_tick_font_size_px(depth_config, track_width_px))
    return DepthAxisFootprint(extra, extra)


__all__ = [
    "DepthAxisSource",
    "DepthAxisStyle",
    "DEPTH_AXIS_SMALL_TICK_SIZE_PX",
    "DEPTH_AXIS_STROKE_WIDTH_PX",
    "DEPTH_AXIS_TICK_SIZE_PX",
    "DepthAxisFootprint",
    "depth_axis_tick_font_size_px",
    "resolve_depth_axis_footprint",
]
