"""Shared circular depth-axis radial sizing helpers."""

from __future__ import annotations

import math

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, NamedTuple, Protocol

from gbdraw.auto_sizes import scalar_axis_tick_font_size
from gbdraw.config.models.objects import ScalarAxisLabelStyle
from gbdraw.core.text import calculate_svg_bbox_dimensions

from .circular import RadialBand
from .circular_ink import CircularInk, CircularInkPolygon, CircularInkRing
from .spatial import Aabb


DEPTH_AXIS_STROKE_WIDTH_PX = 0.8
DEPTH_AXIS_TICK_SIZE_PX = 3.0
DEPTH_AXIS_SMALL_TICK_SIZE_PX = 2.0
# How far a scalar axis tick label starts past the end of its tick.
SCALAR_AXIS_LABEL_GAP_PX = 1.5


class ScalarAxisTick(NamedTuple):
    """One tick of a Circular scalar axis (depth, GC content in percent) as drawn: at ``radius_px`` on the axis
    at the top of the circle, ``length_px`` long to the right, with its ``label`` past it or none."""

    radius_px: float
    length_px: float
    label: str | None


class DepthAxisSource(Protocol):
    """What the radial layout reads of a depth track: its axis, its ticks and their font size."""

    @property
    def show_axis(self) -> bool: ...

    @property
    def show_ticks(self) -> bool: ...

    @property
    def tick_font_size(self) -> float | None: ...

    @property
    def axis_labels(self) -> ScalarAxisLabelStyle: ...


@dataclass(frozen=True)
class DepthAxisStyle:
    """A depth track's axis and ticks as a value: a render stage that replans the
    radial layout takes this instead of the depth configurator, so a depth
    window or scale change does not rebuild it."""

    show_axis: bool
    show_ticks: bool
    tick_font_size: float | None
    axis_labels: ScalarAxisLabelStyle

    @classmethod
    def of(cls, source: DepthAxisSource) -> DepthAxisStyle:
        return cls(
            show_axis=bool(source.show_axis),
            show_ticks=bool(source.show_ticks),
            tick_font_size=source.tick_font_size,
            axis_labels=source.axis_labels,
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
    # Numbered tick labels are centred on their ticks, so the outermost and innermost ones reach half a line past
    # the band; the range label sits inside the band.
    if bool(depth_config.show_ticks) and depth_config.axis_labels == "numbered":
        extra = max(extra, 0.5 * depth_axis_tick_font_size_px(depth_config, track_width_px))
    return DepthAxisFootprint(extra, extra)


# The range label of a Circular scalar axis runs clockwise along the circle inside the band: its alphabetic
# baseline this many label lines inside the band's outer edge (the label's top then sits about a fifth of a line
# inside it), starting past the tick marks right of the axis. The path carries the baseline, not the label's
# centre, because CairoSVG ignores ``dominant-baseline`` on text along a path while browsers honour it.
RANGE_LABEL_BASELINE_INSET_LINES = 0.95
RANGE_LABEL_TOP_INSET_LINES = 0.2
RANGE_LABEL_BOTTOM_INSET_LINES = 1.2


class CircularRangeLabel(NamedTuple):
    """Where a Circular scalar axis draws its range label: ``text`` in ``font_size`` along the circle of
    ``baseline_px`` from ``start_rad`` clockwise of the axis at 12 o'clock, on a path that ends at 3 o'clock."""

    text: str
    font_size: float
    baseline_px: float
    outer_px: float
    start_rad: float

    @property
    def start_deg(self) -> float:
        return math.degrees(self.start_rad)

    @property
    def path_span_deg(self) -> float:
        return 90.0 - self.start_deg

    def bounds(self, font_family: str) -> Aabb:
        """The box the label paints in ``font_family``, cut at 3 o'clock where its path ends."""

        width, _ = calculate_svg_bbox_dimensions(self.text, font_family, self.font_size)
        end = min(self.start_rad + float(width) / self.baseline_px, 0.5 * math.pi)
        top = self.outer_px - RANGE_LABEL_TOP_INSET_LINES * self.font_size
        bottom = max(0.0, self.outer_px - RANGE_LABEL_BOTTOM_INSET_LINES * self.font_size)
        return Aabb(bottom * math.sin(self.start_rad), -top, top * math.sin(end), -bottom * math.cos(end))


def circular_range_label(text: str | None, outer_px: float, *, font_size: float) -> CircularRangeLabel | None:
    """The range label ``text`` of an axis whose band ends at ``outer_px``, starting past the large tick marks and
    their gap; none without a text or when the circle is too small to carry it."""

    if text is None:
        return None
    baseline = float(outer_px) - RANGE_LABEL_BASELINE_INSET_LINES * float(font_size)
    start_px = DEPTH_AXIS_TICK_SIZE_PX + SCALAR_AXIS_LABEL_GAP_PX
    if baseline <= start_px:
        return None
    return CircularRangeLabel(text, float(font_size), baseline, float(outer_px), start_px / baseline)


def scalar_axis_ink(
    baseline_px: float,
    outer_px: float,
    ticks: Sequence[ScalarAxisTick],
    *,
    font_size: float,
    font_family: str,
    range_label: CircularRangeLabel | None = None,
) -> tuple[CircularInk, ...]:
    """What a Circular scalar axis paints: the baseline circle, the axis line up at the top of the circle from
    the baseline to ``outer_px``, its ticks, and their labels (``text-anchor`` start, ``dominant-baseline``
    central, in the family they inherit)."""

    half = 0.5 * DEPTH_AXIS_STROKE_WIDTH_PX
    inks: list[CircularInk] = [
        CircularInkRing(RadialBand(baseline_px - half, baseline_px + half)),
        CircularInkPolygon.segment((0.0, -baseline_px), (0.0, -outer_px), half),
    ]
    for tick in ticks:
        inks.append(CircularInkPolygon.segment((0.0, -tick.radius_px), (tick.length_px, -tick.radius_px), half))
        if tick.label is not None:
            width, _ = calculate_svg_bbox_dimensions(tick.label, font_family, font_size)
            left = tick.length_px + SCALAR_AXIS_LABEL_GAP_PX
            inks.append(
                CircularInkPolygon.box(
                    Aabb(left, -tick.radius_px - 0.5 * font_size, left + width, -tick.radius_px + 0.5 * font_size)
                )
            )
    if range_label is not None:
        inks.append(CircularInkPolygon.box(range_label.bounds(font_family)))
    return tuple(inks)


__all__ = [
    "CircularRangeLabel",
    "circular_range_label",
    "DepthAxisSource",
    "DepthAxisStyle",
    "DEPTH_AXIS_SMALL_TICK_SIZE_PX",
    "DEPTH_AXIS_STROKE_WIDTH_PX",
    "DEPTH_AXIS_TICK_SIZE_PX",
    "SCALAR_AXIS_LABEL_GAP_PX",
    "ScalarAxisTick",
    "DepthAxisFootprint",
    "depth_axis_tick_font_size_px",
    "resolve_depth_axis_footprint",
    "scalar_axis_ink",
]
