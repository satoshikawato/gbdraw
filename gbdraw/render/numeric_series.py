"""One numeric series track as a value that both modes draw: GC content, GC skew and depth.

``NumericSeriesPlot`` says what one series track draws: its samples and baseline rule (``SeriesSamples``), its paint
(``SeriesPaint``, from the configurator's ``series_paint``), the band frame of a GC content track in percent, and its
axis (``ScalarAxis``, whose style is the configurator's ``axis_style``). The three builders below turn a track's
configurator and data frame into that value, and ``numeric_series_plot`` picks the builder of a configurator's kind.
The Circular and Linear drawers (``render/drawers/{circular,linear}/numeric_series.py``) read the value and nothing
else, so what a track paints, what the Circular label placement reserves for it (its ink), the room both modes keep
for its axis and its Legend rows come from one owner.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from pandas import DataFrame

from ..analysis.gc import gc_content_percent_df
from ..configurators import DepthConfigurator, GcContentConfigurator, GcSkewConfigurator
from ..layout.circular_ink import paint_shows
from ..layout.scalar_axis import (
    PERCENT_UNIT,
    ScalarAxis,
    SeriesSamples,
    _depth_axis_bounds,
    _prepare_depth_plot_dataframe,
    band_samples,
    deviation_samples,
    series_baseline_fraction,
)
from ..layout.series_style import ScalarAxisStyle, SeriesPaint

# The render stage that builds each kind's track group; a stage keeps its own name, so a cache entry of one kind
# never serves another.
NumericSeriesStage = Literal["gc_content_content", "gc_skew_content", "depth_content"]
NumericSeriesConfig = GcContentConfigurator | GcSkewConfigurator | DepthConfigurator
# The colour of a series axis, its ticks and its labels, in both modes.
AXIS_STROKE = "#4b5563"


@dataclass(frozen=True)
class SeriesBandFrame:
    """The band of a GC content track in percent: its background (a ring or a rectangle) and its border (two
    circles or a rectangle). Each is drawn whenever the track draws; only those that show count as paint."""

    background_color: str
    background_opacity: float
    border_color: str
    border_width: float

    @property
    def paints_background(self) -> bool:
        return paint_shows(self.background_color, self.background_opacity)

    @property
    def paints_border(self) -> bool:
        return self.border_width > 0 and paint_shows(self.border_color)

    @property
    def border_overhang_px(self) -> float:
        """How far the border paints past each edge of the band: half its width when it paints, else nothing. The
        Circular ink and the Linear slot footprint read it."""

        return 0.5 * self.border_width if self.paints_border else 0.0

    @property
    def paints_edges(self) -> bool:
        """Whether the track paints its band's edges: its background or its border shows. The tick labels next to
        the band keep their margin from it (OV-3116)."""

        return self.paints_background or self.paints_border


@dataclass(frozen=True, eq=False)
class NumericSeriesPlot:
    """What one numeric series track draws, in either mode."""

    # None: the track draws nothing (no windows, or no window with a value).
    samples: SeriesSamples | None
    paint: SeriesPaint
    frame: SeriesBandFrame | None
    # None: no axis (deviation mode, or the axis hidden).
    axis: ScalarAxis | None
    # The dinucleotide the clip id of the part below the baseline names (GC skew); empty without one.
    clip_key: str


def gc_content_band_frame(config: GcContentConfigurator) -> SeriesBandFrame | None:
    """The band frame of a GC content track: in percent mode only."""

    if config.mode != "percent":
        return None
    return SeriesBandFrame(
        str(config.percent_background_color),
        float(config.percent_background_opacity),
        str(config.percent_border_color),
        float(config.percent_border_width),
    )


def _percent_bounds(config: GcContentConfigurator) -> tuple[float, float]:
    return (
        float(0.0 if config.min_percent is None else config.min_percent),
        float(100.0 if config.max_percent is None else config.max_percent),
    )


def _drawn_style(style: ScalarAxisStyle | None) -> ScalarAxisStyle | None:
    return style if style is not None and style.show_axis else None


def gc_content_axis(config: GcContentConfigurator) -> ScalarAxis | None:
    """The axis of a GC content track: in percent mode, when it is shown (``axis_style``); drawn once on a shared
    circle."""

    style = _drawn_style(config.axis_style)
    if style is None:
        return None
    axis_min, axis_max = _percent_bounds(config)
    return ScalarAxis(
        axis_min=axis_min,
        axis_max=axis_max,
        unit=PERCENT_UNIT,
        log_scale=False,
        label_axis_min=True,
        style=style,
        large_tick_interval=config.large_tick_interval,
        small_tick_interval=config.small_tick_interval,
        font_family=config.font_family,
        once_on_shared_circle=True,
    )


def depth_axis(config: DepthConfigurator, depth_df: DataFrame) -> ScalarAxis | None:
    """The axis of a depth track, when it is shown: from the set minimum to the set maximum, else from the
    smallest value (0 for values of 0 or more) to the largest, linear or log-scaled, each number followed by the
    series' ``unit``; each record of a shared circle draws its own. Its lower bound is labelled when it is set or
    below 0."""

    style = _drawn_style(config.axis_style)
    if style is None:
        return None
    axis_min, axis_max = _depth_axis_bounds(depth_df, config.min_depth, config.max_depth)
    return ScalarAxis(
        axis_min=axis_min,
        axis_max=axis_max,
        unit=config.unit,
        log_scale=bool(config.normalize),
        label_axis_min=config.min_depth is not None or axis_min < 0,
        style=style,
        large_tick_interval=config.large_tick_interval,
        small_tick_interval=config.small_tick_interval,
        font_family=config.font_family,
        once_on_shared_circle=False,
    )


def gc_content_series(config: GcContentConfigurator, gc_df: DataFrame) -> NumericSeriesPlot:
    """A GC content track: the deviation of each window's content from the mean, or, in percent mode, the content
    on its axis with the band frame."""

    paint = config.series_paint
    if config.mode != "percent":
        samples = None if gc_df.empty else deviation_samples(
            gc_df[f"{config.dinucleotide} content"], closes_at_record_end=True
        )
        return NumericSeriesPlot(samples, paint, None, None, "")
    axis_min, axis_max = _percent_bounds(config)
    plot_df = gc_content_percent_df(gc_df, dinucleotide=config.dinucleotide, min_percent=axis_min, max_percent=axis_max)
    return NumericSeriesPlot(
        band_samples(plot_df, "value_normalized", closes_at_record_end=True),
        paint,
        gc_content_band_frame(config),
        gc_content_axis(config),
        "",
    )


def gc_skew_series(config: GcSkewConfigurator, gc_df: DataFrame) -> NumericSeriesPlot:
    """A GC skew track: the deviation of each window's skew from the mean, in the high colour above the baseline
    and the low colour below it."""

    samples = None if gc_df.empty else deviation_samples(
        gc_df[f"{config.dinucleotide} skew"], closes_at_record_end=False
    )
    return NumericSeriesPlot(samples, config.series_paint, None, None, config.dinucleotide)


def depth_series(config: DepthConfigurator, depth_df: DataFrame) -> NumericSeriesPlot:
    """A depth track (or a feature density histogram): each window's value as the share of its axis it fills, drawn
    from the series' baseline (``baseline``) as an area or a line (``series_style``)."""

    plot_df = _prepare_depth_plot_dataframe(depth_df, config.min_depth, config.max_depth, normalize=config.normalize)
    axis_min, axis_max = _depth_axis_bounds(depth_df, config.min_depth, config.max_depth)
    mean = float(depth_df["depth"].mean()) if "depth" in depth_df.columns and not depth_df.empty else math.nan
    baseline = series_baseline_fraction(config.baseline, axis_min, axis_max, mean, log_scale=bool(config.normalize))
    return NumericSeriesPlot(
        band_samples(
            plot_df, "depth_normalized", closes_at_record_end=False, baseline=baseline, style=config.series_style
        ),
        config.series_paint,
        None,
        depth_axis(config, depth_df),
        "",
    )


def numeric_series_plot(config: NumericSeriesConfig, frame: DataFrame) -> NumericSeriesPlot:
    """The series a track's configurator draws from ``frame``: the one place that maps a kind to its builder."""

    if isinstance(config, GcContentConfigurator):
        return gc_content_series(config, frame)
    if isinstance(config, GcSkewConfigurator):
        return gc_skew_series(config, frame)
    return depth_series(config, frame)


__all__ = [
    "AXIS_STROKE",
    "NumericSeriesConfig",
    "NumericSeriesPlot",
    "NumericSeriesStage",
    "SeriesBandFrame",
    "depth_axis",
    "depth_series",
    "gc_content_axis",
    "gc_content_band_frame",
    "gc_content_series",
    "gc_skew_series",
    "numeric_series_plot",
]
