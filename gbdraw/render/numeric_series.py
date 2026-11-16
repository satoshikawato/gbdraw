"""One numeric series track as a value that both modes draw: GC content, GC skew and depth.

``NumericSeriesPlot`` says what one series track draws: its samples and baseline rule (``SeriesSamples``), its paint,
the band frame of a GC content track in percent, and its axis. The three builders below turn a track's configurator
and data frame into that value; they are the only code that tells the kinds apart. The Circular and Linear drawers
(``render/drawers/{circular,linear}/numeric_series.py``) read the value and nothing else, so what a track paints, what
the Circular label placement reserves for it (its ink) and what the Linear left room keeps for its axis come from one
value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from pandas import DataFrame

from ..analysis.gc import gc_content_percent_df
from ..configurators import DepthConfigurator, GcContentConfigurator, GcSkewConfigurator
from ..layout.circular_ink import paint_shows
from ..layout.scalar_axis import (
    DEPTH_UNIT,
    PERCENT_UNIT,
    ScalarAxis,
    SeriesSamples,
    _depth_axis_bounds,
    _prepare_depth_plot_dataframe,
    band_samples,
    deviation_samples,
)

DiagramMode = Literal["circular", "linear"]
# The render stage that builds each kind's track group; a stage keeps its own name, so a cache entry of one kind
# never serves another.
NumericSeriesStage = Literal["gc_content_content", "gc_skew_content", "depth_content"]
NumericSeriesConfig = GcContentConfigurator | GcSkewConfigurator | DepthConfigurator
# The colour of a series axis, its ticks and its labels, in both modes.
AXIS_STROKE = "#4b5563"


@dataclass(frozen=True)
class SeriesPaint:
    """How a series path is painted. ``below_fill`` is the fill of the part below the baseline (GC skew), drawn
    by a copy of the path clipped to the band; none for a series of one fill."""

    fill: str
    stroke: str
    stroke_width: float
    fill_opacity: float
    below_fill: str | None
    # The modes in which a stroke that paints nothing (width 0 and colour ``none``) writes no ``stroke-width``;
    # the other modes always write it. GC content leaves it out in both modes, GC skew in Circular only, depth in
    # neither: the old drawers differed, and the output keeps their bytes.
    bare_stroke_modes: frozenset[DiagramMode]

    def writes_stroke_width(self, mode: DiagramMode) -> bool:
        return mode not in self.bare_stroke_modes or self.stroke_width > 0 or self.stroke.strip().lower() != "none"


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
    # The top left corner of a Linear track in its group: GC tracks write it as 0 and depth as 0.0, which svgwrite
    # writes as given, so the output keeps both.
    linear_origin: float


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


def gc_content_axis(config: GcContentConfigurator) -> ScalarAxis | None:
    """The axis of a GC content track: in percent mode, when it is shown; drawn once on a shared circle."""

    if config.mode != "percent" or not bool(config.show_axis):
        return None
    axis_min, axis_max = _percent_bounds(config)
    return ScalarAxis(
        axis_min=axis_min,
        axis_max=axis_max,
        unit=PERCENT_UNIT,
        log_scale=False,
        label_axis_min=True,
        show_ticks=bool(config.show_ticks),
        large_tick_interval=config.large_tick_interval,
        small_tick_interval=config.small_tick_interval,
        tick_font_size=config.tick_font_size,
        axis_labels=config.axis_labels,
        font_family=config.font_family,
        once_on_shared_circle=True,
    )


def depth_axis(config: DepthConfigurator, depth_df: DataFrame) -> ScalarAxis | None:
    """The axis of a depth track, when it is shown: from the set minimum (else 0) to the set maximum (else the
    largest depth), linear or log-scaled; each record of a shared circle draws its own."""

    if not bool(config.show_axis):
        return None
    axis_min, axis_max = _depth_axis_bounds(depth_df, config.min_depth, config.max_depth)
    return ScalarAxis(
        axis_min=axis_min,
        axis_max=axis_max,
        unit=DEPTH_UNIT,
        log_scale=bool(config.normalize),
        label_axis_min=config.min_depth is not None,
        show_ticks=bool(config.show_ticks),
        large_tick_interval=config.large_tick_interval,
        small_tick_interval=config.small_tick_interval,
        tick_font_size=config.tick_font_size,
        axis_labels=config.axis_labels,
        font_family=config.font_family,
        once_on_shared_circle=False,
    )


def gc_content_series(config: GcContentConfigurator, gc_df: DataFrame) -> NumericSeriesPlot:
    """A GC content track: the deviation of each window's content from the mean, or, in percent mode, the content
    on its axis with the band frame."""

    paint = SeriesPaint(
        config.fill_color, config.stroke_color, config.stroke_width, config.fill_opacity, None,
        frozenset({"circular", "linear"}),
    )
    if config.mode != "percent":
        samples = None if gc_df.empty else deviation_samples(
            gc_df[f"{config.dinucleotide} content"], closes_at_record_end=True
        )
        return NumericSeriesPlot(samples, paint, None, None, "", 0)
    axis_min, axis_max = _percent_bounds(config)
    plot_df = gc_content_percent_df(gc_df, dinucleotide=config.dinucleotide, min_percent=axis_min, max_percent=axis_max)
    return NumericSeriesPlot(
        band_samples(plot_df, "value_normalized", closes_at_record_end=True),
        paint,
        gc_content_band_frame(config),
        gc_content_axis(config),
        "",
        0,
    )


def _skew_column(gc_df: DataFrame, dinucleotide: str) -> str:
    """The skew column of a GC frame: the dinucleotide's, else the frame's one skew column."""

    column = f"{dinucleotide} skew"
    if column in gc_df.columns:
        return column
    return [name for name in gc_df.columns if "skew" in name and "cumulative" not in name.lower()][0]


def gc_skew_series(config: GcSkewConfigurator, gc_df: DataFrame) -> NumericSeriesPlot:
    """A GC skew track: the deviation of each window's skew from the mean, in the high colour above the baseline
    and the low colour below it."""

    samples = None if gc_df.empty else deviation_samples(
        gc_df[_skew_column(gc_df, config.dinucleotide)], closes_at_record_end=False
    )
    paint = SeriesPaint(
        config.high_fill_color, config.stroke_color, config.stroke_width, config.fill_opacity, config.low_fill_color,
        frozenset({"circular"}),
    )
    return NumericSeriesPlot(samples, paint, None, None, config.dinucleotide, 0)


def depth_series(config: DepthConfigurator, depth_df: DataFrame) -> NumericSeriesPlot:
    """A depth track (or a feature density histogram): each window's depth as the share of its axis it fills."""

    plot_df = _prepare_depth_plot_dataframe(depth_df, config.min_depth, config.max_depth, normalize=config.normalize)
    paint = SeriesPaint(
        config.fill_color, config.stroke_color, config.stroke_width, config.fill_opacity, None, frozenset()
    )
    return NumericSeriesPlot(
        band_samples(plot_df, "depth_normalized", closes_at_record_end=False),
        paint,
        None,
        depth_axis(config, depth_df),
        "",
        0.0,
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
    "SeriesPaint",
    "depth_axis",
    "depth_series",
    "gc_content_axis",
    "gc_content_band_frame",
    "gc_content_series",
    "gc_skew_series",
    "numeric_series_plot",
]
