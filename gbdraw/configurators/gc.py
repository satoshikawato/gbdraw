#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import copy

from pandas import DataFrame

from gbdraw.config.models import RenderProfile
from gbdraw.config.models.objects import ScalarAxisLabelStyle
from gbdraw.layout.series_style import ScalarAxisStyle, SeriesPaint
from gbdraw.tracks.parsing import SlotWithParams, slot_skew_fill_colors


def _slot_skew_config(skew_config: GcSkewConfigurator, slot: SlotWithParams, dinucleotide: str) -> GcSkewConfigurator:
    """The skew configuration slot ``slot`` draws with: its ``dinucleotide`` and its own fills
    (``slot_skew_fill_colors``)."""

    cloned = copy.copy(skew_config)
    cloned.dinucleotide = str(dinucleotide).upper()
    cloned.high_fill_color, cloned.low_fill_color = slot_skew_fill_colors(
        slot, cloned.high_fill_color, cloned.low_fill_color
    )
    return cloned


class GcContentConfigurator:
    """
    Configurator for GC content calculation parameters.

    This class holds the configuration for computing GC content, including the window size, step size,
    and the specific dinucleotide to consider.

    Attributes:
        window (int): Size of the window used for calculating GC content.
        step (int): Step size for moving the window across the genomic sequence.
        dinucleotide (str): Specific dinucleotide sequence to consider for GC content calculation.
    """

    def __init__(
        self,
        window: int,
        step: int,
        dinucleotide: str,
        profile: RenderProfile,
        default_colors_df: DataFrame,
    ) -> None:
        """
        Initializes the GcContentConfigurator with specified settings.

        Args:
            window (int): Window size for GC content calculation.
            step (int): Step size for the moving window.
            dinucleotide (str): Specific dinucleotide sequence to focus on.
        """
        cfg = profile.config
        self.window: int = window
        self.step: int = step
        self.dinucleotide: str = dinucleotide
        self.high_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "gc_content"]["color"].values[0]
        )
        self.low_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "gc_content"]["color"].values[0]
        )
        self.fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "gc_content"]["color"].values[0]
        )
        self.stroke_color: str = cfg.objects.gc_content.stroke_color
        self.stroke_width: float = cfg.objects.gc_content.stroke_width
        self.fill_opacity: float = cfg.objects.gc_content.fill_opacity
        self.mode: str = cfg.objects.gc_content.mode
        self.min_percent: float | None = cfg.objects.gc_content.min_percent
        self.max_percent: float | None = cfg.objects.gc_content.max_percent
        self.show_axis: bool = cfg.objects.gc_content.show_axis
        self.show_ticks: bool = cfg.objects.gc_content.show_ticks
        self.large_tick_interval: float | None = cfg.objects.gc_content.large_tick_interval
        self.tick_interval: float | None = self.large_tick_interval
        self.small_tick_interval: float | None = cfg.objects.gc_content.small_tick_interval
        self.tick_font_size: float | None = cfg.objects.gc_content.tick_font_size
        # The percent axis labels are measured and drawn in the figure's font (OV-2757).
        self.font_family: str = str(cfg.objects.text.font_family)
        self.axis_labels: ScalarAxisLabelStyle = cfg.objects.gc_content.axis_labels
        self.percent_background_color: str = cfg.objects.gc_content.percent_background_color
        self.percent_background_opacity: float = cfg.objects.gc_content.percent_background_opacity
        self.percent_border_color: str = cfg.objects.gc_content.percent_border_color
        self.percent_border_width: float = cfg.objects.gc_content.percent_border_width

    @property
    def series_paint(self) -> SeriesPaint:
        """How the track's path is painted, as drawn and as its Legend row shows it."""

        return SeriesPaint(
            self.fill_color, self.stroke_color, self.stroke_width, self.fill_opacity, None,
            frozenset({"circular", "linear"}),
        )

    @property
    def axis_style(self) -> ScalarAxisStyle | None:
        """The axis the track draws and keeps room for: in percent mode only."""

        return ScalarAxisStyle.of(self) if self.mode == "percent" else None

class GcSkewConfigurator:
    """
    Configurator for GC skew calculation parameters.

    This class holds the configuration for computing GC skew, including the window size, step size,
    and the specific dinucleotide to consider.

    Attributes:
        window (int): Size of the window used for calculating GC content.
        step (int): Step size for moving the window across the genomic sequence.
        dinucleotide (str): Specific dinucleotide sequence to consider for GC content calculation.
    """

    def __init__(
        self,
        window: int,
        step: int,
        dinucleotide: str,
        profile: RenderProfile,
        default_colors_df: DataFrame,
    ) -> None:
        """
        Initializes the GcSkewConfigurator with specified settings.

        Args:
            window (int): Window size for GC content calculation.
            step (int): Step size for the moving window.
            dinucleotide (str): Specific dinucleotide sequence to focus on.
        """
        cfg = profile.config
        self.window: int = window
        self.step: int = step
        self.dinucleotide: str = dinucleotide
        self.high_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "skew_high"]["color"].values[0]
        )
        self.low_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "skew_low"]["color"].values[0]
        )
        self.stroke_color: str = cfg.objects.gc_skew.stroke_color
        self.stroke_width: float = cfg.objects.gc_skew.stroke_width
        self.fill_opacity: float = cfg.objects.gc_skew.fill_opacity

    @property
    def series_paint(self) -> SeriesPaint:
        """How the track's path is painted, as drawn and as its Legend rows show it: the high fill above the
        baseline and the low fill below it."""

        return SeriesPaint(
            self.high_fill_color, self.stroke_color, self.stroke_width, self.fill_opacity, self.low_fill_color,
            frozenset({"circular"}),
        )

__all__ = [
    "GcContentConfigurator",
    "GcSkewConfigurator",
]
