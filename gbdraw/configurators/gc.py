#!/usr/bin/env python
# coding: utf-8

from typing import Dict

from pandas import DataFrame

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]


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
        config_dict: Dict,
        default_colors_df: DataFrame,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the GcContentConfigurator with specified settings.

        Args:
            window (int): Window size for GC content calculation.
            step (int): Step size for the moving window.
            dinucleotide (str): Specific dinucleotide sequence to focus on.
        """

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
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
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
        self.percent_background_color: str = cfg.objects.gc_content.percent_background_color
        self.percent_background_opacity: float = cfg.objects.gc_content.percent_background_opacity
        self.percent_border_color: str = cfg.objects.gc_content.percent_border_color
        self.percent_border_width: float = cfg.objects.gc_content.percent_border_width
        self.show_gc: bool = cfg.canvas.show_gc


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
        config_dict: Dict,
        default_colors_df: DataFrame,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the GcSkewConfigurator with specified settings.

        Args:
            window (int): Window size for GC content calculation.
            step (int): Step size for the moving window.
            dinucleotide (str): Specific dinucleotide sequence to focus on.
        """

        self.window: int = window
        self.step: int = step
        self.dinucleotide: str = dinucleotide
        self.high_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "skew_high"]["color"].values[0]
        )
        self.low_fill_color: str = (
            default_colors_df[default_colors_df["feature_type"] == "skew_low"]["color"].values[0]
        )
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.stroke_color: str = cfg.objects.gc_skew.stroke_color
        self.stroke_width: float = cfg.objects.gc_skew.stroke_width
        self.fill_opacity: float = cfg.objects.gc_skew.fill_opacity
        self.show_skew: bool = cfg.canvas.show_skew


__all__ = [
    "GcContentConfigurator",
    "GcSkewConfigurator",
]


