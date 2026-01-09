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


