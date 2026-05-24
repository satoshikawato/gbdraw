#!/usr/bin/env python
# coding: utf-8

from typing import Dict

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.core.color import (
    COLLINEAR_ORIENTATION_COLOR_KEYS,
    COLLINEAR_ORIENTATION_MIN_COLOR_KEYS,
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
)


def _default_color(default_colors_df: DataFrame, feature_type: str, fallback: str) -> str:
    matching_rows = default_colors_df[default_colors_df["feature_type"] == feature_type]
    if matching_rows.empty:
        return fallback
    return str(matching_rows["color"].values[0])


class BlastMatchConfigurator:
    """
    Configurator for BLAST match visualization parameters.

    Manages the settings for visualizing BLAST match results, such as thresholds for e-value, bitscore,
    identity, and alignment length, along with a dictionary of sequence lengths.

    Attributes:
        evalue (float): E-value threshold for filtering BLAST matches.
        bitscore (float): Bitscore threshold for filtering BLAST matches.
        identity (float): Identity percentage threshold for filtering BLAST matches.
        alignment_length (int): Minimum alignment length threshold for filtering BLAST matches.
        sequence_length_dict (Dict[str, int]): Dictionary containing the length of each sequence.
    """

    def __init__(
        self,
        evalue: float,
        bitscore: float,
        identity: float,
        alignment_length: int,
        sequence_length_dict: Dict[str, int],
        config_dict: Dict,
        default_colors_df: DataFrame,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the BlastMatchConfigurator with specific thresholds and sequence length data.

        Args:
            evalue (float): E-value threshold for BLAST matches.
            bitscore (float): Bitscore threshold for BLAST matches.
            identity (float): Identity percentage threshold for BLAST matches.
            alignment_length (int): Minimum alignment length threshold for BLAST matches.
            sequence_length_dict (Dict[str, int]): Lengths of the sequences involved in BLAST matches.
        """

        self.evalue: float = evalue
        self.bitscore: float = bitscore
        self.identity: float = identity
        self.alignment_length: int = alignment_length
        self.sequence_length_dict: dict = sequence_length_dict
        self.min_color: str = default_colors_df[default_colors_df["feature_type"] == "pairwise_match_min"][
            "color"
        ].values[0]
        self.max_color: str = default_colors_df[default_colors_df["feature_type"] == "pairwise_match_max"][
            "color"
        ].values[0]
        self.fill_color: str = default_colors_df[default_colors_df["feature_type"] == "pairwise_match"][
            "color"
        ].values[0]
        self.collinearity_orientation_colors: dict[str, str] = {
            orientation: _default_color(
                default_colors_df,
                color_key,
                DEFAULT_COLLINEAR_ORIENTATION_COLORS[orientation],
            )
            for orientation, color_key in COLLINEAR_ORIENTATION_COLOR_KEYS.items()
        }
        self.collinearity_orientation_min_colors: dict[str, str] = {
            orientation: _default_color(
                default_colors_df,
                color_key,
                DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS[orientation],
            )
            for orientation, color_key in COLLINEAR_ORIENTATION_MIN_COLOR_KEYS.items()
        }
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.fill_opacity: float = cfg.objects.blast_match.fill_opacity
        self.stroke_color: str = cfg.objects.blast_match.stroke_color
        self.stroke_width: float = cfg.objects.blast_match.stroke_width
        self.match_style: str = str(cfg.objects.blast_match.style)
        self.curve_tension: float = cfg.objects.blast_match.curve_tension


__all__ = ["BlastMatchConfigurator"]


