#!/usr/bin/env python
# coding: utf-8

from typing import Dict

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]


class BlastMatchConfigurator:
    """
    Configurator for BLAST match visualization parameters.

    Manages the settings for visualizing BLAST match results, such as thresholds for e-value, bitscore,
    and identity, along with a dictionary of sequence lengths.

    Attributes:
        evalue (float): E-value threshold for filtering BLAST matches.
        bitscore (float): Bitscore threshold for filtering BLAST matches.
        identity (float): Identity percentage threshold for filtering BLAST matches.
        sequence_length_dict (Dict[str, int]): Dictionary containing the length of each sequence.
    """

    def __init__(
        self,
        evalue: float,
        bitscore: float,
        identity: float,
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
            sequence_length_dict (Dict[str, int]): Lengths of the sequences involved in BLAST matches.
        """

        self.evalue: float = evalue
        self.bitscore: float = bitscore
        self.identity: float = identity
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
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.fill_opacity: float = cfg.objects.blast_match.fill_opacity
        self.stroke_color: str = cfg.objects.blast_match.stroke_color
        self.stroke_width: float = cfg.objects.blast_match.stroke_width


__all__ = ["BlastMatchConfigurator"]


