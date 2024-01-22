#!/usr/bin/env python
# coding: utf-8
from typing import Dict, Optional
from pandas import DataFrame


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

    def __init__(self, window: int, step: int, dinucleotide: str) -> None:
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


class FeatureDrawingConfigurator:
    """
    Configurator for drawing genomic features.

    This class provides configurations for how genomic features should be visualized, including 
    color mapping and selection of features to display.

    Attributes:
        color_table (Optional[DataFrame]): Table defining color codes for different features.
        default_colors (Optional[DataFrame]): Default color settings for features.
        selected_features_set (str): Identifier for a set of features to be visualized.
    """

    def __init__(self, color_table: Optional[DataFrame], default_colors: Optional[DataFrame], selected_features_set: str) -> None:
        """
        Initializes the FeatureDrawingConfigurator with color settings and feature selection.

        Args:
            color_table (Optional[DataFrame]): Color table for features.
            default_colors (Optional[DataFrame]): Default colors for features.
            selected_features_set (str): Set identifier for selecting features to display.
        """
        self.color_table: Optional[DataFrame] = color_table
        self.default_colors: Optional[DataFrame] = default_colors
        self.selected_features_set: str = selected_features_set


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

    def __init__(self, evalue: float, bitscore: float, identity: float, sequence_length_dict: Dict[str, int]) -> None:
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
