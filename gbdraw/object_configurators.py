#!/usr/bin/env python
# coding: utf-8
from typing import Dict, Optional, List
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

    def __init__(self, window: int, step: int, dinucleotide: str, config_dict: Dict, default_colors_df: DataFrame, show_gc: bool) -> None:
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
        self.high_fill_color: str = default_colors_df[default_colors_df['feature_type'] == 'gc_content']['color'].values[0]
        self.low_fill_color: str = default_colors_df[default_colors_df['feature_type'] == 'gc_content']['color'].values[0]
        self.fill_color: str = default_colors_df[default_colors_df['feature_type'] == 'gc_content']['color'].values[0]
        self.stroke_color: str = config_dict['objects']['gc_content']['stroke_color']
        self.stroke_width: float = config_dict['objects']['gc_content']['stroke_width']
        self.fill_opacity: float = config_dict['objects']['gc_content']['fill_opacity']
        self.show_gc: bool = show_gc
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

    def __init__(self, window: int, step: int, dinucleotide: str, config_dict: Dict, default_colors_df: DataFrame, show_skew: bool) -> None:
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
        self.high_fill_color: str = default_colors_df[default_colors_df['feature_type'] == 'skew_high']['color'].values[0]
        self.low_fill_color: str = default_colors_df[default_colors_df['feature_type'] == 'skew_low']['color'].values[0]
        self.stroke_color: str = config_dict['objects']['gc_skew']['stroke_color']
        self.stroke_width: float = config_dict['objects']['gc_skew']['stroke_width']
        self.fill_opacity: float = config_dict['objects']['gc_skew']['fill_opacity']
        self.show_skew: bool = show_skew
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
    def __init__(self, color_table: Optional[DataFrame], default_colors: DataFrame, selected_features_set: List[str], config_dict: Dict) -> None:
        """
        Initializes the FeatureDrawingConfigurator with color settings and feature selection.

        Args:
            color_table (Optional[DataFrame]): Color table for features.
            default_colors (Optional[DataFrame]): Default colors for features.
            selected_features_set (str): Set identifier for selecting features to display.
        """
        self.color_table: Optional[DataFrame] = color_table
        self.default_colors: DataFrame = default_colors
        self.selected_features_set: List(str) = selected_features_set
        self.block_fill_color: str =  default_colors[default_colors['feature_type'] == 'default']['color'].values[0]
        self.block_stroke_color: str = config_dict['objects']['features']['block_stroke_color']
        self.block_stroke_width: float = config_dict['objects']['features']['block_stroke_width']
        self.line_stroke_color: str = config_dict['objects']['features']['line_stroke_color']
        self.line_stroke_width: float = config_dict['objects']['features']['line_stroke_width']


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

    def __init__(self, evalue: float, bitscore: float, identity: float, sequence_length_dict: Dict[str, int], config_dict: Dict, default_colors_df: DataFrame) -> None:
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
        self.fill_color: str =  default_colors_df[default_colors_df['feature_type'] == 'pairwise_match']['color'].values[0]
        self.fill_opacity: float = config_dict['objects']['blast_match']['fill_opacity']
        self.stroke_color: str = config_dict['objects']['blast_match']['stroke_color']
        self.stroke_width: float = config_dict['objects']['blast_match']['stroke_width']

class LegendDrawingConfigurator:
    def __init__(self, color_table, default_colors, selected_features_set, config_dict, gc_config, skew_config, feature_config, show_gc, show_skew):
        self.color_table = color_table
        self.default_colors = default_colors
        self.selected_features_set = selected_features_set
        self.config_dict = config_dict
        self.show_gc = show_gc
        self.show_skew = show_skew
        self.gc_config = gc_config
        self.skew_config = skew_config
        self.feature_config = feature_config
        self.font_size: str = config_dict['objects']['legends']['font_size']
        self.font_weight: str = config_dict['objects']['legends']['font_weight']
        self.font_family: str = config_dict['objects']['text']['font_family']
        self.text_anchor: str = config_dict['objects']['legends']['text_anchor']
        self.dominant_baseline: str = config_dict['objects']['legends']['dominant_baseline']