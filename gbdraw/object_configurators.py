#!/usr/bin/env python
# coding: utf-8
from typing import Dict, Optional, List
from pandas import DataFrame
from .utility_functions import calculate_bbox_dimensions

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

    def __init__(self, window: int, step: int, dinucleotide: str, config_dict: Dict, default_colors_df: DataFrame) -> None:
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
        self.show_gc: bool = config_dict['canvas']['show_gc']

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

    def __init__(self, window: int, step: int, dinucleotide: str, config_dict: Dict, default_colors_df: DataFrame) -> None:
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
        self.show_skew: bool = config_dict['canvas']['show_skew']
        
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
        self.qualifier_priority: config_dict['labels']['filtering']['qualifier_priority']
        self.blacklist_keywords: List[str] = config_dict['labels']['filtering']['blacklist_keywords']

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
        self.min_color: str = default_colors_df[default_colors_df['feature_type'] == 'pairwise_match_min']['color'].values[0]
        self.max_color: str = default_colors_df[default_colors_df['feature_type'] == 'pairwise_match_max']['color'].values[0]
        self.fill_color: str =  default_colors_df[default_colors_df['feature_type'] == 'pairwise_match']['color'].values[0]
        self.fill_opacity: float = config_dict['objects']['blast_match']['fill_opacity']
        self.stroke_color: str = config_dict['objects']['blast_match']['stroke_color']
        self.stroke_width: float = config_dict['objects']['blast_match']['stroke_width']

class LegendDrawingConfigurator:
    def __init__(self, color_table, default_colors, selected_features_set, config_dict, gc_config, skew_config, feature_config, legend_table=None, blast_config=None) -> None:
        self.color_table = color_table
        self.default_colors = default_colors
        self.selected_features_set = selected_features_set
        self.config_dict = config_dict
        self.show_gc = config_dict['canvas']['show_gc']
        self.show_skew = config_dict['canvas']['show_skew']
        self.gc_config = gc_config
        self.skew_config = skew_config
        self.feature_config = feature_config
        self.blast_config = blast_config
        self.font_size: float = config_dict['objects']['legends']['font_size']
        self.font_weight: str = config_dict['objects']['legends']['font_weight']
        self.font_family: str = config_dict['objects']['text']['font_family']
        self.text_anchor: str = config_dict['objects']['legends']['text_anchor']
        self.color_rect_size: float = config_dict['objects']['legends']['color_rect_size']
        self.dominant_baseline: str = config_dict['objects']['legends']['dominant_baseline']
        self.num_of_columns: int = 1
        self.has_gradient: bool = False
        self.dpi: int = config_dict['png_output']['dpi']
        self.total_feature_legend_width : float = 0
        self.pairwise_legend_width: float = 0
    def calculate_max_bbox_dimensions(self, legend_table):
        longest_key = max(legend_table.keys(), key=len)
        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(longest_key, self.font_family, self.font_size, self.dpi)
        return bbox_width_px, bbox_height_px
    def recalculate_legend_dimensions(self, legend_table, canvas_config):
        line_margin = (24/14) * self.color_rect_size # move to config
        x_margin = (22/14) * self.color_rect_size # move to config
        total_width = canvas_config.fig_width
        for key, properties in legend_table.items():
            if properties.get('type') == 'gradient':
                self.has_gradient = True
        if canvas_config.legend_position == 'top' or canvas_config.legend_position == 'bottom':
            bbox_list = [calculate_bbox_dimensions(item, self.font_family, self.font_size, self.dpi) for item in legend_table]
            print(bbox_list)
            total_feature_legend_width = sum([bbox[0] for bbox in bbox_list]) + x_margin * (len(legend_table) - 1)
            self.pairwise_legend_width = (10 * self.color_rect_size) if self.has_gradient else 0
            if self.has_gradient:
                total_legend_width = total_feature_legend_width + self.pairwise_legend_width
                self.num_of_columns = len(legend_table) - 1
            else:
                total_legend_width = total_feature_legend_width
                self.num_of_columns = len(legend_table)
            if total_legend_width <= total_width:
                self.legend_width = total_legend_width
                self.legend_height = self.color_rect_size + 2 * line_margin
                self.total_feature_legend_width = self.legend_width - self.pairwise_legend_width
                return self
            else:
                # Split into two rows or more. Add one item at a time until width exceeds total_width
                current_width = x_margin
                num_lines = 1
                per_line_item_count = 0
                for item in bbox_list:
                    item_width = item[0]
                    if current_width + item_width + self.pairwise_legend_width + x_margin <= total_width:
                        per_line_item_count += 1
                        current_width += item_width + x_margin + self.pairwise_legend_width
                        if self.has_gradient:
                            current_num_columns = per_line_item_count + 1
                        else:
                            current_num_columns = per_line_item_count
                        if current_num_columns > self.num_of_columns:
                            self.num_of_columns = current_num_columns
                    else:
                        num_lines += 1
                        per_line_item_count = 0
                        current_width = x_margin + item_width + x_margin + self.pairwise_legend_width
                self.legend_width = current_width
                self.total_feature_legend_width = self.legend_width - self.pairwise_legend_width
                self.legend_height = num_lines * (self.color_rect_size + line_margin) + line_margin

                return self
            
        else:
            bbox_width_px, _ = self.calculate_max_bbox_dimensions(legend_table)
            self.legend_width = x_margin + bbox_width_px
            num_lines = 0
            for key, properties in legend_table.items():
                if properties.get('type') == 'gradient':
                    self.has_gradient = True
                else:
                    num_lines += 1
            
            if self.has_gradient:
                num_lines += 2 # グラデーションはタイトルとバーで2行分とみなす
            
            if num_lines > 0:
                self.legend_height = (self.color_rect_size + (num_lines - 1) * line_margin)
            else:
                self.legend_height = 0
            # self.legend_height = (self.color_rect_size + (len(legend_table.keys()) -1) * line_margin)
            return self
        