#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from typing import List, Union, Tuple, Dict, Optional
from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.shapes import Line
from svgwrite.text import Text
from svgwrite.gradients import LinearGradient
from .canvas_generator import LinearCanvasConfigurator
from .canvas_generator import LinearCanvasConfigurator
from .linear_feature_drawer import FeatureDrawer, GcContentDrawer, SkewDrawer, LabelDrawer
from .data_processing import skew_df, prepare_label_list_linear
from .create_feature_objects import create_feature_dict, preprocess_color_tables
from .utility_functions import create_text_element, normalize_position_linear, preprocess_label_filtering, calculate_bbox_dimensions, interpolate_color
from .object_configurators import GcContentConfigurator, FeatureDrawingConfigurator, GcSkewConfigurator
from .circular_path_drawer import generate_text_path
import math
# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


class DefinitionGroup:
    """
    Handles the creation and display of definition details for a SeqRecord in a linear layout.

    Attributes:
        record (SeqRecord): GenBank record containing genomic data.
        title_start_x (float): Horizontal start position for the title text.
        title_start_y (float): Vertical start position for the title text.
        length_start_x (float): Horizontal start position for the length text.
        length_start_y (float): Vertical start position for the length text.
        linear_definition_stroke (float): Stroke width for the text.
        linear_definition_fill (str): Fill color for the text.
        linear_definition_font_size (str): Font size for the text.
        linear_definition_font_weight (str): Font weight for the text.
        linear_definition_font_family (str): Font family for the text.
        linear_text_anchor (str): Text anchor position.
        linear_dominant_baseline (str): Dominant baseline position for the text.
    """

    def __init__(self, record: SeqRecord, config_dict: dict, title_start_x: float = 0, title_start_y: float = -9, length_start_x: float = 0, length_start_y: float = 9) -> None:
        """
        Initializes the DefinitionGroup with the given SeqRecord and configuration.

        Args:
            record (SeqRecord): The genomic record to be used for the definition details.
            config_dict (dict): Configuration dictionary with styling parameters.
            title_start_x (float): Horizontal start position for the title text.
            title_start_y (float): Vertical start position for the title text.
            length_start_x (float): Horizontal start position for the length text.
            length_start_y (float): Vertical start position for the length text.
        """
        self.record: SeqRecord = record
        self.title_start_x: float = title_start_x
        self.title_start_y: float = title_start_y
        self.length_start_x: float = length_start_x
        self.length_start_y: float = length_start_y
        self.definition_bounding_box_width: float = 0
        self.definition_bounding_box_height: float = 0
        self.linear_definition_stroke: float = config_dict['objects']['definition']['linear']['stroke']
        self.linear_definition_fill: str = config_dict['objects']['definition']['linear']['fill']
        self.linear_definition_font_size: str = config_dict['objects']['definition']['linear']['font_size']
        self.linear_definition_font_weight: str = config_dict['objects']['definition']['linear']['font_weight']
        self.linear_definition_font_family: str = config_dict['objects']['text']['font_family']
        self.interval: int = config_dict['objects']['definition']['linear']['interval']
        self.dpi: int = config_dict['canvas']['dpi']
        self.linear_text_anchor: str = config_dict['objects']['definition']['linear']['text_anchor']
        self.linear_dominant_baseline: str = config_dict['objects']['definition']['linear']['dominant_baseline']
        self.get_id_and_length()
        self.name_bounding_box_width, self.name_bounding_box_height = calculate_bbox_dimensions(self.record_name, self.linear_definition_font_family, self.linear_definition_font_size, self.dpi)
        self.length_bounding_box_width, self.length_bounding_box_height = calculate_bbox_dimensions(self.length_label, self.linear_definition_font_family, self.linear_definition_font_size, self.dpi)
        self.calculate_start_coordinates()
        self.definition_group = Group(id=self.track_id)
        self.add_elements_to_group()

    def calculate_start_coordinates(self) -> None:
        max_width = max(self.name_bounding_box_width, self.length_bounding_box_width)
        total_height = self.name_bounding_box_height + self.length_bounding_box_height + self.interval
        self.definition_bounding_box_width = max_width
        self.definition_bounding_box_height = total_height
        self.title_start_x = 0
        self.title_start_y = - (total_height / 2) + (self.name_bounding_box_height / 2)
        self.length_start_x = 0
        self.length_start_y = (total_height / 2) - (self.length_bounding_box_height / 2)


    def get_id_and_length(self) -> None:
        """
        Extracts the ID and length of the SeqRecord for display purposes.
        """
        self.track_id = str(self.record.id)
        self.record_name: str = self.track_id
        self.record_length: int = len(self.record.seq)
        self.length_label: str = "{:,} bp".format(self.record_length)

    def add_elements_to_group(self) -> None:
        """
        Adds the definition elements (like record name and length) to the group.
        """
        self.name_path: Text = create_text_element(self.record_name, self.title_start_x, self.title_start_y,
                                                   self.linear_definition_font_size, self.linear_definition_font_weight, self.linear_definition_font_family, text_anchor=self.linear_text_anchor, dominant_baseline=self.linear_dominant_baseline)
        self.length_path: Text = create_text_element(self.length_label, self.length_start_x, self.length_start_y,
                                                     self.linear_definition_font_size, self.linear_definition_font_weight, self.linear_definition_font_family, text_anchor=self.linear_text_anchor, dominant_baseline=self.linear_dominant_baseline)
        self.definition_group.add(self.name_path)
        self.definition_group.add(self.length_path)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the definition details.

        Returns:
            Group: The SVG group with definition details.
        """
        return self.definition_group


class LengthBarGroup:
    """
    Handles the creation and display of a length scale in a linear layout.
    Supports two styles: 'bar' (a single scale bar) and 'ruler' (a full-axis ruler).
    """

    def __init__(self, fig_width: int, alignment_width: float, longest_genome: int, config_dict: dict, group_id="length_bar") -> None:
        """
        Initializes the LengthBarGroup with the given parameters.
        """
        # --- 1. Load all settings from the 'scale' section of the config ---
        scale_config = config_dict['objects']['scale']
        self.stroke_color: str = scale_config['stroke_color']
        self.stroke_width: float = scale_config['stroke_width']
        self.font_size: str = scale_config['font_size']
        self.font_weight: str = scale_config['font_weight']
        self.font_family: str = config_dict['objects']['text']['font_family']
        self.style: str = scale_config.get('style', 'bar')
        self.manual_interval: Optional[int] = scale_config.get('interval')
        self.scale_group_width: float = 0
        self.scale_group_height: float = 0
        self.dpi =  config_dict['canvas']['dpi']
        # --- 2. Set other properties ---
        self.longest_genome: int = longest_genome
        self.alignment_width: float = alignment_width
        self.group_id: str = group_id
        
        self.scale_group = Group(id=self.group_id)

        # --- 3. Branch to the appropriate setup method based on style ---
        if self.style == 'ruler':
            self.setup_scale_ruler()
        else:
            self.setup_scale_bar()

    def setup_scale_bar(self) -> None:
        """
        Sets up the 'bar' style.
        Uses manual_interval if provided, otherwise calculates an automatic interval.
        """
        if self.manual_interval is not None and self.manual_interval > 0:
            # Use the user-provided interval for the bar length
            self.tick = self.manual_interval
            self.label_text = self._format_tick_label(self.manual_interval, for_bar_style=True)
        else:
            # Automatically determine the best interval and label
            self.define_ticks_by_length_for_bar()
        
        # Configure geometry and add elements to the SVG group
        self.config_bar_geometry()
        self.add_bar_elements_to_group()

    def setup_scale_ruler(self) -> None:
        """
        Sets up the 'ruler' style.
        Uses manual_interval if provided, otherwise falls back to an automatic interval.
        """
        # Determine the tick interval
        if self.manual_interval is not None and self.manual_interval > 0:
            tick_interval = self.manual_interval
        else:
            self.define_ticks_by_length_for_bar() # Use the same auto-logic to get a sensible interval
            tick_interval = self.tick

        if tick_interval <= 0 or self.longest_genome <= 0:
            return 

        # Draw the main axis line
        main_axis = Line(
            start=(0, 0),
            end=(self.alignment_width, 0),
            stroke=self.stroke_color,
            stroke_width=self.stroke_width,
        )
        self.scale_group.add(main_axis)
        scale_ruler_length = self.alignment_width

        # Draw ticks and labels from 0 up to the end of the genome
        position = 0
        first_tick_bbox_width = 0
        last_tick_bbox_width = 0
        max_tick_bbox_height = 0
        while True:
            x_pos = (position / self.longest_genome) * self.alignment_width

            # Draw tick line
            tick_line = Line(
                start=(x_pos, 0-(0.5 * self.stroke_width)), end=(x_pos, 10), 
                stroke=self.stroke_color, stroke_width=self.stroke_width
            )
            self.scale_group.add(tick_line)

            # Draw label
            label_text = self._format_tick_label(position)
            text_element = Text(
                label_text, insert=(x_pos, 15), stroke='none', fill='black',
                font_size=self.font_size, font_weight=self.font_weight,
                font_family=self.font_family, text_anchor="middle",
                dominant_baseline="hanging"
            )
            self.scale_group.add(text_element)
            bbox_width, bbox_height = calculate_bbox_dimensions(label_text, self.font_family, self.font_size, self.dpi)
            max_tick_bbox_height = max(max_tick_bbox_height, bbox_height)
            if first_tick_bbox_width == 0:
                first_tick_bbox_width = bbox_width
            if position >= self.longest_genome:
                last_tick_bbox_width = bbox_width
                self.scale_group_width = (0.5 * first_tick_bbox_width) + scale_ruler_length + (0.5 * last_tick_bbox_width)
                self.scale_group_height: float = (0.5 * self.stroke_width) + 15 + max_tick_bbox_height
                break # Exit after drawing the last tick
            
            position += tick_interval
            # Ensure the final tick is exactly at the end
            if position > self.longest_genome:
                position = self.longest_genome

    def _format_tick_label(self, position: int, for_bar_style: bool = False) -> str:
        """
        Formats a position number into a readable label (e.g., "0.5M", "10kbp").
        """
        if position == 0:
            return "0"

        unit = "bp"
        
        if position >= 1_000_000:
            label = f"{position / 1_000_000:.1f} " + ("M" + unit)
        elif not for_bar_style and position >= 100_000: # For ruler style, use M for values like 500k
            if self.longest_genome >= 1_000_000:
                label = f"{position / 1_000_000:.1f} " + "Mbp"
            else:
                label = f"{position // 1_000} k" + (unit)
        elif position >= 1_000:
            label = f"{position // 1_000} k" + (unit)
        else:
            label = f"{position} {unit}".strip()

        return label

    def define_ticks_by_length_for_bar(self) -> None:
        """
        Automatically determines a sensible tick interval for the 'bar' style
        or as a fallback for the 'ruler' style.
        """
        thresholds = [
            (2000, 100), (20000, 1000), (50000, 5000), (150000, 10000),
            (250000, 50000), (1000000, 100000), (2000000, 200000), 
            (5000000, 500000), (float('inf'), 1000000),
        ]
        for threshold, tick_val in thresholds:
            if self.longest_genome < threshold:
                self.tick = tick_val
                self.label_text = self._format_tick_label(tick_val, for_bar_style=True)
                return

    def config_bar_geometry(self) -> None:
        """Configures the geometry (length and position) for the 'bar' style."""
        if self.longest_genome > 0:
            self.bar_length = self.alignment_width * (self.tick / self.longest_genome)
        else:
            self.bar_length = 0
        self.end_x = self.alignment_width
        self.start_x = self.end_x - self.bar_length
        self.start_y = 0
        self.end_y = 0

    def create_scale_path_linear(self) -> Line:
        """Creates the SVG line element for the 'bar' style."""
        return Line(
            start=(self.start_x, self.start_y), end=(self.end_x, self.end_y),
            stroke=self.stroke_color, stroke_width=self.stroke_width, fill='none'
        )

    def create_scale_text_linear(self) -> Text:
        """Creates the SVG text element for the 'bar' style label."""
        self.tick_bbox_width, self.tick_bbox_height = calculate_bbox_dimensions(self.label_text, self.font_family, self.font_size, self.dpi)
        self.scale_group_width =  self.tick_bbox_width + 10 + self.bar_length
        self.scale_group_height = self.tick_bbox_height
        return Text(
            self.label_text, insert=(self.start_x - 10, self.start_y),
            stroke='none', fill='black', font_size=self.font_size,
            font_weight=self.font_weight, font_family=self.font_family,
            text_anchor="end", dominant_baseline="middle"
        )

    def add_bar_elements_to_group(self) -> None:
        """Adds the line and text elements for the 'bar' style to the SVG group."""
        scale_path = self.create_scale_path_linear()
        scale_text = self.create_scale_text_linear()

        self.scale_group.add(scale_path)
        self.scale_group.add(scale_text)

    def get_group(self) -> Group:
        """Retrieves the final SVG group containing the scale visualization."""
        return self.scale_group


class GcContentGroup:
    """
    Manages the visualization of GC content for a genomic sequence in a linear layout.

    Attributes:
        gb_record (SeqRecord): GenBank record containing genomic data.
        longest_record_len (int): Length of the longest record for normalization.
        track_height (float): Height of the GC content track.
        alignment_width (float): Width of the alignment area.
        gc_config (GcContentConfigurator): Configuration for GC content visualization.
        config_dict (dict): Configuration dictionary with styling parameters.
        start_x (float): Starting x-coordinate for the GC content visualization.
        start_y (float): Starting y-coordinate for the GC content visualization.
    """

    def __init__(self, gb_record: SeqRecord, longest_record_len: int, track_height: float, alignment_width: float, gc_config: GcContentConfigurator, config_dict: dict, start_x: float = 0, start_y: float = 0) -> None:
        """
        Initializes the GcContentGroup with the given parameters and configurations.

        Args:
            gb_record (SeqRecord): GenBank record containing genomic data.
            longest_record_len (int): Length of the longest record for normalization.
            track_height (float): Height of the GC content track.
            alignment_width (float): Width of the alignment area.
            gc_config (GcContentConfigurator): Configuration for GC content visualization.
            config_dict (dict): Configuration dictionary with styling parameters.
            start_x (float): Starting x-coordinate for the GC content visualization.
            start_y (float): Starting y-coordinate for the GC content visualization.
        """
        self.gc_group = Group(id="gc_content")
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.longest_record_len: int = longest_record_len
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.gc_config.window
        self.step: int = self.gc_config.step
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.track_height: float = track_height
        self.alignment_width: float = alignment_width
        self.normalize_length()
        self.generate_gc_df()
        self.add_elements_to_group()

    def normalize_length(self) -> None:
        """
        Normalizes the length of the genomic record relative to the longest record.

        This method calculates the normalization factor based on the length of the current 
        genomic record and the length of the longest record in the dataset. This factor is used 
        to scale the visualization appropriately.
        """
        self.record_len: int = len(self.gb_record.seq)
        self.genome_size_normalization_factor: float = self.record_len / self.longest_record_len

    def generate_gc_df(self) -> None:
        """
        Generates a DataFrame for GC content based on the genomic record.

        This method uses the `skew_df` function to calculate the GC content of the genomic record.
        The resulting DataFrame is used for visualizing the GC content.
        """
        self.gc_df: DataFrame = skew_df(
            self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC content visualization elements to the SVG group.

        This method calls the GcContentDrawer to draw the GC content based on the calculated 
        DataFrame and adds the resulting SVG elements to the group.
        """
        self.gc_group: Group = GcContentDrawer(self.gc_config).draw(self.gc_group, self.gc_df, self.record_len, self.alignment_width, self.genome_size_normalization_factor, self.track_height, self.start_x, self.start_y, self.dinucleotide)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC content visualization.

        Returns:
            Group: The SVG group with GC content visualization elements.
        """
        return self.gc_group

class GcSkewGroup:
    """
    Manages the visualization of GC skew for a genomic sequence in a linear layout.
    (This is the new independent class)
    """

    def __init__(self, gb_record: SeqRecord, longest_record_len: int, track_height: float, alignment_width: float, skew_config: GcSkewConfigurator, config_dict: dict, start_x: float = 0, start_y: float = 0) -> None:
        """
        Initializes the GcSkewGroup with the given parameters and configurations.
        """
        self.skew_group = Group(id="gc_skew")
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.longest_record_len: int = longest_record_len
        self.skew_config: GcSkewConfigurator = skew_config
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.skew_config.window
        self.step: int = self.skew_config.step
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.track_height: float = track_height
        self.alignment_width: float = alignment_width
        self.generate_gc_df()
        self.normalize_length()
        self.add_elements_to_group()

    def normalize_length(self) -> None:
        """
        Normalizes the length of the genomic record relative to the longest record.
        """
        self.record_len: int = len(self.gb_record.seq)
        self.genome_size_normalization_factor: float = self.record_len / self.longest_record_len
    def generate_gc_df(self) -> None:

        self.skew_df: DataFrame = skew_df(
            self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC skew visualization elements to the SVG group.
        """
        self.skew_group: Group = SkewDrawer(self.skew_config).draw(
            self.skew_group, self.skew_df, self.record_len, self.alignment_width,
            self.genome_size_normalization_factor, self.track_height, self.start_x, self.start_y, self.dinucleotide
        )

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC skew visualization.
        """
        return self.skew_group

class SeqRecordGroup:
    """Manages the visualization of a SeqRecord in a linear layout."""

    def __init__(self, gb_record: SeqRecord, canvas_config: LinearCanvasConfigurator, 
                 feature_config: FeatureDrawingConfigurator, config_dict: dict, precalculated_labels: Optional[list] = None) -> None:
        self.gb_record = gb_record
        self.canvas_config = canvas_config
        self.feature_config = feature_config
        self.config_dict = config_dict
        self.precalculated_labels = precalculated_labels
        self.show_labels = self.config_dict['canvas']['show_labels']
        self.label_stroke_color = self.config_dict['labels']['stroke_color']['label_stroke_color']
        self.label_stroke_width = self.config_dict['labels']['stroke_width']['long']
        self.label_filtering = self.config_dict['labels']['filtering']
        self.separate_strands = self.canvas_config.strandedness
        self.resolve_overlaps = self.canvas_config.resolve_overlaps
        self.record_group: Group = self.setup_record_group()
    def draw_linear_axis(self,
                         alignment_width: float,
                         genome_size_normalization_factor: float) -> Line:
        """
        Draws a linear axis for the genomic record visualization.

        Args:
            alignment_width (float): Width of the alignment area.
            genome_size_normalization_factor (float): Normalization factor based on the genome size.

        Returns:
            Line: An SVG line element representing the linear axis.
        """
        bar_length: float = alignment_width * genome_size_normalization_factor
        linear_axis_stroke_color: str = self.config_dict['objects']['axis']['linear']['stroke_color']
        linear_axis_stroke_width: float = self.config_dict['objects']['axis']['linear']['stroke_width']
        start_x: float = 0
        start_y: float = 0
        end_x: float = bar_length
        end_y = 0
        axis_path = Line(
            start=(
                start_x,
                start_y),
            end=(
                end_x,
                end_y),
            stroke=linear_axis_stroke_color,
            stroke_width=linear_axis_stroke_width,
            fill='none')
        return axis_path

    def draw_record(self, feature_dict: dict, record_length: int, cds_height: float, 
                   alignment_width: float, genome_size_normalization_factor: float, 
                   separate_strands: bool, arrow_length: float, group: Group, label_list: list) -> Group:
        """Draws the genomic features onto the provided SVG group."""
        # Draw the axis
        axis_path = self.draw_linear_axis(alignment_width, genome_size_normalization_factor)
        group.add(axis_path)

        # Process labels if enabled
        if self.show_labels:
            for label in label_list:
                if not label["is_embedded"]:
                    line_path = Line(
                        start=(label["middle"], label["feature_middle_y"]),
                        end=(label["middle"], label["middle_y"] + 0.45 * label["height_px"]),
                        stroke=self.label_stroke_color,
                        stroke_width=self.label_stroke_width
                    )
                    group.add(line_path)

        # Draw features
        for feature_object in feature_dict.values():
            feature_strand = feature_object.location[0][2]
            group = FeatureDrawer(self.feature_config).draw(
                feature_object=feature_object,
                group=group,
                genome_length=record_length,
                cds_height=cds_height,
                alignment_width=alignment_width,
                normalization_factor=genome_size_normalization_factor,
                feature_strand=feature_strand,
                separate_strands=separate_strands,
                arrow_length=arrow_length
            )

        # Add labels
        if self.show_labels:
            for label in label_list:
                group = LabelDrawer(self.config_dict).draw(label, group)

        return group

    def setup_record_group(self) -> Group:
        """
        Sets up the SVG group for the SeqRecord visualization.

        Returns:
            Group: The SVG group prepared for SeqRecord visualization.
        """
        alignment_width: float = self.canvas_config.alignment_width
        cds_height: float = self.canvas_config.cds_height
        longest_genome: int = self.canvas_config.longest_genome
        arrow_length: float = self.canvas_config.arrow_length
        color_table: DataFrame | None = self.feature_config.color_table

        separate_strands = self.canvas_config.strandedness
        resolve_overlaps = self.canvas_config.resolve_overlaps
        label_filtering = self.label_filtering  # type: ignore
        track_id = self.gb_record.id
        
        record_group: Group = Group(id=track_id)
        record_group = Group(id=track_id)
        record_length: int = len(self.gb_record.seq)
        genome_size_normalization_factor: float = record_length / longest_genome
        selected_features_set: str = self.feature_config.selected_features_set
        
        default_colors: DataFrame | None = self.feature_config.default_colors
        label_filtering = preprocess_label_filtering(self.label_filtering)
        color_table, default_colors = preprocess_color_tables(color_table, default_colors)
        feature_dict: dict = create_feature_dict(self.gb_record, color_table, selected_features_set, default_colors, separate_strands, resolve_overlaps, label_filtering)
        label_list = []
        if self.show_labels:
            label_list = prepare_label_list_linear(
                feature_dict, record_length, alignment_width,
                genome_size_normalization_factor, cds_height,
                separate_strands, self.config_dict
            )
        
        record_group: Group = self.draw_record(feature_dict, record_length, cds_height, alignment_width, genome_size_normalization_factor, separate_strands, arrow_length, record_group, label_list)
        return record_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the SeqRecord visualization.

        Returns:
            Group: The SVG group with SeqRecord visualization elements.
        """
        return self.record_group


class PairWiseMatchGroup:
    """
    Handles the visualization of pairwise matches (comparisons) in a linear layout.

    Attributes:
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        sequence_length_dict (Dict[str, int]): Dictionary of sequence lengths.
        comparison_df (DataFrame): DataFrame containing comparison data.
        comparison_height (float): Height of the comparison track.
        comparison_count (int): Counter for the number of comparisons.
        config_dict (dict): Configuration dictionary with styling parameters.
    """

    def __init__(self, canvas_config: LinearCanvasConfigurator, sequence_length_dict: dict, comparison_df: DataFrame, actual_comparison_height: float, comparison_count: int, blast_config, records) -> None:
        """
        Initializes the PairWiseMatchGroup with necessary data and configurations.

        Args:
            canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
            sequence_length_dict (Dict[str, int]): Dictionary of sequence lengths.
            comparison_df (DataFrame): DataFrame containing comparison data.
            comparison_height (float): Height of the comparison track.
            comparison_count (int): Counter for the number of comparisons.
            config_dict (dict): Configuration dictionary with styling parameters.
        """
        self.canvas_config: LinearCanvasConfigurator = canvas_config
        self.sequence_length_dict: Dict[str, int] = sequence_length_dict
        self.comparison_df: DataFrame = comparison_df
        self.comparison_height: float = actual_comparison_height
        self.match_fill_color: str = blast_config.fill_color
        self.min_identity: float = blast_config.identity
        self.match_min_color: str = blast_config.min_color
        self.match_max_color: str = blast_config.max_color
        self.match_fill_opacity: float = blast_config.fill_opacity
        self.match_stroke_color: str = blast_config.stroke_color
        self.match_stroke_width: float = blast_config.stroke_width
        self.comparison_count: int = comparison_count
        self.records = records
        self.track_id: str = "comparison" + str(self.comparison_count)
        self.calculate_query_subject_offsets()
        self.match_group = Group(id=self.track_id)
        self.add_elements_to_group()


    def calculate_query_subject_offsets(self) -> Tuple[float, float]:
        """
        Calculates the horizontal offsets for the query and subject sequences.

        This method determines the starting x-coordinates for the query and subject in the alignment, 
        ensuring they are centered if specified in the configuration.

        Args:
            query_id (str): Identifier for the query sequence.
            subject_id (str): Identifier for the subject sequence.

        Returns:
            Tuple[float, float]: The x-coordinate offsets for the query and subject sequences.
        """
        if self.canvas_config.align_center:
            qlen = len(self.records[self.comparison_count-1].seq)
            slen = len(self.records[self.comparison_count].seq)
            self.query_offset_x = (self.canvas_config.longest_genome - qlen) / 2
            self.subject_offset_x = (self.canvas_config.longest_genome - slen) / 2
        else:
            self.query_offset_x = self.subject_offset_x = 0

    def generate_linear_match_path(self, row: DataFrame) -> Path:
        """
        Generates an SVG path for a pairwise match based on the provided row data.

        Args:
            row (DataFrame): A row from the DataFrame containing match information.

        Returns:
            Path: An SVG path element representing the pairwise match.
        """
        query_start: float
        query_end: float
        subject_start: float
        subject_end: float
        query_start_x: float
        query_start_y: float
        query_end_x: float
        query_end_y: float
        subject_start_x: float
        subject_start_y: float
        subject_end_x: float
        subject_end_y: float
        identity_percent = float(row.identity)
        factor = (identity_percent - self.min_identity) / (100 - self.min_identity)
        dynamic_fill_color = interpolate_color(self.match_min_color, self.match_max_color, factor)
        query_start, query_end, subject_start, subject_end = self.calculate_offsets(row)
        query_start_x, query_start_y, query_end_x, query_end_y = self.normalize_positions(
            query_start, query_end, 0)
        subject_start_x, subject_start_y, subject_end_x, subject_end_y = self.normalize_positions(
            subject_start, subject_end, self.comparison_height)

        match_path_desc: str = self.construct_path_description(
            query_start_x, query_start_y, query_end_x, query_end_y, subject_start_x, subject_start_y, subject_end_x, subject_end_y)
        return Path(
            d=match_path_desc,
            fill=dynamic_fill_color,
            fill_opacity=self.match_fill_opacity,
            stroke=self.match_stroke_color,
            stroke_width=self.match_stroke_width)

    def calculate_offsets(self, row: DataFrame) -> tuple[float, float, float, float]:
        """
        Calculates the start and end positions for a pairwise match.

        Args:
            row (DataFrame): A row from the DataFrame containing match information.

        Returns:
            Tuple[float, float, float, float]: Start and end positions for the query and subject matches.
        """
        query_start: float = row.qstart + self.query_offset_x
        query_end: float = row.qend + self.query_offset_x
        subject_start: float = row.sstart + self.subject_offset_x
        subject_end: float = row.send + self.subject_offset_x
        return query_start, query_end, subject_start, subject_end

    def normalize_positions(self, start: float, end: float, y_position: float) -> tuple[float, float, float, float]:
        """
        Normalizes the start and end positions for display on the linear canvas.

        Args:
            start (float): The start position of the match.
            end (float): The end position of the match.
            y_position (float): The y-coordinate for the match.

        Returns:
            Tuple[float, float, float, float]: Normalized start and end x-coordinates and fixed y-coordinates.
        """
        start_x: float = normalize_position_linear(
            start, self.canvas_config.longest_genome, self.canvas_config.alignment_width)
        end_x: float = normalize_position_linear(
            end, self.canvas_config.longest_genome, self.canvas_config.alignment_width)
        return start_x, y_position, end_x, y_position

    def construct_path_description(self, query_start_x: float, query_start_y: float, query_end_x: float, query_end_y: float, subject_start_x: float, subject_start_y: float, subject_end_x: float, subject_end_y: float) -> str:
        """
        Constructs the path description for an SVG path element representing a pairwise match.

        Args:
            query_start_x (float): x-coordinate of the query start.
            query_start_y (float): y-coordinate of the query start.
            query_end_x (float): x-coordinate of the query end.
            query_end_y (float): y-coordinate of the query end.
            subject_start_x (float): x-coordinate of the subject start.
            subject_start_y (float): y-coordinate of the subject start.
            subject_end_x (float): x-coordinate of the subject end.
            subject_end_y (float): y-coordinate of the subject end.

        Returns:
            str: SVG path description string.
        """
        return f"M {query_start_x},{query_start_y}L{query_end_x},{query_end_y} L{subject_end_x},{subject_end_y}L{subject_start_x},{subject_start_y} z"

    def add_elements_to_group(self) -> Group:
        """
        Adds all the pairwise match paths to the SVG group.

        Iterates through the pairwise matches and creates SVG path elements for each, 
        adding them to the group.

        Returns:
            Group: The SVG group with all match paths added.
        """
        for row in self.comparison_df.itertuples():
            match_path: Path = self.generate_linear_match_path(row)
            self.match_group.add(match_path)
        return self.match_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing all the pairwise match visualizations.

        Returns:
            Group: The SVG group with pairwise match elements.
        """
        return self.match_group

class LegendGroup:
    def __init__(self, config_dict, canvas_config, legend_config, legend_table):
        self.legend_group = Group(id="legend")
        self.config_dict = config_dict
        self.canvas_config = canvas_config
        self.legend_config = legend_config
        self.legend_table = legend_table
        self.font_family: str = self.legend_config.font_family
        self.font_weight: str = self.legend_config.font_weight
        self.font_size: float = self.legend_config.font_size
        self.rect_size: float = self.legend_config.color_rect_size
        self.legend_position = self.canvas_config.legend_position
        self.line_height: float = (24/14) * self.rect_size
        self.text_x_offset: float = (22/14) * self.rect_size
        self.num_of_columns = self.legend_config.num_of_columns
        self.has_gradient = self.legend_config.has_gradient
        self.total_feature_legend_width = self.legend_config.total_feature_legend_width
        self.dpi: int = self.config_dict['canvas']['dpi']
        self.legend_width: float = 0
        self.feature_legend_width: float = 0
        self.feature_legend_height: float = 0
        self.pairwise_legend_width = self.legend_config.pairwise_legend_width
        self.pairwise_legend_height: float = 0
        self.legend_height: float = 0
        self.add_elements_to_group()

    def create_rectangle_path_for_legend(self) -> str:
        start_y_top = -self.rect_size / 2
        start_y_bottom = self.rect_size / 2
        return f"M 0,{start_y_top} L {self.rect_size},{start_y_top} L {self.rect_size},{start_y_bottom} L 0,{start_y_bottom} z"

    def add_elements_to_group(self):
        self.feature_legend_height = self.rect_size
        y_offset = self.rect_size / 2
        initial_y_offset = y_offset
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        feature_legend_group = Group(id="feature_legend")
        max_bbox_width = 0
        gradient_y_offset = 0
        grad_bar_width = 0
        current_feature_legend_width = 0
        pairwise_legend_group = Group(id="pairwise_legend")
        if self.num_of_columns > 1:
            if self.has_gradient:
                max_items_per_column = self.num_of_columns - 1
            else:
                max_items_per_column = self.num_of_columns
            current_column = 0
            current_x_offset = 0
            for key, properties in self.legend_table.items():
                if properties['type'] == 'solid':
                    bbox_width, _ = calculate_bbox_dimensions(str(key), self.font_family, self.font_size, self.dpi)
                    max_bbox_width = max(max_bbox_width, bbox_width)
                    if current_feature_legend_width + self.text_x_offset + bbox_width + self.text_x_offset > self.total_feature_legend_width:
                        current_feature_legend_width = 0
                        current_column = 0
                        current_x_offset = 0
                        y_offset += self.line_height
                        self.feature_legend_height += self.line_height
                    rect_path = Path(
                        d=path_desc,
                        fill=properties['fill'],
                        stroke=properties['stroke'],
                        stroke_width=properties['width'])
                    rect_path.translate(current_x_offset, y_offset)
                    feature_legend_group.add(rect_path)
                    current_x_offset +=  self.text_x_offset
                    legend_path = generate_text_path(
                            key, 0, 0, 0, self.font_size, "normal", font, 
                            dominant_baseline='central', text_anchor="start"
                        )
                    legend_path.translate(current_x_offset, y_offset)
                    feature_legend_group.add(legend_path)
                    current_x_offset += (bbox_width + self.text_x_offset)
                    current_feature_legend_width = current_x_offset 
                    self.feature_legend_width = max(current_feature_legend_width, self.feature_legend_width)
                    current_column += 1

                elif properties['type'] == 'gradient':
                    grad_bar_x_start = 0
                    title_x_offset = 0
                    grad_bar_width = self.legend_config.pairwise_legend_width
                    gradient_y_offset = 0
                    gradient_id = f"blast_legend_grad_{abs(hash(properties['min_color'] + properties['max_color']))}"
                    
                    gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                    gradient.add_stop_color(offset="0%", color=properties['min_color'])
                    gradient.add_stop_color(offset="100%", color=properties['max_color'])
                    pairwise_legend_group.add(gradient)
                    title_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='hanging', text_anchor="middle"
                    )
                    title_path_bbox_width, title_path_bbox_height = calculate_bbox_dimensions(str(key), self.font_family, self.font_size, self.dpi)
                    title_x_offset = grad_bar_x_start + (grad_bar_width / 2)
                    title_path.translate(title_x_offset, gradient_y_offset) 
                    pairwise_legend_group.add(title_path) 
                    gradient_y_offset += title_path_bbox_height + (self.rect_size / 2)


                    
                    grad_rect_path_desc = f"M 0,{-self.rect_size / 2} L {grad_bar_width},{-self.rect_size / 2} L {grad_bar_width},{self.rect_size / 2} L 0,{self.rect_size / 2} z"
                    grad_rect = Path(
                        d=grad_rect_path_desc,
                        fill=f"url(#{gradient_id})",
                        stroke=properties['stroke'],
                        stroke_width=properties['width']
                    )
                    grad_rect.translate(grad_bar_x_start, gradient_y_offset) 
                    pairwise_legend_group.add(grad_rect)
                    min_identity = properties.get('min_value', 0)
                    if min_identity == int(min_identity):
                        min_label_text = f"{int(min_identity)}%"
                    else:
                        min_label_text = f"{min_identity}%"
                    label_0 = generate_text_path(
                        min_label_text, 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="start"
                     )
                    label_0.translate(grad_bar_x_start, gradient_y_offset + self.rect_size / 2 + 2) 
                    pairwise_legend_group.add(label_0)
                    _, min_label_bbox_height = calculate_bbox_dimensions(min_label_text, self.font_family, self.font_size, self.dpi)

                    label_100 = generate_text_path(
                        "100%", 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="end"
                    )
                    label_100.translate(grad_bar_x_start + grad_bar_width, gradient_y_offset + self.rect_size / 2 + 2) 
                    pairwise_legend_group.add(label_100)
                    _, max_label_bbox_height = calculate_bbox_dimensions("100%", self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_100)
                    label_bbox_height = max(min_label_bbox_height, max_label_bbox_height)
                    gradient_y_offset += label_bbox_height

                    self.pairwise_legend_width = grad_bar_width
                    self.pairwise_legend_height = gradient_y_offset
            if self.pairwise_legend_width >0:
                self.legend_width = self.feature_legend_width + self.pairwise_legend_width + self.text_x_offset
            else:
                self.legend_width = self.feature_legend_width
            self.legend_height = max(self.feature_legend_height, self.pairwise_legend_height)
            
            if self.feature_legend_height > self.pairwise_legend_height:
                pairwise_legend_group.translate(self.feature_legend_width + self.text_x_offset, (self.feature_legend_height - self.pairwise_legend_height)/2)
            else:
                pairwise_legend_group.translate(self.feature_legend_width  + self.text_x_offset, 0)
                feature_legend_group.translate(0, (self.pairwise_legend_height - self.feature_legend_height)/2)

            self.legend_group.add(pairwise_legend_group)
            self.legend_group.add(feature_legend_group)
        else:
            for key, properties in self.legend_table.items():
                if properties['type'] == 'solid':
                    current_feature_legend_width = 0
                    rect_path = Path(
                        d=path_desc,
                        fill=properties['fill'],
                        stroke=properties['stroke'],
                        stroke_width=properties['width'])
                    rect_path.translate(0, y_offset)
                    feature_legend_group.add(rect_path)
                    bbox_width, bbox_height = calculate_bbox_dimensions(key, self.font_family, self.font_size, self.dpi)
                    max_bbox_width = max(max_bbox_width, bbox_width)
                    legend_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='central', text_anchor="start"
                    )
                    legend_path.translate(self.text_x_offset, y_offset)
                    current_feature_legend_width = self.text_x_offset + bbox_width
                    feature_legend_group.add(legend_path)
                    
                    y_offset += self.line_height 
                    self.feature_legend_width = max(current_feature_legend_width, self.feature_legend_width)
                    
                elif properties['type'] == 'gradient':
                    grad_bar_x_start = 0
                    if self.legend_position == 'top' or self.legend_position == 'bottom':
                        gradient_y_offset = 0
                    else:
                        gradient_y_offset = y_offset
                    grad_bar_width = self.pairwise_legend_width
                    gradient_id = f"blast_legend_grad_{abs(hash(properties['min_color'] + properties['max_color']))}"
                    
                    gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                    gradient.add_stop_color(offset="0%", color=properties['min_color'])
                    gradient.add_stop_color(offset="100%", color=properties['max_color'])
                    pairwise_legend_group.add(gradient)
                    title_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='hanging', text_anchor="middle"
                    )
                    title_path_bbox_width, _ = calculate_bbox_dimensions(key, self.font_family, self.font_size, self.dpi)
                    title_x_offset = grad_bar_x_start + (grad_bar_width / 2)

                    title_path.translate(title_x_offset, gradient_y_offset) 
                    pairwise_legend_group.add(title_path) 
                    gradient_y_offset += self.line_height 


                    
                    grad_rect_path_desc = f"M 0,{-self.rect_size / 2} L {grad_bar_width},{-self.rect_size / 2} L {grad_bar_width},{self.rect_size / 2} L 0,{self.rect_size / 2} z"
                    grad_rect = Path(
                        d=grad_rect_path_desc,
                        fill=f"url(#{gradient_id})",
                        stroke=properties['stroke'],
                        stroke_width=properties['width']
                    )
                    grad_rect.translate(grad_bar_x_start, gradient_y_offset) 
                    pairwise_legend_group.add(grad_rect)
                    min_identity = properties.get('min_value', 0)
                    if min_identity == int(min_identity):
                        min_label_text = f"{int(min_identity)}%"
                    else:
                        min_label_text = f"{min_identity}%"
                    label_0 = generate_text_path(
                        min_label_text, 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="start"
                    )
                    label_0.translate(grad_bar_x_start, gradient_y_offset + self.rect_size / 2 + 2)
                    _, min_label_bbox_height = calculate_bbox_dimensions(min_label_text, self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_0)
                    label_100 = generate_text_path(
                        "100%", 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="end"
                    )
                    label_100.translate(grad_bar_x_start + grad_bar_width, gradient_y_offset + self.rect_size / 2 + 2) 
                    _, max_label_bbox_height = calculate_bbox_dimensions("100%", self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_100)
                    label_bbox_height = max(min_label_bbox_height, max_label_bbox_height)
                    gradient_y_offset += label_bbox_height 
            if self.pairwise_legend_width > 0:
                if self.pairwise_legend_width > self.feature_legend_width:
                    feature_legend_group.translate((self.pairwise_legend_width - self.feature_legend_width)/2, 0)
                    self.legend_group.add(pairwise_legend_group)
                else:
                    # pairwise_legend_group.translate((self.feature_legend_width - self.pairwise_legend_width)/2, 0)
                    self.legend_group.add(pairwise_legend_group)
            self.legend_width = max(self.feature_legend_width, self.pairwise_legend_width)
            self.legend_height = max(gradient_y_offset, y_offset) + initial_y_offset
            self.legend_group.add(feature_legend_group)
        return self.legend_group
    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the figure legends.

        Returns:
            Group: The SVG group with figure legends.
        """
        return self.legend_group