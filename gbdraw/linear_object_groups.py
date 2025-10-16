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
from .canvas_generator import LinearCanvasConfigurator
from .linear_feature_drawer import FeatureDrawer, GcContentDrawer, SkewDrawer, LabelDrawer
from .data_processing import skew_df, prepare_label_list_linear
from .create_feature_objects import create_feature_dict, preprocess_color_tables
from .utility_functions import create_text_element, normalize_position_linear, preprocess_label_filtering, calculate_bbox_dimensions
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
    Handles the creation and display of a length bar in a linear layout.

    This class is responsible for visualizing a scale bar that represents the length 
    of the genomic sequence in comparison to the longest sequence in the dataset.

    Attributes:
        fig_width (int): The width of the figure.
        longest_genome (int): The length of the longest genome in the dataset.
        config_dict (dict): Configuration dictionary with styling parameters.
        group_id (str): Identifier for the SVG group.
    """

    def __init__(self, fig_width: int, alignment_width: float, longest_genome: int, config_dict: dict, group_id="length_bar") -> None:
        """
        Initializes the LengthBarGroup with the given parameters.
        (snip)
        """
        self.length_bar_stroke_color: str = config_dict['objects']['length_bar']['stroke_color']
        self.length_bar_stroke_width: float = config_dict['objects']['length_bar']['stroke_width']
        self.length_bar_font_size: str = config_dict['objects']['length_bar']['font_size']
        self.length_bar_font_weight: str = config_dict['objects']['length_bar']['font_weight']
        self.length_bar_font_family: str = config_dict['objects']['text']['font_family']
        self.longest_genome: int = longest_genome
        self.fig_width: int = fig_width
        self.alignment_width: float = alignment_width # (追加) alignment_widthを保存
        self.group_id: str = group_id
        self.style: str = 'ruler' if config_dict.get('objects', {}).get('length_bar', {}).get('style') else 'default'
        self.define_ticks_by_length()
        self.config_bar()
        self.length_bar_group = Group(id=self.group_id)

        if self.style == 'ruler':
            self.setup_ruler_bar()
        else:
            self.setup_default_bar()
        # self.add_elements_to_group()
    def setup_default_bar(self) -> None:
            """
            Sets up the original, single-tick length bar.
            """
            self.define_ticks_by_length()
            self.config_bar()
            self.add_elements_to_group()
    def setup_ruler_bar(self) -> None:
            """
            Sets up the new ruler-style length bar using the interval from define_ticks_by_length.
            """
            self.define_ticks_by_length()
            tick_interval = self.tick

            if tick_interval == 0:
                return 

            main_axis = Line(
                start=(0, 0),
                end=(self.alignment_width, 0),
                stroke=self.length_bar_stroke_color,
                stroke_width=self.length_bar_stroke_width,
            )
            self.length_bar_group.add(main_axis)

            num_ticks = int(self.longest_genome / tick_interval)
            
            for i in range(num_ticks + 2):
                position = i * tick_interval

                if position > self.longest_genome:
                    position = self.longest_genome
                
                x_pos = (position / self.longest_genome) * self.alignment_width

                tick_line = Line(
                    start=(x_pos, 0),
                    end=(x_pos, 5), 
                    stroke=self.length_bar_stroke_color,
                    stroke_width=self.length_bar_stroke_width,
                )
                self.length_bar_group.add(tick_line)


                label_text = self._format_tick_label(position)
                

                text_element = Text(
                    label_text,
                    insert=(x_pos, 15), 
                    stroke='none',
                    fill='black',
                    font_size=self.length_bar_font_size,
                    font_weight=self.length_bar_font_weight,
                    font_family=self.length_bar_font_family,
                    text_anchor="middle",
                    dominant_baseline="middle"
                )
                self.length_bar_group.add(text_element)


                if position == self.longest_genome:
                    break

    def _calculate_ruler_ticks(self) -> List[int]:
            """
            Calculates the positions for ticks on the ruler based on the genome length.
            """

            magnitude = 10 ** (math.floor(math.log10(self.longest_genome)) -1)
            
            res = self.longest_genome / magnitude
            if res > 5:
                interval = magnitude
            elif res > 2:
                interval = magnitude / 2
            else:
                interval = magnitude / 5

            interval = round(interval / magnitude) * magnitude
            if interval == 0:
                interval = magnitude

            num_ticks = int(self.longest_genome / interval)
            ticks = [i * interval for i in range(num_ticks + 1)]

            if len(ticks) > 1 and (self.longest_genome - ticks[-1]) < interval / 2:
                ticks[-1] = self.longest_genome
            else:
                ticks.append(self.longest_genome)

            return [int(t) for t in ticks]

    def _format_tick_label(self, position: int) -> str:
        """
        Formats the tick position into a human-readable label (e.g., 10k, 1.5M).
        """
        if position == 0:
            return "0"
        if position >= 1_000_000:
            return f"{position / 1_000_000:.1f} Mb"
        if position >= 1_000:
            return f"{position / 1_000:.0f} kb"
        return str(position)
    
    def define_ticks_by_length(self) -> None:
        """
        Defines the scale ticks for the length bar based on the length of the longest genome.

        This method sets the appropriate tick intervals and the label format based on 
        predefined thresholds. It adapts the scale to fit the length of the longest genome.
        """
        # Define a list of (threshold, tick, format) tuples
        thresholds: List[Tuple[Union[float, int], int, str]] = [
            (2000, 100, "{} bp"),
            (20000, 1000, "{} kbp"),
            (50000, 5000, "{} kbp"),
            (150000, 10000, "{} kbp"),
            (250000, 50000, "{} kbp"),
            (1000000, 100000, "{} kbp"),
            (2000000, 200000, "{} kbp"),
            (5000000, 500000, "{} kbp"),
            (float('inf'), 1000000, "{} Mbp"),
        ]
        # Find the appropriate tick and format
        for threshold, tick, label_format in thresholds:
            if self.longest_genome < threshold:
                if self.longest_genome < 5000000:
                    self.tick: int = tick
                    self.label_text: str = label_format.format(
                        int(self.tick / 1000) if self.tick >= 1000 else self.tick)
                else:
                    self.tick: int = tick
                    self.label_text: str = label_format.format(
                        int(self.tick / 1000000) if self.tick >= 1000000 else self.tick)
                break

    def config_bar(self) -> None:
        """
        Configures the length bar dimensions and positions based on the tick scale.
        """
        self.bar_length: float = self.alignment_width * (self.tick / self.longest_genome)
        self.end_x: float = self.alignment_width
        self.start_x: float = self.end_x - self.bar_length
        self.start_y: float = 0
        self.end_y: float = 0

    def create_length_bar_path_linear(self) -> Line:
        """
        Creates the SVG line element for the length bar.

        Returns:
            Line: An SVG line element representing the length bar.
        """
        return Line(
            start=(self.start_x, self.start_y),
            end=(self.end_x, self.end_y),
            stroke=self.length_bar_stroke_color,
            stroke_width=self.length_bar_stroke_width,
            fill='none')

    def create_length_bar_text_linear(self) -> Text:
        """
        Creates the SVG text element for the length bar label.

        Returns:
            Text: An SVG text element representing the label for the length bar.
        """
        return Text(
            self.label_text,
            insert=(self.start_x - 10, self.start_y),
            stroke='none',
            fill='black',
            font_size=self.length_bar_font_size,
            font_weight=self.length_bar_font_weight,
            font_family=self.length_bar_font_family,
            text_anchor="end",
            dominant_baseline="middle")

    def add_elements_to_group(self) -> None:
        """
        Adds the length bar elements (line and text) to the SVG group.

        This method calls the functions to create the line and text elements for the length bar
        and adds them to the group.
        """
        length_bar_path: Line = self.create_length_bar_path_linear()
        length_bar_text_path: Text = self.create_length_bar_text_linear()
        self.length_bar_group.add(length_bar_path)
        self.length_bar_group.add(length_bar_text_path)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the length bar.

        Returns:
            Group: The SVG group with the length bar elements.
        """
        return self.length_bar_group


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
        self.record_group: Group = self.setup_record_group()
        self.separate_strands = self.canvas_config.strandedness
        self.resolve_overlaps = self.canvas_config.resolve_overlaps

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
            if self.precalculated_labels is not None:
                label_list = self.precalculated_labels
            else:
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
        query_start, query_end, subject_start, subject_end = self.calculate_offsets(row)
        query_start_x, query_start_y, query_end_x, query_end_y = self.normalize_positions(
            query_start, query_end, 0)
        subject_start_x, subject_start_y, subject_end_x, subject_end_y = self.normalize_positions(
            subject_start, subject_end, self.comparison_height)

        match_path_desc: str = self.construct_path_description(
            query_start_x, query_start_y, query_end_x, query_end_y, subject_start_x, subject_start_y, subject_end_x, subject_end_y)
        return Path(
            d=match_path_desc,
            fill=self.match_fill_color,
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
        self.font_family: str = self.config_dict['objects']['text']['font_family']
        self.add_elements_to_group()

    def create_rectangle_path_for_legend(self) -> str:
        # Normalize start and end positions
        normalized_start: float = 0
        normalized_end: float = 16
        # Construct the rectangle path
        start_y_top: float
        start_y_bottom: float
        end_y_top: float
        end_y_bottom: float
        start_y_top, start_y_bottom = -8, 8
        end_y_top, end_y_bottom = -8, 8
        rectangle_path: str = f"M {normalized_start},{start_y_top} L {normalized_end},{end_y_top} " f"L {normalized_end},{end_y_bottom} L {normalized_start},{start_y_bottom} z"
        return rectangle_path
    def add_elements_to_group(self):
        count = 0
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        for key in self.legend_table.keys():
            rect_path = Path(
                d=path_desc,
                fill=self.legend_table[key][2],
                stroke=self.legend_table[key][0],
                stroke_width=self.legend_table[key][1])
            rect_path.translate(0, count * 25)
            self.legend_group.add(rect_path)
            legend_path = generate_text_path(key,0, 0, 0, 16, "normal", font, dominant_baseline='central', text_anchor="start")
            legend_path.translate(23, count * 25)
            self.legend_group.add(legend_path)
            count += 1
        return self.legend_group
    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the figure legends.

        Returns:
            Group: The SVG group with figure legends.
        """
        return self.legend_group