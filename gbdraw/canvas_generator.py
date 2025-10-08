#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from svgwrite import Drawing
from typing import Literal
from .utility_functions import determine_length_parameter
# Logging setup
Logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


class CircularCanvasConfigurator:
    """
    Configures the settings for a circular canvas used for genomic data visualization.

    Attributes:
    output_prefix (str): Prefix for the output file.
    total_width (int): Total width of the canvas.
    total_height (int): Total height of the canvas.
    radius (float): Radius of the circular canvas.
    track_ratio (float): Ratio to determine track width.
    show_gc (bool): Flag to display GC content.
    show_skew (bool): Flag to display GC skew.
    strandedness (bool): Flag to display strandedness.

    Methods:
    calculate_dimensions(): Calculates dimensions for the canvas.
    create_svg_canvas(): Creates and returns an SVG canvas for drawing.
    get_track_ids(): Determines the track IDs for visualization.
    """
    def __init__(self, output_prefix: str, config_dict: dict, legend: str, gb_record) -> None:
        """
        Initializes the circular canvas configurator with given settings.

        Args:
        output_prefix (str): Prefix for the output file.
        config_dict (dict): Configuration dictionary with canvas settings.
        show_gc (bool, optional): Flag to display GC content. Defaults to True.
        strandedness (bool, optional): Flag to display strandedness. Defaults to True.
        show_skew (bool, optional): Flag to display GC skew. Defaults to True.
        """
        self.output_prefix: str = output_prefix
        self.config_dict = config_dict
        self.show_labels: bool = self.config_dict['canvas']['show_labels']
        if self.show_labels:
            label_setting = 'with_labels'
        else:
            label_setting = 'without_labels'
        self.default_width: int = self.config_dict['canvas']['circular']['width'][label_setting]
        self.default_height: int = self.config_dict['canvas']['circular']['height']
        self.radius: float = self.config_dict['canvas']['circular']['radius']
        self.track_ratio: float = self.config_dict['canvas']['circular']['track_ratio']
        self.show_gc: bool = self.config_dict['canvas']['show_gc']
        self.show_skew: bool = self.config_dict['canvas']['show_skew']
        self.strandedness: bool = self.config_dict['canvas']['strandedness']
        self.dpi: int = self.config_dict['png_output']['dpi']
        self.length_threshold = self.config_dict['labels']['length_threshold']['circular']
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.track_width = self.config_dict['canvas']['circular']['track_width'][self.length_param]
        self.track_ratio_factors = self.config_dict['canvas']['circular']['track_ratio_factors'][self.length_param]
        self.legend_position: str = legend

        self.calculate_dimensions()
        self.get_track_ids()
    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions and offsets for the circular canvas based on the configuration.
        """
        self.total_height = self.default_height
        self.offset_y: float = self.total_height * 0.5
        if self.legend_position == "left":
            self.total_width = self.default_width * 1.2
            self.offset_x: float = self.default_width * 0.6
        elif self.legend_position == "right":
            self.total_width = self.default_width * 1.2
            self.offset_x: float = self.default_width * 0.5
        else:
            self.total_width = self.default_width
            self.offset_x: float = self.default_width * 0.5
            
        # Create linear canvas
    def recalculate_canvas_dimensions(self, legend_config):
        if self.legend_position == "right":
            self.total_width = self.default_width + (legend_config.legend_width * 1.1)
            self.legend_offset_x = self.default_width + (legend_config.legend_width * 0.05)
            self.legend_offset_y = (self.total_height - legend_config.legend_height) / 2
        elif self.legend_position == "left":
            self.total_width = self.default_width + (legend_config.legend_width * 1.1)
            self.legend_offset_x = legend_config.legend_width * 0.05       
            self.offset_x: float = (self.default_width * 0.5) + (legend_config.legend_width * 1.1)
            self.legend_offset_y = (self.total_height - legend_config.legend_height) / 2
        elif self.legend_position == "upper_left":
            self.legend_offset_x: float =  0.025 * self.total_width
            self.legend_offset_y: float =  0.05 * self.total_height
        elif self.legend_position == "upper_right":
            self.legend_offset_x: float =  0.85 * self.total_width
            self.legend_offset_y: float =  0.05 * self.total_height
        elif self.legend_position == "lower_left":
            self.legend_offset_x: float =  0.025 * self.total_width
            self.legend_offset_y: float =  0.78 * self.total_height
        elif self.legend_position == "lower_right":
            self.legend_offset_x: float =  0.875 * self.total_width
            self.legend_offset_y: float =  0.75 * self.total_height
        elif self.legend_position == "none":
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0
        else:
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0

    def create_svg_canvas(self) -> Drawing:
        """
        Creates and returns an SVG canvas based on the configurator's settings.

        Returns:
        Drawing: An SVG Drawing object.
        """
        return Drawing(
            filename=self.output_prefix + ".svg",
            size=(str(self.total_width) + 'px', str(self.total_height) + 'px'),
            viewBox=('0 0 ' + str(self.total_width) +
                     ' ' + str(self.total_height)),
            debug=False
        )

    def get_track_ids(self) -> None:
        """
        Determines and assigns track IDs for the visualization based on the configurator's settings.
        """
        self.track_ids: dict = {}
        gc_track_id: Literal[2] | None = 2 if self.show_gc or not self.show_skew else None
        skew_track_id: Literal[3, 2] | None = (
            3 if self.show_gc else 2) if self.show_skew else None

        if gc_track_id is not None:
            self.track_ids['gc_track'] = gc_track_id
        if skew_track_id is not None:
            self.track_ids['skew_track'] = skew_track_id


class LinearCanvasConfigurator:
    """
    Configures the settings for a linear canvas used for genomic data visualization.

    Attributes:
    output_prefix (str): Prefix for the output file.
    fig_width (int): Width of the figure.
    vertical_offset (float): Vertical offset for alignment.
    horizontal_offset (float): Horizontal offset for alignment.
    vertical_padding (float): Vertical padding between elements.
    comparison_height (float): Height for comparison tracks.
    canvas_padding (float): Padding around the canvas.
    default_cds_height (float): Default height for coding sequences.
    default_gc_height (float): Default height for GC content.
    show_gc (bool): Flag to display GC content.
    strandedness (bool): Flag to display strandedness.
    num_of_entries (int): Number of entries to visualize.
    align_center (bool): Flag to align content to the center.
    longest_genome (int): Length of the longest genome in the dataset.

    Methods:
    set_gc_height_and_gc_padding(): Sets GC height and padding.
    set_cds_height_and_cds_padding(): Sets CDS height and padding.
    set_arrow_length(): Sets the arrow length for representation.
    calculate_dimensions(): Calculates dimensions for the canvas.
    create_svg_canvas(): Creates and returns an SVG canvas for drawing.
    """

    def __init__(self, num_of_entries: int, longest_genome: int, config_dict: dict, legend: str, output_prefix='out'):
        """
        Initializes the linear canvas configurator with given settings.

        Args:
        num_of_entries (int): Number of entries to visualize.
        longest_genome (int): Length of the longest genome in the dataset.
        config_dict (dict): Configuration dictionary with canvas settings.
        output_prefix (str, optional): Prefix for the output file. Defaults to 'out'.
        show_gc (bool, optional): Flag to display GC content. Defaults to False.
        strandedness (bool, optional): Flag to display strandedness. Defaults to False.
        align_center (bool, optional): Flag to align content to the center. Defaults to False.
        """
        self.output_prefix: str = output_prefix
        self.config_dict = config_dict
        self.longest_genome: int = longest_genome
        self.fig_width: int = self.config_dict['canvas']['linear']['width']
        self.original_vertical_offset: float = self.config_dict['canvas']['linear']['vertical_offset']
        self.vertical_offset: float = self.config_dict['canvas']['linear']['vertical_offset']
        self.horizontal_offset: float = self.config_dict['canvas']['linear']['horizontal_offset']
        self.vertical_padding: float = self.config_dict['canvas']['linear']['vertical_padding']
        self.comparison_height: float = self.config_dict['canvas']['linear']['comparison_height']
        self.canvas_padding: float = self.config_dict['canvas']['linear']['canvas_padding']
        self.length_threshold = self.config_dict['labels']['length_threshold']['linear']
        self.length_param = determine_length_parameter(self.longest_genome, self.length_threshold)
        self.default_cds_height: float = self.config_dict['canvas']['linear']['default_cds_height'][self.length_param]
        self.default_gc_height: float = self.config_dict['canvas']['linear']['default_gc_height']
        self.dpi: int = self.config_dict['png_output']['dpi']
        self.show_gc: bool = self.config_dict['canvas']['show_gc']
        self.show_skew: bool = self.config_dict['canvas']['show_skew']
        self.strandedness: bool = self.config_dict['canvas']['strandedness']
        self.resolve_overlaps: bool = self.config_dict['canvas']['resolve_overlaps']
        self.align_center: bool = self.config_dict['canvas']['linear']['align_center']
        self.show_labels: bool = self.config_dict['canvas']['show_labels']
        self.legend_position = legend
        self.num_of_entries: int = num_of_entries


        self.calculate_dimensions()
        self.set_arrow_length()

    def set_gc_height_and_gc_padding(self) -> None:
        """
        Sets the height and padding for the GC content track based on configuration settings.
        This method adjusts the gc_height and gc_padding attributes.
        """
        if self.show_gc:
            self.gc_height: float = self.default_gc_height
            if self.show_skew:
                self.gc_padding: float = self.gc_height
            else:
                self.gc_padding: float = self.gc_height
        else:
            self.gc_height: float = 0
            self.gc_padding: float = 0
        if self.show_skew:
            self.skew_height: float = self.default_gc_height
            if self.show_gc:
                self.skew_padding: float = self.skew_height
            else:
                self.skew_padding: float = self.skew_height 
        else:
            self.skew_height: float = 0
            self.skew_padding: float = 0

    def set_cds_height_and_cds_padding(self) -> None:
        """
        Sets the height and padding for the coding sequences (CDS) track based on configuration settings.
        This method adjusts the cds_height and cds_padding attributes.
        """
        if self.strandedness:
            self.cds_height: float = self.default_cds_height
            self.cds_padding: float = 0.5 * self.cds_height
        else:
            self.cds_height: float = 0.5 * self.default_cds_height
            self.cds_padding: float = 0.75 * self.cds_height
    def set_arrow_length(self) -> None:
        """
        Sets the length of the arrow used in the representation based on the longest genome.
        This method adjusts the arrow_length attribute.
        """
        self.arrow_length_param = self.config_dict['canvas']['linear']['arrow_length_parameter'][self.length_param]
        self.arrow_length: float = self.arrow_length_param * self.longest_genome
    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions for the linear canvas including the total width and height, 
        considering all the elements and padding. This method updates total_width and total_height attributes.
        """
        
        self.set_gc_height_and_gc_padding()
        self.set_cds_height_and_cds_padding()
        self.add_margin: float | Literal[0] = 2 * self.cds_height if (
            self.show_gc and not self.strandedness) else 0
        self.alignment_width: float = self.fig_width - self.horizontal_offset
        self.total_width = int(self.fig_width + 2 * self.canvas_padding)
        self.total_height = int(2 * self.vertical_offset + (self.cds_height + self.gc_padding) + (self.vertical_padding +
                                self.comparison_height + self.vertical_padding + self.cds_height + self.gc_padding) * (self.num_of_entries - 1))
    def recalculate_canvas_dimensions(self, legend_config):
        """
        Calculates final canvas dimensions and legend offsets, ensuring the legend fits within the canvas.
        """
        def calculate_optimal_legend_y():
            genome_area_top = self.vertical_offset
            genome_area_bottom = self.total_height - self.original_vertical_offset
            genome_area_center_y = genome_area_top + (genome_area_bottom - genome_area_top) / 2
            legend_y = genome_area_center_y - (legend_config.legend_height / 2)

            if legend_y < 0 or (legend_y + legend_config.legend_height) > self.total_height:
                legend_y = (self.total_height - legend_config.legend_height) / 2
            
            return legend_y

        if self.legend_position == "right":
            self.total_width = self.total_width + (legend_config.legend_width * 1.1)
            self.legend_offset_x = self.canvas_padding + self.fig_width + (legend_config.legend_width * 0.05)
            self.legend_offset_y = calculate_optimal_legend_y()

        elif self.legend_position == "left":
            self.total_width = self.total_width + (legend_config.legend_width * 1.1)
            self.legend_offset_x = self.canvas_padding + (legend_config.legend_width * 0.05)       
            self.horizontal_offset: float = self.canvas_padding + (legend_config.legend_width * 1.1)
            self.legend_offset_y = calculate_optimal_legend_y()

        elif self.legend_position == "none":
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0
            
        else:
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0
    def create_svg_canvas(self) -> Drawing:
        """
        Creates and returns an SVG canvas for the linear representation based on the configurator's settings.

        Returns:
        Drawing: An SVG Drawing object representing the linear canvas.
        """
        return Drawing(
            filename=self.output_prefix + ".svg",
            size=(str(self.total_width) + 'px', str(self.total_height) + 'px'),
            viewBox=('0 0 ' + str(self.total_width) +
                     ' ' + str(self.total_height)),
            debug=False
        )
