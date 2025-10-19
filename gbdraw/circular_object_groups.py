#!/usr/bin/env python
# coding: utf-8

import math
from collections import defaultdict
from typing import Optional, Literal, List, Dict
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.shapes import Circle, Line
from .canvas_generator import CircularCanvasConfigurator
from svgwrite.path import Path
from svgwrite.text import Text
from .data_processing import calculate_gc_percent, prepare_label_list
from .utility_functions import parse_mixed_content_text, determine_length_parameter, preprocess_label_filtering
from .create_feature_objects import create_feature_dict, preprocess_color_tables
from .feature_objects import FeatureObject
from .circular_feature_drawer import LabelDrawer, FeatureDrawer, SkewDrawer, GcContentDrawer, DefinitionDrawer
from .circular_path_drawer import generate_text_path, draw_circle_path, generate_circular_tick_paths, generate_circular_tick_labels
from .object_configurators import LegendDrawingConfigurator, GcContentConfigurator, GcSkewConfigurator, FeatureDrawingConfigurator


class GcContentGroup:
    """
    This class is responsible for creating a group for GC content visualization on a circular canvas.

    Attributes:
        gc_group (Group): SVG group for the GC content visualization.
        radius (float): Radius of the circular canvas.
        gc_config (GcContentConfigurator): Configuration settings for the GC content.
        gb_record (SeqRecord): GenBank record containing genomic data.
        config_dict (dict): Configuration dictionary for drawing settings.
        gc_df (DataFrame): DataFrame containing GC content data.
        track_width (float): Width of the GC content track.
        record_len (int): Length of the genomic record.
        norm_factor (float): Normalization factor for scaling GC content values.
    """

    def __init__(self, gb_record: SeqRecord, gc_df: DataFrame, radius: float, track_width: float, gc_config: GcContentConfigurator, config_dict: dict, track_id: str) -> None:
        """
        Initializes the GcContentGroup with necessary parameters and configurations.

        Args:
            gb_record (SeqRecord): GenBank record containing genomic data.
            gc_df (DataFrame): DataFrame containing GC content data.
            radius (float): Radius of the circular canvas.
            track_width (float): Width of the GC content track.
            gc_config (GcContentConfigurator): Configuration settings for the GC content.
            config_dict (dict): Configuration dictionary for drawing settings.
            track_id (str): Identifier for the track.
        """
        self.gc_group = Group(id="gc_content")
        self.radius: float = radius
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.record_len: int = len(self.gb_record.seq)
        self.config_dict: dict = config_dict
        self.gc_df: DataFrame = gc_df
        self.track_width: float = track_width
        self.length_threshold = self.config_dict['labels']['length_threshold']['circular']
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.track_type: str = self.config_dict['canvas']['circular']['track_type']
        self.norm_factor: float = self.config_dict['canvas']['circular']['track_dict'][self.length_param][self.track_type][str(track_id)]
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        """
        Adds GC content visualization elements to the group.
        """
        self.gc_group: Group = GcContentDrawer(self.gc_config).draw(
            self.radius, self.gc_group, self.gc_df, self.record_len, self.track_width, self.norm_factor, self.dinucleotide)

    def get_group(self) -> Group:
        """
        Returns the group with GC content visualization.

        Returns:
            Group: The SVG group containing GC content visualization.
        """
        return self.gc_group


class GcSkewGroup:
    """
    Represents a group for GC skew visualization on a circular genomic plot.

    This class creates a SVG group element that contains visual representations of the 
    genomic GC skew data on a circular plot, using a specified radius and track width.

    Attributes:
        gb_record (SeqRecord): The genomic record from which GC skew data is derived.
        gc_df (DataFrame): DataFrame containing calculated GC skew data.
        radius (float): The radius of the circular plot on which the data will be visualized.
        track_width (float): The width of the track on the circular plot.
        config_dict (dict): Configuration settings for drawing the GC skew.
        record_len (int): The total length of the genomic sequence.
        skew_group (Group): The SVG group element that contains the GC skew visualization.
        norm_factor (float): Normalization factor for scaling the GC skew visualization.
    """

    def __init__(self, gb_record: SeqRecord, gc_df: DataFrame, radius: float, track_width: float, skew_config: GcSkewConfigurator, config_dict: Dict, track_id: str) -> None:
        """
        Constructs the GcSkewGroup object with necessary parameters and configurations.

        Args:
            gb_record (SeqRecord): The genomic record from which GC skew data is derived.
            gc_df (DataFrame): DataFrame containing calculated GC skew data.
            radius (float): The radius of the circular plot on which the data will be visualized.
            track_width (float): The width of the track on the circular plot.
            config_dict (dict): Configuration settings for drawing the GC skew.
            track_id (str): Identifier for the specific track.
        """
        self.gb_record: SeqRecord = gb_record
        self.gc_df: DataFrame = gc_df
        self.radius: float = radius
        self.track_width: float = track_width
        self.skew_config: GcSkewConfigurator = skew_config
        self.record_len: int = len(self.gb_record.seq)
        self.skew_group = Group(id="skew")
        self.config_dict = config_dict
        self.track_type: str = self.config_dict['canvas']['circular']['track_type']
        self.length_threshold = self.config_dict['labels']['length_threshold']['circular']
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.norm_factor: float = self.config_dict['canvas']['circular']['track_dict'][self.length_param][self.track_type][str(track_id)]
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        """
        Adds the visual elements representing the GC skew data to the group.
        """
        self.skew_group: Group = SkewDrawer(self.skew_config).draw(
            self.radius, self.skew_group, self.gc_df, self.record_len, self.track_width, self.norm_factor, self.dinucleotide)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC skew visualization.

        Returns:
            Group: The SVG group with GC skew visualization elements.
        """
        return self.skew_group


class DefinitionGroup:
    """
    Responsible for creating and managing a group for displaying genomic definition information on a circular canvas.

    This class handles the creation of SVG elements that represent various definition information such as 
    species name, strain, GC content, and accession number, and organizes them within an SVG group.

    Attributes:
        gb_record (SeqRecord): The GenBank record containing the genomic data.
        canvas_config (CircularCanvasConfigurator): Configuration settings for the circular canvas.
        species (Optional[str]): The species name, if provided; otherwise, it is extracted from the GenBank record.
        strain (Optional[str]): The strain name, if provided; otherwise, it is extracted from the GenBank record.
        interval (int): The interval between lines in the definition section.
        font_size (int): The font size for the text in the definition section.
        font (str): The font family for the text in the definition section.
        config_dict (dict): Configuration dictionary for drawing settings.
        track_id (str): Identifier for the track, derived from the GenBank record ID.
        definition_group (Group): The SVG group containing the definition elements.
        radius (float): The radius of the circular canvas.
    """

    def __init__(self, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, config_dict: dict, species: Optional[str] = None, strain: Optional[str] = None) -> None:
        """
        Initializes the DefinitionGroup with the necessary data and configurations.

        Args:
            gb_record (SeqRecord): The GenBank record containing genomic data.
            canvas_config (CircularCanvasConfigurator): Configuration settings for the circular canvas.
            config_dict (dict): Configuration dictionary for drawing settings.
            species (Optional[str]): Optionally, the species name to be displayed. If not provided, it will be extracted from the gb_record.
            strain (Optional[str]): Optionally, the strain name to be displayed. If not provided, it will be extracted from the gb_record.
        """
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.species: str | None = species
        self.strain: str | None = strain
        self.replicon: str |None = None
        self.organelle: str | None = None
        self.config_dict: dict = config_dict
        self.interval: int = self.config_dict['objects']['definition']['circular']['interval']
        self.font_size: int = self.config_dict['objects']['definition']['circular']['font_size']
        self.font: str = self.config_dict['objects']['text']['font_family']
        self.track_id: str = str(self.gb_record.id).replace(" ", "_")
        self.definition_group = Group(id=self.track_id)
        self.radius: float = self.canvas_config.radius
        self.calculate_coordinates()
        self.find_organism_name()
        self.add_circular_definitions()

    def calculate_coordinates(self) -> None:
        """
        Calculates the coordinates for placing the definition texts on the circular canvas.

        This method computes the x and y coordinates for positioning the definition texts
        based on the specified radius of the circular canvas.
        """
        self.end_x_1: float = (self.radius) * \
            math.cos(math.radians(360.0 * 0 - 90))
        self.end_x_2: float = (self.radius) * \
            math.cos(math.radians(360.0 * (0.5) - 90))
        self.end_y_1: float = (self.radius) * \
            math.cos(math.radians(360.0 * (0.25) - 90))
        self.end_y_2: float = (self.radius) * \
            math.cos(math.radians(360.0 * (0.75) - 90))
        self.title_x: float = ((self.end_x_1 + self.end_x_2) / 2)
        self.title_y: float = ((self.end_y_1 + self.end_y_2) / 2)

    def find_organism_name(self) -> None:
        """
        Extracts and processes the organism name and strain information from the GenBank record.

        This method navigates through the features of the GenBank record to locate and extract 
        the organism name and strain information. It specifically looks for the 'source' feature 
        type, within which it searches for 'organism', 'isolate', and 'strain' qualifiers. 

        If the species or strain name is explicitly provided to the class during initialization, 
        this method will prioritize that information over what is extracted from the GenBank record. 
        The extracted or provided species and strain names are then further processed to handle 
        mixed content text, preparing them for visualization.

        The method sets the `species_parts` and `strain_parts` attributes of the class, which 
        are lists of dictionaries. Each dictionary in the list represents a segment of the text 
        with associated styling information, aiding in the creation of styled SVG text elements.
        """
        strain_name: str = ""
        record_name: str = ""

        for feature in self.gb_record.features:
            if feature.type == "source":
                if 'organism' in feature.qualifiers.keys():
                    record_name = feature.qualifiers['organism'][0]
                elif 'organism' in self.gb_record.annotations.keys():
                    record_name = str(self.gb_record.annotations['organism'])
                else:
                    pass
                if 'isolate' in feature.qualifiers.keys():
                    strain_name = feature.qualifiers['isolate'][0]
                elif 'strain' in feature.qualifiers.keys():
                    strain_name = feature.qualifiers['strain'][0]
                if 'chromosome' in feature.qualifiers.keys():
                    self.replicon = f"Chromosome {feature.qualifiers['chromosome'][0]}"
                elif 'plasmid' in feature.qualifiers.keys():
                    self.replicon = feature.qualifiers['plasmid'][0]
                if 'organelle' in feature.qualifiers.keys():
                    self.organelle = feature.qualifiers['organelle'][0]
                else:
                    pass
        if self.species:
            self.species_parts: List[Dict[str, str | bool |
                                          None]] = parse_mixed_content_text(self.species)
        else:
            self.species_parts = parse_mixed_content_text(record_name)
        if self.strain:
            self.strain_parts: List[Dict[str, str | bool |
                                         None]] = parse_mixed_content_text(self.strain)
        else:
            self.strain_parts = parse_mixed_content_text(strain_name)
        if self.replicon:
            self.replicon_parts: List[Dict[str, str | bool |
                                            None]] = parse_mixed_content_text(self.replicon)
        else:
            self.replicon_parts = parse_mixed_content_text("")
        if self.organelle:
            self.organelle_parts: List[Dict[str, str | bool |
                                             None]] = parse_mixed_content_text(self.organelle)
        else:
            self.organelle_parts = parse_mixed_content_text("")




    def add_circular_definitions(self) -> None:
        """
        Adds the circular definition texts to the group.

        This method utilizes the DefinitionDrawer to draw texts such as species name, strain,
        GC content, and accession number, and adds them to the definition_group.
        """
        record_length: int = len(self.gb_record.seq)
        accession: str = self.gb_record.id
        gc_percent: float = calculate_gc_percent(self.gb_record.seq)
        self.definition_group: Group = DefinitionDrawer(self.config_dict).draw(
            self.definition_group, self.title_x, self.title_y, self.species_parts, self.strain_parts, self.organelle_parts, self.replicon_parts, gc_percent, accession, record_length)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the definition section.

        Returns:
            Group: The SVG group with the definition section elements.
        """
        return self.definition_group

class LegendGroup:
    def __init__(self, canvas_config, legend_config, legend_table):
        self.legend_group = Group(id="legend")
        self.canvas_config = canvas_config
        self.legend_config = legend_config
        self.legend_table = legend_table
        self.font_family = self.legend_config.font_family
        self.font_size = self.legend_config.font_size
        self.dpi = self.canvas_config.dpi

        self.color_rect_size = self.legend_config.color_rect_size
        self.num_of_lines = len(self.legend_table.keys())
        self.add_elements_to_group()


    def create_rectangle_path_for_legend(self) -> str:
        # Normalize start and end positions
        normalized_start: float = 0
        normalized_end: float = self.color_rect_size
        # Construct the rectangle path
        start_y_top: float
        start_y_bottom: float
        end_y_top: float
        end_y_bottom: float
        start_y_top, start_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        end_y_top, end_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        rectangle_path: str = f"M {normalized_start},{start_y_top} L {normalized_end},{end_y_top} " f"L {normalized_end},{end_y_bottom} L {normalized_start},{start_y_bottom} z"
        return rectangle_path
    def add_elements_to_group(self):
        count = 0
        path_desc = f"M {0},{-0.5 * self.color_rect_size} L {self.legend_config.legend_width},{-0.5 * self.color_rect_size} " f"L {self.legend_config.legend_width},{self.legend_config.legend_height -0.5 * self.color_rect_size} L {0},{self.legend_config.legend_height -0.5 * self.color_rect_size} z"
        rect_path = Path(
                d=path_desc,
                fill="none",
                stroke="none",
                stroke_width=0)
        self.legend_group.add(rect_path)
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        line_margin = (24/14) * self.color_rect_size
        x_margin = (22/14) * self.color_rect_size
        for key, properties in self.legend_table.items():
            if properties['type'] == 'solid':
                rect_path = Path(
                    d=path_desc,
                    fill=properties['fill'],
                    stroke=properties['stroke'],
                    stroke_width=properties['width']
                )
                rect_path.translate(0, count * line_margin)
                self.legend_group.add(rect_path)
                legend_path = generate_text_path(key,0, 0, 0, self.font_size, "normal", font, dominant_baseline='central', text_anchor="start")
                legend_path.translate(x_margin, count * line_margin)
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




class TickGroup:
    """
    Represents a group for tick marks and labels on a circular genomic plot.

    This class creates a SVG group element that contains visual representations of tick marks and labels
    on a circular plot. These ticks help in indicating positions or lengths on the genomic sequence.

    Attributes:
        gb_record (SeqRecord): The genomic record from which the sequence length is derived.
        canvas_config (CircularCanvasConfigurator): Configuration settings for the circular canvas.
        radius (float): The radius of the circular plot.
        tick_group (Group): SVG group containing the tick marks and labels.
        total_len (int): The total length of the genomic sequence.
        config_dict (dict): Configuration settings for drawing the ticks.
        tick_width (float): The width of the tick marks.
        stroke (str): Stroke color for the tick labels.
        fill (str): Fill color for the tick labels.
        font_size (str): Font size for the tick labels.
        font_weight (str): Font weight for the tick labels.
        font_family (str): Font family for the tick labels.
    """

    def __init__(self, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, config_dict: dict) -> None:
        """
        Initializes the TickGroup with necessary parameters and configurations.

        Args:
            gb_record (SeqRecord): The genomic record containing the sequence.
            canvas_config (CircularCanvasConfigurator): Configuration settings for the circular canvas.
            config_dict (dict): Configuration dictionary for drawing settings.
        """
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.radius: float = self.canvas_config.radius
        self.tick_group = Group(id="tick")
        self.total_len: int = len(self.gb_record.seq)
        self.config_dict: dict = config_dict
        ticks_config = self.config_dict['objects']['ticks']
        self.tick_width: float = ticks_config['tick_width']
        self.stroke: str = ticks_config['tick_labels']['stroke']
        self.fill: str = ticks_config['tick_labels']['fill']
        self.font_size: str = ticks_config['tick_labels']['font_size']
        self.font_weight: str = ticks_config['tick_labels']['font_weight']
        self.manual_interval: Optional[int] = config_dict['objects']['scale'].get('interval')
        self.font_family: str = self.config_dict['objects']['text']['font_family']
        self.track_type: str = self.config_dict['canvas']['circular']['track_type']
        self.separate_strands: bool = self.config_dict['canvas']['strandedness']
        self.dpi = self.canvas_config.dpi
        self.set_tick_size()
        self.add_elements_to_group()

    def set_tick_size(self) -> None:
        """
        Determines the size of the tick marks based on the total length of the genomic sequence.

        This method sets the size of large and small ticks, ensuring they are appropriately 
        scaled relative to the length of the genomic sequence. Larger ticks are used for 
        longer sequences, and smaller ticks for shorter sequences.
        """
        if self.manual_interval is not None and self.manual_interval > 0:
            self.tick_large = self.manual_interval
            self.tick_small = self.manual_interval // 10
        else:
            if self.total_len <= 30000:
                tick_large, tick_small = 1000, 100
            elif 30000 < self.total_len <= 50000:
                tick_large, tick_small = 5000, 1000
            elif 50000 < self.total_len <= 150000:
                tick_large, tick_small = 10000, 1000
            elif 150000 < self.total_len <= 1000000:
                tick_large, tick_small = 50000, 10000
            elif 1000000 < self.total_len <= 10000000:
                tick_large, tick_small = 500000, 100000
            else:
                tick_large, tick_small = 1000000, 200000
            self.tick_large: Literal[10000, 50000, 200000, 1000000] = tick_large
            self.tick_small: Literal[1000, 10000, 50000, 200000] = tick_small

    def add_elements_to_group(self) -> Group:
        """
        Adds tick mark elements and labels to the group.

        This method generates paths for tick marks and labels according to the calculated 
        tick sizes and adds them to the tick group. It ensures that the ticks and labels 
        are appropriately positioned and styled according to the configuration settings.

        Returns:
            Group: The updated SVG group with tick marks and labels added.
        """
        ticks_large = list(range(0, self.total_len, self.tick_large))
        size: str = "large"
        tick_paths_large: list[Path] = generate_circular_tick_paths(
            self.radius, self.total_len, size, ticks_large, self.tick_width, self.track_type, self.separate_strands)
        ticks_large_nonzero: list[int] = [x for x in ticks_large if x != 0]
        tick_label_paths_large: list[Text] = generate_circular_tick_labels(
            self.radius, self.total_len, size, ticks_large_nonzero, self.stroke, self.fill, self.font_size, self.font_weight, self.font_family, self.track_type, self.separate_strands, self.dpi)
        for tick_path_large in tick_paths_large:
            self.tick_group.add(tick_path_large)
        for tick_label_path_large in tick_label_paths_large:
            self.tick_group.add(tick_label_path_large)
        return self.tick_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the tick marks and labels.

        Returns:
            Group: The SVG group with tick marks and labels.
        """
        return self.tick_group


class AxisGroup:
    """
    This class is responsible for creating and managing a group for the axis on a circular canvas.

    It generates a visual representation of an axis in the form of a circle and manages it within 
    an SVG group. The axis is drawn based on the specified radius and styled according to the 
    provided stroke color and width.

    Attributes:
        radius (float): The radius of the circular canvas on which the axis is drawn.
        axis_group (Group): The SVG group element containing the axis.
        stroke_color (str): The color used for the stroke of the axis circle.
        stroke_width (float): The width of the stroke for the axis circle.
    """

    def __init__(self, radius: float, config_dict: dict) -> None:
        """
        Initializes the AxisGroup with the necessary radius and configuration settings.

        Args:
            radius (float): The radius of the circular canvas on which the axis is drawn.
            config_dict (dict): A dictionary containing configuration settings, 
                                including stroke color and width for the axis.
        """
        self.radius: float = radius
        self.axis_group = Group(id="Axis")
        self.config_dict = config_dict
        self.stroke_color: str = self.config_dict['objects']['axis']['circular']['stroke_color']
        self.stroke_width: float = self.config_dict['objects']['axis']['circular']['stroke_width']
        self.add_elements_to_group()

    def draw_circular_axis(self) -> None:
        """
        Draws the circular axis based on the specified radius, stroke color, and width.

        This method creates a circle representing the axis and styles it according to the 
        class attributes. The circle is generated using the 'draw_circle_path' utility function.
        """
        self.circular_axis: Circle = draw_circle_path(
            self.radius, self.stroke_color, self.stroke_width)

    def add_elements_to_group(self) -> None:
        """
        Adds the circular axis to the axis group.

        This method is responsible for invoking the drawing of the axis and then adding it to 
        the axis group. It ensures that the axis is correctly positioned and rendered within the group.
        """
        self.draw_circular_axis()
        self.axis_group.add(self.circular_axis)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the axis.

        Returns:
            Group: The SVG group with the axis element.
        """
        return self.axis_group


class SeqRecordGroup:
    """
    Manages the creation and visualization of a SeqRecord (genomic data) on a circular canvas.

    This class is responsible for handling the representation of genomic features, such as genes or other 
    annotated elements, from a SeqRecord object onto a circular plot. It utilizes a configuration for the
    circular canvas and features to properly visualize the genomic data.

    Attributes:
        gb_record (SeqRecord): The GenBank record containing genomic data.
        canvas_config (CircularCanvasConfigurator): Configuration object for the circular canvas.
        feature_config (FeatureDrawingConfigurator): Configuration object for how features should be drawn.
        config_dict (dict): Dictionary containing additional configuration parameters.
        record_group (Group): The SVG group in which the SeqRecord visualization is contained.
    """
    def __init__(self, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, feature_config: FeatureDrawingConfigurator, config_dict: dict) -> None:
        """
        Initializes the SeqRecordGroup with the necessary data and configurations.

        Args:
            gb_record (SeqRecord): The GenBank record containing the genomic data.
            canvas_config (CircularCanvasConfigurator): Configuration settings for the circular canvas.
            feature_config (FeatureDrawingConfigurator): Configuration settings for how genomic features are drawn.
            config_dict (dict): Additional configuration parameters.
        """
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.feature_config: FeatureDrawingConfigurator = feature_config
        self.config_dict: dict = config_dict
        self.length_threshold = self.config_dict['labels']['length_threshold']['circular']
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.font_size: str = self.config_dict['objects']['ticks']['tick_labels']['font_size']
        self.font_family: str = self.config_dict['objects']['text']['font_family']
        self.show_labels = self.config_dict['canvas']['show_labels']
        self.label_stroke_width = self.config_dict['labels']['stroke_width'][self.length_param]
        self.label_stroke_color = self.config_dict['labels']['stroke_color']['label_stroke_color']
        self.label_filtering = self.config_dict['labels']['filtering']
        self.font_size = self.config_dict['objects']['features']['font_size']
        self.dpi =  self.config_dict['canvas']['dpi']
        self.track_type = self.config_dict['canvas']['circular']['track_type']
        self.strandedness = self.config_dict['canvas']['strandedness']
        self.resolve_overlaps = self.config_dict['canvas']['resolve_overlaps']
        self.track_ratio_factors = self.config_dict['canvas']['circular']['track_ratio_factors'][self.length_param]
        self.track_ratio = self.canvas_config.track_ratio
        self.record_group: Group = self.setup_record_group()
    def draw_record(self, feature_dict: Dict[str, FeatureObject], record_length: int, group: Group) -> Group:
        """
        Draws each feature from the genomic data onto the provided SVG group.

        This method iterates over the features in the feature dictionary, drawing each one 
        onto the SVG group using the FeatureDrawer class. It takes into account the length 
        of the record and the radius from the canvas configuration.

        Args:
            feature_dict (Dict[str, FeatureObject]): Dictionary of feature objects to be drawn.
            record_length (int): The length of the genomic sequence in the SeqRecord.
            group (Group): The SVG group where the features will be drawn.

        Returns:
            Group: The SVG group with the drawn features.
        """

        
        if self.show_labels == True:
            label_list = prepare_label_list(feature_dict, record_length, self.canvas_config.radius, self.track_ratio, self.config_dict) 
            for label in label_list:
                if not label["is_embedded"]:
                    group = LabelDrawer(self.config_dict).draw(label, group, record_length, self.canvas_config.radius, self.canvas_config.track_ratio)
                    if label["is_embedded"] == False:
                        line_path = Line(start=(label["middle_x"], label["middle_y"]), end=(label["start_x"], label["start_y"]),
                            stroke=self.label_stroke_color,
                            stroke_width=self.label_stroke_width)
                        group.add(line_path)
                        line_path2 = Line(start=(label["middle_x"], label["middle_y"]), end=(label["feature_middle_x"], label["feature_middle_y"]),
                            stroke=self.label_stroke_color,
                            stroke_width=self.label_stroke_width)
                        group.add(line_path2)
        for feature_object in feature_dict.values():
            group = FeatureDrawer(self.feature_config).draw(feature_object, group, record_length, self.canvas_config.radius, self.canvas_config.track_ratio, self.track_ratio_factors[0], self.track_type, self.strandedness, self.length_param)
        if self.show_labels:
            for label in label_list:
                if label["is_embedded"]:
                    group = LabelDrawer(self.config_dict).draw(label, group, record_length, self.canvas_config.radius, self.canvas_config.track_ratio)                
        return group

    def setup_record_group(self) -> Group:
        """
        Prepares the SVG group for SeqRecord visualization by setting up feature objects.

        This method creates a dictionary of feature objects based on the genomic data and 
        then draws these features onto an SVG group. It handles configuration related to the
        visualization, such as color tables and selected features.

        Returns:
            Group: The SVG group prepared for SeqRecord visualization.
        """
        selected_features_set: str = self.feature_config.selected_features_set
        color_table: Optional[DataFrame] = self.feature_config.color_table
        default_colors: Optional[DataFrame] = self.feature_config.default_colors
        label_filtering = preprocess_label_filtering(self.label_filtering)
        color_table, default_colors = preprocess_color_tables(color_table, default_colors)
        feature_dict: Dict[str, FeatureObject] = create_feature_dict(
            self.gb_record, color_table, selected_features_set, default_colors, self.strandedness, self.resolve_overlaps, label_filtering)
        track_id: str = self.gb_record.id
        record_group = Group(id=track_id)
        record_length: int = len(self.gb_record.seq)

        record_group: Group = self.draw_record(
            feature_dict, record_length, record_group)
        return record_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the SeqRecord visualization.

        Returns:
            Group: The SVG group with the SeqRecord visualization.
        """
        return self.record_group
