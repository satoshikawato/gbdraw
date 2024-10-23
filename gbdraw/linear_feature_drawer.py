#!/usr/bin/env python
# coding: utf-8

from ast import Str
import logging
import sys
from pandas import DataFrame
from svgwrite.path import Path
from svgwrite.container import Group
from typing import Optional, Literal
from .linear_path_drawer import calculate_gc_content_path_desc, create_intron_path_linear, create_arrowhead_path_linear, create_rectangle_path_linear
from .circular_path_drawer import generate_text_path
from .object_configurators import GcSkewConfigurator, FeatureDrawingConfigurator
# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


class FeatureDrawer:
    """
    Responsible for drawing genomic features on a linear canvas.

    This class handles the visualization of genomic features, such as genes or introns,
    by creating SVG path elements based on the provided genomic data.

    Attributes:
        default_feature_color (str): Default fill color for feature blocks.
        default_stroke_color (str): Default stroke color for feature blocks.
        default_stroke_width (float): Default stroke width for feature blocks.
        intron_stroke_color (str): Stroke color for intron lines.
        intron_stroke_width (float): Stroke width for intron lines.
    """

    def __init__(self, feature_config: FeatureDrawingConfigurator) -> None:
        """
        Initializes the FeatureDrawer with default drawing configurations.
        """
        self.default_feature_color: str = feature_config.block_fill_color
        self.default_stroke_color: str = feature_config.block_stroke_color
        self.default_stroke_width: float = feature_config.block_stroke_width
        self.intron_stroke_color: str = feature_config.line_stroke_color
        self.intron_stroke_width: float = feature_config.line_stroke_width

    def draw_path(self, path_data: str, group: Group, fill_color: str, stroke_color_specified: Optional[str] = None, stroke_width_specified: Optional[float] = None) -> None:
        """
        Draws a path on the provided SVG group using the specified style parameters.

        Args:
            path_data (str): SVG path data string.
            group (Group): SVG group to which the path will be added.
            fill_color (str): Fill color for the path.
            stroke_color (Optional[str]): Stroke color for the path. Uses default if None.
            stroke_width (Optional[float]): Stroke width for the path. Uses default if None.
        """
        stroke_color: str = stroke_color_specified if stroke_color_specified is not None else self.default_stroke_color
        stroke_width: float = stroke_width_specified if stroke_width_specified is not None else self.default_stroke_width
        path = Path(
            d=path_data,
            fill=fill_color,
            stroke=stroke_color,
            stroke_width=stroke_width)
        group.add(path)
    def draw_label(self, path_data, group: Group) -> None:
        """
        Draws a path on the provided SVG group using the specified style parameters.

        Args:
            path_data (str): SVG path data string.
            group (Group): SVG group to which the path will be added.
            fill_color (str): Fill color for the path.
            stroke_color (Optional[str]): Stroke color for the path. Uses default if None.
            stroke_width (Optional[float]): Stroke width for the path. Uses default if None.
        """
        group.add(path_data)
    def draw(self, feature_object, group: Group, genome_length: int, cds_height: float, alignment_width: float, normalization_factor: float, strandedness: str, arrow_length: float) -> Group:
        """
        Draws the provided feature object onto the group based on the specified configuration.

        Args:
            feature_object: The feature object containing the genomic feature data.
            group (Group): SVG group where the features will be drawn.
            genome_length (int): The length of the genomic sequence.
            cds_height (float): Height of the coding sequence tracks.
            alignment_width (float): Width of the alignment area.
            normalization_factor (float): Factor for normalizing the size of features.
            strandedness (str): The strandedness ('positive' or 'negative') of the feature.
            arrow_length (float): Length of the arrow for directional features.

        Returns:
            Group: The SVG group with the drawn features.
        """
        gene_paths = FeaturePathGenerator(genome_length, alignment_width, cds_height, normalization_factor, strandedness, arrow_length).generate_linear_gene_path(feature_object)
        for gene_path in gene_paths:
            path_type: str
            path_data: str
            path_type, path_data = gene_path[0], gene_path[1]
            if path_type == "block":
                self.draw_path(path_data, group, fill_color=feature_object.color)
            elif path_type == "line":
                self.draw_path(path_data, group, fill_color="none",
                               stroke_color_specified=self.intron_stroke_color, stroke_width_specified=self.intron_stroke_width)
        return group


class FeaturePathGenerator:
    """
    Generates SVG path data for genomic features on a linear canvas.

    This class is tasked with creating path data for different types of genomic features
    such as genes, introns, and arrowheads based on their coordinates and properties.

    Attributes:
        genome_length (int): Total length of the genomic sequence.
        alignment_width (float): Width of the alignment area.
        cds_height (float): Height of the coding sequence tracks.
        genome_size_normalization_factor (float): Factor for normalizing feature sizes.
        strandedness (str): Strand orientation of the features ('positive' or 'negative').
        arrow_length (float): Length of the arrow for directional features.
    """

    def __init__(self, genome_length: int, alignment_width: float, cds_height: float, genome_size_normalization_factor: float, strandedness: str, arrow_length: float) -> None:
        """
        Initializes the FeaturePathGenerator with configuration settings for the linear canvas.

        Args:
            genome_length (int): Total length of the genomic sequence.
            alignment_width (float): Width of the alignment area.
            cds_height (float): Height of the coding sequence tracks.
            genome_size_normalization_factor (float): Factor for normalizing feature sizes.
            strandedness (str): Strand orientation of the features.
            arrow_length (float): Length of the arrow for directional features.
        """
        self.genome_length: int = genome_length
        self.alignment_width: float = alignment_width
        self.cds_height: float = cds_height
        self.genome_size_normalization_factor: float = genome_size_normalization_factor
        self.strandedness: str = strandedness
        self.arrow_length: float = arrow_length
     
    def generate_linear_gene_path(self, gene_object):
        """
        Generates SVG path data for a given gene object.

        Args:
            gene_object: The gene object containing genomic feature data.

        Returns:
            A list of path data and types for the genomic features in the gene object.
        """

        coords = gene_object.location
        coordinates_paths = []
        for coord in coords:
            coord_dict = {
                'feat_type': coord[0],
                'feat_strand': coord[2],
                'feat_start': coord[3],
                'feat_end': coord[4]}
            feat_type: str = coord_dict['feat_type']
            if feat_type == "line":
                coord_path: list[str] = create_intron_path_linear(
                    coord_dict, self.genome_length, self.alignment_width, self.genome_size_normalization_factor, self.cds_height, self.strandedness)
            elif feat_type == "block":
                if coord[5] and gene_object.is_directional == True:
                    coord_path = create_arrowhead_path_linear(
                        coord_dict, self.arrow_length, self.cds_height, self.strandedness, self.genome_length, self.alignment_width, self.genome_size_normalization_factor)
                else:
                    coord_path = create_rectangle_path_linear(
                        coord_dict, self.genome_length, self.alignment_width, self.genome_size_normalization_factor, self.cds_height, self.strandedness)
            else:
                continue
            coordinates_paths.append(coord_path)

        return coordinates_paths

class LabelDrawer:
    def __init__(self, config_dict: dict) -> None:
        """
        Initializes the DefinitionDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing style settings for the definition section.
        """
        self.config_dict = config_dict

        
    def add_label(self, label_entry, group):
        feature_label_text = label_entry["label_text"]
        middle_x = label_entry["middle_x"]
        middle_y = label_entry["middle_y"]
        label_path = generate_text_path(text = feature_label_text, title_x = middle_x, title_y = middle_y, interval = 0, font_size = label_entry["font_size"], font_weight = 'normal', font = label_entry["font_family"], dominant_baseline = "central", text_anchor = "middle")
        group.add(label_path)
        #feature_path: str = "M " + \
        #str(start_x) + "," + str(start_y) + "L" + \
        #str(end_x) + "," + str(end_y) + " z"
        #group.add(label_path)            
        return group
    def draw(self, label_entry, group):
        group = self.add_label(label_entry, group)
        return group


class GcContentDrawer:
    """
    Handles the drawing of GC content visualization on a linear canvas.

    This class is a static utility for creating SVG paths that represent GC content along a linear genomic plot.

    """
    def __init__(self, gc_config):
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity

    def draw(self, group, gc_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float, start_x: float, start_y: float) -> Group:
        """
        Draws the GC content path on the provided SVG group.

        Args:
            group (Group): SVG group to which the GC content path will be added.
            gc_df (DataFrame): DataFrame containing GC content data.
            record_len (int): Length of the genomic record.
            alignment_width (float): Width of the alignment area.
            genome_size_normalization_factor (float): Normalization factor for scaling the GC content.
            track_height (float): Height of the GC content track.
            gc_path_fill_color (str): Fill color for the GC content path.
            gc_path_stroke_color (str): Stroke color for the GC content path.
            gc_path_fill_opacity (float): Fill opacity for the GC content path.
            start_x (float): Starting x-coordinate for the path.
            start_y (float): Starting y-coordinate for the path.

        Returns:
            Group: The updated SVG group with the GC content path added.
        """
        gc_path_desc: str = calculate_gc_content_path_desc(
            start_x, start_y, gc_df, record_len, alignment_width, genome_size_normalization_factor, track_height)
        gc_path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd")
        group.add(gc_path)
        return group

class SkewDrawer:
    """
    This class is responsible for drawing the GC skew on a circular canvas.

    Attributes:
        skew_high_fill_color (str): Fill color for high GC skew.
        skew_low_fill_color (str): Fill color for low GC skew.
        skew_stroke_color (str): Stroke color for the GC skew path.
        skew_fill_opacity (float): Fill opacity for the GC skew path.
    """

    def __init__(self, skew_config: GcSkewConfigurator) -> None:
        """
        Initializes the SkewDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing color and opacity settings.
        """
        self.skew_high_fill_color: str = skew_config.high_fill_color
        self.skew_low_fill_color: str = skew_config.low_fill_color
        self.skew_stroke_color: str = skew_config.stroke_color
        self.skew_fill_opacity: float = skew_config.fill_opacity

    def draw(self, group, gc_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float, start_x: float, start_y: float) -> Group:
        """
        Draws the GC skew path on the provided SVG group.

        Args:
            radius (float): Radius of the circular canvas.
            group (Group): SVG group to which the skew path will be added.
            gc_df (DataFrame): DataFrame containing GC skew data.
            record_len (int): Length of the genomic record.
            track_width (float): Width of the track on the canvas.
            norm_factor (float): Normalization factor for scaling the skew values.

        Returns:
            Group: The updated SVG group with the GC skew path added.
        """
        skew_desc: str = generate_circular_gc_skew_path_desc(
            radius, gc_df, record_len, track_width, norm_factor)
        circle_desc: str = generate_circle_path_desc(radius, norm_factor)
        circle_path: ClipPath = ClipPath(id='clipper_circle')
        circle_path.add(Path(d=circle_desc, fill="white", stroke='none'))
        skew_high: Path = Path(
            d=skew_desc,
            fill=self.skew_high_fill_color,
            stroke=self.skew_stroke_color,
            fill_opacity=self.skew_fill_opacity,
            fill_rule="evenodd")
        skew_low: Path = Path(
            d=skew_desc,
            fill=self.skew_low_fill_color,
            stroke=self.skew_stroke_color,
            fill_opacity=self.skew_fill_opacity,
            clip_path="url(#clipper_circle)",
            clip_rule="nonzero",
            fill_rule="evenodd")
        group.add(circle_path)
        group.add(skew_high)
        group.add(skew_low)


    
        gc_path_desc: str = calculate_gc_content_path_desc(
            start_x, start_y, gc_df, record_len, alignment_width, genome_size_normalization_factor, track_height)
        gc_path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd")
        group.add(gc_path)
        return group
