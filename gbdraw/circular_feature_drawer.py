#!/usr/bin/env python
# coding: utf-8

import logging
import sys
import math
from typing import Optional, Union, List, Dict
from pandas import DataFrame
from .object_configurators import GcSkewConfigurator, FeatureDrawingConfigurator
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.text import Text
from svgwrite.masking import ClipPath
from .feature_objects import FeatureObject
from .circular_path_drawer import generate_circle_path_desc, generate_circular_gc_skew_path_desc, generate_circular_gc_content_path_desc, generate_circular_rectangle_path, generate_circular_arrowhead_path, generate_circular_intron_path, generate_name_path, generate_text_path


# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


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

    def draw(self, radius: float, group: Group, gc_df: DataFrame, record_len: int, track_width: float, norm_factor: float) -> Group:
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
        return group


class GcContentDrawer:
    """
    This class is responsible for drawing the GC content on a circular canvas.

    Attributes:
        gc_path_fill_color (str): Fill color for the GC content path.
        gc_path_stroke_color (str): Stroke color for the GC content path.
        gc_path_fill_opacity (float): Fill opacity for the GC content path.
    """

    def __init__(self, gc_config: dict) -> None:
        """
        Initializes the GcContentDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing color and opacity settings.
        """
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity
    def draw(self, radius: float, group: Group, gc_df: DataFrame, record_len: int, track_width: float, norm_factor: float) -> Group:
        """
        Draws the GC content path on the provided SVG group.

        Args:
            radius (float): Radius of the circular canvas.
            group (Group): SVG group to which the GC content path will be added.
            gc_df (DataFrame): DataFrame containing GC content data.
            record_len (int): Length of the genomic record.
            track_width (float): Width of the track on the canvas.
            norm_factor (float): Normalization factor for scaling the GC content values.

        Returns:
            Group: The updated SVG group with the GC content path added.
        """
        gc_path_desc: str = generate_circular_gc_content_path_desc(
            radius, record_len, gc_df, track_width, norm_factor)
        gc_path: Path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd")
        group.add(gc_path)
        return group


class FeatureDrawer:
    """
    This class is responsible for drawing genomic features on a circular canvas.

    Attributes:
        default_feature_color (str): Default fill color for feature blocks.
        default_stroke_color (str): Default stroke color for feature blocks.
        default_stroke_width (float): Default stroke width for feature blocks.
        intron_stroke_color (str): Stroke color for intron lines.
        intron_stroke_width (float): Stroke width for intron lines.
    """

    def __init__(self, feature_config: FeatureDrawingConfigurator) -> None:
        """
        Initializes the FeatureDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing color and stroke settings.
        """
        self.default_feature_color: str = feature_config.block_fill_color
        self.default_stroke_color: str = feature_config.block_stroke_color
        self.default_stroke_width: float = feature_config.block_stroke_width
        self.intron_stroke_color: str = feature_config.line_stroke_color
        self.intron_stroke_width: float = feature_config.line_stroke_width

    def draw_path(self, path_data: str, group: Group, fill_color: str, stroke_color_specified: Optional[str] = None, stroke_width_specified: Optional[float] = None) -> None:
        """
        Draws a path with specified attributes on the provided SVG group.

        Args:
            path_data (str): SVG path data string.
            group (Group): SVG group to which the path will be added.
            fill_color (str): Fill color for the path.
            stroke_color_specified (Optional[str]): Optional stroke color for the path.
            stroke_width_specified (Optional[float]): Optional stroke width for the path.
        """
        stroke_color: str = stroke_color_specified if stroke_color_specified is not None else self.default_stroke_color
        stroke_width: float = stroke_width_specified if stroke_width_specified is not None else self.default_stroke_width
        path: Path = Path(
            d=path_data,
            fill=fill_color,
            stroke=stroke_color,
            stroke_width=stroke_width)
        group.add(path)

    def draw(self, feature_object: FeatureObject, group: Group, total_length: int, radius: float, track_ratio: float) -> Group:
        """
        Draws genomic features based on the given FeatureObject.

        Args:
            feature_object (FeatureObject): Object containing feature data.
            group (Group): SVG group to which the features will be added.
            total_length (int): Total length of the genomic record.
            radius (float): Radius of the circular canvas.
            track_ratio (float): Ratio for determining the track width.

        Returns:
            Group: The updated SVG group with the features added.
        """
        gene_paths: list[list[str]] = FeaturePathGenerator(
            radius, total_length, track_ratio).generate_circular_gene_path(feature_object)
        for gene_path in gene_paths:
            path_type: str
            path_data: str
            path_type, path_data = gene_path
            if path_type == "block":
                self.draw_path(path_data, group,
                               fill_color=feature_object.color)

            elif path_type == "line":
                self.draw_path(path_data, group, fill_color="none", stroke_color_specified=self.intron_stroke_color,
                               stroke_width_specified=self.intron_stroke_width)
        return group


class FeaturePathGenerator:
    """
    Generates SVG path data for genomic features on a circular canvas.

    Attributes:
        radius (float): Radius of the circular canvas.
        total_length (int): Total length of the genomic record.
        track_ratio (float): Ratio for determining the track width.
    """

    def __init__(self, radius: float, total_length: int, track_ratio: float) -> None:
        """
        Initializes the FeaturePathGenerator with circular canvas settings.

        Args:
            radius (float): Radius of the circular canvas.
            total_length (int): Total length of the genomic record.
            track_ratio (float): Ratio for determining the track width.
        """
        self.radius: float = radius
        self.total_length: int = total_length
        self.track_ratio: float = track_ratio
        self.set_arrow_length()

    def set_arrow_length(self) -> None:
        """
        Sets the arrow length for feature representation based on the total length of the genomic record.
        The arrow length is calculated using a logistic function to ensure it remains within a reasonable range.
        """
        MIN_ARROW_LENGTH = 30
        MAX_ARROW_LENGTH = 500
        PARAM_A = 3
        PARAM_B = 5
        self.arrow_length: float = (MIN_ARROW_LENGTH + (MAX_ARROW_LENGTH - MIN_ARROW_LENGTH)
                                    * 1 / (1 + math.exp(- PARAM_A * (math.log10(self.total_length) - PARAM_B))))

    def calculate_feature_position_factors_circular(self, strand: str) -> list[float]:
        """
        Calculates position factors for a feature based on its strand orientation on a circular canvas.

        Args:
            strand (str): The strand of the feature ('positive' or 'negative').

        Returns:
            list[float]: A list of factors used to determine the feature's position on the canvas.
        """
        OFFSET = 0.005
        CDS_RATIO: float = self.track_ratio * 0.25
        BASE = 1.0
        factors_positive: list[float] = [
            BASE, BASE + CDS_RATIO * 0.5, BASE + CDS_RATIO]
        factors_negatie: list[float] = [
            BASE - CDS_RATIO, BASE - CDS_RATIO * 0.5, BASE]
        if strand == "positive":
            factors: list[float] = [x + OFFSET for x in factors_positive]
        else:
            factors = [x - OFFSET for x in factors_negatie]
        return factors

    def generate_circular_gene_path(self, feature_object: FeatureObject) -> list[list[str]]:
        """
        Generates SVG path data for each feature in a feature object on a circular canvas.

        Args:
            feature_object (FeatureObject): The feature object containing genomic feature data.

        Returns:
            list[list[str]]: A list containing SVG path data for each feature in the feature object.
        """
        coords: list[List[Union[str, int, bool]]] = feature_object.location
        coordinates_paths: List[List[str]] = []
        for coord in coords:
            coord_path = []
            coord_dict: Dict[str, Union[str, int]] = {
                'feat_type': str(coord[0]),
                'feat_strand': coord[2],
                'feat_start': coord[3],
                'feat_end': coord[4]}
            feat_type: str = str(coord_dict['feat_type'])
            if feat_type == "line":
                coord_path: List[str] = generate_circular_intron_path(
                    self.radius, coord_dict, self.total_length, self.track_ratio)
            elif feat_type == "block":
                if coord[5] and feature_object.is_directional == True:
                    coord_path = generate_circular_arrowhead_path(
                        self.radius, coord_dict, self.total_length, self.arrow_length, self.track_ratio)
                else:
                    coord_path = generate_circular_rectangle_path(
                        self.radius, coord_dict, self.total_length, self.track_ratio)
            coordinates_paths.append(coord_path)
        return coordinates_paths


class DefinitionDrawer:
    """
    This class is responsible for drawing the definition section (title, strain, etc.) on a circular canvas.

    Attributes:
        interval (float): Interval between lines in the definition section.
        fontsize (str): Font size for the text in the definition section.
        font (str): Font family for the text in the definition section.
    """

    def __init__(self, config_dict: dict) -> None:
        """
        Initializes the DefinitionDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing style settings for the definition section.
        """
        self.interval: float = config_dict['objects']['definition']['circular']['interval']
        self.fontsize: str = config_dict['objects']['definition']['circular']['fontsize']
        self.font: str = config_dict['objects']['text']['font_family']

    def draw(self, definition_group: Group, title_x: float, title_y: float, species_parts: list, strain_parts: list, gc_percent: float, accession: str, record_length: int) -> Group:
        """
        Draws the definition section on the provided SVG group.

        Args:
            definition_group (Group): SVG group for the definition section.
            title_x (float): X-coordinate for the title.
            title_y (float): Y-coordinate for the title.
            species_parts (list): List of species name parts.
            strain_parts (list): List of strain name parts.
            gc_percent (float): GC content percentage.
            accession (str): Accession number.
            record_length (int): Length of the genomic record.

        Returns:
            Group: The updated SVG group with the definition section added.
        """
        self.fontsize: str = '{}pt'.format(self.fontsize)
        species_path: Text = generate_name_path(
            species_parts, title_x, title_y, self.interval*-2, self.fontsize, "bold", self.font)
        definition_group.add(species_path)
        strain_path: Text = generate_name_path(
            strain_parts, title_x, title_y, self.interval*-0, self.fontsize, "bold", self.font)
        definition_group.add(strain_path)
        accession_path: Text = generate_text_path(
            accession, title_x, title_y, self.interval*2, self.fontsize, "normal", self.font)
        definition_group.add(accession_path)
        length_text: str = "{:,} bp".format(record_length)
        length_path: Text = generate_text_path(
            length_text, title_x, title_y, self.interval*4, self.fontsize, "normal", self.font)
        definition_group.add(length_path)
        gc_percent_text: str = str(gc_percent) + "% GC"
        gc_percent_path: Text = generate_text_path(
            gc_percent_text, title_x, title_y, self.interval*6, self.fontsize, "normal", self.font)
        definition_group.add(gc_percent_path)
        return definition_group
