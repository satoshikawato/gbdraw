#!/usr/bin/env python
# coding: utf-8

import logging
import sys
import math
from typing import Optional, Union, List, Dict, Literal, Tuple
from pandas import DataFrame
from .object_configurators import GcSkewConfigurator, FeatureDrawingConfigurator
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.text import Text, TSpan, TextPath
from svgwrite.masking import ClipPath
from .feature_objects import FeatureObject
from .circular_path_drawer import get_exon_and_intron_coordinates, calculate_feature_position_factors_circular,  generate_circle_path_desc, generate_circular_gc_skew_path_desc, generate_circular_gc_content_path_desc, generate_circular_rectangle_path, generate_circular_arrowhead_path, generate_circular_intron_path, generate_name_path, generate_text_path
from .utility_functions import determine_length_parameter, calculate_bbox_dimensions, calculate_cds_ratio
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
            stroke_width=stroke_width,
            stroke_linejoin='round',
            stroke_linecap='round',
            stroke_miterlimit=4)
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

    def draw(self, feature_object: FeatureObject, group: Group, total_length: int, radius: float, track_ratio: float, track_ratio_factor, track_type: str, strandedness: bool, length_param) -> Group:
        """
        Draws genomic features based on the given FeatureObject.

        Args:
            feature_object (FeatureObject): Object containing feature data.
            group (Group): SVG group to which the features will be added.
            total_length (int): Total length of the genomic record.
            radius (float): Radius of the circular canvas.
            track_width (float): Ratio for determining the track width.

        Returns:
            Group: The updated SVG group with the features added.
        """
        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        gene_paths = FeaturePathGenerator(
            radius, total_length, track_ratio, cds_ratio, offset, track_type, strandedness).generate_circular_gene_path(feature_object)
        for gene_path in gene_paths:
            if not gene_path[0]:
                continue
            path_type: str
            path_data: str
            path_type, path_data = gene_path[0], gene_path[1]
            if path_type == "block":
                self.draw_path(path_data, group,
                               fill_color=feature_object.color)
            elif path_type == "line":
                self.draw_path(path_data, group, fill_color="none", stroke_color_specified=self.intron_stroke_color,
                               stroke_width_specified=self.intron_stroke_width)
            elif path_type == "label":
                for element in path_data:  # Iterate over SVG elements in the list
                    group.add(element)  # Add each SVG element to the group
        return group

class FeaturePathGenerator:
    """
    Generates SVG path data for genomic features on a circular canvas.

    Attributes:
        radius (float): Radius of the circular canvas.
        total_length (int): Total length of the genomic record.
        track_ratio (float): Ratio for determining the track width.
    """

    def __init__(self, radius: float, total_length: int, track_ratio: float, cds_ratio: float, offset: float, track_type, strandedness) -> None:
        """
        Initializes the FeaturePathGenerator with circular canvas settings.

        Args:
            radius (float): Radius of the circular canvas.
            total_length (int): Total length of the genomic record.
            track_width (float): Ratio for determining the track width.
        """
        self.radius: float = radius
        self.total_length: int = total_length
        self.track_ratio: float = track_ratio
        self.cds_ratio: float = cds_ratio
        self.offset: float = offset
        self.track_type = track_type
        self.strandedness = strandedness
        self.set_arrow_length()

    def set_arrow_length(self) -> None:
        """
        Sets the arrow length for feature representation based on the total length of the genomic record.
        The arrow length is calculated using a logistic function to ensure it remains within a reasonable range.
        """
        MIN_ARROW_LENGTH = 30
        MAX_ARROW_LENGTH = 700
        PARAM_A = 3
        PARAM_B = 5
        self.arrow_length: float = (MIN_ARROW_LENGTH + (MAX_ARROW_LENGTH - MIN_ARROW_LENGTH)
                                    * 1 / (1 + math.exp(- PARAM_A * (math.log10(self.total_length) - PARAM_B))))



    def generate_circular_gene_path(self, feature_object: FeatureObject):
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
                'coord_type': str(coord[0]),
                'coord_strand': coord[2],
                'coord_start': coord[3],
                'coord_end': coord[4]}
            coord_type: str = str(coord_dict['coord_type'])
            if coord_type == "line":
                coord_path: List[str] = generate_circular_intron_path(
                    self.radius, coord_dict, self.total_length, self.track_ratio, self.cds_ratio, self.offset, self.track_type, self.strandedness)
            elif coord_type == "block":
                if coord[5] and feature_object.is_directional == True:
                    coord_path = generate_circular_arrowhead_path(
                        self.radius, coord_dict, self.total_length, self.arrow_length, self.track_ratio, self.cds_ratio, self.offset,  self.track_type, self.strandedness)
                else:
                    coord_path = generate_circular_rectangle_path(
                        self.radius, coord_dict, self.total_length, self.track_ratio, self.cds_ratio, self.offset,  self.track_type, self.strandedness)
            coordinates_paths.append(coord_path)
        return coordinates_paths #, available_tracks


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
        self.font_size: str = config_dict['objects']['definition']['circular']['font_size']
        self.font_family: str = config_dict['objects']['text']['font_family']

    def draw(
        self,
        definition_group: Group,
        title_x: float,
        title_y: float,
        species_parts: list,
        strain_parts: list,
        organelle_parts: list,
        repicon_parts: list,
        gc_percent: float,
        accession: str,
        record_length: int
    ) -> Group:

        lines: list[dict] = []

        if species_parts and species_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": species_parts})

        if strain_parts and strain_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": strain_parts})
        if organelle_parts and organelle_parts[0]["text"] is not None:
                lines.append({"kind": "name", "parts": organelle_parts})
        if repicon_parts and repicon_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": repicon_parts})
        if accession.strip():
            lines.append({"kind": "plain", "text": accession})


        length_text = "{:,} bp".format(record_length)
        lines.append({"kind": "plain", "text": length_text})

        gc_text = f"{gc_percent}% GC"
        lines.append({"kind": "plain", "text": gc_text})

        n = len(lines)
        step = self.interval * 2
        start_offset = -step * (n - 2) / 2      


        for i, line in enumerate(lines):
            y_shift = start_offset + i * step

            if line["kind"] == "name":      # species / strain
                name_path = generate_name_path(
                    line["parts"],
                    title_x,
                    title_y,
                    y_shift,
                    self.font_size,
                    "bold",
                    self.font_family,
                )
                definition_group.add(name_path)

            else:                           # plain text (accession / len / GC)
                text_path = generate_text_path(
                    line["text"],
                    title_x,
                    title_y,
                    y_shift,
                    self.font_size,
                    "normal",
                    self.font_family,
                )
                definition_group.add(text_path)

        return definition_group
class LabelDrawer:
    def __init__(self, config_dict: dict) -> None:
        """
        Initializes the DefinitionDrawer with configuration settings.

        Args:
            config_dict (dict): Configuration dictionary containing style settings for the definition section.
        """
        self.config_dict = config_dict
        self.strandedness: bool = self.config_dict['canvas']['strandedness']
    def set_feature_label_anchor_value(
            self,
            total_len: int,
            tick: float,
            start_x: float,
            is_inner: bool = False
    ) -> Tuple[
            Literal["middle", "start", "end"],
            Literal["text-after-edge", "middle", "hanging"]
    ]:
        """
        Decide `text-anchor` and `dominant-baseline` for a feature label.

        Parameters
        ----------
        total_len : int
            Length of the molecule (bp).
        tick : float
            Genomic coordinate at which the label is centred.
        is_inner : bool, default ``False``
            *False*  → outer-rim labels  
            *True*   → inner-rim labels (i.e. those drawn towards the molecule’s centre).

        Returns
        -------
        anchor_value, baseline_value : str, str
        """
        angle = (360.0 * (tick / total_len)) % 360

        # Same left/right anchor logic for both inner & outer
        if   start_x > 0:
            if is_inner:
                anchor_value = "end"
            else:
                anchor_value = "start"
        else:
            if is_inner:
                anchor_value = "start"
            else:
                anchor_value = "end"

        # Baseline must flip when we draw inside the circle
        if is_inner:
            # Inside: hang text *towards* centre
            if   0   <= angle < 10:    baseline_value = "hanging"
            elif 10  <= angle < 170:   baseline_value = "middle"
            elif 170 <= angle < 190:   baseline_value = "middle"
            elif 190 <= angle < 350:   baseline_value = "middle"
            else:                      baseline_value = "hanging"
        else:
            # Outside: hang text *away* from centre (original gbdraw rule)
            if   0   <= angle < 10:    baseline_value = "text-after-edge"
            elif 10  <= angle < 170:   baseline_value = "middle"
            elif 170 <= angle < 190:   baseline_value = "hanging"
            elif 190 <= angle < 350:   baseline_value = "middle"
            else:                      baseline_value = "text-after-edge"

        return anchor_value, baseline_value
        
    def embed_label(self, group, label, radius, record_length, track_ratio):
        length_param = determine_length_parameter(record_length, self.config_dict['labels']['length_threshold']['circular'])
        track_ratio_factor = self.config_dict['canvas']['circular']['track_ratio_factors'][length_param][0]
        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        factors: list[float] = calculate_feature_position_factors_circular(
        record_length, label["strand"], track_ratio, cds_ratio, offset, self.track_type, self.strandedness)
        angle = 360.0 * (label["middle"] / record_length)
        font_px = float(str(self.font_size).rstrip("ptpx"))
        _, bbox_h = calculate_bbox_dimensions(label["label_text"], self.font_family, font_px, dpi=96)
        center_offset = (bbox_h/4) 
        
        if 0 <= angle < 90:
            param = " 0 0 1 "
            start_x_1: float = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            start_y_1: float = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            end_x: float = (
                radius * factors[1] - center_offset) * math.cos(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            end_y: float = (
                radius * factors[1] - center_offset) * math.sin(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            label_axis_path_desc: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius - center_offset) + "," + str(radius - center_offset) + param + str(end_x) + "," + str(end_y)
        if 90 <= angle < 270:
            param = " 1 0 0 "
            start_x_1: float = (
                radius * factors[1] + center_offset) * math.cos(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            start_y_1: float = (
                radius * factors[1] + center_offset) * math.sin(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            end_x: float = (radius * factors[1] + center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            end_y: float = (radius * factors[1] + center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            label_axis_path_desc: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius + center_offset) + "," + str(radius + center_offset) + param + str(end_x) + "," + str(end_y)
        elif 270 <= angle <= 360:
            param = " 0 0 1 "
            start_x_1: float = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            start_y_1: float = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90))
            end_x: float = (
                radius * factors[1] - center_offset) * math.cos(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            end_y: float = (
                radius * factors[1] - center_offset) * math.sin(math.radians(360.0 * ((label["end"]) / record_length) - 90))
            label_axis_path_desc: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius - center_offset) + "," + str(radius - center_offset) + param + str(end_x) + "," + str(end_y)
        
        label_axis_path = Path(
                d=label_axis_path_desc,
                stroke="none",
                fill="none")
        text_path = Text("") # The text path must go inside a text object. Parameter used here gets ignored
        text_path.add(TextPath(label_axis_path, text=label["label_text"], startOffset="50%", method="align", text_anchor="middle", font_size=self.font_size, font_style='normal',font_weight='normal', font_family=self.font_family, dominant_baseline = "auto"))
        group.add(label_axis_path)
        group.add(text_path)
        return group
    def add_label_on_the_rim(self, group, label, radius, record_length):
            anchor_value, baseline_value = self.set_feature_label_anchor_value(record_length, label["middle"], label["start_x"], label.get("is_inner", False))
            start_x_1 = label["start_x"]
            start_y_1 = label["start_y"]
            label_path = generate_text_path(label["label_text"], start_x_1, start_y_1, interval = 0, font_size = self.font_size, font_weight = 'normal', font = self.font_family, dominant_baseline = baseline_value, text_anchor = anchor_value)
            group.add(label_path)
            return group
    def draw(self, label, group, record_length, radius, track_ratio):
        length_threshold = self.config_dict['labels']['length_threshold']['circular']
        length_param = determine_length_parameter(record_length, length_threshold)
        self.font_size: str = self.config_dict['labels']['font_size'][length_param]
        self.font_family: str = self.config_dict['objects']['text']['font_family']
        self.track_type: str = self.config_dict['canvas']['circular']['track_type']
        self.strandedness: bool = self.config_dict['canvas']['strandedness']
        if label["is_embedded"] == True:
            group = self.embed_label(group, label, radius, record_length, track_ratio)
        else:
            group = self.add_label_on_the_rim(group, label, radius, record_length)
        return group


