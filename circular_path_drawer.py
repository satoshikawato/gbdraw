#!/usr/bin/env python
# coding: utf-8

import logging
import sys
import math
from typing import Literal, Tuple, Dict, Union
from pandas import DataFrame
from svgwrite.shapes import Circle
from svgwrite.path import Path
from svgwrite.text import Text, TSpan, TextPath

from .create_feature_objects import set_arrow_shoulder, get_exon_and_intron_coordinates
from .utility_functions import edit_available_tracks, normalize_position_to_linear_track, calculate_bbox_dimensions

# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def generate_circle_path_desc(radius: float, norm_factor: float) -> str:
    """
    Generates the SVG path description for a circle.

    This function creates a circular path description based on the given radius and normalization factor.
    The path describes a full circle using SVG path syntax.

    Args:
        radius (float): The radius of the circle.
        norm_factor (float): The normalization factor to adjust the radius.

    Returns:
        str: A string representing the SVG path description for the circle.
    """
    norm_radius: float = radius * norm_factor
    circle_radius: float = norm_radius
    circle_start_x: float = circle_radius * \
        math.cos(math.radians(360.0 * 0 - 90))
    circle_start_y: float = circle_radius * \
        math.sin(math.radians(360.0 * 0 - 90))
    circle_desc: str = ''
    start_coordinate: str = "M {} {}".format(circle_start_x, circle_start_y)
    circle_desc += start_coordinate
    n_list: list[str] = []
    for n in range(0, 10, 1):
        deg: float = (n + 1) * 0.1
        n_x: float = circle_radius * math.cos(math.radians(360.0 * deg - 90))
        n_y: float = circle_radius * math.sin(math.radians(360.0 * deg - 90))
        n_list.append(
            "A {} {} 0 0 1 {} {}".format(
                circle_radius,
                circle_radius,
                n_x,
                n_y))
    circle_desc += ''.join(n_list)
    circle_desc += "z"
    return circle_desc


def generate_circular_gc_content_path_desc(radius: float, record_len: int, gc_df: DataFrame, track_width: float, norm_factor: float) -> str:
    """
    Generates the SVG path description for circular GC content representation.

    This function creates a path description that visually represents GC content variation along a circular genome plot.
    The path is adjusted based on the GC content data, radius, track width, and normalization factor.

    Args:
        radius (float): The radius of the circular plot.
        record_len (int): The length of the genomic sequence.
        gc_df (DataFrame): DataFrame containing GC content data.
        track_width (float): The width of the GC content track.
        norm_factor (float): The normalization factor to adjust the path.

    Returns:
        str: A string representing the SVG path description for the GC content.
    """
    norm_radius: float = radius * norm_factor
    coodinates_list: list[str] = []
    skew_start_x: float = norm_radius * \
        math.cos(math.radians(360.0 * (0 / record_len) - 90))
    skew_start_y: float = norm_radius * \
        math.sin(math.radians(360.0 * (0 / record_len) - 90))
    skew_start_position: str = "M{} {}".format(skew_start_x, skew_start_y)
    coodinates_list.append(skew_start_position)
    column: str = 'GC content'
    mean = float(gc_df[column].mean())
    max_diff = float((gc_df[column] - mean).abs().max())
    for index, row in gc_df.iterrows():
        value = float(row[column])
        diff: float = (value - mean)
        radius_of_coordinate: float = (
            norm_radius + (0.5 * track_width * (diff / max_diff)))
        x_corrdinate: float = radius_of_coordinate * \
            math.cos(math.radians(360.0 * (index / record_len) - 90)
                     )  # type: ignore
        y_corrdinate: float = radius_of_coordinate * \
            math.sin(math.radians(360.0 * (index / record_len) - 90)
                     )  # type: ignore
        corrdinate: str = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
        coodinates_list.append(corrdinate)

    gc_desc: str = "{}".format(''.join(coodinates_list))
    gc_desc += "z"
    circle_desc: str = generate_circle_path_desc(radius, norm_factor)
    gc_desc += circle_desc
    return gc_desc


def generate_circular_gc_skew_path_desc(radius: float, df: DataFrame, total_len: int, track_width: float, norm_factor: float) -> str:
    """
    Generates the SVG path description for circular GC skew representation.

    This function creates a path description that visually represents GC skew variation along a circular genome plot.
    The path is adjusted based on the GC skew data, radius, track width, and normalization factor.

    Args:
        radius (float): The radius of the circular plot.
        df (DataFrame): DataFrame containing GC skew data.
        total_len (int): The total length of the genomic sequence.
        track_width (float): The width of the GC skew track.
        norm_factor (float): The normalization factor to adjust the path.

    Returns:
        str: A string representing the SVG path description for the GC skew.
    """
    norm_radius: float = radius * norm_factor
    skew_desc_list: list[str] = []
    skew_start_x: float = norm_radius * \
        math.cos(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_y: float = norm_radius * \
        math.sin(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_position: str = "M{} {}".format(skew_start_x, skew_start_y)
    skew_desc_list.append(skew_start_position)
    column: str = 'GC skew'
    mean = float(df[column].mean())
    max_diff = float((df[column] - mean).abs().max())
    for index, row in df.iterrows():
        value = float(row[column])
        diff: float = (value - mean)
        radius_of_coordinate: float = (
            norm_radius + (0.5 * track_width * (diff / max_diff)))
        x_corrdinate: float = radius_of_coordinate * \
            math.cos(math.radians(360.0 * (index / total_len) - 90)
                     )  # type: ignore
        y_corrdinate: float = radius_of_coordinate * \
            math.sin(math.radians(360.0 * (index / total_len) - 90)
                     )  # type: ignore
        corrdinate: str = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
        skew_desc_list.append(corrdinate)
    skew_desc: str = "{}".format(''.join(skew_desc_list))
    skew_desc += "z"
    circle_desc: str = generate_circle_path_desc(radius, norm_factor)
    skew_desc += circle_desc
    return skew_desc

def calculate_cds_ratio(track_ratio, seq_length):
    if seq_length < 50000:
        cds_ratio = track_ratio * 0.50
        offset = 0.01
    else:
        cds_ratio = track_ratio * 0.25
        offset = 0.005
    return cds_ratio, offset
def calculate_feature_position_factors_circular(total_length, strand: str, track_ratio: float, track_type="tuckin", strandedness=True) -> list[float]:
    """
    Calculates position factors for a feature based on its strand orientation on a circular canvas.

    Args:
        strand (str): The strand of the feature ('positive' or 'negative').
        track_ratio (float): Ratio to determine the track width.

    Returns:
        list[float]: A list of factors used to determine the feature's position on the canvas.
    """
    CDS_RATIO, OFFSET = calculate_cds_ratio(track_ratio, total_length)
    BASE: float = 1.0
    if strandedness == True:
        if track_type == "middle":
            factors_positive: list[float] = [
                BASE, BASE + CDS_RATIO * 0.5, BASE + CDS_RATIO]
            factors_negative: list[float] = [
                BASE - CDS_RATIO, BASE - CDS_RATIO * 0.5, BASE]
        elif track_type == "spreadout":
            factors_positive: list[float] = [
                BASE + CDS_RATIO * 1.4, BASE + CDS_RATIO * 1.9, BASE + CDS_RATIO * 2.4]
            factors_negative: list[float] = [
                BASE + CDS_RATIO * 0.4, BASE + CDS_RATIO * 0.9, BASE + CDS_RATIO * 1.4]
        elif track_type == "tuckin":
            factors_positive: list[float] = [
                BASE - CDS_RATIO * 1.7, BASE - CDS_RATIO * 1.2, BASE - CDS_RATIO * 0.7]
            factors_negative: list[float] = [
                BASE - CDS_RATIO * 2.7, BASE - CDS_RATIO * 2.2, BASE - CDS_RATIO * 1.7]
        else:
            factors_positive: list[float] = [
                BASE, BASE + CDS_RATIO * 0.5, BASE + CDS_RATIO]
            factors_negative: list[float] = [
                BASE - CDS_RATIO, BASE - CDS_RATIO * 0.5, BASE]
        if strand == "positive":
            factors: list[float] = [x + OFFSET for x in factors_positive]
        else:
            factors = [x - OFFSET for x in factors_negative]
    else:
        # No strand separation: use the same three radii for both strands,
        # chosen by track_type.
        if track_type == "middle":
            base_factors = [BASE - CDS_RATIO * 0.5, BASE, BASE + CDS_RATIO * 0.5]
        elif track_type == "spreadout":
            base_factors = [
                BASE + CDS_RATIO * 0.4,
                BASE + CDS_RATIO * 0.9,
                BASE + CDS_RATIO * 1.4
            ]
        elif track_type == "tuckin":
            base_factors = [
                BASE - CDS_RATIO * 1.7,
                BASE - CDS_RATIO * 1.2,
                BASE - CDS_RATIO * 0.7
            ]
        else:
            # fallback to middle-style
            base_factors = [BASE - CDS_RATIO * 0.5, BASE, BASE + CDS_RATIO * 0.5]

        # Apply the tiny OFFSET inward/outward
        if strand == "positive":
            return [x for x in base_factors]
        else:
            return [x for x in base_factors]   
    return factors


def generate_circular_intron_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, track_ratio: float, track_type: str, strandedness: bool) -> list[str]:
    """
    Generates the SVG path description for an intron feature on a circular canvas.

    Args:
        radius (float): Radius of the circular canvas.
        coord_dict (Dict[str, Union[str, int]]): Dictionary with coordinates and feature information.
        total_length (int): Total length of the genomic sequence.
        track_ratio (float): Ratio for determining the track width.

    Returns:
        list[str]: SVG path description for the intron feature.
    """
    coord_strand: str = str(coord_dict['coord_strand'])
    coord_start: int = int(coord_dict['coord_start'])
    coord_end: int = int(coord_dict['coord_end'])
    # I admit they're magic numbers
    strand_dict: Dict[str, str] = {
        "positive": " 0 0 0 ", "negative": " 0 0 1 "}
    param: str = strand_dict[coord_strand]
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, track_type, strandedness)
    start_x_1: float = (radius * factors[1]) * math.cos(
        math.radians(360.0 * ((coord_start) / total_length) - 90))
    start_y_1: float = (radius * factors[1]) * math.sin(
        math.radians(360.0 * ((coord_start) / total_length) - 90))
    end_x_1: float = (radius * factors[1]) * \
        math.cos(math.radians(360.0 * ((coord_end) / total_length) - 90))
    end_y_1: float = (radius * factors[1]) * \
        math.sin(math.radians(360.0 * ((coord_end) / total_length) - 90))
    feature_path: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(
        # +" z"
        radius) + "," + str(radius) + param + str(end_x_1) + "," + str(end_y_1)
    return ["line", feature_path]


def generate_circular_arrowhead_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, cds_arrow_length: float, track_ratio: float, track_type: str, strandedness: bool) -> list[str]:
    """
    Generates the SVG path description for an arrowhead feature on a circular canvas.

    Args:
        radius (float): Radius of the circular canvas.
        coord_dict (Dict[str, Union[str, int]]): Dictionary with coordinates and feature information.
        total_length (int): Total length of the genomic sequence.
        cds_arrow_length (float): Length of the coding sequence arrow.
        track_ratio (float): Ratio for determining the track width.

    Returns:
        list[str]: SVG path description for the arrowhead feature.
    """
    coord_strand: str = str(coord_dict['coord_strand'])
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, track_type, strandedness)
    coord_start = int(coord_dict["coord_start"])
    coord_end = int(coord_dict["coord_end"])
    coord_len: float = abs(
        int(int(coord_dict["coord_end"]) - int(coord_dict["coord_start"])))
    coord_strand = str(coord_dict["coord_strand"])
    arrow_strand_dict: Dict[str, Tuple[int, int, str, str]] = {
        "positive": (
            coord_start,
            coord_end,
            " 0 0 1 ",
            " 0 0 0 "),
        "negative": (
            coord_end,
            coord_start,
            " 0 0 0 ",
            " 0 0 1 ")}
    arrow_start, arrow_end, param_1, param_2 = arrow_strand_dict[coord_strand]
    shoulder: float = set_arrow_shoulder(
        coord_strand, arrow_end, cds_arrow_length)
    if abs(coord_len) < cds_arrow_length:
        point_x: float = (
            radius * factors[1]) * math.cos(math.radians(360.0 * (arrow_end / total_length) - 90))
        point_y: float = (
            radius * factors[1]) * math.sin(math.radians(360.0 * (arrow_end / total_length) - 90))
        start_x_1: float = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_1: float = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_x_2: float = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_2: float = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        feature_path: str = "M " + str(start_x_1) + "," + str(start_y_1) + " L" + str(
            point_x) + "," + str(point_y) + " L" + str(start_x_2) + "," + str(start_y_2) + " z"
    else:
        point_x = (
            radius * factors[1]) * math.cos(math.radians(360.0 * (arrow_end / total_length) - 90))
        point_y = (
            radius * factors[1]) * math.sin(math.radians(360.0 * (arrow_end / total_length) - 90))
        start_x_1 = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_1 = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_x_2 = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_2 = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        end_x_1: float = (
            radius * factors[0]) * math.cos(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_y_1: float = (
            radius * factors[0]) * math.sin(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_x_2: float = (
            radius * factors[2]) * math.cos(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_y_2: float = (
            radius * factors[2]) * math.sin(math.radians(360.0 * ((shoulder) / total_length) - 90))
        feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius) + "," + str(radius) + param_1 + str(end_x_1) + "," + str(end_y_1) + " L" + str(
            point_x) + "," + str(point_y) + " L" + str(end_x_2) + "," + str(end_y_2) + "A" + str(radius) + "," + str(radius) + param_2 + str(start_x_2) + "," + str(start_y_2) + " z"
    return ["block", feature_path]


def generate_circular_rectangle_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, track_ratio: float, track_type: str, strandedness: bool) -> list[str]:
    """
    Generates the SVG path description for a rectangular feature on a circular canvas.

    Args:
        radius (float): Radius of the circular canvas.
        coord_dict (Dict[str, Union[str, int]]): Dictionary with coordinates and feature information.
        total_length (int): Total length of the genomic sequence.
        track_ratio (float): Ratio for determining the track width.

    Returns:
        list[str]: SVG path description for the rectangular feature.
    """
    coord_strand: str = str(coord_dict['coord_strand'])
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, track_type, strandedness)
    coord_start: int = int(coord_dict['coord_start'])
    coord_end: int = int(coord_dict['coord_end'])
    rect_strand_dict: dict[str, Tuple[str, str]] = {
        "positive": (
            " 0 0 1 ", " 0 0 0 "), "negative": (
            " 0 0 0 ", " 0 0 1 ")}
    param_1: str
    param_2: str
    param_1, param_2 = rect_strand_dict[coord_strand]
    start_x_1: float = (
        radius * factors[0]) * math.cos(math.radians(360.0 * (coord_start / total_length) - 90))
    start_y_1: float = (
        radius * factors[0]) * math.sin(math.radians(360.0 * (coord_start / total_length) - 90))
    start_x_2: float = (
        radius * factors[2]) * math.cos(math.radians(360.0 * (coord_start / total_length) - 90))
    start_y_2: float = (
        radius * factors[2]) * math.sin(math.radians(360.0 * (coord_start / total_length) - 90))
    end_x_1: float = (radius * factors[0]) * \
        math.cos(math.radians(360.0 * ((coord_end) / total_length) - 90))
    end_y_1: float = (radius * factors[0]) * \
        math.sin(math.radians(360.0 * ((coord_end) / total_length) - 90))
    end_x_2: float = (radius * factors[2]) * \
        math.cos(math.radians(360.0 * ((coord_end) / total_length) - 90))
    end_y_2: float = (radius * factors[2]) * \
        math.sin(math.radians(360.0 * ((coord_end) / total_length) - 90))
    feature_path: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius) + "," + str(radius) + param_1 + str(end_x_1) + "," + str(
        end_y_1) + " L" + str(end_x_2) + "," + str(end_y_2) + "A" + str(radius) + "," + str(radius) + param_2 + str(start_x_2) + "," + str(start_y_2) + " z"
    return ["block", feature_path]


def draw_circle_path(radius: float, stroke_color: str, stroke_width: float) -> Circle:
    """
    Draws a circle path for the circular canvas.

    Args:
        radius (float): Radius of the circle.
        stroke_color (str): Stroke color of the circle.
        stroke_width (float): Stroke width of the circle.

    Returns:
        Circle: An SVG circle element.
    """
    circle_path = Circle(
        center=(
            0,
            0),
        r=radius,  # type: ignore
        stroke=stroke_color,
        stroke_width=stroke_width,
        fill='none')
    return circle_path


def generate_circular_tick_paths(radius: float, total_len: int, size: str, ticks: list, tick_width: float, track_type: str) -> list[Path]:
    """
    Generates SVG path descriptions for tick marks on a circular canvas.

    Args:
        radius (float): Radius of the circular canvas.
        total_len (int): Total length of the genomic sequence.
        size (str): Size of the ticks ('small' or 'large').
        ticks (list): List of positions for the ticks.
        tick_width (float): Width of the tick marks.

    Returns:
        list[Path]: List of SVG path elements for the tick marks.
    """
    tick_paths_list: list[Path] = []
    if track_type == "middle":
        ratio: dict[str, list[float]] = {
            'small': [0.915, 0.93], 'large': [0.91, 0.93]}
    elif track_type == "spreadout":
        ratio: dict[str, list[float]] = {
            'small': [0.985, 1.0], 'large': [0.98, 1.0]}
    elif track_type == "tuckin":
        ratio: dict[str, list[float]] = {
            'small': [0.845, 0.86], 'large': [0.84, 0.86]}
    else:
        ratio: dict[str, list[float]] = {
            'small': [1.06, 1.08], 'large': [0.98, 1.0]}
    prox: float
    dist: float
    prox, dist = ratio[size]
    for tick in ticks:
        prox_x: float = (radius * prox) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        prox_y: float = (radius * prox) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        dist_x: float = (radius * dist) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        dist_y: float = (radius * dist) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        tick_path_desc: str = "M " + \
            str(prox_x) + "," + str(prox_y) + " L" + \
            str(dist_x) + "," + str(dist_y) + " z"
        tick_path = Path(
            d=tick_path_desc,
            stroke='gray',
            stroke_width=tick_width)
        tick_paths_list.append(tick_path)
    return tick_paths_list


# Draw tick labels on a circle

# generate_circular_tick_labels(radius: float, total_len: int, size: str, ticks: list, stroke: str, fill: str, font_size: float, font_weight: str, font_family: str, track_type: str) -> list[Text]:
def generate_circular_tick_labels(radius: float, total_len: int, size: str, ticks: list, stroke: str, fill: str, font_size: float, font_weight: str, font_family: str, track_type: str, dpi: int) -> list[Text]:
    tick_label_paths_list: list[Text] = []
    if track_type == "middle":
        ratio: dict[str, list[float]] = {
            'small': [0.89, 1.10], 'large': [0.89, 1.13]}
    elif track_type == "spreadout":
        ratio: dict[str, list[float]] = {
            'small': [0.96, 1.21], 'large': [0.96, 1.24]}
    elif track_type == "tuckin":
        ratio: dict[str, list[float]] = {
            'small': [0.81, 1.10], 'large': [0.81, 1.13]}
    else:
        ratio: dict[str, list[float]] = {
            'small': [0.98, 1.0], 'large': [0.98, 1.0]}
    prox: float
    dist: float
    prox, dist = ratio[size]
    for tick in ticks:
        anchor_value, baseline_value = set_tick_label_anchor_value(
            total_len, tick)
            
        #factors: list[float] = calculate_feature_position_factors_circular(
        #record_length, label["strand"], track_ratio, self.track_type)
        angle = 360.0 * (tick / total_len)
        label_text: str = str(int(tick / 1000)) + " kbp"
        
        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(label_text, font_family, font_size, dpi)
        label_as_feature_length = total_len * bbox_width_px/(2*math.pi*radius)
        label_start = tick - (label_as_feature_length/2)
        label_end = tick + (label_as_feature_length/2)
        if 0 <= angle < 90:
            param = " 0 0 1 "
            start_x_1: float = (radius * prox) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90))
            start_y_1: float = (radius * prox) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90))
            end_x: float = (
                radius * prox) * math.cos(math.radians(360.0 * (label_end / total_len) - 90))
            end_y: float = (
                radius * prox) * math.sin(math.radians(360.0 * (label_end / total_len) - 90))
        if 90 <= angle < 270:
            param = " 1 0 0 "
            start_x_1: float = (
                radius * prox) * math.cos(math.radians(360.0 * (label_end / total_len) - 90))
            start_y_1: float = (
                radius * prox) * math.sin(math.radians(360.0 * (label_end / total_len) - 90))
            end_x: float = (radius * prox) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90))
            end_y: float = (radius * prox) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90))
        elif 270 <= angle <= 360:
            param = " 0 0 1 "
            start_x_1: float = (radius * prox) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90))
            start_y_1: float = (radius * prox) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90))
            end_x: float = (
                radius * prox) * math.cos(math.radians(360.0 * ((label_end) / total_len) - 90))
            end_y: float = (
                radius * prox) * math.sin(math.radians(360.0 * ((label_end) / total_len) - 90))
        label_axis_path_desc: str = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius) + "," + str(radius) + param + str(end_x) + "," + str(end_y)
        label_axis_path = Path(
                d=label_axis_path_desc,
                stroke="none",
                fill="none")
        text_path = Text("") # The text path must go inside a text object. Parameter used here gets ignored
        text_path.add(TextPath(label_axis_path, text=label_text, startOffset="50%", method="align", text_anchor="middle", font_size=font_size, font_style='normal',font_weight='normal', font_family=font_family, dominant_baseline = "middle"))
        tick_label_paths_list.append(label_axis_path)
        tick_label_paths_list.append(text_path)
    return tick_label_paths_list


def set_tick_label_anchor_value(total_len: int, tick: float) -> tuple[Literal['middle', 'start', 'end'], Literal['text-after-edge', 'middle', 'hanging']]:
    """
    Determines the anchor and baseline values for tick labels based on their position.

    Args:
        total_len (int): Total length of the genomic sequence.
        tick (float): The position of the tick on the genomic sequence.

    Returns:
        tuple[Literal['middle', 'start', 'end'], Literal['text-after-edge', 'middle', 'hanging']]: 
        A tuple containing the text anchor and dominant baseline values for the tick label.
    """
    angle: float = (360.0 * (tick / total_len))
    if 0 <= angle < 45:
        anchor_value, baseline_value = "middle", "text-after-edge"
    elif 45 <= angle < 155:
        anchor_value, baseline_value = "start", "middle"
    elif 155 <= angle < 205:
        anchor_value, baseline_value = "middle", "hanging"
    elif 205 <= angle < 315:
        anchor_value, baseline_value = "end", "middle"
    elif 315 <= angle < 360:
        anchor_value, baseline_value = "middle", "text-after-edge"
    else:
        raise ValueError("Abnormal angle: verify the ticks and total length")
    return anchor_value, baseline_value



def generate_name_path(name_parts: list, title_x: float, title_y: float, interval: float, font_size: str, font_weight: str, font_family: str) -> Text:
    """
    Generates an SVG text element for a name or title on the circular canvas.

    Args:
        name_parts (list): List of parts of the name or title, each with potential styling.
        title_x (float): X-coordinate for the title.
        title_y (float): Y-coordinate for the title.
        interval (float): Interval for line spacing.
        font_size (str): Font size for the text.
        font_weight (str): Font weight for the text.
        font_family (str): Font family for the text.

    Returns:
        Text: An SVG text element for the name or title.
    """
    name_path = Text("", insert=(title_x, (title_y + (interval))),
                     stroke='none',
                     fill='black',
                     font_size=font_size,
                     font_weight=font_weight,
                     font_family=font_family,
                     text_anchor="middle")
    for part in name_parts:
        tspan = TSpan(
            part['text'], font_style='italic' if part['italic'] else 'normal')
        name_path.add(tspan)
    return name_path


def generate_text_path(text: str, title_x: float, title_y: float, interval: float, font_size: float, font_weight: str, font: str, dominant_baseline: str = "middle", text_anchor: str = "middle") -> Text:
    """
    Generates an SVG text element for general text on the circular canvas.

    Args:
        text (str): The text content to be displayed.
        title_x (float): X-coordinate for the text.
        title_y (float): Y-coordinate for the text.
        interval (float): Interval for line spacing.
        fontsize (str): Font size for the text.
        fontweight (str): Font weight for the text.
        font (str): Font family for the text.

    Returns:
        Text: An SVG text element for the provided text.
    """
    text_path = Text(text,
                     insert=(title_x, (title_y + (interval))),
                     stroke='none',
                     fill='black',
                     font_size=font_size,
                     font_style='normal',
                     font_weight=font_weight,
                     font_family=font,
                     dominant_baseline=dominant_baseline,
                     text_anchor=text_anchor)
    return text_path


