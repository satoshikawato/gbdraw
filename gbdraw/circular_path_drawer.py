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
    if seq_length < 25000:
        cds_ratio = track_ratio * 0.50
    elif 25000 <= seq_length < 50000:
        cds_ratio = track_ratio * 0.35
    else:
        cds_ratio = track_ratio * 0.25
    return cds_ratio
def calculate_feature_position_factors_circular(total_length, strand: str, track_ratio: float, track_type="tuckin") -> list[float]:
    """
    Calculates position factors for a feature based on its strand orientation on a circular canvas.

    Args:
        strand (str): The strand of the feature ('positive' or 'negative').
        track_ratio (float): Ratio to determine the track width.

    Returns:
        list[float]: A list of factors used to determine the feature's position on the canvas.
    """
    OFFSET: float = 0.005
    CDS_RATIO = calculate_cds_ratio(track_ratio, total_length)
    BASE: float = 1.0
    if track_type == "default":
        factors_positive: list[float] = [
            BASE, BASE + CDS_RATIO * 0.5, BASE + CDS_RATIO]
        factors_negatie: list[float] = [
            BASE - CDS_RATIO, BASE - CDS_RATIO * 0.5, BASE]
    elif track_type == "spreadout":
        factors_positive: list[float] = [
            BASE + CDS_RATIO * 1.4, BASE + CDS_RATIO * 1.9, BASE + CDS_RATIO * 2.4]
        factors_negatie: list[float] = [
            BASE + CDS_RATIO * 0.4, BASE + CDS_RATIO * 0.9, BASE + CDS_RATIO * 1.4]
    elif track_type == "tuckin":
        factors_positive: list[float] = [
            BASE - CDS_RATIO * 1.7, BASE - CDS_RATIO * 1.2, BASE - CDS_RATIO * 0.7]
        factors_negatie: list[float] = [
            BASE - CDS_RATIO * 2.7, BASE - CDS_RATIO * 2.2, BASE - CDS_RATIO * 1.7]
    else:
        factors_positive: list[float] = [
            BASE, BASE + CDS_RATIO * 0.5, BASE + CDS_RATIO]
        factors_negatie: list[float] = [
            BASE - CDS_RATIO, BASE - CDS_RATIO * 0.5, BASE]
    if strand == "positive":
        factors: list[float] = [x + OFFSET for x in factors_positive]
    else:
        factors = [x - OFFSET for x in factors_negatie]
    return factors


def generate_circular_intron_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, track_ratio: float) -> list[str]:
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
        total_length, coord_strand, track_ratio)
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


def generate_circular_arrowhead_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, cds_arrow_length: float, track_ratio: float):
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
        total_length, coord_strand, track_ratio)
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


def generate_circular_rectangle_path(radius: float, coord_dict: Dict[str, Union[str, int]], total_length: int, track_ratio: float) -> list[str]:
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
        total_length, coord_strand, track_ratio)
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


def generate_circular_tick_paths(radius: float, total_len: int, size: str, ticks: list, tick_width: float, track_type = "tuckin") -> list[Path]:
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
    if track_type == "default":
        ratio: dict[str, list[float]] = {
            'small': [1.06, 1.08], 'large': [1.06, 1.11]}
    elif track_type == "spreadout":
        ratio: dict[str, list[float]] = {
            'small': [1.13, 1.15], 'large': [1.13, 1.18]}
    elif track_type == "tuckin":
        ratio: dict[str, list[float]] = {
            'small': [0.985, 1.00], 'large': [0.98, 1.00]}
    else:
        ratio: dict[str, list[float]] = {
            'small': [1.06, 1.08], 'large': [1.06, 1.11]}
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


def generate_circular_tick_labels(radius: float, total_len: int, size: str, ticks: list, stroke: str, fill: str, font_size: str, font_weight: str, font_family: str, track_type = "tuckin") -> list[Text]:
    """
    Generates SVG text elements for tick labels on a circular canvas.

    Args:
        radius (float): Radius of the circular canvas.
        total_len (int): Total length of the genomic sequence.
        size (str): Size of the ticks ('small' or 'large').
        ticks (list): List of positions for the tick labels.
        stroke (str): Stroke color for the text.
        fill (str): Fill color for the text.
        font_size (str): Font size of the text.
        font_weight (str): Font weight of the text.
        font_family (str): Font family of the text.

    Returns:
        list[Text]: List of SVG text elements for the tick labels.
    """
    tick_label_paths_list: list[Text] = []
    if track_type == "default":
        ratio: dict[str, list[float]] = {
            'small': [1.10, 1.11], 'large': [1.12, 1.14]}
    elif track_type == "spreadout":
        ratio: dict[str, list[float]] = {
            'small': [1.19, 1.21], 'large': [1.19, 1.24]}
    elif track_type == "tuckin":
        ratio: dict[str, list[float]] = {
            'small': [1.04, 1.06], 'large': [1.04, 1.09]}
    else:
        ratio: dict[str, list[float]] = {
            'small': [1.10, 1.11], 'large': [1.12, 1.14]}
    prox: float
    dist: float
    prox, dist = ratio[size]
    for tick in ticks:
        anchor_value, baseline_value = set_tick_label_anchor_value(
            total_len, tick)
        label_text: str = str(int(tick / 1000)) + " kbp"
        tixk_label_x: float = (radius * prox) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        tick_label_y: float = (radius * prox) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        label_path = Text(
            label_text,
            insert=(
                tixk_label_x,
                tick_label_y),
            stroke=stroke,
            fill=fill,
            font_size=font_size,
            font_weight=font_weight,
            font_family=font_family,
            text_anchor=anchor_value,
            dominant_baseline=baseline_value)
        tick_label_paths_list.append(label_path)
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



def generate_name_path(name_parts: list, title_x: float, title_y: float, interval: float, fontsize: str, fontweight: str, font: str) -> Text:
    """
    Generates an SVG text element for a name or title on the circular canvas.

    Args:
        name_parts (list): List of parts of the name or title, each with potential styling.
        title_x (float): X-coordinate for the title.
        title_y (float): Y-coordinate for the title.
        interval (float): Interval for line spacing.
        fontsize (str): Font size for the text.
        fontweight (str): Font weight for the text.
        font (str): Font family for the text.

    Returns:
        Text: An SVG text element for the name or title.
    """
    name_path = Text("", insert=(title_x, (title_y + (interval))),
                     stroke='none',
                     fill='black',
                     font_size=fontsize,
                     font_weight=fontweight,
                     font_family=font,
                     text_anchor="middle")
    for part in name_parts:
        tspan = TSpan(
            part['text'], font_style='italic' if part['italic'] else 'normal')
        name_path.add(tspan)
    return name_path


def generate_text_path(text: str, title_x: float, title_y: float, interval: float, fontsize: str, fontweight: str, font: str, dominant_baseline: str = "auto", text_anchor: str = "middle") -> Text:
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
                     font_size=fontsize,
                     font_style='normal',
                     font_weight=fontweight,
                     font_family=font,
                     dominant_baseline=dominant_baseline,
                     text_anchor=text_anchor)
    return text_path

def generate_circular_label_path(radius: float, total_length: int, track_ratio: float, feature_object, available_tracks):
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
    coordinates = feature_object.coordinates
    location: list[Tuple[str, str, str, int, int, bool]
                   ] = get_exon_and_intron_coordinates(coordinates, total_length)
    feat_strand: str = str(feature_object.coordinates)
    strand_dict = {
        "positive": (
            " 0 0 1 ", " 0 0 0 "), "negative": (
            " 0 0 0 ", " 0 0 1 ")}
    param_1, param_2 = strand_dict[feat_strand]
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, feat_strand, track_ratio)
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])
    feat_middle = (feat_start + feat_end)/2
    feat_len: float = abs(
        int(int(coord_dict["feat_end"]) - int(coord_dict["feat_start"])))
    feat_strand = str(coord_dict["feat_strand"])
    feature_lengh_in_pixels = (2*math.pi*radius) * (feat_len)/total_length
    feature_label_text = feature_object.label_text
    if feature_label_text == '':
        return [], available_tracks
    else:       
        bbox_width_px, _ = calculate_bbox_dimensions(feature_label_text, "Liberation Sans", 8, 96)
    bbox_as_feature_length = total_length * bbox_width_px/(2*math.pi*radius)
    bbox_start = feat_middle - (bbox_as_feature_length/2)
    bbox_end = feat_middle + (bbox_as_feature_length/2)
    if bbox_width_px <= abs(feature_lengh_in_pixels):
        if (0.25*total_length) < feat_middle < (0.75 *total_length):
            end_x: float = (radius * factors[1]) * math.cos(
                math.radians(360.0 * (bbox_start / total_length) - 90))
            end_y: float = (radius * factors[1]) * math.sin(
                math.radians(360.0 * (bbox_start / total_length) - 90))
            start_x: float = (
                radius * factors[1]) * math.cos(math.radians(360.0 * ((bbox_end) / total_length) - 90))
            start_y: float = (
                radius * factors[1]) * math.sin(math.radians(360.0 * ((bbox_end) / total_length) - 90))  
        else:
            start_x: float = (radius * factors[1]) * math.cos(
                math.radians(360.0 * (bbox_start / total_length) - 90))
            start_y: float = (radius * factors[1]) * math.sin(
                math.radians(360.0 * (bbox_start / total_length) - 90))
            end_x: float = (
                radius * factors[1]) * math.cos(math.radians(360.0 * ((bbox_end) / total_length) - 90))
            end_y: float = (
                radius * factors[1]) * math.sin(math.radians(360.0 * ((bbox_end) / total_length) - 90))
        label_axis_path_desc: str = "M " + str(start_x) + "," + str(start_y) + "A" + str(radius) + param_1 + "," + str(radius) + param_2 + str(end_x) + "," + str(end_y)
        label_axis_path = Path(
                d=label_axis_path_desc,
                stroke="none",
                fill="none")
        text_path = Text("") # The text path must go inside a text object. Parameter used here gets ignored
        text_path.add(TextPath(label_axis_path, text=feature_label_text, startOffset="50%", method="align", text_anchor="middle", font_size='8pt',font_style='normal',font_weight='normal',font_family='Liberation Sans', dominant_baseline = "central"))
        path_list = []
        path_list.append(label_axis_path)
        path_list.append(text_path)
    else:
        anchor_value, baseline_value = set_feature_label_anchor_value(total_length, feat_middle)
        path_list = []
        start_x_1: float = (1.1 * radius) * math.cos(
                    math.radians(360.0 * (feat_middle / total_length) - 90))
        start_y_1: float = (1.2 * radius) * math.sin(
            math.radians(360.0 * (feat_middle / total_length) - 90))
        label_path = generate_text_path(feature_label_text, start_x_1, start_y_1, interval = 0, fontsize = '8pt', fontweight = 'normal', font = "Liberation Sans", dominant_baseline = baseline_value, text_anchor = anchor_value)
        path_list.append(label_path)
    return ["label", path_list], available_tracks
''' 
'''   
    

