#!/usr/bin/env python
# coding: utf-8

import logging
import sys

from pandas import DataFrame

from .create_feature_objects import set_arrow_shoulder
from .utility_functions import normalize_position_to_linear_track
# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def create_intron_path_linear(coord_dict: dict, genome_length: int, alignment_width: float, genome_size_normalization_factor: float, cds_height: float, strandedness: str) -> list[str]:
    """
    Creates a linear SVG path for an intron feature.

    Args:
        coord_dict (dict): Dictionary containing coordinates and other data of the feature.
        genome_length (int): Total length of the genome.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for genome size.
        cds_height (float): Height of the coding sequence tracks.
        strandedness (str): Strand orientation of the feature.

    Returns:
        list[str]: A list containing the type of the path and the path data.
    """
    feat_start: int = coord_dict['feat_start']
    feat_end: int = coord_dict['feat_end']
    feat_strand: str = coord_dict["feat_strand"]
    factors: list[float] = calculate_feature_position_factors_linear(
        feat_strand, strandedness)
    normalized_start: float = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end: float = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)
    start_x: float = normalized_start
    start_y: float = cds_height * factors[1]
    end_x: float = normalized_end
    end_y: float = cds_height * factors[1]
    feature_path: str = "M " + \
        str(start_x) + "," + str(start_y) + "L" + \
        str(end_x) + "," + str(end_y) + " z"
    return ["line", feature_path]


def calculate_feature_position_factors_linear(strand: str, strandedness: str) -> list[float]:
    """
    Calculates position factors for linear feature representation based on strand and strandedness.

    Args:
        strand (str): Strand of the feature ('positive' or 'negative').
        strandedness (str): Strandedness property of the genomic data.

    Returns:
        list[float]: A list of position factors to determine the vertical placement of the feature.
    """
    if strandedness:
        OFFSET = 0.1
        HEIGHT = 0.45
        factors_positive: list[float] = [-(HEIGHT + OFFSET), -
                                         ((HEIGHT + 2 * OFFSET) / 2), -OFFSET]
        factors_negative: list[float] = [
            OFFSET, ((HEIGHT + 2 * OFFSET) / 2), (HEIGHT + OFFSET)]
        factors: list[float] = factors_positive if strand == "positive" else factors_negative
    else:
        factors = [-0.5, 0, 0.5]
    return factors


def normalize_feature_positions(coord_dict: dict, genome_length: int, alignment_width: float, genome_size_normalization_factor: float) -> tuple[float, float]:
    """
    Normalizes the feature positions for linear display.

    Args:
        coord_dict (dict): Dictionary containing feature coordinates.
        genome_length (int): Total length of the genome.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.

    Returns:
        tuple[float, float]: Normalized start and end positions of the feature.
    """
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])
    normalized_start: float = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end: float = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)
    return normalized_start, normalized_end


def calculate_normalized_arrow_length(arrow_length: float, genome_length: int, alignment_width: float, genome_size_normalization_factor: float) -> float:
    """
    Calculates the normalized length of an arrow for feature representation.

    Args:
        arrow_length (float): The specified arrow length.
        genome_length (int): Total length of the genome.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.

    Returns:
        float: Normalized arrow length.
    """
    return alignment_width * (arrow_length / genome_length) * genome_size_normalization_factor


def get_arrow_strand_positions(normalized_start: float, normalized_end: float) -> dict[str, list[float]]:
    """
    Determines the start and end positions of an arrow based on strand orientation.

    Args:
        normalized_start (float): Normalized start position of the feature.
        normalized_end (float): Normalized end position of the feature.

    Returns:
        dict[str, list[float]]: Dictionary with strand keys and corresponding position values.
    """
    return {
        "positive": [normalized_start, normalized_end],
        "negative": [normalized_end, normalized_start]
    }


def create_arrowhead_path_linear(coord_dict: dict, arrow_length: float, cds_height: float, strandedness: str, genome_length: int, alignment_width: float, genome_size_normalization_factor: float) -> list[str]:
    """
    Creates a linear SVG path for an arrowhead feature.

    Args:
        coord_dict (dict): Dictionary containing coordinates and other data of the feature.
        arrow_length (float): Length of the arrowhead.
        cds_height (float): Height of the coding sequence tracks.
        strandedness (str): Strand orientation of the feature.
        genome_length (int): Total length of the genome.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.

    Returns:
        list[str]: A list containing the type of the path and the path data.
    """
    normalized_start: float
    normalized_end: float
    normalized_start, normalized_end = normalize_feature_positions(
        coord_dict, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_feat_len: float = normalized_end - normalized_start
    normalized_arrow_length: float = calculate_normalized_arrow_length(
        arrow_length, genome_length, alignment_width, genome_size_normalization_factor)

    arrow_strand_dict: dict[str, list[float]] = {
        "positive": [normalized_start, normalized_end],
        "negative": [normalized_end, normalized_start]
    }
    arrow_start: float
    arrow_end: float
    arrow_start, arrow_end = arrow_strand_dict[coord_dict["feat_strand"]]
    shoulder: float = set_arrow_shoulder(
        coord_dict["feat_strand"], arrow_end, normalized_arrow_length)

    factors: list[float] = calculate_feature_position_factors_linear(
        coord_dict["feat_strand"], strandedness)
    return construct_arrowhead_path(arrow_start, arrow_end, shoulder, factors, normalized_feat_len, normalized_arrow_length, cds_height)


def construct_arrowhead_path(arrow_start: float, arrow_end: float, shoulder: float, factors: list[float], normalized_feat_len: float, normalized_arrow_length: float, cds_height: float) -> list[str]:
    """
    Constructs the SVG path description for an arrowhead feature.

    Args:
        arrow_start (float): Normalized start position of the arrow.
        arrow_end (float): Normalized end position of the arrow.
        shoulder (float): Position of the arrow shoulder.
        factors (list[float]): Position factors for the feature.
        normalized_feat_len (float): Normalized length of the feature.
        normalized_arrow_length (float): Normalized length of the arrow.
        cds_height (float): Height of the coding sequence tracks.

    Returns:
        list[str]: A list containing the type of the path and the path data.
    """
    point_x: float
    point_y: float
    start_x_1: float
    start_y_1: float
    start_x_2: float
    start_y_2: float
    point_x, point_y = arrow_end, cds_height * factors[1]
    start_x_1, start_y_1 = arrow_start, cds_height * factors[0]
    start_x_2, start_y_2 = arrow_start, cds_height * factors[2]

    if abs(normalized_feat_len) < normalized_arrow_length:
        # If the feature is shorter than the arrow length, make the entire feature an arrowhead
        feature_path: str = f"M {start_x_1},{start_y_1} L {
            point_x},{point_y} L {start_x_2},{start_y_2} z"
    else:
        # If the feature is longer than the arrow length, construct a path with an arrowhead at the end
        end_x_1: float
        end_y_1: float
        end_x_2: float
        end_y_2: float
        end_x_1, end_y_1 = shoulder, cds_height * factors[0]
        end_x_2, end_y_2 = shoulder, cds_height * factors[2]
        feature_path = f"M {start_x_1},{start_y_1} L {end_x_1},{end_y_1} L {
            point_x},{point_y} L {end_x_2},{end_y_2} L {start_x_2},{start_y_2} z"

    return ["block", feature_path]


def create_rectangle_path_linear(coord_dict: dict, genome_length: int, alignment_width: float, genome_size_normalization_factor: float, cds_height: float, strandedness: str) -> list[str]:
    """
    Creates a linear SVG path for a rectangular feature.

    Args:
        coord_dict (dict): Dictionary containing coordinates and other data of the feature.
        genome_length (int): Total length of the genome.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.
        cds_height (float): Height of the coding sequence tracks.
        strandedness (str): Strand orientation of the feature.

    Returns:
        list[str]: A list containing the type of the path and the path data.
    """
    feat_start: float
    feat_end: float
    feat_start, feat_end = coord_dict['feat_start'], coord_dict['feat_end']
    factors: list[float] = calculate_feature_position_factors_linear(
        coord_dict["feat_strand"], strandedness)

    # Normalize start and end positions
    normalized_start: float = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end: float = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)

    # Construct the rectangle path
    start_y_top: float
    start_y_bottom: float
    end_y_top: float
    end_y_bottom: float
    start_y_top, start_y_bottom = cds_height * \
        factors[0], cds_height * factors[2]
    end_y_top, end_y_bottom = start_y_top, start_y_bottom
    feature_path: str = f"M {normalized_start},{start_y_top} L {normalized_end},{end_y_top} " \
        f"L {normalized_end},{end_y_bottom} L {
            normalized_start},{start_y_bottom} z"

    return ["block", feature_path]


def calculate_corrdinate(index: int, value: float, mean: float, max_diff: float, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float) -> str:
    """
    Calculates the coordinate for a point in the GC content path.

    Args:
        index (int): Position index in the genome.
        value (float): GC content value at the index.
        mean (float): Mean value of the GC content.
        max_diff (float): Maximum difference from the mean in the GC content.
        record_len (int): Length of the genomic record.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.
        track_height (float): Height of the GC content track.

    Returns:
        str: Coordinate string for the GC content path.
    """
    diff: float = (value - mean)
    x_corrdinate: float = normalize_position_to_linear_track(
        index, record_len, alignment_width, genome_size_normalization_factor)
    y_corrdinate: float = - (0.5 * track_height * (diff / max_diff))
    corrdinate: str = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
    return corrdinate, x_corrdinate


def calculate_gc_content_path_desc(start_x: float, start_y: float, gc_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float) -> str:
    """
    Calculates the SVG path description for GC content visualization.

    Args:
        start_x (float): Starting x-coordinate for the path.
        start_y (float): Starting y-coordinate for the path.
        gc_df (DataFrame): DataFrame containing GC content data.
        record_len (int): Length of the genomic record.
        alignment_width (float): Width of the alignment area.
        genome_size_normalization_factor (float): Normalization factor for the genome size.
        track_height (float): Height of the GC content track.

    Returns:
        str: SVG path description for the GC content.
    """
    coodinates_list: list = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column = 'GC content'
    mean = float(gc_df[column].mean())
    max_diff = float((gc_df[column] - mean).abs().max())
    for index, row in gc_df.iterrows():
        value = float(row[column])
        corrdinate, x_corrdinate = calculate_corrdinate(
            index, value, mean, max_diff, record_len, alignment_width, genome_size_normalization_factor, track_height)
        coodinates_list.append(corrdinate)
    penultimate_coordinate: str = "L{} {}".format(str(x_corrdinate), str(start_y))
    coodinates_list.append(penultimate_coordinate) 
    end_coordinate: str = "L{} {}".format(str(start_x), str(start_y))
    coodinates_list.append(end_coordinate)
    gc_content_desc: str = "{}".format(''.join(coodinates_list))
    gc_content_desc += "z"
    return gc_content_desc

def generate_gc_skew_path_desc(start_x: float, start_y: float, gc_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float) -> str:
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