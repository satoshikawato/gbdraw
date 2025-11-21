#!/usr/bin/env python
# coding: utf-8

import logging
from pandas import DataFrame
from .circular_path_drawer import generate_text_path
from .create_feature_objects import set_arrow_shoulder
from .utility_functions import  normalize_position_to_linear_track, calculate_bbox_dimensions, get_label_text

# Logging setup
logger = logging.getLogger(__name__)

def calculate_feature_position_factors_linear(strand: str, track_id: int, separate_strands: bool) -> list[float]:
    """
    Calculates position factors for linear feature representation based on strand and strandedness.
    
    Args:
        strand (str): Feature strand ('positive' or 'negative')
        track_id (int): Track number for the feature
        separate_strands (bool): Whether to separate strands
        
    Returns:
        list[float]: Three position factors [top, middle, bottom] for feature placement
    """

    
    # Constants for strand-separated layout
    INITIAL_OFFSET = 0.1           # Base offset from axis
    TRACK_SPACING = 0.5            # Vertical space per track
    FEATURE_HEIGHT = 0.5           # Height of feature
    FEATURE_OFFSET = 0.25          # Half of feature height


        
    def calculate_base_position(track_num: int) -> float:
        """Calculate base y-position for track"""
        return -INITIAL_OFFSET - track_num * TRACK_SPACING
    
    def calculate_track_offset(track_num: int) -> float:
        """Calculate cumulative offset based on track number"""
        if separate_strands:
            return track_num * -1.5 * INITIAL_OFFSET if strand == "positive" else (track_num-1) * -1.5 *INITIAL_OFFSET
        else:
            return track_num * - 3* INITIAL_OFFSET if strand == "positive" else (track_num) * - 3* INITIAL_OFFSET

    # Calculate base position and offsets
    base_pos = calculate_base_position(track_id)
    track_offset = calculate_track_offset(track_id)
    
    if not separate_strands:
        # Simple stacking without strand separation - use positive numbering for all features
        track = track_id # Convert negative track numbers to positive
        return [-track - FEATURE_HEIGHT + track_offset,    # Top
                -track - 0 + track_offset,      # Middle
                -track + FEATURE_HEIGHT + track_offset]    # Bottom

    # Return [top, middle, bottom] positions
    return [
        base_pos - FEATURE_HEIGHT + track_offset,      # Top
        base_pos - FEATURE_OFFSET + track_offset,      # Middle
        base_pos + track_offset                        # Bottom
    ]

def create_intron_path_linear(coord_dict: dict, genome_length: int, 
                            alignment_width: float, genome_size_normalization_factor: float, 
                            cds_height: float, feature_strand: str,
                            separate_strands: bool, feature_track_id: int) -> list[str]:
    """
    Creates a linear SVG path for an intron feature.
    """
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    
    # Get position factors
    factors = calculate_feature_position_factors_linear(
        strand=feature_strand,
        track_id=feature_track_id,
        separate_strands=separate_strands
    )
    
    # Normalize positions
    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)
    
    # Calculate y-position (middle)
    y_position = cds_height * factors[1]
    
    # Create path
    feature_path = f"M {normalized_start},{y_position} L{normalized_end},{y_position} z"
    
    return ["line", feature_path]

def create_rectangle_path_linear(coord_dict: dict, genome_length: int, 
                               alignment_width: float, genome_size_normalization_factor: float, 
                               cds_height: float, feature_strand: str,
                               separate_strands: bool, feature_track_id: int) -> list[str]:
    """Creates a linear SVG path for a rectangular feature."""
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    
    # Get position factors
    factors = calculate_feature_position_factors_linear(
        strand=feature_strand,
        track_id=feature_track_id,
        separate_strands=separate_strands
    )
    
    # Normalize positions
    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)
    
    # Calculate y-positions
    start_y_top = cds_height * factors[0]
    start_y_bottom = cds_height * factors[2]
    
    # Create path
    feature_path = f"M {normalized_start},{start_y_top} " \
                  f"L {normalized_end},{start_y_top} " \
                  f"L {normalized_end},{start_y_bottom} " \
                  f"L {normalized_start},{start_y_bottom} z"
    
    return ["block", feature_path]

def create_arrowhead_path_linear(coord_dict: dict, arrow_length: float, 
                               cds_height: float, feature_strand: str, 
                               genome_length: int, alignment_width: float, 
                               genome_size_normalization_factor: float,
                               separate_strands: bool, feature_track_id: int = 0) -> list[str]:
    """Creates a linear SVG path for an arrowhead feature."""
    # Normalize positions
    normalized_start, normalized_end = normalize_feature_positions(
        coord_dict, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_feat_len = normalized_end - normalized_start
    
    # Calculate arrow properties
    normalized_arrow_length = calculate_normalized_arrow_length(
        arrow_length, genome_length, alignment_width, genome_size_normalization_factor)
    
    # Get arrow positions based on strand
    arrow_start, arrow_end = get_arrow_strand_positions(normalized_start, normalized_end)[coord_dict["feat_strand"]]
    shoulder = set_arrow_shoulder(coord_dict["feat_strand"], arrow_end, normalized_arrow_length)
    
    # Get position factors
    factors = calculate_feature_position_factors_linear(
        strand=feature_strand,
        track_id=feature_track_id,
        separate_strands=separate_strands
    )
    
    return construct_arrowhead_path(
        arrow_start, arrow_end, shoulder, factors,
        normalized_feat_len, normalized_arrow_length, cds_height
    )

def get_arrow_strand_positions(normalized_start: float, normalized_end: float) -> dict[str, list[float]]:
    """Get start and end positions for arrow based on strand."""
    return {
        "positive": [normalized_start, normalized_end],
        "negative": [normalized_end, normalized_start]
    }

def calculate_normalized_arrow_length(arrow_length: float, genome_length: int, 
                                   alignment_width: float, 
                                   genome_size_normalization_factor: float) -> float:
    """Calculate normalized arrow length."""
    return alignment_width * (arrow_length / genome_length) * genome_size_normalization_factor

def normalize_feature_positions(coord_dict: dict, genome_length: int, 
                             alignment_width: float,
                             genome_size_normalization_factor: float) -> tuple[float, float]:
    """Normalize feature start and end positions."""
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])
    
    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor)
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor)
    
    return normalized_start, normalized_end

def construct_arrowhead_path(arrow_start: float, arrow_end: float, shoulder: float,
                           factors: list[float], normalized_feat_len: float,
                           normalized_arrow_length: float, cds_height: float) -> list[str]:
    """Construct the SVG path for an arrowhead."""
    # Calculate points
    point_x, point_y = arrow_end, cds_height * factors[1]
    start_x_1, start_y_1 = arrow_start, cds_height * factors[0]
    start_x_2, start_y_2 = arrow_start, cds_height * factors[2]
    
    if abs(normalized_feat_len) < normalized_arrow_length:
        # Short feature - make entire feature an arrowhead
        feature_path = f"M {start_x_1},{start_y_1} " \
                      f"L {point_x},{point_y} " \
                      f"L {start_x_2},{start_y_2} z"
    else:
        # Normal feature with arrowhead
        end_x_1, end_y_1 = shoulder, cds_height * factors[0]
        end_x_2, end_y_2 = shoulder, cds_height * factors[2]
        feature_path = f"M {start_x_1},{start_y_1} " \
                      f"L {end_x_1},{end_y_1} " \
                      f"L {point_x},{point_y} " \
                      f"L {end_x_2},{end_y_2} " \
                      f"L {start_x_2},{start_y_2} z"
    
    return ["block", feature_path]


def edit_available_tracks(available_tracks, bbox_start, bbox_end):
    track_found =  False
    for track in available_tracks.keys():
        track_end = available_tracks[track][1]
        available_start = track_end + 1
        if bbox_start >= available_start:
            track_factor = list(available_tracks).index(track) + 1
            available_tracks[track][1] = bbox_end
            track_found = True
            return available_tracks, track_factor
    if track_found ==  False:
        new_track_id = "track_{}".format((len(available_tracks.keys()) + 1))
        available_tracks[new_track_id] = [bbox_start, bbox_end]
        track_factor = len(available_tracks.keys())
    return available_tracks, track_factor



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


def calculate_gc_content_path_desc(start_x: float, start_y: float, gc_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float, dinucleotide: str) -> str:
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
    column: str = f'{dinucleotide} content'
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

def calculate_gc_skew_path_desc(start_x: float, start_y: float, skew_df: DataFrame, record_len: int, alignment_width: float, genome_size_normalization_factor: float, track_height: float) -> str:
    """
    Generates the SVG path description for linear GC skew representation.
    (This is the new function to be added)
    """
    coodinates_list: list = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column = [col for col in skew_df.columns if 'skew' in col and 'cumulative' not in col.lower()][0]
    mean = float(skew_df[column].mean())
    max_diff = float((skew_df[column] - mean).abs().max()) if float((skew_df[column] - mean).abs().max()) > 0 else 1.0

    for index, row in skew_df.iterrows():
        value = float(row[column])
        corrdinate, x_corrdinate = calculate_corrdinate(
            index, value, mean, max_diff, record_len, alignment_width, genome_size_normalization_factor, track_height)
        coodinates_list.append(corrdinate)

    penultimate_coordinate: str = "L{} {}".format(str(x_corrdinate), str(start_y))
    coodinates_list.append(penultimate_coordinate)
    end_coordinate: str = "L{} {}".format(str(start_x), str(start_y))
    coodinates_list.append(end_coordinate)
    gc_skew_desc: str = "{}".format(''.join(coodinates_list))
    gc_skew_desc += "z"
    return gc_skew_desc