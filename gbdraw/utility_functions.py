#!/usr/bin/env python
# coding: utf-8
import os
import math
from Bio.SeqRecord import SeqRecord
from svgwrite.text import Text
from pandas import DataFrame
import xml.etree.ElementTree as ET
from typing import List, Dict, Union, Literal
from .find_font_files import get_text_bbox_size_pixels, get_font_dict
from .feature_objects import FeatureObject, GeneObject, RepeatObject

def calculate_bbox_dimensions(text, font_family, font_size, dpi):
    fonts = [font.strip("'") for font in font_family.split(', ')]
    primary_font_family = fonts[0]
    font_file_dict = get_font_dict(fonts, ["Regular"]) 
    font_path = font_file_dict[primary_font_family]["Regular"]
    bbox_width_px, bbox_height_px = get_text_bbox_size_pixels(font_path, text, font_size, dpi)
    return bbox_width_px, bbox_height_px

def normalize_position_linear(position: float, longest_genome: int, alignment_width: float) -> float:
    """
    Normalizes a genomic position for linear representation.

    Args:
        position (float): The position in the genome to be normalized.
        longest_genome (int): The length of the longest genome in the dataset.
        alignment_width (float): The width of the linear alignment area.

    Returns:
        float: Normalized position as a proportion of the alignment width.
    """
    return alignment_width * (position / longest_genome)


def create_dict_for_sequence_lengths(records: list[SeqRecord]) -> dict[str, int]:
    """
    Creates a dictionary mapping sequence IDs to their lengths.

    Args:
        records (list[SeqRecord]): A list of SeqRecord objects.

    Returns:
        dict[str, int]: A dictionary where keys are sequence IDs and values are sequence lengths.
    """
    return {record.id: len(record.seq) for record in records}


def normalize_position_to_linear_track(
        position: int,
        genome_length: float,
        alignment_width: float,
        genome_size_normalization_factor: float) -> float:
    """
    Normalizes a genomic position to a linear track.

    Args:
        position (int): Genomic position to normalize.
        genome_length (float): Total length of the genome.
        alignment_width (float): Width of the alignment track.
        genome_size_normalization_factor (float): Normalization factor based on the genome size.

    Returns:
        float: Normalized position on the linear track.
    """
    normalized_position: float = alignment_width * \
        (position / genome_length) * genome_size_normalization_factor
    return normalized_position


def create_text_element(text: str, x: float, y: float, font_size: str, font_weight: str, font_family: str) -> Text:
    """
    Creates an SVG text element.

    Args:
        text (str): Text content for the element.
        x (float): x-coordinate for the text placement.
        y (float): y-coordinate for the text placement.
        font_size (str): Font size for the text.
        font_weight (str): Font weight for the text.
        font_family (str): Font family for the text.

    Returns:
        Text: An SVG text element.
    """
    return Text(
        text,
        insert=(x, y),
        stroke='none',
        fill='black',
        font_size=font_size,
        font_weight=font_weight,
        font_family=font_family,
        text_anchor="middle",
        dominant_baseline="middle")


def parse_mixed_content_text(input_text: str) -> List[Dict[str, Union[str, bool, None]]]:
    """
    Parses text containing mixed regular and italic content.

    This function is particularly useful for parsing biological names where certain parts 
    need to be italicized.

    Args:
        input_text (str): Text string containing mixed content.

    Returns:
        List[Dict[str, Union[str, bool, None]]]: A list of text parts, each with the text content 
        and a boolean indicating whether it is italicized.
    """
    parts: List[Dict[str, Union[str, bool, None]]] = []
    try:
        wrapped_text: str = f"<root>{input_text}</root>"
        root: ET.Element = ET.fromstring(wrapped_text)
        if list(root):  # Check if the root element has any child elements
            for element in root:
                if element.tag == 'i':
                    parts.append({'text': element.text, 'italic': True})
                else:
                    parts.append({'text': element.text, 'italic': False})
                if element.tail:
                    parts.append({'text': element.tail, 'italic': False})
        else:
            # If there are no child elements, treat the entire string as regular text
            parts.append({'text': root.text, 'italic': False})
    except ET.ParseError:
        # Fallback for non-XML content
        parts.append({'text': input_text, 'italic': False})
    return parts


def suppress_gc_content_and_skew(suppress_gc: bool, suppress_skew: bool) -> tuple[bool, bool]:
    """
    Determines whether to show or suppress GC content and GC skew.

    Args:
        suppress_gc (bool): Flag to suppress GC content if True.
        suppress_skew (bool): Flag to suppress GC skew if True.

    Returns:
        tuple[bool, bool]: Tuple indicating whether to show GC content and GC skew.
    """
    show_gc, show_skew = True, True
    if suppress_gc == True:
        show_gc = False
    if suppress_skew == True:
        show_skew = False
    return show_gc, show_skew

def update_config_value(config_dict, path, value):
    """ Helper function to update nested dictionary """
    keys = path.split('.')
    for key in keys[:-1]:
        config_dict = config_dict.setdefault(key, {})
    config_dict[keys[-1]] = value


def modify_config_dict(config_dict, 
                       block_stroke_width=None, 
                       block_stroke_color=None, 
                       line_stroke_color=None, 
                       line_stroke_width=None, 
                       gc_stroke_color=None,
                       label_font_size=None, 
                       show_gc=None, 
                       show_skew=None, 
                       show_labels=None, 
                       align_center=None, 
                       cicular_width_with_labels=None, 
                       track_type=None, 
                       strandedness=None, 
                       resolve_overlaps=None, 
                       allow_inner_labels=None,
                       label_radius_offset=None,
                       label_blacklist=None,
                       qualifier_priority=None,
                       outer_label_x_radius_offset=None,
                       outer_label_y_radius_offset=None,
                       inner_label_x_radius_offset=None,
                       inner_label_y_radius_offset=None) -> dict:
    # Mapping of parameter names to their paths in the config_dict
    label_font_size_circular_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['long']
    label_font_size_circular_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['short']
    label_font_size_linear_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['long']
    label_font_size_linear_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['short']

    if label_blacklist:
        if os.path.isfile(label_blacklist):
            with open(label_blacklist, 'r') as f:
                keywords = [line.strip() for line in f if line.strip()]
            update_config_value(config_dict, 'labels.filtering.blacklist_keywords', keywords)
        else:
            update_config_value(config_dict, 'labels.filtering.blacklist_keywords', [k.strip() for k in label_blacklist.split(',')])

    # NEW: ラベル優先順位の処理を追加
    if qualifier_priority is not None and isinstance(qualifier_priority, DataFrame):
        # DataFrameを行ごとに処理し、{feature_type: [priority1, priority2, ...]} の辞書を作成
        priority_dict = {
            row['feature_type']: [p.strip() for p in row['priorities'].split(',')]
            for _, row in qualifier_priority.iterrows()
        }
        update_config_value(config_dict, 'labels.filtering.qualifier_priority', priority_dict)

    param_paths = {
        'block_stroke_width': 'objects.features.block_stroke_width',
        'block_stroke_color': 'objects.features.block_stroke_color',
        'line_stroke_color': 'objects.features.line_stroke_color',
        'line_stroke_width': 'objects.features.line_stroke_width',
        'gc_stroke_color': 'objects.gc_content.stroke_color',
        'label_font_size_circular_long': 'labels.font_size.long',
        'label_font_size_circular_short': 'labels.font_size.short',
        'label_font_size_linear_long': 'labels.font_size.linear.long',
        'label_font_size_linear_short': 'labels.font_size.linear.short',
        'strandedness': 'canvas.strandedness',
        'show_gc': 'canvas.show_gc',
        'show_skew': 'canvas.show_skew',
        'show_labels': 'canvas.show_labels',
        'align_center': 'canvas.linear.align_center',
        'track_type': 'canvas.circular.track_type',
        'resolve_overlaps': 'canvas.resolve_overlaps',
        'allow_inner_labels': 'canvas.circular.allow_inner_labels',
        'outer_label_x_radius_offset': 'labels.unified_adjustment.outer_labels.x_radius_offset',
        'outer_label_y_radius_offset': 'labels.unified_adjustment.outer_labels.y_radius_offset',
        'inner_label_x_radius_offset': 'labels.unified_adjustment.inner_labels.x_radius_offset',
        'inner_label_y_radius_offset': 'labels.unified_adjustment.inner_labels.y_radius_offset',
    }
    # Update the config_dict for each specified parameter
    for param, path in param_paths.items():
        value = locals()[param]
        if value is not None:
            update_config_value(config_dict, path, value)
    return config_dict

def determine_output_file_prefix(gb_records, output_prefix, record_count, accession):
    """
    Determines the output file prefix based on the number of records, output prefix, record count, and accession.

    Args:
        gb_records (list): List of GenBank records.
        output_prefix (str): The prefix provided for the output file.
        record_count (int): The current record count in a loop of multiple records.
        accession (str): The accession number of the current record.

    Returns:
        str: The determined output file prefix.
    """
    if len(gb_records) > 1 and output_prefix is not None:
        return "{}_{}".format(output_prefix, record_count)
    elif len(gb_records) == 1 and output_prefix is not None:
        return output_prefix
    else:
        return accession
    
def check_feature_presence(records: Union[List[SeqRecord], SeqRecord], features_list: List[str]) -> dict:
    # If records is a single SeqRecord, wrap it in a list
    if isinstance(records, SeqRecord):
        records = [records]

    features_to_be_checked = set(features_list)  # Use a set for efficient lookup
    features_present = []  # Initialize all as False

    for record in records:
        for feature in record.features:
            if feature.type in features_to_be_checked:
                features_present.append(feature.type)
                # Once a feature type is found, no need to check further for the same type
                features_to_be_checked.remove(feature.type)
                if not features_to_be_checked:  # Break early if all selected features are found
                    break
    return features_present

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


def get_label_text(seq_feature, filtering_config) -> str:
    """
    Extracts a label for a feature based on a priority list of qualifiers
    determined by the feature type (from config.toml), and filters it against a blacklist.
    """
    feature_type = seq_feature.type
    priority_config = filtering_config['qualifier_priority']
    blacklist = filtering_config['blacklist_keywords']

    # Determine which priority list to use based on the feature's type string.
    # This logic corresponds to how GeneObject, RepeatObject, etc., are defined.
    if feature_type in ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'gene']:
        priority_list = priority_config['gene']
    elif feature_type == 'repeat_region':
        priority_list = priority_config['repeat']
    else:
        priority_list = priority_config['feature']

    text = ''
    for priority in priority_list:
        if hasattr(seq_feature, priority) and seq_feature.__getattribute__(priority):
            potential_text = seq_feature.__getattribute__(priority)
            if isinstance(potential_text, list):
                potential_text = ', '.join(text)  # Join list items with a comma
            if not any(keyword in potential_text.lower() for keyword in blacklist):
                text = potential_text
                break  # Stop at the first valid label found
    return text # Return an empty string if no suitable label is found

def get_coordinates_of_longest_segment(feature_object):
    coords: list[List[Union[str, int, bool]]] = feature_object.coordinates
    for coord in coords:
        exon_strand: str = get_strand(exon_line.strand)  # type: ignore
    for coord in coords:
        coord_dict: Dict[str, Union[str, int]] = {
            'coord_type': str(coord[0]),
            'coord_strand': coord[2],
            'coord_start': coord[3],
            'coord_end': coord[4]}
        coord_type: str = str(coord_dict['coord_type'])
        

# Function to convert degrees to radians
def deg_to_rad(deg):
    return deg * math.pi / 180.0

# Function to calculate the x and y coordinates from the center
def calculate_coordinates(center_x, center_y, radius, angle_degrees):
    angle_radians = deg_to_rad(angle_degrees)
    x = center_x + radius * math.cos(angle_radians)
    y = center_y + radius * math.sin(angle_radians)
    return x, y


def determine_length_parameter(record_length: int, length_threshold: int) -> str:
    if record_length < length_threshold:
        return "short"
    else:
        return "long"
    

def calculate_cds_ratio(track_ratio, length_param, track_ratio_factor):
    if length_param == "short":
        cds_ratio = float(track_ratio * track_ratio_factor)
        offset = float(0.01)
    else:
        cds_ratio = float(track_ratio * track_ratio_factor)
        offset = float(0.005)
    return cds_ratio, offset