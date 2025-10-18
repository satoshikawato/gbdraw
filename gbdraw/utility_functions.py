#!/usr/bin/env python
# coding: utf-8
import os
import math
import logging
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from svgwrite.text import Text
import pandas as pd
from pandas import DataFrame
import xml.etree.ElementTree as ET
from typing import List, Dict, Union, Literal
from .find_font_files import get_text_bbox_size_pixels, get_font_dict
from .feature_objects import FeatureObject, GeneObject, RepeatObject
from .file_processing import read_filter_list_file 
import functools 
import math

logger = logging.getLogger(__name__)


def interpolate_color(color_min: str, color_max: str, factor: float) -> str:
    """
    Calculates an intermediate color between two hex colors.

    Args:
        color_min (str): The starting hex color string (e.g., "#RRGGBB").
        color_max (str): The ending hex color string (e.g., "#RRGGBB").
        factor (float): The interpolation factor, from 0.0 to 1.0.

    Returns:
        str: The calculated intermediate hex color string.
    """
    # Convert hex colors to RGB components
    r_min, g_min, b_min = int(color_min[1:3], 16), int(color_min[3:5], 16), int(color_min[5:7], 16)
    r_max, g_max, b_max = int(color_max[1:3], 16), int(color_max[3:5], 16), int(color_max[5:7], 16)

    # Interpolate each RGB component
    r = int(r_min + (r_max - r_min) * factor)
    g = int(g_min + (g_max - g_min) * factor)
    b = int(b_min + (b_max - b_min) * factor)

    # Convert back to hex color string
    return f"#{r:02x}{g:02x}{b:02x}"

@functools.lru_cache(maxsize=4096)
def calculate_bbox_dimensions(text, font_family, font_size, dpi):
    fonts = [font.strip("'") for font in font_family.split(', ')]
    primary_font_family = fonts[0]
    font_file_dict = get_font_dict(tuple(fonts), ("Regular",)) 
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


def create_text_element(text: str, x: float, y: float, font_size: str, font_weight: str, font_family: str, text_anchor: str = "middle", dominant_baseline:str = "middle") -> Text:
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
        text_anchor=text_anchor,
        dominant_baseline=dominant_baseline)


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
                       circular_axis_stroke_color=None, 
                       circular_axis_stroke_width=None,
                       linear_axis_stroke_color=None,
                       linear_axis_stroke_width=None, 
                       line_stroke_color=None, 
                       line_stroke_width=None, 
                       gc_stroke_color=None,
                       linear_definition_font_size=None,
                       circular_definition_font_size=None,
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
                       label_whitelist=None,
                       qualifier_priority=None,
                       outer_label_x_radius_offset=None,
                       outer_label_y_radius_offset=None,
                       inner_label_x_radius_offset=None,
                       inner_label_y_radius_offset=None,
                       comparison_height=None,
                       font_family=None,
                       default_cds_height=None,
                       gc_height=None,
                       scale_style=None,
                       scale_stroke_color=None,
                       scale_stroke_width=None,
                       scale_font_size=None,
                       scale_interval=None,
                       blast_color_min=None,
                       blast_color_max=None)-> dict:
    # Mapping of parameter names to their paths in the config_dict
    label_font_size_circular_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['long']
    label_font_size_circular_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['short']
    label_font_size_linear_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['long']
    label_font_size_linear_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['short']
    
    circular_definition_font_interval = None
    if default_cds_height is not None:
        default_cds_height_short = default_cds_height
        default_cds_height_long = default_cds_height
    else:
        default_cds_height_short = config_dict['canvas']['linear']['default_cds_height']['short']
        default_cds_height_long = config_dict['canvas']['linear']['default_cds_height']['long']
    if circular_definition_font_size is not None:
        circular_definition_font_interval = float(circular_definition_font_size) + 2
    # Process label_blacklist only if the argument was explicitly passed
    if label_blacklist is not None:
        if label_blacklist == "":
            update_config_value(config_dict, 'labels.filtering.blacklist_keywords', [])
        else:
            is_safe_path = False
            try:
                safe_dir = os.path.realpath(os.getcwd())
                requested_path = os.path.realpath(os.path.join(safe_dir, label_blacklist))
                if os.path.commonpath([requested_path, safe_dir]) == safe_dir:
                    is_safe_path = True
            except TypeError:
                is_safe_path = False

            if is_safe_path and os.path.isfile(requested_path):
                with open(requested_path, 'r') as f:
                    keywords = [line.strip() for line in f if line.strip()]
                update_config_value(config_dict, 'labels.filtering.blacklist_keywords', keywords)
            else:
                if not is_safe_path and os.path.exists(label_blacklist):
                     logger.warning(
                        f"Security Warning: Path '{label_blacklist}' is outside the current working directory. "
                        "It will be treated as a comma-separated string."
                    )
                update_config_value(config_dict, 'labels.filtering.blacklist_keywords', [k.strip() for k in label_blacklist.split(',')])
    
    if label_whitelist is not None:
        if label_whitelist == "":
             update_config_value(config_dict, 'labels.filtering.whitelist_df', None)
        else:
            whitelist_df = read_filter_list_file(label_whitelist)
            update_config_value(config_dict, 'labels.filtering.whitelist_df', whitelist_df)
    else:
        update_config_value(config_dict, 'labels.filtering.whitelist_df', None)
    if qualifier_priority is not None and isinstance(qualifier_priority, DataFrame):
        priority_dict = {
            row['feature_type']: [p.strip() for p in row['priorities'].split(',')]
            for _, row in qualifier_priority.iterrows()
        }
        update_config_value(config_dict, 'labels.filtering.qualifier_priority', priority_dict)

    param_paths = {
        'block_stroke_width': 'objects.features.block_stroke_width',
        'block_stroke_color': 'objects.features.block_stroke_color',
        'circular_axis_stroke_color': 'objects.axis.circular.stroke_color',
        'circular_axis_stroke_width': 'objects.axis.circular.stroke_width',
        'line_stroke_color': 'objects.features.line_stroke_color',
        'line_stroke_width': 'objects.features.line_stroke_width',
        'gc_stroke_color': 'objects.gc_content.stroke_color',
        'linear_axis_stroke_color': 'objects.axis.linear.stroke_color',
        'linear_axis_stroke_width': 'objects.axis.linear.stroke_width',
        'linear_definition_font_size': 'objects.definition.linear.font_size',
        'circular_definition_font_size': 'objects.definition.circular.font_size',
        'circular_definition_font_interval': 'objects.definition.circular.interval',
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
        'comparison_height': 'canvas.linear.comparison_height',
        'font_family': 'objects.text.font_family',
        'default_cds_height_long': 'canvas.linear.default_cds_height.long',
        'default_cds_height_short': 'canvas.linear.default_cds_height.short',
        'gc_height': 'canvas.linear.default_gc_height',
        'scale_style': 'objects.scale.style',
        'scale_stroke_color': 'objects.scale.stroke_color',
        'scale_stroke_width': 'objects.scale.stroke_width',
        'scale_font_size': 'objects.scale.font_size',
        'scale_interval': 'objects.scale.interval',
        'blast_color_min': 'objects.blast_match.min_color',
        'blast_color_max': 'objects.blast_match.max_color'
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

def preprocess_label_filtering(label_filtering: dict):
    whitelist_df = label_filtering.get('whitelist_df')
    if whitelist_df is not None and not whitelist_df.empty:
        whitelist_map = {}
        for row in whitelist_df.itertuples(index=False):
            whitelist_map.setdefault(row.feature_type, {}).setdefault(row.qualifier, set()).add(row.keyword)
    else:
        whitelist_map = None

    priority_df = label_filtering.get('qualifier_priority_df')
    if priority_df is not None and not priority_df.empty:
        priority_map = {}
        for row in priority_df.itertuples(index=False):
            priority_map[row.feature_type] = [p.strip() for p in row.priorities.split(',')]
    else:
        priority_map = {}

    label_filtering['whitelist_map'] = whitelist_map
    label_filtering['priority_map'] = priority_map
    return label_filtering

def get_label_text(feature: SeqFeature, label_filtering: dict) -> str:
    feature_type = feature.type
    qualifiers = feature.qualifiers
    whitelist_map = label_filtering.get('whitelist_map')
    blacklist = label_filtering.get('blacklist_keywords', [])
    priority_map = label_filtering.get('priority_map', {})

    # --- Step 1: whitelist check ---
    if whitelist_map:
        is_eligible = False
        rules = whitelist_map.get(feature_type, {})
        for key, values in qualifiers.items():
            if key in rules and any(v in rules[key] for v in values):
                is_eligible = True
                break
        if not is_eligible:
            return ""

    # --- Step 2: priority-based label extraction ---
    priority_list = priority_map.get(feature_type, ['product', 'gene', 'locus_tag',
                                                    'protein_id', 'old_locus_tag', 'note'])
    final_label = ""
    for key in priority_list:
        if key in qualifiers and qualifiers[key]:
            final_label = qualifiers[key][0]
            break

    # --- Step 3: blacklist fallback ---
    if not whitelist_map and any(bl in final_label.lower() for bl in blacklist):
        return ""

    return final_label

def get_coordinates_of_longest_segment(feature_object):
    coords: list[List[Union[str, int, bool]]] = feature_object.coordinates
    if not coords:
        return None, -1
        
    longest_segment_info = None
    max_length = -1

    for coord in coords:
        try:
            # Assuming coord format is [type, strand, start, end, ...]
            start, end = int(coord[2]), int(coord[3])
            length = abs(end - start)
            if length > max_length:
                max_length = length
                longest_segment_info = coord
        except (IndexError, TypeError, ValueError):
            # Skip malformed coordinate entries
            continue
            
    return longest_segment_info, max_length
        

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