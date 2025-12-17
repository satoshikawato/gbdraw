#!/usr/bin/env python
# coding: utf-8
import os
import math
import sys
import logging
import functools
import xml.etree.ElementTree as ET
from typing import Optional, List, Dict, Union
from importlib import resources

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from svgwrite.text import Text
import pandas as pd
from pandas import DataFrame

# PILではなく、fontToolsを使用
from fontTools.ttLib import TTFont

from .feature_objects import FeatureObject, GeneObject, RepeatObject

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
#  Logic ported from find_font_files.py to remove dependency
# ------------------------------------------------------------------
def get_text_bbox_size_pixels(font_path, text, font_size, dpi):
    """
    Directly parses the font file using fontTools to calculate text dimensions.
    This logic was originally in find_font_files.py.
    """
    try:
        font = TTFont(font_path)
    except Exception as e:
        logger.warning(f"Failed to load font file {font_path}: {e}")
        # Fallback approximation
        return len(str(text)) * float(font_size) * 0.6, float(font_size)

    hmtx = font['hmtx']
    cmap = font['cmap']
    head = font['head']
    units_per_em = head.unitsPerEm
    
    # Get cmap (Format 4 or 12 preferred)
    t = cmap.getcmap(3, 1).cmap
    
    total_width = 0
    ymaxes = []
    ymins = []
    
    text_str = str(text)
    if not text_str:
        return 0.0, 0.0

    rsb_previous = 0
    
    for i, char in enumerate(text_str):
        char_code = ord(char)
        if char_code in t:
            glyph_index = t[char_code]
        else:
            glyph_index = 0
            
        try:
            # Try to get vertical metrics from glyf table if available
            if 'glyf' in font:
                g = font['glyf'][glyph_index]
                ymax = g.yMax if hasattr(g, 'yMax') else 0
                ymin = g.yMin if hasattr(g, 'yMin') else 0
                xmax = g.xMax if hasattr(g, 'xMax') else 0
            else:
                # Fallback for OTF/CFF fonts without glyf table
                ymax = head.yMax
                ymin = head.yMin
                xmax = 0 
        except Exception:
            ymax, ymin, xmax = 0, 0, 0

        ymaxes.append(ymax)
        ymins.append(ymin)
        
        # Calculate horizontal advance
        try:
            advance_width, lsb = hmtx[glyph_index]
        except KeyError:
             advance_width, lsb = 1000, 0

        if xmax == 0:
             advance_width -= rsb_previous
        
        if lsb > 0:
            total_width -= lsb
        else:
             total_width += lsb
             
        rsb = (xmax - advance_width)
        total_width += advance_width
        
        if rsb > 0:
             total_width -= rsb
        else:
             total_width += rsb
             
        if i == 0:
             total_width += lsb
             
        rsb_previous = rsb

    # Formula: pixels = (design_units * font_size * dpi) / (72 * units_per_em)
    scale_factor = (float(font_size) * dpi) / (72 * units_per_em)
    
    text_width_pixels = total_width * scale_factor
    
    if ymaxes and ymins:
        max_y = max(ymaxes)
        min_y = min(ymins)
        text_height_pixels = (max_y + abs(min_y)) * scale_factor
    else:
        text_height_pixels = float(font_size)

    return text_width_pixels, text_height_pixels


# ------------------------------------------------------------------
#  Main BBox Calculation Function
# ------------------------------------------------------------------
@functools.lru_cache(maxsize=4096)
def calculate_bbox_dimensions(text, font_family, font_size, dpi):
    """
    Calculates bounding box dimensions using bundled font files in gbdraw package.
    Locates the font in gbdraw.data.fonts using importlib.resources.
    """
    target_font_path = None
    
    try:
        # Access the bundled fonts directory
        # Assumes structure: gbdraw/data/fonts/*.ttf
        font_dir = resources.files('gbdraw.data').joinpath('fonts')
        
        # List available TTF/OTF files
        available_fonts = [f for f in font_dir.iterdir() if f.name.lower().endswith(('.ttf', '.otf'))]
        
        if not available_fonts:
            # raise FileNotFoundError("No font files found in gbdraw.data.fonts")
            pass # Fallback to approximation below

        else:
            # Simple matching logic
            requested = font_family.split(',')[0].strip().lower().replace(" ", "")
            
            for f in available_fonts:
                if requested in f.name.lower():
                    target_font_path = str(f)
                    break
            
            # Fallback to the first available font
            if not target_font_path:
                target_font_path = str(available_fonts[0])

            # Calculate using fontTools
            return get_text_bbox_size_pixels(target_font_path, text, font_size, dpi)

    except Exception as e:
        logger.debug(f"Font calculation failed ({e}). Using approximation.")
    
    # Fallback Approximation
    try:
        f_size = float(font_size)
    except:
        f_size = 12.0
    return len(str(text)) * f_size * 0.6, f_size


# ------------------------------------------------------------------
#  Original Utility Functions
# ------------------------------------------------------------------

def read_qualifier_priority_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = ['feature_type', 'priorities']

    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines='error',  
            engine='python'
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in qualifier priority file '{filepath}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Qualifier priority file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        sys.exit(1)

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    logger.info(f"Successfully loaded qualifier priority from {filepath}")
    return df

def read_filter_list_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = ['feature_type', 'qualifier', 'keyword']

    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=required_cols,
            dtype=str,
            comment='#',
            on_bad_lines='error',
            engine='python'
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in filter list file '{filepath}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Filter list file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        sys.exit(1)

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    logger.info(f"Successfully loaded filter list from {filepath}")
    return df

def interpolate_color(color_min: str, color_max: str, factor: float) -> str:
    r_min, g_min, b_min = int(color_min[1:3], 16), int(color_min[3:5], 16), int(color_min[5:7], 16)
    r_max, g_max, b_max = int(color_max[1:3], 16), int(color_max[3:5], 16), int(color_max[5:7], 16)

    r = int(r_min + (r_max - r_min) * factor)
    g = int(g_min + (g_max - g_min) * factor)
    b = int(b_min + (b_max - b_min) * factor)

    return f"#{r:02x}{g:02x}{b:02x}"


def normalize_position_linear(position: float, longest_genome: int, alignment_width: float) -> float:
    return alignment_width * (position / longest_genome)


def create_dict_for_sequence_lengths(records: list[SeqRecord]) -> dict[str, int]:
    return {record.id: len(record.seq) for record in records}


def normalize_position_to_linear_track(
        position: int,
        genome_length: float,
        alignment_width: float,
        genome_size_normalization_factor: float) -> float:
    normalized_position: float = alignment_width * \
        (position / genome_length) * genome_size_normalization_factor
    return normalized_position


def create_text_element(text: str, x: float, y: float, font_size: str, font_weight: str, font_family: str, text_anchor: str = "middle", dominant_baseline:str = "middle") -> Text:
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
    parts: List[Dict[str, Union[str, bool, None]]] = []
    try:
        wrapped_text: str = f"<root>{input_text}</root>"
        root: ET.Element = ET.fromstring(wrapped_text)
        if list(root):
            for element in root:
                if element.tag == 'i':
                    parts.append({'text': element.text, 'italic': True})
                else:
                    parts.append({'text': element.text, 'italic': False})
                if element.tail:
                    parts.append({'text': element.tail, 'italic': False})
        else:
            parts.append({'text': root.text, 'italic': False})
    except ET.ParseError:
        parts.append({'text': input_text, 'italic': False})
    return parts


def suppress_gc_content_and_skew(suppress_gc: bool, suppress_skew: bool) -> tuple[bool, bool]:
    show_gc, show_skew = True, True
    if suppress_gc == True:
        show_gc = False
    if suppress_skew == True:
        show_skew = False
    return show_gc, show_skew


def update_config_value(config_dict, path, value):
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
                       blast_color_max=None,
                       legend_box_size=None,
                       legend_font_size=None,
                       normalize_length=None
                       )-> dict:
    
    label_font_size_circular_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['long']
    label_font_size_circular_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['short']
    label_font_size_linear_long = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['long']
    label_font_size_linear_short = label_font_size if label_font_size is not None else config_dict['labels']['font_size']['linear']['short']
        
    if linear_definition_font_size is not None:
        linear_definition_font_size_short = linear_definition_font_size
        linear_definition_font_size_long = linear_definition_font_size
    else:
        linear_definition_font_size_short = config_dict['objects']['definition']['linear']['font_size']['short']
        linear_definition_font_size_long = config_dict['objects']['definition']['linear']['font_size']['long']

    if default_cds_height is not None:
        default_cds_height_short = default_cds_height
        default_cds_height_long = default_cds_height
    else:
        default_cds_height_short = config_dict['canvas']['linear']['default_cds_height']['short']
        default_cds_height_long = config_dict['canvas']['linear']['default_cds_height']['long']
        
    if circular_definition_font_size is not None:
        circular_definition_font_interval = float(circular_definition_font_size) + 2
    else:
        circular_definition_font_interval = None

    if legend_box_size is not None:
        legend_box_size_short = legend_box_size
        legend_box_size_long = legend_box_size
    else:
        legend_box_size_short = config_dict['objects']['legends']['color_rect_size']['short']
        legend_box_size_long = config_dict['objects']['legends']['color_rect_size']['long']

    if legend_font_size is not None:
        legend_font_size_short = legend_font_size
        legend_font_size_long = legend_font_size
    else:
        legend_font_size_short = config_dict['objects']['legends']['font_size']['short']
        legend_font_size_long = config_dict['objects']['legends']['font_size']['long']

    if scale_font_size is not None:
        scale_font_size_short = scale_font_size
        scale_font_size_long = scale_font_size
    else:
        scale_font_size_short = config_dict['objects']['scale']['font_size']['short']
        scale_font_size_long = config_dict['objects']['scale']['font_size']['long']       

    if linear_axis_stroke_width is not None:
        linear_axis_stroke_width_short = linear_axis_stroke_width
        linear_axis_stroke_width_long = linear_axis_stroke_width
    else:
        linear_axis_stroke_width_short = config_dict['objects']['axis']['linear']['stroke_width']['short']
        linear_axis_stroke_width_long = config_dict['objects']['axis']['linear']['stroke_width']['long']

    if line_stroke_width is not None:
        line_stroke_width_short = line_stroke_width
        line_stroke_width_long = line_stroke_width
    else:
        line_stroke_width_short = config_dict['objects']['features']['line_stroke_width']['short']
        line_stroke_width_long = config_dict['objects']['features']['line_stroke_width']['long']

    if line_stroke_color is not None:
        line_stroke_color = line_stroke_color
    else:
        line_stroke_color = config_dict['objects']['features']['line_stroke_color']

    if block_stroke_color is not None:
        block_stroke_color = block_stroke_color
    else:
        config_dict['objects']['features']['block_stroke_color']

    if block_stroke_width is not None:
        block_stroke_width_short = block_stroke_width
        block_stroke_width_long = block_stroke_width
    else:
        block_stroke_width_short = config_dict['objects']['features']['block_stroke_width']['short']
        block_stroke_width_long = config_dict['objects']['features']['block_stroke_width']['long']        

    if circular_axis_stroke_width is not None:
        circular_axis_stroke_width_short = circular_axis_stroke_width
        circular_axis_stroke_width_long = circular_axis_stroke_width        
    else:
        circular_axis_stroke_width_short = config_dict['objects']['axis']['circular']['stroke_width']['short']
        circular_axis_stroke_width_long = config_dict['objects']['axis']['circular']['stroke_width']['long']

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
        'block_stroke_color': 'objects.features.block_stroke_color',
        'block_stroke_width_short': 'objects.features.block_stroke_width.short',
        'block_stroke_width_long': 'objects.features.block_stroke_width.long',
        'circular_axis_stroke_color': 'objects.axis.circular.stroke_color',
        'circular_axis_stroke_width_short': 'objects.axis.circular.stroke_width.short',
        'circular_axis_stroke_width_long': 'objects.axis.circular.stroke_width.long',
        'line_stroke_color': 'objects.features.line_stroke_color',
        'line_stroke_width_short': 'objects.features.line_stroke_width.short',
        'line_stroke_width_long': 'objects.features.line_stroke_width.long',
        'gc_stroke_color': 'objects.gc_content.stroke_color',
        'linear_axis_stroke_color': 'objects.axis.linear.stroke_color',
        'linear_axis_stroke_width_short': 'objects.axis.linear.stroke_width.short',
        'linear_axis_stroke_width_long': 'objects.axis.linear.stroke_width.long',
        'linear_definition_font_size_short': 'objects.definition.linear.font_size.short',
        'linear_definition_font_size_long': 'objects.definition.linear.font_size.long',
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
        'scale_font_size_short': 'objects.scale.font_size.short',
        'scale_font_size_long': 'objects.scale.font_size.long',
        'scale_interval': 'objects.scale.interval',
        'blast_color_min': 'objects.blast_match.min_color',
        'blast_color_max': 'objects.blast_match.max_color',
        'legend_box_size_short': 'objects.legends.color_rect_size.short',
        'legend_box_size_long': 'objects.legends.color_rect_size.long',
        'legend_font_size_short': 'objects.legends.font_size.short',
        'legend_font_size_long': 'objects.legends.font_size.long',
        'normalize_length': 'canvas.linear.normalize_length'
    }
    for param, path in param_paths.items():
        value = locals()[param]
        if value is not None:
            update_config_value(config_dict, path, value)

    return config_dict

def determine_output_file_prefix(gb_records, output_prefix, record_count, accession):
    if len(gb_records) > 1 and output_prefix is not None:
        return "{}_{}".format(output_prefix, record_count)
    elif len(gb_records) == 1 and output_prefix is not None:
        return output_prefix
    else:
        return accession
    
def check_feature_presence(records: Union[List[SeqRecord], SeqRecord], features_list: List[str]) -> dict:
    if isinstance(records, SeqRecord):
        records = [records]

    features_to_be_checked = set(features_list)
    features_present = []

    for record in records:
        for feature in record.features:
            if feature.type in features_to_be_checked:
                features_present.append(feature.type)
                features_to_be_checked.remove(feature.type)
                if not features_to_be_checked:
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

    if whitelist_map:
        is_eligible = False
        rules = whitelist_map.get(feature_type, {})
        for key, values in qualifiers.items():
            if key in rules and any(v in rules[key] for v in values):
                is_eligible = True
                break
        if not is_eligible:
            return ""

    priority_list = priority_map.get(feature_type, ['product', 'gene', 'locus_tag',
                                                    'protein_id', 'old_locus_tag', 'note'])
    final_label = ""
    for key in priority_list:
        if key in qualifiers and qualifiers[key]:
            final_label = qualifiers[key][0]
            break

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
            start, end = int(coord[2]), int(coord[3])
            length = abs(end - start)
            if length > max_length:
                max_length = length
                longest_segment_info = coord
        except (IndexError, TypeError, ValueError):
            continue
            
    return longest_segment_info, max_length
        

def deg_to_rad(deg):
    return deg * math.pi / 180.0

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

def calculate_angle_degrees(center_x, center_y, x, y, middle, start_angle, end_angle, total_length, x_radius, y_radius, normalize):
    """Calculate the angle in degrees from a point (x, y) relative to the center of an ellipse."""
    x_normalized = (x - center_x) / x_radius
    y_normalized = (y - center_y) / y_radius
    
    angle_radians = math.atan2(y_normalized, x_normalized)
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees

def place_labels_on_arc_fc(labels: list[dict],center_x: float,center_y: float,x_radius: float,y_radius: float,start_angle: float,end_angle: float,total_length: int) -> list[dict]:
    def check_overlap(label1, label2, total_length, margin):
        return y_overlap(label1, label2, total_length, margin) and x_overlap(label1, label2, minimum_margin=2)
    
    if not labels:
        return []

    rearranged_labels = []
    labels = sort_labels(labels)
    current_angle = -75
    increment = 0.1

    for i, label in enumerate(labels):
        if i == 0:
            label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, current_angle)
            rearranged_labels.append(label)
        else:
            new_angle = (current_angle + increment)
            if new_angle <-75: 
                new_angle = -75
            elif -75 <= new_angle <85:
                if label['middle'] > (total_length / 2) or i >= len(labels) * (2/3):
                    new_angle = 85
                else:
                    new_agle = new_angle
            elif 85 < new_angle < 90:
                new_angle = 90
            else:
                new_angle = new_angle
            if new_angle > 269:
                new_angle = 269
            label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, new_angle)
            
            # Note: calculate_coordinates only takes 4 args in the original def above, but place_labels logic 
            # implies simple circular calculation. Using the simplified version defined earlier.
            
            while rearranged_labels and check_overlap(label, rearranged_labels[-1], total_length, 0.1):
                new_angle += 0.01
                label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, new_angle)
                
            rearranged_labels.append(label)
            current_angle = new_angle

    return rearranged_labels

def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def sort_labels(labels):
    return sorted(labels, key=lambda x: x['middle'])

def y_overlap(label1, label2, total_len, minimum_margin):
    # Adjusted to consider absolute values for y coordinates
    label1_start_y = label1["start_y"]
    label2_start_y = label2["start_y"]
    label1_angle = (360.0 * (label1["middle"] / total_len)) % 360
    label2_angle = (360.0 * (label2["middle"] / total_len)) % 360
    
    # Helper to calculate max/min Y based on angle
    def get_y_bounds(lbl, angle, start_y):
        h = lbl["height_px"]
        if lbl.get("is_inner", False) == False:
            if 0 <= angle < 10: return start_y, start_y - h - 0.5 * minimum_margin
            elif 10 <= angle < 170: return start_y + 0.5 * h + 0.5 * minimum_margin, start_y - 0.5 * h - 0.5 * minimum_margin
            elif 170 <= angle < 190: return start_y + h + 0.5 * minimum_margin, start_y - 0.5 * minimum_margin
            elif 190 <= angle < 350: return start_y + 0.5 * h + 0.5 * minimum_margin, start_y - 0.5 * h - 0.5 * minimum_margin
            else: return start_y, start_y - h - 0.5 * minimum_margin
        else:
            if 0 <= angle < 10: return start_y + h + 0.5 * minimum_margin, start_y
            elif 10 <= angle < 170: return start_y + 0.5 * h + 0.5 * minimum_margin, start_y - 0.5 * h - 0.5 * minimum_margin
            elif 170 <= angle < 190: return start_y + 0.5 * h + 0.5 * minimum_margin, start_y - 0.5 * h - 0.5 * minimum_margin
            elif 190 <= angle < 350: return start_y + 0.5 * h + 0.5 * minimum_margin, start_y - 0.5 * h - 0.5 * minimum_margin
            else: return start_y, start_y - h - 0.5 * minimum_margin

    max_y1, min_y1 = get_y_bounds(label1, label1_angle, label1_start_y)
    max_y2, min_y2 = get_y_bounds(label2, label2_angle, label2_start_y)

    if min_y1 < min_y2:
        return max_y1 >= min_y2
    else:
        return max_y2 >= min_y1

def x_overlap(label1, label2, minimum_margin=1.0):
    def get_x_bounds(lbl):
        w = lbl["width_px"]
        sx = lbl["start_x"]
        if lbl.get("is_inner", False) == False:
            if sx > 0: return sx + w + 0.5 * minimum_margin, sx - 0.5 * minimum_margin
            else: return sx + 0.5 * minimum_margin, sx - w - 0.5 * minimum_margin
        else:
            if sx > 0: return sx + 0.5 * minimum_margin, sx - w - 0.5 * minimum_margin
            else: return sx + w + 0.5 * minimum_margin, sx - 0.5 * minimum_margin

    max_x1, min_x1 = get_x_bounds(label1)
    max_x2, min_x2 = get_x_bounds(label2)

    if min_x1 < min_x2:
        return max_x1 >= min_x2 or max_x1 >= max_x2
    else:
        return max_x2 >= min_x1 or max_x2 >= max_x1

def improved_label_placement_fc(labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle, y_margin=0.1, max_iterations=10000):
    def move_label(label, angle):
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def calculate_angle_of_three_points(x1, y1, x2, y2, x3, y3):
        v1 = (x1 - x2, y1 - y2)
        v2 = (x3 - x2, y3 - y2)
        dot_product = v1[0] * v2[0] + v1[1] * v2[1]
        mag1 = math.sqrt(v1[0]**2 + v1[1]**2)
        mag2 = math.sqrt(v2[0]**2 + v2[1]**2)
        try:
            cos_angle = dot_product / (mag1 * mag2)
            angle = math.acos(max(-1, min(1, cos_angle)))
            return math.degrees(angle)
        except ZeroDivisionError:
            return 0

    def check_overlap(label1, label2, total_length):
        return y_overlap(label1, label2, total_length, y_margin) and x_overlap(label1, label2, minimum_margin=1)

    labels = sort_labels(labels)
    for iteration in range(max_iterations):
        changes_made = False
        for i, label in enumerate(reversed(labels)):
            reverse_i = len(labels) -1 - i
            normalize=False
            current_angle = calculate_angle_degrees(center_x, center_y, label['start_x'], label['start_y'], label['middle'], start_angle, end_angle, total_length, x_radius, y_radius, normalize=normalize)

            current_score = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, label['start_x'], label['start_y'])
            
            overlaps_prev = False
            overlaps_next = False

            if i == 0:
                overlaps_prev = check_overlap(label, labels[reverse_i-1], total_length)
                overlaps_next = check_overlap(label, labels[0], total_length)
            elif 0 < i < len(labels) - 1:
                overlaps_prev = check_overlap(labels[reverse_i-1], label, total_length)
                overlaps_next = check_overlap(label, labels[reverse_i+1], total_length)
            elif i == len(labels)- 1 :
                overlaps_prev = check_overlap(label, labels[-1], total_length)
                overlaps_next = check_overlap(label, labels[reverse_i+1], total_length)  
            
            if overlaps_prev and overlaps_next:
                continue

            if overlaps_prev:
                direction = 1 
            elif overlaps_next:
                direction = -1
            else:
                test_angle_plus = (current_angle + 1)
                test_x_plus, test_y_plus = move_label(label, test_angle_plus)
                score_plus = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_plus, test_y_plus)
                
                test_angle_minus = (current_angle - 1)
                test_x_minus, test_y_minus = move_label(label, test_angle_minus)
                score_minus = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_minus, test_y_minus)
                direction = 1 if abs(score_plus) < abs(score_minus) else -1

            creates_new_overlap = False
            while True:
                new_angle = (current_angle + direction * 0.05)
                new_x, new_y = move_label(label, new_angle)
                new_score = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, new_x, new_y)
                
                label_copy = label.copy()
                label_copy['start_x'], label_copy['start_y'] = new_x, new_y               
                if i == 0:
                    if overlaps_prev: creates_new_overlap = (check_overlap(label_copy, labels[0], total_length))
                    elif overlaps_next: creates_new_overlap = (check_overlap(label_copy, labels[reverse_i - 1], total_length))
                    else: creates_new_overlap = (check_overlap(label_copy, labels[reverse_i - 1], total_length)) or (check_overlap(label_copy, labels[0], total_length))
                elif 0 < i < len(labels)-1:
                    prev_label = labels[reverse_i - 1]
                    next_label = labels[reverse_i + 1]
                    if overlaps_prev: creates_new_overlap = (check_overlap(label_copy, next_label, total_length))
                    elif overlaps_next: creates_new_overlap = (check_overlap(label_copy, prev_label, total_length))
                    else: creates_new_overlap = (check_overlap(label_copy, prev_label, total_length)) or (check_overlap(label_copy, next_label, total_length))
                elif i == len(labels) -1:
                    if overlaps_prev: creates_new_overlap = (check_overlap(label_copy, labels[reverse_i + 1], total_length))
                    elif overlaps_next: creates_new_overlap = (check_overlap(label_copy, labels[-1], total_length))
                    else: creates_new_overlap = (check_overlap(label_copy, labels[-1], total_length)) or (check_overlap(label_copy, labels[reverse_i +1], total_length))
                
                if (abs(new_score) <= abs(current_score)) and not creates_new_overlap:
                    label['start_x'], label['start_y'] = new_x, new_y
                    current_angle = new_angle
                    current_score = new_score
                    changes_made = True
                else:
                    break

        if not changes_made:
            break

    return labels


def rearrange_labels_fc(labels, feature_radius, total_length, genome_len, config_dict, strands, is_outer):
    track_type = config_dict['canvas']['circular']['track_type']
    
    if is_outer:
        offset_config = config_dict['labels']['unified_adjustment']['outer_labels']
        x_radius_factor = config_dict['labels']['arc_x_radius_factor'][track_type][strands][genome_len]
        y_radius_factor = config_dict['labels']['arc_y_radius_factor'][track_type][strands][genome_len]
        default_center_x = config_dict['labels']['arc_center_x'][track_type][genome_len]
    else:
        offset_config = config_dict['labels']['unified_adjustment']['inner_labels']
        x_radius_factor = config_dict['labels']['inner_arc_x_radius_factor'][track_type][strands][genome_len]
        y_radius_factor = config_dict['labels']['inner_arc_y_radius_factor'][track_type][strands][genome_len]
        default_center_x = config_dict['labels']['inner_arc_center_x'][track_type][genome_len]

    x_radius_offset = offset_config['x_radius_offset']
    y_radius_offset = offset_config['y_radius_offset']

    if is_outer:
        x_radius = feature_radius * x_radius_factor * x_radius_offset
        y_radius = feature_radius * y_radius_factor * y_radius_offset
    else:
        x_radius = feature_radius * x_radius_factor * (2 - x_radius_offset)
        y_radius = feature_radius * y_radius_factor * (2 - y_radius_offset)

    center_y = 0
    center_x = default_center_x
    start_angle = 0
    end_angle = 360

    labels = sorted(labels, key=lambda x: x['middle'])    
    labels = place_labels_on_arc_fc(labels, center_x, center_y, x_radius, y_radius, start_angle, end_angle, total_length)
    labels = improved_label_placement_fc(labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle)
    
    return labels


def prepare_label_list(feature_dict, total_length, radius, track_ratio, config_dict):
    embedded_labels = []
    outer_labels = []
    inner_labels = []
    
    label_filtering = config_dict['labels']['filtering']
    label_filtering = preprocess_label_filtering(label_filtering)
    length_threshold = config_dict['labels']['length_threshold']['circular']
    length_param = determine_length_parameter(total_length, length_threshold)
    track_type = config_dict['canvas']['circular']['track_type']
    strandedness = config_dict['canvas']['strandedness']
    
    strands = "separate" if strandedness else "single"
    allow_inner_labels = config_dict['canvas']['circular']['allow_inner_labels']
    radius_factor = config_dict['labels']['radius_factor'][track_type][strands][length_param] 
    inner_radius_factor = config_dict['labels']['inner_radius_factor'][track_type][strands][length_param]  

    font_family = config_dict['objects']['text']['font_family']
    font_size: str = config_dict['labels']['font_size'][length_param]
    interval = config_dict['canvas']['dpi']

    track_ratio_factor = config_dict['canvas']['circular']['track_ratio_factors'][length_param][0]
    
    # Pre-fetch for bbox calculation
    from .circular_path_drawer import calculate_feature_position_factors_circular
    from .create_feature_objects import get_strand
    
    for feature_object in feature_dict.values():
        feature_label_text = get_label_text(feature_object, label_filtering)
        if feature_label_text == '':
            continue
        else:      
            label_entry = dict()
            longest_segment_length = 0
            is_embedded = False
            label_middle = 0
            coordinate_strand: str = "undefined"
            feature_location_list = feature_object.location
            list_of_coordinates = feature_object.coordinates
            feature_object.strand = feature_location_list[0][2]

            feature_location_count=0
            longest_segment_start = 0
            longest_segment_end = 0
            longeset_segment_middle = 0

            for coordinate in list_of_coordinates:
                if feature_location_list[feature_location_count][0] == "line":
                    feature_location_count += 1
                    continue
                else:
                    coordinate_start = int(coordinate.start)
                    coordinate_end = int(coordinate.end)
                    coordinate_strand = get_strand(coordinate.strand)
                    interval_length = abs(int(coordinate_end - coordinate_start ) + 1)
                    interval_middle = int(coordinate_end + coordinate_start) / 2
                    feature_location_count += 1
                    if interval_length > longest_segment_length:
                        longest_segment_start = coordinate_start
                        longest_segment_end = coordinate_end
                        longeset_segment_middle = interval_middle
                        longest_segment_length = interval_length
            
            cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
            factors: list[float] = calculate_feature_position_factors_circular(total_length, coordinate_strand, track_ratio, cds_ratio, offset, track_type, strandedness)
            
            bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
            
            label_middle = longeset_segment_middle
            label_as_feature_length = total_length * (1.1 * bbox_width_px)/(2*math.pi*radius)
            label_start = label_middle - (label_as_feature_length/2)
            label_end = label_middle + (label_as_feature_length/2)
            feature_middle_x: float = (radius * factors[1]) * math.cos(math.radians(360.0 * ((label_middle) / total_length) - 90))
            feature_middle_y: float = (radius * factors[1]) * math.sin(math.radians(360.0 * ((label_middle) / total_length) - 90))
            
            if feature_object.strand == "positive" or allow_inner_labels == False:
                middle_x: float = (radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y: float = (radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            else:
                middle_x: float = (inner_radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y: float = (inner_radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            
            if label_start > longest_segment_start and label_end < longest_segment_end: 
                is_embedded = True
            else:
                is_embedded = False
                
            label_entry["label_text"] = feature_label_text
            label_entry["middle"] = label_middle
            label_entry["start"] = label_start
            label_entry["end"] = label_end
            label_entry["middle_x"] = middle_x
            label_entry["middle_y"] = middle_y
            label_entry["feature_middle_x"] = feature_middle_x
            label_entry["feature_middle_y"] = feature_middle_y                        
            label_entry["width_px"] = bbox_width_px
            label_entry["height_px"] = bbox_height_px
            label_entry["strand"] = coordinate_strand
            label_entry["is_embedded"] = is_embedded
            label_entry["font_size"] = font_size
            label_entry["font_family"] = font_family

            if is_embedded:
                embedded_labels.append(label_entry)
            else:
                if label_entry["middle"] > (total_length / 2):
                    if feature_object.strand == "positive" or allow_inner_labels == False:
                        label_entry["is_inner"] = False
                        outer_labels.append(label_entry)
                    else:
                        label_entry["is_inner"] = True
                        inner_labels.append(label_entry)
                else:
                    if feature_object.strand == "positive" or allow_inner_labels == False:
                        label_entry["is_inner"] = False
                        outer_labels.append(label_entry)
                    else:
                        label_entry["is_inner"] = True
                        inner_labels.append(label_entry)

    outer_labels_rearranged = rearrange_labels_fc(outer_labels, radius, total_length, length_param, config_dict, strands, is_outer=True)
    inner_labels_rearranged = []
    if allow_inner_labels:
        inner_labels_rearranged = rearrange_labels_fc(inner_labels, radius, total_length, length_param, config_dict, strands, is_outer=False)
    
    label_list_fc = embedded_labels + outer_labels_rearranged + inner_labels_rearranged
    return label_list_fc

def check_label_overlap(label1, label2):
   return not (label1["end"] < label2["start"] or label2["end"] < label1["start"])

def find_lowest_available_track(track_dict, label):
   track_num = 1
   while True:
       track_id = f"track_{track_num}"
       has_overlap = False
       if track_id in track_dict:
           for existing_label in track_dict[track_id]:
               if check_label_overlap(label, existing_label):
                   has_overlap = True
                   break
       
       if not has_overlap:
           return track_num
       track_num += 1

def prepare_label_list_linear(feature_dict, genome_length, alignment_width, 
                           genome_size_normalization_factor, cds_height, 
                           strandedness, config_dict):
   
   # Pre-fetch needed functions
   from .linear_path_drawer import calculate_feature_position_factors_linear
   from .create_feature_objects import get_strand
   from collections import defaultdict
   
   embedded_labels = []
   external_labels = []
   track_dict = defaultdict(list)
   feature_track_positions = {} 
   
   length_threshold = config_dict['labels']['length_threshold']['circular']
   length_param = determine_length_parameter(genome_length, length_threshold)
   font_family = config_dict['objects']['text']['font_family']
   font_size = config_dict['labels']['font_size']['linear'][length_param]
   interval = config_dict['canvas']['dpi']
   label_filtering = config_dict['labels']['filtering']

   max_feature_track = 0
   for feature_object in feature_dict.values():
       if feature_object.feature_track_id > max_feature_track:
           max_feature_track = feature_object.feature_track_id
   
   top_factors = calculate_feature_position_factors_linear(
       strand='positive', 
       track_id=max_feature_track, 
       separate_strands=strandedness
   )
   top_feature_y_limit = cds_height * top_factors[0]

   for feature_id, feature_object in feature_dict.items():
       if len(feature_object.coordinates) == 0:
           continue
       coordinate = feature_object.coordinates[0]
       strand = get_strand(coordinate.strand)
       feature_track_id = feature_object.feature_track_id
       factors = calculate_feature_position_factors_linear(strand, feature_track_id, strandedness)
       track_y_position = cds_height * factors[1]
       feature_track_positions[feature_id] = track_y_position
   
   max_bbox_height = 0
   for feature_id, feature_object in reversed(list(feature_dict.items())):
       feature_label_text = get_label_text(feature_object, label_filtering)
       feature_track_id = feature_object.feature_track_id
       if not feature_label_text:
           continue
           
       bbox_width_px, bbox_height_px = calculate_bbox_dimensions(
           feature_label_text, font_family, font_size, interval)
       if bbox_height_px > max_bbox_height:
            max_bbox_height = bbox_height_px

       longest_segment_length = 0
       label_middle = 0
       coordinate_strand = None
       factors = None
       longest_segment_start = 0
       longest_segment_end = 0
       
       for coordinate in feature_object.coordinates:
           coordinate_strand = get_strand(coordinate.strand)
           factors = calculate_feature_position_factors_linear(coordinate_strand, feature_track_id, strandedness)
           start = int(coordinate.start)
           end = int(coordinate.end)
           segment_length = abs(end - start + 1)
           if segment_length > longest_segment_length:
               longest_segment_start = start
               longest_segment_end = end
               longest_segment_length = segment_length
               label_middle = (end + start) / 2
               
       normalized_start = normalize_position_to_linear_track(
           longest_segment_start, genome_length, alignment_width, genome_size_normalization_factor)
       normalized_end = normalize_position_to_linear_track(
           longest_segment_end, genome_length, alignment_width, genome_size_normalization_factor)
       longest_segment_length_in_pixels = abs(normalized_end - normalized_start) + 1
       normalized_middle = (normalized_start + normalized_end) / 2
       bbox_start = normalized_middle - (bbox_width_px / 2)
       bbox_end = normalized_middle + (bbox_width_px / 2)
       
       feature_y = feature_track_positions.get(feature_id, cds_height * factors[1])
       
       label_entry = {
           "label_text": feature_label_text,
           "middle": normalized_middle,
           "start": bbox_start,
           "end": bbox_end,
           "middle_x": normalized_middle,
           "width_px": bbox_width_px,
           "height_px": bbox_height_px,
           "strand": coordinate_strand,
           "feature_middle_y": feature_y,
           "font_size": font_size,
           "font_family": font_family
       }

       if bbox_width_px < longest_segment_length_in_pixels:
           label_entry.update({
               "middle_y": feature_y,
               "is_embedded": True,
               "track_id": "track_0"
           })
           track_dict["track_0"].append(label_entry)
       else:
           label_entry.update({
               "middle_y": 0,
               "is_embedded": False
           })
           best_track = find_lowest_available_track(track_dict, label_entry)
           label_entry["track_id"] = f"track_{best_track}"
           track_dict[f"track_{best_track}"].append(label_entry)

   if "track_0" in track_dict:
       for label in track_dict["track_0"]:
           embedded_labels.append(label)
   
   track_height = max_bbox_height * 1.1 
   for track_id in sorted(track_dict.keys()):
       if track_id == "track_0":
           continue
       track_labels = track_dict[track_id]
       track_num = int(track_id.split("_")[1])
       for label in track_labels:
           label["middle_y"] = top_feature_y_limit - (track_height * track_num)
           external_labels.append(label)

   return embedded_labels + external_labels