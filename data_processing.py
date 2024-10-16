#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import math
from collections import defaultdict
from Bio.SeqFeature import SimpleLocation
from pandas import DataFrame
from typing import Generator, Any, Dict, List, Optional
from Bio.SeqRecord import SeqRecord
from .create_feature_objects import get_strand
from .utility_functions import calculate_bbox_dimensions, get_label_text

def skew_df(record: SeqRecord, window: int, step: int, nt: str) -> DataFrame:
    """
    Calculates dinucleotide skew and content in a DNA sequence, returning a DataFrame.

    Args:
        record (SeqRecord): BioPython SeqRecord object containing the DNA sequence.
        window (int): Window size for calculating skew and content.
        step (int): Step size for the sliding window.
        nt (str): Dinucleotide pair for skew calculation (e.g., 'GC').

    Returns:
        DataFrame: A pandas DataFrame containing columns for dinucleotide content,
                   skew, and cumulative skew, each indexed by sequence position.

    This function calculates the skew and content of a specified dinucleotide pair across
    a DNA sequence. The skew is calculated for each window of the sequence, and both the
    cumulative and individual skews are normalized and returned along with the content in a DataFrame.
    """
    nt_list = list(nt)
    nt_1: str = nt_list[0]
    nt_2: str = nt_list[1]
    skew_sum: float = 0
    skew_dict: dict[int, float] = {}
    content_dict: dict[int, float] = {}
    skew_cumulative_dict: dict[int, float] = {}
    seq: str = record.seq.upper()
    for start, seq_part in sliding_window(seq, window, step):
        skew: float = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        dinucleotide_content: float = (seq_part.count(
            nt_1) + seq_part.count(nt_2)) / len(seq_part)
        content_dict[start] = dinucleotide_content
        skew_dict[start] = skew
        skew_sum = (skew_sum + skew)
        skew_cumulative_dict[start] = (skew_sum)
    max_skew_abs: float = abs(max(skew_dict.values(), key=abs))
    max_skew_cumulative_abs: float = abs(
        max(skew_cumulative_dict.values(), key=abs))
    factor: float = (max_skew_abs / max_skew_cumulative_abs)
    skew_cumulative_dict.update((x, (y * factor))
                                for x, y in skew_cumulative_dict.items())
    content_legend: str = "{} content".format(nt)
    skew_legend: str = "{} skew".format(nt)
    cumulative_skew_legend: str = "Cumulative {} skew, normalized".format(nt)
    df = pd.DataFrame({content_legend: pd.Series(content_dict), skew_legend: pd.Series(
        skew_dict), cumulative_skew_legend: pd.Series(skew_cumulative_dict)})
    return df


def calculate_dinucleotide_skew(seq: str, base1: str, base2: str) -> float:
    """
    Calculate the skew of two nucleotides in a DNA sequence.

    The skew is calculated as (Count(Base1) - Count(Base2)) / (Count(Base1) + Count(Base2)).
    This metric helps identify biases in the occurrence of specific nucleotides.

    Args:
        seq (str): DNA sequence to analyze.
        base1 (str): First nucleotide for skew calculation.
        base2 (str): Second nucleotide for skew calculation.

    Returns:
        float: Skew value for the two nucleotides in the given sequence.
    """
    base1_count: int = seq.count(base1)
    base2_count: int = seq.count(base2)
    skew: float = (base1_count - base2_count) / (base1_count + base2_count)
    return skew


def calculate_gc_percent(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        float: The GC content as a percentage of the total sequence length.
    """
    g_count: int = sequence.count("G")
    c_count: int = sequence.count("C")
    gc_percent: float = round(100 * ((g_count + c_count) / len(sequence)), 2)
    return gc_percent


def sliding_window(seq: str, window: int, step: int) -> Generator[tuple[int, str], Any, None]:
    """
    Generate a sequence of substrings from a given sequence using a sliding window approach.

    This generator function iterates over a sequence and yields substrings of a specified length (window),
    moving a certain number of bases each time (step). It is designed to handle circular sequences, so when
    the window extends beyond the end of the sequence, it wraps around to the beginning.

    Args:
        seq (str): The sequence to be processed.
        window (int): The length of the window (substring) to be extracted.
        step (int): The number of bases to move the window each time.

    Yields:
        tuple: A tuple containing the start position of the window and the extracted substring.

    For each iteration, the function calculates the end position of the window. If the end extends beyond the
    sequence length, it wraps around, combining the end of the sequence with the beginning to form a circular
    sequence. This approach is particularly useful for genomes with circular topology.
    """
    for start in range(0, len(seq), step):
        end: int = start + window
        if end > len(seq):
            overhang_length: int = (end - len(seq))
            # assuming circular sequence
            out_seq: str = seq[start:len(seq)] + seq[0:overhang_length]
        else:
            out_seq = seq[start:end]
        yield start, out_seq

def prepare_legend_table(gc_config, skew_config, feature_config, features_present):
    legend_table = dict()
    color_table: Optional[DataFrame] = feature_config.color_table
    default_colors: DataFrame = feature_config.default_colors
    features_present: List[str] = features_present
    block_stroke_color: str = feature_config.block_stroke_color
    block_stroke_width: float = feature_config.block_stroke_width
    show_gc = gc_config.show_gc
    gc_stroke_color: str = gc_config.stroke_color
    gc_stroke_width: float = gc_config.stroke_width
    gc_high_fill_color: str = gc_config.high_fill_color
    gc_low_fill_color: str = gc_config.low_fill_color
    show_skew = skew_config.show_skew
    skew_high_fill_color: str = skew_config.high_fill_color
    skew_low_fill_color: str = skew_config.low_fill_color
    skew_stroke_color: str = skew_config.stroke_color
    skew_stroke_width: float = skew_config.stroke_width
    feature_specific_colors = dict()
    if color_table is not None and not color_table.empty:
        for _, row in color_table.iterrows():
            feature_type = row['feature_type']
            if feature_type not in feature_specific_colors:
                feature_specific_colors[feature_type] = []
            feature_specific_colors[feature_type].append((row['caption'], row['color']))
    for selected_feature in features_present:
        if selected_feature in feature_specific_colors.keys():
            for entry in feature_specific_colors[selected_feature]:
                specific_caption = entry[0]
                specific_fill_color = entry[1]
                legend_table[specific_caption] = (block_stroke_color, block_stroke_width, specific_fill_color)
            if selected_feature == 'CDS':
                new_selected_key_name = 'other proteins'
            else:
                new_selected_key_name = f'other {selected_feature}s'
            feature_fill_color = default_colors[default_colors['feature_type'] == selected_feature]['color'].values[0]
            legend_table[new_selected_key_name] = (block_stroke_color, block_stroke_width, feature_fill_color)
        else:
            feature_fill_color = default_colors[default_colors['feature_type'] == selected_feature]['color'].values[0]
            legend_table[selected_feature] = (block_stroke_color, block_stroke_width, feature_fill_color)        
    if show_gc:
        if gc_high_fill_color == gc_low_fill_color:
            legend_table['GC content'] = (gc_stroke_color, gc_stroke_width, gc_high_fill_color)
        else:
            legend_table['GC content (+)'] = (gc_stroke_color, gc_stroke_width, gc_high_fill_color)
            legend_table['GC content (-)'] = (gc_stroke_color, gc_stroke_width, gc_low_fill_color)
    if show_skew:
        if skew_high_fill_color == skew_low_fill_color:
            legend_table['GC skew'] = (skew_stroke_color, skew_stroke_width, skew_high_fill_color)
        else:
            legend_table['GC skew (+)'] = (skew_stroke_color, skew_stroke_width, skew_high_fill_color)
            legend_table['GC dkew (-)'] = (skew_stroke_color, skew_stroke_width, skew_low_fill_color)
    return legend_table



def y_overlap(label1, label2, minimum_margin):
    # Adjusted to consider absolute values for y coordinates
    label1_start_y = label1["start_y"]
    label2_start_y = label2["start_y"]
    if label1_start_y < label2_start_y:
        if label1_start_y + 0.5 * label1["height_px"] + minimum_margin > (label2_start_y - 0.5 * label2["height_px"]):
            return True
        else:
            return False
    else:
        if (label2_start_y + 0.5 * label2["height_px"] + minimum_margin) > (label1_start_y - 0.5 * label1["height_px"] ):
            return True
        else:
            return False

def x_overlap(label1, label2):

    # Adjusted to directly return the evaluated condition
    min_x1 = label1["start_x"]
    max_x1 = label1["start_x"] + label1["width_px"]
    min_x2 = label2["start_x"]
    max_x2 = label2["start_x"] + label2["width_px"]
    if min_x1 < min_x2:
        if (max_x1 - min_x2) >0:
            return True
        else:
            return False
    else:
        if (max_x2 - min_x1) >0:
            return True
        else:
            return False
        

def calculate_angle_degrees(center_x, center_y, x, y, middle, start_angle, end_angle, total_length):
    """Calculate the angle in degrees from a point (x, y) relative to a center (center_x, center_y)."""
    angle_radians = math.atan2(y - center_y, x - center_x)
    angle_degrees = math.degrees(angle_radians)
    if x >= 0:
        angle_degrees = angle_degrees
    else:
        if y >= 0:
            angle_degrees = angle_degrees
        else:
            angle_degrees = 360 + angle_degrees
    return angle_degrees

def calculate_coordinates(center_x, center_y, radius, angle_degrees, middle, total_length):
    angle_radians = math.radians(angle_degrees)
    y = center_y + radius * math.sin(angle_radians)
    x = center_x + radius * math.cos(angle_radians)
    
    # Adjust x-coordinate for left side labels
    if middle >= total_length/2:
        x = center_x - abs(radius * math.cos(angle_radians))
    
    return x, y

def calculate_angle_for_y(center_y, radius, y):
    """Calculate the angle in degrees for a given y-coordinate on a circle."""
    # Ensure the y-coordinate is within the circle's bounds to avoid math domain errors
    if center_y - radius <= y <= center_y + radius:
        # Calculate the arcsine of the y-coordinate
        angle_radians = math.asin((y - center_y) / radius)
        angle_degrees = math.degrees(angle_radians)
        return angle_degrees
    else:
        return None  # Indicates the y-coordinate is outside the circle's bounds

def place_labels_on_arc(labels, center_x, center_y, radius, start_angle, end_angle, total_length):
    start_y = center_y + radius * math.sin(math.radians(start_angle))
    end_y = center_y + radius * math.sin(math.radians(end_angle))
    total_y_range = end_y - start_y
    y_increment = total_y_range / (len(labels) - 1)

    for i, label in enumerate(labels):
        y = start_y + i * y_increment
        angle_degrees = calculate_angle_for_y(center_y, radius, y)
        # print(label["middle"], angle_degrees)
        if angle_degrees is not None:
            label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, radius, angle_degrees, label['middle'], total_length)
            start_angle = 0
            end_angle = 0
            degree = calculate_angle_degrees(center_x, center_y, label['start_x'], label['start_y'], label["middle"], start_angle, end_angle, total_length)
            # print(label['middle'], label['start_x'], label['start_y'], degree)
        else:
            print(f"Y-coordinate {y} is outside the circle's bounds for label {label}.")

    return labels

def get_shifted_label_and_distance(label, current_angle, step, center_x, center_y, arc_radius, total_length):
    x1 = label["middle_x"]
    y1 = label["middle_y"]
    new_label = label
    new_angle = current_angle + step
    new_start_x, new_start_y = calculate_coordinates(center_x, center_y, arc_radius, new_angle, label['middle'], total_length)
    new_label["start_x"] = new_start_x
    new_label["start_y"] = new_start_y
    new_distance = math.sqrt((new_start_x - x1) ** 2 + (new_start_y - y1) ** 2)
    return new_label, new_distance
def sort_labels(labels):
    return sorted(labels, key=lambda x: x['middle'])

def shift_label_group(labels, start_index, end_index, shift_amount, center_x, center_y, radius, total_length):
    for i in range(start_index, end_index + 1):
        label = labels[i]
        current_angle = calculate_angle_degrees(center_x, center_y, label["start_x"], label["start_y"], label["middle"], 0, 360, total_length)
        new_angle = current_angle + shift_amount
        label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, radius, new_angle, label['middle'], total_length)
    return labels

def iteratively_refine_labels(labels, center_x, center_y, arc_radius, minimum_margin, start_angle, end_angle, total_length, num_iterations=1):
    labels = sort_labels(labels)
    for _ in range(num_iterations):
        i = 0
        while i < len(labels) - 1:
            current_label = labels[i]
            next_label = labels[i + 1]
            if y_overlap(current_label, next_label, minimum_margin):
                # Find the group of overlapping labels
                group_end = i + 1
                while group_end < len(labels) - 1 and y_overlap(labels[group_end], labels[group_end + 1], minimum_margin):
                    group_end += 1
                
                # Calculate the total shift needed
                total_shift = sum(label["height_px"] for label in labels[i:group_end+1]) + minimum_margin * (group_end - i)
                shift_per_label = total_shift / (group_end - i + 1)
                
                # Shift the group of labels
                shift_amount = math.degrees(shift_per_label / arc_radius)
                labels = shift_label_group(labels, i, group_end, shift_amount, center_x, center_y, arc_radius, total_length)
                
                i = group_end + 1
            else:
                i += 1
    return labels

def shift_label_positions(labels, center_x, center_y, arc_radius, start_angle, end_angle, total_length):
    minimum_margin = 5
    start_y = center_y + arc_radius * math.sin(math.radians(start_angle))
    end_y = center_y + arc_radius * math.sin(math.radians(end_angle))
    arc_height = abs(end_y - start_y)
    total_bbox_height = sum(label["height_px"] for label in labels)
    total_margin = minimum_margin * (len(labels)-1)
    total_label_height = total_bbox_height + total_margin
    if arc_height >= total_label_height:
        labels = iteratively_refine_labels(labels, center_x, center_y, arc_radius, minimum_margin, start_angle, end_angle, total_length, num_iterations=1000)
    else:
        print("Labels are placed too tight")
    return labels

def rearrange_labels(labels, radius, total_length, is_right):
    arc_radius = radius * 1.75
    center_y = 0
    if is_right:
        center_x = -220
        start_angle = -40
        end_angle = 40
    else:
        center_x = 220
        start_angle = 140
        end_angle = 220       
    labels = place_labels_on_arc(labels, center_x, center_y, arc_radius, start_angle, end_angle, total_length)
    labels = shift_label_positions(labels, center_x, center_y, arc_radius, start_angle, end_angle, total_length)
    return labels

def prepare_label_list(feature_dict, total_length, radius, font, fontsize, interval):
    embedded_labels = []
    left_labels = []
    right_labels = []
    label_list = []

    for feature_object in feature_dict.values():
        feature_label_text = get_label_text(feature_object)
        if feature_label_text == '':
            continue
        else:      
            label_entry = dict()
            longest_segment_length = 0
            is_embedded = False
            label_middle = 0
            coordinate_strand: str = "undefined"
            list_of_coordinates: List[SimpleLocation] = feature_object.coordinates
            for coordinate in list_of_coordinates:
                coordinate_start = int(coordinate.start)
                coordinate_end = int(coordinate.end)
                coordinate_strand = get_strand(coordinate.strand)
                interval_length = int(coordinate_end - coordinate_start + 1)
                interval_middle = int(coordinate_end + coordinate_start) / 2
                if interval_length > longest_segment_length:
                    longest_segment_length = interval_length
                    label_middle = interval_middle
            longest_segment_lengh_in_pixels = (2*math.pi*radius) * (longest_segment_length)/total_length
            bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font, fontsize, interval)
            label_as_feature_length = total_length * bbox_width_px/(2*math.pi*radius)
            label_start = label_middle - (label_as_feature_length/2)
            label_end = label_middle + (label_as_feature_length/2)
            middle_x: float = (1.05 * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
            middle_y: float = (1.05 * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))  
            if bbox_width_px < longest_segment_lengh_in_pixels:
                is_embedded = True
            else:
                is_embedded = False
            label_entry["label_text"] = feature_label_text
            label_entry["middle"] = label_middle
            label_entry["start"] = label_start
            label_entry["end"] = label_end
            label_entry["middle_x"] = middle_x
            label_entry["middle_y"] = middle_y            
            label_entry["width_px"] = bbox_width_px
            label_entry["height_px"] = bbox_height_px
            label_entry["strand"] = coordinate_strand
            label_entry["is_embedded"] = is_embedded
            if is_embedded:
                embedded_labels.append(label_entry)
            else:
                if label_entry["middle"] > (total_length / 2):
                    left_labels.append(label_entry)
                else:
                    right_labels.append(label_entry)
    right_labels_rearranged = rearrange_labels(right_labels, radius, total_length, is_right=True)
    left_labels_rearranged = rearrange_labels(left_labels, radius, total_length, is_right=False)
    label_list = embedded_labels + right_labels_rearranged + left_labels_rearranged
    return label_list