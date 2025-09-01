#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import random
import math
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from Bio.SeqFeature import SimpleLocation
from pandas import DataFrame
from typing import Generator, Any, Dict, List, Optional
from Bio.SeqRecord import SeqRecord
from .create_feature_objects import get_strand
from .utility_functions import calculate_bbox_dimensions, get_label_text , determine_length_parameter, calculate_cds_ratio
from .circular_path_drawer import calculate_feature_position_factors_circular
from .linear_path_drawer import calculate_feature_position_factors_linear, normalize_position_to_linear_track, generate_text_path



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
            matching_rows = default_colors[default_colors['feature_type'] == selected_feature]
            if not matching_rows.empty:
                feature_fill_color = default_colors[default_colors['feature_type'] == selected_feature]['color'].values[0]
            else:
               feature_fill_color = default_colors[default_colors['feature_type'] == "default"]['color'].values[0]
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
            legend_table['GC skew (-)'] = (skew_stroke_color, skew_stroke_width, skew_low_fill_color)
    return legend_table

def y_overlap(label1, label2, total_len, minimum_margin):
    # Adjusted to consider absolute values for y coordinates
    label1_start_y = label1["start_y"]
    label2_start_y = label2["start_y"]
    label1_angle = (360.0 * (label1["middle"] / total_len)) % 360
    label2_angle = (360.0 * (label2["middle"] / total_len)) % 360
    if label1["is_inner"] == False:
        if   0   <= label1_angle < 10:    # baseline_value = "text-after-edge"
            max_y1 = label1_start_y 
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin 
        elif 10  <= label1_angle < 170:   # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 170 <= label1_angle < 190:   # baseline_value = "hanging"
            max_y1 = label1_start_y + 1.0 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * minimum_margin 
        elif 190 <= label1_angle < 350:   # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        else:                      # baseline_value = "text-after-edge"
            max_y1 = label1_start_y 
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin
    else:
        if   0   <= label1_angle < 10:    # baseline_value = "hanging"
            max_y1 = label1_start_y + 1.0 * label1["height_px"] + 0.5 * minimum_margin 
            min_y1 = label1_start_y 
        elif 10  <= label1_angle < 170:   # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 170 <= label1_angle < 190:   # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 190 <= label1_angle < 350:   # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        else:                      # baseline_value = "hanging"
            max_y1 = label1_start_y 
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin
    if label2["is_inner"] == False:
        if   0   <= label2_angle < 10:    # baseline_value = "text-after-edge"
            max_y2 = label2_start_y 
            min_y2 = label2_start_y - 1.0 * label2["height_px"] - 0.5 * minimum_margin
        elif 10  <= label2_angle < 170:   # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 170 <= label2_angle < 190:   # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * minimum_margin
        elif 190 <= label2_angle < 350:   # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        else:                      # baseline_value = "text-after-edge"
            max_y2 = label2_start_y
            min_y2 = label2_start_y - 1.0 * label2["height_px"] - 0.5 * minimum_margin
    else:
        if   0   <= label2_angle < 10:    # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y
        elif 10  <= label2_angle < 170:   # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 170 <= label2_angle < 190:   # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 190 <= label2_angle < 350:   # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        else:                      # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y
    if min_y1 < min_y2:
        if max_y1 >= min_y2:
            return True
        else:
            return False
    else:
        if max_y2 >= min_y1:
            return True
        else:
            return False

            
def x_overlap(label1, label2, minimum_margin=1.0):
    # Adjusted to directly return the evaluated condition

    if label1["is_inner"] == False:
        if label1["start_x"] > 0:
            max_x1 = label1["start_x"] + label1["width_px"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - 0.5 * minimum_margin
        else:
            max_x1 = label1["start_x"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - label1["width_px"] - 0.5 * minimum_margin
    else:
        if label1["start_x"] > 0:
            max_x1 = label1["start_x"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - label1["width_px"] - 0.5 * minimum_margin
        else:
            max_x1 = label1["start_x"] + label1["width_px"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - 0.5 * minimum_margin
    if label2["is_inner"] == False:
        if label2["start_x"] > 0:
            max_x2 = label2["start_x"] + label2["width_px"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - 0.5 * minimum_margin
        else:
            max_x2 = label2["start_x"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - label2["width_px"] - 0.5 * minimum_margin
    else:
        if label2["start_x"] > 0:
            max_x2 = label2["start_x"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - label2["width_px"] - 0.5 * minimum_margin
        else:
            max_x2 = label2["start_x"] + label2["width_px"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - 0.5 * minimum_margin
    if min_x1 < min_x2:
        if max_x1 >= min_x2:
            return True
        # when nested also true
        elif max_x1 >= max_x2:
            return True
        else:
            return False
    else:
        if max_x2 >= min_x1:
            return True
        # when nested also true
        elif max_x2 >= max_x1:
            return True
        else:
            return False
        

def calculate_angle_degrees(center_x, center_y, x, y, middle, start_angle, end_angle, total_length, x_radius, y_radius, normalize):
    """Calculate the angle in degrees from a point (x, y) relative to the center of an ellipse."""
    # Normalize coordinates to a unit circle
    x_normalized = (x - center_x) / x_radius
    y_normalized = (y - center_y) / y_radius
    
    angle_radians = math.atan2(y_normalized, x_normalized)
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees

def calculate_coordinates(center_x, center_y, x_radius, y_radius, angle_degrees, middle, total_length):
    angle_radians = math.radians(angle_degrees)
    y = center_y + y_radius * math.sin(angle_radians)
    x = center_x + x_radius * math.cos(angle_radians)

    return x, y

def calculate_angle_for_y(center_y, y_radius, y):
    """Calculate the angle in degrees for a given y-coordinate on an ellipse."""
    if center_y - y_radius <= y <= center_y + y_radius:
        # Calculate the arcsine of the normalized y-coordinate
        angle_radians = math.asin((y - center_y) / y_radius)
        angle_degrees = math.degrees(angle_radians)
        return angle_degrees
    else:
        return None  # Indicates the y-coordinate is outside the ellipse's bounds


def place_labels_on_arc_fc(labels: list[dict],center_x: float,center_y: float,x_radius: float,y_radius: float,start_angle: float,end_angle: float,total_length: int) -> list[dict]:
    def calculate_angle(x, y, origin_x, origin_y):
        return math.degrees(math.atan2((y - origin_y), (x - origin_x))) % 360 

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
        cos_angle = dot_product / (mag1 * mag2)
        angle = math.acos(max(-1, min(1, cos_angle)))
        return math.degrees(angle)

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
            label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, current_angle, label['middle'], total_length)
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
            label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, new_angle, label['middle'], total_length)
            # if the new position overlaps with the previous label, adjust the angle
            while rearranged_labels and check_overlap(label, rearranged_labels[-1], total_length, 0.1):
                new_angle += 0.01
                label['start_x'], label['start_y'] = calculate_coordinates(center_x, center_y, x_radius, y_radius, new_angle, label['middle'], total_length)
                
            rearranged_labels.append(label)
            current_angle = new_angle

    return rearranged_labels

def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def sort_labels(labels):
    return sorted(labels, key=lambda x: x['middle'])

def improved_label_placement_fc(labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle, y_margin=0.1, max_iterations=10000):
    def calculate_angle(x, y, origin_x, origin_y):
        return math.degrees(math.atan2((y - origin_y), (x - origin_x))) % 360 

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
        cos_angle = dot_product / (mag1 * mag2)
        angle = math.acos(max(-1, min(1, cos_angle)))
        return math.degrees(angle)

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
            # Check overlaps with neighbors
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
            # Determine movement direction
            if overlaps_prev:
                direction = 1  # Move towards end angle
            elif overlaps_next:
                direction = -1  # Move towards start angle
            else:
                # If no overlap, determine direction based on which way reduces the score
                test_angle_plus = (current_angle + 1)
                test_x_plus, test_y_plus = move_label(label, test_angle_plus)
                score_plus = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_plus, test_y_plus)
                
                test_angle_minus = (current_angle - 1)
                test_x_minus, test_y_minus = move_label(label, test_angle_minus)
                score_minus = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_minus, test_y_minus)
                direction = 1 if abs(score_plus) < abs(score_minus) else -1

            # Move label
            creates_new_overlap = False
            while True:
                new_angle = (current_angle + direction * 0.05)
                new_x, new_y = move_label(label, new_angle)
                new_score = calculate_angle_of_three_points(label["feature_middle_x"], label["feature_middle_y"], 0, 0, new_x, new_y)
                # Check if this move would create overlap with neighbors
                label_copy = label.copy()
                label_copy['start_x'], label_copy['start_y'] = new_x, new_y               
                if i == 0:
                    if overlaps_prev:
                        creates_new_overlap = (check_overlap(label_copy, labels[0], total_length))
                    elif overlaps_next:
                        creates_new_overlap = (check_overlap(label_copy, labels[reverse_i - 1], total_length))
                    else: 
                        creates_new_overlap = (check_overlap(label_copy, labels[reverse_i - 1], total_length)) or (check_overlap(label_copy, labels[0], total_length))
                elif 0 < i < len(labels)-1:
                    prev_label = labels[reverse_i - 1]
                    next_label = labels[reverse_i + 1]
                    if overlaps_prev:
                        creates_new_overlap = (check_overlap(label_copy, next_label, total_length))
                    elif overlaps_next:
                        creates_new_overlap = (check_overlap(label_copy, prev_label, total_length))
                    else:
                        creates_new_overlap = (check_overlap(label_copy, prev_label, total_length)) or (check_overlap(label_copy, next_label, total_length))
                elif i == len(labels) -1:
                    if overlaps_prev:
                        creates_new_overlap = (check_overlap(label_copy, labels[reverse_i + 1], total_length))
                    elif overlaps_next:
                        creates_new_overlap = (check_overlap(label_copy, labels[-1], total_length))
                    else:
                        creates_new_overlap = (check_overlap(label_copy, labels[-1], total_length)) or (check_overlap(label_copy, labels[reverse_i +1], total_length))
                if (abs(new_score) <= abs(current_score)) and not creates_new_overlap:
                    label['start_x'], label['start_y'] = new_x, new_y
                    current_angle = new_angle
                    current_score = new_score
                    changes_made = True
                else:
                    break  # Can't move further without issues or no improvement
        if not changes_made:
            break  # No changes were made in this iteration, so we can stop

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
    # Initial placement of labels
    labels = place_labels_on_arc_fc(labels, center_x, center_y, x_radius, y_radius, start_angle, end_angle, total_length)
    
    # Apply improved label placement
    labels = improved_label_placement_fc(labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle)
    
    return labels


def prepare_label_list(feature_dict, total_length, radius, track_ratio, config_dict):
    embedded_labels = []
    left_labels = []
    right_labels = []
    left_inner_labels = []
    right_inner_labels = []
    outer_labels_rearranged = []
    inner_labels_rearranged = []
    outer_labels = []
    inner_labels = []
    label_list = []
    label_filtering = config_dict['labels']['filtering']
    length_threshold = config_dict['labels']['length_threshold']['circular']
    length_param = determine_length_parameter(total_length, length_threshold)
    track_type = config_dict['canvas']['circular']['track_type']
    strandedness = config_dict['canvas']['strandedness']
    label_radius_offset = config_dict['labels']['unified_adjustment']

    strands = "separate" if strandedness else "single"
    allow_inner_labels = config_dict['canvas']['circular']['allow_inner_labels']
    radius_factor = config_dict['labels']['radius_factor'][track_type][strands][length_param] 
    inner_radius_factor = config_dict['labels']['inner_radius_factor'][track_type][strands][length_param]  

    font_family = config_dict['objects']['text']['font_family']
    font_size: str = config_dict['labels']['font_size'][length_param]
    interval = config_dict['canvas']['dpi']

    track_ratio_factor = config_dict['canvas']['circular']['track_ratio_factors'][length_param][0]
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
            longest_segment_length_in_pixels = (2*math.pi*radius_factor*radius) * (longest_segment_length)/total_length
            bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
            label_middle = longeset_segment_middle
            label_as_feature_length = total_length * bbox_width_px/(2*math.pi*radius)
            label_start = label_middle - (label_as_feature_length/2)
            label_end = label_middle + (label_as_feature_length/2)
            feature_middle_x: float = (radius * factors[1]) * math.cos(math.radians(360.0 * ((label_middle) / total_length) - 90))
            feature_middle_y: float = (radius * factors[1]) * math.sin(math.radians(360.0 * ((label_middle) / total_length) - 90))
            if feature_object.strand == "positive" or allow_inner_labels == False:
                middle_x: float = (radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90)) # 1.05?
                middle_y: float = (radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90)) # 1.05?  
            else:
                middle_x: float = (inner_radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y: float = (inner_radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            if bbox_width_px  < longest_segment_length_in_pixels *1.05:
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
    if allow_inner_labels:
        inner_labels_rearranged = rearrange_labels_fc(inner_labels, radius, total_length, length_param, config_dict, strands, is_outer=False)
    
    label_list_fc = embedded_labels + outer_labels_rearranged + inner_labels_rearranged
    return label_list_fc

def check_label_overlap(label1, label2):
   """Check if two labels overlap horizontally"""
   return not (label1["end"] < label2["start"] or label2["end"] < label1["start"])

def find_lowest_available_track(track_dict, label):
   """Find the lowest track (closest to track_1) where the label can be placed without overlap"""
   track_num = 1
   while True:
       track_id = f"track_{track_num}"
       # Check if this track has overlaps
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
   """
   Prepares a list of labels for linear genome visualization with proper track organization.
   """
   embedded_labels = []
   external_labels = []
   track_dict = defaultdict(list)
   feature_track_positions = {}  # Store feature track positions
   
   # Get configuration values
   length_threshold = config_dict['labels']['length_threshold']['circular']
   length_param = determine_length_parameter(genome_length, length_threshold)
   font_family = config_dict['objects']['text']['font_family']
   font_size = config_dict['labels']['font_size']['linear'][length_param]
   interval = config_dict['canvas']['dpi']
   label_filtering = config_dict['labels']['filtering']
   # First pass: Calculate feature track positions
   for feature_id, feature_object in feature_dict.items():
       if len(feature_object.coordinates) == 0:
           continue
           
       # Get feature track info
       coordinate = feature_object.coordinates[0]  # Use first coordinate for strand
       strand = get_strand(coordinate.strand)
       feature_track_id = feature_object.feature_track_id
       
       # Calculate track position using the same logic as for features
       factors = calculate_feature_position_factors_linear(strand, feature_track_id, strandedness)
       
       track_y_position = cds_height * factors[1]  # Use middle factor
       
       # Store track position for this feature
       feature_track_positions[feature_id] = track_y_position
   
   # Second pass: Process labels
   for feature_id, feature_object in reversed(list(feature_dict.items())):
       feature_label_text = get_label_text(feature_object, label_filtering)
       feature_track_id = feature_object.feature_track_id
       if not feature_label_text:
           continue
           
       # Calculate label dimensions and positions
       bbox_width_px, bbox_height_px = calculate_bbox_dimensions(
           feature_label_text, font_family, font_size, interval)
       
       # Find the longest segment and its middle point
       longest_segment_length = 0
       label_middle = 0
       coordinate_strand = None
       factors = None
       longest_segment_start = 0
       longest_segment_end = 0
       
       for coordinate in feature_object.coordinates:
           coordinate_strand = get_strand(coordinate.strand)
           factors = calculate_feature_position_factors_linear(coordinate_strand, strandedness, feature_track_id)
           start = int(coordinate.start)
           end = int(coordinate.end)
           segment_length = abs(end - start + 1)
           if segment_length > longest_segment_length:
               longest_segment_start = start
               longest_segment_end = end
               longest_segment_length = segment_length
               label_middle = (end + start) / 2
               
       # Calculate normalized positions
       normalized_start = normalize_position_to_linear_track(
           longest_segment_start, genome_length, alignment_width, genome_size_normalization_factor)
       normalized_end = normalize_position_to_linear_track(
           longest_segment_end, genome_length, alignment_width, genome_size_normalization_factor)
       longest_segment_length_in_pixels = abs(normalized_end - normalized_start) + 1
       normalized_middle = (normalized_start + normalized_end) / 2
       bbox_start = normalized_middle - (bbox_width_px / 2)
       bbox_end = normalized_middle + (bbox_width_px / 2)
       
       # Get actual feature track position
       feature_y = feature_track_positions.get(feature_id, cds_height * factors[1])
       
       # Create base label entry
       label_entry = {
           "label_text": feature_label_text,
           "middle": normalized_middle,
           "start": bbox_start,
           "end": bbox_end,
           "middle_x": normalized_middle,
           "width_px": bbox_width_px,
           "height_px": bbox_height_px,
           "strand": coordinate_strand,
           "feature_middle_y": feature_y,  # Use actual feature position
           "font_size": font_size,
           "font_family": font_family
       }

       # Determine if label should be embedded
       if bbox_width_px < longest_segment_length_in_pixels:
           label_entry.update({
               "middle_y": feature_y,  # Use actual feature position for embedded labels
               "is_embedded": True,
               "track_id": "track_0"
           })
           track_dict["track_0"].append(label_entry)
       else:
           # For external labels: find lowest possible track
           label_entry.update({
               "middle_y": 0,
               "is_embedded": False
           })
           
           # Find lowest track where label can be placed without overlaps
           best_track = find_lowest_available_track(track_dict, label_entry)
           label_entry["track_id"] = f"track_{best_track}"
           track_dict[f"track_{best_track}"].append(label_entry)

   # Process embedded labels
   if "track_0" in track_dict:
       for label in track_dict["track_0"]:
           embedded_labels.append(label)
   
   # Adjust track height based on strandedness
   if strandedness:
       track_height = 0.275 * cds_height  # Reduced height for separate strands mode
   else:
       track_height = 0.50 * cds_height  # Original height for default mode
   for track_id in sorted(track_dict.keys()):
       if track_id == "track_0":
           continue
           
       track_labels = track_dict[track_id]
       track_num = int(track_id.split("_")[1])
       
       for label in track_labels:
           # Compact vertical positioning for external labels
           label["middle_y"] = (-0.75 * cds_height - (track_height * track_num))
           external_labels.append(label)

   return embedded_labels + external_labels

