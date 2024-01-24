#!/usr/bin/env python
# coding: utf-8


import pandas as pd
from collections import defaultdict
from pandas import DataFrame
from typing import Generator, Any, Dict, List, Optional
from Bio.SeqRecord import SeqRecord


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