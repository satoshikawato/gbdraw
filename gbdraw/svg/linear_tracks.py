#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame

from ..layout.linear_coords import normalize_position_to_linear_track


def calculate_corrdinate(
    index: int,
    value: float,
    mean: float,
    max_diff: float,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
):
    diff: float = value - mean
    x_corrdinate: float = normalize_position_to_linear_track(
        index, record_len, alignment_width, genome_size_normalization_factor
    )
    # If the series is completely flat, avoid division by zero and draw a flat line.
    if max_diff == 0:
        y_corrdinate = 0.0
    else:
        y_corrdinate = -(0.5 * track_height * (diff / max_diff))
    corrdinate: str = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
    return corrdinate, x_corrdinate


def calculate_gc_content_path_desc(
    start_x: float,
    start_y: float,
    gc_df: DataFrame,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
    dinucleotide: str,
) -> str:
    coodinates_list: list[str] = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column: str = f"{dinucleotide} content"
    mean = float(gc_df[column].mean())
    max_diff = float((gc_df[column] - mean).abs().max())
    for index, row in gc_df.iterrows():
        value = float(row[column])
        corrdinate, x_corrdinate = calculate_corrdinate(
            index, value, mean, max_diff, record_len, alignment_width, genome_size_normalization_factor, track_height
        )
        coodinates_list.append(corrdinate)
    penultimate_coordinate: str = "L{} {}".format(str(x_corrdinate), str(start_y))
    coodinates_list.append(penultimate_coordinate)
    end_coordinate: str = "L{} {}".format(str(start_x), str(start_y))
    coodinates_list.append(end_coordinate)
    gc_content_desc: str = "{}".format("".join(coodinates_list))
    gc_content_desc += "z"
    return gc_content_desc


def calculate_gc_skew_path_desc(
    start_x: float,
    start_y: float,
    skew_df: DataFrame,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
) -> str:
    coodinates_list: list[str] = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column = [col for col in skew_df.columns if "skew" in col and "cumulative" not in col.lower()][0]
    mean = float(skew_df[column].mean())
    max_diff = float((skew_df[column] - mean).abs().max()) if float((skew_df[column] - mean).abs().max()) > 0 else 1.0

    for index, row in skew_df.iterrows():
        value = float(row[column])
        corrdinate, x_corrdinate = calculate_corrdinate(
            index, value, mean, max_diff, record_len, alignment_width, genome_size_normalization_factor, track_height
        )
        coodinates_list.append(corrdinate)

    penultimate_coordinate: str = "L{} {}".format(str(x_corrdinate), str(start_y))
    coodinates_list.append(penultimate_coordinate)
    end_coordinate: str = "L{} {}".format(str(start_x), str(start_y))
    coodinates_list.append(end_coordinate)
    gc_skew_desc: str = "{}".format("".join(coodinates_list))
    gc_skew_desc += "z"
    return gc_skew_desc


__all__ = ["calculate_corrdinate", "calculate_gc_content_path_desc", "calculate_gc_skew_path_desc"]


