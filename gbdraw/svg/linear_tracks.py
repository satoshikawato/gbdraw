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


def calculate_depth_path_desc(
    start_x: float,
    start_y: float,
    depth_df: DataFrame,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
) -> str:
    """Return a filled linear area path for binned depth coverage."""

    if depth_df.empty or record_len <= 0 or track_height <= 0:
        return ""

    baseline_y = float(start_y) + float(track_height)
    x_values: list[float] = []
    y_values: list[float] = []
    for _, row in depth_df.iterrows():
        position = int(row["position"])
        value = max(0.0, min(1.0, float(row["depth_normalized"])))
        x_value = normalize_position_to_linear_track(
            position, record_len, alignment_width, genome_size_normalization_factor
        )
        y_value = baseline_y - (float(track_height) * value)
        x_values.append(float(x_value))
        y_values.append(float(y_value))

    if not x_values:
        return ""

    final_x = normalize_position_to_linear_track(
        record_len, record_len, alignment_width, genome_size_normalization_factor
    )
    path_segments = [f"M{start_x} {baseline_y}", f"L{x_values[0]} {y_values[0]}"]
    path_segments.extend(
        f"L{x_value} {y_value}" for x_value, y_value in zip(x_values[1:], y_values[1:])
    )
    path_segments.append(f"L{final_x} {y_values[-1]}")
    path_segments.append(f"L{final_x} {baseline_y}")
    path_segments.append(f"L{start_x} {baseline_y}z")
    return "".join(path_segments)


__all__ = [
    "calculate_corrdinate",
    "calculate_depth_path_desc",
    "calculate_gc_content_path_desc",
    "calculate_gc_skew_path_desc",
]


