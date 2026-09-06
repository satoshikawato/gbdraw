#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame

from ..layout.linear_coords import normalize_position_to_linear_track


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from ..layout.scalar_axis import project_scalar_samples

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
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    if gc_df.empty:
        return ""
    coodinates_list: list[str] = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column: str = f"{dinucleotide} content"
    mean = float(gc_df[column].mean())
    max_diff = float((gc_df[column] - mean).abs().max())
    samples = gc_df[column].items()
    if record_transform is not None and record_transform.start_coordinate is not None:
        samples = ((p.position, p.value) for segment in project_scalar_samples(
            gc_df.index, gc_df[column], record_transform) for p in segment.points)
    x_corrdinate = start_x
    for index, value in samples:
        value = float(value)
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
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    if skew_df.empty:
        return ""
    coodinates_list: list[str] = []
    start_position: str = "M{} {}".format(start_x, start_y)
    coodinates_list.append(start_position)
    column = [col for col in skew_df.columns if "skew" in col and "cumulative" not in col.lower()][0]
    mean = float(skew_df[column].mean())
    max_diff = float((skew_df[column] - mean).abs().max()) if float((skew_df[column] - mean).abs().max()) > 0 else 1.0

    samples = skew_df[column].items()
    if record_transform is not None and record_transform.start_coordinate is not None:
        samples = ((p.position, p.value) for segment in project_scalar_samples(
            skew_df.index, skew_df[column], record_transform) for p in segment.points)
    x_corrdinate = start_x
    for index, value in samples:
        value = float(value)
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


def calculate_linear_scalar_area_path_desc(
    start_x: float,
    start_y: float,
    scalar_df: DataFrame,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
    *,
    value_column: str = "value_normalized",
    position_column: str = "position",
    source_positions: bool = False,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return a filled linear area path for normalized scalar values."""

    if scalar_df.empty or record_len <= 0 or track_height <= 0:
        return ""

    baseline_y = float(start_y) + float(track_height)
    x_values: list[float] = []
    y_values: list[float] = []
    samples = zip(scalar_df[position_column], scalar_df[value_column], strict=True)
    projected = record_transform is not None and record_transform.start_coordinate is not None
    if projected:
        samples = ((p.position, p.value) for segment in project_scalar_samples(
            scalar_df[position_column], scalar_df[value_column], record_transform,
            source_positions=source_positions) for p in segment.points)
    for position, value in samples:
        position = float(position) if projected else int(position)
        value = max(0.0, min(1.0, float(value)))
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


def calculate_depth_path_desc(
    start_x: float,
    start_y: float,
    depth_df: DataFrame,
    record_len: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    track_height: float,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return a filled linear area path for binned depth coverage."""

    return calculate_linear_scalar_area_path_desc(
        start_x,
        start_y,
        depth_df,
        record_len,
        alignment_width,
        genome_size_normalization_factor,
        track_height,
        value_column="depth_normalized",
        source_positions=True,
        record_transform=record_transform,
    )


__all__ = [
    "calculate_corrdinate",
    "calculate_depth_path_desc",
    "calculate_gc_content_path_desc",
    "calculate_gc_skew_path_desc",
    "calculate_linear_scalar_area_path_desc",
]


