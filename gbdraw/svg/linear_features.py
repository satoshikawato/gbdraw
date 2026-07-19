#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Tuple

from .arrows import set_arrow_shoulder
from ..layout.linear import calculate_feature_position_factors_linear
from ..layout.linear_coords import normalize_position_to_linear_track


FeatureYPositions = Tuple[float, float, float]


def _resolve_feature_y_positions(
    *,
    cds_height: float,
    feature_strand: str,
    separate_strands: bool,
    feature_track_id: int,
    track_layout: str,
    track_axis_gap: float | None,
    feature_y_positions: FeatureYPositions | None,
) -> FeatureYPositions:
    """Resolve legacy factors only when no measured lane was supplied."""

    if feature_y_positions is not None:
        top_y, middle_y, bottom_y = feature_y_positions
        return float(top_y), float(middle_y), float(bottom_y)
    axis_gap_factor = (
        (float(track_axis_gap) / float(cds_height))
        if (track_axis_gap is not None and float(cds_height) > 0.0)
        else None
    )
    factors = calculate_feature_position_factors_linear(
        strand=feature_strand,
        track_id=feature_track_id,
        separate_strands=separate_strands,
        track_layout=track_layout,
        axis_gap_factor=axis_gap_factor,
    )
    height = float(cds_height)
    return (
        height * float(factors[0]),
        height * float(factors[1]),
        height * float(factors[2]),
    )


def create_intron_path_linear(
    coord_dict: dict,
    genome_length: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    cds_height: float,
    feature_strand: str,
    separate_strands: bool,
    feature_track_id: int,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
) -> list[str]:
    """
    Creates a linear SVG path for an intron feature.
    """
    feat_start = coord_dict["feat_start"]
    feat_end = coord_dict["feat_end"]
    _top_y, middle_y, _bottom_y = _resolve_feature_y_positions(
        cds_height=cds_height,
        feature_strand=feature_strand,
        separate_strands=separate_strands,
        feature_track_id=feature_track_id,
        track_layout=track_layout,
        track_axis_gap=track_axis_gap,
        feature_y_positions=feature_y_positions,
    )

    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor
    )

    feature_path = f"M {normalized_start},{middle_y} L{normalized_end},{middle_y} z"
    return ["line", feature_path]


def create_rectangle_path_linear(
    coord_dict: dict,
    genome_length: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    cds_height: float,
    feature_strand: str,
    separate_strands: bool,
    feature_track_id: int,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
) -> list[str]:
    """Creates a linear SVG path for a rectangular feature."""
    feat_start = coord_dict["feat_start"]
    feat_end = coord_dict["feat_end"]
    start_y_top, _middle_y, start_y_bottom = _resolve_feature_y_positions(
        cds_height=cds_height,
        feature_strand=feature_strand,
        separate_strands=separate_strands,
        feature_track_id=feature_track_id,
        track_layout=track_layout,
        track_axis_gap=track_axis_gap,
        feature_y_positions=feature_y_positions,
    )

    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor
    )

    feature_path = (
        f"M {normalized_start},{start_y_top} "
        f"L {normalized_end},{start_y_top} "
        f"L {normalized_end},{start_y_bottom} "
        f"L {normalized_start},{start_y_bottom} z"
    )
    return ["block", feature_path]


def get_arrow_strand_positions(normalized_start: float, normalized_end: float) -> dict[str, list[float]]:
    """Get start and end positions for arrow based on strand."""
    return {"positive": [normalized_start, normalized_end], "negative": [normalized_end, normalized_start]}


def calculate_normalized_arrow_length(
    arrow_length: float, genome_length: int, alignment_width: float, genome_size_normalization_factor: float
) -> float:
    return alignment_width * (arrow_length / genome_length) * genome_size_normalization_factor


def normalize_feature_positions(
    coord_dict: dict, genome_length: int, alignment_width: float, genome_size_normalization_factor: float
) -> Tuple[float, float]:
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])

    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width, genome_size_normalization_factor
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width, genome_size_normalization_factor
    )

    return normalized_start, normalized_end


def construct_arrowhead_path(
    arrow_start: float,
    arrow_end: float,
    shoulder: float,
    factors: list[float] | None,
    normalized_feat_len: float,
    normalized_arrow_length: float,
    cds_height: float,
    feature_y_positions: FeatureYPositions | None = None,
) -> list[str]:
    """Construct the SVG path for an arrowhead."""
    if feature_y_positions is None:
        if factors is None:
            raise ValueError("factors or feature_y_positions is required")
        top_y, middle_y, bottom_y = (
            cds_height * factors[0],
            cds_height * factors[1],
            cds_height * factors[2],
        )
    else:
        top_y, middle_y, bottom_y = feature_y_positions
    point_x, point_y = arrow_end, middle_y
    start_x_1, start_y_1 = arrow_start, top_y
    start_x_2, start_y_2 = arrow_start, bottom_y

    if abs(normalized_feat_len) < normalized_arrow_length:
        feature_path = f"M {start_x_1},{start_y_1} L {point_x},{point_y} L {start_x_2},{start_y_2} z"
    else:
        end_x_1, end_y_1 = shoulder, top_y
        end_x_2, end_y_2 = shoulder, bottom_y
        feature_path = (
            f"M {start_x_1},{start_y_1} "
            f"L {end_x_1},{end_y_1} "
            f"L {point_x},{point_y} "
            f"L {end_x_2},{end_y_2} "
            f"L {start_x_2},{start_y_2} z"
        )

    return ["block", feature_path]


def create_arrowhead_path_linear(
    coord_dict: dict,
    arrow_length: float,
    cds_height: float,
    feature_strand: str,
    genome_length: int,
    alignment_width: float,
    genome_size_normalization_factor: float,
    separate_strands: bool,
    feature_track_id: int = 0,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
) -> list[str]:
    """Creates a linear SVG path for an arrowhead feature."""
    normalized_start, normalized_end = normalize_feature_positions(
        coord_dict, genome_length, alignment_width, genome_size_normalization_factor
    )
    normalized_feat_len = normalized_end - normalized_start

    normalized_arrow_length = calculate_normalized_arrow_length(
        arrow_length, genome_length, alignment_width, genome_size_normalization_factor
    )

    arrow_start, arrow_end = get_arrow_strand_positions(normalized_start, normalized_end)[coord_dict["feat_strand"]]
    shoulder = set_arrow_shoulder(coord_dict["feat_strand"], arrow_end, normalized_arrow_length)
    resolved_y_positions = _resolve_feature_y_positions(
        cds_height=cds_height,
        feature_strand=feature_strand,
        separate_strands=separate_strands,
        feature_track_id=feature_track_id,
        track_layout=track_layout,
        track_axis_gap=track_axis_gap,
        feature_y_positions=feature_y_positions,
    )

    return construct_arrowhead_path(
        arrow_start,
        arrow_end,
        shoulder,
        None,
        normalized_feat_len,
        normalized_arrow_length,
        cds_height,
        feature_y_positions=resolved_y_positions,
    )


__all__ = [
    "calculate_normalized_arrow_length",
    "construct_arrowhead_path",
    "create_arrowhead_path_linear",
    "create_intron_path_linear",
    "create_rectangle_path_linear",
    "get_arrow_strand_positions",
    "normalize_feature_positions",
]


