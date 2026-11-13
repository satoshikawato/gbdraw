#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Tuple

from .arrows import (
    ArrowHeadLengthRatio,
    LinearBlockOutline,
    block_path_with_open_caps,
    calculate_normalized_arrow_length,
    linear_block_outline,
)
from ..layout.linear import calculate_feature_position_factors_linear
from ..layout.linear_coords import normalize_position_to_linear_track
from .path_data import path_point

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
        feat_start, genome_length, alignment_width
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width
    )

    feature_path = f"M {path_point(normalized_start, middle_y)} L{path_point(normalized_end, middle_y)} z"
    return ["line", feature_path]


def create_rectangle_path_linear(
    coord_dict: dict,
    genome_length: int,
    alignment_width: float,
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
        feat_start, genome_length, alignment_width
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width
    )

    feature_path = (
        f"M {path_point(normalized_start, start_y_top)} "
        f"L {path_point(normalized_end, start_y_top)} "
        f"L {path_point(normalized_end, start_y_bottom)} "
        f"L {path_point(normalized_start, start_y_bottom)} z"
    )
    return block_path_with_open_caps(
        feature_path, open_tail=bool(coord_dict.get("open_start")),
        end_cap=f"L {path_point(normalized_end, start_y_bottom)}" if coord_dict.get("open_end") else None,
        tail_cap=f"L {path_point(normalized_start, start_y_top)}",
    )


def normalize_feature_positions(
    coord_dict: dict, genome_length: int, alignment_width: float
) -> Tuple[float, float]:
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])

    normalized_start = normalize_position_to_linear_track(
        feat_start, genome_length, alignment_width
    )
    normalized_end = normalize_position_to_linear_track(
        feat_end, genome_length, alignment_width
    )

    return normalized_start, normalized_end


def _arrow_path(outline: LinearBlockOutline, full_width: bool) -> str:
    """The path text of an arrow outline: a triangle, the legacy five vertices at full shaft width, else
    seven with the narrowed shaft."""

    tail_x, shoulder_x, tip_x = outline.tail_x, outline.shoulder_x, outline.tip_x
    top_y, middle_y, bottom_y = outline.top_y, outline.middle_y, outline.bottom_y
    if outline.triangle:
        return f"M {path_point(tail_x, top_y)} L {path_point(tip_x, middle_y)} L {path_point(tail_x, bottom_y)} z"
    if full_width:
        return (
            f"M {path_point(tail_x, top_y)} "
            f"L {path_point(shoulder_x, top_y)} "
            f"L {path_point(tip_x, middle_y)} "
            f"L {path_point(shoulder_x, bottom_y)} "
            f"L {path_point(tail_x, bottom_y)} z"
        )
    return (
        f"M {path_point(tail_x, outline.shaft_top_y)} "
        f"L {path_point(shoulder_x, outline.shaft_top_y)} "
        f"L {path_point(shoulder_x, top_y)} "
        f"L {path_point(tip_x, middle_y)} "
        f"L {path_point(shoulder_x, bottom_y)} "
        f"L {path_point(shoulder_x, outline.shaft_bottom_y)} "
        f"L {path_point(tail_x, outline.shaft_bottom_y)} z"
    )


def create_arrow_path_linear(
    coord_dict: dict,
    arrow_length: float,
    cds_height: float,
    feature_strand: str,
    genome_length: int,
    alignment_width: float,
    separate_strands: bool,
    feature_track_id: int = 0,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
    head_length_ratio: ArrowHeadLengthRatio = "auto",
    shaft_width_ratio: float = 1.0,
) -> list[str]:
    """Create a Linear arrow from its outline (``linear_block_outline``), the legacy path at full shaft
    width."""
    normalized_start, normalized_end = normalize_feature_positions(
        coord_dict, genome_length, alignment_width
    )
    resolved_y_positions = _resolve_feature_y_positions(
        cds_height=cds_height,
        feature_strand=feature_strand,
        separate_strands=separate_strands,
        feature_track_id=feature_track_id,
        track_layout=track_layout,
        track_axis_gap=track_axis_gap,
        feature_y_positions=feature_y_positions,
    )
    outline = linear_block_outline(
        "arrow",
        normalized_start,
        normalized_end,
        coord_dict["feat_strand"],
        resolved_y_positions,
        calculate_normalized_arrow_length(arrow_length, genome_length, alignment_width),
        head_length_ratio,
        shaft_width_ratio,
    )
    return block_path_with_open_caps(
        _arrow_path(outline, float(shaft_width_ratio) == 1.0),
        open_tail=bool(coord_dict.get("open_start" if coord_dict["feat_strand"] == "positive" else "open_end")),
    )


def create_arrowhead_path_linear(
    coord_dict: dict,
    arrow_length: float,
    cds_height: float,
    feature_strand: str,
    genome_length: int,
    alignment_width: float,
    separate_strands: bool,
    feature_track_id: int = 0,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
    head_length_ratio: ArrowHeadLengthRatio = "auto",
) -> list[str]:
    """Compatibility wrapper for the legacy five-vertex Linear arrow."""
    return create_arrow_path_linear(
        coord_dict=coord_dict,
        arrow_length=arrow_length,
        cds_height=cds_height,
        feature_strand=feature_strand,
        genome_length=genome_length,
        alignment_width=alignment_width,
        separate_strands=separate_strands,
        feature_track_id=feature_track_id,
        track_layout=track_layout,
        track_axis_gap=track_axis_gap,
        feature_y_positions=feature_y_positions,
        head_length_ratio=head_length_ratio,
        shaft_width_ratio=1.0,
    )


def create_arrow_shaft_path_linear(
    coord_dict: dict,
    genome_length: int,
    alignment_width: float,
    cds_height: float,
    feature_strand: str,
    separate_strands: bool,
    feature_track_id: int,
    shaft_width_ratio: float,
    track_layout: str = "middle",
    track_axis_gap: float | None = None,
    feature_y_positions: FeatureYPositions | None = None,
) -> list[str]:
    """Create a reduced-width shaft block for a multipart arrow."""
    normalized_start, normalized_end = normalize_feature_positions(
        coord_dict, genome_length, alignment_width
    )
    outline = linear_block_outline(
        "shaft",
        normalized_start,
        normalized_end,
        coord_dict["feat_strand"],
        _resolve_feature_y_positions(
            cds_height=cds_height,
            feature_strand=feature_strand,
            separate_strands=separate_strands,
            feature_track_id=feature_track_id,
            track_layout=track_layout,
            track_axis_gap=track_axis_gap,
            feature_y_positions=feature_y_positions,
        ),
        0.0,
        "auto",
        shaft_width_ratio,
    )
    shaft_top_y, shaft_bottom_y = outline.shaft_top_y, outline.shaft_bottom_y
    feature_path = (
        f"M {path_point(normalized_start, shaft_top_y)} "
        f"L {path_point(normalized_end, shaft_top_y)} "
        f"L {path_point(normalized_end, shaft_bottom_y)} "
        f"L {path_point(normalized_start, shaft_bottom_y)} z"
    )
    return block_path_with_open_caps(
        feature_path, open_tail=bool(coord_dict.get("open_start")),
        end_cap=f"L {path_point(normalized_end, shaft_bottom_y)}" if coord_dict.get("open_end") else None,
        tail_cap=f"L {path_point(normalized_start, shaft_top_y)}",
    )


__all__ = [
    "create_arrow_path_linear",
    "create_arrow_shaft_path_linear",
    "create_arrowhead_path_linear",
    "create_intron_path_linear",
    "create_rectangle_path_linear",
    "normalize_feature_positions",
]

