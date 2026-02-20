#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations


def _legacy_middle_factors(strand: str, track_id: int, separate_strands: bool) -> list[float]:
    """Legacy factor calculation used for backward-compatible 'middle' layout."""
    initial_offset = 0.1
    track_spacing = 0.5
    feature_height = 0.5
    feature_offset = 0.25

    def calculate_base_position(track_num: int) -> float:
        return -initial_offset - track_num * track_spacing

    def calculate_track_offset(track_num: int) -> float:
        if separate_strands:
            return (
                track_num * -1.5 * initial_offset
                if strand == "positive"
                else (track_num - 1) * -1.5 * initial_offset
            )
        return track_num * -3 * initial_offset

    if not separate_strands:
        base_pos = calculate_base_position(track_id)
        track_offset = calculate_track_offset(track_id)
        track = track_id
        return [
            -track - feature_height + track_offset,
            -track + track_offset,
            -track + feature_height + track_offset,
        ]

    # Keep middle layout vertically symmetric around axis in separated-strand mode.
    # Positive track 0 and negative track -1 should be mirrored with equal axis gaps.
    half_height = feature_height / 2.0
    track_step = 0.65
    axis_to_center = 0.35
    if strand == "negative":
        negative_index = max(0, abs(int(track_id)) - 1)
        middle = axis_to_center + (track_step * negative_index)
    else:
        positive_index = max(0, int(track_id))
        middle = -(axis_to_center + (track_step * positive_index))
    return [middle - half_height, middle, middle + half_height]


def _resolve_axis_gap_factor(
    *,
    separate_strands: bool,
    axis_gap_factor: float | None,
) -> float:
    """Resolve axis-to-feature edge gap in factor units (scaled by cds_height)."""
    if axis_gap_factor is not None:
        return max(0.0, float(axis_gap_factor))
    # Auto defaults tuned for above/below layouts to keep labels and features legible.
    if separate_strands:
        return 0.25
    return 0.30


def calculate_feature_position_factors_linear(
    strand: str,
    track_id: int,
    separate_strands: bool,
    track_layout: str = "middle",
    axis_gap_factor: float | None = None,
) -> list[float]:
    """
    Calculates feature position factors for linear track layouts.

    Returns:
        list[float]: Three position factors [top, middle, bottom] for feature placement.
    """
    layout = str(track_layout).strip().lower()
    if layout in {"spreadout", "above"}:
        layout = "above"
    elif layout in {"tuckin", "below"}:
        layout = "below"
    else:
        layout = "middle"

    if layout == "middle":
        return _legacy_middle_factors(strand, track_id, separate_strands)

    resolved_axis_gap = _resolve_axis_gap_factor(
        separate_strands=separate_strands,
        axis_gap_factor=axis_gap_factor,
    )

    if not separate_strands:
        track_index = max(0, abs(int(track_id)))
        track_spacing = 1.3
        half_height = 0.5
        middle = resolved_axis_gap + half_height + (track_spacing * track_index)
        if layout == "above":
            middle = -middle
        return [middle - half_height, middle, middle + half_height]

    # Separate-strands above/below modes keep positive strand above negative strand
    # while shifting both bands to the selected side of the axis.
    track_spacing = 0.65
    half_height = 0.25
    strand_band_gap = 0.6
    positive_index = max(0, int(track_id))
    negative_index = max(0, abs(int(track_id)) - 1)

    if layout == "above":
        if strand == "negative":
            middle = -(resolved_axis_gap + half_height + (track_spacing * negative_index))
        else:
            middle = -(resolved_axis_gap + half_height + strand_band_gap + (track_spacing * positive_index))
    else:  # below
        if strand == "negative":
            middle = resolved_axis_gap + half_height + strand_band_gap + (track_spacing * negative_index)
        else:
            middle = resolved_axis_gap + half_height + (track_spacing * positive_index)

    return [middle - half_height, middle, middle + half_height]


__all__ = ["calculate_feature_position_factors_linear"]


