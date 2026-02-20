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

    base_pos = calculate_base_position(track_id)
    track_offset = calculate_track_offset(track_id)

    if not separate_strands:
        track = track_id
        return [
            -track - feature_height + track_offset,
            -track + track_offset,
            -track + feature_height + track_offset,
        ]

    return [
        base_pos - feature_height + track_offset,
        base_pos - feature_offset + track_offset,
        base_pos + track_offset,
    ]


def calculate_feature_position_factors_linear(
    strand: str,
    track_id: int,
    separate_strands: bool,
    track_layout: str = "middle",
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

    if not separate_strands:
        track_index = max(0, abs(int(track_id)))
        track_spacing = 1.3
        half_height = 0.5
        axis_gap = 0.1
        middle = axis_gap + half_height + (track_spacing * track_index)
        if layout == "above":
            middle = -middle
        return [middle - half_height, middle, middle + half_height]

    # Separate-strands above/below modes keep positive strand above negative strand
    # while shifting both bands to the selected side of the axis.
    track_spacing = 0.65
    half_height = 0.25
    axis_gap = 0.1
    strand_band_gap = 0.6
    positive_index = max(0, int(track_id))
    negative_index = max(0, abs(int(track_id)) - 1)

    if layout == "above":
        if strand == "negative":
            middle = -(axis_gap + half_height + (track_spacing * negative_index))
        else:
            middle = -(axis_gap + half_height + strand_band_gap + (track_spacing * positive_index))
    else:  # below
        if strand == "negative":
            middle = axis_gap + half_height + strand_band_gap + (track_spacing * negative_index)
        else:
            middle = axis_gap + half_height + (track_spacing * positive_index)

    return [middle - half_height, middle, middle + half_height]


__all__ = ["calculate_feature_position_factors_linear"]


