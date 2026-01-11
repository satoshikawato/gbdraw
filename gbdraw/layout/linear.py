#!/usr/bin/env python
# coding: utf-8

def calculate_feature_position_factors_linear(
    strand: str, track_id: int, separate_strands: bool
) -> list[float]:
    """
    Calculates position factors for linear feature representation based on strand and strandedness.

    Returns:
        list[float]: Three position factors [top, middle, bottom] for feature placement.
    """
    # Constants for strand-separated layout
    INITIAL_OFFSET = 0.1  # Base offset from axis
    TRACK_SPACING = 0.5  # Vertical space per track
    FEATURE_HEIGHT = 0.5  # Height of feature
    FEATURE_OFFSET = 0.25  # Half of feature height

    def calculate_base_position(track_num: int) -> float:
        return -INITIAL_OFFSET - track_num * TRACK_SPACING

    def calculate_track_offset(track_num: int) -> float:
        if separate_strands:
            return (
                track_num * -1.5 * INITIAL_OFFSET
                if strand == "positive"
                else (track_num - 1) * -1.5 * INITIAL_OFFSET
            )
        else:
            return (
                track_num * -3 * INITIAL_OFFSET
                if strand == "positive"
                else track_num * -3 * INITIAL_OFFSET
            )

    base_pos = calculate_base_position(track_id)
    track_offset = calculate_track_offset(track_id)

    if not separate_strands:
        track = track_id
        return [
            -track - FEATURE_HEIGHT + track_offset,  # Top
            -track - 0 + track_offset,  # Middle
            -track + FEATURE_HEIGHT + track_offset,  # Bottom
        ]

    return [
        base_pos - FEATURE_HEIGHT + track_offset,  # Top
        base_pos - FEATURE_OFFSET + track_offset,  # Middle
        base_pos + track_offset,  # Bottom
    ]


__all__ = ["calculate_feature_position_factors_linear"]


