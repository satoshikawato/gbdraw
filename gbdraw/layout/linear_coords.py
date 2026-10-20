#!/usr/bin/env python
# coding: utf-8

def normalize_position_to_linear_track(
    position: int,
    genome_length: float,
    alignment_width: float,
    genome_size_normalization_factor: float,
) -> float:
    normalized_position: float = (
        alignment_width * (position / genome_length) * genome_size_normalization_factor
    )
    return normalized_position


__all__ = ["normalize_position_to_linear_track"]


