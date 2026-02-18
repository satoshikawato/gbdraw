#!/usr/bin/env python
# coding: utf-8

import math


def set_arrow_shoulder(feat_strand: str, arrow_end: float, cds_arrow_length: float) -> float:
    """
    Calculates the shoulder position for an arrow representing a genomic feature.
    """
    if feat_strand == "positive":
        shoulder = float(arrow_end - cds_arrow_length)
    else:
        shoulder = float(arrow_end + cds_arrow_length)
    return shoulder


def calculate_circular_arrow_length(total_length: int) -> float:
    """Return circular arrowhead length (bp) derived from genome size."""
    min_arrow_length = 30.0
    max_arrow_length = 700.0
    param_a = 3.0
    param_b = 5.0
    return min_arrow_length + (max_arrow_length - min_arrow_length) * (
        1.0 / (1.0 + math.exp(-param_a * (math.log10(total_length) - param_b)))
    )


__all__ = ["calculate_circular_arrow_length", "set_arrow_shoulder"]


