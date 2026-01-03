#!/usr/bin/env python
# coding: utf-8

def calculate_feature_position_factors_circular(
    total_length: int,
    strand: str,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str = "tuckin",
    strandedness: bool = True,
) -> list[float]:
    """
    Calculates position factors for a feature based on its strand orientation on a circular canvas.
    """
    BASE: float = 1.0
    cds_ratio = float(cds_ratio)
    offset = float(offset)
    if strandedness is True:
        if track_type == "middle":
            factors_positive: list[float] = [BASE, BASE + cds_ratio * 0.5, BASE + cds_ratio]
            factors_negative: list[float] = [BASE - cds_ratio, BASE - cds_ratio * 0.5, BASE]
        elif track_type == "spreadout":
            factors_positive = [BASE + cds_ratio * 1.4, BASE + cds_ratio * 1.9, BASE + cds_ratio * 2.4]
            factors_negative = [BASE + cds_ratio * 0.4, BASE + cds_ratio * 0.9, BASE + cds_ratio * 1.4]
        elif track_type == "tuckin":
            factors_positive = [BASE - cds_ratio * 1.5, BASE - cds_ratio * 1.0, BASE - cds_ratio * 0.5]
            factors_negative = [BASE - cds_ratio * 2.5, BASE - cds_ratio * 2.0, BASE - cds_ratio * 1.5]
        else:
            factors_positive = [BASE, BASE + cds_ratio * 0.5, BASE + cds_ratio]
            factors_negative = [BASE - cds_ratio, BASE - cds_ratio * 0.5, BASE]
        if strand == "positive":
            factors: list[float] = [x + offset for x in factors_positive]
        else:
            factors = [x - offset for x in factors_negative]
    else:
        # No strand separation: use the same three radii for both strands, chosen by track_type.
        if track_type == "middle":
            base_factors = [BASE - cds_ratio * 0.5, BASE, BASE + cds_ratio * 0.5]
        elif track_type == "spreadout":
            base_factors = [BASE + cds_ratio * 0.4, BASE + cds_ratio * 0.9, BASE + cds_ratio * 1.4]
        elif track_type == "tuckin":
            base_factors = [BASE - cds_ratio * 1.7, BASE - cds_ratio * 1.2, BASE - cds_ratio * 0.7]
        else:
            base_factors = [BASE - cds_ratio * 0.5, BASE, BASE + cds_ratio * 0.5]

        return [x for x in base_factors]

    return factors


__all__ = ["calculate_feature_position_factors_circular"]


