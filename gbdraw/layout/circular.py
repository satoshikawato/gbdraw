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
    track_id: int = 0,
) -> list[float]:
    """
    Calculates position factors for a feature based on its strand orientation on a circular canvas.

    The factors determine the inner, middle, and outer radii of the feature arc relative to
    the base radius. When resolve_overlaps is enabled with strandedness=False, track_id is
    used to offset features to avoid visual overlap.

    Args:
        total_length: Total length of the genome (not used in calculation but kept for API consistency)
        strand: "positive", "negative", or "undefined"
        track_ratio: Base track ratio from config
        cds_ratio: Calculated CDS ratio for feature height
        offset: Base offset value
        track_type: "tuckin", "middle", or "spreadout"
        strandedness: Whether to separate strands
        track_id: Track number for overlap resolution. Only effective when strandedness=False.
                  - spreadout: displaced outward
                  - middle: positive/undefined displaced outward, negative displaced inward
                  - tuckin: displaced inward

    Returns:
        List of three floats [inner_factor, middle_factor, outer_factor] used to calculate
        the actual radii by multiplying with the base radius.
    """
    BASE: float = 1.0
    cds_ratio = float(cds_ratio)
    offset = float(offset)

    lane_width = cds_ratio
    lane_spacing = 0.01
    lane_step = lane_width + lane_spacing

    def factors_from_center(center: float) -> list[float]:
        return [
            center - (0.5 * lane_width),
            center,
            center + (0.5 * lane_width),
        ]

    if strandedness is True:
        if track_type == "middle":
            center_positive = BASE + (0.5 * lane_spacing) + (0.5 * lane_width)
            center_negative = BASE - (0.5 * lane_spacing) - (0.5 * lane_width)
        elif track_type == "spreadout":
            center_positive = BASE + lane_spacing + (0.5 * lane_width)
            center_negative = center_positive + lane_step
        elif track_type == "tuckin":
            center_positive = BASE - lane_spacing - (0.5 * lane_width)
            center_negative = center_positive - lane_step
        else:
            center_positive = BASE + (0.5 * lane_spacing) + (0.5 * lane_width)
            center_negative = BASE - (0.5 * lane_spacing) - (0.5 * lane_width)

        if strand == "positive":
            factors = factors_from_center(center_positive + offset)
        else:
            factors = factors_from_center(center_negative - offset)

        return factors

    else:
        if track_type == "middle":
            if track_id != 0:
                if int(track_id) < 0 or strand == "negative":
                    center = BASE - (abs(track_id) * lane_step)
                else:
                    center = BASE + (abs(track_id) * lane_step)
            else:
                center = BASE
        elif track_type == "spreadout":
            center = BASE + lane_spacing + (0.5 * lane_width) + (abs(track_id) * lane_step)
        elif track_type == "tuckin":
            center = BASE - lane_spacing - (0.5 * lane_width) - (abs(track_id) * lane_step)
        else:
            center = BASE

        return factors_from_center(center)


__all__ = ["calculate_feature_position_factors_circular"]
