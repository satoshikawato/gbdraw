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
                  - middle: displaced outward
                  - tuckin: displaced inward
    
    Returns:
        List of three floats [inner_factor, middle_factor, outer_factor] used to calculate
        the actual radii by multiplying with the base radius.
    """
    BASE: float = 1.0
    cds_ratio = float(cds_ratio)
    offset = float(offset)
    
    # Calculate track offset for resolve_overlaps
    # Each track is spaced by cds_ratio * TRACK_SPACING_MULTIPLIER
    TRACK_SPACING_MULTIPLIER = 1.2
    track_offset = abs(track_id) * cds_ratio * TRACK_SPACING_MULTIPLIER
    
    if strandedness is True:
        # With strand separation: resolve_overlaps is NOT supported
        # track_id is ignored
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
        
        return factors
    
    else:
        # No strand separation: resolve_overlaps is supported
        if track_type == "middle":
            base_factors = [BASE - cds_ratio * 0.5, BASE, BASE + cds_ratio * 0.5]
            # resolve_overlaps: push overlapping features OUTWARD
            if track_id != 0:
                base_factors = [x + track_offset for x in base_factors]
        elif track_type == "spreadout":
            base_factors = [BASE + cds_ratio * 0.4, BASE + cds_ratio * 0.9, BASE + cds_ratio * 1.4]
            # resolve_overlaps: push overlapping features OUTWARD
            if track_id != 0:
                base_factors = [x + track_offset for x in base_factors]
        elif track_type == "tuckin":
            base_factors = [BASE - cds_ratio * 1.7, BASE - cds_ratio * 1.2, BASE - cds_ratio * 0.7]
            # resolve_overlaps: push overlapping features INWARD (toward center)
            if track_id != 0:
                base_factors = [x - track_offset for x in base_factors]
        else:
            base_factors = [BASE - cds_ratio * 0.5, BASE, BASE + cds_ratio * 0.5]

        return [x for x in base_factors]


__all__ = ["calculate_feature_position_factors_circular"]
