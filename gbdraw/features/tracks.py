#!/usr/bin/env python
# coding: utf-8

from typing import Dict, List, Optional, Tuple

from .objects import FeatureObject


def get_feature_ends(feature, genome_length: Optional[int] = None) -> Tuple[int, int, str]:
    """
    Get feature start/end positions and strand.
    
    For origin-spanning features on circular genomes, returns start > end
    to indicate the feature crosses the origin.
    
    Args:
        feature: FeatureObject with location/coordinates
        genome_length: Total genome length (required for origin-spanning detection)
    
    Returns:
        (start, end, strand) tuple. If start > end, feature spans the origin.
    """
    strand = feature.location[0].strand if feature.location else "undefined"
    
    if hasattr(feature, "coordinates") and feature.coordinates:
        parts = list(feature.coordinates)
    else:
        parts = list(feature.location)
    
    if not parts:
        return 1, 1, strand
    
    starts = [max(1, int(part.start)) for part in parts]
    ends = [max(1, int(part.end)) for part in parts]

    # Origin-spanning features touch both the left and right genome boundaries.
    # Avoid using part ordering for detection because negative-strand multipart features
    # can be listed in descending genomic order without crossing the origin.
    if len(parts) > 1 and genome_length:
        genome_end = int(genome_length)
        touches_left_boundary = min(starts) <= 1
        touches_right_boundary = max(ends) >= genome_end
        if touches_left_boundary and touches_right_boundary:
            # Represent origin-spanning features as start > end.
            start = max(starts)
            end = min(ends)
            return start, end, strand

    # Normal case: use min/max
    start = min(starts)
    end = max(ends)
    
    return start, end, strand


def calculate_feature_metrics(feature, genome_length: Optional[int] = None) -> Tuple[int, int]:
    """
    Calculate total span and occupied length with corrected calculations.
    
    Handles origin-spanning features correctly when genome_length is provided.

    Returns:
        (total_span, occupied_length)
    """
    start, end, _ = get_feature_ends(feature, genome_length)
    
    # Calculate total span, accounting for origin-spanning
    if start > end and genome_length:
        # Origin-spanning: spans from start to genome_length, then from 1 to end
        total_span = (genome_length - start) + end
    else:
        total_span = abs(end - start)

    occupied_length = 0
    for part in feature.location:
        if part.kind == "block":
            occupied_length += abs(part.end - part.start)

    occupied_length = min(occupied_length, total_span) if total_span > 0 else occupied_length
    return total_span, occupied_length


def check_feature_overlap(
    a: dict, b: dict, separate_strands: bool, genome_length: Optional[int] = None
) -> bool:
    """
    Check if two features overlap on a circular genome.
    
    Handles origin-spanning features when genome_length is provided.
    Origin-spanning features have start > end.
    
    Args:
        a: Feature dict with 'start', 'end', 'strand' keys
        b: Feature dict with 'start', 'end', 'strand' keys
        separate_strands: If True, features on different strands don't overlap
        genome_length: Total genome length (required for origin-spanning detection)
    
    Returns:
        True if features overlap, False otherwise
    """
    if separate_strands and a["strand"] != b["strand"]:
        return False
    
    a_start, a_end = a["start"], a["end"]
    b_start, b_end = b["start"], b["end"]
    
    # Check if either feature spans the origin
    a_spans_origin = a_start > a_end
    b_spans_origin = b_start > b_end
    
    if not a_spans_origin and not b_spans_origin:
        # Simple linear overlap check
        return not (a_end < b_start or a_start > b_end)
    
    if genome_length is None:
        # Can't determine accurately without genome_length, assume overlap
        return True
    
    # Handle origin-spanning cases
    def overlaps_circular(s1: int, e1: int, s2: int, e2: int, length: int) -> bool:
        """Check overlap for circular genome coordinates."""
        span1_origin = s1 > e1
        span2_origin = s2 > e2
        
        if not span1_origin and not span2_origin:
            # Neither spans origin
            return not (e1 < s2 or s1 > e2)
        
        if span1_origin and span2_origin:
            # Both span origin - they definitely overlap (both cover position 1 and length)
            return True
        
        if span1_origin:
            # Feature 1 spans origin: covers [s1, length] and [1, e1]
            # Feature 2 is normal: covers [s2, e2]
            # Overlap if: s2 <= length and s2 >= s1 (overlaps right part)
            #         or: e2 >= 1 and e2 <= e1 (overlaps left part)
            #         or: s2 <= e1 (feature 2 is within left part)
            #         or: e2 >= s1 (feature 2 is within right part)
            return (s2 >= s1) or (e2 <= e1) or (s2 <= e1 and e2 >= 1) or (e2 >= s1 and s2 <= length)
        
        # span2_origin is True
        # Feature 2 spans origin, feature 1 is normal
        return (s1 >= s2) or (e1 <= e2) or (s1 <= e2 and e1 >= 1) or (e1 >= s2 and s1 <= length)
    
    return overlaps_circular(a_start, a_end, b_start, b_end, genome_length)


def find_best_track(
    feature: dict,
    track_dict: Dict[str, List[dict]],
    separate_strands: bool,
    resolve_overlaps: bool,
    genome_length: Optional[int] = None,
    max_track: int = 100,
) -> int:
    """
    Find the best track for a feature, avoiding overlaps if resolve_overlaps is True.
    
    Args:
        feature: Feature metrics dict
        track_dict: Dict of existing features per track
        separate_strands: Whether to separate positive/negative strands
        resolve_overlaps: Whether to resolve overlapping features
        genome_length: Total genome length (for origin-spanning detection)
        max_track: Maximum number of tracks to consider
    
    Returns:
        Track number (positive for positive strand, negative for negative strand)
    """
    if not separate_strands:
        track_nums = [0] if not resolve_overlaps else list(range(0, max_track))
    else:
        if not resolve_overlaps:
            track_nums = [0] if feature["strand"] == "positive" else [-1]
        else:
            if feature["strand"] == "positive":
                track_nums = list(range(0, max_track))
            else:
                track_nums = list(range(-1, -max_track - 1, -1))

    if resolve_overlaps:
        for tn in track_nums:
            key = f"track_{abs(tn)}"
            if key not in track_dict or not track_dict[key]:
                return tn
            has_overlap = False
            for existing in track_dict[key]:
                if check_feature_overlap(feature, existing, separate_strands, genome_length):
                    has_overlap = True
                    break
            if not has_overlap:
                return tn

    return track_nums[0]


def arrange_feature_tracks(
    feature_dict: Dict[str, FeatureObject],
    separate_strands: bool,
    resolve_overlaps: bool,
    genome_length: Optional[int] = None,
) -> Dict[str, FeatureObject]:
    """
    Arrange features in tracks with improved strand handling and track assignment.
    
    Args:
        feature_dict: Dict of feature_id -> FeatureObject
        separate_strands: Whether to separate positive/negative strands
        resolve_overlaps: Whether to resolve overlapping features by assigning to different tracks
        genome_length: Total genome length (for origin-spanning feature detection)
    
    Returns:
        Updated feature_dict with feature_track_id set on each feature
    """
    feature_metrics = {}
    for feat_id, feature in feature_dict.items():
        total_span, occupied_length = calculate_feature_metrics(feature, genome_length)
        start, end, strand = get_feature_ends(feature, genome_length)

        occupation_ratio = occupied_length / total_span if total_span > 0 else 0

        feature_metrics[feat_id] = {
            "id": feat_id,
            "start": start,
            "end": end,
            "strand": strand,
            "total_span": total_span,
            "occupied_length": occupied_length,
            "occupation_ratio": occupation_ratio,
        }

    def sort_key(item):
        metrics = item[1]
        if separate_strands:
            return (
                0 if metrics["strand"] == "positive" else 1,
                -metrics["occupied_length"],
                metrics["start"],
            )
        else:
            return (
                -metrics["occupied_length"],
                metrics["start"],
            )

    sorted_features = sorted(feature_metrics.items(), key=sort_key)

    pos_tracks: Dict[str, List[dict]] = {}
    neg_tracks: Dict[str, List[dict]] | None = {} if separate_strands else None

    for feat_id, feat_metrics in sorted_features:
        if separate_strands:
            track_dict = neg_tracks if feat_metrics["strand"] == "negative" else pos_tracks  # type: ignore[assignment]
        else:
            track_dict = pos_tracks

        track_num = find_best_track(
            feat_metrics, track_dict, separate_strands, resolve_overlaps, genome_length
        )
        track_id = f"track_{abs(track_num)}"

        if track_id not in track_dict:
            track_dict[track_id] = []
        track_dict[track_id].append(feat_metrics)

        feature_dict[feat_id].feature_track_id = track_num

    return feature_dict


__all__ = [
    "arrange_feature_tracks",
    "calculate_feature_metrics",
    "check_feature_overlap",
    "find_best_track",
    "get_feature_ends",
]
