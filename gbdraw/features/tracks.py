#!/usr/bin/env python
# coding: utf-8

from typing import Dict, List, NoReturn, Optional, Tuple

from ..exceptions import ValidationError
from ..layout.spatial import IntervalIndex, split_circular_interval
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


def feature_overlap_bp(
    a: dict, b: dict, *, separate_strands: bool = False,
    genome_length: Optional[int] = None,
) -> int:
    """Length of the shared outer envelope, using existing [1, length) wrap pieces."""
    if separate_strands and _strand_pool(a["strand"]) != _strand_pool(b["strand"]):
        return 0
    if genome_length is None and (a["start"] > a["end"] or b["start"] > b["end"]):
        raise ValidationError("Origin-spanning overlap length requires genome_length.")

    def pieces(feature):
        start, end = feature["start"], feature["end"]
        return [(start, end)] if start <= end else [(start, int(genome_length)), (1, end)]

    return sum(
        max(0, min(a_end, b_end) - max(a_start, b_start))
        for a_start, a_end in pieces(a)
        for b_start, b_end in pieces(b)
    )


def _strand_pool(strand: str) -> str:
    # The allocator's existing nominal track is 0 for positive, -1 otherwise.
    return "positive" if strand == "positive" else "negative"


def features_conflict(
    a: dict, b: dict, *, tolerance_bp: int = 0,
    separate_strands: bool = False, genome_length: Optional[int] = None,
) -> bool:
    if separate_strands and _strand_pool(a["strand"]) != _strand_pool(b["strand"]):
        return False
    if genome_length is None and (a["start"] > a["end"] or b["start"] > b["end"]):
        if tolerance_bp:
            raise ValidationError("Origin-spanning overlap tolerance requires genome_length.")
        return True  # Preserve the existing conservative no-length default.
    return feature_overlap_bp(
        a, b, separate_strands=separate_strands, genome_length=genome_length,
    ) > tolerance_bp


def check_feature_overlap(
    a: dict, b: dict, separate_strands: bool, genome_length: Optional[int] = None,
    tolerance_bp: int = 0,
) -> bool:
    """Compatibility boolean entry; all allocation uses the same final predicate."""
    return features_conflict(
        a, b, tolerance_bp=tolerance_bp,
        separate_strands=separate_strands, genome_length=genome_length,
    )


def _raise_feature_track_limit(feature: dict, max_track: int) -> NoReturn:
    feature_id = feature.get("id", "<unknown>")
    raise ValidationError(
        f"Unable to place feature {feature_id!r} without overlap: "
        f"all {max_track} feature lanes are occupied."
    )


def find_best_track(
    feature: dict,
    track_dict: Dict[str, List[dict]],
    separate_strands: bool,
    resolve_overlaps: bool,
    genome_length: Optional[int] = None,
    max_track: int = 100,
    tolerance_bp: int = 0,
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
                if check_feature_overlap(feature, existing, separate_strands, genome_length, tolerance_bp):
                    has_overlap = True
                    break
            if not has_overlap:
                return tn

    if resolve_overlaps:
        _raise_feature_track_limit(feature, max_track)
    return track_nums[0]


def _feature_interval_bucket_size(genome_length: Optional[int]) -> float:
    if genome_length is None or int(genome_length) <= 0:
        return 1000.0
    return max(1.0, float(genome_length) / 512.0)


def _feature_candidates_for_track(
    feature: dict,
    track_key: str,
    track_dict: Dict[str, List[dict]],
    track_indexes: dict[str, IntervalIndex],
    feature_by_id: dict[str, dict],
    genome_length: Optional[int],
) -> list[dict]:
    existing_features = track_dict.get(track_key, [])
    if not existing_features:
        return []

    index = track_indexes.get(track_key)
    if index is None:
        return existing_features

    seen: set[str] = set()
    candidate_ids: list[str] = []
    for interval in split_circular_interval(feature["start"], feature["end"], genome_length):
        for candidate_id in index.query(interval):
            feature_id = str(candidate_id)
            if feature_id in seen:
                continue
            seen.add(feature_id)
            candidate_ids.append(feature_id)
    return [feature_by_id[feature_id] for feature_id in candidate_ids]


def _insert_feature_track_index(
    track_indexes: dict[str, IntervalIndex],
    track_key: str,
    feature: dict,
    genome_length: Optional[int],
    bucket_size: float,
) -> None:
    index = track_indexes.setdefault(track_key, IntervalIndex(bucket_size=bucket_size))
    for interval in split_circular_interval(feature["start"], feature["end"], genome_length):
        index.insert(str(feature["id"]), interval)


def _find_best_track_indexed(
    feature: dict,
    track_dict: Dict[str, List[dict]],
    track_indexes: dict[str, IntervalIndex],
    feature_by_id: dict[str, dict],
    separate_strands: bool,
    resolve_overlaps: bool,
    genome_length: Optional[int] = None,
    max_track: int = 100,
    tolerance_bp: int = 0,
) -> int:
    """Find a feature track using index candidates and exact overlap checks."""
    if not separate_strands:
        track_nums = [0] if not resolve_overlaps else list(range(0, max_track))
    else:
        if not resolve_overlaps:
            track_nums = [0] if feature["strand"] == "positive" else [-1]
        elif feature["strand"] == "positive":
            track_nums = list(range(0, max_track))
        else:
            track_nums = list(range(-1, -max_track - 1, -1))

    if resolve_overlaps:
        for tn in track_nums:
            key = f"track_{abs(tn)}"
            if key not in track_dict or not track_dict[key]:
                return tn
            candidates = _feature_candidates_for_track(
                feature,
                key,
                track_dict,
                track_indexes,
                feature_by_id,
                genome_length,
            )
            if not any(check_feature_overlap(feature, existing, separate_strands, genome_length, tolerance_bp) for existing in candidates):
                return tn

    if resolve_overlaps:
        _raise_feature_track_limit(feature, max_track)
    return track_nums[0]




def _find_best_track_split_overlaps_by_strand_indexed(
    feature: dict,
    center_track: List[dict],
    center_track_index: IntervalIndex,
    outer_tracks: Dict[str, List[dict]],
    outer_track_indexes: dict[str, IntervalIndex],
    inner_tracks: Dict[str, List[dict]],
    inner_track_indexes: dict[str, IntervalIndex],
    feature_by_id: dict[str, dict],
    genome_length: Optional[int] = None,
    max_track: int = 100,
    tolerance_bp: int = 0,
) -> int:
    """Find split inner/outer track using index candidates and exact checks."""
    center_candidates = _feature_candidates_for_track(
        feature,
        "track_0",
        {"track_0": center_track},
        {"track_0": center_track_index},
        feature_by_id,
        genome_length,
    )
    if not any(check_feature_overlap(feature, existing, False, genome_length, tolerance_bp) for existing in center_candidates):
        return 0

    is_negative = feature["strand"] == "negative"
    track_dict = inner_tracks if is_negative else outer_tracks
    track_indexes = inner_track_indexes if is_negative else outer_track_indexes
    sign = -1 if is_negative else 1

    for track_index in range(1, max_track):
        key = f"track_{track_index}"
        if key not in track_dict or not track_dict[key]:
            return sign * track_index
        candidates = _feature_candidates_for_track(
            feature,
            key,
            track_dict,
            track_indexes,
            feature_by_id,
            genome_length,
        )
        if not any(check_feature_overlap(feature, existing, False, genome_length, tolerance_bp) for existing in candidates):
            return sign * track_index

    _raise_feature_track_limit(feature, max_track)


def arrange_feature_tracks(
    feature_dict: Dict[str, FeatureObject],
    separate_strands: bool,
    resolve_overlaps: bool,
    split_overlaps_by_strand: bool = False,
    genome_length: Optional[int] = None,
    *,
    fixed_tracks: dict[str, int] | None = None,
    tolerance_bp: int = 0,
) -> Dict[str, FeatureObject]:
    """
    Arrange features in tracks with improved strand handling and track assignment.

    Args:
        feature_dict: Dict of feature_id -> FeatureObject
        separate_strands: Whether to separate positive/negative strands
        resolve_overlaps: Whether to resolve overlapping features by assigning to different tracks
        split_overlaps_by_strand: When True (and only when separate_strands is False),
            uses shared center track 0 plus strand-specific displacement pools.
            Intended for circular middle resolve_overlaps behavior.
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

    # Occupancy is partitioned by physical lane, never by a raw strand third pool.
    split = bool(split_overlaps_by_strand) and not separate_strands
    pos_tracks: Dict[str, List[dict]] = {}
    neg_tracks: Dict[str, List[dict]] = {}
    pos_indexes: dict[str, IntervalIndex] = {}
    neg_indexes: dict[str, IntervalIndex] = {}
    center_track: List[dict] = []
    center_index = IntervalIndex(bucket_size=_feature_interval_bucket_size(genome_length))
    bucket_size = _feature_interval_bucket_size(genome_length)

    def occupancy(track_num):
        if split and track_num == 0:
            return {"track_0": center_track}, {"track_0": center_index}, "track_0"
        if track_num < 0:
            return neg_tracks, neg_indexes, f"track_{abs(track_num)}"
        return pos_tracks, pos_indexes, f"track_{track_num}"

    def reserve(feat_id, track_num, *, fixed=False):
        metrics = feature_metrics[feat_id]
        tracks, indexes, key = occupancy(track_num)
        if fixed:
            candidates = _feature_candidates_for_track(
                metrics, key, tracks, indexes, feature_metrics, genome_length,
            )
            if any(features_conflict(
                metrics, other, tolerance_bp=tolerance_bp, genome_length=genome_length,
            ) for other in candidates):
                raise ValidationError(
                    f"Fixed feature placement conflict for {feat_id!r} on lane {track_num}."
                )
        tracks.setdefault(key, []).append(metrics)
        _insert_feature_track_index(indexes, key, metrics, genome_length, bucket_size)
        feature_dict[feat_id].feature_track_id = track_num

    fixed_tracks = fixed_tracks or {}
    for feat_id, track_num in fixed_tracks.items():
        reserve(feat_id, track_num, fixed=True)

    for feat_id, feat_metrics in sorted_features:
        if feat_id in fixed_tracks:
            continue
        if split and resolve_overlaps:
            track_num = _find_best_track_split_overlaps_by_strand_indexed(
                feat_metrics, center_track, center_index,
                pos_tracks, pos_indexes, neg_tracks, neg_indexes,
                feature_metrics, genome_length=genome_length, tolerance_bp=tolerance_bp,
            )
        else:
            negative = separate_strands and _strand_pool(feat_metrics["strand"]) == "negative"
            track_num = _find_best_track_indexed(
                feat_metrics,
                neg_tracks if negative else pos_tracks,
                neg_indexes if negative else pos_indexes,
                feature_metrics, separate_strands, resolve_overlaps,
                genome_length, tolerance_bp=tolerance_bp,
            )
        reserve(feat_id, track_num)
    return feature_dict


__all__ = [
    "arrange_feature_tracks",
    "calculate_feature_metrics",
    "check_feature_overlap",
    "feature_overlap_bp",
    "features_conflict",
    "find_best_track",
    "get_feature_ends",
]
