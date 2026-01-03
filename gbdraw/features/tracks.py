#!/usr/bin/env python
# coding: utf-8

from typing import Dict, List, Tuple

from ..feature_objects import FeatureObject


def get_feature_ends(feature):
    strand = feature.location[0].strand if feature.location else "undefined"
    if hasattr(feature, "coordinates") and feature.coordinates:
        start = min(part.start for part in feature.coordinates)
        if start < 1:
            start = 1
        end = max(part.end for part in feature.coordinates)
        if end < 1:
            end = 1
    else:
        start = min(part.start for part in feature.location)
        end = max(part.end for part in feature.location)
        if start < 1:
            start = 1
        if end < 1:
            end = 1
    return start, end, strand


def calculate_feature_metrics(feature) -> Tuple[int, int]:
    """
    Calculate total span and occupied length with corrected calculations.

    Returns:
        (total_span, occupied_length)
    """
    if hasattr(feature, "coordinates") and feature.coordinates:
        start = min(part.start for part in feature.coordinates)
        if start < 1:
            start = 1
        end = max(part.end for part in feature.coordinates)
        if end < 1:
            end = 1
    else:
        start = min(part.start for part in feature.location)
        end = max(part.end for part in feature.location)
        if start < 1:
            start = 1
        if end < 1:
            end = 1

    total_span = abs(end - start)

    occupied_length = 0
    for part in feature.location:
        if part.kind == "block":
            occupied_length += abs(part.end - part.start)

    occupied_length = min(occupied_length, total_span)
    return total_span, occupied_length


def check_feature_overlap(a: dict, b: dict, separate_strands: bool) -> bool:
    if separate_strands and a["strand"] != b["strand"]:
        return False
    return not (a["end"] < b["start"] or a["start"] > b["end"])


def find_best_track(
    feature: dict,
    track_dict: Dict[str, List[dict]],
    separate_strands: bool,
    resolve_overlaps: bool,
    max_track: int = 100,
) -> int:
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
            for existing in track_dict[key]:
                if check_feature_overlap(feature, existing, separate_strands):
                    break
            else:
                return tn

    return track_nums[0]


def arrange_feature_tracks(
    feature_dict: Dict[str, FeatureObject], separate_strands: bool, resolve_overlaps: bool
) -> Dict[str, FeatureObject]:
    """
    Arrange features in tracks with improved strand handling and track assignment.
    """
    feature_metrics = {}
    for feat_id, feature in feature_dict.items():
        total_span, occupied_length = calculate_feature_metrics(feature)
        start, end, strand = get_feature_ends(feature)

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

        track_num = find_best_track(feat_metrics, track_dict, separate_strands, resolve_overlaps)
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


