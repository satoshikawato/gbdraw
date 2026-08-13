#!/usr/bin/env python
# coding: utf-8

from typing import List, Union

from Bio.SeqRecord import SeqRecord

from ..features.visibility import should_render_feature
from ..features.instance_identity import FeatureInstanceIdentityPlan
from ..features.semantic_selectors import FeatureSemanticSelectorContext


def create_dict_for_sequence_lengths(records: list[SeqRecord]) -> dict[str, int]:
    return {record.id: len(record.seq) for record in records}


def determine_length_parameter(record_length: int, length_threshold: int) -> str:
    if record_length < length_threshold:
        return "short"
    return "long"


def check_feature_presence(
    records: Union[List[SeqRecord], SeqRecord],
    features_list: List[str],
    feature_visibility_rules=None,
    specific_color_rules=None,
    feature_instance_identity_plan: FeatureInstanceIdentityPlan | None = None,
    feature_semantic_selector_context: FeatureSemanticSelectorContext | None = None,
) -> list[str]:
    if isinstance(records, SeqRecord):
        records = [records]

    features_present: list[str] = []
    seen_feature_types: set[str] = set()

    for record in records:
        for feature in record.features:
            feature_instance_identity = (
                feature_instance_identity_plan.identity_for_feature(feature)
                if feature_instance_identity_plan is not None
                else None
            )
            if not should_render_feature(
                feature,
                features_list,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=specific_color_rules,
                feature_instance_identity=feature_instance_identity,
                feature_semantic_selector_context=feature_semantic_selector_context,
            ):
                continue
            if feature.type in seen_feature_types:
                continue
            seen_feature_types.add(feature.type)
            features_present.append(feature.type)
    return features_present


def get_coordinates_of_longest_segment(feature_object):
    coords = feature_object.coordinates
    if not coords:
        return None, -1

    longest_segment_info = None
    max_length = -1

    for coord in coords:
        try:
            start, end = int(coord[2]), int(coord[3])
            length = abs(end - start)
            if length > max_length:
                max_length = length
                longest_segment_info = coord
        except (IndexError, TypeError, ValueError):
            continue

    return longest_segment_info, max_length


__all__ = [
    "check_feature_presence",
    "create_dict_for_sequence_lengths",
    "determine_length_parameter",
    "get_coordinates_of_longest_segment",
]
