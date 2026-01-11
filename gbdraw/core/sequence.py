#!/usr/bin/env python
# coding: utf-8

from typing import List, Union

from Bio.SeqRecord import SeqRecord


def create_dict_for_sequence_lengths(records: list[SeqRecord]) -> dict[str, int]:
    return {record.id: len(record.seq) for record in records}


def determine_length_parameter(record_length: int, length_threshold: int) -> str:
    if record_length < length_threshold:
        return "short"
    return "long"


def determine_output_file_prefix(gb_records, output_prefix, record_count, accession):
    if len(gb_records) > 1 and output_prefix is not None:
        return "{}_{}".format(output_prefix, record_count)
    elif len(gb_records) == 1 and output_prefix is not None:
        return output_prefix
    else:
        return accession


def check_feature_presence(
    records: Union[List[SeqRecord], SeqRecord], features_list: List[str]
) -> dict:
    if isinstance(records, SeqRecord):
        records = [records]

    features_to_be_checked = set(features_list)
    features_present: list[str] = []

    for record in records:
        for feature in record.features:
            if feature.type in features_to_be_checked:
                features_present.append(feature.type)
                features_to_be_checked.remove(feature.type)
                if not features_to_be_checked:
                    break
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
    "determine_output_file_prefix",
    "get_coordinates_of_longest_segment",
]


