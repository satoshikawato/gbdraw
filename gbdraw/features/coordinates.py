#!/usr/bin/env python
# coding: utf-8

import logging
from typing import List, Tuple

from Bio.SeqFeature import SimpleLocation

from .objects import FeatureLocation, FeatureLocationPart, Strand


logger = logging.getLogger(__name__)


def get_exon_coordinate(
    exon_line: SimpleLocation, previous_exon_count: int, last_or_not: bool
) -> Tuple[int, FeatureLocationPart]:
    """
    Calculates the coordinate of an exon.

    Returns:
        (exon_count, exon_coordinate)
    """
    exon_strand: Strand = get_strand(exon_line.strand)  # type: ignore
    exon_count: int = previous_exon_count + 1
    exon_id: str = str(exon_count).zfill(3)
    exon_start = int(exon_line.start)  # type: ignore
    if exon_start < 1:
        exon_start = 1

    exon_end = int(exon_line.end)  # type: ignore
    if exon_end < 1:
        exon_end = 1
    exon_coordinate = FeatureLocationPart("block", exon_id, exon_strand, exon_start, exon_end, last_or_not)
    return exon_count, exon_coordinate


def get_intron_coordinate(
    previous_exon: FeatureLocationPart,
    current_exon: FeatureLocationPart,
    previous_intron_count: int,
    genome_length: int,
) -> tuple[int, list[FeatureLocationPart]]:
    """
    Generates coordinates for introns based on adjacent exons.
    Returns an empty list if exons overlap or if coordinates would be invalid
    (e.g., RNA editing frameshift features with coordinate fallback).
    """
    intron_count: int = previous_intron_count + 1
    intron_id: str = str(intron_count).zfill(3)
    intron_strand: Strand = previous_exon.strand
    previous_exon_start: int = previous_exon.start
    previous_exon_end: int = previous_exon.end
    current_exon_start: int = current_exon.start
    current_exon_end: int = current_exon.end
    intron_parts: List[FeatureLocationPart] = []

    # Check for overlapping exons or invalid intron coordinates
    # This handles cases like RNA editing frameshift (e.g., join(1850..2340,2339..3026))
    # where the second exon starts before or at the first exon's end
    if intron_strand == "positive" and previous_exon_end >= current_exon_start:
        # Exons overlap or are adjacent with no gap - no intron to draw
        return intron_count, intron_parts
    elif intron_strand == "negative" and current_exon_end >= previous_exon_start:
        # Exons overlap or are adjacent with no gap - no intron to draw
        return intron_count, intron_parts

    if (intron_strand == "positive" and previous_exon_end > current_exon_start + 1) or (
        intron_strand == "negative" and previous_exon_end + 1 < current_exon_start
    ):
        if intron_strand == "positive":
            intron_start_1: int = min(int(previous_exon_end) + 1, genome_length)
            intron_end_1: int = genome_length
            intron_coordinate_1 = FeatureLocationPart(
                "line", intron_id, intron_strand, intron_start_1, intron_end_1, False
            )
            intron_parts.append(intron_coordinate_1)
            intron_start_2: int = 0
            intron_end_2: int = int(current_exon_start) - 1
            intron_coordinate_2 = FeatureLocationPart(
                "line", intron_id, intron_strand, intron_start_2, intron_end_2, False
            )
            intron_parts.append(intron_coordinate_2)
        elif intron_strand == "negative":
            intron_start_1: int = int(previous_exon_start) - 1
            intron_end_1: int = 0
            intron_coordinate_1 = FeatureLocationPart(
                "line", intron_id, intron_strand, intron_start_1, intron_end_1, False
            )
            intron_parts.append(intron_coordinate_1)
            intron_start_2: int = genome_length
            intron_end_2: int = min(int(current_exon_end) + 1, genome_length)
            intron_coordinate_2 = FeatureLocationPart(
                "line", intron_id, intron_strand, intron_start_2, intron_end_2, False
            )
            intron_parts.append(intron_coordinate_2)
    else:
        intron_start: int = 0
        intron_end: int = 0
        if intron_strand == "positive":
            intron_start = int(previous_exon_end) + 1
            intron_end = int(current_exon_start) - 1
        elif intron_strand == "negative":
            intron_start = int(current_exon_end) + 1
            intron_end = int(previous_exon_start) - 1
        else:
            logger.error(f"Undefined intron_strand: {intron_strand}")
        intron_coordinate = FeatureLocationPart(
            "line", intron_id, intron_strand, intron_start, intron_end, False
        )
        intron_parts.append(intron_coordinate)
    return intron_count, intron_parts


def get_exon_and_intron_coordinates(
    exons: List[SimpleLocation], genome_length: int, is_trans_spliced: bool = False
) -> FeatureLocation:
    """
    Processes a list of exons to determine the coordinates of both exons and introns.
    Handles trans-spliced features by drawing only exons without introns.
    """
    coordinates_list: FeatureLocation = []
    exon_count: int = 0
    intron_count: int = 0
    num_exons: int = len(exons)

    if num_exons == 0:
        return []

    if num_exons == 1:
        exon_count, exon_coord = get_exon_coordinate(exons[0], exon_count, True)
        coordinates_list.append(exon_coord)

    elif num_exons > 1 and not is_trans_spliced:
        for i, exon in enumerate(exons):
            is_last_exon = i == num_exons - 1
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last_exon)

            if i > 0:
                previous_exon_coord = coordinates_list[-1]
                intron_count, intron_coords = get_intron_coordinate(
                    previous_exon_coord, exon_coord, intron_count, genome_length
                )
                coordinates_list.extend(intron_coords)

            coordinates_list.append(exon_coord)

    elif num_exons > 1 and is_trans_spliced:
        for i, exon in enumerate(exons):
            is_last = i == num_exons - 1
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last)
            coordinates_list.append(exon_coord)

    return coordinates_list


def get_coordinate(coord: FeatureLocationPart) -> Tuple[str, Strand, int, int]:
    """
    Extracts and returns basic feature coordinates.
    """
    feat_type: str = coord.kind
    feat_strand: Strand = coord.strand
    feat_start: int = coord.start
    feat_end: int = coord.end
    return feat_type, feat_strand, feat_start, feat_end


def get_strand(strand_value: int) -> Strand:
    """
    Converts strand value to a string representation.
    """
    if strand_value == 1:
        strand: str = "positive"
    elif strand_value == -1:
        strand = "negative"
    else:
        strand = "undefined"
    return strand


__all__ = [
    "get_coordinate",
    "get_exon_and_intron_coordinates",
    "get_exon_coordinate",
    "get_intron_coordinate",
    "get_strand",
]


