#!/usr/bin/env python
# coding: utf-8
import re
import sys
import logging
from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from typing import List, Dict, Tuple
from .feature_objects import GeneObject, RepeatObject, FeatureObject
from .utility_functions import get_label_text
# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def get_color(feature: SeqFeature, color_table: DataFrame, default_colors: DataFrame) -> str:
    """
    Determines the color for a genomic feature based on a color table and default colors.

    Args:
        feature (SeqFeature): Genomic feature object containing qualifiers.
        color_table (DataFrame): Table mapping feature types and qualifiers to specific colors.
        default_colors (DataFrame): Table containing default colors for feature types.

    Returns:
        str: Hex color code for the feature.

    This function searches the color table for a color matching the feature's type and qualifiers.
    If no specific match is found, it uses the default color for the feature's type. A fallback color
    is returned if the feature type is not in the default color table.
    """
    # Check each row in color_table against feature's qualifiers
    if isinstance(color_table, DataFrame):
        for _, row in color_table.iterrows():
            if row['feature_type'] != feature.type:
                continue
            qualifier_key: str = row['qualifier_key']
            qualifier_value: str = row['value']
            if qualifier_key in feature.qualifiers:
                for value in feature.qualifiers[qualifier_key]:
                    if re.search(qualifier_value, value, re.IGNORECASE):
                        return row['color']

    # Default color if no specific match is found
    default_color_row = default_colors[default_colors['feature_type']
                                       == feature.type]
    if not default_color_row.empty:
        return default_color_row.iloc[0]['color']

    return "#d3d3d3"  # Fallback color


def create_repeat_object(repeat_id: str, feature: SeqFeature, color_table: DataFrame, default_colors: DataFrame, genome_length: int, label_filtering) -> RepeatObject:
    """
    Creates a RepeatObject representing a repeat region in a genome.

    Args:
        repeat_id (str): Identifier for the repeat region.
        feature (SeqFeature): BioPython SeqFeature object representing the repeat region.
        color_table (DataFrame): Table mapping feature types and qualifiers to specific colors.
        default_colors (DataFrame): Table containing default colors for feature types.
        genome_length (int): Length of the genome sequence.

    Returns:
        RepeatObject: An object encapsulating data about the repeat region.

    This function constructs a RepeatObject by extracting and processing information from the feature,
    including its location, directional status, color, and other properties such as repeat family and type.
    """
    coordinates = feature.location.parts
    is_directional: bool = False
    rpt_family: str = feature.qualifiers.get('rpt_family', ["undefined"])[0]
    rpt_type: str = feature.qualifiers.get('rpt_type', ["undefined"])[0]
    note: str = feature.qualifiers.get('note', [""])
    location: list[Tuple[str, str, str, int, int, bool]
                   ] = get_exon_and_intron_coordinates(coordinates, genome_length)
    color: str = get_color(feature, color_table, default_colors)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering)

    repeat_object = RepeatObject(
        repeat_id, location, is_directional, color, note, rpt_family, rpt_type, label_text, coordinates, feature_type)
    return repeat_object


def create_feature_object(feature_id: str, feature: SeqFeature, color_table: DataFrame, default_colors: DataFrame, genome_length: int, label_filtering) -> FeatureObject:
    """
    Creates a FeatureObject representing a generic genomic feature.

    Args:
        feature_id (str): Identifier for the feature.
        feature (SeqFeature): BioPython SeqFeature object representing the genomic feature.
        color_table (DataFrame): Table mapping feature types and qualifiers to specific colors.
        default_colors (DataFrame): Table containing default colors for feature types.
        genome_length (int): Length of the genome sequence.

    Returns:
        FeatureObject: An object encapsulating data about the genomic feature.

    This function constructs a FeatureObject by extracting information from the feature,
    including its location, directional status, color, and additional notes.
    """
    coordinates = feature.location.parts
    is_directional: bool = False
    note: str = feature.qualifiers.get('note', [""])
    location: list[Tuple[str, str, str, int, int, bool]
                   ] = get_exon_and_intron_coordinates(coordinates, genome_length)
    color: str = get_color(feature, color_table, default_colors)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering)
    feature_object = FeatureObject(
        feature_id, location, is_directional, color, note, label_text, coordinates, feature_type)
    return feature_object


def create_gene_object(feature_id: str, feature: SeqFeature, color_table: DataFrame, default_colors: DataFrame, genome_length: int, label_filtering) -> GeneObject:
    """
    Creates a GeneObject representing a gene in a genome.

    Args:
        feature_id (str): Identifier for the gene.
        feature (SeqFeature): BioPython SeqFeature object representing the gene.
        color_table (DataFrame): Table mapping feature types and qualifiers to specific colors.
        default_colors (DataFrame): Table containing default colors for feature types.
        genome_length (int): Length of the genome sequence.

    Returns:
        GeneObject: An object encapsulating data about the gene.

    This function constructs a GeneObject by extracting and processing information from the feature,
    including its location, directional status, color, and other properties such as product and gene biotype.
    """
    is_directional: bool = True
    coordinates: List[SimpleLocation] = feature.location.parts
    note: str = feature.qualifiers.get('note', [""])[0]
    product: str = feature.qualifiers.get('product', [""])[0]
    gene: str = feature.qualifiers.get('gene', [""])[0]
    is_trans_spliced = 'trans_splicing' in feature.qualifiers
    feature_type = feature.type
    location: list[Tuple[str, str, str, int, int, bool]] = get_exon_and_intron_coordinates(coordinates, genome_length, is_trans_spliced)
    color: str = get_color(feature, color_table, default_colors)
    label_text = get_label_text(feature, label_filtering)
    gene_object = GeneObject(
        feature_id, location, is_directional, color, note, product, feature.type, gene, label_text, coordinates, feature_type)
    return gene_object


def get_feature_ends(feature):
    strand = feature.location[0][2]
    if hasattr(feature, 'coordinates') and feature.coordinates:
        start = min(part.start for part in feature.coordinates)
        if start < 1:
            start = 1
        end = max(part.end for part in feature.coordinates)
        if end < 1:
            end = 1
    else:
        start = min(part[3] for part in feature.location)
        end = max(part[4] for part in feature.location)
        if start < 1:
            start = 1
        if end < 1:
            end = 1
    return start, end, strand


def calculate_feature_metrics(feature) -> Tuple[int, int]:
    """
    Calculate total span and occupied length with corrected calculations
    
    Returns:
        Tuple[int, int]: (total_span, occupied_length)
    """
    # Get the full span first
    if hasattr(feature, 'coordinates') and feature.coordinates:
        start = min(part.start for part in feature.coordinates)
        if start < 1:
            start = 1
        end = max(part.end for part in feature.coordinates)
        if end < 1:
            end = 1
    else:
        start = min(part[3] for part in feature.location)
        end = max(part[4] for part in feature.location)
        if start < 1:
            start = 1
        if end < 1:
            end = 1

    
    total_span = abs(end - start)
    
    # Calculate occupied length (only blocks/exons)
    occupied_length = 0
    if hasattr(feature, 'coordinates') and feature.coordinates:
        for part in feature.location:
            if part[0] == "block":
                occupied_length += abs(part[4] - part[3])
    else:
        for part in feature.location:
            if part[0] == "block":
                occupied_length += abs(part[4] - part[3])
    
    # Ensure occupied length cannot exceed total span
    occupied_length = min(occupied_length, total_span)
    
    return total_span, occupied_length

def check_feature_overlap(feature1: dict, feature2: dict, separate_strands: bool = False) -> bool:
    """
    Check if two features overlap, considering strands only when separate_strands is True
    
    Args:
        feature1, feature2 (dict): Feature dictionaries containing start, end and strand positions
        separate_strands (bool): If True, features on different strands don't count as overlapping
        
    Returns:
        bool: True if features overlap (considering strand settings), False otherwise
    """
    # Only consider strand differences if we're separating strands
    if separate_strands and feature1["strand"] != feature2["strand"]:
        return False
    
    # Normalize coordinates to ensure start < end
    def normalize_coords(start: int, end: int) -> tuple[int, int]:
        start = max(1, start)  # Ensure coordinates are at least 1
        end = max(1, end)
        return (min(start, end), max(start, end))
    
    f1_start, f1_end = normalize_coords(feature1["start"], feature1["end"])
    f2_start, f2_end = normalize_coords(feature2["start"], feature2["end"])
    
    # Check for any overlap including nesting
    overlap_conditions = [
        # Regular overlap from either side
        (f1_start <= f2_start <= f1_end),
        (f1_start <= f2_end <= f1_end),
        (f2_start <= f1_start <= f2_end),
        (f2_start <= f1_end <= f2_end),
        
        # Complete nesting
        (f1_start <= f2_start and f1_end >= f2_end),
        (f2_start <= f1_start and f2_end >= f1_end)
    ]
    
    return any(overlap_conditions)
def check_feature_overlap(a: dict, b: dict, separate_strands: bool) -> bool:

    if separate_strands and a["strand"] != b["strand"]:
        return False
    return not (a["end"] < b["start"] or a["start"] > b["end"])

def find_best_track(
    feature: dict,
    track_dict: Dict[str, List[dict]],
    separate_strands: bool,
    resolve_overlaps: bool,
    max_track: int = 100
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
                track_nums = list(range(-1, -max_track-1, -1))

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

def arrange_feature_tracks(feature_dict: Dict[str, FeatureObject], separate_strands: bool, resolve_overlaps: bool) -> Dict[str, FeatureObject]:
    """
    Arrange features in tracks with improved strand handling and track assignment
    """
    # Calculate metrics for all features
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
            "occupation_ratio": occupation_ratio
        }

    # Sorting strategy - when not separating strands, ignore strand in sorting
    def sort_key(item):
        metrics = item[1]
        if separate_strands:
            return (
                0 if metrics["strand"] == "positive" else 1,  # Group by strand
                -metrics["occupied_length"],   # Larger features first
                metrics["start"]               # Left to right within strand
            )
        else:
            return (
                -metrics["occupied_length"],   # Larger features first
                metrics["start"]               # Left to right
            )

    sorted_features = sorted(feature_metrics.items(), key=sort_key)
    
    # Track dictionaries - only use pos_tracks when not separating strands
    pos_tracks = {}
    neg_tracks = {} if separate_strands else None
    
    # Place features
    for feat_id, feat_metrics in sorted_features:
        if separate_strands:
            # Use separate track dictionaries for each strand
            track_dict = neg_tracks if feat_metrics["strand"] == "negative" else pos_tracks
        else:
            # Use single track dictionary for all features
            track_dict = pos_tracks
        
        track_num = find_best_track(feat_metrics, track_dict, separate_strands, resolve_overlaps)
        track_id = f"track_{abs(track_num)}"  # Use absolute track number for dictionary key
        
        if track_id not in track_dict:
            track_dict[track_id] = []
        track_dict[track_id].append(feat_metrics)
        
        feature_dict[feat_id].feature_track_id = track_num
    
    return feature_dict

def create_feature_dict(gb_record: SeqRecord, color_table: DataFrame, selected_features_set: List[str], default_colors: DataFrame, separate_strands: bool, resolve_overlaps:bool, label_filtering) -> Dict[str, FeatureObject]:
    """
    Creates a dictionary mapping feature IDs to FeatureObjects from a GenBank record.

    Args:
        gb_record (SeqRecord): BioPython SeqRecord object representing a genomic sequence.
        color_table (DataFrame): Table mapping feature types and qualifiers to specific colors.
        selected_features_set (List[str]): Set of feature types to include.
        default_colors (DataFrame): Table containing default colors for feature types.

    Returns:
        Dict[str, FeatureObject]: A dictionary where keys are feature IDs and values are FeatureObjects.

    This function iterates over features in the GenBank record, creating appropriate FeatureObjects (including
    GeneObject and RepeatObject) based on the feature type and adding them to the dictionary.
    """
    feature_dict: dict = {}
    locus_count: int = 0
    repeat_count: int = 0
    feature_count: int = 0
    genome_length: int = len(gb_record.seq)
    separate_strands: bool = separate_strands
    resolve_overlaps: bool = resolve_overlaps
    label_filtering = label_filtering
    for feature in gb_record.features:
        if feature.type not in selected_features_set:
            continue
        else:
            if (feature.type == 'CDS') or (
                    feature.type == 'rRNA') or (feature.type == 'tRNA') or (feature.type == 'tmRNA') or  (feature.type == 'ncRNA') or (feature.type == 'misc_RNA'):
                locus_count: int = locus_count + 1
                locus_id: str = "gene_" + str(locus_count).zfill(9)
                gene_object: GeneObject = create_gene_object(
                    locus_id, feature, color_table, default_colors, genome_length, label_filtering)
                feature_dict[locus_id] = gene_object
            elif feature.type == 'repeat_region':
                repeat_count: int = repeat_count + 1
                repeat_id: str = "crt_" + str(repeat_count).zfill(9)
                repeat_object: RepeatObject = create_repeat_object(
                    repeat_id, feature, color_table, default_colors, genome_length, label_filtering)
                feature_dict[repeat_id] = repeat_object
            else:
                feature_count: int = feature_count + 1
                feature_id: str = "feature_" + str(feature_count).zfill(9)
                feature_object: FeatureObject = create_feature_object(
                    feature_id, feature, color_table, default_colors, genome_length, label_filtering)
                feature_dict[feature_id] = feature_object
    feature_dict = arrange_feature_tracks(feature_dict, separate_strands, resolve_overlaps)
    if locus_count == 0:
        logger.warning(f"WARNING: No genes were found in {gb_record.id}. Are you sure the GenBank file is in the correct format?")
    return feature_dict


def get_exon_coordinate(exon_line: SimpleLocation, previous_exon_count: int, last_or_not: bool) -> Tuple[int, Tuple[str, str, str, int, int, bool]]:
    """
    Calculates the coordinate of an exon.

    Args:
        exon_line (SimpleLocation): The location of the exon within the genome.
        previous_exon_count (int): Count of exons processed before this one.
        last_or_not (bool): Indicator if this is the last exon.

    Returns:
        Tuple[int, Tuple[str, str, str, int, int, bool]]: A tuple containing the updated exon count and 
        a tuple with exon details (type, id, strand, start, end, and if it's the last exon).

    This function generates a unique ID for the exon and calculates its start and end coordinates, along with 
    the strand information. It returns these details in a structured format.
    """
    exon_strand: str = get_strand(exon_line.strand)  # type: ignore
    exon_count: int = previous_exon_count + 1
    exon_id: str = str(exon_count).zfill(3)
    exon_start = int(exon_line.start)  # type: ignore
    if exon_start < 1:
        exon_start = 1

    exon_end = int(exon_line.end)  # type: ignore
    if exon_end < 1:
        exon_end = 1
    exon_coordinate: Tuple[str, str, str, int, int, bool] = (
        "block",
        exon_id,
        exon_strand,
        exon_start,
        exon_end,
        last_or_not)
    return exon_count, exon_coordinate


def get_intron_coordinate(previous_exon: Tuple[str, str, str, int, int, bool], current_exon: Tuple[str, str, str, int, int, bool], previous_intron_count: int, genome_length: int) -> tuple[int, list[Tuple[str, str, str, int, int, bool]]]:
    """
    Generates coordinates for introns based on adjacent exons.

    Args:
        previous_exon (Tuple[str, str, str, int, int, bool]): Details of the previous exon.
        current_exon (Tuple[str, str, str, int, int, bool]): Details of the current exon.
        previous_intron_count (int): Count of introns processed before this one.
        genome_length (int): Total length of the genome.

    Returns:
        tuple[int, list[Tuple[str, str, str, int, int, bool]]]: A tuple containing the updated intron count 
        and a list of tuples with intron details (type, id, strand, start, end, and a boolean flag).

    This function calculates the start and end coordinates of introns based on the end of the previous exon 
    and the start of the current exon. Handles cases for different strand orientations and circular genomes.
    """
    intron_count: int = previous_intron_count + 1
    intron_id: str = str(intron_count).zfill(3)
    intron_strand: str = previous_exon[2]
    previous_exon_start: int = previous_exon[3]
    previous_exon_end: int = previous_exon[4]
    current_exon_start: int = current_exon[3]
    current_exon_end: int = current_exon[4]
    intron_parts: List[Tuple[str, str, str, int, int, bool]] = []
    if (previous_exon_end > current_exon_start + 1 and intron_strand == "positive") or (previous_exon_end + 1 < current_exon_start and intron_strand == "negative"):
        if intron_strand == "positive":
            intron_start_1: int = int(previous_exon_end) + 1
            intron_end_1: int = genome_length
            intron_coordinate_1: Tuple[str, str, str, int, int, bool] = (
                "line", intron_id, intron_strand, intron_start_1, intron_end_1, False)
            intron_parts.append(intron_coordinate_1)
            intron_start_2: int = 0
            intron_end_2: int = int(current_exon_start) - 1
            intron_coordinate_2 = (
                "line", intron_id, intron_strand, intron_start_2, intron_end_2, False)
            intron_parts.append(intron_coordinate_2)
        elif intron_strand == "negative":
            intron_start_1: int = int(current_exon_end) + 1
            intron_end_1: int = 0
            intron_coordinate_1: Tuple[str, str, str, int, int, bool] = (
                "line", intron_id, intron_strand, intron_start_1, intron_end_1, False)
            intron_parts.append(intron_coordinate_1)
            intron_start_2: int = genome_length
            intron_end_2: int = int(current_exon_end) + 1
            intron_coordinate_2: Tuple[str, str, str, int, int, bool] = (
                "line", intron_id, intron_strand, intron_start_2, intron_end_2, False)
            intron_parts.append(intron_coordinate_2)
    else:
        intron_start: int = 0
        intron_end: int = 0
        if intron_strand == "positive":
            intron_start: int = int(previous_exon_end) + 1
            intron_end: int = int(current_exon_start) - 1
        elif intron_strand == "negative":
            intron_start = int(current_exon_end) + 1
            intron_end = int(previous_exon_start) - 1
        else:
            # Handle the case where intron_strand is neither "positive" nor "negative"
            logger.error(f"Undefined intron_strand: {intron_strand}")
        intron_coordinate: Tuple[str, str, str, int, int, bool] = (
            "line",
            intron_id,
            intron_strand,
            intron_start,
            intron_end,
            False)
        intron_parts.append(intron_coordinate)
    return intron_count, intron_parts


def get_exon_and_intron_coordinates(exons: List[SimpleLocation], genome_length: int, is_trans_spliced: bool = False) -> list[Tuple[str, str, str, int, int, bool]]:
    """
    Processes a list of exons to determine the coordinates of both exons and introns.
    Handles trans-spliced features by drawing only exons without introns.
    """
    coordinates_list: list[Tuple[str, str, str, int, int, bool]] = []
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
            is_last_exon = (i == num_exons - 1)
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last_exon)

            if i > 0:
                previous_exon_coord = coordinates_list[-1]
                intron_count, intron_coords = get_intron_coordinate(
                    previous_exon_coord, exon_coord, intron_count, genome_length)
                coordinates_list.extend(intron_coords)
            
            coordinates_list.append(exon_coord)


    elif num_exons > 1 and is_trans_spliced:

        for i, exon in enumerate(exons):
            is_last = (i == num_exons - 1)
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last)
            coordinates_list.append(exon_coord)
            
    return coordinates_list


def set_arrow_shoulder(feat_strand: str, arrow_end: float, cds_arrow_length: float) -> float:
    """
    Calculates the shoulder position for an arrow representing a genomic feature.

    Args:
        feat_strand (str): The strand of the feature ('positive' or 'negative').
        arrow_end (float): The end position of the feature.
        cds_arrow_length (float): The length of the arrow representing the feature.

    Returns:
        float: The calculated shoulder position of the arrow.

    This function is used in the visualization of genomic features, particularly for determining where 
    the arrow (representing a feature like a gene) should start or end on the genome diagram.
    """
    if feat_strand == "positive":
        shoulder = float(arrow_end - cds_arrow_length)
    else:
        shoulder = float(arrow_end + cds_arrow_length)
    return shoulder


def get_coordinate(coord) -> Tuple[str, str, int, int]:
    """
    Extracts and returns basic feature coordinates.

    Args:
        coord (Tuple[str, str, int, int]): A tuple containing feature type, strand, start, and end positions.

    Returns:
        Tuple[str, str, int, int]: The same tuple passed as an argument.

    A utility function primarily used for unpacking and repacking feature coordinates. It's helpful for clarity 
    and consistency in handling feature data.
    """
    feat_type: str = coord[0]
    feat_strand: str = coord[2]
    feat_start: int = coord[3]
    feat_end: int = coord[4]
    return feat_type, feat_strand, feat_start, feat_end


def get_strand(strand_value: int) -> str:
    """
    Converts strand value to a string representation.

    Args:
        strand_value (int): Numeric representation of the strand (1 for positive, -1 for negative).

    Returns:
        str: String representation of the strand ('positive', 'negative', or 'undefined').

    This function translates the numeric strand representation used in bioinformatics into a more 
    human-readable string format.
    """
    if strand_value == 1:
        strand: str = "positive"
    elif strand_value == -1:
        strand: str = "negative"
    else:
        strand: str = "undefined"
    return strand
