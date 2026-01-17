#!/usr/bin/env python
# coding: utf-8

from typing import List, Dict, Optional, Set, Tuple

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from .objects import GeneObject, RepeatObject, FeatureObject
from ..labels.filtering import get_label_text
from .colors import get_color, get_color_with_info
from .coordinates import get_exon_and_intron_coordinates
from .tracks import arrange_feature_tracks


def create_repeat_object(
    repeat_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    genome_length: int,
    label_filtering,
    record_id: Optional[str] = None,
    color: Optional[str] = None,
) -> RepeatObject:
    """
    Creates a RepeatObject representing a repeat region in a genome.
    """
    coordinates = feature.location.parts
    is_directional: bool = False
    rpt_family: str = feature.qualifiers.get("rpt_family", ["undefined"])[0]
    rpt_type: str = feature.qualifiers.get("rpt_type", ["undefined"])[0]
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(coordinates, genome_length)
    if color is None:
        color = get_color(feature, color_table, default_colors, record_id=record_id)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering)

    repeat_object = RepeatObject(
        repeat_id,
        location,
        is_directional,
        color,
        note,
        rpt_family,
        rpt_type,
        label_text,
        coordinates,
        feature_type,
        qualifiers=feature.qualifiers,
        record_id=record_id,
    )
    return repeat_object


def create_feature_object(
    feature_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    genome_length: int,
    label_filtering,
    record_id: Optional[str] = None,
    color: Optional[str] = None,
) -> FeatureObject:
    """
    Creates a FeatureObject representing a generic genomic feature.
    """
    coordinates = feature.location.parts
    is_directional: bool = False
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(coordinates, genome_length)
    if color is None:
        color = get_color(feature, color_table, default_colors, record_id=record_id)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering)

    feature_object = FeatureObject(
        feature_id,
        location,
        is_directional,
        color,
        note,
        label_text,
        coordinates,
        feature_type,
        qualifiers=feature.qualifiers,
        record_id=record_id,
    )
    return feature_object


def create_gene_object(
    feature_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    genome_length: int,
    label_filtering,
    record_id: Optional[str] = None,
    color: Optional[str] = None,
) -> GeneObject:
    """
    Creates a GeneObject representing a gene in a genome.
    """
    is_directional: bool = True
    coordinates: List[SimpleLocation] = feature.location.parts
    note: str = feature.qualifiers.get("note", [""])[0]
    product: str = feature.qualifiers.get("product", [""])[0]
    gene: str = feature.qualifiers.get("gene", [""])[0]
    is_trans_spliced = "trans_splicing" in feature.qualifiers
    feature_type = feature.type
    location = get_exon_and_intron_coordinates(coordinates, genome_length, is_trans_spliced)
    if color is None:
        color = get_color(feature, color_table, default_colors, record_id=record_id)
    label_text = get_label_text(feature, label_filtering)

    gene_object = GeneObject(
        feature_id,
        location,
        is_directional,
        color,
        note,
        product,
        feature.type,
        gene,
        label_text,
        coordinates,
        feature_type,
        qualifiers=feature.qualifiers,
        record_id=record_id,
    )
    return gene_object


def _create_feature_dict_core(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
) -> Tuple[Dict[str, FeatureObject], Set[Tuple[str, str]], Set[str]]:
    """
    Creates a dictionary mapping feature IDs to FeatureObjects from a GenBank record.

    NOTE: `color_table` / `default_colors` are expected to be preprocessed maps
    produced by `gbdraw.features.colors.preprocess_color_tables`.

    Returns:
        Tuple of (feature_dict, used_color_rules, default_used_features) where:
        - used_color_rules: set of (caption, color) tuples for rules that actually matched
        - default_used_features: set of feature types that fell back to default color
    """
    feature_dict: Dict[str, FeatureObject] = {}
    used_color_rules: Set[Tuple[str, str]] = set()
    default_used_features: Set[str] = set()
    locus_count: int = 0
    repeat_count: int = 0
    feature_count: int = 0
    genome_length: int = len(gb_record.seq)

    for feature in gb_record.features:
        if feature.type not in selected_features_set:
            continue

        # Track which color rules are actually used
        color, caption = get_color_with_info(feature, color_table, default_colors, record_id=gb_record.id)
        if caption is None:
            default_used_features.add(feature.type)
        elif caption:
            used_color_rules.add((caption, color))

        if feature.type in {"CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"}:
            locus_count += 1
            locus_id: str = "gene_" + str(locus_count).zfill(9)
            gene_object: GeneObject = create_gene_object(
                locus_id,
                feature,
                color_table,
                default_colors,
                genome_length,
                label_filtering,
                record_id=gb_record.id,
                color=color,
            )
            feature_dict[locus_id] = gene_object
        elif feature.type == "repeat_region":
            repeat_count += 1
            repeat_id: str = "crt_" + str(repeat_count).zfill(9)
            repeat_object: RepeatObject = create_repeat_object(
                repeat_id,
                feature,
                color_table,
                default_colors,
                genome_length,
                label_filtering,
                record_id=gb_record.id,
                color=color,
            )
            feature_dict[repeat_id] = repeat_object
        else:
            feature_count += 1
            feature_id: str = "feature_" + str(feature_count).zfill(9)
            feature_object: FeatureObject = create_feature_object(
                feature_id,
                feature,
                color_table,
                default_colors,
                genome_length,
                label_filtering,
                record_id=gb_record.id,
                color=color,
            )
            feature_dict[feature_id] = feature_object

    feature_dict = arrange_feature_tracks(
        feature_dict, separate_strands, resolve_overlaps, genome_length
    )
    return feature_dict, used_color_rules, default_used_features


def create_feature_dict(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
) -> Tuple[Dict[str, FeatureObject], Set[Tuple[str, str]]]:
    """
    Creates a dictionary mapping feature IDs to FeatureObjects from a GenBank record.

    NOTE: `color_table` / `default_colors` are expected to be preprocessed maps
    produced by `gbdraw.features.colors.preprocess_color_tables`.

    Returns:
        Tuple of (feature_dict, used_color_rules) where used_color_rules is a set
        of (caption, color) tuples for rules that actually matched features.
    """
    feature_dict, used_color_rules, _ = _create_feature_dict_core(
        gb_record,
        color_table,
        selected_features_set,
        default_colors,
        separate_strands,
        resolve_overlaps,
        label_filtering,
    )
    return feature_dict, used_color_rules


def create_feature_dict_with_color_usage(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
) -> Tuple[Dict[str, FeatureObject], Set[Tuple[str, str]], Set[str]]:
    """
    Variant of create_feature_dict that also returns default-used feature types.
    """
    return _create_feature_dict_core(
        gb_record,
        color_table,
        selected_features_set,
        default_colors,
        separate_strands,
        resolve_overlaps,
        label_filtering,
    )


__all__ = [
    "create_feature_dict",
    "create_feature_dict_with_color_usage",
    "create_feature_object",
    "create_gene_object",
    "create_repeat_object",
]


