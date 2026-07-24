#!/usr/bin/env python
# coding: utf-8

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from .objects import GeneObject, RepeatObject, FeatureObject
from .visibility import should_render_feature
from ..labels.filtering import get_label_text
from .colors import get_color, get_color_with_info
from .coordinates import get_exon_and_intron_coordinates
from .shapes import (
    DEFAULT_DIRECTIONAL_FEATURE_TYPES,
    FeatureRendering,
    default_feature_rendering,
    normalize_feature_shape_overrides,
)
from .tracks import arrange_feature_tracks


def create_repeat_object(
    repeat_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    genome_length: int,
    label_filtering,
    is_directional: bool,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
) -> RepeatObject:
    """
    Creates a RepeatObject representing a repeat region in a genome.
    """
    coordinates = feature.location.parts
    rpt_family: str = feature.qualifiers.get("rpt_family", ["undefined"])[0]
    rpt_type: str = feature.qualifiers.get("rpt_type", ["undefined"])[0]
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(coordinates, genome_length)
    color: str = get_color(feature, color_table, default_colors, record_id=record_id)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering, record_id=record_id) if compute_label_text else ""

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
    is_directional: bool,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
) -> FeatureObject:
    """
    Creates a FeatureObject representing a generic genomic feature.
    """
    coordinates = feature.location.parts
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(coordinates, genome_length)
    color: str = get_color(feature, color_table, default_colors, record_id=record_id)
    feature_type = feature.type
    label_text = get_label_text(feature, label_filtering, record_id=record_id) if compute_label_text else ""

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
    is_directional: bool,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
) -> GeneObject:
    """
    Creates a GeneObject representing a gene in a genome.
    """
    coordinates: List[SimpleLocation] = feature.location.parts
    note: str = feature.qualifiers.get("note", [""])[0]
    product: str = feature.qualifiers.get("product", [""])[0]
    gene: str = feature.qualifiers.get("gene", [""])[0]
    is_trans_spliced = "trans_splicing" in feature.qualifiers
    feature_type = feature.type
    location = get_exon_and_intron_coordinates(coordinates, genome_length, is_trans_spliced)
    color: str = get_color(feature, color_table, default_colors, record_id=record_id)
    label_text = get_label_text(feature, label_filtering, record_id=record_id) if compute_label_text else ""

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


@dataclass(frozen=True)
class FeatureBuildResult:
    """Visible features partitioned by their resolved rendering layer."""

    foreground_features: dict[str, FeatureObject]
    underlay_features: tuple[FeatureObject, ...]
    used_color_rules: frozenset[tuple[str, str]]


def _build_feature_layers(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
    rendering_resolver: Callable[[str], FeatureRendering],
    split_overlaps_by_strand: bool = False,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    compute_label_text: bool = True,
) -> FeatureBuildResult:
    foreground_features: Dict[str, FeatureObject] = {}
    underlay_features: list[FeatureObject] = []
    used_color_rules: Set[Tuple[str, str]] = set()
    locus_count: int = 0
    repeat_count: int = 0
    feature_count: int = 0
    genome_length: int = len(gb_record.seq)

    for feature in gb_record.features:
        if not should_render_feature(
            feature,
            selected_features_set,
            feature_visibility_rules=feature_visibility_rules,
            record_id=gb_record.id,
            specific_color_rules=color_table,
        ):
            continue

        # Track which color rules are actually used
        color, caption = get_color_with_info(feature, color_table, default_colors, record_id=gb_record.id)
        if caption:
            used_color_rules.add((caption, color))
        rendering = rendering_resolver(str(feature.type))
        is_directional = rendering == "arrow"
        include_label = compute_label_text and rendering != "underlay"

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
                is_directional,
                record_id=gb_record.id,
                compute_label_text=include_label,
            )
            feature_id, feature_object = locus_id, gene_object
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
                is_directional,
                record_id=gb_record.id,
                compute_label_text=include_label,
            )
            feature_id, feature_object = repeat_id, repeat_object
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
                is_directional,
                record_id=gb_record.id,
                compute_label_text=include_label,
            )
        if rendering == "underlay":
            underlay_features.append(feature_object)
        else:
            foreground_features[feature_id] = feature_object

    foreground_features = arrange_feature_tracks(
        foreground_features,
        separate_strands,
        resolve_overlaps,
        split_overlaps_by_strand=split_overlaps_by_strand,
        genome_length=genome_length,
    )
    return FeatureBuildResult(
        foreground_features=foreground_features,
        underlay_features=tuple(underlay_features),
        used_color_rules=frozenset(used_color_rules),
    )


def create_feature_layers(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
    split_overlaps_by_strand: bool = False,
    feature_shapes: Mapping[str, str] | None = None,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    compute_label_text: bool = True,
) -> FeatureBuildResult:
    """Build visible features using the current three-value rendering contract."""

    normalized_shapes = normalize_feature_shape_overrides(feature_shapes)
    return _build_feature_layers(
        gb_record,
        color_table,
        selected_features_set,
        default_colors,
        separate_strands,
        resolve_overlaps,
        label_filtering,
        rendering_resolver=lambda feature_type: normalized_shapes.get(
            feature_type,
            default_feature_rendering(feature_type),
        ),
        split_overlaps_by_strand=split_overlaps_by_strand,
        feature_visibility_rules=feature_visibility_rules,
        compute_label_text=compute_label_text,
    )


def create_feature_dict(
    gb_record: SeqRecord,
    color_table,
    selected_features_set: List[str],
    default_colors,
    separate_strands: bool,
    resolve_overlaps: bool,
    label_filtering,
    split_overlaps_by_strand: bool = False,
    directional_feature_types: Optional[Set[str]] = None,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    compute_label_text: bool = True,
) -> Tuple[Dict[str, FeatureObject], Set[Tuple[str, str]]]:
    """Build the legacy foreground-only feature dictionary.

    This compatibility entry point deliberately has no underlay input. All visible
    features, including ``repeat_region``, remain foreground objects. Main diagram
    assembly uses :func:`create_feature_layers` instead.
    """

    directional_types = (
        {str(feature_type) for feature_type in directional_feature_types}
        if directional_feature_types is not None
        else set(DEFAULT_DIRECTIONAL_FEATURE_TYPES)
    )
    result = _build_feature_layers(
        gb_record,
        color_table,
        selected_features_set,
        default_colors,
        separate_strands,
        resolve_overlaps,
        label_filtering,
        rendering_resolver=lambda feature_type: (
            "arrow" if feature_type in directional_types else "rectangle"
        ),
        split_overlaps_by_strand=split_overlaps_by_strand,
        feature_visibility_rules=feature_visibility_rules,
        compute_label_text=compute_label_text,
    )
    return result.foreground_features, set(result.used_color_rules)


__all__ = [
    "FeatureBuildResult",
    "create_feature_dict",
    "create_feature_layers",
    "create_feature_object",
    "create_gene_object",
    "create_repeat_object",
]
