#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, NamedTuple, Optional, Set, Tuple

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from ..core.record_metadata import (
    location_over_origin,
    _feature_source_index_map,
)
from .ids import compute_source_feature_hash, feature_instances, shared_source_feature_hashes
from .objects import GeneObject, RepeatObject, FeatureObject
from .overrides import feature_override_lookup
from .visibility import should_render_feature
from ..labels.filtering import get_label_text
from .colors import get_color, get_color_with_info
from .coordinates import get_exon_and_intron_coordinates, origin_join_order, project_feature_parts
from ..layout.record_coordinates import RecordDisplayTransform
from ..layout.record_transform import spans_may_cross_origin, unset_record_transform
from .shapes import (
    DEFAULT_DIRECTIONAL_FEATURE_TYPES,
    FeatureGlyph,
    FeatureRendering,
    default_feature_rendering,
    normalize_feature_shape_overrides,
)
from gbdraw.exceptions import ValidationError
from .tracks import arrange_feature_tracks
if TYPE_CHECKING:
    from .overrides import ResolvedFeatureOverride
    from .placement import FeaturePlacementSlot, ResolvedRecordFeatureInputs


class DrawnParts(NamedTuple):
    """The parts a feature object draws and the record length over whose origin they may run (None when its spans
    may not cross the origin)."""

    parts: Sequence[SimpleLocation]
    origin_length: int | None
    # The same parts in the order they read over the origin when the view joins them there into one block
    # (``origin_join_order``), else None; the feature's location is built in this order, its coordinates keep ``parts``.
    over_origin: Sequence[SimpleLocation] | None = None

    @property
    def in_drawing_order(self) -> list[SimpleLocation]:
        return list(self.over_origin if self.over_origin is not None else self.parts)


@dataclass(frozen=True)
class FeatureLayerMode:
    """What a diagram mode changes in its feature layers."""

    # Circular: overlapping features of the two strands take tracks apart.
    split_overlaps_by_strand: bool = False
    # Linear: a part past the end of a record whose spans may cross its origin is drawn as the join over the origin
    # (``location_over_origin``, OV-3724); Circular draws it by its angles.
    past_end_as_join: bool = False
    # Circular: two blocks that meet over the origin of a record drawn as a circle are drawn as one block there
    # (``origin_join_order``).
    joins_over_origin: bool = False


def create_repeat_object(
    repeat_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    label_filtering,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
    glyph_kind: FeatureGlyph | None = None,
    feature_override: ResolvedFeatureOverride | None = None,
    *,
    instance: Optional[str],
    drawn: DrawnParts,
) -> RepeatObject:
    """
    Creates a RepeatObject representing a repeat region in a genome, drawn as ``drawn``.
    """
    coordinates = list(drawn.parts)
    rpt_family: str = feature.qualifiers.get("rpt_family", ["undefined"])[0]
    rpt_type: str = feature.qualifiers.get("rpt_type", ["undefined"])[0]
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(drawn.in_drawing_order, origin_length=drawn.origin_length)
    color: str = get_color(feature, color_table, default_colors, record_id=record_id, instance=instance)
    feature_type = feature.type
    label_text = (
        get_label_text(feature, label_filtering, record_id=record_id, feature_override=feature_override)
        if compute_label_text else ""
    )

    repeat_object = RepeatObject(
        repeat_id,
        location,
        None,
        color,
        note,
        rpt_family,
        rpt_type,
        label_text,
        coordinates,
        feature_type,
        qualifiers=feature.qualifiers,
        record_id=record_id,
        glyph_kind=glyph_kind,
    )
    repeat_object.feature_override = feature_override
    return repeat_object


def create_feature_object(
    feature_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    label_filtering,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
    glyph_kind: FeatureGlyph | None = None,
    feature_override: ResolvedFeatureOverride | None = None,
    *,
    instance: Optional[str],
    drawn: DrawnParts,
) -> FeatureObject:
    """
    Creates a FeatureObject representing a generic genomic feature, drawn as ``drawn``.
    """
    coordinates = list(drawn.parts)
    note: str = feature.qualifiers.get("note", [""])[0]
    location = get_exon_and_intron_coordinates(drawn.in_drawing_order, origin_length=drawn.origin_length)
    color: str = get_color(feature, color_table, default_colors, record_id=record_id, instance=instance)
    feature_type = feature.type
    label_text = (
        get_label_text(feature, label_filtering, record_id=record_id, feature_override=feature_override)
        if compute_label_text else ""
    )

    feature_object = FeatureObject(
        feature_id,
        location,
        None,
        color,
        note,
        label_text,
        coordinates,
        feature_type,
        qualifiers=feature.qualifiers,
        record_id=record_id,
        glyph_kind=glyph_kind,
    )
    feature_object.feature_override = feature_override
    return feature_object


def create_gene_object(
    feature_id: str,
    feature: SeqFeature,
    color_table,
    default_colors,
    label_filtering,
    record_id: Optional[str] = None,
    compute_label_text: bool = True,
    glyph_kind: FeatureGlyph | None = None,
    feature_override: ResolvedFeatureOverride | None = None,
    *,
    instance: Optional[str],
    drawn: DrawnParts,
) -> GeneObject:
    """
    Creates a GeneObject representing a gene in a genome, drawn as ``drawn``.
    """
    coordinates: List[SimpleLocation] = list(drawn.parts)
    note: str = feature.qualifiers.get("note", [""])[0]
    product: str = feature.qualifiers.get("product", [""])[0]
    gene: str = feature.qualifiers.get("gene", [""])[0]
    is_trans_spliced = "trans_splicing" in feature.qualifiers
    feature_type = feature.type
    location = get_exon_and_intron_coordinates(
        drawn.in_drawing_order, is_trans_spliced, origin_length=drawn.origin_length
    )
    color: str = get_color(feature, color_table, default_colors, record_id=record_id, instance=instance)
    label_text = (
        get_label_text(feature, label_filtering, record_id=record_id, feature_override=feature_override)
        if compute_label_text else ""
    )

    gene_object = GeneObject(
        feature_id,
        location,
        None,
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
        glyph_kind=glyph_kind,
    )
    gene_object.feature_override = feature_override
    return gene_object


class _FeatureGeometry(NamedTuple):
    """What one feature draws in the view a layer is built for (``_feature_geometry``)."""

    drawn: DrawnParts
    # The parts the display projection draws, in drawing order; None when the record is not projected.
    projected: list[SimpleLocation] | None
    # Its two blocks are drawn as one block over the origin in this view (``FeatureObject.joins_over_origin``).
    joins_over_origin: bool


def _feature_geometry(
    feature: SeqFeature, genome_length: int, origin_length: int | None, layer_mode: FeatureLayerMode, *, projects: bool
) -> _FeatureGeometry:
    """The parts ``feature`` draws on its record, decided once per feature and view.

    A part past the end of a record whose spans may cross its origin is the join over the origin it is (OV-3724): a
    Linear axis and a display projection draw that join; the feature keeps its location, and Circular draws the
    part by its angles. A view that joins blocks over the origin draws the parts the contact rule joins
    (``origin_join_order``) in the order that reads over it; the drawers and the label anchors read the decision.
    """

    def over_origin(parts: Sequence[SimpleLocation]) -> list[SimpleLocation] | None:
        return origin_join_order(parts, origin_length) if layer_mode.joins_over_origin else None

    joined_parts = location_over_origin(feature.location, genome_length).parts if origin_length is not None else None
    parts = joined_parts if joined_parts is not None and layer_mode.past_end_as_join else feature.location.parts
    drawn = DrawnParts(parts, origin_length, over_origin(parts))
    if not projects:
        return _FeatureGeometry(drawn, None, drawn.over_origin is not None)
    projected = joined_parts if joined_parts is not None else feature.location.parts
    projected_over_origin = over_origin(projected)
    return _FeatureGeometry(
        drawn,
        projected_over_origin if projected_over_origin is not None else list(projected),
        projected_over_origin is not None,
    )


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
    layer_mode: FeatureLayerMode = FeatureLayerMode(),
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    compute_label_text: bool = True,
    record_transform: RecordDisplayTransform | None = None,
    record_features: ResolvedRecordFeatureInputs | None = None,
    placement_slot: FeaturePlacementSlot | None = None,
    feature_overlap_tolerance_bp: int = 0,
) -> FeatureBuildResult:
    foreground_features: Dict[str, FeatureObject] = {}
    underlay_features: list[FeatureObject] = []
    used_color_rules: Set[Tuple[str, str]] = set()
    locus_count: int = 0
    repeat_count: int = 0
    feature_count: int = 0
    genome_length: int = len(gb_record.seq)
    # Connector lines over the origin are drawn only where a span may cross it, in every view.
    origin_length = genome_length if spans_may_cross_origin(
        gb_record, record_transform if record_transform is not None else unset_record_transform(gb_record)
    ) else None
    source_indexes = _feature_source_index_map(gb_record.features)
    shared_hashes = shared_source_feature_hashes(gb_record)
    instance_of = feature_instances(gb_record)
    override_of = feature_override_lookup(
        gb_record, record_features.overrides if record_features is not None else None
    )
    # A canceled Worker request stops between features: this build's state is local (D7).
    from gbdraw.api.prepared import cancel_checkpoint

    for feature in gb_record.features:
        cancel_checkpoint()
        override = override_of(feature)
        instance = instance_of(feature)
        if not should_render_feature(
            feature,
            selected_features_set,
            feature_visibility_rules=feature_visibility_rules,
            record_id=gb_record.id,
            specific_color_rules=color_table,
            feature_override=override,
            instance=instance,
        ):
            continue

        # Track which color rules are actually used
        color, caption = get_color_with_info(
            feature, color_table, default_colors, record_id=gb_record.id, instance=instance
        )
        if caption:
            used_color_rules.add((caption, color))
        rendering = rendering_resolver(str(feature.type))
        glyph_kind: FeatureGlyph = (
            "rectangle" if rendering == "underlay" else rendering
        )
        include_label = compute_label_text and rendering != "underlay"
        geometry = _feature_geometry(
            feature, genome_length, origin_length, layer_mode,
            projects=record_transform is not None and record_transform.projects,
        )

        feature_id: str
        feature_object: FeatureObject
        if feature.type in {"CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"}:
            locus_count += 1
            locus_id: str = "gene_" + str(locus_count).zfill(9)
            gene_object: GeneObject = create_gene_object(
                locus_id,
                feature,
                color_table,
                default_colors,
                label_filtering,
                record_id=gb_record.id,
                compute_label_text=include_label,
                glyph_kind=glyph_kind,
                feature_override=override,
                instance=instance,
                drawn=geometry.drawn,
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
                label_filtering,
                record_id=gb_record.id,
                compute_label_text=include_label,
                glyph_kind=glyph_kind,
                feature_override=override,
                instance=instance,
                drawn=geometry.drawn,
            )
            feature_id, feature_object = repeat_id, repeat_object
        else:
            feature_count += 1
            feature_id = "feature_" + str(feature_count).zfill(9)
            feature_object = create_feature_object(
                feature_id,
                feature,
                color_table,
                default_colors,
                label_filtering,
                record_id=gb_record.id,
                compute_label_text=include_label,
                glyph_kind=glyph_kind,
                feature_override=override,
                instance=instance,
                drawn=geometry.drawn,
            )
        feature_object.source_feature_index = source_indexes[id(feature)]
        feature_object.feature_hash = compute_source_feature_hash(feature, gb_record)
        if feature_object.feature_hash in shared_hashes:
            feature_object.rendered_instance_index = feature_object.source_feature_index
        if geometry.projected is not None and record_transform is not None:
            feature_object.display_parts = project_feature_parts(
                geometry.projected, record_transform,
                origin_length=origin_length,
                is_trans_spliced="trans_splicing" in feature.qualifiers,
            )
        feature_object.joins_over_origin = geometry.joins_over_origin
        if rendering == "underlay":
            underlay_features.append(feature_object)
        else:
            foreground_features[feature_id] = feature_object

    if placement_slot is not None:
        from .placement import plan_feature_placements

        plan_feature_placements(
            foreground_features, slot=placement_slot, record_features=record_features,
            resolve_overlaps=resolve_overlaps, tolerance_bp=feature_overlap_tolerance_bp,
            genome_length=genome_length,
        )
    else:
        if record_features is not None and record_features.placements:
            raise ValidationError("Feature placement requires resolved slot geometry.")
        arrange_feature_tracks(
            foreground_features, separate_strands, resolve_overlaps,
            split_overlaps_by_strand=layer_mode.split_overlaps_by_strand,
            genome_length=genome_length, tolerance_bp=feature_overlap_tolerance_bp,
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
    layer_mode: FeatureLayerMode = FeatureLayerMode(),
    feature_shapes: Mapping[str, str] | None = None,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    compute_label_text: bool = True,
    record_transform: RecordDisplayTransform | None = None,
    record_features: ResolvedRecordFeatureInputs | None = None,
    placement_slot: FeaturePlacementSlot | None = None,
    feature_overlap_tolerance_bp: int = 0,
) -> FeatureBuildResult:
    """Build visible features using the current rendering contract; ``layer_mode`` holds what the diagram mode
    changes (``FeatureLayerMode``)."""

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
        layer_mode=layer_mode,
        feature_visibility_rules=feature_visibility_rules,
        compute_label_text=compute_label_text,
        record_transform=record_transform,
        record_features=record_features,
        placement_slot=placement_slot,
        feature_overlap_tolerance_bp=feature_overlap_tolerance_bp,
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
        layer_mode=FeatureLayerMode(split_overlaps_by_strand=split_overlaps_by_strand),
        feature_visibility_rules=feature_visibility_rules,
        compute_label_text=compute_label_text,
    )
    return result.foreground_features, set(result.used_color_rules)


__all__ = [
    "FeatureBuildResult",
    "create_feature_dict",
    "create_feature_layers",
    "DrawnParts",
    "FeatureLayerMode",
    "create_feature_object",
    "create_gene_object",
    "create_repeat_object",
]
