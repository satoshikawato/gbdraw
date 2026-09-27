"""Resolve typed record inputs once while retaining their source provenance."""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Hashable, Literal, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.analysis.protein_colinearity import (
    OrthogroupGraphResult,
    OrthogroupMember,
    OrthogroupResult,
)
from gbdraw.exceptions import ValidationError
from gbdraw.core.record_metadata import _iter_source_features, _read_coord_map, _source_feature_index
from gbdraw.io.cli_tables import (
    read_comparisons_table,
    read_conservation_table,
    read_records_table,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.genome import load_gbks, load_gff_fasta
from gbdraw.io.record_select import (
    RecordSelector,
    parse_record_selector,
    reverse_records,
)
from gbdraw.io.regions import RegionSpec, apply_region_specs, parse_region_spec
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_placement import resolve_record_row_positions
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentEvidenceEdge,
    AmbiguousAlignmentRecord,
    SimilarityAlignmentCandidate,
    SimilarityAlignmentPlan,
    resolve_similarity_alignment,
)
from gbdraw.linear_comparison import LinearComparison

from .options import (
    CircularDiagramOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from gbdraw.features.source import SourceFeatureIdentity, build_source_feature_catalog
from gbdraw.features.ids import compute_feature_hash_from_location_parts

from .prepared import (
    ParsedRecordInputs,
    PreparedResourceIdentity,
    get_or_build_parsed_source,
    prepared_resource_identity,
)
from .requests import (
    CircularBatchOutputPolicy,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    RecordCardinality,
    RecordCollectionOptions,
    RecordDisplayOptions,
    RecordInput,
    RecordInputSource,
    RecordPresentation,
    RenderOutputRequest,
)

logger = logging.getLogger(__name__)

GenBankLoader = Callable[..., list[SeqRecord]]
GffFastaLoader = Callable[..., list[SeqRecord]]


@dataclass(frozen=True)
class ResolvedRecordDisplay:
    """Transient display intent resolved from complete source facts."""

    detected_topology: Literal["circular", "linear", "unknown"]
    is_circular: bool
    start_coordinate: int | None
    current_start_coordinate: int
    orientation_step: Literal[-1, 1]


def resolve_record_display(
    options: RecordDisplayOptions,
    *,
    source_length: int,
    detected_topology: Literal["circular", "linear", "unknown"],
    source_base: int,
    source_step: Literal[-1, 1],
    has_input_region: bool = False,
    has_collection_region: bool = False,
    is_cropped: bool = False,
) -> ResolvedRecordDisplay:
    """Validate supplied source facts without loading or transforming records.

    The caller supplies the full source length/topology, the existing affine
    map from core.record_metadata, and crop applicability for this record.
    Neither source facts nor crop state are inferred from a materialized length.
    current_start_coordinate is the existing affine start before extra rotation.
    """
    if not isinstance(options, RecordDisplayOptions):
        raise ValidationError("Record display has an unsupported type.")
    if detected_topology not in ("circular", "linear", "unknown"):
        raise ValidationError("Detected topology must be circular, linear, or unknown.")
    crop_flags = (has_input_region, has_collection_region, is_cropped)
    if not all(isinstance(value, bool) for value in crop_flags):
        raise ValidationError("Resolved crop flags must be booleans.")
    if options.start_coordinate is not None and any(crop_flags):
        raise ValidationError("An explicit display start cannot be combined with a crop.")
    is_circular = (
        detected_topology == "circular"
        if options.is_circular is None
        else options.is_circular
    )
    # The transform owns length, affine-domain, orientation and anchor validation.
    RecordDisplayTransform(
        source_length, source_base, source_step, options.start_coordinate, is_circular,
    )
    return ResolvedRecordDisplay(
        detected_topology=detected_topology,
        is_circular=is_circular,
        start_coordinate=options.start_coordinate,
        current_start_coordinate=source_base,
        orientation_step=source_step,
    )


@dataclass(frozen=True)
class ResolvedRecordProvenance:
    """Stable source coordinates for one displayed record."""

    resolved_index: int
    input_index: int
    source_record_index: int
    source_record_count: int
    source_record_id: str
    source_kind: Literal["genbank", "gff_fasta", "memory"]
    source_paths: tuple[str, ...]
    record_key: str
    cardinality: RecordCardinality
    selector: RecordSelector | None
    region: RegionSpec | None
    presentation: RecordPresentation
    display: RecordDisplayOptions = field(default_factory=RecordDisplayOptions)
    source_length: int | None = None
    detected_topology: Literal["circular", "linear", "unknown"] = "unknown"
    is_cropped: bool = False
    has_collection_region: bool = False
    resolved_display: ResolvedRecordDisplay | None = None
    source_feature_catalog: tuple[SourceFeatureIdentity, ...] | None = field(default=None, repr=False)


@dataclass(frozen=True)
class ResolvedRecordCollection:
    """Flattened records and one provenance entry per displayed record."""

    records: tuple[SeqRecord, ...]
    provenance: tuple[ResolvedRecordProvenance, ...]
    displays: tuple[ResolvedRecordDisplay, ...] = field(init=False)
    transforms: tuple[RecordDisplayTransform, ...] = field(init=False)

    def __post_init__(self) -> None:
        if not self.records:
            raise ValidationError("A resolved record collection cannot be empty.")
        if len(self.records) != len(self.provenance):
            raise ValidationError(
                "Resolved record provenance must align with displayed records."
            )
        provenance = []
        displays = []
        transforms = []
        for record, item in zip(self.records, self.provenance, strict=True):
            base, step = _read_coord_map(record)
            cropped = item.is_cropped or bool(record.annotations.get("gbdraw_region_applied"))
            length = item.source_length
            if length is None and not cropped:
                length = len(record)
            try:
                display = resolve_record_display(
                    item.display,
                    source_length=length,
                    detected_topology=item.detected_topology,
                    source_base=base,
                    source_step=step,
                    has_input_region=item.region is not None,
                    has_collection_region=item.has_collection_region,
                    is_cropped=cropped,
                )
                transform = RecordDisplayTransform(
                    length, base, step, display.start_coordinate, display.is_circular,
                )
            except ValidationError as exc:
                raise ValidationError(
                    f"Record {item.record_key!r} (input {item.input_index + 1}, "
                    f"source {item.source_paths or item.source_kind}, "
                    f"selector {item.selector}, record {item.source_record_id!r}): {exc}"
                ) from exc
            displays.append(display)
            transforms.append(transform)
            provenance.append(replace(item, resolved_display=display))
        object.__setattr__(self, "provenance", tuple(provenance))
        object.__setattr__(self, "displays", tuple(displays))
        object.__setattr__(self, "transforms", tuple(transforms))


def project_source_bound_comparisons(
    options: LinearDiagramOptions, collection: ResolvedRecordCollection,
) -> LinearDiagramOptions:
    """Reproject existing source-bound evidence into the requested record views.

    View hashes must match this source/crop's current or opposite orientation.
    Source identities, scores and block membership remain unchanged. Coordinates
    are record-local; the existing renderer alone applies a circular display cut.
    """
    frames = [comparison.matches for comparison in options.linear_comparisons or ()]
    frames.extend(options.protein_comparisons or ())
    required = {f"{role}_{name}" for role in ("query", "subject")
                for name in ("feature_index", "feature_svg_id")}
    if not any(not frame.empty and required <= set(frame.columns) for frame in frames):
        return options
    bindings = []
    for record, provenance in zip(collection.records, collection.provenance, strict=True):
        current, opposite = {}, {}
        for ordinal, feature in enumerate(_iter_source_features(record.features)):
            if feature.location is None:
                continue
            index = _source_feature_index(feature)
            index = ordinal if index is None else index
            for target, location in ((current, feature.location),
                                     (opposite, feature.location._flip(len(record)))):
                target[index] = compute_feature_hash_from_location_parts(
                    feature.type,
                    [(int(part.start), int(part.end), part.strand) for part in location.parts],
                    record_id=record.id,
                )
        source = {feature.source_feature_index: feature.stable_feature_id
                  for feature in provenance.source_feature_catalog or ()}
        bindings.append((source, current, opposite, len(record)))

    def project(frame: pd.DataFrame, query_index: int, subject_index: int) -> pd.DataFrame:
        required = {f"{role}_{name}" for role in ("query", "subject")
                    for name in ("feature_index", "feature_svg_id")}
        if frame.empty or not required <= set(frame.columns):
            return frame
        updated = None
        for row_index, row in frame.iterrows():
            flips = []
            for role, record_index, prefix in (("query", query_index, "q"), ("subject", subject_index, "s")):
                if not 0 <= record_index < len(bindings):
                    raise ValidationError("Comparison record endpoint is outside the resolved collection.")
                source, current, opposite, length = bindings[record_index]
                try:
                    indices = [int(value) for value in str(row[f"{role}_feature_index"]).split(";")]
                except ValueError:
                    raise ValidationError("Comparison source feature indexes must be integers.") from None
                source_ids = str(row[f"{role}_feature_svg_id"]).split(";")
                view_ids = str(row.get(f"{role}_view_feature_svg_id", row[f"{role}_feature_svg_id"])).split(";")
                if not len(indices) == len(source_ids) == len(view_ids):
                    raise ValidationError("Comparison source feature binding has inconsistent coverage.")
                if any(source.get(index) != source_id for index, source_id in zip(indices, source_ids, strict=True)):
                    raise ValidationError("Comparison source feature index conflicts with its source feature ID.")
                current_ids, opposite_ids = ([views.get(index) for index in indices] for views in (current, opposite))
                flip = view_ids != current_ids
                if flip and view_ids != opposite_ids:
                    raise ValidationError("Comparison view feature IDs do not match the current source/crop binding.")
                flips.append(flip)
                if not flip:
                    continue
                if updated is None:
                    updated = frame.copy(deep=True)
                for coordinate_field in (f"{prefix}start", f"{prefix}end"):
                    coordinate = row[coordinate_field]
                    if not 1 <= coordinate <= length or int(coordinate) != coordinate:
                        raise ValidationError("Comparison view coordinate must be a finite genomic base within its record.")
                    updated.at[row_index, coordinate_field] = length + 1 - int(coordinate)
                updated.at[row_index, f"{role}_view_feature_svg_id"] = ";".join(current_ids)
            if updated is not None and flips[0] != flips[1] and "collinearity_orientation" in frame.columns:
                orientation = row["collinearity_orientation"]
                if orientation in ("plus", "minus"):
                    updated.at[row_index, "collinearity_orientation"] = "minus" if orientation == "plus" else "plus"
        return frame if updated is None else updated

    explicit = tuple(replace(comparison, matches=project(comparison.matches,
        comparison.query_record_index, comparison.subject_record_index))
        for comparison in options.linear_comparisons or ())
    protein = None if options.protein_comparisons is None else tuple(
        project(frame, index, index + 1) for index, frame in enumerate(options.protein_comparisons))
    changed = any(before.matches is not after.matches for before, after in
                  zip(options.linear_comparisons or (), explicit, strict=True))
    changed = changed or (protein is not None and any(before is not after for before, after in
                         zip(options.protein_comparisons or (), protein, strict=True)))
    return replace(options, linear_comparisons=explicit, protein_comparisons=protein) if changed else options


def project_similarity_alignment_centers(
    collection: ResolvedRecordCollection,
    plan: SimilarityAlignmentPlan | None,
) -> tuple[float | None, ...]:
    """Project selected anchor centers from already resolved record orientations."""

    if plan is None:
        return tuple(None for _ in collection.records)
    record_keys = tuple(item.record_key for item in collection.provenance)
    plan.validate_record_coverage(record_keys)
    decisions = {decision.record_key: decision for decision in plan.records}
    return tuple(
        _project_alignment_anchor_center(
            decision.anchor,
            provenance=provenance,
            transform=transform,
            displayed_length=len(record),
        )
        if decision.status is not AlignmentDecisionStatus.SKIPPED
        else None
        for record, provenance, transform, decision in zip(
            collection.records,
            collection.provenance,
            collection.transforms,
            (decisions[key] for key in record_keys),
            strict=True,
        )
    )


def _project_alignment_anchor_center(
    anchor: AlignmentAnchorIdentity | None,
    *,
    provenance: ResolvedRecordProvenance,
    transform: RecordDisplayTransform,
    displayed_length: int,
) -> float:
    feature = _alignment_source_feature(anchor, provenance)
    center = (
        min(part[0] for part in feature.location_parts)
        + max(part[1] for part in feature.location_parts)
    ) / 2.0
    projected = transform.source_position_to_display_offset(center)
    if not 0.0 <= projected <= float(displayed_length):
        raise ValidationError(
            "Similarity alignment anchor is outside the current crop for record "
            f"{provenance.record_key!r}."
        )
    return projected



def _alignment_source_feature(
    anchor: AlignmentAnchorIdentity | None,
    provenance: ResolvedRecordProvenance,
) -> SourceFeatureIdentity:
    if anchor is None:
        raise ValidationError(
            f"Similarity alignment record {provenance.record_key!r} has no anchor."
        )
    catalog = provenance.source_feature_catalog
    if catalog is None:
        raise ValidationError(
            f"Similarity alignment record {provenance.record_key!r} has no source feature catalog."
        )
    if anchor.source_feature_index is not None:
        matches = [
            feature
            for feature in catalog
            if feature.source_feature_index == anchor.source_feature_index
        ]
    elif anchor.stable_feature_svg_id is not None:
        matches = [
            feature
            for feature in catalog
            if feature.stable_feature_id == anchor.stable_feature_svg_id
        ]
    else:
        matches = [
            feature
            for feature in catalog
            if anchor.biological_feature_id
            in {feature.biological_feature_id, feature.stable_feature_id}
        ]
    if (
        len(matches) == 1
        and anchor.stable_feature_svg_id is not None
        and matches[0].stable_feature_id != anchor.stable_feature_svg_id
    ):
        matches = []
    if len(matches) != 1:
        raise ValidationError(
            "Similarity alignment anchor must resolve to exactly one source feature "
            f"for record {provenance.record_key!r}."
        )
    return matches[0]


@dataclass(frozen=True)
class AlignmentAnchorDisplayFact:
    """Source-bound anchor geometry in one already resolved display transform."""

    source_start: int
    source_end: int
    source_strand: int | None
    displayed_strand: int | None
    display_center: float | None


def project_similarity_alignment_anchor_fact(
    collection: ResolvedRecordCollection,
    anchor: AlignmentAnchorIdentity,
) -> AlignmentAnchorDisplayFact:
    """Use the same source identity and center projection as the renderer."""

    matches = [i for i, item in enumerate(collection.provenance)
               if item.record_key == anchor.record_key]
    if len(matches) != 1:
        raise ValidationError("Alignment anchor record coverage is invalid.")
    index = matches[0]
    provenance = collection.provenance[index]
    feature = _alignment_source_feature(anchor, provenance)
    if anchor.biological_feature_id not in {feature.biological_feature_id, feature.stable_feature_id}:
        raise ValidationError("Alignment anchor biological identity changed.")
    strands = {part[2] for part in feature.location_parts}
    strand = next(iter(strands)) if len(strands) == 1 and strands <= {-1, 1} else None
    try:
        center = _project_alignment_anchor_center(
            anchor, provenance=provenance, transform=collection.transforms[index],
            displayed_length=len(collection.records[index]),
        )
    except ValidationError:
        center = None  # An exact source anchor outside the crop is unusable.
    return AlignmentAnchorDisplayFact(
        min(part[0] for part in feature.location_parts),
        max(part[1] for part in feature.location_parts),
        strand, None if strand is None else strand * collection.transforms[index].source_step, center,
    )

def resolve_cli_similarity_alignment_plan(
    collection: ResolvedRecordCollection,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    *,
    exact_reference: str,
) -> SimilarityAlignmentPlan:
    """Resolve the strict non-interactive CLI adapter through the shared resolver."""

    target = str(exact_reference).strip()
    if not target or "\0" in target:
        raise ValidationError(
            "--align_orthogroup_feature requires a non-empty exact feature/protein ID."
        )
    if orthogroups is None:
        raise ValidationError(
            "--align_orthogroup_feature requires orthogroup metadata from the requested analysis."
        )
    if target in orthogroups.orthogroups:
        raise ValidationError(
            "--align_orthogroup_feature accepts an exact feature/protein ID, not "
            f"Similarity Group ID {target!r}."
        )
    record_keys = tuple(item.record_key for item in collection.provenance)
    candidate_rows: list[tuple[OrthogroupMember, SimilarityAlignmentCandidate]] = []
    exact_matches: list[tuple[OrthogroupMember, SimilarityAlignmentCandidate]] = []
    for members in orthogroups.orthogroups.values():
        for member in members:
            if member.record_index < 0 or member.record_index >= len(collection.records):
                continue
            candidate = _similarity_candidate_from_orthogroup_member(
                collection,
                member,
            )
            candidate_rows.append((member, candidate))
            catalog_alias = candidate.anchor.biological_feature_id
            aliases = {
                str(value)
                for value in (
                    member.protein_id,
                    member.source_protein_id,
                    member.feature_svg_id,
                    catalog_alias,
                    candidate.anchor.stable_feature_svg_id,
                )
                if value
            }
            if target in aliases:
                exact_matches.append((member, candidate))
    if not exact_matches:
        raise ValidationError(
            "--align_orthogroup_feature did not match an exact feature/protein ID."
        )
    if len(exact_matches) != 1:
        matches = ", ".join(
            sorted(
                f"{candidate.anchor.record_key}:{member.protein_id}"
                for member, candidate in exact_matches
            )
        )
        raise ValidationError(
            "--align_orthogroup_feature is ambiguous; use one exact candidate ID "
            f"from: {matches}."
        )
    reference_member, reference = exact_matches[0]
    group_id = reference_member.orthogroup_id
    group_candidates = tuple(
        candidate
        for member, candidate in candidate_rows
        if member.orthogroup_id == group_id
    )
    identity_by_protein = {
        member.protein_id: candidate.anchor
        for member, candidate in candidate_rows
        if member.orthogroup_id == group_id
    }
    edges = tuple(
        AlignmentEvidenceEdge(
            group_id=group_id,
            query=identity_by_protein[edge.query_protein_id],
            subject=identity_by_protein[edge.subject_protein_id],
            edge_kind=str(edge.edge_kind),
        )
        for edge in orthogroups.ortholog_edges_by_orthogroup_id.get(group_id, ())
        if edge.query_protein_id in identity_by_protein
        and edge.subject_protein_id in identity_by_protein
    )
    resolution = resolve_similarity_alignment(
        record_keys=record_keys,
        group_id=group_id,
        reference=reference.anchor,
        candidates=group_candidates,
        edges=edges,
    )
    if resolution.ambiguities:
        details = "; ".join(
            _cli_ambiguity_message(ambiguity, candidate_rows)
            for ambiguity in resolution.ambiguities
        )
        raise ValidationError(
            "--align_orthogroup_feature cannot choose among multiple candidates; "
            f"select an exact candidate for each record: {details}."
        )
    return resolution.require_plan()


def _similarity_candidate_from_orthogroup_member(
    collection: ResolvedRecordCollection,
    member: OrthogroupMember,
) -> SimilarityAlignmentCandidate:
    provenance = collection.provenance[member.record_index]
    catalog = provenance.source_feature_catalog or ()
    matches = [
        feature
        for feature in catalog
        if member.feature_svg_id is not None
        and feature.stable_feature_id == member.feature_svg_id
    ]
    if not matches:
        matches = [
            feature
            for feature in catalog
            if feature.source_feature_index == member.feature_index
        ]
    identity_is_unique = len(matches) == 1
    feature = matches[0] if identity_is_unique else None
    anchor = AlignmentAnchorIdentity(
        record_key=provenance.record_key,
        biological_feature_id=(
            feature.biological_feature_id
            if feature is not None
            else str(member.feature_svg_id or member.protein_id)
        ),
        source_feature_index=(
            feature.source_feature_index if feature is not None else member.feature_index
        ),
        stable_feature_svg_id=(
            feature.stable_feature_id
            if feature is not None
            else member.feature_svg_id
        ),
    )
    display_center: float | None = None
    if feature is not None:
        source_center = (
            min(part[0] for part in feature.location_parts)
            + max(part[1] for part in feature.location_parts)
        ) / 2.0
        projected = collection.transforms[
            member.record_index
        ].source_position_to_display_offset(source_center)
        if 0.0 <= projected <= float(len(collection.records[member.record_index])):
            display_center = projected
    return SimilarityAlignmentCandidate(
        group_id=member.orthogroup_id,
        anchor=anchor,
        displayed_strand=member.strand if member.strand in (-1, 1) else None,
        center_mappable=display_center is not None,
        display_center=display_center,
        identity_is_unique=identity_is_unique,
        representative=member.representative,
        role=str(member.role),
        source_start=member.start,
        source_end=member.end,
        display_name=member.gene or member.product or member.label or member.protein_id,
    )


def _cli_ambiguity_message(
    ambiguity: AmbiguousAlignmentRecord,
    rows: Sequence[tuple[OrthogroupMember, SimilarityAlignmentCandidate]],
) -> str:
    protein_by_key = {
        candidate.anchor.canonical_key: member.protein_id
        for member, candidate in rows
    }
    candidate_ids = ", ".join(
        protein_by_key.get(candidate.anchor.canonical_key, candidate.anchor.biological_feature_id)
        for candidate in ambiguity.candidates
    )
    return f"record {ambiguity.record_key!r} candidates [{candidate_ids}]"


def _detected_topology(record: SeqRecord, source_kind: str):
    value = str(record.annotations.get("topology", "")).strip().lower()
    return value if source_kind != "gff_fasta" and value in {"circular", "linear"} else "unknown"


@dataclass(frozen=True)
class RecordInputManifest:
    """Unresolved typed record inputs derived from CLI path syntax."""

    records: tuple[RecordInput, ...]
    record_options: RecordCollectionOptions
    source_paths: tuple[str, ...]
    multi_record_positions: tuple[str, ...] = ()


def _expanded_cli_track_values(
    values: Sequence[object] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[object | None]:
    items = list(values or ())
    if not items:
        return [None] * track_count
    if len(items) == 1:
        return items * track_count
    if len(items) != track_count:
        raise ValidationError(
            f"{field_name} count must be one or equal to the number of "
            f"depth tracks ({track_count})."
        )
    return items


def _optional_positive_float(value: object | None, *, field_name: str) -> float | None:
    if value is None or str(value).strip().lower() in {
        "",
        "auto",
        "none",
        "null",
        "-",
    }:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{field_name} values must be numeric or auto.") from exc
    if numeric <= 0:
        raise ValidationError(f"{field_name} values must be > 0.")
    return numeric


def depth_track_inputs_from_cli(
    groups: Sequence[Sequence[str]] | None,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    heights: Sequence[object] | None = None,
    large_tick_intervals: Sequence[object] | None = None,
    small_tick_intervals: Sequence[object] | None = None,
    tick_font_sizes: Sequence[object] | None = None,
) -> tuple[DepthTrackInput, ...] | None:
    """Translate logical CLI depth groups without knowing record cardinality."""

    if not groups:
        return None
    track_count = len(groups)
    label_values = _expanded_cli_track_values(
        labels,
        track_count=track_count,
        field_name="depth_track_labels",
    )
    color_values = _expanded_cli_track_values(
        colors,
        track_count=track_count,
        field_name="depth_track_colors",
    )
    height_values = _expanded_cli_track_values(
        heights,
        track_count=track_count,
        field_name="depth_track_heights",
    )
    large_values = _expanded_cli_track_values(
        large_tick_intervals,
        track_count=track_count,
        field_name="depth_track_large_tick_intervals",
    )
    small_values = _expanded_cli_track_values(
        small_tick_intervals,
        track_count=track_count,
        field_name="depth_track_small_tick_intervals",
    )
    font_values = _expanded_cli_track_values(
        tick_font_sizes,
        track_count=track_count,
        field_name="depth_track_tick_font_sizes",
    )
    inputs: list[DepthTrackInput] = []
    for index, group in enumerate(groups):
        sources = tuple(
            None
            if str(path).strip().lower() in {"", "-", "none", "null"}
            else str(path)
            for path in group
        )
        if not sources or not any(source is not None for source in sources):
            raise ValidationError(
                f"--depth_track #{index + 1} must include at least one file."
            )
        source: object = sources[0] if len(sources) == 1 else sources
        inputs.append(
            DepthTrackInput(
                source=source,
                label=(
                    str(label_values[index]).strip()
                    if label_values[index] is not None
                    else None
                ),
                color=(
                    str(color_values[index]).strip()
                    if color_values[index] is not None
                    else None
                ),
                height=_optional_positive_float(
                    height_values[index],
                    field_name="depth_track_heights",
                ),
                large_tick_interval=_optional_positive_float(
                    large_values[index],
                    field_name="depth_track_large_tick_intervals",
                ),
                small_tick_interval=_optional_positive_float(
                    small_values[index],
                    field_name="depth_track_small_tick_intervals",
                ),
                tick_font_size=_optional_positive_float(
                    font_values[index],
                    field_name="depth_track_tick_font_sizes",
                ),
            )
        )
    return tuple(inputs)


def _source_key(source: RecordInputSource) -> tuple[object, ...]:
    if isinstance(source, GenBankInputSource):
        return ("genbank", str(source.path))
    if isinstance(source, GffFastaInputSource):
        return ("gff_fasta", str(source.gff_path), str(source.fasta_path))
    if isinstance(source, InMemoryRecordSource):
        return ("memory", id(source.record), id(source.source_feature_catalog))
    raise ValidationError("Unsupported record input source.")


def _source_details(
    source: RecordInputSource,
) -> tuple[Literal["genbank", "gff_fasta", "memory"], tuple[str, ...]]:
    if isinstance(source, GenBankInputSource):
        return "genbank", (str(source.path),)
    if isinstance(source, GffFastaInputSource):
        return "gff_fasta", (str(source.gff_path), str(source.fasta_path))
    if isinstance(source, InMemoryRecordSource):
        return "memory", ()
    raise ValidationError("Unsupported record input source.")


def _load_source_records(
    source: RecordInputSource,
    *,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    genbank_loader: GenBankLoader,
    gff_loader: GffFastaLoader,
) -> ParsedRecordInputs:
    def load() -> ParsedRecordInputs:
        gff_catalogs: list[tuple[SourceFeatureIdentity, ...]] = []
        if isinstance(source, GenBankInputSource):
            records = genbank_loader([str(source.path)])
        elif isinstance(source, GffFastaInputSource):
            records = gff_loader(
                [str(source.gff_path)],
                [str(source.fasta_path)],
                selected_features_set=gff_candidate_features,
                keep_all_features=gff_keep_all_features,
                source_feature_catalogs=gff_catalogs,
            )
        elif isinstance(source, InMemoryRecordSource):
            records = [source.record]
        else:  # pragma: no cover - RecordInput validates this union.
            raise ValidationError("Unsupported record input source.")
        if not records:
            raise ValidationError("A record input source resolved to no records.")
        catalogs = (
            (source.source_feature_catalog,)
            if isinstance(source, InMemoryRecordSource) and source.source_feature_catalog is not None
            else tuple(gff_catalogs) if isinstance(source, GffFastaInputSource)
            else tuple(build_source_feature_catalog(record) for record in records)
        )
        return ParsedRecordInputs(tuple(records), catalogs)

    cache_spec = _prepared_source_cache_spec(
        source,
        gff_candidate_features=gff_candidate_features,
        gff_keep_all_features=gff_keep_all_features,
    )
    if cache_spec is None:
        return load()
    key, identities = cache_spec
    return get_or_build_parsed_source(
        key,
        identities,
        load,
        publish=lambda value: bool(value),
    )


def _prepared_source_cache_spec(
    source: RecordInputSource,
    *,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
) -> tuple[Hashable, frozenset[PreparedResourceIdentity]] | None:
    """Build a path-independent parsed-source cache key for a Web resource."""

    if isinstance(source, GenBankInputSource):
        identity = prepared_resource_identity(source.path)
        if identity is None:
            return None
        return ("parsed-source-v1", "genbank", identity), frozenset({identity})
    if isinstance(source, GffFastaInputSource):
        gff_identity = prepared_resource_identity(source.gff_path)
        fasta_identity = prepared_resource_identity(source.fasta_path)
        if gff_identity is None or fasta_identity is None:
            return None
        identities = frozenset({gff_identity, fasta_identity})
        return (
            (
                "parsed-source-v1",
                "gff-fasta",
                gff_identity,
                fasta_identity,
                tuple(sorted(gff_candidate_features or ())),
                bool(gff_keep_all_features),
            ),
            identities,
        )
    return None


def _selector_from_region(region: RegionSpec | None) -> RecordSelector | None:
    if region is None:
        return None
    if region.record_id is not None:
        return RecordSelector(
            raw=region.record_id,
            record_id=region.record_id,
            record_index=None,
        )
    if region.record_index is not None:
        return RecordSelector(
            raw=f"#{region.record_index + 1}",
            record_id=None,
            record_index=region.record_index,
        )
    return None


def _selected_source_indexes(
    records: Sequence[SeqRecord],
    selector: RecordSelector | None,
) -> list[int]:
    if selector is None:
        return list(range(len(records)))
    if selector.record_index is not None:
        index = selector.record_index
        if index < 0 or index >= len(records):
            raise ValidationError(
                f"Record selector {selector.label()} is out of range "
                f"(loaded {len(records)} record(s))."
            )
        return [index]
    record_id = selector.record_id or ""
    matches = [
        index for index, record in enumerate(records) if record.id == record_id
    ]
    if not matches:
        raise ValidationError(
            f"Record selector '{record_id}' did not match any record ID."
        )
    if len(matches) > 1:
        raise ValidationError(
            f"Record selector '{record_id}' matched multiple records. "
            "Use #index to disambiguate."
        )
    return matches


def _cardinality_indexes(
    indexes: Sequence[int],
    *,
    cardinality: RecordCardinality,
    input_index: int,
) -> tuple[int, ...]:
    selected = tuple(indexes)
    if cardinality is RecordCardinality.EXACTLY_ONE:
        if len(selected) != 1:
            raise ValidationError(
                f"RecordInput #{input_index + 1} requires exactly one record; "
                f"resolved {len(selected)}. Add a selector, use "
                "RecordCardinality.FIRST, or explicitly use RecordCardinality.ALL. "
                "For per-record CLI settings, use --records_table."
            )
        return selected
    if not selected:
        raise ValidationError(
            f"RecordInput #{input_index + 1} resolved no records."
        )
    if cardinality is RecordCardinality.FIRST:
        return selected[:1]
    if cardinality is RecordCardinality.ALL:
        return selected
    raise ValidationError("Unsupported record cardinality.")


def _unqualified_region(region: RegionSpec) -> RegionSpec:
    return replace(
        region,
        file_selector=None,
        record_id=None,
        record_index=None,
    )


def _apply_presentation(
    record: SeqRecord,
    presentation: RecordPresentation,
    *,
    record_key: str,
) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    if presentation.label:
        record.annotations["gbdraw_record_label"] = presentation.label
    if presentation.subtitle:
        record.annotations["gbdraw_record_subtitle"] = presentation.subtitle
    record.annotations["gbdraw_record_key"] = record_key


def _apply_provenance_annotations(
    record: SeqRecord,
    provenance: ResolvedRecordProvenance,
) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations.update(
        {
            "gbdraw_input_index": provenance.input_index,
            "gbdraw_source_record_index": provenance.source_record_index,
            "gbdraw_source_record_count": provenance.source_record_count,
            "gbdraw_source_record_id": provenance.source_record_id,
            "gbdraw_record_key": provenance.record_key,
            "gbdraw_record_cardinality": provenance.cardinality.value,
            "gbdraw_source_kind": provenance.source_kind,
            "gbdraw_source_paths": provenance.source_paths,
        }
    )
    if provenance.source_length is not None:
        record.annotations["gbdraw_source_length"] = provenance.source_length
    if provenance.source_paths:
        record.annotations["gbdraw_source_file"] = provenance.source_paths[0]
        record.annotations["gbdraw_source_basename"] = Path(
            provenance.source_paths[0]
        ).name


def _apply_collection_options(
    records: list[SeqRecord],
    provenance: Sequence[ResolvedRecordProvenance],
    options: RecordCollectionOptions,
) -> list[SeqRecord]:
    if options.regions:
        try:
            records = apply_region_specs(records, options.regions, log=logger)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
    for field_name, annotation_key in (
        ("labels", "gbdraw_record_label"),
        ("subtitles", "gbdraw_record_subtitle"),
    ):
        values = getattr(options, field_name)
        if len(values) > len(records):
            logger.warning(
                "WARNING: More record %s were provided than records resolved; "
                "extra values will be ignored.",
                field_name,
            )
        for index, value in enumerate(values[: len(records)]):
            if value:
                records[index].annotations[annotation_key] = value
    for record, item in zip(records, provenance, strict=True):
        _apply_provenance_annotations(record, item)
    return records


def resolve_record_inputs(
    record_inputs: Sequence[RecordInput],
    *,
    record_options: RecordCollectionOptions | None = None,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    genbank_loader: GenBankLoader = load_gbks,
    gff_loader: GffFastaLoader = load_gff_fasta,
) -> ResolvedRecordCollection:
    """Load each unique source once, then apply typed selection and transforms."""

    inputs = tuple(record_inputs)
    if not inputs:
        raise ValidationError("A request requires at least one RecordInput.")
    cache: dict[tuple[object, ...], ParsedRecordInputs] = {}
    records: list[SeqRecord] = []
    provenance: list[ResolvedRecordProvenance] = []
    for input_index, record_input in enumerate(inputs):
        key = _source_key(record_input.source)
        parsed = cache.get(key)
        if parsed is None:
            parsed = _load_source_records(
                record_input.source,
                gff_candidate_features=gff_candidate_features,
                gff_keep_all_features=gff_keep_all_features,
                genbank_loader=genbank_loader,
                gff_loader=gff_loader,
            )
            cache[key] = parsed
        raw_records = parsed.records
        selector = record_input.selector or _selector_from_region(record_input.region)
        source_indexes = _cardinality_indexes(
            _selected_source_indexes(raw_records, selector),
            cardinality=record_input.cardinality,
            input_index=input_index,
        )
        source_kind, source_paths = _source_details(record_input.source)
        expands = len(source_indexes) > 1
        for source_record_index in source_indexes:
            source_record = raw_records[source_record_index]
            source_cropped = bool(source_record.annotations.get("gbdraw_region_applied"))
            record = copy.deepcopy(source_record)
            record = reverse_records(
                (record,),
                record_input.presentation.reverse_complement,
                log=logger,
            )[0]
            if record_input.region is not None:
                try:
                    record = apply_region_specs(
                        (record,),
                        (_unqualified_region(record_input.region),),
                        log=logger,
                    )[0]
                except ValueError as exc:
                    raise ValidationError(str(exc)) from exc
            base_key = record_input.record_key or f"record-{input_index + 1}"
            record_key = (
                f"{base_key}:{source_record_index + 1}"
                if expands
                else base_key
            )
            item = ResolvedRecordProvenance(
                resolved_index=len(records),
                input_index=input_index,
                source_record_index=source_record_index,
                source_record_count=len(raw_records),
                source_record_id=str(raw_records[source_record_index].id),
                source_kind=source_kind,
                source_paths=source_paths,
                record_key=record_key,
                cardinality=record_input.cardinality,
                selector=selector,
                region=record_input.region,
                presentation=record_input.presentation,
                display=record_input.display,
                source_length=(source_record.annotations.get("gbdraw_source_length")
                               if source_cropped else len(source_record)),
                detected_topology=_detected_topology(source_record, source_kind),
                is_cropped=source_cropped,
                source_feature_catalog=parsed.source_feature_catalogs[source_record_index],
            )
            _apply_presentation(
                record,
                record_input.presentation,
                record_key=record_key,
            )
            _apply_provenance_annotations(record, item)
            records.append(record)
            provenance.append(item)
    before_collection = records
    records = _apply_collection_options(
        records,
        provenance,
        record_options or RecordCollectionOptions(),
    )
    provenance = [replace(item, has_collection_region=record is not before)
                  for item, record, before in zip(provenance, records, before_collection, strict=True)]
    return ResolvedRecordCollection(tuple(records), tuple(provenance))


def _normalized_cli_values(
    values: Sequence[object] | None,
    *,
    count: int,
    field_name: str,
    default: object,
) -> list[object]:
    normalized = list(values or ())
    if len(normalized) > count:
        raise ValidationError(
            f"Too many {field_name} values (expected at most {count})."
        )
    normalized.extend(default for _ in range(count - len(normalized)))
    return normalized


def _parse_cli_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", "", "none", "null", "-"}:
        return False
    raise ValidationError(f"Invalid reverse_complement value: {value}")


def record_input_manifest_from_paths(
    *,
    gbk_paths: Sequence[str] | None = None,
    gff_paths: Sequence[str] | None = None,
    fasta_paths: Sequence[str] | None = None,
    cardinalities: Sequence[RecordCardinality] | RecordCardinality,
    selectors: Sequence[str] | None = None,
    reverse_flags: Sequence[object] | None = None,
    labels: Sequence[str] | None = None,
    subtitles: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
) -> RecordInputManifest:
    """Translate CLI path lists into unresolved typed record inputs."""

    gbks = tuple(gbk_paths or ())
    gffs = tuple(gff_paths or ())
    fastas = tuple(fasta_paths or ())
    if bool(gbks) == bool(gffs):
        raise ValidationError(
            "Specify either GenBank paths or GFF3/FASTA paths."
        )
    if gffs and len(gffs) != len(fastas):
        raise ValidationError("GFF3 and FASTA path counts must match.")
    count = len(gbks) if gbks else len(gffs)
    if isinstance(cardinalities, RecordCardinality):
        cardinality_values = [cardinalities] * count
    else:
        cardinality_values = list(cardinalities)
        if len(cardinality_values) != count:
            raise ValidationError(
                "Record cardinality count must match the number of input sources."
            )
    if not all(
        isinstance(cardinality, RecordCardinality)
        for cardinality in cardinality_values
    ):
        raise ValidationError(
            "Record cardinalities must be RecordCardinality values."
        )
    selector_values = _normalized_cli_values(
        selectors,
        count=count,
        field_name="record_id",
        default="",
    )
    reverse_values = _normalized_cli_values(
        reverse_flags,
        count=count,
        field_name="reverse_complement",
        default=False,
    )
    records: list[RecordInput] = []
    source_paths: list[str] = []
    for index in range(count):
        source: RecordInputSource
        if gbks:
            source = GenBankInputSource(gbks[index])
            source_paths.append(gbks[index])
        else:
            source = GffFastaInputSource(gffs[index], fastas[index])
            source_paths.append(gffs[index])
        selector = parse_record_selector(str(selector_values[index] or ""))
        cardinality = cardinality_values[index]
        if selector is not None:
            cardinality = RecordCardinality.EXACTLY_ONE
        records.append(
            RecordInput(
                source=source,
                cardinality=cardinality,
                selector=selector,
                presentation=RecordPresentation(
                    reverse_complement=_parse_cli_bool(reverse_values[index]),
                ),
                record_key=f"record-{index + 1}",
            )
        )
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(
            regions=tuple(parse_region_spec(value) for value in (regions or ())),
            labels=tuple(labels or ()),
            subtitles=tuple(subtitles or ()),
        ),
        source_paths=tuple(source_paths),
    )


def record_input_manifest_from_table(path: str) -> RecordInputManifest:
    """Translate one records table into exact-one unresolved row inputs."""

    table = read_records_table(path)
    records: list[RecordInput] = []
    source_paths: list[str] = []
    for index, row in enumerate(table.rows):
        if table.input_kind == "gbk":
            source: RecordInputSource = GenBankInputSource(row.gbk)
            source_paths.append(row.gbk)
        else:
            source = GffFastaInputSource(row.gff, row.fasta)
            source_paths.append(row.gff)
        records.append(
            RecordInput(
                source=source,
                cardinality=RecordCardinality.EXACTLY_ONE,
                selector=parse_record_selector(row.record_id),
                display=RecordDisplayOptions(
                    is_circular=None if row.topology is None else row.topology == "circular",
                    start_coordinate=row.display_start,
                ),
                region=parse_region_spec(row.region) if row.region else None,
                presentation=RecordPresentation(
                    label=row.record_label or None,
                    subtitle=row.record_subtitle or None,
                    reverse_complement=row.reverse_complement,
                    grid_row=row.row,
                    grid_column=row.column,
                ),
                record_key=f"record-{index + 1}",
            )
        )
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(),
        source_paths=tuple(source_paths),
        multi_record_positions=tuple(table.multi_record_positions()),
    )


def resolve_circular_options(
    options: CircularDiagramOptions,
) -> CircularDiagramOptions:
    """Resolve an optional Circular comparison manifest once."""

    if options.conservation_table_file is None:
        return options
    table = read_conservation_table(options.conservation_table_file)
    return replace(
        options,
        conservation_table_file=None,
        conservation_blast_files=tuple(table.conservation_blast_files),
        conservation_fasta_files=(
            tuple(table.comparison_fasta_files)
            if table.comparison_fasta_files is not None
            else None
        ),
        conservation_labels=(
            tuple(table.labels) if table.labels is not None else None
        ),
        conservation_colors=(
            tuple(table.colors) if table.colors is not None else None
        ),
    )


def _resolve_linear_comparison_selector(
    records: Sequence[SeqRecord],
    selector: str,
    *,
    table_path: str,
    row_number: int,
    column: str,
) -> int:
    if selector.startswith("#"):
        try:
            index = int(selector[1:]) - 1
        except ValueError as exc:
            raise ValidationError(
                f"{table_path}: row {row_number}, column {column!r}: "
                f"invalid record selector {selector!r}."
            ) from exc
        if 0 <= index < len(records):
            return index
    else:
        matches = [
            index
            for index, record in enumerate(records)
            if str(record.id) == selector
        ]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValidationError(
                f"{table_path}: row {row_number}, column {column!r}: "
                f"selector {selector!r} matched multiple record IDs; use #index."
            )
    raise ValidationError(
        f"{table_path}: row {row_number}, column {column!r}: unresolved "
        f"record selector {selector!r}."
    )


def resolve_linear_options(
    options: LinearDiagramOptions,
    *,
    records: Sequence[SeqRecord],
    layout: LinearMultiRecordOptions | None,
) -> LinearDiagramOptions:
    """Resolve an optional comparison manifest against displayed records."""

    positions = layout.multi_record_positions if layout is not None else None
    _ordered, rows_by_record = resolve_record_row_positions(records, positions)
    multi_record_rows = len(set(rows_by_record)) < len(records)
    if multi_record_rows and options.blast_files:
        raise ValidationError(
            "-b/--blast is ambiguous when a Linear row contains multiple records; "
            "use a comparison table with explicit query and subject selectors."
        )
    if options.comparison_table_file is None:
        return options
    table = read_comparisons_table(options.comparison_table_file)
    comparisons: list[LinearComparison] = []
    for row in table.rows:
        query_index = _resolve_linear_comparison_selector(
            records,
            row.query,
            table_path=table.table_path,
            row_number=row.row_number,
            column="query",
        )
        subject_index = _resolve_linear_comparison_selector(
            records,
            row.subject,
            table_path=table.table_path,
            row_number=row.row_number,
            column="subject",
        )
        query_row = int(rows_by_record[query_index])
        subject_row = int(rows_by_record[subject_index])
        if query_index == subject_index:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                "query and subject resolved to the same record."
            )
        if query_row == subject_row or abs(query_row - subject_row) != 1:
            topology = (
                "different rows"
                if query_row == subject_row
                else "adjacent rows"
            )
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                f"query row {query_row + 1} and subject row "
                f"{subject_row + 1} must be in {topology}."
            )
        try:
            matches = pd.read_csv(
                row.blast,
                sep="\t",
                comment="#",
                names=COMPARISON_COLUMNS,
            )
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'blast': "
                f"could not parse {row.blast}."
            ) from exc
        comparisons.append(
            LinearComparison(query_index, subject_index, matches)
        )
    return replace(
        options,
        comparison_table_file=None,
        linear_comparisons=tuple(comparisons),
    )


def resolve_circular_batch_outputs(
    policy: CircularBatchOutputPolicy,
    records: Sequence[SeqRecord],
) -> tuple[RenderOutputRequest, ...]:
    """Resolve one output request per displayed record after expansion."""

    prefix = str(policy.output_prefix) if policy.output_prefix is not None else None
    if prefix is None:
        raw_prefixes: list[str] = []
        used: set[str] = set()
        for record in records:
            base = resolve_implicit_record_output_prefix(record.id)
            candidate = base
            suffix = 2
            while candidate in used:
                candidate = f"{base}_{suffix}"
                suffix += 1
            used.add(candidate)
            raw_prefixes.append(candidate)
    elif len(records) == 1:
        raw_prefixes = [prefix]
    else:
        raw_prefixes = [
            f"{prefix}_{index}" for index in range(1, len(records) + 1)
        ]
    outputs: list[RenderOutputRequest] = []
    for raw_prefix in raw_prefixes:
        path = Path(raw_prefix)
        outputs.append(
            RenderOutputRequest(
                output_prefix=path.name,
                output_directory=(
                    path.parent if path.parent != Path(".") else None
                ),
                formats=policy.formats,
                overwrite=policy.overwrite,
                interactive_metadata_policy=policy.interactive_metadata_policy,
            )
        )
    return tuple(outputs)


def resolve_implicit_record_output_prefix(record_id: object) -> str:
    """Require a record-derived output name to be one filename component."""

    raw_record_id = str(record_id)
    try:
        output = RenderOutputRequest(output_prefix=raw_record_id)
    except ValidationError as exc:
        raise ValidationError(
            f"Record ID {raw_record_id!r} cannot be used as an implicit output "
            "filename prefix. Specify an explicit output prefix."
        ) from exc
    return output.output_prefix


__all__ = [
    "RecordInputManifest",
    "ResolvedRecordCollection",
    "ResolvedRecordProvenance",
    "depth_track_inputs_from_cli",
    "record_input_manifest_from_paths",
    "record_input_manifest_from_table",
    "resolve_circular_batch_outputs",
    "resolve_circular_options",
    "resolve_implicit_record_output_prefix",
    "resolve_linear_options",
    "resolve_record_inputs",
]
