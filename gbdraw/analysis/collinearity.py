#!/usr/bin/env python
# coding: utf-8

"""Native gbdraw collinear block calling."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from heapq import heappop, heappush
from typing import Literal, Mapping, NamedTuple, Sequence, SupportsInt, cast, get_args

import pandas as pd
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.analysis.origin_arcs import (
    LinkSide,
    arc_link_endpoints,
    cover_arcs,
    directed_arc,
    rows_inside_records,
)
from gbdraw.analysis.collinearity_units import (
    CollinearityUnit,
    CollinearityUnitIndex,
    CollinearityUnitMode,
    build_collinearity_unit_index,
)
from gbdraw.analysis.protein_colinearity import (
    CdsProtein,
    LosatpCacheManager,
    LosatpRunner,
    OrthogroupMembershipMode,
    OrthogroupResult,
    OrthogroupGraphResult,
    _OrthogroupEdgeIndex,
    _edge_metadata_for_protein_pair,
    _genomic_link_coordinates,
    _orthogroup_member_counts,
    ProteinExtractionResult,
    _consecutive_display_filenames,
    _search_losatp_record_pairs,
    _select_member_candidate_hits_per_query,
    extract_cds_proteins,
    filter_protein_hits_by_thresholds,
    normalize_orthogroup_membership_mode,
    select_reciprocal_best_hit_edges,
    select_rbh_orthogroup_edges_from_directional_hits,
)
from gbdraw.comparisons.losat_jobs import losatp_job_specs
from gbdraw.contract_shapes import Spellings, read_spelling
from gbdraw.core.numeric import compensated_sum
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS

CollinearityOrientation = Literal["plus", "minus"]
CollinearityBlockKind = Literal["syntenic", "cluster", "singleton"]
CollinearityColorMode = Literal["average_identity", "orientation", "orientation_identity"]
CollinearityAnchorMode = Literal["all", "rbh"]
CollinearityMergeOrientation = Literal["strand", "order"]
CollinearitySearchScope = Literal["adjacent", "all"]

COLLINEARITY_METADATA_COLUMNS = (
    "collinearity_block_id",
    "collinearity_block_kind",
    "collinearity_orientation",
    "collinearity_block_score",
    "collinearity_anchor_index",
    "collinearity_anchor_count",
    "collinearity_color_mode",
    "orthogroup_id",
    "query_unit_id",
    "subject_unit_id",
    "query_unit_kind",
    "subject_unit_kind",
    "query_locus_id",
    "subject_locus_id",
    "query_display_name",
    "subject_display_name",
    "query_protein_id",
    "subject_protein_id",
    "query_feature_svg_id",
    "subject_feature_svg_id",
    "query_view_feature_svg_id",
    "subject_view_feature_svg_id",
    "query_feature_index",
    "subject_feature_index",
    "rbh_orthogroup_id",
    "ortholog_path_id",
    "edge_kind",
    "render_role",
    "query_orthogroup_representative",
    "subject_orthogroup_representative",
    "query_orthogroup_member_count",
    "subject_orthogroup_member_count",
)
COLLINEARITY_COMPARISON_COLUMNS = tuple(COMPARISON_COLUMNS) + COLLINEARITY_METADATA_COLUMNS
COLLINEARITY_COLOR_MODES: tuple[CollinearityColorMode, ...] = get_args(CollinearityColorMode)
COLLINEARITY_ANCHOR_MODES: tuple[CollinearityAnchorMode, ...] = get_args(CollinearityAnchorMode)
COLLINEARITY_SEARCH_SCOPES: tuple[CollinearitySearchScope, ...] = get_args(CollinearitySearchScope)
COLLINEARITY_MERGE_ORIENTATIONS: tuple[CollinearityMergeOrientation, ...] = get_args(CollinearityMergeOrientation)
# The other spellings each normalizer below reads (gbdraw.contract_shapes.read_spelling);
# the request description names them for the Web's Session Load.
COLLINEARITY_COLOR_MODE_SPELLINGS: Spellings = {
    "aliases": {"identity": "average_identity"}, "foldHyphens": True, "blank": "orientation",
}
COLLINEARITY_ANCHOR_MODE_SPELLINGS: Spellings = {
    "aliases": {
        "all_hits": "all",
        "raw": "all",
        "top_n": "all",
        "topn": "all",
        "reciprocal_best": "rbh",
        "strict_rbh": "rbh",
    },
    "foldHyphens": True,
    "blank": "rbh",
}
COLLINEARITY_SEARCH_SCOPE_SPELLINGS: Spellings = {"foldHyphens": True, "blank": "adjacent"}
# "either" was a third value that behaved as "strand" (OV-2301).
COLLINEARITY_MERGE_ORIENTATION_SPELLINGS: Spellings = {"aliases": {"either": "strand"}, "blank": "strand"}
# The anchor mode "one_to_one" (forward mutual best hits) was removed.
# A stored request or Web draft that holds it, or one of its spellings that
# requests before schema 9 hold, is read as "rbh", the default, and its reader
# reports the replacement.
RETIRED_COLLINEARITY_ANCHOR_MODES: Mapping[str, CollinearityAnchorMode] = {
    spelling: "rbh"
    for spelling in ("one_to_one", "one2one", "one_to_ones", "mutual_best", "top1", "top_1")
}


@dataclass(frozen=True)
class CollinearityAnchor:
    query_protein_id: str
    subject_protein_id: str
    query_record_index: int
    subject_record_index: int
    query_order: int
    subject_order: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    identity: float
    evalue: float
    bitscore: float
    alignment_length: int
    query_feature_svg_id: str
    subject_feature_svg_id: str
    source: str
    query_unit_id: str
    subject_unit_id: str
    query_unit_kind: Literal["locus", "cds"]
    subject_unit_kind: Literal["locus", "cds"]
    query_locus_id: str | None
    subject_locus_id: str | None
    query_display_name: str
    subject_display_name: str
    query_strand: int | None = None
    subject_strand: int | None = None
    orthogroup_id: str = ""
    rbh_orthogroup_id: str = ""
    ortholog_path_id: str = ""
    edge_kind: str = ""
    render_role: str = ""
    query_orthogroup_representative: bool = False
    subject_orthogroup_representative: bool = False
    query_orthogroup_member_count: int = 0
    subject_orthogroup_member_count: int = 0
    query_view_feature_svg_id: str = ""
    subject_view_feature_svg_id: str = ""
    query_feature_index: int | None = None
    subject_feature_index: int | None = None


@dataclass(frozen=True)
class CollinearityBlock:
    block_id: str
    query_record_index: int
    subject_record_index: int
    orientation: CollinearityOrientation
    score: float
    anchors: tuple[CollinearityAnchor, ...]
    # The block caller writes only "cluster" and "singleton" and leaves
    # block_evalue and label None (OV-2300). The fields and "syntenic" stay
    # because saved collinearity-result resources carry them and their reader
    # requires every field.
    kind: CollinearityBlockKind = "syntenic"
    block_evalue: float | None = None
    label: str | None = None


@dataclass(frozen=True)
class CollinearityResult:
    blocks: tuple[CollinearityBlock, ...]
    unblocked_anchors: tuple[CollinearityAnchor, ...] = ()
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None


@dataclass(frozen=True)
class LosslessCollinearityParameters:
    min_anchors: int = 1
    max_unit_gap: int = 0
    max_diagonal_drift: int = 0
    max_conflicts: int = 1
    merge_orientation: CollinearityMergeOrientation = "strand"

    def __post_init__(self) -> None:
        # Every reader of the parameters reads the old spelling "either" as "strand".
        normalized = read_spelling(
            self.merge_orientation,
            COLLINEARITY_MERGE_ORIENTATIONS,
            COLLINEARITY_MERGE_ORIENTATION_SPELLINGS,
        )
        if normalized is not None:
            object.__setattr__(self, "merge_orientation", normalized)

    def validate(self) -> None:
        if int(self.min_anchors) <= 0:
            raise ValidationError("collinear_min_anchors must be > 0")
        if int(self.max_unit_gap) < 0:
            raise ValidationError("collinear_max_unit_gap must be >= 0")
        if int(self.max_diagonal_drift) < 0:
            raise ValidationError("collinear_max_diagonal_drift must be >= 0")
        if int(self.max_conflicts) < 0:
            raise ValidationError("collinear_max_conflicts must be >= 0")
        if str(self.merge_orientation) not in COLLINEARITY_MERGE_ORIENTATIONS:
            raise ValidationError(
                "collinear_merge_orientation must be one of: "
                + ", ".join(COLLINEARITY_MERGE_ORIENTATIONS),
                diagnostic={"code": "INPUT_INVALID", "field": "collinear_merge_orientation", "reason": "COLLINEAR_MERGE_ORIENTATION"},
            )


def normalize_collinearity_color_mode(mode: str | None) -> CollinearityColorMode:
    normalized = read_spelling(mode, COLLINEARITY_COLOR_MODES, COLLINEARITY_COLOR_MODE_SPELLINGS)
    if normalized is None:
        raise ValidationError(
            "collinear_color_mode must be one of: "
            + ", ".join(COLLINEARITY_COLOR_MODES),
            diagnostic={"code": "INPUT_INVALID", "field": "collinear_color_mode", "reason": "COLLINEAR_COLOR_MODE"},
        )
    return cast(CollinearityColorMode, normalized)


def normalize_collinearity_anchor_mode(mode: str | None) -> CollinearityAnchorMode:
    normalized = read_spelling(mode, COLLINEARITY_ANCHOR_MODES, COLLINEARITY_ANCHOR_MODE_SPELLINGS)
    if normalized is None:
        raise ValidationError(
            "collinear_anchor_mode must be one of: "
            + ", ".join(COLLINEARITY_ANCHOR_MODES),
            diagnostic={"code": "INPUT_INVALID", "field": "collinear_anchor_mode", "reason": "COLLINEAR_ANCHOR_MODE"},
        )
    return cast(CollinearityAnchorMode, normalized)


def normalize_collinearity_search_scope(scope: str | None) -> CollinearitySearchScope:
    normalized = read_spelling(scope, COLLINEARITY_SEARCH_SCOPES, COLLINEARITY_SEARCH_SCOPE_SPELLINGS)
    if normalized is None:
        raise ValidationError(
            "collinear_search_scope must be one of: "
            + ", ".join(COLLINEARITY_SEARCH_SCOPES),
            diagnostic={"code": "INPUT_INVALID", "field": "collinear_search_scope", "reason": "ADJACENT_ALL"},
        )
    return cast(CollinearitySearchScope, normalized)


def iter_collinearity_search_pairs(
    record_count: int,
    *,
    scope: CollinearitySearchScope | str,
) -> tuple[tuple[int, int], ...]:
    """Return the unordered record pairs that ``scope`` searches for Collinear evidence.

    A view of the LOSAT job plan (:func:`gbdraw.comparisons.losat_jobs.losatp_job_specs`)
    for consecutive records with similarity-group inference on, the only case
    in which Search between applies.
    """

    specs = losatp_job_specs(
        "collinear",
        record_count=record_count,
        infer_orthogroups=True,
        search_scope=normalize_collinearity_search_scope(str(scope)),
    )
    return tuple(dict.fromkeys(
        cast(tuple[int, int], tuple(sorted(pair))) for pair in specs if pair[0] != pair[1]
    ))


def _float_from_row(row: object, column: str, default: float = 0.0) -> float:
    try:
        return float(getattr(row, column))
    except (TypeError, ValueError, AttributeError):
        return float(default)


def _int_from_row(row: object, column: str, default: int = 0) -> int:
    try:
        return int(getattr(row, column))
    except (TypeError, ValueError, AttributeError):
        return int(default)


def _unit_genomic_link_coordinates(unit: CollinearityUnit) -> tuple[int, int]:
    return arc_link_endpoints((int(unit.start), int(unit.end)), unit.strand != -1)


def _normalized_strand(value: object) -> int | None:
    try:
        # Values that int() rejects raise TypeError or ValueError, handled below.
        strand = int(cast(SupportsInt, value))
    except (TypeError, ValueError):
        return None
    return strand if strand in {-1, 1} else None


def _view_feature_svg_id(protein: CdsProtein) -> str:
    return str(
        protein.view_feature_svg_id
        or protein.feature_svg_id
        or ""
    )


def _anchor_strand_orientation(anchor: CollinearityAnchor) -> CollinearityOrientation | None:
    query_strand = _normalized_strand(anchor.query_strand)
    subject_strand = _normalized_strand(anchor.subject_strand)
    if query_strand is None or subject_strand is None:
        return None
    return "plus" if query_strand == subject_strand else "minus"


def _path_sort_key(
    anchor: CollinearityAnchor,
    orientation: CollinearityOrientation,
) -> tuple[int, int, str, str, str]:
    return (
        int(anchor.query_order),
        _orientation_subject_value(anchor, orientation),
        str(anchor.orthogroup_id),
        str(anchor.query_protein_id),
        str(anchor.subject_protein_id),
    )


def _path_sorted_anchors(
    anchors: Sequence[CollinearityAnchor],
    orientation: CollinearityOrientation,
) -> tuple[CollinearityAnchor, ...]:
    return tuple(sorted(anchors, key=lambda anchor: _path_sort_key(anchor, orientation)))


def _orientation_subject_value(
    anchor: CollinearityAnchor,
    orientation: CollinearityOrientation,
) -> int:
    subject_order = int(anchor.subject_order)
    return -subject_order if orientation == "minus" else subject_order


def _orientation_diagonal(
    anchor: CollinearityAnchor,
    orientation: CollinearityOrientation,
) -> int:
    if orientation == "minus":
        return int(anchor.query_order) + int(anchor.subject_order)
    return int(anchor.subject_order) - int(anchor.query_order)


def _path_gap(
    previous: CollinearityAnchor,
    current: CollinearityAnchor,
    orientation: CollinearityOrientation,
) -> tuple[int, int]:
    query_gap = int(current.query_order) - int(previous.query_order) - 1
    if orientation == "minus":
        subject_gap = int(previous.subject_order) - int(current.subject_order) - 1
    else:
        subject_gap = int(current.subject_order) - int(previous.subject_order) - 1
    return query_gap, subject_gap


def _anchors_are_path_compatible(
    previous: CollinearityAnchor,
    current: CollinearityAnchor,
    *,
    orientation: CollinearityOrientation,
    max_gap: int,
    max_diagonal_drift: int,
) -> bool:
    query_gap, subject_gap = _path_gap(previous, current, orientation)
    if query_gap < 0 or subject_gap < 0:
        return False
    if query_gap > int(max_gap) or subject_gap > int(max_gap):
        return False
    diagonal_drift = abs(
        _orientation_diagonal(current, orientation)
        - _orientation_diagonal(previous, orientation)
    )
    return diagonal_drift <= int(max_diagonal_drift)


def _anchor_sort_key(anchor: CollinearityAnchor) -> tuple[int, int, int, int, str, str, str]:
    return (
        int(anchor.query_record_index),
        int(anchor.subject_record_index),
        int(anchor.query_order),
        int(anchor.subject_order),
        str(anchor.orthogroup_id),
        str(anchor.query_protein_id),
        str(anchor.subject_protein_id),
    )


def _protein_order_by_id(protein_map: Mapping[str, CdsProtein]) -> dict[str, int]:
    grouped: dict[int, list[CdsProtein]] = {}
    for protein in protein_map.values():
        grouped.setdefault(int(protein.record_index), []).append(protein)
    order_by_id: dict[str, int] = {}
    for proteins in grouped.values():
        for order, protein in enumerate(
            sorted(
                proteins,
                key=lambda item: (
                    int(item.start),
                    int(item.end),
                    int(item.feature_index),
                    str(item.protein_id),
                ),
            )
        ):
            order_by_id[str(protein.protein_id)] = order
    return order_by_id


def _orthogroup_member_counts_by_record(
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
) -> dict[tuple[str, int], int]:
    counts: dict[tuple[str, int], int] = {}
    if orthogroups is None:
        return counts
    for orthogroup_id, members in orthogroups.orthogroups.items():
        for record_index, count in _orthogroup_member_counts(members).items():
            key = (str(orthogroup_id), record_index)
            counts[key] = counts.get(key, 0) + count
    return counts


def _orthogroup_edge_metadata_for_anchor(
    query_id: str,
    subject_id: str,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    edge_indexes: dict[str, _OrthogroupEdgeIndex],
) -> dict[str, object]:
    metadata = {
        "orthogroup_id": "",
        "rbh_orthogroup_id": "",
        "ortholog_path_id": "",
        "edge_kind": "",
        "render_role": "",
        "query_orthogroup_representative": False,
        "subject_orthogroup_representative": False,
    }
    if orthogroups is None:
        return metadata
    query_member = orthogroups.member_by_protein_id.get(query_id)
    subject_member = orthogroups.member_by_protein_id.get(subject_id)
    if query_member is not None:
        metadata["query_orthogroup_representative"] = bool(query_member.representative)
    if subject_member is not None:
        metadata["subject_orthogroup_representative"] = bool(subject_member.representative)
    orthogroup_id = ""
    if query_member is not None and subject_member is not None:
        if query_member.orthogroup_id == subject_member.orthogroup_id:
            orthogroup_id = str(query_member.orthogroup_id)
    metadata["orthogroup_id"] = orthogroup_id
    metadata.update(_edge_metadata_for_protein_pair(
        orthogroups, orthogroup_id, query_id, subject_id, edge_indexes,
    ))
    return metadata


def _protein_locus_metadata(protein: CdsProtein) -> str | None:
    for attr in ("locus_tag", "gene_id", "old_locus_tag", "gene"):
        text = str(getattr(protein, attr, None) or "").strip()
        if text:
            return text
    return None


def _lossless_anchor_from_edge_row(
    row: object,
    *,
    query_record_index: int,
    subject_record_index: int,
    protein_map: Mapping[str, CdsProtein],
    order_by_id: Mapping[str, int],
    unit_index: CollinearityUnitIndex | None,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    member_counts_by_record: Mapping[tuple[str, int], int],
    edge_indexes: dict[str, _OrthogroupEdgeIndex],
) -> CollinearityAnchor | None:
    query_id = str(getattr(row, "query"))
    subject_id = str(getattr(row, "subject"))
    query_protein = protein_map.get(query_id)
    subject_protein = protein_map.get(subject_id)
    if query_protein is None or subject_protein is None:
        return None
    query_unit = (
        unit_index.unit_by_protein_id.get(query_id)
        if unit_index is not None
        else None
    )
    subject_unit = (
        unit_index.unit_by_protein_id.get(subject_id)
        if unit_index is not None
        else None
    )
    query_representative = (
        protein_map.get(query_unit.representative_protein_id, query_protein)
        if query_unit is not None
        else query_protein
    )
    subject_representative = (
        protein_map.get(subject_unit.representative_protein_id, subject_protein)
        if subject_unit is not None
        else subject_protein
    )
    edge_metadata = _orthogroup_edge_metadata_for_anchor(
        query_id,
        subject_id,
        orthogroups,
        edge_indexes,
    )
    orthogroup_id = str(edge_metadata["orthogroup_id"])
    qstart, qend = (
        _unit_genomic_link_coordinates(query_unit)
        if query_unit is not None
        else _genomic_link_coordinates(query_protein)
    )
    sstart, send = (
        _unit_genomic_link_coordinates(subject_unit)
        if subject_unit is not None
        else _genomic_link_coordinates(subject_protein)
    )
    return CollinearityAnchor(
        query_protein_id=query_id,
        subject_protein_id=subject_id,
        query_record_index=int(query_record_index),
        subject_record_index=int(subject_record_index),
        query_order=int(
            query_unit.order
            if query_unit is not None
            else order_by_id.get(query_id, query_protein.feature_index)
        ),
        subject_order=int(
            subject_unit.order
            if subject_unit is not None
            else order_by_id.get(subject_id, subject_protein.feature_index)
        ),
        query_start=qstart,
        query_end=qend,
        subject_start=sstart,
        subject_end=send,
        query_strand=_normalized_strand(
            query_unit.strand if query_unit is not None else query_protein.strand
        ),
        subject_strand=_normalized_strand(
            subject_unit.strand if subject_unit is not None else subject_protein.strand
        ),
        identity=_float_from_row(row, "identity", 0.0),
        evalue=_float_from_row(row, "evalue", 1.0),
        bitscore=_float_from_row(row, "bitscore", 0.0),
        alignment_length=_int_from_row(row, "alignment_length", 0),
        query_feature_svg_id=str(query_representative.feature_svg_id or ""),
        subject_feature_svg_id=str(subject_representative.feature_svg_id or ""),
        query_view_feature_svg_id=_view_feature_svg_id(query_representative),
        subject_view_feature_svg_id=_view_feature_svg_id(subject_representative),
        query_feature_index=query_representative.feature_index,
        subject_feature_index=subject_representative.feature_index,
        source="orthogroup_display_edge" if orthogroups is not None else "protein_hit",
        query_unit_id=query_unit.unit_id if query_unit is not None else query_id,
        subject_unit_id=subject_unit.unit_id if subject_unit is not None else subject_id,
        query_unit_kind=query_unit.unit_kind if query_unit is not None else "cds",
        subject_unit_kind=subject_unit.unit_kind if subject_unit is not None else "cds",
        query_locus_id=(
            query_unit.locus_id
            if query_unit is not None
            else _protein_locus_metadata(query_protein)
        ),
        subject_locus_id=(
            subject_unit.locus_id
            if subject_unit is not None
            else _protein_locus_metadata(subject_protein)
        ),
        query_display_name=(
            query_unit.display_name
            if query_unit is not None
            else str(query_protein.label or query_id)
        ),
        subject_display_name=(
            subject_unit.display_name
            if subject_unit is not None
            else str(subject_protein.label or subject_id)
        ),
        orthogroup_id=orthogroup_id,
        rbh_orthogroup_id=str(edge_metadata["rbh_orthogroup_id"]),
        ortholog_path_id=str(edge_metadata["ortholog_path_id"]),
        edge_kind=str(edge_metadata["edge_kind"]),
        render_role=str(edge_metadata["render_role"]),
        query_orthogroup_representative=bool(edge_metadata["query_orthogroup_representative"]),
        subject_orthogroup_representative=bool(edge_metadata["subject_orthogroup_representative"]),
        query_orthogroup_member_count=int(
            member_counts_by_record.get((orthogroup_id, int(query_record_index)), 0)
        ),
        subject_orthogroup_member_count=int(
            member_counts_by_record.get((orthogroup_id, int(subject_record_index)), 0)
        ),
    )


def orthogroup_edges_to_lossless_collinearity_anchors(
    adjacent_edges_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    *,
    unit_index: CollinearityUnitIndex | None = None,
) -> tuple[CollinearityAnchor, ...]:
    """Convert adjacent Orthogroup display edges one-to-one into anchors."""

    order_by_id = _protein_order_by_id(protein_map)
    member_counts_by_record = _orthogroup_member_counts_by_record(orthogroups)
    edge_indexes: dict[str, _OrthogroupEdgeIndex] = {}
    anchors: list[CollinearityAnchor] = []
    missing_ids: set[str] = set()
    for query_record_index, subject_record_index in sorted(adjacent_edges_by_pair):
        edges = adjacent_edges_by_pair[(query_record_index, subject_record_index)]
        if edges is None or edges.empty:
            continue
        for row in edges.itertuples(index=False):
            query_id = str(getattr(row, "query"))
            subject_id = str(getattr(row, "subject"))
            if query_id not in protein_map:
                missing_ids.add(query_id)
            if subject_id not in protein_map:
                missing_ids.add(subject_id)
            anchor = _lossless_anchor_from_edge_row(
                row,
                query_record_index=int(query_record_index),
                subject_record_index=int(subject_record_index),
                protein_map=protein_map,
                order_by_id=order_by_id,
                unit_index=unit_index,
                orthogroups=orthogroups,
                member_counts_by_record=member_counts_by_record,
                edge_indexes=edge_indexes,
            )
            if anchor is not None:
                anchors.append(anchor)
    if missing_ids:
        raise ParseError(
            "Orthogroup display edge contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )
    if unit_index is not None:
        strongest_by_unit_pair: dict[
            tuple[int, int, str, str], CollinearityAnchor
        ] = {}
        for anchor in anchors:
            key = (
                int(anchor.query_record_index),
                int(anchor.subject_record_index),
                str(anchor.query_unit_id),
                str(anchor.subject_unit_id),
            )
            current = strongest_by_unit_pair.get(key)
            rank = (
                -float(anchor.bitscore),
                float(anchor.evalue),
                -float(anchor.identity),
                -int(anchor.alignment_length),
                str(anchor.query_protein_id),
                str(anchor.subject_protein_id),
            )
            current_rank = (
                (
                    -float(current.bitscore),
                    float(current.evalue),
                    -float(current.identity),
                    -int(current.alignment_length),
                    str(current.query_protein_id),
                    str(current.subject_protein_id),
                )
                if current is not None
                else None
            )
            if current_rank is None or rank < current_rank:
                strongest_by_unit_pair[key] = anchor
        anchors = list(strongest_by_unit_pair.values())
    return tuple(sorted(anchors, key=_anchor_sort_key))


def _resolve_lossless_params(
    params: LosslessCollinearityParameters | None,
) -> LosslessCollinearityParameters:
    resolved = params or LosslessCollinearityParameters()
    resolved.validate()
    return resolved


def _lossless_anchor_supports_orientation(
    anchor: CollinearityAnchor,
    orientation: CollinearityOrientation,
    params: LosslessCollinearityParameters,
) -> bool:
    if params.merge_orientation == "order":
        return True
    strand_orientation = _anchor_strand_orientation(anchor)
    return strand_orientation is None or strand_orientation == orientation


def _lossless_anchors_are_compatible(
    previous: CollinearityAnchor,
    current: CollinearityAnchor,
    *,
    orientation: CollinearityOrientation,
    params: LosslessCollinearityParameters,
) -> bool:
    if int(previous.query_record_index) != int(current.query_record_index):
        return False
    if int(previous.subject_record_index) != int(current.subject_record_index):
        return False
    if not _lossless_anchor_supports_orientation(previous, orientation, params):
        return False
    if not _lossless_anchor_supports_orientation(current, orientation, params):
        return False
    return _anchors_are_path_compatible(
        previous,
        current,
        orientation=orientation,
        max_gap=int(params.max_unit_gap),
        max_diagonal_drift=int(params.max_diagonal_drift),
    )


def _lossless_singleton_orientation(anchor: CollinearityAnchor) -> CollinearityOrientation:
    return _anchor_strand_orientation(anchor) or "plus"


def _anchor_bitscore_total(anchors: Sequence[CollinearityAnchor]) -> float:
    """Total anchor bitscore, with the same bits on every Python (OV-1237)."""

    return compensated_sum(float(anchor.bitscore) for anchor in anchors)


def _lossless_block_from_anchors(
    *,
    block_id: str,
    pair: tuple[int, int],
    orientation: CollinearityOrientation,
    anchors: Sequence[CollinearityAnchor],
) -> CollinearityBlock:
    ordered = _path_sorted_anchors(anchors, orientation)
    return CollinearityBlock(
        block_id=block_id,
        query_record_index=int(pair[0]),
        subject_record_index=int(pair[1]),
        orientation=orientation,
        score=_anchor_bitscore_total(ordered),
        anchors=ordered,
        kind="cluster" if len(ordered) > 1 else "singleton",
    )


class _PairOrderIndex(NamedTuple):
    """Every anchor of one record pair: the query orders that hold an anchor,
    ascending, and the sorted subject orders at each of them."""

    query_orders: list[int]
    subject_orders: list[list[int]]


def _pair_order_index(anchors: Sequence[CollinearityAnchor]) -> _PairOrderIndex:
    by_query: dict[int, list[int]] = {}
    for anchor in anchors:
        by_query.setdefault(int(anchor.query_order), []).append(int(anchor.subject_order))
    query_orders = sorted(by_query)
    return _PairOrderIndex(query_orders, [sorted(by_query[order]) for order in query_orders])


def _lossless_step_conflicts_exceed(
    previous: CollinearityAnchor,
    current: CollinearityAnchor,
    index: _PairOrderIndex,
    max_conflicts: int,
) -> bool:
    """Whether more than ``max_conflicts`` anchors of the pair lie strictly
    inside the box between two chain neighbours, in both order intervals.

    Only the query orders that hold an anchor are visited, so a wide box costs
    the anchors near it, not its width.
    """

    query_low, query_high = sorted((int(previous.query_order), int(current.query_order)))
    subject_low, subject_high = sorted((int(previous.subject_order), int(current.subject_order)))
    if query_high - query_low < 2 or subject_high - subject_low < 2:
        return False
    conflicts = 0
    for position in range(
        bisect_right(index.query_orders, query_low), bisect_left(index.query_orders, query_high)
    ):
        subject_orders = index.subject_orders[position]
        conflicts += bisect_left(subject_orders, subject_high) - bisect_right(subject_orders, subject_low)
        if conflicts > max_conflicts:
            return True
    return False


def _oriented_subjects_between(
    index: _PairOrderIndex,
    query_order: int,
    orientation: CollinearityOrientation,
    low: int,
    high: int,
) -> list[int]:
    """Oriented subject orders strictly between ``low`` and ``high`` of every
    anchor of the pair at one query order."""

    position = bisect_left(index.query_orders, query_order)
    if position == len(index.query_orders) or index.query_orders[position] != query_order:
        return []
    subject_orders = index.subject_orders[position]
    if orientation == "minus":
        start = bisect_right(subject_orders, -high)
        end = bisect_left(subject_orders, -low)
        return [-subject for subject in subject_orders[start:end]]
    return subject_orders[bisect_right(subject_orders, low):bisect_left(subject_orders, high)]

_CHAIN_ORIENTATIONS: tuple[CollinearityOrientation, ...] = ("plus", "minus")


@dataclass
class _OrientationChains:
    """Longest-chain state of one record pair in one orientation.

    Nodes are the pair's anchors that support the orientation, ranked by the
    path sort key; ``predecessors[rank]`` lists, in rank order, the nodes that
    may precede that node in a chain. ``count``/``score``/``previous`` hold the
    best chain ending at each node: most anchors, then the highest summed
    bitscore, then the earliest predecessor rank.
    """

    orientation: CollinearityOrientation
    anchor_indexes: list[int]
    predecessors: list[list[int]]
    count: list[int]
    score: list[float]
    previous: list[int]
    followers: list[set[int]]


def _orientation_chains(
    anchors: Sequence[CollinearityAnchor],
    orientation: CollinearityOrientation,
    params: LosslessCollinearityParameters,
    order_index: _PairOrderIndex,
) -> _OrientationChains:
    anchor_indexes = sorted(
        (
            index
            for index, anchor in enumerate(anchors)
            if _lossless_anchor_supports_orientation(anchor, orientation, params)
        ),
        key=lambda index: _path_sort_key(anchors[index], orientation),
    )
    max_gap = int(params.max_unit_gap)
    max_conflicts = int(params.max_conflicts)
    max_drift = int(params.max_diagonal_drift)
    # Query orders placed so far, ascending, and at each one its oriented
    # subject orders and ranks, both in rank order.
    placed_queries: list[int] = []
    placed: list[tuple[list[int], list[int]]] = []
    predecessors: list[list[int]] = []
    for rank, anchor_index in enumerate(anchor_indexes):
        anchor = anchors[anchor_index]
        query_order = int(anchor.query_order)
        subject_value = _orientation_subject_value(anchor, orientation)
        found: list[int] = []
        # A compatible predecessor lies at most max_gap + 1 units before the
        # node in both orders. Visit only the placed query orders in that
        # window, nearest first, and bisect their oriented subject orders, so
        # the cost is the anchors in the window, not the width of max_gap.
        # A predecessor dx units back has dy >= dx - max_drift, so the anchors
        # of nearer units whose oriented subject order lies above
        # subject_value - dx + max_drift are inside the box of every
        # predecessor from dx back; once more than max_conflicts are, the
        # farther units cannot hold one.
        enclosed = 0
        inside: list[int] = []  # oriented subject orders, negated (max-heap)
        for column in range(
            bisect_left(placed_queries, query_order) - 1,
            bisect_left(placed_queries, query_order - max_gap - 1) - 1,
            -1,
        ):
            previous_query = placed_queries[column]
            floor = subject_value - (query_order - previous_query) + max_drift
            while inside and -inside[0] > floor:
                heappop(inside)
                enclosed += 1
            if enclosed > max_conflicts:
                break
            subject_values, ranks = placed[column]
            for position in range(
                bisect_left(subject_values, subject_value - max_gap - 1),
                bisect_left(subject_values, subject_value),
            ):
                previous_rank = ranks[position]
                previous = anchors[anchor_indexes[previous_rank]]
                if _lossless_anchors_are_compatible(
                    previous, anchor, orientation=orientation, params=params
                ) and not _lossless_step_conflicts_exceed(
                    previous, anchor, order_index, max_conflicts
                ):
                    found.append(previous_rank)
            for oriented in _oriented_subjects_between(
                order_index, previous_query, orientation,
                subject_value - max_gap - 1 + max_drift, subject_value,
            ):
                heappush(inside, -oriented)
        found.sort()
        predecessors.append(found)
        if not placed_queries or placed_queries[-1] != query_order:
            placed_queries.append(query_order)
            placed.append(([], []))
        bucket = placed[-1]
        bucket[0].append(subject_value)
        bucket[1].append(rank)
    size = len(anchor_indexes)
    return _OrientationChains(
        orientation=orientation,
        anchor_indexes=anchor_indexes,
        predecessors=predecessors,
        count=[0] * size,
        score=[0.0] * size,
        previous=[-1] * size,
        followers=[set() for _ in range(size)],
    )


def _update_best_chain_end(
    chains: _OrientationChains,
    rank: int,
    anchors: Sequence[CollinearityAnchor],
    alive: Sequence[bool],
) -> None:
    best, best_count, best_score = -1, 0, 0.0
    for previous_rank in chains.predecessors[rank]:
        if not alive[chains.anchor_indexes[previous_rank]]:
            continue
        count = chains.count[previous_rank]
        score = chains.score[previous_rank]
        if count > best_count or (count == best_count and score > best_score):
            best, best_count, best_score = previous_rank, count, score
    old = chains.previous[rank]
    if old >= 0:
        chains.followers[old].discard(rank)
    if best >= 0:
        chains.followers[best].add(rank)
    chains.previous[rank] = best
    chains.count[rank] = best_count + 1
    chains.score[rank] = best_score + float(anchors[chains.anchor_indexes[rank]].bitscore)


def _lossless_chains_for_pair(
    pair: tuple[int, int],
    anchors: Sequence[CollinearityAnchor],
    params: LosslessCollinearityParameters,
) -> list[CollinearityBlock]:
    """Split one record pair's anchors into longest compatible chains.

    Each orientation's anchors, in query order, form a DAG whose edges join
    path-compatible neighbours with at most ``max_conflicts`` other anchors in
    the box between them. The longest chain over both orientations is taken
    out and the rest are chained again until no chain of two anchors remains;
    every anchor left is a singleton block.
    """

    order_index = _pair_order_index(anchors)
    alive = [True] * len(anchors)
    graphs = [
        _orientation_chains(anchors, orientation, params, order_index)
        for orientation in _CHAIN_ORIENTATIONS
    ]
    rank_by_anchor = [
        {anchor_index: rank for rank, anchor_index in enumerate(graph.anchor_indexes)}
        for graph in graphs
    ]
    # Chain ends by (most anchors, highest score, plus before minus, earliest rank);
    # an entry is stale once its node is taken or its best chain changed.
    heap: list[tuple[int, float, int, int]] = []
    for graph_index, graph in enumerate(graphs):
        for rank in range(len(graph.anchor_indexes)):
            _update_best_chain_end(graph, rank, anchors, alive)
            if graph.count[rank] > 1:
                heappush(heap, (-graph.count[rank], -graph.score[rank], graph_index, rank))

    blocks: list[CollinearityBlock] = []
    while heap:
        negative_count, negative_score, graph_index, rank = heappop(heap)
        graph = graphs[graph_index]
        if (
            not alive[graph.anchor_indexes[rank]]
            or graph.count[rank] != -negative_count
            or graph.score[rank] != -negative_score
        ):
            continue
        taken: list[int] = []
        while rank >= 0:
            taken.append(graph.anchor_indexes[rank])
            rank = graph.previous[rank]
        for anchor_index in taken:
            alive[anchor_index] = False
        blocks.append(
            _lossless_block_from_anchors(
                block_id=f"pending_{len(blocks) + 1:04d}",
                pair=pair,
                orientation=graph.orientation,
                anchors=[anchors[anchor_index] for anchor_index in taken],
            )
        )
        for other_index, other in enumerate(graphs):
            # Only chains that ran through a taken anchor can change.
            stale: set[int] = set()
            pending = [
                rank_by_anchor[other_index][anchor_index]
                for anchor_index in taken
                if anchor_index in rank_by_anchor[other_index]
            ]
            while pending:
                for follower in other.followers[pending.pop()]:
                    if follower not in stale and alive[other.anchor_indexes[follower]]:
                        stale.add(follower)
                        pending.append(follower)
            for stale_rank in sorted(stale):
                _update_best_chain_end(other, stale_rank, anchors, alive)
                if other.count[stale_rank] > 1:
                    heappush(
                        heap,
                        (-other.count[stale_rank], -other.score[stale_rank], other_index, stale_rank),
                    )

    for anchor_index in sorted(
        (index for index, is_alive in enumerate(alive) if is_alive),
        key=lambda index: _anchor_sort_key(anchors[index]),
    ):
        anchor = anchors[anchor_index]
        blocks.append(
            _lossless_block_from_anchors(
                block_id=f"pending_{len(blocks) + 1:04d}",
                pair=pair,
                orientation=_lossless_singleton_orientation(anchor),
                anchors=[anchor],
            )
        )
    return blocks


def _filter_blocks_by_min_anchors(
    blocks: Sequence[CollinearityBlock],
    *,
    min_anchors: int,
) -> tuple[tuple[CollinearityBlock, ...], tuple[CollinearityAnchor, ...]]:
    threshold = int(min_anchors)
    if threshold <= 0:
        raise ValidationError("collinear_min_anchors must be > 0")
    kept: list[CollinearityBlock] = []
    dropped_anchors: list[CollinearityAnchor] = []
    for block in blocks:
        if len(block.anchors) >= threshold:
            kept.append(block)
        else:
            dropped_anchors.extend(block.anchors)
    return tuple(kept), tuple(sorted(dropped_anchors, key=_anchor_sort_key))


def _final_block_sort_key(block: CollinearityBlock) -> tuple[int, int, int, int, int, str]:
    anchors = _path_sorted_anchors(block.anchors, block.orientation)
    first = anchors[0] if anchors else None
    kind_rank = 1 if block.kind == "singleton" else 0
    return (
        int(block.query_record_index),
        int(block.subject_record_index),
        int(first.query_order if first is not None else 0),
        _orientation_subject_value(first, block.orientation) if first is not None else 0,
        kind_rank,
        str(block.block_id),
    )


def _renumber_lossless_blocks(blocks: Sequence[CollinearityBlock]) -> tuple[CollinearityBlock, ...]:
    renumbered: list[CollinearityBlock] = []
    for index, block in enumerate(sorted(blocks, key=_final_block_sort_key), start=1):
        renumbered.append(
            _lossless_block_from_anchors(
                block_id=f"block_{index:04d}",
                pair=(block.query_record_index, block.subject_record_index),
                orientation=block.orientation,
                anchors=block.anchors,
            )
        )
    return tuple(renumbered)


def cluster_lossless_collinearity_anchors(
    anchors: Sequence[CollinearityAnchor],
    *,
    params: LosslessCollinearityParameters | None = None,
) -> CollinearityResult:
    """Cluster every input anchor, preserving singleton links."""

    resolved_params = params or LosslessCollinearityParameters()
    resolved_params.validate()
    anchors_by_pair: dict[tuple[int, int], list[CollinearityAnchor]] = {}
    for anchor in anchors:
        anchors_by_pair.setdefault(
            (int(anchor.query_record_index), int(anchor.subject_record_index)),
            [],
        ).append(anchor)

    blocks: list[CollinearityBlock] = []
    for pair in sorted(anchors_by_pair):
        blocks.extend(
            _lossless_chains_for_pair(pair, anchors_by_pair[pair], resolved_params)
        )
    filtered_blocks, unblocked_anchors = _filter_blocks_by_min_anchors(
        blocks,
        min_anchors=resolved_params.min_anchors,
    )
    return CollinearityResult(
        blocks=_renumber_lossless_blocks(filtered_blocks),
        unblocked_anchors=unblocked_anchors,
    )


def call_collinearity_blocks(
    anchors: Sequence[CollinearityAnchor],
    *,
    params: LosslessCollinearityParameters | None = None,
) -> CollinearityResult:
    """Cluster collinearity anchors with the lossless parameter contract."""

    return cluster_lossless_collinearity_anchors(anchors, params=params)


def _empty_hits() -> DataFrame:
    return pd.DataFrame(columns=COMPARISON_COLUMNS)


def _directional_hit_tables_from_sequence(
    hits_by_pair: Sequence[DataFrame],
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> dict[tuple[int, int], DataFrame]:
    tables: dict[tuple[int, int], DataFrame] = {}
    for pair_index, hits in enumerate(hits_by_pair):
        tables[(pair_index, pair_index + 1)] = hits if hits is not None else _empty_hits()
    if reverse_hits_by_pair is not None:
        for pair_index, hits in enumerate(reverse_hits_by_pair):
            tables[(pair_index + 1, pair_index)] = hits if hits is not None else _empty_hits()
    return tables


def _normalize_hit_tables(
    hits_by_pair: Sequence[DataFrame] | Mapping[tuple[int, int], DataFrame],
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> dict[tuple[int, int], DataFrame]:
    if isinstance(hits_by_pair, Mapping):
        return {
            (int(query_index), int(subject_index)): hits if hits is not None else _empty_hits()
            for (query_index, subject_index), hits in hits_by_pair.items()
        }
    return _directional_hit_tables_from_sequence(hits_by_pair, reverse_hits_by_pair)


def _filter_hit_tables_by_search_scope(
    hit_tables: Mapping[tuple[int, int], DataFrame],
    *,
    record_count: int,
    scope: CollinearitySearchScope | str,
    comparison_pairs: Sequence[tuple[int, int]] | None = None,
) -> dict[tuple[int, int], DataFrame]:
    """Keep local semantic evidence and scope only cross-record evidence."""
    normalized_scope = normalize_collinearity_search_scope(str(scope))
    if normalized_scope == "all":
        return dict(hit_tables)
    allowed_pairs = (
        {
            tuple(sorted((int(query_index), int(subject_index))))
            for query_index, subject_index in comparison_pairs
        }
        if comparison_pairs is not None
        else set(iter_collinearity_search_pairs(record_count, scope=normalized_scope))
    )
    return {
        (query_index, subject_index): hits
        for (query_index, subject_index), hits in hit_tables.items()
        if query_index == subject_index
        or (min(query_index, subject_index), max(query_index, subject_index)) in allowed_pairs
    }


def _select_orthogroup_edges(
    hit_tables: Mapping[tuple[int, int], DataFrame],
    *,
    edge_mode: CollinearityAnchorMode | str,
) -> dict[tuple[int, int], DataFrame]:
    normalized_edge_mode = normalize_collinearity_anchor_mode(str(edge_mode))
    selected: dict[tuple[int, int], DataFrame] = {}
    unordered_pairs = sorted(
        {
            (min(query_index, subject_index), max(query_index, subject_index))
            for query_index, subject_index in hit_tables
            if query_index != subject_index
        }
    )
    for query_index, subject_index in unordered_pairs:
        forward_hits = hit_tables.get((query_index, subject_index), _empty_hits())
        if normalized_edge_mode == "rbh":
            reverse_hits = hit_tables.get((subject_index, query_index), _empty_hits())
            selected[(query_index, subject_index)] = select_reciprocal_best_hit_edges(
                forward_hits,
                reverse_hits,
            )
        else:
            selected[(query_index, subject_index)] = forward_hits.copy()
    return selected


def build_orthogroup_collinearity_blocks_from_hits(
    hits_by_pair: Sequence[DataFrame] | Mapping[tuple[int, int], DataFrame],
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    params: LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    edge_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    orthogroup_member_max_hits: int | None = None,
    infer_orthogroups: bool,
    max_paralog_links_per_orthogroup: int = 2,
    comparison_pairs: Sequence[tuple[int, int]] | None = None,
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> CollinearityResult:
    """Call lossless blocks from filtered hits, optionally inferring orthogroups."""

    if not isinstance(infer_orthogroups, bool):
        raise ValidationError("infer_orthogroups must be a boolean")
    lossless_params = _resolve_lossless_params(params)
    unit_index = build_collinearity_unit_index(
        extraction,
        records=records,
        mode=unit_mode,
    )
    directional_tables = _normalize_hit_tables(hits_by_pair, reverse_hits_by_pair)
    normalized_edge_mode = normalize_collinearity_anchor_mode(str(edge_mode))
    normalized_search_scope = normalize_collinearity_search_scope(str(search_scope))
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if orthogroup_member_max_hits is not None and int(orthogroup_member_max_hits) <= 0:
        raise ValidationError("orthogroup_member_max_hits must be > 0 or None")
    resolved_max_paralog_links = int(max_paralog_links_per_orthogroup)
    if resolved_max_paralog_links <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    record_count = len(records) if records is not None else len(extraction.proteins_by_record)
    if comparison_pairs is None:
        normalized_comparison_pairs = tuple(
            (record_index, record_index + 1)
            for record_index in range(max(0, record_count - 1))
        )
    else:
        normalized_pairs: list[tuple[int, int]] = []
        for raw_pair in comparison_pairs:
            if len(raw_pair) != 2:
                raise ValidationError("comparison_pairs must contain record-index pairs.")
            query_index, subject_index = sorted((int(raw_pair[0]), int(raw_pair[1])))
            if query_index < 0 or subject_index >= record_count or query_index == subject_index:
                raise ValidationError("comparison_pairs contains an invalid record-index pair.")
            normalized_pairs.append((query_index, subject_index))
        if len(set(normalized_pairs)) != len(normalized_pairs):
            raise ValidationError("comparison_pairs must not contain duplicates.")
        normalized_comparison_pairs = tuple(normalized_pairs)
    semantic_hit_tables = _filter_hit_tables_by_search_scope(
        directional_tables,
        record_count=record_count,
        scope=normalized_search_scope,
        comparison_pairs=(
            normalized_comparison_pairs
            if normalized_search_scope == "adjacent" and comparison_pairs is not None
            else None
        ),
    )
    geometry_hit_tables = {
        pair: hits
        for pair, hits in semantic_hit_tables.items()
        if pair[0] != pair[1]
    }
    orthogroups = None
    if infer_orthogroups:
        edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
            semantic_hit_tables,
            extraction.protein_map,
            record_count=record_count,
            orthogroup_membership_mode=normalized_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_related_edges_per_orthogroup=resolved_max_paralog_links,
            comparison_pairs=(),
        )
        orthogroups = edge_selection.orthogroups
        if normalized_edge_mode == "rbh":
            edge_tables = {
                pair: edges
                for pair, edges in edge_selection.all_edges_by_pair.items()
                if pair in geometry_hit_tables or pair[::-1] in geometry_hit_tables
            }
        else:
            edge_tables = _select_orthogroup_edges(
                geometry_hit_tables, edge_mode=normalized_edge_mode,
            )
    else:
        edge_tables = _select_orthogroup_edges(
            {
                pair: _select_member_candidate_hits_per_query(hits, max_hits=orthogroup_member_max_hits)
                for pair, hits in geometry_hit_tables.items()
            },
            edge_mode=normalized_edge_mode,
        )
    comparison_edges_by_pair = {
        pair: edge_tables.get(pair, _empty_hits())
        for pair in normalized_comparison_pairs
    }
    anchors = orthogroup_edges_to_lossless_collinearity_anchors(
        comparison_edges_by_pair,
        extraction.protein_map,
        orthogroups,
        unit_index=unit_index,
    )
    result = cluster_lossless_collinearity_anchors(anchors, params=lossless_params)
    return CollinearityResult(
        blocks=result.blocks,
        unblocked_anchors=result.unblocked_anchors,
        orthogroups=orthogroups,
    )


def _validate_collinearity_extraction(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
) -> None:
    empty_record_ids = [
        str(records[index].id)
        for index, proteins in enumerate(extraction.proteins_by_record)
        if not proteins
    ]
    if empty_record_ids:
        raise ValidationError(
            "protein_blastp_mode='collinear' requires at least one CDS protein in each record; "
            "no CDS proteins were found in: "
            + ", ".join(empty_record_ids)
        )


def build_orthogroup_collinearity_blocks(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    losatp_threads: int | None = None,
    candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    orthogroup_member_max_hits: int | None = None,
    infer_orthogroups: bool,
    max_paralog_links_per_orthogroup: int = 2,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    params: LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    edge_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    comparison_pairs: Sequence[tuple[int, int]] | None = None,
    runner: LosatpRunner | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    cache_filenames: Sequence[str] | None = None,
) -> CollinearityResult:
    """Run scoped LOSATP searches and call blocks, with optional orthogroup inference."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='collinear' requires at least two records")
    if losatp_threads is not None and int(losatp_threads) <= 0:
        raise ValidationError("losatp_threads must be > 0 or None")
    if not isinstance(infer_orthogroups, bool):
        raise ValidationError("infer_orthogroups must be a boolean")
    lossless_params = _resolve_lossless_params(params)
    normalized_edge_mode = normalize_collinearity_anchor_mode(str(edge_mode))
    normalized_search_scope = normalize_collinearity_search_scope(str(search_scope))
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if orthogroup_member_max_hits is not None and int(orthogroup_member_max_hits) <= 0:
        raise ValidationError("orthogroup_member_max_hits must be > 0 or None")
    resolved_max_paralog_links = int(max_paralog_links_per_orthogroup)
    if resolved_max_paralog_links <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    extraction = protein_extraction or extract_cds_proteins(
        records,
        feature_visibility_rules=feature_visibility_rules,
    )
    _validate_collinearity_extraction(records, extraction)

    specs = losatp_job_specs(
        "collinear",
        record_count=len(records),
        pairs=(
            tuple(
                cast(tuple[int, int], tuple(sorted((int(query_index), int(subject_index)))))
                for query_index, subject_index in comparison_pairs
            )
            if comparison_pairs is not None
            else None
        ),
        infer_orthogroups=infer_orthogroups,
        search_scope=normalized_search_scope,
    )
    hits_by_pair = _search_losatp_record_pairs(
        records,
        extraction,
        specs,
        losatp_bin=losatp_bin,
        ncbi_blastp_bin=ncbi_blastp_bin,
        losatp_threads=losatp_threads,
        candidate_limit=candidate_limit,
        max_hsps_per_subject=None,
        runner=runner,
        losatp_cache=losatp_cache,
        display_filenames={
            pair: filename
            for pair, filename in _consecutive_display_filenames(
                len(records), cache_filenames
            ).items()
            if pair in specs
        },
    )
    directional_tables: dict[tuple[int, int], DataFrame] = {
        pair: filter_protein_hits_by_thresholds(
            hits_by_pair[pair],
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        for pair in specs
    }

    return build_orthogroup_collinearity_blocks_from_hits(
        directional_tables,
        extraction,
        records=records,
        params=lossless_params,
        unit_mode=unit_mode,
        edge_mode=normalized_edge_mode,
        search_scope=normalized_search_scope,
        orthogroup_membership_mode=normalized_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        infer_orthogroups=infer_orthogroups,
        max_paralog_links_per_orthogroup=resolved_max_paralog_links,
        comparison_pairs=comparison_pairs,
    )


def build_collinearity_blocks_from_hits(
    hits_by_pair: Sequence[DataFrame] | Mapping[tuple[int, int], DataFrame],
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    params: LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    anchor_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    orthogroup_member_max_hits: int | None = None,
    infer_orthogroups: bool,
    max_paralog_links_per_orthogroup: int = 2,
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> CollinearityResult:
    """Call blocks from filtered directional hits, with optional orthogroup inference."""

    return build_orthogroup_collinearity_blocks_from_hits(
        hits_by_pair,
        extraction,
        records=records,
        params=params,
        unit_mode=unit_mode,
        edge_mode=anchor_mode,
        search_scope=search_scope,
        orthogroup_membership_mode=orthogroup_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        infer_orthogroups=infer_orthogroups,
        max_paralog_links_per_orthogroup=max_paralog_links_per_orthogroup,
        reverse_hits_by_pair=reverse_hits_by_pair,
    )


def _record_ids(
    records: Sequence[SeqRecord] | None,
    record_ids: Sequence[str] | None,
) -> list[str]:
    if record_ids is not None:
        return [str(record_id) for record_id in record_ids]
    if records is not None:
        return [str(record.id) for record in records]
    return []


def _block_link_side(
    anchors: Sequence[CollinearityAnchor],
    role: Literal["query", "subject"],
    *,
    forward: bool,
    record_length: int | None,
) -> LinkSide:
    """One side of a block: the arc that covers its anchors, as directed endpoints.

    An anchor's endpoints disagree with its strand when its unit crosses the
    record origin (``gbdraw.analysis.origin_arcs``); the block then covers the
    origin too. Without the record length no anchor crosses it.
    """

    arcs = [
        directed_arc(
            int(getattr(anchor, f"{role}_start")),
            int(getattr(anchor, f"{role}_end")),
            getattr(anchor, f"{role}_strand"),
            record_length,
        )[0]
        for anchor in anchors
    ]
    # No cover: the arcs, one across the origin (so the length is known), cover the record.
    cover = cover_arcs(arcs) or (0, int(record_length or 0))
    start, end = arc_link_endpoints(cover, forward)
    return LinkSide(start, end, cover[0] > cover[1], int(record_length or 0))


def _weighted_identity(anchors: Sequence[CollinearityAnchor]) -> float:
    total_length = sum(max(0, int(anchor.alignment_length)) for anchor in anchors)
    if total_length > 0:
        return compensated_sum(
            float(anchor.identity) * max(0, int(anchor.alignment_length))
            for anchor in anchors
        ) / total_length
    return compensated_sum(float(anchor.identity) for anchor in anchors) / len(anchors)


def _joined_anchor_values(anchors: Sequence[CollinearityAnchor], attr: str) -> str:
    values: list[str] = []
    seen: set[str] = set()
    for anchor in anchors:
        value = getattr(anchor, attr)
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            values.append(text)
    return ";".join(values)


def _joined_anchor_endpoint_values(
    anchors: Sequence[CollinearityAnchor],
    role: Literal["query", "subject"],
) -> tuple[str, str, str]:
    """Return tuple-aligned stable IDs, source indexes, and view IDs."""

    values: list[tuple[str, int, str]] = []
    seen: set[tuple[str, int, str]] = set()
    for anchor in anchors:
        stable_id = str(
            getattr(anchor, f"{role}_feature_svg_id", "") or ""
        ).strip()
        view_id = str(
            getattr(anchor, f"{role}_view_feature_svg_id", "")
            or stable_id
        ).strip()
        raw_feature_index = getattr(anchor, f"{role}_feature_index", None)
        try:
            feature_index = int(cast(SupportsInt, raw_feature_index))
        except (TypeError, ValueError):
            feature_index = -1
        value = (stable_id, feature_index, view_id)
        if value in seen:
            continue
        seen.add(value)
        values.append(value)
    return (
        ";".join(value[0] for value in values),
        (
            ";".join(str(value[1]) for value in values)
            if all(value[1] >= 0 for value in values)
            else ""
        ),
        ";".join(value[2] for value in values),
    )


def convert_collinearity_blocks_to_pair_comparisons(
    result: CollinearityResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    record_ids: Sequence[str] | None = None,
    color_mode: CollinearityColorMode | str = "orientation",
    search_scope: CollinearitySearchScope | str | None = None,
) -> dict[tuple[int, int], DataFrame]:
    """Convert accepted blocks into comparison frames keyed by record endpoints.

    The native TSV keeps anchor-level details, but the diagram should display a
    collinear block as one continuous span ribbon.  The existing comparison
    renderer can draw that block ribbon when we pass it one span row per block.
    The rows name the scope of their similarity groups (``search_scope``) only
    when an anchor belongs to a group: blocks built without inferred groups
    have none to name (OV-2302).
    """

    normalized_color_mode = normalize_collinearity_color_mode(str(color_mode))
    resolved_record_ids = _record_ids(records, record_ids)
    record_count = len(resolved_record_ids)
    if records is not None:
        record_count = len(records)
    else:
        for block in result.blocks:
            record_count = max(record_count, int(block.query_record_index) + 1, int(block.subject_record_index) + 1)
    if record_count < 2:
        return {}

    record_lengths = [len(record) for record in records] if records is not None else None
    rows_by_pair: dict[tuple[int, int], list[dict[str, object]]] = {}
    for block in result.blocks:
        pair = (int(block.query_record_index), int(block.subject_record_index))
        if pair[0] < 0 or pair[1] >= record_count or pair[0] >= pair[1]:
            raise ValidationError("Collinearity block contains invalid record endpoints.")
        if not block.anchors:
            continue
        first_anchor = block.anchors[0]
        query_name = (
            resolved_record_ids[first_anchor.query_record_index]
            if first_anchor.query_record_index < len(resolved_record_ids)
            else str(first_anchor.query_record_index + 1)
        )
        subject_name = (
            resolved_record_ids[first_anchor.subject_record_index]
            if first_anchor.subject_record_index < len(resolved_record_ids)
            else str(first_anchor.subject_record_index + 1)
        )
        query_side = _block_link_side(
            block.anchors, "query", forward=True,
            record_length=record_lengths[pair[0]] if record_lengths is not None else None,
        )
        subject_side = _block_link_side(
            block.anchors, "subject", forward=block.orientation != "minus",
            record_length=record_lengths[pair[1]] if record_lengths is not None else None,
        )
        (
            query_feature_svg_ids,
            query_feature_indexes,
            query_view_feature_svg_ids,
        ) = _joined_anchor_endpoint_values(block.anchors, "query")
        (
            subject_feature_svg_ids,
            subject_feature_indexes,
            subject_view_feature_svg_ids,
        ) = _joined_anchor_endpoint_values(block.anchors, "subject")
        rows_by_pair.setdefault(pair, []).extend(rows_inside_records(
            {
                "query": query_name,
                "subject": subject_name,
                "identity": float(_weighted_identity(block.anchors)),
                "alignment_length": sum(max(0, int(anchor.alignment_length)) for anchor in block.anchors),
                "mismatches": 0,
                "gap_opens": 0,
                "qstart": query_side.start,
                "qend": query_side.end,
                "sstart": subject_side.start,
                "send": subject_side.end,
                "evalue": min(float(anchor.evalue) for anchor in block.anchors),
                "bitscore": _anchor_bitscore_total(block.anchors),
                "collinearity_block_id": block.block_id,
                "collinearity_block_kind": block.kind,
                "collinearity_orientation": block.orientation,
                "collinearity_block_score": float(block.score),
                "collinearity_anchor_index": 1,
                "collinearity_anchor_count": len(block.anchors),
                "collinearity_color_mode": normalized_color_mode,
                "orthogroup_id": _joined_anchor_values(block.anchors, "orthogroup_id"),
                "query_unit_id": _joined_anchor_values(block.anchors, "query_unit_id"),
                "subject_unit_id": _joined_anchor_values(block.anchors, "subject_unit_id"),
                "query_unit_kind": first_anchor.query_unit_kind,
                "subject_unit_kind": first_anchor.subject_unit_kind,
                "query_locus_id": _joined_anchor_values(block.anchors, "query_locus_id"),
                "subject_locus_id": _joined_anchor_values(block.anchors, "subject_locus_id"),
                "query_display_name": _joined_anchor_values(block.anchors, "query_display_name"),
                "subject_display_name": _joined_anchor_values(block.anchors, "subject_display_name"),
                "query_protein_id": _joined_anchor_values(block.anchors, "query_protein_id"),
                "subject_protein_id": _joined_anchor_values(block.anchors, "subject_protein_id"),
                "query_feature_svg_id": query_feature_svg_ids,
                "subject_feature_svg_id": subject_feature_svg_ids,
                "query_feature_index": query_feature_indexes,
                "subject_feature_index": subject_feature_indexes,
                "query_view_feature_svg_id": query_view_feature_svg_ids,
                "subject_view_feature_svg_id": subject_view_feature_svg_ids,
                "rbh_orthogroup_id": _joined_anchor_values(block.anchors, "rbh_orthogroup_id"),
                "ortholog_path_id": _joined_anchor_values(block.anchors, "ortholog_path_id"),
                "edge_kind": _joined_anchor_values(block.anchors, "edge_kind") or "rbh",
                "render_role": _joined_anchor_values(block.anchors, "render_role") or "block_anchor",
                "query_orthogroup_representative": _joined_anchor_values(block.anchors, "query_orthogroup_representative"),
                "subject_orthogroup_representative": _joined_anchor_values(block.anchors, "subject_orthogroup_representative"),
                "query_orthogroup_member_count": _joined_anchor_values(block.anchors, "query_orthogroup_member_count"),
                "subject_orthogroup_member_count": _joined_anchor_values(block.anchors, "subject_orthogroup_member_count"),
            },
            query_side,
            subject_side,
        ))

    comparisons = {
        pair: pd.DataFrame.from_records(rows, columns=COLLINEARITY_COMPARISON_COLUMNS)
        for pair, rows in sorted(rows_by_pair.items())
    }
    has_groups = any(anchor.orthogroup_id for block in result.blocks for anchor in block.anchors)
    if search_scope is not None and has_groups:
        normalized_search_scope = normalize_collinearity_search_scope(
            str(search_scope)
        )
        presentation_scope = (
            "global_collinear"
            if normalized_search_scope == "all"
            else "adjacent_local"
        )
        group_kind = (
            "orthogroup"
            if presentation_scope == "global_collinear"
            else "collinear_gene_group"
        )
        for comparison in comparisons.values():
            comparison["group_kind"] = group_kind
            comparison["group_scope"] = presentation_scope
            comparison["collinear_group_scope"] = presentation_scope
    return comparisons


def convert_collinearity_blocks_to_comparisons(
    result: CollinearityResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    record_ids: Sequence[str] | None = None,
    color_mode: CollinearityColorMode | str = "orientation",
    search_scope: CollinearitySearchScope | str | None = None,
) -> list[DataFrame]:
    """Convert accepted blocks into legacy adjacent-record comparison frames."""

    pair_comparisons = convert_collinearity_blocks_to_pair_comparisons(
        result,
        records=records,
        record_ids=record_ids,
        color_mode=color_mode,
        search_scope=search_scope,
    )
    record_count = len(records) if records is not None else len(_record_ids(records, record_ids))
    if records is None:
        for query_index, subject_index in pair_comparisons:
            record_count = max(record_count, query_index + 1, subject_index + 1)
    non_adjacent = [pair for pair in pair_comparisons if pair[1] != pair[0] + 1]
    if non_adjacent:
        raise ValidationError("Collinearity rendering supports adjacent record pairs only.")
    return [
        pair_comparisons.get(
            (record_index, record_index + 1),
            pd.DataFrame(columns=COLLINEARITY_COMPARISON_COLUMNS),
        )
        for record_index in range(max(0, record_count - 1))
    ]


__all__ = [
    "COLLINEARITY_ANCHOR_MODES",
    "COLLINEARITY_COLOR_MODES",
    "COLLINEARITY_COMPARISON_COLUMNS",
    "COLLINEARITY_MERGE_ORIENTATIONS",
    "COLLINEARITY_METADATA_COLUMNS",
    "COLLINEARITY_SEARCH_SCOPES",
    "CollinearityAnchor",
    "CollinearityAnchorMode",
    "CollinearityBlock",
    "CollinearityBlockKind",
    "CollinearityColorMode",
    "CollinearityMergeOrientation",
    "CollinearityResult",
    "CollinearitySearchScope",
    "LosslessCollinearityParameters",
    "build_collinearity_blocks_from_hits",
    "build_orthogroup_collinearity_blocks",
    "build_orthogroup_collinearity_blocks_from_hits",
    "call_collinearity_blocks",
    "cluster_lossless_collinearity_anchors",
    "convert_collinearity_blocks_to_comparisons",
    "convert_collinearity_blocks_to_pair_comparisons",
    "iter_collinearity_search_pairs",
    "normalize_collinearity_anchor_mode",
    "normalize_collinearity_color_mode",
    "normalize_collinearity_search_scope",
    "orthogroup_edges_to_lossless_collinearity_anchors",
]
