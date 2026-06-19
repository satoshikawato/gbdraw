#!/usr/bin/env python
# coding: utf-8

"""Native gbdraw collinear block calling."""

from __future__ import annotations

from dataclasses import dataclass
import math
import sys
from typing import Literal, Mapping, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity_units import (
    CollinearityUnit,
    CollinearityUnitIndex,
    CollinearityUnitMode,
)
from gbdraw.analysis.protein_colinearity import (
    CdsProtein,
    LosatpRunner,
    OrthogroupMember,
    OrthogroupMembershipMode,
    OrthogroupResult,
    ProteinExtractionResult,
    _run_losatp_search,
    build_orthogroups_from_protein_hits,
    extract_cds_proteins,
    filter_protein_hits_by_thresholds,
    normalize_orthogroup_membership_mode,
    proteins_to_fasta,
    select_reciprocal_best_hit_edges,
    select_reciprocal_best_hits,
    select_rbh_orthogroup_edges_from_directional_hits,
)
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS

CollinearityOrientation = Literal["plus", "minus"]
CollinearityBlockKind = Literal["syntenic", "cluster", "singleton"]
CollinearityScoreMode = Literal["constant", "bitscore"]
CollinearityColorMode = Literal["average_identity", "orientation", "orientation_identity"]
CollinearityAnchorMode = Literal["all", "one_to_one", "rbh"]
CollinearitySearchScope = Literal["adjacent", "all"]

COLLINEARITY_METADATA_COLUMNS = (
    "collinearity_block_id",
    "collinearity_block_kind",
    "collinearity_orientation",
    "collinearity_block_score",
    "collinearity_block_evalue",
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
COLLINEARITY_COLOR_MODES = ("average_identity", "orientation", "orientation_identity")
COLLINEARITY_ANCHOR_MODES = ("all", "one_to_one", "rbh")
COLLINEARITY_SEARCH_SCOPES = ("adjacent", "all")


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


@dataclass(frozen=True)
class CollinearityBlock:
    block_id: str
    query_record_index: int
    subject_record_index: int
    orientation: CollinearityOrientation
    score: float
    anchors: tuple[CollinearityAnchor, ...]
    kind: CollinearityBlockKind = "syntenic"
    block_evalue: float | None = None
    label: str | None = None


@dataclass(frozen=True)
class CollinearityResult:
    blocks: tuple[CollinearityBlock, ...]
    unblocked_anchors: tuple[CollinearityAnchor, ...] = ()
    orthogroups: OrthogroupResult | None = None


@dataclass(frozen=True)
class CollinearityParameters:
    min_anchors: int = 1
    max_gene_gap: int = 25
    block_merge_gap: int = 50
    singleton_merge_gap: int = 25
    max_diagonal_drift: int = 25
    max_conflicts_in_merge_gap: int = 1
    max_paralog_links_per_orthogroup: int = 2
    gap_penalty: float = 1.0
    nearby_duplicate_window: int = 0
    score_mode: CollinearityScoreMode = "constant"
    constant_anchor_score: float = 50.0
    min_block_score: float | None = None
    block_evalue: float | None = None

    def validate(self) -> None:
        if int(self.min_anchors) <= 0:
            raise ValidationError("collinear_min_anchors must be > 0")
        if int(self.max_gene_gap) < 0:
            raise ValidationError("collinear_max_gene_gap must be >= 0")
        if int(self.block_merge_gap) < 0:
            raise ValidationError("collinear_block_merge_gap must be >= 0")
        if int(self.singleton_merge_gap) < 0:
            raise ValidationError("collinear_singleton_merge_gap must be >= 0")
        if int(self.max_diagonal_drift) < 0:
            raise ValidationError("collinear_max_diagonal_drift must be >= 0")
        if int(self.max_conflicts_in_merge_gap) < 0:
            raise ValidationError("collinear_max_conflicts_in_merge_gap must be >= 0")
        if int(self.max_paralog_links_per_orthogroup) <= 0:
            raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
        if float(self.gap_penalty) < 0:
            raise ValidationError("collinear_gap_penalty must be >= 0")
        if int(self.nearby_duplicate_window) < 0:
            raise ValidationError("collinear_nearby_duplicate_window must be >= 0")
        if str(self.score_mode) not in {"constant", "bitscore"}:
            raise ValidationError("collinear_score_mode must be one of: constant, bitscore")
        if float(self.constant_anchor_score) <= 0:
            raise ValidationError("collinear_constant_anchor_score must be > 0")
        if self.block_evalue is not None:
            block_evalue = float(self.block_evalue)
            if not math.isfinite(block_evalue) or block_evalue < 0:
                raise ValidationError("collinear_block_evalue must be a finite value >= 0 or None")

    def effective_min_block_score(self) -> float:
        if self.min_block_score is not None:
            return float(self.min_block_score)
        if self.score_mode == "constant":
            return float(self.constant_anchor_score) * int(self.min_anchors)
        return 0.0


@dataclass(frozen=True)
class LosslessCollinearityParameters:
    min_anchors: int = 1
    max_unit_gap: int = 0
    max_diagonal_drift: int = 0
    max_conflicts: int = 0
    merge_orientation: Literal["strand", "order", "either"] = "either"

    def validate(self) -> None:
        if int(self.min_anchors) <= 0:
            raise ValidationError("collinear_min_anchors must be > 0")
        if int(self.max_unit_gap) < 0:
            raise ValidationError("collinear_max_unit_gap must be >= 0")
        if int(self.max_diagonal_drift) < 0:
            raise ValidationError("collinear_max_diagonal_drift must be >= 0")
        if int(self.max_conflicts) < 0:
            raise ValidationError("collinear_max_conflicts must be >= 0")
        if str(self.merge_orientation) not in {"strand", "order", "either"}:
            raise ValidationError("collinear_merge_orientation must be one of: strand, order, either")


@dataclass(frozen=True)
class _Chain:
    orientation: CollinearityOrientation
    score: float
    anchors: tuple[CollinearityAnchor, ...]


def normalize_collinearity_color_mode(mode: str | None) -> CollinearityColorMode:
    normalized = str(mode or "orientation").strip().lower().replace("-", "_")
    if normalized == "identity":
        normalized = "average_identity"
    if normalized not in COLLINEARITY_COLOR_MODES:
        raise ValidationError(
            "collinear_color_mode must be one of: "
            + ", ".join(COLLINEARITY_COLOR_MODES)
        )
    return normalized  # type: ignore[return-value]


def normalize_collinearity_anchor_mode(mode: str | None) -> CollinearityAnchorMode:
    normalized = str(mode or "rbh").strip().lower().replace("-", "_")
    aliases = {
        "all_hits": "all",
        "raw": "all",
        "top_n": "all",
        "topn": "all",
        "one2one": "one_to_one",
        "one_to_ones": "one_to_one",
        "mutual_best": "one_to_one",
        "top1": "one_to_one",
        "top_1": "one_to_one",
        "reciprocal_best": "rbh",
        "strict_rbh": "rbh",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in COLLINEARITY_ANCHOR_MODES:
        raise ValidationError(
            "collinear_anchor_mode must be one of: "
            + ", ".join(COLLINEARITY_ANCHOR_MODES)
        )
    return normalized  # type: ignore[return-value]


def normalize_collinearity_search_scope(scope: str | None) -> CollinearitySearchScope:
    normalized = str(scope or "adjacent").strip().lower().replace("-", "_")
    if normalized not in COLLINEARITY_SEARCH_SCOPES:
        raise ValidationError(
            "collinear_search_scope must be one of: "
            + ", ".join(COLLINEARITY_SEARCH_SCOPES)
        )
    return normalized  # type: ignore[return-value]


def iter_collinearity_search_pairs(
    record_count: int,
    *,
    scope: CollinearitySearchScope | str,
) -> tuple[tuple[int, int], ...]:
    """Return unordered record pairs searched for Collinear evidence."""

    normalized_scope = normalize_collinearity_search_scope(str(scope))
    count = max(0, int(record_count))
    if normalized_scope == "adjacent":
        return tuple((index, index + 1) for index in range(max(0, count - 1)))
    return tuple(
        (query_index, subject_index)
        for query_index in range(count)
        for subject_index in range(query_index + 1, count)
    )


def select_collinearity_anchor_hits(
    hits: DataFrame,
    *,
    anchor_mode: CollinearityAnchorMode | str = "rbh",
    reverse_hits: DataFrame | None = None,
) -> DataFrame:
    """Reduce filtered protein hits to the anchors used for block calling."""

    normalized_anchor_mode = normalize_collinearity_anchor_mode(str(anchor_mode))
    if normalized_anchor_mode == "all":
        return hits.copy()
    if normalized_anchor_mode == "one_to_one":
        return select_reciprocal_best_hits(hits)
    if reverse_hits is None:
        raise ValidationError("collinear_anchor_mode='rbh' requires reverse-direction hits.")
    return select_reciprocal_best_hit_edges(hits, reverse_hits)


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


def _anchor_strength_key(anchor: CollinearityAnchor) -> tuple[float, float, float, int, str, str]:
    return (
        float(anchor.evalue),
        -float(anchor.bitscore),
        -float(anchor.identity),
        -int(anchor.alignment_length),
        str(anchor.query_protein_id),
        str(anchor.subject_protein_id),
    )


def _anchor_score(anchor: CollinearityAnchor, params: CollinearityParameters) -> float:
    if params.score_mode == "bitscore":
        return float(anchor.bitscore)
    return float(params.constant_anchor_score)


def _genomic_link_coordinates(protein: CdsProtein) -> tuple[int, int]:
    if protein.strand == -1:
        return int(protein.end), int(protein.start) + 1
    return int(protein.start) + 1, int(protein.end)


def _normalized_strand(value: object) -> int | None:
    try:
        strand = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return strand if strand in {-1, 1} else None


def _anchor_strand_orientation(anchor: CollinearityAnchor) -> CollinearityOrientation | None:
    query_strand = _normalized_strand(anchor.query_strand)
    subject_strand = _normalized_strand(anchor.subject_strand)
    if query_strand is None or subject_strand is None:
        return None
    return "plus" if query_strand == subject_strand else "minus"


def _chain_strand_support(chain: _Chain) -> tuple[int, int]:
    supported = 0
    conflicted = 0
    for anchor in chain.anchors:
        strand_orientation = _anchor_strand_orientation(anchor)
        if strand_orientation is None:
            continue
        if strand_orientation == chain.orientation:
            supported += 1
        else:
            conflicted += 1
    return supported, conflicted


def _representative_proteins(
    units: Sequence[CollinearityUnit],
    protein_map: Mapping[str, CdsProtein],
) -> list[CdsProtein]:
    proteins: list[CdsProtein] = []
    missing_ids: list[str] = []
    for unit in units:
        protein = protein_map.get(unit.representative_protein_id)
        if protein is None:
            missing_ids.append(unit.representative_protein_id)
        else:
            proteins.append(protein)
    if missing_ids:
        raise ValidationError(
            "collinearity unit representatives are missing from the protein map: "
            + ", ".join(sorted(missing_ids)[:10])
        )
    return proteins


def protein_hits_to_collinearity_anchors(
    hits: DataFrame,
    *,
    query_units: Mapping[str, CollinearityUnit],
    subject_units: Mapping[str, CollinearityUnit],
    query_protein_map: Mapping[str, CdsProtein],
    subject_protein_map: Mapping[str, CdsProtein],
    source: str = "losatp",
) -> tuple[CollinearityAnchor, ...]:
    """Map retained protein hits to collinearity unit anchors."""

    if hits.empty:
        return ()
    missing_ids: set[str] = set()
    anchors: list[CollinearityAnchor] = []

    for row in hits.itertuples(index=False):
        query_id = str(row.query)
        subject_id = str(row.subject)
        query_unit = query_units.get(query_id)
        subject_unit = subject_units.get(subject_id)
        query_protein = query_protein_map.get(query_unit.representative_protein_id) if query_unit else None
        subject_protein = subject_protein_map.get(subject_unit.representative_protein_id) if subject_unit else None
        if query_unit is None:
            missing_ids.add(query_id)
        if subject_unit is None:
            missing_ids.add(subject_id)
        if query_unit is None or subject_unit is None or query_protein is None or subject_protein is None:
            continue
        if query_unit.unit_id == subject_unit.unit_id:
            continue
        if query_unit.record_index == subject_unit.record_index:
            continue
        qstart, qend = _genomic_link_coordinates(query_protein)
        sstart, send = _genomic_link_coordinates(subject_protein)
        anchors.append(
            CollinearityAnchor(
                query_protein_id=query_protein.protein_id,
                subject_protein_id=subject_protein.protein_id,
                query_record_index=int(query_unit.record_index),
                subject_record_index=int(subject_unit.record_index),
                query_order=int(query_unit.order),
                subject_order=int(subject_unit.order),
                query_start=qstart,
                query_end=qend,
                subject_start=sstart,
                subject_end=send,
                query_strand=_normalized_strand(query_protein.strand),
                subject_strand=_normalized_strand(subject_protein.strand),
                identity=_float_from_row(row, "identity"),
                evalue=_float_from_row(row, "evalue"),
                bitscore=_float_from_row(row, "bitscore"),
                alignment_length=_int_from_row(row, "alignment_length"),
                query_feature_svg_id=query_unit.representative_feature_svg_id,
                subject_feature_svg_id=subject_unit.representative_feature_svg_id,
                source=source,
                query_unit_id=query_unit.unit_id,
                subject_unit_id=subject_unit.unit_id,
                query_unit_kind=query_unit.unit_kind,
                subject_unit_kind=subject_unit.unit_kind,
                query_locus_id=query_unit.locus_id,
                subject_locus_id=subject_unit.locus_id,
                query_display_name=query_unit.display_name,
                subject_display_name=subject_unit.display_name,
                orthogroup_id=str(getattr(row, "orthogroup_id", "") or ""),
            )
        )

    if missing_ids:
        raise ParseError(
            "LOSATP blastp output contains protein IDs that do not resolve to collinearity units: "
            + ", ".join(sorted(missing_ids)[:20])
        )
    return tuple(anchors)


def deduplicate_unit_pair_anchors(
    anchors: Sequence[CollinearityAnchor],
) -> tuple[CollinearityAnchor, ...]:
    """Keep the strongest anchor for each query-subject unit pair."""

    best_by_pair: dict[tuple[str, str], CollinearityAnchor] = {}
    for anchor in anchors:
        key = (anchor.query_unit_id, anchor.subject_unit_id)
        current = best_by_pair.get(key)
        if current is None or _anchor_strength_key(anchor) < _anchor_strength_key(current):
            best_by_pair[key] = anchor
    return tuple(
        sorted(
            best_by_pair.values(),
            key=lambda anchor: (
                anchor.query_record_index,
                anchor.subject_record_index,
                anchor.query_order,
                anchor.subject_order,
                anchor.query_protein_id,
                anchor.subject_protein_id,
            ),
        )
    )


def _collapse_duplicate_pass(
    anchors: Sequence[CollinearityAnchor],
    *,
    fixed_attr: str,
    nearby_attr: str,
    window: int,
) -> tuple[CollinearityAnchor, ...]:
    if window <= 0 or len(anchors) <= 1:
        return tuple(anchors)
    sorted_anchors = sorted(
        anchors,
        key=lambda anchor: (
            getattr(anchor, fixed_attr),
            getattr(anchor, nearby_attr),
            anchor.query_order,
            anchor.subject_order,
            anchor.query_protein_id,
            anchor.subject_protein_id,
        ),
    )
    kept: list[CollinearityAnchor] = []
    cluster: list[CollinearityAnchor] = []
    cluster_fixed: int | None = None
    cluster_max_nearby: int | None = None

    def flush() -> None:
        nonlocal cluster
        if cluster:
            kept.append(min(cluster, key=_anchor_strength_key))
            cluster = []

    for anchor in sorted_anchors:
        fixed = int(getattr(anchor, fixed_attr))
        nearby = int(getattr(anchor, nearby_attr))
        starts_new = (
            cluster_fixed is None
            or fixed != cluster_fixed
            or cluster_max_nearby is None
            or nearby - cluster_max_nearby > int(window)
        )
        if starts_new:
            flush()
            cluster_fixed = fixed
            cluster_max_nearby = nearby
        else:
            cluster_max_nearby = max(cluster_max_nearby, nearby)
        cluster.append(anchor)
    flush()
    return tuple(kept)


def collapse_nearby_duplicate_anchors(
    anchors: Sequence[CollinearityAnchor],
    *,
    window: int,
) -> tuple[CollinearityAnchor, ...]:
    """Collapse tandem/repetitive nearby unit anchors in both dimensions."""

    if int(window) <= 0:
        return tuple(anchors)
    query_collapsed = _collapse_duplicate_pass(
        anchors,
        fixed_attr="query_order",
        nearby_attr="subject_order",
        window=int(window),
    )
    subject_collapsed = _collapse_duplicate_pass(
        query_collapsed,
        fixed_attr="subject_order",
        nearby_attr="query_order",
        window=int(window),
    )
    return tuple(
        sorted(
            subject_collapsed,
            key=lambda anchor: (
                anchor.query_record_index,
                anchor.subject_record_index,
                anchor.query_order,
                anchor.subject_order,
                anchor.query_protein_id,
                anchor.subject_protein_id,
            ),
        )
    )


def _anchor_identity(anchor: CollinearityAnchor) -> tuple[str, str, int, int, str]:
    return (
        str(anchor.query_unit_id),
        str(anchor.subject_unit_id),
        int(anchor.query_record_index),
        int(anchor.subject_record_index),
        str(anchor.orthogroup_id),
    )


def _path_sorted_anchors(
    anchors: Sequence[CollinearityAnchor],
    orientation: CollinearityOrientation,
) -> tuple[CollinearityAnchor, ...]:
    if orientation == "minus":
        return tuple(
            sorted(
                anchors,
                key=lambda anchor: (
                    int(anchor.query_order),
                    -int(anchor.subject_order),
                    str(anchor.orthogroup_id),
                    str(anchor.query_protein_id),
                    str(anchor.subject_protein_id),
                ),
            )
        )
    return tuple(
        sorted(
            anchors,
            key=lambda anchor: (
                int(anchor.query_order),
                int(anchor.subject_order),
                str(anchor.orthogroup_id),
                str(anchor.query_protein_id),
                str(anchor.subject_protein_id),
            ),
        )
    )


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


def _chain_density(chain: _Chain) -> float:
    anchors = _path_sorted_anchors(chain.anchors, chain.orientation)
    if not anchors:
        return 0.0
    query_span = max(anchor.query_order for anchor in anchors) - min(anchor.query_order for anchor in anchors) + 1
    subject_span = max(anchor.subject_order for anchor in anchors) - min(anchor.subject_order for anchor in anchors) + 1
    return len(anchors) / max(1, max(int(query_span), int(subject_span)))


def _chain_diagonal_drift(chain: _Chain) -> int:
    values = [
        _orientation_diagonal(anchor, chain.orientation)
        for anchor in chain.anchors
    ]
    if not values:
        return 0
    return max(values) - min(values)


def _orthogroup_score(anchors: Sequence[CollinearityAnchor], params: CollinearityParameters) -> float:
    return float(sum(_anchor_score(anchor, params) for anchor in anchors))


def _candidate_run_rank(chain: _Chain) -> tuple[float, float, int, float, int, int, str, str, int]:
    anchors = _path_sorted_anchors(chain.anchors, chain.orientation)
    first = anchors[0]
    orientation_rank = 0 if chain.orientation == "plus" else 1
    return (
        -len(anchors),
        -_chain_density(chain),
        _chain_diagonal_drift(chain),
        -float(chain.score),
        int(first.query_order),
        _orientation_subject_value(first, chain.orientation),
        str(first.orthogroup_id),
        str(first.query_protein_id),
        orientation_rank,
    )


def _candidate_runs_for_orientation(
    anchors: Sequence[CollinearityAnchor],
    *,
    orientation: CollinearityOrientation,
    params: CollinearityParameters,
) -> list[_Chain]:
    ordered = sorted(
        anchors,
        key=lambda anchor: (
            int(anchor.query_order),
            _orientation_subject_value(anchor, orientation),
            str(anchor.orthogroup_id),
            str(anchor.query_protein_id),
            str(anchor.subject_protein_id),
        ),
    )
    runs: list[_Chain] = []
    seen_paths: set[tuple[tuple[str, str, int, int, str], ...]] = set()
    for start_index, start_anchor in enumerate(ordered):
        run = [start_anchor]
        last = start_anchor
        for candidate in ordered[start_index + 1 :]:
            if _anchors_are_path_compatible(
                last,
                candidate,
                orientation=orientation,
                max_gap=int(params.max_gene_gap),
                max_diagonal_drift=int(params.max_diagonal_drift),
            ):
                run.append(candidate)
                last = candidate
        path = _path_sorted_anchors(run, orientation)
        path_key = tuple(_anchor_identity(anchor) for anchor in path)
        if path_key in seen_paths:
            continue
        seen_paths.add(path_key)
        runs.append(
            _Chain(
                orientation=orientation,
                score=_orthogroup_score(path, params),
                anchors=path,
            )
        )
    return runs


def _resolve_overlapping_runs(
    runs: Sequence[_Chain],
    *,
    params: CollinearityParameters,
) -> tuple[_Chain, ...]:
    minimum_syntenic_anchors = max(2, int(params.min_anchors))
    accepted: list[_Chain] = []
    consumed: set[tuple[str, str, int, int, str]] = set()
    for run in sorted(runs, key=_candidate_run_rank):
        if len(run.anchors) < minimum_syntenic_anchors:
            continue
        run_ids = {_anchor_identity(anchor) for anchor in run.anchors}
        if run_ids & consumed:
            continue
        accepted.append(run)
        consumed.update(run_ids)
    return tuple(accepted)


def _block_from_anchors(
    *,
    block_id: str,
    pair: tuple[int, int],
    orientation: CollinearityOrientation,
    kind: CollinearityBlockKind,
    anchors: Sequence[CollinearityAnchor],
    params: CollinearityParameters,
) -> CollinearityBlock:
    path = _path_sorted_anchors(anchors, orientation)
    return CollinearityBlock(
        block_id=block_id,
        query_record_index=int(pair[0]),
        subject_record_index=int(pair[1]),
        orientation=orientation,
        kind=kind,
        score=_orthogroup_score(path, params),
        anchors=path,
        block_evalue=None,
    )


def _conflicts_between_blocks(
    left: CollinearityBlock,
    right: CollinearityBlock,
    anchors: Sequence[CollinearityAnchor],
) -> int:
    left_path = _path_sorted_anchors(left.anchors, left.orientation)
    right_path = _path_sorted_anchors(right.anchors, right.orientation)
    if not left_path or not right_path:
        return 0
    left_end = left_path[-1]
    right_start = right_path[0]
    query_min = min(int(left_end.query_order), int(right_start.query_order))
    query_max = max(int(left_end.query_order), int(right_start.query_order))
    subject_min = min(int(left_end.subject_order), int(right_start.subject_order))
    subject_max = max(int(left_end.subject_order), int(right_start.subject_order))
    block_ids = {_anchor_identity(anchor) for anchor in (*left.anchors, *right.anchors)}
    return sum(
        1
        for anchor in anchors
        if _anchor_identity(anchor) not in block_ids
        and query_min < int(anchor.query_order) < query_max
        and subject_min < int(anchor.subject_order) < subject_max
    )


def _blocks_can_merge(
    left: CollinearityBlock,
    right: CollinearityBlock,
    *,
    anchors: Sequence[CollinearityAnchor],
    params: CollinearityParameters,
) -> bool:
    if left.kind != "syntenic" or right.kind != "syntenic":
        return False
    if left.orientation != right.orientation:
        return False
    if int(left.query_record_index) != int(right.query_record_index):
        return False
    if int(left.subject_record_index) != int(right.subject_record_index):
        return False
    left_path = _path_sorted_anchors(left.anchors, left.orientation)
    right_path = _path_sorted_anchors(right.anchors, right.orientation)
    if not left_path or not right_path:
        return False
    if not _anchors_are_path_compatible(
        left_path[-1],
        right_path[0],
        orientation=left.orientation,
        max_gap=int(params.block_merge_gap),
        max_diagonal_drift=int(params.max_diagonal_drift),
    ):
        return False
    conflict_count = _conflicts_between_blocks(left, right, anchors)
    return conflict_count <= int(params.max_conflicts_in_merge_gap)


def _merge_syntenic_blocks(
    blocks: Sequence[CollinearityBlock],
    *,
    anchors: Sequence[CollinearityAnchor],
    params: CollinearityParameters,
) -> tuple[CollinearityBlock, ...]:
    merged: list[CollinearityBlock] = []
    for block in sorted(blocks, key=_final_block_sort_key):
        if not merged:
            merged.append(block)
            continue
        previous = merged[-1]
        if _blocks_can_merge(previous, block, anchors=anchors, params=params):
            merged[-1] = _block_from_anchors(
                block_id=previous.block_id,
                pair=(previous.query_record_index, previous.subject_record_index),
                orientation=previous.orientation,
                kind="syntenic",
                anchors=(*previous.anchors, *block.anchors),
                params=params,
            )
        else:
            merged.append(block)
    return tuple(merged)


def _anchor_singleton_orientation(anchor: CollinearityAnchor) -> CollinearityOrientation:
    return _anchor_strand_orientation(anchor) or "plus"


def _anchor_can_absorb_into_block(
    anchor: CollinearityAnchor,
    block: CollinearityBlock,
    *,
    params: CollinearityParameters,
) -> bool:
    if block.kind != "syntenic":
        return False
    if int(anchor.query_record_index) != int(block.query_record_index):
        return False
    if int(anchor.subject_record_index) != int(block.subject_record_index):
        return False
    anchor_orientation = _anchor_strand_orientation(anchor)
    if anchor_orientation is not None and anchor_orientation != block.orientation:
        return False
    path = list(_path_sorted_anchors((*block.anchors, anchor), block.orientation))
    try:
        index = path.index(anchor)
    except ValueError:
        return False
    if index > 0 and not _anchors_are_path_compatible(
        path[index - 1],
        anchor,
        orientation=block.orientation,
        max_gap=int(params.singleton_merge_gap),
        max_diagonal_drift=int(params.max_diagonal_drift),
    ):
        return False
    if index + 1 < len(path) and not _anchors_are_path_compatible(
        anchor,
        path[index + 1],
        orientation=block.orientation,
        max_gap=int(params.singleton_merge_gap),
        max_diagonal_drift=int(params.max_diagonal_drift),
    ):
        return False
    return True


def _absorption_rank(anchor: CollinearityAnchor, block: CollinearityBlock) -> tuple[int, int, int, str]:
    distances = [
        max(
            abs(int(anchor.query_order) - int(block_anchor.query_order)),
            abs(int(anchor.subject_order) - int(block_anchor.subject_order)),
        )
        for block_anchor in block.anchors
    ]
    diagonal_distances = [
        abs(
            _orientation_diagonal(anchor, block.orientation)
            - _orientation_diagonal(block_anchor, block.orientation)
        )
        for block_anchor in block.anchors
    ]
    return (
        min(distances) if distances else 0,
        min(diagonal_distances) if diagonal_distances else 0,
        int(block.query_record_index),
        str(block.block_id),
    )


def _absorb_singleton_anchors(
    blocks: Sequence[CollinearityBlock],
    singletons: Sequence[CollinearityAnchor],
    *,
    params: CollinearityParameters,
) -> tuple[tuple[CollinearityBlock, ...], tuple[CollinearityAnchor, ...]]:
    updated = list(blocks)
    remaining: list[CollinearityAnchor] = []
    for anchor in sorted(singletons, key=_anchor_sort_key):
        candidates = [
            (index, block)
            for index, block in enumerate(updated)
            if _anchor_can_absorb_into_block(anchor, block, params=params)
        ]
        if not candidates:
            remaining.append(anchor)
            continue
        best_index, best_block = min(candidates, key=lambda item: _absorption_rank(anchor, item[1]))
        updated[best_index] = _block_from_anchors(
            block_id=best_block.block_id,
            pair=(best_block.query_record_index, best_block.subject_record_index),
            orientation=best_block.orientation,
            kind="syntenic",
            anchors=(*best_block.anchors, anchor),
            params=params,
        )
    return tuple(updated), tuple(remaining)


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
    orthogroups: OrthogroupResult,
) -> dict[tuple[str, int], int]:
    counts: dict[tuple[str, int], int] = {}
    for orthogroup_id, members in orthogroups.orthogroups.items():
        for member in members:
            key = (str(orthogroup_id), int(member.record_index))
            counts[key] = counts.get(key, 0) + 1
    return counts


def _orthogroup_id_for_edge(
    query_id: str,
    subject_id: str,
    orthogroups: OrthogroupResult,
) -> str:
    query_member = orthogroups.member_by_protein_id.get(query_id)
    subject_member = orthogroups.member_by_protein_id.get(subject_id)
    if query_member is None or subject_member is None:
        return ""
    if query_member.orthogroup_id != subject_member.orthogroup_id:
        return ""
    return str(query_member.orthogroup_id)


def _orthogroup_edge_metadata_for_anchor(
    query_id: str,
    subject_id: str,
    orthogroup_id: str,
    orthogroups: OrthogroupResult,
) -> dict[str, object]:
    metadata = {
        "rbh_orthogroup_id": "",
        "ortholog_path_id": "",
        "edge_kind": "",
        "render_role": "",
        "query_orthogroup_representative": False,
        "subject_orthogroup_representative": False,
    }
    query_member = orthogroups.member_by_protein_id.get(query_id)
    subject_member = orthogroups.member_by_protein_id.get(subject_id)
    if query_member is not None:
        metadata["query_orthogroup_representative"] = bool(query_member.representative)
    if subject_member is not None:
        metadata["subject_orthogroup_representative"] = bool(subject_member.representative)
    if not orthogroup_id:
        return metadata
    candidate_edges = [
        *orthogroups.ortholog_edges_by_orthogroup_id.get(orthogroup_id, ()),
        *orthogroups.related_edges_by_orthogroup_id.get(orthogroup_id, ()),
    ]
    for edge in candidate_edges:
        if (
            edge.query_protein_id == query_id
            and edge.subject_protein_id == subject_id
        ) or (
            edge.query_protein_id == subject_id
            and edge.subject_protein_id == query_id
        ):
            source_group = str(edge.source_rbh_orthogroup_id or "")
            target_group = str(edge.target_rbh_orthogroup_id or "")
            if source_group and target_group and source_group != target_group:
                rbh_group = f"{source_group};{target_group}"
            else:
                rbh_group = source_group or target_group
            metadata.update(
                {
                    "rbh_orthogroup_id": rbh_group,
                    "ortholog_path_id": str(edge.path_id or ""),
                    "edge_kind": edge.edge_kind,
                    "render_role": edge.render_role,
                }
            )
            return metadata
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
    orthogroups: OrthogroupResult,
    member_counts_by_record: Mapping[tuple[str, int], int],
) -> CollinearityAnchor | None:
    query_id = str(getattr(row, "query"))
    subject_id = str(getattr(row, "subject"))
    query_protein = protein_map.get(query_id)
    subject_protein = protein_map.get(subject_id)
    if query_protein is None or subject_protein is None:
        return None
    orthogroup_id = _orthogroup_id_for_edge(query_id, subject_id, orthogroups)
    edge_metadata = _orthogroup_edge_metadata_for_anchor(
        query_id,
        subject_id,
        orthogroup_id,
        orthogroups,
    )
    qstart, qend = _genomic_link_coordinates(query_protein)
    sstart, send = _genomic_link_coordinates(subject_protein)
    return CollinearityAnchor(
        query_protein_id=query_id,
        subject_protein_id=subject_id,
        query_record_index=int(query_record_index),
        subject_record_index=int(subject_record_index),
        query_order=int(order_by_id.get(query_id, query_protein.feature_index)),
        subject_order=int(order_by_id.get(subject_id, subject_protein.feature_index)),
        query_start=qstart,
        query_end=qend,
        subject_start=sstart,
        subject_end=send,
        query_strand=_normalized_strand(query_protein.strand),
        subject_strand=_normalized_strand(subject_protein.strand),
        identity=_float_from_row(row, "identity", 0.0),
        evalue=_float_from_row(row, "evalue", 1.0),
        bitscore=_float_from_row(row, "bitscore", 0.0),
        alignment_length=_int_from_row(row, "alignment_length", 0),
        query_feature_svg_id=str(query_protein.feature_svg_id or ""),
        subject_feature_svg_id=str(subject_protein.feature_svg_id or ""),
        source="orthogroup_display_edge",
        query_unit_id=query_id,
        subject_unit_id=subject_id,
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=_protein_locus_metadata(query_protein),
        subject_locus_id=_protein_locus_metadata(subject_protein),
        query_display_name=str(query_protein.label or query_id),
        subject_display_name=str(subject_protein.label or subject_id),
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
    orthogroups: OrthogroupResult,
) -> tuple[CollinearityAnchor, ...]:
    """Convert adjacent Orthogroup display edges one-to-one into anchors."""

    order_by_id = _protein_order_by_id(protein_map)
    member_counts_by_record = _orthogroup_member_counts_by_record(orthogroups)
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
                orthogroups=orthogroups,
                member_counts_by_record=member_counts_by_record,
            )
            if anchor is not None:
                anchors.append(anchor)
    if missing_ids:
        raise ParseError(
            "Orthogroup display edge contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )
    return tuple(sorted(anchors, key=_anchor_sort_key))


def _lossless_params_from_legacy(
    params: CollinearityParameters | LosslessCollinearityParameters | None,
) -> LosslessCollinearityParameters:
    if params is None:
        resolved = LosslessCollinearityParameters()
    elif isinstance(params, LosslessCollinearityParameters):
        resolved = params
    else:
        params.validate()
        resolved = LosslessCollinearityParameters(
            min_anchors=int(params.min_anchors),
            max_unit_gap=int(params.max_gene_gap),
            max_diagonal_drift=int(params.max_diagonal_drift),
            max_conflicts=int(params.max_conflicts_in_merge_gap),
        )
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


def _lossless_boundary_orientation(
    previous: CollinearityAnchor,
    current: CollinearityAnchor,
    params: LosslessCollinearityParameters,
) -> CollinearityOrientation:
    if params.merge_orientation in {"strand", "either"}:
        previous_strand = _anchor_strand_orientation(previous)
        current_strand = _anchor_strand_orientation(current)
        if previous_strand is not None and current_strand in {None, previous_strand}:
            return previous_strand
        if current_strand is not None and previous_strand is None:
            return current_strand
    return "minus" if int(current.subject_order) < int(previous.subject_order) else "plus"


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


def _lossless_block_score(anchors: Sequence[CollinearityAnchor]) -> float:
    return float(sum(float(anchor.bitscore) for anchor in anchors))


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
        score=_lossless_block_score(ordered),
        anchors=ordered,
        kind="cluster" if len(ordered) > 1 else "singleton",
        block_evalue=None,
    )


def _lossless_initial_clusters_for_pair(
    pair: tuple[int, int],
    anchors: Sequence[CollinearityAnchor],
    params: LosslessCollinearityParameters,
) -> list[CollinearityBlock]:
    ordered = sorted(
        anchors,
        key=lambda anchor: (
            int(anchor.query_order),
            int(anchor.subject_order),
            str(anchor.orthogroup_id),
            str(anchor.query_protein_id),
            str(anchor.subject_protein_id),
        ),
    )
    blocks: list[CollinearityBlock] = []
    current: list[CollinearityAnchor] = []
    current_orientation: CollinearityOrientation | None = None

    def flush() -> None:
        nonlocal current, current_orientation
        if not current:
            return
        orientation = current_orientation or _lossless_singleton_orientation(current[0])
        blocks.append(
            _lossless_block_from_anchors(
                block_id=f"pending_{len(blocks) + 1:04d}",
                pair=pair,
                orientation=orientation,
                anchors=current,
            )
        )
        current = []
        current_orientation = None

    for anchor in ordered:
        if not current:
            current = [anchor]
            current_orientation = None
            continue
        previous = current[-1]
        orientation = current_orientation or _lossless_boundary_orientation(previous, anchor, params)
        if _lossless_anchors_are_compatible(
            previous,
            anchor,
            orientation=orientation,
            params=params,
        ):
            current.append(anchor)
            current_orientation = orientation
            continue
        flush()
        current = [anchor]
    flush()
    return blocks


def _lossless_clusters_can_merge(
    left: CollinearityBlock,
    right: CollinearityBlock,
    params: LosslessCollinearityParameters,
) -> bool:
    if left.kind == "singleton" and right.kind == "singleton":
        return False
    if left.orientation != right.orientation:
        return False
    left_path = _path_sorted_anchors(left.anchors, left.orientation)
    right_path = _path_sorted_anchors(right.anchors, right.orientation)
    if not left_path or not right_path:
        return False
    return _lossless_anchors_are_compatible(
        left_path[-1],
        right_path[0],
        orientation=left.orientation,
        params=params,
    )


def _merge_lossless_clusters(
    blocks: Sequence[CollinearityBlock],
    params: LosslessCollinearityParameters,
) -> tuple[CollinearityBlock, ...]:
    merged: list[CollinearityBlock] = []
    for block in sorted(blocks, key=_final_block_sort_key):
        if not merged or not _lossless_clusters_can_merge(merged[-1], block, params):
            merged.append(block)
            continue
        previous = merged[-1]
        merged[-1] = _lossless_block_from_anchors(
            block_id=previous.block_id,
            pair=(previous.query_record_index, previous.subject_record_index),
            orientation=previous.orientation,
            anchors=(*previous.anchors, *block.anchors),
        )
    return tuple(merged)


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
        pair_blocks = _lossless_initial_clusters_for_pair(
            pair,
            anchors_by_pair[pair],
            resolved_params,
        )
        blocks.extend(_merge_lossless_clusters(pair_blocks, resolved_params))
    filtered_blocks, unblocked_anchors = _filter_blocks_by_min_anchors(
        blocks,
        min_anchors=resolved_params.min_anchors,
    )
    return CollinearityResult(
        blocks=_renumber_lossless_blocks(filtered_blocks),
        unblocked_anchors=unblocked_anchors,
    )


def _renumber_blocks(blocks: Sequence[CollinearityBlock], params: CollinearityParameters) -> tuple[CollinearityBlock, ...]:
    renumbered: list[CollinearityBlock] = []
    for index, block in enumerate(sorted(blocks, key=_final_block_sort_key), start=1):
        renumbered.append(
            _block_from_anchors(
                block_id=f"block_{index:04d}",
                pair=(block.query_record_index, block.subject_record_index),
                orientation=block.orientation,
                kind=block.kind,
                anchors=block.anchors,
                params=params,
            )
        )
    return tuple(renumbered)


def _chain_rank(chain: _Chain) -> tuple[float, int, int, int, int, int, str, str, int]:
    first = chain.anchors[0]
    orientation_rank = 0 if chain.orientation == "plus" else 1
    strand_supported, strand_conflicted = _chain_strand_support(chain)
    return (
        -float(chain.score),
        -len(chain.anchors),
        -strand_supported,
        strand_conflicted,
        int(first.query_order),
        int(first.subject_order),
        str(first.query_protein_id),
        str(first.subject_protein_id),
        orientation_rank,
    )


def _chain_path_anchors(chain: _Chain) -> tuple[CollinearityAnchor, ...]:
    if chain.orientation == "minus":
        return tuple(
            sorted(
                chain.anchors,
                key=lambda anchor: (
                    int(anchor.query_order),
                    -int(anchor.subject_order),
                    str(anchor.query_protein_id),
                    str(anchor.subject_protein_id),
                ),
            )
        )
    return tuple(
        sorted(
            chain.anchors,
            key=lambda anchor: (
                int(anchor.query_order),
                int(anchor.subject_order),
                str(anchor.query_protein_id),
                str(anchor.subject_protein_id),
            ),
        )
    )


def _log_permutation(n: int, m: int) -> float:
    if m < 0 or n < m:
        return math.inf
    if m == 0:
        return 0.0
    return math.lgamma(n + 1) - math.lgamma(n - m + 1)


def _anchor_genomic_position(anchor: CollinearityAnchor, start_attr: str, end_attr: str) -> float:
    start = float(getattr(anchor, start_attr))
    end = float(getattr(anchor, end_attr))
    return min(start, end)


def calculate_collinearity_block_evalue(
    chain: _Chain,
    available_anchors: Sequence[CollinearityAnchor],
) -> float:
    """Calculate a MCScanX-style block expectation value in genomic coordinate space."""

    path_anchors = _chain_path_anchors(chain)
    m = len(path_anchors)
    if m == 0:
        return math.inf

    query_positions = [
        _anchor_genomic_position(anchor, "query_start", "query_end")
        for anchor in path_anchors
    ]
    subject_positions = [
        _anchor_genomic_position(anchor, "subject_start", "subject_end")
        for anchor in path_anchors
    ]
    query_min = min(query_positions)
    query_max = max(query_positions)
    subject_min = min(subject_positions)
    subject_max = max(subject_positions)
    region_query_length = max(1.0, query_max - query_min)
    region_subject_length = max(1.0, subject_max - subject_min)
    n = sum(
        1
        for anchor in available_anchors
        if query_min <= _anchor_genomic_position(anchor, "query_start", "query_end") <= query_max
        and subject_min <= _anchor_genomic_position(anchor, "subject_start", "subject_end") <= subject_max
    )
    if n < m:
        return math.inf

    log_evalue = (
        math.log(2.0)
        + _log_permutation(n, m)
        - (m - 1) * (math.log(region_query_length) + math.log(region_subject_length))
    )
    for previous, current in zip(path_anchors, path_anchors[1:]):
        query_distance = max(
            1.0,
            abs(
                _anchor_genomic_position(current, "query_start", "query_end")
                - _anchor_genomic_position(previous, "query_start", "query_end")
            ),
        )
        subject_distance = max(
            1.0,
            abs(
                _anchor_genomic_position(current, "subject_start", "subject_end")
                - _anchor_genomic_position(previous, "subject_start", "subject_end")
            ),
        )
        log_evalue += math.log(query_distance) + math.log(subject_distance)

    if log_evalue <= math.log(sys.float_info.min):
        return 0.0
    if log_evalue >= math.log(sys.float_info.max):
        return math.inf
    return float(math.exp(log_evalue))


def _is_chain_accepted(
    chain: _Chain,
    params: CollinearityParameters,
    *,
    block_evalue: float,
) -> bool:
    if len(chain.anchors) < int(params.min_anchors):
        return False
    if float(chain.score) < params.effective_min_block_score():
        return False
    if params.block_evalue is not None and float(block_evalue) > float(params.block_evalue):
        return False
    return True


def _find_best_chain(
    anchors: Sequence[CollinearityAnchor],
    *,
    orientation: CollinearityOrientation,
    params: CollinearityParameters,
) -> _Chain | None:
    if not anchors:
        return None
    max_subject_order = max(int(anchor.subject_order) for anchor in anchors)

    def transformed_subject_order(anchor: CollinearityAnchor) -> int:
        if orientation == "plus":
            return int(anchor.subject_order)
        return max_subject_order - int(anchor.subject_order)

    ordered = sorted(
        anchors,
        key=lambda anchor: (
            int(anchor.query_order),
            transformed_subject_order(anchor),
            str(anchor.query_protein_id),
            str(anchor.subject_protein_id),
        ),
    )
    path_scores: list[float] = []
    path_indices: list[tuple[int, ...]] = []

    for j, anchor_j in enumerate(ordered):
        best_score = _anchor_score(anchor_j, params)
        best_path = (j,)
        x_j = int(anchor_j.query_order)
        y_j = transformed_subject_order(anchor_j)
        for i in range(j):
            anchor_i = ordered[i]
            x_i = int(anchor_i.query_order)
            y_i = transformed_subject_order(anchor_i)
            dx = x_j - x_i - 1
            dy = y_j - y_i - 1
            if dx < 0 or dy < 0:
                continue
            if dx > int(params.max_gene_gap) or dy > int(params.max_gene_gap):
                continue
            candidate_score = (
                path_scores[i]
                + _anchor_score(anchor_j, params)
                - max(dx, dy) * float(params.gap_penalty)
            )
            candidate_path = path_indices[i] + (j,)
            candidate_chain = _Chain(
                orientation=orientation,
                score=candidate_score,
                anchors=tuple(ordered[index] for index in candidate_path),
            )
            current_chain = _Chain(
                orientation=orientation,
                score=best_score,
                anchors=tuple(ordered[index] for index in best_path),
            )
            if _chain_rank(candidate_chain) < _chain_rank(current_chain):
                best_score = candidate_score
                best_path = candidate_path
        path_scores.append(best_score)
        path_indices.append(best_path)

    best_chain: _Chain | None = None
    for score, indices in zip(path_scores, path_indices):
        chain = _Chain(
            orientation=orientation,
            score=score,
            anchors=tuple(ordered[index] for index in indices),
        )
        if best_chain is None or _chain_rank(chain) < _chain_rank(best_chain):
            best_chain = chain
    return best_chain


def call_collinearity_blocks(
    anchors: Sequence[CollinearityAnchor],
    *,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
) -> CollinearityResult:
    """Call orthogroup-order syntenic blocks and singleton links."""

    params = params or CollinearityParameters()
    params.validate()

    deduplicated = deduplicate_unit_pair_anchors(anchors)
    collapsed = collapse_nearby_duplicate_anchors(
        deduplicated,
        window=int(params.nearby_duplicate_window),
    )
    anchors_by_pair: dict[tuple[int, int], list[CollinearityAnchor]] = {}
    for anchor in collapsed:
        anchors_by_pair.setdefault(
            (int(anchor.query_record_index), int(anchor.subject_record_index)),
            [],
        ).append(anchor)

    blocks: list[CollinearityBlock] = []

    for pair_key in sorted(anchors_by_pair):
        pair_anchors = tuple(sorted(anchors_by_pair[pair_key], key=_anchor_sort_key))
        candidate_runs = [
            *_candidate_runs_for_orientation(pair_anchors, orientation="plus", params=params),
            *_candidate_runs_for_orientation(pair_anchors, orientation="minus", params=params),
        ]
        syntenic_runs = _resolve_overlapping_runs(candidate_runs, params=params)
        syntenic_blocks = [
            _block_from_anchors(
                block_id=f"pending_{index:04d}",
                pair=pair_key,
                orientation=run.orientation,
                kind="syntenic",
                anchors=run.anchors,
                params=params,
            )
            for index, run in enumerate(syntenic_runs, start=1)
        ]
        merged_blocks = _merge_syntenic_blocks(
            syntenic_blocks,
            anchors=pair_anchors,
            params=params,
        )
        consumed_ids = {
            _anchor_identity(anchor)
            for block in merged_blocks
            for anchor in block.anchors
        }
        singleton_candidates = [
            anchor
            for anchor in pair_anchors
            if _anchor_identity(anchor) not in consumed_ids
        ]
        absorbed_blocks, _ = _absorb_singleton_anchors(
            merged_blocks,
            singleton_candidates,
            params=params,
        )
        blocks.extend(absorbed_blocks)
        absorbed_ids = {
            _anchor_identity(anchor)
            for block in absorbed_blocks
            for anchor in block.anchors
        }
        emit_candidates = [
            anchor
            for anchor in pair_anchors
            if _anchor_identity(anchor) not in absorbed_ids
        ]
        for singleton_index, anchor in enumerate(sorted(emit_candidates, key=_anchor_sort_key), start=1):
            orientation = _anchor_singleton_orientation(anchor)
            blocks.append(
                _block_from_anchors(
                    block_id=f"singleton_{pair_key[0] + 1:04d}_{singleton_index:04d}",
                    pair=pair_key,
                    orientation=orientation,
                    kind="singleton",
                    anchors=(anchor,),
                    params=params,
                )
            )

    filtered_blocks, filtered_unblocked = _filter_blocks_by_min_anchors(
        blocks,
        min_anchors=params.min_anchors,
    )

    return CollinearityResult(
        blocks=_renumber_blocks(filtered_blocks, params),
        unblocked_anchors=filtered_unblocked,
    )


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
) -> dict[tuple[int, int], DataFrame]:
    normalized_scope = normalize_collinearity_search_scope(str(scope))
    if normalized_scope == "all":
        return dict(hit_tables)
    allowed_pairs = set(
        iter_collinearity_search_pairs(record_count, scope=normalized_scope)
    )
    return {
        (query_index, subject_index): hits
        for (query_index, subject_index), hits in hit_tables.items()
        if query_index != subject_index
        and (min(query_index, subject_index), max(query_index, subject_index)) in allowed_pairs
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
        elif normalized_edge_mode == "one_to_one":
            selected[(query_index, subject_index)] = select_reciprocal_best_hits(forward_hits)
        else:
            selected[(query_index, subject_index)] = forward_hits.copy()
    return selected


def _evidence_lookup(
    edge_tables: Mapping[tuple[int, int], DataFrame],
) -> dict[tuple[int, int, str, str], object]:
    lookup: dict[tuple[int, int, str, str], object] = {}
    for (query_index, subject_index), hits in edge_tables.items():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            key = (
                int(query_index),
                int(subject_index),
                str(row.query),
                str(row.subject),
            )
            current = lookup.get(key)
            if current is None or _anchor_strength_key_from_row(row) < _anchor_strength_key_from_row(current):
                lookup[key] = row
    return lookup


def _anchor_strength_key_from_row(row: object) -> tuple[float, float, float, int, str, str]:
    return (
        _float_from_row(row, "evalue", 1.0),
        -_float_from_row(row, "bitscore", 0.0),
        -_float_from_row(row, "identity", 0.0),
        -_int_from_row(row, "alignment_length", 0),
        str(getattr(row, "query", "")),
        str(getattr(row, "subject", "")),
    )


def _members_by_record(
    members: Sequence[OrthogroupMember],
) -> dict[int, list[OrthogroupMember]]:
    grouped: dict[int, list[OrthogroupMember]] = {}
    for member in members:
        grouped.setdefault(int(member.record_index), []).append(member)
    for record_members in grouped.values():
        record_members.sort(
            key=lambda member: (
                int(member.start),
                int(member.end),
                int(member.feature_index),
                str(member.protein_id),
            )
        )
    return grouped


def _member_unit_pairs(
    members: Sequence[OrthogroupMember],
    unit_index: CollinearityUnitIndex,
) -> list[tuple[OrthogroupMember, CollinearityUnit]]:
    pairs: list[tuple[OrthogroupMember, CollinearityUnit]] = []
    seen_unit_ids: set[str] = set()
    for member in members:
        unit = unit_index.unit_by_protein_id.get(str(member.protein_id))
        if unit is None or unit.unit_id in seen_unit_ids:
            continue
        seen_unit_ids.add(unit.unit_id)
        pairs.append((member, unit))
    return pairs


def _candidate_pair_rank(
    query_member: OrthogroupMember,
    query_unit: CollinearityUnit,
    subject_member: OrthogroupMember,
    subject_unit: CollinearityUnit,
    evidence: object | None,
) -> tuple[bool, bool, bool, bool, float, float, float, int, int, int, str, str]:
    return (
        not (bool(query_member.representative) and bool(subject_member.representative)),
        not bool(query_member.representative),
        not bool(subject_member.representative),
        evidence is None,
        _float_from_row(evidence, "evalue", 1.0),
        -_float_from_row(evidence, "bitscore", 0.0),
        -_float_from_row(evidence, "identity", 0.0),
        -_int_from_row(evidence, "alignment_length", 0),
        abs(int(query_unit.order) - int(subject_unit.order)),
        int(query_unit.order),
        str(query_member.protein_id),
        str(subject_member.protein_id),
    )


def _make_orthogroup_anchor(
    *,
    orthogroup_id: str,
    query_member: OrthogroupMember,
    query_unit: CollinearityUnit,
    subject_member: OrthogroupMember,
    subject_unit: CollinearityUnit,
    query_protein: CdsProtein,
    subject_protein: CdsProtein,
    evidence: object | None,
    query_member_count: int,
    subject_member_count: int,
) -> CollinearityAnchor:
    qstart, qend = _genomic_link_coordinates(query_protein)
    sstart, send = _genomic_link_coordinates(subject_protein)
    return CollinearityAnchor(
        query_protein_id=str(query_member.protein_id),
        subject_protein_id=str(subject_member.protein_id),
        query_record_index=int(query_unit.record_index),
        subject_record_index=int(subject_unit.record_index),
        query_order=int(query_unit.order),
        subject_order=int(subject_unit.order),
        query_start=qstart,
        query_end=qend,
        subject_start=sstart,
        subject_end=send,
        query_strand=_normalized_strand(query_protein.strand),
        subject_strand=_normalized_strand(subject_protein.strand),
        identity=_float_from_row(evidence, "identity", 0.0),
        evalue=_float_from_row(evidence, "evalue", 1.0),
        bitscore=_float_from_row(evidence, "bitscore", 0.0),
        alignment_length=_int_from_row(evidence, "alignment_length", 0),
        query_feature_svg_id=str(query_protein.feature_svg_id or query_unit.representative_feature_svg_id),
        subject_feature_svg_id=str(subject_protein.feature_svg_id or subject_unit.representative_feature_svg_id),
        source="orthogroup",
        query_unit_id=query_unit.unit_id,
        subject_unit_id=subject_unit.unit_id,
        query_unit_kind=query_unit.unit_kind,
        subject_unit_kind=subject_unit.unit_kind,
        query_locus_id=query_unit.locus_id,
        subject_locus_id=subject_unit.locus_id,
        query_display_name=query_unit.display_name,
        subject_display_name=subject_unit.display_name,
        orthogroup_id=str(orthogroup_id),
        rbh_orthogroup_id=str(orthogroup_id),
        edge_kind="rbh" if query_member.representative and subject_member.representative else "coortholog",
        render_role="display_edge",
        query_orthogroup_representative=bool(query_member.representative),
        subject_orthogroup_representative=bool(subject_member.representative),
        query_orthogroup_member_count=int(query_member_count),
        subject_orthogroup_member_count=int(subject_member_count),
    )


def build_orthogroup_collinearity_anchors(
    orthogroups: OrthogroupResult,
    extraction: ProteinExtractionResult,
    unit_index: CollinearityUnitIndex,
    edge_tables: Mapping[tuple[int, int], DataFrame],
    *,
    max_paralog_links_per_orthogroup: int = 2,
) -> tuple[CollinearityAnchor, ...]:
    """Build adjacent-record anchors from orthogroup membership."""

    evidence_by_pair = _evidence_lookup(edge_tables)
    anchors: list[CollinearityAnchor] = []
    record_count = len(extraction.proteins_by_record)
    for orthogroup_id in sorted(orthogroups.orthogroups):
        members_by_record = _members_by_record(orthogroups.orthogroups[orthogroup_id])
        for query_index in range(record_count - 1):
            subject_index = query_index + 1
            query_members = members_by_record.get(query_index, [])
            subject_members = members_by_record.get(subject_index, [])
            if not query_members or not subject_members:
                continue
            query_pairs = _member_unit_pairs(query_members, unit_index)
            subject_pairs = _member_unit_pairs(subject_members, unit_index)
            candidates: list[tuple[tuple[object, ...], OrthogroupMember, CollinearityUnit, OrthogroupMember, CollinearityUnit, object | None]] = []
            for query_member, query_unit in query_pairs:
                for subject_member, subject_unit in subject_pairs:
                    evidence = evidence_by_pair.get(
                        (
                            query_index,
                            subject_index,
                            str(query_member.protein_id),
                            str(subject_member.protein_id),
                        )
                    )
                    candidates.append(
                        (
                            _candidate_pair_rank(
                                query_member,
                                query_unit,
                                subject_member,
                                subject_unit,
                                evidence,
                            ),
                            query_member,
                            query_unit,
                            subject_member,
                            subject_unit,
                            evidence,
                        )
                    )
            used_query_units: set[str] = set()
            used_subject_units: set[str] = set()
            selected_count = 0
            for _rank, query_member, query_unit, subject_member, subject_unit, evidence in sorted(candidates, key=lambda item: item[0]):
                if query_unit.unit_id in used_query_units or subject_unit.unit_id in used_subject_units:
                    continue
                query_protein = extraction.protein_map.get(str(query_member.protein_id))
                subject_protein = extraction.protein_map.get(str(subject_member.protein_id))
                if query_protein is None or subject_protein is None:
                    continue
                anchors.append(
                    _make_orthogroup_anchor(
                        orthogroup_id=orthogroup_id,
                        query_member=query_member,
                        query_unit=query_unit,
                        subject_member=subject_member,
                        subject_unit=subject_unit,
                        query_protein=query_protein,
                        subject_protein=subject_protein,
                        evidence=evidence,
                        query_member_count=len(query_members),
                        subject_member_count=len(subject_members),
                    )
                )
                used_query_units.add(query_unit.unit_id)
                used_subject_units.add(subject_unit.unit_id)
                selected_count += 1
                if selected_count >= int(max_paralog_links_per_orthogroup):
                    break
    return tuple(sorted(anchors, key=_anchor_sort_key))


def build_orthogroup_collinearity_blocks_from_hits(
    hits_by_pair: Sequence[DataFrame] | Mapping[tuple[int, int], DataFrame],
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    edge_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "rbh",
    orthogroup_member_max_hits: int = 5,
    max_paralog_links_per_orthogroup: int = 2,
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> CollinearityResult:
    """Call lossless Orthogroup-sourced collinearity blocks from filtered hits."""

    lossless_params = _lossless_params_from_legacy(params)
    directional_tables = _normalize_hit_tables(hits_by_pair, reverse_hits_by_pair)
    normalized_edge_mode = normalize_collinearity_anchor_mode(str(edge_mode))
    normalized_search_scope = normalize_collinearity_search_scope(str(search_scope))
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if int(orthogroup_member_max_hits) <= 0:
        raise ValidationError("orthogroup_member_max_hits must be > 0")
    resolved_max_paralog_links = int(max_paralog_links_per_orthogroup)
    if (
        isinstance(params, CollinearityParameters)
        and resolved_max_paralog_links == 2
    ):
        resolved_max_paralog_links = int(params.max_paralog_links_per_orthogroup)
    if resolved_max_paralog_links <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    record_count = len(records) if records is not None else len(extraction.proteins_by_record)
    directional_tables = _filter_hit_tables_by_search_scope(
        directional_tables,
        record_count=record_count,
        scope=normalized_search_scope,
    )
    if normalized_edge_mode == "rbh":
        edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
            directional_tables,
            extraction.protein_map,
            record_count=record_count,
            orthogroup_membership_mode=normalized_membership_mode,
            orthogroup_member_max_hits=int(orthogroup_member_max_hits),
            max_related_edges_per_orthogroup=resolved_max_paralog_links,
        )
        orthogroups = edge_selection.orthogroups
        adjacent_edges_by_pair = edge_selection.adjacent_anchor_edges_by_pair
    else:
        edge_tables = _select_orthogroup_edges(directional_tables, edge_mode=normalized_edge_mode)
        if normalized_membership_mode == "rbh":
            orthogroups = build_orthogroups_from_protein_hits(
                tuple(edge_tables.values()),
                extraction.protein_map,
                include_singletons=False,
            )
        else:
            seed_selection = select_rbh_orthogroup_edges_from_directional_hits(
                directional_tables,
                extraction.protein_map,
                record_count=record_count,
                orthogroup_membership_mode=normalized_membership_mode,
                orthogroup_member_max_hits=int(orthogroup_member_max_hits),
                max_related_edges_per_orthogroup=resolved_max_paralog_links,
            )
            orthogroups = seed_selection.orthogroups
        adjacent_edges_by_pair = {
            pair: table
            for pair, table in edge_tables.items()
            if int(pair[1]) == int(pair[0]) + 1
        }
        for query_index in range(max(0, record_count - 1)):
            adjacent_edges_by_pair.setdefault((query_index, query_index + 1), _empty_hits())
    anchors = orthogroup_edges_to_lossless_collinearity_anchors(
        adjacent_edges_by_pair,
        extraction.protein_map,
        orthogroups,
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
    losatp_threads: int | None = None,
    candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "rbh",
    orthogroup_member_max_hits: int = 5,
    max_paralog_links_per_orthogroup: int = 2,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    edge_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    runner: LosatpRunner | None = None,
) -> CollinearityResult:
    """Run scoped LOSATP blastp evidence searches, infer orthogroups, and call blocks."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='collinear' requires at least two records")
    if losatp_threads is not None and int(losatp_threads) <= 0:
        raise ValidationError("losatp_threads must be > 0 or None")
    lossless_params = _lossless_params_from_legacy(params)
    normalized_edge_mode = normalize_collinearity_anchor_mode(str(edge_mode))
    normalized_search_scope = normalize_collinearity_search_scope(str(search_scope))
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if int(orthogroup_member_max_hits) <= 0:
        raise ValidationError("orthogroup_member_max_hits must be > 0")
    resolved_max_paralog_links = int(max_paralog_links_per_orthogroup)
    if (
        isinstance(params, CollinearityParameters)
        and resolved_max_paralog_links == 2
    ):
        resolved_max_paralog_links = int(params.max_paralog_links_per_orthogroup)
    if resolved_max_paralog_links <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    extraction = extract_cds_proteins(records)
    _validate_collinearity_extraction(records, extraction)

    search_candidate_limit = (
        candidate_limit
        if candidate_limit is not None
        else (
            None
            if normalized_edge_mode == "all"
            else (1 if normalized_membership_mode == "rbh" else int(orthogroup_member_max_hits))
        )
    )
    directional_tables: dict[tuple[int, int], DataFrame] = {}
    for query_index, subject_index in iter_collinearity_search_pairs(
        len(records),
        scope=normalized_search_scope,
    ):
        query_fasta = proteins_to_fasta(extraction.proteins_by_record[query_index])
        subject_fasta = proteins_to_fasta(extraction.proteins_by_record[subject_index])
        forward_hits = _run_losatp_search(
            query_fasta,
            subject_fasta,
            losatp_bin=losatp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=search_candidate_limit,
            runner=runner,
        )
        directional_tables[(query_index, subject_index)] = filter_protein_hits_by_thresholds(
            forward_hits,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        if normalized_edge_mode == "rbh" or normalized_membership_mode != "rbh":
            reverse_hits = _run_losatp_search(
                subject_fasta,
                query_fasta,
                losatp_bin=losatp_bin,
                losatp_threads=losatp_threads,
                candidate_limit=search_candidate_limit,
                runner=runner,
            )
            directional_tables[(subject_index, query_index)] = filter_protein_hits_by_thresholds(
                reverse_hits,
                evalue=evalue,
                bitscore=bitscore,
                identity=identity,
                alignment_length=alignment_length,
            )

    return build_orthogroup_collinearity_blocks_from_hits(
        directional_tables,
        extraction,
        records=records,
        params=lossless_params,
        unit_mode=unit_mode,
        edge_mode=normalized_edge_mode,
        search_scope=normalized_search_scope,
        orthogroup_membership_mode=normalized_membership_mode,
        orthogroup_member_max_hits=int(orthogroup_member_max_hits),
        max_paralog_links_per_orthogroup=resolved_max_paralog_links,
    )


def build_native_collinearity_blocks(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    losatp_threads: int | None = None,
    candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "rbh",
    orthogroup_member_max_hits: int = 5,
    max_paralog_links_per_orthogroup: int = 2,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    anchor_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    runner: LosatpRunner | None = None,
) -> CollinearityResult:
    """Compatibility entry point for the orthogroup-first collinearity caller."""

    return build_orthogroup_collinearity_blocks(
        records,
        losatp_bin=losatp_bin,
        losatp_threads=losatp_threads,
        candidate_limit=candidate_limit,
        orthogroup_membership_mode=orthogroup_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        max_paralog_links_per_orthogroup=max_paralog_links_per_orthogroup,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        params=params,
        unit_mode=unit_mode,
        edge_mode=anchor_mode,
        search_scope=search_scope,
        runner=runner,
    )


def build_collinearity_blocks_from_hits(
    hits_by_pair: Sequence[DataFrame] | Mapping[tuple[int, int], DataFrame],
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    anchor_mode: CollinearityAnchorMode | str = "rbh",
    search_scope: CollinearitySearchScope | str = "adjacent",
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "rbh",
    orthogroup_member_max_hits: int = 5,
    max_paralog_links_per_orthogroup: int = 2,
    reverse_hits_by_pair: Sequence[DataFrame] | None = None,
) -> CollinearityResult:
    """Call orthogroup-first blocks from already filtered directional hit tables."""

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


def _coordinate_span(anchors: Sequence[CollinearityAnchor], start_attr: str, end_attr: str) -> tuple[int, int]:
    starts: list[int] = []
    ends: list[int] = []
    for anchor in anchors:
        a = int(getattr(anchor, start_attr))
        b = int(getattr(anchor, end_attr))
        starts.append(min(a, b))
        ends.append(max(a, b))
    return min(starts), max(ends)


def _weighted_identity(anchors: Sequence[CollinearityAnchor]) -> float:
    total_length = sum(max(0, int(anchor.alignment_length)) for anchor in anchors)
    if total_length > 0:
        return sum(
            float(anchor.identity) * max(0, int(anchor.alignment_length))
            for anchor in anchors
        ) / total_length
    return sum(float(anchor.identity) for anchor in anchors) / len(anchors)


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


def convert_collinearity_blocks_to_comparisons(
    result: CollinearityResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    record_ids: Sequence[str] | None = None,
    color_mode: CollinearityColorMode | str = "orientation",
) -> list[DataFrame]:
    """Convert accepted blocks into existing linear comparison rows.

    The native TSV keeps anchor-level details, but the diagram should display a
    collinear block as one continuous span ribbon.  The existing comparison
    renderer can draw that block ribbon when we pass it one span row per block.
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
        return []

    rows_by_pair: list[list[dict[str, object]]] = [[] for _ in range(record_count - 1)]
    for block in result.blocks:
        pair_index = int(block.query_record_index)
        if int(block.subject_record_index) != pair_index + 1:
            raise ValidationError("Collinearity rendering supports adjacent record pairs only.")
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
        query_start, query_end = _coordinate_span(block.anchors, "query_start", "query_end")
        subject_min, subject_max = _coordinate_span(block.anchors, "subject_start", "subject_end")
        if block.orientation == "minus":
            subject_start, subject_end = subject_max, subject_min
        else:
            subject_start, subject_end = subject_min, subject_max
        rows_by_pair[pair_index].append(
            {
                "query": query_name,
                "subject": subject_name,
                "identity": float(_weighted_identity(block.anchors)),
                "alignment_length": sum(max(0, int(anchor.alignment_length)) for anchor in block.anchors),
                "mismatches": 0,
                "gap_opens": 0,
                "qstart": query_start,
                "qend": query_end,
                "sstart": subject_start,
                "send": subject_end,
                "evalue": min(float(anchor.evalue) for anchor in block.anchors),
                "bitscore": sum(float(anchor.bitscore) for anchor in block.anchors),
                "collinearity_block_id": block.block_id,
                "collinearity_block_kind": block.kind,
                "collinearity_orientation": block.orientation,
                "collinearity_block_score": float(block.score),
                "collinearity_block_evalue": block.block_evalue,
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
                "query_feature_svg_id": _joined_anchor_values(block.anchors, "query_feature_svg_id"),
                "subject_feature_svg_id": _joined_anchor_values(block.anchors, "subject_feature_svg_id"),
                "rbh_orthogroup_id": _joined_anchor_values(block.anchors, "rbh_orthogroup_id"),
                "ortholog_path_id": _joined_anchor_values(block.anchors, "ortholog_path_id"),
                "edge_kind": _joined_anchor_values(block.anchors, "edge_kind") or "rbh",
                "render_role": _joined_anchor_values(block.anchors, "render_role") or "block_anchor",
                "query_orthogroup_representative": _joined_anchor_values(block.anchors, "query_orthogroup_representative"),
                "subject_orthogroup_representative": _joined_anchor_values(block.anchors, "subject_orthogroup_representative"),
                "query_orthogroup_member_count": _joined_anchor_values(block.anchors, "query_orthogroup_member_count"),
                "subject_orthogroup_member_count": _joined_anchor_values(block.anchors, "subject_orthogroup_member_count"),
            }
        )

    return [
        pd.DataFrame.from_records(rows, columns=COLLINEARITY_COMPARISON_COLUMNS)
        for rows in rows_by_pair
    ]


__all__ = [
    "COLLINEARITY_ANCHOR_MODES",
    "COLLINEARITY_COLOR_MODES",
    "COLLINEARITY_COMPARISON_COLUMNS",
    "COLLINEARITY_METADATA_COLUMNS",
    "COLLINEARITY_SEARCH_SCOPES",
    "CollinearityAnchor",
    "CollinearityAnchorMode",
    "CollinearityBlock",
    "CollinearityBlockKind",
    "CollinearityColorMode",
    "CollinearityParameters",
    "CollinearityResult",
    "CollinearitySearchScope",
    "LosslessCollinearityParameters",
    "build_collinearity_blocks_from_hits",
    "build_native_collinearity_blocks",
    "build_orthogroup_collinearity_anchors",
    "build_orthogroup_collinearity_blocks",
    "build_orthogroup_collinearity_blocks_from_hits",
    "call_collinearity_blocks",
    "calculate_collinearity_block_evalue",
    "cluster_lossless_collinearity_anchors",
    "collapse_nearby_duplicate_anchors",
    "convert_collinearity_blocks_to_comparisons",
    "deduplicate_unit_pair_anchors",
    "iter_collinearity_search_pairs",
    "normalize_collinearity_anchor_mode",
    "normalize_collinearity_color_mode",
    "normalize_collinearity_search_scope",
    "orthogroup_edges_to_lossless_collinearity_anchors",
    "protein_hits_to_collinearity_anchors",
    "select_collinearity_anchor_hits",
]
