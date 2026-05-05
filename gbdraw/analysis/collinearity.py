#!/usr/bin/env python
# coding: utf-8

"""Native gbdraw collinear block calling."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity_units import (
    CollinearityUnit,
    CollinearityUnitMode,
    build_collinearity_unit_index,
)
from gbdraw.analysis.protein_colinearity import (
    CdsProtein,
    LosatpRunner,
    ProteinExtractionResult,
    _run_losatp_search,
    extract_cds_proteins,
    filter_protein_hits_by_thresholds,
    proteins_to_fasta,
)
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS

CollinearityOrientation = Literal["plus", "minus"]
CollinearityScoreMode = Literal["constant", "bitscore"]
CollinearityColorMode = Literal["average_identity", "orientation"]

COLLINEARITY_METADATA_COLUMNS = (
    "collinearity_block_id",
    "collinearity_orientation",
    "collinearity_block_score",
    "collinearity_anchor_index",
    "collinearity_anchor_count",
    "collinearity_color_mode",
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
)
COLLINEARITY_COMPARISON_COLUMNS = tuple(COMPARISON_COLUMNS) + COLLINEARITY_METADATA_COLUMNS
COLLINEARITY_COLOR_MODES = ("average_identity", "orientation")


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


@dataclass(frozen=True)
class CollinearityBlock:
    block_id: str
    query_record_index: int
    subject_record_index: int
    orientation: CollinearityOrientation
    score: float
    anchors: tuple[CollinearityAnchor, ...]
    label: str | None = None


@dataclass(frozen=True)
class CollinearityResult:
    blocks: tuple[CollinearityBlock, ...]
    unblocked_anchors: tuple[CollinearityAnchor, ...] = ()


@dataclass(frozen=True)
class CollinearityParameters:
    min_anchors: int = 5
    max_gene_gap: int = 25
    gap_penalty: float = 1.0
    nearby_duplicate_window: int = 5
    score_mode: CollinearityScoreMode = "constant"
    constant_anchor_score: float = 50.0
    min_block_score: float | None = None

    def validate(self) -> None:
        if int(self.min_anchors) <= 0:
            raise ValidationError("collinear_min_anchors must be > 0")
        if int(self.max_gene_gap) < 0:
            raise ValidationError("collinear_max_gene_gap must be >= 0")
        if float(self.gap_penalty) < 0:
            raise ValidationError("collinear_gap_penalty must be >= 0")
        if int(self.nearby_duplicate_window) < 0:
            raise ValidationError("collinear_nearby_duplicate_window must be >= 0")
        if str(self.score_mode) not in {"constant", "bitscore"}:
            raise ValidationError("collinear_score_mode must be one of: constant, bitscore")
        if float(self.constant_anchor_score) <= 0:
            raise ValidationError("collinear_constant_anchor_score must be > 0")

    def effective_min_block_score(self) -> float:
        if self.min_block_score is not None:
            return float(self.min_block_score)
        if self.score_mode == "constant":
            return float(self.constant_anchor_score) * int(self.min_anchors)
        return 0.0


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


def _chain_rank(chain: _Chain) -> tuple[float, int, int, int, str, str, int]:
    first = chain.anchors[0]
    orientation_rank = 0 if chain.orientation == "plus" else 1
    return (
        -float(chain.score),
        -len(chain.anchors),
        int(first.query_order),
        int(first.subject_order),
        str(first.query_protein_id),
        str(first.subject_protein_id),
        orientation_rank,
    )


def _is_chain_accepted(chain: _Chain, params: CollinearityParameters) -> bool:
    return (
        len(chain.anchors) >= int(params.min_anchors)
        and float(chain.score) >= params.effective_min_block_score()
    )


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
    params: CollinearityParameters | None = None,
) -> CollinearityResult:
    """Call plus/minus collinear blocks from unit-order anchors."""

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
    consumed_ids: set[tuple[str, str, int, int]] = set()
    rejected_anchors: set[CollinearityAnchor] = set(collapsed)
    block_counter = 0

    for pair_key in sorted(anchors_by_pair):
        available = list(anchors_by_pair[pair_key])
        while available:
            candidate_chains = [
                chain
                for chain in (
                    _find_best_chain(available, orientation="plus", params=params),
                    _find_best_chain(available, orientation="minus", params=params),
                )
                if chain is not None and _is_chain_accepted(chain, params)
            ]
            if not candidate_chains:
                break
            best_chain = min(candidate_chains, key=_chain_rank)
            block_counter += 1
            block_id = f"block_{block_counter:04d}"
            chain_anchor_ids = {
                (
                    anchor.query_unit_id,
                    anchor.subject_unit_id,
                    int(anchor.query_record_index),
                    int(anchor.subject_record_index),
                )
                for anchor in best_chain.anchors
            }
            consumed_ids.update(chain_anchor_ids)
            for anchor in best_chain.anchors:
                rejected_anchors.discard(anchor)
            blocks.append(
                CollinearityBlock(
                    block_id=block_id,
                    query_record_index=pair_key[0],
                    subject_record_index=pair_key[1],
                    orientation=best_chain.orientation,
                    score=float(best_chain.score),
                    anchors=tuple(
                        sorted(
                            best_chain.anchors,
                            key=lambda anchor: (
                                int(anchor.query_order),
                                int(anchor.subject_order),
                                str(anchor.query_protein_id),
                                str(anchor.subject_protein_id),
                            ),
                        )
                    ),
                )
            )
            available = [
                anchor
                for anchor in available
                if (
                    anchor.query_unit_id,
                    anchor.subject_unit_id,
                    int(anchor.query_record_index),
                    int(anchor.subject_record_index),
                )
                not in consumed_ids
            ]

    return CollinearityResult(
        blocks=tuple(blocks),
        unblocked_anchors=tuple(
            sorted(
                rejected_anchors,
                key=lambda anchor: (
                    anchor.query_record_index,
                    anchor.subject_record_index,
                    anchor.query_order,
                    anchor.subject_order,
                    anchor.query_protein_id,
                    anchor.subject_protein_id,
                ),
            )
        ),
    )


def build_native_collinearity_blocks(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    candidate_limit: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    params: CollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
    runner: LosatpRunner | None = None,
) -> CollinearityResult:
    """Run LOSATP blastp on adjacent records and call native blocks."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='collinear' requires at least two records")
    params = params or CollinearityParameters()
    params.validate()
    extraction = extract_cds_proteins(records)
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
    unit_index = build_collinearity_unit_index(
        extraction,
        records=records,
        mode=unit_mode,
    )
    all_anchors: list[CollinearityAnchor] = []
    for record_index in range(len(records) - 1):
        query_units = unit_index.units_by_record[record_index]
        subject_units = unit_index.units_by_record[record_index + 1]
        query_representatives = _representative_proteins(query_units, extraction.protein_map)
        subject_representatives = _representative_proteins(subject_units, extraction.protein_map)
        protein_hits = _run_losatp_search(
            proteins_to_fasta(query_representatives),
            proteins_to_fasta(subject_representatives),
            losatp_bin=losatp_bin,
            candidate_limit=candidate_limit,
            runner=runner,
        )
        filtered_hits = filter_protein_hits_by_thresholds(
            protein_hits,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        all_anchors.extend(
            protein_hits_to_collinearity_anchors(
                filtered_hits,
                query_units=unit_index.unit_by_protein_id,
                subject_units=unit_index.unit_by_protein_id,
                query_protein_map=extraction.protein_map,
                subject_protein_map=extraction.protein_map,
                source="losatp",
            )
        )
    return call_collinearity_blocks(all_anchors, params=params)


def build_collinearity_blocks_from_hits(
    hits_by_pair: Sequence[DataFrame],
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    params: CollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
) -> CollinearityResult:
    """Call native blocks from already filtered adjacent-pair protein hit tables."""

    unit_index = build_collinearity_unit_index(
        extraction,
        records=records,
        mode=unit_mode,
    )
    anchors: list[CollinearityAnchor] = []
    for pair_index, hits in enumerate(hits_by_pair):
        if hits is None or hits.empty:
            continue
        anchors.extend(
            protein_hits_to_collinearity_anchors(
                hits,
                query_units=unit_index.unit_by_protein_id,
                subject_units=unit_index.unit_by_protein_id,
                query_protein_map=extraction.protein_map,
                subject_protein_map=extraction.protein_map,
                source="losatp",
            )
        )
    return call_collinearity_blocks(anchors, params=params)


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
                "collinearity_orientation": block.orientation,
                "collinearity_block_score": float(block.score),
                "collinearity_anchor_index": 1,
                "collinearity_anchor_count": len(block.anchors),
                "collinearity_color_mode": normalized_color_mode,
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
            }
        )

    return [
        pd.DataFrame.from_records(rows, columns=COLLINEARITY_COMPARISON_COLUMNS)
        for rows in rows_by_pair
    ]


__all__ = [
    "COLLINEARITY_COLOR_MODES",
    "COLLINEARITY_COMPARISON_COLUMNS",
    "COLLINEARITY_METADATA_COLUMNS",
    "CollinearityAnchor",
    "CollinearityBlock",
    "CollinearityColorMode",
    "CollinearityParameters",
    "CollinearityResult",
    "build_collinearity_blocks_from_hits",
    "build_native_collinearity_blocks",
    "call_collinearity_blocks",
    "collapse_nearby_duplicate_anchors",
    "convert_collinearity_blocks_to_comparisons",
    "deduplicate_unit_pair_anchors",
    "normalize_collinearity_color_mode",
    "protein_hits_to_collinearity_anchors",
]
