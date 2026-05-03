#!/usr/bin/env python
# coding: utf-8

"""Protein-level colinearity helpers for linear pairwise comparisons."""

from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO
import logging
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Callable, Literal, Mapping, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from Bio.SeqFeature import SeqFeature  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.features.colors import compute_feature_hash
from gbdraw.io.comparisons import COMPARISON_COLUMNS

logger = logging.getLogger(__name__)

_VALID_PROTEIN_RE = re.compile(r"^[A-Z]+$")
_VALID_FASTA_ID_RE = re.compile(r"^\S+$")
_NUMERIC_COMPARISON_COLUMNS = COMPARISON_COLUMNS[2:]
LOSATP_METADATA_COLUMNS = (
    "query_protein_id",
    "subject_protein_id",
    "query_source_protein_id",
    "subject_source_protein_id",
    "query_record_index",
    "subject_record_index",
    "query_feature_index",
    "subject_feature_index",
    "query_feature_svg_id",
    "subject_feature_svg_id",
    "orthogroup_id",
    "query_orthogroup_representative",
    "subject_orthogroup_representative",
)
LOSATP_COMPARISON_COLUMNS = tuple(COMPARISON_COLUMNS) + LOSATP_METADATA_COLUMNS
PROTEIN_BLASTP_MODES = ("none", "pairwise", "orthogroup")
ProteinBlastpMode = Literal["none", "pairwise", "orthogroup"]


@dataclass(frozen=True)
class CdsProtein:
    """A protein sequence and its genomic CDS coordinates."""

    protein_id: str
    record_index: int
    feature_index: int
    record_id: str
    start: int
    end: int
    strand: int | None
    label: str
    protein_length: int
    sequence: str
    source_protein_id: str | None = None
    feature_svg_id: str | None = None
    gene: str | None = None
    product: str | None = None
    note: str | None = None


@dataclass(frozen=True)
class ProteinExtractionResult:
    """CDS protein extraction output grouped for pairwise LOSATP runs."""

    proteins_by_record: list[list[CdsProtein]]
    protein_map: dict[str, CdsProtein]


@dataclass(frozen=True)
class OrthogroupMember:
    """A CDS/protein member assigned to an orthogroup."""

    orthogroup_id: str
    protein_id: str
    record_index: int
    feature_index: int
    record_id: str
    label: str
    start: int
    end: int
    strand: int | None
    feature_svg_id: str | None
    source_protein_id: str | None
    gene: str | None = None
    product: str | None = None
    note: str | None = None
    representative: bool = False


@dataclass(frozen=True)
class OrthogroupNameCandidate:
    """Annotation-derived display-name candidate for an orthogroup."""

    text: str
    source: str
    member_count: int
    record_coverage_count: int
    representative_count: int
    score: float


@dataclass(frozen=True)
class OrthogroupResult:
    """Connected-component orthogroups and per-protein lookup metadata."""

    orthogroups: dict[str, list[OrthogroupMember]]
    member_by_protein_id: dict[str, OrthogroupMember]
    names_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    descriptions_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    name_candidates_by_orthogroup_id: dict[str, list[OrthogroupNameCandidate]] = field(default_factory=dict)
    confidence_by_orthogroup_id: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class ProteinBlastpResult:
    """LOSATP blastp display comparisons plus optional orthogroup metadata."""

    comparisons: list[DataFrame]
    orthogroups: OrthogroupResult | None = None


@dataclass(frozen=True)
class _CdsProteinCandidate:
    protein_id: str
    source_protein_id: str | None
    record_index: int
    feature_index: int
    record_id: str
    start: int
    end: int
    strand: int | None
    label: str
    protein_length: int
    sequence: str
    feature_svg_id: str | None
    gene: str | None
    product: str | None
    note: str | None


LosatpRunner = Callable[[str, str], DataFrame]


def _first_qualifier(feature: SeqFeature, key: str) -> str | None:
    value = feature.qualifiers.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        value = value[0]
    text = str(value).strip()
    return text or None


def _first_qualifier_from_any(feature: SeqFeature, keys: Sequence[str]) -> str | None:
    for key in keys:
        value = _first_qualifier(feature, key)
        if value:
            return value
    return None


def _fasta_safe_protein_id(protein_id: str | None) -> str | None:
    if protein_id is None:
        return None
    return protein_id if _VALID_FASTA_ID_RE.fullmatch(protein_id) else None


def _unique_protein_id(
    preferred_id: str | None,
    fallback_id: str,
    used_ids: set[str],
) -> str:
    if preferred_id is not None and preferred_id not in used_ids:
        return preferred_id
    if fallback_id not in used_ids:
        return fallback_id
    suffix = 2
    while f"{fallback_id}_{suffix}" in used_ids:
        suffix += 1
    return f"{fallback_id}_{suffix}"


def _clean_protein_sequence(sequence: object | None) -> str | None:
    if sequence is None:
        return None
    protein = "".join(str(sequence).split()).upper()
    while protein.endswith("*"):
        protein = protein[:-1]
    if not protein or "*" in protein:
        return None
    if not _VALID_PROTEIN_RE.fullmatch(protein):
        return None
    return protein


def _translation_table(feature: SeqFeature) -> int | str:
    raw_table = _first_qualifier(feature, "transl_table")
    if raw_table is None:
        return 1
    try:
        return int(raw_table)
    except ValueError:
        return raw_table


def _codon_start(feature: SeqFeature) -> int:
    raw_start = _first_qualifier(feature, "codon_start")
    if raw_start is None:
        return 1
    try:
        codon_start = int(raw_start)
    except ValueError:
        return 1
    return codon_start if codon_start in {1, 2, 3} else 1


def _translate_cds_feature(record: SeqRecord, feature: SeqFeature) -> str | None:
    try:
        nucleotide_sequence = feature.extract(record.seq)
        offset = _codon_start(feature) - 1
        if offset:
            nucleotide_sequence = nucleotide_sequence[offset:]
        protein = nucleotide_sequence.translate(table=_translation_table(feature), to_stop=False)
    except Exception as exc:
        logger.debug(
            "Skipping CDS feature %s on %s: translation failed: %s",
            getattr(feature, "id", ""),
            record.id,
            exc,
        )
        return None
    return _clean_protein_sequence(protein)


def _cds_span(feature: SeqFeature) -> tuple[int, int, int | None] | None:
    location = feature.location
    if location is None:
        return None
    start = int(location.start)
    end = int(location.end)
    if start < 0 or end <= start:
        return None
    strand = location.strand if location.strand in {-1, 1} else None
    return start, end, strand


def extract_cds_proteins(
    records: Sequence[SeqRecord],
    *,
    record_index_offset: int = 0,
    prefer_source_ids: bool = True,
) -> ProteinExtractionResult:
    """Extract CDS proteins using protein IDs where possible and genomic spans.

    Coordinates in the returned metadata are Python/Biopython-style 0-based start
    and 0-based exclusive end values.
    """

    candidate_records: list[list[_CdsProteinCandidate]] = []
    source_id_counts: dict[str, int] = {}
    protein_map: dict[str, CdsProtein] = {}

    record_index_offset = int(record_index_offset)

    for record_index, record in enumerate(records):
        global_record_index = record_index + record_index_offset
        record_candidates: list[_CdsProteinCandidate] = []
        cds_count = 0
        for feature_index, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue

            span = _cds_span(feature)
            if span is None:
                continue
            start, end, strand = span

            protein_sequence = _clean_protein_sequence(
                _first_qualifier(feature, "translation")
            )
            if protein_sequence is None:
                protein_sequence = _translate_cds_feature(record, feature)
            if protein_sequence is None:
                continue

            cds_count += 1
            synthetic_protein_id = f"gbd_r{global_record_index + 1:04d}_cds{cds_count:06d}"
            source_protein_id = _first_qualifier(feature, "protein_id")
            fasta_safe_source_id = _fasta_safe_protein_id(source_protein_id)
            if fasta_safe_source_id is not None:
                source_id_counts[fasta_safe_source_id] = source_id_counts.get(fasta_safe_source_id, 0) + 1
            gene = _first_qualifier(feature, "gene")
            product = _first_qualifier(feature, "product")
            note = _first_qualifier(feature, "note")
            label = (
                _first_qualifier_from_any(
                    feature, ("locus_tag", "protein_id", "gene", "product")
                )
                or synthetic_protein_id
            )
            record_candidates.append(
                _CdsProteinCandidate(
                    protein_id=synthetic_protein_id,
                    source_protein_id=source_protein_id,
                    record_index=global_record_index,
                    feature_index=feature_index,
                    record_id=str(record.id),
                    start=start,
                    end=end,
                    strand=strand,
                    label=label,
                    protein_length=len(protein_sequence),
                    sequence=protein_sequence,
                    feature_svg_id=compute_feature_hash(feature, record_id=record.id),
                    gene=gene,
                    product=product,
                    note=note,
                )
            )

        candidate_records.append(record_candidates)

    proteins_by_record: list[list[CdsProtein]] = []
    used_ids: set[str] = set()
    for record_candidates in candidate_records:
        record_proteins: list[CdsProtein] = []
        for candidate in record_candidates:
            fasta_safe_source_id = _fasta_safe_protein_id(candidate.source_protein_id)
            preferred_id = (
                fasta_safe_source_id
                if prefer_source_ids
                and fasta_safe_source_id is not None
                and source_id_counts.get(fasta_safe_source_id, 0) == 1
                else candidate.protein_id
            )
            protein_id = _unique_protein_id(preferred_id, candidate.protein_id, used_ids)
            cds_protein = CdsProtein(
                protein_id=protein_id,
                record_index=candidate.record_index,
                feature_index=candidate.feature_index,
                record_id=candidate.record_id,
                start=candidate.start,
                end=candidate.end,
                strand=candidate.strand,
                label=candidate.label,
                protein_length=candidate.protein_length,
                sequence=candidate.sequence,
                source_protein_id=candidate.source_protein_id,
                feature_svg_id=candidate.feature_svg_id,
                gene=candidate.gene,
                product=candidate.product,
                note=candidate.note,
            )
            record_proteins.append(cds_protein)
            protein_map[protein_id] = cds_protein
            used_ids.add(protein_id)

        proteins_by_record.append(record_proteins)

    return ProteinExtractionResult(
        proteins_by_record=proteins_by_record,
        protein_map=protein_map,
    )


def proteins_to_fasta(proteins: Sequence[CdsProtein]) -> str:
    """Return FASTA text for CDS proteins."""

    lines: list[str] = []
    for protein in proteins:
        lines.append(f">{protein.protein_id} {protein.label}")
        sequence = protein.sequence
        lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    return "\n".join(lines) + ("\n" if lines else "")


def _coerce_outfmt6_numeric_columns(df: DataFrame) -> DataFrame:
    coerced = df.copy()
    for column in _NUMERIC_COMPARISON_COLUMNS:
        coerced[column] = pd.to_numeric(coerced[column], errors="coerce")
    if coerced[list(_NUMERIC_COMPARISON_COLUMNS)].isna().any(axis=None):
        raise ParseError("LOSATP blastp output contains non-numeric outfmt 6 fields.")
    return coerced


def parse_losatp_outfmt6(text: str) -> DataFrame:
    """Parse LOSATP blastp outfmt 6 text into the standard comparison columns."""

    data_lines = [
        line
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not data_lines:
        return pd.DataFrame(columns=COMPARISON_COLUMNS)

    try:
        df = pd.read_csv(
            StringIO("\n".join(data_lines)),
            sep="\t",
            names=COMPARISON_COLUMNS,
        )
    except Exception as exc:
        raise ParseError(f"Failed to parse LOSATP blastp output: {exc}") from exc

    if len(df.columns) != len(COMPARISON_COLUMNS):
        raise ParseError("LOSATP blastp output does not match outfmt 6 columns.")
    return _coerce_outfmt6_numeric_columns(df)


def _validate_max_hits(max_hits: int, *, option_name: str = "protein_blastp_max_hits") -> None:
    if int(max_hits) <= 0:
        raise ValidationError(f"{option_name} must be > 0")


def _validate_candidate_limit(candidate_limit: int | None) -> None:
    if candidate_limit is not None and int(candidate_limit) <= 0:
        raise ValidationError("protein_blastp_candidate_limit must be > 0 or None")


def normalize_protein_blastp_mode(mode: str | None) -> ProteinBlastpMode:
    """Return a validated LOSATP blastp mode."""

    normalized = str(mode or "none").strip().lower()
    if normalized not in PROTEIN_BLASTP_MODES:
        raise ValidationError(
            "protein_blastp_mode must be one of: "
            + ", ".join(PROTEIN_BLASTP_MODES)
        )
    return normalized  # type: ignore[return-value]


def _validate_comparison_columns(hits: DataFrame) -> None:
    missing_columns = set(COMPARISON_COLUMNS).difference(hits.columns)
    if missing_columns:
        raise ParseError(
            "LOSATP blastp output is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )


def filter_protein_hits_by_thresholds(
    hits: DataFrame,
    *,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
) -> DataFrame:
    """Keep LOSATP blastp hits that pass visible pairwise-match thresholds."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    try:
        evalue_threshold = float(evalue)
        bitscore_threshold = float(bitscore)
        identity_threshold = float(identity)
        alignment_length_threshold = int(alignment_length)
    except (TypeError, ValueError) as exc:
        raise ValidationError("protein colinearity thresholds must be numeric") from exc
    if alignment_length_threshold < 0:
        raise ValidationError("alignment_length must be >= 0")

    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    return coerced_hits.loc[
        (coerced_hits["evalue"] <= evalue_threshold)
        & (coerced_hits["bitscore"] >= bitscore_threshold)
        & (coerced_hits["identity"] >= identity_threshold)
        & (coerced_hits["alignment_length"] >= alignment_length_threshold)
    ].reset_index(drop=True)


def _sort_hits_for_query_best(hits: DataFrame) -> DataFrame:
    return hits.sort_values(
        ["query", "bitscore", "evalue", "identity", "alignment_length", "subject"],
        ascending=[True, False, True, False, False, True],
        kind="mergesort",
    )


def _sort_hits_for_subject_best(hits: DataFrame) -> DataFrame:
    return hits.sort_values(
        ["subject", "bitscore", "evalue", "identity", "alignment_length", "query"],
        ascending=[True, False, True, False, False, True],
        kind="mergesort",
    )


def select_top_hits_per_query(
    hits: DataFrame,
    *,
    max_hits: int,
) -> DataFrame:
    """Keep the strongest distinct subject protein hits per query protein."""

    _validate_max_hits(max_hits)
    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    sorted_hits = _sort_hits_for_query_best(_coerce_outfmt6_numeric_columns(hits))
    sorted_hits = sorted_hits.drop_duplicates(["query", "subject"], keep="first")
    return (
        sorted_hits.groupby("query", group_keys=False, sort=False)
        .head(int(max_hits))
        .reset_index(drop=True)
    )


def select_reciprocal_best_hits(hits: DataFrame) -> DataFrame:
    """Keep query-subject pairs that are mutual best hits in one adjacent-pair table."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    query_best = (
        _sort_hits_for_query_best(coerced_hits)
        .drop_duplicates(["query", "subject"], keep="first")
        .drop_duplicates("query", keep="first")
    )
    subject_best = (
        _sort_hits_for_subject_best(coerced_hits)
        .drop_duplicates(["query", "subject"], keep="first")
        .drop_duplicates("subject", keep="first")
    )
    reciprocal_pairs = {
        (str(row.query), str(row.subject))
        for row in subject_best.itertuples(index=False)
    }
    return query_best.loc[
        [
            (str(row.query), str(row.subject)) in reciprocal_pairs
            for row in query_best.itertuples(index=False)
        ]
    ].reset_index(drop=True)


def select_best_hits_per_query(
    hits: DataFrame,
) -> DataFrame:
    """Choose the deterministic best subject for each query protein."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    return (
        _sort_hits_for_query_best(_coerce_outfmt6_numeric_columns(hits))
        .drop_duplicates(["query", "subject"], keep="first")
        .drop_duplicates("query", keep="first")
        .reset_index(drop=True)
    )


def select_reciprocal_best_hit_edges(
    forward_hits: DataFrame,
    reverse_hits: DataFrame,
) -> DataFrame:
    """Keep RBH edges from directional best-hit tables.

    Returned rows use the forward orientation: query record i, subject record j.
    """

    if forward_hits.empty or reverse_hits.empty:
        return forward_hits.iloc[0:0].copy()

    forward_best = select_best_hits_per_query(forward_hits)
    reverse_best = select_best_hits_per_query(reverse_hits)
    reverse_best_by_query = {
        str(row.query): str(row.subject)
        for row in reverse_best.itertuples(index=False)
    }
    keep_mask = [
        reverse_best_by_query.get(str(row.subject)) == str(row.query)
        for row in forward_best.itertuples(index=False)
    ]
    return forward_best.loc[keep_mask].reset_index(drop=True)


def cap_hits_per_query(
    hits: DataFrame,
    *,
    max_hits: int = 5,
    distinct_subjects: bool = True,
) -> DataFrame:
    """Keep the strongest LOSATP hits per query protein."""

    _validate_max_hits(max_hits)
    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    sorted_hits = _sort_hits_for_query_best(_coerce_outfmt6_numeric_columns(hits))
    if distinct_subjects:
        sorted_hits = sorted_hits.drop_duplicates(["query", "subject"], keep="first")
    return (
        sorted_hits.groupby("query", group_keys=False, sort=False)
        .head(int(max_hits))
        .reset_index(drop=True)
    )


class _UnionFind:
    def __init__(self) -> None:
        self._parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        if item not in self._parent:
            self._parent[item] = item

    def find(self, item: str) -> str:
        self.add(item)
        parent = self._parent[item]
        if parent != item:
            parent = self.find(parent)
            self._parent[item] = parent
        return parent

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if right_root < left_root:
            left_root, right_root = right_root, left_root
        self._parent[right_root] = left_root


def _row_float(row: object, column: str, default: float = 0.0) -> float:
    try:
        return float(getattr(row, column))
    except (TypeError, ValueError, AttributeError):
        return float(default)


def _member_rank_from_row(row: object) -> tuple[float, float, float, float]:
    return (
        _row_float(row, "bitscore", 0.0),
        _row_float(row, "evalue", float("inf")),
        _row_float(row, "identity", 0.0),
        _row_float(row, "alignment_length", 0.0),
    )


def _is_better_member_rank(
    candidate: tuple[float, float, float, float],
    current: tuple[float, float, float, float] | None,
) -> bool:
    if current is None:
        return True
    return (
        candidate[0] > current[0]
        or (candidate[0] == current[0] and candidate[1] < current[1])
        or (
            candidate[0] == current[0]
            and candidate[1] == current[1]
            and candidate[2] > current[2]
        )
        or (
            candidate[0] == current[0]
            and candidate[1] == current[1]
            and candidate[2] == current[2]
            and candidate[3] > current[3]
        )
    )


def _protein_sort_key(protein: CdsProtein) -> tuple[int, int, int, str]:
    return (
        int(protein.record_index),
        int(protein.start),
        int(protein.end),
        str(protein.protein_id),
    )


_ORTHOGROUP_NAME_SOURCE_WEIGHTS = {
    "product": 80,
    "gene": 55,
    "note": 35,
    "label": 10,
}
_ORTHOGROUP_NAME_SOURCE_ORDER = ("product", "gene", "note", "label")
_ORTHOGROUP_NAME_SOURCE_RANK = {
    source: index for index, source in enumerate(_ORTHOGROUP_NAME_SOURCE_ORDER)
}
_NOTE_PREFIX_RE = re.compile(r"^(?:product|gene|note)\s*:\s*", re.IGNORECASE)
_DUF_CANDIDATE_RE = re.compile(r"\bDUF\d+\b.*\bdomain-containing protein\b", re.IGNORECASE)
_GENERIC_ORTHOGROUP_NAME_CANDIDATES = {
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
    "predicted protein",
    "putative protein",
}
_WEAK_ORTHOGROUP_NAME_CANDIDATES = {
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
    "predicted protein",
    "putative protein",
}


def _normalize_orthogroup_name_text(
    text: str | None,
    source: str,
) -> tuple[str, str] | None:
    if text is None:
        return None
    normalized = re.sub(r"\s+", " ", str(text)).strip()
    normalized = _NOTE_PREFIX_RE.sub("", normalized).strip()
    normalized = normalized.strip(" \t\r\n.,;:()[]{}")
    if not normalized:
        return None
    compact_length = len(re.sub(r"[^0-9A-Za-z]+", "", normalized))
    if source != "gene" and compact_length < 4:
        return None
    key = normalized.casefold()
    if key in _GENERIC_ORTHOGROUP_NAME_CANDIDATES:
        return None
    return key, normalized


def _is_weak_orthogroup_name_candidate(text: str, source: str) -> bool:
    normalized = re.sub(r"\s+", " ", str(text)).strip().casefold()
    if source == "label":
        return True
    if normalized in _WEAK_ORTHOGROUP_NAME_CANDIDATES:
        return True
    return bool(_DUF_CANDIDATE_RE.search(text))


def _member_annotation_value(member: OrthogroupMember, source: str) -> str | None:
    if source == "product":
        return member.product
    if source == "gene":
        return member.gene
    if source == "note":
        return member.note
    if source == "label":
        return member.label
    return None


def _iter_orthogroup_name_candidates(
    members: Sequence[OrthogroupMember],
    *,
    include_label: bool,
) -> list[OrthogroupNameCandidate]:
    candidate_sources = _ORTHOGROUP_NAME_SOURCE_ORDER if include_label else _ORTHOGROUP_NAME_SOURCE_ORDER[:-1]
    accumulators: dict[str, dict[str, object]] = {}

    for member in members:
        member_id = str(member.protein_id)
        for source in candidate_sources:
            normalized = _normalize_orthogroup_name_text(
                _member_annotation_value(member, source),
                source,
            )
            if normalized is None:
                continue
            key, display_text = normalized
            source_weight = _ORTHOGROUP_NAME_SOURCE_WEIGHTS[source]
            source_rank = _ORTHOGROUP_NAME_SOURCE_RANK[source]
            accumulator = accumulators.setdefault(
                key,
                {
                    "text": display_text,
                    "source": source,
                    "source_weight": source_weight,
                    "source_rank": source_rank,
                    "members": set(),
                    "records": set(),
                    "representatives": set(),
                    "weak": _is_weak_orthogroup_name_candidate(display_text, source),
                },
            )
            if (
                source_weight > int(accumulator["source_weight"])
                or (
                    source_weight == int(accumulator["source_weight"])
                    and source_rank < int(accumulator["source_rank"])
                )
            ):
                accumulator["text"] = display_text
                accumulator["source"] = source
                accumulator["source_weight"] = source_weight
                accumulator["source_rank"] = source_rank
            if not _is_weak_orthogroup_name_candidate(display_text, source):
                accumulator["weak"] = False
            members_set = accumulator["members"]
            records_set = accumulator["records"]
            representatives_set = accumulator["representatives"]
            if isinstance(members_set, set):
                members_set.add(member_id)
            if isinstance(records_set, set):
                records_set.add(int(member.record_index))
            if member.representative and isinstance(representatives_set, set):
                representatives_set.add(member_id)

    candidates: list[OrthogroupNameCandidate] = []
    for accumulator in accumulators.values():
        members_set = accumulator["members"]
        records_set = accumulator["records"]
        representatives_set = accumulator["representatives"]
        member_count = len(members_set) if isinstance(members_set, set) else 0
        record_coverage_count = len(records_set) if isinstance(records_set, set) else 0
        representative_count = len(representatives_set) if isinstance(representatives_set, set) else 0
        source_weight = int(accumulator["source_weight"])
        weak_penalty = 45 if bool(accumulator["weak"]) else 0
        score = (
            source_weight
            + member_count * 10
            + record_coverage_count * 20
            + representative_count * 5
            - weak_penalty
        )
        candidates.append(
            OrthogroupNameCandidate(
                text=str(accumulator["text"]),
                source=str(accumulator["source"]),
                member_count=member_count,
                record_coverage_count=record_coverage_count,
                representative_count=representative_count,
                score=float(score),
            )
        )

    return sorted(
        candidates,
        key=lambda candidate: (
            -candidate.score,
            _ORTHOGROUP_NAME_SOURCE_RANK.get(candidate.source, 99),
            -candidate.record_coverage_count,
            -candidate.member_count,
            candidate.text.casefold(),
        ),
    )


def _orthogroup_name_confidence(
    candidate: OrthogroupNameCandidate | None,
) -> str:
    if candidate is None:
        return "none"
    if _is_weak_orthogroup_name_candidate(candidate.text, candidate.source):
        return "low"
    if candidate.member_count <= 1:
        return "low"
    if candidate.record_coverage_count >= 2:
        return "high"
    return "medium"


def _orthogroup_description(
    candidate: OrthogroupNameCandidate | None,
    confidence: str,
    total_record_coverage: int,
) -> str:
    if candidate is None or confidence in {"none", "low"}:
        return "No informative product/gene/note consensus was found."
    record_word = "record" if total_record_coverage == 1 else "records"
    return (
        f"Suggested from {candidate.source} annotations in "
        f"{candidate.record_coverage_count} of {total_record_coverage} {record_word}."
    )


def _build_orthogroup_name_metadata(
    orthogroups: Mapping[str, Sequence[OrthogroupMember]],
) -> tuple[
    dict[str, str],
    dict[str, str],
    dict[str, list[OrthogroupNameCandidate]],
    dict[str, str],
]:
    names_by_orthogroup_id: dict[str, str] = {}
    descriptions_by_orthogroup_id: dict[str, str] = {}
    name_candidates_by_orthogroup_id: dict[str, list[OrthogroupNameCandidate]] = {}
    confidence_by_orthogroup_id: dict[str, str] = {}

    for orthogroup_id, members in orthogroups.items():
        annotation_candidates = _iter_orthogroup_name_candidates(members, include_label=False)
        candidates = annotation_candidates or _iter_orthogroup_name_candidates(members, include_label=True)
        top_candidate = candidates[0] if candidates else None
        total_record_coverage = len({int(member.record_index) for member in members})
        confidence = _orthogroup_name_confidence(top_candidate)

        if top_candidate is not None:
            names_by_orthogroup_id[orthogroup_id] = top_candidate.text
        descriptions_by_orthogroup_id[orthogroup_id] = _orthogroup_description(
            top_candidate,
            confidence,
            total_record_coverage,
        )
        name_candidates_by_orthogroup_id[orthogroup_id] = candidates
        confidence_by_orthogroup_id[orthogroup_id] = confidence

    return (
        names_by_orthogroup_id,
        descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id,
    )


def build_orthogroups_from_protein_hits(
    hits_by_pair: Sequence[DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> OrthogroupResult:
    """Build connected-component orthogroups from retained LOSATP protein hits."""

    union_find = _UnionFind()
    edge_nodes: list[tuple[str, str]] = []
    member_ranks: dict[str, tuple[float, float, float, float]] = {}
    missing_ids: set[str] = set()

    for hits in hits_by_pair:
        if hits is None or hits.empty:
            continue
        missing_columns = {"query", "subject"}.difference(hits.columns)
        if missing_columns:
            raise ParseError(
                "LOSATP blastp output is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )
        for row in hits.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            if query_id not in protein_map:
                missing_ids.add(query_id)
            if subject_id not in protein_map:
                missing_ids.add(subject_id)
            if query_id not in protein_map or subject_id not in protein_map:
                continue
            if protein_map[query_id].record_index == protein_map[subject_id].record_index:
                continue
            union_find.union(query_id, subject_id)
            edge_nodes.append((query_id, subject_id))
            rank = _member_rank_from_row(row)
            if _is_better_member_rank(rank, member_ranks.get(query_id)):
                member_ranks[query_id] = rank
            if _is_better_member_rank(rank, member_ranks.get(subject_id)):
                member_ranks[subject_id] = rank

    if missing_ids:
        raise ParseError(
            "LOSATP blastp output contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )

    components: dict[str, set[str]] = {}
    for query_id, subject_id in edge_nodes:
        components.setdefault(union_find.find(query_id), set()).update({query_id, subject_id})

    sorted_components = sorted(
        components.values(),
        key=lambda member_ids: min(_protein_sort_key(protein_map[member_id]) for member_id in member_ids),
    )

    orthogroups: dict[str, list[OrthogroupMember]] = {}
    member_by_protein_id: dict[str, OrthogroupMember] = {}
    for component_index, member_ids in enumerate(sorted_components, start=1):
        orthogroup_id = f"og_{component_index}"
        member_ids_sorted = sorted(member_ids, key=lambda member_id: _protein_sort_key(protein_map[member_id]))

        representative_ids: set[str] = set()
        members_by_record: dict[int, list[str]] = {}
        for member_id in member_ids_sorted:
            protein = protein_map[member_id]
            members_by_record.setdefault(int(protein.record_index), []).append(member_id)

        for record_member_ids in members_by_record.values():
            representative_id = min(
                record_member_ids,
                key=lambda member_id: (
                    -member_ranks.get(member_id, (0.0, float("inf"), 0.0, 0.0))[0],
                    member_ranks.get(member_id, (0.0, float("inf"), 0.0, 0.0))[1],
                    -member_ranks.get(member_id, (0.0, float("inf"), 0.0, 0.0))[2],
                    -member_ranks.get(member_id, (0.0, float("inf"), 0.0, 0.0))[3],
                    str(member_id),
                ),
            )
            representative_ids.add(representative_id)

        group_members: list[OrthogroupMember] = []
        for member_id in member_ids_sorted:
            protein = protein_map[member_id]
            member = OrthogroupMember(
                orthogroup_id=orthogroup_id,
                protein_id=protein.protein_id,
                record_index=protein.record_index,
                feature_index=protein.feature_index,
                record_id=protein.record_id,
                label=protein.label,
                start=protein.start,
                end=protein.end,
                strand=protein.strand,
                feature_svg_id=protein.feature_svg_id,
                source_protein_id=protein.source_protein_id,
                gene=protein.gene,
                product=protein.product,
                note=protein.note,
                representative=member_id in representative_ids,
            )
            group_members.append(member)
            member_by_protein_id[member_id] = member
        orthogroups[orthogroup_id] = group_members

    (
        names_by_orthogroup_id,
        descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id,
    ) = _build_orthogroup_name_metadata(orthogroups)

    return OrthogroupResult(
        orthogroups=orthogroups,
        member_by_protein_id=member_by_protein_id,
        names_by_orthogroup_id=names_by_orthogroup_id,
        descriptions_by_orthogroup_id=descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id=name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id=confidence_by_orthogroup_id,
    )


def _genomic_link_coordinates(protein: CdsProtein) -> tuple[int, int]:
    if protein.strand == -1:
        return protein.end, protein.start + 1
    return protein.start + 1, protein.end


def convert_protein_hits_to_genomic_links(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
    orthogroups: OrthogroupResult | None = None,
) -> DataFrame:
    """Convert protein hit rows to genomic-coordinate comparison rows."""

    return convert_pair_protein_hits_to_genomic_links(
        hits,
        protein_map,
        protein_map,
        orthogroups=orthogroups,
    )


def _feature_svg_id(protein: CdsProtein) -> str:
    return str(protein.feature_svg_id or "")


def _protein_metadata_value(value: object | None) -> object:
    return "" if value is None else value


def convert_pair_protein_hits_to_genomic_links(
    hits: DataFrame,
    query_protein_map: Mapping[str, CdsProtein],
    subject_protein_map: Mapping[str, CdsProtein],
    orthogroups: OrthogroupResult | None = None,
) -> DataFrame:
    """Convert pairwise protein hit rows using separate query and subject maps."""

    if hits.empty:
        return pd.DataFrame(columns=LOSATP_COMPARISON_COLUMNS)

    rows: list[dict[str, object]] = []
    missing_ids: set[str] = set()
    for row in hits.itertuples(index=False):
        query_id = str(row.query)
        subject_id = str(row.subject)
        query_protein = query_protein_map.get(query_id)
        subject_protein = subject_protein_map.get(subject_id)
        if query_protein is None:
            missing_ids.add(query_id)
        if subject_protein is None:
            missing_ids.add(subject_id)
        if query_protein is None or subject_protein is None:
            continue

        qstart, qend = _genomic_link_coordinates(query_protein)
        sstart, send = _genomic_link_coordinates(subject_protein)
        query_member = (
            orthogroups.member_by_protein_id.get(query_id)
            if orthogroups is not None
            else None
        )
        subject_member = (
            orthogroups.member_by_protein_id.get(subject_id)
            if orthogroups is not None
            else None
        )
        orthogroup_id = ""
        if query_member is not None and subject_member is not None:
            if query_member.orthogroup_id == subject_member.orthogroup_id:
                orthogroup_id = query_member.orthogroup_id
        rows.append(
            {
                "query": query_protein.record_id,
                "subject": subject_protein.record_id,
                "identity": row.identity,
                "alignment_length": row.alignment_length,
                "mismatches": row.mismatches,
                "gap_opens": row.gap_opens,
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "evalue": row.evalue,
                "bitscore": row.bitscore,
                "query_protein_id": query_protein.protein_id,
                "subject_protein_id": subject_protein.protein_id,
                "query_source_protein_id": _protein_metadata_value(query_protein.source_protein_id),
                "subject_source_protein_id": _protein_metadata_value(subject_protein.source_protein_id),
                "query_record_index": query_protein.record_index,
                "subject_record_index": subject_protein.record_index,
                "query_feature_index": query_protein.feature_index,
                "subject_feature_index": subject_protein.feature_index,
                "query_feature_svg_id": _feature_svg_id(query_protein),
                "subject_feature_svg_id": _feature_svg_id(subject_protein),
                "orthogroup_id": orthogroup_id,
                "query_orthogroup_representative": (
                    bool(query_member.representative) if query_member is not None else False
                ),
                "subject_orthogroup_representative": (
                    bool(subject_member.representative) if subject_member is not None else False
                ),
            }
        )

    if missing_ids:
        raise ParseError(
            "LOSATP blastp output contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )

    return pd.DataFrame.from_records(rows, columns=LOSATP_COMPARISON_COLUMNS)


def run_losatp_blastp(
    query_fasta: str,
    subject_fasta: str,
    *,
    losatp_bin: str = "losat",
    max_hits: int | None = None,
) -> DataFrame:
    """Run external LOSATP blastp and parse outfmt 6 output."""

    if max_hits is not None:
        _validate_max_hits(max_hits, option_name="protein_blastp_candidate_limit")
    with tempfile.TemporaryDirectory(prefix="gbdraw_losatp_") as temp_dir:
        temp_path = Path(temp_dir)
        query_path = temp_path / "query.faa"
        subject_path = temp_path / "subject.faa"
        query_path.write_text(query_fasta, encoding="utf-8")
        subject_path.write_text(subject_fasta, encoding="utf-8")

        command = [
            losatp_bin,
            "blastp",
            "-query",
            str(query_path),
            "-subject",
            str(subject_path),
            "-outfmt",
            "6",
            "-max_hsps_per_subject",
            "1",
        ]
        if max_hits is not None:
            command.extend(["-max_target_seqs", str(int(max_hits))])
        logger.info("INFO: Running LOSATP blastp for protein colinearity.")
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError as exc:
            raise ValidationError(f"LOSATP binary not found: {losatp_bin}") from exc

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        detail = f": {stderr}" if stderr else ""
        raise ValidationError(f"LOSATP blastp failed with exit code {completed.returncode}{detail}")
    return parse_losatp_outfmt6(completed.stdout)


def _validate_extraction_has_proteins(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    *,
    option_name: str,
) -> None:
    empty_record_ids = [
        str(records[index].id)
        for index, proteins in enumerate(extraction.proteins_by_record)
        if not proteins
    ]
    if empty_record_ids:
        raise ValidationError(
            f"{option_name} requires at least one CDS protein in each record; "
            "no CDS proteins were found in: "
            + ", ".join(empty_record_ids)
        )


def _run_losatp_search(
    query_fasta: str,
    subject_fasta: str,
    *,
    losatp_bin: str,
    candidate_limit: int | None,
    runner: LosatpRunner | None,
) -> DataFrame:
    if runner is not None:
        return runner(query_fasta, subject_fasta)
    return run_losatp_blastp(
        query_fasta,
        subject_fasta,
        losatp_bin=losatp_bin,
        max_hits=candidate_limit,
    )


def build_pairwise_protein_blastp_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    max_hits: int = 5,
    candidate_limit: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
) -> ProteinBlastpResult:
    """Generate adjacent-record LOSATP blastp display comparisons."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='pairwise' requires at least two records")
    _validate_max_hits(max_hits)
    _validate_candidate_limit(candidate_limit)
    if int(alignment_length) < 0:
        raise ValidationError("alignment_length must be >= 0")

    extraction = extract_cds_proteins(records)
    _validate_extraction_has_proteins(
        records,
        extraction,
        option_name="protein_blastp_mode='pairwise'",
    )

    comparisons: list[DataFrame] = []
    for record_index in range(len(records) - 1):
        query_fasta = proteins_to_fasta(extraction.proteins_by_record[record_index])
        subject_fasta = proteins_to_fasta(extraction.proteins_by_record[record_index + 1])
        protein_hits = _run_losatp_search(
            query_fasta,
            subject_fasta,
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
        display_hits = select_top_hits_per_query(filtered_hits, max_hits=max_hits)
        comparisons.append(
            convert_protein_hits_to_genomic_links(
                display_hits,
                extraction.protein_map,
                orthogroups=None,
            )
        )
    return ProteinBlastpResult(comparisons=comparisons, orthogroups=None)


def build_rbh_orthogroup_protein_blastp_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    candidate_limit: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
) -> ProteinBlastpResult:
    """Infer all-vs-all RBH orthogroups and return adjacent RBH display links."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='orthogroup' requires at least two records")
    _validate_candidate_limit(candidate_limit)
    if int(alignment_length) < 0:
        raise ValidationError("alignment_length must be >= 0")

    extraction = extract_cds_proteins(records)
    _validate_extraction_has_proteins(
        records,
        extraction,
        option_name="protein_blastp_mode='orthogroup'",
    )

    adjacent_rbh_edges_by_pair: list[DataFrame] = [
        pd.DataFrame(columns=COMPARISON_COLUMNS)
        for _ in range(len(records) - 1)
    ]
    orthogroup_edges: list[DataFrame] = []

    for query_index in range(len(records)):
        for subject_index in range(query_index + 1, len(records)):
            query_fasta = proteins_to_fasta(extraction.proteins_by_record[query_index])
            subject_fasta = proteins_to_fasta(extraction.proteins_by_record[subject_index])

            forward_hits = _run_losatp_search(
                query_fasta,
                subject_fasta,
                losatp_bin=losatp_bin,
                candidate_limit=candidate_limit,
                runner=runner,
            )
            reverse_hits = _run_losatp_search(
                subject_fasta,
                query_fasta,
                losatp_bin=losatp_bin,
                candidate_limit=candidate_limit,
                runner=runner,
            )
            filtered_forward_hits = filter_protein_hits_by_thresholds(
                forward_hits,
                evalue=evalue,
                bitscore=bitscore,
                identity=identity,
                alignment_length=alignment_length,
            )
            filtered_reverse_hits = filter_protein_hits_by_thresholds(
                reverse_hits,
                evalue=evalue,
                bitscore=bitscore,
                identity=identity,
                alignment_length=alignment_length,
            )
            rbh_edges = select_reciprocal_best_hit_edges(
                filtered_forward_hits,
                filtered_reverse_hits,
            )
            orthogroup_edges.append(rbh_edges)
            if subject_index == query_index + 1:
                adjacent_rbh_edges_by_pair[query_index] = rbh_edges

    orthogroups = build_orthogroups_from_protein_hits(
        orthogroup_edges,
        extraction.protein_map,
    )
    comparisons = [
        convert_protein_hits_to_genomic_links(
            rbh_edges,
            extraction.protein_map,
            orthogroups=orthogroups,
        )
        for rbh_edges in adjacent_rbh_edges_by_pair
    ]
    return ProteinBlastpResult(comparisons=comparisons, orthogroups=orthogroups)


def build_protein_colinearity_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    max_hits: int = 5,
    candidate_limit: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
) -> list[DataFrame]:
    """Compatibility helper returning pairwise LOSATP blastp display comparisons."""

    return build_pairwise_protein_blastp_comparisons(
        records,
        losatp_bin=losatp_bin,
        max_hits=max_hits,
        candidate_limit=candidate_limit,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        runner=runner,
    ).comparisons


__all__ = [
    "CdsProtein",
    "LOSATP_COMPARISON_COLUMNS",
    "LOSATP_METADATA_COLUMNS",
    "PROTEIN_BLASTP_MODES",
    "OrthogroupMember",
    "OrthogroupNameCandidate",
    "OrthogroupResult",
    "ProteinBlastpMode",
    "ProteinBlastpResult",
    "ProteinExtractionResult",
    "build_orthogroups_from_protein_hits",
    "build_pairwise_protein_blastp_comparisons",
    "build_protein_colinearity_comparisons",
    "build_rbh_orthogroup_protein_blastp_comparisons",
    "cap_hits_per_query",
    "convert_pair_protein_hits_to_genomic_links",
    "convert_protein_hits_to_genomic_links",
    "extract_cds_proteins",
    "filter_protein_hits_by_thresholds",
    "normalize_protein_blastp_mode",
    "parse_losatp_outfmt6",
    "proteins_to_fasta",
    "run_losatp_blastp",
    "select_best_hits_per_query",
    "select_reciprocal_best_hit_edges",
    "select_reciprocal_best_hits",
    "select_top_hits_per_query",
]
