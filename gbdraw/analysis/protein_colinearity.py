#!/usr/bin/env python
# coding: utf-8

"""Protein-level colinearity helpers for linear pairwise comparisons."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from io import StringIO
import logging
import math
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
    "rbh_orthogroup_id",
    "ortholog_path_id",
    "edge_kind",
    "render_role",
    "query_orthogroup_representative",
    "subject_orthogroup_representative",
    "query_orthogroup_member_count",
    "subject_orthogroup_member_count",
    "query_orthogroup_role",
    "subject_orthogroup_role",
    "query_orthogroup_confidence",
    "subject_orthogroup_confidence",
    "query_orthogroup_assignment_reason",
    "subject_orthogroup_assignment_reason",
)
LOSATP_COMPARISON_COLUMNS = tuple(COMPARISON_COLUMNS) + LOSATP_METADATA_COLUMNS
PROTEIN_BLASTP_MODES = ("none", "pairwise", "orthogroup", "collinear")
ProteinBlastpMode = Literal["none", "pairwise", "orthogroup", "collinear"]
ORTHOGROUP_INFERENCE_VERSION = "anchor_core_v1"
ORTHOGROUP_MEMBERSHIP_MODES = (ORTHOGROUP_INFERENCE_VERSION,)
_LEGACY_ORTHOGROUP_MEMBERSHIP_MODES = ("rbh", "family_merge", "distribution_split")
OrthogroupMembershipMode = Literal["anchor_core_v1"]
OrthogroupMemberRole = Literal["anchor", "coortholog", "inparalog", "low_confidence"]
OrthogroupMemberConfidence = Literal["high", "medium", "low"]
OrthologEdgeKind = Literal[
    "rbh",
    "coortholog",
    "same_record_inparalog",
    "related_homolog",
    "ambiguous_paralog",
    "weak_bridge",
    "domain_only",
    "related_paralog",
]
OrthologRenderRole = Literal["block_anchor", "display_edge"]
ORTHOLOG_DISPLAY_EDGE_KIND_RANK = {
    "rbh": 0,
    "coortholog": 1,
    "same_record_inparalog": 2,
    "related_homolog": 3,
    "ambiguous_paralog": 4,
    "weak_bridge": 5,
    "domain_only": 6,
    "related_paralog": 7,
}
_NEAR_RECIPROCAL_MIN_RATIO = 0.85
_CORE_BRIDGE_MIN_RATIO = 0.85
_INPARALOG_MIN_SAME_RECORD_RATIO_TO_LOCAL_ANCHOR = 0.75
_MIN_BEST_SECOND_CORE_RATIO = 1.25
_MIN_MEMBERSHIP_MIN_COVERAGE = 0.30
_DOMAIN_ONLY_MAX_MIN_COVERAGE = 0.25
_LOW_CONFIDENCE_MIN_BEST_SECOND_CORE_RATIO = 1.10
_ORTHOGROUP_SPLIT_MIN_SEPARATION_RATIO = 1.35
_ORTHOGROUP_SPLIT_MAX_CONDUCTANCE = 0.20
_ORTHOGROUP_SPLIT_MIN_STABILITY = 2
_ORTHOGROUP_SPLIT_MIN_CHILD_SIZE = 2
_ORTHOGROUP_SPLIT_TOP_FRACTION = 0.05
_ORTHOGROUP_SPLIT_MIN_COVERAGE = 0.30
_ORTHOGROUP_SPLIT_MIN_FIT_HITS = 12
_ORTHOGROUP_SPLIT_EPSILON = 1e-12


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
    locus_tag: str | None = None
    gene_id: str | None = None
    old_locus_tag: str | None = None
    db_xref: tuple[str, ...] = ()
    gff_id: str | None = None
    parent_ids: tuple[str, ...] = ()
    gene_parent_id: str | None = None
    feature_type: str = "CDS"
    feature_hash_start: int | None = None
    feature_hash_end: int | None = None
    feature_hash_strand: int | None = None


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
    role: OrthogroupMemberRole = "anchor"
    confidence: OrthogroupMemberConfidence = "high"
    assignment_reason: str = ""
    supporting_edges: tuple[str, ...] = ()
    best_core_support: float = 0.0
    second_best_core_support: float = 0.0


@dataclass(frozen=True)
class OrthologEdge:
    """Selected evidence edge describing relationships inside or near an orthogroup."""

    orthogroup_id: str
    source_rbh_orthogroup_id: str | None
    target_rbh_orthogroup_id: str | None
    query_protein_id: str
    subject_protein_id: str
    query_record_index: int
    subject_record_index: int
    edge_kind: OrthologEdgeKind
    render_role: OrthologRenderRole
    path_id: str | None
    identity: float
    evalue: float
    bitscore: float
    alignment_length: int


@dataclass(frozen=True)
class OrthologPath:
    """Traceable ortholog/co-ortholog path inside one broad orthogroup."""

    orthogroup_id: str
    path_id: str
    protein_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    shared_protein_ids: tuple[str, ...] = ()


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
    rbh_orthogroups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    ortholog_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)
    ortholog_paths_by_orthogroup_id: dict[str, tuple[OrthologPath, ...]] = field(default_factory=dict)
    related_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)


@dataclass(frozen=True)
class ProteinBlastpResult:
    """LOSATP blastp display comparisons plus optional orthogroup metadata."""

    comparisons: list[DataFrame]
    orthogroups: OrthogroupResult | None = None


@dataclass(frozen=True)
class OrthogroupEdgeSelectionResult:
    """Orthogroup selection plus separate adjacent anchor and display edges."""

    orthogroups: OrthogroupResult
    all_edges_by_pair: dict[tuple[int, int], DataFrame]
    adjacent_anchor_edges_by_pair: dict[tuple[int, int], DataFrame]
    adjacent_display_edges_by_pair: dict[tuple[int, int], DataFrame]


@dataclass(frozen=True)
class EvidenceIndex:
    """Best directional evidence rows keyed by record pair and protein IDs."""

    best_rows_by_edge: dict[tuple[int, int, str, str], object]

    def get(
        self,
        query_record_index: int,
        subject_record_index: int,
        query_protein_id: str,
        subject_protein_id: str,
    ) -> object | None:
        return self.best_rows_by_edge.get(
            (
                int(query_record_index),
                int(subject_record_index),
                str(query_protein_id),
                str(subject_protein_id),
            )
        )

    def get_reverse(
        self,
        query_record_index: int,
        subject_record_index: int,
        query_protein_id: str,
        subject_protein_id: str,
    ) -> object | None:
        return self.get(
            int(subject_record_index),
            int(query_record_index),
            str(subject_protein_id),
            str(query_protein_id),
        )

    def strength_key(self, row: object) -> tuple[float, float, float, int, str, str]:
        return _hit_strength_key_from_row(row)


@dataclass(frozen=True)
class _NormalizedProteinHit:
    query_id: str
    subject_id: str
    query_record_index: int
    subject_record_index: int
    bitscore: float
    evalue: float
    identity: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_coverage: float
    subject_coverage: float
    min_coverage: float
    normalized_score: float
    record_pair_normalization_source: str
    domain_only: bool


@dataclass(frozen=True)
class _LocalThreshold:
    protein_id: str
    score: float
    source: str
    rbnh_count: int
    accepted_edge_count: int


@dataclass(frozen=True)
class _SplitDecision:
    accepted: bool
    child_components: tuple[tuple[str, ...], ...]
    threshold: float | None
    separation_ratio: float
    conductance: float
    stability: int
    reason: str


@dataclass(frozen=True)
class _DistributionEvidenceEdge:
    query_id: str
    subject_id: str
    row: object
    score: float
    query_score: float | None
    subject_score: float | None
    is_anchor: bool
    fallback_normalized: bool


@dataclass(frozen=True)
class _AnchorCoreEvidenceEdge:
    query_id: str
    subject_id: str
    row: object
    score: float
    reverse_score: float
    edge_kind: Literal["rbh", "coortholog"]


@dataclass(frozen=True)
class _CoreSupportCandidate:
    group_id: str
    support: float
    same_record_score: float
    cross_record_score: float
    evidence_row: object
    evidence_query_id: str
    evidence_subject_id: str
    relation_kind: OrthologEdgeKind
    min_coverage: float
    domain_only: bool
    high_confidence_pass: bool
    low_confidence_pass: bool
    reason: str


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
    locus_tag: str | None
    gene_id: str | None
    old_locus_tag: str | None
    db_xref: tuple[str, ...]
    gff_id: str | None
    parent_ids: tuple[str, ...]
    gene_parent_id: str | None
    feature_type: str
    feature_hash_start: int | None
    feature_hash_end: int | None
    feature_hash_strand: int | None


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


def _qualifier_values(feature: SeqFeature, key: str) -> tuple[str, ...]:
    value = feature.qualifiers.get(key)
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        value = (value,)
    values: list[str] = []
    for item in value:
        text = str(item).strip()
        if text:
            values.append(text)
    return tuple(values)


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


def _feature_identifier(feature: SeqFeature) -> str | None:
    feature_id = _first_qualifier(feature, "ID")
    if feature_id:
        return feature_id
    raw_id = str(getattr(feature, "id", "") or "").strip()
    if raw_id and raw_id != "<unknown id>":
        return raw_id
    return None


def _iter_record_features(
    record: SeqRecord,
) -> list[tuple[int, SeqFeature, tuple[SeqFeature, ...]]]:
    items: list[tuple[int, SeqFeature, tuple[SeqFeature, ...]]] = []

    def walk(feature: SeqFeature, ancestors: tuple[SeqFeature, ...]) -> None:
        items.append((len(items), feature, ancestors))
        sub_features = getattr(feature, "sub_features", None) or []
        for sub_feature in sub_features:
            walk(sub_feature, ancestors + (feature,))

    for feature in record.features:
        walk(feature, ())
    return items


def _build_parent_graph(
    feature_items: Sequence[tuple[int, SeqFeature, tuple[SeqFeature, ...]]],
) -> dict[str, tuple[str, tuple[str, ...]]]:
    graph: dict[str, tuple[str, tuple[str, ...]]] = {}
    for _feature_index, feature, ancestors in feature_items:
        feature_id = _feature_identifier(feature)
        if not feature_id:
            continue
        parent_ids = _qualifier_values(feature, "Parent")
        if not parent_ids and ancestors:
            ancestor_id = _feature_identifier(ancestors[-1])
            if ancestor_id:
                parent_ids = (ancestor_id,)
        graph[feature_id] = (str(feature.type), tuple(parent_ids))
    return graph


def _resolve_gene_parent_id(
    feature: SeqFeature,
    ancestors: Sequence[SeqFeature],
    parent_graph: Mapping[str, tuple[str, tuple[str, ...]]],
) -> str | None:
    for ancestor in reversed(tuple(ancestors)):
        if str(ancestor.type).lower() == "gene":
            ancestor_id = _feature_identifier(ancestor)
            if ancestor_id:
                return ancestor_id

    start_ids = list(_qualifier_values(feature, "Parent"))
    if not start_ids and ancestors:
        ancestor_id = _feature_identifier(ancestors[-1])
        if ancestor_id:
            start_ids.append(ancestor_id)

    seen: set[str] = set()

    def resolve(feature_id: str) -> str | None:
        if feature_id in seen:
            return None
        seen.add(feature_id)
        feature_type, parent_ids = parent_graph.get(feature_id, ("", ()))
        if str(feature_type).lower() == "gene":
            return feature_id
        for parent_id in parent_ids:
            resolved = resolve(parent_id)
            if resolved:
                return resolved
        return None

    for parent_id in start_ids:
        resolved = resolve(parent_id)
        if resolved:
            return resolved
    return None


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


def _feature_hash_inputs(feature: SeqFeature) -> tuple[int | None, int | None, int | None]:
    location = feature.location
    if location is None:
        return None, None, None
    if hasattr(location, "parts") and location.parts:
        location = location.parts[0]
    strand = location.strand if location.strand in {-1, 1} else None
    return int(location.start), int(location.end), strand


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
        feature_items = _iter_record_features(record)
        parent_graph = _build_parent_graph(feature_items)
        for feature_index, feature, ancestors in feature_items:
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
            hash_start, hash_end, hash_strand = _feature_hash_inputs(feature)
            synthetic_protein_id = f"gbd_r{global_record_index + 1:04d}_cds{cds_count:06d}"
            source_protein_id = _first_qualifier(feature, "protein_id")
            fasta_safe_source_id = _fasta_safe_protein_id(source_protein_id)
            if fasta_safe_source_id is not None:
                source_id_counts[fasta_safe_source_id] = source_id_counts.get(fasta_safe_source_id, 0) + 1
            gene = _first_qualifier(feature, "gene")
            product = _first_qualifier(feature, "product")
            note = _first_qualifier(feature, "note")
            locus_tag = _first_qualifier(feature, "locus_tag")
            gene_id = _first_qualifier(feature, "gene_id")
            old_locus_tag = _first_qualifier(feature, "old_locus_tag")
            db_xref = _qualifier_values(feature, "db_xref")
            gff_id = _feature_identifier(feature)
            parent_ids = _qualifier_values(feature, "Parent")
            if not parent_ids and ancestors:
                ancestor_id = _feature_identifier(ancestors[-1])
                if ancestor_id:
                    parent_ids = (ancestor_id,)
            gene_parent_id = _resolve_gene_parent_id(feature, ancestors, parent_graph)
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
                    locus_tag=locus_tag,
                    gene_id=gene_id,
                    old_locus_tag=old_locus_tag,
                    db_xref=db_xref,
                    gff_id=gff_id,
                    parent_ids=parent_ids,
                    gene_parent_id=gene_parent_id,
                    feature_type=str(feature.type),
                    feature_hash_start=hash_start,
                    feature_hash_end=hash_end,
                    feature_hash_strand=hash_strand,
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
                locus_tag=candidate.locus_tag,
                gene_id=candidate.gene_id,
                old_locus_tag=candidate.old_locus_tag,
                db_xref=candidate.db_xref,
                gff_id=candidate.gff_id,
                parent_ids=candidate.parent_ids,
                gene_parent_id=candidate.gene_parent_id,
                feature_type=candidate.feature_type,
                feature_hash_start=candidate.feature_hash_start,
                feature_hash_end=candidate.feature_hash_end,
                feature_hash_strand=candidate.feature_hash_strand,
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


def _validate_losatp_threads(threads: int | None) -> None:
    if threads is not None and int(threads) <= 0:
        raise ValidationError("losatp_threads must be > 0 or None")


def normalize_protein_blastp_mode(mode: str | None) -> ProteinBlastpMode:
    """Return a validated LOSATP blastp mode."""

    normalized = str(mode or "none").strip().lower()
    if normalized not in PROTEIN_BLASTP_MODES:
        raise ValidationError(
            "protein_blastp_mode must be one of: "
            + ", ".join(PROTEIN_BLASTP_MODES)
        )
    return normalized  # type: ignore[return-value]


def normalize_orthogroup_membership_mode(mode: str | None) -> OrthogroupMembershipMode:
    """Return the active orthogroup inference model.

    Legacy user/session values are accepted as aliases so old configurations
    load without reactivating the removed inference modes.
    """

    normalized = str(mode or ORTHOGROUP_INFERENCE_VERSION).strip().lower().replace("-", "_")
    aliases = {
        "legacy": ORTHOGROUP_INFERENCE_VERSION,
        "rbh": ORTHOGROUP_INFERENCE_VERSION,
        "rbh_only": ORTHOGROUP_INFERENCE_VERSION,
        "merge": ORTHOGROUP_INFERENCE_VERSION,
        "family": ORTHOGROUP_INFERENCE_VERSION,
        "family_merge": ORTHOGROUP_INFERENCE_VERSION,
        "local_split": ORTHOGROUP_INFERENCE_VERSION,
        "density_split": ORTHOGROUP_INFERENCE_VERSION,
        "outparalog_split": ORTHOGROUP_INFERENCE_VERSION,
        "distribution_split": ORTHOGROUP_INFERENCE_VERSION,
        "orthogroups": ORTHOGROUP_INFERENCE_VERSION,
        "anchor_core": ORTHOGROUP_INFERENCE_VERSION,
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in ORTHOGROUP_MEMBERSHIP_MODES:
        raise ValidationError(
            "orthogroup_membership_mode must be "
            f"{ORTHOGROUP_INFERENCE_VERSION} or a legacy alias: "
            + ", ".join(_LEGACY_ORTHOGROUP_MEMBERSHIP_MODES)
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


def _stable_quantile(values: Sequence[float], fraction: float) -> float:
    cleaned = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not cleaned:
        return 0.0
    if len(cleaned) == 1:
        return cleaned[0]
    bounded_fraction = min(1.0, max(0.0, float(fraction)))
    position = bounded_fraction * (len(cleaned) - 1)
    lower_index = int(math.floor(position))
    upper_index = int(math.ceil(position))
    if lower_index == upper_index:
        return cleaned[lower_index]
    lower = cleaned[lower_index]
    upper = cleaned[upper_index]
    return lower + (upper - lower) * (position - lower_index)


def _stable_median(values: Sequence[float]) -> float:
    return _stable_quantile(values, 0.5)


def _normalized_score_from_row(row: object) -> float:
    score = _row_float(row, "normalized_score", 0.0)
    return score if math.isfinite(score) and score > 0.0 else 0.0


def _normalized_hit_rank_from_row(row: object) -> tuple[float, float, float, int, str, str]:
    return (
        -_normalized_score_from_row(row),
        _row_float(row, "evalue", float("inf")),
        -_row_float(row, "identity", 0.0),
        -int(_row_float(row, "alignment_length", 0.0)),
        str(getattr(row, "query", "")),
        str(getattr(row, "subject", "")),
    )


def _sort_hits_for_normalized_query_best(hits: DataFrame) -> DataFrame:
    return hits.sort_values(
        ["query", "normalized_score", "evalue", "identity", "alignment_length", "subject"],
        ascending=[True, False, True, False, False, True],
        kind="mergesort",
    )


def _select_normalized_fit_rows(
    normalized_hits: DataFrame,
    *,
    top_fraction: float,
) -> list[tuple[float, float]]:
    if normalized_hits.empty:
        return []
    sorted_hits = normalized_hits.sort_values(
        ["length_product", "query", "subject"],
        ascending=[True, True, True],
        kind="mergesort",
    )
    bin_size = max(4, int(math.ceil(len(sorted_hits) / 8)))
    fit_points: list[tuple[float, float]] = []
    for start_index in range(0, len(sorted_hits), bin_size):
        bin_hits = sorted_hits.iloc[start_index : start_index + bin_size]
        top_count = max(1, int(math.ceil(len(bin_hits) * float(top_fraction))))
        top_hits = bin_hits.sort_values(
            ["bitscore", "query", "subject"],
            ascending=[False, True, True],
            kind="mergesort",
        ).head(top_count)
        for row in top_hits.itertuples(index=False):
            length_product = _row_float(row, "length_product", 0.0)
            bitscore = _row_float(row, "bitscore", 0.0)
            if length_product <= 0.0 or bitscore <= 0.0:
                continue
            fit_points.append((math.log10(length_product), math.log10(bitscore)))
    return fit_points


def _fit_expected_bitscore_model(
    normalized_hits: DataFrame,
    *,
    top_fraction: float,
) -> tuple[float, float] | None:
    if len(normalized_hits) < _ORTHOGROUP_SPLIT_MIN_FIT_HITS:
        return None
    fit_points = _select_normalized_fit_rows(
        normalized_hits,
        top_fraction=top_fraction,
    )
    if len(fit_points) < 3:
        return None
    x_values = [point[0] for point in fit_points]
    y_values = [point[1] for point in fit_points]
    if len(set(round(value, 9) for value in x_values)) < 2:
        return None
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    denominator = sum((x_value - x_mean) ** 2 for x_value in x_values)
    if denominator <= _ORTHOGROUP_SPLIT_EPSILON:
        return None
    slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x_values, y_values)
    ) / denominator
    intercept = y_mean - slope * x_mean
    if not (math.isfinite(slope) and math.isfinite(intercept)):
        return None
    return slope, intercept


def _empty_normalized_hit_table(columns: Sequence[str]) -> DataFrame:
    extra_columns = (
        "query_length",
        "subject_length",
        "length_product",
        "query_coverage",
        "subject_coverage",
        "min_coverage",
        "normalized_score",
        "normalization_fallback",
        "record_pair_normalization_source",
        "domain_only",
    )
    return pd.DataFrame(columns=tuple(dict.fromkeys([*columns, *extra_columns])))


def _normalize_directional_hit_table(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
    *,
    top_fraction: float = _ORTHOGROUP_SPLIT_TOP_FRACTION,
    min_coverage: float = _ORTHOGROUP_SPLIT_MIN_COVERAGE,
) -> DataFrame:
    if hits is None or hits.empty:
        columns = tuple(hits.columns) if hits is not None else tuple(COMPARISON_COLUMNS)
        return _empty_normalized_hit_table(columns)
    _validate_comparison_columns(hits)
    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    rows: list[dict[str, object]] = []
    for row in coerced_hits.itertuples(index=False):
        query_id = str(row.query)
        subject_id = str(row.subject)
        query_protein = protein_map.get(query_id)
        subject_protein = protein_map.get(subject_id)
        if query_protein is None or subject_protein is None:
            continue
        query_length = int(query_protein.protein_length)
        subject_length = int(subject_protein.protein_length)
        bitscore = _row_float(row, "bitscore", 0.0)
        alignment_length = int(_row_float(row, "alignment_length", 0.0))
        if query_length <= 0 or subject_length <= 0 or bitscore <= 0.0 or alignment_length <= 0:
            continue
        query_coverage = min(1.0, float(alignment_length) / float(query_length))
        subject_coverage = min(1.0, float(alignment_length) / float(subject_length))
        min_hit_coverage = min(query_coverage, subject_coverage)
        if min_hit_coverage < float(min_coverage):
            continue
        length_product = float(query_length * subject_length)
        record = dict(row._asdict())
        record.update(
            {
                "query_length": query_length,
                "subject_length": subject_length,
                "length_product": length_product,
                "query_coverage": query_coverage,
                "subject_coverage": subject_coverage,
                "min_coverage": min_hit_coverage,
            }
        )
        rows.append(record)

    if not rows:
        return _empty_normalized_hit_table(tuple(coerced_hits.columns))

    normalized_hits = pd.DataFrame.from_records(rows)
    model = _fit_expected_bitscore_model(
        normalized_hits,
        top_fraction=top_fraction,
    )
    fallback = model is None
    normalized_scores: list[float] = []
    for row in normalized_hits.itertuples(index=False):
        bitscore = _row_float(row, "bitscore", 0.0)
        length_product = _row_float(row, "length_product", 0.0)
        if fallback:
            denominator = math.sqrt(max(length_product, _ORTHOGROUP_SPLIT_EPSILON))
            normalized_scores.append(bitscore / denominator)
            continue
        slope, intercept = model
        expected_bitscore = 10 ** (slope * math.log10(length_product) + intercept)
        if not math.isfinite(expected_bitscore) or expected_bitscore <= 0.0:
            denominator = math.sqrt(max(length_product, _ORTHOGROUP_SPLIT_EPSILON))
            normalized_scores.append(bitscore / denominator)
        else:
            normalized_scores.append(bitscore / expected_bitscore)
    normalized_hits["normalized_score"] = normalized_scores
    normalized_hits["normalization_fallback"] = bool(fallback)
    normalized_hits["record_pair_normalization_source"] = (
        "sqrt_length_fallback" if fallback else "fit"
    )
    normalized_hits["domain_only"] = normalized_hits["min_coverage"].astype(float) <= _DOMAIN_ONLY_MAX_MIN_COVERAGE
    return normalized_hits.reset_index(drop=True)


def _normalize_directional_hit_tables(
    directional_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    min_coverage: float = _ORTHOGROUP_SPLIT_MIN_COVERAGE,
) -> dict[tuple[int, int], DataFrame]:
    return {
        (int(pair[0]), int(pair[1])): _normalize_directional_hit_table(
            hits,
            protein_map,
            min_coverage=float(min_coverage),
        )
        for pair, hits in directional_tables.items()
    }


def _cap_hits_per_query_by_normalized_score(
    hits: DataFrame,
    *,
    max_hits: int,
    distinct_subjects: bool = True,
) -> DataFrame:
    _validate_max_hits(max_hits, option_name="orthogroup_member_max_hits")
    if hits.empty:
        return hits.copy()
    sorted_hits = _sort_hits_for_normalized_query_best(hits)
    if distinct_subjects:
        sorted_hits = sorted_hits.drop_duplicates(["query", "subject"], keep="first")
    return (
        sorted_hits.groupby("query", group_keys=False, sort=False)
        .head(int(max_hits))
        .reset_index(drop=True)
    )


def _select_normalized_best_hits_per_query(hits: DataFrame) -> DataFrame:
    if hits.empty:
        return hits.copy()
    return (
        _sort_hits_for_normalized_query_best(hits)
        .drop_duplicates(["query", "subject"], keep="first")
        .drop_duplicates("query", keep="first")
        .reset_index(drop=True)
    )


def _select_reciprocal_best_normalized_hit_edges(
    forward_hits: DataFrame,
    reverse_hits: DataFrame,
) -> DataFrame:
    if forward_hits.empty or reverse_hits.empty:
        return forward_hits.iloc[0:0].copy()
    forward_best = _select_normalized_best_hits_per_query(forward_hits)
    reverse_best = _select_normalized_best_hits_per_query(reverse_hits)
    reverse_best_by_query = {
        str(row.query): str(row.subject)
        for row in reverse_best.itertuples(index=False)
    }
    keep_mask = [
        reverse_best_by_query.get(str(row.subject)) == str(row.query)
        for row in forward_best.itertuples(index=False)
    ]
    return forward_best.loc[keep_mask].reset_index(drop=True)


def _select_rbnh_edge_tables(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
) -> dict[tuple[int, int], DataFrame]:
    unordered_pairs = sorted(
        {
            (min(int(query_index), int(subject_index)), max(int(query_index), int(subject_index)))
            for query_index, subject_index in normalized_tables
            if int(query_index) != int(subject_index)
        }
    )
    rbnh_edge_tables: dict[tuple[int, int], DataFrame] = {}
    for query_index, subject_index in unordered_pairs:
        forward_hits = normalized_tables.get((query_index, subject_index), _empty_normalized_hit_table(COMPARISON_COLUMNS))
        reverse_hits = normalized_tables.get((subject_index, query_index), _empty_normalized_hit_table(COMPARISON_COLUMNS))
        rbnh_edge_tables[(query_index, subject_index)] = _select_reciprocal_best_normalized_hit_edges(
            forward_hits,
            reverse_hits,
        )
    return rbnh_edge_tables


def _comparison_columns_only(hits: DataFrame) -> DataFrame:
    if hits.empty:
        return pd.DataFrame(columns=COMPARISON_COLUMNS)
    return hits.loc[:, list(COMPARISON_COLUMNS)].copy()


def _normalized_hit_from_row(
    row: object,
    protein_map: Mapping[str, CdsProtein],
) -> _NormalizedProteinHit | None:
    query_id = str(getattr(row, "query", ""))
    subject_id = str(getattr(row, "subject", ""))
    query_protein = protein_map.get(query_id)
    subject_protein = protein_map.get(subject_id)
    if query_protein is None or subject_protein is None:
        return None
    return _NormalizedProteinHit(
        query_id=query_id,
        subject_id=subject_id,
        query_record_index=int(query_protein.record_index),
        subject_record_index=int(subject_protein.record_index),
        bitscore=_row_float(row, "bitscore", 0.0),
        evalue=_row_float(row, "evalue", float("inf")),
        identity=_row_float(row, "identity", 0.0),
        alignment_length=int(_row_float(row, "alignment_length", 0.0)),
        query_length=int(_row_float(row, "query_length", query_protein.protein_length)),
        subject_length=int(_row_float(row, "subject_length", subject_protein.protein_length)),
        query_coverage=_row_float(row, "query_coverage", 0.0),
        subject_coverage=_row_float(row, "subject_coverage", 0.0),
        min_coverage=_row_float(row, "min_coverage", 0.0),
        normalized_score=_normalized_score_from_row(row),
        record_pair_normalization_source=str(
            getattr(row, "record_pair_normalization_source", "")
            or ("sqrt_length_fallback" if bool(getattr(row, "normalization_fallback", False)) else "fit")
        ),
        domain_only=bool(getattr(row, "domain_only", False)),
    )


def _canonical_edge_key_from_ids(
    query_id: str,
    subject_id: str,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[str, str] | None:
    if query_id not in protein_map or subject_id not in protein_map:
        return None
    if protein_map[query_id].record_index == protein_map[subject_id].record_index:
        return None
    return _canonical_edge_endpoint_ids(query_id, subject_id, protein_map)


def _anchor_pairs_from_edge_tables(
    rbnh_edge_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> set[tuple[str, str]]:
    anchor_pairs: set[tuple[str, str]] = set()
    for hits in rbnh_edge_tables.values():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            key = _canonical_edge_key_from_ids(str(row.query), str(row.subject), protein_map)
            if key is not None:
                anchor_pairs.add(key)
    return anchor_pairs


def _collect_distribution_evidence_edges(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    anchor_pairs: set[tuple[str, str]],
    max_hits: int,
) -> tuple[_DistributionEvidenceEdge, ...]:
    edge_data: dict[tuple[str, str], dict[str, object]] = {}
    for hits in normalized_tables.values():
        if hits is None or hits.empty:
            continue
        capped_hits = _cap_hits_per_query_by_normalized_score(
            hits,
            max_hits=int(max_hits),
            distinct_subjects=True,
        )
        for row in capped_hits.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            key = _canonical_edge_key_from_ids(query_id, subject_id, protein_map)
            if key is None:
                continue
            score = _normalized_score_from_row(row)
            if score <= 0.0:
                continue
            record = edge_data.setdefault(
                key,
                {
                    "row": row,
                    "rank": _normalized_hit_rank_from_row(row),
                    "score": score,
                    "endpoint_scores": {},
                    "fallback": bool(getattr(row, "normalization_fallback", False)),
                },
            )
            endpoint_scores = record["endpoint_scores"]
            if isinstance(endpoint_scores, dict):
                endpoint_scores[query_id] = max(float(endpoint_scores.get(query_id, 0.0)), score)
            if score > float(record["score"]):
                record["score"] = score
            record["fallback"] = bool(record["fallback"]) and bool(getattr(row, "normalization_fallback", False))
            rank = _normalized_hit_rank_from_row(row)
            if rank < record["rank"]:
                record["row"] = row
                record["rank"] = rank

    evidence_edges: list[_DistributionEvidenceEdge] = []
    for (query_id, subject_id), record in edge_data.items():
        endpoint_scores = record["endpoint_scores"]
        query_score = None
        subject_score = None
        if isinstance(endpoint_scores, dict):
            query_score = endpoint_scores.get(query_id)
            subject_score = endpoint_scores.get(subject_id)
        evidence_edges.append(
            _DistributionEvidenceEdge(
                query_id=query_id,
                subject_id=subject_id,
                row=record["row"],
                score=float(record["score"]),
                query_score=float(query_score) if query_score is not None else None,
                subject_score=float(subject_score) if subject_score is not None else None,
                is_anchor=(query_id, subject_id) in anchor_pairs,
                fallback_normalized=bool(record["fallback"]),
            )
        )
    return tuple(sorted(evidence_edges, key=_distribution_edge_sort_key))


def _distribution_edge_sort_key(edge: _DistributionEvidenceEdge) -> tuple[float, float, float, int, str, str]:
    return (
        -float(edge.score),
        _row_float(edge.row, "evalue", float("inf")),
        -_row_float(edge.row, "identity", 0.0),
        -int(_row_float(edge.row, "alignment_length", 0.0)),
        edge.query_id,
        edge.subject_id,
    )


def _score_valley_threshold(scores: Sequence[float], rbnh_floor: float) -> float | None:
    unique_scores = sorted(
        {float(score) for score in scores if math.isfinite(float(score)) and float(score) > 0.0},
        reverse=True,
    )
    if len(unique_scores) < 3 or rbnh_floor <= 0.0:
        return None
    best_ratio = 0.0
    best_threshold: float | None = None
    for high_score, low_score in zip(unique_scores, unique_scores[1:]):
        if low_score <= 0.0:
            continue
        ratio = high_score / low_score
        if (
            ratio >= _ORTHOGROUP_SPLIT_MIN_SEPARATION_RATIO
            and high_score >= rbnh_floor * 0.5
            and low_score < rbnh_floor * 0.9
            and ratio > best_ratio
        ):
            best_ratio = ratio
            best_threshold = (high_score + low_score) / 2.0
    return best_threshold


def _derive_local_thresholds(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
    rbnh_edge_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> dict[str, _LocalThreshold]:
    scores_by_protein: dict[str, list[float]] = {}
    rbnh_scores_by_protein: dict[str, list[float]] = {}
    for hits in normalized_tables.values():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            hit = _normalized_hit_from_row(row, protein_map)
            if hit is None or hit.query_record_index == hit.subject_record_index:
                continue
            if hit.normalized_score <= 0.0:
                continue
            scores_by_protein.setdefault(hit.query_id, []).append(hit.normalized_score)

    for hits in rbnh_edge_tables.values():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            hit = _normalized_hit_from_row(row, protein_map)
            if hit is None or hit.normalized_score <= 0.0:
                continue
            for protein_id in (hit.query_id, hit.subject_id):
                scores_by_protein.setdefault(protein_id, []).append(hit.normalized_score)
                rbnh_scores_by_protein.setdefault(protein_id, []).append(hit.normalized_score)

    thresholds: dict[str, _LocalThreshold] = {}
    for protein_id, scores in sorted(scores_by_protein.items()):
        positive_scores = [score for score in scores if score > 0.0 and math.isfinite(score)]
        if not positive_scores:
            continue
        rbnh_scores = [
            score
            for score in rbnh_scores_by_protein.get(protein_id, [])
            if score > 0.0 and math.isfinite(score)
        ]
        if rbnh_scores:
            rbnh_floor = _stable_quantile(rbnh_scores, 0.25)
            threshold = rbnh_floor * 0.70
            source = "rbnh"
            valley_threshold = _score_valley_threshold(positive_scores, rbnh_floor)
            if valley_threshold is not None and valley_threshold > threshold:
                threshold = valley_threshold
                source = "valley"
        else:
            threshold = max(positive_scores) * 0.85
            source = "fallback"
        thresholds[protein_id] = _LocalThreshold(
            protein_id=protein_id,
            score=float(threshold),
            source=source,
            rbnh_count=len(rbnh_scores),
            accepted_edge_count=0,
        )
    return thresholds


def _edge_passes_local_thresholds(
    edge: _DistributionEvidenceEdge,
    thresholds: Mapping[str, _LocalThreshold],
) -> bool:
    if edge.is_anchor:
        return True
    query_threshold = thresholds.get(edge.query_id)
    subject_threshold = thresholds.get(edge.subject_id)
    query_passes = (
        query_threshold is not None
        and edge.query_score is not None
        and float(edge.query_score) >= float(query_threshold.score)
    )
    subject_passes = (
        subject_threshold is not None
        and edge.subject_score is not None
        and float(edge.subject_score) >= float(subject_threshold.score)
    )
    if query_threshold is not None and subject_threshold is not None:
        return query_passes and subject_passes
    if query_threshold is not None:
        return query_passes
    if subject_threshold is not None:
        return subject_passes
    return False


def _select_local_threshold_edges(
    evidence_edges: Sequence[_DistributionEvidenceEdge],
    thresholds: Mapping[str, _LocalThreshold],
) -> tuple[tuple[_DistributionEvidenceEdge, ...], dict[str, _LocalThreshold]]:
    accepted_edges: list[_DistributionEvidenceEdge] = []
    accepted_counts: dict[str, int] = {}
    for edge in evidence_edges:
        if not _edge_passes_local_thresholds(edge, thresholds):
            continue
        accepted_edges.append(edge)
        accepted_counts[edge.query_id] = accepted_counts.get(edge.query_id, 0) + 1
        accepted_counts[edge.subject_id] = accepted_counts.get(edge.subject_id, 0) + 1

    updated_thresholds = {
        protein_id: replace(
            threshold,
            accepted_edge_count=accepted_counts.get(protein_id, 0),
        )
        for protein_id, threshold in thresholds.items()
    }
    for threshold in updated_thresholds.values():
        logger.debug(
            "orthogroup threshold: protein=%s source=%s score=%.6g accepted_edges=%s",
            threshold.protein_id,
            threshold.source,
            threshold.score,
            threshold.accepted_edge_count,
        )
    return tuple(sorted(accepted_edges, key=_distribution_edge_sort_key)), updated_thresholds


def _components_from_distribution_edges(
    member_ids: set[str],
    edges: Sequence[_DistributionEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[str, ...], ...]:
    union_find = _UnionFind()
    for member_id in member_ids:
        union_find.add(member_id)
    for edge in edges:
        if edge.query_id in member_ids and edge.subject_id in member_ids:
            union_find.union(edge.query_id, edge.subject_id)
    components: dict[str, set[str]] = {}
    for member_id in member_ids:
        components.setdefault(union_find.find(member_id), set()).add(member_id)
    return tuple(
        tuple(sorted(component, key=lambda protein_id: _protein_sort_key(protein_map[protein_id])))
        for component in sorted(
            components.values(),
            key=lambda component: min(_protein_sort_key(protein_map[protein_id]) for protein_id in component),
        )
    )


def _partition_key(child_components: Sequence[Sequence[str]]) -> tuple[tuple[str, ...], ...]:
    return tuple(sorted(tuple(sorted(component)) for component in child_components if component))


def _score_split_partition(
    child_components: Sequence[Sequence[str]],
    edges: Sequence[_DistributionEvidenceEdge],
) -> tuple[float, float]:
    component_by_protein: dict[str, int] = {}
    for component_index, component in enumerate(child_components):
        for protein_id in component:
            component_by_protein[str(protein_id)] = component_index
    internal_scores_by_component: dict[int, list[float]] = {}
    cross_scores: list[float] = []
    internal_weight = 0.0
    cross_weight = 0.0
    for edge in edges:
        query_component = component_by_protein.get(edge.query_id)
        subject_component = component_by_protein.get(edge.subject_id)
        if query_component is None or subject_component is None:
            continue
        if query_component == subject_component:
            internal_scores_by_component.setdefault(query_component, []).append(float(edge.score))
            internal_weight += float(edge.score)
        else:
            cross_scores.append(float(edge.score))
            cross_weight += float(edge.score)

    within_medians = [
        _stable_median(scores)
        for component_index, scores in internal_scores_by_component.items()
        if len(child_components[component_index]) >= _ORTHOGROUP_SPLIT_MIN_CHILD_SIZE and scores
    ]
    within_median_min = min(within_medians) if within_medians else 0.0
    cross_p95 = _stable_quantile(cross_scores, 0.95) if cross_scores else _ORTHOGROUP_SPLIT_EPSILON
    separation_ratio = within_median_min / max(cross_p95, _ORTHOGROUP_SPLIT_EPSILON)
    conductance = cross_weight / max(cross_weight + internal_weight, _ORTHOGROUP_SPLIT_EPSILON)
    return separation_ratio, conductance


def _child_has_strong_internal_evidence(
    child_component: Sequence[str],
    edges: Sequence[_DistributionEvidenceEdge],
) -> bool:
    child_ids = {str(protein_id) for protein_id in child_component}
    if len(child_ids) <= 1:
        return True
    return any(
        edge.query_id in child_ids
        and edge.subject_id in child_ids
        and (edge.is_anchor or not edge.fallback_normalized)
        for edge in edges
    )


def _evaluate_split_partition(
    parent_component: Sequence[str],
    child_components: Sequence[Sequence[str]],
    edges: Sequence[_DistributionEvidenceEdge],
    *,
    threshold: float | None,
    stability: int,
    reason: str,
) -> _SplitDecision:
    normalized_children = tuple(
        tuple(child)
        for child in child_components
        if child
    )
    parent_key = _partition_key((tuple(parent_component),))
    child_key = _partition_key(normalized_children)
    if len(normalized_children) < 2 or child_key == parent_key:
        return _SplitDecision(
            accepted=False,
            child_components=normalized_children,
            threshold=threshold,
            separation_ratio=0.0,
            conductance=1.0,
            stability=stability,
            reason="not_split",
        )
    large_children = [
        child
        for child in normalized_children
        if len(child) >= _ORTHOGROUP_SPLIT_MIN_CHILD_SIZE
    ]
    if len(large_children) < 2:
        return _SplitDecision(
            accepted=False,
            child_components=normalized_children,
            threshold=threshold,
            separation_ratio=0.0,
            conductance=1.0,
            stability=stability,
            reason="small_child",
        )
    if stability < _ORTHOGROUP_SPLIT_MIN_STABILITY:
        separation_ratio, conductance = _score_split_partition(normalized_children, edges)
        return _SplitDecision(
            accepted=False,
            child_components=normalized_children,
            threshold=threshold,
            separation_ratio=separation_ratio,
            conductance=conductance,
            stability=stability,
            reason="unstable",
        )
    if not all(_child_has_strong_internal_evidence(child, edges) for child in large_children):
        separation_ratio, conductance = _score_split_partition(normalized_children, edges)
        return _SplitDecision(
            accepted=False,
            child_components=normalized_children,
            threshold=threshold,
            separation_ratio=separation_ratio,
            conductance=conductance,
            stability=stability,
            reason="fallback_only_child",
        )
    separation_ratio, conductance = _score_split_partition(normalized_children, edges)
    accepted = (
        separation_ratio >= _ORTHOGROUP_SPLIT_MIN_SEPARATION_RATIO
        and conductance <= _ORTHOGROUP_SPLIT_MAX_CONDUCTANCE
    )
    return _SplitDecision(
        accepted=accepted,
        child_components=normalized_children,
        threshold=threshold,
        separation_ratio=separation_ratio,
        conductance=conductance,
        stability=stability,
        reason=reason if accepted else "weak_separation",
    )


def _threshold_sweep_split_decision(
    member_ids: set[str],
    edges: Sequence[_DistributionEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> _SplitDecision:
    unique_scores = sorted(
        {float(edge.score) for edge in edges if edge.query_id in member_ids and edge.subject_id in member_ids},
        reverse=True,
    )
    if len(unique_scores) < 2:
        return _SplitDecision(False, (), None, 0.0, 1.0, 0, "no_sweep_levels")

    partition_counts: dict[tuple[tuple[str, ...], ...], tuple[tuple[tuple[str, ...], ...], int, float]] = {}
    for threshold in unique_scores[:-1]:
        kept_edges = [
            edge
            for edge in edges
            if edge.query_id in member_ids
            and edge.subject_id in member_ids
            and float(edge.score) >= float(threshold)
        ]
        child_components = _components_from_distribution_edges(member_ids, kept_edges, protein_map)
        key = _partition_key(child_components)
        if len(key) < 2:
            continue
        current = partition_counts.get(key)
        if current is None:
            partition_counts[key] = (child_components, 1, float(threshold))
        else:
            components, count, first_threshold = current
            partition_counts[key] = (components, count + 1, first_threshold)

    best_decision = _SplitDecision(False, (), None, 0.0, 1.0, 0, "no_candidate")
    for child_components, stability, threshold in partition_counts.values():
        decision = _evaluate_split_partition(
            tuple(sorted(member_ids)),
            child_components,
            edges,
            threshold=threshold,
            stability=stability,
            reason="threshold_sweep",
        )
        if not decision.accepted:
            continue
        if (
            not best_decision.accepted
            or decision.separation_ratio > best_decision.separation_ratio
            or (
                decision.separation_ratio == best_decision.separation_ratio
                and decision.conductance < best_decision.conductance
            )
        ):
            best_decision = decision
    return best_decision


def _split_component_recursively(
    member_ids: set[str],
    *,
    all_edges: Sequence[_DistributionEvidenceEdge],
    accepted_edges: Sequence[_DistributionEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[str, ...], ...]:
    if len(member_ids) <= 2:
        return (
            tuple(sorted(member_ids, key=lambda protein_id: _protein_sort_key(protein_map[protein_id]))),
        )
    internal_edges = [
        edge
        for edge in all_edges
        if edge.query_id in member_ids and edge.subject_id in member_ids
    ]
    internal_accepted_edges = [
        edge
        for edge in accepted_edges
        if edge.query_id in member_ids and edge.subject_id in member_ids
    ]
    local_components = _components_from_distribution_edges(
        member_ids,
        internal_accepted_edges,
        protein_map,
    )
    local_decision = _evaluate_split_partition(
        tuple(sorted(member_ids)),
        local_components,
        internal_edges,
        threshold=None,
        stability=_ORTHOGROUP_SPLIT_MIN_STABILITY,
        reason="local_threshold",
    )
    decision = local_decision
    if not decision.accepted:
        decision = _threshold_sweep_split_decision(member_ids, internal_edges, protein_map)

    if not decision.accepted:
        return (
            tuple(sorted(member_ids, key=lambda protein_id: _protein_sort_key(protein_map[protein_id]))),
        )

    logger.debug(
        "orthogroup split: members=%s accepted=true children=%s separation=%.6g conductance=%.6g",
        len(member_ids),
        len(decision.child_components),
        decision.separation_ratio,
        decision.conductance,
    )
    split_components: list[tuple[str, ...]] = []
    for child_component in decision.child_components:
        split_components.extend(
            _split_component_recursively(
                set(child_component),
                all_edges=internal_edges,
                accepted_edges=internal_accepted_edges,
                protein_map=protein_map,
            )
        )
    return tuple(split_components)


def _representative_ids_for_distribution_groups(
    group_member_ids: Mapping[str, set[str]],
    evidence_edges: Sequence[_DistributionEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> dict[str, set[str]]:
    group_by_protein = {
        protein_id: group_id
        for group_id, member_ids in group_member_ids.items()
        for protein_id in member_ids
    }
    member_ranks: dict[str, tuple[float, float, float, float]] = {}
    for edge in evidence_edges:
        if group_by_protein.get(edge.query_id) != group_by_protein.get(edge.subject_id):
            continue
        rank = _member_rank_from_row(edge.row)
        if _is_better_member_rank(rank, member_ranks.get(edge.query_id)):
            member_ranks[edge.query_id] = rank
        if _is_better_member_rank(rank, member_ranks.get(edge.subject_id)):
            member_ranks[edge.subject_id] = rank

    representative_ids_by_group: dict[str, set[str]] = {}
    for group_id, member_ids in group_member_ids.items():
        members_by_record: dict[int, list[str]] = {}
        for member_id in member_ids:
            protein = protein_map.get(member_id)
            if protein is None:
                continue
            members_by_record.setdefault(int(protein.record_index), []).append(member_id)
        representative_ids_by_group[group_id] = set()
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
            representative_ids_by_group[group_id].add(representative_id)
    return representative_ids_by_group


def _build_distribution_split_orthogroups_from_normalized_tables(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
    rbnh_edge_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    include_singletons: bool,
    member_max_hits: int,
) -> OrthogroupResult:
    anchor_pairs = _anchor_pairs_from_edge_tables(rbnh_edge_tables, protein_map)
    evidence_edges = _collect_distribution_evidence_edges(
        normalized_tables,
        protein_map,
        anchor_pairs=anchor_pairs,
        max_hits=int(member_max_hits),
    )
    thresholds = _derive_local_thresholds(
        normalized_tables,
        rbnh_edge_tables,
        protein_map,
    )
    accepted_edges, _updated_thresholds = _select_local_threshold_edges(
        evidence_edges,
        thresholds,
    )

    broad_member_ids: set[str] = set()
    for edge in evidence_edges:
        broad_member_ids.update((edge.query_id, edge.subject_id))
    if include_singletons:
        broad_member_ids.update(str(protein_id) for protein_id in protein_map)

    broad_components = _components_from_distribution_edges(
        broad_member_ids,
        evidence_edges,
        protein_map,
    )
    final_components: list[tuple[str, ...]] = []
    for broad_component in broad_components:
        final_components.extend(
            _split_component_recursively(
                set(broad_component),
                all_edges=evidence_edges,
                accepted_edges=accepted_edges,
                protein_map=protein_map,
            )
        )
    if include_singletons:
        assigned_ids = {
            protein_id
            for component in final_components
            for protein_id in component
        }
        for protein_id in sorted(
            set(protein_map).difference(assigned_ids),
            key=lambda member_id: _protein_sort_key(protein_map[member_id]),
        ):
            final_components.append((protein_id,))

    final_components = sorted(
        {
            tuple(sorted(component, key=lambda protein_id: _protein_sort_key(protein_map[protein_id])))
            for component in final_components
            if component
        },
        key=lambda component: min(_protein_sort_key(protein_map[member_id]) for member_id in component),
    )
    group_member_ids: dict[str, set[str]] = {
        f"og_{component_index}": set(component)
        for component_index, component in enumerate(final_components, start=1)
    }
    group_order = list(group_member_ids)
    representative_ids_by_group = _representative_ids_for_distribution_groups(
        group_member_ids,
        [edge for edge in accepted_edges if edge.is_anchor or _edge_passes_local_thresholds(edge, thresholds)],
        protein_map,
    )
    group_by_protein = {
        protein_id: group_id
        for group_id, member_ids in group_member_ids.items()
        for protein_id in member_ids
    }

    ortholog_edges_by_group: dict[str, list[OrthologEdge]] = {}
    rbh_members_by_group: dict[str, set[str]] = {
        group_id: set()
        for group_id in group_member_ids
    }
    accepted_keys = {
        (edge.query_id, edge.subject_id)
        for edge in accepted_edges
    }
    for edge in sorted(evidence_edges, key=_distribution_edge_sort_key):
        group_id = group_by_protein.get(edge.query_id)
        if group_id is None or group_by_protein.get(edge.subject_id) != group_id:
            continue
        if not edge.is_anchor and (edge.query_id, edge.subject_id) not in accepted_keys:
            continue
        if edge.is_anchor:
            rbh_members_by_group.setdefault(group_id, set()).update((edge.query_id, edge.subject_id))
        existing_edge_ids = {
            (item.query_protein_id, item.subject_protein_id, item.edge_kind)
            for item in ortholog_edges_by_group.get(group_id, [])
        }
        edge_kind: OrthologEdgeKind = "rbh" if edge.is_anchor else "coortholog"
        if (
            (edge.query_id, edge.subject_id, edge_kind) in existing_edge_ids
            or (edge.subject_id, edge.query_id, edge_kind) in existing_edge_ids
        ):
            continue
        ortholog_edges_by_group.setdefault(group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=group_id,
                source_rbh_orthogroup_id=group_id,
                target_rbh_orthogroup_id=group_id,
                query_id=edge.query_id,
                subject_id=edge.subject_id,
                row=edge.row,
                protein_map=protein_map,
                edge_kind=edge_kind,
                render_role="block_anchor" if edge.is_anchor else "display_edge",
            )
        )

    for group_id, member_ids in group_member_ids.items():
        if not rbh_members_by_group.get(group_id):
            rbh_members_by_group[group_id] = set(member_ids)

    updated_edges_by_group, paths_by_group = _build_ortholog_paths(
        {
            group_id: sorted(
                edges,
                key=lambda edge: (
                    edge.query_record_index,
                    edge.subject_record_index,
                    edge.query_protein_id,
                    edge.subject_protein_id,
                    edge.edge_kind,
                ),
            )
            for group_id, edges in ortholog_edges_by_group.items()
        },
        protein_map,
    )

    return _orthogroup_result_from_member_ids(
        group_member_ids,
        representative_ids_by_group,
        protein_map,
        group_order=group_order,
        rbh_orthogroups={
            group_id: tuple(
                sorted(member_ids, key=lambda member_id: _protein_sort_key(protein_map[member_id]))
            )
            for group_id, member_ids in rbh_members_by_group.items()
        },
        ortholog_edges_by_orthogroup_id=updated_edges_by_group,
        ortholog_paths_by_orthogroup_id=paths_by_group,
        related_edges_by_orthogroup_id={},
    )


def _select_distribution_split_orthogroup_edges_from_directional_hits(
    directional_hits_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    record_count: int | None,
    include_singletons: bool,
    orthogroup_member_max_hits: int,
    max_related_edges_per_orthogroup: int,
) -> OrthogroupEdgeSelectionResult:
    if int(max_related_edges_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    normalized_tables = _normalize_directional_hit_tables(
        directional_hits_by_pair,
        protein_map,
    )
    rbnh_edge_tables = _select_rbnh_edge_tables(normalized_tables)
    if record_count is None:
        record_count = 0
        for query_index, subject_index in normalized_tables:
            record_count = max(record_count, int(query_index) + 1, int(subject_index) + 1)

    all_edges_by_pair = {
        pair: _comparison_columns_only(table)
        for pair, table in rbnh_edge_tables.items()
    }
    adjacent_anchor_edges_by_pair = {
        pair: _comparison_columns_only(table)
        for pair, table in rbnh_edge_tables.items()
        if int(pair[1]) == int(pair[0]) + 1
    }
    adjacent_candidate_edges_by_pair = {
        (query_index, query_index + 1): _comparison_columns_only(
            directional_hits_by_pair.get(
                (query_index, query_index + 1),
                _empty_comparison_hits(),
            )
        )
        for query_index in range(max(0, int(record_count) - 1))
    }
    for query_index in range(max(0, int(record_count) - 1)):
        adjacent_anchor_edges_by_pair.setdefault(
            (query_index, query_index + 1),
            _empty_comparison_hits(),
        )

    orthogroups = _build_distribution_split_orthogroups_from_normalized_tables(
        normalized_tables,
        rbnh_edge_tables,
        protein_map,
        include_singletons=include_singletons,
        member_max_hits=int(orthogroup_member_max_hits),
    )
    adjacent_display_edges_by_pair = _build_adjacent_display_edges_by_pair(
        adjacent_anchor_edges_by_pair,
        orthogroups,
        record_count=int(record_count),
        max_display_edges_per_orthogroup=int(max_related_edges_per_orthogroup),
        adjacent_candidate_edges_by_pair=adjacent_candidate_edges_by_pair,
    )
    return OrthogroupEdgeSelectionResult(
        orthogroups=orthogroups,
        all_edges_by_pair=all_edges_by_pair,
        adjacent_anchor_edges_by_pair=adjacent_anchor_edges_by_pair,
        adjacent_display_edges_by_pair=adjacent_display_edges_by_pair,
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
    "gene": 80,
    "product": 55,
    "note": 35,
    "label": 10,
}
_ORTHOGROUP_NAME_SOURCE_ORDER = ("gene", "product", "note", "label")
_ORTHOGROUP_NAME_SOURCE_RANK = {
    source: index for index, source in enumerate(_ORTHOGROUP_NAME_SOURCE_ORDER)
}
_NOTE_PREFIX_RE = re.compile(r"^(?:product|gene|note)\s*:\s*", re.IGNORECASE)
_DUF_CANDIDATE_RE = re.compile(r"\bDUF\d+\b.*\bdomain-containing protein\b", re.IGNORECASE)
_ANNOTATION_PROVENANCE_NOTE_RE = re.compile(
    r"^derived by automated computational analysis\b",
    re.IGNORECASE,
)
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
    if source == "note" and _ANNOTATION_PROVENANCE_NOTE_RE.search(normalized):
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

    return sorted(candidates, key=_orthogroup_name_candidate_sort_key)


def _orthogroup_name_candidate_sort_key(
    candidate: OrthogroupNameCandidate,
) -> tuple[int, float, int, int, int, str]:
    gene_priority = (
        candidate.source == "gene"
        and not _is_weak_orthogroup_name_candidate(candidate.text, candidate.source)
        and (
            candidate.record_coverage_count >= 2
            or candidate.representative_count > 0
        )
    )
    return (
        0 if gene_priority else 1,
        -candidate.score,
        _ORTHOGROUP_NAME_SOURCE_RANK.get(candidate.source, 99),
        -candidate.record_coverage_count,
        -candidate.member_count,
        candidate.text.casefold(),
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

        if top_candidate is not None and confidence not in {"none", "low"}:
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


def _hit_strength_key_from_row(row: object) -> tuple[float, float, float, int, str, str]:
    return (
        _row_float(row, "evalue", float("inf")),
        -_row_float(row, "bitscore", 0.0),
        -_row_float(row, "identity", 0.0),
        -int(_row_float(row, "alignment_length", 0.0)),
        str(getattr(row, "query", "")),
        str(getattr(row, "subject", "")),
    )


def build_pair_evidence_index(
    directional_tables: Mapping[tuple[int, int], DataFrame],
) -> EvidenceIndex:
    """Build best directional evidence lookups from filtered hit tables."""

    best_rows_by_edge: dict[tuple[int, int, str, str], object] = {}
    for (query_index, subject_index), hits in directional_tables.items():
        if hits is None or hits.empty:
            continue
        missing_columns = {"query", "subject"}.difference(hits.columns)
        if missing_columns:
            raise ParseError(
                "LOSATP blastp output is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )
        coerced_hits = _coerce_outfmt6_numeric_columns(hits)
        for row in coerced_hits.itertuples(index=False):
            key = (
                int(query_index),
                int(subject_index),
                str(row.query),
                str(row.subject),
            )
            current = best_rows_by_edge.get(key)
            if current is None or _hit_strength_key_from_row(row) < _hit_strength_key_from_row(current):
                best_rows_by_edge[key] = row
    return EvidenceIndex(best_rows_by_edge=best_rows_by_edge)


def _edge_id(edge: OrthologEdge) -> str:
    return (
        f"{edge.orthogroup_id}:"
        f"{edge.query_record_index}:{edge.query_protein_id}->"
        f"{edge.subject_record_index}:{edge.subject_protein_id}:"
        f"{edge.edge_kind}"
    )


def _make_orthogroup_member(
    orthogroup_id: str,
    protein: CdsProtein,
    *,
    representative: bool,
    role: OrthogroupMemberRole = "anchor",
    confidence: OrthogroupMemberConfidence = "high",
    assignment_reason: str = "",
    supporting_edges: Sequence[str] = (),
    best_core_support: float = 0.0,
    second_best_core_support: float = 0.0,
) -> OrthogroupMember:
    return OrthogroupMember(
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
        representative=representative,
        role=role,
        confidence=confidence,
        assignment_reason=assignment_reason,
        supporting_edges=tuple(str(edge_id) for edge_id in supporting_edges),
        best_core_support=float(best_core_support),
        second_best_core_support=float(second_best_core_support),
    )


def _member_ids_for_result(orthogroups: OrthogroupResult) -> dict[str, tuple[str, ...]]:
    return {
        orthogroup_id: tuple(str(member.protein_id) for member in members)
        for orthogroup_id, members in orthogroups.orthogroups.items()
    }


def _copy_orthogroup_result_with_metadata(
    result: OrthogroupResult,
    *,
    rbh_orthogroups: Mapping[str, Sequence[str]] | None = None,
    ortholog_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]] | None = None,
    ortholog_paths_by_orthogroup_id: Mapping[str, Sequence[OrthologPath]] | None = None,
    related_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]] | None = None,
) -> OrthogroupResult:
    return replace(
        result,
        rbh_orthogroups={
            str(key): tuple(str(protein_id) for protein_id in value)
            for key, value in (rbh_orthogroups or {}).items()
        },
        ortholog_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in (ortholog_edges_by_orthogroup_id or {}).items()
        },
        ortholog_paths_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in (ortholog_paths_by_orthogroup_id or {}).items()
        },
        related_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in (related_edges_by_orthogroup_id or {}).items()
        },
    )


def _orthogroup_result_from_member_ids(
    group_member_ids: Mapping[str, set[str]],
    representative_ids_by_group: Mapping[str, set[str]],
    protein_map: Mapping[str, CdsProtein],
    *,
    group_order: Sequence[str],
    rbh_orthogroups: Mapping[str, Sequence[str]],
    ortholog_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]],
    ortholog_paths_by_orthogroup_id: Mapping[str, Sequence[OrthologPath]],
    related_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]],
    member_roles_by_protein_id: Mapping[str, OrthogroupMemberRole] | None = None,
    member_confidence_by_protein_id: Mapping[str, OrthogroupMemberConfidence] | None = None,
    assignment_reason_by_protein_id: Mapping[str, str] | None = None,
    supporting_edges_by_protein_id: Mapping[str, Sequence[str]] | None = None,
    best_core_support_by_protein_id: Mapping[str, float] | None = None,
    second_best_core_support_by_protein_id: Mapping[str, float] | None = None,
) -> OrthogroupResult:
    orthogroups: dict[str, list[OrthogroupMember]] = {}
    member_by_protein_id: dict[str, OrthogroupMember] = {}
    ordered_group_ids = [
        group_id
        for group_id in group_order
        if group_id in group_member_ids and group_member_ids[group_id]
    ]
    remaining_group_ids = sorted(
        set(group_member_ids).difference(ordered_group_ids),
        key=lambda group_id: min(
            _protein_sort_key(protein_map[member_id])
            for member_id in group_member_ids[group_id]
            if member_id in protein_map
        ),
    )
    for group_id in [*ordered_group_ids, *remaining_group_ids]:
        member_ids = sorted(
            group_member_ids[group_id],
            key=lambda member_id: _protein_sort_key(protein_map[member_id]),
        )
        representative_ids = representative_ids_by_group.get(group_id, set())
        members: list[OrthogroupMember] = []
        for member_id in member_ids:
            protein = protein_map[member_id]
            member = _make_orthogroup_member(
                group_id,
                protein,
                representative=member_id in representative_ids,
                role=(member_roles_by_protein_id or {}).get(member_id, "anchor"),
                confidence=(member_confidence_by_protein_id or {}).get(member_id, "high"),
                assignment_reason=(assignment_reason_by_protein_id or {}).get(member_id, ""),
                supporting_edges=(supporting_edges_by_protein_id or {}).get(member_id, ()),
                best_core_support=(best_core_support_by_protein_id or {}).get(member_id, 0.0),
                second_best_core_support=(second_best_core_support_by_protein_id or {}).get(member_id, 0.0),
            )
            members.append(member)
            member_by_protein_id[member_id] = member
        orthogroups[group_id] = members

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
        rbh_orthogroups={
            str(key): tuple(str(protein_id) for protein_id in value)
            for key, value in rbh_orthogroups.items()
        },
        ortholog_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in ortholog_edges_by_orthogroup_id.items()
        },
        ortholog_paths_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in ortholog_paths_by_orthogroup_id.items()
        },
        related_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in related_edges_by_orthogroup_id.items()
        },
    )


def _canonical_edge_endpoint_ids(
    query_id: str,
    subject_id: str,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[str, str]:
    query_protein = protein_map[query_id]
    subject_protein = protein_map[subject_id]
    if int(query_protein.record_index) < int(subject_protein.record_index):
        return query_id, subject_id
    if int(query_protein.record_index) > int(subject_protein.record_index):
        return subject_id, query_id
    return tuple(sorted((query_id, subject_id)))  # type: ignore[return-value]


def _make_ortholog_edge(
    *,
    orthogroup_id: str,
    source_rbh_orthogroup_id: str | None,
    target_rbh_orthogroup_id: str | None,
    query_id: str,
    subject_id: str,
    row: object,
    protein_map: Mapping[str, CdsProtein],
    edge_kind: OrthologEdgeKind,
    render_role: OrthologRenderRole,
    path_id: str | None = None,
) -> OrthologEdge:
    query_protein = protein_map[query_id]
    subject_protein = protein_map[subject_id]
    return OrthologEdge(
        orthogroup_id=str(orthogroup_id),
        source_rbh_orthogroup_id=source_rbh_orthogroup_id,
        target_rbh_orthogroup_id=target_rbh_orthogroup_id,
        query_protein_id=str(query_id),
        subject_protein_id=str(subject_id),
        query_record_index=int(query_protein.record_index),
        subject_record_index=int(subject_protein.record_index),
        edge_kind=edge_kind,
        render_role=render_role,
        path_id=path_id,
        identity=_row_float(row, "identity", 0.0),
        evalue=_row_float(row, "evalue", 1.0),
        bitscore=_row_float(row, "bitscore", 0.0),
        alignment_length=int(_row_float(row, "alignment_length", 0.0)),
    )


def _collect_ranked_membership_edges(
    directional_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    max_hits: int,
) -> list[tuple[tuple[float, float, float, int, str, str], str, str, object]]:
    _validate_max_hits(max_hits, option_name="orthogroup_member_max_hits")
    best_by_pair: dict[tuple[str, str], tuple[tuple[float, float, float, int, str, str], str, str, object]] = {}
    for hits in directional_tables.values():
        if hits is None or hits.empty:
            continue
        capped_hits = cap_hits_per_query(hits, max_hits=int(max_hits), distinct_subjects=True)
        for row in capped_hits.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            if query_id not in protein_map or subject_id not in protein_map:
                continue
            if protein_map[query_id].record_index == protein_map[subject_id].record_index:
                continue
            canonical_query_id, canonical_subject_id = _canonical_edge_endpoint_ids(
                query_id,
                subject_id,
                protein_map,
            )
            key = (canonical_query_id, canonical_subject_id)
            rank = _hit_strength_key_from_row(row)
            current = best_by_pair.get(key)
            if current is None or rank < current[0]:
                best_by_pair[key] = (rank, canonical_query_id, canonical_subject_id, row)
    return sorted(best_by_pair.values(), key=lambda item: item[0])


def _rbh_edges_from_edge_tables(
    edge_tables: Mapping[tuple[int, int], DataFrame],
    seed_orthogroups: OrthogroupResult,
    protein_map: Mapping[str, CdsProtein],
) -> dict[str, list[OrthologEdge]]:
    edges_by_group: dict[str, list[OrthologEdge]] = {}
    for hits in edge_tables.values():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            if query_id not in protein_map or subject_id not in protein_map:
                continue
            canonical_query_id, canonical_subject_id = _canonical_edge_endpoint_ids(
                query_id,
                subject_id,
                protein_map,
            )
            query_member = seed_orthogroups.member_by_protein_id.get(canonical_query_id)
            subject_member = seed_orthogroups.member_by_protein_id.get(canonical_subject_id)
            if query_member is None or subject_member is None:
                continue
            if query_member.orthogroup_id != subject_member.orthogroup_id:
                continue
            group_id = str(query_member.orthogroup_id)
            edges_by_group.setdefault(group_id, []).append(
                _make_ortholog_edge(
                    orthogroup_id=group_id,
                    source_rbh_orthogroup_id=group_id,
                    target_rbh_orthogroup_id=group_id,
                    query_id=canonical_query_id,
                    subject_id=canonical_subject_id,
                    row=row,
                    protein_map=protein_map,
                    edge_kind="rbh",
                    render_role="block_anchor",
                )
            )
    return {
        group_id: sorted(edges, key=lambda edge: (edge.query_record_index, edge.subject_record_index, edge.query_protein_id, edge.subject_protein_id))
        for group_id, edges in edges_by_group.items()
    }


def _path_sort_key(
    protein_ids: Sequence[str],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[int, int, int, str], ...]:
    return tuple(_protein_sort_key(protein_map[protein_id]) for protein_id in protein_ids)


def _build_ortholog_paths(
    edges_by_group: Mapping[str, Sequence[OrthologEdge]],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[dict[str, tuple[OrthologEdge, ...]], dict[str, tuple[OrthologPath, ...]]]:
    updated_edges_by_group: dict[str, tuple[OrthologEdge, ...]] = {}
    paths_by_group: dict[str, tuple[OrthologPath, ...]] = {}
    for group_id, edges in edges_by_group.items():
        path_edges = [
            edge
            for edge in edges
            if edge.edge_kind in {"rbh", "coortholog"}
            and edge.query_protein_id in protein_map
            and edge.subject_protein_id in protein_map
        ]
        if not path_edges:
            updated_edges_by_group[group_id] = tuple(edges)
            paths_by_group[group_id] = ()
            continue
        outgoing: dict[str, list[OrthologEdge]] = {}
        incoming: dict[str, list[OrthologEdge]] = {}
        for edge in path_edges:
            outgoing.setdefault(edge.query_protein_id, []).append(edge)
            incoming.setdefault(edge.subject_protein_id, []).append(edge)
        for edge_list in outgoing.values():
            edge_list.sort(
                key=lambda edge: (
                    edge.subject_record_index,
                    _protein_sort_key(protein_map[edge.subject_protein_id]),
                    _edge_id(edge),
                )
            )
        nodes = set(outgoing).union(incoming)
        start_nodes = [
            node
            for node in nodes
            if node not in incoming
        ] or list(nodes)
        start_nodes.sort(key=lambda protein_id: _protein_sort_key(protein_map[protein_id]))

        raw_paths: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

        def walk(node: str, protein_path: tuple[str, ...], edge_path: tuple[str, ...]) -> None:
            next_edges = outgoing.get(node, [])
            if not next_edges:
                if edge_path:
                    raw_paths.append((protein_path, edge_path))
                return
            for edge in next_edges:
                if edge.subject_protein_id in protein_path:
                    if edge_path:
                        raw_paths.append((protein_path, edge_path))
                    continue
                walk(
                    edge.subject_protein_id,
                    (*protein_path, edge.subject_protein_id),
                    (*edge_path, _edge_id(edge)),
                )

        for start_node in start_nodes:
            walk(start_node, (start_node,), ())

        deduped: dict[tuple[str, ...], tuple[str, ...]] = {}
        for protein_path, edge_path in raw_paths:
            current = deduped.get(protein_path)
            if current is None or edge_path < current:
                deduped[protein_path] = edge_path
        sorted_paths = sorted(
            deduped.items(),
            key=lambda item: (_path_sort_key(item[0], protein_map), item[1]),
        )
        protein_path_counts: dict[str, int] = {}
        for protein_path, _edge_path in sorted_paths:
            for protein_id in set(protein_path):
                protein_path_counts[protein_id] = protein_path_counts.get(protein_id, 0) + 1

        edge_path_id: dict[str, str] = {}
        paths: list[OrthologPath] = []
        for path_index, (protein_path, edge_path) in enumerate(sorted_paths, start=1):
            path_id = f"{group_id}.path_{path_index}"
            for edge_id in edge_path:
                edge_path_id.setdefault(edge_id, path_id)
            shared_protein_ids = tuple(
                sorted(
                    (
                        protein_id
                        for protein_id in protein_path
                        if protein_path_counts.get(protein_id, 0) > 1
                    ),
                    key=lambda protein_id: _protein_sort_key(protein_map[protein_id]),
                )
            )
            paths.append(
                OrthologPath(
                    orthogroup_id=group_id,
                    path_id=path_id,
                    protein_ids=tuple(protein_path),
                    edge_ids=tuple(edge_path),
                    shared_protein_ids=shared_protein_ids,
                )
            )

        updated_edges_by_group[group_id] = tuple(
            replace(edge, path_id=edge_path_id.get(_edge_id(edge), edge.path_id))
            for edge in edges
        )
        paths_by_group[group_id] = tuple(paths)
    return updated_edges_by_group, paths_by_group


def _anchor_core_hit_rank(
    row: object,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[float, float, float, float, int, int, str, int, str]:
    query_id = str(getattr(row, "query", ""))
    subject_id = str(getattr(row, "subject", ""))
    query_record_index = int(protein_map[query_id].record_index) if query_id in protein_map else 0
    subject_record_index = int(protein_map[subject_id].record_index) if subject_id in protein_map else 0
    return (
        -_normalized_score_from_row(row),
        _row_float(row, "evalue", float("inf")),
        -_row_float(row, "min_coverage", 0.0),
        -_row_float(row, "identity", 0.0),
        -int(_row_float(row, "alignment_length", 0.0)),
        query_record_index,
        query_id,
        subject_record_index,
        subject_id,
    )


def _row_min_coverage(row: object) -> float:
    return _row_float(row, "min_coverage", 0.0)


def _row_domain_only(row: object) -> bool:
    return bool(getattr(row, "domain_only", False)) or _row_min_coverage(row) <= _DOMAIN_ONLY_MAX_MIN_COVERAGE


def _row_supports_membership(row: object) -> bool:
    return _row_min_coverage(row) >= _MIN_MEMBERSHIP_MIN_COVERAGE and not _row_domain_only(row)


def _dedupe_anchor_core_directional_rows(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> dict[tuple[str, str], object]:
    best_by_direction: dict[tuple[str, str], object] = {}
    for hits in normalized_tables.values():
        if hits is None or hits.empty:
            continue
        for row in hits.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            if query_id == subject_id:
                continue
            if query_id not in protein_map or subject_id not in protein_map:
                continue
            if _normalized_score_from_row(row) <= 0.0:
                continue
            key = (query_id, subject_id)
            current = best_by_direction.get(key)
            if current is None or _anchor_core_hit_rank(row, protein_map) < _anchor_core_hit_rank(current, protein_map):
                best_by_direction[key] = row
    return best_by_direction


def _best_rows_by_query_target_record(
    best_by_direction: Mapping[tuple[str, str], object],
    protein_map: Mapping[str, CdsProtein],
) -> dict[tuple[str, int], list[object]]:
    rows_by_query_record: dict[tuple[str, int], list[object]] = {}
    for (query_id, subject_id), row in best_by_direction.items():
        query_protein = protein_map.get(query_id)
        subject_protein = protein_map.get(subject_id)
        if query_protein is None or subject_protein is None:
            continue
        if int(query_protein.record_index) == int(subject_protein.record_index):
            continue
        rows_by_query_record.setdefault((query_id, int(subject_protein.record_index)), []).append(row)
    for rows in rows_by_query_record.values():
        rows.sort(key=lambda row: _anchor_core_hit_rank(row, protein_map))
    return rows_by_query_record


def _comparison_record_for_ids(row: object, query_id: str, subject_id: str) -> dict[str, object]:
    record = {column: getattr(row, column) for column in COMPARISON_COLUMNS}
    original_query = str(getattr(row, "query", ""))
    original_subject = str(getattr(row, "subject", ""))
    if original_query == subject_id and original_subject == query_id:
        record["qstart"], record["sstart"] = record["sstart"], record["qstart"]
        record["qend"], record["send"] = record["send"], record["qend"]
    record["query"] = query_id
    record["subject"] = subject_id
    return record


def _anchor_core_edge_sort_key(
    edge: _AnchorCoreEvidenceEdge,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[int, int, int, str], tuple[int, int, int, str], int, float, float, str, str]:
    return (
        _protein_sort_key(protein_map[edge.query_id]),
        _protein_sort_key(protein_map[edge.subject_id]),
        0 if edge.edge_kind == "rbh" else 1,
        -max(float(edge.score), float(edge.reverse_score)),
        _row_float(edge.row, "evalue", float("inf")),
        edge.query_id,
        edge.subject_id,
    )


def _select_anchor_core_edges(
    best_by_direction: Mapping[tuple[str, str], object],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[_AnchorCoreEvidenceEdge, ...]:
    rows_by_query_record = _best_rows_by_query_target_record(best_by_direction, protein_map)
    best_subject_by_query_record = {
        key: str(rows[0].subject)
        for key, rows in rows_by_query_record.items()
        if rows
    }
    best_score_by_query_record = {
        key: _normalized_score_from_row(rows[0])
        for key, rows in rows_by_query_record.items()
        if rows
    }
    selected_by_pair: dict[tuple[str, str], _AnchorCoreEvidenceEdge] = {}
    for (query_id, subject_id), row in sorted(
        best_by_direction.items(),
        key=lambda item: _anchor_core_hit_rank(item[1], protein_map),
    ):
        query_protein = protein_map.get(query_id)
        subject_protein = protein_map.get(subject_id)
        if query_protein is None or subject_protein is None:
            continue
        query_record_index = int(query_protein.record_index)
        subject_record_index = int(subject_protein.record_index)
        if query_record_index == subject_record_index:
            continue
        if not _row_supports_membership(row):
            continue
        reverse_row = best_by_direction.get((subject_id, query_id))
        if reverse_row is None or not _row_supports_membership(reverse_row):
            continue
        score = _normalized_score_from_row(row)
        reverse_score = _normalized_score_from_row(reverse_row)
        query_best = best_score_by_query_record.get((query_id, subject_record_index), 0.0)
        subject_best = best_score_by_query_record.get((subject_id, query_record_index), 0.0)
        if query_best <= 0.0 or subject_best <= 0.0:
            continue
        if score < _NEAR_RECIPROCAL_MIN_RATIO * query_best:
            continue
        if reverse_score < _NEAR_RECIPROCAL_MIN_RATIO * subject_best:
            continue
        strict_rbh = (
            best_subject_by_query_record.get((query_id, subject_record_index)) == subject_id
            and best_subject_by_query_record.get((subject_id, query_record_index)) == query_id
        )
        canonical_query_id, canonical_subject_id = _canonical_edge_endpoint_ids(
            query_id,
            subject_id,
            protein_map,
        )
        canonical_row = row if canonical_query_id == query_id else reverse_row
        canonical_score = score if canonical_query_id == query_id else reverse_score
        canonical_reverse_score = reverse_score if canonical_query_id == query_id else score
        candidate = _AnchorCoreEvidenceEdge(
            query_id=canonical_query_id,
            subject_id=canonical_subject_id,
            row=canonical_row,
            score=float(canonical_score),
            reverse_score=float(canonical_reverse_score),
            edge_kind="rbh" if strict_rbh else "coortholog",
        )
        key = (canonical_query_id, canonical_subject_id)
        current = selected_by_pair.get(key)
        if current is None:
            selected_by_pair[key] = candidate
            continue
        current_rank = _anchor_core_edge_sort_key(current, protein_map)
        candidate_rank = _anchor_core_edge_sort_key(candidate, protein_map)
        if candidate_rank < current_rank:
            selected_by_pair[key] = candidate
    return tuple(
        sorted(
            selected_by_pair.values(),
            key=lambda edge: _anchor_core_edge_sort_key(edge, protein_map),
        )
    )


def _derive_anchor_core_thresholds(
    best_by_direction: Mapping[tuple[str, str], object],
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> dict[str, _LocalThreshold]:
    cross_scores_by_protein: dict[str, list[float]] = {}
    for (query_id, subject_id), row in best_by_direction.items():
        query_protein = protein_map.get(query_id)
        subject_protein = protein_map.get(subject_id)
        if query_protein is None or subject_protein is None:
            continue
        if int(query_protein.record_index) == int(subject_protein.record_index):
            continue
        score = _normalized_score_from_row(row)
        if score > 0.0 and math.isfinite(score):
            cross_scores_by_protein.setdefault(query_id, []).append(score)

    anchor_scores_by_protein: dict[str, list[float]] = {}
    for edge in anchor_edges:
        if edge.score > 0.0:
            anchor_scores_by_protein.setdefault(edge.query_id, []).append(float(edge.score))
        if edge.reverse_score > 0.0:
            anchor_scores_by_protein.setdefault(edge.subject_id, []).append(float(edge.reverse_score))

    thresholds: dict[str, _LocalThreshold] = {}
    for protein_id in sorted(set(cross_scores_by_protein).union(anchor_scores_by_protein)):
        anchor_scores = [
            score
            for score in anchor_scores_by_protein.get(protein_id, [])
            if score > 0.0 and math.isfinite(score)
        ]
        if anchor_scores:
            anchor_floor = min(anchor_scores)
            thresholds[protein_id] = _LocalThreshold(
                protein_id=protein_id,
                score=float(anchor_floor * _INPARALOG_MIN_SAME_RECORD_RATIO_TO_LOCAL_ANCHOR),
                source="anchor",
                rbnh_count=len(anchor_scores),
                accepted_edge_count=len(anchor_scores),
            )
            continue
        cross_scores = [
            score
            for score in cross_scores_by_protein.get(protein_id, [])
            if score > 0.0 and math.isfinite(score)
        ]
        if not cross_scores:
            continue
        fallback_floor = max(cross_scores)
        thresholds[protein_id] = _LocalThreshold(
            protein_id=protein_id,
            score=float(fallback_floor * _INPARALOG_MIN_SAME_RECORD_RATIO_TO_LOCAL_ANCHOR),
            source="fallback",
            rbnh_count=0,
            accepted_edge_count=0,
        )
    return thresholds


def _best_evidence_between_protein_and_members(
    protein_id: str,
    member_ids: Sequence[str],
    best_by_direction: Mapping[tuple[str, str], object],
    protein_map: Mapping[str, CdsProtein],
    *,
    same_record: bool,
) -> tuple[float, object | None, str, str]:
    protein = protein_map[protein_id]
    best: tuple[float, object | None, str, str] = (0.0, None, "", "")
    member_set = set(member_ids)
    for member_id in member_set:
        if member_id == protein_id or member_id not in protein_map:
            continue
        member = protein_map[member_id]
        is_same_record = int(protein.record_index) == int(member.record_index)
        if is_same_record != bool(same_record):
            continue
        for query_id, subject_id in ((protein_id, member_id), (member_id, protein_id)):
            row = best_by_direction.get((query_id, subject_id))
            if row is None:
                continue
            score = _normalized_score_from_row(row)
            if score <= 0.0:
                continue
            current_row = best[1]
            if (
                current_row is None
                or score > best[0]
                or (
                    score == best[0]
                    and _anchor_core_hit_rank(row, protein_map) < _anchor_core_hit_rank(current_row, protein_map)
                )
            ):
                best = (float(score), row, query_id, subject_id)
    return best


def _support_gap_ratio(best_support: float, second_support: float) -> float:
    if best_support <= 0.0:
        return 0.0
    if second_support <= _ORTHOGROUP_SPLIT_EPSILON:
        return float("inf")
    return float(best_support) / max(float(second_support), _ORTHOGROUP_SPLIT_EPSILON)


def _build_core_support_candidate(
    protein_id: str,
    group_id: str,
    member_ids: Sequence[str],
    best_by_direction: Mapping[tuple[str, str], object],
    thresholds: Mapping[str, _LocalThreshold],
    protein_map: Mapping[str, CdsProtein],
) -> _CoreSupportCandidate | None:
    same_score, same_row, same_query_id, same_subject_id = _best_evidence_between_protein_and_members(
        protein_id,
        member_ids,
        best_by_direction,
        protein_map,
        same_record=True,
    )
    cross_score, cross_row, cross_query_id, cross_subject_id = _best_evidence_between_protein_and_members(
        protein_id,
        member_ids,
        best_by_direction,
        protein_map,
        same_record=False,
    )
    if same_row is None and cross_row is None:
        return None

    support = float(cross_score) + 0.5 * float(same_score)
    if support <= 0.0:
        return None

    candidate_threshold = thresholds.get(protein_id)
    same_member_id = ""
    if same_row is not None:
        same_member_id = same_subject_id if same_query_id == protein_id else same_query_id
    same_member_threshold = thresholds.get(same_member_id)
    same_pass = (
        same_row is not None
        and _row_supports_membership(same_row)
        and (
            (
                candidate_threshold is not None
                and same_score >= float(candidate_threshold.score)
            )
            or (
                same_member_threshold is not None
                and same_score >= float(same_member_threshold.score)
            )
        )
    )
    cross_pass = (
        cross_row is not None
        and _row_supports_membership(cross_row)
        and candidate_threshold is not None
        and cross_score >= float(candidate_threshold.score)
    )

    if same_pass or same_score >= cross_score:
        evidence_row = same_row
        evidence_query_id = same_query_id
        evidence_subject_id = same_subject_id
        relation_kind: OrthologEdgeKind = "same_record_inparalog"
    else:
        evidence_row = cross_row
        evidence_query_id = cross_query_id
        evidence_subject_id = cross_subject_id
        relation_kind = "coortholog"
    if evidence_row is None:
        return None

    min_coverage = _row_min_coverage(evidence_row)
    domain_only = _row_domain_only(evidence_row)
    high_confidence_pass = bool(same_pass or cross_pass) and min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    low_confidence_pass = min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    if same_pass:
        reason = "same-record support passed local anchor threshold"
    elif cross_pass:
        reason = "cross-record support passed local fallback threshold"
    elif domain_only:
        reason = "support is domain-only"
    elif min_coverage < _MIN_MEMBERSHIP_MIN_COVERAGE:
        reason = "support is below membership coverage"
    else:
        reason = "support is separated but below high-confidence threshold"
    return _CoreSupportCandidate(
        group_id=group_id,
        support=support,
        same_record_score=float(same_score),
        cross_record_score=float(cross_score),
        evidence_row=evidence_row,
        evidence_query_id=evidence_query_id,
        evidence_subject_id=evidence_subject_id,
        relation_kind=relation_kind,
        min_coverage=float(min_coverage),
        domain_only=bool(domain_only),
        high_confidence_pass=high_confidence_pass,
        low_confidence_pass=low_confidence_pass,
        reason=reason,
    )


def _core_support_sort_key(candidate: _CoreSupportCandidate) -> tuple[float, str]:
    return (-float(candidate.support), candidate.group_id)


def _append_limited_related_edge(
    related_edges_by_group: dict[str, list[OrthologEdge]],
    *,
    group_id: str,
    max_related_edges_per_orthogroup: int,
    edge: OrthologEdge,
) -> None:
    edges = related_edges_by_group.setdefault(group_id, [])
    if len(edges) >= int(max_related_edges_per_orthogroup):
        return
    existing = {
        (item.query_protein_id, item.subject_protein_id, item.edge_kind)
        for item in edges
    }
    key = (edge.query_protein_id, edge.subject_protein_id, edge.edge_kind)
    reverse_key = (edge.subject_protein_id, edge.query_protein_id, edge.edge_kind)
    if key in existing or reverse_key in existing:
        return
    edges.append(edge)


def _anchor_core_edge_tables(
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> dict[tuple[int, int], DataFrame]:
    rows_by_pair: dict[tuple[int, int], list[dict[str, object]]] = {}
    for edge in anchor_edges:
        query_record_index = int(protein_map[edge.query_id].record_index)
        subject_record_index = int(protein_map[edge.subject_id].record_index)
        if query_record_index == subject_record_index:
            continue
        pair = (min(query_record_index, subject_record_index), max(query_record_index, subject_record_index))
        query_id, subject_id = edge.query_id, edge.subject_id
        if query_record_index > subject_record_index:
            query_id, subject_id = subject_id, query_id
        rows_by_pair.setdefault(pair, []).append(_comparison_record_for_ids(edge.row, query_id, subject_id))
    return {
        pair: pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)
        for pair, rows in rows_by_pair.items()
    }


def _build_anchor_core_orthogroups(
    best_by_direction: Mapping[tuple[str, str], object],
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
    *,
    include_singletons: bool,
    max_related_edges_per_orthogroup: int,
) -> OrthogroupResult:
    thresholds = _derive_anchor_core_thresholds(best_by_direction, anchor_edges, protein_map)
    union_find = _UnionFind()
    member_ranks: dict[str, tuple[float, float, float, float]] = {}
    member_roles: dict[str, OrthogroupMemberRole] = {}
    member_confidence: dict[str, OrthogroupMemberConfidence] = {}
    assignment_reasons: dict[str, str] = {}
    supporting_edges: dict[str, tuple[str, ...]] = {}
    best_core_support: dict[str, float] = {}
    second_core_support: dict[str, float] = {}
    ortholog_edges_by_group: dict[str, list[OrthologEdge]] = {}

    for edge in anchor_edges:
        union_find.union(edge.query_id, edge.subject_id)
        rank = _member_rank_from_row(edge.row)
        if _is_better_member_rank(rank, member_ranks.get(edge.query_id)):
            member_ranks[edge.query_id] = rank
        if _is_better_member_rank(rank, member_ranks.get(edge.subject_id)):
            member_ranks[edge.subject_id] = rank
        for protein_id in (edge.query_id, edge.subject_id):
            current_role = member_roles.get(protein_id)
            if edge.edge_kind == "rbh" or current_role is None:
                member_roles[protein_id] = "anchor" if edge.edge_kind == "rbh" else "coortholog"
            member_confidence[protein_id] = "high"
            assignment_reasons[protein_id] = (
                "cross-record reciprocal best anchor"
                if edge.edge_kind == "rbh"
                else "cross-record near-reciprocal anchor"
            )
            supporting_edges[protein_id] = (*supporting_edges.get(protein_id, ()), f"{edge.query_id}->{edge.subject_id}:{edge.edge_kind}")

    components: dict[str, set[str]] = {}
    for edge in anchor_edges:
        root = union_find.find(edge.query_id)
        components.setdefault(root, set()).update((edge.query_id, edge.subject_id))
    if include_singletons:
        for protein_id in protein_map:
            root = union_find.find(str(protein_id))
            components.setdefault(root, set()).add(str(protein_id))

    sorted_components = sorted(
        (member_ids for member_ids in components.values() if member_ids),
        key=lambda member_ids: min(_protein_sort_key(protein_map[member_id]) for member_id in member_ids),
    )
    group_member_ids: dict[str, set[str]] = {
        f"og_{index}": set(member_ids)
        for index, member_ids in enumerate(sorted_components, start=1)
    }
    group_order = list(group_member_ids)
    group_by_protein: dict[str, str] = {
        protein_id: group_id
        for group_id, member_ids in group_member_ids.items()
        for protein_id in member_ids
    }
    representative_ids_by_group: dict[str, set[str]] = {}
    for group_id, member_ids in group_member_ids.items():
        members_by_record: dict[int, list[str]] = {}
        for member_id in member_ids:
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
            representative_ids_by_group.setdefault(group_id, set()).add(representative_id)

    for edge in anchor_edges:
        group_id = group_by_protein.get(edge.query_id)
        if group_id is None:
            continue
        ortholog_edges_by_group.setdefault(group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=group_id,
                source_rbh_orthogroup_id=group_id,
                target_rbh_orthogroup_id=group_id,
                query_id=edge.query_id,
                subject_id=edge.subject_id,
                row=edge.row,
                protein_map=protein_map,
                edge_kind=edge.edge_kind,
                render_role="block_anchor",
            )
        )

    core_member_snapshot = {
        group_id: set(member_ids)
        for group_id, member_ids in group_member_ids.items()
    }
    related_edges_by_group: dict[str, list[OrthologEdge]] = {}
    for protein_id in sorted(protein_map, key=lambda item: _protein_sort_key(protein_map[item])):
        if protein_id in group_by_protein:
            continue
        candidates = [
            candidate
            for group_id, member_ids in core_member_snapshot.items()
            if (
                candidate := _build_core_support_candidate(
                    protein_id,
                    group_id,
                    sorted(member_ids, key=lambda item: _protein_sort_key(protein_map[item])),
                    best_by_direction,
                    thresholds,
                    protein_map,
                )
            )
            is not None
        ]
        if not candidates:
            continue
        candidates.sort(key=_core_support_sort_key)
        best = candidates[0]
        second = candidates[1] if len(candidates) > 1 else None
        second_support = float(second.support) if second is not None else 0.0
        gap_ratio = _support_gap_ratio(best.support, second_support)
        assigned_role: OrthogroupMemberRole | None = None
        assigned_confidence: OrthogroupMemberConfidence | None = None
        assigned_reason = best.reason
        if best.high_confidence_pass and gap_ratio >= _MIN_BEST_SECOND_CORE_RATIO:
            assigned_role = "inparalog" if best.relation_kind == "same_record_inparalog" else "coortholog"
            assigned_confidence = "high"
        elif best.low_confidence_pass and gap_ratio >= _LOW_CONFIDENCE_MIN_BEST_SECOND_CORE_RATIO:
            assigned_role = "low_confidence"
            assigned_confidence = "low"
            assigned_reason = f"low-confidence assignment: {best.reason}"

        if assigned_role is None or assigned_confidence is None:
            if best.domain_only:
                unassigned_relation_kind: OrthologEdgeKind = "domain_only"
            elif second is not None and best.low_confidence_pass:
                unassigned_relation_kind = "ambiguous_paralog"
            else:
                unassigned_relation_kind = best.relation_kind
            if unassigned_relation_kind in {"ambiguous_paralog", "domain_only"}:
                _append_limited_related_edge(
                    related_edges_by_group,
                    group_id=best.group_id,
                    max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
                    edge=_make_ortholog_edge(
                        orthogroup_id=best.group_id,
                        source_rbh_orthogroup_id=best.group_id,
                        target_rbh_orthogroup_id=None,
                        query_id=best.evidence_query_id,
                        subject_id=best.evidence_subject_id,
                        row=best.evidence_row,
                        protein_map=protein_map,
                        edge_kind=unassigned_relation_kind,
                        render_role="display_edge",
                    ),
                )
            continue

        group_member_ids.setdefault(best.group_id, set()).add(protein_id)
        group_by_protein[protein_id] = best.group_id
        member_roles[protein_id] = assigned_role
        member_confidence[protein_id] = assigned_confidence
        assignment_reasons[protein_id] = assigned_reason
        best_core_support[protein_id] = float(best.support)
        second_core_support[protein_id] = float(second_support)
        support_edge_id = f"{best.evidence_query_id}->{best.evidence_subject_id}:{best.relation_kind}"
        supporting_edges[protein_id] = (support_edge_id,)
        ortholog_edges_by_group.setdefault(best.group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=best.group_id,
                source_rbh_orthogroup_id=best.group_id,
                target_rbh_orthogroup_id=best.group_id,
                query_id=best.evidence_query_id,
                subject_id=best.evidence_subject_id,
                row=best.evidence_row,
                protein_map=protein_map,
                edge_kind=best.relation_kind,
                render_role="display_edge",
            )
        )

    membership_edge_keys = {
        tuple(sorted((edge.query_protein_id, edge.subject_protein_id))) + (edge.edge_kind,)
        for edges in ortholog_edges_by_group.values()
        for edge in edges
    }
    seen_related_pairs: set[tuple[str, str, OrthologEdgeKind]] = set()
    for (query_id, subject_id), row in sorted(
        best_by_direction.items(),
        key=lambda item: _anchor_core_hit_rank(item[1], protein_map),
    ):
        if query_id not in protein_map or subject_id not in protein_map or query_id == subject_id:
            continue
        canonical_query_id, canonical_subject_id = _canonical_edge_endpoint_ids(query_id, subject_id, protein_map)
        query_group = group_by_protein.get(canonical_query_id)
        subject_group = group_by_protein.get(canonical_subject_id)
        if query_group is None and subject_group is None:
            continue
        relation_kind: OrthologEdgeKind
        if _row_domain_only(row):
            relation_kind = "domain_only"
        elif query_group is not None and subject_group is not None and query_group != subject_group:
            relation_kind = "weak_bridge"
        elif query_group == subject_group:
            query_role = member_roles.get(canonical_query_id)
            subject_role = member_roles.get(canonical_subject_id)
            if "inparalog" in {query_role, subject_role}:
                relation_kind = "same_record_inparalog"
            else:
                relation_kind = "related_homolog"
        else:
            relation_kind = "related_homolog"
        edge_key = tuple(sorted((canonical_query_id, canonical_subject_id))) + (relation_kind,)
        if edge_key in membership_edge_keys or edge_key in seen_related_pairs:
            continue
        seen_related_pairs.add(edge_key)
        for group_id in dict.fromkeys(group for group in (query_group, subject_group) if group is not None):
            _append_limited_related_edge(
                related_edges_by_group,
                group_id=group_id,
                max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
                edge=_make_ortholog_edge(
                    orthogroup_id=group_id,
                    source_rbh_orthogroup_id=query_group,
                    target_rbh_orthogroup_id=subject_group,
                    query_id=canonical_query_id,
                    subject_id=canonical_subject_id,
                    row=row,
                    protein_map=protein_map,
                    edge_kind=relation_kind,
                    render_role="display_edge",
                ),
            )

    updated_edges_by_group, paths_by_group = _build_ortholog_paths(
        {
            group_id: sorted(
                edges,
                key=lambda edge: (
                    edge.query_record_index,
                    edge.subject_record_index,
                    edge.query_protein_id,
                    edge.subject_protein_id,
                    edge.edge_kind,
                ),
            )
            for group_id, edges in ortholog_edges_by_group.items()
        },
        protein_map,
    )
    anchor_members_by_group = {
        group_id: [
            member_id
            for member_id in member_ids
            if member_roles.get(member_id) == "anchor"
        ]
        for group_id, member_ids in group_member_ids.items()
    }
    rbh_orthogroups = {
        group_id: tuple(
            sorted(
                anchor_members_by_group.get(group_id) or list(member_ids),
                key=lambda member_id: _protein_sort_key(protein_map[member_id]),
            )
        )
        for group_id, member_ids in group_member_ids.items()
    }
    return _orthogroup_result_from_member_ids(
        group_member_ids,
        representative_ids_by_group,
        protein_map,
        group_order=group_order,
        rbh_orthogroups=rbh_orthogroups,
        ortholog_edges_by_orthogroup_id=updated_edges_by_group,
        ortholog_paths_by_orthogroup_id=paths_by_group,
        related_edges_by_orthogroup_id={
            group_id: tuple(
                sorted(
                    edges,
                    key=lambda edge: (
                        ORTHOLOG_DISPLAY_EDGE_KIND_RANK.get(str(edge.edge_kind), 99),
                        edge.query_record_index,
                        edge.subject_record_index,
                        edge.query_protein_id,
                        edge.subject_protein_id,
                    ),
                )
            )
            for group_id, edges in related_edges_by_group.items()
        },
        member_roles_by_protein_id=member_roles,
        member_confidence_by_protein_id=member_confidence,
        assignment_reason_by_protein_id=assignment_reasons,
        supporting_edges_by_protein_id=supporting_edges,
        best_core_support_by_protein_id=best_core_support,
        second_best_core_support_by_protein_id=second_core_support,
    )


def _select_anchor_core_orthogroup_edges_from_directional_hits(
    directional_hits_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    record_count: int | None,
    include_singletons: bool,
    max_related_edges_per_orthogroup: int,
) -> OrthogroupEdgeSelectionResult:
    if int(max_related_edges_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    normalized_tables = _normalize_directional_hit_tables(
        directional_hits_by_pair,
        protein_map,
        min_coverage=0.0,
    )
    best_by_direction = _dedupe_anchor_core_directional_rows(normalized_tables, protein_map)
    anchor_edges = _select_anchor_core_edges(best_by_direction, protein_map)
    anchor_edge_tables = _anchor_core_edge_tables(anchor_edges, protein_map)
    if record_count is None:
        record_count = 0
        for query_index, subject_index in directional_hits_by_pair:
            record_count = max(record_count, int(query_index) + 1, int(subject_index) + 1)

    all_edges_by_pair = {
        pair: _comparison_columns_only(table)
        for pair, table in anchor_edge_tables.items()
    }
    adjacent_anchor_edges_by_pair = {
        pair: _comparison_columns_only(table)
        for pair, table in anchor_edge_tables.items()
        if int(pair[1]) == int(pair[0]) + 1
    }
    adjacent_candidate_edges_by_pair = {
        (query_index, query_index + 1): _comparison_columns_only(
            directional_hits_by_pair.get(
                (query_index, query_index + 1),
                _empty_comparison_hits(),
            )
        )
        for query_index in range(max(0, int(record_count) - 1))
    }
    for query_index in range(max(0, int(record_count) - 1)):
        adjacent_anchor_edges_by_pair.setdefault(
            (query_index, query_index + 1),
            _empty_comparison_hits(),
        )

    orthogroups = _build_anchor_core_orthogroups(
        best_by_direction,
        anchor_edges,
        protein_map,
        include_singletons=include_singletons,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )
    adjacent_display_edges_by_pair = _build_adjacent_display_edges_by_pair(
        adjacent_anchor_edges_by_pair,
        orthogroups,
        record_count=int(record_count),
        max_display_edges_per_orthogroup=int(max_related_edges_per_orthogroup),
        adjacent_candidate_edges_by_pair=adjacent_candidate_edges_by_pair,
    )
    return OrthogroupEdgeSelectionResult(
        orthogroups=orthogroups,
        all_edges_by_pair=all_edges_by_pair,
        adjacent_anchor_edges_by_pair=adjacent_anchor_edges_by_pair,
        adjacent_display_edges_by_pair=adjacent_display_edges_by_pair,
    )


def expand_orthogroup_membership_from_evidence(
    seed_orthogroups: OrthogroupResult,
    rbh_edge_tables: Mapping[tuple[int, int], DataFrame],
    directional_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    membership_mode: OrthogroupMembershipMode | str = ORTHOGROUP_INFERENCE_VERSION,
    member_max_hits: int = 5,
    max_related_edges_per_orthogroup: int = 2,
) -> OrthogroupResult:
    """Expand RBH seed orthogroups with strong family-level evidence."""

    normalized_mode = normalize_orthogroup_membership_mode(str(membership_mode))
    _validate_max_hits(member_max_hits, option_name="orthogroup_member_max_hits")
    if int(max_related_edges_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")

    if normalized_mode == ORTHOGROUP_INFERENCE_VERSION:
        normalized_tables = _normalize_directional_hit_tables(
            directional_tables,
            protein_map,
            min_coverage=0.0,
        )
        best_by_direction = _dedupe_anchor_core_directional_rows(normalized_tables, protein_map)
        anchor_edges = _select_anchor_core_edges(best_by_direction, protein_map)
        return _build_anchor_core_orthogroups(
            best_by_direction,
            anchor_edges,
            protein_map,
            include_singletons=set(seed_orthogroups.member_by_protein_id) == set(protein_map),
            max_related_edges_per_orthogroup=int(max_related_edges_per_orthogroup),
        )

    if normalized_mode == "distribution_split":
        normalized_tables = _normalize_directional_hit_tables(
            directional_tables,
            protein_map,
        )
        rbnh_edge_tables = _select_rbnh_edge_tables(normalized_tables)
        return _build_distribution_split_orthogroups_from_normalized_tables(
            normalized_tables,
            rbnh_edge_tables,
            protein_map,
            include_singletons=set(seed_orthogroups.member_by_protein_id) == set(protein_map),
            member_max_hits=int(member_max_hits),
        )

    rbh_orthogroups = _member_ids_for_result(seed_orthogroups)
    group_order = list(seed_orthogroups.orthogroups)
    group_member_ids: dict[str, set[str]] = {
        group_id: {str(member.protein_id) for member in members}
        for group_id, members in seed_orthogroups.orthogroups.items()
    }
    representative_ids_by_group: dict[str, set[str]] = {
        group_id: {
            str(member.protein_id)
            for member in members
            if bool(member.representative)
        }
        for group_id, members in seed_orthogroups.orthogroups.items()
    }
    seed_group_by_protein: dict[str, str] = {
        str(member.protein_id): group_id
        for group_id, members in seed_orthogroups.orthogroups.items()
        for member in members
    }
    final_group_by_protein: dict[str, str] = dict(seed_group_by_protein)

    rbh_edges_by_group = _rbh_edges_from_edge_tables(
        rbh_edge_tables,
        seed_orthogroups,
        protein_map,
    )
    if normalized_mode == "rbh":
        updated_edges_by_group, paths_by_group = _build_ortholog_paths(
            rbh_edges_by_group,
            protein_map,
        )
        return _copy_orthogroup_result_with_metadata(
            seed_orthogroups,
            rbh_orthogroups=rbh_orthogroups,
            ortholog_edges_by_orthogroup_id=updated_edges_by_group,
            ortholog_paths_by_orthogroup_id=paths_by_group,
            related_edges_by_orthogroup_id={},
        )

    related_edges_by_group: dict[str, list[OrthologEdge]] = {}
    ortholog_edges_by_group: dict[str, list[OrthologEdge]] = {
        group_id: list(edges)
        for group_id, edges in rbh_edges_by_group.items()
    }

    group_union = _UnionFind()
    for group_id in group_member_ids:
        group_union.add(group_id)

    def current_group(group_id: str) -> str:
        root = group_union.find(group_id)
        if root == group_id:
            return group_id
        return root

    def add_related_edge(
        query_id: str,
        subject_id: str,
        row: object,
        left_group: str | None,
        right_group: str | None,
    ) -> None:
        candidate_groups = list(dict.fromkeys(group for group in (left_group, right_group) if group is not None))
        for group_id in candidate_groups:
            if len(related_edges_by_group.get(group_id, [])) >= int(max_related_edges_per_orthogroup):
                continue
            related_edges_by_group.setdefault(group_id, []).append(
                _make_ortholog_edge(
                    orthogroup_id=group_id,
                    source_rbh_orthogroup_id=left_group,
                    target_rbh_orthogroup_id=right_group,
                    query_id=query_id,
                    subject_id=subject_id,
                    row=row,
                    protein_map=protein_map,
                    edge_kind="related_paralog",
                    render_role="display_edge",
                )
            )

    ranked_edges = _collect_ranked_membership_edges(
        directional_tables,
        protein_map,
        max_hits=int(member_max_hits),
    )
    for _rank, query_id, subject_id, row in ranked_edges:
        query_seed_group = seed_group_by_protein.get(query_id)
        subject_seed_group = seed_group_by_protein.get(subject_id)
        query_group = final_group_by_protein.get(query_id)
        subject_group = final_group_by_protein.get(subject_id)

        if query_seed_group and subject_seed_group:
            left_root = current_group(query_seed_group)
            right_root = current_group(subject_seed_group)
            if left_root != right_root:
                group_union.union(left_root, right_root)
                merged_root = current_group(left_root)
                merged_other = right_root if merged_root == left_root else left_root
                group_member_ids.setdefault(merged_root, set()).update(
                    group_member_ids.pop(merged_other, set())
                )
                representative_ids_by_group.setdefault(merged_root, set()).update(
                    representative_ids_by_group.pop(merged_other, set())
                )
                ortholog_edges_by_group.setdefault(merged_root, []).extend(
                    ortholog_edges_by_group.pop(merged_other, [])
                )
                related_edges_by_group.setdefault(merged_root, []).extend(
                    related_edges_by_group.pop(merged_other, [])
                )
                for protein_id, group_id in list(final_group_by_protein.items()):
                    if current_group(group_id) == merged_root:
                        final_group_by_protein[protein_id] = merged_root
                query_group = final_group_by_protein.get(query_id)
                subject_group = final_group_by_protein.get(subject_id)

        if query_group and subject_group:
            query_group = current_group(query_group)
            subject_group = current_group(subject_group)
            if query_group == subject_group:
                if query_seed_group and subject_seed_group and query_seed_group != subject_seed_group:
                    add_related_edge(query_id, subject_id, row, query_group, query_group)
                    continue
                existing_edge_ids = {
                    (edge.query_protein_id, edge.subject_protein_id, edge.edge_kind)
                    for edge in ortholog_edges_by_group.get(query_group, [])
                }
                edge_key = (query_id, subject_id, "coortholog")
                reverse_edge_key = (subject_id, query_id, "coortholog")
                rbh_key = (query_id, subject_id, "rbh")
                reverse_rbh_key = (subject_id, query_id, "rbh")
                if (
                    edge_key not in existing_edge_ids
                    and reverse_edge_key not in existing_edge_ids
                    and rbh_key not in existing_edge_ids
                    and reverse_rbh_key not in existing_edge_ids
                    and (query_seed_group != subject_seed_group or query_seed_group is None)
                ):
                    ortholog_edges_by_group.setdefault(query_group, []).append(
                        _make_ortholog_edge(
                            orthogroup_id=query_group,
                            source_rbh_orthogroup_id=query_seed_group,
                            target_rbh_orthogroup_id=subject_seed_group,
                            query_id=query_id,
                            subject_id=subject_id,
                            row=row,
                            protein_map=protein_map,
                            edge_kind="coortholog",
                            render_role="display_edge",
                        )
                    )
                continue
            add_related_edge(query_id, subject_id, row, query_group, subject_group)
            continue

        if query_group is None and subject_group is None:
            continue
        anchor_group = current_group(query_group or subject_group or "")
        unassigned_id = subject_id if query_group is not None else query_id
        group_member_ids.setdefault(anchor_group, set()).add(unassigned_id)
        final_group_by_protein[unassigned_id] = anchor_group
        ortholog_edges_by_group.setdefault(anchor_group, []).append(
            _make_ortholog_edge(
                orthogroup_id=anchor_group,
                source_rbh_orthogroup_id=query_seed_group,
                target_rbh_orthogroup_id=subject_seed_group,
                query_id=query_id,
                subject_id=subject_id,
                row=row,
                protein_map=protein_map,
                edge_kind="coortholog",
                render_role="display_edge",
            )
        )

    updated_edges_by_group, paths_by_group = _build_ortholog_paths(
        {
            group_id: sorted(
                edges,
                key=lambda edge: (
                    edge.query_record_index,
                    edge.subject_record_index,
                    edge.query_protein_id,
                    edge.subject_protein_id,
                    edge.edge_kind,
                ),
            )
            for group_id, edges in ortholog_edges_by_group.items()
        },
        protein_map,
    )

    return _orthogroup_result_from_member_ids(
        group_member_ids,
        representative_ids_by_group,
        protein_map,
        group_order=group_order,
        rbh_orthogroups=rbh_orthogroups,
        ortholog_edges_by_orthogroup_id=updated_edges_by_group,
        ortholog_paths_by_orthogroup_id=paths_by_group,
        related_edges_by_orthogroup_id={
            group_id: tuple(
                sorted(
                    edges[: int(max_related_edges_per_orthogroup)],
                    key=lambda edge: (
                        edge.query_record_index,
                        edge.subject_record_index,
                        edge.query_protein_id,
                        edge.subject_protein_id,
                    ),
                )
            )
            for group_id, edges in related_edges_by_group.items()
        },
    )


def build_orthogroups_from_protein_hits(
    hits_by_pair: Sequence[DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    include_singletons: bool = False,
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
    if include_singletons:
        for protein_id in protein_map:
            root = union_find.find(str(protein_id))
            components.setdefault(root, set()).add(str(protein_id))

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


def _orthogroup_member_count(
    orthogroups: OrthogroupResult | None,
    orthogroup_id: str,
    record_index: int,
) -> int:
    if orthogroups is None or not orthogroup_id:
        return 0
    return sum(
        1
        for member in orthogroups.orthogroups.get(orthogroup_id, [])
        if int(member.record_index) == int(record_index)
    )


def _edge_metadata_for_protein_pair(
    orthogroups: OrthogroupResult | None,
    orthogroup_id: str,
    query_id: str,
    subject_id: str,
) -> dict[str, object]:
    metadata = {
        "rbh_orthogroup_id": "",
        "ortholog_path_id": "",
        "edge_kind": "",
        "render_role": "",
    }
    if orthogroups is None or not orthogroup_id:
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
        edge_metadata = _edge_metadata_for_protein_pair(
            orthogroups,
            orthogroup_id,
            query_id,
            subject_id,
        )
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
                "rbh_orthogroup_id": edge_metadata["rbh_orthogroup_id"],
                "ortholog_path_id": edge_metadata["ortholog_path_id"],
                "edge_kind": edge_metadata["edge_kind"],
                "render_role": edge_metadata["render_role"],
                "query_orthogroup_representative": (
                    bool(query_member.representative) if query_member is not None else False
                ),
                "subject_orthogroup_representative": (
                    bool(subject_member.representative) if subject_member is not None else False
                ),
                "query_orthogroup_member_count": _orthogroup_member_count(
                    orthogroups,
                    orthogroup_id,
                    query_protein.record_index,
                ),
                "subject_orthogroup_member_count": _orthogroup_member_count(
                    orthogroups,
                    orthogroup_id,
                    subject_protein.record_index,
                ),
                "query_orthogroup_role": (
                    str(query_member.role) if query_member is not None else ""
                ),
                "subject_orthogroup_role": (
                    str(subject_member.role) if subject_member is not None else ""
                ),
                "query_orthogroup_confidence": (
                    str(query_member.confidence) if query_member is not None else ""
                ),
                "subject_orthogroup_confidence": (
                    str(subject_member.confidence) if subject_member is not None else ""
                ),
                "query_orthogroup_assignment_reason": (
                    str(query_member.assignment_reason) if query_member is not None else ""
                ),
                "subject_orthogroup_assignment_reason": (
                    str(subject_member.assignment_reason) if subject_member is not None else ""
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
    threads: int | None = None,
) -> DataFrame:
    """Run external LOSATP blastp and parse outfmt 6 output."""

    if max_hits is not None:
        _validate_max_hits(max_hits, option_name="protein_blastp_candidate_limit")
    _validate_losatp_threads(threads)
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
        if threads is not None:
            command.extend(["--num-threads", str(int(threads))])
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
    losatp_threads: int | None,
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
        threads=losatp_threads,
    )


def _empty_comparison_hits() -> DataFrame:
    return pd.DataFrame(columns=COMPARISON_COLUMNS)


def _ortholog_edge_display_rank(
    edge: OrthologEdge,
    orthogroups: OrthogroupResult,
) -> tuple[object, ...]:
    query_member = orthogroups.member_by_protein_id.get(edge.query_protein_id)
    subject_member = orthogroups.member_by_protein_id.get(edge.subject_protein_id)
    both_representative = (
        query_member is not None
        and subject_member is not None
        and bool(query_member.representative)
        and bool(subject_member.representative)
    )
    return (
        ORTHOLOG_DISPLAY_EDGE_KIND_RANK.get(str(edge.edge_kind), 99),
        0 if both_representative else 1,
        float(edge.evalue),
        -float(edge.bitscore),
        -float(edge.identity),
        -int(edge.alignment_length),
        int(edge.query_record_index),
        int(edge.subject_record_index),
        str(edge.query_protein_id),
        str(edge.subject_protein_id),
        str(edge.edge_kind),
    )


def _display_pair_orthogroup_id(
    query_id: str,
    subject_id: str,
    orthogroups: OrthogroupResult,
) -> str | None:
    query_member = orthogroups.member_by_protein_id.get(str(query_id))
    subject_member = orthogroups.member_by_protein_id.get(str(subject_id))
    if query_member is None or subject_member is None:
        return None
    if query_member.orthogroup_id != subject_member.orthogroup_id:
        return None
    return str(query_member.orthogroup_id)


def _comparison_row_from_ortholog_edge_for_pair(
    edge: OrthologEdge,
    pair: tuple[int, int],
) -> dict[str, object] | None:
    query_index, subject_index = int(pair[0]), int(pair[1])
    if (
        int(edge.query_record_index) == query_index
        and int(edge.subject_record_index) == subject_index
    ):
        query_id = edge.query_protein_id
        subject_id = edge.subject_protein_id
    elif (
        int(edge.query_record_index) == subject_index
        and int(edge.subject_record_index) == query_index
    ):
        query_id = edge.subject_protein_id
        subject_id = edge.query_protein_id
    else:
        return None

    alignment_length = max(0, int(edge.alignment_length))
    return {
        "query": str(query_id),
        "subject": str(subject_id),
        "identity": float(edge.identity),
        "alignment_length": alignment_length,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": 1,
        "qend": alignment_length,
        "sstart": 1,
        "send": alignment_length,
        "evalue": float(edge.evalue),
        "bitscore": float(edge.bitscore),
    }


def _comparison_row_key(row: Mapping[str, object]) -> tuple[str, str]:
    return str(row["query"]), str(row["subject"])


def _comparison_row_display_rank(
    row: Mapping[str, object],
    orthogroups: OrthogroupResult,
    *,
    edge_kind_rank: int,
) -> tuple[object, ...]:
    query_id, subject_id = _comparison_row_key(row)
    query_member = orthogroups.member_by_protein_id.get(query_id)
    subject_member = orthogroups.member_by_protein_id.get(subject_id)
    both_representative = (
        query_member is not None
        and subject_member is not None
        and bool(query_member.representative)
        and bool(subject_member.representative)
    )
    return (
        int(edge_kind_rank),
        0 if both_representative else 1,
        float(row.get("evalue", float("inf"))),
        -float(row.get("bitscore", 0.0)),
        -float(row.get("identity", 0.0)),
        -int(float(row.get("alignment_length", 0.0))),
        int(query_member.record_index) if query_member is not None else 0,
        int(subject_member.record_index) if subject_member is not None else 0,
        query_id,
        subject_id,
        "direct_orthogroup",
    )


def _comparison_row_from_hit_row(row: object) -> dict[str, object]:
    return {
        column: getattr(row, column)
        for column in COMPARISON_COLUMNS
    }


def _has_uncovered_display_alternative(
    row_key: tuple[str, str],
    *,
    current_kind_rank: int,
    candidate_subjects_by_query: Mapping[str, set[str]],
    candidate_queries_by_subject: Mapping[str, set[str]],
    candidate_kind_rank_by_pair: Mapping[tuple[str, str], int],
    covered_queries: set[str],
    covered_subjects: set[str],
    existing_pairs: set[tuple[str, str]],
) -> bool:
    query_id, subject_id = row_key
    if query_id not in covered_queries and subject_id in covered_subjects:
        for alternative_subject_id in candidate_subjects_by_query.get(query_id, set()):
            alternative_key = (query_id, alternative_subject_id)
            if alternative_subject_id == subject_id:
                continue
            if alternative_subject_id in covered_subjects:
                continue
            if alternative_key in existing_pairs or (alternative_key[1], alternative_key[0]) in existing_pairs:
                continue
            if candidate_kind_rank_by_pair.get(alternative_key, 99) > current_kind_rank:
                continue
            return True
    if query_id in covered_queries and subject_id not in covered_subjects:
        for alternative_query_id in candidate_queries_by_subject.get(subject_id, set()):
            alternative_key = (alternative_query_id, subject_id)
            if alternative_query_id == query_id:
                continue
            if alternative_query_id in covered_queries:
                continue
            if alternative_key in existing_pairs or (alternative_key[1], alternative_key[0]) in existing_pairs:
                continue
            if candidate_kind_rank_by_pair.get(alternative_key, 99) > current_kind_rank:
                continue
            return True
    return False


def _add_display_row(
    row: dict[str, object],
    *,
    existing_pairs: set[tuple[str, str]],
    covered_queries: set[str],
    covered_subjects: set[str],
    added_rows: list[dict[str, object]],
) -> None:
    row_key = _comparison_row_key(row)
    existing_pairs.add(row_key)
    existing_pairs.add((row_key[1], row_key[0]))
    covered_queries.add(row_key[0])
    covered_subjects.add(row_key[1])
    added_rows.append(row)


def _build_adjacent_display_edges_by_pair(
    adjacent_anchor_edges_by_pair: Mapping[tuple[int, int], DataFrame],
    orthogroups: OrthogroupResult,
    *,
    record_count: int,
    max_display_edges_per_orthogroup: int,
    adjacent_candidate_edges_by_pair: Mapping[tuple[int, int], DataFrame] | None = None,
) -> dict[tuple[int, int], DataFrame]:
    """Overlay bounded non-anchor orthogroup edges for adjacent display."""

    cap = int(max_display_edges_per_orthogroup)
    if cap <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")

    display_edges_by_pair: dict[tuple[int, int], DataFrame] = {
        (int(pair[0]), int(pair[1])): table.copy()
        for pair, table in adjacent_anchor_edges_by_pair.items()
    }
    for query_index in range(max(0, int(record_count) - 1)):
        display_edges_by_pair.setdefault(
            (query_index, query_index + 1),
            _empty_comparison_hits(),
        )

    all_secondary_edges_by_group = {
        orthogroup_id: tuple(
            edge
            for edge in (
                *orthogroups.ortholog_edges_by_orthogroup_id.get(orthogroup_id, ()),
                *orthogroups.related_edges_by_orthogroup_id.get(orthogroup_id, ()),
            )
            if edge.render_role == "display_edge"
        )
        for orthogroup_id in orthogroups.orthogroups
    }

    for pair, anchor_table in list(display_edges_by_pair.items()):
        existing_pairs = {
            (str(row.query), str(row.subject))
            for row in anchor_table.itertuples(index=False)
        }
        existing_pairs.update((subject_id, query_id) for query_id, subject_id in list(existing_pairs))
        # Secondary edges should reveal extra paralog endpoints, not redraw
        # cross-links between endpoints already covered by stronger edges.
        covered_query_by_group: dict[str, set[str]] = {}
        covered_subject_by_group: dict[str, set[str]] = {}
        for row in anchor_table.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            orthogroup_id = _display_pair_orthogroup_id(
                query_id,
                subject_id,
                orthogroups,
            )
            if orthogroup_id is None:
                continue
            covered_query_by_group.setdefault(orthogroup_id, set()).add(query_id)
            covered_subject_by_group.setdefault(orthogroup_id, set()).add(subject_id)
        added_rows: list[dict[str, object]] = []
        direct_candidates_by_group: dict[str, list[tuple[tuple[object, ...], dict[str, object]]]] = {}
        direct_candidate_table = (
            adjacent_candidate_edges_by_pair.get(pair)
            if adjacent_candidate_edges_by_pair is not None
            else None
        )
        if direct_candidate_table is not None and not direct_candidate_table.empty:
            direct_rows = _sort_hits_for_query_best(
                _coerce_outfmt6_numeric_columns(
                    direct_candidate_table.loc[:, list(COMPARISON_COLUMNS)]
                )
            ).drop_duplicates(["query", "subject"], keep="first")
            for row in direct_rows.itertuples(index=False):
                query_id = str(row.query)
                subject_id = str(row.subject)
                orthogroup_id = _display_pair_orthogroup_id(
                    query_id,
                    subject_id,
                    orthogroups,
                )
                if orthogroup_id is None:
                    continue
                row_record = _comparison_row_from_hit_row(row)
                direct_candidates_by_group.setdefault(orthogroup_id, []).append(
                    (
                        _comparison_row_display_rank(
                            row_record,
                            orthogroups,
                            edge_kind_rank=ORTHOLOG_DISPLAY_EDGE_KIND_RANK["coortholog"],
                        ),
                        row_record,
                    )
                )

        for orthogroup_id in sorted(set(all_secondary_edges_by_group).union(direct_candidates_by_group)):
            candidates: list[tuple[tuple[object, ...], dict[str, object]]] = []
            candidates.extend(direct_candidates_by_group.get(orthogroup_id, []))
            for edge in all_secondary_edges_by_group[orthogroup_id]:
                row = _comparison_row_from_ortholog_edge_for_pair(edge, pair)
                if row is None:
                    continue
                row_key = (str(row["query"]), str(row["subject"]))
                if row_key in existing_pairs or (row_key[1], row_key[0]) in existing_pairs:
                    continue
                candidates.append((_ortholog_edge_display_rank(edge, orthogroups), row))
            added_for_group = 0
            covered_queries = covered_query_by_group.setdefault(orthogroup_id, set())
            covered_subjects = covered_subject_by_group.setdefault(orthogroup_id, set())
            sorted_candidates = sorted(candidates, key=lambda item: item[0])
            candidate_subjects_by_query: dict[str, set[str]] = {}
            candidate_queries_by_subject: dict[str, set[str]] = {}
            candidate_kind_rank_by_pair: dict[tuple[str, str], int] = {}
            for rank, row in sorted_candidates:
                query_id, subject_id = _comparison_row_key(row)
                candidate_subjects_by_query.setdefault(query_id, set()).add(subject_id)
                candidate_queries_by_subject.setdefault(subject_id, set()).add(query_id)
                candidate_kind_rank_by_pair[(query_id, subject_id)] = int(rank[0])

            # Display selection is stricter than membership expansion: prefer
            # same-quality links that introduce new endpoints on both records.
            for require_both_uncovered in (True, False):
                for rank, row in sorted_candidates:
                    if added_for_group >= cap:
                        break
                    row_key = _comparison_row_key(row)
                    if row_key in existing_pairs or (row_key[1], row_key[0]) in existing_pairs:
                        continue
                    query_covered = row_key[0] in covered_queries
                    subject_covered = row_key[1] in covered_subjects
                    if query_covered and subject_covered:
                        continue
                    if require_both_uncovered:
                        if query_covered or subject_covered:
                            continue
                    else:
                        if not (query_covered or subject_covered):
                            continue
                        if _has_uncovered_display_alternative(
                            row_key,
                            current_kind_rank=int(rank[0]),
                            candidate_subjects_by_query=candidate_subjects_by_query,
                            candidate_queries_by_subject=candidate_queries_by_subject,
                            candidate_kind_rank_by_pair=candidate_kind_rank_by_pair,
                            covered_queries=covered_queries,
                            covered_subjects=covered_subjects,
                            existing_pairs=existing_pairs,
                        ):
                            continue
                    _add_display_row(
                        row,
                        existing_pairs=existing_pairs,
                        covered_queries=covered_queries,
                        covered_subjects=covered_subjects,
                        added_rows=added_rows,
                    )
                    added_for_group += 1
                if added_for_group >= cap:
                    break
        if added_rows:
            added_frame = pd.DataFrame.from_records(added_rows, columns=COMPARISON_COLUMNS)
            if anchor_table.empty:
                display_edges_by_pair[pair] = added_frame.reset_index(drop=True)
            else:
                display_edges_by_pair[pair] = pd.concat(
                    [
                        anchor_table.loc[:, list(COMPARISON_COLUMNS)],
                        added_frame,
                    ],
                    ignore_index=True,
                )
    return display_edges_by_pair


def select_rbh_orthogroup_edges_from_directional_hits(
    directional_hits_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    record_count: int | None = None,
    include_singletons: bool = False,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = ORTHOGROUP_INFERENCE_VERSION,
    orthogroup_member_max_hits: int = 5,
    max_related_edges_per_orthogroup: int = 2,
) -> OrthogroupEdgeSelectionResult:
    """Select RBH anchors and bounded expanded edges used by Orthogroup mode.

    Input tables are already threshold-filtered LOSATP outfmt6 rows keyed by
    directional record pair. Anchor and display edges are keyed by forward
    adjacent pairs, e.g. ``(0, 1)``. Expanded membership modes keep RBH
    anchors separate from additional non-anchor display edges.
    """

    normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    return _select_anchor_core_orthogroup_edges_from_directional_hits(
        directional_hits_by_pair,
        protein_map,
        record_count=record_count,
        include_singletons=include_singletons,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )

    unordered_pairs = sorted(
        {
            (min(int(query_index), int(subject_index)), max(int(query_index), int(subject_index)))
            for query_index, subject_index in directional_hits_by_pair
            if int(query_index) != int(subject_index)
        }
    )
    all_edges_by_pair: dict[tuple[int, int], DataFrame] = {}
    adjacent_anchor_edges_by_pair: dict[tuple[int, int], DataFrame] = {}
    for query_index, subject_index in unordered_pairs:
        forward_hits = directional_hits_by_pair.get((query_index, subject_index), _empty_comparison_hits())
        reverse_hits = directional_hits_by_pair.get((subject_index, query_index), _empty_comparison_hits())
        rbh_edges = select_reciprocal_best_hit_edges(forward_hits, reverse_hits)
        pair = (query_index, subject_index)
        all_edges_by_pair[pair] = rbh_edges
        if subject_index == query_index + 1:
            adjacent_anchor_edges_by_pair[pair] = rbh_edges

    if record_count is None:
        record_count = 0
        for query_index, subject_index in unordered_pairs:
            record_count = max(record_count, query_index + 1, subject_index + 1)
    for query_index in range(max(0, int(record_count) - 1)):
        adjacent_anchor_edges_by_pair.setdefault(
            (query_index, query_index + 1),
            _empty_comparison_hits(),
        )
    adjacent_candidate_edges_by_pair = {
        (query_index, query_index + 1): _comparison_columns_only(
            directional_hits_by_pair.get(
                (query_index, query_index + 1),
                _empty_comparison_hits(),
            )
        )
        for query_index in range(max(0, int(record_count) - 1))
    }

    seed_orthogroups = build_orthogroups_from_protein_hits(
        tuple(all_edges_by_pair.values()),
        protein_map,
        include_singletons=include_singletons,
    )
    orthogroups = expand_orthogroup_membership_from_evidence(
        seed_orthogroups,
        all_edges_by_pair,
        directional_hits_by_pair,
        protein_map,
        membership_mode=normalized_membership_mode,
        member_max_hits=orthogroup_member_max_hits,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )
    adjacent_display_edges_by_pair = _build_adjacent_display_edges_by_pair(
        adjacent_anchor_edges_by_pair,
        orthogroups,
        record_count=int(record_count),
        max_display_edges_per_orthogroup=int(max_related_edges_per_orthogroup),
        adjacent_candidate_edges_by_pair=adjacent_candidate_edges_by_pair,
    )
    return OrthogroupEdgeSelectionResult(
        orthogroups=orthogroups,
        all_edges_by_pair=all_edges_by_pair,
        adjacent_anchor_edges_by_pair=adjacent_anchor_edges_by_pair,
        adjacent_display_edges_by_pair=adjacent_display_edges_by_pair,
    )


def build_pairwise_protein_blastp_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    losatp_threads: int | None = None,
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
    _validate_losatp_threads(losatp_threads)
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
            losatp_threads=losatp_threads,
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
    losatp_threads: int | None = None,
    candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = ORTHOGROUP_INFERENCE_VERSION,
    orthogroup_member_max_hits: int = 5,
    max_related_edges_per_orthogroup: int = 2,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
) -> ProteinBlastpResult:
    """Infer all-vs-all RBH-seeded orthogroups and return adjacent display links."""

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='orthogroup' requires at least two records")
    _validate_losatp_threads(losatp_threads)
    _validate_candidate_limit(candidate_limit)
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    _validate_max_hits(orthogroup_member_max_hits, option_name="orthogroup_member_max_hits")
    if int(max_related_edges_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    if int(alignment_length) < 0:
        raise ValidationError("alignment_length must be >= 0")

    extraction = extract_cds_proteins(records)
    _validate_extraction_has_proteins(
        records,
        extraction,
        option_name="protein_blastp_mode='orthogroup'",
    )

    search_candidate_limit = candidate_limit
    directional_hits_by_pair: dict[tuple[int, int], DataFrame] = {}
    for query_index in range(len(records)):
        for subject_index in range(query_index, len(records)):
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
            filtered_forward_hits = filter_protein_hits_by_thresholds(
                forward_hits,
                evalue=evalue,
                bitscore=bitscore,
                identity=identity,
                alignment_length=alignment_length,
            )
            directional_hits_by_pair[(query_index, subject_index)] = filtered_forward_hits
            if query_index == subject_index:
                continue
            reverse_hits = _run_losatp_search(
                subject_fasta,
                query_fasta,
                losatp_bin=losatp_bin,
                losatp_threads=losatp_threads,
                candidate_limit=search_candidate_limit,
                runner=runner,
            )
            filtered_reverse_hits = filter_protein_hits_by_thresholds(
                reverse_hits,
                evalue=evalue,
                bitscore=bitscore,
                identity=identity,
                alignment_length=alignment_length,
            )
            directional_hits_by_pair[(subject_index, query_index)] = filtered_reverse_hits

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits_by_pair,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode=normalized_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )
    comparisons = [
        convert_protein_hits_to_genomic_links(
            edge_selection.adjacent_display_edges_by_pair.get(
                (record_index, record_index + 1),
                _empty_comparison_hits(),
            ),
            extraction.protein_map,
            orthogroups=edge_selection.orthogroups,
        )
        for record_index in range(len(records) - 1)
    ]
    return ProteinBlastpResult(comparisons=comparisons, orthogroups=edge_selection.orthogroups)


def build_protein_colinearity_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    losatp_threads: int | None = None,
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
        losatp_threads=losatp_threads,
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
    "ORTHOGROUP_INFERENCE_VERSION",
    "ORTHOGROUP_MEMBERSHIP_MODES",
    "PROTEIN_BLASTP_MODES",
    "EvidenceIndex",
    "OrthogroupEdgeSelectionResult",
    "OrthogroupMember",
    "OrthogroupMemberConfidence",
    "OrthogroupMemberRole",
    "OrthogroupMembershipMode",
    "OrthogroupNameCandidate",
    "OrthogroupResult",
    "OrthologEdge",
    "OrthologEdgeKind",
    "OrthologPath",
    "OrthologRenderRole",
    "ProteinBlastpMode",
    "ProteinBlastpResult",
    "ProteinExtractionResult",
    "build_pair_evidence_index",
    "build_orthogroups_from_protein_hits",
    "build_pairwise_protein_blastp_comparisons",
    "build_protein_colinearity_comparisons",
    "build_rbh_orthogroup_protein_blastp_comparisons",
    "cap_hits_per_query",
    "convert_pair_protein_hits_to_genomic_links",
    "convert_protein_hits_to_genomic_links",
    "extract_cds_proteins",
    "filter_protein_hits_by_thresholds",
    "expand_orthogroup_membership_from_evidence",
    "normalize_orthogroup_membership_mode",
    "normalize_protein_blastp_mode",
    "parse_losatp_outfmt6",
    "proteins_to_fasta",
    "run_losatp_blastp",
    "select_best_hits_per_query",
    "select_rbh_orthogroup_edges_from_directional_hits",
    "select_reciprocal_best_hit_edges",
    "select_reciprocal_best_hits",
    "select_top_hits_per_query",
]
