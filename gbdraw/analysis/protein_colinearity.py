#!/usr/bin/env python
# coding: utf-8

"""Protein-level colinearity helpers for linear pairwise comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import logging
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Callable, Mapping, Sequence

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
    start: int
    end: int
    strand: int | None
    feature_svg_id: str | None
    source_protein_id: str | None
    representative: bool = False


@dataclass(frozen=True)
class OrthogroupResult:
    """Connected-component orthogroups and per-protein lookup metadata."""

    orthogroups: dict[str, list[OrthogroupMember]]
    member_by_protein_id: dict[str, OrthogroupMember]


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


def _validate_max_hits(max_hits: int) -> None:
    if int(max_hits) <= 0:
        raise ValidationError("losatp_max_hits must be > 0")


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

    missing_columns = set(COMPARISON_COLUMNS).difference(hits.columns)
    if missing_columns:
        raise ParseError(
            "LOSATP blastp output is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    sorted_hits = hits.sort_values(
        ["query", "bitscore", "evalue", "identity", "alignment_length"],
        ascending=[True, False, True, False, False],
        kind="mergesort",
    )
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


def _member_rank_from_row(row: object) -> tuple[float, float, float]:
    return (
        _row_float(row, "bitscore", 0.0),
        _row_float(row, "evalue", float("inf")),
        _row_float(row, "identity", 0.0),
    )


def _is_better_member_rank(
    candidate: tuple[float, float, float],
    current: tuple[float, float, float] | None,
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
    )


def _protein_sort_key(protein: CdsProtein) -> tuple[int, int, int, str]:
    return (
        int(protein.record_index),
        int(protein.start),
        int(protein.end),
        str(protein.protein_id),
    )


def build_orthogroups_from_protein_hits(
    hits_by_pair: Sequence[DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> OrthogroupResult:
    """Build connected-component orthogroups from retained LOSATP protein hits."""

    union_find = _UnionFind()
    edge_nodes: list[tuple[str, str]] = []
    member_ranks: dict[str, tuple[float, float, float]] = {}
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
                    -member_ranks.get(member_id, (0.0, float("inf"), 0.0))[0],
                    member_ranks.get(member_id, (0.0, float("inf"), 0.0))[1],
                    -member_ranks.get(member_id, (0.0, float("inf"), 0.0))[2],
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
                start=protein.start,
                end=protein.end,
                strand=protein.strand,
                feature_svg_id=protein.feature_svg_id,
                source_protein_id=protein.source_protein_id,
                representative=member_id in representative_ids,
            )
            group_members.append(member)
            member_by_protein_id[member_id] = member
        orthogroups[orthogroup_id] = group_members

    return OrthogroupResult(
        orthogroups=orthogroups,
        member_by_protein_id=member_by_protein_id,
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
        elif query_member is not None:
            orthogroup_id = query_member.orthogroup_id
        elif subject_member is not None:
            orthogroup_id = subject_member.orthogroup_id
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
    max_hits: int = 5,
) -> DataFrame:
    """Run external LOSATP blastp and parse outfmt 6 output."""

    _validate_max_hits(max_hits)
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
            "-max_target_seqs",
            str(int(max_hits)),
            "-max_hsps_per_subject",
            "1",
        ]
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


def build_protein_colinearity_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    max_hits: int = 5,
    runner: LosatpRunner | None = None,
) -> list[DataFrame]:
    """Generate adjacent-record genomic comparison DataFrames with LOSATP blastp."""

    if len(records) < 2:
        raise ValidationError("protein_colinearity requires at least two records")
    _validate_max_hits(max_hits)

    extraction = extract_cds_proteins(records)
    empty_record_ids = [
        str(records[index].id)
        for index, proteins in enumerate(extraction.proteins_by_record)
        if not proteins
    ]
    if empty_record_ids:
        raise ValidationError(
            "protein_colinearity requires at least one CDS protein in each record; "
            "no CDS proteins were found in: "
            + ", ".join(empty_record_ids)
        )

    capped_hits_by_pair: list[DataFrame] = []
    for record_index in range(len(records) - 1):
        query_fasta = proteins_to_fasta(extraction.proteins_by_record[record_index])
        subject_fasta = proteins_to_fasta(extraction.proteins_by_record[record_index + 1])
        protein_hits = (
            runner(query_fasta, subject_fasta)
            if runner is not None
            else run_losatp_blastp(
                query_fasta,
                subject_fasta,
                losatp_bin=losatp_bin,
                max_hits=max_hits,
            )
        )
        capped_hits = cap_hits_per_query(protein_hits, max_hits=max_hits)
        capped_hits_by_pair.append(capped_hits)

    orthogroups = build_orthogroups_from_protein_hits(
        capped_hits_by_pair,
        extraction.protein_map,
    )
    comparisons = [
        convert_protein_hits_to_genomic_links(
            capped_hits,
            extraction.protein_map,
            orthogroups=orthogroups,
        )
        for capped_hits in capped_hits_by_pair
    ]
    return comparisons


__all__ = [
    "CdsProtein",
    "LOSATP_COMPARISON_COLUMNS",
    "LOSATP_METADATA_COLUMNS",
    "OrthogroupMember",
    "OrthogroupResult",
    "ProteinExtractionResult",
    "build_orthogroups_from_protein_hits",
    "build_protein_colinearity_comparisons",
    "cap_hits_per_query",
    "convert_pair_protein_hits_to_genomic_links",
    "convert_protein_hits_to_genomic_links",
    "extract_cds_proteins",
    "parse_losatp_outfmt6",
    "proteins_to_fasta",
    "run_losatp_blastp",
]
