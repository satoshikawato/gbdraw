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
from gbdraw.io.comparisons import COMPARISON_COLUMNS

logger = logging.getLogger(__name__)

_VALID_PROTEIN_RE = re.compile(r"^[A-Z]+$")
_VALID_FASTA_ID_RE = re.compile(r"^\S+$")
_NUMERIC_COMPARISON_COLUMNS = COMPARISON_COLUMNS[2:]


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


@dataclass(frozen=True)
class ProteinExtractionResult:
    """CDS protein extraction output grouped for pairwise LOSATP runs."""

    proteins_by_record: list[list[CdsProtein]]
    protein_map: dict[str, CdsProtein]


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
                if fasta_safe_source_id is not None
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


def _genomic_link_coordinates(protein: CdsProtein) -> tuple[int, int]:
    if protein.strand == -1:
        return protein.end, protein.start + 1
    return protein.start + 1, protein.end


def convert_protein_hits_to_genomic_links(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
) -> DataFrame:
    """Convert protein hit rows to genomic-coordinate comparison rows."""

    return convert_pair_protein_hits_to_genomic_links(hits, protein_map, protein_map)


def convert_pair_protein_hits_to_genomic_links(
    hits: DataFrame,
    query_protein_map: Mapping[str, CdsProtein],
    subject_protein_map: Mapping[str, CdsProtein],
) -> DataFrame:
    """Convert pairwise protein hit rows using separate query and subject maps."""

    if hits.empty:
        return pd.DataFrame(columns=COMPARISON_COLUMNS)

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
            }
        )

    if missing_ids:
        raise ParseError(
            "LOSATP blastp output contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )

    return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)


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

    comparisons: list[DataFrame] = []
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
        comparisons.append(
            convert_protein_hits_to_genomic_links(capped_hits, extraction.protein_map)
        )
    return comparisons


__all__ = [
    "CdsProtein",
    "ProteinExtractionResult",
    "build_protein_colinearity_comparisons",
    "cap_hits_per_query",
    "convert_pair_protein_hits_to_genomic_links",
    "convert_protein_hits_to_genomic_links",
    "extract_cds_proteins",
    "parse_losatp_outfmt6",
    "proteins_to_fasta",
    "run_losatp_blastp",
]
