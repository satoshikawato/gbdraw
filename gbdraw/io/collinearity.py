#!/usr/bin/env python
# coding: utf-8

"""Native gbdraw collinearity TSV import/export."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import csv
import math
from typing import Mapping, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    LosslessCollinearityParameters,
    _filter_blocks_by_min_anchors,
)
from gbdraw.analysis.collinearity_units import (
    CollinearityUnit,
    CollinearityUnitIndex,
    CollinearityUnitMode,
    build_collinearity_unit_index,
)
from gbdraw.analysis.protein_colinearity import CdsProtein, extract_cds_proteins
from gbdraw.exceptions import ParseError, ValidationError

NATIVE_COLLINEARITY_REQUIRED_COLUMNS = (
    "block_id",
    "query_record",
    "query_unit",
    "subject_record",
    "subject_unit",
    "orientation",
)
NATIVE_COLLINEARITY_WRITER_COLUMNS = (
    "block_id",
    "block_kind",
    "anchor_index",
    "orthogroup_id",
    "query_orthogroup_member_count",
    "subject_orthogroup_member_count",
    "query_record",
    "query_unit",
    "query_unit_kind",
    "query_locus_id",
    "query_protein_id",
    "query_display_name",
    "query_start",
    "query_end",
    "query_strand",
    "subject_record",
    "subject_unit",
    "subject_unit_kind",
    "subject_locus_id",
    "subject_protein_id",
    "subject_display_name",
    "subject_start",
    "subject_end",
    "subject_strand",
    "orientation",
    "identity",
    "evalue",
    "bitscore",
    "alignment_length",
    "score",
    "block_evalue",
)


def _read_text(source: str | Path) -> str:
    path = Path(str(source))
    try:
        if "\n" not in str(source) and path.exists():
            return path.read_text(encoding="utf-8")
    except OSError:
        pass
    return str(source)


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text == "."


def _optional_float(value: object, default: float) -> float:
    if _is_missing(value):
        return float(default)
    try:
        return float(str(value).strip())
    except ValueError as exc:
        raise ParseError(f"Invalid numeric value in native collinearity TSV: {value}") from exc


def _optional_nullable_float(value: object) -> float | None:
    if _is_missing(value):
        return None
    try:
        parsed = float(str(value).strip())
    except ValueError as exc:
        raise ParseError(f"Invalid numeric value in native collinearity TSV: {value}") from exc
    if not math.isfinite(parsed) or parsed < 0:
        raise ParseError(f"Invalid numeric value in native collinearity TSV: {value}")
    return parsed


def _optional_int(value: object, default: int) -> int:
    if _is_missing(value):
        return int(default)
    try:
        return int(float(str(value).strip()))
    except ValueError as exc:
        raise ParseError(f"Invalid integer value in native collinearity TSV: {value}") from exc


def _record_name_candidates(record: SeqRecord) -> set[str]:
    values = {str(record.id)}
    name = str(getattr(record, "name", "") or "").strip()
    if name:
        values.add(name)
    description = str(getattr(record, "description", "") or "").strip()
    if description:
        values.add(description)
    return {value for value in values if value}


def _resolve_record_selector(records: Sequence[SeqRecord], selector: object) -> int:
    text = str(selector or "").strip()
    if not text:
        raise ParseError("Record selector is empty in native collinearity TSV.")
    if text.startswith("#"):
        index_text = text[1:].strip()
        if not index_text.isdigit():
            raise ParseError(f"Invalid record selector '{text}'. Use #<number>, record ID, or record name.")
        index = int(index_text) - 1
        if index < 0 or index >= len(records):
            raise ParseError(f"Record selector '{text}' is out of range for {len(records)} records.")
        return index
    matches = [
        index
        for index, record in enumerate(records)
        if text in _record_name_candidates(record)
    ]
    if not matches:
        raise ParseError(f"Record selector '{text}' did not match any record.")
    if len(matches) > 1:
        raise ParseError(f"Record selector '{text}' matched multiple records; use #index.")
    return matches[0]


def _resolve_unit_selector(
    unit_index: CollinearityUnitIndex,
    record_index: int,
    selector: object,
) -> CollinearityUnit:
    text = str(selector or "").strip()
    if not text:
        raise ParseError("Unit selector is empty in native collinearity TSV.")
    if record_index >= len(unit_index.aliases_by_record):
        raise ParseError(f"Record index #{record_index + 1} has no collinearity units.")
    if text in unit_index.ambiguous_aliases_by_record[record_index]:
        raise ParseError(
            f"Unit selector '{text}' is ambiguous within record #{record_index + 1}."
        )
    unit_id = unit_index.aliases_by_record[record_index].get(text)
    if not unit_id:
        raise ParseError(
            f"Unit selector '{text}' did not match any collinearity unit in record #{record_index + 1}."
        )
    return unit_index.unit_by_id[unit_id]


def _protein_coordinates(unit: CollinearityUnit, protein_map: Mapping[str, CdsProtein]) -> tuple[int, int]:
    protein = protein_map.get(unit.representative_protein_id)
    if protein is None:
        if unit.strand == -1:
            return int(unit.end), int(unit.start) + 1
        return int(unit.start) + 1, int(unit.end)
    if protein.strand == -1:
        return int(protein.end), int(protein.start) + 1
    return int(protein.start) + 1, int(protein.end)


def _normalize_header(row: Mapping[str, str]) -> dict[str, str]:
    normalized = {str(key).strip(): value for key, value in row.items() if key is not None}
    if "query_unit" not in normalized and "query_cds" in normalized:
        normalized["query_unit"] = normalized["query_cds"]
    if "subject_unit" not in normalized and "subject_cds" in normalized:
        normalized["subject_unit"] = normalized["subject_cds"]
    return normalized


def _normalize_block_kind(value: object) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return "syntenic"
    if text not in {"syntenic", "cluster", "singleton"}:
        raise ParseError(f"Invalid native collinearity block_kind '{text}'.")
    return text


def _parse_rows(text: str) -> list[dict[str, str]]:
    lines = [
        line
        for line in str(text).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise ParseError("Native collinearity TSV is empty.")
    reader = csv.DictReader(lines, delimiter="\t")
    if not reader.fieldnames:
        raise ParseError("Native collinearity TSV must include a header row.")
    rows = [_normalize_header(row) for row in reader]
    if not rows:
        return []
    missing = set(NATIVE_COLLINEARITY_REQUIRED_COLUMNS).difference(rows[0].keys())
    if missing:
        raise ParseError(
            "Native collinearity TSV is missing required columns: "
            + ", ".join(sorted(missing))
        )
    return rows


def parse_native_collinearity_tsv(
    source: str | Path,
    records: Sequence[SeqRecord],
    *,
    params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    unit_mode: CollinearityUnitMode | str = "auto",
) -> CollinearityResult:
    """Parse and validate native headered collinearity TSV."""

    if len(records) < 2:
        raise ValidationError("Native collinearity TSV import requires at least two records.")
    params = params or LosslessCollinearityParameters()
    params.validate()
    min_anchors = int(params.min_anchors)
    default_anchor_score = float(getattr(params, "constant_anchor_score", 50.0))
    rows = _parse_rows(_read_text(source))
    extraction = extract_cds_proteins(records)
    unit_index = build_collinearity_unit_index(
        extraction,
        records=records,
        mode=unit_mode,
    )

    grouped: dict[str, list[tuple[int, CollinearityAnchor, float]]] = {}
    block_record_pairs: dict[str, tuple[int, int]] = {}
    block_orientations: dict[str, str] = {}
    block_kinds: dict[str, str] = {}
    block_evalues: dict[str, float | None] = {}
    for row_number, row in enumerate(rows, start=2):
        block_id = str(row.get("block_id", "")).strip()
        if not block_id:
            raise ParseError(f"Missing block_id on native collinearity TSV row {row_number}.")
        query_record_index = _resolve_record_selector(records, row.get("query_record"))
        subject_record_index = _resolve_record_selector(records, row.get("subject_record"))
        if subject_record_index != query_record_index + 1:
            raise ValidationError(
                f"Native collinearity TSV block '{block_id}' references non-adjacent records; "
                "MVP rendering supports query #N to subject #N+1 only."
            )
        pair = (query_record_index, subject_record_index)
        existing_pair = block_record_pairs.setdefault(block_id, pair)
        if existing_pair != pair:
            raise ValidationError(f"Native collinearity TSV block '{block_id}' has mixed record pairs.")

        orientation = str(row.get("orientation", "")).strip().lower()
        if orientation not in {"plus", "minus"}:
            raise ParseError(f"Invalid orientation '{orientation}' in block '{block_id}'.")
        existing_orientation = block_orientations.setdefault(block_id, orientation)
        if existing_orientation != orientation:
            raise ValidationError(f"Native collinearity TSV block '{block_id}' has mixed orientations.")
        block_kind = _normalize_block_kind(row.get("block_kind"))
        existing_kind = block_kinds.setdefault(block_id, block_kind)
        if existing_kind != block_kind:
            raise ValidationError(f"Native collinearity TSV block '{block_id}' has mixed block_kind values.")
        block_evalue = _optional_nullable_float(row.get("block_evalue"))
        existing_block_evalue = block_evalues.setdefault(block_id, block_evalue)
        if existing_block_evalue != block_evalue:
            raise ValidationError(f"Native collinearity TSV block '{block_id}' has conflicting block_evalue values.")

        query_unit = _resolve_unit_selector(unit_index, query_record_index, row.get("query_unit"))
        subject_unit = _resolve_unit_selector(unit_index, subject_record_index, row.get("subject_unit"))
        default_qstart, default_qend = _protein_coordinates(query_unit, extraction.protein_map)
        default_sstart, default_send = _protein_coordinates(subject_unit, extraction.protein_map)
        qstart = _optional_int(row.get("query_start"), default_qstart)
        qend = _optional_int(row.get("query_end"), default_qend)
        sstart = _optional_int(row.get("subject_start"), default_sstart)
        send = _optional_int(row.get("subject_end"), default_send)
        score = _optional_float(row.get("score"), default_anchor_score)
        anchor = CollinearityAnchor(
            query_protein_id=query_unit.representative_protein_id,
            subject_protein_id=subject_unit.representative_protein_id,
            query_record_index=query_record_index,
            subject_record_index=subject_record_index,
            query_order=query_unit.order,
            subject_order=subject_unit.order,
            query_start=qstart,
            query_end=qend,
            subject_start=sstart,
            subject_end=send,
            query_strand=query_unit.strand,
            subject_strand=subject_unit.strand,
            identity=_optional_float(row.get("identity"), 0.0),
            evalue=_optional_float(row.get("evalue"), 0.0),
            bitscore=_optional_float(row.get("bitscore"), score),
            alignment_length=_optional_int(row.get("alignment_length"), 0),
            query_feature_svg_id=query_unit.representative_feature_svg_id,
            subject_feature_svg_id=subject_unit.representative_feature_svg_id,
            source="native_tsv",
            query_unit_id=query_unit.unit_id,
            subject_unit_id=subject_unit.unit_id,
            query_unit_kind=query_unit.unit_kind,
            subject_unit_kind=subject_unit.unit_kind,
            query_locus_id=query_unit.locus_id,
            subject_locus_id=subject_unit.locus_id,
            query_display_name=query_unit.display_name,
            subject_display_name=subject_unit.display_name,
            orthogroup_id=str(row.get("orthogroup_id", "") or "").strip(),
            query_orthogroup_member_count=_optional_int(row.get("query_orthogroup_member_count"), 0),
            subject_orthogroup_member_count=_optional_int(row.get("subject_orthogroup_member_count"), 0),
        )
        anchor_index = _optional_int(row.get("anchor_index"), len(grouped.get(block_id, [])) + 1)
        grouped.setdefault(block_id, []).append((anchor_index, anchor, score))

    blocks: list[CollinearityBlock] = []
    for block_id in sorted(grouped):
        ordered_entries = sorted(
            grouped[block_id],
            key=lambda item: (
                item[0],
                item[1].query_order,
                item[1].subject_order,
                item[1].query_protein_id,
                item[1].subject_protein_id,
            ),
        )
        anchors = tuple(item[1] for item in ordered_entries)
        block_kind = block_kinds.get(block_id, "syntenic")
        if block_kind == "singleton" and len(anchors) != 1:
            raise ValidationError(
                f"Native collinearity TSV singleton block '{block_id}' must contain exactly one anchor."
            )
        score = float(sum(item[2] for item in ordered_entries))
        pair = block_record_pairs[block_id]
        blocks.append(
            CollinearityBlock(
                block_id=block_id,
                query_record_index=pair[0],
                subject_record_index=pair[1],
                orientation=block_orientations[block_id],  # type: ignore[arg-type]
                score=score,
                kind=block_kind,  # type: ignore[arg-type]
                block_evalue=block_evalues.get(block_id),
                anchors=anchors,
            )
        )
    filtered_blocks, unblocked_anchors = _filter_blocks_by_min_anchors(
        blocks,
        min_anchors=min_anchors,
    )
    return CollinearityResult(blocks=filtered_blocks, unblocked_anchors=unblocked_anchors)


def write_native_collinearity_tsv(result: CollinearityResult) -> str:
    """Serialize native collinearity blocks as headered TSV text."""

    handle = StringIO()
    writer = csv.DictWriter(
        handle,
        fieldnames=list(NATIVE_COLLINEARITY_WRITER_COLUMNS),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for block in result.blocks:
        for anchor_index, anchor in enumerate(block.anchors, start=1):
            writer.writerow(
                {
                    "block_id": block.block_id,
                    "block_kind": block.kind,
                    "anchor_index": anchor_index,
                    "orthogroup_id": anchor.orthogroup_id,
                    "query_orthogroup_member_count": anchor.query_orthogroup_member_count,
                    "subject_orthogroup_member_count": anchor.subject_orthogroup_member_count,
                    "query_record": f"#{anchor.query_record_index + 1}",
                    "query_unit": anchor.query_unit_id,
                    "query_unit_kind": anchor.query_unit_kind,
                    "query_locus_id": anchor.query_locus_id or "",
                    "query_protein_id": anchor.query_protein_id,
                    "query_display_name": anchor.query_display_name,
                    "query_start": anchor.query_start,
                    "query_end": anchor.query_end,
                    "query_strand": anchor.query_strand or "",
                    "subject_record": f"#{anchor.subject_record_index + 1}",
                    "subject_unit": anchor.subject_unit_id,
                    "subject_unit_kind": anchor.subject_unit_kind,
                    "subject_locus_id": anchor.subject_locus_id or "",
                    "subject_protein_id": anchor.subject_protein_id,
                    "subject_display_name": anchor.subject_display_name,
                    "subject_start": anchor.subject_start,
                    "subject_end": anchor.subject_end,
                    "subject_strand": anchor.subject_strand or "",
                    "orientation": block.orientation,
                    "identity": anchor.identity,
                    "evalue": anchor.evalue,
                    "bitscore": anchor.bitscore,
                    "alignment_length": anchor.alignment_length,
                    "score": (
                        anchor.bitscore
                        if pd.notna(anchor.bitscore)
                        else ""
                    ),
                    "block_evalue": (
                        block.block_evalue
                        if block.block_evalue is not None
                        else "."
                    ),
                }
            )
    return handle.getvalue()


def write_native_collinearity_tsv_file(path: str | Path, result: CollinearityResult) -> None:
    Path(path).write_text(write_native_collinearity_tsv(result), encoding="utf-8")


__all__ = [
    "NATIVE_COLLINEARITY_REQUIRED_COLUMNS",
    "NATIVE_COLLINEARITY_WRITER_COLUMNS",
    "parse_native_collinearity_tsv",
    "write_native_collinearity_tsv",
    "write_native_collinearity_tsv_file",
]
