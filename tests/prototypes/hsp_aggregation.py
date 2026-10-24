"""S03 comparison oracle: pre-change aggregation from 9a4f7e29.

The pair loop and its row helpers are frozen here; the empty-table owner is the
unchanged production helper. This module must never be imported by production.
"""
from __future__ import annotations

import math
import sys
from typing import Mapping, Sequence, SupportsFloat, cast

import pandas as pd
from pandas import DataFrame

from gbdraw.analysis.protein_colinearity import CdsProtein, _empty_normalized_hit_table
from gbdraw.io.comparisons import COMPARISON_COLUMNS


def _row_float(row: object, column: str, default: float = 0.0) -> float:
    try:
        return float(getattr(row, column))
    except (TypeError, ValueError, AttributeError):
        return float(default)


def _coverage_coordinate(value: object) -> int | None:
    try:
        coordinate = float(cast(SupportsFloat, value))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(coordinate):
        return None
    return int(coordinate)


def _coverage_interval_from_hsp(start: object, end: object, length: int) -> tuple[int, int] | None:
    if int(length) <= 0:
        return None
    start_coordinate = _coverage_coordinate(start)
    end_coordinate = _coverage_coordinate(end)
    if start_coordinate is None or end_coordinate is None:
        return None
    low = min(start_coordinate, end_coordinate)
    high = max(start_coordinate, end_coordinate)
    clamped_low = max(1, low)
    clamped_high = min(int(length), high)
    if clamped_low > clamped_high:
        return None
    return clamped_low, clamped_high


def _merge_coverage_intervals(intervals: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    cleaned = sorted((int(start), int(end)) for start, end in intervals if int(start) <= int(end))
    if not cleaned:
        return ()
    merged: list[tuple[int, int]] = []
    current_start, current_end = cleaned[0]
    for start, end in cleaned[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
            continue
        merged.append((current_start, current_end))
        current_start, current_end = start, end
    merged.append((current_start, current_end))
    return tuple(merged)


def _covered_length(intervals: Sequence[tuple[int, int]]) -> int:
    return sum(max(0, int(end) - int(start) + 1) for start, end in _merge_coverage_intervals(intervals))


def _raw_hsp_coordinate_rank(row: object, column: str) -> int:
    coordinate = _coverage_coordinate(getattr(row, column, None))
    return coordinate if coordinate is not None else sys.maxsize


def _raw_hsp_representative_rank(row: object, row_index: int) -> tuple[float, float, float, int, int, int, int, int, int]:
    return (
        -_row_float(row, "bitscore", 0.0),
        _row_float(row, "evalue", float("inf")),
        -_row_float(row, "identity", 0.0),
        -int(_row_float(row, "alignment_length", 0.0)),
        _raw_hsp_coordinate_rank(row, "qstart"),
        _raw_hsp_coordinate_rank(row, "qend"),
        _raw_hsp_coordinate_rank(row, "sstart"),
        _raw_hsp_coordinate_rank(row, "send"),
        int(row_index),
    )


def aggregate_hsps_baseline(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
) -> DataFrame:
    if hits is None or hits.empty:
        columns = tuple(hits.columns) if hits is not None else tuple(COMPARISON_COLUMNS)
        return _empty_normalized_hit_table(columns)

    rows: list[dict[str, object]] = []
    for (query_id, subject_id), group in hits.groupby(["query", "subject"], sort=False):
        query_protein = protein_map.get(str(query_id))
        subject_protein = protein_map.get(str(subject_id))
        if query_protein is None or subject_protein is None:
            continue
        query_length = int(query_protein.protein_length)
        subject_length = int(subject_protein.protein_length)
        if query_length <= 0 or subject_length <= 0:
            continue
        hsp_rows = list(group.itertuples(index=False))
        if not hsp_rows:
            continue
        representative_index, representative_row = min(
            enumerate(hsp_rows),
            key=lambda item: _raw_hsp_representative_rank(item[1], item[0]),
        )
        bitscore = _row_float(representative_row, "bitscore", 0.0)
        representative_alignment_length = int(_row_float(representative_row, "alignment_length", 0.0))
        if bitscore <= 0.0 or representative_alignment_length <= 0:
            continue

        query_intervals: list[tuple[int, int]] = []
        subject_intervals: list[tuple[int, int]] = []
        total_hsp_alignment_length = 0
        for hsp_row in hsp_rows:
            alignment_length = _row_float(hsp_row, "alignment_length", 0.0)
            if math.isfinite(alignment_length) and alignment_length > 0:
                total_hsp_alignment_length += int(alignment_length)
            query_interval = _coverage_interval_from_hsp(
                getattr(hsp_row, "qstart", None),
                getattr(hsp_row, "qend", None),
                query_length,
            )
            if query_interval is not None:
                query_intervals.append(query_interval)
            subject_interval = _coverage_interval_from_hsp(
                getattr(hsp_row, "sstart", None),
                getattr(hsp_row, "send", None),
                subject_length,
            )
            if subject_interval is not None:
                subject_intervals.append(subject_interval)

        query_covered_length = _covered_length(query_intervals)
        subject_covered_length = _covered_length(subject_intervals)
        query_coverage = min(1.0, float(query_covered_length) / float(query_length))
        subject_coverage = min(1.0, float(subject_covered_length) / float(subject_length))
        min_hit_coverage = min(query_coverage, subject_coverage)

        record = dict(representative_row._asdict())
        record.update(
            {
                "query_length": query_length,
                "subject_length": subject_length,
                "length_product": float(query_length * subject_length),
                "hsp_count": int(len(hsp_rows)),
                "query_covered_length": int(query_covered_length),
                "subject_covered_length": int(subject_covered_length),
                "query_coverage": query_coverage,
                "subject_coverage": subject_coverage,
                "min_coverage": min_hit_coverage,
                "representative_alignment_length": int(representative_alignment_length),
                "total_hsp_alignment_length": int(total_hsp_alignment_length),
                "coverage_source": "hsp_union",
            }
        )
        rows.append(record)

    if not rows:
        return _empty_normalized_hit_table(tuple(hits.columns))
    return pd.DataFrame.from_records(rows)
