"""S03 comparison oracle: pre-change aggregation from 9a4f7e29.

Only the removed pair loop is frozen here. Rank/interval/empty-table owners are
unchanged production helpers; this module must never be imported by production.
"""
from __future__ import annotations

import math
from typing import Mapping

import pandas as pd
from pandas import DataFrame

from gbdraw.analysis.protein_colinearity import (
    CdsProtein, _empty_normalized_hit_table, _raw_hsp_representative_rank,
    _row_float, _coverage_interval_from_hsp, _covered_length,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS


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
