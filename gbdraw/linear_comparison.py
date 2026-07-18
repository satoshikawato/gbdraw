"""Normalized comparison contracts for Linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError


@dataclass(frozen=True)
class LinearComparison:
    """A comparison result with explicit input-record endpoints."""

    query_record_index: int
    subject_record_index: int
    matches: DataFrame

    def __post_init__(self) -> None:
        for name in ("query_record_index", "subject_record_index"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValidationError(f"{name} must be a non-negative integer.")
        if self.query_record_index == self.subject_record_index:
            raise ValidationError("A Linear comparison must connect two different records.")
        if not isinstance(self.matches, DataFrame):
            raise ValidationError("LinearComparison.matches must be a pandas DataFrame.")


def merge_linear_comparisons(
    comparisons: Sequence[LinearComparison],
) -> tuple[LinearComparison, ...]:
    """Merge multiple sources for the same directed endpoint pair."""

    grouped: dict[tuple[int, int], list[DataFrame]] = {}
    order: list[tuple[int, int]] = []
    for comparison in comparisons:
        if not isinstance(comparison, LinearComparison):
            raise ValidationError("linear_comparisons must contain LinearComparison values.")
        key = (comparison.query_record_index, comparison.subject_record_index)
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(comparison.matches)
    return tuple(
        LinearComparison(
            query_record_index=query_index,
            subject_record_index=subject_index,
            matches=(
                frames[0].reset_index(drop=True)
                if len(frames) == 1
                else pd.concat(frames, ignore_index=True)
            ),
        )
        for (query_index, subject_index) in order
        for frames in (grouped[(query_index, subject_index)],)
    )


def validate_linear_comparison_topology(
    comparisons: Sequence[LinearComparison],
    rows_by_record: Sequence[int],
) -> None:
    """Reject same-row, same-record, and non-adjacent-row comparison edges."""

    record_count = len(rows_by_record)
    for comparison in comparisons:
        query_index = comparison.query_record_index
        subject_index = comparison.subject_record_index
        if query_index >= record_count or subject_index >= record_count:
            raise ValidationError(
                "Linear comparison endpoint is outside the loaded records: "
                f"query=#{query_index + 1}, subject=#{subject_index + 1}, "
                f"record_count={record_count}."
            )
        query_row = int(rows_by_record[query_index])
        subject_row = int(rows_by_record[subject_index])
        if query_row == subject_row:
            raise ValidationError(
                "Linear comparison endpoints must be in different rows: "
                f"query=#{query_index + 1} row={query_row + 1}, "
                f"subject=#{subject_index + 1} row={subject_row + 1}."
            )
        if abs(query_row - subject_row) != 1:
            raise ValidationError(
                "Linear comparison endpoints must be in adjacent rows: "
                f"query=#{query_index + 1} row={query_row + 1}, "
                f"subject=#{subject_index + 1} row={subject_row + 1}."
            )


__all__ = [
    "LinearComparison",
    "merge_linear_comparisons",
    "validate_linear_comparison_topology",
]
