"""Normalized comparison contracts for Linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.layout.record_coordinates import RecordDisplayTransform, SourceInterval, alignment_cut_breakpoints


def project_match_endpoints(
    endpoints: tuple[int, int, int, int],
    transforms: tuple[RecordDisplayTransform, RecordDisplayTransform],
) -> tuple[tuple[float, float, float, float], ...]:
    """Split existing record-local, directed inclusive endpoints at common t cuts.

    The source HSP remains untouched. Fractional opposite-side endpoints describe
    endpoint-linear geometry, including gapped hits, and are never sequence data.
    Both unset transforms retain the historical endpoint geometry exactly.
    """
    if all(transform.start_coordinate is None for transform in transforms):
        return (endpoints,)
    projected = []
    source_spans = []
    for (start, end), transform in zip((endpoints[:2], endpoints[2:]), transforms, strict=True):
        if any(isinstance(value, bool) or int(value) != value for value in (start, end)):
            raise ValidationError("Comparison endpoints must be integer bases.")
        start, end = int(start), int(end)
        if not 1 <= min(start, end) <= max(start, end) <= transform.length:
            raise ValidationError("Comparison endpoints must lie within their record.")
        strand = 1 if start <= end else -1
        span = SourceInterval(min(start, end) - 1, max(start, end), strand)
        if transform.start_coordinate is None:
            directed = (span.start, span.end) if strand == 1 else (span.end, span.start)
            projected.append(((0.0, 1.0, *directed),))
            continue
        fragments = transform.project_local_parts((span,))
        source_span = SourceInterval(
            min(part.source_start for part in fragments),
            max(part.source_end for part in fragments),
            fragments[0].orientation * transform.source_step,
        )
        source_spans.append((transform, source_span))
        intervals = []
        traversed = 0
        length = span.end - span.start
        for part in fragments:
            width = part.display_end - part.display_start
            directed = ((part.display_start, part.display_end) if part.orientation == 1
                        else (part.display_end, part.display_start))
            intervals.append((traversed / length, (traversed + width) / length, *directed))
            traversed += width
        projected.append(tuple(intervals))
    cuts = (0.0, *alignment_cut_breakpoints(*source_spans), 1.0)
    result = []
    for left, right in zip(cuts, cuts[1:]):
        coordinates = []
        for intervals in projected:
            lo, hi, start, end = next(part for part in intervals if part[0] <= (left + right) / 2 < part[1])
            coordinates.extend(start + (end - start) * (t - lo) / (hi - lo) for t in (left, right))
        result.append(tuple(coordinates))
    return tuple(result)


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
