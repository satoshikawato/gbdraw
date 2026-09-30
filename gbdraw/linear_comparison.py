"""Normalized comparison contracts for Linear diagrams."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Literal, Sequence

import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ComparisonIdentityError, ValidationError
from gbdraw.layout.record_coordinates import RecordDisplayTransform, SourceInterval, alignment_cut_breakpoints

logger = logging.getLogger(__name__)


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


_VERSION_SUFFIX = re.compile(r"\.[0-9]+$")
_FEATURE_BINDING_COLUMNS = frozenset(
    f"{role}_{name}" for role in ("query", "subject") for name in ("feature_index", "feature_svg_id")
)


def _record_id_index(records: Sequence[object]) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    """Index record IDs and names exactly and without a version suffix."""

    exact: dict[str, list[int]] = {}
    loose: dict[str, list[int]] = {}
    for index, record in enumerate(records):
        for attribute in ("id", "name"):
            text = str(getattr(record, attribute, "") or "").strip()
            if not text or text == "<unknown name>":
                continue
            for keys, key in ((exact, text), (loose, _VERSION_SUFFIX.sub("", text))):
                if index not in keys.setdefault(key, []):
                    keys[key].append(index)
    return exact, loose


def _resolve_table_id(
    value: str,
    endpoint: int,
    index: tuple[dict[str, list[int]], dict[str, list[int]]],
) -> tuple[Literal["match", "conflict", "unknown"], int | None]:
    """Match exact IDs before version-tolerant IDs so distinct versions stay distinct."""

    exact, loose = index
    for named in (exact.get(value, []), loose.get(_VERSION_SUFFIX.sub("", value), [])):
        if endpoint in named:
            return "match", endpoint
        if named:
            return "conflict", named[0]
    return "unknown", None


def validate_linear_comparison_record_ids(
    comparisons: Sequence[LinearComparison],
    records: Sequence[object],
) -> None:
    """Check table query/subject IDs against each comparison's endpoint records.

    A row that names the opposite endpoint or another displayed record raises
    ComparisonIdentityError. IDs that name no displayed record keep positional
    placement with a warning. Version suffixes (``.1``) are tolerated. Tables
    whose rows are bound to source features are checked by feature identity
    instead. Every Linear comparison input (``-b`` files, comparison tables,
    typed and Web requests) reaches this check beside the topology check.
    """

    id_index = _record_id_index(records)
    for comparison in comparisons:
        frame = comparison.matches
        if frame.empty or _FEATURE_BINDING_COLUMNS <= set(frame.columns):
            continue
        endpoints = {
            "query": comparison.query_record_index,
            "subject": comparison.subject_record_index,
        }
        pair = (
            f"query record #{endpoints['query'] + 1} "
            f"{str(getattr(records[endpoints['query']], 'id', '')).strip()!r} and subject record "
            f"#{endpoints['subject'] + 1} {str(getattr(records[endpoints['subject']], 'id', '')).strip()!r}"
        )
        unknown_rows = pd.Series(False, index=frame.index)
        unknown_ids: list[str] = []
        for role, endpoint in endpoints.items():
            if role not in frame.columns:
                continue
            table_ids = frame[role].astype(str).str.strip()
            unknown_role_ids: list[str] = []
            for value in table_ids.unique():
                status, named = _resolve_table_id(str(value), endpoint, id_index)
                if status == "unknown":
                    unknown_role_ids.append(str(value))
                elif status == "conflict" and named is not None:
                    other = endpoints["subject" if role == "query" else "query"]
                    target = (
                        f"the {'subject' if role == 'query' else 'query'} record"
                        if named == other
                        else f"record #{named + 1} {str(getattr(records[named], 'id', '')).strip()!r}"
                    )
                    raise ComparisonIdentityError(
                        f"Comparison between {pair}: the {role} column names {str(value)!r}, "
                        f"which is {target}, not the {role} record. Swap the query and subject "
                        "columns of the table or assign the table to the matching record pair.",
                        reason="RECORD_ID",
                    )
            unknown_rows |= table_ids.isin(unknown_role_ids)
            unknown_ids.extend(value for value in unknown_role_ids if value not in unknown_ids)
        if unknown_ids:
            logger.warning(
                "WARNING: Comparison between %s: %d row(s) use sequence IDs that match no "
                "displayed record (%s); these rows are drawn on the records assigned by position.",
                pair,
                int(unknown_rows.sum()),
                ", ".join(repr(value) for value in unknown_ids[:3]),
            )


__all__ = [
    "LinearComparison",
    "merge_linear_comparisons",
    "validate_linear_comparison_record_ids",
    "validate_linear_comparison_topology",
]
