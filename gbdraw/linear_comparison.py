"""Normalized comparison contracts for Linear diagrams."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Literal, Sequence, cast

import pandas as pd
from pandas import DataFrame

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.exceptions import ComparisonIdentityError, ValidationError
from gbdraw.layout.record_coordinates import (
    DisplayFragment,
    RecordDisplayTransform,
    SourceInterval,
    alignment_cut_breakpoints,
)

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

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
    projected: list[tuple[tuple[float, float, int, int], ...]] = []
    source_spans = []
    for (start, end), transform in zip((endpoints[:2], endpoints[2:]), transforms, strict=True):
        if any(isinstance(value, bool) or int(value) != value for value in (start, end)):
            raise ValidationError("Comparison endpoints must be integer bases.")
        start, end = int(start), int(end)
        if not 1 <= min(start, end) <= max(start, end) <= transform.length:
            raise ValidationError("Comparison endpoints must lie within their record.")
        strand: Literal[-1, 1] = 1 if start <= end else -1
        span = SourceInterval(min(start, end) - 1, max(start, end), strand)
        if transform.start_coordinate is None:
            directed = (span.start, span.end) if strand == 1 else (span.end, span.start)
            projected.append(((0.0, 1.0, *directed),))
            continue
        # A transform with a start coordinate always returns display fragments.
        fragments = cast("tuple[DisplayFragment, ...]", transform.project_local_parts((span,)))
        source_span = SourceInterval(
            min(part.source_start for part in fragments),
            max(part.source_end for part in fragments),
            cast("Literal[-1, 1]", fragments[0].orientation * transform.source_step),
        )
        source_spans.append((transform, source_span))
        intervals: list[tuple[float, float, int, int]] = []
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
    result: list[tuple[float, float, float, float]] = []
    for left, right in zip(cuts, cuts[1:]):
        coordinates: list[float] = []
        for endpoint_intervals in projected:
            lo, hi, start, end = next(
                part for part in endpoint_intervals if part[0] <= (left + right) / 2 < part[1]
            )
            coordinates.extend(start + (end - start) * (t - lo) / (hi - lo) for t in (left, right))
        result.append(cast("tuple[float, float, float, float]", tuple(coordinates)))
    return tuple(result)


@dataclass(frozen=True)
class LinearComparison:
    """A comparison result with explicit input-record endpoints.

    ``search_frame_text`` is the raw BLAST outfmt 6 text in the search frame
    when the planner produced ``matches`` from a LOSAT search or the Session
    decoder read them from a ``nucleotideBlast`` resource; the Session persists
    it as a ``nucleotideBlast`` resource, as the Web does.
    """

    query_record_index: int
    subject_record_index: int
    matches: DataFrame
    search_frame_text: str | None = field(default=None, repr=False, compare=False)

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


def _rows_without_feature_binding(frame: DataFrame) -> DataFrame:
    """Return the rows that carry no source-feature binding.

    Sources for one endpoint pair are merged before this check, so uploaded
    rows can share a frame with bound LOSATP or saved protein rows.
    """

    if not _FEATURE_BINDING_COLUMNS <= set(frame.columns):
        return frame
    binding = frame.loc[:, sorted(_FEATURE_BINDING_COLUMNS)]
    present = binding.notna() & binding.astype(str).apply(lambda column: column.str.strip() != "")
    return frame.loc[~present.any(axis=1)]


def project_search_frame_comparisons(
    comparisons: Sequence[LinearComparison],
    records: Sequence[SeqRecord],
) -> tuple[LinearComparison, ...]:
    """Project comparison table rows from the search frame into each record (PD-OI-073).

    Every comparison table (``-b`` files, comparison tables, uploaded and
    generated LOSAT rows) is read in the search frame: the selected and cropped
    record, 1-based, in the source strand. A reverse-complemented record maps a
    coordinate x to L + 1 - x. Rows bound to source features keep the view
    projection of the planner. A row outside 1..L of its record is rejected
    instead of drawn beyond the record. This is the one owner of the projection;
    persisted requests and Sessions keep the search frame.
    """

    projected: list[LinearComparison] = []
    for comparison in comparisons:
        frame = _rows_without_feature_binding(comparison.matches)
        if frame.empty:
            projected.append(comparison)
            continue
        for role, record_index, columns in (
            ("query", comparison.query_record_index, ("qstart", "qend")),
            ("subject", comparison.subject_record_index, ("sstart", "send")),
        ):
            record = records[record_index]
            length = len(record)
            values = frame.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce")
            outside = ~((values >= 1) & (values <= length) & (values == values.round())).all(axis=1)
            if outside.any():
                first = frame.loc[outside].iloc[0]
                raise ValidationError(
                    f"Comparison between query record #{comparison.query_record_index + 1} and subject "
                    f"record #{comparison.subject_record_index + 1}: {int(outside.sum())} row(s) have "
                    f"{role} coordinates outside 1..{length} of record "
                    f"{str(getattr(record, 'id', '')).strip()!r} (for example "
                    f"{first[columns[0]]}..{first[columns[1]]}). Comparison tables use coordinates "
                    "of the selected and cropped record in the source strand.",
                    diagnostic={
                        "code": "COMPARISON_INPUT",
                        "reason": "SEARCH_FRAME",
                        "column": 7 if role == "query" else 9,  # qstart or sstart
                    },
                )
        projected.append(reverse_unbound_endpoint_rows(
            comparison,
            *(_endpoint_frame(records[index]) for index in (comparison.query_record_index, comparison.subject_record_index)),
        ))
    return tuple(projected)


def _endpoint_frame(record: SeqRecord) -> tuple[int, bool]:
    return len(record), _read_coord_map(record)[1] == -1


def reverse_unbound_endpoint_rows(
    comparison: LinearComparison,
    query: tuple[int, bool],
    subject: tuple[int, bool],
) -> LinearComparison:
    """Map the unbound rows of each reversed endpoint ``(L, True)`` x -> L + 1 - x.

    The map is its own inverse: it projects search-frame rows onto a reversed
    record and converts rows of a reversed record back to its search frame.
    """

    frame = _rows_without_feature_binding(comparison.matches)
    if frame.empty or not (query[1] or subject[1]):
        return comparison
    updated = comparison.matches.copy()
    for (length, reverse), columns in ((query, ("qstart", "qend")), (subject, ("sstart", "send"))):
        if reverse:
            for column in columns:
                updated.loc[frame.index, column] = int(length) + 1 - pd.to_numeric(frame[column]).astype(int)
    return replace(comparison, matches=updated)


def reverse_endpoint_table_text(
    text: str,
    query: tuple[int, bool],
    subject: tuple[int, bool],
) -> str:
    """Rewrite one nucleotide table between a record's frame and its reverse complement.

    ``query`` and ``subject`` are ``(L, reversed)`` of each endpoint record;
    rows of a reversed endpoint map x -> L + 1 - x, the map of
    :func:`project_search_frame_comparisons`, which is its own inverse. The Web
    Session reader converts origin/main Sessions (version 42 and older, rows
    stored after the reverse complement) once at Load; the Session writer
    converts rows into the frame of a reverse-complemented record it persists
    as a sequence.
    """

    from io import StringIO

    from gbdraw.io.comparisons import read_comparison_table

    frame = read_comparison_table(StringIO(str(text or "")), label="comparison table")
    for (length, reverse), columns in ((query, ("qstart", "qend")), (subject, ("sstart", "send"))):
        values = frame.loc[:, list(columns)].astype(int)
        if reverse and not ((values >= 1) & (values <= int(length))).all(axis=None):
            raise ValidationError(
                f"Comparison rows lie outside 1..{int(length)} of their reverse-complemented record.",
                diagnostic={"code": "COMPARISON_INPUT", "reason": "SEARCH_FRAME"},
            )
    converted = reverse_unbound_endpoint_rows(LinearComparison(0, 1, frame), query, subject).matches
    return converted.to_csv(sep="\t", header=False, index=False, lineterminator="\n")


COMPARISON_RECORD_ID_UNMATCHED = "comparison_record_id_unmatched"


@dataclass(frozen=True)
class ComparisonRecordIdWarning:
    """Table rows whose sequence IDs name no displayed record (PD-OI-074).

    The rows stay on the endpoint records assigned by position. The CLI logs
    ``message``; the Web adapter carries it in the Run metadata because browser
    execution discards logging.
    """

    query_record_index: int
    subject_record_index: int
    query_record_id: str
    subject_record_id: str
    row_count: int
    example_ids: tuple[str, ...]
    message: str
    code: str = COMPARISON_RECORD_ID_UNMATCHED


def validate_linear_comparison_record_ids(
    comparisons: Sequence[LinearComparison],
    records: Sequence[object],
) -> tuple[ComparisonRecordIdWarning, ...]:
    """Check table query/subject IDs against each comparison's endpoint records.

    A row that names the opposite endpoint or another displayed record raises
    ComparisonIdentityError. IDs that name no displayed record keep positional
    placement; each affected endpoint pair is logged and returned as a warning.
    Version suffixes (``.1``) are tolerated. Rows bound to source features are
    checked by feature identity instead. Every Linear comparison input (``-b``
    files, comparison tables, typed and Web requests) reaches this check beside
    the topology check.
    """

    warnings: list[ComparisonRecordIdWarning] = []
    id_index = _record_id_index(records)
    for comparison in comparisons:
        frame = _rows_without_feature_binding(comparison.matches)
        if frame.empty:
            continue
        endpoints = {
            "query": comparison.query_record_index,
            "subject": comparison.subject_record_index,
        }
        record_ids = {
            role: str(getattr(records[endpoint], "id", "")).strip()
            for role, endpoint in endpoints.items()
        }
        pair = (
            f"query record #{endpoints['query'] + 1} {record_ids['query']!r} and subject record "
            f"#{endpoints['subject'] + 1} {record_ids['subject']!r}"
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
                        diagnostic={
                            "code": "COMPARISON_INPUT",
                            "reason": "RECORD_ID",
                            "column": 1 if role == "query" else 2,  # BLAST outfmt 6 column
                        },
                    )
            unknown_rows |= table_ids.isin(unknown_role_ids)
            unknown_ids.extend(value for value in unknown_role_ids if value not in unknown_ids)
        if unknown_ids:
            row_count = int(unknown_rows.sum())
            example_ids = tuple(unknown_ids[:3])
            message = (
                f"Comparison between {pair}: {row_count} row(s) use sequence IDs that match no "
                f"displayed record ({', '.join(repr(value) for value in example_ids)}); "
                "these rows are drawn on the records assigned by position."
            )
            logger.warning("WARNING: %s", message)
            warnings.append(
                ComparisonRecordIdWarning(
                    query_record_index=endpoints["query"],
                    subject_record_index=endpoints["subject"],
                    query_record_id=record_ids["query"],
                    subject_record_id=record_ids["subject"],
                    row_count=row_count,
                    example_ids=example_ids,
                    message=message,
                )
            )
    return tuple(warnings)


__all__ = [
    "COMPARISON_RECORD_ID_UNMATCHED",
    "ComparisonRecordIdWarning",
    "LinearComparison",
    "merge_linear_comparisons",
    "project_search_frame_comparisons",
    "validate_linear_comparison_record_ids",
    "validate_linear_comparison_topology",
]
