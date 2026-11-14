"""What one diagram read from its inputs, and what it drew of them: the input summary.

A render that draws nothing from a file it read, or from a region annotation set it was
given, must say why. This module is the
one owner of that report for every surface: the CLI logs :func:`input_summary_lines`
as INFO and each :class:`InputWarning` as WARNING, and the Web shows the same text,
which the Worker carries as :func:`input_summary_payload` beside its Results.

The feature counts come from the feature-drawing decision
(:func:`gbdraw.features.visibility.should_render_feature`, with the inputs the
drawing reads) and from the GFF3 type filter of record planning; the comparison
counts come from the threshold tests of
:func:`gbdraw.io.comparisons.filter_comparison_dataframe`.
"""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal, TypedDict

from Bio.SeqRecord import SeqRecord

from gbdraw.auto_sizes import DrawnExtent
from gbdraw.config.models.canvas import LinearWrapLength
from gbdraw.core.record_metadata import _drawn_source_span
from gbdraw.features.ids import feature_instances
from gbdraw.features.overrides import feature_override_lookup
from gbdraw.features.placement import ResolvedRecordFeatureInputs
from gbdraw.features.selector_values import get_feature_type
from gbdraw.features.visibility import should_render_feature
from gbdraw.io.comparison_sequences import RepeatedRecordIds
from gbdraw.io.comparisons import ComparisonThresholdName, ComparisonThresholdTally
from gbdraw.layout.linear_wrap import LinearWrapSummary
from gbdraw.linear_comparison import ComparisonTableTally

from .record_planning import ResolvedRecordProvenance

logger = logging.getLogger(__name__)

InputWarningCode = Literal[
    "no_drawable_features",
    "comparison_rows_filtered_out",
    "comparison_table_empty",
    "annotation_set_not_drawn",
    "comparison_record_id_repeated",
]
NO_DRAWABLE_FEATURES: InputWarningCode = "no_drawable_features"
COMPARISON_ROWS_FILTERED_OUT: InputWarningCode = "comparison_rows_filtered_out"
COMPARISON_TABLE_EMPTY: InputWarningCode = "comparison_table_empty"
ANNOTATION_SET_NOT_DRAWN: InputWarningCode = "annotation_set_not_drawn"
COMPARISON_RECORD_ID_REPEATED: InputWarningCode = "comparison_record_id_repeated"
# Warnings of inputs drawn in full: the Web does not head them "Part of the input was not drawn".
INPUT_WARNINGS_OF_DRAWN_INPUT: frozenset[InputWarningCode] = frozenset({COMPARISON_RECORD_ID_REPEATED})

_THRESHOLD_WORDING: dict[ComparisonThresholdName, tuple[str, str, str]] = {
    # name: (relation of a left-out row, threshold name, unit)
    "evalue": ("above", "e-value", ""),
    "bitscore": ("below", "bitscore", ""),
    "identity": ("below", "identity", "%"),
    "alignment_length": ("below", "alignment length", ""),
}


@dataclass(frozen=True)
class FeatureTypeTally:
    """The features of one type in one drawn record: drawn, and left out by reason."""

    feature_type: str
    drawn: int
    # Its type is not among the feature types to draw.
    type_not_drawn: int
    # Its type is drawn, but a visibility rule or a feature edit hides it.
    hidden: int
    # Left out when the GFF3 file was read (types not drawn), counted over the
    # whole source record: the load filters before a region crops the record.
    not_read: int


@dataclass(frozen=True)
class RecordInputSummary:
    """One drawn record: where it was read, which part is drawn, and its feature tallies."""

    record_index: int
    record_id: str
    source_name: str
    # The files it was read from, as the request named them: with ``source_record_index`` and ``region``, the record's
    # identity, which a batch's log tells apart (two files of one base name in two folders are two inputs).
    source_paths: tuple[str, ...]
    source_record_index: int
    source_record_count: int
    drawn_length: int
    source_length: int
    # The drawn source bases, 1-based and inclusive, when the record is cropped; a crop across the origin starts after
    # it ends.
    region: tuple[int, int] | None
    features: tuple[FeatureTypeTally, ...]
    # The region of the source record the request asked for (1-based, inclusive, as asked: its end may run past the
    # record), ``None`` when it asked for none. The Web compares it with the record's card.
    requested_region: tuple[int, int] | None = None
    # False when the track stack has no features track (``slots_draw_features``): no feature is drawn, and
    # ``features`` is empty.
    features_track: bool = True
    # The shortest record one circle draws, when this record is shorter and so not drawn
    # (``gbdraw.layout.record_units.left_off_circle``); ``features`` is then empty.
    left_off_circle_under_bp: int | None = None

    @property
    def drawn_features(self) -> int:
        return sum(item.drawn for item in self.features)


@dataclass(frozen=True)
class InputWarning:
    """A reason the drawing shows less than its inputs hold, or names part of them ambiguously; the run still draws."""

    code: InputWarningCode
    message: str
    # The input it concerns, by identity: a record's place in its files, a comparison table's file and record pair,
    # a region annotation set's ID, a comparison sequence file's path. A batch logs one warning of a code per input.
    subject: tuple[str, ...]


# A warning of one input (``InputWarning.code``, ``InputWarning.subject``): what a batch logs once.
InputWarningKey = tuple[InputWarningCode, tuple[str, ...]]


@dataclass(frozen=True)
class InputSummary:
    """What one Result read and drew."""

    records: tuple[RecordInputSummary, ...]
    comparisons: tuple[ComparisonTableTally, ...] = ()
    # The region annotation sets with annotations that the explicit track stack draws nothing of
    # (``gbdraw.annotations.planning.annotation_sets_left_out``), in the order the sets were given.
    annotation_sets_left_out: tuple[str, ...] = ()
    # The comparison genome files the render read whose record IDs repeat (OV-2954).
    repeated_comparison_ids: tuple[RepeatedRecordIds, ...] = ()


class InputWarningPayload(TypedDict):
    code: InputWarningCode
    message: str


class LinearWrapPayload(TypedDict):
    """The rows a wrapped Linear Result drew (:class:`LinearWrapSummary`); the Web reads it as sent."""

    rowLength: int
    rowCount: int
    # What was asked for: ``"auto"`` or a length in bp, never ``"none"`` (nothing is cut then).
    requested: LinearWrapLength
    recordRows: list[int]
    # The bases on the longest drawn axis, which sizes and ticks follow (DrawnExtent.drawn_axis_length).
    drawnAxisLength: int


class RegionPayload(TypedDict):
    start: int
    end: int


class RecordPayload(TypedDict):
    """One drawn record (:class:`RecordInputSummary`): its source record's bases and the region of the source the
    request asked for. Its drawn bases follow from them by the crop rule (``region_slice``), which the Web
    export check applies in Python (D-70)."""

    sourceLength: int
    region: RegionPayload | None


class InputSummaryPayload(TypedDict):
    """The input summary of one Web Result (``metadata.inputSummaries``); not saved in a Session."""

    resultIndex: int
    resultName: str
    lines: list[str]
    warnings: list[InputWarningPayload]
    # Each drawn record, in drawing order: the Web tells from them whether the records its cards ask for now are the
    # ones this Result drew (and a wrapped Result's export check sends them back).
    records: list[RecordPayload]
    # ``None`` when the Result cut no record (a Circular Result, or a Linear one drawn on one row each).
    linearWrap: LinearWrapPayload | None


def _record_feature_tallies(
    record: SeqRecord,
    *,
    selected_features: Sequence[str],
    feature_visibility_rules: list[dict[str, object]] | None,
    specific_color_rules: Mapping[str, object] | None,
    record_features: ResolvedRecordFeatureInputs | None,
    types_left_out_at_load: Sequence[tuple[str, int]],
) -> tuple[FeatureTypeTally, ...]:
    selected = set(selected_features)
    drawn: Counter[str] = Counter()
    type_not_drawn: Counter[str] = Counter()
    not_read: Counter[str] = Counter(dict(types_left_out_at_load))
    hidden: Counter[str] = Counter()
    override_of = feature_override_lookup(record, record_features.overrides if record_features else None)
    instance_of = feature_instances(record)
    for feature in record.features:
        feature_type = get_feature_type(feature)
        if should_render_feature(
            feature,
            selected,
            feature_visibility_rules=feature_visibility_rules,
            record_id=record.id,
            specific_color_rules=specific_color_rules,
            feature_override=override_of(feature),
            instance=instance_of(feature),
        ):
            drawn[feature_type] += 1
        elif feature_type in selected:
            hidden[feature_type] += 1
        else:
            type_not_drawn[feature_type] += 1
    return tuple(
        FeatureTypeTally(name, drawn[name], type_not_drawn[name], hidden[name], not_read[name])
        for name in sorted({*drawn, *type_not_drawn, *hidden, *not_read}, key=lambda name: (name.lower(), name))
    )


def _drawn_region(record: SeqRecord, source_length: int) -> tuple[int, int] | None:
    return None if len(record.seq) >= source_length else _drawn_source_span(record)


def _region_text(region: tuple[int, int]) -> str:
    """A drawn region in the summary's words; a region whose start is after its end runs across the origin."""

    start, end = region
    return f"{start:,}..{end:,}{' across the origin' if start > end else ''}"


def summarize_inputs(
    records: Sequence[SeqRecord],
    provenance: Sequence[ResolvedRecordProvenance],
    *,
    selected_features: Sequence[str],
    feature_visibility_rules: list[dict[str, object]] | None,
    specific_color_rules: Mapping[str, object] | None,
    record_features: Sequence[ResolvedRecordFeatureInputs] = (),
    comparisons: Sequence[ComparisonTableTally] = (),
    draws_features: bool = True,
    annotation_sets_left_out: Sequence[str] = (),
    left_off_circle: Mapping[int, int] | None = None,
) -> InputSummary:
    """The input summary of the drawn ``records`` (aligned with their ``provenance``).

    ``record_features`` is empty or aligned with ``records``; the other feature
    arguments are the ones the drawing reads. ``draws_features`` is false when
    the track stack has no features track: no feature is counted as drawn.
    ``annotation_sets_left_out`` names the region annotation sets the track
    stack draws nothing of. ``left_off_circle`` maps the index of each record one circle leaves out to the
    shortest record it draws.
    """

    summaries = []
    left_off = left_off_circle or {}
    for index, (record, item) in enumerate(zip(records, provenance, strict=True)):
        source_length = item.source_length if item.source_length is not None else len(record.seq)
        under_bp = left_off.get(index)
        summaries.append(RecordInputSummary(
            record_index=index,
            record_id=str(record.id),
            source_name=item.source_name,
            source_paths=item.source_paths,
            source_record_index=item.source_record_index,
            source_record_count=item.source_record_count,
            drawn_length=len(record.seq),
            source_length=source_length,
            region=_drawn_region(record, source_length),
            requested_region=None if item.region is None else (item.region.start, item.region.end),
            features=_record_feature_tallies(
                record,
                selected_features=selected_features,
                feature_visibility_rules=feature_visibility_rules,
                specific_color_rules=specific_color_rules,
                record_features=record_features[index] if record_features else None,
                types_left_out_at_load=item.feature_types_left_out_at_load,
            ) if draws_features and under_bp is None else (),
            features_track=draws_features,
            left_off_circle_under_bp=under_bp,
        ))
    return InputSummary(
        records=tuple(summaries),
        comparisons=tuple(comparisons),
        annotation_sets_left_out=tuple(annotation_sets_left_out),
    )


def _count(value: int, noun: str, plural: str | None = None) -> str:
    return f"{value:,} {noun}" if value == 1 else f"{value:,} {plural or noun + 's'}"


def _type_counts(counts: Sequence[tuple[str, int]]) -> str:
    return ", ".join(f"{value:,} {name}" for name, value in counts)


def _types_not_drawn(record: RecordInputSummary) -> list[tuple[str, int]]:
    """The features of types not drawn, within the drawn part of the record.

    The GFF3 types left out at load count over the whole record, so a cropped
    record names them apart (``_types_not_read_in_record``).
    """

    whole = 0 if record.region is not None else 1
    counts = [(item.feature_type, item.type_not_drawn + whole * item.not_read) for item in record.features]
    return [(name, value) for name, value in counts if value]


def _types_not_read_in_record(record: RecordInputSummary) -> list[tuple[str, int]]:
    if record.region is None:
        return []
    return [(item.feature_type, item.not_read) for item in record.features if item.not_read]


def _record_subject(record: RecordInputSummary) -> tuple[str, ...]:
    region = "" if record.region is None else _region_text(record.region)
    return (str(record.source_record_index), region, *record.source_paths)


def _record_label(record: RecordInputSummary) -> str:
    return f"Record {record.record_id} ({record.source_name})"


def _record_line(record: RecordInputSummary) -> str:
    place = f"record {record.source_record_index + 1:,} of {record.source_record_count:,}"
    extent = (
        f"{_region_text(record.region)} of {record.source_length:,} bp"
        if record.region is not None
        else f"{record.drawn_length:,} bp"
    )
    if record.left_off_circle_under_bp is not None:
        return (
            f"{_record_label(record)}, {place}, {extent}: not drawn, shorter than "
            f"{record.left_off_circle_under_bp:,} bp, the shortest record drawn on one circle."
        )
    if not record.features_track:
        return f"{_record_label(record)}, {place}, {extent}: drew no features (no enabled features track slot)."
    drawn = [(item.feature_type, item.drawn) for item in record.features if item.drawn]
    line = f"{_record_label(record)}, {place}, {extent}: drew {_count(record.drawn_features, 'feature')}"
    if drawn:
        line += f" ({_type_counts(drawn)})"
    not_drawn = _types_not_drawn(record)
    hidden = [(item.feature_type, item.hidden) for item in record.features if item.hidden]
    if not_drawn:
        line += f"; left out {_type_counts(not_drawn)} (type{'s' if len(not_drawn) > 1 else ''} not drawn)"
    if not_read := _types_not_read_in_record(record):
        line += f"; the whole record has {_type_counts(not_read)} of type{'s' if len(not_read) > 1 else ''} not drawn"
    if hidden:
        line += f"; hid {_type_counts(hidden)} (visibility rule or feature edit)"
    return line + "."


def _comparison_label(table: ComparisonTableTally, records: Sequence[RecordInputSummary]) -> str:
    def endpoint(index: int) -> str:
        record_id = records[index].record_id if index < len(records) else ""
        return f"record #{index + 1} {record_id}".rstrip()

    pair = f"{endpoint(table.query_record_index)} and {endpoint(table.subject_record_index)}"
    if table.protein_matches:
        return f"Protein comparison of {pair}"
    if table.source_name:
        return f"Comparison table {table.source_name} ({pair})"
    return f"Comparison of {pair}"


def _threshold_text(item: ComparisonThresholdTally) -> str:
    relation, name, unit = _THRESHOLD_WORDING[item.name]
    return (
        f"{_count(item.rows_left_out, 'row')} {relation} {name} {item.threshold:g}{unit} "
        f"(best {item.best_value:g}{unit})"
    )


def _comparison_line(table: ComparisonTableTally, records: Sequence[RecordInputSummary]) -> str:
    tally = table.tally
    if table.protein_matches:
        # The rows the protein analysis left out are not counted (OV-2389).
        return f"{_comparison_label(table, records)}: drew {_count(tally.rows_kept, 'protein match', 'protein matches')}."
    line = f"{_comparison_label(table, records)}: kept {tally.rows_kept:,} of {_count(tally.rows_read, 'row')}"
    if tally.thresholds:
        line += "; " + "; ".join(_threshold_text(item) for item in tally.thresholds)
    return line + "."


def input_summary_lines(summary: InputSummary) -> tuple[str, ...]:
    """One line per drawn record and per comparison table (the CLI's INFO lines)."""

    return (
        *(_record_line(record) for record in summary.records),
        *(_comparison_line(table, summary.records) for table in summary.comparisons),
    )


def _no_feature_reason(record: RecordInputSummary) -> str:
    hidden = [(item.feature_type, item.hidden) for item in record.features if item.hidden]
    not_drawn = _types_not_drawn(record)
    if record.region is not None and not hidden and not not_drawn:
        return f"its region {_region_text(record.region)} holds no feature of a type drawn"
    if not record.features:
        return "it has no features"
    reasons = []
    if not_drawn:
        reasons.append(
            f"its {_type_counts(not_drawn)} features are of types not drawn; add a type to the "
            "feature types to draw (-k/--features) to draw it"
        )
    if hidden:
        reasons.append(f"its {_type_counts(hidden)} features are hidden by a visibility rule or a feature edit")
    return "; ".join(reasons)


def input_warnings(summary: InputSummary) -> tuple[InputWarning, ...]:
    """The reasons the drawing shows nothing of a record, of a comparison table, or of a region annotation set, and
    the comparison files whose hover and span export name some records ambiguously: one per input and code, so one
    file given twice warns once."""

    warnings = [
        InputWarning(
            NO_DRAWABLE_FEATURES,
            f"{_record_label(record)} has no feature to draw: {_no_feature_reason(record)}.",
            _record_subject(record),
        )
        for record in summary.records
        # A stack without a features track draws none by request, and one circle leaves a short record out (its
        # line says so); the inputs are not the reason.
        if record.features_track and record.left_off_circle_under_bp is None and record.drawn_features == 0
    ]
    for table in summary.comparisons:
        tally = table.tally
        label = _comparison_label(table, summary.records)
        # A Linear comparison: a table's file, or a record pair (a Linear run draws one Result).
        table_subject = (table.source_name or "", str(table.query_record_index), str(table.subject_record_index))
        if table.protein_matches and tally.rows_kept == 0:
            warnings.append(InputWarning(COMPARISON_TABLE_EMPTY, f"{label} has no protein match to draw.", table_subject))
        elif tally.rows_read == 0:
            warnings.append(
                InputWarning(COMPARISON_TABLE_EMPTY, f"{label} has no rows, so no match is drawn.", table_subject)
            )
        elif tally.rows_kept == 0:
            warnings.append(InputWarning(
                COMPARISON_ROWS_FILTERED_OUT,
                f"{label}: the thresholds leave out all {_count(tally.rows_read, 'row')}, so no match is drawn: "
                + "; ".join(_threshold_text(item) for item in tally.thresholds) + ".",
                table_subject,
            ))
    warnings.extend(
        InputWarning(
            ANNOTATION_SET_NOT_DRAWN,
            f"Region annotation set {set_id} is not drawn: no annotations track slot names it. "
            f"Add an annotations track slot for it (set_id={set_id}) to draw it.",
            (set_id,),
        )
        for set_id in summary.annotation_sets_left_out
    )
    warnings.extend(
        InputWarning(COMPARISON_RECORD_ID_REPEATED, item.message, (item.source_path,))
        for item in summary.repeated_comparison_ids
    )
    unique: dict[InputWarningKey, InputWarning] = {}
    for warning in warnings:
        unique.setdefault((warning.code, warning.subject), warning)
    return tuple(unique.values())


def log_input_warnings(summary: InputSummary, *, warned: set[InputWarningKey] | None = None) -> None:
    """Write the summary's warnings as WARNING lines (the Python API's ``draw_circular`` and ``draw_linear``).

    ``warned`` holds the warnings the earlier Results of a batch wrote, by code and
    input: a warning of an input they share (one comparison file every Result
    reads) is written once, and the new ones join it. Two inputs whose warnings
    read alike (two files of one base name) are each written.
    """

    for warning in input_warnings(summary):
        if warned is not None:
            key: InputWarningKey = (warning.code, warning.subject)
            if key in warned:
                continue
            warned.add(key)
        logger.warning("%s", warning.message)


def log_input_summary(summary: InputSummary, *, warned: set[InputWarningKey] | None = None) -> None:
    """Write the summary as CLI INFO lines and its warnings as WARNING lines, each of a batch once (``warned``)."""

    for line in input_summary_lines(summary):
        logger.info("%s", line)
    log_input_warnings(summary, warned=warned)


def with_source_names(summary: InputSummary, names: Mapping[str, str]) -> InputSummary:
    """The summary with each source file named by ``names`` (file name read -> file name shown).

    The Web reads each file from a staged copy named by its content; the page
    knows the name of the file it read.
    """

    return replace(
        summary,
        records=tuple(
            replace(record, source_name=names.get(record.source_name, record.source_name))
            for record in summary.records
        ),
        comparisons=tuple(
            replace(table, source_name=names.get(table.source_name, table.source_name))
            if table.source_name is not None else table
            for table in summary.comparisons
        ),
        repeated_comparison_ids=tuple(
            replace(item, source_name=names.get(item.source_name, item.source_name))
            for item in summary.repeated_comparison_ids
        ),
    )


def _linear_wrap_payload(wrap: LinearWrapSummary) -> LinearWrapPayload:
    return {
        "rowLength": wrap.row_length,
        "rowCount": wrap.row_count,
        "requested": wrap.requested,
        "recordRows": list(wrap.record_rows),
        "drawnAxisLength": DrawnExtent("linear", wrap.record_lengths, wrap.row_length).drawn_axis_length,
    }


def _record_payload(record: RecordInputSummary) -> RecordPayload:
    region = record.requested_region
    return {
        "sourceLength": record.source_length,
        "region": None if region is None else {"start": region[0], "end": region[1]},
    }


def input_summary_payload(
    summary: InputSummary,
    *,
    result_index: int,
    result_name: str,
    linear_wrap: LinearWrapSummary | None = None,
) -> InputSummaryPayload:
    """The Web's copy of the summary of one Result: the same lines and warnings the CLI logs, and the rows a wrapped
    Linear Result drew (``linear_wrap``, the CLI logs its message)."""

    return {
        "resultIndex": result_index,
        "resultName": result_name,
        "lines": list(input_summary_lines(summary)),
        "warnings": [{"code": item.code, "message": item.message} for item in input_warnings(summary)],
        "records": [_record_payload(record) for record in summary.records],
        "linearWrap": None if linear_wrap is None else _linear_wrap_payload(linear_wrap),
    }


__all__ = [
    "ANNOTATION_SET_NOT_DRAWN",
    "COMPARISON_RECORD_ID_REPEATED",
    "COMPARISON_ROWS_FILTERED_OUT",
    "COMPARISON_TABLE_EMPTY",
    "FeatureTypeTally",
    "INPUT_WARNINGS_OF_DRAWN_INPUT",
    "InputSummary",
    "InputSummaryPayload",
    "InputWarning",
    "InputWarningCode",
    "InputWarningKey",
    "InputWarningPayload",
    "LinearWrapPayload",
    "NO_DRAWABLE_FEATURES",
    "RecordInputSummary",
    "RecordPayload",
    "RegionPayload",
    "input_summary_lines",
    "input_summary_payload",
    "input_warnings",
    "log_input_summary",
    "log_input_warnings",
    "summarize_inputs",
    "with_source_names",
]
