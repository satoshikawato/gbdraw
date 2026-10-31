"""What one diagram read from its inputs, and what it drew of them: the input summary.

A render that draws nothing from a file it read must say why. This module is the
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

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.features.ids import feature_instances
from gbdraw.features.overrides import feature_override_lookup
from gbdraw.features.placement import ResolvedRecordFeatureInputs
from gbdraw.features.selector_values import get_feature_type
from gbdraw.features.visibility import should_render_feature
from gbdraw.io.comparisons import ComparisonThresholdName, ComparisonThresholdTally
from gbdraw.linear_comparison import ComparisonTableTally

from .record_planning import ResolvedRecordProvenance

logger = logging.getLogger(__name__)

InputWarningCode = Literal[
    "no_drawable_features",
    "comparison_rows_filtered_out",
    "comparison_table_empty",
]
NO_DRAWABLE_FEATURES: InputWarningCode = "no_drawable_features"
COMPARISON_ROWS_FILTERED_OUT: InputWarningCode = "comparison_rows_filtered_out"
COMPARISON_TABLE_EMPTY: InputWarningCode = "comparison_table_empty"

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


@dataclass(frozen=True)
class RecordInputSummary:
    """One drawn record: where it was read, which part is drawn, and its feature tallies."""

    record_index: int
    record_id: str
    source_name: str
    source_record_index: int
    source_record_count: int
    drawn_length: int
    source_length: int
    # The drawn source bases, 1-based and inclusive, when the record is cropped.
    region: tuple[int, int] | None
    features: tuple[FeatureTypeTally, ...]

    @property
    def drawn_features(self) -> int:
        return sum(item.drawn for item in self.features)


@dataclass(frozen=True)
class InputWarning:
    """A reason the drawing shows less than its inputs hold; the run still draws."""

    code: InputWarningCode
    message: str


@dataclass(frozen=True)
class InputSummary:
    """What one Result read and drew."""

    records: tuple[RecordInputSummary, ...]
    comparisons: tuple[ComparisonTableTally, ...] = ()


class InputWarningPayload(TypedDict):
    code: InputWarningCode
    message: str


class InputSummaryPayload(TypedDict):
    """The input summary of one Web Result (``metadata.inputSummaries``); not saved in a Session."""

    resultIndex: int
    resultName: str
    lines: list[str]
    warnings: list[InputWarningPayload]


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
    type_not_drawn: Counter[str] = Counter(dict(types_left_out_at_load))
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
        FeatureTypeTally(name, drawn[name], type_not_drawn[name], hidden[name])
        for name in sorted({*drawn, *type_not_drawn, *hidden}, key=lambda name: (name.lower(), name))
    )


def _drawn_region(record: SeqRecord, source_length: int) -> tuple[int, int] | None:
    length = len(record.seq)
    if length >= source_length:
        return None
    base, step = _read_coord_map(record)
    end = base + step * (length - 1)
    return (min(base, end), max(base, end))


def summarize_inputs(
    records: Sequence[SeqRecord],
    provenance: Sequence[ResolvedRecordProvenance],
    *,
    selected_features: Sequence[str],
    feature_visibility_rules: list[dict[str, object]] | None,
    specific_color_rules: Mapping[str, object] | None,
    record_features: Sequence[ResolvedRecordFeatureInputs] = (),
    comparisons: Sequence[ComparisonTableTally] = (),
) -> InputSummary:
    """The input summary of the drawn ``records`` (aligned with their ``provenance``).

    ``record_features`` is empty or aligned with ``records``; the other feature
    arguments are the ones the drawing reads.
    """

    summaries = []
    for index, (record, item) in enumerate(zip(records, provenance, strict=True)):
        source_length = item.source_length if item.source_length is not None else len(record.seq)
        summaries.append(RecordInputSummary(
            record_index=index,
            record_id=str(record.id),
            source_name=item.source_name,
            source_record_index=item.source_record_index,
            source_record_count=item.source_record_count,
            drawn_length=len(record.seq),
            source_length=source_length,
            region=_drawn_region(record, source_length),
            features=_record_feature_tallies(
                record,
                selected_features=selected_features,
                feature_visibility_rules=feature_visibility_rules,
                specific_color_rules=specific_color_rules,
                record_features=record_features[index] if record_features else None,
                types_left_out_at_load=item.feature_types_left_out_at_load,
            ),
        ))
    return InputSummary(records=tuple(summaries), comparisons=tuple(comparisons))


def _count(value: int, noun: str) -> str:
    return f"{value:,} {noun}" if value == 1 else f"{value:,} {noun}s"


def _type_counts(counts: Sequence[tuple[str, int]]) -> str:
    return ", ".join(f"{value:,} {name}" for name, value in counts)


def _record_label(record: RecordInputSummary) -> str:
    return f"Record {record.record_id} ({record.source_name})"


def _record_line(record: RecordInputSummary) -> str:
    place = f"record {record.source_record_index + 1:,} of {record.source_record_count:,}"
    extent = (
        f"{record.region[0]:,}..{record.region[1]:,} of {record.source_length:,} bp"
        if record.region is not None
        else f"{record.drawn_length:,} bp"
    )
    drawn = [(item.feature_type, item.drawn) for item in record.features if item.drawn]
    line = f"{_record_label(record)}, {place}, {extent}: drew {_count(record.drawn_features, 'feature')}"
    if drawn:
        line += f" ({_type_counts(drawn)})"
    not_drawn = [(item.feature_type, item.type_not_drawn) for item in record.features if item.type_not_drawn]
    hidden = [(item.feature_type, item.hidden) for item in record.features if item.hidden]
    if not_drawn:
        line += f"; left out {_type_counts(not_drawn)} (type{'s' if len(not_drawn) > 1 else ''} not drawn)"
    if hidden:
        line += f"; hid {_type_counts(hidden)} (visibility rule or feature edit)"
    return line + "."


def _comparison_label(table: ComparisonTableTally, records: Sequence[RecordInputSummary]) -> str:
    def endpoint(index: int) -> str:
        record_id = records[index].record_id if index < len(records) else ""
        return f"record #{index + 1} {record_id}".rstrip()

    pair = f"{endpoint(table.query_record_index)} and {endpoint(table.subject_record_index)}"
    if table.protein_hits:
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
    if table.protein_hits:
        # The hits the thresholds leave out are not counted (OV-2389).
        return f"{_comparison_label(table, records)}: drew {_count(tally.rows_kept, 'protein hit')} passing the thresholds."
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
    if record.region is not None and not hidden:
        # The types left out when a GFF3 file is read are counted over the whole record.
        return f"its region {record.region[0]:,}..{record.region[1]:,} holds no feature of a type drawn"
    if not record.features:
        return "it has no features"
    reasons = []
    not_drawn = [(item.feature_type, item.type_not_drawn) for item in record.features if item.type_not_drawn]
    if not_drawn:
        reasons.append(
            f"its {_type_counts(not_drawn)} features are of types not drawn; add a type to the "
            "feature types to draw (-k/--features) to draw it"
        )
    if hidden:
        reasons.append(f"its {_type_counts(hidden)} features are hidden by a visibility rule or a feature edit")
    return "; ".join(reasons)


def input_warnings(summary: InputSummary) -> tuple[InputWarning, ...]:
    """The reasons the drawing shows nothing of a record or of a comparison table."""

    warnings = [
        InputWarning(
            NO_DRAWABLE_FEATURES,
            f"{_record_label(record)} has no feature to draw: {_no_feature_reason(record)}.",
        )
        for record in summary.records
        if record.drawn_features == 0
    ]
    for table in summary.comparisons:
        tally = table.tally
        label = _comparison_label(table, summary.records)
        if table.protein_hits and tally.rows_kept == 0:
            warnings.append(InputWarning(
                COMPARISON_TABLE_EMPTY, f"{label} has no protein hit that passes the thresholds, so no match is drawn."
            ))
        elif tally.rows_read == 0:
            warnings.append(InputWarning(COMPARISON_TABLE_EMPTY, f"{label} has no rows, so no match is drawn."))
        elif tally.rows_kept == 0:
            warnings.append(InputWarning(
                COMPARISON_ROWS_FILTERED_OUT,
                f"{label}: the thresholds leave out all {_count(tally.rows_read, 'row')}, so no match is drawn: "
                + "; ".join(_threshold_text(item) for item in tally.thresholds) + ".",
            ))
    return tuple(warnings)


def log_input_summary(summary: InputSummary) -> None:
    """Write the summary as CLI INFO lines and its warnings as WARNING lines."""

    for line in input_summary_lines(summary):
        logger.info("INFO: %s", line)
    for warning in input_warnings(summary):
        logger.warning("WARNING: %s", warning.message)


def with_source_names(summary: InputSummary, names: Mapping[str, str]) -> InputSummary:
    """The summary with each source file named by ``names`` (file name read -> file name shown).

    The Web reads each file from a staged copy named by its content; the page
    knows the name of the file it read.
    """

    return InputSummary(
        records=tuple(
            replace(record, source_name=names.get(record.source_name, record.source_name))
            for record in summary.records
        ),
        comparisons=tuple(
            replace(table, source_name=names.get(table.source_name, table.source_name))
            if table.source_name is not None else table
            for table in summary.comparisons
        ),
    )


def input_summary_payload(summary: InputSummary, *, result_index: int, result_name: str) -> InputSummaryPayload:
    """The Web's copy of the summary of one Result: the same lines and warnings the CLI logs."""

    return {
        "resultIndex": result_index,
        "resultName": result_name,
        "lines": list(input_summary_lines(summary)),
        "warnings": [{"code": item.code, "message": item.message} for item in input_warnings(summary)],
    }


__all__ = [
    "COMPARISON_ROWS_FILTERED_OUT",
    "COMPARISON_TABLE_EMPTY",
    "FeatureTypeTally",
    "InputSummary",
    "InputSummaryPayload",
    "InputWarning",
    "InputWarningCode",
    "InputWarningPayload",
    "NO_DRAWABLE_FEATURES",
    "RecordInputSummary",
    "input_summary_lines",
    "input_summary_payload",
    "input_warnings",
    "log_input_summary",
    "summarize_inputs",
    "with_source_names",
]
