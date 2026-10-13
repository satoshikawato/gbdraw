"""The label fit report as the Web reads it, and the rows it computes later.

Python owns every value and every word of the report: the Web dialog shows
the rows' titles, the settings they change, the records' counts, and the
canvas sizes as given, and nothing in JavaScript places labels or searches
values (LP-D43, LP-D48). Rows a bounded search did not reach arrive
``not_computed`` (OD-F14 for a drawn figure, the font scan bound when nothing
was drawn); :class:`LabelFitRowSearches` keeps the last such figure or failure
of each mode so the Web can compute one row through :func:`label_fit_row`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict, cast, get_args

from svgwrite import Drawing

from gbdraw.api.diagram import label_fit_report, label_fit_row
from gbdraw.exceptions import LabelsDoNotFitError, ValidationError
from gbdraw.labels.fit import (
    ALL_OUTSIDE_SCOPE,
    FONT_SIZE_PATHS,
    PLACEMENT_PATH,
    RADIAL_PLACEMENT,
    ROW_TITLES,
    SCOPE_PATH,
    SPACING_PATH,
    X_OFFSET_PATH,
    Y_OFFSET_PATH,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelFitRowStatus,
    LabelFitValue,
    LabelKind,
    label_fit_record_counts,
    label_kind_words,
)

LABEL_FIT_PAYLOAD_SCHEMA: Literal[1] = 1

# The Web settings a row changes: the name and unit of each control, by the
# config path the request mapping projects it to. Paths that one control sets
# together (the font size classes) share a name and show once.
_SETTING_NAMES: Mapping[str, tuple[str, str]] = {
    X_OFFSET_PATH: ("Outer X Offset", ""),
    Y_OFFSET_PATH: ("Outer Y Offset", ""),
    **{path: ("Label Font Size", "pt") for path in FONT_SIZE_PATHS},
    SPACING_PATH: ("Circular Label Spacing", "px"),
    SCOPE_PATH: ("Label Mode", ""),
    PLACEMENT_PATH: ("Circular Label Placement", ""),
}
# A word value as the Web control names its option.
_VALUE_WORDS: Mapping[tuple[str, str], str] = {
    (SCOPE_PATH, ALL_OUTSIDE_SCOPE): "Out",
    (PLACEMENT_PATH, RADIAL_PLACEMENT): "Radial",
}
_ROW_KINDS: frozenset[str] = frozenset(get_args(LabelFitRowKind))
_STALE_REASON = "Generate again to compute this choice."
_UNKNOWN_CONFIG_PATH = {"code": "INPUT_INVALID", "field": "configOverrides", "reason": "UNKNOWN_CONFIG_PATH"}
_UNKNOWN_ROW_KIND = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


class LabelFitSettingPayload(TypedDict):
    name: str
    unit: str
    value: float | str


class LabelFitCanvasSizePayload(TypedDict):
    width: float
    height: float


class LabelFitRowPayload(TypedDict):
    kind: LabelFitRowKind
    title: str
    status: LabelFitRowStatus
    values: dict[str, float | str] | None
    settings: list[LabelFitSettingPayload]
    canvasSize: LabelFitCanvasSizePayload | None
    reason: str | None


class LabelFitRecordPayload(TypedDict):
    recordId: str
    labelKind: LabelKind
    labelKindTitle: str
    labelCount: int
    fontSize: float
    overlapPairs: int
    crossingPairs: int
    definitionOverlaps: int
    refused: bool
    counts: str


class LabelFitReportPayload(TypedDict):
    schema: Literal[1]
    drawn: bool
    overlapPairs: int
    crossingPairs: int
    definitionOverlaps: int
    records: list[LabelFitRecordPayload]
    rows: list[LabelFitRowPayload]


class LabelFitFailurePayload(TypedDict):
    """The report of a render that drew nothing, and the token of its rows to compute."""

    report: LabelFitReportPayload
    rowToken: str | None


def _setting(path: str, value: LabelFitValue) -> LabelFitSettingPayload:
    if path not in _SETTING_NAMES or (isinstance(value, str) and (path, value) not in _VALUE_WORDS):
        raise ValidationError(
            f"The label fit row sets an unknown config path or value: {path}.", diagnostic=_UNKNOWN_CONFIG_PATH
        )
    name, unit = _SETTING_NAMES[path]
    shown = _VALUE_WORDS[(path, value)] if isinstance(value, str) else float(value)
    return {"name": name, "unit": unit, "value": shown}


def _settings(values: Mapping[str, LabelFitValue]) -> list[LabelFitSettingPayload]:
    settings: dict[str, LabelFitSettingPayload] = {}
    for path, value in values.items():
        setting = _setting(path, value)
        settings.setdefault(setting["name"], setting)
    return list(settings.values())


def label_fit_row_payload(row: LabelFitRow) -> LabelFitRowPayload:
    """One row with its title, the settings it changes, and its canvas size."""

    values = dict(row.values) if row.values is not None else None
    return {
        "kind": row.kind,
        "title": ROW_TITLES[row.kind],
        "status": row.status,
        "values": (
            {path: value if isinstance(value, str) else float(value) for path, value in values.items()}
            if values
            else None
        ),
        "settings": _settings(values) if values else [],
        "canvasSize": (
            {"width": float(row.canvas_size[0]), "height": float(row.canvas_size[1])}
            if row.canvas_size is not None
            else None
        ),
        "reason": row.reason,
    }


def _record_payload(record: LabelFitRecord) -> LabelFitRecordPayload:
    return {
        "recordId": record.record_id,
        "labelKind": record.label_kind,
        "labelKindTitle": f"{label_kind_words(record)} labels",
        "labelCount": record.label_count,
        "fontSize": float(record.font_size),
        "overlapPairs": record.fit.overlap_pairs,
        "crossingPairs": record.fit.crossing_pairs,
        "definitionOverlaps": record.fit.definition_overlaps,
        "refused": record.refused,
        # A refused radial record has no counts; the reason stands in (OV-699).
        "counts": label_fit_record_counts(record),
    }


def label_fit_report_payload(report: LabelFitReport) -> LabelFitReportPayload:
    """The bounded JSON form of a report: identifiers, counts, values, and Python's words."""

    return {
        "schema": LABEL_FIT_PAYLOAD_SCHEMA,
        "drawn": report.drawn,
        "overlapPairs": report.overlap_pairs,
        "crossingPairs": report.crossing_pairs,
        "definitionOverlaps": report.definition_overlaps,
        "records": [_record_payload(record) for record in report.records],
        "rows": [label_fit_row_payload(row) for row in report.rows],
    }


def label_fit_failure(error: BaseException) -> LabelsDoNotFitError | None:
    """The labels-do-not-fit error that stopped a render, if ``error`` is or wraps one."""

    seen: set[int] = set()
    current: BaseException | None = error
    while current is not None and id(current) not in seen:
        if isinstance(current, LabelsDoNotFitError):
            return current
        seen.add(id(current))
        current = current.__cause__ or current.__context__
    return None


def _has_rows_to_compute(report: LabelFitReport) -> bool:
    return any(row.status == "not_computed" for row in report.rows)


class LabelFitRowSearches:
    """The last render of each mode whose report has rows to compute.

    Bounded: one render per mode, its drawn figures or the failure that drew
    nothing. A newer render of a mode releases the older one, and its tokens
    then answer with a ``not_computed`` row.
    """

    def __init__(self) -> None:
        self._generation = 0
        self._renders: dict[str, tuple[int, dict[int, Drawing | LabelsDoNotFitError]]] = {}

    def _keep(self, mode: str, sources: Sequence[Drawing | LabelsDoNotFitError | None]) -> list[str | None]:
        self._generation += 1
        kept = {index: source for index, source in enumerate(sources) if source is not None}
        self._renders[mode] = (self._generation, kept)
        return [f"{mode}:{self._generation}:{index}" if index in kept else None for index in range(len(sources))]

    def retain(self, mode: str, drawings: Sequence[Drawing]) -> list[str | None]:
        """Keep the drawings whose report has a row not computed; one token per drawing."""

        def computable(drawing: Drawing) -> Drawing | None:
            report = label_fit_report(drawing)
            return drawing if report is not None and _has_rows_to_compute(report) else None

        return self._keep(mode, [computable(drawing) for drawing in drawings])

    def retain_failure(self, mode: str, error: LabelsDoNotFitError) -> str | None:
        """Keep a failure whose report has a row not computed; its token."""

        computable = error.row_search is not None and _has_rows_to_compute(error.report)
        return self._keep(mode, [error if computable else None])[0]

    def failure_payload(self, error: LabelsDoNotFitError) -> LabelFitFailurePayload:
        """The failure's report and the token of its kept rows (by identity), if kept."""

        token = next(
            (
                f"{mode}:{generation}:{index}"
                for mode, (generation, kept) in self._renders.items()
                for index, source in kept.items()
                if source is error
            ),
            None,
        )
        return {"report": label_fit_report_payload(error.report), "rowToken": token}

    def row(self, token: str, kind: str) -> LabelFitRowPayload:
        """Compute one row of a kept figure or failure (no time limit)."""

        if kind not in _ROW_KINDS:
            raise ValidationError(f"Unknown label fit row kind: {kind}.", diagnostic=_UNKNOWN_ROW_KIND)
        row_kind = cast(LabelFitRowKind, kind)
        mode, _, rest = token.partition(":")
        generation_text, _, index_text = rest.partition(":")
        render = self._renders.get(mode)
        source = (
            render[1].get(int(index_text))
            if render is not None
            and generation_text.isdigit()
            and index_text.isdigit()
            and render[0] == int(generation_text)
            else None
        )
        row = label_fit_row(source, row_kind) if source is not None else None
        if row is None:
            return label_fit_row_payload(LabelFitRow(row_kind, "not_computed", reason=_STALE_REASON))
        return label_fit_row_payload(row)


__all__ = [
    "LABEL_FIT_PAYLOAD_SCHEMA",
    "LabelFitFailurePayload",
    "LabelFitRecordPayload",
    "LabelFitReportPayload",
    "LabelFitRowPayload",
    "LabelFitRowSearches",
    "label_fit_failure",
    "label_fit_report_payload",
    "label_fit_row_payload",
]
