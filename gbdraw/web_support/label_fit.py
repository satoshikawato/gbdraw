"""The label fit report as the Web reads it, and the rows it computes later.

Python owns every value and every word of the report: the Web dialog shows
the rows' titles, the settings they change, the records' counts, and the
canvas sizes as given, and nothing in JavaScript places labels or searches
values (LP-D43, LP-D48). Rows a bounded search did not reach arrive
``not_computed`` (OD-F14 for a drawn figure, the font scan bound when nothing
was drawn); :class:`LabelFitRowSearches` keeps the last such figure or failure
of each mode so the Web can compute one row through :func:`label_fit_row`,
searched for at most ``WEB_ROW_SEARCH_SECONDS`` (LP-D60).
"""

from __future__ import annotations

import json
import time
import warnings
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from typing import Literal, TextIO, TypeAlias, TypedDict, cast, get_args

from svgwrite import Drawing

from gbdraw.api.diagram import label_fit_report, label_fit_row
from gbdraw.core.numeric import parse_ascii_nonnegative_integer
from gbdraw.exceptions import LabelsDoNotFitError, LabelsMayNotFitWarning, ValidationError
from gbdraw.labels.fit import (
    ALL_OUTSIDE_SCOPE,
    CENTER_RESERVED_RADIUS_PATH,
    DEFINITION_FONT_PATH,
    FONT_SIZE_PATHS,
    PLACEMENT_PATH,
    RADIAL_PLACEMENT,
    ROW_TITLES,
    SCOPE_PATH,
    SEPARATE_FEATURES_PATH,
    SPACING_PATH,
    X_OFFSET_PATH,
    Y_OFFSET_PATH,
    HorizontalLabelKind,
    LabelCapacityRemedy,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelFitRowStatus,
    LabelFitValue,
    LabelKind,
    label_fit_notice,
    label_fit_record_counts,
    label_kind_title,
    row_note,
    shown_capacity_ratio,
    split_row_option_values,
)

LABEL_FIT_PAYLOAD_SCHEMA: Literal[1] = 1
# How long the Web searches one row on demand (LP-D60): Apply and Cancel wait
# at most this long for the row being computed.
WEB_ROW_SEARCH_SECONDS = 30.0
LabelFitMode: TypeAlias = Literal["circular", "linear"]
_MODES: frozenset[str] = frozenset(get_args(LabelFitMode))

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
    # Record Labels > Default font size: the Circular centre definition's font (D-12).
    DEFINITION_FONT_PATH: ("Default font size", "pt"),
    # Advanced > Center Reserved Radius: "narrow the inside tracks" (PD-OI-097).
    CENTER_RESERVED_RADIUS_PATH: ("Center Reserved Radius", "px"),
}
# A control whose name differs by mode: the Circular and the Linear checkbox
# of `canvas.resolve_overlaps`.
_MODE_SETTING_NAMES: Mapping[tuple[LabelFitMode, str], tuple[str, str]] = {
    ("circular", SEPARATE_FEATURES_PATH): ("Resolve Overlaps", ""),
    ("linear", SEPARATE_FEATURES_PATH): ("Resolve Feature Overlaps", ""),
}
# A word or switch value as the Web control names its option.
_VALUE_WORDS: Mapping[tuple[str, str | bool], str] = {
    (SCOPE_PATH, ALL_OUTSIDE_SCOPE): "Out",
    (PLACEMENT_PATH, RADIAL_PLACEMENT): "Radial",
    (SEPARATE_FEATURES_PATH, True): "On",
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


class LabelFitTrackOptionsPayload(TypedDict, total=False):
    """The canonical render request's Circular track options (``diagramOptions.tracks``) a row sets, by their
    request names: "narrow the inside tracks" sets ``centerReservedRadius`` (PD-OI-097)."""

    centerReservedRadius: float


class LabelFitRowPayload(TypedDict):
    kind: LabelFitRowKind
    title: str
    status: LabelFitRowStatus
    # The config overrides a row sets, by config path; set when offered.
    values: dict[str, bool | float | str] | None
    trackOptions: LabelFitTrackOptionsPayload
    settings: list[LabelFitSettingPayload]
    canvasSize: LabelFitCanvasSizePayload | None
    reason: str | None
    note: str | None


class LabelFitRecordPayload(TypedDict):
    recordId: str
    # The record's drawn index (``data-gbdraw-record-index``); None for a Linear record.
    recordIndex: int | None
    labelKind: LabelKind
    labelKindTitle: str
    labelCount: int
    fontSize: float
    overlapPairs: int
    crossingPairs: int
    definitionOverlaps: int
    regionOverlaps: int
    refused: bool
    # The room the labels need as a multiple of the room around the map, as shown (``shown_capacity_ratio``), when
    # they were not placed (CAPACITY-STOP); the Web words such a record from its values, not from ``counts``.
    capacityRatio: float | None
    counts: str


class LabelCapacityNoticePayload(TypedDict):
    """One side of a record whose labels may not fit (CAPACITY-WARN), as the Web words it while the render runs."""

    recordId: str
    labelKind: HorizontalLabelKind
    labelCount: int
    fontSize: float
    # As shown (``shown_capacity_ratio``).
    capacityRatio: float


class LabelFitReportPayload(TypedDict):
    schema: Literal[1]
    drawn: bool
    overlapPairs: int
    crossingPairs: int
    definitionOverlaps: int
    regionOverlaps: int
    notice: list[str]
    records: list[LabelFitRecordPayload]
    rows: list[LabelFitRowPayload]
    # The ways to make room for records that were not placed (CAPACITY-STOP), or when no row has a value (NO-ROW-WAYS,
    # listed once every row is known not to fit); the Web words each with its controls.
    remedies: list[LabelCapacityRemedy]


class LabelFitFailurePayload(TypedDict):
    """The report of a render that drew nothing, and the token of its rows to compute."""

    report: LabelFitReportPayload
    rowToken: str | None


def label_fit_mode(value: object) -> LabelFitMode:
    """The diagram mode of a Web render, whose controls a row's settings name."""

    if not isinstance(value, str) or value not in _MODES:
        raise ValidationError(f"Unknown label fit mode: {value!r}.", diagnostic=_UNKNOWN_ROW_KIND)
    return cast(LabelFitMode, value)


def _setting(path: str, value: LabelFitValue, mode: LabelFitMode) -> LabelFitSettingPayload:
    named = _MODE_SETTING_NAMES.get((mode, path)) or _SETTING_NAMES.get(path)
    word = isinstance(value, (str, bool))
    if named is None or (word and (path, value) not in _VALUE_WORDS):
        raise ValidationError(
            f"The label fit row sets an unknown config path or value: {path}.", diagnostic=_UNKNOWN_CONFIG_PATH
        )
    name, unit = named
    shown = _VALUE_WORDS[(path, value)] if isinstance(value, (str, bool)) else float(value)
    return {"name": name, "unit": unit, "value": shown}


def _settings(values: Mapping[str, LabelFitValue], mode: LabelFitMode) -> list[LabelFitSettingPayload]:
    settings: dict[str, LabelFitSettingPayload] = {}
    for path, value in values.items():
        setting = _setting(path, value, mode)
        settings.setdefault(setting["name"], setting)
    return list(settings.values())


def label_fit_row_payload(row: LabelFitRow, mode: LabelFitMode) -> LabelFitRowPayload:
    """One row with its title, the settings it changes (as the controls of ``mode`` name them), and its canvas size."""

    values = dict(row.values) if row.values is not None else None
    # Drawing again takes the row's settings as config overrides and track options (``split_row_option_values``).
    config, center_reserved_radius = split_row_option_values(values) if values else ({}, None)
    track_options: LabelFitTrackOptionsPayload = {}
    if center_reserved_radius is not None:
        track_options["centerReservedRadius"] = center_reserved_radius
    return {
        "kind": row.kind,
        "title": ROW_TITLES[row.kind],
        "status": row.status,
        "values": (
            {path: value if isinstance(value, (str, bool)) else float(value) for path, value in config.items()}
            if values
            else None
        ),
        "trackOptions": track_options,
        "settings": _settings(values, mode) if values else [],
        "canvasSize": (
            {"width": float(row.canvas_size[0]), "height": float(row.canvas_size[1])}
            if row.canvas_size is not None
            else None
        ),
        "reason": row.reason,
        # What a row that draws with the notice leaves (LP-D57).
        "note": row_note(row),
    }


def _record_payload(record: LabelFitRecord) -> LabelFitRecordPayload:
    return {
        "recordId": record.record_id,
        "recordIndex": record.record_index,
        "labelKind": record.label_kind,
        "labelKindTitle": label_kind_title(record),
        "labelCount": record.label_count,
        "fontSize": float(record.font_size),
        "overlapPairs": record.fit.overlap_pairs,
        "crossingPairs": record.fit.crossing_pairs,
        "definitionOverlaps": record.fit.definition_overlaps,
        "regionOverlaps": record.fit.region_overlaps,
        "refused": record.refused,
        "capacityRatio": shown_capacity_ratio(record.capacity) if record.capacity is not None else None,
        # A refused radial record has no counts; the reason stands in (OV-699).
        "counts": label_fit_record_counts(record),
    }


def label_fit_report_payload(report: LabelFitReport, mode: LabelFitMode) -> LabelFitReportPayload:
    """The bounded JSON form of a report: identifiers, counts, values, and Python's words."""

    return {
        "schema": LABEL_FIT_PAYLOAD_SCHEMA,
        "drawn": report.drawn,
        "overlapPairs": report.overlap_pairs,
        "crossingPairs": report.crossing_pairs,
        "definitionOverlaps": report.definition_overlaps,
        "regionOverlaps": report.region_overlaps,
        "notice": list(label_fit_notice(report)),
        "records": [_record_payload(record) for record in report.records],
        "rows": [label_fit_row_payload(row, mode) for row in report.rows],
        "remedies": list(report.room_remedies),
    }


def label_capacity_notice_payload(warning: LabelsMayNotFitWarning) -> LabelCapacityNoticePayload:
    """The values of a ``LabelsMayNotFitWarning`` the Web words in the UI language (OV-2179)."""

    return {
        "recordId": warning.record_id,
        "labelKind": warning.label_kind,
        "labelCount": warning.capacity.label_count,
        "fontSize": float(warning.capacity.font_size),
        "capacityRatio": shown_capacity_ratio(warning.capacity),
    }


@contextmanager
def label_capacity_notices(notice: Callable[[str], object]) -> Iterator[None]:
    """Send each ``LabelsMayNotFitWarning`` of a render to ``notice`` as the JSON text of its
    :class:`LabelCapacityNoticePayload`, when it is warned (CAPACITY-WARN): before the labels are placed, so the Web
    shows it while the render runs."""

    with warnings.catch_warnings():
        warnings.simplefilter("always", LabelsMayNotFitWarning)
        show_other = warnings.showwarning

        def show(
            message: Warning | str,
            category: type[Warning],
            filename: str,
            lineno: int,
            file: TextIO | None = None,
            line: str | None = None,
        ) -> None:
            if isinstance(message, LabelsMayNotFitWarning):
                notice(json.dumps(label_capacity_notice_payload(message)))
                return
            show_other(message, category, filename, lineno, file, line)

        warnings.showwarning = show
        yield


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


_Kept: TypeAlias = tuple[int, dict[int, Drawing | LabelsDoNotFitError]]


class LabelFitRowSearches:
    """The renders of each mode whose report has rows to compute.

    Bounded: at most two per mode, the last drawn figures (the shown Result)
    and the last failure that drew nothing. A newer drawn render of a mode
    releases both, a newer failure only the older failure (the shown Result
    stays), and released tokens answer with a ``not_computed`` row.
    """

    def __init__(self, *, row_search_seconds: float = WEB_ROW_SEARCH_SECONDS) -> None:
        self._row_search_seconds = row_search_seconds
        self._generation = 0
        self._drawn: dict[LabelFitMode, _Kept] = {}
        self._failed: dict[LabelFitMode, _Kept] = {}

    def _keep(
        self,
        kept_by_mode: dict[LabelFitMode, _Kept],
        mode: LabelFitMode,
        sources: Sequence[Drawing | LabelsDoNotFitError | None],
    ) -> list[str | None]:
        self._generation += 1
        kept = {index: source for index, source in enumerate(sources) if source is not None}
        kept_by_mode[mode] = (self._generation, kept)
        return [f"{mode}:{self._generation}:{index}" if index in kept else None for index in range(len(sources))]

    def _kept(self) -> Iterator[tuple[LabelFitMode, _Kept]]:
        yield from self._drawn.items()
        yield from self._failed.items()

    def retain(self, mode: LabelFitMode, drawings: Sequence[Drawing]) -> list[str | None]:
        """Keep the drawings whose report has a row not computed; one token per drawing."""

        def computable(drawing: Drawing) -> Drawing | None:
            report = label_fit_report(drawing)
            return drawing if report is not None and _has_rows_to_compute(report) else None

        self._failed.pop(mode, None)
        return self._keep(self._drawn, mode, [computable(drawing) for drawing in drawings])

    def retain_failure(self, mode: LabelFitMode, error: LabelsDoNotFitError) -> str | None:
        """Keep a failure whose report has a row not computed; its token."""

        computable = error.row_search is not None and _has_rows_to_compute(error.report)
        return self._keep(self._failed, mode, [error if computable else None])[0]

    def failure_payload(self, error: LabelsDoNotFitError, mode: LabelFitMode) -> LabelFitFailurePayload:
        """The failure's report and the token of its kept rows (by identity), if kept."""

        token = next(
            (
                f"{kept_mode}:{generation}:{index}"
                for kept_mode, (generation, kept) in self._kept()
                for index, source in kept.items()
                if source is error
            ),
            None,
        )
        return {"report": label_fit_report_payload(error.report, mode), "rowToken": token}

    def row(self, token: str, kind: str) -> LabelFitRowPayload:
        """Compute one row of a kept figure or failure, searched for at most the row bound."""

        if kind not in _ROW_KINDS:
            raise ValidationError(f"Unknown label fit row kind: {kind}.", diagnostic=_UNKNOWN_ROW_KIND)
        row_kind = cast(LabelFitRowKind, kind)
        mode_text, _, rest = token.partition(":")
        mode = label_fit_mode(mode_text)
        generation_text, _, index_text = rest.partition(":")
        wanted_generation = parse_ascii_nonnegative_integer(generation_text)
        wanted_index = parse_ascii_nonnegative_integer(index_text)
        source = (
            next(
                (
                    kept.get(wanted_index)
                    for kept_mode, (generation, kept) in self._kept()
                    if kept_mode == mode and generation == wanted_generation
                ),
                None,
            )
            if wanted_generation is not None and wanted_index is not None
            else None
        )
        deadline = time.monotonic() + self._row_search_seconds
        row = label_fit_row(source, row_kind, deadline=deadline) if source is not None else None
        if row is None:
            return label_fit_row_payload(LabelFitRow(row_kind, "not_computed", reason=_STALE_REASON), mode)
        if row.status == "not_computed":
            # The bound stopped the search: no value, which is neither "does not fit" nor stale.
            reason = f"No value found within {self._row_search_seconds:g} seconds."
            return label_fit_row_payload(LabelFitRow(row_kind, "not_computed", reason=reason), mode)
        return label_fit_row_payload(row, mode)


__all__ = [
    "LABEL_FIT_PAYLOAD_SCHEMA",
    "WEB_ROW_SEARCH_SECONDS",
    "LabelCapacityNoticePayload",
    "LabelFitFailurePayload",
    "LabelFitMode",
    "LabelFitRecordPayload",
    "LabelFitReportPayload",
    "LabelFitRowPayload",
    "LabelFitRowSearches",
    "label_capacity_notice_payload",
    "label_capacity_notices",
    "label_fit_failure",
    "label_fit_mode",
    "label_fit_report_payload",
    "label_fit_row_payload",
]
