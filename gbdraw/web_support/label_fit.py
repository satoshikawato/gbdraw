"""The label fit report as the Web reads it, and the rows it computes later.

Python owns every value and every word of the report: the Web dialog shows
the rows' titles, the settings they change, and the canvas sizes as given,
and nothing in JavaScript places labels or searches values (LP-D43, LP-D48).
A figure drawn with crossing leaders comes with rows that were not searched
(OD-F14); :class:`LabelFitRowSearches` keeps the last such figures of each
mode so the Web can compute one row through :func:`label_fit_row`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict, cast, get_args

from svgwrite import Drawing

from gbdraw.api.diagram import label_fit_report, label_fit_row
from gbdraw.exceptions import LabelsDoNotFitError, ValidationError
from gbdraw.labels.fit import (
    FONT_SIZE_PATHS,
    ROW_TITLES,
    SPACING_PATH,
    X_OFFSET_PATH,
    Y_OFFSET_PATH,
)
from gbdraw.labels.fit import (
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelFitRowStatus,
    LabelKind,
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
}
_LABEL_KIND_TITLES: Mapping[LabelKind, str] = {"circular_outer": "outer labels"}
_ROW_KINDS: frozenset[str] = frozenset(get_args(LabelFitRowKind))
_STALE_REASON = "Generate again to compute this choice."
_UNKNOWN_CONFIG_PATH = {"code": "INPUT_INVALID", "field": "configOverrides", "reason": "UNKNOWN_CONFIG_PATH"}
_UNKNOWN_ROW_KIND = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


class LabelFitSettingPayload(TypedDict):
    name: str
    unit: str
    value: float


class LabelFitCanvasSizePayload(TypedDict):
    width: float
    height: float


class LabelFitRowPayload(TypedDict):
    kind: LabelFitRowKind
    title: str
    status: LabelFitRowStatus
    values: dict[str, float] | None
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


class LabelFitReportPayload(TypedDict):
    schema: Literal[1]
    drawn: bool
    overlapPairs: int
    crossingPairs: int
    records: list[LabelFitRecordPayload]
    rows: list[LabelFitRowPayload]


def _settings(values: Mapping[str, float]) -> list[LabelFitSettingPayload]:
    settings: dict[str, LabelFitSettingPayload] = {}
    for path, value in values.items():
        if path not in _SETTING_NAMES:
            raise ValidationError(
                f"The label fit row sets an unknown config path: {path}.", diagnostic=_UNKNOWN_CONFIG_PATH
            )
        name, unit = _SETTING_NAMES[path]
        settings.setdefault(name, {"name": name, "unit": unit, "value": float(value)})
    return list(settings.values())


def label_fit_row_payload(row: LabelFitRow) -> LabelFitRowPayload:
    """One row with its title, the settings it changes, and its canvas size."""

    values = dict(row.values) if row.values is not None else None
    return {
        "kind": row.kind,
        "title": ROW_TITLES[row.kind],
        "status": row.status,
        "values": {path: float(value) for path, value in values.items()} if values else None,
        "settings": _settings(values) if values else [],
        "canvasSize": (
            {"width": float(row.canvas_size[0]), "height": float(row.canvas_size[1])}
            if row.canvas_size is not None
            else None
        ),
        "reason": row.reason,
    }


def label_fit_report_payload(report: LabelFitReport) -> LabelFitReportPayload:
    """The bounded JSON form of a report: identifiers, counts, values, and Python's words."""

    return {
        "schema": LABEL_FIT_PAYLOAD_SCHEMA,
        "drawn": report.drawn,
        "overlapPairs": report.overlap_pairs,
        "crossingPairs": report.crossing_pairs,
        "records": [
            {
                "recordId": record.record_id,
                "labelKind": record.label_kind,
                "labelKindTitle": _LABEL_KIND_TITLES[record.label_kind],
                "labelCount": record.label_count,
                "fontSize": float(record.font_size),
                "overlapPairs": record.fit.overlap_pairs,
                "crossingPairs": record.fit.crossing_pairs,
            }
            for record in report.records
        ],
        "rows": [label_fit_row_payload(row) for row in report.rows],
    }


def label_fit_error_payload(error: BaseException) -> LabelFitReportPayload | None:
    """The report of overlapping labels that stopped a render, if ``error`` is one."""

    seen: set[int] = set()
    current: BaseException | None = error
    while current is not None and id(current) not in seen:
        if isinstance(current, LabelsDoNotFitError):
            return label_fit_report_payload(current.report)
        seen.add(id(current))
        current = current.__cause__ or current.__context__
    return None


class LabelFitRowSearches:
    """The figures of the last render of each mode whose report has rows to compute.

    Bounded: one render per mode. A newer render of a mode releases the older
    figures, and their tokens then answer with a ``not_computed`` row.
    """

    def __init__(self) -> None:
        self._generation = 0
        self._renders: dict[str, tuple[int, dict[int, Drawing]]] = {}

    def retain(self, mode: str, drawings: Sequence[Drawing]) -> list[str | None]:
        """Keep the drawings whose report has a row not computed; one token per drawing."""

        self._generation += 1
        kept: dict[int, Drawing] = {}
        tokens: list[str | None] = []
        for index, drawing in enumerate(drawings):
            report = label_fit_report(drawing)
            if report is not None and any(row.status == "not_computed" for row in report.rows):
                kept[index] = drawing
                tokens.append(f"{mode}:{self._generation}:{index}")
            else:
                tokens.append(None)
        self._renders[mode] = (self._generation, kept)
        return tokens

    def row(self, token: str, kind: str) -> LabelFitRowPayload:
        """Compute one row of a kept figure (no time limit)."""

        if kind not in _ROW_KINDS:
            raise ValidationError(f"Unknown label fit row kind: {kind}.", diagnostic=_UNKNOWN_ROW_KIND)
        row_kind = cast(LabelFitRowKind, kind)
        mode, _, rest = token.partition(":")
        generation_text, _, index_text = rest.partition(":")
        render = self._renders.get(mode)
        drawing = (
            render[1].get(int(index_text))
            if render is not None
            and generation_text.isdigit()
            and index_text.isdigit()
            and render[0] == int(generation_text)
            else None
        )
        row = label_fit_row(drawing, row_kind) if drawing is not None else None
        if row is None:
            return label_fit_row_payload(LabelFitRow(row_kind, "not_computed", reason=_STALE_REASON))
        return label_fit_row_payload(row)


__all__ = [
    "LABEL_FIT_PAYLOAD_SCHEMA",
    "LabelFitRecordPayload",
    "LabelFitReportPayload",
    "LabelFitRowPayload",
    "LabelFitRowSearches",
    "label_fit_error_payload",
    "label_fit_report_payload",
    "label_fit_row_payload",
]
