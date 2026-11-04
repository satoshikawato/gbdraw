"""The Web export's size check at the Worker boundary (W-26, D-04).

``gbdraw.render.export_limits`` owns the PNG and PDF limits and the advice. This function only reads the browser's
request: the SVG size of the figure, the PNG DPI, and the wrap length leaf of the exported Result's committed request,
as stored, with each record's source length and requested region. The config validator reads the leaf, the drawing's
crop rule (``check_start_end_coords``) gives each record's drawn bases, and whether a record is cut and the rows follow
from them here, by ``gbdraw.layout.linear_wrap`` as the drawing did, so a Result without an input summary (a loaded
Session, whatever wrote its leaf, a cropped record included) is checked as the CLI checks it (D-52, D-55, D-70). A
figure over the limit raises the ``EXPORT_LIMIT`` error.

An SVG has no limit; the check of an SVG of a Result drawn to a requested size answers the root ``width`` and
``height`` its file states, by ``export_root_size`` (``gbdraw/render/document.py``, the CLI's owner) from the
committed request's size leaves as stored and the exported root's px size. A Result drawn at Auto is not checked.
"""

from __future__ import annotations

import json
from typing import Literal, TypedDict

from gbdraw.api.input_summary import RecordPayload
from gbdraw.config.models.canvas import (
    FIGURE_HEIGHT_PATH,
    FIGURE_WIDTH_PATH,
    figure_length_from_config,
    linear_wrap_length_from_config,
)
from gbdraw.core.figure_length import RequestedFigureSize, parse_figure_length
from gbdraw.crop_genbank import check_start_end_coords
from gbdraw.layout.linear_wrap import LinearWrapSummary, resolve_wrap_length
from gbdraw.render.document import export_root_size
from gbdraw.render.export_limits import check_export_limit

_FORMATS: dict[str, Literal["png", "pdf"]] = {"png": "png", "pdf": "pdf"}


class _WrapPayload(TypedDict):
    """The exported Result's wrap: its committed request's ``canvas.linear.wrap_length`` leaf as stored (any value the
    config validator reads, ``"5kb"`` and ``None`` included) and each record it drew, at least one, as its input
    summary names it (source length and requested region)."""

    requested: str | int | float | None
    records: list[RecordPayload]


class _FigureSizePayload(TypedDict):
    """The exported Result's committed request leaves ``canvas.figure_width`` and ``canvas.figure_height`` as
    stored (``None`` for an Auto side; any other value is read by the config validator, which refuses it)."""

    width: object
    height: object


class _ExportLimitRequired(TypedDict):
    format: str
    width: float
    height: float
    dpi: float | None
    linearWrap: _WrapPayload | None


class ExportLimitRequest(_ExportLimitRequired, total=False):
    """``format`` ``png``, ``pdf`` or ``svg``; ``width`` and ``height`` in SVG pixels; ``dpi`` the PNG's (``None``
    otherwise); ``linearWrap`` the exported Result's wrap, ``None`` when it was drawn with wrap off (or is no Linear
    Result); ``figureSize`` the Result's requested size, which an SVG check names (PNG and PDF send none)."""

    figureSize: _FigureSizePayload | None


class ExportRootSize(TypedDict):
    """The root ``width`` and ``height`` an exported SVG states (``export_root_size``)."""

    width: str
    height: str


def _requested_size(payload: _FigureSizePayload | None) -> RequestedFigureSize:
    """The size the leaves ask for, read by the config validator as a Generate reads them."""

    if payload is None:
        return RequestedFigureSize()
    return RequestedFigureSize(
        parse_figure_length(figure_length_from_config(payload["width"], path=FIGURE_WIDTH_PATH)),
        parse_figure_length(figure_length_from_config(payload["height"], path=FIGURE_HEIGHT_PATH)),
    )


def _svg_root(request: ExportLimitRequest) -> ExportRootSize:
    root = export_root_size(
        _requested_size(request.get("figureSize")), float(request["width"]), float(request["height"])
    )
    if root is None:
        raise ValueError("An SVG export check names the requested figure size of a sized Result.")
    return {"width": root["width"], "height": root["height"]}


def _drawn_length(record: RecordPayload, ordinal: int) -> int:
    """The bases the drawing draws of a record: the whole source, or its region as the drawing crops it."""

    length = int(record["sourceLength"])
    region = record["region"]
    if region is None:
        return length
    start_0, end_0 = check_start_end_coords(length, int(region["start"]), int(region["end"]), f"#{ordinal}")
    return end_0 - start_0


def _wrap(payload: _WrapPayload | None) -> LinearWrapSummary | None:
    """What the Result drew with its wrap: ``None`` when no record is cut, as ``resolve_wrap`` decides."""

    if payload is None:
        return None
    requested = linear_wrap_length_from_config(payload["requested"])
    lengths = tuple(_drawn_length(record, ordinal) for ordinal, record in enumerate(payload["records"], start=1))
    row_length = resolve_wrap_length(requested, lengths)
    return None if row_length is None else LinearWrapSummary(row_length, requested, lengths)


def check_web_export_limit_json(request_json: str) -> str:
    """``{"ok": true}`` when the figure fits its format; otherwise ``check_export_limit``'s ``ExportError``.

    For an SVG, ``{"ok": true, "root": {"width", "height"}}``: the root size its file states.
    """

    request: ExportLimitRequest = json.loads(str(request_json))
    if request["format"] == "svg":
        return json.dumps({"ok": True, "root": _svg_root(request)})
    fmt = _FORMATS.get(str(request["format"]))
    if fmt is None:
        raise ValueError(f"Unknown export format: {request['format']}")
    dpi = request["dpi"]
    if fmt == "png" and dpi is None:
        raise ValueError("A PNG export names its DPI.")
    check_export_limit(
        fmt,
        float(request["width"]),
        float(request["height"]),
        wrap=_wrap(request["linearWrap"]),
        # The Web PNG scales by its DPI; a PDF page never depends on it.
        dpi=float(dpi) if fmt == "png" and dpi is not None else None,
    )
    return json.dumps({"ok": True})


__all__ = ["ExportLimitRequest", "ExportRootSize", "check_web_export_limit_json"]
