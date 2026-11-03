"""The Web export's size check at the Worker boundary (W-26, D-04).

``gbdraw.render.export_limits`` owns the PNG and PDF limits and the advice. This function only reads the browser's
request: the SVG size of the figure, the PNG DPI, and the wrap length leaf of the exported Result's committed request,
as stored, with its records' lengths. The config validator reads the leaf, and whether a record is cut and the rows
follow from them here, by ``gbdraw.layout.linear_wrap`` as the drawing did, so a Result without an input summary (a
loaded Session, whatever wrote its leaf) is checked as the CLI checks it (D-52, D-55). A figure over the limit raises
the ``EXPORT_LIMIT`` error.
"""

from __future__ import annotations

import json
from typing import Literal, TypedDict

from gbdraw.config.models.canvas import linear_wrap_length_from_config
from gbdraw.layout.linear_wrap import LinearWrapSummary, resolve_wrap_length
from gbdraw.render.export_limits import check_export_limit

_FORMATS: dict[str, Literal["png", "pdf"]] = {"png": "png", "pdf": "pdf"}


class _WrapPayload(TypedDict):
    """The exported Result's wrap: its committed request's ``canvas.linear.wrap_length`` leaf as stored (any value the
    config validator reads, ``"5kb"`` and ``None`` included) and the bases of each record it drew, at least one."""

    requested: str | int | float | None
    recordLengths: list[int]


class ExportLimitRequest(TypedDict):
    """``format`` ``png`` or ``pdf``; ``width`` and ``height`` in SVG pixels; ``dpi`` the PNG's (``None`` for a
    PDF); ``linearWrap`` the exported Result's wrap, ``None`` when it was drawn with wrap off (or is no Linear
    Result)."""

    format: str
    width: float
    height: float
    dpi: float | None
    linearWrap: _WrapPayload | None


def _wrap(payload: _WrapPayload | None) -> LinearWrapSummary | None:
    """What the Result drew with its wrap: ``None`` when no record is cut, as ``resolve_wrap`` decides."""

    if payload is None:
        return None
    requested = linear_wrap_length_from_config(payload["requested"])
    lengths = tuple(int(length) for length in payload["recordLengths"])
    row_length = resolve_wrap_length(requested, lengths)
    return None if row_length is None else LinearWrapSummary(row_length, requested, lengths)


def check_web_export_limit_json(request_json: str) -> str:
    """``{"ok": true}`` when the figure fits its format; otherwise ``check_export_limit``'s ``ExportError``."""

    request: ExportLimitRequest = json.loads(str(request_json))
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


__all__ = ["ExportLimitRequest", "check_web_export_limit_json"]
