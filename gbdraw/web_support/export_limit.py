"""The Web export's size check at the Worker boundary (W-26, D-04).

``gbdraw.render.export_limits`` owns the PNG and PDF limits and the advice. This function only reads the browser's
request: the SVG size of the figure, the PNG DPI, and what a wrapped Result drew (the input summary's
``linearWrap``, as Python sent it). A figure over the limit raises the ``EXPORT_LIMIT`` error.
"""

from __future__ import annotations

import json
from typing import Literal, TypedDict

from gbdraw.layout.linear_wrap import LinearWrapSummary
from gbdraw.render.export_limits import check_export_limit

_FORMATS: dict[str, Literal["png", "pdf"]] = {"png": "png", "pdf": "pdf"}


class _WrapPayload(TypedDict):
    """The fields of the input summary's ``linearWrap`` (``gbdraw.api.input_summary.LinearWrapPayload``) the check
    reads; the rows follow from them."""

    rowLength: int
    requested: int | Literal["auto"]
    recordLengths: list[int]


class ExportLimitRequest(TypedDict):
    """``format`` ``png`` or ``pdf``; ``width`` and ``height`` in SVG pixels; ``dpi`` the PNG's (``None`` for a
    PDF); ``linearWrap`` what the exported Result drew, ``None`` when it cut no record."""

    format: str
    width: float
    height: float
    dpi: float | None
    linearWrap: _WrapPayload | None


def _wrap(payload: _WrapPayload | None) -> LinearWrapSummary | None:
    if payload is None:
        return None
    requested = payload["requested"]
    return LinearWrapSummary(
        int(payload["rowLength"]),
        "auto" if requested == "auto" else int(requested),
        tuple(int(length) for length in payload["recordLengths"]),
    )


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
