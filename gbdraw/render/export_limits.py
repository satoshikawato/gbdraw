"""The size limits of the PNG and PDF formats, checked before an export starts.

A PNG side over ``PNG_MAX_PIXELS`` cannot be drawn (cairo refuses the surface). A PDF page side over
``PDF_MAX_PAGE_PT`` is written by cairo but is over what PDF viewers accept, so only a wrapped figure, new in
this format, is refused for it; an unwrapped PDF stays as it was. SVG, interactive SVG, EPS and PS have no limit.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import floor

from svgwrite import Drawing

from gbdraw.exceptions import ExportError
from gbdraw.layout.linear_wrap import LinearWrapSummary, longer_row_lengths, wrap_length_option
from gbdraw.render.document import drawing_size

PNG_MAX_PIXELS = 32_767
PDF_MAX_PAGE_PT = 14_400
CSS_PX_PER_INCH = 96.0
PDF_PT_PER_SVG_PX = 0.75

_EXPORT_WRAP_ATTR = "_gbdraw_export_wrap"


def bind_export_wrap(drawing: Drawing, wrap: LinearWrapSummary | None) -> None:
    """Keep what a wrapped figure drew with its drawing, for the check at export."""

    setattr(drawing, _EXPORT_WRAP_ATTR, wrap)


def _export_wrap(drawing: Drawing) -> LinearWrapSummary | None:
    wrap = getattr(drawing, _EXPORT_WRAP_ATTR, None)
    return wrap if isinstance(wrap, LinearWrapSummary) else None


def _row_length_advice(wrap: LinearWrapSummary, width: float, height: float, limit: int) -> tuple[int, float] | None:
    """The next longer 1-2-5 row length that fits and its estimated height, from the drawn height per row: the rows
    share the height. ``None`` when no longer row length would fit."""

    if width <= limit:
        drawn_rows = wrap.row_count
        for length in longer_row_lengths(wrap.row_length, max(wrap.record_lengths)):
            estimate = height * wrap.rows_at(length) / drawn_rows
            if estimate <= limit:
                return length, estimate
    return None


def check_export_limit(
    fmt: str,
    width_px: float,
    height_px: float,
    *,
    wrap: LinearWrapSummary | None = None,
    dpi: float | None = None,
) -> None:
    """Raise ``ExportError`` when a figure of this SVG size is over ``fmt``'s limit.

    ``wrap`` is what a wrapped figure drew, and it decides the advice. ``dpi`` is given only by a caller whose PNG
    scales by it (the Web app); gbdraw's CairoSVG PNG is one pixel per SVG pixel whatever its DPI, and a PDF page
    never depends on it. The error's diagnostic carries the same facts for the Web: the limit and the size in the
    format's unit, and the advised row length and DPI when one fits.
    """

    if fmt == "png":
        scale = (dpi if dpi is not None else CSS_PX_PER_INCH) / CSS_PX_PER_INCH
        limit, unit = PNG_MAX_PIXELS, "px"
        # CairoSVG rounds a PNG's size to whole pixels.
        over = max(round(width_px * scale), round(height_px * scale)) > limit
    elif fmt == "pdf" and wrap is not None:
        scale = PDF_PT_PER_SVG_PX
        limit, unit = PDF_MAX_PAGE_PT, "pt"
        over = max(width_px, height_px) * scale > limit
    else:
        return
    if not over:
        return
    width, height = width_px * scale, height_px * scale
    diagnostic: dict[str, object] = {
        "code": "EXPORT_LIMIT",
        "reason": "PNG_SIZE" if fmt == "png" else "PDF_SIZE",
        "limit": limit,
        "width": round(width),
        "height": round(height),
    }
    message = (
        f"{fmt.upper()} is limited to {limit:,} {unit} per side, and this figure is {width:,.0f} x {height:,.0f} {unit}."
    )
    if wrap is not None:
        advice = _row_length_advice(wrap, width, height, limit)
        if advice is None:
            message += " No longer row length would fit (estimated from the drawn height)."
        else:
            length, estimate = advice
            diagnostic["rowLength"] = length
            message += (
                f" A longer row length should fit: {wrap_length_option(length)} is about {estimate:,.0f} {unit} "
                "tall (estimated from the drawn height)."
            )
    if dpi is not None and fmt == "png":
        lower_dpi = floor(limit * CSS_PX_PER_INCH / max(width_px, height_px))
        if lower_dpi >= 1:
            diagnostic["dpi"] = lower_dpi
            message += f" A DPI of {lower_dpi} or lower would fit."
    raise ExportError(f"{message} SVG has no size limit (-f svg).", diagnostic=diagnostic)


def check_drawing_export_limits(drawing: Drawing, formats: Sequence[str]) -> None:
    """Check each requested format against a drawing, before anything is converted or written."""

    width, height = drawing_size(drawing)
    wrap = _export_wrap(drawing)
    for fmt in formats:
        check_export_limit(fmt, width, height, wrap=wrap)


__all__ = [
    "PDF_MAX_PAGE_PT",
    "PNG_MAX_PIXELS",
    "bind_export_wrap",
    "check_drawing_export_limits",
    "check_export_limit",
]
