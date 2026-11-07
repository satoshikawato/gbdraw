"""The size limits of the PNG and PDF formats, checked before an export starts, and the default PNG DPI.

A PNG side over ``PNG_MAX_PIXELS`` cannot be drawn (cairo refuses the surface). A PDF page side over
``PDF_MAX_PAGE_PT`` is written by cairo but is over what PDF viewers accept, so only a wrapped figure, new in
this format, is refused for it; an unwrapped PDF stays as it was. SVG, interactive SVG, EPS and PS have no limit.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import floor, isfinite

from svgwrite import Drawing

from gbdraw.core.figure_length import CSS_PX_PER_INCH
from gbdraw.exceptions import ExportError, ValidationError
from gbdraw.layout.linear_wrap import LinearWrapSummary, longer_row_lengths, wrap_length_option
from gbdraw.render.document import drawing_size

PNG_MAX_PIXELS = 32_767
# The DPI of a PNG whose caller names none (the CLI without --dpi, Python without dpi=, a Session without a
# readable ui.downloadDpi), and of the Web app's DPI menu: the Web reads it generated (DEFAULT_PNG_DPI).
DEFAULT_PNG_DPI: int = 300
PDF_MAX_PAGE_PT = 14_400
PDF_PT_PER_SVG_PX = 0.75

_EXPORT_WRAP_ATTR = "_gbdraw_export_wrap"


def bind_export_wrap(drawing: Drawing, wrap: LinearWrapSummary | None) -> None:
    """Keep what a wrapped figure drew with its drawing, for the check at export."""

    setattr(drawing, _EXPORT_WRAP_ATTR, wrap)


def _export_wrap(drawing: Drawing) -> LinearWrapSummary | None:
    wrap = getattr(drawing, _EXPORT_WRAP_ATTR, None)
    return wrap if isinstance(wrap, LinearWrapSummary) else None


def checked_png_dpi(dpi: float | None) -> float | None:
    """A PNG DPI as an export takes it: ``None`` (:data:`DEFAULT_PNG_DPI`) or a positive finite number."""

    if dpi is not None and (
        isinstance(dpi, bool) or not isinstance(dpi, (int, float)) or not isfinite(dpi) or dpi <= 0
    ):
        raise ValidationError(
            f"dpi must be a positive number; got {dpi!r}.",
            diagnostic={"code": "INPUT_INVALID", "reason": "POSITIVE"},
        )
    return dpi


def png_export_dpi(dpi: float | None) -> float:
    """The DPI a PNG is drawn at: ``dpi``, else :data:`DEFAULT_PNG_DPI`. Every PNG writer and check reads it here."""

    return DEFAULT_PNG_DPI if dpi is None else dpi


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

    ``wrap`` is what a wrapped figure drew, and it decides the advice. ``dpi`` is the PNG's DPI (``dpi / 96`` pixels
    per SVG pixel: the Web app's, the CLI ``--dpi`` and the ``dpi`` of the Python exports); ``None`` is
    :data:`DEFAULT_PNG_DPI`, which the PNG refusal names with how to set another. A PDF page never depends on it.
    The error's diagnostic carries the same facts for the Web: the limit and the size in the format's unit, and the
    advised row length and DPI when one fits.
    """

    png_dpi = png_export_dpi(dpi)
    if fmt == "png":
        scale = png_dpi / CSS_PX_PER_INCH
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
    message = f"{fmt.upper()} is limited to {limit:,} {unit} per side, and this figure is {width:,.0f} x {height:,.0f} {unit}"
    if fmt == "png":
        message += f" at {png_dpi:g} dpi" + (" (the default PNG DPI)" if dpi is None else "")
    message += "."
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
    if fmt == "png":
        lower_dpi = floor(limit * CSS_PX_PER_INCH / max(width_px, height_px))
        if lower_dpi >= 1:
            diagnostic["dpi"] = lower_dpi
            message += f" A DPI of {lower_dpi} or lower would fit"
            # The caller's DPI and the Session's (gbdraw.session.noting_session_png_dpi) say how they are set.
            message += "; --dpi on the command line or dpi= in Python sets it." if dpi is None else "."
    raise ExportError(f"{message} SVG has no size limit (-f svg).", diagnostic=diagnostic)


def check_drawing_export_limits(drawing: Drawing, formats: Sequence[str], *, dpi: float | None = None) -> None:
    """Check each requested format against a drawing, before anything is converted or written (``dpi``: the PNG's)."""

    width, height = drawing_size(drawing)
    wrap = _export_wrap(drawing)
    for fmt in formats:
        check_export_limit(fmt, width, height, wrap=wrap, dpi=dpi)


__all__ = [
    "DEFAULT_PNG_DPI",
    "PDF_MAX_PAGE_PT",
    "PNG_MAX_PIXELS",
    "bind_export_wrap",
    "check_drawing_export_limits",
    "check_export_limit",
    "checked_png_dpi",
    "png_export_dpi",
]
