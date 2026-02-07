"""Render/export helpers for the public API layer."""

from __future__ import annotations

import logging
import os
import sys
from typing import Iterable, Sequence

from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.render import export as _export  # type: ignore[reportMissingImports]

logger = logging.getLogger(__name__)

_ACCEPTED_FORMATS = {"svg", "png", "eps", "ps", "pdf"}


def _normalize_formats(formats: Sequence[str] | str | None) -> list[str]:
    if formats is None:
        return ["svg"]
    if isinstance(formats, str):
        return _export.parse_formats(formats)
    out: list[str] = []
    for fmt in formats:
        if fmt is None:
            continue
        value = str(fmt).strip().lower()
        if not value:
            continue
        if value not in _ACCEPTED_FORMATS:
            logger.warning("WARNING: Unaccepted/unrecognized output file format: %s", value)
            continue
        if value not in out:
            out.append(value)
    if not out:
        logger.warning("WARNING: No valid output file format was specified; generate a PNG file.")
        return ["png"]
    return out


def _resolve_base_prefix(canvas: Drawing, output_prefix: str | None, output_dir: str | None) -> str:
    base_from_canvas = os.path.splitext(getattr(canvas, "filename", "") or "out.svg")[0]
    base_prefix = output_prefix or base_from_canvas or "out"
    base_prefix = os.path.splitext(base_prefix)[0]

    if output_dir:
        if output_prefix and os.path.dirname(output_prefix) not in {"", "."}:
            raise ValidationError(
                "output_dir cannot be combined with a path-like output_prefix; choose one."
            )
        os.makedirs(output_dir, exist_ok=True)
        base_prefix = os.path.join(output_dir, os.path.basename(base_prefix))

    return base_prefix


def _ensure_overwrite_ok(paths: Iterable[str], overwrite: bool) -> None:
    if overwrite:
        return
    existing = [path for path in paths if os.path.exists(path)]
    if existing:
        raise ValidationError(
            "Output file(s) already exist: " + ", ".join(existing) + ". Use overwrite=True to replace."
        )


def save_figure_to(
    canvas: Drawing,
    formats: Sequence[str] | str | None,
    *,
    output_dir: str | None = None,
    output_prefix: str | None = None,
    overwrite: bool = False,
) -> list[str]:
    """Save a figure to an explicit output directory/prefix.

    This always writes an SVG, then optionally converts to other formats using CairoSVG
    (when available), mirroring the CLI behavior.
    """

    fmt_list = _normalize_formats(formats)
    base_prefix = _resolve_base_prefix(canvas, output_prefix, output_dir)
    svg_filename = f"{base_prefix}.svg"

    # Always save SVG (base format).
    output_paths = [svg_filename]
    for fmt in fmt_list:
        if fmt == "svg":
            continue
        output_paths.append(f"{base_prefix}.{fmt}")

    _ensure_overwrite_ok(output_paths, overwrite)

    canvas.saveas(svg_filename)
    logger.info("Generated SVG: %s", svg_filename)

    formats_to_process = [f for f in fmt_list if f != "svg"]
    if not formats_to_process:
        return output_paths

    if "pyodide" in sys.modules:
        logger.info("Running in WebAssembly: Image conversion will be handled by the browser.")
        return output_paths

    if not _export.CAIROSVG_AVAILABLE:
        missing_formats = ", ".join([f.upper() for f in formats_to_process])
        logger.warning(
            "Skipping generation of: %s\n   CairoSVG is not installed.\n   To enable PNG/PDF support, run: pip install gbdraw[export]",
            missing_formats,
        )
        return output_paths

    try:
        svg_string = canvas.tostring()
        for fmt in formats_to_process:
            out_file = f"{base_prefix}.{fmt}"
            if fmt == "png":
                _export.cairosvg.svg2png(bytestring=svg_string.encode("utf-8"), write_to=out_file)
            elif fmt == "pdf":
                _export.cairosvg.svg2pdf(bytestring=svg_string.encode("utf-8"), write_to=out_file)
            elif fmt == "ps":
                _export.cairosvg.svg2ps(bytestring=svg_string.encode("utf-8"), write_to=out_file)
            elif fmt == "eps":
                _export.cairosvg.svg2ps(bytestring=svg_string.encode("utf-8"), write_to=out_file)
            logger.info("Generated %s: %s", fmt.upper(), out_file)
    except Exception as exc:
        logger.error("Failed to generate images using CairoSVG: %s", exc)

    return output_paths


def render_to_bytes(canvas: Drawing, fmt: str) -> bytes:
    """Render a canvas to bytes (SVG always; PNG/PDF/PS/EPS require CairoSVG)."""

    fmt_norm = str(fmt).strip().lower().lstrip(".")
    if not fmt_norm:
        raise ValidationError("Format must be specified (e.g., 'svg', 'png').")

    if fmt_norm == "svg":
        return canvas.tostring().encode("utf-8")

    if "pyodide" in sys.modules:
        raise ValidationError("Binary export is not available under WebAssembly (pyodide).")

    if not _export.CAIROSVG_AVAILABLE:
        raise ValidationError(
            "CairoSVG is not installed. Install with: pip install gbdraw[export]"
        )

    svg_string = canvas.tostring().encode("utf-8")
    if fmt_norm == "png":
        return _export.cairosvg.svg2png(bytestring=svg_string)
    if fmt_norm == "pdf":
        return _export.cairosvg.svg2pdf(bytestring=svg_string)
    if fmt_norm in {"ps", "eps"}:
        return _export.cairosvg.svg2ps(bytestring=svg_string)

    raise ValidationError(f"Unsupported format: {fmt}")


__all__ = ["parse_formats", "render_to_bytes", "save_figure", "save_figure_to"]

# Re-exported for convenience.
parse_formats = _export.parse_formats
save_figure = _export.save_figure
