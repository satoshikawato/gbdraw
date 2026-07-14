"""Render/export helpers for the public API layer."""

from __future__ import annotations

import logging
import os
import sys
from typing import Iterable, Sequence

from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ExportError, GbdrawError, ValidationError  # type: ignore[reportMissingImports]
from gbdraw.render import export as _export  # type: ignore[reportMissingImports]
from gbdraw.render.formats import (
    CAIROSVG_FORMATS,
    INTERACTIVE_SVG_FORMAT,
    SVG_FORMAT,
    classify_formats,
    normalize_format_sequence,
    normalize_format_token,
    resolve_format_output_path,
    resolve_output_paths,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

logger = logging.getLogger(__name__)


def _normalize_formats(formats: Sequence[str] | str | None) -> list[str]:
    return normalize_format_sequence(formats, logger=logger)


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
        for path in paths:
            if not os.path.exists(path):
                continue
            try:
                os.remove(path)
            except OSError as exc:
                raise ExportError(f"Could not replace existing output file: {path}") from exc
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
    interactive_context: InteractiveSvgContext | None = None,
) -> list[str]:
    """Save a figure to an explicit output directory/prefix.

    This always writes an SVG, then optionally converts to other formats using
    CairoSVG. Unlike the CLI-oriented :func:`save_figure`, explicitly requested
    formats are strict: failure to generate any one of them raises an exception.
    """

    fmt_list = _normalize_formats(formats)
    base_prefix = _resolve_base_prefix(canvas, output_prefix, output_dir)
    svg_filename = resolve_format_output_path(base_prefix, SVG_FORMAT)

    output_paths = resolve_output_paths(base_prefix, fmt_list, include_base_svg=True)

    _ensure_overwrite_ok(output_paths, overwrite)

    try:
        canvas.saveas(svg_filename)
    except Exception as exc:
        raise ExportError(f"Failed to generate SVG: {exc}") from exc
    if not os.path.isfile(svg_filename):
        raise ExportError(f"SVG export did not create the requested file: {svg_filename}")
    logger.info("Generated SVG: %s", svg_filename)

    classification = classify_formats(fmt_list)
    svg_source = canvas.tostring()
    if classification.interactive:
        interactive_filename = resolve_format_output_path(base_prefix, INTERACTIVE_SVG_FORMAT)
        try:
            interactive_svg = enrich_svg(svg_source, context=interactive_context)
        except GbdrawError:
            raise
        except Exception as exc:
            raise GbdrawError(f"Interactive SVG export failed: {exc}") from exc
        try:
            with open(interactive_filename, "w", encoding="utf-8") as handle:
                handle.write(interactive_svg)
        except OSError as exc:
            raise ExportError(f"Failed to write interactive SVG: {exc}") from exc
        if not os.path.isfile(interactive_filename):
            raise ExportError(
                f"Interactive SVG export did not create the requested file: {interactive_filename}"
            )
        logger.info("Generated interactive SVG: %s", interactive_filename)

    formats_to_process = list(classification.cairosvg)
    if not formats_to_process:
        return output_paths

    if "pyodide" in sys.modules:
        raise ValidationError(
            "Binary file export is not available under WebAssembly (pyodide); "
            "browser-side conversion does not create local output paths."
        )

    try:
        cairosvg_module = _export.get_cairosvg()
    except ImportError as exc:
        missing_formats = ", ".join([f.upper() for f in formats_to_process])
        raise ValidationError(
            f"Cannot generate {missing_formats}: CairoSVG is not installed. "
            "Install with: pip install gbdraw[export]"
        ) from exc

    svg_bytes = svg_source.encode("utf-8")
    for fmt in formats_to_process:
        out_file = resolve_format_output_path(base_prefix, fmt)
        try:
            if fmt not in CAIROSVG_FORMATS:
                continue
            if fmt == "png":
                cairosvg_module.svg2png(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "pdf":
                cairosvg_module.svg2pdf(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "ps":
                cairosvg_module.svg2ps(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "eps":
                cairosvg_module.svg2ps(bytestring=svg_bytes, write_to=out_file)
        except Exception as exc:
            raise ExportError(f"Failed to generate {fmt.upper()}: {exc}") from exc
        if not os.path.isfile(out_file):
            raise ExportError(
                f"{fmt.upper()} export did not create the requested file: {out_file}"
            )
        logger.info("Generated %s: %s", fmt.upper(), out_file)

    generated_paths = [path for path in output_paths if os.path.isfile(path)]
    if len(generated_paths) != len(output_paths):
        missing = [path for path in output_paths if path not in generated_paths]
        raise ExportError("Export completed without creating: " + ", ".join(missing))
    return generated_paths


def render_to_bytes(
    canvas: Drawing,
    fmt: str,
    *,
    interactive_context: InteractiveSvgContext | None = None,
) -> bytes:
    """Render a canvas to bytes (SVG always; PNG/PDF/PS/EPS require CairoSVG)."""

    fmt_norm = normalize_format_token(fmt)
    if not fmt_norm:
        raise ValidationError("Format must be specified (e.g., 'svg', 'png').")

    if fmt_norm == "svg":
        return canvas.tostring().encode("utf-8")
    if fmt_norm == INTERACTIVE_SVG_FORMAT:
        return enrich_svg(canvas.tostring(), context=interactive_context).encode("utf-8")
    if fmt_norm not in CAIROSVG_FORMATS:
        raise ValidationError(f"Unsupported format: {fmt}")

    if "pyodide" in sys.modules:
        raise ValidationError("Binary export is not available under WebAssembly (pyodide).")

    try:
        cairosvg_module = _export.get_cairosvg()
    except ImportError as exc:
        raise ValidationError(
            "CairoSVG is not installed. Install with: pip install gbdraw[export]"
        ) from exc

    svg_string = canvas.tostring().encode("utf-8")
    try:
        if fmt_norm == "png":
            rendered = cairosvg_module.svg2png(bytestring=svg_string)
        elif fmt_norm == "pdf":
            rendered = cairosvg_module.svg2pdf(bytestring=svg_string)
        elif fmt_norm in {"ps", "eps"}:
            rendered = cairosvg_module.svg2ps(bytestring=svg_string)
        else:
            raise ValidationError(f"Unsupported format: {fmt}")
    except GbdrawError:
        raise
    except Exception as exc:
        raise ExportError(f"Failed to render {fmt_norm.upper()} bytes: {exc}") from exc
    if not isinstance(rendered, (bytes, bytearray)):
        raise ExportError(
            f"{fmt_norm.upper()} byte export returned no binary payload."
        )
    return bytes(rendered)


__all__ = ["parse_formats", "render_to_bytes", "save_figure", "save_figure_to"]

# Re-exported for convenience.
parse_formats = _export.parse_formats
save_figure = _export.save_figure
