"""Render/export helpers for the public API layer."""

from __future__ import annotations

import logging
import os
import shutil
import sys
import tempfile
from pathlib import Path
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
from gbdraw.render.output_paths import preflight_output_paths

logger = logging.getLogger(__name__)


def _normalize_formats(formats: Sequence[str] | str | None) -> list[str]:
    return normalize_format_sequence(formats, logger=logger)


def _late_output_collision_message(path: str, *, overwrite: bool) -> str:
    if overwrite:
        return (
            f"Output target appeared during export: {path}. "
            "Refusing to replace it."
        )
    return (
        f"Output file already exists: {path}. "
        "Use overwrite=True to replace it."
    )


def _resolve_base_prefix(canvas: Drawing, output_prefix: str | None, output_dir: str | None) -> str:
    base_from_canvas = os.path.splitext(getattr(canvas, "filename", "") or "out.svg")[0]
    base_prefix = output_prefix or base_from_canvas or "out"

    if output_dir:
        if output_prefix and os.path.dirname(output_prefix) not in {"", "."}:
            raise ValidationError(
                "output_dir cannot be combined with a path-like output_prefix; choose one."
            )
        base_prefix = os.path.join(output_dir, os.path.basename(base_prefix))

    return base_prefix


def _ensure_overwrite_ok(paths: Iterable[str | Path], overwrite: bool) -> None:
    output_paths = preflight_output_paths(paths, overwrite=overwrite)
    if not overwrite:
        return
    for path in output_paths:
        if not path.exists():
            continue
        try:
            path.unlink()
        except OSError as exc:
            raise ExportError(
                f"Could not replace existing output file: {path}"
            ) from exc


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
    CairoSVG. Unlike the deprecated
    :func:`gbdraw.render.export.save_figure`, explicitly requested formats are
    strict: failure to generate any one of them raises an exception.
    """

    fmt_list = _normalize_formats(formats)
    base_prefix = _resolve_base_prefix(canvas, output_prefix, output_dir)
    svg_filename = resolve_format_output_path(base_prefix, SVG_FORMAT)

    output_paths = resolve_output_paths(base_prefix, fmt_list, include_base_svg=True)

    _ensure_overwrite_ok(output_paths, overwrite)
    try:
        Path(base_prefix).parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValidationError(
            f"Could not prepare output directory: {Path(base_prefix).parent}"
        ) from exc

    try:
        canvas.filename = svg_filename
        with open(svg_filename, "x", encoding="utf-8") as svg_file:
            canvas.write(svg_file)
    except FileExistsError as exc:
        raise ValidationError(
            _late_output_collision_message(
                svg_filename,
                overwrite=overwrite,
            )
        ) from exc
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
            interactive_svg = enrich_svg(
                svg_source,
                context=interactive_context,
                result_name=os.path.basename(interactive_filename),
            )
        except GbdrawError:
            raise
        except Exception as exc:
            raise GbdrawError(f"Interactive SVG export failed: {exc}") from exc
        try:
            with open(interactive_filename, "x", encoding="utf-8") as handle:
                handle.write(interactive_svg)
        except FileExistsError as exc:
            raise ValidationError(
                _late_output_collision_message(
                    interactive_filename,
                    overwrite=overwrite,
                )
            ) from exc
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
            with tempfile.SpooledTemporaryFile(
                max_size=8 * 1024 * 1024,
                mode="w+b",
            ) as staged_output:
                if fmt == "png":
                    cairosvg_module.svg2png(
                        bytestring=svg_bytes,
                        write_to=staged_output,
                    )
                elif fmt == "pdf":
                    cairosvg_module.svg2pdf(
                        bytestring=svg_bytes,
                        write_to=staged_output,
                    )
                elif fmt in {"ps", "eps"}:
                    cairosvg_module.svg2ps(
                        bytestring=svg_bytes,
                        write_to=staged_output,
                    )
                staged_output.seek(0)
                with open(out_file, "xb") as output_file:
                    shutil.copyfileobj(staged_output, output_file)
        except FileExistsError as exc:
            raise ValidationError(
                _late_output_collision_message(
                    out_file,
                    overwrite=overwrite,
                )
            ) from exc
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


__all__ = ["preflight_output_paths", "render_to_bytes", "save_figure_to"]
