"""Render/export helpers for the public API layer."""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Sequence

from svgwrite import Drawing

from gbdraw.exceptions import ExportError, GbdrawError, ValidationError
from gbdraw.render import export as _export
from gbdraw.render.document import write_svg_document, write_text
from gbdraw.render.formats import (
    INTERACTIVE_SVG_FORMAT,
    SVG_FORMAT,
    classify_formats,
    parse_output_formats,
    resolve_format_output_path,
    resolve_output_paths,
)
from gbdraw.render.font_embedding import FontEmbedding, embed_font_css, svg_font_embedding
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg
from gbdraw.render.output_paths import (
    commit_staged_outputs,
    preflight_output_paths,
    staged_output_files,
)

logger = logging.getLogger(__name__)


def _font_embedding(svg_source: str, embed_fonts: bool) -> FontEmbedding:
    """The fonts the SVG outputs of one drawing embed (D-05), computed once."""
    if not embed_fonts:
        return FontEmbedding(css="", refused_families=())
    try:
        embedding = svg_font_embedding(svg_source)
    except Exception as exc:
        raise ExportError(
            f"Could not embed the fonts in the SVG: {exc}. Turn font embedding off to write it without them "
            "(--no_font_embedding on the command line, embed_fonts=False in Python)."
        ) from exc
    if embedding.refused_families:
        logger.warning(
            "Font %s is not embedded in the SVG because its license (OS/2 fsType) does not allow "
            "embedding; the SVG shows it only where it is installed.",
            ", ".join(repr(family) for family in embedding.refused_families),
        )
    return embedding


def _normalize_formats(formats: Sequence[str] | str | None) -> list[str]:
    return [SVG_FORMAT] if formats is None else list(parse_output_formats(formats))


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


def _cairosvg_module(formats: Sequence[str]) -> ModuleType | None:
    """The CairoSVG module for ``formats``, checked before anything is written."""
    if not formats:
        return None
    if "pyodide" in sys.modules:
        raise ValidationError(
            "Binary file export is not available under WebAssembly (pyodide); "
            "browser-side conversion does not create local output paths."
        )
    try:
        return _export.get_cairosvg()
    except ImportError as exc:
        missing_formats = ", ".join([f.upper() for f in formats])
        raise ValidationError(
            f"Cannot generate {missing_formats}: CairoSVG is not installed. "
            "Install with: pip install gbdraw[export]"
        ) from exc


def save_figure_to(
    canvas: Drawing,
    formats: Sequence[str] | str | None,
    *,
    output_dir: str | None = None,
    output_prefix: str | None = None,
    overwrite: bool = False,
    interactive_context: InteractiveSvgContext | None = None,
    embed_fonts: bool = True,
) -> list[str]:
    """Save a figure to an explicit output directory/prefix.

    This always writes an SVG, then optionally converts to other formats using
    CairoSVG. Unlike the deprecated
    :func:`gbdraw.render.export.save_figure`, explicitly requested formats are
    strict: failure to generate any one of them raises an exception.
    Every format is generated into staged files beside its output before any
    output is written, so a failure changes no existing file and leaves no new
    one, with or without ``overwrite``.
    With ``embed_fonts``, the SVG and interactive SVG carry a subset of each
    font they use that gbdraw does not bundle. PNG, PDF, EPS and PS draw with
    the registered and installed fonts whatever ``embed_fonts``; PDF, EPS and
    PS draw a font whose license forbids embedding it in the nearest
    Liberation font (``gbdraw.render.font_embedding``).
    """

    fmt_list = _normalize_formats(formats)
    base_prefix = _resolve_base_prefix(canvas, output_prefix, output_dir)
    svg_filename = resolve_format_output_path(base_prefix, SVG_FORMAT)
    output_paths = resolve_output_paths(base_prefix, fmt_list, include_base_svg=True)
    preflight_output_paths(output_paths, overwrite=overwrite)
    classification = classify_formats(fmt_list)
    cairosvg_module = _cairosvg_module(classification.cairosvg)
    try:
        Path(base_prefix).parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValidationError(
            f"Could not prepare output directory: {Path(base_prefix).parent}"
        ) from exc

    labels = {svg_filename: "SVG"}
    with staged_output_files(output_paths) as staged:
        try:
            canvas.filename = svg_filename
            document = write_svg_document(canvas)
            embedding = _font_embedding(document.body, embed_fonts)
            with open(staged[Path(svg_filename)], "w", encoding="utf-8") as svg_file:
                svg_file.write(document.prolog)
                write_text(svg_file, embed_font_css(document.body, embedding))
        except GbdrawError:
            raise
        except Exception as exc:
            raise ExportError(f"Failed to generate SVG: {exc}") from exc

        if classification.interactive:
            interactive_filename = resolve_format_output_path(base_prefix, INTERACTIVE_SVG_FORMAT)
            labels[interactive_filename] = "interactive SVG"
            try:
                interactive_svg = embed_font_css(
                    enrich_svg(
                        document.body,
                        context=interactive_context,
                        result_name=os.path.basename(interactive_filename),
                    ),
                    embedding,
                )
            except GbdrawError:
                raise
            except Exception as exc:
                raise GbdrawError(f"Interactive SVG export failed: {exc}") from exc
            try:
                with open(staged[Path(interactive_filename)], "w", encoding="utf-8") as handle:
                    write_text(handle, interactive_svg)
            except OSError as exc:
                raise ExportError(f"Failed to write interactive SVG: {exc}") from exc
            # Written: CairoSVG below does not hold it beside its own tree.
            del interactive_svg

        if cairosvg_module is not None:
            svg_bytes = document.body.encode("utf-8")
            # Nothing reads the text again: CairoSVG holds the bytes and its own tree.
            del document
            for fmt in classification.cairosvg:
                out_file = resolve_format_output_path(base_prefix, fmt)
                labels[out_file] = fmt.upper()
                try:
                    with open(staged[Path(out_file)], "wb") as output_file:
                        _export.convert_svg_with_cairosvg(
                            svg_bytes,
                            fmt,
                            write_to=output_file,
                            cairosvg_module=cairosvg_module,
                        )
                except Exception as exc:
                    raise ExportError(f"Failed to generate {fmt.upper()}: {exc}") from exc

        commit_staged_outputs(staged, overwrite=overwrite)
    for path in output_paths:
        logger.info("Generated %s: %s", labels[path], path)
    return output_paths


def render_to_bytes(
    canvas: Drawing,
    fmt: str,
    *,
    interactive_context: InteractiveSvgContext | None = None,
    embed_fonts: bool = True,
) -> bytes:
    """Render a canvas to bytes (SVG always; PNG/PDF/PS/EPS require CairoSVG).

    ``embed_fonts`` as in :func:`save_figure_to`: SVG and interactive SVG carry
    the custom fonts they use.
    """

    fmt_norm = parse_output_formats((fmt,))[0]
    if fmt_norm in {SVG_FORMAT, INTERACTIVE_SVG_FORMAT}:
        svg_source = write_svg_document(canvas).body
        embedding = _font_embedding(svg_source, embed_fonts)
        if fmt_norm == INTERACTIVE_SVG_FORMAT:
            svg_source = enrich_svg(svg_source, context=interactive_context)
        return embed_font_css(svg_source, embedding).encode("utf-8")

    if "pyodide" in sys.modules:
        raise ValidationError("Binary export is not available under WebAssembly (pyodide).")

    try:
        cairosvg_module = _export.get_cairosvg()
    except ImportError as exc:
        raise ValidationError(
            "CairoSVG is not installed. Install with: pip install gbdraw[export]"
        ) from exc

    svg_string = write_svg_document(canvas).body.encode("utf-8")
    try:
        rendered = _export.convert_svg_with_cairosvg(
            svg_string,
            fmt_norm,
            cairosvg_module=cairosvg_module,
        )
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
