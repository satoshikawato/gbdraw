#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import functools
import importlib
import os
import re
import sys
import logging
import warnings
import xml.etree.ElementTree as ET
from types import ModuleType
from typing import Any, List

from svgwrite import Drawing

from gbdraw.core.text import dominant_baseline_shift_em, get_font_vertical_metrics
from gbdraw.exceptions import GbdrawError
from gbdraw.render.formats import (
    CAIROSVG_FORMATS,
    INTERACTIVE_SVG_FORMAT,
    SVG_FORMAT,
    classify_formats,
    parse_format_string,
    resolve_format_output_path,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

logger = logging.getLogger(__name__)

_cairosvg_module: ModuleType | None = None


def _load_cairosvg() -> ModuleType | None:
    global _cairosvg_module
    if _cairosvg_module is not None:
        return _cairosvg_module
    try:
        _cairosvg_module = importlib.import_module("cairosvg")
    except (ImportError, OSError):
        return None
    return _cairosvg_module


def has_cairosvg() -> bool:
    return _load_cairosvg() is not None


def get_cairosvg() -> ModuleType:
    cairosvg_module = _load_cairosvg()
    if cairosvg_module is None:
        raise ImportError("CairoSVG is not installed. Install with: pip install gbdraw[export]")
    return cairosvg_module


_SVG_NAMESPACE = "http://www.w3.org/2000/svg"
_CAIROSVG_CONVERTERS = {"png": "svg2png", "pdf": "svg2pdf", "ps": "svg2ps", "eps": "svg2ps"}
# CairoSVG 2.x ignores dominant-baseline on textPath and places these values on
# plain text with other metrics than browsers (hanging and middle use the
# text-before-edge and central formulas; the others are ignored).
_CAIROSVG_MISPLACED_TEXT_BASELINES = frozenset({"hanging", "middle", "mathematical", "ideographic"})
_TEXT_POSITION_ATTRIBUTES = ("y", "dy", "dominant-baseline")
_NUMBER_RE = re.compile(r"^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*(px|pt)?\s*$")


def _local_name(tag: object) -> str:
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def _font_size_px(value: str | None, inherited: float | None) -> float | None:
    if value is None:
        return inherited
    match = _NUMBER_RE.match(value)
    if match is None:
        return None
    size = float(match.group(1))
    return size * 96.0 / 72.0 if match.group(2) == "pt" else size


def _shifted_dy(dy: str | None, shift: float) -> str | None:
    if dy is None or not dy.strip():
        return f"{shift:.4f}"
    first, _, rest = dy.strip().replace(",", " ").partition(" ")
    match = _NUMBER_RE.match(first)
    if match is None or match.group(2) == "pt":
        return None
    return f"{float(match.group(1)) + shift:.4f} {rest}".strip()


def _apply_baseline_shift(element: ET.Element, baseline: str, style: dict[str, Any]) -> bool:
    if style["font_size"] is None:
        return False
    metrics = get_font_vertical_metrics(style["font_family"], style["font_weight"], style["font_style"])
    shift_em = dominant_baseline_shift_em(baseline, metrics) if metrics is not None else None
    if not shift_em:
        return False
    dy = _shifted_dy(element.get("dy"), shift_em * style["font_size"])
    if dy is None:
        return False
    element.set("dy", dy)
    element.set("dominant-baseline", "alphabetic")
    return True


def _rewrite_text_baselines(
    element: ET.Element,
    inherited: dict[str, Any],
    *,
    in_text: bool = False,
) -> bool:
    if "dominant-baseline" in (element.get("style") or ""):
        return False
    style = dict(inherited)
    for attribute in ("dominant-baseline", "font-family", "font-weight", "font-style"):
        if element.get(attribute) is not None:
            style[attribute.replace("-", "_")] = element.get(attribute)
    style["font_size"] = _font_size_px(element.get("font-size"), style["font_size"])
    tag = _local_name(element.tag)
    baseline = str(style["dominant_baseline"] or "auto").strip().lower()
    if tag == "textPath":
        return _apply_baseline_shift(element, baseline, style)
    if tag == "text" and not in_text and baseline in _CAIROSVG_MISPLACED_TEXT_BASELINES:
        descendants = list(element.iter())[1:]
        if not any(
            _local_name(child.tag) == "textPath" or any(child.get(name) is not None for name in _TEXT_POSITION_ATTRIBUTES)
            for child in descendants
        ):
            return _apply_baseline_shift(element, baseline, style)
    changed = False
    for child in element:
        changed = _rewrite_text_baselines(child, style, in_text=in_text or tag == "text") or changed
    return changed


@functools.lru_cache(maxsize=1)
def prepare_svg_for_cairosvg(svg_bytes: bytes) -> bytes:
    """Return SVG bytes whose text CairoSVG places where browsers draw it.

    CairoSVG ignores ``dominant-baseline`` on ``textPath`` (circular tick
    labels) and uses other metrics than browsers for some values on plain
    text. Those elements get the browser baseline offset as ``dy``, computed
    from the packaged font metrics, and ``dominant-baseline="alphabetic"``.
    Values CairoSVG already places like browsers are left unchanged, and the
    input is returned unchanged when nothing needs rewriting.
    """
    if b"dominant-baseline" not in svg_bytes:
        return svg_bytes
    try:
        root = ET.fromstring(svg_bytes)
    except ET.ParseError:
        return svg_bytes
    inherited: dict[str, Any] = {
        "dominant_baseline": "auto",
        "font_family": "sans-serif",
        "font_weight": "normal",
        "font_style": "normal",
        "font_size": 16.0,
    }
    if not _rewrite_text_baselines(root, inherited):
        return svg_bytes
    try:
        return ET.tostring(root, encoding="utf-8", xml_declaration=True, default_namespace=_SVG_NAMESPACE)
    except ValueError:
        return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def convert_svg_with_cairosvg(
    svg_source: bytes | str,
    fmt: str,
    *,
    write_to: Any = None,
    cairosvg_module: ModuleType | None = None,
    **options: Any,
) -> Any:
    """Convert SVG to PNG, PDF, PS, or EPS with CairoSVG.

    Every CairoSVG conversion goes through this function, so raster and
    vector exports get :func:`prepare_svg_for_cairosvg`. ``options`` are passed
    to CairoSVG (for example ``url``, ``output_width``, ``background_color``).
    Returns the converted bytes when ``write_to`` is not given.
    """
    converter_name = _CAIROSVG_CONVERTERS.get(str(fmt).lower())
    if converter_name is None:
        raise ValueError(f"CairoSVG cannot convert to format: {fmt}")
    module = cairosvg_module if cairosvg_module is not None else get_cairosvg()
    source = svg_source.encode("utf-8") if isinstance(svg_source, str) else bytes(svg_source)
    kwargs = dict(options)
    kwargs["bytestring"] = prepare_svg_for_cairosvg(source)
    if write_to is not None:
        kwargs["write_to"] = write_to
    return getattr(module, converter_name)(**kwargs)


def parse_formats(out_formats: str) -> list[str]:
    return parse_format_string(out_formats, logger=logger)


def save_figure(
    canvas: Drawing,
    list_of_formats: List[str],
    *,
    interactive_context: InteractiveSvgContext | None = None,
) -> None:
    """Save a figure using the deprecated warning-and-skip export contract.

    Use :func:`gbdraw.api.save_figure_to` for strict file export or
    :func:`gbdraw.api.render_to_bytes` for in-memory output.
    """
    warnings.warn(
        "gbdraw.render.export.save_figure() is deprecated and will be removed "
        "in gbdraw 0.16; use gbdraw.api.save_figure_to() or "
        "gbdraw.api.render_to_bytes().",
        DeprecationWarning,
        stacklevel=2,
    )
    base_filename = os.path.splitext(canvas.filename)[0]
    svg_filename = resolve_format_output_path(base_filename, SVG_FORMAT)

    # 1. Always save SVG
    canvas.saveas(svg_filename)
    logger.info(f"Generated SVG: {svg_filename}")

    classification = classify_formats(list_of_formats)
    svg_source = canvas.tostring()

    if classification.interactive:
        interactive_filename = resolve_format_output_path(
            base_filename,
            INTERACTIVE_SVG_FORMAT,
        )
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
        with open(interactive_filename, "w", encoding="utf-8") as handle:
            handle.write(interactive_svg)
        logger.info(f"Generated interactive SVG: {interactive_filename}")

    formats_to_process = list(classification.cairosvg)
    if not formats_to_process:
        return

    # --- WebAssembly (Pyodide) Check ---
    if "pyodide" in sys.modules:
        if formats_to_process:
            logger.info(
                "Running in WebAssembly: Image conversion will be handled by the browser."
            )
        return

    # --- CLI Conversion Logic (CairoSVG Only) ---
    try:
        cairosvg_module = get_cairosvg()
    except ImportError:
        # CairoSVG not available; warn user about skipped formats
        missing_formats = ", ".join([f.upper() for f in formats_to_process])
        logger.warning(
            f"⚠️  Skipping generation of: {missing_formats}\n"
            f"   CairoSVG is not installed.\n"
            f"   To enable PNG/PDF support, run: pip install gbdraw[export]\n"
        )
        return

    try:
        svg_bytes = svg_source.encode("utf-8")
        for fmt in formats_to_process:
            if fmt not in CAIROSVG_FORMATS:
                continue
            out_file = resolve_format_output_path(base_filename, fmt)

            convert_svg_with_cairosvg(
                svg_bytes,
                fmt,
                write_to=out_file,
                cairosvg_module=cairosvg_module,
            )

            logger.info(f"Generated {fmt.upper()}: {out_file}")
    except Exception as e:
        logger.error(f"Failed to generate images using CairoSVG: {e}")


__all__ = [
    "convert_svg_with_cairosvg",
    "get_cairosvg",
    "has_cairosvg",
    "parse_formats",
    "prepare_svg_for_cairosvg",
    "save_figure",
]


