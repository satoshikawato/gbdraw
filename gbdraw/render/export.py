#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import importlib
import os
import sys
import logging
from types import ModuleType
from typing import List

from svgwrite import Drawing

logger = logging.getLogger(__name__)

_cairosvg_module: ModuleType | None = None


class _LazyCairoSvgAvailability:
    """Bool-like compatibility proxy that resolves CairoSVG only on demand."""

    def __bool__(self) -> bool:
        return has_cairosvg()

    def __repr__(self) -> str:
        return str(has_cairosvg())


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


CAIROSVG_AVAILABLE = _LazyCairoSvgAvailability()


def parse_formats(out_formats: str) -> list[str]:
    list_of_formats: list[str] = [fmt.strip().lower() for fmt in out_formats.split(",")]
    accepted_formats: list[str] = ["svg", "png", "eps", "ps", "pdf"]

    accepted_list_of_formats: list[str] = []

    for fmt in list_of_formats:
        if fmt in accepted_formats:
            accepted_list_of_formats.append(fmt)
        else:
            logger.warning(f"WARNING: Unaccepted/unrecognized output file format: {fmt}")

    # If no valid formats are found, default to 'png'
    if not accepted_list_of_formats:
        logger.warning(
            "WARNING: No valid output file format was specified; generate a PNG file."
        )
        accepted_list_of_formats = ["png"]

    # Remove duplicates if exist
    accepted_list_of_formats = list(set(accepted_list_of_formats))
    return accepted_list_of_formats


def save_figure(canvas: Drawing, list_of_formats: List[str]) -> None:
    """
    Saves the rendered figure.
    Logic:
      1. Always save SVG (base format).
      2. If Pyodide (Web), skip conversion (JS handles it).
      3. If CLI, try CairoSVG. If not available, warn and skip.
    """
    base_filename = os.path.splitext(canvas.filename)[0]
    svg_filename = f"{base_filename}.svg"

    # 1. Always save SVG
    canvas.saveas(svg_filename)
    logger.info(f"Generated SVG: {svg_filename}")

    # 2. Prepare other formats (excluding SVG)
    formats_to_process = [f for f in list_of_formats if f != "svg"]

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
        svg_bytes = canvas.tostring().encode("utf-8")
        for fmt in formats_to_process:
            out_file = f"{base_filename}.{fmt}"

            if fmt == "png":
                cairosvg_module.svg2png(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "pdf":
                cairosvg_module.svg2pdf(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "ps":
                cairosvg_module.svg2ps(bytestring=svg_bytes, write_to=out_file)
            elif fmt == "eps":
                cairosvg_module.svg2ps(bytestring=svg_bytes, write_to=out_file)

            logger.info(f"Generated {fmt.upper()}: {out_file}")
    except Exception as e:
        logger.error(f"Failed to generate images using CairoSVG: {e}")


__all__ = [
    "CAIROSVG_AVAILABLE",
    "get_cairosvg",
    "has_cairosvg",
    "parse_formats",
    "save_figure",
]


