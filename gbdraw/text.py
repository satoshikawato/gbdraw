#!/usr/bin/env python
# coding: utf-8

import functools
import logging
import xml.etree.ElementTree as ET
from importlib import resources
from typing import Dict, List, Optional, Union

from fontTools.ttLib import TTFont
from svgwrite.text import Text

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
#  Logic ported from find_font_files.py to remove dependency
# ------------------------------------------------------------------
def get_text_bbox_size_pixels(font_path, text, font_size, dpi):
    """
    Directly parses the font file using fontTools to calculate text dimensions.
    This logic was originally in find_font_files.py.
    """
    try:
        font = TTFont(font_path)
    except Exception as e:
        logger.warning(f"Failed to load font file {font_path}: {e}")
        # Fallback approximation
        return len(str(text)) * float(font_size) * 0.6, float(font_size)

    hmtx = font["hmtx"]
    cmap = font["cmap"]
    head = font["head"]
    units_per_em = head.unitsPerEm

    # Get cmap (Format 4 or 12 preferred)
    t = cmap.getcmap(3, 1).cmap

    total_width = 0
    ymaxes = []
    ymins = []

    text_str = str(text)
    if not text_str:
        return 0.0, 0.0

    rsb_previous = 0

    for i, char in enumerate(text_str):
        char_code = ord(char)
        if char_code in t:
            glyph_index = t[char_code]
        else:
            glyph_index = 0

        try:
            # Try to get vertical metrics from glyf table if available
            if "glyf" in font:
                g = font["glyf"][glyph_index]
                ymax = g.yMax if hasattr(g, "yMax") else 0
                ymin = g.yMin if hasattr(g, "yMin") else 0
                xmax = g.xMax if hasattr(g, "xMax") else 0
            else:
                # Fallback for OTF/CFF fonts without glyf table
                ymax = head.yMax
                ymin = head.yMin
                xmax = 0
        except Exception:
            ymax, ymin, xmax = 0, 0, 0

        ymaxes.append(ymax)
        ymins.append(ymin)

        # Calculate horizontal advance
        try:
            advance_width, lsb = hmtx[glyph_index]
        except KeyError:
            advance_width, lsb = 1000, 0

        if xmax == 0:
            advance_width -= rsb_previous

        if lsb > 0:
            total_width -= lsb
        else:
            total_width += lsb

        rsb = xmax - advance_width
        total_width += advance_width

        if rsb > 0:
            total_width -= rsb
        else:
            total_width += rsb

        if i == 0:
            total_width += lsb

        rsb_previous = rsb

    # Formula: pixels = (design_units * font_size * dpi) / (72 * units_per_em)
    scale_factor = (float(font_size) * dpi) / (72 * units_per_em)

    text_width_pixels = total_width * scale_factor

    if ymaxes and ymins:
        max_y = max(ymaxes)
        min_y = min(ymins)
        text_height_pixels = (max_y + abs(min_y)) * scale_factor
    else:
        text_height_pixels = float(font_size)

    return text_width_pixels, text_height_pixels


# ------------------------------------------------------------------
#  Main BBox Calculation Function
# ------------------------------------------------------------------
@functools.lru_cache(maxsize=4096)
def calculate_bbox_dimensions(text, font_family, font_size, dpi):
    """
    Calculates bounding box dimensions using bundled font files in gbdraw package.
    Locates the font in gbdraw.data.fonts using importlib.resources.
    """
    target_font_path = None

    try:
        # Access the bundled fonts directory
        # Assumes structure: gbdraw/data/fonts/*.ttf
        font_dir = resources.files("gbdraw.data").joinpath("fonts")

        # List available TTF/OTF files
        available_fonts = [
            f for f in font_dir.iterdir() if f.name.lower().endswith((".ttf", ".otf"))
        ]

        if not available_fonts:
            # raise FileNotFoundError("No font files found in gbdraw.data.fonts")
            pass  # Fallback to approximation below

        else:
            # Simple matching logic
            requested = font_family.split(",")[0].strip().lower().replace(" ", "")

            for f in available_fonts:
                if requested in f.name.lower():
                    target_font_path = str(f)
                    break

            # Fallback to the first available font
            if not target_font_path:
                target_font_path = str(available_fonts[0])

            # Calculate using fontTools
            return get_text_bbox_size_pixels(target_font_path, text, font_size, dpi)

    except Exception as e:
        logger.debug(f"Font calculation failed ({e}). Using approximation.")

    # Fallback Approximation
    try:
        f_size = float(font_size)
    except Exception:
        f_size = 12.0
    return len(str(text)) * f_size * 0.6, f_size


def create_text_element(
    text: str,
    x: float,
    y: float,
    font_size: str,
    font_weight: str,
    font_family: str,
    text_anchor: str = "middle",
    dominant_baseline: str = "middle",
) -> Text:
    return Text(
        text,
        insert=(x, y),
        stroke="none",
        fill="black",
        font_size=font_size,
        font_weight=font_weight,
        font_family=font_family,
        text_anchor=text_anchor,
        dominant_baseline=dominant_baseline,
    )


def parse_mixed_content_text(input_text: str) -> List[Dict[str, Union[str, bool, None]]]:
    parts: List[Dict[str, Union[str, bool, None]]] = []
    try:
        wrapped_text: str = f"<root>{input_text}</root>"
        root: ET.Element = ET.fromstring(wrapped_text)
        if list(root):
            for element in root:
                if element.tag == "i":
                    parts.append({"text": element.text, "italic": True})
                else:
                    parts.append({"text": element.text, "italic": False})
                if element.tail:
                    parts.append({"text": element.tail, "italic": False})
        else:
            parts.append({"text": root.text, "italic": False})
    except ET.ParseError:
        parts.append({"text": input_text, "italic": False})
    return parts


__all__ = [
    "calculate_bbox_dimensions",
    "create_text_element",
    "get_text_bbox_size_pixels",
    "parse_mixed_content_text",
]


