#!/usr/bin/env python
# coding: utf-8

import functools
import logging
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from importlib import resources
from typing import Dict, List, Mapping, Optional, Union

from fontTools.ttLib import TTFont
from svgwrite.text import Text

logger = logging.getLogger(__name__)

# Global caches for font objects to avoid repeated file loading
_font_cache: Dict[str, TTFont] = {}
_font_path_cache: Dict[str, Optional[str]] = {}

_BUNDLED_FONT_FAMILY_ALIASES = {
    "liberationsans": "LiberationSans",
    "arial": "LiberationSans",
    "helvetica": "LiberationSans",
    "nimbussansl": "LiberationSans",
    "sans-serif": "LiberationSans",
    "sansserif": "LiberationSans",
    "liberationserif": "LiberationSerif",
    "timesnewroman": "LiberationSerif",
    "times": "LiberationSerif",
    "serif": "LiberationSerif",
    "liberationmono": "LiberationMono",
    "couriernew": "LiberationMono",
    "courier": "LiberationMono",
    "monospace": "LiberationMono",
}


# ------------------------------------------------------------------
#  Kerning Support (CLI only - Pyodide uses Canvas API)
# ------------------------------------------------------------------
def _get_kerning_value(font, left_glyph_index, right_glyph_index):
    """
    Get kerning adjustment value for a glyph pair from kern or GPOS tables.

    Args:
        font: TTFont object
        left_glyph_index: Glyph index of the left glyph
        right_glyph_index: Glyph index of the right glyph

    Returns:
        int: Kerning adjustment value in font design units, or 0 if not found
    """
    if left_glyph_index is None or right_glyph_index is None:
        return 0

    try:
        # Try kern table first (simpler, faster, used in older fonts)
        if "kern" in font:
            kern_table = font["kern"]
            if hasattr(kern_table, "kernTables"):
                for subtable in kern_table.kernTables:
                    if subtable.coverage & 1:
                        pair = (left_glyph_index, right_glyph_index)
                        if hasattr(subtable, "kernTable") and pair in subtable.kernTable:
                            return subtable.kernTable[pair]
    except Exception as e:
        logger.debug(f"Error reading kern table: {e}")

    try:
        # Try GPOS table (modern OpenType fonts)
        if "GPOS" in font:
            gpos_table = font["GPOS"]
            if hasattr(gpos_table, "table") and gpos_table.table:
                lookup_list = gpos_table.table.LookupList
                if lookup_list:
                    for lookup in lookup_list.Lookup:
                        if lookup.LookupType == 2:
                            for subtable in lookup.SubTable:
                                if hasattr(subtable, "PairSets") and hasattr(subtable, "Coverage"):
                                    coverage = subtable.Coverage
                                    if hasattr(coverage, "glyphs") and left_glyph_index in coverage.glyphs:
                                        try:
                                            left_idx = coverage.glyphs.index(left_glyph_index)
                                            if left_idx < len(subtable.PairSets):
                                                pair_set = subtable.PairSets[left_idx]
                                                if pair_set:
                                                    for pvr in pair_set:
                                                        if hasattr(pvr, "SecondGlyph") and pvr.SecondGlyph == right_glyph_index:
                                                            if hasattr(pvr, "Value1") and hasattr(pvr.Value1, "XAdvance"):
                                                                return pvr.Value1.XAdvance
                                        except (ValueError, IndexError, AttributeError):
                                            pass
                                elif hasattr(subtable, "ClassDef1") and hasattr(subtable, "ClassDef2"):
                                    class_def1 = subtable.ClassDef1
                                    class_def2 = subtable.ClassDef2
                                    class1 = class_def1.get(left_glyph_index, 0) if hasattr(class_def1, "get") else 0
                                    class2 = class_def2.get(right_glyph_index, 0) if hasattr(class_def2, "get") else 0
                                    if hasattr(subtable, "Class1Record") and class1 < len(subtable.Class1Record):
                                        class1_record = subtable.Class1Record[class1]
                                        if hasattr(class1_record, "Class2Record") and class2 < len(class1_record.Class2Record):
                                            class2_record = class1_record.Class2Record[class2]
                                            if hasattr(class2_record, "Value1") and hasattr(class2_record.Value1, "XAdvance"):
                                                return class2_record.Value1.XAdvance
    except Exception as e:
        logger.debug(f"Error reading GPOS table: {e}")

    return 0


# ------------------------------------------------------------------
#  Font Loading with Caching
# ------------------------------------------------------------------
def _get_cached_font(font_path: str) -> Optional[TTFont]:
    """Load a font file with caching to avoid repeated disk I/O."""
    if font_path in _font_cache:
        return _font_cache[font_path]

    try:
        font = TTFont(font_path)
        _font_cache[font_path] = font
        return font
    except Exception as e:
        logger.warning(f"Failed to load font file {font_path}: {e}")
        return None


def _normalize_font_family_candidate(value: str) -> str:
    return value.strip().strip("\"'").strip()


def _comparable_font_name(value: str) -> str:
    return _normalize_font_family_candidate(value).lower().replace(" ", "")


def _normalize_font_weight(font_weight: str | int | float | None) -> str:
    value = "normal" if font_weight is None else str(font_weight).strip().lower()
    if value in {"bold", "bolder"}:
        return "bold"
    try:
        return "bold" if int(float(value)) >= 600 else "normal"
    except (TypeError, ValueError):
        return "normal"


def _normalize_font_style(font_style: str | None) -> str:
    value = "normal" if font_style is None else str(font_style).strip().lower()
    return "italic" if value in {"italic", "oblique"} else "normal"


def _font_style_suffix(font_weight: str | int | float | None, font_style: str | None) -> str:
    weight = _normalize_font_weight(font_weight)
    style = _normalize_font_style(font_style)
    if weight == "bold" and style == "italic":
        return "BoldItalic"
    if weight == "bold":
        return "Bold"
    if style == "italic":
        return "Italic"
    return "Regular"


def _bundled_font_dirs():
    data_dir = resources.files("gbdraw.data")
    return (data_dir.joinpath("fonts"), data_dir)


def _available_bundled_fonts():
    available_fonts = []
    seen = set()
    for font_dir in _bundled_font_dirs():
        try:
            if not font_dir.is_dir():
                continue
            for font_file in font_dir.iterdir():
                if not font_file.name.lower().endswith((".ttf", ".otf")):
                    continue
                font_path = str(font_file)
                if font_path in seen:
                    continue
                seen.add(font_path)
                available_fonts.append(font_file)
        except Exception as e:
            logger.debug(f"Bundled font directory lookup failed for {font_dir!s}: {e}")
    return available_fonts


def _font_family_prefix(family: str) -> str | None:
    comparable_family = _comparable_font_name(family)
    return _BUNDLED_FONT_FAMILY_ALIASES.get(comparable_family)


def _select_bundled_font_variant(
    available_fonts,
    *,
    prefix: str,
    font_weight: str | int | float | None,
    font_style: str | None,
) -> Optional[str]:
    desired_suffix = _font_style_suffix(font_weight, font_style)
    desired_stem = f"{prefix}-{desired_suffix}".lower()
    regular_stem = f"{prefix}-Regular".lower()

    by_stem = {os.path.splitext(font_file.name)[0].lower(): font_file for font_file in available_fonts}
    if desired_stem in by_stem:
        return str(by_stem[desired_stem])
    if regular_stem in by_stem:
        logger.debug(
            "Bundled font variant %s was not found; falling back to %s-Regular.",
            desired_stem,
            prefix,
        )
        return str(by_stem[regular_stem])
    return None


def _resolve_font_path(
    font_family: str,
    *,
    allow_system: bool = False,
    font_weight: str | int | float | None = "normal",
    font_style: str | None = "normal",
) -> Optional[str]:
    """Resolve font family name to a font file path with caching."""
    cache_key = (
        f"{font_family}\0system={int(allow_system)}"
        f"\0weight={_normalize_font_weight(font_weight)}"
        f"\0style={_normalize_font_style(font_style)}"
    )
    if cache_key in _font_path_cache:
        return _font_path_cache[cache_key]

    families = [
        family
        for family in (_normalize_font_family_candidate(part) for part in str(font_family).split(","))
        if family
    ]
    target_font_path = None
    try:
        available_fonts = _available_bundled_fonts()

        if available_fonts:
            for family in families:
                prefix = _font_family_prefix(family)
                if prefix is None:
                    continue
                target_font_path = _select_bundled_font_variant(
                    available_fonts,
                    prefix=prefix,
                    font_weight=font_weight,
                    font_style=font_style,
                )
                if target_font_path:
                    break

            requested_names = [_comparable_font_name(family) for family in families]

            if target_font_path is None:
                for requested in requested_names:
                    for f in available_fonts:
                        if requested in _comparable_font_name(f.name):
                            prefix = os.path.splitext(f.name)[0].split("-", 1)[0]
                            target_font_path = _select_bundled_font_variant(
                                available_fonts,
                                prefix=prefix,
                                font_weight=font_weight,
                                font_style=font_style,
                            ) or str(f)
                            break
                    if target_font_path:
                        break

            if not target_font_path:
                target_font_path = _select_bundled_font_variant(
                    available_fonts,
                    prefix="LiberationSans",
                    font_weight=font_weight,
                    font_style=font_style,
                ) or str(available_fonts[0])
    except Exception as e:
        logger.debug(f"Font path resolution failed ({e}).")

    if allow_system and target_font_path is None:
        fc_match = shutil.which("fc-match")
        if fc_match:
            for family in families:
                if _comparable_font_name(family) in {"sans-serif", "serif", "monospace", "cursive", "fantasy"}:
                    continue
                try:
                    result = subprocess.run(
                        [fc_match, "-f", "%{file}", family],
                        check=False,
                        capture_output=True,
                        text=True,
                        timeout=2,
                    )
                except Exception as e:
                    logger.debug(f"Fontconfig lookup failed for {family!r}: {e}")
                    continue
                candidate = result.stdout.strip()
                if result.returncode == 0 and candidate and os.path.exists(candidate):
                    target_font_path = candidate
                    break

    _font_path_cache[cache_key] = target_font_path
    return target_font_path


# ------------------------------------------------------------------
#  Logic ported from find_font_files.py to remove dependency
# ------------------------------------------------------------------
def get_text_bbox_size_pixels(font_path, text, font_size, dpi, *, svg_units: bool = False):
    """
    Directly parses the font file using fontTools to calculate text dimensions.
    This logic was originally in find_font_files.py.

    The width calculation accounts for kerning adjustments from the font's kern
    table (legacy fonts) or GPOS table (modern OpenType fonts), providing
    accurate bounding box dimensions that match the actual rendered text.

    Args:
        font_path (str): Path to the font file (e.g., .ttf, .otf).
        text (str): The text string for which to calculate the size.
        font_size (int or float): The font size (in points).
        dpi (int): Dots Per Inch (DPI). Used to convert font units into pixels.
    Returns:
        tuple[float, float]:
            text_width_pixels (float): The calculated width of the text in pixels,
                including kerning adjustments.
            text_height_pixels (float): The calculated height of the text in pixels.
    """
    font = _get_cached_font(font_path)
    if font is None:
        # Fallback approximation
        return len(str(text)) * float(font_size) * 0.6, float(font_size)

    hmtx = font["hmtx"]
    cmap = font["cmap"]
    head = font["head"]
    units_per_em = head.unitsPerEm

    # Get cmap (Format 4 or 12 preferred)
    t = cmap.getcmap(3, 1).cmap

    total_width = 0
    total_advance_width = 0
    ymaxes = []
    ymins = []

    text_str = str(text)
    if not text_str:
        return 0.0, 0.0

    rsb_previous = 0
    previous_glyph_index = None

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

        # Apply kerning adjustment between consecutive glyphs
        if previous_glyph_index is not None:
            kerning_value = _get_kerning_value(font, previous_glyph_index, glyph_index)
            total_width += kerning_value
            total_advance_width += kerning_value

        # Calculate horizontal advance
        try:
            advance_width, lsb = hmtx[glyph_index]
        except KeyError:
            advance_width, lsb = 1000, 0

        total_advance_width += advance_width

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
        previous_glyph_index = glyph_index

    if svg_units:
        # SVG font-size values are user units (CSS px in the generated SVG),
        # not typographic points.
        scale_factor = float(font_size) / units_per_em
    else:
        # Legacy layout calculations treat configured font sizes as points.
        scale_factor = (float(font_size) * dpi) / (72 * units_per_em)

    text_width_pixels = (total_advance_width if svg_units else total_width) * scale_factor

    if ymaxes and ymins:
        max_y = max(ymaxes)
        min_y = min(ymins)
        text_height_pixels = (max_y + abs(min_y)) * scale_factor
    else:
        text_height_pixels = float(font_size)

    return text_width_pixels, text_height_pixels


# ------------------------------------------------------------------
#  Browser Canvas Text Measurement (for Pyodide)
# ------------------------------------------------------------------
def _measure_text_canvas(
    text: str,
    font_family: str,
    font_size: float,
    font_weight: str | int | float | None = "normal",
    font_style: str | None = "normal",
) -> tuple:
    """
    Measure text dimensions using browser's Canvas API.
    This is much faster than fontTools in Pyodide because it uses
    the browser's native text rendering engine.

    Args:
        text: The text to measure
        font_family: Font family name (e.g., "Liberation Sans")
        font_size: Font size in pixels

    Returns:
        tuple: (width, height) in pixels
    """
    try:
        from js import document  # type: ignore[import-not-found]

        canvas = document.createElement("canvas")
        ctx = canvas.getContext("2d")

        # Map bundled font names to web-safe equivalents
        # Liberation Sans -> Arial (metrics are nearly identical)
        # Liberation Serif -> Times New Roman
        # Liberation Mono -> Courier New
        web_font = font_family
        font_lower = font_family.lower()
        if "liberation" in font_lower:
            if "sans" in font_lower:
                web_font = "Arial, sans-serif"
            elif "serif" in font_lower:
                web_font = "Times New Roman, serif"
            elif "mono" in font_lower:
                web_font = "Courier New, monospace"

        ctx.font = (
            f"{_normalize_font_style(font_style)} "
            f"{_normalize_font_weight(font_weight)} "
            f"{font_size}px {web_font}"
        )
        metrics = ctx.measureText(str(text))

        # Width from measureText
        width = metrics.width

        # Height estimation: actualBoundingBoxAscent + actualBoundingBoxDescent
        # These are available in modern browsers
        if hasattr(metrics, "actualBoundingBoxAscent") and hasattr(
            metrics, "actualBoundingBoxDescent"
        ):
            height = metrics.actualBoundingBoxAscent + metrics.actualBoundingBoxDescent
        else:
            # Fallback: approximate height as font_size * 1.2
            height = font_size * 1.2

        return float(width), float(height)
    except Exception as e:
        logger.debug(f"Canvas text measurement failed: {e}")
        # Fallback approximation
        return len(str(text)) * float(font_size) * 0.6, float(font_size)


# ------------------------------------------------------------------
#  Main BBox Calculation Function
# ------------------------------------------------------------------
@functools.lru_cache(maxsize=4096)
def calculate_bbox_dimensions(
    text,
    font_family,
    font_size,
    dpi,
    font_weight: str | int | float | None = "normal",
    font_style: str | None = "normal",
):
    """
    Calculates bounding box dimensions using bundled font files in gbdraw package.
    In Pyodide/browser environments, uses the native Canvas API for fast measurement.
    In native Python, uses fontTools for accurate measurement with kerning.
    """
    # Use browser's Canvas API in Pyodide (much faster, uses native text engine)
    if "pyodide" in sys.modules:
        return _measure_text_canvas(
            text,
            font_family,
            float(font_size),
            font_weight=font_weight,
            font_style=font_style,
        )

    # Native Python: Use fontTools with kerning support
    target_font_path = _resolve_font_path(
        font_family,
        font_weight=font_weight,
        font_style=font_style,
    )

    if target_font_path:
        # Calculate using fontTools (font object is cached internally)
        return get_text_bbox_size_pixels(target_font_path, text, font_size, dpi)

    # Fallback Approximation
    try:
        f_size = float(font_size)
    except Exception:
        f_size = 12.0
    return len(str(text)) * f_size * 0.6, f_size


@functools.lru_cache(maxsize=4096)
def calculate_svg_bbox_dimensions(
    text,
    font_family,
    font_size,
    dpi,
    font_weight: str | int | float | None = "normal",
    font_style: str | None = "normal",
):
    """
    Calculates text dimensions in SVG user units for width-sensitive layout.

    This is separate from ``calculate_bbox_dimensions`` because older layout
    reserves rely on its point-style sizing. Linear definition columns need the
    actual rendered SVG width so their right edge can sit close to the record
    axis.
    """
    if "pyodide" in sys.modules:
        return _measure_text_canvas(
            text,
            font_family,
            float(font_size),
            font_weight=font_weight,
            font_style=font_style,
        )

    target_font_path = _resolve_font_path(
        font_family,
        allow_system=True,
        font_weight=font_weight,
        font_style=font_style,
    )
    if target_font_path:
        return get_text_bbox_size_pixels(
            target_font_path,
            text,
            font_size,
            dpi,
            svg_units=True,
        )

    try:
        f_size = float(font_size)
    except Exception:
        f_size = 12.0

    estimated_width = 0.0
    for char in str(text):
        if char.isspace():
            estimated_width += 0.28 * f_size
        elif char in "ilI.,:;|![]()'`":
            estimated_width += 0.28 * f_size
        elif char in "mwMW@#%&":
            estimated_width += 0.82 * f_size
        elif char.isupper():
            estimated_width += 0.62 * f_size
        elif char.isdigit():
            estimated_width += 0.55 * f_size
        else:
            estimated_width += 0.48 * f_size
    return estimated_width, f_size


def create_text_element(
    text: str,
    x: float,
    y: float,
    font_size: str | float,
    font_weight: str,
    font_family: str,
    text_anchor: str = "middle",
    dominant_baseline: str = "middle",
    *,
    stroke: str = "none",
    fill: str = "black",
    extra_attrs: Mapping[str, str] | None = None,
) -> Text:
    text_el = Text(
        text,
        insert=(x, y),
        stroke=stroke,
        fill=fill,
        font_size=font_size,
        font_weight=font_weight,
        font_family=font_family,
        text_anchor=text_anchor,
        dominant_baseline=dominant_baseline,
        debug=False if extra_attrs else True,
    )
    for key, value in (extra_attrs or {}).items():
        text_el.attribs[str(key)] = str(value)
    return text_el


def parse_mixed_content_text(input_text: str) -> List[Dict[str, Union[str, bool, None]]]:
    parts: List[Dict[str, Union[str, bool, None]]] = []
    try:
        wrapped_text: str = f"<root>{input_text}</root>"
        root: ET.Element = ET.fromstring(wrapped_text)
        if list(root):
            if root.text is not None:
                parts.append({"text": root.text, "italic": False})
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
    "calculate_svg_bbox_dimensions",
    "create_text_element",
    "get_text_bbox_size_pixels",
    "parse_mixed_content_text",
]


