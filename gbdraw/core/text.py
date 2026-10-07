#!/usr/bin/env python
# coding: utf-8

import functools
import logging
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
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
#  Kerning Support
# ------------------------------------------------------------------
# Pair kerning follows HarfBuzz, which Chromium and Firefox shape text with
# (hb-ot-shape.cc, ``hb_ot_shape_planner_t::compile``): when the font's GPOS
# table has a ``kern`` feature, its pair-positioning lookups give the kerning
# and the legacy ``kern`` table is ignored; otherwise the first horizontal
# ``kern`` subtable that lists the pair gives it. The pairs are looked up in
# string order, without skipping marks; scripts are not told apart, so the
# lookups of every ``kern`` feature apply (their coverage is per script in
# the bundled fonts).


@dataclass(frozen=True)
class _GlyphPairSubtable:
    """GPOS PairPos format 1: explicit first/second glyph pairs."""

    values: Mapping[tuple[str, str], int]

    def match(self, left: str, right: str) -> Optional[int]:
        return self.values.get((left, right))


@dataclass(frozen=True)
class _ClassPairSubtable:
    """GPOS PairPos format 2: first-glyph coverage and two class definitions."""

    first_glyphs: frozenset[str]
    first_classes: Mapping[str, int]
    second_classes: Mapping[str, int]
    values: tuple[tuple[int, ...], ...]

    def match(self, left: str, right: str) -> Optional[int]:
        if left not in self.first_glyphs:
            return None
        first_class = self.first_classes.get(left, 0)
        second_class = self.second_classes.get(right, 0)
        if first_class >= len(self.values) or second_class >= len(self.values[first_class]):
            return None
        return self.values[first_class][second_class]


_PairSubtable = Union[_GlyphPairSubtable, _ClassPairSubtable]


@dataclass(frozen=True)
class _FontPairKerning:
    """Pair kerning of one font, in font design units."""

    gpos_lookups: Optional[tuple[tuple[_PairSubtable, ...], ...]]
    legacy_pairs: Mapping[tuple[str, str], int]

    def value(self, left: str, right: str) -> int:
        if self.gpos_lookups is None:
            return self.legacy_pairs.get((left, right), 0)
        total = 0
        for subtables in self.gpos_lookups:
            # Within a lookup the first subtable that applies to the pair wins;
            # every lookup of the feature applies.
            for subtable in subtables:
                adjustment = subtable.match(left, right)
                if adjustment is not None:
                    total += adjustment
                    break
        return total


_pair_kerning_cache: Dict[int, tuple[TTFont, _FontPairKerning]] = {}


def _x_advance(value_record) -> int:
    return int(getattr(value_record, "XAdvance", 0) or 0) if value_record is not None else 0


def _pair_subtable(subtable) -> Optional[_PairSubtable]:
    if getattr(subtable, "Format", None) == 1:
        values: Dict[tuple[str, str], int] = {}
        for first, pair_set in zip(subtable.Coverage.glyphs, subtable.PairSet):
            for record in pair_set.PairValueRecord:
                values.setdefault((first, record.SecondGlyph), _x_advance(record.Value1))
        return _GlyphPairSubtable(values)
    if getattr(subtable, "Format", None) == 2:
        return _ClassPairSubtable(
            first_glyphs=frozenset(subtable.Coverage.glyphs),
            first_classes=dict(subtable.ClassDef1.classDefs) if subtable.ClassDef1 else {},
            second_classes=dict(subtable.ClassDef2.classDefs) if subtable.ClassDef2 else {},
            values=tuple(
                tuple(_x_advance(class2.Value1) for class2 in class1.Class2Record)
                for class1 in subtable.Class1Record
            ),
        )
    return None


def _gpos_kern_lookups(font) -> Optional[tuple[tuple[_PairSubtable, ...], ...]]:
    """Return the pair lookups of GPOS ``kern`` features, or None without one."""
    if "GPOS" not in font:
        return None
    table = font["GPOS"].table
    if table is None or table.FeatureList is None or table.LookupList is None:
        return None
    lookup_indices = sorted(
        {
            int(index)
            for record in table.FeatureList.FeatureRecord
            if record.FeatureTag == "kern"
            for index in record.Feature.LookupListIndex
        }
    )
    if not lookup_indices:
        return None
    lookups: List[tuple[_PairSubtable, ...]] = []
    for index in lookup_indices:
        lookup = table.LookupList.Lookup[index]
        raw_subtables = list(lookup.SubTable)
        if lookup.LookupType == 9:
            raw_subtables = [
                subtable.ExtSubTable
                for subtable in raw_subtables
                if subtable.ExtensionLookupType == 2
            ]
        elif lookup.LookupType != 2:
            continue
        subtables = tuple(
            parsed
            for parsed in (_pair_subtable(subtable) for subtable in raw_subtables)
            if parsed is not None
        )
        if subtables:
            lookups.append(subtables)
    return tuple(lookups)


def _legacy_kern_pairs(font) -> Dict[tuple[str, str], int]:
    pairs: Dict[tuple[str, str], int] = {}
    if "kern" not in font:
        return pairs
    for subtable in getattr(font["kern"], "kernTables", ()):
        if subtable.coverage & 1 and hasattr(subtable, "kernTable"):
            for pair, value in subtable.kernTable.items():
                pairs.setdefault(pair, int(value))
    return pairs


def _font_pair_kerning(font) -> _FontPairKerning:
    cached = _pair_kerning_cache.get(id(font))
    if cached is not None and cached[0] is font:
        return cached[1]
    try:
        gpos_lookups = _gpos_kern_lookups(font)
    except Exception as e:
        logger.debug(f"Error reading GPOS table: {e}")
        gpos_lookups = None
    legacy_pairs: Dict[tuple[str, str], int] = {}
    if gpos_lookups is None:
        try:
            legacy_pairs = _legacy_kern_pairs(font)
        except Exception as e:
            logger.debug(f"Error reading kern table: {e}")
    kerning = _FontPairKerning(gpos_lookups=gpos_lookups, legacy_pairs=legacy_pairs)
    _pair_kerning_cache[id(font)] = (font, kerning)
    return kerning


def _get_kerning_value(font, left_glyph, right_glyph) -> int:
    """
    Get the kerning adjustment for a glyph pair from GPOS or the kern table.

    Args:
        font: TTFont object
        left_glyph: Glyph name of the left glyph (as the cmap returns it)
        right_glyph: Glyph name of the right glyph

    Returns:
        int: Kerning adjustment value in font design units, or 0 if not found
    """
    if left_glyph is None or right_glyph is None:
        return 0
    return _font_pair_kerning(font).value(left_glyph, right_glyph)


def _class_pair_candidates(subtable: _ClassPairSubtable, glyph_order) -> set[tuple[str, str]]:
    seconds_by_class: Dict[int, List[str]] = {}
    for glyph in glyph_order:
        seconds_by_class.setdefault(subtable.second_classes.get(glyph, 0), []).append(glyph)
    pairs: set[tuple[str, str]] = set()
    for first in subtable.first_glyphs:
        first_class = subtable.first_classes.get(first, 0)
        row = subtable.values[first_class] if first_class < len(subtable.values) else ()
        for second_class, adjustment in enumerate(row):
            if adjustment:
                pairs.update((first, second) for second in seconds_by_class.get(second_class, ()))
    return pairs


def font_pair_kerning_table(font) -> Dict[tuple[str, str], int]:
    """Return every glyph pair that ``get_text_bbox_size_pixels`` kerns in ``font``.

    Keys are glyph-name pairs in string order; values are the non-zero
    adjustments in font design units, sorted by key.
    """
    kerning = _font_pair_kerning(font)
    candidates: set[tuple[str, str]] = set(kerning.legacy_pairs)
    for subtables in kerning.gpos_lookups or ():
        for subtable in subtables:
            if isinstance(subtable, _GlyphPairSubtable):
                candidates.update(subtable.values)
            else:
                candidates.update(_class_pair_candidates(subtable, font.getGlyphOrder()))
    table: Dict[tuple[str, str], int] = {}
    for left, right in sorted(candidates):
        adjustment = kerning.value(left, right)
        if adjustment:
            table[(left, right)] = adjustment
    return table


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

    The width calculation accounts for pair kerning: the GPOS ``kern`` feature
    when the font has one, else the legacy ``kern`` table (see
    ``_get_kerning_value``), as browsers shaping with HarfBuzz apply it.

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
    Uses bundled fonts for reproducible measurement across native and browser runtimes.
    """
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


@dataclass(frozen=True)
class FontVerticalMetrics:
    """Vertical font metrics as fractions of the font size (em units)."""

    ascent: float
    descent: float
    x_height: float


@functools.lru_cache(maxsize=64)
def get_font_vertical_metrics(
    font_family: str,
    font_weight: str | int | float | None = "normal",
    font_style: str | None = "normal",
) -> Optional[FontVerticalMetrics]:
    """Return the ascent, descent, and x-height of the font used for layout.

    The font is resolved like :func:`calculate_bbox_dimensions`. Ascent and
    descent follow FreeType, which browsers on Linux use: the OS/2 typo
    metrics when USE_TYPO_METRICS is set, otherwise the hhea metrics.
    """
    font_path = _resolve_font_path(
        font_family,
        font_weight=font_weight,
        font_style=font_style,
    )
    font = _get_cached_font(font_path) if font_path else None
    if font is None:
        return None
    units_per_em = float(font["head"].unitsPerEm)
    os2 = font["OS/2"] if "OS/2" in font else None
    hhea = font["hhea"] if "hhea" in font else None
    if os2 is not None and os2.fsSelection & (1 << 7):
        ascent, descent = os2.sTypoAscender, -os2.sTypoDescender
    elif hhea is not None and (hhea.ascent or hhea.descent):
        ascent, descent = hhea.ascent, -hhea.descent
    elif os2 is not None:
        ascent, descent = os2.sTypoAscender, -os2.sTypoDescender
    else:
        ascent, descent = font["head"].yMax, -font["head"].yMin
    x_height = float(getattr(os2, "sxHeight", 0) or 0)
    if x_height <= 0 and "glyf" in font:
        glyph_name = font.getBestCmap().get(ord("x"))
        if glyph_name is not None:
            x_height = float(getattr(font["glyf"][glyph_name], "yMax", 0) or 0)
    if x_height <= 0:
        x_height = ascent / 2.0
    return FontVerticalMetrics(
        ascent=ascent / units_per_em,
        descent=descent / units_per_em,
        x_height=x_height / units_per_em,
    )


def dominant_baseline_shift_em(dominant_baseline: str, metrics: FontVerticalMetrics) -> Optional[float]:
    """Return how far browsers move horizontal glyphs down for ``dominant-baseline``.

    The result is in em units, relative to the alphabetic baseline; positive
    values move glyphs down (toward the glyph bottom). Values follow Chromium
    for fonts without a BASE table. Unknown values return ``None``.
    """
    value = str(dominant_baseline or "").strip().lower()
    if value in {"auto", "alphabetic"}:
        return 0.0
    if value == "middle":
        return metrics.x_height / 2.0
    if value == "central":
        return (metrics.ascent - metrics.descent) / 2.0
    if value == "hanging":
        return 0.8 * metrics.ascent
    if value == "mathematical":
        return metrics.ascent / 2.0
    if value == "text-before-edge":
        return metrics.ascent
    if value in {"text-after-edge", "ideographic"}:
        return -metrics.descent
    return None


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
    "FontVerticalMetrics",
    "calculate_bbox_dimensions",
    "calculate_svg_bbox_dimensions",
    "create_text_element",
    "dominant_baseline_shift_em",
    "font_pair_kerning_table",
    "get_font_vertical_metrics",
    "get_text_bbox_size_pixels",
    "parse_mixed_content_text",
]


