#!/usr/bin/env python
# coding: utf-8

import functools
import logging
import math
import re
import weakref
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from collections.abc import Iterable, Sequence
from typing import Dict, List, Literal, Mapping, Optional, TypedDict, TypeVar, Union

from fontTools.pens.basePen import BasePen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.ttLib import TTFont

from gbdraw.core.fonts import FontFace, FontWeight, load_font, resolve_font_face
from gbdraw.core.numeric import compensated_sum
from gbdraw.svg.elements import Text

logger = logging.getLogger(__name__)


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


# Keyed by the font object and held weakly: an entry lives as long as its
# font, so a face ``gbdraw.core.fonts`` forgets (unregistered or replaced)
# takes its kerning with it, and a new font is a new key even at the address
# of a freed one. The kerning holds no reference back to its font.
_pair_kerning_cache: "weakref.WeakKeyDictionary[TTFont, _FontPairKerning]" = weakref.WeakKeyDictionary()


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
    cached = _pair_kerning_cache.get(font)
    if cached is not None:
        return cached
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
    _pair_kerning_cache[font] = kerning
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
    """Return every glyph pair that ``_face_text_box`` kerns in ``font``.

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
#  Logic ported from find_font_files.py to remove dependency
# ------------------------------------------------------------------
_SVG_SPACE_RUN = re.compile(" +")


def collapse_svg_white_space(text: str) -> str:
    """Collapse white space as an SVG renderer draws ``text``.

    gbdraw writes text without ``xml:space``, so renderers turn tabs and line
    breaks into spaces and runs of spaces into one. Leading and trailing
    spaces are kept: callers also measure the styled fragments (``tspan``) of
    one text, and a renderer drops spaces only at the ends of the whole text.
    """
    spaced = text.replace("\r\n", " ").translate({0x09: 0x20, 0x0A: 0x20, 0x0D: 0x20})
    return _SVG_SPACE_RUN.sub(" ", spaced)


def _glyph_bounds(font: TTFont, glyph_name: str) -> tuple[float, float, float]:
    """``(yMax, yMin, xMax)`` of a glyph; from the outline for CFF fonts."""
    if "glyf" in font:
        glyph = font["glyf"][glyph_name]
        return (
            getattr(glyph, "yMax", 0),
            getattr(glyph, "yMin", 0),
            getattr(glyph, "xMax", 0),
        )
    glyph_set = font.getGlyphSet()
    pen = BoundsPen(glyph_set)
    glyph_set[glyph_name].draw(pen)
    if pen.bounds is None:
        return 0, 0, 0
    _x_min, y_min, x_max, y_max = pen.bounds
    return y_max, y_min, x_max


def _measured_cmap(font: TTFont) -> Mapping[int, str]:
    """The Windows BMP cmap; the best available one for a font without it."""
    bmp_cmap = font["cmap"].getcmap(3, 1)
    return bmp_cmap.cmap if bmp_cmap is not None else (font.getBestCmap() or {})


class FontMetricsTable(TypedDict):
    """The facts text measurement reads from one face, in font units.

    The Web's port of the Legend layout (gbdraw/web/js/services/legend-layout.js)
    measures with them: ``cmap`` is [code point, glyph index, ...] in ascending
    code points, ``glyphs`` is [advance, lsb, xMax, yMax, yMin, ...] by glyph
    index, ``kerning`` is [left glyph, right glyph, value, ...] for the pairs
    ``font_pair_kerning_table`` applies between glyphs a character reaches, and
    ``lineAscent``/``lineDescent`` are the line box ``drawn_text_box`` measures
    heights with (``_line_box_units``).
    """

    unitsPerEm: int
    lineAscent: int
    lineDescent: int
    cmap: list[int]
    glyphs: list[float]
    kerning: list[int]


def _font_unit(value: float) -> float:
    return value if isinstance(value, float) and not value.is_integer() else int(value)


def _glyph_facts(font: TTFont, glyph_name: str) -> list[float]:
    # The expressions of _face_text_box, read once per glyph.
    try:
        ymax, ymin, xmax = _glyph_bounds(font, glyph_name)
    except Exception:
        ymax, ymin, xmax = 0, 0, 0
    try:
        advance, lsb = font["hmtx"][glyph_name]
    except KeyError:
        advance, lsb = font["head"].unitsPerEm, 0
    return [_font_unit(value) for value in (advance, lsb, xmax, ymax, ymin)]


def font_metrics_table(face: FontFace) -> FontMetricsTable:
    """The measurement facts of one face (``FontMetricsTable``).

    tools/generate_legend_font_metrics.py writes them for the bundled faces;
    the Web builds them for a registered font (``FontFace.measurement_name``
    names the face, as a Result's ``legendReflow.fontFile`` does).
    """
    font = load_font(face)
    if font is None:
        raise ValueError(f"Cannot read font file {face.path}")
    order = font.getGlyphOrder()
    index = {glyph: position for position, glyph in enumerate(order)}
    cmap = _measured_cmap(font)
    mapped = set(cmap.values())
    kerning: list[int] = []
    for (left, right), value in font_pair_kerning_table(font).items():
        if value and left in mapped and right in mapped:
            kerning.extend([index[left], index[right], int(value)])
    line_ascent, line_descent = _line_box_units(font)
    return {
        "unitsPerEm": int(font["head"].unitsPerEm),
        "lineAscent": int(line_ascent),
        "lineDescent": int(line_descent),
        "cmap": [value for code_point in sorted(cmap) for value in (code_point, index[cmap[code_point]])],
        "glyphs": [value for glyph in order for value in _glyph_facts(font, glyph)],
        "kerning": kerning,
    }


@functools.lru_cache(maxsize=8192)
def _face_text_box(face: FontFace, text, font_size):
    """The text box of ``text`` in one face, in SVG user units (px at ``font_size``).

    The width is the advance with pair kerning (the GPOS ``kern`` feature when
    the font has one, else the legacy ``kern`` table; see
    ``_get_kerning_value``), as browsers shaping with HarfBuzz apply it. White
    space is collapsed as renderers draw it (``collapse_svg_white_space``). A
    character the face lacks counts as an em box (advance and ink one em wide,
    ``hhea`` ascent to descent high), as a fallback font draws it. The height
    is the tallest glyph's top to the lowest glyph's bottom.
    """
    font = load_font(face)
    if font is None:
        # Fallback approximation
        return len(str(text)) * float(font_size) * 0.6, float(font_size)

    hmtx = font["hmtx"]
    head = font["head"]
    units_per_em = head.unitsPerEm

    t = _measured_cmap(font)

    total_advance_width = 0
    ymaxes = []
    ymins = []

    text_str = collapse_svg_white_space(str(text))
    if not text_str:
        return 0.0, 0.0

    previous_glyph_name: Optional[str] = None

    for char in text_str:
        # None for a character the face lacks.
        glyph_name = t.get(ord(char))

        if glyph_name is None:
            # The face lacks the character, so a renderer draws it from a
            # fallback font, typically one em wide (CJK): an em box.
            ymax, ymin = font["hhea"].ascent, font["hhea"].descent
        else:
            try:
                ymax, ymin, _ = _glyph_bounds(font, glyph_name)
            except Exception:
                ymax, ymin = 0, 0

        ymaxes.append(ymax)
        ymins.append(ymin)

        # Apply kerning adjustment between consecutive glyphs
        if previous_glyph_name is not None:
            kerning_value = _get_kerning_value(font, previous_glyph_name, glyph_name)
            total_advance_width += kerning_value

        # Calculate horizontal advance
        if glyph_name is None:
            advance_width = units_per_em
        else:
            try:
                advance_width = hmtx[glyph_name][0]
            except KeyError:
                advance_width = units_per_em

        total_advance_width += advance_width

        previous_glyph_name = glyph_name

    # SVG font-size values are user units (CSS px in the generated SVG).
    scale_factor = float(font_size) / units_per_em
    text_width_pixels = total_advance_width * scale_factor
    text_height_pixels = (max(ymaxes) + abs(min(ymins))) * scale_factor

    return text_width_pixels, text_height_pixels


def measurement_font_file(
    font_family: str,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
) -> Optional[str]:
    """Return the name of the face ``calculate_svg_bbox_dimensions`` measures with.

    For example ``"LiberationSans-Regular"``, or ``user-<hash>-<style>`` for a
    registered font (``FontFace.measurement_name``); None when no bundled font
    exists and the measurement falls back to its approximation.
    """
    face = resolve_font_face(font_family, font_weight, font_style)
    return face.measurement_name if face is not None else None


def calculate_svg_bbox_dimensions(
    text,
    font_family,
    font_size,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
):
    """
    Calculates text dimensions in SVG user units: the advance width and the
    glyph-extent height, with the font ``font_family`` resolves to
    (``gbdraw.core.fonts.resolve_font_face``): a bundled Liberation face, a
    registered font file, or the file fontconfig finds for the family.

    Without a readable font the width is a per-character estimate and the
    height is ``font_size``.
    """
    face = resolve_font_face(font_family, font_weight, font_style)
    if face is not None:
        return _face_text_box(face, text, font_size)

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


def measure_mixed_text_line(
    parts: Iterable[tuple[str, str]],
    font_family: str,
    font_size: float,
    *,
    font_weight: FontWeight = "normal",
) -> tuple[float, float]:
    """Width and height of one text line drawn as styled tspans.

    ``parts`` are ``(text, font_style)`` pairs in drawing order. The width is
    the painted SVG width of the parts added up
    (``calculate_svg_bbox_dimensions``, italic parts measured italic). The
    height is the tallest part's height. No parts measure ``(0.0, 0.0)``.
    """
    widths: list[float] = []
    heights: list[float] = []
    for text, font_style in parts:
        width, height = calculate_svg_bbox_dimensions(
            text, font_family, font_size, font_weight=font_weight, font_style=font_style
        )
        widths.append(float(width))
        heights.append(float(height))
    return compensated_sum(widths), max(heights, default=0.0)


@dataclass(frozen=True)
class TextInk:
    """Where a text's letters are painted, in em, from the start of its advance on the alphabetic baseline.

    ``left`` and ``right`` bound the ink horizontally, ``advance`` is the
    width ``text-anchor`` aligns, and ``top`` and ``bottom`` bound the ink
    vertically (positive above the baseline, as in the font).
    """

    left: float
    right: float
    advance: float
    top: float
    bottom: float


def text_ink(
    text: str,
    font_family: str,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
) -> Optional[TextInk]:
    """The painted extent of ``text`` in the font ``font_family`` resolves to (the ink box, LP-D25).

    Kerning and missing characters are measured as in ``calculate_svg_bbox_dimensions``;
    None when the font cannot be read or nothing is painted.
    """
    face = resolve_font_face(font_family, font_weight, font_style)
    return _face_text_ink(face, collapse_svg_white_space(str(text))) if face is not None else None


@functools.lru_cache(maxsize=8192)
def _face_text_ink(face: FontFace, text: str) -> Optional[TextInk]:
    font = load_font(face)
    if font is None or not text:
        return None
    units_per_em = float(font["head"].unitsPerEm)
    cmap = _measured_cmap(font)
    origin = 0.0
    left = right = math.inf
    top = bottom = 0.0
    previous: Optional[str] = None
    for char in text:
        glyph_name = cmap.get(ord(char))
        if previous is not None:
            origin += _get_kerning_value(font, previous, glyph_name)
        if glyph_name is None:
            # A fallback font paints an em box (``_face_text_box``).
            advance, ink_left, ink_right = units_per_em, 0.0, units_per_em
            ink_top, ink_bottom = float(font["hhea"].ascent), float(font["hhea"].descent)
        else:
            advance, ink_left, ink_right, ink_top, ink_bottom = _glyph_facts(font, glyph_name)
        if ink_right > ink_left or ink_top > ink_bottom:
            left = min(left, origin + ink_left)
            right = origin + ink_right if math.isinf(right) else max(right, origin + ink_right)
            top = max(top, ink_top)
            bottom = min(bottom, ink_bottom)
        origin += advance
        previous = glyph_name
    if math.isinf(left):
        return None
    return TextInk(
        left=left / units_per_em,
        right=right / units_per_em,
        advance=origin / units_per_em,
        top=top / units_per_em,
        bottom=bottom / units_per_em,
    )


# Where ``text-anchor`` puts the anchor point along a text's advance.
TEXT_ANCHOR_SHARE: Mapping[str, float] = {"start": 0.0, "middle": 0.5, "end": 1.0}


def text_ink_offsets(
    text: str,
    font_family: str,
    font_size: float,
    *,
    dominant_baseline: str,
    text_anchor: str,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
) -> Optional[tuple[float, float, float, float]]:
    """The box ``text`` paints at ``font_size``, in px from the point it is drawn at.

    Returns (left, right, top, bottom), y down, for ``dominant_baseline`` and
    ``text_anchor``, in the face of ``font_weight`` and ``font_style``; None
    when the font cannot be measured.
    """

    return text_runs_ink_offsets(
        ((text, font_style),),
        font_family,
        font_size,
        dominant_baseline=dominant_baseline,
        text_anchor=text_anchor,
        font_weight=font_weight,
    )


# The style a run carries through ``collapse_runs_white_space``.
_RunStyle = TypeVar("_RunStyle")


def collapse_runs_white_space(runs: Sequence[tuple[str, _RunStyle]]) -> list[tuple[str, _RunStyle]]:
    """``runs`` (``(text, font_style)``) with white space collapsed as a renderer draws one text of ``tspan``s.

    Each run collapses as ``collapse_svg_white_space`` does, and a run's leading
    space after a space (ending an earlier run) is not drawn: browsers collapse
    white space across ``tspan`` boundaries. Nor are the spaces at the start and
    the end of the whole text.
    """

    collapsed: list[tuple[str, _RunStyle]] = []
    after_space = True
    for run, font_style in runs:
        text = collapse_svg_white_space(run)
        if after_space and text.startswith(" "):
            text = text[1:]
        if text:
            after_space = text.endswith(" ")
        collapsed.append((text, font_style))
    # The trailing space of the whole text: that of its last run with any text left.
    for index in range(len(collapsed) - 1, -1, -1):
        text, font_style = collapsed[index]
        if text:
            collapsed[index] = (text.rstrip(" "), font_style)
            break
    return collapsed


def text_runs_ink_offsets(
    runs: Sequence[tuple[str, str | None]],
    font_family: str,
    font_size: float,
    *,
    dominant_baseline: str,
    text_anchor: str,
    font_weight: FontWeight = "normal",
) -> Optional[tuple[float, float, float, float]]:
    """``text_ink_offsets`` of one text drawn as styled runs (``tspan``), ``runs`` being ``(text, font_style)``.

    Each run is measured in its own face and starts at the end of the previous
    run's advance, white space collapsed across runs (``collapse_runs_white_space``);
    ``text_anchor`` aligns the runs' total advance. A run that paints nothing
    only advances.
    """

    metrics = get_font_vertical_metrics(font_family, font_weight)
    shift = None if metrics is None else dominant_baseline_shift_em(dominant_baseline, metrics)
    if shift is None:
        return None
    size = float(font_size)
    advance = 0.0
    left = right = top = bottom = math.nan
    for text, font_style in collapse_runs_white_space(runs):
        ink = text_ink(text, font_family, font_weight, font_style)
        if ink is None:
            if resolve_font_face(font_family, font_weight, font_style) is None:
                return None
            advance += float(calculate_svg_bbox_dimensions(text, font_family, size, font_weight, font_style)[0]) / size
            continue
        run_left, run_right = advance + ink.left, advance + ink.right
        left = run_left if math.isnan(left) else min(left, run_left)
        right = run_right if math.isnan(right) else max(right, run_right)
        top = ink.top if math.isnan(top) else max(top, ink.top)
        bottom = ink.bottom if math.isnan(bottom) else min(bottom, ink.bottom)
        advance += ink.advance
    if math.isnan(left):
        return None
    origin = -TEXT_ANCHOR_SHARE[text_anchor] * advance * size
    return (
        origin + left * size,
        origin + right * size,
        (shift - top) * size,
        (shift - bottom) * size,
    )


@dataclass(frozen=True)
class GlyphOutline:
    """One glyph of a text as it is painted at its font size, in px: ``middle`` is the middle of its advance from the
    start of the text's advance, and ``contours`` its outline as closed polygons (curves flattened), x from that
    middle and y up from the alphabetic baseline. A character the face lacks is an em box (``_face_text_box``)."""

    middle: float
    contours: tuple[tuple[tuple[float, float], ...], ...]


class _FlatPen(BasePen):
    """A glyph's contours as polygons, each curve as straight segments."""

    def __init__(self, glyph_set: Mapping[str, object]) -> None:
        super().__init__(glyph_set)
        self.contours: list[list[tuple[float, float]]] = []

    def _moveTo(self, point: tuple[float, float]) -> None:
        self.contours.append([point])

    def _lineTo(self, point: tuple[float, float]) -> None:
        self.contours[-1].append(point)

    def _curveToOne(self, p1: tuple[float, float], p2: tuple[float, float], p3: tuple[float, float]) -> None:
        (x0, y0) = self.contours[-1][-1]
        for step in range(1, _CURVE_STEPS + 1):
            t = step / _CURVE_STEPS
            a, b, c, d = (1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t * t, t**3
            self.contours[-1].append(
                (a * x0 + b * p1[0] + c * p2[0] + d * p3[0], a * y0 + b * p1[1] + c * p2[1] + d * p3[1])
            )

    def _qCurveToOne(self, p1: tuple[float, float], p2: tuple[float, float]) -> None:
        (x0, y0) = self.contours[-1][-1]
        for step in range(1, _CURVE_STEPS + 1):
            t = step / _CURVE_STEPS
            a, b, c = (1 - t) ** 2, 2 * (1 - t) * t, t * t
            self.contours[-1].append((a * x0 + b * p1[0] + c * p2[0], a * y0 + b * p1[1] + c * p2[1]))


# Straight segments per curve of a glyph outline: a quarter circle of an em-high letter at 60 px strays less than
# 0.1 px from its chords.
_CURVE_STEPS = 8


@functools.lru_cache(maxsize=4096)
def _face_glyph_contours(face: FontFace, glyph_name: str) -> tuple[tuple[tuple[float, float], ...], ...]:
    font = load_font(face)
    if font is None:
        return ()
    glyph_set = font.getGlyphSet()
    pen = _FlatPen(glyph_set)
    glyph_set[glyph_name].draw(pen)
    return tuple(tuple(contour) for contour in pen.contours if len(contour) > 2)


@functools.lru_cache(maxsize=2048)
def text_runs_glyph_outlines(
    runs: tuple[tuple[str, str | None], ...],
    font_family: str,
    font_size: float,
    font_weight: FontWeight = "normal",
) -> tuple[GlyphOutline, ...]:
    """The painted glyphs of one text drawn as styled runs (``tspan``), ``runs`` being ``(text, font_style)``, laid
    out as ``text_runs_ink_offsets`` measures them: each run in its own face, kerned within it, white space collapsed
    across runs. Empty when a face cannot be read."""

    size = float(font_size)
    glyphs: list[GlyphOutline] = []
    origin = 0.0
    for text, font_style in collapse_runs_white_space(runs):
        face = resolve_font_face(font_family, font_weight, font_style)
        font = None if face is None else load_font(face)
        if face is None or font is None:
            return ()
        scale = size / float(font["head"].unitsPerEm)
        cmap = _measured_cmap(font)
        previous: Optional[str] = None
        for char in text:
            glyph_name = cmap.get(ord(char))
            if previous is not None:
                origin += _get_kerning_value(font, previous, glyph_name) * scale
            if glyph_name is None:
                em = float(font["head"].unitsPerEm)
                top, bottom = float(font["hhea"].ascent), float(font["hhea"].descent)
                advance = em * scale
                contours: tuple[tuple[tuple[float, float], ...], ...] = (((0.0, bottom), (em, bottom), (em, top), (0.0, top)),)
            else:
                advance = float(font["hmtx"][glyph_name][0]) * scale
                contours = _face_glyph_contours(face, glyph_name)
            if contours:
                half = advance / 2
                glyphs.append(
                    GlyphOutline(
                        origin + half,
                        tuple(tuple((x * scale - half, y * scale) for x, y in contour) for contour in contours),
                    )
                )
            origin += advance
            previous = glyph_name
    return tuple(glyphs)


@dataclass(frozen=True)
class FontVerticalMetrics:
    """Vertical font metrics as fractions of the font size (em units)."""

    ascent: float
    descent: float
    x_height: float


def get_font_vertical_metrics(
    font_family: str,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
) -> Optional[FontVerticalMetrics]:
    """Return the ascent, descent, and x-height of the font used for layout.

    The font is resolved like :func:`calculate_svg_bbox_dimensions`. Ascent and
    descent follow FreeType, which browsers on Linux use: the OS/2 typo
    metrics when USE_TYPO_METRICS is set, otherwise the hhea metrics.
    """
    face = resolve_font_face(font_family, font_weight, font_style)
    return _face_vertical_metrics(face) if face is not None else None


def _line_box_units(font: TTFont) -> tuple[float, float]:
    """The ascent and descent (both positive, font units) FreeType gives a line of the face."""
    os2 = font["OS/2"] if "OS/2" in font else None
    hhea = font["hhea"] if "hhea" in font else None
    if os2 is not None and os2.fsSelection & (1 << 7):
        return os2.sTypoAscender, -os2.sTypoDescender
    if hhea is not None and (hhea.ascent or hhea.descent):
        return hhea.ascent, -hhea.descent
    if os2 is not None:
        return os2.sTypoAscender, -os2.sTypoDescender
    return font["head"].yMax, -font["head"].yMin


@functools.lru_cache(maxsize=64)
def _face_vertical_metrics(face: FontFace) -> Optional[FontVerticalMetrics]:
    font = load_font(face)
    if font is None:
        return None
    units_per_em = float(font["head"].unitsPerEm)
    os2 = font["OS/2"] if "OS/2" in font else None
    ascent, descent = _line_box_units(font)
    x_height = float(getattr(os2, "sxHeight", 0) or 0)
    if x_height <= 0:
        glyph_name = (font.getBestCmap() or {}).get(ord("x"))
        if glyph_name is not None:
            x_height = float(_glyph_bounds(font, glyph_name)[0] or 0)
    if x_height <= 0:
        x_height = ascent / 2.0
    return FontVerticalMetrics(
        ascent=ascent / units_per_em,
        descent=descent / units_per_em,
        x_height=x_height / units_per_em,
    )


def drawn_text_box(
    text: str,
    font_family: str,
    font_size: float,
    font_weight: FontWeight = "normal",
    font_style: str | None = "normal",
) -> tuple[float, float]:
    """The box ``text`` is drawn in at ``font_size`` px: (width, height) in px.

    The width is the text's advance (``calculate_svg_bbox_dimensions``, what a
    browser's ``getComputedTextLength`` reports) or, when a last glyph paints
    past its advance, the right edge of its ink (``text_ink``), so the box
    contains what is painted. The height is the face's line box (ascent plus
    descent, ``get_font_vertical_metrics``), the box a browser reports for a
    text element, whatever characters the text has. Text that collapses to
    nothing measures ``(0.0, 0.0)``; without a readable font both values are
    ``calculate_svg_bbox_dimensions``' estimate. Layout reserves exactly this
    box and adds its own named margins (OV-165).
    """

    width, height = calculate_svg_bbox_dimensions(text, font_family, font_size, font_weight, font_style)
    metrics = get_font_vertical_metrics(font_family, font_weight, font_style)
    if metrics is None or not collapse_svg_white_space(str(text)):
        return float(width), float(height)
    size = float(font_size)
    ink = text_ink(text, font_family, font_weight, font_style)
    if ink is not None:
        width = max(width, ink.right * size)
    return float(width), (metrics.ascent + metrics.descent) * size


def drawn_runs_box(
    runs: Sequence[tuple[str, str]],
    font_family: str,
    font_size: float,
) -> tuple[float, float]:
    """``drawn_text_box`` of one text drawn as styled runs (``tspan``), ``runs`` being ``(text, font_style)``.

    Each run is measured in its own face and starts at the end of the previous
    run's advance, white space collapsed across runs
    (``collapse_runs_white_space``). The width is the last
    advance or, when a run's box ends past it (an italic glyph painting past
    its advance), that box's right edge; the height is the tallest run's line
    box. One run measures as ``drawn_text_box``. Unlike
    ``measure_mixed_text_line`` (advance width, glyph height), the box
    contains the painted ink and has the line height a browser reports, so a
    text measures the same whether it is one run or several.
    """

    if len(runs) == 1:
        text, font_style = runs[0]
        return drawn_text_box(text, font_family, font_size, "normal", font_style)
    right = origin = height = 0.0
    for text, font_style in collapse_runs_white_space(runs):
        width, run_height = drawn_text_box(text, font_family, font_size, "normal", font_style)
        right = max(right, origin + width)
        origin += float(calculate_svg_bbox_dimensions(text, font_family, font_size, "normal", font_style)[0])
        height = max(height, run_height)
    return max(right, origin), height


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


TextRunStyle = Literal["normal", "italic"]
TextRun = tuple[str, TextRunStyle]

# The markup user text drawn in a Legend caption or a region label may use (OV-2760, OV-2808): non-empty
# ``<i>...</i>`` runs, lower case and not nested, with no other ``<``, ``>`` or ``&`` and no character XML rewrites
# or refuses (a carriage return, the other control characters but tab and line feed, a lone surrogate, U+FFFE,
# U+FFFF). In this subset the XML parse and a plain split give the same runs, so the Web
# (services/legend-layout.js ``legendCaptionRuns``) needs no XML parser.
_MARKUP_CHARACTER = r"[^<>&\x00-\x08\x0b-\x1f\ud800-\udfff\ufffe\uffff]"
_ITALIC_MARKUP = re.compile(rf"{_MARKUP_CHARACTER}*(?:<i>{_MARKUP_CHARACTER}+</i>{_MARKUP_CHARACTER}*)+")


def italic_markup_runs(text: str) -> tuple[TextRun, ...]:
    """The styled runs user text (a Legend caption, a region label) is drawn and measured as, in drawing order.

    A text in the italic subset gives the runs ``parse_mixed_content_text`` (the owner of ``<i>`` markup in
    Species, record names and plot titles) reads from it; any other text is one literal upright run.
    """

    if _ITALIC_MARKUP.fullmatch(text) is None:
        return ((text, "normal"),)
    return tuple(
        (str(part["text"]), "italic" if part["italic"] else "normal")
        for part in parse_mixed_content_text(text)
        if part["text"]
    )


__all__ = [
    "FontMetricsTable",
    "FontVerticalMetrics",
    "GlyphOutline",
    "TextInk",
    "TextRun",
    "TextRunStyle",
    "calculate_svg_bbox_dimensions",
    "collapse_runs_white_space",
    "collapse_svg_white_space",
    "create_text_element",
    "dominant_baseline_shift_em",
    "drawn_runs_box",
    "drawn_text_box",
    "font_metrics_table",
    "font_pair_kerning_table",
    "get_font_vertical_metrics",
    "italic_markup_runs",
    "measure_mixed_text_line",
    "measurement_font_file",
    "parse_mixed_content_text",
    "text_ink",
    "text_ink_offsets",
    "text_runs_glyph_outlines",
    "text_runs_ink_offsets",
]


