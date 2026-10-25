"""White space and missing characters are measured as renderers draw them (OV-165, part 1).

gbdraw writes text without ``xml:space``, so renderers turn tabs and line
breaks into spaces and runs of spaces into one; no-break spaces are kept.
Spaces at the ends of a styled fragment (``tspan``) are kept, because a
renderer drops spaces only at the ends of the whole text. A character the
bundled face lacks is drawn from a fallback font one em wide (Chromium draws a
CJK character exactly one em wide), so it counts as an em box.
Expected values come from the font tables; the width pins were read in
Chromium (``getComputedTextLength`` at ``font-size`` = units per em, with the
bundled face loaded through ``@font-face``).
"""

from __future__ import annotations

import os
import uuid
from pathlib import Path

import pytest
from fontTools.ttLib import TTFont

from gbdraw.core.fonts import FontFace
from gbdraw.core.text import (
    _face_text_box,
    calculate_svg_bbox_dimensions,
    collapse_svg_white_space,
    font_pair_kerning_table,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
FACES = sorted(path.stem for path in DATA_DIR.glob("Liberation*.ttf"))


def _path(face: str) -> str:
    return str(DATA_DIR / f"{face}.ttf")


def _box(font_file: str, text: str, size: float) -> tuple[float, float]:
    """The measured (width, height) of ``text`` in one font file at ``size`` px."""
    return _face_text_box(FontFace(font_file), text, size)


def _box_units(face: str, text: str) -> tuple[float, float]:
    """Width and height in font units (font-size = units per em)."""
    return _box(_path(face), text, 2048)


def _advance_units(face: str, text: str) -> float:
    """Width in font units (font-size = units per em)."""
    return _box_units(face, text)[0]


@pytest.mark.parametrize(
    ("raw", "drawn"),
    [
        ("a  b", "a b"),
        ("a\tb", "a b"),
        ("a\nb", "a b"),
        ("a\r\nb", "a b"),
        ("a \t\n b", "a b"),
        ("  a  ", " a "),
        ("a  b", "a  b"),
    ],
)
def test_white_space_collapses_as_drawn_and_keeps_no_break_spaces(raw: str, drawn: str) -> None:
    assert collapse_svg_white_space(raw) == drawn


@pytest.mark.parametrize("face", FACES)
def test_repeated_spaces_tabs_and_line_breaks_measure_as_one_space(face: str) -> None:
    one_space = _advance_units(face, "a b")
    for variant in ("a  b", "a\tb", "a\nb", "a\r\nb", "a \t b"):
        assert _advance_units(face, variant) == one_space, repr(variant)


@pytest.mark.parametrize("face", FACES)
def test_no_break_spaces_are_kept_and_counted(face: str) -> None:
    font = TTFont(_path(face))
    cmap = font["cmap"].getcmap(3, 1).cmap
    nbsp = cmap[0xA0]
    kerning = font_pair_kerning_table(font)

    # Two no-break spaces are drawn as two; two spaces as one.
    assert _advance_units(face, "a  b") == _advance_units(face, "a b")
    assert _advance_units(face, "a  b") - _advance_units(face, "a b") == (
        font["hmtx"][nbsp][0] + kerning.get((nbsp, nbsp), 0)
    )


@pytest.mark.parametrize("face", FACES)
def test_fragment_keeps_one_boundary_space(face: str) -> None:
    # Styled fragments are measured one by one, so their end spaces stay;
    # only runs collapse (the renderer collapses them across the whole text).
    font = TTFont(_path(face))
    cmap = font["cmap"].getcmap(3, 1).cmap
    n, space = cmap[ord("n")], cmap[ord(" ")]
    kerning = font_pair_kerning_table(font)
    space_width = font["hmtx"][space][0] + kerning.get((n, space), 0)

    assert _advance_units(face, "Plain  ") == _advance_units(face, "Plain ")
    assert _advance_units(face, "Plain ") == _advance_units(face, "Plain") + space_width
    assert _advance_units(face, "  tail") == _advance_units(face, " tail")


@pytest.mark.parametrize("face", FACES)
def test_missing_character_counts_as_an_em_box(face: str) -> None:
    font = TTFont(_path(face))
    em = font["head"].unitsPerEm
    height = font["hhea"].ascent - font["hhea"].descent

    assert _box_units(face, "中") == (em, height)
    assert _advance_units(face, "中") == em
    assert _advance_units(face, "a中b") == _advance_units(face, "ab") + em
    assert _advance_units(face, "日本語") == 3 * em


# Chromium (Playwright) getComputedTextLength at font-size 2048 px, bundled face via @font-face.
CHROMIUM_TEXT_LENGTH = {
    "LiberationSans-Regular": {"a  b": 2847, "a\tb": 2847, "中": 2048, "a中b": 4326},
    "LiberationSerif-Regular": {"a  b": 2445, "a\tb": 2445, "中": 2048, "a中b": 3981},
    "LiberationMono-Regular": {"a  b": 3687, "a\tb": 3687, "中": 2048, "a中b": 4506},
}


@pytest.mark.parametrize(
    ("face", "text"),
    [(face, text) for face, values in CHROMIUM_TEXT_LENGTH.items() for text in values],
)
def test_width_matches_chromium_text_length(face: str, text: str) -> None:
    assert _advance_units(face, text) == CHROMIUM_TEXT_LENGTH[face][text]


# OV-680: a family fontconfig finds is measured with its own file, even when a
# bundled generic alias follows it in the list; fontconfig's fallback for a
# family it does not have is not used.


def _renamed_face(source: str, family: str, target: Path) -> str:
    font = TTFont(_path(source))
    name_table = font["name"]
    name_table.names = [record for record in name_table.names if record.nameID not in {1, 4, 6, 16, 17}]
    for name_id, value in ((1, family), (4, family), (6, family.replace(" ", "")), (16, family)):
        name_table.setName(value, name_id, 3, 1, 0x409)
    font.save(str(target))
    return str(target)


def _fake_fc_match(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, family: str, font_file: str) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "fc-match"
    script.write_text(f"#!/bin/sh\nprintf '%s\\t%s\\t0' '{family}' '{font_file}'\n")
    script.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")


def test_ov680_a_named_system_family_before_a_generic_is_measured_with_its_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    family = f"Gbdraw OV680 {uuid.uuid4().hex[:8]}"
    font_file = _renamed_face("LiberationSerif-Regular", family, tmp_path / "ov680.ttf")
    _fake_fc_match(tmp_path, monkeypatch, family, font_file)
    text = "Hypothetical protein"

    width = calculate_svg_bbox_dimensions(text, f"'{family}', sans-serif", 14)

    assert width == _box(font_file, text, 14)
    assert width != _box(_path("LiberationSans-Regular"), text, 14)


def test_ov680_fontconfig_fallback_for_another_family_is_not_used(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    family = f"Gbdraw OV680 {uuid.uuid4().hex[:8]}"
    font_file = _renamed_face("LiberationSerif-Regular", "Some Other Family", tmp_path / "other.ttf")
    _fake_fc_match(tmp_path, monkeypatch, "Some Other Family", font_file)
    text = "Hypothetical protein"

    assert calculate_svg_bbox_dimensions(text, f"'{family}', serif", 14) == _box(
        _path("LiberationSerif-Regular"), text, 14
    )


# OV-681: a CFF (.otf) face is measured from its outlines, as the same glyphs
# in a TrueType face are. OV-682: a face without the Windows BMP cmap (3, 1)
# is measured with its best cmap instead of raising.


def _box_font(path: Path, *, cff: bool) -> str:
    from fontTools.fontBuilder import FontBuilder
    from fontTools.pens.t2CharStringPen import T2CharStringPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen

    builder = FontBuilder(1000, isTTF=not cff)
    order = [".notdef", "A", "V"]
    builder.setupGlyphOrder(order)
    builder.setupCharacterMap({ord("A"): "A", ord("V"): "V"})
    boxes = {".notdef": (50, 0, 450, 700), "A": (40, -120, 610, 730), "V": (10, 0, 560, 680)}
    advances = {".notdef": 500, "A": 650, "V": 570}
    outlines = {}
    for name, (x0, y0, x1, y1) in boxes.items():
        pen = T2CharStringPen(advances[name], None) if cff else TTGlyphPen(None)
        pen.moveTo((x0, y0))
        pen.lineTo((x0, y1))
        pen.lineTo((x1, y1))
        pen.lineTo((x1, y0))
        pen.closePath()
        outlines[name] = pen.getCharString() if cff else pen.glyph()
    if cff:
        builder.setupCFF("GbdrawBox", {"FullName": "Gbdraw Box"}, outlines, {})
    else:
        builder.setupGlyf(outlines)
    builder.setupHorizontalMetrics({name: (advances[name], boxes[name][0]) for name in order})
    builder.setupHorizontalHeader(ascent=800, descent=-200)
    builder.setupNameTable({"familyName": "Gbdraw Box", "styleName": "Regular"})
    builder.setupOS2()
    builder.setupPost()
    builder.save(str(path))
    return str(path)


def test_ov681_a_cff_face_is_measured_from_its_outlines(tmp_path: Path) -> None:
    truetype = _box_font(tmp_path / "box.ttf", cff=False)
    cff = _box_font(tmp_path / "box.otf", cff=True)

    for text in ("A", "AV", "VAV"):
        assert _box(cff, text, 1000) == _box(truetype, text, 1000)
        assert _box(cff, text, 14) == _box(truetype, text, 14)
    assert _box(cff, "A", 1000)[1] == 850


def test_ov682_a_face_without_the_bmp_cmap_is_measured_with_its_best_cmap(tmp_path: Path) -> None:
    font = TTFont(_path("LiberationSans-Regular"))
    font["cmap"].tables = [table for table in font["cmap"].tables if (table.platformID, table.platEncID) != (3, 1)]
    target = tmp_path / "no-bmp-cmap.ttf"
    font.save(str(target))
    text = "Hypothetical protein"

    assert _box(str(target), text, 14) == _box(_path("LiberationSans-Regular"), text, 14)
