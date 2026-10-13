"""White space and missing characters are measured as renderers draw them (OV-165, part 1).

gbdraw writes text without ``xml:space``, so renderers turn tabs and line
breaks into spaces and runs of spaces into one; no-break spaces are kept.
Spaces at the ends of a styled fragment (``tspan``) are kept, because a
renderer drops spaces only at the ends of the whole text. A character the
bundled face lacks is drawn from a fallback font one em wide (Chromium draws a
CJK character exactly one em wide), so it counts as an em box. Interior spaces
still add nothing to the point-mode width (deferred, see OV-165 part 2).
Expected values come from the font tables; the SVG-mode pins were read in
Chromium (``getComputedTextLength`` at ``font-size`` = units per em, with the
bundled face loaded through ``@font-face``).
"""

from __future__ import annotations

import os
import uuid
from pathlib import Path

import pytest
from fontTools.ttLib import TTFont

from gbdraw.core.text import (
    calculate_bbox_dimensions,
    calculate_svg_bbox_dimensions,
    collapse_svg_white_space,
    font_pair_kerning_table,
    get_text_bbox_size_pixels,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
FACES = sorted(path.stem for path in DATA_DIR.glob("Liberation*.ttf"))


def _path(face: str) -> str:
    return str(DATA_DIR / f"{face}.ttf")


def _svg_units(face: str, text: str) -> float:
    """SVG-mode width in font units (font-size = units per em)."""
    return get_text_bbox_size_pixels(_path(face), text, 2048, 96, svg_units=True)[0]


def _point_units(face: str, text: str) -> tuple[float, float]:
    """Point-mode width and height in font units (scale 1)."""
    return get_text_bbox_size_pixels(_path(face), text, 72, 2048)


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
    one_space = _svg_units(face, "a b")
    for variant in ("a  b", "a\tb", "a\nb", "a\r\nb", "a \t b"):
        assert _svg_units(face, variant) == one_space, repr(variant)
        assert _point_units(face, variant) == _point_units(face, "a b"), repr(variant)


@pytest.mark.parametrize("face", FACES)
def test_no_break_spaces_are_kept_and_counted(face: str) -> None:
    font = TTFont(_path(face))
    cmap = font["cmap"].getcmap(3, 1).cmap
    nbsp = cmap[0xA0]
    kerning = font_pair_kerning_table(font)

    # Two no-break spaces are drawn as two; two spaces as one.
    assert _svg_units(face, "a  b") == _svg_units(face, "a b")
    assert _svg_units(face, "a  b") - _svg_units(face, "a b") == (
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

    assert _svg_units(face, "Plain  ") == _svg_units(face, "Plain ")
    assert _svg_units(face, "Plain ") == _svg_units(face, "Plain") + space_width
    assert _svg_units(face, "  tail") == _svg_units(face, " tail")


@pytest.mark.parametrize("face", FACES)
def test_missing_character_counts_as_an_em_box(face: str) -> None:
    font = TTFont(_path(face))
    em = font["head"].unitsPerEm
    height = font["hhea"].ascent - font["hhea"].descent

    assert _point_units(face, "中") == (em, height)
    assert _svg_units(face, "中") == em
    assert _point_units(face, "a中b")[0] == _point_units(face, "ab")[0] + em
    assert _point_units(face, "日本語")[0] == 3 * em


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
def test_svg_mode_matches_chromium_text_length(face: str, text: str) -> None:
    assert _svg_units(face, text) == CHROMIUM_TEXT_LENGTH[face][text]


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

    svg_width = calculate_svg_bbox_dimensions(text, f"'{family}', sans-serif", 14, 96)
    point_width = calculate_bbox_dimensions(text, f"'{family}', sans-serif", 14, 72)

    assert svg_width == get_text_bbox_size_pixels(font_file, text, 14, 96, svg_units=True)
    assert point_width == get_text_bbox_size_pixels(font_file, text, 14, 72)
    assert svg_width != get_text_bbox_size_pixels(_path("LiberationSans-Regular"), text, 14, 96, svg_units=True)


def test_ov680_fontconfig_fallback_for_another_family_is_not_used(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    family = f"Gbdraw OV680 {uuid.uuid4().hex[:8]}"
    font_file = _renamed_face("LiberationSerif-Regular", "Some Other Family", tmp_path / "other.ttf")
    _fake_fc_match(tmp_path, monkeypatch, "Some Other Family", font_file)
    text = "Hypothetical protein"

    assert calculate_svg_bbox_dimensions(text, f"'{family}', serif", 14, 96) == get_text_bbox_size_pixels(
        _path("LiberationSerif-Regular"), text, 14, 96, svg_units=True
    )
