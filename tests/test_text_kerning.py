"""Pair kerning in text measurement follows the bundled fonts' GPOS ``kern``.

The expected adjustments were read with HarfBuzz 11.2.1, an independent
shaper (``hb-shape --output-format=json`` with and without ``--features=-kern``,
sum of ``ax``), and are checked against values read directly from the font
tables below.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString
from fontTools.ttLib import TTFont

from gbdraw.core.text import (
    _get_kerning_value,
    font_pair_kerning_table,
    get_text_bbox_size_pixels,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
FACES = sorted(path.stem for path in DATA_DIR.glob("Liberation*.ttf"))

# hb-shape 11.2.1, font design units (2048 per em).
HARFBUZZ_PAIR_KERNING = {
    "LiberationMono-Bold": {"AV": 0, "To": 0, "LT": 0, "Yo": 0, "אל": 0},
    "LiberationMono-BoldItalic": {"AV": 0, "To": 0, "LT": 0, "Yo": 0, "אל": 0},
    "LiberationMono-Italic": {"AV": 0, "To": 0, "LT": 0, "Yo": 0, "אל": 0},
    "LiberationMono-Regular": {"AV": 0, "To": 0, "LT": 0, "Yo": 0, "אל": 0},
    "LiberationSans-Bold": {"AV": -152, "To": -152, "LT": -152, "Yo": -152, "אל": -41},
    "LiberationSans-BoldItalic": {"AV": -152, "To": -76, "LT": -152, "Yo": -76, "אל": -41},
    "LiberationSans-Italic": {"AV": -113, "To": -188, "LT": -152, "Yo": -113, "אל": -41},
    "LiberationSans-Regular": {"AV": -152, "To": -227, "LT": -152, "Yo": -188, "אל": -41},
    "LiberationSerif-Bold": {"AV": -264, "To": -188, "LT": -188, "Yo": -227, "אל": 0},
    "LiberationSerif-BoldItalic": {"AV": -152, "To": -188, "LT": -37, "Yo": -227, "אל": 0},
    "LiberationSerif-Italic": {"AV": -102, "To": -188, "LT": -41, "Yo": -188, "אל": 0},
    "LiberationSerif-Regular": {"AV": -264, "To": -143, "LT": -188, "Yo": -205, "אל": 0},
}

PAIRS = [(face, pair) for face in FACES for pair in HARFBUZZ_PAIR_KERNING[face]]


def _font(face: str) -> TTFont:
    return TTFont(str(DATA_DIR / f"{face}.ttf"))


def _glyphs(font: TTFont, text: str) -> list[str]:
    cmap = font["cmap"].getcmap(3, 1).cmap
    return [cmap[ord(char)] for char in text]


def _table_gpos_kern(font: TTFont) -> dict[tuple[str, str], int]:
    """Read the GPOS ``kern`` feature pairs straight from the tables."""
    table = font["GPOS"].table
    indices = {
        index
        for record in table.FeatureList.FeatureRecord
        if record.FeatureTag == "kern"
        for index in record.Feature.LookupListIndex
    }
    pairs: dict[tuple[str, str], int] = {}
    for index in sorted(indices):
        lookup = table.LookupList.Lookup[index]
        assert lookup.LookupType == 2
        for subtable in lookup.SubTable:
            assert subtable.Format == 1
            for first, pair_set in zip(subtable.Coverage.glyphs, subtable.PairSet):
                for record in pair_set.PairValueRecord:
                    value = record.Value1.XAdvance
                    if value:
                        pairs[(first, record.SecondGlyph)] = value
    return pairs


def _table_legacy_kern(font: TTFont) -> dict[tuple[str, str], int]:
    if "kern" not in font:
        return {}
    (subtable,) = font["kern"].kernTables
    return {pair: value for pair, value in subtable.kernTable.items() if value}


def test_bundled_faces_are_the_twelve_liberation_faces() -> None:
    assert FACES == sorted(HARFBUZZ_PAIR_KERNING)


@pytest.mark.parametrize(("face", "pair"), PAIRS)
def test_pair_kerning_matches_harfbuzz(face: str, pair: str) -> None:
    font = _font(face)
    left, right = _glyphs(font, pair)
    expected = HARFBUZZ_PAIR_KERNING[face][pair]

    assert _get_kerning_value(font, left, right) == expected
    assert _table_gpos_kern(font).get((left, right), 0) == expected


@pytest.mark.parametrize(("face", "pair"), PAIRS)
def test_pair_width_applies_the_kerning(face: str, pair: str) -> None:
    font = _font(face)
    path = str(DATA_DIR / f"{face}.ttf")
    left, right = _glyphs(font, pair)
    kerning = HARFBUZZ_PAIR_KERNING[face][pair]
    units_per_em = font["head"].unitsPerEm
    (left_advance, left_lsb), (right_advance, right_lsb) = font["hmtx"][left], font["hmtx"][right]
    left_xmax, right_xmax = font["glyf"][left].xMax, font["glyf"][right].xMax

    # SVG units: advances plus the pair adjustment, scaled by size / em.
    svg_width, _ = get_text_bbox_size_pixels(path, pair, 14, 96, svg_units=True)
    assert svg_width == pytest.approx((left_advance + right_advance + kerning) * 14 / units_per_em)

    # Point units: the tight ink width of the two glyphs plus the adjustment.
    tight_units = (
        left_lsb - abs(left_lsb) + left_advance - abs(left_xmax - left_advance)
        + kerning
        - abs(right_lsb) + right_advance - abs(right_xmax - right_advance)
    )
    width, _ = get_text_bbox_size_pixels(path, pair, 14, 96)
    assert width == pytest.approx(tight_units * 14 * 96 / (72 * units_per_em))


@pytest.mark.parametrize("face", FACES)
def test_kerning_table_is_the_gpos_kern_feature(face: str) -> None:
    font = _font(face)
    gpos_pairs = _table_gpos_kern(font)

    assert font_pair_kerning_table(font) == gpos_pairs
    # Latin, Greek and Cyrillic pairs equal the legacy table pair for pair;
    # GPOS adds only the Hebrew lookup of the Sans faces.
    legacy_pairs = _table_legacy_kern(font)
    assert {pair: gpos_pairs[pair] for pair in legacy_pairs} == legacy_pairs
    extra = set(gpos_pairs) - set(legacy_pairs)
    assert len(extra) == (1107 if face.startswith("LiberationSans") else 0)


def _reloaded(font: TTFont, tmp_path: Path, name: str) -> str:
    path = tmp_path / f"{name}.ttf"
    buffer = io.BytesIO()
    font.save(buffer)
    path.write_bytes(buffer.getvalue())
    return str(path)


def test_gpos_kern_takes_precedence_over_the_kern_table(tmp_path: Path) -> None:
    font = _font("LiberationSans-Regular")
    font["kern"].kernTables[0].kernTable[("A", "V")] = -500
    both = TTFont(_reloaded(font, tmp_path, "both"))
    assert _get_kerning_value(both, "A", "V") == -152

    for record in both["GPOS"].table.FeatureList.FeatureRecord:
        if record.FeatureTag == "kern":
            record.FeatureTag = "zzzz"
    legacy_only = TTFont(_reloaded(both, tmp_path, "legacy-only"))
    assert _get_kerning_value(legacy_only, "A", "V") == -500
    assert _get_kerning_value(legacy_only, "T", "o") == -227


def test_class_pair_kerning_applies_format_2(tmp_path: Path) -> None:
    font = _font("LiberationMono-Regular")
    addOpenTypeFeaturesFromString(
        font,
        "@LEFT = [A T]; @RIGHT = [V o];\n"
        "feature kern { pos A Y -30; pos @LEFT @RIGHT -80; } kern;\n",
        tables=["GPOS"],
    )
    path = _reloaded(font, tmp_path, "class-kern")
    kerned = TTFont(path)
    formats = {
        subtable.Format
        for lookup in kerned["GPOS"].table.LookupList.Lookup
        for subtable in lookup.SubTable
    }
    assert formats == {1, 2}

    assert _get_kerning_value(kerned, "A", "Y") == -30
    assert _get_kerning_value(kerned, "A", "V") == -80
    assert _get_kerning_value(kerned, "T", "o") == -80
    assert _get_kerning_value(kerned, "T", "A") == 0
    assert _get_kerning_value(kerned, "V", "A") == 0
    assert font_pair_kerning_table(kerned) == {
        ("A", "V"): -80,
        ("A", "Y"): -30,
        ("A", "o"): -80,
        ("T", "V"): -80,
        ("T", "o"): -80,
    }
    advance = kerned["hmtx"]["A"][0]
    width, _ = get_text_bbox_size_pixels(path, "TAVo", 10, 96, svg_units=True)
    assert width == pytest.approx((4 * advance - 80) * 10 / kerned["head"].unitsPerEm)


def test_missing_glyph_pairs_are_not_kerned() -> None:
    font = _font("LiberationSans-Regular")
    assert _get_kerning_value(font, 0, "V") == 0
    assert _get_kerning_value(font, "A", None) == 0
