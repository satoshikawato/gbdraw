"""Custom fonts in exported files: SVG embedding, --font_file, PDF TrueType (D-05, D-07)."""

from __future__ import annotations

import base64
import io
import json
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
import svgwrite
from fontTools.ttLib import TTCollection, TTFont

import gbdraw.api.diagram as diagram_module
from gbdraw.api.render import render_to_bytes
from gbdraw.circular import circular_main
from gbdraw.cli_utils.session import render_main
from gbdraw.core import fonts
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.render import font_embedding
from gbdraw.web_support import export_fonts
from tests.utils.session_drawing_view import session_view

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
FONT_FACE = re.compile(r"@font-face\s*\{([^}]*)\}")


def _renamed_font(source: str, family: str, target: Path, *, weight: int = 400, fs_type: int = 0) -> Path:
    font = TTFont(str(DATA_DIR / f"{source}.ttf"))
    name_table = font["name"]
    name_table.names = [record for record in name_table.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
    style = "Bold" if weight >= 600 else "Regular"
    for name_id, value in ((1, family), (2, style), (4, f"{family} {style}"), (6, family.replace(" ", "")), (16, family)):
        name_table.setName(value, name_id, 3, 1, 0x409)
    font["OS/2"].usWeightClass = weight
    font["OS/2"].fsSelection = (font["OS/2"].fsSelection & ~0x60) | (0x20 if weight >= 600 else 0x40)
    font["OS/2"].fsType = fs_type
    font.save(str(target))
    return target


def _cff_font(path: Path, family: str) -> Path:
    """An OpenType (CFF) font with a curved glyph, built in place."""
    from fontTools.fontBuilder import FontBuilder
    from fontTools.pens.t2CharStringPen import T2CharStringPen

    builder = FontBuilder(1000, isTTF=False)
    order = [".notdef", "space", "O", "V"]
    builder.setupGlyphOrder(order)
    builder.setupCharacterMap({0x20: "space", ord("O"): "O", ord("V"): "V"})
    advances = {".notdef": 500, "space": 250, "O": 720, "V": 570}
    charstrings = {}
    for name in order:
        pen = T2CharStringPen(advances[name], None)
        if name == "O":
            pen.moveTo((360, 0))
            pen.curveTo((560, 0), (680, 160), (680, 350))
            pen.curveTo((680, 540), (560, 700), (360, 700))
            pen.curveTo((160, 700), (40, 540), (40, 350))
            pen.curveTo((40, 160), (160, 0), (360, 0))
            pen.closePath()
        elif name != "space":
            pen.moveTo((20, 0))
            pen.lineTo((20, 680))
            pen.lineTo((advances[name] - 20, 680))
            pen.lineTo((advances[name] - 20, 0))
            pen.closePath()
        charstrings[name] = pen.getCharString()
    builder.setupCFF(family.replace(" ", ""), {"FullName": family}, charstrings, {})
    builder.setupHorizontalMetrics({name: (advances[name], 20) for name in order})
    builder.setupHorizontalHeader(ascent=800, descent=-200)
    builder.setupNameTable({"familyName": family, "styleName": "Regular"})
    builder.setupOS2(usWeightClass=400, fsType=0)
    builder.setupPost()
    builder.save(str(path))
    return path


def _genbank(path: Path) -> Path:
    sequence = "".join(f"{index * 60 + 1:>9} {'acgt' * 15}\n" for index in range(6))
    path.write_text(
        "LOCUS       font_test                360 bp    DNA     circular UNA 01-JAN-2000\n"
        "DEFINITION  Font test.\nACCESSION   font_test\nVERSION     font_test\nKEYWORDS    .\n"
        "SOURCE      synthetic construct\n  ORGANISM  synthetic construct\n            .\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..360\n"
        "     CDS             41..125\n"
        '                     /gene="fontA"\n'
        '                     /product="hypothetical protein"\n'
        f"ORIGIN\n{sequence}//\n"
    )
    return path


def _svg(*texts: tuple[str, str, str]) -> str:
    """An SVG with one <text> per (font-family, font-weight, text)."""
    drawing = svgwrite.Drawing(size=("200px", "100px"))
    for index, (family, weight, value) in enumerate(texts):
        drawing.add(drawing.text(value, insert=(0, 20 * (index + 1)), font_family=family, font_weight=weight))
    return drawing.tostring()


def _faces(css: str) -> list[dict[str, str]]:
    faces = []
    for body in FONT_FACE.findall(css):
        declarations = dict(re.findall(r"(font-family|font-weight|font-style):\s*([^;]+?)\s*;", body))
        source = re.search(r"src:\s*(url\([^)]*\)[^;]*);", body)
        assert source, body
        declarations["src"] = source.group(1)
        faces.append(declarations)
    return faces


def _embedded_font(face: dict[str, str]) -> tuple[str, TTFont]:
    match = re.match(r'url\(data:(font/\w+);base64,([A-Za-z0-9+/=]+)\)', face["src"])
    assert match, face["src"]
    return match.group(1), TTFont(io.BytesIO(base64.b64decode(match.group(2))))


@pytest.fixture(autouse=True)
def isolated_fonts(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(fonts, "_fontconfig_cache", {})
    monkeypatch.setattr(fonts, "_platform_index", {}, raising=False)
    monkeypatch.setattr(diagram_module, "_warned_font_families", set(), raising=False)
    fonts.clear_registered_fonts()
    clear_face_urls = getattr(getattr(font_embedding, "_face_data_url_for", None), "cache_clear", None)
    if clear_face_urls is not None:
        clear_face_urls()
    yield
    fonts.clear_registered_fonts()


# ---- the embedding owner ----


# ---- SVGs that need no embedding are not parsed (L2-04) ----

_BUNDLED_ONLY = (
    ("'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif", "normal", "Hi"),
    ("serif", "bold", "Bold"),
    ("'Courier New', monospace", "normal", "Mono"),
)


def _refuse_parse(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*_args: object, **_kwargs: object) -> ET.Element:
        raise AssertionError("the SVG was parsed")

    monkeypatch.setattr(font_embedding.ET, "fromstring", refuse)


def test_an_svg_drawn_only_in_bundled_faces_is_not_parsed(monkeypatch: pytest.MonkeyPatch) -> None:
    svg = _svg(*_BUNDLED_ONLY)
    _refuse_parse(monkeypatch)

    assert font_embedding.svg_font_embedding(svg) == font_embedding.FontEmbedding(css="", refused_families=())
    assert font_embedding.fonts_for_cairosvg(svg.encode("utf-8"), "pdf") == svg.encode("utf-8")


def _custom_font_svgs(family: str) -> dict[str, str]:
    """SVGs that draw text in ``family`` through each spelling of a font-family value."""
    encoded = "&#" + str(ord(family[0])) + ";" + family[1:]
    text = _svg(("'Liberation Sans', sans-serif", "normal", "plain"))
    return {
        "attribute": _svg((f"'{family}', serif", "normal", "Hello")),
        "after an unknown family": _svg((f"'Gbdraw Not Installed', '{family}', serif", "normal", "Hello")),
        "bold only": _svg((f"'{family}', serif", "bold", "Hello")),
        "style declaration": text.replace("<text ", f"<text style=\"font-size: 9px; font-family: '{family}', serif\" ", 1)
        .replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1),
        "single-quoted attribute": text.replace("<text ", f"<text font-family='\"{family}\", serif' ", 1).replace(
            'font-family="\'Liberation Sans\', sans-serif" ', "", 1
        ),
        "character reference": text.replace("<text ", f'<text font-family="\'{encoded}\', serif" ', 1).replace(
            'font-family="\'Liberation Sans\', sans-serif" ', "", 1
        ),
        "style name as a character reference": text.replace(
            "<text ", f"<text style=\"font&#45;family: '{family}', serif\" ", 1
        ).replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1),
        "style colon as a character reference": text.replace(
            "<text ", f"<text style=\"font-family&#58; '{family}', serif\" ", 1
        ).replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1),
        "DTD attribute default": text.replace(
            "<svg ", f"<!DOCTYPE svg [<!ATTLIST text font-family CDATA \"'{family}'\">]><svg ", 1
        ).replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1),
        "entity in a style value": text.replace("<text ", "<text style=\"font-family: &fam;\" ", 1).replace(
            "<svg ", f"<!DOCTYPE svg [<!ENTITY fam \"'{family}'\">]><svg ", 1
        ).replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1),
        "inherited from a group": text.replace("<text ", f"<g font-family=\"'{family}'\"><text ", 1)
        .replace('font-family="\'Liberation Sans\', sans-serif" ', "", 1)
        .replace("</text>", "</text></g>", 1),
    }


@pytest.mark.parametrize("source", ["registered", "fontconfig"])
def test_an_svg_that_may_draw_a_custom_face_is_embedded_as_the_parse_finds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    family = "Gbdraw Embed Prefilter"
    path = _renamed_font("LiberationSerif-Regular", family, tmp_path / "p.ttf")
    if source == "registered":
        fonts.register_font_file(path)
    else:
        real_fontconfig_face = fonts._fontconfig_face

        def installed(name: str, weight: int, italic: bool) -> tuple[fonts.FontFace | None, str]:
            if fonts._comparable_name(name) == fonts._comparable_name(family):
                return fonts.FontFace(str(path)), "answered"
            return real_fontconfig_face(name, weight, italic)

        monkeypatch.setattr(fonts, "_fontconfig_face", installed)
        fonts._clear_resolution_caches()

    for spelling, svg in _custom_font_svgs(family).items():
        assert font_embedding._may_use_custom_face(svg), spelling
        assert font_embedding._may_use_custom_face(svg.encode("utf-8")), spelling
        embedding = font_embedding.svg_font_embedding(svg)
        assert [face["font-family"] for face in _faces(embedding.css)] == [f"'{family}'"], spelling
    for svg in (_svg(*_BUNDLED_ONLY), _svg(("'Gbdraw Not Installed', serif", "normal", "x"))):
        assert not font_embedding._may_use_custom_face(svg)


def test_a_registered_truetype_font_is_embedded_as_a_subset_of_its_used_characters(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Embed Serif", tmp_path / "s.ttf"))

    embedding = font_embedding.svg_font_embedding(_svg(("'Gbdraw Embed Serif', serif", "normal", "Hello you")))

    assert embedding.refused_families == ()
    (face,) = _faces(embedding.css)
    assert face["font-family"] == "'Gbdraw Embed Serif'"
    assert (face["font-weight"], face["font-style"]) == ("400", "normal")
    mime, font = _embedded_font(face)
    assert mime == "font/ttf" and "glyf" in font
    assert set(font.getBestCmap()) == {ord(char) for char in "Hello you"}
    assert len(font.getGlyphOrder()) == len(set("Hello you")) + 1  # and .notdef


def test_a_second_export_of_the_same_characters_subsets_the_face_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # OV-1568: SVG and PNG exports of one Result embed the same subset.
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Embed Twice", tmp_path / "t.ttf"))
    subsets: list[set[str]] = []
    real = font_embedding._subset

    def counting(font: TTFont, characters: set[str]) -> None:
        subsets.append(set(characters))
        real(font, characters)

    monkeypatch.setattr(font_embedding, "_subset", counting)
    svg = _svg(("'Gbdraw Embed Twice', serif", "normal", "Hello"))
    first = font_embedding.svg_font_embedding(svg)
    assert font_embedding.svg_font_embedding(svg) == first
    assert font_embedding.svg_font_embedding(svg, kind="raster") == first
    assert subsets == [set("Hello")]
    # Other characters are another subset.
    font_embedding.svg_font_embedding(_svg(("'Gbdraw Embed Twice', serif", "normal", "Help")))
    assert subsets == [set("Hello"), set("Help")]


def test_an_opentype_cff_font_keeps_its_flavour(tmp_path: Path) -> None:
    fonts.register_font_file(_cff_font(tmp_path / "c.otf", "Gbdraw Embed Cff"))

    embedding = font_embedding.svg_font_embedding(_svg(("'Gbdraw Embed Cff', sans-serif", "normal", "OV")))

    (face,) = _faces(embedding.css)
    mime, font = _embedded_font(face)
    assert mime == "font/otf" and "CFF " in font and "glyf" not in font
    assert set(font.getBestCmap()) == {ord("O"), ord("V")}


def test_each_used_face_is_embedded_once_with_its_own_weight(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSans-Regular", "Gbdraw Two Faces", tmp_path / "r.ttf"))
    fonts.register_font_file(
        _renamed_font("LiberationSans-Bold", "Gbdraw Two Faces", tmp_path / "b.ttf", weight=700)
    )
    family = "'Gbdraw Two Faces', sans-serif"

    embedding = font_embedding.svg_font_embedding(
        _svg((family, "normal", "ab"), (family, "bold", "bc"), (family, "normal", "a"))
    )

    faces = sorted(_faces(embedding.css), key=lambda face: face["font-weight"])
    assert [face["font-weight"] for face in faces] == ["400", "700"]
    assert set(_embedded_font(faces[0])[1].getBestCmap()) == {ord("a"), ord("b")}
    assert set(_embedded_font(faces[1])[1].getBestCmap()) == {ord("b"), ord("c")}


def test_bundled_fonts_are_never_embedded() -> None:
    svg = _svg(
        ("'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif", "bold", "Bundled"),
        ("'Times New Roman', serif", "normal", "Alias"),
        ("monospace", "normal", "Generic"),
    )

    embedding = font_embedding.svg_font_embedding(svg)

    assert embedding == font_embedding.FontEmbedding(css="", refused_families=())
    assert font_embedding.embed_font_css(svg, embedding) == svg


def test_a_font_whose_license_forbids_embedding_is_refused_by_family(tmp_path: Path) -> None:
    fonts.register_font_file(
        _renamed_font("LiberationSans-Regular", "Gbdraw Restricted", tmp_path / "x.ttf", fs_type=0x0002)
    )
    fonts.register_font_file(
        _renamed_font("LiberationSerif-Regular", "Gbdraw Printable", tmp_path / "p.ttf", fs_type=0x0002 | 0x0004)
    )

    embedding = font_embedding.svg_font_embedding(
        _svg(("'Gbdraw Restricted', sans-serif", "normal", "No"), ("'Gbdraw Printable', serif", "normal", "Yes"))
    )

    assert embedding.refused_families == ("Gbdraw Restricted",)
    assert [face["font-family"] for face in _faces(embedding.css)] == ["'Gbdraw Printable'"]


def test_a_font_that_forbids_subsetting_is_embedded_whole(tmp_path: Path) -> None:
    # D-09: fsType bit 0x0100 (no subsetting) with embedding allowed.
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Whole Serif", tmp_path / "w.ttf", fs_type=0x0100)
    fonts.register_font_file(path)

    embedding = font_embedding.svg_font_embedding(_svg(("'Gbdraw Whole Serif', serif", "normal", "Hi")))

    assert embedding.refused_families == ()
    (face,) = _faces(embedding.css)
    _mime, font = _embedded_font(face)
    assert font.getGlyphOrder() == TTFont(str(path)).getGlyphOrder()
    assert set(font.getBestCmap()) == set(TTFont(str(path)).getBestCmap())


def test_the_css_goes_into_a_defs_of_its_own_and_comes_out_again(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw <Odd> & Co", tmp_path / "o.ttf"))
    svg = _svg(("\"Gbdraw <Odd> & Co\", serif", "normal", "Text"))

    embedded = font_embedding.embed_font_css(svg, font_embedding.svg_font_embedding(svg))

    root = ET.fromstring(embedded)
    styles = root.findall("{http://www.w3.org/2000/svg}defs/{http://www.w3.org/2000/svg}style")
    assert len(styles) == 1 and "@font-face" in (styles[0].text or "")
    assert "Gbdraw <Odd> & Co" in (styles[0].text or "")
    assert font_embedding.remove_embedded_font_css(embedded) == svg


def test_embedding_is_deterministic(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Same Bytes", tmp_path / "d.ttf"))
    svg = _svg(("'Gbdraw Same Bytes', serif", "normal", "Twice"))

    assert font_embedding.svg_font_embedding(svg) == font_embedding.svg_font_embedding(svg)


def test_the_api_writers_embed_by_default_and_can_turn_it_off(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Api Serif", tmp_path / "a.ttf"))
    drawing = svgwrite.Drawing(size=("200px", "100px"))
    drawing.add(drawing.text("Api", insert=(0, 20), font_family="'Gbdraw Api Serif', serif"))

    embedded = render_to_bytes(drawing, "svg").decode("utf-8")
    plain = render_to_bytes(drawing, "svg", embed_fonts=False).decode("utf-8")

    assert embedded.count("@font-face") == 1
    assert plain == drawing.tostring()
    assert font_embedding.remove_embedded_font_css(embedded) == plain


# ---- PDF ----


def test_a_cff_font_converts_to_truetype_outlines_with_the_same_advances(tmp_path: Path) -> None:
    path = _cff_font(tmp_path / "c.otf", "Gbdraw Pdf Cff")
    registered = fonts.register_font_file(path)

    converted = TTFont(io.BytesIO(font_embedding.truetype_bytes_for_pdf(registered.face)))

    source = TTFont(str(path))
    assert "glyf" in converted and "CFF " not in converted
    assert converted.sfntVersion == "\x00\x01\x00\x00"
    assert {name: converted["hmtx"][name][0] for name in source.getGlyphOrder()} == {
        name: source["hmtx"][name][0] for name in source.getGlyphOrder()
    }
    assert converted["glyf"]["O"].numberOfContours == 1
    assert font_embedding.truetype_bytes_for_pdf(path) == font_embedding.truetype_bytes_for_pdf(registered.face)


def test_a_truetype_font_comes_back_as_its_own_outlines(tmp_path: Path) -> None:
    path = _renamed_font("LiberationSans-Regular", "Gbdraw Pdf Ttf", tmp_path / "t.ttf")

    converted = TTFont(io.BytesIO(font_embedding.truetype_bytes_for_pdf(path)))

    assert "glyf" in converted and len(converted.getGlyphOrder()) == len(TTFont(str(path)).getGlyphOrder())


@pytest.mark.parametrize("kind", ["cff", "ttf"])
def test_a_pdf_face_without_a_postscript_name_gets_one(tmp_path: Path, kind: str) -> None:
    # jsPDF reads name ID 6 (the PostScript name) of every font it embeds and
    # fails the export without one; _cff_font has none, like some real fonts.
    if kind == "cff":
        path = _cff_font(tmp_path / "c.otf", "Gbdraw Pdf No Ps")
    else:
        path = _renamed_font("LiberationSans-Regular", "Gbdraw Pdf No Ps", tmp_path / "t.ttf")
        font = TTFont(str(path))
        font["name"].names = [record for record in font["name"].names if record.nameID != 6]
        font.save(str(path))

    converted = TTFont(io.BytesIO(font_embedding.truetype_bytes_for_pdf(path)))

    assert converted["name"].getDebugName(6) == "GbdrawPdfNoPs-Regular"


def test_a_font_that_cannot_convert_is_a_font_file_error(tmp_path: Path) -> None:
    broken = tmp_path / "broken.otf"
    broken.write_bytes(b"OTTO" + b"\x00" * 64)

    with pytest.raises(ValidationError) as caught:
        font_embedding.truetype_bytes_for_pdf(broken)

    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FONT_FILE"}


def _pdf_faces(choice: font_embedding.PdfFonts) -> list[tuple[str, str]]:
    return [(face.family, face.file) for face in choice.faces]


def test_pdf_fonts_name_the_bundled_face_of_each_run_as_jspdf_registers_it() -> None:
    choice = font_embedding.pdf_fonts(
        [
            ("'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif", True, False, "Sans"),
            ("serif", False, True, "Serif"),
            ("'Courier New', monospace", True, True, "Mono"),
        ],
        known=("LiberationSerif-Italic.ttf",),
    )

    assert _pdf_faces(choice) == [
        ("LiberationSans", "LiberationSans-Bold.ttf"),
        ("LiberationSerif", "LiberationSerif-Italic.ttf"),
        ("LiberationMono", "LiberationMono-BoldItalic.ttf"),
    ]
    # Only the files the caller does not hold yet, byte for byte the bundled files.
    assert set(choice.fonts) == {"LiberationSans-Bold.ttf", "LiberationMono-BoldItalic.ttf"}
    assert choice.fonts["LiberationSans-Bold.ttf"] == (DATA_DIR / "LiberationSans-Bold.ttf").read_bytes()
    assert choice.substitutions == ()


def test_pdf_fonts_draw_a_registered_cff_face_with_truetype_outlines_sent_once(tmp_path: Path) -> None:
    registered = fonts.register_font_file(_cff_font(tmp_path / "c.otf", "Gbdraw Pdf Run"))
    name = registered.face.measurement_name
    regular = ("'Gbdraw Pdf Run', sans-serif", False, False, "OV")

    choice = font_embedding.pdf_fonts([regular, ("'Gbdraw Pdf Run', sans-serif", True, False, "V")], known=())

    # The one face draws both runs (CSS matching), under one jsPDF family.
    (face,) = set(choice.faces)
    assert _pdf_faces(choice) == [(name.rsplit("-", 1)[0], face.file)] * 2
    assert face.file.startswith(f"{name}-") and face.file.endswith(".ttf")
    assert choice.fonts == {face.file: font_embedding.truetype_bytes_for_pdf(registered.face, "OV")}
    assert "glyf" in TTFont(io.BytesIO(choice.fonts[face.file]))
    assert font_embedding.pdf_fonts([regular], known=(face.file,)).fonts == {}


def test_pdf_fonts_substitute_the_nearest_liberation_face_and_name_the_font(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fonts.register_font_file(
        _renamed_font("LiberationSans-Regular", "Gbdraw Pdf Restricted", tmp_path / "x.ttf", fs_type=0x0002)
    )
    broken = fonts.register_font_file(_renamed_font("LiberationMono-Regular", "Gbdraw Pdf Broken", tmp_path / "b.ttf"))
    real = font_embedding._truetype_bytes

    def failing(face: fonts.FontFace, characters: frozenset[str] | None) -> bytes:
        if face == broken.face:
            raise ValidationError("cannot convert", diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"})
        return real(face, characters)

    monkeypatch.setattr(font_embedding, "_truetype_bytes", failing)

    choice = font_embedding.pdf_fonts(
        [("'Gbdraw Pdf Restricted', serif", False, False, "R"), ("'Gbdraw Pdf Broken', monospace", True, False, "B")],
        known=(),
    )

    assert _pdf_faces(choice) == [
        ("LiberationSerif", "LiberationSerif-Regular.ttf"),
        ("LiberationMono", "LiberationMono-Bold.ttf"),
    ]
    assert choice.substitutions == (
        font_embedding.PdfFontSubstitution("Gbdraw Pdf Restricted", "Liberation Serif", "license"),
        font_embedding.PdfFontSubstitution("Gbdraw Pdf Broken", "Liberation Mono", "font-file"),
    )


def test_the_web_export_helpers_reply_in_json(tmp_path: Path) -> None:
    registered = fonts.register_font_file(
        _renamed_font("LiberationSerif-Regular", "Gbdraw Web Export", tmp_path / "w.ttf")
    )
    svg = _svg(("'Gbdraw Web Export', serif", "normal", "Web"))

    reply = json.loads(export_fonts.svg_font_embedding_json(svg))
    assert reply == {
        "block": font_embedding.embedded_fonts_block(font_embedding.svg_font_embedding(svg)), "refusedFamilies": []
    }

    name = registered.face.measurement_name
    pdf = json.loads(export_fonts.pdf_fonts_json(
        json.dumps([{"fontFamily": "'Gbdraw Web Export', serif", "bold": False, "italic": False, "text": "Web"}]), "[]"
    ))
    [face] = pdf["faces"]
    assert face["family"] == name.rsplit("-", 1)[0] and face["file"].startswith(f"{name}-")
    assert base64.b64decode(pdf["fonts"][face["file"]]) == font_embedding.truetype_bytes_for_pdf(registered.face, "Web")
    assert pdf["substitutions"] == []


# ---- CLI ----


def _box_font(target: Path, family: str) -> Path:
    """Liberation Serif with every glyph drawn as its bounding box: the same
    metrics, unmistakable ink."""
    from fontTools.pens.boundsPen import BoundsPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen

    path = _renamed_font("LiberationSerif-Regular", family, target)
    font = TTFont(str(path))
    glyphs = font.getGlyphSet()
    boxes = {}
    for name in font.getGlyphOrder():
        bounds = BoundsPen(glyphs)
        glyphs[name].draw(bounds)
        boxes[name] = bounds.bounds
    for name, box in boxes.items():
        pen = TTGlyphPen(None)
        if box:
            x0, y0, x1, y1 = box
            pen.moveTo((x0, y0))
            pen.lineTo((x0, y1))
            pen.lineTo((x1, y1))
            pen.lineTo((x1, y0))
            pen.closePath()
        font["glyf"][name] = pen.glyph()
    font.save(str(path))
    return path


def test_a_font_file_that_is_not_installed_draws_in_cairo_raster_output(tmp_path: Path) -> None:
    # Receipt 〔書き出し〕: PNG and PDF draw the chosen font. CairoSVG finds
    # fonts through fontconfig, so a registered file must be visible there.
    pytest.importorskip("cairosvg")
    from PIL import Image

    from gbdraw.render import export as export_module

    # D-13: a PNG draws the chosen font whatever its license.
    path = _box_font(tmp_path / "box.ttf", "Gbdraw Raster Box")
    font = TTFont(str(path))
    font["OS/2"].fsType = 0x0002
    font.save(str(path))
    fonts.register_font_file_faces(path)

    def ink(font_family: str) -> int:
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80">'
            '<rect width="400" height="80" fill="white"/>'
            f'<text x="10" y="50" font-size="40" font-family="{font_family}">Hello gbdraw</text></svg>'
        )
        image = Image.open(io.BytesIO(export_module.convert_svg_with_cairosvg(svg, "png"))).convert("L")
        return sum(1 for value in image.getdata() if value < 64)

    assert ink("'Gbdraw Raster Box', serif") > 1.5 * ink("'Liberation Serif', serif")




def _font_faces_in(svg_path: Path) -> list[dict[str, str]]:
    return _faces(svg_path.read_text(encoding="utf-8"))


@pytest.mark.parametrize("main", [circular_main, linear_main], ids=["circular", "linear"])
def test_cli_font_file_draws_with_its_family_and_embeds_it(tmp_path: Path, main) -> None:
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw File Serif", tmp_path / "f.ttf")
    gbk = _genbank(tmp_path / "in.gb")

    main(["--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg,interactive_svg", "--font_file", str(font)])

    svg = tmp_path / "out.svg"
    assert set(re.findall(r'font-family="([^"]*)"', svg.read_text())) == {"'Gbdraw File Serif', serif"}
    faces = _font_faces_in(svg)
    assert [face["font-family"] for face in faces] == ["'Gbdraw File Serif'"]
    ET.parse(svg)
    interactive = tmp_path / "out.interactive.svg"
    assert [face["font-family"] for face in _font_faces_in(interactive)] == ["'Gbdraw File Serif'"]
    ET.parse(interactive)


def test_cli_no_font_embedding_writes_the_svg_without_fonts(tmp_path: Path) -> None:
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw Off Serif", tmp_path / "f.ttf")
    gbk = _genbank(tmp_path / "in.gb")
    base = ["--gbk", str(gbk), "-f", "svg", "--font_file", str(font)]

    circular_main([*base, "-o", str(tmp_path / "on")])
    circular_main([*base, "-o", str(tmp_path / "off"), "--no_font_embedding"])

    off = (tmp_path / "off.svg").read_text(encoding="utf-8")
    on = (tmp_path / "on.svg").read_text(encoding="utf-8")
    assert "@font-face" not in off
    assert font_embedding.remove_embedded_font_css(on) == off


def test_cli_logs_the_families_it_may_not_embed(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    font = _renamed_font("LiberationSans-Regular", "Gbdraw Locked Sans", tmp_path / "l.ttf", fs_type=0x0002)
    gbk = _genbank(tmp_path / "in.gb")

    with caplog.at_level(logging.WARNING):
        circular_main(["--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg", "--font_file", str(font)])

    assert "Gbdraw Locked Sans" in caplog.text and "not embedded" in caplog.text
    assert "@font-face" not in (tmp_path / "out.svg").read_text()


def test_cli_registers_every_face_of_a_collection(tmp_path: Path) -> None:
    regular = _renamed_font("LiberationSans-Regular", "Gbdraw Collection", tmp_path / "r.ttf")
    bold = _renamed_font("LiberationSans-Bold", "Gbdraw Collection", tmp_path / "b.ttf", weight=700)
    collection = TTCollection()
    collection.fonts = [TTFont(str(regular)), TTFont(str(bold))]
    ttc = tmp_path / "both.ttc"
    collection.save(str(ttc))
    gbk = _genbank(tmp_path / "in.gb")

    circular_main(["--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg", "--font_file", str(ttc)])

    registered = sorted((entry.family, entry.weight, entry.face.index) for entry in fonts.registered_fonts())
    assert registered == [("Gbdraw Collection", 400, 0), ("Gbdraw Collection", 700, 1)]
    assert "'Gbdraw Collection', sans-serif" in (tmp_path / "out.svg").read_text()


def test_cli_font_family_names_one_of_several_font_files(tmp_path: Path) -> None:
    first = _renamed_font("LiberationSans-Regular", "Gbdraw First", tmp_path / "1.ttf")
    second = _renamed_font("LiberationMono-Regular", "Gbdraw Second", tmp_path / "2.ttf")
    gbk = _genbank(tmp_path / "in.gb")

    circular_main([
        "--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg",
        "--font_file", str(first), "--font_file", str(second), "--font_family", "Gbdraw Second",
    ])

    assert set(re.findall(r'font-family="([^"]*)"', (tmp_path / "out.svg").read_text())) == {
        "'Gbdraw Second', monospace"
    }


def test_cli_missing_font_file_is_a_font_file_error(tmp_path: Path) -> None:
    gbk = _genbank(tmp_path / "in.gb")

    with pytest.raises(ValidationError) as caught:
        circular_main([
            "--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg",
            "--font_file", str(tmp_path / "missing.ttf"),
        ])

    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FONT_FILE"}
    assert not (tmp_path / "out.svg").exists()


def test_cli_session_keeps_the_family_name_and_a_font_file_draws_it_again(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw Session Serif", tmp_path / "s.ttf")
    gbk = _genbank(tmp_path / "in.gb")
    session_path = tmp_path / "font.gbdraw-session.json"

    with caplog.at_level(logging.WARNING):
        circular_main([
            "--gbk", str(gbk), "-o", str(tmp_path / "font"), "-f", "svg",
            "--font_file", str(font), "--no_font_embedding",
            "--save_session", "--session_output", str(session_path),
        ])
    saved = session_path.read_text(encoding="utf-8")

    # D-02: the name only; no font file, its path, or embedded font data.
    assert "The Session stores the font name 'Gbdraw Session Serif' but not the font file" in caplog.text
    assert "'Gbdraw Session Serif', serif" in json.dumps(session_view(json.loads(saved))["renderRequest"])
    assert str(font) not in saved and "--font_file" not in saved and "--no_font_embedding" not in saved
    assert "@font-face" not in saved

    fonts.clear_registered_fonts()
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        circular_main(["--session", str(session_path), "-o", str(tmp_path / "substituted"), "-f", "svg"])
    assert "Font 'Gbdraw Session Serif' is not available" in caplog.text

    circular_main([
        "--session", str(session_path), "-o", str(tmp_path / "replay"), "-f", "svg",
        "--font_file", str(font),
    ])
    replay = (tmp_path / "replay.svg").read_text(encoding="utf-8")
    assert [face["font-family"] for face in _faces(replay)] == ["'Gbdraw Session Serif'"]
    assert font_embedding.remove_embedded_font_css(replay) == (tmp_path / "font.svg").read_text(encoding="utf-8")

    circular_main([
        "--session", str(session_path), "-o", str(tmp_path / "plain"), "-f", "svg",
        "--font_file", str(font), "--no_font_embedding",
    ])
    assert (tmp_path / "plain.svg").read_text(encoding="utf-8") == (tmp_path / "font.svg").read_text(encoding="utf-8")


def test_a_saved_session_keeps_no_embedded_font_in_its_preview(tmp_path: Path) -> None:
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw Preview Serif", tmp_path / "s.ttf")
    gbk = _genbank(tmp_path / "in.gb")
    session_path = tmp_path / "embedded.gbdraw-session.json"

    circular_main([
        "--gbk", str(gbk), "-o", str(tmp_path / "embedded"), "-f", "svg",
        "--font_file", str(font), "--save_session", "--session_output", str(session_path),
    ])

    assert "@font-face" in (tmp_path / "embedded.svg").read_text(encoding="utf-8")
    assert "@font-face" not in session_path.read_text(encoding="utf-8")


def test_ov688_a_session_saved_again_warns_that_it_keeps_only_the_font_name(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw Resave Serif", tmp_path / "s.ttf")
    gbk = _genbank(tmp_path / "in.gb")
    first = tmp_path / "first.gbdraw-session.json"
    circular_main([
        "--gbk", str(gbk), "-o", str(tmp_path / "first"), "-f", "svg",
        "--font_file", str(font), "--save_session", "--session_output", str(first),
    ])
    caplog.clear()

    with caplog.at_level(logging.WARNING):
        circular_main([
            "--session", str(first), "-o", str(tmp_path / "second"), "-f", "svg", "--font_file", str(font),
            "--session_output", str(tmp_path / "second.gbdraw-session.json"),
        ])

    assert "The Session stores the font name 'Gbdraw Resave Serif' but not the font file" in caplog.text
    assert "@font-face" not in (tmp_path / "second.gbdraw-session.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("kind, reason", [("svg", "FONT_EMBEDDING"), ("raster", "FONT_RASTER")])
def test_a_web_embedding_failure_names_the_opt_out_only_for_svg(
    monkeypatch: pytest.MonkeyPatch, kind: str, reason: str
) -> None:
    # F2: the opt-out does not change PNG, so a PNG failure has its own guidance.
    def broken(_svg: str, _kind: str) -> font_embedding.FontEmbedding:
        raise RuntimeError("subsetter failed")

    monkeypatch.setattr(export_fonts, "svg_font_embedding", broken)

    with pytest.raises(ValidationError) as caught:
        export_fonts.svg_font_embedding_json("<svg/>", kind)

    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "reason": reason}


def test_f2_a_raster_embeds_the_whole_face_when_subsetting_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The raster SVG never leaves the page, so the whole face may go in it.
    fonts.register_font_file(_renamed_font("LiberationSans-Regular", "Gbdraw Raster Whole", tmp_path / "w.ttf"))
    svg = _svg(("'Gbdraw Raster Whole', sans-serif", "normal", "Ab"))

    def broken(_font: TTFont, _characters: set[str]) -> None:
        raise RuntimeError("subsetter failed")

    monkeypatch.setattr(font_embedding, "_subset", broken)

    (face,) = _faces(font_embedding.svg_font_embedding(svg, kind="raster").css)
    assert {ord("A"), ord("b"), ord("z")} <= set(_embedded_font(face)[1].getBestCmap())
    with pytest.raises(RuntimeError, match="subsetter failed"):
        font_embedding.svg_font_embedding(svg, kind="svg")


def test_gbdraw_render_takes_font_files_and_no_font_embedding(tmp_path: Path) -> None:
    # D-09: `gbdraw render --session` replays a Session that names a custom font.
    font = _renamed_font("LiberationSerif-Regular", "Gbdraw Render Serif", tmp_path / "s.ttf")
    gbk = _genbank(tmp_path / "in.gb")
    session_path = tmp_path / "font.gbdraw-session.json"
    circular_main([
        "--gbk", str(gbk), "-o", str(tmp_path / "font"), "-f", "svg",
        "--font_file", str(font), "--no_font_embedding",
        "--save_session", "--session_output", str(session_path),
    ])
    fonts.clear_registered_fonts()

    render_main(["--session", str(session_path), "-o", str(tmp_path / "embedded"), "-f", "svg", "--font_file", str(font)])
    render_main([
        "--session", str(session_path), "-o", str(tmp_path / "plain"), "-f", "svg",
        "--font_file", str(font), "--no_font_embedding",
    ])

    embedded = (tmp_path / "embedded.svg").read_text(encoding="utf-8")
    assert [face["font-family"] for face in _faces(embedded)] == ["'Gbdraw Render Serif'"]
    plain = (tmp_path / "plain.svg").read_text(encoding="utf-8")
    assert plain == (tmp_path / "font.svg").read_text(encoding="utf-8")
    assert font_embedding.remove_embedded_font_css(embedded) == plain


def test_a_font_file_cairo_cannot_see_is_named_in_a_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    import ctypes
    import ctypes.util

    def missing(_name: str) -> None:
        raise OSError("no fontconfig")

    pytest.importorskip("cairosvg")
    from gbdraw.render import export as export_module

    monkeypatch.setattr(ctypes, "CDLL", missing)
    monkeypatch.setattr(ctypes.util, "find_library", lambda _name: None)
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw No Fc Serif", tmp_path / "n.ttf")
    fonts.register_font_file_faces(path)
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="40">'
        "<text x=\"5\" y=\"30\" font-family=\"'Gbdraw No Fc Serif', serif\">Fc</text></svg>"
    )

    # M2: at the CairoSVG export, once per file.
    with caplog.at_level(logging.WARNING, logger="gbdraw.core.fonts"):
        export_module.convert_svg_with_cairosvg(svg, "png")
        export_module.convert_svg_with_cairosvg(svg, "pdf")

    message = "PNG, PDF, EPS and PS output draw 'Gbdraw No Fc Serif' with an installed font"
    assert caplog.text.count(message) == 1


# ---- PR-2 review fixes (D-13): one owner decides per output kind ----

# fsType values and what each output kind may do with the face (D-13).
_FS_TYPE_CASES = {
    "installable": (0x0000, "subset", "subset", "subset"),
    "restricted": (0x0002, "refused", "refused", "subset"),
    "restricted-printable": (0x0002 | 0x0004, "subset", "subset", "subset"),
    "bitmap-only": (0x0200, "refused", "refused", "subset"),
    "no-subsetting": (0x0100, "whole", "refused", "subset"),
}


@pytest.mark.parametrize("case", sorted(_FS_TYPE_CASES))
def test_one_owner_decides_how_each_output_kind_carries_a_face(tmp_path: Path, case: str) -> None:
    fs_type, svg, document, raster = _FS_TYPE_CASES[case]
    font = TTFont(str(_renamed_font("LiberationSans-Regular", "Gbdraw Kind", tmp_path / "k.ttf", fs_type=fs_type)))

    assert [font_embedding.face_embedding(font, kind) for kind in ("svg", "document", "raster")] == [
        svg, document, raster
    ]


def test_s1_a_bitmap_only_font_is_refused_in_the_svg(tmp_path: Path) -> None:
    fonts.register_font_file(
        _renamed_font("LiberationSans-Regular", "Gbdraw Bitmap Only", tmp_path / "m.ttf", fs_type=0x0200)
    )

    embedding = font_embedding.svg_font_embedding(_svg(("'Gbdraw Bitmap Only', sans-serif", "normal", "No")))

    assert embedding == font_embedding.FontEmbedding(css="", refused_families=("Gbdraw Bitmap Only",))


def test_d13_a_raster_carries_every_custom_face_whatever_its_license(tmp_path: Path) -> None:
    fonts.register_font_file(
        _renamed_font("LiberationSans-Regular", "Gbdraw Raster Locked", tmp_path / "r.ttf", fs_type=0x0002)
    )

    embedding = font_embedding.svg_font_embedding(
        _svg(("'Gbdraw Raster Locked', sans-serif", "normal", "Ab")), kind="raster"
    )

    assert embedding.refused_families == ()
    (face,) = _faces(embedding.css)
    assert set(_embedded_font(face)[1].getBestCmap()) == {ord("A"), ord("b")}


@pytest.mark.parametrize("fs_type", [0x0002, 0x0200, 0x0100], ids=["restricted", "bitmap-only", "no-subsetting"])
def test_s1_s2_a_pdf_draws_a_face_it_may_not_carry_in_the_nearest_liberation_face(
    tmp_path: Path, fs_type: int
) -> None:
    fonts.register_font_file(
        _renamed_font("LiberationSerif-Regular", "Gbdraw Pdf Refused", tmp_path / "x.ttf", fs_type=fs_type)
    )

    choice = font_embedding.pdf_fonts([("'Gbdraw Pdf Refused', serif", False, False, "Text")], known=())

    assert _pdf_faces(choice) == [("LiberationSerif", "LiberationSerif-Regular.ttf")]
    assert choice.substitutions == (
        font_embedding.PdfFontSubstitution("Gbdraw Pdf Refused", "Liberation Serif", "license"),
    )


def test_s3_a_pdf_carries_a_subset_of_the_characters_each_custom_face_draws(tmp_path: Path) -> None:
    registered = fonts.register_font_file(_cff_font(tmp_path / "c.otf", "Gbdraw Pdf Subset"))
    family = "'Gbdraw Pdf Subset', sans-serif"

    choice = font_embedding.pdf_fonts([(family, False, False, "O"), (family, True, False, "V ")], known=())

    (file,) = set(face.file for face in choice.faces)
    font = TTFont(io.BytesIO(choice.fonts[file]))
    assert "glyf" in font and set(font.getBestCmap()) == {ord("O"), ord("V"), ord(" ")}
    assert file.startswith(registered.face.measurement_name)
    # Other characters are another file; the same characters are known.
    other = font_embedding.pdf_fonts([(family, False, False, "OOO")], known=(file,))
    assert other.faces[0].file != file and set(other.fonts) == {other.faces[0].file}
    assert font_embedding.pdf_fonts([(family, False, False, "VO ")], known=(file,)).fonts == {}


def test_s6_the_web_helper_returns_the_block_embed_font_css_inserts(tmp_path: Path) -> None:
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Web Block", tmp_path / "w.ttf"))
    svg = _svg(("'Gbdraw Web Block', serif", "normal", "Block"))

    reply = json.loads(export_fonts.svg_font_embedding_json(svg))

    embedding = font_embedding.svg_font_embedding(svg)
    block = font_embedding.embedded_fonts_block(embedding)
    assert reply == {"block": block, "refusedFamilies": []}
    assert block.startswith('<defs id="gbdraw-embedded-fonts">') and "@font-face" in block
    assert font_embedding.embed_font_css(svg, embedding).replace(block, "", 1) == svg
    assert json.loads(export_fonts.svg_font_embedding_json(svg, "raster"))["block"] == block


def test_m1_an_embedding_failure_names_the_cli_and_python_opt_outs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gbdraw.api.render as render_module
    from gbdraw.exceptions import ExportError

    def broken(_svg: str) -> font_embedding.FontEmbedding:
        raise RuntimeError("subsetter failed")

    monkeypatch.setattr(render_module, "svg_font_embedding", broken)
    drawing = svgwrite.Drawing(size=("10px", "10px"))

    with pytest.raises(ExportError) as caught:
        render_to_bytes(drawing, "svg")

    assert "--no_font_embedding" in str(caught.value) and "embed_fonts=False" in str(caught.value)


def _pdf_with(tmp_path: Path, family: str, fs_type: int, caplog: pytest.LogCaptureFixture) -> bytes:
    pytest.importorskip("cairosvg")
    from gbdraw.render import export as export_module

    fonts.register_font_file(
        _renamed_font("LiberationSerif-Regular", family, tmp_path / "f.ttf", fs_type=fs_type)
    )
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="60">'
        f'<text x="5" y="40" font-size="20" font-family="\'{family}\', serif">Hello</text></svg>'
    )
    with caplog.at_level(logging.WARNING):
        return export_module.convert_svg_with_cairosvg(svg, "pdf")


@pytest.mark.parametrize("fs_type", [0x0002, 0x0200, 0x0100], ids=["restricted", "bitmap-only", "no-subsetting"])
def test_b1_a_cairo_pdf_never_carries_a_face_whose_license_forbids_it(
    tmp_path: Path, fs_type: int, caplog: pytest.LogCaptureFixture
) -> None:
    family = f"Gbdraw Cairo Refused {fs_type:04x}"

    pdf = _pdf_with(tmp_path, family, fs_type, caplog)

    assert family.replace(" ", "").encode() not in pdf
    assert b"LiberationSerif" in pdf
    assert f"Font '{family}' is not embedded in the PDF" in caplog.text
    assert "Its text is drawn in 'Liberation Serif'" in caplog.text


def test_b1_a_cairo_pdf_carries_a_face_it_may_embed(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    family = "Gbdraw Cairo Allowed"

    pdf = _pdf_with(tmp_path, family, 0x0000, caplog)

    assert family.replace(" ", "").encode() in pdf
    assert "not embedded" not in caplog.text


def test_m2_registering_a_font_file_shows_nothing_to_fontconfig(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    import ctypes

    def refuse(_name: str) -> None:
        raise AssertionError("fontconfig is touched at registration")

    monkeypatch.setattr(ctypes, "CDLL", refuse)
    with caplog.at_level(logging.WARNING):
        fonts.register_font_file(_renamed_font("LiberationSerif-Regular", "Gbdraw Svg Only", tmp_path / "s.ttf"))

    assert caplog.text == ""


def test_ov744_a_font_with_a_reserved_unicode_range_bit_still_subsets(tmp_path: Path) -> None:
    # GNU Unifont sets OS/2 Unicode range bit 123, which fontTools cannot
    # recalculate when it prunes the ranges of a subset.
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Range Bit", tmp_path / "r.ttf")
    font = TTFont(str(path))
    font["OS/2"].ulUnicodeRange4 |= 1 << (123 - 96)
    font.save(str(path))
    registered = fonts.register_font_file(path)

    embedding = font_embedding.svg_font_embedding(_svg(("'Gbdraw Range Bit', serif", "normal", "Ok")))

    assert [face["font-family"] for face in _faces(embedding.css)] == ["'Gbdraw Range Bit'"]
    assert "glyf" in TTFont(io.BytesIO(font_embedding.truetype_bytes_for_pdf(registered.face, "Ok")))


def test_bytes_in_another_declared_encoding_are_read_as_declared(tmp_path: Path) -> None:
    family = "Gbdraw Caf\u00e9 Prefilter"
    fonts.register_font_file(_renamed_font("LiberationSerif-Regular", family, tmp_path / "e.ttf"))
    svg = _svg((f"'{family}', serif", "normal", "Hello"))
    declared = '<?xml version="1.0" encoding="ISO-8859-1"?>' + svg
    assert font_embedding._may_use_custom_face(declared.encode("latin-1"))
    assert font_embedding._may_use_custom_face(svg.encode("utf-16"))


def test_css_text_with_a_reference_does_not_force_the_parse() -> None:
    svg = _svg(*_BUNDLED_ONLY).replace("</svg>", "<style>a &gt; b { color: red }</style></svg>")
    assert not font_embedding._may_use_custom_face(svg)
    assert not font_embedding._may_use_custom_face(svg.encode("utf-8"))
