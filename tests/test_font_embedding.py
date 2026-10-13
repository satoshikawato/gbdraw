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
from gbdraw.core import fonts
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.render import font_embedding

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
    yield
    fonts.clear_registered_fonts()


# ---- the embedding owner ----


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


def test_a_font_that_cannot_convert_is_a_font_file_error(tmp_path: Path) -> None:
    broken = tmp_path / "broken.otf"
    broken.write_bytes(b"OTTO" + b"\x00" * 64)

    with pytest.raises(ValidationError) as caught:
        font_embedding.truetype_bytes_for_pdf(broken)

    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FONT_FILE"}


# ---- CLI ----


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
    assert "'Gbdraw Session Serif', serif" in json.dumps(json.loads(saved)["renderRequest"])
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
