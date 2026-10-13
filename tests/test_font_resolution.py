"""Font availability, resolution order, substitution, and the --font_family flag.

The owner is ``gbdraw.core.fonts`` (custom-fonts D-01, D-03): a family resolves
to a bundled Liberation face or alias, a registered font file, or the file
fontconfig finds for exactly that family; a drawing names only families that
resolve, and the nearest Liberation family replaces missing ones with a warning.
fontconfig is replaced by a fake ``fc-match`` so the tests do not depend on the
fonts of the machine.
"""

from __future__ import annotations

import base64
import json
import logging
import os
import re
from pathlib import Path

import pytest
from fontTools.ttLib import TTFont

from gbdraw.api.config import load_default_config
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main
from gbdraw.core import fonts
from gbdraw.core.text import (
    calculate_bbox_dimensions,
    calculate_svg_bbox_dimensions,
    font_metrics_table,
    get_text_bbox_size_pixels,
    measurement_font_file,
)
from gbdraw.session_io import session_to_cli_args

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
DEFAULT_FAMILY = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
TEXT = "Hypothetical protein"


def _bundled(stem: str) -> str:
    return str(DATA_DIR / f"{stem}.ttf")


def _renamed_font(source: str, family: str, target: Path, *, bold: bool = False) -> Path:
    font = TTFont(_bundled(source))
    name_table = font["name"]
    name_table.names = [record for record in name_table.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
    style = "Bold" if bold else "Regular"
    for name_id, value in ((1, family), (2, style), (4, f"{family} {style}"), (6, family.replace(" ", "")), (16, family)):
        name_table.setName(value, name_id, 3, 1, 0x409)
    font["OS/2"].usWeightClass = 700 if bold else 400
    font.save(str(target))
    return target


@pytest.fixture(autouse=True)
def isolated_fonts(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(fonts, "_fontconfig_cache", {})
    fonts.clear_registered_fonts()
    yield
    fonts.clear_registered_fonts()


@pytest.fixture
def fake_fc_match(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Install an ``fc-match`` that answers ``family<TAB>file<TAB>0`` per requested family."""

    answers: dict[str, str] = {}
    bin_dir = tmp_path / "fake-bin"
    bin_dir.mkdir()
    table = tmp_path / "fc-answers.json"
    script = bin_dir / "fc-match"
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import json, sys\n"
        f"answers = json.load(open({str(table)!r}))\n"
        "requested = sys.argv[-1].split(':', 1)[0].replace('\\\\', '')\n"
        "family, path = answers.get(requested, answers.get('*', ['DejaVu Sans', '/nonexistent/DejaVuSans.ttf']))\n"
        "sys.stdout.write(f'{family}\\t{path}\\t0')\n"
    )
    script.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")

    def answer(requested: str, family: str, path: Path | str) -> None:
        answers[requested] = [family, str(path)]
        table.write_text(json.dumps(answers))

    table.write_text("{}")
    return answer


# ---- registry and measurement ----


def test_a_registered_font_is_measured_with_its_own_file(tmp_path: Path) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "a.ttf")

    registered = fonts.register_font_file(path)

    assert (registered.family, registered.weight, registered.italic) == ("Gbdraw Test Serif", 400, False)
    assert registered.generic == "serif"
    family = "'Gbdraw Test Serif', sans-serif"
    assert calculate_bbox_dimensions(TEXT, family, 14, 72) == get_text_bbox_size_pixels(path, TEXT, 14, 72)
    assert calculate_svg_bbox_dimensions(TEXT, family, 14, 96) == get_text_bbox_size_pixels(
        path, TEXT, 14, 96, svg_units=True
    )
    assert calculate_bbox_dimensions(TEXT, family, 14, 72) != calculate_bbox_dimensions(TEXT, "sans-serif", 14, 72)


def test_replacing_the_bytes_of_a_family_changes_its_measurement(tmp_path: Path) -> None:
    """The Web keeps one Worker: the same family can later name other bytes."""

    path = tmp_path / "user-font.ttf"
    _renamed_font("LiberationSerif-Regular", "Gbdraw Test Fam", path)
    fonts.register_font_file(path)
    serif_width = calculate_bbox_dimensions(TEXT, "'Gbdraw Test Fam'", 14, 72)
    serif_name = measurement_font_file("'Gbdraw Test Fam'")

    _renamed_font("LiberationMono-Regular", "Gbdraw Test Fam", path)
    fonts.register_font_file(path)

    assert calculate_bbox_dimensions(TEXT, "'Gbdraw Test Fam'", 14, 72) == get_text_bbox_size_pixels(
        _bundled("LiberationMono-Regular"), TEXT, 14, 72
    )
    assert calculate_bbox_dimensions(TEXT, "'Gbdraw Test Fam'", 14, 72) != serif_width
    assert measurement_font_file("'Gbdraw Test Fam'") != serif_name


def test_unregistering_returns_the_family_to_substitution(tmp_path: Path) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "a.ttf")
    fonts.register_font_file(path)
    assert not fonts.resolve_font_family("'Gbdraw Test Serif', serif").substituted

    removed = fonts.unregister_font_file(path)

    assert [font.family for font in removed] == ["Gbdraw Test Serif"]
    assert calculate_bbox_dimensions(TEXT, "'Gbdraw Test Serif'", 14, 72) == get_text_bbox_size_pixels(
        _bundled("LiberationSans-Regular"), TEXT, 14, 72
    )
    assert fonts.resolve_font_family("'Gbdraw Test Serif', serif").substituted


def test_a_registered_family_uses_its_regular_face_for_a_missing_variant(tmp_path: Path) -> None:
    regular = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "r.ttf")
    bold = _renamed_font("LiberationSerif-Bold", "Gbdraw Test Serif", tmp_path / "b.ttf", bold=True)
    fonts.register_font_file(regular)
    fonts.register_font_file(bold)

    assert fonts.resolve_font_face("'Gbdraw Test Serif'", "bold").path == str(bold)
    assert fonts.resolve_font_face("'Gbdraw Test Serif'", "normal", "italic").path == str(regular)


def test_measurement_names_bundled_faces_by_stem_and_registered_faces_by_their_bytes(tmp_path: Path) -> None:
    assert measurement_font_file(DEFAULT_FAMILY) == "LiberationSans-Regular"
    assert measurement_font_file(DEFAULT_FAMILY, "bold", "italic") == "LiberationSans-BoldItalic"
    first = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "first.ttf")
    second = tmp_path / "second-copy.ttf"
    second.write_bytes(first.read_bytes())

    name = measurement_font_file(fonts.register_font_file(first).family)
    fonts.clear_registered_fonts()
    fonts.register_font_file(second)

    assert re.fullmatch(r"user-[0-9a-f]{16}-Regular", name or "")
    assert measurement_font_file("'Gbdraw Test Serif'") == name


def test_the_metrics_table_of_a_bundled_face_is_the_generated_one() -> None:
    face = fonts.resolve_font_face("serif")
    table = font_metrics_table(face)
    generated = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/utils/legend-font-metrics.generated.js").read_text()

    assert face.measurement_name == "LiberationSerif-Regular"
    assert f"unitsPerEm: {table['unitsPerEm']}," in generated
    assert len(table["glyphs"]) == 5 * len(TTFont(face.path).getGlyphOrder())


def test_a_file_that_is_not_a_font_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "not-a-font.ttf"
    path.write_bytes(b"not a font")

    with pytest.raises(Exception, match="Cannot read font file"):
        fonts.register_font_file(path)


# ---- resolution order ----


def test_a_named_family_before_a_later_generic_wins(tmp_path: Path, fake_fc_match) -> None:
    path = _renamed_font("LiberationMono-Regular", "Gbdraw System Mono", tmp_path / "m.ttf")
    fake_fc_match("Gbdraw System Mono", "Gbdraw System Mono", path)

    face = fonts.resolve_font_face("'Gbdraw System Mono', Arial, sans-serif")

    assert face is not None and face.path == str(path)


def test_fontconfig_fallback_of_another_family_is_rejected(tmp_path: Path, fake_fc_match) -> None:
    path = _renamed_font("LiberationMono-Regular", "DejaVu Sans", tmp_path / "d.ttf")
    fake_fc_match("*", "DejaVu Sans", path)

    face = fonts.resolve_font_face("'Comic Sans MS', serif")

    assert face is not None and face.measurement_name == "LiberationSerif-Regular"


def test_bundled_aliases_never_ask_fontconfig(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fonts.shutil, "which", lambda *_args: pytest.fail("fc-match was asked"))

    assert fonts.resolve_font_face(DEFAULT_FAMILY, "bold").measurement_name == "LiberationSans-Bold"
    assert fonts.resolve_font_family(DEFAULT_FAMILY).effective == DEFAULT_FAMILY


def test_without_fontconfig_a_system_family_is_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fonts.shutil, "which", lambda *_args: None)

    assert fonts.resolve_font_family("'DejaVu Sans', sans-serif").missing == ("DejaVu Sans",)


# ---- substitution ----


@pytest.mark.parametrize(
    ("requested", "effective"),
    [
        (DEFAULT_FAMILY, DEFAULT_FAMILY),
        ("Arial", "Arial"),
        ("'Comic Sans MS', serif", "'Liberation Serif', serif"),
        ("'Comic Sans MS', monospace", "'Liberation Mono', monospace"),
        ("'Comic Sans MS', cursive", "'Liberation Sans', cursive"),
        ("'Comic Sans MS'", "'Liberation Sans', sans-serif"),
        ("'Fira Code'", "'Liberation Mono', monospace"),
        ("'Noto Serif CJK JP'", "'Liberation Serif', serif"),
        ("'Hiragino Mincho ProN'", "'Liberation Serif', serif"),
        ("'Noto Sans'", "'Liberation Sans', sans-serif"),
        ("'Comic Sans MS', 'Times New Roman', serif", "'Times New Roman', serif"),
    ],
)
def test_missing_families_are_replaced_with_the_nearest_liberation_family(
    requested: str, effective: str, fake_fc_match
) -> None:
    resolution = fonts.resolve_font_family(requested)

    assert resolution.effective == effective
    assert resolution.substituted == (requested != effective)


def test_the_substitution_message_names_the_missing_font_and_its_replacement(fake_fc_match) -> None:
    resolution = fonts.resolve_font_family("'Comic Sans MS', serif")

    assert resolution.missing == ("Comic Sans MS",)
    assert resolution.message == (
        "Font 'Comic Sans MS' is not available; text uses 'Liberation Serif' instead "
        "(font-family: 'Liberation Serif', serif)."
    )
    assert fonts.resolve_font_family(DEFAULT_FAMILY).message is None


def test_the_default_config_font_family_comes_out_unchanged() -> None:
    requested = load_default_config()["objects"]["text"]["font_family"]

    assert requested == DEFAULT_FAMILY
    assert fonts.resolve_font_family(requested).effective == requested


# ---- generic family of a file ----


@pytest.mark.parametrize(
    ("stem", "generic"),
    [
        ("LiberationSans-Regular", "sans-serif"),
        ("LiberationSerif-Regular", "serif"),
        ("LiberationMono-Regular", "monospace"),
    ],
)
def test_the_generic_family_of_a_font_file(stem: str, generic: str) -> None:
    assert fonts.font_generic_family(_bundled(stem)) == generic


def test_one_family_name_gets_its_generic(tmp_path: Path, fake_fc_match) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw System Serif", tmp_path / "s.ttf")
    fake_fc_match("Gbdraw System Serif", "Gbdraw System Serif", path)

    assert fonts.complete_font_family("Gbdraw System Serif") == "'Gbdraw System Serif', serif"
    assert fonts.complete_font_family("Courier New") == "'Courier New', monospace"
    assert fonts.complete_font_family("serif") == "serif"
    assert fonts.complete_font_family("Missing Family") == "Missing Family"
    assert fonts.complete_font_family("A, B") == "A, B"
    assert fonts.font_family_list_for("O'Brien Sans", "sans-serif") == "\"O'Brien Sans\", sans-serif"


# ---- CLI and Session ----


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


def _font_families(svg: Path) -> set[str]:
    return set(re.findall(r'font-family="([^"]*)"', svg.read_text()))


def test_cli_font_family_draws_and_measures_with_the_font(tmp_path: Path, fake_fc_match) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw System Serif", tmp_path / "s.ttf")
    fake_fc_match("Gbdraw System Serif", "Gbdraw System Serif", path)
    gbk = _genbank(tmp_path / "in.gb")
    session_path = tmp_path / "font.gbdraw-session.json"

    circular_main([
        "--gbk", str(gbk), "-o", str(tmp_path / "font"), "-f", "svg",
        "--font_family", "Gbdraw System Serif",
        "--save_session", "--session_output", str(session_path),
    ])

    assert _font_families(tmp_path / "font.svg") == {"'Gbdraw System Serif', serif"}
    assert "'Gbdraw System Serif', serif" in json.dumps(json.loads(session_path.read_text())["renderRequest"])
    circular_main(["--session", str(session_path), "-o", str(tmp_path / "replay"), "-f", "svg"])
    assert (tmp_path / "replay.svg").read_text() == (tmp_path / "font.svg").read_text()


@pytest.mark.parametrize("main", [circular_main, linear_main], ids=["circular", "linear"])
def test_cli_with_a_missing_font_warns_and_substitutes(
    tmp_path: Path, fake_fc_match, caplog: pytest.LogCaptureFixture, main
) -> None:
    gbk = _genbank(tmp_path / "in.gb")

    with caplog.at_level(logging.WARNING):
        main([
            "--gbk", str(gbk), "-o", str(tmp_path / "missing"), "-f", "svg",
            "--font_family", "'Gbdraw Missing Serif', serif",
        ])

    assert "Font 'Gbdraw Missing Serif' is not available" in caplog.text
    families = _font_families(tmp_path / "missing.svg")
    assert families and all("Gbdraw Missing Serif" not in family for family in families)
    assert any("Liberation Serif" in family for family in families)


def test_a_gui_session_font_family_becomes_the_cli_flag(tmp_path: Path) -> None:
    session = {
        "format": "gbdraw-session",
        "version": 30,
        "createdAt": "2026-10-13T00:00:00Z",
        "config": {"form": {"prefix": "out"}, "adv": {"font_family": "'Noto Sans', sans-serif"}},
        "ui": {"mode": "circular", "cInputType": "gb", "lInputType": "gb"},
        "files": {
            "c_gb": {
                "name": "input.gb",
                "type": "application/octet-stream",
                "size": 17,
                "lastModified": 0,
                "data": base64.b64encode(b"LOCUS       TEST\n").decode("ascii"),
            }
        },
    }

    spec = session_to_cli_args(session, mode="circular", temp_dir=tmp_path, output_override=None, format_override=None)

    assert spec.args[spec.args.index("--font_family") + 1] == "'Noto Sans', sans-serif"
