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
import gc
import json
import logging
import os
import re
import subprocess
import weakref
from pathlib import Path

import pytest
from fontTools.ttLib import TTCollection, TTFont

from gbdraw.api.config import load_default_config
import gbdraw.api.diagram as diagram_module
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main
from gbdraw.core import fonts
from gbdraw.core.text import (
    calculate_bbox_dimensions,
    calculate_svg_bbox_dimensions,
    font_metrics_table,
    font_pair_kerning_table,
    get_text_bbox_size_pixels,
    measurement_font_file,
)
from gbdraw.render.font_embedding import remove_embedded_font_css
from gbdraw.session_io import session_to_cli_args

DATA_DIR = Path(__file__).resolve().parents[1] / "gbdraw" / "data"
DEFAULT_FAMILY = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
TEXT = "Hypothetical protein"


def _bundled(stem: str) -> str:
    return str(DATA_DIR / f"{stem}.ttf")


def _renamed(source: str, family: str, *, weight: int = 400) -> TTFont:
    font = TTFont(_bundled(source))
    name_table = font["name"]
    name_table.names = [record for record in name_table.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
    style = f"W{weight}"
    for name_id, value in ((1, family), (2, style), (4, f"{family} {style}"), (6, family.replace(" ", "")), (16, family)):
        name_table.setName(value, name_id, 3, 1, 0x409)
    font["OS/2"].usWeightClass = weight
    return font


def _renamed_font(source: str, family: str, target: Path, *, weight: int = 400) -> Path:
    _renamed(source, family, weight=weight).save(str(target))
    return target


@pytest.fixture(autouse=True)
def isolated_fonts(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(fonts, "_fontconfig_cache", {})
    # No platform font folders unless a test names them.
    monkeypatch.setattr(fonts, "_platform_index", {}, raising=False)
    monkeypatch.setattr(diagram_module, "_warned_font_families", set(), raising=False)
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


def test_custom_font_families_among_counts_a_registered_alias_name() -> None:
    # D-10: a loaded alias name is custom; an alias that was not loaded is not.
    assert fonts.custom_font_families_among("'Times New Roman', 'DejaVu Serif', serif", ()) == ("DejaVu Serif",)
    assert fonts.custom_font_families_among(
        "'Times New Roman', 'DejaVu Serif', serif", ("times new  roman",)
    ) == ("Times New Roman", "DejaVu Serif")
    assert fonts.custom_font_families_among("serif", ("serif",)) == ()


def test_a_registered_font_wins_over_the_bundled_alias_of_its_name(tmp_path: Path) -> None:
    # Review R2: a registered "Times New Roman" is measured with its own bytes.
    assert fonts.resolve_font_face("'Times New Roman', serif").measurement_name == "LiberationSerif-Regular"
    assert fonts.custom_font_families("'Times New Roman', serif") == ()
    path = _renamed_font("LiberationMono-Regular", "Times New Roman", tmp_path / "tnr.ttf")

    registered = fonts.register_font_file(path)

    assert fonts.resolve_font_face("'Times New Roman', serif") == registered.face
    assert calculate_bbox_dimensions(TEXT, "'Times New Roman', serif", 14, 72) == get_text_bbox_size_pixels(
        path, TEXT, 14, 72
    )
    # Generic keywords keep their Liberation face, and other aliases keep today's mapping.
    assert fonts.resolve_font_face("serif").measurement_name == "LiberationSerif-Regular"
    assert fonts.resolve_font_face("Arial").measurement_name == "LiberationSans-Regular"
    # A Session elsewhere draws it with Liberation Serif, so it is a custom font (D-02).
    assert fonts.custom_font_families("'Times New Roman', serif") == ("Times New Roman",)
    fonts.unregister_font_file(path)
    assert fonts.resolve_font_face("'Times New Roman', serif").measurement_name == "LiberationSerif-Regular"


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


@pytest.mark.parametrize("order", ["light first", "black first"])
def test_registered_faces_are_chosen_by_css_weight_matching(tmp_path: Path, order: str) -> None:
    """Local Font Access returns every face of a family, in any order (review F1)."""

    faces = {
        300: _renamed_font("LiberationSans-Regular", "Gbdraw Weights", tmp_path / "w300.ttf", weight=300),
        400: _renamed_font("LiberationSerif-Regular", "Gbdraw Weights", tmp_path / "w400.ttf", weight=400),
        700: _renamed_font("LiberationSerif-Bold", "Gbdraw Weights", tmp_path / "w700.ttf", weight=700),
        900: _renamed_font("LiberationMono-Bold", "Gbdraw Weights", tmp_path / "w900.ttf", weight=900),
    }
    weights = [300, 400, 700, 900] if order == "light first" else [900, 700, 400, 300]
    for weight in weights:
        fonts.register_font_file(faces[weight])

    picked = {
        requested: Path(fonts.resolve_font_face("'Gbdraw Weights'", requested).path).name
        for requested in ("normal", "bold", 200, 300, 500, 600, 800, 900)
    }

    assert picked == {
        "normal": "w400.ttf", "bold": "w700.ttf", 200: "w300.ttf", 300: "w300.ttf",
        500: "w400.ttf", 600: "w700.ttf", 800: "w900.ttf", 900: "w900.ttf",
    }
    assert calculate_bbox_dimensions(TEXT, "'Gbdraw Weights'", 14, 72) == get_text_bbox_size_pixels(
        faces[400], TEXT, 14, 72
    )


def test_a_registered_family_uses_its_regular_face_for_a_missing_variant(tmp_path: Path) -> None:
    regular = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "r.ttf")
    bold = _renamed_font("LiberationSerif-Bold", "Gbdraw Test Serif", tmp_path / "b.ttf", weight=700)
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


def test_an_unregistered_face_is_not_read_again_from_its_path(tmp_path: Path) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "a.ttf")
    face = fonts.register_font_file(path).face
    fonts.unregister_font_file(path)
    _renamed_font("LiberationMono-Regular", "Gbdraw Test Serif", path)

    with pytest.raises(ValueError, match="no longer registered"):
        fonts.load_font(face)


@pytest.mark.parametrize("forget", ["unregister", "replace"])
def test_a_forgotten_face_keeps_no_kerning_or_font(tmp_path: Path, forget: str) -> None:
    """Unregistering a face, or registering another face of its family,
    weight and style, lets its font go: its kerning is not kept for it."""
    path = _renamed_font("LiberationSans-Regular", "Gbdraw Test Sans", tmp_path / "a.ttf")
    face = fonts.register_font_file(path).face
    font = fonts.load_font(face)
    assert font_pair_kerning_table(font)
    forgotten = weakref.ref(font)
    del font

    if forget == "unregister":
        fonts.unregister_font_file(path)
    else:
        fonts.register_font_file(_renamed_font("LiberationSans-Bold", "Gbdraw Test Sans", tmp_path / "b.ttf"))
    gc.collect()

    assert forgotten() is None


def test_a_font_at_a_reused_address_gets_its_own_kerning() -> None:
    """The kerning cache is keyed by the font object: a font allocated where a
    freed font was never gets that font's kerning."""
    data_dir = Path(fonts.__file__).resolve().parents[1] / "data"
    sans = TTFont(str(data_dir / "LiberationSans-Regular.ttf"))
    assert font_pair_kerning_table(sans)
    freed_id = id(sans)
    del sans
    gc.collect()

    candidates = []
    for _attempt in range(64):
        mono = TTFont(str(data_dir / "LiberationMono-Regular.ttf"))
        if id(mono) == freed_id:
            break
        candidates.append(mono)

    assert font_pair_kerning_table(mono) == {}


def test_a_registered_font_has_the_metrics_table_of_its_bytes(tmp_path: Path) -> None:
    """The Web measures a registered font's Legend with this table; the bundled
    tables are checked byte for byte by tests/test_legend_layout_vectors.py."""

    bundled = fonts.resolve_font_face("serif")
    registered = fonts.register_font_file(
        _renamed_font("LiberationSerif-Regular", "Gbdraw Test Serif", tmp_path / "a.ttf")
    ).face
    table = font_metrics_table(registered)

    assert registered.measurement_name != bundled.measurement_name
    assert {key: value for key, value in table.items() if key != "cmap"} == {
        key: value for key, value in font_metrics_table(bundled).items() if key != "cmap"
    }
    assert table["cmap"] == font_metrics_table(bundled)["cmap"]
    assert len(table["glyphs"]) == 5 * len(TTFont(bundled.path).getGlyphOrder())


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


def _without_font_names(svg: Path) -> str:
    """The SVG without the font: its embedded subset, font-family and the recorded face."""
    text = re.sub(r'font-family="[^"]*"', 'font-family=""', remove_embedded_font_css(svg.read_text()))
    return re.sub(r"fontFile&quot;:&quot;[^&]*&quot;", "fontFile&quot;:&quot;&quot;", text)


def test_cli_font_family_draws_and_measures_with_the_font(
    tmp_path: Path, fake_fc_match, caplog: pytest.LogCaptureFixture
) -> None:
    path = _renamed_font("LiberationMono-Regular", "Gbdraw System Mono", tmp_path / "m.ttf")
    fake_fc_match("Gbdraw System Mono", "Gbdraw System Mono", path)
    gbk = _genbank(tmp_path / "in.gb")
    session_path = tmp_path / "font.gbdraw-session.json"

    with caplog.at_level(logging.WARNING):
        circular_main([
            "--gbk", str(gbk), "-o", str(tmp_path / "font"), "-f", "svg",
            "--font_family", "Gbdraw System Mono",
            "--save_session", "--session_output", str(session_path),
        ])
    for prefix, family in (("same-metrics", "'Liberation Mono', monospace"), ("default", None)):
        circular_main(
            ["--gbk", str(gbk), "-o", str(tmp_path / prefix), "-f", "svg"]
            + (["--font_family", family] if family else [])
        )

    assert _font_families(tmp_path / "font.svg") == {"'Gbdraw System Mono', monospace"}
    # Laid out with the font's own metrics: the same as Liberation Mono, whose
    # outlines it has, and not as the default Liberation Sans.
    assert _without_font_names(tmp_path / "font.svg") == _without_font_names(tmp_path / "same-metrics.svg")
    assert _without_font_names(tmp_path / "font.svg") != _without_font_names(tmp_path / "default.svg")
    # D-02: the Session keeps the name only, and says so.
    assert "'Gbdraw System Mono', monospace" in json.dumps(json.loads(session_path.read_text())["renderRequest"])
    assert "The Session stores the font name 'Gbdraw System Mono' but not the font file" in caplog.text
    circular_main(["--session", str(session_path), "-o", str(tmp_path / "replay"), "-f", "svg"])
    assert (tmp_path / "replay.svg").read_text() == (tmp_path / "font.svg").read_text()


def test_a_session_with_bundled_fonts_says_nothing_about_fonts(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    gbk = _genbank(tmp_path / "in.gb")

    with caplog.at_level(logging.WARNING):
        circular_main(["--gbk", str(gbk), "-o", str(tmp_path / "plain"), "-f", "svg", "--save_session"])

    assert "font" not in caplog.text.lower()


def test_the_substitution_is_warned_once_per_family_and_run(
    tmp_path: Path, fake_fc_match, caplog: pytest.LogCaptureFixture
) -> None:
    gbk = _genbank(tmp_path / "in.gb")

    with caplog.at_level(logging.WARNING):
        for prefix in ("first", "second"):
            circular_main([
                "--gbk", str(gbk), "-o", str(tmp_path / prefix), "-f", "svg",
                "--font_family", "'Gbdraw Missing Once', serif",
            ])

    assert caplog.text.count("Font 'Gbdraw Missing Once' is not available") == 1


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


# ---- fontconfig failures, encodings, variable fonts (review F4-F6) ----


def test_a_fontconfig_timeout_is_reported_and_tried_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _renamed_font("LiberationSerif-Regular", "Gbdraw Slow Serif", tmp_path / "s.ttf")
    calls: list[float] = []

    def run(*_args, timeout, **_kwargs):
        calls.append(timeout)
        if len(calls) == 1:
            raise subprocess.TimeoutExpired("fc-match", timeout)
        return subprocess.CompletedProcess([], 0, stdout=f"Gbdraw Slow Serif\t{path}\t0", stderr="")

    monkeypatch.setattr(fonts.shutil, "which", lambda *_args: "/usr/bin/fc-match")
    monkeypatch.setattr(fonts.subprocess, "run", run)

    first = fonts.resolve_font_family("'Gbdraw Slow Serif', serif")
    second = fonts.resolve_font_family("'Gbdraw Slow Serif', serif")

    assert first.unchecked == ("Gbdraw Slow Serif",)
    assert first.message.startswith(
        "Font 'Gbdraw Slow Serif' could not be checked because fontconfig did not answer; text uses"
    )
    assert not second.substituted
    assert calls[0] >= 10


def test_fontconfig_output_is_read_as_utf8(tmp_path: Path, fake_fc_match) -> None:
    family = "Gbdraw 日本語 Gothic"
    path = _renamed_font("LiberationSans-Regular", family, tmp_path / "日本語.ttf")
    fake_fc_match(family, family, path)

    face = fonts.resolve_font_face(f"'{family}'")

    assert face is not None and face.path == str(path)


def _variable_font(path: Path) -> Path:
    """A one-glyph variable font: "A" is 500 units wide at wght 400 and 700 at 900."""

    from fontTools.fontBuilder import FontBuilder
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    from fontTools.ttLib.tables.TupleVariation import TupleVariation

    builder = FontBuilder(1000, isTTF=True)
    builder.setupGlyphOrder([".notdef", "A"])
    builder.setupCharacterMap({ord("A"): "A"})
    glyphs = {}
    for name in (".notdef", "A"):
        pen = TTGlyphPen(None)
        pen.moveTo((50, 0))
        pen.lineTo((50, 700))
        pen.lineTo((450, 700))
        pen.lineTo((450, 0))
        pen.closePath()
        glyphs[name] = pen.glyph()
    builder.setupGlyf(glyphs)
    builder.setupHorizontalMetrics({".notdef": (500, 50), "A": (500, 50)})
    builder.setupHorizontalHeader(ascent=800, descent=-200)
    builder.setupNameTable({"familyName": "Gbdraw Variable", "styleName": "Regular"})
    builder.setupOS2()
    builder.setupPost()
    builder.setupFvar(
        axes=[("wght", 100, 400, 900, "Weight")],
        instances=[{"location": {"wght": 900}, "stylename": "Black"}],
    )
    wider = [(0, 0), (0, 0), (200, 0), (200, 0), (0, 0), (200, 0), (0, 0), (0, 0)]
    builder.setupGvar({"A": [TupleVariation({"wght": (0, 1.0, 1.0)}, wider)], ".notdef": []})
    builder.save(str(path))
    return path


def test_a_fontconfig_named_instance_is_measured_as_that_instance(tmp_path: Path, fake_fc_match) -> None:
    path = _variable_font(tmp_path / "variable.ttf")
    fake_fc_match("Gbdraw Variable", "Gbdraw Variable", path)
    # fontconfig's index: named instance 1 in the high 16 bits, face 0.
    answers = json.loads((tmp_path / "fc-answers.json").read_text())
    answers["Gbdraw Variable"] = ["Gbdraw Variable", f"{path}\t{1 << 16}"]
    (tmp_path / "fc-answers.json").write_text(json.dumps(answers))

    face = fonts.resolve_font_face("'Gbdraw Variable'")

    assert (face.index, face.instance) == (0, 1)
    assert calculate_svg_bbox_dimensions("A", "'Gbdraw Variable'", 1000, 96)[0] == 700
    assert get_text_bbox_size_pixels(path, "A", 1000, 96, svg_units=True)[0] == 500


# ---- platform font folders without fontconfig (review F7) ----


def test_without_fontconfig_a_family_is_found_in_the_platform_font_folders(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    folder = tmp_path / "Fonts" / "nested"
    folder.mkdir(parents=True)
    regular = _renamed_font("LiberationSerif-Regular", "Gbdraw Folder Serif", folder / "r.ttf")
    bold = _renamed_font("LiberationSerif-Bold", "Gbdraw Folder Serif", folder / "b.ttf", weight=700)
    (folder / "broken.ttf").write_bytes(b"not a font")
    collection = TTCollection()
    collection.fonts = [
        _renamed("LiberationSans-Regular", "Gbdraw Pair"),
        _renamed("LiberationMono-Regular", "Gbdraw Pair Mono"),
    ]
    collection.save(str(folder / "pair.ttc"))
    indexed: list[tuple[str, ...]] = []
    real_index = fonts._index_font_dirs
    monkeypatch.setattr(fonts.shutil, "which", lambda *_args: None)
    monkeypatch.setattr(fonts, "_platform_index", None)
    monkeypatch.setattr(fonts, "_platform_font_dirs", lambda: (str(tmp_path / "Fonts"),))
    monkeypatch.setattr(fonts, "_index_font_dirs", lambda dirs: indexed.append(tuple(dirs)) or real_index(dirs))

    assert not fonts.resolve_font_family("'Gbdraw Folder Serif', serif").substituted
    assert fonts.resolve_font_face("'Gbdraw Folder Serif'", "bold").path == str(bold)
    assert fonts.resolve_font_face("'Gbdraw Folder Serif'").path == str(regular)
    assert fonts.resolve_font_face("'Gbdraw Pair Mono'") == fonts.FontFace(str(folder / "pair.ttc"), 1)
    assert fonts.complete_font_family("Gbdraw Pair Mono") == "'Gbdraw Pair Mono', monospace"
    assert fonts.resolve_font_family("'Not In The Folder', serif").substituted
    assert len(indexed) == 1


def test_pyodide_has_no_platform_font_folders(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fonts.sys, "platform", "emscripten")

    assert fonts._platform_font_dirs() == ()


def test_the_platform_font_folders_of_each_system(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(fonts.os.path, "isdir", lambda _path: True)
    monkeypatch.setenv("WINDIR", "C:\\Windows")
    monkeypatch.setenv("LOCALAPPDATA", "C:\\Users\\me\\AppData\\Local")
    monkeypatch.setattr(fonts.sys, "platform", "win32")
    windows = fonts._platform_font_dirs()
    monkeypatch.setattr(fonts.sys, "platform", "darwin")
    mac = fonts._platform_font_dirs()

    assert windows[0].endswith("Fonts") and "Microsoft" in windows[1]
    assert mac[:2] == ("/System/Library/Fonts", "/Library/Fonts")


def test_ov743_no_warning_when_the_first_family_is_drawn() -> None:
    resolution = fonts.resolve_font_family("'Liberation Serif', 'Gbdraw Never Installed', serif")

    assert resolution.replaced == ()
    assert resolution.message is None
    assert resolution.effective == "'Liberation Serif', serif"
    replaced = fonts.resolve_font_family("'Gbdraw Never Installed', 'Gbdraw Gone Too', serif")
    assert replaced.replaced == ("Gbdraw Never Installed", "Gbdraw Gone Too")
    assert replaced.message is not None

