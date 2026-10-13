"""The Web's font registration helper (custom fonts PR-2, review R3 and R7)."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from gbdraw.core.fonts import clear_registered_fonts, font_family_list_for, registered_fonts

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def helpers():
    source = (REPO_ROOT / "gbdraw/web/js/app/python-helpers.js").read_text()
    namespace: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


@pytest.fixture(autouse=True)
def _clean_registry():
    clear_registered_fonts()
    yield
    clear_registered_fonts()


def _register(helpers, entries):
    return json.loads(helpers["register_web_fonts_json"](json.dumps(entries)))


def test_the_reply_names_the_stored_list_python_writes(helpers, tmp_path):
    path = tmp_path / "serif.ttf"
    shutil.copy(REPO_ROOT / "gbdraw/data/LiberationSerif-Regular.ttf", path)
    reply = _register(helpers, [{"identity": "a" * 64, "path": str(path)}])
    [font] = reply["fonts"]
    assert font["fontFamily"] == font_family_list_for(font["family"], font["generic"])


def test_a_call_registers_all_of_its_files_or_none(helpers, tmp_path):
    good = tmp_path / "serif.ttf"
    shutil.copy(REPO_ROOT / "gbdraw/data/LiberationSerif-Regular.ttf", good)
    broken = tmp_path / "broken.ttf"
    broken.write_bytes(b"not a font")
    reply = _register(helpers, [
        {"identity": "a" * 64, "path": str(good)},
        {"identity": "b" * 64, "path": str(broken)},
    ])
    assert reply["error"]["code"] == "INPUT_INVALID"
    assert registered_fonts() == ()


def test_s5_a_failed_call_leaves_earlier_registrations_and_replaced_faces_as_they_were(helpers, tmp_path):
    from fontTools.ttLib import TTFont

    first = tmp_path / "v1.ttf"
    shutil.copy(REPO_ROOT / "gbdraw/data/LiberationSerif-Regular.ttf", first)
    second = tmp_path / "v2.ttf"
    font = TTFont(str(first))
    font["head"].fontRevision = 9.0  # other bytes, same family, weight and style
    font.save(str(second))
    broken = tmp_path / "broken.ttf"
    broken.write_bytes(b"not a font")

    assert "fonts" in _register(helpers, [{"identity": "a" * 64, "path": str(first)}])
    before = registered_fonts()
    reply = _register(helpers, [
        {"identity": "b" * 64, "path": str(second)},
        {"identity": "c" * 64, "path": str(broken)},
    ])

    assert reply["error"]["code"] == "INPUT_INVALID"
    assert registered_fonts() == before
    assert before[0].face.path == str(first)
    from gbdraw.core.fonts import face_source_bytes
    assert face_source_bytes(before[0].face) == first.read_bytes()
