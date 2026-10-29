"""The packaged palette ``herbal_garden`` was spelled ``hearbal_garden`` through
0.14.2. The former name still reads as the palette everywhere a palette name is
looked up; writers and every list show only ``herbal_garden``."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

import gbdraw.circular as circular_cli_module
from gbdraw import FeatureOptions
from gbdraw.api import ColorOptions
from gbdraw.io.colors import PALETTE_NAME_ALIASES, packaged_color_palettes, require_palette_name

REPO_ROOT = Path(__file__).resolve().parents[1]
HUMAN = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
OLD, NEW = "hearbal_garden", "herbal_garden"
RELEASED_SESSION = (
    REPO_ROOT / "tests" / "fixtures" / "sessions" / "HmmtDNA-hearbal-garden-cli.v48.gbdraw-session.json.gz"
)


def test_the_palette_is_listed_by_its_current_name_only() -> None:
    palettes = packaged_color_palettes()
    assert NEW in list(palettes)
    assert OLD not in list(palettes)
    assert PALETTE_NAME_ALIASES == {OLD: NEW}
    # A lookup by the former name finds the same palette.
    assert palettes[OLD] == palettes[NEW]
    assert palettes.get(OLD) == palettes[NEW]
    assert require_palette_name(OLD) == OLD
    assert FeatureOptions(palette=OLD).palette == OLD
    # Every request reaches its palette through ColorOptions, which keeps the current name.
    assert ColorOptions(default_colors_palette=OLD).default_colors_palette == NEW
    assert ColorOptions(default_colors_palette="nope").default_colors_palette == "nope"
    with pytest.raises(Exception, match=r"choose from: .*herbal_garden") as excinfo:
        require_palette_name("nope")
    assert OLD not in str(excinfo.value)


@pytest.mark.parametrize(
    "path",
    [
        "gbdraw/data/color_palettes.toml",
        "gbdraw/web/gallery/palettes/palettes.json",
        "examples/color_palette_examples.md",
    ],
)
def test_lists_of_palettes_show_the_current_name(path: str) -> None:
    source = (REPO_ROOT / path).read_text(encoding="utf-8")
    assert NEW in source
    assert OLD not in source


def _render(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str, palette: str) -> Path:
    directory = tmp_path / name
    directory.mkdir()
    monkeypatch.chdir(directory)
    circular_cli_module.circular_main(
        ["--gbk", str(HUMAN), "-p", palette, "-o", "out", "-f", "svg", "--save_session"]
    )
    return directory


def test_the_former_cli_name_renders_like_the_current_one_and_saves_the_current_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    new = _render(tmp_path, monkeypatch, "new", NEW)
    old = _render(tmp_path, monkeypatch, "old", OLD)
    assert (old / "out.svg").read_bytes() == (new / "out.svg").read_bytes()
    session = json.loads((old / "out.gbdraw-session.json").read_text(encoding="utf-8"))
    assert session["renderRequest"]["diagramOptions"]["colors"]["defaultColorsPalette"] == NEW
    # The invocation record keeps the command line as typed.
    assert session["cliInvocation"]["args"][session["cliInvocation"]["args"].index("-p") + 1] == OLD


def test_a_released_session_that_names_the_former_palette_replays_like_the_current_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """gbdraw 0.14.2 wrote this Session for ``-p hearbal_garden`` (its provenance.json)."""

    text = gzip.decompress(RELEASED_SESSION.read_bytes()).decode("utf-8")
    stored = json.loads(text)
    assert stored["renderRequest"]["diagramOptions"]["colors"]["defaultColorsPalette"] == OLD
    svgs = {}
    for palette in (OLD, NEW, "default"):
        directory = tmp_path / palette
        directory.mkdir()
        session = directory / "in.gbdraw-session.json"
        session.write_text(text.replace(f'"{OLD}"', f'"{palette}"'), encoding="utf-8")
        monkeypatch.chdir(directory)
        circular_cli_module.circular_main(
            ["--session", str(session), "-o", "replay", "-f", "svg", "--session_output", "resaved.json"]
        )
        svgs[palette] = (directory / "replay.svg").read_bytes()
    assert svgs[OLD] == svgs[NEW]
    assert svgs[NEW] != svgs["default"]
    resaved = json.loads((tmp_path / OLD / "resaved.json").read_text(encoding="utf-8"))
    assert resaved["renderRequest"]["diagramOptions"]["colors"]["defaultColorsPalette"] == NEW
