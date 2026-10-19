"""OV-1244 (Owner Q-07 a): an unknown output format or palette name is refused
where a user enters it, with an error that names the valid choices; a released
legacy Session that stores one still replays with the old fallback."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw import FeatureOptions
from gbdraw.api import RenderOutputRequest, save_figure_to
from gbdraw.exceptions import ValidationError
from gbdraw.render.formats import ACCEPTED_FORMATS

REPO_ROOT = Path(__file__).resolve().parents[1]
HUMAN = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
# The 0.13.0 CLI wrote this Session (scale-interval.provenance.json); its args replay.
RELEASED_CLI_LINEAR = (
    REPO_ROOT / "tests" / "fixtures" / "sessions" / "scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz"
)
CLI_MODULES = {"circular": circular_cli_module, "linear": linear_cli_module}
FORMAT_CHOICES = ", ".join(ACCEPTED_FORMATS)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    ("option", "value", "named"),
    [
        ("-f", "foo", f"choose from: {FORMAT_CHOICES}"),
        ("-f", "svg,foo", f"choose from: {FORMAT_CHOICES}"),
        ("-p", "nope", "choose from: default, mint,"),
    ],
)
def test_cli_refuses_an_unknown_format_or_palette(mode: str, option: str, value: str, named: str) -> None:
    with pytest.raises(ValidationError, match=repr(value.split(",")[-1])) as excinfo:
        CLI_MODULES[mode]._get_args(["--gbk", str(HUMAN), option, value])
    assert named in str(excinfo.value)
    assert excinfo.value.diagnostic == {"code": "INPUT_INVALID"}


def test_cli_run_with_an_unknown_format_writes_nothing(tmp_path: Path, gbdraw_runner) -> None:
    returncode, output, _ = gbdraw_runner.run(
        "circular", [HUMAN], "unknown_format", tmp_path, extra_args=["-f", "foo"]
    )
    assert returncode != 0
    assert f"ERROR: Unknown output format 'foo'; choose from: {FORMAT_CHOICES}." in output
    assert not list(tmp_path.glob("unknown_format*"))


def test_python_options_refuse_an_unknown_palette_or_format(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match=r"Unknown palette 'nope'; choose from: default, ") as excinfo:
        FeatureOptions(palette="nope")
    assert excinfo.value.diagnostic == {"code": "INPUT_INVALID"}
    assert FeatureOptions(palette="mint").palette == "mint"
    with pytest.raises(ValidationError, match=f"choose from: {FORMAT_CHOICES}"):
        RenderOutputRequest(output_directory=tmp_path, formats=("svg", "foo"))
    with pytest.raises(ValidationError, match=f"Unknown output format 'foo'; choose from: {FORMAT_CHOICES}"):
        save_figure_to(object(), "foo", output_prefix=str(tmp_path / "out"))  # type: ignore[arg-type]
    assert not list(tmp_path.iterdir())


def _session(path: Path) -> dict[str, Any]:
    return json.loads(gzip.decompress(path.read_bytes()))


def test_a_released_legacy_session_with_an_unknown_palette_or_format_still_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    svgs = {}
    for name, extra in (("stored", ["-p", "nope", "-f", "foo,svg"]), ("default", ["-f", "svg"])):
        session = _session(RELEASED_CLI_LINEAR)
        assert session["version"] == 30
        args = session["cliInvocation"]["args"]
        del args[args.index("-f") : args.index("-f") + 2]
        args[args.index("-o") + 1] = name
        args.extend(extra)
        path = tmp_path / f"{name}.gbdraw-session.json"
        path.write_text(json.dumps(session), encoding="utf-8")
        linear_cli_module.linear_main(["--session", str(path)])
        svgs[name] = (tmp_path / f"{name}.svg").read_text(encoding="utf-8").replace(name, "OUT")
    assert svgs["stored"] == svgs["default"]


def test_a_stored_request_with_an_unknown_palette_still_renders(tmp_path: Path) -> None:
    from gbdraw.api import load_session_document, materialize_session, render_session

    fixture = REPO_ROOT / "tests" / "fixtures" / "sessions" / "HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz"
    svgs = []
    for palette in ("nope", "default"):
        session = _session(fixture)
        colors = session["renderRequest"]["diagramOptions"]["colors"]
        assert colors["defaultColorsPalette"] == "default"
        colors["defaultColorsPalette"] = palette
        output_directory = tmp_path / palette
        with materialize_session(load_session_document(session), output_directory=output_directory) as materialized:
            render_session(materialized)
        svgs.append(sorted(path.read_bytes() for path in output_directory.glob("*.svg")))
    assert svgs[0] and svgs[0] == svgs[1]


@pytest.mark.parametrize("command", ["render", "linear"])
def test_a_session_format_override_refuses_an_unknown_format(command: str) -> None:
    from gbdraw.cli_utils.session import parse_render_args, parse_session_pre_args

    args = ["--session", str(RELEASED_CLI_LINEAR), "-f", "svg,foo"]
    with pytest.raises(ValidationError, match=f"Unknown output format 'foo'; choose from: {FORMAT_CHOICES}"):
        parse_render_args(args) if command == "render" else parse_session_pre_args(args, mode="linear")
    request = parse_render_args(["--session", str(RELEASED_CLI_LINEAR), "-f", "SVG,interactive-svg"])
    assert request.format == "svg,interactive_svg"
