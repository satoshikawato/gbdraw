"""OV-1282 (Owner Q-13 a): a Session that stores a Linear definition line font
size <= 0 loads in the CLI and the Python API as Auto, with the Web's
lost-value notice; a size <= 0 entered as an option is still refused
(PD-OI-081)."""

from __future__ import annotations

import gzip
import json
import logging
import re
from pathlib import Path
from typing import Any

import pytest

import gbdraw.linear as linear_cli_module
from gbdraw.api import LinearDiagramOptions, load_session_document, materialize_session, render_session
from gbdraw.exceptions import ValidationError

REPO_ROOT = Path(__file__).resolve().parents[1]
SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
# Web Save (configOverrides), CLI (full config), and a released 0.13.0 Web
# Session 30 without a request (its draft replays through CLI arguments).
WEB_OVERRIDES = SESSIONS / "lambda_basic_linear.v44-schema8.gbdraw-session.json.gz"
CLI_CONFIG = SESSIONS / "cli-linear-tables.v44.gbdraw-session.json.gz"
LEGACY_WEB_DRAFT = SESSIONS / "hash-color-rule-linear-rc.v30.gbdraw-session.json.gz"
SUBTITLE = "objects.definition.linear.line_styles.subtitle.font_size"
NAME = "objects.definition.linear.line_styles.name.font_size"


def _session(path: Path) -> dict[str, Any]:
    return json.loads(gzip.decompress(path.read_bytes()))


def _with_sizes(fixture: Path, sizes: dict[str, object]) -> dict[str, Any]:
    session = _session(fixture)
    if fixture == LEGACY_WEB_DRAFT:
        styles = session["config"]["adv"]["linear_definition_line_styles"]
        for path, size in sizes.items():
            styles[path.split(".")[-2]]["font_size"] = size
        return session
    options = session["renderRequest"]["diagramOptions"]
    for path, size in sizes.items():
        if fixture == CLI_CONFIG:
            options["config"]["objects"]["definition"]["linear"]["line_styles"][path.split(".")[-2]]["font_size"] = size
        else:
            options["configOverrides"][path] = size
    return session


def _cli_svg(session: dict[str, Any], tmp_path: Path, name: str) -> bytes:
    path = tmp_path / f"{name}.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    linear_cli_module.linear_main(["--session", str(path), "-o", name, "-f", "svg"])
    return (tmp_path / f"{name}.svg").read_bytes()


def _python_svgs(session: dict[str, Any], tmp_path: Path, name: str) -> list[bytes]:
    output_directory = tmp_path / name
    with materialize_session(load_session_document(session), output_directory=output_directory) as materialized:
        render_session(materialized)
    return sorted(path.read_bytes() for path in output_directory.glob("*.svg"))


@pytest.mark.parametrize(
    ("fixture", "where"),
    [
        (WEB_OVERRIDES, 'renderRequest.diagramOptions.configOverrides["{path}"]'),
        (CLI_CONFIG, "renderRequest.diagramOptions.config.{path}"),
        (LEGACY_WEB_DRAFT, "config.adv.linear_definition_line_styles.{kind}.font_size"),
    ],
    ids=["web-overrides", "cli-config", "legacy-web-draft"],
)
def test_a_stored_size_of_zero_or_less_draws_auto_with_the_web_notice(
    fixture: Path, where: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.chdir(tmp_path)
    stored = _with_sizes(fixture, {SUBTITLE: 0, NAME: -1.5})
    auto = _with_sizes(fixture, {SUBTITLE: None, NAME: None})
    caplog.set_level(logging.WARNING, logger="gbdraw.definition_line_styles")
    assert _cli_svg(stored, tmp_path, "stored") == _cli_svg(auto, tmp_path, "auto")
    notices = [record.getMessage() for record in caplog.records if "cannot be used" in record.getMessage()]
    expected = [
        f'WARNING: Settings: the value "{value}" ({where.format(path=path, kind=path.split(".")[-2])}) '
        "cannot be used; its default is used."
        for path, value in ((NAME, "-1.5"), (SUBTITLE, "0"))
    ]
    assert sorted(notices) == sorted(expected)
    if fixture != LEGACY_WEB_DRAFT:  # Sessions 27-30 have no typed request.
        assert _python_svgs(stored, tmp_path, "py-stored") == _python_svgs(auto, tmp_path, "py-auto")


def test_an_entered_size_of_zero_or_less_is_still_refused(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(ValidationError, match="font sizes must be finite numbers greater than zero"):
        LinearDiagramOptions(config_overrides={SUBTITLE: 0})
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(
            ["--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"), "--definition_line_style", "subtitle:size=0"]
        )
    assert "positive finite number" in capsys.readouterr().err


def test_the_notice_keeps_the_web_load_wording() -> None:
    wording = (REPO_ROOT / "gbdraw" / "web" / "js" / "utils" / "load-repairs.js").read_text(encoding="utf-8")
    assert re.search(r"`Settings: \$\{stored === undefined \? `the value of \$\{path\}` : `the value \"\$\{stored\}\" \(\$\{path\}\)`\} cannot be used; `", wording)
    assert "defaulted: 'its default is used'" in wording
