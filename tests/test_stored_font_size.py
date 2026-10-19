"""OV-1282 (Owner Q-13 a): a Session that stores a font size <= 0 loads in the
CLI and the Python API as Auto, with the Web's lost-value notice; a size <= 0
entered as an option is still refused (PD-OI-081)."""

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
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.root import _typed_leaves, style_leaf_domain
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import CanonicalRequestDecodingError, decode_canonical_request

REPO_ROOT = Path(__file__).resolve().parents[1]
SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
# Web Save (configOverrides), CLI (full config), and a released 0.13.0 Web
# Session 30 without a request (its draft replays through CLI arguments).
WEB_OVERRIDES = SESSIONS / "lambda_basic_linear.v44-schema8.gbdraw-session.json.gz"
CLI_CONFIG = SESSIONS / "cli-linear-tables.v44.gbdraw-session.json.gz"
LEGACY_WEB_DRAFT = SESSIONS / "hash-color-rule-linear-rc.v30.gbdraw-session.json.gz"
CIRCULAR_WEB_OVERRIDES = SESSIONS / "HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz"
PACKAGED = load_config_toml("gbdraw.data", "config.toml")
# Every font size the typed configuration constrains (a POSITIVE style leaf).
FONT_SIZE_LEAVES = tuple(
    path
    for path, _ in _typed_leaves(GbdrawConfig.from_dict(PACKAGED))
    if (style_leaf_domain(path) or ("",))[0] == "POSITIVE"
)
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
    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
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


def _packaged(path: str) -> object:
    tree: object = PACKAGED
    for key in path.split("."):
        tree = tree.get(key) if isinstance(tree, dict) else None
    return tree


def _decoded(fixture: Path, change, tmp_path: Path) -> object:
    session = _session(fixture)
    change(session["renderRequest"]["diagramOptions"])
    with materialize_session(load_session_document(_session(fixture)), output_directory=tmp_path) as materialized:
        try:
            options = decode_canonical_request(
                session["renderRequest"], resource_paths=materialized.resource_paths, output_directory=tmp_path
            ).options
            return options.config, options.config_overrides
        except CanonicalRequestDecodingError as error:
            return str(error)


def _set_override(path: str, value: object):
    def change(options: dict[str, Any]) -> None:
        options.setdefault("configOverrides", {})[path] = value
    return change


def _drop_override(path: str):
    def change(options: dict[str, Any]) -> None:
        options.setdefault("configOverrides", {}).pop(path, None)
    return change


def _set_config(path: str, value: object):
    def change(options: dict[str, Any]) -> None:
        *parents, leaf = path.split(".")
        tree = options["config"]
        for key in parents:
            tree = tree.setdefault(key, {})
        tree[leaf] = value
    return change


def test_every_typed_font_size_leaf_is_covered() -> None:
    assert len(FONT_SIZE_LEAVES) == 22
    assert {"objects.legends.font_size.long", "labels.font_size.long", "objects.scale.font_size.long",
            "objects.definition.linear.font_size.long", "objects.gc_content.tick_font_size"} <= set(FONT_SIZE_LEAVES)


@pytest.mark.parametrize("path", FONT_SIZE_LEAVES)
@pytest.mark.parametrize("fixture", [WEB_OVERRIDES, CIRCULAR_WEB_OVERRIDES], ids=["linear", "circular"])
def test_a_stored_override_font_size_of_zero_or_less_is_dropped_with_the_web_notice(
    path: str, fixture: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
    # Also for a setting of the other mode, which a stored size > 0 makes refused.
    stored = _decoded(fixture, _set_override(path, 0), tmp_path)
    assert not isinstance(stored, str), stored
    assert stored == _decoded(fixture, _drop_override(path), tmp_path)
    assert f'Settings: the value "0" (renderRequest.diagramOptions.configOverrides["{path}"]) cannot be used' in caplog.text


@pytest.mark.parametrize("path", FONT_SIZE_LEAVES)
def test_a_stored_config_font_size_of_zero_or_less_reads_the_packaged_default(
    path: str, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
    stored = _decoded(CLI_CONFIG, _set_config(path, -2), tmp_path)
    assert not isinstance(stored, str), stored
    assert stored == _decoded(CLI_CONFIG, _set_config(path, _packaged(path)), tmp_path)
    assert f'Settings: the value "-2" (renderRequest.diagramOptions.config.{path}) cannot be used' in caplog.text


def test_a_legacy_web_draft_font_size_of_zero_or_less_replays_as_auto(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
    sizes = ("legend_font_size", "scale_font_size", "label_font_size", "def_font_size")
    stored, auto = _session(LEGACY_WEB_DRAFT), _session(LEGACY_WEB_DRAFT)
    for key in sizes:
        stored["config"]["adv"][key] = 0
        auto["config"]["adv"].pop(key, None)
    assert _cli_svg(stored, tmp_path, "stored") == _cli_svg(auto, tmp_path, "auto")
    for key in sizes:
        assert f'Settings: the value "0" (config.adv.{key}) cannot be used' in caplog.text
