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
from gbdraw.api.session_compat import decode_session_request
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
            options = decode_session_request(
                session["version"], session["renderRequest"],
                resource_paths=materialized.resource_paths, output_directory=tmp_path,
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


def test_a_fresh_request_still_refuses_a_font_size_of_zero_or_less(tmp_path: Path) -> None:
    """A Web Generate request is an input, not a stored Session (PD-OI-081)."""
    session = _session(WEB_OVERRIDES)
    session["renderRequest"]["diagramOptions"]["configOverrides"]["objects.legends.font_size.long"] = 0
    with materialize_session(load_session_document(_session(WEB_OVERRIDES)), output_directory=tmp_path) as materialized:
        with pytest.raises(CanonicalRequestDecodingError, match="font sizes must be finite numbers greater than zero"):
            decode_canonical_request(
                session["renderRequest"], resource_paths=materialized.resource_paths, output_directory=tmp_path
            )


# Written by the released 0.13.0 CLI (tag 0.13.0, 8be9e3a5): `gbdraw linear --gbk F1.gb
# --label_font_size 0` and `gbdraw circular --gbk F1.gb --label_font_size 0
# --legend_font_size 0`, each with --session_output; both drew. Provenance:
# tests/fixtures/sessions/legacy-font-size.provenance.json.
RELEASED_CLI = {
    "linear": SESSIONS / "legacy-font-size-zero-linear-cli.v30.gbdraw-session.json.gz",
    "circular": SESSIONS / "legacy-font-size-zero-circular-cli.v30.gbdraw-session.json.gz",
}


def _replay(mode: str, session: dict[str, Any], tmp_path: Path, name: str) -> bytes:
    import gbdraw.circular as circular_cli_module

    path = tmp_path / f"{name}.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    main = linear_cli_module.linear_main if mode == "linear" else circular_cli_module.circular_main
    main(["--session", str(path), "-o", name, "-f", "svg"])
    return (tmp_path / f"{name}.svg").read_bytes()


def _without_font_size_args(session: dict[str, Any]) -> dict[str, Any]:
    args = session["cliInvocation"]["args"]
    kept: list[str] = []
    skip = False
    for token in args:
        if skip:
            skip = False
        elif token.endswith("font_size"):
            skip = True
        else:
            kept.append(token)
    session["cliInvocation"]["args"] = kept
    return session


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_a_released_cli_session_with_a_font_size_of_zero_replays_as_auto(
    mode: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
    stored = _session(RELEASED_CLI[mode])
    assert stored["version"] == 30
    options = [token for token in stored["cliInvocation"]["args"] if token.endswith("font_size")]
    assert options and _replay(mode, stored, tmp_path, "stored") == _replay(
        mode, _without_font_size_args(_session(RELEASED_CLI[mode])), tmp_path, "auto"
    )
    for option in options:
        assert f'Settings: the value "0" (cliInvocation.args {option}) cannot be used' in caplog.text


def test_a_released_cli_session_with_a_font_size_still_invalid_gives_the_cli_error(tmp_path: Path) -> None:
    import os
    import subprocess
    import sys

    session = _session(RELEASED_CLI["linear"])
    args = session["cliInvocation"]["args"]
    args[args.index("--label_font_size") + 1] = "nan"
    path = tmp_path / "stored.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "linear", "--session", str(path), "-o", "out", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
    )
    assert result.returncode != 0
    assert "Traceback" not in result.stderr
    assert "label_font_size" in result.stderr
    assert not (tmp_path / "out.svg").exists()


def test_a_stored_unset_font_size_stays_unset_without_a_notice(caplog: pytest.LogCaptureFixture) -> None:
    from gbdraw.session_request_codec import with_stored_font_sizes_as_auto

    caplog.set_level(logging.WARNING, logger="gbdraw.config.models.root")
    unset = {"diagramOptions": {
        "config": {"objects": {"gc_content": {"tick_font_size": None}}, "labels": {"font_size": {"long": None}}},
        "configOverrides": {"objects.legends.font_size.long": None},
    }}
    assert with_stored_font_sizes_as_auto(unset) == unset
    assert "cannot be used" not in caplog.text
    zero = {"diagramOptions": {"config": {"labels": {"font_size": {"long": 0}}}}}
    assert with_stored_font_sizes_as_auto(zero)["diagramOptions"]["config"]["labels"]["font_size"]["long"] == 8
    assert "cannot be used" in caplog.text


def test_the_packaged_config_is_read_only_when_a_size_is_replaced(monkeypatch: pytest.MonkeyPatch) -> None:
    import gbdraw.session_request_codec as codec

    calls: list[tuple[str, str]] = []
    real = codec.load_config_toml

    def counting(directory: str, name: str) -> dict:
        calls.append((directory, name))
        return real(directory, name)

    monkeypatch.setattr(codec, "load_config_toml", counting)
    codec.with_stored_font_sizes_as_auto(_session(CLI_CONFIG)["renderRequest"])
    assert calls == []
    two_zeros = {"diagramOptions": {"config": {"labels": {"font_size": {"long": 0, "short": -1}}}}}
    codec.with_stored_font_sizes_as_auto(two_zeros)
    assert calls == [("gbdraw.data", "config.toml")]
