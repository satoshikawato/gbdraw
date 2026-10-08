"""D-04: a scale interval must be > 0 where a user enters one; a stored value <= 0 still draws the automatic interval."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw import CircularOptions, LinearOptions, draw_circular, draw_linear, read_genbank
from gbdraw.api import (
    CircularDiagramOptions,
    LinearDiagramOptions,
    load_session_document,
    materialize_session,
    render_session,
)
from gbdraw.api.config import load_default_config
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.web_support.error_adapter import serialize_web_error

REPO_ROOT = Path(__file__).resolve().parents[1]
HUMAN = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
# Real writers, see scale-interval.provenance.json: main fe6861f0 CLI (diagramOptions.config),
# main fe6861f0 typed API (diagramOptions.configOverrides), and the 0.13.0 CLI (legacy cliInvocation).
MAIN_CLI_CIRCULAR_ZERO = SESSIONS / "scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz"
MAIN_API_LINEAR_NEGATIVE = SESSIONS / "scale-interval-negative-linear-api.v44.gbdraw-session.json.gz"
RELEASED_CLI_LINEAR_NEGATIVE = SESSIONS / "scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz"
RELEASED_CLI_CIRCULAR_ZERO = SESSIONS / "scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz"
CLI_MODULES = {"circular": circular_cli_module, "linear": linear_cli_module}
EXPECTED_DIAGNOSTIC = {
    "code": "INPUT_INVALID",
    "configPath": "objects.scale.interval",
    "reason": "POSITIVE_INTEGER_OR_AUTO",
}


def _session(path: Path) -> dict[str, Any]:
    return json.loads(gzip.decompress(path.read_bytes()))


def _without_stored_interval(session: dict[str, Any]) -> dict[str, Any]:
    stored = json.loads(json.dumps(session))
    options = stored["renderRequest"]["diagramOptions"]
    if options.get("config") is not None:
        options["config"]["objects"]["scale"]["interval"] = None
    if options.get("configOverrides") is not None:
        options["configOverrides"]["objects.scale.interval"] = None
    return stored


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("value", ["0", "-5"])
def test_cli_rejects_a_non_positive_scale_interval(
    mode: str, value: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        CLI_MODULES[mode]._get_args(["--gbk", str(HUMAN), "--scale_interval", value])
    assert excinfo.value.code == 2
    assert "--scale_interval must be > 0" in capsys.readouterr().err


@pytest.mark.parametrize("options_type", [CircularDiagramOptions, LinearDiagramOptions])
@pytest.mark.parametrize("value", [0, -5])
def test_typed_api_rejects_a_non_positive_scale_interval(options_type: type, value: int) -> None:
    with pytest.raises(ValidationError, match=r"'objects\.scale\.interval'") as excinfo:
        options_type(config_overrides={"objects.scale.interval": value})
    assert excinfo.value.diagnostic == EXPECTED_DIAGNOSTIC
    config = load_default_config()
    config["objects"]["scale"]["interval"] = value
    with pytest.raises(ValidationError, match=r"'objects\.scale\.interval'") as excinfo:
        options_type(config=config)
    assert excinfo.value.diagnostic == EXPECTED_DIAGNOSTIC
    assert options_type(config_overrides={"objects.scale.interval": 1}).config_overrides


@pytest.mark.parametrize(
    ("draw", "options_type"), [(draw_circular, CircularOptions), (draw_linear, LinearOptions)]
)
def test_python_api_rejects_a_non_positive_scale_interval(draw: Any, options_type: type) -> None:
    records = read_genbank(HUMAN)
    for value in (0, -5):
        with pytest.raises(ValidationError, match=r"'objects\.scale\.interval'") as excinfo:
            draw(records, options=options_type(config_overrides={"objects.scale.interval": value}))
        assert excinfo.value.diagnostic == EXPECTED_DIAGNOSTIC
        config = load_default_config()
        config["objects"]["scale"]["interval"] = value
        with pytest.raises(ValidationError, match=r"'objects\.scale\.interval'") as excinfo:
            draw(records, options=options_type(config=config))
        assert excinfo.value.diagnostic == EXPECTED_DIAGNOSTIC


def test_a_web_generate_request_rejects_what_a_session_load_reads_as_automatic(tmp_path: Path) -> None:
    # Web Generate decodes its request without the Session read rules.
    with materialize_session(
        load_session_document(_session(MAIN_API_LINEAR_NEGATIVE)), output_directory=tmp_path
    ) as materialized:
        with pytest.raises(Exception, match=r"objects\.scale\.interval") as excinfo:
            decode_canonical_request(
                materialized.document._data["renderRequest"],
                resource_paths=materialized.resource_paths,
                output_directory=tmp_path,
            )
    web_error = serialize_web_error(excinfo.value, operation="render", stage="decode")
    assert web_error["code"] == "INPUT_INVALID"
    assert web_error["context"] == {"reason": "POSITIVE_INTEGER_OR_AUTO", "configPath": "objects.scale.interval"}


def _python_replay_svg(session: dict[str, Any], tmp_path: Path, name: str) -> bytes:
    output_directory = tmp_path / name
    with materialize_session(load_session_document(session), output_directory=output_directory) as materialized:
        render_session(materialized)
    (svg,) = output_directory.glob("*.svg")
    return svg.read_bytes()


def _cli_replay_svg(session: dict[str, Any], tmp_path: Path, name: str) -> bytes:
    path = tmp_path / f"{name}.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    mode = (session.get("renderRequest") or session["cliInvocation"])["mode"]
    main = getattr(CLI_MODULES[mode], f"{mode}_main")
    main(["--session", str(path), "-o", name, "-f", "svg"])
    return (tmp_path / f"{name}.svg").read_bytes()


@pytest.mark.parametrize("fixture", [MAIN_CLI_CIRCULAR_ZERO, MAIN_API_LINEAR_NEGATIVE], ids=["cli-config", "api-overrides"])
def test_a_main_session_with_a_non_positive_interval_draws_the_automatic_interval(
    fixture: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    session = _session(fixture)
    automatic = _without_stored_interval(session)
    assert _python_replay_svg(session, tmp_path, "stored") == _python_replay_svg(automatic, tmp_path, "automatic")
    assert _cli_replay_svg(session, tmp_path, "cli-stored") == _cli_replay_svg(automatic, tmp_path, "cli-automatic")


@pytest.mark.parametrize(
    ("fixture", "stored"),
    [(RELEASED_CLI_LINEAR_NEGATIVE, "-5"), (RELEASED_CLI_CIRCULAR_ZERO, "0")],
    ids=["linear-negative", "circular-zero"],
)
def test_a_released_legacy_cli_session_with_a_non_positive_interval_still_replays(
    fixture: Path, stored: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    session = _session(fixture)
    assert session["version"] == 30
    automatic = json.loads(json.dumps(session))
    args = automatic["cliInvocation"]["args"]
    index = args.index("--scale_interval")
    assert args[index + 1] == stored
    del args[index : index + 2]
    assert _cli_replay_svg(session, tmp_path, "stored") == _cli_replay_svg(automatic, tmp_path, "automatic")
