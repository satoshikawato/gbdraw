"""D-04: a scale interval must be > 0 where a user types it; a stored value <= 0 still draws the automatic interval."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw import CircularOptions, LinearOptions
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    render_request,
)
from gbdraw.api.requests import GenBankInputSource
from gbdraw.exceptions import ValidationError

REPO_ROOT = Path(__file__).resolve().parents[1]
HUMAN = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
LEGACY_LINEAR_SESSION = REPO_ROOT / "tests" / "fixtures" / "sessions" / "cli-linear-protein.v30.gbdraw-session.json.gz"
CLI_MODULES = {"circular": circular_cli_module, "linear": linear_cli_module}


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("value", ["0", "-5"])
def test_cli_rejects_a_non_positive_scale_interval(
    mode: str, value: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        CLI_MODULES[mode]._get_args(["--gbk", str(HUMAN), "--scale_interval", value])
    assert excinfo.value.code == 2
    assert "--scale_interval must be > 0" in capsys.readouterr().err


@pytest.mark.parametrize("options_type", [CircularOptions, LinearOptions])
@pytest.mark.parametrize("value", [0, -5])
def test_python_api_rejects_a_non_positive_scale_interval(options_type: type, value: int) -> None:
    with pytest.raises(ValidationError, match=r"'objects\.scale\.interval'; expected a positive integer") as excinfo:
        options_type(config_overrides={"objects.scale.interval": value})
    assert excinfo.value.diagnostic == {
        "code": "INPUT_INVALID",
        "configPath": "objects.scale.interval",
        "reason": "POSITIVE_INTEGER_OR_AUTO",
    }
    assert options_type(config_overrides={"objects.scale.interval": 1}).config_overrides


def _typed_svg(mode: str, tmp_path: Path, interval: int | None) -> bytes:
    overrides = {} if interval is None else {"objects.scale.interval": interval}
    request_type, options_type = (
        (CircularDiagramRequest, CircularDiagramOptions)
        if mode == "circular"
        else (LinearDiagramRequest, LinearDiagramOptions)
    )
    prefix = f"{mode}-{interval}"
    render_request(
        request_type(
            records=(RecordInput(source=GenBankInputSource(path=str(HUMAN))),),
            options=options_type(config_overrides=overrides),
            output=RenderOutputRequest(output_prefix=prefix, output_directory=tmp_path, formats=("svg",)),
        )
    )
    return (tmp_path / f"{prefix}.svg").read_bytes()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_stored_non_positive_interval_in_a_request_draws_the_automatic_interval(
    mode: str, tmp_path: Path
) -> None:
    # The typed request is also what a Session decodes to, so it keeps accepting <= 0.
    automatic = _typed_svg(mode, tmp_path, None)
    assert _typed_svg(mode, tmp_path, 0) == automatic
    assert _typed_svg(mode, tmp_path, -5) == automatic


def test_a_legacy_cli_session_with_a_non_positive_interval_still_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    session = json.loads(gzip.decompress(LEGACY_LINEAR_SESSION.read_bytes()))
    monkeypatch.chdir(tmp_path)
    outputs = {}
    for name, extra in (("automatic", []), ("zero", ["--scale_interval", "0"]), ("negative", ["--scale_interval", "-5"])):
        stored = json.loads(json.dumps(session))
        args = stored["cliInvocation"]["args"]
        index = args.index("-o")
        stored["cliInvocation"]["args"] = [*args[:index], *extra, *args[index:]]
        path = tmp_path / f"{name}.gbdraw-session.json"
        path.write_text(json.dumps(stored), encoding="utf-8")
        linear_cli_module.linear_main(["--session", str(path), "-o", name, "-f", "svg"])
        outputs[name] = (tmp_path / f"{name}.svg").read_bytes()
    assert outputs["zero"] == outputs["automatic"]
    assert outputs["negative"] == outputs["automatic"]
