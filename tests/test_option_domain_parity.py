"""Shared option-domain vectors: the CLI and the typed request reject alike."""

from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path

import pytest

from gbdraw import cli
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.session_request_codec import (
    CanonicalRequestDecodingError,
    decode_canonical_request,
    encode_canonical_request,
)
from gbdraw.web_support.error_adapter import FIELDS, serialize_web_error

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = json.loads(
    (REPO_ROOT / "tests" / "fixtures" / "option_domain_vectors.json").read_text(
        encoding="utf-8"
    )
)
CLI_RECORD = REPO_ROOT / VECTORS["cliRecord"]


def _cases(kind: str, *, cli_only: bool = False):
    for vector in VECTORS[kind]:
        if cli_only and "cli" not in vector:
            continue
        for mode in vector["modes"]:
            yield pytest.param(mode, vector, id=f"{vector['id']}-{mode}")


def _request_value(value: object) -> object:
    return float("nan") if value == "NaN" else value


def _decode_with_options(tmp_path: Path, mode: str, options: dict) -> None:
    source = tmp_path / "record.gbk"
    source.write_text(CLI_RECORD.read_text(encoding="utf-8"), encoding="utf-8")
    request_type = CircularDiagramRequest if mode == "circular" else LinearDiagramRequest
    encoded = encode_canonical_request(
        request_type(records=(RecordInput(source=GenBankInputSource(source)),))
    )
    payload = copy.deepcopy(encoded.payload)
    diagram_options = payload["diagramOptions"]
    for key, value in options.items():
        if key == "configOverrides":
            merged = dict(diagram_options.get("configOverrides") or {})
            merged.update(value)
            diagram_options[key] = merged
        else:
            diagram_options[key] = _request_value(value)
    resource_paths = {
        resource.resource_id: resource.source_path
        for resource in encoded.resources
        if resource.source_path is not None
    }
    decode_canonical_request(
        payload,
        resource_paths=resource_paths,
        output_directory=tmp_path / "output",
    )


def _run_cli(monkeypatch, capsys, tmp_path: Path, mode: str, extra: list[str]):
    argv = [
        "gbdraw",
        mode,
        "--gbk",
        str(CLI_RECORD),
        "-o",
        str(tmp_path / "out"),
        "-f",
        "svg",
        *extra,
    ]
    monkeypatch.setattr(sys, "argv", argv)
    try:
        cli.main()
    except SystemExit as exit_:
        code = exit_.code
    else:
        code = 0
    return code, capsys.readouterr().err


@pytest.mark.parametrize(("mode", "vector"), list(_cases("invalid")))
def test_typed_request_rejects_vector_with_shared_reason(tmp_path, mode, vector):
    with pytest.raises(CanonicalRequestDecodingError) as caught:
        _decode_with_options(tmp_path, mode, vector["request"])
    expected = {
        key: vector[key]
        for key in ("code", "field", "reason", "configPath")
        if key in vector
    }
    # (ii) the producer's diagnostic in the decoding error chain.
    cause: BaseException | None = caught.value
    while cause is not None and getattr(cause, "diagnostic", None) is None:
        cause = cause.__cause__
    assert cause is not None and cause.diagnostic == expected
    # (iii) the Web payload keeps the published identifiers of that meaning.
    payload = serialize_web_error(caught.value, operation="generate", stage="request-validation")
    assert payload["code"] == vector["code"]
    assert payload["context"] == {
        key: value
        for key, value in expected.items()
        if key != "code" and (key != "field" or value in FIELDS)
    }


@pytest.mark.parametrize(("mode", "vector"), list(_cases("accepted")))
def test_typed_request_accepts_vector(tmp_path, mode, vector):
    _decode_with_options(tmp_path, mode, vector["request"])


@pytest.mark.parametrize(("mode", "vector"), list(_cases("invalid", cli_only=True)))
def test_cli_rejects_vector_without_traceback(monkeypatch, capsys, tmp_path, mode, vector):
    code, stderr = _run_cli(monkeypatch, capsys, tmp_path, mode, vector["cli"])
    assert code not in (0, None), stderr
    assert "Traceback" not in stderr
    assert not (tmp_path / "out.svg").exists()
    if vector.get("cliRejectedBy") == "argparse":
        assert "error:" in stderr
    else:
        # The typed owner's message names the rejected field or setting.
        assert "ERROR:" in stderr or "error:" in stderr
        path = vector.get("configPath", "")
        named = vector.get("field") or ("font_size" if "font_size" in path else "stroke_width")
        assert named in stderr, stderr


@pytest.mark.parametrize(("mode", "vector"), list(_cases("accepted", cli_only=True)))
def test_cli_accepts_vector(monkeypatch, capsys, tmp_path, mode, vector):
    code, stderr = _run_cli(monkeypatch, capsys, tmp_path, mode, vector["cli"])
    assert code in (0, None), stderr
    assert (tmp_path / "out.svg").is_file()


def test_vectors_cover_each_validated_request_field():
    covered = {vector.get("field") for vector in VECTORS["invalid"]}
    assert {
        "window",
        "step",
        "depth_window",
        "depth_step",
        "dinucleotide",
        "evalue",
        "bitscore",
        "identity",
        "alignment_length",
        "plot_title_font_size",
    } <= covered
    assert all(
        not math.isnan(value) if isinstance(value, float) else True
        for vector in VECTORS["invalid"]
        for value in vector["request"].values()
    )
