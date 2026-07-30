from __future__ import annotations

from pathlib import Path
from typing import Callable

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.cli_utils.session as cli_session
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.circular import _get_args as get_circular_args
from gbdraw.circular import circular_main
from gbdraw.cli_utils.session import parse_session_pre_args, strip_session_output_args
from gbdraw.exceptions import ValidationError
from gbdraw.linear import _get_args as get_linear_args
from gbdraw.linear import linear_main
from gbdraw.session import (
    build_session_document,
    save_session_document,
    with_request_output,
)


def _write_genbank(path: Path) -> None:
    record = SeqRecord(
        Seq("ATGC" * 25),
        id="overwrite-test",
        annotations={"molecule_type": "DNA"},
    )
    SeqIO.write(record, path, "genbank")


@pytest.mark.parametrize("get_args", [get_circular_args, get_linear_args])
def test_cli_overwrite_flag_is_explicit(get_args: Callable) -> None:
    assert get_args(["--gbk", "input.gb"]).overwrite is False
    assert get_args(["--gbk", "input.gb", "--overwrite"]).overwrite is True


@pytest.mark.parametrize("run_cli", [circular_main, linear_main])
def test_cli_refuses_existing_output_without_overwrite(
    run_cli: Callable[[list[str]], None],
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    output_path = output_prefix.with_suffix(".svg")
    _write_genbank(input_path)
    output_path.write_text("keep this file", encoding="utf-8")

    with pytest.raises(ValidationError, match="already exist"):
        run_cli(
            [
                "--gbk",
                str(input_path),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )

    assert output_path.read_text(encoding="utf-8") == "keep this file"


@pytest.mark.parametrize("run_cli", [circular_main, linear_main])
def test_cli_overwrite_replaces_existing_output(
    run_cli: Callable[[list[str]], None],
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    output_path = output_prefix.with_suffix(".svg")
    _write_genbank(input_path)
    output_path.write_text("stale output", encoding="utf-8")

    run_cli(
        [
            "--gbk",
            str(input_path),
            "--output",
            str(output_prefix),
            "--format",
            "svg",
            "--overwrite",
            "--no-gc",
            "--no-skew",
            "--legend",
            "none",
        ]
    )

    assert output_path.read_text(encoding="utf-8").startswith("<?xml")


@pytest.mark.parametrize("run_cli", [circular_main, linear_main])
def test_cli_rejects_existing_session_sidecar_before_rendering(
    run_cli: Callable[[list[str]], None],
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    output_path = output_prefix.with_suffix(".svg")
    sidecar_path = tmp_path / "saved.gbdraw-session.json"
    _write_genbank(input_path)
    sidecar_path.write_text("keep this session", encoding="utf-8")

    with pytest.raises(ValidationError, match="Session output already exists"):
        run_cli(
            [
                "--gbk",
                str(input_path),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(sidecar_path),
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )

    assert sidecar_path.read_text(encoding="utf-8") == "keep this session"
    assert not output_path.exists()


@pytest.mark.parametrize(
    ("run_cli", "request_type"),
    (
        (circular_main, CircularDiagramRequest),
        (linear_main, LinearDiagramRequest),
    ),
)
def test_canonical_replay_rejects_existing_sidecar_before_rendering(
    run_cli: Callable[[list[str]], None],
    request_type: type[CircularDiagramRequest] | type[LinearDiagramRequest],
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGC" * 25),
        id="replay-test",
        annotations={"molecule_type": "DNA"},
    )
    session_path = tmp_path / "source.gbdraw-session.json"
    output_prefix = tmp_path / "replayed"
    output_path = output_prefix.with_suffix(".svg")
    sidecar_path = tmp_path / "saved.gbdraw-session.json"
    request = request_type(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(output_prefix="stored"),
    )
    save_session_document(session_path, request)
    sidecar_path.write_text("keep this session", encoding="utf-8")

    with pytest.raises(ValidationError, match="Session output already exists"):
        run_cli(
            [
                "--session",
                str(session_path),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(sidecar_path),
            ]
        )

    assert sidecar_path.read_text(encoding="utf-8") == "keep this session"
    assert not output_path.exists()


def test_circular_implicit_sidecar_collision_is_rejected_before_rendering(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_path = tmp_path / "overwrite-test.svg"
    sidecar_path = tmp_path / "overwrite-test.gbdraw-session.json"
    _write_genbank(input_path)
    sidecar_path.write_text("keep this session", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValidationError, match="Session output already exists"):
        circular_main(
            [
                "--gbk",
                str(input_path),
                "--format",
                "svg",
                "--save_session",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )

    assert sidecar_path.read_text(encoding="utf-8") == "keep this session"
    assert not output_path.exists()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_session_replay_requires_current_overwrite_permission(mode: str) -> None:
    default_request = parse_session_pre_args(
        ["--session", "diagram.gbdraw-session.json"],
        mode=mode,
    )
    overwrite_request = parse_session_pre_args(
        ["--session", "diagram.gbdraw-session.json", "--overwrite"],
        mode=mode,
    )

    assert default_request is not None
    assert default_request.overwrite is False
    assert overwrite_request is not None
    assert overwrite_request.overwrite is True


def test_session_output_override_can_disable_stored_overwrite() -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(overwrite=True),
    )

    updated = with_request_output(request, overwrite=False)

    assert request.output.overwrite is True
    assert updated.output.overwrite is False


@pytest.mark.parametrize("overwrite", [False, True])
def test_canonical_session_replay_uses_current_overwrite_permission(
    overwrite: bool,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    stored_request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(overwrite=True),
    )
    session = build_session_document(stored_request).to_dict()
    captured: dict[str, bool] = {}

    def fake_render(request, *, session_artifacts=None):
        captured["overwrite"] = request.output.overwrite
        assert session_artifacts is not None
        return object()

    monkeypatch.setattr(cli_session, "_render_request", fake_render)

    assert cli_session.render_canonical_session_if_present(
        session,
        mode="circular",
        output_override=str(tmp_path / "replayed"),
        format_override="svg",
        overwrite=overwrite,
        save_session=False,
        session_output=None,
    )
    assert captured["overwrite"] is overwrite


def test_saved_cli_invocation_does_not_persist_overwrite_permission() -> None:
    assert strip_session_output_args(
        [
            "--gbk",
            "input.gb",
            "--overwrite",
            "--save_session",
            "--session_output",
            "diagram.gbdraw-session.json",
        ]
    ) == ["--gbk", "input.gb"]


def test_web_wrapper_explicitly_replaces_ephemeral_outputs() -> None:
    helper_path = (
        Path(__file__).resolve().parents[1]
        / "gbdraw"
        / "web"
        / "js"
        / "app"
        / "python-helpers.js"
    )

    assert 'full_args = args + ["-f", "svg", "--overwrite"]' in helper_path.read_text(
        encoding="utf-8"
    )
