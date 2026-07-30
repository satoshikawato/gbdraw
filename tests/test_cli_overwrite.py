from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Callable

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.circular as circular_module
import gbdraw.cli_utils.session as cli_session
import gbdraw.session as session_module
from gbdraw.api.requests import (
    CircularBatchRequest,
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
    SessionFormatError,
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


@pytest.mark.parametrize("run_cli", [circular_main, linear_main])
@pytest.mark.parametrize("overwrite", [False, True])
def test_plain_cli_rejects_diagram_and_sidecar_path_collision_before_rendering(
    run_cli: Callable[[list[str]], None],
    overwrite: bool,
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    shared_path = output_prefix.with_suffix(".svg")
    _write_genbank(input_path)
    shared_path.write_text("keep this file", encoding="utf-8")
    args = [
        "--gbk",
        str(input_path),
        "--output",
        str(output_prefix),
        "--format",
        "svg",
        "--session_output",
        str(shared_path),
        "--no-gc",
        "--no-skew",
        "--legend",
        "none",
    ]
    if overwrite:
        args.append("--overwrite")

    with pytest.raises(ValidationError, match="collides with diagram output"):
        run_cli(args)

    assert shared_path.read_text(encoding="utf-8") == "keep this file"


def test_circular_sidecar_collision_is_rejected_before_diagram_build(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    _write_genbank(input_path)

    def unexpected_build(*_args, **_kwargs):
        raise AssertionError("diagram build started before sidecar preflight")

    monkeypatch.setattr(
        circular_module,
        "build_request_plan_diagram",
        unexpected_build,
    )

    with pytest.raises(ValidationError, match="collides with diagram output"):
        circular_main(
            [
                "--gbk",
                str(input_path),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(output_prefix.with_suffix(".svg")),
                "--overwrite",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )


def test_circular_batch_output_collision_with_sidecar_is_rejected_before_build(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    first_input = tmp_path / "first.gb"
    second_input = tmp_path / "second.gb"
    output_prefix = tmp_path / "batch"
    occupied_output = tmp_path / "batch_2.svg"
    _write_genbank(first_input)
    _write_genbank(second_input)
    occupied_output.write_text("keep this file", encoding="utf-8")

    def unexpected_build(*_args, **_kwargs):
        raise AssertionError("diagram build started before batch output preflight")

    monkeypatch.setattr(
        circular_module,
        "build_request_plan_diagram",
        unexpected_build,
    )

    with pytest.raises(ValidationError, match="already exist"):
        circular_main(
            [
                "--gbk",
                str(first_input),
                str(second_input),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(tmp_path / "sidecar.json"),
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )

    assert not (tmp_path / "batch_1.svg").exists()
    assert occupied_output.read_text(encoding="utf-8") == "keep this file"


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


@pytest.mark.parametrize(
    ("run_cli", "request_type"),
    (
        (circular_main, CircularDiagramRequest),
        (linear_main, LinearDiagramRequest),
    ),
)
@pytest.mark.parametrize("overwrite", [False, True])
def test_canonical_replay_rejects_diagram_and_sidecar_path_collision(
    run_cli: Callable[[list[str]], None],
    request_type: type[CircularDiagramRequest] | type[LinearDiagramRequest],
    overwrite: bool,
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGC" * 25),
        id="replay-test",
        annotations={"molecule_type": "DNA"},
    )
    session_path = tmp_path / "source.gbdraw-session.json"
    output_prefix = tmp_path / "replayed"
    shared_path = output_prefix.with_suffix(".svg")
    save_session_document(
        session_path,
        request_type(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            output=RenderOutputRequest(output_prefix="stored"),
        ),
    )
    shared_path.write_text("keep this file", encoding="utf-8")
    args = [
        "--session",
        str(session_path),
        "--output",
        str(output_prefix),
        "--format",
        "svg",
        "--session_output",
        str(shared_path),
    ]
    if overwrite:
        args.append("--overwrite")

    with pytest.raises(ValidationError, match="collides with diagram output"):
        run_cli(args)

    assert shared_path.read_text(encoding="utf-8") == "keep this file"


@pytest.mark.parametrize("overwrite", [False, True])
def test_canonical_circular_batch_rejects_sidecar_collision_with_any_item(
    overwrite: bool,
    tmp_path: Path,
) -> None:
    records = tuple(
        SeqRecord(
            Seq("ATGC" * 25),
            id=record_id,
            annotations={"molecule_type": "DNA"},
        )
        for record_id in ("first", "second")
    )
    session_path = tmp_path / "batch.gbdraw-session.json"
    save_session_document(
        session_path,
        CircularBatchRequest(
            records=tuple(
                RecordInput(source=InMemoryRecordSource(record))
                for record in records
            ),
            outputs=(
                RenderOutputRequest(output_prefix="stored-1"),
                RenderOutputRequest(output_prefix="stored-2"),
            ),
        ),
    )
    output_prefix = tmp_path / "batch-replay"
    first_path = tmp_path / "batch-replay_1.svg"
    shared_path = tmp_path / "batch-replay_2.svg"
    shared_path.write_text("keep this file", encoding="utf-8")
    args = [
        "--session",
        str(session_path),
        "--output",
        str(output_prefix),
        "--format",
        "svg",
        "--session_output",
        str(shared_path),
    ]
    if overwrite:
        args.append("--overwrite")

    with pytest.raises(ValidationError, match="collides with diagram output"):
        circular_main(args)

    assert not first_path.exists()
    assert shared_path.read_text(encoding="utf-8") == "keep this file"


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


def test_circular_explicit_output_implicit_sidecar_preflights_actual_path(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "explicit-output"
    output_path = output_prefix.with_suffix(".svg")
    sidecar_path = Path(f"{output_prefix}.gbdraw-session.json")
    _write_genbank(input_path)
    sidecar_path.write_text("keep this session", encoding="utf-8")

    with pytest.raises(ValidationError, match="Session output already exists"):
        circular_main(
            [
                "--gbk",
                str(input_path),
                "--output",
                str(output_prefix),
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


@pytest.mark.parametrize(
    "target_kind",
    ("directory", "fifo", "file-parent", "dangling-symlink"),
)
def test_circular_rejects_unwritable_sidecar_target_before_diagram_build(
    target_kind: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.gb"
    output_prefix = tmp_path / "diagram"
    _write_genbank(input_path)
    if target_kind == "directory":
        sidecar_path = tmp_path / "sidecar-target"
        sidecar_path.mkdir()
        expected = "not replaceable files"
    elif target_kind == "fifo":
        if not hasattr(os, "mkfifo"):
            pytest.skip("FIFOs are unavailable on this platform")
        sidecar_path = tmp_path / "sidecar-target"
        os.mkfifo(sidecar_path)
        expected = "not replaceable files"
    elif target_kind == "file-parent":
        blocked_parent = tmp_path / "blocked-parent"
        blocked_parent.write_text("not a directory", encoding="utf-8")
        sidecar_path = blocked_parent / "sidecar.json"
        expected = "parent path.*not directories"
    else:
        sidecar_path = tmp_path / "sidecar-target"
        try:
            sidecar_path.symlink_to(tmp_path / "missing-sidecar-target")
        except OSError as exc:
            pytest.skip(f"file symlinks are unavailable: {exc}")
        expected = "not replaceable files"

    def unexpected_build(*_args, **_kwargs):
        raise AssertionError("diagram build started before sidecar target preflight")

    monkeypatch.setattr(
        circular_module,
        "build_request_plan_diagram",
        unexpected_build,
    )

    with pytest.raises(ValidationError, match=expected):
        circular_main(
            [
                "--gbk",
                str(input_path),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(sidecar_path),
                "--overwrite",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ]
        )

    assert not output_prefix.with_suffix(".svg").exists()
    if target_kind == "dangling-symlink":
        assert sidecar_path.is_symlink()
        assert not sidecar_path.exists()


def test_save_session_document_refuses_dangling_symlink_without_overwrite(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "saved.gbdraw-session.json"
    missing_target = tmp_path / "missing-session-target"
    try:
        output_path.symlink_to(missing_target)
    except OSError as exc:
        pytest.skip(f"file symlinks are unavailable: {exc}")
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(
                    SeqRecord(
                        Seq("ATGC"),
                        id="record",
                        annotations={"molecule_type": "DNA"},
                    )
                )
            ),
        ),
    )

    with pytest.raises(SessionFormatError, match="not replaceable files"):
        save_session_document(output_path, request)

    assert output_path.is_symlink()
    assert not missing_target.exists()


@pytest.mark.parametrize(
    "target_kind",
    ("directory", "fifo", "dangling-symlink"),
)
def test_save_session_document_rejects_nonreplaceable_target_with_overwrite(
    target_kind: str,
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "saved.gbdraw-session.json"
    if target_kind == "directory":
        output_path.mkdir()
    elif target_kind == "fifo":
        if not hasattr(os, "mkfifo"):
            pytest.skip("FIFOs are unavailable on this platform")
        os.mkfifo(output_path)
    else:
        try:
            output_path.symlink_to(tmp_path / "missing-session-target")
        except OSError as exc:
            pytest.skip(f"file symlinks are unavailable: {exc}")
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(
                    SeqRecord(
                        Seq("ATGC"),
                        id="record",
                        annotations={"molecule_type": "DNA"},
                    )
                )
            ),
        ),
    )

    with pytest.raises(SessionFormatError, match="not replaceable files"):
        save_session_document(output_path, request, overwrite=True)

    if target_kind == "directory":
        assert output_path.is_dir()
    elif target_kind == "fifo":
        assert output_path.exists()
    else:
        assert output_path.is_symlink()


@pytest.mark.parametrize("name", ("side\ncar.json", "side\x00car.json"))
def test_session_sidecar_rejects_control_characters(
    name: str,
    tmp_path: Path,
) -> None:
    with pytest.raises(
        ValidationError,
        match="must not contain ASCII control characters",
    ):
        cli_session.preflight_session_sidecar_if_requested(
            save_session=True,
            session_output=str(tmp_path / name),
            output_prefix=None,
        )


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

    def fake_render(request, *, session_document=None):
        captured["overwrite"] = request.output.overwrite
        assert session_document is not None
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


def test_legacy_canonical_sidecar_saves_rendered_request_and_migrated_adjunct(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    stored_request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(output_prefix="stored"),
    )
    session = build_session_document(stored_request).to_dict()
    session["version"] = 33
    session["config"] = {
        "adv": {
            "depth_tick_interval": 10,
            "depth_tracks": [{"tick_interval": 5}],
        },
        "losat": {"blastp": {"collinearMaxGeneGap": 2}},
    }
    session["ui"] = {"mode": "linear"}
    captured: dict[str, object] = {}

    def fake_render(request, *, session_document=None):
        assert session_document is not None
        rendered_request = replace(
            request,
            output=replace(request.output, output_prefix="adapter-owned"),
        )
        captured["rendered_request"] = rendered_request
        return SimpleNamespace(
            request=rendered_request,
            drawing=object(),
            output_paths=(),
            interactive_context=None,
        )

    def fake_save(path, request, **kwargs):
        captured["saved_path"] = path
        captured["saved_request"] = request
        captured["adjunct"] = kwargs["adjunct"]

    monkeypatch.setattr(cli_session, "_render_request", fake_render)
    monkeypatch.setattr(session_module, "save_session_document", fake_save)
    sidecar_path = tmp_path / "saved.gbdraw-session.json"

    assert cli_session.render_canonical_session_if_present(
        session,
        mode="linear",
        output_override=str(tmp_path / "replayed"),
        format_override="svg",
        overwrite=False,
        save_session=True,
        session_output=str(sidecar_path),
    )

    assert captured["saved_path"] == sidecar_path
    assert captured["saved_request"] is captured["rendered_request"]
    adjunct = captured["adjunct"]
    assert isinstance(adjunct, dict)
    config = adjunct["config"]
    assert config["adv"]["depth_large_tick_interval"] == 10
    assert config["adv"]["depth_tracks"] == [{"large_tick_interval": 5}]
    assert config["losat"]["blastp"]["collinearMaxUnitGap"] == 2
    assert "depth_tick_interval" in session["config"]["adv"]
    assert "collinearMaxGeneGap" in session["config"]["losat"]["blastp"]


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


def test_web_render_isolates_ephemeral_outputs_without_overwrite_permission() -> None:
    worker_path = (
        Path(__file__).resolve().parents[1]
        / "gbdraw"
        / "web"
        / "js"
        / "workers"
        / "diagram-generation-worker.js"
    )
    request_path = (
        Path(__file__).resolve().parents[1]
        / "gbdraw"
        / "web"
        / "js"
        / "services"
        / "session-request.js"
    )

    worker_source = worker_path.read_text(encoding="utf-8")
    request_source = request_path.read_text(encoding="utf-8")

    assert "const workspace = `/gbdraw-web-render-${Number(requestId) || 0}`" in worker_source
    assert "pyodide.FS.analyzePath(workspace).exists" in worker_source
    assert "Diagram render workspace cleanup invariant failed" in worker_source
    assert "overwrite: false" in request_source
