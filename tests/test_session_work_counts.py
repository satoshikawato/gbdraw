"""Work counts of Session load, render and write (perf 0.14.x, PERF-P 6, 8, 10).

Each test counts calls of one piece of work during one operation and asserts
the exact number, so a second validation, copy or extraction of the same
object fails here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Callable

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.api.request_render  # noqa: F401 - imported so its names are patched
import gbdraw.api.session_compat  # noqa: F401
import gbdraw.circular  # noqa: F401
import gbdraw.cli_utils.session  # noqa: F401
import gbdraw.linear  # noqa: F401
import gbdraw.session  # noqa: F401
from gbdraw import cli
from gbdraw.api import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.session import save_session_document
from gbdraw.session_io import validate_session, write_session_json

REPO_ROOT = Path(__file__).resolve().parents[1]
GALLERY_SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
CIRCULAR_SESSION = GALLERY_SESSIONS / "HmmtDNA_basic_circular.gbdraw-session.json"


def _count_calls(
    monkeypatch: pytest.MonkeyPatch,
    function: Callable[..., Any],
) -> list[tuple[Any, ...]]:
    """Patch every gbdraw module name bound to ``function`` with a counter."""

    calls: list[tuple[Any, ...]] = []

    def counted(*args: Any, **kwargs: Any) -> Any:
        calls.append(args)
        return function(*args, **kwargs)

    for module in list(sys.modules.values()):
        if not getattr(module, "__name__", "").startswith("gbdraw"):
            continue
        for name, value in list(vars(module).items()):
            if value is function:
                monkeypatch.setattr(module, name, counted)
    return calls


def _run_cli(*args: str) -> None:
    argv = sys.argv
    sys.argv = ["gbdraw", *args]
    try:
        cli.main()
    finally:
        sys.argv = argv


def _small_request(tmp_path: Path) -> CircularDiagramRequest:
    record = SeqRecord(Seq("ATGC" * 50), id="record", annotations={"molecule_type": "DNA"})
    return CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
    )


def test_session_writers_validate_each_payload_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    validations = _count_calls(monkeypatch, validate_session)

    # A built document was validated when it was built; its write does not
    # validate it again.
    save_session_document(tmp_path / "built.gbdraw-session.json", _small_request(tmp_path))
    assert len(validations) == 1

    # A payload from outside is validated by the writer.
    validations.clear()
    payload = json.loads((tmp_path / "built.gbdraw-session.json").read_text(encoding="utf-8"))
    write_session_json(tmp_path / "external.gbdraw-session.json", payload)
    assert len(validations) == 1
    assert (tmp_path / "external.gbdraw-session.json").read_bytes() == (
        tmp_path / "built.gbdraw-session.json"
    ).read_bytes()


def test_cli_save_session_validates_the_built_sidecar_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    validations = _count_calls(monkeypatch, validate_session)
    _run_cli(
        "circular",
        "--gbk",
        str(REPO_ROOT / "examples" / "MellatMJNV.gb"),
        "-o",
        str(tmp_path / "saved"),
        "-f",
        "svg",
        "--save_session",
    )
    # The canonical document inside the sidecar, then the sidecar that
    # build_session_json returns; the write validates neither again.
    assert len(validations) == 2
    assert (tmp_path / "saved.gbdraw-session.json").is_file()
