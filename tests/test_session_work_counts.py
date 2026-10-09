"""Work counts of Session load, render and write (perf 0.14.x, PERF-P 6, 8, 10).

Each test counts calls of one piece of work during one operation and asserts
the exact number, so a second validation, copy or extraction of the same
object fails here.
"""

from __future__ import annotations

import base64
import copy
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.api.request_render  # noqa: F401 - imported so its names are patched
import gbdraw.api.session_compat  # noqa: F401
import gbdraw.circular  # noqa: F401
import gbdraw.cli_utils.session  # noqa: F401
import gbdraw.linear  # noqa: F401
import gbdraw.session  # noqa: F401
from gbdraw import cli
from gbdraw.api.request_render import _extract_linear_request_proteins
from gbdraw.api import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.session import save_session_document
from gbdraw.session_io import SESSION_FORMAT, validate_session, write_session_json

REPO_ROOT = Path(__file__).resolve().parents[1]
GALLERY_SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
CIRCULAR_SESSION = GALLERY_SESSIONS / "HmmtDNA_basic_circular.gbdraw-session.json"
# Five records with current LOSATP raw entries and a protein identity manifest.
PROTEIN_SESSION = GALLERY_SESSIONS / "BGC0000708-BGC0000713.gbdraw-session.json"


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


def _count_deepcopy_callers(monkeypatch: pytest.MonkeyPatch) -> Counter[str]:
    """Count ``copy.deepcopy`` calls made directly by gbdraw functions."""

    original = copy.deepcopy
    callers: Counter[str] = Counter()

    def counted(value: Any, memo: Any = None, _nil: Any = []) -> Any:
        frame = sys._getframe(1)
        module = str(frame.f_globals.get("__name__", ""))
        if module.startswith("gbdraw"):
            callers[f"{module}.{frame.f_code.co_name}"] += 1
        return original(value, memo)

    monkeypatch.setattr(copy, "deepcopy", counted)
    return callers


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


def test_cli_session_render_validates_and_copies_the_session_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    validations = _count_calls(monkeypatch, validate_session)
    copies = _count_deepcopy_callers(monkeypatch)
    _run_cli("circular", "--session", str(CIRCULAR_SESSION), "-o", str(tmp_path / "replayed"), "-f", "svg")
    # The Session is validated once, when it is loaded, and the loaded document
    # is not copied: the render reads it in place.
    assert len(validations) == 1
    assert copies == Counter(
        {
            "gbdraw.api.session_compat.canonical_payload_for_session_decode": 1,
            # CurrentRequestArtifacts detaches the Session's identity manifest.
            "gbdraw.api.request_render.__post_init__": 1,
        }
    )

    # Writing the Session again validates the loaded and the written document.
    validations.clear()
    copies.clear()
    _run_cli(
        "circular",
        "--session",
        str(CIRCULAR_SESSION),
        "-o",
        str(tmp_path / "resaved"),
        "-f",
        "svg",
        "--session_output",
        str(tmp_path / "resaved.gbdraw-session.json"),
    )
    assert len(validations) == 2
    assert copies == Counter(
        {
            "gbdraw.api.session_compat.canonical_payload_for_session_decode": 1,
            "gbdraw.api.request_render.__post_init__": 1,
            # The sidecar's source state, its adjunct, and the written document.
            "gbdraw.session.to_dict": 1,
            "gbdraw.session._build_session_document_from_resolved_request": 1,
            "gbdraw.session.__post_init__": 1,
        }
    )


def _legacy_session(tmp_path: Path, mode: str) -> Path:
    """A Session 30 with one GenBank file: the CLI replays it as CLI arguments."""

    record = SeqRecord(Seq("ATGC" * 90), id="legacy", annotations={"molecule_type": "DNA"})
    genbank = tmp_path / "legacy.gb"
    SeqIO.write([record], genbank, "genbank")
    content = genbank.read_bytes()
    embedded = {
        "name": "legacy.gb",
        "type": "application/octet-stream",
        "size": len(content),
        "lastModified": 0,
        "data": base64.b64encode(content).decode("ascii"),
    }
    session = {
        "format": SESSION_FORMAT,
        "version": 30,
        "createdAt": "2026-06-22T00:00:00Z",
        "config": {"form": {"prefix": "out"}, "adv": {}},
        "ui": {"mode": mode, "cInputType": "gb", "lInputType": "gb"},
        "files": {"c_gb": embedded} if mode == "circular" else {"linearSeqs": [{"gb": embedded}]},
    }
    path = tmp_path / "legacy.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    return path


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_legacy_session_replay_reads_the_loaded_payload(
    mode: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    session_path = _legacy_session(tmp_path, mode)
    validations = _count_calls(monkeypatch, validate_session)
    copies = _count_deepcopy_callers(monkeypatch)
    _run_cli(mode, "--session", str(session_path), "-o", str(tmp_path / "replayed"), "-f", "svg")
    assert (tmp_path / "replayed.svg").exists()
    # Loading validates the Session, and so does the public session_to_cli_args;
    # the linear run also validates the source Session once more when it renders.
    assert len(validations) == (2 if mode == "circular" else 3)
    # The replay reads the loaded payload; only the run copies its own input.
    expected = Counter({f"gbdraw.{mode}.run_{mode}_from_namespace": 1})
    if mode == "linear":
        # The linear run's CurrentRequestArtifacts detaches the identity manifest.
        expected["gbdraw.api.request_render.__post_init__"] = 1
    assert copies == expected


def test_cli_linear_session_extracts_its_proteins_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    extractions = _count_calls(monkeypatch, _extract_linear_request_proteins)
    _run_cli("linear", "--session", str(PROTEIN_SESSION), "-o", str(tmp_path / "replayed"), "-f", "svg")
    # The Session adapter extracts the proteins to check the saved LOSATP
    # entries; the build of the same plan reuses that extraction.
    assert len(extractions) == 1
