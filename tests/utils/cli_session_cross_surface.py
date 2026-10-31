"""CLI Session cross-surface matrix (D-01, OV-221).

A Session that the CLI writes draws the CLI's figure when the CLI, the Python
API, or the Web app loads it, and when one of them saves it again. The cases
live in ``tests/fixtures/cli_session_cross_surface/cases.json``; a case's
``legacySessions`` (OV-269) are its Sessions written by older CLI writers on
``main``, loaded in place of the current CLI Session. A case names the
``drawings`` it compares by ID (default: the one drawing its ``mode`` names);
a case with ``session`` and no ``args`` starts from that Session 49 fixture, and
the CLI figure of each of its drawings is ``gbdraw render --drawing ID``. pytest checks the CLI and
Python surfaces, and ``tests/web/session-cli-compatibility.test.mjs``
the Web app through ``tests/web/helpers/cli-session-cross-surface.py``.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from gbdraw.api import (
    load_session_document,
    materialize_session,
    render_session,
    save_session_document,
    session_to_request,
)
from tests.utils.svg_compare import SVGComparisonResult, compare_svgs

ROOT = Path(__file__).resolve().parents[2]
CASES_PATH = ROOT / "tests" / "fixtures" / "cli_session_cross_surface" / "cases.json"


@dataclass(frozen=True)
class CaseFiles:
    """The CLI figure of one case and the Sessions each writer saved for it."""

    svgs: dict[str, Path]
    cli_session: Path
    cli_resave: Path
    python_resaves: dict[str, Path]


def load_cases() -> list[dict]:
    """Every case, then each legacy Session of a case as a case of its own."""

    cases = json.loads(CASES_PATH.read_text(encoding="utf-8"))["cases"]
    return cases + [
        {**case, "id": f"{case['id']}@{name}", "session": fixture}
        for case in cases
        for name, fixture in case.get("legacySessions", {}).items()
    ]


def case_drawings(case: dict) -> list[str]:
    """The IDs of the drawings a case compares: its ``drawings``, else the one drawing of its mode."""

    return case.get("drawings") or [case["mode"]]


def _absolute(arg: str) -> str:
    return str(ROOT / arg) if arg.startswith("tests/") else arg


def run_cli(mode: str, args: list[str], output_prefix: Path, *extra: str) -> Path:
    """Run one CLI drawing and return its SVG."""

    result = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", mode, *map(_absolute, args),
         "-o", str(output_prefix), "-f", "svg", *extra],
        capture_output=True, text=True, cwd=output_prefix.parent,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
    )
    if result.returncode:
        raise AssertionError(f"gbdraw {mode} failed:\n{result.stdout}{result.stderr}")
    return output_prefix.with_suffix(".svg")


def replay_cli(mode: str, session: Path, output_prefix: Path, *extra: str, drawing: str | None = None) -> Path:
    """Replay a Session on the CLI: the mode command, or ``render --drawing ID`` for one drawing of it."""

    if drawing:
        return run_cli("render", ["--session", str(session), "--drawing", drawing], output_prefix, *extra)
    return run_cli(mode, ["--session", str(session)], output_prefix, *extra)


def render_python(session: Path, output_directory: Path, drawing: str | None = None) -> Path:
    """Draw ``drawing`` of a Session, by default the drawing it opens in Web."""

    document = load_session_document(session)
    with materialize_session(document, output_directory=output_directory) as materialized:
        result = render_session(materialized, drawing=drawing or document.active_drawing_id)
    return next(path for path in result.output_paths if path.suffix == ".svg")


def resave_python(session: Path, path: Path, drawing: str | None = None) -> Path:
    """Save one drawing of a Session again from Python (``drawing`` is needed when it holds several)."""

    with materialize_session(load_session_document(session), output_directory=path.parent / "python-resave-out") as materialized:
        save_session_document(path, session_to_request(materialized, drawing=drawing))
    return path


def write_case(case: dict, directory: Path) -> CaseFiles:
    """Draw a case on the CLI and save it again from the CLI and from Python.

    A legacy case loads its stored Session instead of the one the CLI writes
    now; a session-sourced case loads its Session 49 fixture and draws each of
    its drawings with ``gbdraw render --drawing ID``.
    """

    directory.mkdir(parents=True, exist_ok=True)
    cli_session = directory / "cli.gbdraw-session.json"
    cli_resave = directory / "cli-resave.gbdraw-session.json"
    stored = case.get("session")
    if "args" in case:
        svg = run_cli(case["mode"], case["args"], directory / "cli",
                      *(() if stored else ("--session_output", str(cli_session))))
        svgs = {case_drawings(case)[0]: svg}
    if stored:
        cli_session.write_bytes(gzip.decompress((ROOT / stored).read_bytes()))
    if "args" in case:
        replay_cli(case["mode"], cli_session, directory / "cli-replay", "--session_output", str(cli_resave))
    else:
        svgs = {drawing: replay_cli("render", cli_session, directory / f"cli-{drawing}", drawing=drawing)
                for drawing in case_drawings(case)}
        run_cli("render", ["--session", str(cli_session)], directory / "cli-replay", "--save_session")
        (directory / "cli-replay.gbdraw-session.json").rename(cli_resave)
    python_resaves = {drawing: resave_python(cli_session, directory / f"python-resave-{drawing}.gbdraw-session.json",
                                             drawing if "args" not in case else None)
                      for drawing in svgs}
    return CaseFiles(svgs, cli_session, cli_resave, python_resaves)


def check(expected: Path, actual: Path) -> SVGComparisonResult:
    return compare_svgs(expected, actual)
