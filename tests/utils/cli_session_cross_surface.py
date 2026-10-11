"""CLI Session cross-surface matrix (D-01, OV-221).

A Session that the CLI writes draws the CLI's figure when the CLI, the Python
API, or the Web app loads it, and when one of them saves it again. The cases
live in ``tests/fixtures/cli_session_cross_surface/cases.json``; a case's
``legacySessions`` (OV-269) are its Sessions written by older CLI writers on
``main``, loaded in place of the current CLI Session. pytest checks the CLI and
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

    svg: Path
    cli_session: Path
    cli_resave: Path
    python_resave: Path


def load_cases() -> list[dict]:
    """Every case, then each legacy Session of a case as a case of its own."""

    cases = json.loads(CASES_PATH.read_text(encoding="utf-8"))["cases"]
    return cases + [
        {**case, "id": f"{case['id']}@{name}", "session": fixture}
        for case in cases
        for name, fixture in case.get("legacySessions", {}).items()
    ]


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


def replay_cli(mode: str, session: Path, output_prefix: Path, *extra: str) -> Path:
    return run_cli(mode, ["--session", str(session)], output_prefix, *extra)


def render_python(session: Path, output_directory: Path) -> Path:
    with materialize_session(load_session_document(session), output_directory=output_directory) as materialized:
        result = render_session(materialized)
    return next(path for path in result.output_paths if path.suffix == ".svg")


def resave_python(session: Path, path: Path) -> Path:
    with materialize_session(load_session_document(session), output_directory=path.parent / "python-resave-out") as materialized:
        save_session_document(path, session_to_request(materialized))
    return path


def write_case(case: dict, directory: Path) -> CaseFiles:
    """Draw a case on the CLI and save it again from the CLI and from Python.

    A legacy case loads its stored Session instead of the one the CLI writes now.
    """

    directory.mkdir(parents=True, exist_ok=True)
    cli_session = directory / "cli.gbdraw-session.json"
    legacy = case.get("session")
    svg = run_cli(case["mode"], case["args"], directory / "cli",
                  *(() if legacy else ("--session_output", str(cli_session))))
    if legacy:
        cli_session.write_bytes(gzip.decompress((ROOT / legacy).read_bytes()))
    cli_resave = directory / "cli-resave.gbdraw-session.json"
    replay_cli(case["mode"], cli_session, directory / "cli-replay", "--session_output", str(cli_resave))
    return CaseFiles(svg, cli_session, cli_resave, resave_python(cli_session, directory / "python-resave.gbdraw-session.json"))


def check(expected: Path, actual: Path) -> SVGComparisonResult:
    return compare_svgs(expected, actual)
