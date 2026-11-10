"""The in-process ``gbdraw`` runner answers as a fresh ``python -m gbdraw.cli`` does and leaves the process as it was."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

from gbdraw.core import fonts
from tests.utils.cli_in_process import run_gbdraw

LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"


@pytest.mark.parametrize(
    "args",
    [
        pytest.param(["unknown"], id="unknown-command"),
        pytest.param(["linear", "--no_such_option"], id="parser-error"),
        pytest.param(["linear", "--gbk", str(LAMBDA), "-o", "out", "-f", "svg"], id="linear-figure"),
        pytest.param(["circular", "--gbk", str(LAMBDA), "-o", "out", "-f", "svg", "--labels"], id="circular-figure"),
        pytest.param(["linear", "--session", "missing.gbdraw-session.json", "-o", "out"], id="session-error"),
    ],
)
def test_in_process_run_answers_as_a_fresh_process(args: list[str], tmp_path: Path) -> None:
    process_dir = tmp_path / "process"
    in_process_dir = tmp_path / "in_process"
    process_dir.mkdir()
    in_process_dir.mkdir()
    process = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", *args], cwd=process_dir, capture_output=True, text=True, check=False
    )

    result = run_gbdraw(args, cwd=in_process_dir)

    assert (result.returncode, result.stdout, result.stderr) == (process.returncode, process.stdout, process.stderr)
    assert sorted(path.name for path in in_process_dir.iterdir()) == sorted(path.name for path in process_dir.iterdir())
    for path in process_dir.iterdir():
        assert (in_process_dir / path.name).read_bytes() == path.read_bytes(), path.name


def test_in_process_run_leaves_the_process_as_it_found_it(tmp_path: Path) -> None:
    root = logging.getLogger()
    before = (
        list(root.handlers), root.level, list(sys.argv), os.getcwd(), list(warnings.filters), warnings.showwarning,
        fonts.registered_fonts(),
    )

    result = run_gbdraw(["linear", "--gbk", str(LAMBDA), "-o", "out", "-f", "svg"], cwd=tmp_path)

    assert result.returncode == 0, result.stderr
    assert "INFO:" in result.stderr
    after = (
        list(root.handlers), root.level, list(sys.argv), os.getcwd(), list(warnings.filters), warnings.showwarning,
        fonts.registered_fonts(),
    )
    assert after == before
