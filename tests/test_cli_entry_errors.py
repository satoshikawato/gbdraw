"""How the ``gbdraw`` command reports help, its version, and input errors.

Each case runs the command in a subprocess, as a user does, and checks what it
prints and its exit status (release-0140 stabilization: OV-1239, OV-1240,
OV-1241, OV-1251, OV-1252).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

TEST_INPUTS = Path(__file__).parent / "test_inputs"


def _gbdraw(*args: str, cwd: Path, timeout: int = 120) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", *args],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=cwd,
    )


@pytest.mark.parametrize("flag", ["--help", "-h"])
def test_gui_help_prints_usage_without_starting_the_server(flag: str, tmp_path: Path) -> None:
    """OV-1239: ``gbdraw gui --help`` started the server and the browser."""

    result = _gbdraw("gui", flag, cwd=tmp_path, timeout=30)
    assert result.returncode == 0, result.stderr
    assert "usage: gbdraw gui" in result.stdout
    assert "Local server started" not in result.stdout


def test_gui_refuses_unknown_arguments(tmp_path: Path) -> None:
    result = _gbdraw("gui", "--port", "8000", cwd=tmp_path, timeout=30)
    assert result.returncode == 2
    assert "unrecognized arguments: --port 8000" in result.stderr
    assert "Local server started" not in result.stdout


def test_version_is_a_top_level_flag_only(tmp_path: Path) -> None:
    """OV-1239: ``--version`` anywhere in the arguments printed the version and exited 0."""

    top = _gbdraw("--version", cwd=tmp_path)
    assert top.returncode == 0 and top.stdout.startswith("gbdraw version ")

    sub = _gbdraw("circular", "--gbk", "missing.gb", "--version", cwd=tmp_path)
    assert sub.returncode != 0
    assert "gbdraw version" not in sub.stdout
    assert "unrecognized arguments: --version" in sub.stderr
