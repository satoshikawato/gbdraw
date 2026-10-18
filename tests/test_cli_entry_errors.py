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


def test_cli_existing_output_error_names_the_overwrite_option(tmp_path: Path) -> None:
    """OV-1240: the CLI said "Use overwrite=True", the Python spelling, instead of ``--overwrite``."""

    gbk = str(TEST_INPUTS / "NC_001416.gb")
    first = _gbdraw("circular", "--gbk", gbk, "-o", "out", "--save_session", cwd=tmp_path, timeout=300)
    assert first.returncode == 0, first.stderr
    again = _gbdraw("circular", "--gbk", gbk, "-o", "out", cwd=tmp_path, timeout=300)
    assert again.returncode == 1
    assert "ERROR: Output file already exists: " in again.stderr
    assert "Use --overwrite to replace it." in again.stderr
    assert "overwrite=True" not in again.stderr
    sessions = sorted(tmp_path.glob("*.gbdraw-session.json*"))
    assert sessions, sorted(path.name for path in tmp_path.iterdir())
    replay = _gbdraw("render", "--session", str(sessions[0]), "-o", "replay",
                     "--session_output", sessions[0].name, cwd=tmp_path, timeout=300)
    assert replay.returncode == 1
    assert "ERROR: Session output already exists: " in replay.stderr
    assert "Use --overwrite to replace it." in replay.stderr
    assert "overwrite=True" not in replay.stderr


def test_python_existing_output_errors_share_one_wording(tmp_path: Path) -> None:
    """OV-1240: one owner words the refusal; Python callers read ``overwrite=True``."""

    from gbdraw.render.output_paths import OutputExistsError, preflight_output_paths

    first = tmp_path / "a.svg"
    second = tmp_path / "b.svg"
    first.write_text("x")
    second.write_text("x")
    with pytest.raises(OutputExistsError) as one:
        preflight_output_paths((first,), overwrite=False)
    assert str(one.value) == f"Output file already exists: {first}. Use overwrite=True to replace it."
    assert one.value.message("--overwrite") == f"Output file already exists: {first}. Use --overwrite to replace it."
    with pytest.raises(OutputExistsError) as two:
        preflight_output_paths((first, second), overwrite=False)
    assert str(two.value) == f"Output files already exist: {first}, {second}. Use overwrite=True to replace them."
