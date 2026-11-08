"""Where the gbdraw command writes its log: every line to standard error, one level label each."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

from gbdraw.config.toml import load_config_toml

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ROOT / "tests" / "fixtures" / "input_variants"
QUERY = VARIANTS / "prokka_locus_tag_only.gb"
SUBJECT = VARIANTS / "rast_seed_dbxref.gb"
ROWS_BELOW_IDENTITY = VARIANTS / "blast_rows_below_identity.tsv"

_LABEL = re.compile(r"^(INFO|WARNING|ERROR): ")
_TWO_LABELS = re.compile(r"^(INFO|WARNING|ERROR): +(INFO|WARNING|ERROR):")


def _python(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=cwd, capture_output=True, text=True, env={"PYTHONPATH": str(ROOT), "PATH": ""}, check=False,
    )


def _gbdraw(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return _python("-m", "gbdraw.cli", *args, cwd=cwd)


def _assert_log_on_stderr_only(run: subprocess.CompletedProcess[str]) -> list[str]:
    assert run.returncode == 0, run.stderr
    assert run.stdout == ""
    lines = run.stderr.splitlines()
    assert lines
    assert all(_LABEL.match(line) for line in lines), lines
    assert not any(_TWO_LABELS.match(line) for line in lines), lines
    return lines


def test_a_circular_run_writes_its_log_to_stderr_and_nothing_to_stdout(tmp_path: Path) -> None:
    run = _gbdraw("circular", "--gbk", str(QUERY), "-o", "out", "-f", "svg", cwd=tmp_path)

    lines = _assert_log_on_stderr_only(run)
    assert "INFO: Generated SVG: out.svg" in lines


def test_a_linear_run_writes_its_warning_to_stderr_and_reads_the_config_once(tmp_path: Path) -> None:
    run = _gbdraw(
        "linear", "--gbk", str(QUERY), str(SUBJECT), "-b", str(ROWS_BELOW_IDENTITY), "--identity", "70",
        "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    lines = _assert_log_on_stderr_only(run)
    assert any(
        line.startswith("WARNING: Comparison table blast_rows_below_identity.tsv") and "leave out all 5 rows" in line
        for line in lines
    ), lines
    assert [line for line in lines if line.startswith("INFO: Loading config file: ")] == [
        f"INFO: Loading config file: {(ROOT / 'gbdraw' / 'data' / 'config.toml').resolve()}"
    ]
    assert "INFO: Generated SVG: out.svg" in lines


@pytest.mark.parametrize(
    ("module", "inputs"),
    [("gbdraw.circular", [QUERY]), ("gbdraw.linear", [QUERY, SUBJECT])],
)
def test_running_a_mode_module_directly_logs_like_the_command(tmp_path: Path, module: str, inputs: list[Path]) -> None:
    run = _python("-m", module, "--gbk", *(str(path) for path in inputs), "-o", "out", "-f", "svg", cwd=tmp_path)

    lines = _assert_log_on_stderr_only(run)
    assert "INFO: Generated SVG: out.svg" in lines


def test_importing_the_mode_modules_leaves_logging_alone(tmp_path: Path) -> None:
    run = _python(
        "-c",
        "import logging; import gbdraw.circular, gbdraw.linear, gbdraw.cli; "
        "root = logging.getLogger(); print(len(root.handlers), logging.getLevelName(root.level))",
        cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    assert run.stdout.split() == ["0", "WARNING"]


def test_each_config_load_returns_its_own_dict() -> None:
    first = load_config_toml("gbdraw.data", "config.toml")
    first["canvas"]["linear"]["track_layout"] = "changed"
    first.setdefault("labels", {})["added"] = True

    second = load_config_toml("gbdraw.data", "config.toml")

    assert second["canvas"]["linear"]["track_layout"] != "changed"
    assert "added" not in second["labels"]


def test_gui_status_lines_go_to_stderr(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    import socketserver
    import webbrowser

    from gbdraw import cli

    def _interrupt(self: socketserver.BaseServer, *args: object, **kwargs: object) -> None:
        raise KeyboardInterrupt

    opened: list[str] = []
    monkeypatch.setattr(webbrowser, "open", opened.append)
    monkeypatch.setattr(socketserver.TCPServer, "serve_forever", _interrupt)
    monkeypatch.setattr(socketserver.TCPServer, "shutdown", lambda self: None)

    cli.start_local_server(str(tmp_path))

    captured = capsys.readouterr()
    assert captured.out == ""
    assert f"Local server started at: {opened[0]}" in captured.err
    assert "Stopping server..." in captured.err


def test_an_argument_error_carries_one_label(tmp_path: Path) -> None:
    run = _gbdraw("circular", "--gff", str(tmp_path / "genome.gff"), cwd=tmp_path)

    assert run.returncode == 2
    assert run.stdout == ""
    assert "error: --gff requires --fasta." in run.stderr
    assert "Error:" not in run.stderr
