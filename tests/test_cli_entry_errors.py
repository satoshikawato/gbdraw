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


def test_cli_typed_option_error_is_not_its_own_cause(tmp_path: Path) -> None:
    """A rejected Circular option was re-raised from itself, so walking its causes never ended."""

    gbk = str(TEST_INPUTS / "NC_001416.gb")
    sub = _gbdraw("circular", "--gbk", gbk, "-o", "out", "--window", "0", cwd=tmp_path, timeout=120)
    assert sub.returncode == 1
    assert "window must be a positive integer or None." in sub.stderr
    assert "Traceback" not in sub.stderr


def test_parser_reason_is_in_the_raised_error(tmp_path: Path) -> None:
    """OV-1241: the parser's reason went to the log only; API and Web callers lost it."""

    import gbdraw
    from gbdraw.exceptions import ParseError

    cut = tmp_path / "cut.gb"
    cut.write_bytes((TEST_INPUTS / "NC_001416.gb").read_bytes()[:3000])
    with pytest.raises(ParseError, match="Premature end of line during sequence data"):
        gbdraw.read_genbank(cut)
    gff = tmp_path / "bad.gff"
    gff.write_text("##gff-version 3\nseq1\tsrc\tCDS\tabc\t10\t.\t+\t0\tID=a\n", encoding="utf-8")
    fasta = tmp_path / "a.fa"
    fasta.write_text(">seq1\nACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(ParseError, match="invalid literal for int"):
        gbdraw.read_gff(gff, fasta)
    sub = _gbdraw("circular", "--gbk", str(cut), "-o", "out", cwd=tmp_path)
    assert sub.returncode == 1
    assert sub.stderr.count("Premature end of line during sequence data") == 1, sub.stderr


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


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("missing_gbk", "File does not exist or is not accessible: nonexist.gb"),
        ("default_colors_extra_column", "Malformed line in default colors file 'd3.tsv' at line 1: expected 2 columns."),
        ("missing_color_table", "Color table file not found:"),
        ("color_table_six_columns", "Malformed line in color table 't6.tsv' at line 1: expected 5 columns."),
        ("whitelist_blank_keyword", "Missing values in 'wl.tsv' at line 1. Missing columns: ['keyword']."),
    ],
)
def test_cli_input_error_prints_once(case: str, message: str, tmp_path: Path) -> None:
    """OV-1241: the reader logged the error and the CLI printed it again."""

    (tmp_path / "d3.tsv").write_text("CDS\t#ff0000\textra\n", encoding="utf-8")
    (tmp_path / "t6.tsv").write_text("CDS\tproduct\tx\t#ff0000\tcap\textra\n", encoding="utf-8")
    (tmp_path / "wl.tsv").write_text("CDS\tproduct\t\n", encoding="utf-8")
    gbk = str(TEST_INPUTS / "NC_001416.gb")
    argv = {
        "missing_gbk": ["--gbk", "nonexist.gb"],
        "default_colors_extra_column": ["--gbk", gbk, "-d", "d3.tsv"],
        "missing_color_table": ["--gbk", gbk, "-t", "nope.tsv"],
        "color_table_six_columns": ["--gbk", gbk, "-t", "t6.tsv"],
        "whitelist_blank_keyword": ["--gbk", gbk, "--label_whitelist", "wl.tsv"],
    }[case]
    result = _gbdraw("circular", *argv, "-o", "out", cwd=tmp_path)
    assert result.returncode == 1
    output = result.stdout + result.stderr
    assert output.count(message) == 1, output
    assert output.count("ERROR:") == 1, output


def test_no_reader_logs_an_error_it_raises() -> None:
    """OV-1241: an error a gbdraw function raises is reported once, by its caller.

    A ``logger.error`` (or ``critical``/``exception``) directly before a
    ``raise`` printed the same error twice on the CLI.
    """

    import ast

    import gbdraw

    offenders = []
    package = Path(gbdraw.__file__).parent
    for path in sorted(package.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            for field in ("body", "orelse", "finalbody"):
                body = getattr(node, field, None)
                if not isinstance(body, list):
                    continue
                for statement, following in zip(body, body[1:]):
                    call = statement.value if isinstance(statement, ast.Expr) else None
                    if (
                        isinstance(call, ast.Call)
                        and isinstance(call.func, ast.Attribute)
                        and call.func.attr in {"error", "critical", "exception"}
                        # a logger (``logger``, ``logging``), not ``parser.error``
                        and isinstance(call.func.value, ast.Name)
                        and "log" in call.func.value.id.lower()
                        and isinstance(following, ast.Raise)
                    ):
                        offenders.append(f"{path.relative_to(package)}:{statement.lineno}")
    assert offenders == []


GALLERY_SESSION = Path(__file__).parents[1] / "gbdraw" / "web" / "gallery" / "sessions" / "HmmtDNA_basic_circular.gbdraw-session.json"


@pytest.mark.parametrize("damage", ["truncated_gzip", "corrupt_gzip", "invalid_utf8"])
def test_cli_reports_an_unreadable_session_file_without_a_traceback(damage: str, tmp_path: Path) -> None:
    """OV-1251: a truncated or corrupt gzip Session, or text that is not UTF-8, ended in a traceback."""

    import gzip

    text = GALLERY_SESSION.read_bytes()
    compressed = gzip.compress(text)
    target = tmp_path / ("bad.json" if damage == "invalid_utf8" else "bad.json.gz")
    target.write_bytes({
        "truncated_gzip": compressed[: len(compressed) // 2],
        "corrupt_gzip": compressed[:20] + bytes(64) + compressed[84:],
        "invalid_utf8": b'{"format": "\xff\xfe"}',
    }[damage])
    result = _gbdraw("render", "--session", str(target), "-o", "out", cwd=tmp_path)
    assert result.returncode == 1
    assert "Traceback" not in result.stderr
    assert f"ERROR: Not a valid JSON session file: {target}" in result.stderr


def test_python_session_loaders_share_the_unreadable_file_error(tmp_path: Path) -> None:
    from gbdraw.exceptions import ValidationError
    from gbdraw.session import SessionFormatError, load_session_document
    from gbdraw.session_io import load_session

    target = tmp_path / "bad.json.gz"
    target.write_bytes(b"\x1f\x8b" + bytes(10))
    with pytest.raises(ValidationError, match="Not a valid JSON session file"):
        load_session(target)
    with pytest.raises(SessionFormatError, match="Not a valid JSON session file"):
        load_session_document(target)


def test_session_shape_errors_name_the_file_version(tmp_path: Path) -> None:
    """OV-1252: a version-48 file was told "Session version 40 requires a results array"."""

    import json

    from gbdraw.session_io import CURRENT_SESSION_VERSION

    payload = json.loads(GALLERY_SESSION.read_text(encoding="utf-8"))
    assert payload["version"] == CURRENT_SESSION_VERSION
    payload["results"] = "x"
    target = tmp_path / "badresults.json"
    target.write_text(json.dumps(payload), encoding="utf-8")
    result = _gbdraw("render", "--session", str(target), "-o", "out", cwd=tmp_path)
    assert result.returncode == 1
    assert f"Session version {CURRENT_SESSION_VERSION} requires a results array." in result.stderr
    assert "version 40" not in result.stderr
