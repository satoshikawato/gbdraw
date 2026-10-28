"""A requested PNG, PDF, EPS or PS that CairoSVG cannot make fails the run.

Without CairoSVG the command line stops before it draws, prints the formats
it cannot make and the install command on stderr, exits non-zero, and writes
no file, as the Python API does. The loader is patched, so the tests do not
depend on CairoSVG being absent: CI installs the ``export`` extra.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from svgwrite import Drawing

from gbdraw import cli
from gbdraw.api import render as api_render
from gbdraw.exceptions import ValidationError
from gbdraw.render import export as export_module

GENBANK = Path(__file__).resolve().parent / "test_inputs" / "HmmtDNA.gbk"
INSTALL_HINT = "pip install gbdraw[export]"


@pytest.fixture
def without_cairosvg(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(export_module, "_load_cairosvg", lambda: None)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    ("formats", "missing"),
    [("png", "PNG"), ("svg,png", "PNG"), ("svg,pdf,eps", "PDF, EPS")],
)
def test_cli_exits_non_zero_when_cairosvg_cannot_make_a_requested_format(
    without_cairosvg: None,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    mode: str,
    formats: str,
    missing: str,
) -> None:
    prefix = tmp_path / "figure"
    monkeypatch.setattr(
        sys, "argv", ["gbdraw", mode, "--gbk", str(GENBANK), "-o", str(prefix), "-f", formats]
    )

    with pytest.raises(SystemExit) as raised:
        cli.main()

    assert raised.value.code == 1
    stderr = capsys.readouterr().err
    assert f"ERROR: Cannot generate {missing}: CairoSVG is not installed." in stderr
    assert INSTALL_HINT in stderr
    # The run stops before it draws: no SVG, as the Python API writes none.
    assert sorted(tmp_path.iterdir()) == []


def test_cli_and_python_api_report_missing_cairosvg_with_one_message(
    without_cairosvg: None,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["gbdraw", "circular", "--gbk", str(GENBANK), "-o", str(tmp_path / "cli"), "-f", "svg,png"],
    )
    with pytest.raises(SystemExit):
        cli.main()
    cli_message = capsys.readouterr().err.strip().removeprefix("ERROR: ")

    drawing = Drawing(filename=str(tmp_path / "api.svg"))
    with pytest.raises(ValidationError) as saved:
        api_render.save_figure_to(drawing, "svg,png", output_dir=str(tmp_path), output_prefix="api")
    with pytest.raises(ValidationError) as rendered:
        api_render.render_to_bytes(drawing, "png")

    assert str(saved.value) == cli_message
    assert str(rendered.value) == cli_message
    assert sorted(tmp_path.iterdir()) == []
