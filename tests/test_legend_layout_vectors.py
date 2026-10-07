"""Python's side of the shared Legend layout vectors.

``tests/fixtures/legend_layout_vectors.json`` holds Python's Legend
measurement, layout, local bounds and composition for shared inputs; the Web
port (``gbdraw/web/js/services/legend-layout.js``) must equal them exactly
(``tests/web/legend-layout-vectors.test.mjs``). These tests keep the stored
values equal to what Python computes now, and the browser font table equal to
the bundled fonts, so a change on the Python side fails until the vectors and
the port follow it.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from gbdraw.core.text import calculate_bbox_dimensions

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "legend_layout_vectors.json"


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_generated_legend_font_metrics_match_bundled_fonts() -> None:
    result = _run("tools/generate_legend_font_metrics.py", "--check")
    assert result.returncode == 0, result.stdout + result.stderr


def test_legend_layout_vectors_match_python() -> None:
    result = _run("tools/generate_legend_layout_vectors.py", "--check")
    assert result.returncode == 0, result.stdout + result.stderr


def test_legend_layout_vectors_cover_gallery_edits_and_branches() -> None:
    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))
    gallery = [case for case in vectors["layout"] if not case["name"].startswith("synthetic")]
    assert sum(case["edit"] == "none" for case in gallery) == 10
    assert sum(case["edit"] != "none" for case in gallery) >= 80
    assert {case["mode"] for case in vectors["layout"]} == {"circular", "linear"}
    assert len(vectors["measurement"]) >= 200


@pytest.mark.parametrize("index", [0, 40, 100, 200])
def test_measurement_vectors_are_calculate_bbox_dimensions(index: int) -> None:
    vector = json.loads(VECTORS.read_text(encoding="utf-8"))["measurement"][index]
    width, height = calculate_bbox_dimensions(
        vector["text"], vector["fontFamily"], vector["fontSize"], vector["dpi"]
    )
    assert vector["expected"] == {"width": width, "height": height}
