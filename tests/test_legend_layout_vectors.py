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
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pytest

from gbdraw.core.text import drawn_text_box
from gbdraw.core.numeric import compensated_sum
from gbdraw.legend.metrics import legend_caption_box
from gbdraw.session_container import SESSION_DOCUMENT_ENTRY
from gbdraw.session_drawings import checked_drawing_views

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "legend_layout_vectors.json"
GALLERY_SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
SVG = "{http://www.w3.org/2000/svg}"


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


def test_legend_layout_vectors_follow_the_gallery_and_python() -> None:
    # The Gallery Legends are read again from the Sessions, so a Gallery refresh
    # that changes a Legend fails here, not late in a Playwright shard.
    result = _run("tools/generate_legend_layout_vectors.py", "--check", "--from-gallery")
    assert result.returncode == 0, result.stdout + result.stderr


def _translation(element: ET.Element) -> tuple[float, float]:
    pairs = re.findall(r"translate\(\s*([^,\s)]+)\s*[, ]\s*([^)\s]+)?\s*\)", element.get("transform", ""))
    return (sum(float(x) for x, _ in pairs), sum(float(y or 0) for _, y in pairs))


def _drawn_gallery_legends() -> dict[str, dict[str, object]]:
    """Each Gallery Result's Legend as Python drew it: its local box and the offset of its gradient group."""
    drawn: dict[str, dict[str, object]] = {}
    for path in sorted(GALLERY_SESSIONS.glob("*.gbdraw")):
        with zipfile.ZipFile(path) as archive:
            session = json.loads(archive.read(SESSION_DOCUMENT_ENTRY))
        results = [result for view in checked_drawing_views(session) for result in view.get("results") or []]
        for index, result in enumerate(results):
            root = ET.fromstring(result["content"].encode())
            ids = {element.get("id"): element for element in root.iter() if element.get("id")}
            if "legend" not in ids:
                continue
            box = json.loads(root.get("data-gbdraw-composition"))["legend"]["localBounds"]
            gradient = ids.get("pairwise_legend_h") if "legend_horizontal" in ids else ids.get("conservation_identity_legend")
            drawn[f"{path.name.removesuffix('.gbdraw')}#{index}"] = {
                "localBounds": [box["x"], box["y"], box["x"] + box["width"], box["y"] + box["height"]],
                "gradientOffset": None if gradient is None else _translation(gradient),
            }
    return drawn


def test_gallery_legend_vectors_are_the_drawn_gallery_legends() -> None:
    # Every Gallery Legend has a vector, and Python's layout of the rows read
    # from it is the Legend Python drew (OV-3840: similarity-ring rows read as
    # gradient rows laid the E. coli rings Legend out 11.4 px wider).
    drawn = _drawn_gallery_legends()
    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))
    stored = {
        case["name"]: case["expected"] for case in vectors["layout"]
        if case["edit"] == "none" and not case["name"].startswith("synthetic")
    }
    assert sorted(stored) == sorted(drawn)
    for name, legend in drawn.items():
        expected = stored[name]
        assert expected["localBounds"] == pytest.approx(legend["localBounds"], rel=0, abs=1e-9), name
        if legend["gradientOffset"] is not None:
            linear = "activeOrientation" in expected
            source = expected["horizontal"] if linear else expected
            assert (source["gradientX"], source["gradientY"]) == legend["gradientOffset"], name


def test_legend_layout_vectors_cover_gallery_edits_and_branches() -> None:
    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))
    gallery = [case for case in vectors["layout"] if not case["name"].startswith("synthetic")]
    assert sum(case["edit"] != "none" for case in gallery) >= 80
    assert {case["mode"] for case in vectors["layout"]} == {"circular", "linear"}
    assert len(vectors["measurement"]) >= 200


@pytest.mark.parametrize("index", [0, 40, 100, 200])
def test_measurement_vectors_are_legend_caption_box(index: int) -> None:
    vector = json.loads(VECTORS.read_text(encoding="utf-8"))["measurement"][index]
    width, height = legend_caption_box(vector["text"], vector["fontFamily"], vector["fontSize"])
    assert vector["expected"] == {"width": width, "height": height}


def test_literal_measurement_vectors_are_drawn_text_box() -> None:
    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))["measurement"]
    literal = [vector for vector in vectors if "<i>" not in vector["text"]]
    assert len(literal) >= 456
    for vector in literal[::37]:
        width, height = drawn_text_box(vector["text"], vector["fontFamily"], vector["fontSize"])
        assert vector["expected"] == {"width": width, "height": height}, vector["text"]


def _naive_sum(values: list[float]) -> float:
    total = 0.0
    for value in values:
        total += value
    return total


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([0.1] * 10, 1.0),
        ([1e100, 1.0, -1e100], 1.0),
        ([131.05357142857142, 131.05357142857142, 108.93452380952381], 371.04166666666663),
    ],
)
def test_compensated_sum_is_cpython_312_sum_on_every_python(values: list[float], expected: float) -> None:
    # The Circular Legend sums row widths; Python 3.10 and 3.11 add naively,
    # 3.12+ compensate, and the Web port compensates. The helper gives the
    # 3.12+ doubles on every version, so the layout does not depend on it.
    assert compensated_sum(values) == expected
    if sys.version_info >= (3, 12):
        assert sum(values) == expected


def test_compensated_sum_differs_from_naive_addition() -> None:
    values = [0.1] * 10
    assert _naive_sum(values) == 0.9999999999999999
    assert compensated_sum(values) == 1.0


def test_compensated_sum_edges() -> None:
    assert compensated_sum([]) == 0.0
    assert str(compensated_sum([-0.0])) == "0.0"
    assert compensated_sum([float("inf"), 1.0]) == float("inf")
    assert compensated_sum(iter([2.5, 1])) == 3.5
