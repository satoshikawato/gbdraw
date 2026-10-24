"""The Legend reserves the width its captions are drawn at (OV-165 part 2, OV-1740).

A Legend caption is drawn as SVG text: the browser advances it by the sum of
its glyph advances in CSS px. The Legend layout sizes its text column, its box
and so the canvas from that drawn advance (``drawn_text_box``), in every face,
with or without interior spaces.
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from gbdraw.circular import circular_main
from gbdraw.configurators.legend import _circular_legend_local_bounds, _linear_legend_local_bounds
from gbdraw.core.text import calculate_svg_bbox_dimensions, text_ink
from gbdraw.legend.circular_layout import build_circular_legend_layout
from gbdraw.legend.linear_layout import build_linear_legend_layout

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "legend_layout_vectors.json"
SVG = "{http://www.w3.org/2000/svg}"
FAMILIES = (
    "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif",
    "Times New Roman, serif",
    "Courier New, monospace",
)
# Captions with 0 to 5 interior spaces.
CAPTIONS = (
    "CDS",
    "GC content",
    "GC skew (-)",
    "type IV pilin protein",
    "a b c d e",
    "i l i l i l",
)
FONT_SIZE = 16.0
RECT = 20.0


def _solid(fill: str = "#54bcf8") -> dict[str, object]:
    return {"type": "solid", "fill": fill, "stroke": "none", "width": 0}


def _table(captions: tuple[str, ...]) -> dict[str, dict[str, object]]:
    return {caption: _solid() for caption in captions}


def _drawn(text: str, family: str, size: float = FONT_SIZE) -> float:
    return float(calculate_svg_bbox_dimensions(text, family, size, 96)[0])


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("side", ["right", "bottom"])
def test_circular_legend_box_contains_every_drawn_caption(family: str, side: str) -> None:
    layout = build_circular_legend_layout(
        _table(CAPTIONS),
        legend_position=side,
        canvas_width=2000.0,
        font_family=family,
        font_size=FONT_SIZE,
        color_rect_size=RECT,
    )
    bounds = _circular_legend_local_bounds(layout, color_rect_size=RECT)
    for entry in layout.solid_entries:
        right = entry.text_x + _drawn(str(entry.key), family)
        assert right <= bounds.max_x + 1e-9, (entry.key, right, bounds.max_x)


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("side", ["right", "bottom"])
def test_linear_legend_box_contains_every_drawn_caption(family: str, side: str) -> None:
    layout = build_linear_legend_layout(
        _table(CAPTIONS),
        legend_position=side,
        canvas_width=2000.0,
        font_family=family,
        font_size=FONT_SIZE,
        color_rect_size=RECT,
    )
    bounds = _linear_legend_local_bounds(layout, color_rect_size=RECT)
    active = layout.active
    for entry in active.feature.entries:
        right = active.feature_x + entry.text_x + _drawn(str(entry.key), family)
        assert right <= bounds.max_x + 1e-9, (entry.key, right, bounds.max_x)


def _corpus() -> list[tuple[str, str, float]]:
    """Every Legend caption of the shared vectors (the Gallery Legends and their edits) in each face."""
    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))
    captions = {row["key"] for case in vectors["layout"] for row in case["rows"]}
    captions |= {case["text"] for case in vectors["measurement"] if case["text"].strip()}
    sizes = {float(case["options"]["fontSize"]) for case in vectors["layout"]}
    return [(caption, family, size) for caption in sorted(captions) for family in FAMILIES for size in sorted(sizes)]


def test_legend_reserves_at_least_the_drawn_ink_of_every_corpus_caption() -> None:
    misses = []
    for caption, family, size in _corpus():
        layout = build_circular_legend_layout(
            {caption: _solid()},
            legend_position="right",
            canvas_width=2000.0,
            font_family=family,
            font_size=size,
            color_rect_size=RECT,
        )
        reserved = layout.solid_entries[0].text_width
        ink = text_ink(caption, family)
        if ink is not None and reserved < ink.right * size - 1e-9:
            misses.append((caption, family, size, reserved, ink.right * size))
    assert not misses, f"{len(misses)} captions drawn past their reserve, e.g. {misses[:5]}"


def _translate(element: ET.Element) -> tuple[float, float]:
    match = re.search(r"translate\(\s*([-\d.e+]+)\s*[, ]\s*([-\d.e+]+)?", element.get("transform", ""))
    return (float(match.group(1)), float(match.group(2) or 0)) if match else (0.0, 0.0)


def test_mono_legend_stays_in_its_box_and_the_canvas(tmp_path: Path) -> None:
    """OV-1740: `GC skew (-)` in a monospace face ended 18.1 px past the Legend box and 2.1 px past the canvas."""
    family = "Courier New, monospace"
    circular_main([
        "--gbk", str(REPO_ROOT / "examples" / "MeenMJNV.gb"), "-f", "svg", "-l", "right",
        "--font_family", family, "--no_font_embedding", "-o", str(tmp_path / "mono"),
    ])
    root = ET.parse(tmp_path / "mono.svg").getroot()
    canvas_width = float(root.get("viewBox", "").split()[2])
    composition = json.loads(root.get("data-gbdraw-composition", "{}"))
    local = composition["legend"]["localBounds"]
    legend = next(element for element in root.iter() if element.get("data-gbdraw-composition-role") == "legend")
    legend_x, _ = _translate(legend)
    captions = 0
    for row in legend.iter(SVG + "g"):
        if row.get("data-legend-key") is None:
            continue
        for label in row.findall(SVG + "text"):
            text_x, _ = _translate(label)
            right = text_x + _drawn(label.text or "", label.get("font-family", ""), float(label.get("font-size", "0")))
            assert right <= local["x"] + local["width"] + 1e-6, (label.text, right, local)
            assert legend_x + right <= canvas_width + 1e-6, (label.text, legend_x + right, canvas_width)
            captions += 1
    assert captions == 4
