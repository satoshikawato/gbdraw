"""The drawing has no background: the SVG paints nothing behind the figure, and a PNG keeps that transparency.

The Web PNG download keeps it transparent too unless White background is on
(``tests/web/presentation-options.playwright.spec.js``). PDF, EPS and PS pages have no background either (CairoSVG
paints none).
"""

from __future__ import annotations

import io
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.api import render_to_bytes
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.render.export import has_cairosvg

HMMT = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
SVG = "{http://www.w3.org/2000/svg}"


def _record():
    return next(SeqIO.parse(str(HMMT), "genbank"))


def _canvas_rects(svg_text: str) -> list[ET.Element]:
    """The rectangles that cover the whole canvas: a painted background."""

    root = ET.fromstring(svg_text)
    width, height = (float(str(root.get(name)).removesuffix("px")) for name in ("width", "height"))
    return [
        rect
        for rect in root.iter(f"{SVG}rect")
        if float(rect.get("width", "0").removesuffix("px") or 0) >= width
        and float(rect.get("height", "0").removesuffix("px") or 0) >= height
    ]


@pytest.mark.parametrize("draw", [build_circular_diagram, lambda record: build_linear_diagram([record])])
def test_the_svg_paints_no_background(draw) -> None:
    assert _canvas_rects(draw(_record()).tostring()) == []


@pytest.mark.skipif(not has_cairosvg(), reason="CairoSVG is not installed; the PNG export needs it")
def test_the_png_keeps_the_transparent_background() -> None:
    from PIL import Image

    image = Image.open(io.BytesIO(render_to_bytes(build_circular_diagram(_record()), "png")))
    assert image.mode == "RGBA"
    assert image.getpixel((0, 0))[3] == 0
