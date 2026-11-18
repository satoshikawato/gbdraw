"""A Circular label is drawn inside its feature only when its glyphs stay inside the feature's drawn outline
(OV-3649): the embed test read the label's length against the feature's span alone, so a label near the end of its
feature ran into the arrowhead or touched the narrowed outline beside it.

The checks read the drawn SVG alone: each embedded label's glyph box from its text path and the font's ink
(``text_ink``), and its feature's outline from the feature path, less half the path's stroke.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from gbdraw.core.text import text_ink
from tests.utils.cli_in_process import run_gbdraw

INPUTS = Path(__file__).parent / "test_inputs"
SESSIONS = Path(__file__).parent.parent / "gbdraw" / "web" / "gallery" / "sessions"
QUALIFIER_PRIORITY = Path(__file__).parent.parent / "examples" / "qualifier_priority.tsv"
SVG = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}"
_NUMBER = r"-?\d+(?:\.\d+)?(?:e-?\d+)?"
Point = tuple[float, float]

# Feature ids (``data-gbdraw-feature-id``) of the human mitochondrial genome (HmmtDNA.gbk, NC_012920.1).
HMMTDNA_ATP8 = "f31dbd2ba"
HMMTDNA_ND4L = "f63cda8a2"
HMMTDNA_ND6 = "f6b2a293e"


def _draw(tmp_path: Path, gbk: str, *options: str) -> ET.Element:
    return _draw_args(tmp_path, "--gbk", str(INPUTS / gbk), *options)


def _draw_args(tmp_path: Path, *args: str) -> ET.Element:
    result = run_gbdraw(["circular", *args, "-o", "out", "-f", "svg"], cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return ET.parse(tmp_path / "out.svg").getroot()


def _gene_labels(tmp_path: Path) -> Path:
    path = tmp_path / "cds_gene_qualifier_priority.tsv"
    path.write_text("CDS\tgene\n", encoding="utf-8")
    return path


def _outline(d: str) -> list[Point]:
    """A feature path as a polygon: its moves and lines, and each arc sampled along its circle."""

    points: list[Point] = []
    tokens = re.findall(rf"[MLAZmlaz]|{_NUMBER}", d)
    command = ""
    index = 0
    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
            continue
        if command == "A":
            sweep = int(float(tokens[index + 4]))
            x, y = float(tokens[index + 5]), float(tokens[index + 6])
            x0, y0 = points[-1]
            a0 = math.atan2(y0, x0)
            turn = (math.atan2(y, x) - a0) % math.tau if sweep else -((a0 - math.atan2(y, x)) % math.tau)
            r0, r1 = math.hypot(x0, y0), math.hypot(x, y)
            for k in range(1, 64):
                radius = r0 + (r1 - r0) * k / 64
                points.append((radius * math.cos(a0 + turn * k / 64), radius * math.sin(a0 + turn * k / 64)))
            points.append((x, y))
            index += 7
        else:
            points.append((float(tokens[index]), float(tokens[index + 1])))
            index += 2
    return points


def _signed_gap(point: Point, polygon: list[Point]) -> float:
    """The distance from ``point`` to the polygon's outline: positive inside it, negative outside (even-odd rule)."""

    x, y = point
    edges = list(zip(polygon, polygon[1:] + polygon[:1]))
    crossings = sum(1 for (x1, y1), (x2, y2) in edges if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1))
    gap = math.inf
    for (x1, y1), (x2, y2) in edges:
        dx, dy = x2 - x1, y2 - y1
        t = max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / ((dx * dx + dy * dy) or 1.0)))
        gap = min(gap, math.hypot(x - x1 - t * dx, y - y1 - t * dy))
    return gap if crossings % 2 else -gap


def _glyph_corners(text_path: ET.Element, arc: str) -> list[Point]:
    """The corners of the glyph box an embedded label paints along its text path (``startOffset="50%"``,
    ``text-anchor="middle"``): a clockwise path turns the glyphs' tops outward, an anticlockwise one inward."""

    # The arc may be drawn in several pieces, all on one circle and turning one way.
    values = [float(value) for value in re.findall(_NUMBER, arc)]
    x0, y0, radius, sweep, x1, y1 = values[0], values[1], values[2], int(values[6]), values[-2], values[-1]
    a0, a1 = math.atan2(y0, x0), math.atan2(y1, x1)
    direction = 1 if sweep else -1
    span = (a1 - a0) % math.tau if sweep else (a0 - a1) % math.tau
    size = float(text_path.get("font-size", "nan"))
    ink = text_ink(text_path.text or "", text_path.get("font-family", ""), "normal", text_path.get("font-style", "normal"))
    assert ink is not None
    middle = a0 + direction * span / 2
    band = (radius + ink.bottom * size, radius + ink.top * size) if direction == 1 else (
        radius - ink.top * size, radius - ink.bottom * size
    )
    corners = []
    for along in ((ink.left - ink.advance / 2) * size, (ink.right - ink.advance / 2) * size):
        angle = middle + direction * along / radius
        corners += [(r * math.cos(angle), r * math.sin(angle)) for r in band]
    return corners


def _embedded_glyph_gaps(root: ET.Element) -> dict[str, float]:
    """For each embedded label, by its feature id: how far its glyph box stays inside its feature's outline, past half
    the outline's stroke (negative: the glyphs touch or cross the stroke)."""

    arcs = {path.get("id"): path.get("d", "") for path in root.iter(f"{SVG}path") if path.get("id")}
    outlines: dict[str, list[tuple[list[Point], float]]] = {}
    for path in root.iter(f"{SVG}path"):
        feature_id = path.get("data-gbdraw-feature-id")
        if feature_id and path.get("fill", "none") != "none":
            outlines.setdefault(feature_id, []).append((_outline(path.get("d", "")), float(path.get("stroke-width") or 0)))
    gaps = {}
    for text in root.iter(f"{SVG}text"):
        text_path = text.find(f"{SVG}textPath")
        feature_id = text.get("data-label-feature-id")
        if text_path is None or not feature_id:
            continue
        corners = _glyph_corners(text_path, arcs[(text_path.get(f"{XLINK}href") or text_path.get("href", ""))[1:]])
        gaps[feature_id] = max(
            min(_signed_gap(corner, polygon) for corner in corners) - 0.5 * stroke
            for polygon, stroke in outlines[feature_id]
        )
    return gaps


@pytest.mark.parametrize(
    ("gbk", "options"),
    [
        pytest.param(
            # ATP8 runs into its arrowhead once its length is read at the lane centre without the outline (H-CLI-11).
            "HmmtDNA.gbk",
            ("-k", "CDS,rRNA,tRNA", "--feature_shape", "CDS=arrow", "--arrow_shaft_width_ratio", "0.72",
             "--resolve_overlaps", "--separate_strands", "--track_type", "spreadout", "--labels", "both",
             "--block_stroke_width", "1.5"),
            id="human-mitochondrion-narrow-shafts",
        ),
        pytest.param(
            # ndhH touched its arrowhead's stroke; petA crosses its head when read at the lane centre alone.
            "NC_001879.gbk",
            ("-k", "CDS,rRNA,tRNA", "--qualifier_priority", str(QUALIFIER_PRIORITY), "--separate_strands",
             "--track_type", "tuckin", "--labels", "both", "--label_placement", "radial", "--block_stroke_width", "1"),
            id="tobacco-plastome",
        ),
    ],
)
def test_embedded_labels_keep_their_glyphs_inside_the_drawn_outline(
    tmp_path: Path, gbk: str, options: tuple[str, ...]
) -> None:
    """OV-3649: every embedded label's glyphs clear its feature's stroke, beside arrowheads too."""

    extra = ("--qualifier_priority", str(_gene_labels(tmp_path))) if gbk == "HmmtDNA.gbk" else ()
    gaps = _embedded_glyph_gaps(_draw(tmp_path, gbk, *extra, *options))

    assert len(gaps) >= 10
    assert {feature_id: round(gap, 2) for feature_id, gap in gaps.items() if gap < 0.0} == {}


def test_the_tobacco_gallery_figure_keeps_embedded_glyphs_inside_the_drawn_outline(tmp_path: Path) -> None:
    """OV-3649: in the Gallery tobacco plastome ndhD touched its arrowhead's stroke, and psbC ran into its head once
    its length was read at the lane centre without the outline."""

    gaps = _embedded_glyph_gaps(_draw_args(tmp_path, "--session", str(SESSIONS / "tobacco-chloroplast.gbdraw")))

    assert len(gaps) >= 10
    assert {feature_id: round(gap, 2) for feature_id, gap in gaps.items() if gap < 0.0} == {}


def test_a_label_whose_glyphs_clear_its_outline_is_drawn_inside_an_outer_lane(tmp_path: Path) -> None:
    """OV-3649: the embed length was read at the axis radius as well as the lane centre, so on a lane outside the
    axis ND4L went outside although its glyphs clear its outline; it is read at the lane centre alone now."""

    gaps = _embedded_glyph_gaps(
        _draw(
            tmp_path, "HmmtDNA.gbk", "--qualifier_priority", str(_gene_labels(tmp_path)), "--separate_strands",
            "--track_type", "middle", "--labels", "out",
        )
    )

    assert HMMTDNA_ND4L in gaps
    assert gaps[HMMTDNA_ND4L] > 0.0


def test_a_label_whose_glyphs_clear_its_outline_is_drawn_inside_a_lane_inside_the_axis(tmp_path: Path) -> None:
    """OV-3649: on the minus-strand lane of ``--separate_strands``, inside the axis, "ND6 prot" fits inside ND6's
    drawn outline; the 1.1x span test at the lane centre put it outside although its glyphs clear the outline."""

    label_table = tmp_path / "labels.tsv"
    label_table.write_text(
        "record_id\tfeature_type\tqualifier\tvalue\tlabel_text\n*\t*\tlabel\t^NADH dehydrogenase subunit 6$\tND6 prot\n",
        encoding="utf-8",
    )
    root = _draw(
        tmp_path, "HmmtDNA.gbk", "--label_table", str(label_table), "--separate_strands", "--track_type", "tuckin",
        "--labels", "out",
    )
    labels = {
        text.get("data-label-feature-id"): "".join(text.itertext())
        for text in root.iter(f"{SVG}text")
        if text.get("data-label-feature-id")
    }
    gaps = _embedded_glyph_gaps(root)

    assert labels[HMMTDNA_ND6] == "ND6 prot"
    assert HMMTDNA_ND6 in gaps
    assert gaps[HMMTDNA_ND6] > 0.0
    assert min(gaps.values()) > 0.0
