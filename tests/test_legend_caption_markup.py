"""Legend captions written with ``<i>`` draw those words in italics (OV-2760).

A Legend caption (a color rule's Legend name, a comparison ring's label, a
gradient title) that uses the ``<i>...</i>`` markup Species, record names and
plot titles accept is drawn as one ``<text>`` with an italic ``<tspan>`` per
italic run, and the Legend reserves the width of the drawn runs. Any other
caption is drawn and measured as literal text, as before.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from svgwrite import Drawing

from gbdraw.circular import circular_main
from gbdraw.core.text import calculate_svg_bbox_dimensions, drawn_text_box
from gbdraw.legend.metrics import legend_caption_box, legend_caption_runs
from gbdraw.legend.circular_layout import build_circular_legend_layout
from gbdraw.legend.linear_layout import build_linear_legend_layout
from gbdraw.linear import linear_main
from gbdraw.render.groups.circular.legend import LegendGroup as CircularLegendGroup
from gbdraw.render.groups.linear.legend import LegendGroup as LinearLegendGroup

from tests.test_annotation_gradient_legend import (
    PAIRWISE_ROW,
    RING_ROW,
    _circular_measurement,
    _linear_measurement,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SVG = "{http://www.w3.org/2000/svg}"
FAMILY = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
SERIF = "Times New Roman, serif"
SIZE = 16.0
DANIO = "<i>Danio rerio</i> (zebrafish)"


@pytest.mark.parametrize(
    ("caption", "runs"),
    [
        ("CDS", (("CDS", "normal"),)),
        ("", (("", "normal"),)),
        (DANIO, (("Danio rerio", "italic"), (" (zebrafish)", "normal"))),
        ("<i>E. coli</i> and <i>V. cholerae</i>", (("E. coli", "italic"), (" and ", "normal"), ("V. cholerae", "italic"))),
        ("<i>Homo sapiens</i>", (("Homo sapiens", "italic"),)),
        # Outside the subset: literal text, as before.
        ("<b>x</b>", (("<b>x</b>", "normal"),)),
        ("a & <i>b</i>", (("a & <i>b</i>", "normal"),)),
        ("<i>a &amp; b</i>", (("<i>a &amp; b</i>", "normal"),)),
        ("<i></i>x", (("<i></i>x", "normal"),)),
        ("<i><i>a</i></i>", (("<i><i>a</i></i>", "normal"),)),
        ("<I>a</I>", (("<I>a</I>", "normal"),)),
        ("<i>a", (("<i>a", "normal"),)),
        ("a</i>", (("a</i>", "normal"),)),
        ("1 < 2", (("1 < 2", "normal"),)),
    ],
)
def test_caption_runs_follow_the_italic_subset(caption: str, runs: tuple[tuple[str, str], ...]) -> None:
    assert legend_caption_runs(caption) == runs


@pytest.mark.parametrize("caption", ["CDS", "GC skew (-)", "<b>x</b>", "a & <i>b</i>", "  spaced  ", "", "日本語"])
def test_a_literal_caption_measures_as_drawn_text_box(caption: str) -> None:
    assert legend_caption_box(caption, FAMILY, SIZE) == drawn_text_box(caption, FAMILY, SIZE)


def _run_sum(caption: str, family: str = FAMILY) -> float:
    """The drawn width of the caption's runs, each in its own face, each starting at the previous run's advance."""

    right = origin = 0.0
    for text, style in legend_caption_runs(caption):
        right = max(right, origin + drawn_text_box(text, family, SIZE, "normal", style)[0])
        origin += float(calculate_svg_bbox_dimensions(text, family, SIZE, "normal", style)[0])
    return max(right, origin)


@pytest.mark.parametrize("family", [FAMILY, SERIF])
def test_a_markup_caption_measures_its_runs_in_their_faces(family: str) -> None:
    width, height = legend_caption_box(DANIO, family, SIZE)
    assert width == _run_sum(DANIO, family)
    assert height == drawn_text_box("Danio rerio", family, SIZE, "normal", "italic")[1]


def test_italic_runs_are_measured_in_the_italic_face() -> None:
    # Liberation Serif Italic advances differ from its Regular ones (Liberation Sans shares them).
    assert legend_caption_box(DANIO, SERIF, SIZE)[0] != drawn_text_box("Danio rerio (zebrafish)", SERIF, SIZE)[0]


def _svg(group: object) -> ET.Element:
    drawing = Drawing(debug=False)
    drawing.add(group)  # type: ignore[arg-type]
    return ET.fromstring(drawing.tostring())


def _caption_of(row: ET.Element) -> ET.Element:
    text = row.find(f"{SVG}text")
    assert text is not None
    return text


def _assert_italic_caption(text: ET.Element) -> None:
    spans = text.findall(f"{SVG}tspan")
    assert [(span.text, span.get("font-style")) for span in spans] == [
        ("Danio rerio", "italic"), (" (zebrafish)", None)
    ]
    assert "<i>" not in "".join(text.itertext())


def _rows(root: ET.Element, name: str) -> list[ET.Element]:
    return [
        row for row in root.iter(f"{SVG}g")
        if (row.get("data-legend-name") or row.get("data-legend-key")) == name
    ]


def _color_table(tmp_path: Path) -> Path:
    table = tmp_path / "rules.tsv"
    table.write_text(f"CDS\tproduct\tcytochrome c oxidase subunit I\t#e7298a\t{DANIO}\n", encoding="utf-8")
    return table


def test_circular_cli_draws_a_color_rule_legend_name_in_italics(tmp_path: Path) -> None:
    circular_main([
        "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"), "-t", str(_color_table(tmp_path)),
        "-f", "svg", "--no_font_embedding", "-o", str(tmp_path / "circular"),
    ])
    root = ET.parse(tmp_path / "circular.svg").getroot()
    rows = _rows(root, DANIO)
    assert len(rows) == 1
    _assert_italic_caption(_caption_of(rows[0]))
    assert "&lt;i&gt;" not in (tmp_path / "circular.svg").read_text(encoding="utf-8").split('data-legend')[0]


def test_linear_cli_draws_a_color_rule_legend_name_in_italics(tmp_path: Path) -> None:
    linear_main([
        "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"), "-t", str(_color_table(tmp_path)),
        "-f", "svg", "--no_font_embedding", "-o", str(tmp_path / "linear"),
    ])
    root = ET.parse(tmp_path / "linear.svg").getroot()
    rows = _rows(root, DANIO)
    assert rows
    for row in rows:
        _assert_italic_caption(_caption_of(row))


def _solid(fill: str = "#54bcf8") -> dict[str, object]:
    return {"type": "solid", "fill": fill, "stroke": "none", "width": 0}


@pytest.mark.parametrize("side", ["right", "bottom"])
def test_circular_layout_reserves_the_drawn_runs(side: str) -> None:
    layout = build_circular_legend_layout(
        {"CDS": _solid(), DANIO: _solid("#e7298a")}, legend_position=side, canvas_width=2000.0,
        font_family=FAMILY, font_size=SIZE, color_rect_size=20.0,
    )
    widths = {entry.key: entry.text_width for entry in layout.solid_entries}
    assert widths[DANIO] == _run_sum(DANIO)
    assert widths["CDS"] == drawn_text_box("CDS", FAMILY, SIZE)[0]


def test_linear_layout_reserves_the_drawn_runs() -> None:
    plain = build_linear_legend_layout(
        {"Danio rerio (zebrafish)": _solid()}, legend_position="right", canvas_width=2000.0,
        font_family=FAMILY, font_size=SIZE, color_rect_size=20.0,
    )
    italic = build_linear_legend_layout(
        {DANIO: _solid()}, legend_position="right", canvas_width=2000.0,
        font_family=FAMILY, font_size=SIZE, color_rect_size=20.0,
    )
    plain_width = drawn_text_box("Danio rerio (zebrafish)", FAMILY, SIZE)[0]
    assert italic.vertical.feature.width - plain.vertical.feature.width == pytest.approx(_run_sum(DANIO) - plain_width)


def test_circular_ring_and_gradient_captions_draw_in_italics() -> None:
    table = {DANIO: {**RING_ROW}, "<i>Danio rerio</i> (zebrafish) identity": {**PAIRWISE_ROW}}
    measurement = _circular_measurement(table)
    root = _svg(CircularLegendGroup(measurement, table).get_group())
    ring = _rows(root, DANIO)
    assert len(ring) == 1
    _assert_italic_caption(_caption_of(ring[0]))
    gradient = _rows(root, "<i>Danio rerio</i> (zebrafish) identity")
    spans = _caption_of(gradient[0]).findall(f"{SVG}tspan")
    assert [(span.text, span.get("font-style")) for span in spans] == [
        ("Danio rerio", "italic"), (" (zebrafish) identity", None)
    ]
    layout = measurement.circular_layout
    assert layout is not None and layout.gradient is not None
    assert layout.gradient.width >= layout.gradient.rows[0].text_x + _run_sum(DANIO)


def test_linear_gradient_title_draws_in_italics() -> None:
    table = {DANIO: {**PAIRWISE_ROW}}
    root = _svg(LinearLegendGroup(_linear_measurement(table), table).get_group())
    rows = _rows(root, DANIO)
    assert rows
    for row in rows:
        _assert_italic_caption(_caption_of(row))


def test_white_space_collapses_across_runs_as_drawn() -> None:
    # A renderer draws one space where a run's space meets the next run's space.
    assert legend_caption_box("<i>a </i> b", SERIF, SIZE) == legend_caption_box("<i>a </i>b", SERIF, SIZE)
    assert legend_caption_box("<i>a</i>  b", SERIF, SIZE) == legend_caption_box("<i>a</i> b", SERIF, SIZE)
