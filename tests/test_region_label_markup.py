"""Region labels draw ``<i>`` markup as italics and measure the drawn runs (OV-2808).

A region label is user text like a Legend caption (OV-2760): ``italic_markup_runs`` reads its runs, the label's ink
(``label_ink_offsets``, which the lane room, placement and canvas read) measures those runs, and the drawers write
one ``tspan`` per run. Text outside the italic subset stays one literal run, also when a curved label is cut.
"""

from __future__ import annotations

from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    CoordinateSpan,
    LinearDiagramOptions,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.annotations.models import ResolvedRegionAnnotation
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.exceptions import LabelLeadersCrossWarning
from gbdraw.render.drawers.circular.annotations import circular_annotation_label_extent

SANS = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
SERIF = "'Times New Roman', 'Liberation Serif', serif"
MARKUP = "<i>stx2</i> phage"
SVG = "{http://www.w3.org/2000/svg}"


def _region(label: str, orientation: str = "auto", end: int = 300) -> RegionAnnotation:
    return RegionAnnotation(
        "r", CoordinateSpan(None, 100, end), label=label, style=RegionAnnotationStyle(label_orientation=orientation)
    )


def _resolved(label: str) -> ResolvedRegionAnnotation:
    return ResolvedRegionAnnotation(
        "r", "regions", 0, ((100, 300),), 200.0, 200, label, "bracket", None, RegionAnnotationStyle()
    )


def _circular_root(annotation: RegionAnnotation, length: int = 1000) -> ElementTree.Element:
    options = CircularDiagramOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),)),
        tracks=CircularRequestTrackOptions(
            circular_track_slots=("features:features", "regions:annotations@side=outside,set_id=regions,w=20px")
        ),
    )
    record = SeqRecord(Seq("A" * length), id="r1", name="r1")
    return ElementTree.fromstring(build_circular_diagram(record, options=options).tostring())


def _linear_root(annotation: RegionAnnotation) -> ElementTree.Element:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    options = LinearDiagramOptions(annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),)))
    return ElementTree.fromstring(build_linear_diagram([record], options=options).tostring())


def _label_element(root: ElementTree.Element) -> ElementTree.Element:
    """The ``<text>`` (or its ``<textPath>``) drawing the region's label."""

    (group,) = [element for element in root.iter() if element.attrib.get("data-gbdraw-annotation-id") == "r"]
    (text,) = [element for element in group.iter() if element.tag == f"{SVG}text"]
    paths = [element for element in text if element.tag == f"{SVG}textPath"]
    return paths[0] if paths else text


def _drawn_runs(element: ElementTree.Element) -> list[tuple[str, str]]:
    """``element``'s runs as drawn: its own text upright, then each ``tspan`` in its ``font-style``."""

    runs = [(element.text, "normal")] if element.text else []
    for child in element:
        assert child.tag == f"{SVG}tspan"
        runs.append((child.text or "", child.attrib.get("font-style", "normal")))
    return runs


@pytest.mark.parametrize("orientation", ["auto", "horizontal", "tangent", "radial"])
def test_circular_region_label_draws_italic_runs(orientation: str) -> None:
    label = _label_element(_circular_root(_region(MARKUP, orientation)))

    assert _drawn_runs(label) == [("stx2", "italic"), (" phage", "normal")]


def test_linear_region_label_draws_italic_runs() -> None:
    label = _label_element(_linear_root(_region(MARKUP, "horizontal")))

    assert _drawn_runs(label) == [("stx2", "italic"), (" phage", "normal")]


@pytest.mark.parametrize("text", ["Region", "a & <i>b</i>", "<b>x</b>", "<i>a</i"])
def test_region_label_outside_the_subset_stays_literal(text: str) -> None:
    for root in (_circular_root(_region(text)), _linear_root(_region(text, "horizontal"))):
        assert _drawn_runs(_label_element(root)) == [(text, "normal")]


@pytest.mark.parametrize("family", [SANS, SERIF])
@pytest.mark.parametrize("orientation", ["auto", "horizontal"])
def test_region_label_room_measures_the_drawn_runs(family: str, orientation: str) -> None:
    # The lane packer's measure of a label: the bp its letters reach about its place.
    def reach(label: str) -> float:
        extent = circular_annotation_label_extent(
            _resolved(label),
            RegionAnnotationStyle(label_orientation=orientation),
            12.0,
            200.0,
            4.0,
            font_family=family,
            record_length=1000,
            lanes_from_px=300.0,
            outward=True,
        )
        return extent.high_bp - extent.low_bp

    drawn, plain, literal = reach(MARKUP), reach("stx2 phage"), reach("[i]stx2[/i] phage")
    assert drawn == pytest.approx(plain, rel=0.1)
    assert drawn < 0.9 * literal


def test_cut_curved_region_label_keeps_italic_runs() -> None:
    text = "<i>" + "stx2 " * 300 + "</i>phage"
    with pytest.warns(LabelLeadersCrossWarning, match="region label cut"):
        label = _label_element(_circular_root(_region(text, "arc", end=200)))
    runs = _drawn_runs(label)

    assert label.tag == f"{SVG}textPath"
    # The cut keeps the italic start and ends with the ellipsis in that run.
    assert len(runs) == 1 and runs[0][1] == "italic" and runs[0][0].startswith("stx2 stx2")
    assert runs[0][0].endswith("…")


def test_cut_curved_literal_label_stays_literal() -> None:
    # The text is outside the subset (``&``); its cut start "<i>a</i> …" would be inside it, but stays literal.
    text = "<i>a</i> " + "x " * 1000 + "&"
    with pytest.warns(LabelLeadersCrossWarning, match="region label cut"):
        label = _label_element(_circular_root(_region(text, "arc", end=200)))
    runs = _drawn_runs(label)

    assert len(runs) == 1 and runs[0][1] == "normal" and runs[0][0].startswith("<i>a</i> x x")
    assert runs[0][0].endswith("…")
