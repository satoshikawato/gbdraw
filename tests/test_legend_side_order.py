"""The Legend lists track rows and Annotations entries by side, in both modes (OV-3187, Owner D-62).

The outer side (outside, above) first, then the rows on the Axis, then the inner side (inside, below); the written
order within a side. An annotation overlay row is listed right after its anchor, and the overlays on one anchor in
their written order. This is the order the Web stack holds, so a CLI stack written out of side order draws the
Legend the Web draws for it.
"""

from __future__ import annotations

from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    CoordinateSpan,
    LinearDiagramOptions,
    LinearRequestTrackOptions,
    RegionAnnotation,
)
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.io.record_select import parse_record_selector

_SETS = {"a": "Alpha set", "b": "Beta set", "r": "Rho set"}


def _record(length: int) -> SeqRecord:
    record = SeqRecord(Seq("ACGGT" * (length // 5)), id="r1", name="r1")
    record.features = [SeqFeature(FeatureLocation(100, 900, strand=1), type="CDS")]
    return record


def _annotations() -> AnnotationOptions:
    return AnnotationOptions(
        sets=tuple(
            AnnotationSet(
                set_id,
                (RegionAnnotation(f"{set_id}1", CoordinateSpan(parse_record_selector("#1"), 200, 700),
                                  mark="highlight", legend_label=caption),),
            )
            for set_id, caption in _SETS.items()
        )
    )


def _legend_texts(drawing) -> list[str]:
    root = ElementTree.fromstring(drawing.tostring())
    legend = next(element for element in root.iter() if element.attrib.get("id") == "legend")
    return [text for element in legend.iter() if element.tag.endswith("text") for text in [("".join(element.itertext())).strip()] if text]


def _circular_legend(*slots: str) -> list[str]:
    return _legend_texts(build_circular_diagram(
        _record(20000),
        options=CircularDiagramOptions(annotations=_annotations(), tracks=CircularRequestTrackOptions(circular_track_slots=slots)),
    ))


def _linear_legend(*slots: str) -> list[str]:
    return _legend_texts(build_linear_diagram(
        [_record(5000)],
        options=LinearDiagramOptions(annotations=_annotations(), tracks=LinearRequestTrackOptions(linear_track_slots=slots)),
    ))


def _order(texts: list[str], names: tuple[str, ...]) -> list[str]:
    """The first entry of each of ``names`` (a name or its ``(+)``/``(-)`` rows), in the Legend's order."""

    listed: list[str] = []
    for text in texts:
        name = text.removesuffix(" (+)").removesuffix(" (-)")
        if name in names and name not in listed:
            listed.append(name)
    return listed


def test_circular_track_rows_follow_the_side_order() -> None:
    texts = _circular_legend("features:features", "gc:content", "sk:skew@side=outside")
    assert _order(texts, ("GC content", "GC skew")) == ["GC skew", "GC content"]


def test_linear_annotation_entries_follow_the_side_order() -> None:
    texts = _linear_legend("features:features", "b:annotations@set_id=b,side=below", "r:annotations@set_id=r,side=above")
    assert _order(texts, tuple(_SETS.values())) == ["Rho set", "Beta set"]


@pytest.mark.parametrize(
    ("legend", "slots", "expected"),
    [
        # An overlay written before its anchor is listed right after it, among the rows on the anchor's side.
        (_circular_legend, ("a:annotations@set_id=a,side=overlay,anchor_slot=features,z=1", "features:features",
                            "r:annotations@set_id=r,side=outside", "b:annotations@set_id=b,side=inside"),
         ["Rho set", "Alpha set", "Beta set"]),
        (_linear_legend, ("a:annotations@set_id=a,side=overlay,anchor_slot=features,z=1", "features:features",
                          "r:annotations@set_id=r,side=above", "b:annotations@set_id=b,side=below"),
         ["Rho set", "Alpha set", "Beta set"]),
        # Two overlays on one anchor keep their written order, right after the anchor.
        (_circular_legend, ("b:annotations@set_id=b,side=overlay,anchor_slot=gc,z=1", "features:features",
                            "r:annotations@set_id=r,side=inside", "gc:content",
                            "a:annotations@set_id=a,side=overlay,anchor_slot=gc,z=2"),
         ["Rho set", "Beta set", "Alpha set"]),
        (_linear_legend, ("b:annotations@set_id=b,side=overlay,anchor_slot=gc,z=1", "features:features",
                          "r:annotations@set_id=r,side=above", "gc:content@side=below",
                          "a:annotations@set_id=a,side=overlay,anchor_slot=gc,z=2"),
         ["Rho set", "Beta set", "Alpha set"]),
    ],
)
def test_overlay_entries_follow_their_anchor(legend, slots, expected) -> None:
    assert _order(legend(*slots), tuple(_SETS.values())) == expected


def _draw_circular_default(annotations: AnnotationOptions):
    return build_circular_diagram(_record(20000), options=CircularDiagramOptions(annotations=annotations))


def _draw_linear_default(annotations: AnnotationOptions):
    return build_linear_diagram([_record(5000)], options=LinearDiagramOptions(annotations=annotations))


@pytest.mark.parametrize("draw", [_draw_circular_default, _draw_linear_default], ids=["circular", "linear"])
def test_default_figure_lists_a_highlight_after_a_later_sets_row(draw) -> None:
    # No track stack: the default planner draws set A's highlights over the features and set B's band outside them,
    # so B's row is listed before A's highlights (Owner 2026-11-10, after seeing this change in default figures).
    span = parse_record_selector("#1")
    annotations = AnnotationOptions(sets=(
        AnnotationSet("A", (RegionAnnotation("a1", CoordinateSpan(span, 200, 700), mark="highlight",
                                             legend_label="A highlight"),)),
        AnnotationSet("B", (RegionAnnotation("b1", CoordinateSpan(span, 1200, 1700), mark="band",
                                             legend_label="B band"),)),
    ))
    assert _order(_legend_texts(draw(annotations)), ("A highlight", "B band")) == ["B band", "A highlight"]
