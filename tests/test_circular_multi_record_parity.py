"""A one-record Multi-Record Canvas equals the single-record Circular path.

The Multi-Record Canvas owns placement only (R8). Its slot geometry, legend
and center definition come from the same owners as the single-record path, and
"no depth input" (None) is distinct from "depth input without a cell for this
record" ([]).
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
    CoordinateSpan,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram

REPO_ROOT = Path(__file__).resolve().parents[1]
HMMT = REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
SVG_NS = {"svg": "http://www.w3.org/2000/svg"}

# The Web defaults for a fresh Circular Generate (Multi-Record Canvas on).
WEB_DEFAULT_OVERRIDES = {
    "canvas.circular.track_type": "tuckin",
    "canvas.strandedness": True,
    "canvas.show_gc": True,
    "canvas.show_skew": True,
    "labels.circular.scope": "none",
}


def _hmmt_record() -> SeqRecord:
    return next(SeqIO.parse(str(HMMT), "genbank"))


def _synthetic_record(record_id: str = "rec1", length: int = 2400) -> SeqRecord:
    record = SeqRecord(Seq("ATGCGC" * (length // 6)), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(FeatureLocation(0, length, strand=1), type="source",
                   qualifiers={"organism": ["Synthetic organism"]}),
        SeqFeature(FeatureLocation(100, 400, strand=1), type="CDS",
                   qualifiers={"product": ["alpha"]}),
        SeqFeature(FeatureLocation(700, 1000, strand=-1), type="CDS",
                   qualifiers={"product": ["beta"]}),
        SeqFeature(FeatureLocation(1300, 1370, strand=1), type="tRNA",
                   qualifiers={"product": ["tRNA-Ala"]}),
    ]
    return record


def _constant_depth_table(reference: str, length: int, depth: float = 20.0) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [reference] * length,
            "position": list(range(1, length + 1)),
            "depth": [depth] * length,
        }
    )


def _hit(subject: str, send: int) -> tuple[object, ...]:
    return ("query1", subject, 90.0, 20, 0, 0, 1, 20, 1, send, 1e-20, 100.0)


def _slot_geometry(drawing) -> list[tuple[str, str, float, float]]:
    geometry = getattr(drawing, "_gbdraw_track_slot_geometry")
    return [
        (
            str(slot["slotId"]),
            str(slot["renderer"]),
            round(float(slot["widthPx"]), 3),
            round(float(slot["radiusFactor"]), 5),
        )
        for record in geometry["records"]
        for slot in record["slots"]
    ]


def _legend_rows(drawing) -> list[tuple[str, str]]:
    root = ET.fromstring(drawing.tostring())
    gradients = {
        element.get("id"): tuple(
            stop.get("stop-color") for stop in element.findall("svg:stop", SVG_NS)
        )
        for element in root.iter(f"{{{SVG_NS['svg']}}}linearGradient")
    }
    rows: list[tuple[str, str]] = []
    for entry in root.findall(".//svg:g[@data-legend-key]", SVG_NS):
        fill = next(
            (
                str(path.get("fill"))
                for path in entry.iter(f"{{{SVG_NS['svg']}}}path")
                if path.get("fill") not in (None, "none")
            ),
            "none",
        )
        if fill.startswith("url(#"):
            fill = "gradient:" + ",".join(gradients.get(fill[5:-1], ()))
        rows.append((str(entry.get("data-legend-key")), fill.lower()))
    return rows


def _definition_lines(drawing) -> list[str]:
    root = ET.fromstring(drawing.tostring())
    group = root.find(".//svg:g[@data-gbdraw-role='record-definition']", SVG_NS)
    assert group is not None
    return ["".join(text.itertext()) for text in group.findall("svg:text", SVG_NS)]


def _build_both(record: SeqRecord, options: CircularDiagramOptions):
    single = build_circular_diagram(record, options=options)
    grid = build_circular_multi_diagram(
        [record],
        options=options,
        layout=CircularMultiRecordOptions(multi_record_positions=["#1@1"]),
    )
    return single, grid


def _web_default_options(**kwargs) -> CircularDiagramOptions:
    overrides = dict(WEB_DEFAULT_OVERRIDES)
    overrides.update(kwargs.pop("config_overrides", {}))
    return CircularDiagramOptions(
        config_overrides=overrides,
        output=CircularOutputOptions(legend="left"),
        **kwargs,
    )


def _parity_cases() -> list:
    record = _synthetic_record()
    slots = (
        "features:features",
        "ticks:ticks",
        "gc_content:dinucleotide_content@legend_label=MY GC",
        "gc_skew:dinucleotide_skew@positive_color=#ff0000,negative_color=#0000ff",
        "at_skew:dinucleotide_skew@nt=AT,w=20px",
    )
    annotation = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "a1",
                        CoordinateSpan(None, 200, 900),
                        mark="bracket",
                        legend_label="MY REGION",
                        style=RegionAnnotationStyle(stroke="#00aa00"),
                    ),
                ),
            ),
        )
    )
    return [
        pytest.param(_hmmt_record(), _web_default_options(), id="web-default-hmmt"),
        pytest.param(
            record,
            _web_default_options(
                tracks=CircularRequestTrackOptions(circular_track_slots=slots)
            ),
            id="custom-slots",
        ),
        pytest.param(record, _web_default_options(annotations=annotation), id="annotation-legend-label"),
        pytest.param(
            record,
            _web_default_options(
                depth_track_tables=[[_constant_depth_table("rec1", len(record.seq))]],
                depth_track_labels=["Sample A"],
                tracks=CircularRequestTrackOptions(
                    circular_track_slots=(
                        "features:features",
                        "ticks:ticks",
                        "depth:depth@track_index=0,legend_label=MY DEPTH",
                        "gc_content:dinucleotide_content",
                    )
                ),
                window=100,
                step=100,
                depth_window=100,
                depth_step=100,
            ),
            id="depth-slot-legend-label",
        ),
        pytest.param(
            record,
            _web_default_options(
                conservation_dataframes=[
                    pd.DataFrame([_hit("rec1", len(record.seq))], columns=COMPARISON_COLUMNS)
                ],
                conservation_reference="subject",
                conservation_labels=["barcode13"],
                conservation_colors=["#E15759"],
            ),
            id="colored-conservation",
        ),
    ]


@pytest.mark.circular
@pytest.mark.parametrize(("record", "options"), _parity_cases())
def test_one_record_multi_record_canvas_matches_single_record_path(record, options) -> None:
    single, grid = _build_both(record, options)

    assert _slot_geometry(grid) == _slot_geometry(single)
    assert _legend_rows(grid) == _legend_rows(single)
    assert _definition_lines(grid) == _definition_lines(single)


@pytest.mark.circular
def test_multi_record_canvas_without_depth_input_reserves_no_depth_slot() -> None:
    records = [_synthetic_record("rec1"), _synthetic_record("rec2")]
    drawing = build_circular_multi_diagram(records, options=_web_default_options())

    geometry = getattr(drawing, "_gbdraw_track_slot_geometry")
    assert len(geometry["records"]) == 2
    for record_geometry in geometry["records"]:
        assert "depth" not in {slot["slotId"] for slot in record_geometry["slots"]}


@pytest.mark.circular
def test_web_default_one_record_gc_widths_match_single_record() -> None:
    single, grid = _build_both(_hmmt_record(), _web_default_options())

    def widths(drawing) -> dict[str, float]:
        return {
            slot_id: width
            for slot_id, _renderer, width, _radius in _slot_geometry(drawing)
            if slot_id in {"gc_content", "gc_skew"}
        }

    assert widths(grid) == widths(single)
    assert widths(grid)["gc_content"] == pytest.approx(74.1, abs=0.05)


@pytest.mark.circular
def test_multi_record_canvas_legend_includes_slot_rows_for_every_record() -> None:
    records = [_synthetic_record("rec1"), _synthetic_record("rec2")]
    drawing = build_circular_multi_diagram(
        records,
        options=_web_default_options(
            tracks=CircularRequestTrackOptions(
                circular_track_slots=(
                    "features:features",
                    "ticks:ticks",
                    "gc_content:dinucleotide_content@legend_label=MY GC",
                    "at_skew:dinucleotide_skew@nt=AT,w=20px",
                )
            )
        ),
    )

    keys = [key for key, _fill in _legend_rows(drawing)]
    assert "MY GC" in keys
    assert "AT skew (+)" in keys and "AT skew (-)" in keys
    assert "GC content" not in keys
