"""The axis geometry each record's track-slot metadata carries (Select region).

Python writes, per drawn record, where its axis lies in the frame of the
record's drawn group and which source base each drawn base shows
(``gbdraw.render.track_slot_metadata``). The Web's pointer converter
(``gbdraw/web/js/services/region-axis.js``) reads the same vectors
(``tests/web/region-axis.test.mjs``): each probe is a point in the record's
frame and the display index and source base Python says it shows.

Rewrite the vectors with ``python tests/test_region_axis_geometry.py`` after an
intended change.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    CircularMultiRecordOptions,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
)
from gbdraw.io.regions import parse_region_spec

INPUTS = Path(__file__).parent / "test_inputs"
LAMBDA = INPUTS / "NC_001416.gb"
MITO = INPUTS / "HmmtDNA.gbk"
VECTORS = Path(__file__).parent / "fixtures" / "region_axis_geometry_vectors.json"
SVG = "{http://www.w3.org/2000/svg}"


def _linear(*records: RecordInput, overrides: dict[str, object] | None = None) -> Any:
    return plan_request(
        LinearDiagramRequest(records=records, options=LinearDiagramOptions(config_overrides=overrides or {}))
    ).build().drawing


def _circular(*records: RecordInput, grid: bool = False) -> Any:
    return plan_request(
        CircularDiagramRequest(records=records, layout=CircularMultiRecordOptions() if grid else None)
    ).build()


def _drawings() -> dict[str, Any]:
    return {
        "linear-single": _linear(RecordInput(GenBankInputSource(str(LAMBDA)), record_key="lambda")),
        "linear-wrap-reverse": _linear(
            RecordInput(
                GenBankInputSource(str(LAMBDA)),
                record_key="lambda",
                presentation=RecordPresentation(reverse_complement=True),
            ),
            overrides={"canvas.linear.wrap_length": 20_000},
        ),
        "linear-two-records-crop-reverse": _linear(
            RecordInput(
                GenBankInputSource(str(LAMBDA)),
                record_key="crop",
                region=parse_region_spec("NC_001416.1:1001-9000:rc"),
            ),
            RecordInput(GenBankInputSource(str(LAMBDA)), record_key="whole"),
        ),
        "circular-display-start": _circular(
            RecordInput(
                GenBankInputSource(str(MITO)),
                record_key="mito",
                display=RecordDisplayOptions(start_coordinate=16_000),
            )
        ),
        "circular-multi-record-canvas": _circular(
            RecordInput(GenBankInputSource(str(MITO)), record_key="mito"),
            RecordInput(GenBankInputSource(str(LAMBDA)), record_key="lambda", display=RecordDisplayOptions(is_circular=True)),
            grid=True,
        ),
    }


def _source_base(source: dict[str, int], index: int) -> int:
    """The payload's documented map from a display index to a source base."""

    return (source["anchorBase"] - 1 + source["step"] * index) % source["sourceLength"] + 1


def _probe_indexes(drawn_length: int) -> list[int]:
    return sorted({0, 1, drawn_length // 3, drawn_length // 2, drawn_length - 1})


def _probes(axis: dict[str, Any]) -> list[dict[str, Any]]:
    """Points at the middle of a few drawn bases, with the base each shows."""

    source = axis["source"]
    probes: list[dict[str, Any]] = []
    for index in _probe_indexes(source["drawnLength"]):
        if axis["kind"] == "circular":
            radius = (axis["bandInnerPx"] + axis["bandOuterPx"]) / 2
            angle = math.radians(axis["startAngleDeg"] + axis["degreesPerBp"] * (index + 0.5))
            point = [radius * math.cos(angle), radius * math.sin(angle)]
        else:
            row = next(row for row in axis["rows"] if row["displayStart"] <= index < row["displayEnd"])
            px_per_bp = (row["xEndPx"] - row["xStartPx"]) / (row["displayEnd"] - row["displayStart"])
            point = [row["xStartPx"] + px_per_bp * (index - row["displayStart"] + 0.5), row["yPx"]]
        probes.append(
            {"point": [round(value, 6) for value in point], "displayIndex": index, "sourceBase": _source_base(source, index)}
        )
    return probes


def compute_vectors() -> list[dict[str, Any]]:
    cases = []
    for name, drawing in _drawings().items():
        records = drawing._gbdraw_track_slot_geometry["records"]
        cases.append(
            {
                "name": name,
                "records": [
                    {
                        "recordIndex": record["recordIndex"],
                        "recordKey": record["recordKey"],
                        "axis": record["axis"],
                        "probes": _probes(record["axis"]),
                    }
                    for record in records
                ],
            }
        )
    return cases


def _case(name: str) -> list[dict[str, Any]]:
    return next(case["records"] for case in compute_vectors() if case["name"] == name)


def test_axis_vectors_match_python() -> None:
    stored = json.loads(VECTORS.read_text(encoding="utf-8"))["cases"]
    assert stored == compute_vectors(), (
        "Python's axis geometry changed; rewrite the vectors with `python tests/test_region_axis_geometry.py`"
    )


def test_source_bases_are_the_input_file_coordinates() -> None:
    """Domain facts: the source base each end of the drawn axis shows."""

    def ends(name: str) -> list[tuple[str, int, int]]:
        return [
            (
                record["recordKey"],
                _source_base(record["axis"]["source"], 0),
                _source_base(record["axis"]["source"], record["axis"]["source"]["drawnLength"] - 1),
            )
            for record in _case(name)
        ]

    assert ends("linear-single") == [("lambda", 1, 48_502)]
    # A reverse-complemented record reads from its last base.
    assert ends("linear-wrap-reverse") == [("lambda", 48_502, 1)]
    assert ends("linear-two-records-crop-reverse") == [("crop", 9_000, 1_001), ("whole", 1, 48_502)]
    # Display start 16,000: the drawn axis wraps through the origin.
    assert ends("circular-display-start") == [("mito", 16_000, 15_999)]
    mito = _case("circular-display-start")[0]["axis"]["source"]
    assert _source_base(mito, 569) == 16_569 and _source_base(mito, 570) == 1


def _translate(element: ElementTree.Element) -> tuple[float, float]:
    total = [0.0, 0.0]
    for x, y in re.findall(r"translate\(\s*([-0-9.e]+)[ ,]+([-0-9.e]+)\s*\)", element.get("transform", "")):
        total[0] += float(x)
        total[1] += float(y)
    return total[0], total[1]


def test_linear_rows_lie_on_the_drawn_axis_lines() -> None:
    drawing = _drawings()["linear-wrap-reverse"]
    root = ElementTree.fromstring(drawing.tostring())
    rows = drawing._gbdraw_track_slot_geometry["records"][0]["axis"]["rows"]
    drawn = []
    for group in root.iter(f"{SVG}g"):
        if "data-gbdraw-wrap-row" not in group.attrib:
            continue
        line = next(group.iter(f"{SVG}line"), None)
        if line is not None and line.get("y1") == "0":
            drawn.append((_translate(group)[1], float(line.get("x2", "nan"))))
    assert sorted(set(drawn)) == [(row["yPx"], row["xEndPx"]) for row in rows]
    assert [(row["displayStart"], row["displayEnd"]) for row in rows] == [(0, 20_000), (20_000, 40_000), (40_000, 48_502)]


def test_circular_axis_radius_is_the_drawn_axis_circle() -> None:
    drawing = _drawings()["circular-multi-record-canvas"]
    root = ElementTree.fromstring(drawing.tostring())
    radii = sorted(
        float(circle.get("r", "nan"))
        for circle in root.iter(f"{SVG}circle")
        if circle.get("stroke") == "gray" and circle.get("cx") == "0"
    )
    records = drawing._gbdraw_track_slot_geometry["records"]
    assert sorted(record["axis"]["radiusPx"] for record in records) == radii
    assert [record["recordKey"] for record in records] == ["mito", "lambda"]
    for record in records:
        axis = record["axis"]
        assert axis["bandInnerPx"] < axis["bandOuterPx"]
        assert axis["startAngleDeg"] == -90.0


if __name__ == "__main__":
    VECTORS.write_text(json.dumps({"cases": compute_vectors()}, indent=1) + "\n", encoding="utf-8")
