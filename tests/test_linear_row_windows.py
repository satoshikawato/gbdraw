"""A Linear row of record windows: one name, a numbered list, edge and gap marks (PD-OI-122), the row text kept off
the comparison ribbons (OV-2712) and the scale bar in a figure with several records in a row (OV-2714)."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import copy
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import render_request
from gbdraw.api.assembly_inputs import LinearAssemblyInputs
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import LinearMultiRecordOptions
from gbdraw.api.record_planning import resolve_record_inputs
from gbdraw.api.requests import (
    GenBankInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
)
from gbdraw.features.placement import FeaturePlacementOverride, FeaturePlacementTarget
from gbdraw.features.source import capture_source_features
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.regions import RecordWindow
from gbdraw.layout.linear_row_windows import (
    LinearRowWindows,
    omitted_length_text,
    plan_row_windows,
    window_span_text,
)
from gbdraw.layout.record_units import RecordUnit, record_units
from gbdraw.linear_comparison import LinearComparison

REPO = Path(__file__).resolve().parents[1]
MITO = REPO / "tests" / "test_inputs" / "HmmtDNA.gbk"
SVG = "{http://www.w3.org/2000/svg}"


def _unit(index: int, span: tuple[int, int], *, key: str = "a", reverse: bool = False, label: str = "",
          circular: bool = False, length: int = 10_000, region: bool = True,
          record_id: str = "R1") -> tuple[RecordUnit, str]:
    return RecordUnit(index, key, record_id, span, reverse, length, circular, region, None), label


def _plan(units: list[tuple[RecordUnit, str]]) -> LinearRowWindows | None:
    """The row's plan from its units and their record labels."""
    return plan_row_windows([unit for unit, _ in units], [label for _, label in units])


# --- the owner: which rows are windowed, the name, the list, the marks ---------------------------------------


def test_a_row_of_whole_records_or_of_one_window_is_not_a_windowed_row() -> None:
    assert _plan([_unit(0, (1, 10_000), region=False), _unit(1, (1, 10_000), key="b", region=False)]) is None
    assert _plan([_unit(0, (100, 900))]) is None


def test_the_row_name_is_the_first_label_and_a_different_label_goes_to_its_line() -> None:
    plan = _plan(
        [
            _unit(0, (1, 500), key="c17", record_id="C17", reverse=True),
            _unit(1, (10, 900), key="c34", record_id="C34", label="<i>B. vel</i> SK19"),
            _unit(2, (1, 300), key="c39", record_id="C39", label="SK19 draft"),
            _unit(3, (5, 40), key="c32", record_id="C32", label="<i>B. vel</i> SK19"),
        ]
    )
    assert plan is not None
    assert plan.name_record_index == 1
    assert [(line.number, line.label, line.text) for line in plan.lines] == [
        (1, "", "C17:1-500 (−)"),
        (2, "", "C34:10-900"),
        (3, "SK19 draft", "C39:1-300"),
        (4, "", "C32:5-40"),
    ]
    # Windows of different records: no "//".
    assert plan.gaps == ()


def test_edges_are_bars_at_the_ends_of_a_linear_record_and_dashes_elsewhere() -> None:
    plan = _plan(
        [
            _unit(0, (1, 2_000)),
            _unit(1, (5_001, 10_000)),
            _unit(2, (1, 2_000), key="m", circular=True),
            _unit(3, (9_000, 10_000), key="r", reverse=True),
        ]
    )
    assert plan is not None
    # A reverse window's left edge is its last base.
    assert plan.edges == (("bar", "dash"), ("dash", "bar"), ("dash", "dash"), ("bar", "dash"))


@pytest.mark.parametrize(
    ("first", "second", "reverse", "omitted"),
    [
        ((1_000, 2_000), (3_541_001, 3_600_000), False, 3_539_000),
        ((1_000, 2_000), (2_001, 3_000), False, None),  # touching
        ((1_000, 2_000), (1_500, 3_000), False, None),  # overlapping
        ((5_000, 6_000), (1_000, 2_000), False, None),  # backwards
        ((5_000, 6_000), (1_000, 2_000), True, 2_999),  # both reversed: drawn left to right, source runs down
    ],
)
def test_the_gap_between_two_windows_of_one_record(first, second, reverse, omitted) -> None:
    plan = _plan(
        [_unit(0, first, reverse=reverse, length=4_000_000), _unit(1, second, reverse=reverse, length=4_000_000)]
    )
    assert plan is not None
    [gap] = plan.gaps
    assert gap.omitted_bp == omitted


def test_windows_in_opposite_directions_or_not_next_to_each_other_have_no_length() -> None:
    opposite = _plan([_unit(0, (1, 100)), _unit(1, (500, 900), reverse=True)])
    assert opposite is not None and [gap.omitted_bp for gap in opposite.gaps] == [None]
    # The third window of the record lies in the interval between the first two.
    apart = _plan([_unit(0, (1, 100)), _unit(1, (900, 1_000)), _unit(2, (400, 500))])
    assert apart is not None and [gap.omitted_bp for gap in apart.gaps] == [None, None]


def test_a_gap_across_the_origin_of_a_circular_record_counts_round_it() -> None:
    plan = _plan(
        [_unit(0, (9_000, 9_500), circular=True), _unit(1, (9_900, 200), circular=True)]
    )
    assert plan is not None and [gap.omitted_bp for gap in plan.gaps] == [399]


def test_the_texts_of_spans_and_omitted_lengths() -> None:
    assert window_span_text(1_906_335, 1_994_004) == "1,906,335-1,994,004"
    assert window_span_text(65_000, 4_000) == "65,000-4,000"
    assert [omitted_length_text(n) for n in (1_535_000, 12_345, 850)] == ["1.54 Mb", "12.3 kb", "850 bp"]



def test_the_unit_of_an_omitted_length_is_chosen_after_rounding() -> None:
    assert [omitted_length_text(n) for n in (999, 1_000, 999_949, 999_950, 999_999, 1_000_000)] == [
        "999 bp", "1.0 kb", "999.9 kb", "1.00 Mb", "1.00 Mb", "1.00 Mb",
    ]


@pytest.mark.parametrize(
    ("first", "second", "text"),
    [
        ((16_000, 500), (2_000, 3_000), "1.5 kb"),  # the first window crosses the origin
        ((15_000, 16_000), (1_000, 2_000), "1.6 kb"),  # the interval left out crosses the origin
    ],
)
def test_on_a_circular_record_the_interval_runs_forward_through_the_origin(first, second, text) -> None:
    plan = _plan(
        [_unit(0, first, circular=True, length=16_569), _unit(1, second, circular=True, length=16_569)]
    )
    assert plan is not None and [gap.text for gap in plan.gaps] == [text]


def test_on_a_linear_record_windows_that_run_backwards_have_no_length() -> None:
    plan = _plan([_unit(0, (15_000, 16_000), length=16_569), _unit(1, (1_000, 2_000), length=16_569)])
    assert plan is not None and [gap.text for gap in plan.gaps] == [""]


# --- the unit list: what each drawn record shows, read from the record planner's output ------------------------


def test_the_unit_list_reads_each_drawn_record_as_its_record_and_window() -> None:
    """One unit per drawn record: the record it draws (shared by its windows), the window span or the whole record,
    the orientation, the record's length and topology, and the window's place among its record's windows."""

    resolved = resolve_record_inputs(
        (
            RecordInput(
                InMemoryRecordSource(_linear_record()),
                windows=(RecordWindow(1, 2_000), RecordWindow(5_001, 8_000, True)),
                record_key="lin",
            ),
            RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(16_000, 500),), record_key="mito"),
            RecordInput(InMemoryRecordSource(_linear_record()), record_key="whole"),
        ),
        gff_candidate_features=None,
        gff_keep_all_features=True,
    )
    units = record_units(resolved.records, None)
    assert [
        (u.record_index, u.record_key, u.record_id, u.span, u.reverse, u.source_length, u.circular, u.is_window,
         u.window_index)
        for u in units
    ] == [
        (0, "lin", "LIN", (1, 2_000), False, 10_000, False, True, 0),
        (1, "lin", "LIN", (5_001, 8_000), True, 10_000, False, True, 1),
        (2, "mito", "NC_012920.1", (16_000, 500), False, 16_569, True, True, None),
        (3, "whole", "LIN", (1, 10_000), False, 10_000, False, False, None),
    ]
    assert [unit.length for unit in units] == [2_000, 3_000, 1_070, 10_000]
    # A reverse unit draws its last base first.
    assert units[1].edge_bases() == (8_000, 5_001)


# --- the drawing ------------------------------------------------------------------------------------------------


def _linear_record() -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * 2_500), id="LIN", name="LIN")
    record.annotations.update({"molecule_type": "DNA", "topology": "linear"})
    record.features = [
        SeqFeature(FeatureLocation(start, start + 600, strand=1), type="CDS", qualifiers={"locus_tag": [f"t{start}"]})
        for start in range(100, 9_000, 1_000)
    ]
    return record


def _windowed_figure(tmp_path: Path) -> ET.Element:
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                InMemoryRecordSource(_linear_record()),
                windows=(RecordWindow(1, 2_000, column=1), RecordWindow(5_001, 8_000, column=2)),
                record_key="lin",
                presentation=RecordPresentation(label="<i>Linear</i> one", subtitle="sub", grid_row=1),
            ),
            RecordInput(
                GenBankInputSource(MITO),
                windows=(RecordWindow(100, 900, column=1), RecordWindow(5_000, 5_600, True, column=2)),
                record_key="mito",
                presentation=RecordPresentation(grid_row=2),
            ),
        ),
        layout=LinearMultiRecordOptions(),
    )
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    return ET.fromstring(svg.read_text(encoding="utf-8"))


def _texts(group: ET.Element) -> list[str]:
    return ["".join(text.itertext()) for text in group.iter(f"{SVG}text")]


def test_a_windowed_row_shows_one_name_and_a_numbered_list_and_no_text_above_its_windows(tmp_path: Path) -> None:
    root = _windowed_figure(tmp_path)
    headings = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "record-definition-row"]
    lines = [
        ["".join(text.itertext()) for text in g.iter(f"{SVG}text") if text.get("data-definition-line-kind")]
        for g in headings
    ]
    assert lines == [
        ["Linear one", "LIN:1-2,000", "LIN:5,001-8,000"],
        ["NC_012920.1:100-900", "NC_012920.1:5,000-5,600 (−)"],
    ]
    for heading in headings:
        numbers = [g for g in heading.iter(f"{SVG}g") if g.get("data-gbdraw-row-mark") == "number"]
        assert [_texts(g) for g in numbers] == [["1"], ["2"]]
    above = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "record-definition"]
    assert above and all(_texts(g) == [] for g in above)


def test_a_windowed_row_draws_window_numbers_edge_marks_and_gap_marks(tmp_path: Path) -> None:
    root = _windowed_figure(tmp_path)
    marks = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "row-window-marks"]
    assert len(marks) == 2
    kinds = [
        [element.get("data-gbdraw-row-mark") for element in group if element.get("data-gbdraw-row-mark")]
        for group in marks
    ]
    # Linear record: a bar at base 1, dashes at the cuts inside it; "//" with 3,000 bp left out.
    assert kinds[0].count("bar") == 1 and kinds[0].count("dash") == 3
    # Circular record: dashes only; opposite directions draw "//" without a length.
    assert kinds[1].count("bar") == 0 and kinds[1].count("dash") == 4
    assert [_texts(g) for g in marks[0] if g.get("data-gbdraw-row-mark") == "gap"] == [["3.0 kb"]]
    assert [_texts(g) for g in marks[1] if g.get("data-gbdraw-row-mark") == "gap"] == [[]]
    assert [len([n for n in g if n.get("data-gbdraw-row-mark") == "number"]) for g in marks] == [2, 2]
    # "//" breaks the axis: drawn with the axis's stroke colour and width, as the dashed axis ends are (Owner 2026-11-09).
    slashes = [
        line for g in marks for gap in g if gap.get("data-gbdraw-row-mark") == "gap" for line in gap.iter(f"{SVG}line")
    ]
    dashes = [line for g in marks for line in g if line.get("data-gbdraw-row-mark") == "dash"]
    axis_strokes = {(line.get("stroke"), line.get("stroke-width")) for line in dashes}
    assert len(slashes) == 4 and len(axis_strokes) == 1
    assert {(line.get("stroke"), line.get("stroke-width")) for line in slashes} == axis_strokes
    # The white box behind the omitted length lies under "//", so it never hides the strokes' round ends.
    [labelled] = [gap for gap in marks[0] if gap.get("data-gbdraw-row-mark") == "gap"]
    assert [element.tag.removeprefix(SVG) for element in labelled] == ["rect", "line", "line", "text"]


def _comparison(query: int, subject: int) -> LinearComparison:
    row = [f"r{query + 1}", f"r{subject + 1}", 90.0, 100, 0, 0, 10, 900, 20, 910, 1e-20, 200]
    return LinearComparison(query, subject, pd.DataFrame([row], columns=COMPARISON_COLUMNS))


def _multi_record_figure() -> ET.Element:
    records = [SeqRecord(Seq("A" * 1000), id=f"r{index}") for index in range(1, 4)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
        record.features.append(SeqFeature(FeatureLocation(10, 900, strand=1), type="CDS"))
    canvas = assemble_linear_diagram_from_records(
        records,
        LinearAssemblyInputs(
            cfg=apply_config_overrides(
                None, {"labels.linear.scope": "none", "canvas.show_gc": False, "canvas.show_skew": False}
            ),
            linear_comparisons=[_comparison(0, 2), _comparison(1, 2)],
            legend="none",
        ),
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1", "#3@2")),
    )
    return ET.fromstring(canvas.tostring())


def _translate_y(element: ET.Element) -> float:
    return sum(float(y) for _x, y in re.findall(r"translate\(([-\d.]+),([-\d.]+)\)", element.get("transform", "")))


def test_text_above_the_records_of_a_row_stays_off_the_comparison_ribbons() -> None:
    root = _multi_record_figure()
    ribbon_bottom = max(
        _translate_y(group) + max(float(y) for y in re.findall(r"[ ,L]([-\d.]+)(?=[ L z]|$)", path.get("d", "")))
        for group in root.iter(f"{SVG}g")
        if group.get("data-subject-row") == "1"
        for path in group.iter(f"{SVG}path")
    )
    [header] = [
        g for g in root.iter(f"{SVG}g")
        if g.get("data-gbdraw-role") == "record-definition" and g.get("data-gbdraw-record-index") == "2"
    ]
    font_sizes = [float(text.get("font-size", "0")) for text in header.iter(f"{SVG}text")]
    header_top = _translate_y(header) + min(float(text.get("y", "0")) for text in header.iter(f"{SVG}text")) - (
        0.5 * max(font_sizes)
    )
    assert _texts(header)
    assert ribbon_bottom <= header_top


def test_a_figure_with_several_records_in_a_row_draws_the_scale_bar() -> None:
    root = _multi_record_figure()
    assert [g for g in root.iter(f"{SVG}g") if g.get("id") == "length_bar"]


def test_a_ruler_scale_in_a_figure_with_several_records_in_a_row_spans_one_record() -> None:
    # A ruler across the whole row would number the second record's bases on from the first (0..2 kb for two 1 kb
    # records); it spans the longest record, as in a figure with one record per row.
    records = [SeqRecord(Seq("A" * 1000), id=f"r{index}") for index in range(1, 3)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
        record.features.append(SeqFeature(FeatureLocation(10, 900, strand=1), type="CDS"))
    canvas = assemble_linear_diagram_from_records(
        records,
        LinearAssemblyInputs(
            cfg=apply_config_overrides(None, {"objects.scale.style": "ruler", "labels.linear.scope": "none"}),
            legend="none",
        ),
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")),
    )
    root = ET.fromstring(canvas.tostring())
    [bar] = [g for g in root.iter(f"{SVG}g") if g.get("id") == "length_bar"]
    widest = max(record["sequenceWidthPx"] for record in canvas._gbdraw_track_slot_geometry["records"])
    ruler_end = max(float(line.get("x2", "0")) for line in bar.iter(f"{SVG}line"))
    assert ruler_end == pytest.approx(widest)



def _mito_windows_figure(tmp_path: Path, windows: tuple[RecordWindow, ...]) -> ET.Element:
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                GenBankInputSource(MITO),
                windows=windows,
                record_key="mito",
                presentation=RecordPresentation(grid_row=1),
            ),
        ),
        layout=LinearMultiRecordOptions(),
    )
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    return ET.fromstring(svg.read_text(encoding="utf-8"))


def test_a_window_across_the_origin_is_listed_with_its_span_and_counts_the_gap_round_the_origin(tmp_path: Path) -> None:
    root = _mito_windows_figure(
        tmp_path, (RecordWindow(16_000, 500, column=1), RecordWindow(2_000, 3_000, column=2))
    )
    [heading] = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "record-definition-row"]
    assert [
        "".join(text.itertext()) for text in heading.iter(f"{SVG}text") if text.get("data-definition-line-kind")
    ] == ["NC_012920.1:16,000-500", "NC_012920.1:2,000-3,000"]
    [marks] = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "row-window-marks"]
    assert [e.get("data-gbdraw-row-mark") for e in marks if e.get("data-gbdraw-row-mark") in ("bar", "dash")] == [
        "dash", "dash", "dash", "dash",
    ]
    assert [_texts(g) for g in marks if g.get("data-gbdraw-row-mark") == "gap"] == [["1.5 kb"]]


def _record(record_id: str, length: int = 1_000) -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * (length // 4)), id=record_id, name=record_id)
    record.annotations.update({"molecule_type": "DNA", "topology": "linear"})
    record.features = [SeqFeature(FeatureLocation(10, 900, strand=1), type="CDS")]
    return record


def _window_input(record_id: str, column: int, row: int, label: str | None = None) -> RecordInput:
    return RecordInput(
        InMemoryRecordSource(_record(record_id)),
        windows=(RecordWindow(1, 950),),
        record_key=record_id.lower(),
        presentation=RecordPresentation(label=label, grid_row=row, grid_column=column),
    )


def _render(tmp_path: Path, request: LinearDiagramRequest) -> ET.Element:
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    return ET.fromstring(svg.read_text(encoding="utf-8"))


def test_the_row_heading_carries_the_identity_of_the_record_whose_label_it_shows(tmp_path: Path) -> None:
    root = _render(
        tmp_path,
        LinearDiagramRequest(
            records=(_window_input("A", 1, 1), _window_input("B", 2, 1, label="Bname")),
            layout=LinearMultiRecordOptions(),
        ),
    )
    [heading] = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "record-definition-row"]
    assert "Bname" in _texts(heading)
    # A click on the name opens B's card, and an edit there writes B's label.
    assert (heading.get("data-gbdraw-record-index"), heading.get("data-gbdraw-record-id")) == ("1", "B")


def test_ribbons_end_above_the_window_numbers_of_the_row_they_reach(tmp_path: Path) -> None:
    comparison = LinearComparison(
        0, 1, pd.DataFrame([["T", "A", 90.0, 700, 0, 0, 100, 800, 100, 800, 1e-20, 200]], columns=COMPARISON_COLUMNS)
    )
    second = LinearComparison(
        0, 2, pd.DataFrame([["T", "B", 90.0, 700, 0, 0, 100, 800, 100, 800, 1e-20, 200]], columns=COMPARISON_COLUMNS)
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                InMemoryRecordSource(_record("T", 2_000)),
                record_key="t",
                presentation=RecordPresentation(grid_row=1, grid_column=1),
            ),
            _window_input("A", 1, 2),
            _window_input("B", 2, 2),
        ),
        layout=LinearMultiRecordOptions(),
    )
    root = _render(
        tmp_path,
        replace(request, options=replace(request.options, linear_comparisons=(comparison, second))),
    )
    ribbon_bottom = max(
        _translate_y(group) + max(float(y) for y in re.findall(r"[ ,L]([-\d.]+)(?=[ L z]|$)", path.get("d", "")))
        for group in root.iter(f"{SVG}g")
        if group.get("data-subject-row") == "1"
        for path in group.iter(f"{SVG}path")
    )
    [marks] = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "row-window-marks"]
    circles = [
        _translate_y(marks) + float(circle.get("cy", "0")) - float(circle.get("r", "0"))
        for number in marks
        if number.get("data-gbdraw-row-mark") == "number"
        for circle in number.iter(f"{SVG}circle")
    ]
    assert len(circles) == 2
    assert ribbon_bottom <= min(circles)


def _translate(element: ET.Element) -> tuple[float, float]:
    pairs = re.findall(r"translate\(([-\d.]+),([-\d.]+)\)", element.get("transform", ""))
    return sum(float(x) for x, _y in pairs), sum(float(y) for _x, y in pairs)


def _path_points(group: ET.Element, path: ET.Element) -> list[tuple[float, float]]:
    dx, dy = _translate(group)
    return [(dx + float(x), dy + float(y)) for x, y in re.findall(r"([-\d.]+),([-\d.]+)", path.get("d", ""))]


def _stranded_record(record_id: str, length: int, features: list[tuple[int, int, int]]) -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * (length // 4)), id=record_id, name=record_id)
    record.annotations.update({"molecule_type": "DNA", "topology": "linear"})
    record.features = [SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS") for start, end, strand in features]
    return record


def _match(
    query: int, subject: int, query_id: str, subject_id: str, query_feature: str = ""
) -> LinearComparison:
    """One match of bases 100-400 of each side; ``query_feature`` binds its query end to that drawn feature."""

    row = [query_id, subject_id, 90.0, 300, 0, 0, 100, 400, 100, 400, 1e-20, 200]
    matches = pd.DataFrame([row], columns=COMPARISON_COLUMNS)
    if query_feature:
        matches["query_view_feature_svg_id"] = [query_feature]
    return LinearComparison(query, subject, matches)


# Drawn units of ``_row_between_ribbons_figure``: T, the row's windows C1 (reversed), C2a, C2b, and B.
_INTO_ROW = ((0, 1, "T", "C1"), (0, 2, "T", "C2"), (0, 3, "T", "C2"))
_OUT_OF_ROW = ((1, 4, "C1", "B"), (2, 4, "C2", "B"), (3, 4, "C2", "B"))


def _row_between_ribbons_figure(
    tmp_path: Path,
    comparisons: tuple[LinearComparison, ...] = tuple(_match(*pair) for pair in (*_INTO_ROW, *_OUT_OF_ROW)),
    overrides: dict[str, object] | None = None,
    feature_placements: tuple[FeaturePlacementOverride, ...] = (),
) -> ET.Element:
    """Strands drawn apart; a windowed middle row between two whole records, ribbons above and below it. The row
    holds a reversed window that reaches its contig's start (C1, every feature on the lower strand), then two
    windows of one record (C2): the first without a lower-strand feature, the second without an upper one."""

    request = LinearDiagramRequest(
        records=(
            RecordInput(
                InMemoryRecordSource(_stranded_record("T", 3_000, [(100, 600, 1), (1_200, 1_700, -1), (2_000, 2_600, 1)])),
                record_key="t",
                presentation=RecordPresentation(grid_row=1),
            ),
            RecordInput(
                InMemoryRecordSource(_stranded_record("C1", 1_000, [(100, 400, 1), (500, 900, 1)])),
                windows=(RecordWindow(1, 950, True),),
                record_key="c1",
                presentation=RecordPresentation(label="Row", grid_row=2),
            ),
            RecordInput(
                InMemoryRecordSource(_stranded_record("C2", 4_000, [(100, 600, 1), (2_100, 2_600, -1)])),
                windows=(RecordWindow(1, 1_000), RecordWindow(2_001, 3_000)),
                record_key="c2",
                presentation=RecordPresentation(grid_row=2),
            ),
            RecordInput(
                InMemoryRecordSource(_stranded_record("B", 3_000, [(100, 600, 1), (1_200, 1_700, -1)])),
                record_key="b",
                presentation=RecordPresentation(grid_row=3),
            ),
        ),
        layout=LinearMultiRecordOptions(),
    )
    options = replace(
        request.options,
        linear_comparisons=comparisons,
        config_overrides={
            "canvas.strandedness": True, "canvas.show_gc": False, "canvas.show_skew": False, **(overrides or {})
        },
        feature_placements=feature_placements,
    )
    return _render(tmp_path, replace(request, options=options))


def test_window_numbers_stand_at_the_upper_left_of_each_window_clear_of_ribbons_and_features(tmp_path: Path) -> None:
    """OV-3374: a window without an upper-strand feature put its number at axis height and drew the "//" length
    there, over the ribbons; the ribbons to it ended inside the row's upper lane, and those from a window without a
    lower-strand feature started inside its lower lane. The windows of a row share one top and one bottom: every
    number and length stands on the row's top, the ribbons end above them and start below the row's lower lane."""

    root = _row_between_ribbons_figure(tmp_path)
    [marks] = [g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-role") == "row-window-marks"]
    marks_x, marks_y = _translate(marks)
    circles = [
        (marks_x + float(circle.get("cx", "0")), marks_y + float(circle.get("cy", "0")), float(circle.get("r", "0")))
        for number in marks
        if number.get("data-gbdraw-row-mark") == "number"
        for circle in number.iter(f"{SVG}circle")
    ]
    [gap] = [g for g in marks if g.get("data-gbdraw-row-mark") == "gap"]
    [box] = gap.iter(f"{SVG}rect")
    [length] = gap.iter(f"{SVG}text")
    assert _texts(gap) == ["1.0 kb"]
    assert len(circles) == 3
    # One line for the row's numbers and its omitted length.
    assert len({round(cy, 3) for _cx, cy, _r in circles} | {round(marks_y + float(length.get("y", "0")), 3)}) == 1

    windows = {
        index: g
        for g in root.iter(f"{SVG}g")
        if g.get("data-gbdraw-slot-renderer") == "features" and (index := g.get("data-record-index")) in {"1", "2", "3"}
    }
    axis_y = {round(_translate(g)[1], 3) for g in windows.values()}
    assert len(axis_y) == 1
    features = [
        _path_points(g, path) for g in windows.values() for path in g.iter(f"{SVG}path") if path.get("data-gbdraw-feature-id")
    ]
    row_top = min(y for points in features for _x, y in points)
    row_bottom = max(y for points in features for _x, y in points)
    # Each number stands at its own window's start: left of it, right of the window before it.
    spans = [
        (_translate(windows[index])[0], _translate(windows[index])[0] + float(next(windows[index].iter(f"{SVG}line")).get("x2", "0")))
        for index in ("1", "2", "3")
    ]
    for number, (cx, _cy, r) in enumerate(circles):
        assert cx + r <= spans[number][0]
        assert number == 0 or cx - r >= spans[number - 1][1]
    for cx, cy, r in circles:
        # Above the axis, on the row's upper edge, and in no feature's box.
        assert row_top - r <= cy <= row_top < min(axis_y)
        for points in features:
            left, right = min(x for x, _y in points), max(x for x, _y in points)
            top, bottom = min(y for _x, y in points), max(y for _x, y in points)
            assert cx + r <= left or cx - r >= right or cy + r <= top or cy - r >= bottom

    def ribbon_points(attribute: str) -> list[list[tuple[float, float]]]:
        return [
            _path_points(g, path)
            for g in root.iter(f"{SVG}g")
            if g.get(attribute) == "1"
            for path in g.iter(f"{SVG}path")
        ]

    into_row = ribbon_points("data-subject-row")
    out_of_row = ribbon_points("data-query-row")
    assert len(into_row) == 3 and len(out_of_row) == 3
    # Every ribbon into the row ends at one y, above its numbers and the length's box.
    ends = {round(max(y for _x, y in points), 3) for points in into_row}
    assert len(ends) == 1
    assert max(ends) <= min(cy - r for _cx, cy, r in circles)
    assert max(ends) <= marks_y + float(box.get("y", "0"))
    # Every ribbon out of the row starts at one y, below the row's lower lane.
    starts = {round(min(y for _x, y in points), 3) for points in out_of_row}
    assert len(starts) == 1 and min(starts) >= row_bottom


def _ribbon_heights(root: ET.Element, subject_row: str) -> list[float]:
    return [
        max(y for _x, y in points) - min(y for _x, y in points)
        for g in root.iter(f"{SVG}g")
        if g.get("data-subject-row") == subject_row
        for path in g.iter(f"{SVG}path")
        for points in [_path_points(g, path)]
    ]


def test_ribbons_from_a_shallow_window_keep_the_comparison_height(tmp_path: Path) -> None:
    """The rows are spaced from the row's shared edge, where the ribbons attach: ribbons from the one window of the
    row without a lower-strand feature (C2a) are as tall as the set comparison height, not shorter by the depth of
    the row's other windows."""

    root = _row_between_ribbons_figure(
        tmp_path,
        comparisons=(*(_match(*pair) for pair in _INTO_ROW), _match(*_OUT_OF_ROW[1])),
        overrides={"canvas.linear.comparison_height": 200},
    )
    heights = _ribbon_heights(root, "2")
    assert heights and all(height == pytest.approx(200.0, abs=0.01) for height in heights)


def test_ribbons_bound_to_a_moved_feature_start_on_the_row_edge(tmp_path: Path) -> None:
    """Strands drawn on the axis. The one feature of window C2a is moved to the lane above the axis, so the ribbons
    of that window follow its features; a ribbon bound to the feature still starts on the row's lower edge, with
    the ribbons of the row's other windows, not under that lane (across the axis and the row's lower lane)."""

    source = _stranded_record("C2", 4_000, [(100, 600, 1), (2_100, 2_600, -1)])
    moved = FeaturePlacementOverride(
        "c2",
        capture_source_features(copy.deepcopy(source)).entries[0].biological_feature_id,
        FeaturePlacementTarget("lane", "above", 1),
    )
    on_axis = {"canvas.strandedness": False}
    plain = _row_between_ribbons_figure(tmp_path / "plain", overrides=on_axis, feature_placements=(moved,))
    [feature] = [
        path.get("data-gbdraw-stable-feature-id")
        for g in plain.iter(f"{SVG}g")
        if g.get("data-gbdraw-slot-renderer") == "features" and g.get("data-record-index") == "2"
        for path in g.iter(f"{SVG}path")
        if path.get("data-gbdraw-stable-feature-id")
    ]
    bound = (_match(1, 4, "C1", "B"), _match(2, 4, "C2", "B", query_feature=feature))
    root = _row_between_ribbons_figure(
        tmp_path / "bound",
        comparisons=(*(_match(*pair) for pair in _INTO_ROW), *bound),
        overrides=on_axis,
        feature_placements=(moved,),
    )
    starts = {
        g.get("data-query-record-index"): round(min(y for path in g.iter(f"{SVG}path") for _x, y in _path_points(g, path)), 3)
        for g in root.iter(f"{SVG}g")
        if g.get("data-query-row") == "1"
    }
    assert set(starts) == {"1", "2"}
    assert starts["2"] == starts["1"]


def _axis_gap_under_row_of_windows(tmp_path: Path, shallow_first: bool) -> float:
    """The axis gap from a row of two windows (C1 with upper-strand features only, C2 with lower-strand ones) to a
    short record S below it that lies under the first window only; strands drawn apart, no comparisons."""

    shallow = RecordInput(
        InMemoryRecordSource(_stranded_record("C1", 2_000, [(100, 400, 1), (500, 900, 1)])),
        windows=(RecordWindow(1, 1_000),),
        record_key="c1",
        presentation=RecordPresentation(label="Row", grid_row=1),
    )
    deep = RecordInput(
        InMemoryRecordSource(_stranded_record("C2", 4_000, [(100, 600, -1), (700, 900, -1)])),
        windows=(RecordWindow(1, 1_000),),
        record_key="c2",
        presentation=RecordPresentation(label="Row", grid_row=1),
    )
    short = RecordInput(
        InMemoryRecordSource(_stranded_record("S", 600, [(100, 500, 1)])),
        record_key="s",
        presentation=RecordPresentation(grid_row=2),
    )
    request = LinearDiagramRequest(
        records=(shallow, deep, short) if shallow_first else (deep, shallow, short),
        layout=LinearMultiRecordOptions(),
    )
    overrides = {"canvas.strandedness": True, "canvas.show_gc": False, "canvas.show_skew": False}
    root = _render(tmp_path, replace(request, options=replace(request.options, config_overrides=overrides)))
    axes = {
        g.get("data-gbdraw-record-id"): _translate(g)[1]
        for g in root.iter(f"{SVG}g")
        if g.get("data-gbdraw-slot-renderer") == "features"
    }
    assert axes["C1"] == axes["C2"]
    return axes["S"] - axes["C1"]


def test_a_row_of_windows_has_one_lower_edge_for_the_row_below(tmp_path: Path) -> None:
    """A row of windows has one lower edge, as it has one line of numbers and one edge for its ribbons: a record
    under its shallow window stands as far below it as under its deep window, with or without ribbons."""

    under_shallow = _axis_gap_under_row_of_windows(tmp_path / "shallow", shallow_first=True)
    under_deep = _axis_gap_under_row_of_windows(tmp_path / "deep", shallow_first=False)
    assert under_shallow == pytest.approx(under_deep, abs=1e-6)
