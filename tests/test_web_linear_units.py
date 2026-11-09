"""The drawn units and unit pairs Python plans for the Web's comparison searches (OV-2738, D-01, D-37).

The Web searches the units of this plan and writes their indexes into the request, which Python reads as drawn-unit
indexes. The plan is the one the render and the CLI search use: the units are the resolved records, and the pairs of
a record pair are the pairs a comparisons-table ``source=losat`` row searches.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.api import GenBankInputSource
from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions, LosatSearchOptions
from gbdraw.api.record_planning import resolve_linear_options
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import LinearDiagramRequest, RecordInput
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import RecordWindow
from gbdraw.web_support.linear_units import plan_linear_units, plan_linear_units_json

REPO = Path(__file__).resolve().parents[1]
CASES = json.loads((REPO / "tests" / "fixtures" / "linear_unit_plan_cases.json").read_text(encoding="utf-8"))["cases"]
MITO = REPO / "tests" / "test_inputs" / "HmmtDNA.gbk"


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_python_plans_the_units_and_unit_pairs_of_the_shared_cases(case: dict) -> None:
    assert plan_linear_units(case["records"], case["pairs"]) == case["plan"]
    assert json.loads(plan_linear_units_json(json.dumps(case["records"]), json.dumps(case["pairs"]))) == case["plan"]


def test_the_web_plan_is_the_search_and_the_drawing_python_runs_for_the_same_records(tmp_path: Path) -> None:
    """Parity: the units are the drawn units of the request (keys in drawing order), and the unit pairs of a record
    pair are the pairs a ``source=losat`` comparisons-table row searches (CLI ``--losat``)."""

    windows = (RecordWindow(100, 900), RecordWindow(5000, 5600, True))
    request = LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(MITO), windows=windows, record_key="mito"),
            RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(1, 500),), record_key="other"),
        )
    )
    collection = resolve_request_records(request)
    table = tmp_path / "comparisons.tsv"
    table.write_text("blast\tquery\tsubject\tsource\n\t#1\t#2\tlosat\n", encoding="utf-8")
    searched = resolve_linear_options(
        LinearDiagramOptions(comparison_table_file=str(table), losat_search=LosatSearchOptions(program="losatn")),
        records=collection.records,
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1", "#3@2")),
    ).losat_search
    assert searched is not None and searched.pairs is not None

    plan = plan_linear_units(
        [
            {"recordKey": "mito", "windows": [
                {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None}
                for w in windows
            ]},
            {"recordKey": "other", "windows": [{"start": 1, "end": 500, "reverseComplement": False, "column": None}]},
        ],
        [{"queryRecordKey": "mito", "subjectRecordKey": "other", "source": "losat"}],
    )
    assert [unit["unitKey"] for unit in plan["units"]] == [item.record_key for item in collection.provenance]
    assert [unit["windowOf"] for unit in plan["units"]] == [
        item.window.window_of if item.window is not None else None for item in collection.provenance
    ]
    assert [(pair["queryUnit"], pair["subjectUnit"]) for pair in plan["pairs"]] == list(searched.pairs) == [
        (0, 2),
        (1, 2),
    ]


def test_an_uploaded_comparison_of_a_record_drawn_in_windows_is_refused() -> None:
    """A BLAST TSV gives the coordinates of one drawn unit: the Web must search such a pair (as the CLI refuses a
    comparisons-table file for it)."""

    case = CASES[0]
    upload = [{**case["pairs"][0], "source": "upload"}]
    with pytest.raises(ValidationError, match="drawn in several windows") as caught:
        plan_linear_units(case["records"], upload)
    assert caught.value.diagnostic == {"code": "COMPARISON_INPUT", "reason": "WINDOWED_FILE_COMPARISON"}


@pytest.mark.parametrize(
    ("records", "pairs"),
    [
        ([{"recordKey": "a", "windows": [{"start": 0, "end": 5}]}], []),
        ([{"recordKey": "a", "windows": []}, {"recordKey": "a", "windows": []}], []),
        ([{"recordKey": "a", "windows": []}], [{"queryRecordKey": "a", "subjectRecordKey": "z", "source": "losat"}]),
        ([{"recordKey": "a", "windows": []}], [{"queryRecordKey": "a", "subjectRecordKey": "a", "source": "table"}]),
    ],
    ids=["bad-window", "record-twice", "unknown-record", "unknown-source"],
)
def test_a_malformed_plan_request_is_a_helper_protocol_error(records: list, pairs: list) -> None:
    with pytest.raises(ValidationError) as caught:
        plan_linear_units(records, pairs)
    assert caught.value.diagnostic == {"code": "HELPER_PROTOCOL"}


def test_the_cli_and_the_web_write_one_protein_identity_for_each_window_unit() -> None:
    """A record analysis names the span its unit's search reads (``searchRegion``), on the CLI as on the Web, so a
    Session saved by either reuses the other's protein results; a record drawn whole names none (PD-OI-122)."""

    from gbdraw.api.request_render import _extract_linear_request_proteins, plan_linear_request
    from tests.test_origin_crossing_comparisons import _web_helpers

    windows = (RecordWindow(3000, 6000), RecordWindow(7000, 12000, True), RecordWindow(14000, 4300))
    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(MITO), windows=windows, record_key="mito"),
        RecordInput(GenBankInputSource(MITO), record_key="other"),
    ))
    plan = plan_linear_request(request)
    cli = _extract_linear_request_proteins(plan.request, plan.records, plan.inputs, plan.provenance).identity_manifest
    assert cli is not None
    units = plan_linear_units(
        [{"recordKey": "mito", "windows": [
            {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None} for w in windows
        ]}, {"recordKey": "other", "windows": []}],
        [],
    )["units"]
    helpers = _web_helpers()
    for index, unit in enumerate(units):
        window = unit["window"]
        web = json.loads(str(helpers["extract_cds_protein_fasta"](
            str(MITO), "genbank", None, None if window is None else window["searchRegion"], None, "0", index,
            unit["unitKey"], None, None, unit["windowOf"],
        )))
        assert "error" not in web, web.get("error")
        analysis = cli["recordAnalyses"][cli["recordInstances"][unit["unitKey"]]["recordAnalysisId"]]
        assert analysis["region"] == (None if window is None else window["searchRegion"]), unit["unitKey"]
        assert web["record_analysis_id"] == cli["recordInstances"][unit["unitKey"]]["recordAnalysisId"], unit["unitKey"]
        assert web["runtime_binding_hash"] == cli["recordInstances"][unit["unitKey"]]["runtimeBindingHash"]
        assert web["display_binding_hash"] == cli["recordInstances"][unit["unitKey"]]["displayBindingHash"]


def _two_record_genbank(tmp_path: Path) -> Path:
    """HmmtDNA twice, the second under its own ID: a file whose records a selector or an ALL input picks."""

    import copy

    from Bio import SeqIO

    first = SeqIO.read(MITO, "genbank")
    second = copy.deepcopy(first)
    second.id = second.name = "SECOND"
    path = tmp_path / "two.gb"
    SeqIO.write([first, second], path, "genbank")
    return path


@pytest.mark.parametrize(
    "case",
    ["reversed-window", "reversed-record", "selector", "all-split"],
)
def test_the_cli_and_the_web_write_one_record_analysis_and_display_binding(tmp_path: Path, case: str) -> None:
    """A protein record analysis is the Web's on the CLI: its display metadata is in the search frame (the crop before
    its reverse complement), and its selector names the record of the file the Web card names (OV-2855, OV-2856). A
    Session saved by either surface then reuses the other's derived protein results."""

    from gbdraw.api.request_render import _extract_linear_request_proteins, plan_linear_request
    from gbdraw.api.requests import RecordCardinality, RecordPresentation
    from gbdraw.io.record_select import parse_record_selector
    from tests.test_origin_crossing_comparisons import _web_helpers

    path = MITO
    # Each case: the CLI record input, then the Web card's helper arguments for each unit (path, region, selector).
    if case == "reversed-window":
        record = RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(7000, 12000, True),), record_key="mito")
        web_units = [("mito", "7000-12000", None)]
    elif case == "reversed-record":
        record = RecordInput(
            GenBankInputSource(MITO), presentation=RecordPresentation(reverse_complement=True), record_key="mito",
        )
        web_units = [("mito", None, None)]
    elif case == "selector":
        path = _two_record_genbank(tmp_path)
        record = RecordInput(GenBankInputSource(path), selector=parse_record_selector("#2"), record_key="b")
        web_units = [("b", None, "#2")]
    else:
        # An ALL input expands to `<key>:<n>`; the Web splits it into one card per record, selected by `#n`.
        path = _two_record_genbank(tmp_path)
        record = RecordInput(GenBankInputSource(path), cardinality=RecordCardinality.ALL, record_key="all")
        web_units = [("all:1", None, "#1"), ("all:2", None, "#2")]
    plan = plan_linear_request(LinearDiagramRequest(records=(record,)))
    cli = _extract_linear_request_proteins(plan.request, plan.records, plan.inputs, plan.provenance).identity_manifest
    assert cli is not None
    helpers = _web_helpers()
    for index, (unit_key, region, selector) in enumerate(web_units):
        web = json.loads(str(helpers["extract_cds_protein_fasta"](
            str(path), "genbank", None, region, selector, "0", index, unit_key, None, None, None,
        )))
        assert "error" not in web, web.get("error")
        web_binding = web["identity_manifest"]["recordInstances"][unit_key]
        cli_binding = cli["recordInstances"][unit_key]
        assert cli["recordAnalyses"][cli_binding["recordAnalysisId"]] == (
            web["identity_manifest"]["recordAnalyses"][web_binding["recordAnalysisId"]]
        ), unit_key
        assert cli_binding == web_binding, unit_key


def _card(key, windows, length, circular, row=None, column=None, reverse=False):
    return {
        "recordKey": key, "windows": [
            {"start": start, "end": end, "reverseComplement": rc, "column": col} for start, end, rc, col in windows
        ],
        "gridRow": row, "gridColumn": column, "reverseComplement": reverse, "sourceLength": length, "circular": circular,
    }


def test_card_window_plan_numbers_and_omitted_lengths_follow_the_drawing() -> None:
    """The card's circled numbers and "bp left out" come from the row owners the drawing uses (PD-OI-124, OV-2992):
    two windows of a circular record that overlap leave nothing out; a third window covering the interval leaves
    nothing out; after Up/Down the row order (not the list order) pairs the neighbours."""
    from gbdraw.web_support.linear_units import plan_card_windows, plan_card_windows_json

    plain = plan_card_windows([_card("a", [(1, 100, False, None), (301, 400, False, None)], 10_000, False),
                               _card("b", [], 5_000, False)])
    assert plain == {
        "windows": [{"recordKey": "a", "window": 0, "number": 1}, {"recordKey": "a", "window": 1, "number": 2}],
        "gaps": [{"recordKey": "a", "leftWindow": 0, "rightWindow": 1, "omittedBp": 200}],
    }
    overlap = plan_card_windows([_card("a", [(1, 600, False, None), (500, 900, False, None)], 1_000, True)])
    assert overlap["gaps"] == []
    covered = plan_card_windows([_card("a", [(1, 100, False, None), (301, 400, False, None), (50, 350, False, None)],
                                       10_000, False)])
    assert covered["gaps"] == []
    # Arranged in rows, Up on the third window: the row draws windows 1, 3, 2, and numbers them so.
    moved = plan_card_windows([_card("a", [(1, 100, False, 1), (201, 300, False, 3), (401, 500, False, 2)],
                                     10_000, False, row=1)])
    assert [(entry["window"], entry["number"]) for entry in moved["windows"]] == [(0, 1), (2, 2), (1, 3)]
    assert moved["gaps"] == []
    # A whole record placed first in the row numbers the windows after it, as the drawing does.
    shared = plan_card_windows([_card("w", [], 5_000, False, row=1, column=1),
                                _card("a", [(1, 100, False, 2), (301, 400, False, 3)], 10_000, False, row=1)])
    assert [entry["number"] for entry in shared["windows"]] == [2, 3]
    assert json.loads(plan_card_windows_json(json.dumps([_card("a", [(1, 100, False, None)], 500, False)]))) == {
        "windows": [], "gaps": []
    }


def test_card_window_plan_refuses_a_record_without_its_length() -> None:
    from gbdraw.web_support.linear_units import plan_card_windows

    record = _card("a", [(1, 100, False, None), (301, 400, False, None)], 10_000, False)
    record["sourceLength"] = None
    with pytest.raises(ValidationError):
        plan_card_windows([record])
