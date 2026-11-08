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
