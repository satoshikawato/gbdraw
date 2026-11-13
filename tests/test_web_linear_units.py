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
import gbdraw.api.request_render as request_render_module
from gbdraw.api.request_render import LinearRequestPlan, plan_request, resolve_request_records
from gbdraw.api.requests import LinearDiagramRequest, RecordInput, RecordPresentation
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import RecordWindow
from gbdraw.layout.linear_multi_record import grid_display_pairs
from gbdraw.layout.record_placement import resolve_record_row_positions
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
            {"recordKey": "mito", "gridRow": None, "windows": [
                {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None}
                for w in windows
            ]},
            {
                "recordKey": "other", "gridRow": None,
                "windows": [{"start": 1, "end": 500, "reverseComplement": False, "column": None}],
            },
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


_ONE_WINDOW = (RecordWindow(100, 900),)
_TWO_WINDOWS = (RecordWindow(100, 900), RecordWindow(5000, 5600))


@pytest.mark.parametrize(
    ("mito_windows", "rows", "grid"),
    [
        (_ONE_WINDOW, (1, 1, 2), True),
        (_ONE_WINDOW, (2, 1, 1), True),
        (_ONE_WINDOW, (3, 1, 2), True),
        (_TWO_WINDOWS, (None, None, None), True),
        (_ONE_WINDOW, (None, None, None), False),
        (_ONE_WINDOW, (1, 2, 3), False),
    ],
    ids=["two-above-one", "one-above-two", "rows-out-of-order", "two-windows-share-a-row", "no-rows", "one-per-row"],
)
def test_the_web_display_pairs_are_the_pairs_the_drawing_and_the_derived_identity_read(
    mito_windows: tuple[RecordWindow, ...], rows: tuple, grid: bool
) -> None:
    """Parity (OV-3248): the unit pairs the Web's unit plan says the drawing displays are the pairs the Python drawing
    displays for the same request (``grid_display_pairs`` over its rows), and they decide the derived LOSATP
    identity's ``explicitDisplayPairs``: a grid arranged in the Web (rows, no columns), rows out of record order, and
    a record drawn in two windows, which share its row. One record per row in record order has no grid pairs."""

    request = LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(MITO), windows=mito_windows, record_key="mito",
                        presentation=RecordPresentation(grid_row=rows[0])),
            RecordInput(GenBankInputSource(MITO), record_key="other", presentation=RecordPresentation(grid_row=rows[1])),
            RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(1, 500),), record_key="third",
                        presentation=RecordPresentation(grid_row=rows[2])),
        ),
        layout=LinearMultiRecordOptions(),
    )
    plan = plan_request(request)
    assert isinstance(plan, LinearRequestPlan)
    drawn_rows = resolve_record_row_positions(
        plan.records, plan.layout.multi_record_positions if plan.layout is not None else None
    )[1]
    web = plan_linear_units(
        [
            {"recordKey": key, "gridRow": row, "windows": [
                {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None}
                for w in record_windows
            ]}
            for key, row, record_windows in (
                ("mito", rows[0], mito_windows), ("other", rows[1], ()), ("third", rows[2], (RecordWindow(1, 500),))
            )
        ],
        [],
    )
    drawn = grid_display_pairs(drawn_rows)
    assert (drawn is not None) is grid
    assert web["gridPairs"] == (
        None if drawn is None else [{"queryUnit": query, "subjectUnit": subject} for query, subject in drawn]
    )
    settings = request_render_module._request_derived_settings(plan.request, "collinear", drawn_rows)
    assert settings.explicit_display_pairs is grid


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
        ([{"recordKey": "a", "windows": [{"start": 0, "end": 5}], "gridRow": None}], []),
        ([{"recordKey": "a", "windows": [], "gridRow": None}, {"recordKey": "a", "windows": [], "gridRow": None}], []),
        (
            [{"recordKey": "a", "windows": [], "gridRow": None}],
            [{"queryRecordKey": "a", "subjectRecordKey": "z", "source": "losat"}],
        ),
        (
            [{"recordKey": "a", "windows": [], "gridRow": None}],
            [{"queryRecordKey": "a", "subjectRecordKey": "a", "source": "table"}],
        ),
        ([{"recordKey": "a", "windows": []}], []),
        ([{"recordKey": "a", "windows": [], "gridRow": 0}], []),
        ([{"recordKey": "a", "windows": [], "gridRow": 1}, {"recordKey": "b", "windows": [], "gridRow": None}], []),
    ],
    ids=["bad-window", "record-twice", "unknown-record", "unknown-source", "no-grid-row", "bad-grid-row", "partial-rows"],
)
def test_a_malformed_plan_request_is_a_helper_protocol_error(records: list, pairs: list) -> None:
    with pytest.raises(ValidationError) as caught:
        plan_linear_units(records, pairs)
    assert caught.value.diagnostic == {"code": "HELPER_PROTOCOL"}


@pytest.mark.parametrize(
    ("record", "message"),
    [({"recordKey": "a", "windows": []}, "has no gridRow"), ({"recordKey": "a", "windows": [], "gridRow": 0}, "invalid gridRow")],
    ids=["missing", "zero"],
)
def test_a_plan_request_names_a_missing_grid_row_apart_from_an_invalid_one(record: dict, message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        plan_linear_units([record], [])


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
        [{"recordKey": "mito", "gridRow": None, "windows": [
            {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None} for w in windows
        ]}, {"recordKey": "other", "gridRow": None, "windows": []}],
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


def test_the_cli_and_the_web_search_one_side_fasta_for_a_reversed_window_of_a_multi_record_side(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Parity (OV-3457): a LOSATP side of several records of one File, one of them a window drawn reversed, is the
    same FASTA on the CLI as on the Web: each record's proteins in its search frame (the crop before its reverse
    complement, as the Web helper extracts them), so the CLI's batches carry the searchContext keys the Web stored."""

    import pandas as pd

    import gbdraw.analysis.protein_colinearity as protein_colinearity
    from gbdraw.api.request_render import _extract_linear_request_proteins, plan_linear_request
    from gbdraw.comparisons.losat_jobs import losat_record_uids
    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    from tests.test_origin_crossing_comparisons import _web_helpers

    other = _two_record_genbank(tmp_path)
    windows = (RecordWindow(3000, 6000), RecordWindow(7000, 12000, True))
    plan = plan_linear_request(LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(MITO), windows=windows, record_key="mito"),
        RecordInput(GenBankInputSource(other), selector=parse_record_selector("#2"), record_key="other"),
    )))
    extraction = _extract_linear_request_proteins(plan.request, plan.records, plan.inputs, plan.provenance)
    units = plan_linear_units(
        [{"recordKey": "mito", "gridRow": None, "windows": [
            {"start": w.start, "end": w.end, "reverseComplement": w.reverse_complement, "column": None} for w in windows
        ]}, {"recordKey": "other", "gridRow": None, "windows": []}],
        [],
    )["units"]
    helpers = _web_helpers()
    web_fasta: dict[str, str] = {}
    for index, unit in enumerate(units):
        window = unit["window"]
        path, selector = (MITO, None) if unit["windowOf"] else (other, "#2")
        web = json.loads(str(helpers["extract_cds_protein_fasta"](
            str(path), "genbank", None, None if window is None else window["searchRegion"], selector, "0", index,
            unit["unitKey"], None, None, unit["windowOf"],
        )))
        assert "error" not in web, web.get("error")
        web_fasta[unit["unitKey"]] = web["fasta"]

    prepare = protein_colinearity.prepare_losat_batches
    seen: list[tuple[object, object]] = []

    def capture(jobs, *, uids, record_fasta, protein=False):
        cli = prepare(jobs, uids=uids, record_fasta=record_fasta, protein=protein)
        web = prepare(jobs, uids=uids, record_fasta=lambda index: web_fasta[str(uids[index])], protein=protein)
        seen.append((cli, web))
        return cli

    monkeypatch.setattr(protein_colinearity, "prepare_losat_batches", capture)
    uids = losat_record_uids(plan.records)
    pairs = [(uids.index(query), uids.index("other")) for query in uids if query != "other"]
    protein_colinearity.build_pairwise_protein_blastp_comparisons(
        plan.records, protein_extraction=extraction, pairs=pairs,
        runner=lambda _query, _subject: pd.DataFrame(columns=COMPARISON_COLUMNS),
    )
    assert seen
    for cli_batches, web_batches in seen:
        assert any(len(batch.query.indexes) > 1 for batch in cli_batches)
        for cli_batch, web_batch in zip(cli_batches, web_batches, strict=True):
            assert cli_batch.query.fasta == web_batch.query.fasta
            assert cli_batch.subject.fasta == web_batch.subject.fasta
            assert cli_batch.search_context == web_batch.search_context


def _card(key, windows, length, circular, row=None, column=None, reverse=False, drawn=1):
    return {
        "recordKey": key, "windows": [
            {"start": start, "end": end, "reverseComplement": rc, "column": col} for start, end, rc, col in windows
        ],
        "gridRow": row, "gridColumn": column, "reverseComplement": reverse, "sourceLength": length, "circular": circular,
        "drawnRecords": drawn,
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
    # A record with windows draws one record.
    with pytest.raises(ValidationError):
        plan_card_windows([_card("a", [(1, 100, False, None), (301, 400, False, None)], 10_000, False, drawn=2)])


def test_card_window_plan_counts_each_record_a_file_card_draws_whole() -> None:
    """A file card of two records drawn whole, no length known, before a windowed card in one row: the drawing draws
    both records first, so it numbers the windows 3 and 4; in rows of their own the windows are 1 and 2."""
    from gbdraw.web_support.linear_units import plan_card_windows

    windowed = _card("a", [(1, 100, False, None), (301, 400, False, None)], 10_000, False, row=1)
    shared = plan_card_windows([_card("m", [], None, False, row=1, drawn=2), windowed])
    assert [entry["number"] for entry in shared["windows"]] == [3, 4]
    apart = plan_card_windows([_card("m", [], None, False, drawn=2), {**windowed, "gridRow": None}])
    assert [entry["number"] for entry in apart["windows"]] == [1, 2]
    # A placed file card between the windows takes a column per record (OV-3072): the windows are 1 and 4.
    placed = plan_card_windows([
        _card("a", [(1, 100, False, 1), (301, 400, False, 3)], 10_000, False, row=1),
        _card("m", [], None, False, row=1, column=2, drawn=2),
    ])
    assert [entry["number"] for entry in placed["windows"]] == [1, 4]


def test_plan_record_places_numbers_a_placed_file_drawn_whole_and_moves_the_later_columns_of_its_row() -> None:
    """OV-3072: the columns a Source recipe's records table writes are the drawing's (``numbered_places``)."""
    import json

    from gbdraw.web_support.linear_units import plan_record_places, plan_record_places_json

    records = [
        {"gridRow": 1, "gridColumn": 1, "windowColumns": [], "drawnRecords": 3},
        {"gridRow": 1, "gridColumn": None, "windowColumns": [2, 3], "drawnRecords": 1},
        {"gridRow": 2, "gridColumn": 1, "windowColumns": [], "drawnRecords": 1},
        {"gridRow": 1, "gridColumn": 4, "windowColumns": [None], "drawnRecords": 1},
    ]
    assert plan_record_places(records) == [
        {"columns": [1, 2, 3]}, {"columns": [4, 5]}, {"columns": [1]}, {"columns": [6]},
    ]
    assert json.loads(plan_record_places_json(json.dumps(records[2:3]))) == [{"columns": [1]}]
    with pytest.raises(ValidationError):
        plan_record_places([{"gridRow": 1, "gridColumn": 1, "windowColumns": [2, 3], "drawnRecords": 2}])
