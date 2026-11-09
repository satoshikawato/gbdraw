"""Find in other files in the Web (PD-OI-123): the Worker helpers over the shown drawing's request.

The drawing is a Linear request with one record per File (a Linear File card),
keyed as the Web keys its cards. Searches run with the fake LOSAT of
``test_locate_command``; the plan the helpers give equals the plan ``gbdraw
locate`` writes for the same Files, and the new drawing is derived from the
source drawing as New drawing from a region derives one.
"""

from __future__ import annotations

import json
import re
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordCardinality, RecordInput
from gbdraw.api import render_request
from gbdraw.comparisons.locate import LocatePlan, LocateWindowChoice, choose_locate_windows, locate_file_counts
from gbdraw.comparisons.losat_runtime import run_losat_search
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request
from gbdraw.web_support.locate import LocateOptionsInput, LocateSearchInput, locate_web_json
from tests.test_locate_command import _fake_losat, _genomes, _interleaved, _locate_args, _run

OPTIONS: LocateOptionsInput = {"evalue": 1e-2, "bitscore": 50.0, "identity": 0.0, "alignmentLength": 0}


def _source(paths: dict[str, Path], *names: str) -> tuple[dict[str, Any], dict[str, str]]:
    # A card of a File of several records draws them all, as the Web's request writes it.
    request = LinearDiagramRequest(records=[
        RecordInput(GenBankInputSource(paths[name]), record_key=name, cardinality=RecordCardinality.ALL) for name in names
    ])
    encoded = encode_canonical_request(request)
    return encoded.payload, {resource.resource_id: str(resource.source_path) for resource in encoded.resources}


def _call(step: str, source: tuple[dict[str, Any], dict[str, str]], tmp_path: Path, payload: object) -> Any:
    """One step of the Worker operation, through its JSON boundary."""

    request, resources = source
    return json.loads(locate_web_json(
        json.dumps(request), json.dumps(resources), str(tmp_path / "out"), json.dumps({"step": step, step: payload})
    ))


def _windows(*spans: tuple[int, int]) -> list[dict[str, Any]]:
    return [{"start": start, "end": end, "reverseComplement": False} for start, end in spans]


def test_the_reference_counts_its_genes_or_names_the_rule_a_window_breaks(tmp_path: Path) -> None:
    source = _source(_genomes(tmp_path), "ref", "a")

    ready = _call("reference", source, tmp_path, {"recordKey": "ref", "windows": _windows((1501, 8100))})
    assert ready == {"status": "ready", "fileKey": "ref", "recordId": "REF", "genes": 6}
    refusals = {
        "WINDOW_WITHOUT_CDS": _windows((100, 900)),
        "WINDOWS_OVERLAP": _windows((1501, 4500), (4000, 8100)),
        "WINDOW_PAST_END": _windows((1501, 99_000)),
        "WINDOW_ACROSS_LINEAR_ORIGIN": _windows((8100, 1501)),
    }
    for reason, windows in refusals.items():
        reply = _call("reference", source, tmp_path, {"recordKey": "ref", "windows": windows})
        assert (reply["status"], reply["reason"]) == ("refused", reason)
    assert _call("reference", source, tmp_path, {"recordKey": "ref", "windows": _windows((1501, 4500), (4000, 8100))})[
        "start"
    ] == 4000


def test_the_files_list_every_record_of_every_file_of_the_drawing(tmp_path: Path) -> None:
    source = _source(_genomes(tmp_path), "ref", "b")
    request, resources = source
    reply = json.loads(locate_web_json(json.dumps(request), json.dumps(resources), str(tmp_path / "out"), json.dumps({"step": "files"})))
    assert reply == {"files": [
        {"fileKey": "ref", "records": [{"recordKey": "ref", "recordId": "REF", "length": 12_000, "circular": False}]},
        {"fileKey": "b", "records": [
            {"recordKey": "b:1", "recordId": "B1", "length": 5_500, "circular": False},
            {"recordKey": "b:2", "recordId": "B2", "length": 8_500, "circular": False},
        ]},
    ]}


def _search(paths: dict[str, Path], tmp_path: Path, targets: tuple[str, ...] = ("a", "b", "c")) -> tuple[Any, LocateSearchInput, Any]:
    losat = _fake_losat(tmp_path)
    source = _source(paths, "ref", *targets)
    search: LocateSearchInput = {
        "reference": {"recordKey": "ref", "windows": _windows((1501, 8100))},
        "targets": list(targets),
        "program": "tlosatn",
    }
    jobs = _call("search", source, tmp_path, search)
    reply = _call("plan", source, tmp_path, {
        "search": search,
        "outcomes": [
            {"fileKey": job["fileKey"], "text": run_losat_search("tlosatn", job["queryFasta"], job["subjectFasta"], losat_bin=losat), "error": None}
            for job in jobs["jobs"]
        ],
        "options": OPTIONS,
    })
    return source, search, reply


def _keyed(plan: dict[str, Any]) -> list[Any]:
    """Each File's windows and gene counts without the keys a surface gives its Files, records and proteins."""

    return [
        (
            [
                (item["record_id"], item["window"], item["selected"],
                 [{key: value for key, value in gene.items() if key != "protein_id"} for gene in item["genes"]])
                for item in file["windows"]
            ],
            {key: len(value) for key, value in file["coverage"].items()},
        )
        for file in plan["files"]
    ]


def test_the_plan_of_the_web_search_is_the_plan_the_command_writes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = _genomes(tmp_path)
    prefix = tmp_path / "cli"
    losat = _fake_losat(tmp_path)
    assert _run([*_locate_args(paths, losat, prefix, "a", "b", "c"), "--bitscore", "50", "--evalue", "0.01", "--identity", "0",
                 "--alignment_length", "0"], monkeypatch) == 0
    written = json.loads(Path(str(prefix) + ".plan.json").read_text(encoding="utf-8"))

    _source_, _search_, reply = _search(paths, tmp_path)

    assert [file["file_key"] for file in reply["plan"]["files"]] == ["a", "b", "c"]
    assert _keyed(reply["plan"]) == _keyed(written)
    # The drawing's keys: the reference record and each target record by the key the drawing gives it.
    assert reply["plan"]["reference"]["record_key"] == "ref"
    assert [w["record_key"] for w in reply["plan"]["files"][1]["windows"]] == ["b:1", "b:2"]
    assert [(pair["query"]["recordKey"], pair["subject"]["recordKey"]) for pair in reply["pairs"]] == [
        ("ref", "a"), ("a", "b:1"), ("a", "b:2"),
    ]
    # The name of each reference gene, by its protein ID: its locus tag here.
    assert [reply["geneNames"][gene["cds"]["protein_id"]] for gene in reply["plan"]["reference"]["genes"]] == [
        f"g{index}" for index in range(1, 7)
    ]
    # The CDS a window edge cuts is named too: the proposal names it.
    assert [reply["geneNames"][cds["protein_id"]] for cds in reply["plan"]["reference"]["cut_genes"]] == ["g7"]


def test_the_plan_reads_back_from_its_json(tmp_path: Path) -> None:
    plan = _search(_genomes(tmp_path), tmp_path)[2]["plan"]
    assert LocatePlan.from_json(plan).to_json() == plan
    with pytest.raises(Exception):
        LocatePlan.from_json({**plan, "files": [{**plan["files"][0], "windows": "none"}]})


def test_kept_windows_and_edits_change_only_those_windows_and_the_coverage(tmp_path: Path) -> None:
    plan = LocatePlan.from_json(_search(_genomes(tmp_path), tmp_path)[2]["plan"])
    b = plan.files[1]
    first, second = b.windows
    assert b.coverage.absent == ()

    # B2 left out: B's genes g4..g6 are no longer found, the plan says so, and nothing else changes.
    chosen = choose_locate_windows(plan, [
        *(LocateWindowChoice(file.file_key, w.window.column or 0, w.window.start, w.window.end, w.window.reverse_complement)
          for file in (plan.files[0],) for w in file.windows),
        LocateWindowChoice("b", first.window.column or 0, first.window.start, first.window.end, False),
    ])
    assert chosen.files[0] == plan.files[0]
    kept_b = chosen.files[1]
    assert [w.record_id for w in kept_b.windows] == ["B1"]
    assert [w.record_id for w in kept_b.candidates] == ["B2"] and kept_b.candidates[0].reason is None
    assert len(kept_b.coverage.absent) == 3
    assert ("GENES_ABSENT", (("count", 3),)) in [(w.code, w.args) for w in kept_b.warnings]
    assert chosen.files[2].windows == ()  # C had none; nothing chosen there either

    # An edit moves only that window: its span, orientation and contig ends follow; its genes stay as found.
    edited = choose_locate_windows(plan, [
        LocateWindowChoice("b", first.window.column or 0, 200, first.window.end, True),
        LocateWindowChoice("b", second.window.column or 0, second.window.start, second.window.end, False),
    ])
    window = edited.files[1].windows[0]
    assert (window.window.start, window.window.reverse_complement, window.touches_start) == (200, True, False)
    assert window.genes == first.genes
    assert edited.files[1].windows[1] == second

    with pytest.raises(Exception) as refused:
        choose_locate_windows(plan, [LocateWindowChoice("b", first.window.column or 0, 3000, 100, False)])
    assert (refused.value.reason, refused.value.place) == ("WINDOW_ACROSS_LINEAR_ORIGIN", ("b", first.window.column))


def test_the_new_drawing_draws_the_kept_windows_from_the_source_drawing(tmp_path: Path) -> None:
    paths = _genomes(tmp_path)
    source, _search_input, reply = _search(paths, tmp_path)
    plan = reply["plan"]
    kept = [
        {"fileKey": file["file_key"], "column": w["window"]["column"], "start": w["window"]["start"],
         "end": w["window"]["end"], "reverseComplement": w["window"]["reverse_complement"]}
        for file in plan["files"] for w in file["windows"]
    ]
    labels = {"ref": "ref.gb", "a": "a.gb", "b": "b.gb", "c": "c.gb"}
    derived = _call("drawing", source, tmp_path, {"plan": plan, "choices": kept, "labels": labels, "adaptSizes": True})

    assert derived["status"] == "derived"
    assert derived["pairs"] == reply["pairs"]
    drawing = derived["drawing"]
    assert (drawing["status"], drawing["mode"]) == ("derived", "linear")
    request = drawing["renderRequest"]
    # Row 1 the reference, then one row per File with a kept window, labelled with the File's name.
    assert [(record["presentation"]["gridRow"], record["presentation"].get("label")) for record in request["records"]] == [
        (1, "ref.gb"), (2, "a.gb"), (3, "b.gb"), (3, None),
    ]
    # The new drawing reads the source drawing's Files: each resource is one the source drawing has.
    assert set(drawing["resourceSources"].values()) <= set(source[1])
    decoded = decode_canonical_request(
        request, resource_paths={key: source[1][value] for key, value in drawing["resourceSources"].items()},
        output_directory=tmp_path / "out",
    )
    assert isinstance(decoded, LinearDiagramRequest)
    assert decoded.options.losat_search is not None and len(decoded.options.losat_search.pairs) == 3

    # Leaving B2 out leaves its record out of the drawing and out of the pairs.
    fewer = _call("drawing", source, tmp_path, {
        "plan": plan, "choices": kept[:-1], "labels": labels, "adaptSizes": True,
    })
    assert len(fewer["drawing"]["renderRequest"]["records"]) == 3
    assert len(fewer["pairs"]) == 2

    # An edit the plan cannot take names its window.
    refused = _call("drawing", source, tmp_path, {
        "plan": plan, "choices": [{**kept[2], "start": 3000, "end": 100}], "labels": labels, "adaptSizes": True,
    })
    assert refused == {
        "status": "refused", "reason": "WINDOW_ACROSS_LINEAR_ORIGIN", "start": 3000, "end": 100,
        "fileKey": "b", "column": kept[2]["column"],
    }


def test_a_span_picked_in_a_window_of_a_record_drawn_in_several_is_a_reference_of_that_record(tmp_path: Path) -> None:
    # The Region popup names the picked span by the drawn unit's geometry entry: its record's key and the span in
    # source coordinates. The dialog lists the record under the same key and takes the span as the reference.
    from gbdraw.api.request_render import plan_linear_request
    from gbdraw.api.requests import GenBankInputSource as Source
    from gbdraw.io.regions import RecordWindow

    mito = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    request = LinearDiagramRequest(records=(
        RecordInput(Source(mito), windows=(RecordWindow(100, 3000), RecordWindow(5000, 9000)), record_key="mito"),
        RecordInput(Source(mito), record_key="other"),
    ))
    geometry = plan_linear_request(request).build().drawing._gbdraw_track_slot_geometry["records"]
    units = [(entry["recordKey"], entry.get("windowIndex")) for entry in geometry]
    assert units[:2] == [("mito", 0), ("mito", 1)]
    encoded = encode_canonical_request(request)
    source = (encoded.payload, {resource.resource_id: str(resource.source_path) for resource in encoded.resources})

    files = _call("files", source, tmp_path, {})
    assert [record["recordKey"] for record in files["files"][0]["records"]] == ["mito"]
    ready = _call("reference", source, tmp_path, {"recordKey": units[1][0], "windows": _windows((5200, 8800))})
    assert (ready["status"], ready["fileKey"], ready["recordId"]) == ("ready", "mito", "NC_012920.1")


def _kept(plan: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"fileKey": file["file_key"], "column": w["window"]["column"], "start": w["window"]["start"],
         "end": w["window"]["end"], "reverseComplement": w["window"]["reverse_complement"]}
        for file in plan["files"] for w in file["windows"]
    ]


def test_the_proposal_counts_are_the_report_counts_and_number_the_windows_as_the_drawing_does(tmp_path: Path) -> None:
    # PD-OI-123: the circled number of a window is the one it has in the drawing; Python gives it, with the genes
    # found and the alignments the thresholds removed, to the dialog and to the report of `gbdraw locate`
    # (`locate_file_counts`). D's windows run D1, D2, D1 in the order of the reference (PD-OI-122: by column).
    paths = _interleaved(tmp_path)
    source, _search_input, reply = _search(paths, tmp_path, ("d", "c"))
    plan = reply["plan"]
    derived = _call("drawing", source, tmp_path, {
        "plan": plan, "choices": _kept(plan), "labels": {"ref": "ref.gb", "d": "d.gb", "c": "c.gb"}, "adaptSizes": True,
    })
    d, c = derived["counts"]
    assert [(w["record_id"], w["window"]["column"]) for w in plan["files"][0]["windows"]] == [("D1", 1), ("D2", 2), ("D1", 3)]
    assert d == {
        "fileKey": "d", "windowNumbers": [{"column": 1, "number": 1}, {"column": 2, "number": 2}, {"column": 3, "number": 3}],
        "genesFound": 6, "genesTotal": 6, "alignmentsRemoved": 0,
    }
    assert (c["fileKey"], c["windowNumbers"], c["genesFound"]) == ("c", [], 0)
    chosen = LocatePlan.from_json(derived["plan"])
    assert derived["counts"] == [
        {"fileKey": item.file_key, "windowNumbers": [{"column": a, "number": b} for a, b in item.window_numbers],
         "genesFound": item.genes_found, "genesTotal": item.genes_total, "alignmentsRemoved": item.alignments_removed}
        for item in locate_file_counts(chosen)
    ]

    # The new drawing draws D's row in that order: its window list reads D1, D2, D1.
    drawing = derived["drawing"]
    request = decode_canonical_request(
        drawing["renderRequest"], resource_paths={key: source[1][value] for key, value in drawing["resourceSources"].items()},
        output_directory=tmp_path / "out",
    )
    assert isinstance(request, LinearDiagramRequest)
    request = replace(request, options=replace(request.options, losat_search=None))
    svg = render_request(request).drawing.tostring()
    spans = [f"{w['record_id']}:{w['window']['start']:,}-{w['window']['end']:,}" for w in plan["files"][0]["windows"]]
    assert re.findall(r">(D[12]:[^<]*)<", svg) == spans


def _numbers(counts: dict[str, Any]) -> list[tuple[int, int]]:
    return [(item["column"], item["number"]) for item in counts["windowNumbers"]]


def _drawn_list(derived: dict[str, Any], source: tuple[dict[str, Any], dict[str, str]], tmp_path: Path) -> list[str]:
    """The window list of the new drawing's windowed row, as the SVG writes it (numbers 1.. in that order)."""

    drawing = derived["drawing"]
    request = decode_canonical_request(
        drawing["renderRequest"], resource_paths={key: source[1][value] for key, value in drawing["resourceSources"].items()},
        output_directory=tmp_path / "out",
    )
    assert isinstance(request, LinearDiagramRequest)
    request = replace(request, options=replace(request.options, losat_search=None))
    return re.findall(r">([DE][12]:[^<]*)<", render_request(request).drawing.tostring())


def test_an_unticked_first_window_leaves_the_others_numbered_as_the_drawing_numbers_them(tmp_path: Path) -> None:
    # The number of a window is its place in the drawn row, found by its column in the plan: with D1 ① left out,
    # D2 is drawn first and numbered 1, the second D1 window 2.
    paths = _interleaved(tmp_path)
    source, _search_input, reply = _search(paths, tmp_path, ("d",))
    plan = reply["plan"]
    labels = {"ref": "ref.gb", "d": "d.gb"}
    derived = _call("drawing", source, tmp_path, {
        "plan": plan, "choices": [item for item in _kept(plan) if (item["fileKey"], item["column"]) != ("d", 1)],
        "labels": labels, "adaptSizes": True,
    })
    [d] = derived["counts"]
    assert _numbers(d) == [(2, 1), (3, 2)]
    windows = {w["window"]["column"]: w for w in plan["files"][0]["windows"]}
    assert _drawn_list(derived, source, tmp_path) == [
        f"{windows[column]['record_id']}:{windows[column]['window']['start']:,}-{windows[column]['window']['end']:,}"
        for column, _number in _numbers(d)
    ]


def test_a_kept_window_after_one_not_selected_is_numbered_1_in_the_proposal(tmp_path: Path) -> None:
    from tests.test_locate_command import _upstream_copy

    paths = _upstream_copy(tmp_path)
    source, _search_input, reply = _search(paths, tmp_path, ("e",))
    plan = reply["plan"]
    derived = _call("drawing", source, tmp_path, {
        "plan": plan, "choices": _kept(plan), "labels": {"ref": "ref.gb", "e": "e.gb"}, "adaptSizes": True,
    })
    [e] = derived["counts"]
    assert _numbers(e) == [(2, 1)]


def _split(paths: dict[str, Path], name: str, *record_ids: str) -> list[RecordInput]:
    """A card the request splits into one record per record of its File (as the Web writes a card whose records
    differ in a per-record setting): records that share one source."""

    from gbdraw.io.record_select import parse_record_selector

    return [
        RecordInput(GenBankInputSource(paths[name]), selector=parse_record_selector(record_id), record_key=f"{name}:{index}")
        for index, record_id in enumerate(record_ids, start=1)
    ]


def _encoded(records: list[RecordInput]) -> tuple[dict[str, Any], dict[str, str]]:
    encoded = encode_canonical_request(LinearDiagramRequest(records=records))
    return encoded.payload, {resource.resource_id: str(resource.source_path) for resource in encoded.resources}


def test_records_of_one_source_are_one_file_searched_once(tmp_path: Path) -> None:
    # Files are the drawing's sources: D's two request records are one File with both drawn record keys, searched
    # once and listed once.
    paths = _interleaved(tmp_path)
    source = _encoded([
        RecordInput(GenBankInputSource(paths["ref"]), record_key="ref", cardinality=RecordCardinality.ALL),
        *_split(paths, "d", "D1", "D2"),
    ])
    files = _call("files", source, tmp_path, {})
    assert [(file["fileKey"], [record["recordKey"] for record in file["records"]]) for file in files["files"]] == [
        ("ref", ["ref"]), ("d:1", ["d:1", "d:2"]),
    ]
    search: LocateSearchInput = {
        "reference": {"recordKey": "ref", "windows": _windows((1501, 8100))}, "targets": ["d:1"], "program": "tlosatn",
    }
    jobs = _call("search", source, tmp_path, search)["jobs"]
    assert [job["fileKey"] for job in jobs] == ["d:1"]
    with pytest.raises(Exception, match="other Files of the drawing"):
        _call("search", source, tmp_path, {**search, "targets": ["d:1", "d:2"]})


def test_a_record_that_shares_the_reference_source_is_not_another_file(tmp_path: Path) -> None:
    # The reference's File is not searched (PD-OI-123): a record of the same source is part of it, not a target.
    paths = _interleaved(tmp_path)
    source = _encoded([*_split(paths, "d", "D1", "D2"), RecordInput(GenBankInputSource(paths["ref"]), record_key="ref")])
    files = _call("files", source, tmp_path, {})
    assert [file["fileKey"] for file in files["files"]] == ["d:1", "ref"]
    ready = _call("reference", source, tmp_path, {"recordKey": "d:2", "windows": _windows((2001, 3300))})
    assert (ready["status"], ready["fileKey"], ready["recordId"]) == ("ready", "d:1", "D2")
    search: LocateSearchInput = {
        "reference": {"recordKey": "d:2", "windows": _windows((2001, 3300))}, "targets": ["d:1"], "program": "tlosatn",
    }
    with pytest.raises(Exception, match="other Files of the drawing"):
        _call("search", source, tmp_path, search)
    assert [job["fileKey"] for job in _call("search", source, tmp_path, {**search, "targets": ["ref"]})["jobs"]] == ["ref"]
