"""Find in other files in the Web (PD-OI-123): the Worker helpers over the shown drawing's request.

The drawing is a Linear request with one record per File (a Linear File card),
keyed as the Web keys its cards. Searches run with the fake LOSAT of
``test_locate_command``; the plan the helpers give equals the plan ``gbdraw
locate`` writes for the same Files, and the new drawing is derived from the
source drawing as New drawing from a region derives one.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordCardinality, RecordInput
from gbdraw.comparisons.locate import LocatePlan, LocateWindowChoice, choose_locate_windows
from gbdraw.comparisons.losat_runtime import run_losat_search
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request
from gbdraw.web_support.locate import (
    LocateOptionsInput,
    LocateSearchInput,
    locate_web_drawing_json,
    locate_web_plan_json,
    locate_web_reference_json,
    locate_web_search_plan_json,
)
from tests.test_locate_command import _fake_losat, _genomes, _locate_args, _run

OPTIONS: LocateOptionsInput = {
    "evalue": 1e-2, "bitscore": 50.0, "identity": 0.0, "alignmentLength": 0, "flank": 500, "mergeDistance": 10_000,
    "maxWindows": 4,
}


def _source(paths: dict[str, Path], *names: str) -> tuple[dict[str, Any], dict[str, str]]:
    # A card of a File of several records draws them all, as the Web's request writes it.
    request = LinearDiagramRequest(records=[
        RecordInput(GenBankInputSource(paths[name]), record_key=name, cardinality=RecordCardinality.ALL) for name in names
    ])
    encoded = encode_canonical_request(request)
    return encoded.payload, {resource.resource_id: str(resource.source_path) for resource in encoded.resources}


def _call(helper: Any, source: tuple[dict[str, Any], dict[str, str]], tmp_path: Path, payload: object) -> Any:
    request, resources = source
    return json.loads(helper(json.dumps(request), json.dumps(resources), str(tmp_path / "out"), json.dumps(payload)))


def _windows(*spans: tuple[int, int]) -> list[dict[str, Any]]:
    return [{"start": start, "end": end, "reverseComplement": False} for start, end in spans]


def test_the_reference_counts_its_genes_or_names_the_rule_a_window_breaks(tmp_path: Path) -> None:
    source = _source(_genomes(tmp_path), "ref", "a")

    ready = _call(locate_web_reference_json, source, tmp_path, {"recordKey": "ref", "windows": _windows((1501, 8100))})
    assert ready == {"status": "ready", "fileKey": "ref", "recordId": "REF", "genes": 6, "cutGenes": ["g7"]}
    refusals = {
        "WINDOW_WITHOUT_CDS": _windows((100, 900)),
        "WINDOWS_OVERLAP": _windows((1501, 4500), (4000, 8100)),
        "WINDOW_PAST_END": _windows((1501, 99_000)),
        "WINDOW_ACROSS_LINEAR_ORIGIN": _windows((8100, 1501)),
    }
    for reason, windows in refusals.items():
        reply = _call(locate_web_reference_json, source, tmp_path, {"recordKey": "ref", "windows": windows})
        assert (reply["status"], reply["reason"]) == ("refused", reason)
    assert _call(locate_web_reference_json, source, tmp_path, {"recordKey": "ref", "windows": _windows((1501, 4500), (4000, 8100))})[
        "start"
    ] == 4000


def _search(paths: dict[str, Path], tmp_path: Path) -> tuple[Any, LocateSearchInput, Any]:
    losat = _fake_losat(tmp_path)
    source = _source(paths, "ref", "a", "b", "c")
    search: LocateSearchInput = {
        "reference": {"recordKey": "ref", "windows": _windows((1501, 8100))},
        "targets": ["a", "b", "c"],
        "program": "tlosatn",
    }
    jobs = _call(locate_web_search_plan_json, source, tmp_path, search)
    reply = _call(locate_web_plan_json, source, tmp_path, {
        **search,
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
    assert _run([*_locate_args(paths, losat, prefix, "a", "b", "c"), "--flank", "500", "--max_windows", "4"], monkeypatch) == 0
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
    derived = _call(locate_web_drawing_json, source, tmp_path, {"plan": plan, "choices": kept, "labels": labels, "adaptSizes": True})

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
    fewer = _call(locate_web_drawing_json, source, tmp_path, {
        "plan": plan, "choices": kept[:-1], "labels": labels, "adaptSizes": True,
    })
    assert len(fewer["drawing"]["renderRequest"]["records"]) == 3
    assert len(fewer["pairs"]) == 2

    # An edit the plan cannot take names its window.
    refused = _call(locate_web_drawing_json, source, tmp_path, {
        "plan": plan, "choices": [{**kept[2], "start": 3000, "end": 100}], "labels": labels, "adaptSizes": True,
    })
    assert refused == {
        "status": "refused", "reason": "WINDOW_ACROSS_LINEAR_ORIGIN", "start": 3000, "end": 100,
        "fileKey": "b", "column": kept[2]["column"],
    }
