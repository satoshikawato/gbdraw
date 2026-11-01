"""The Worker helper behind the Web's "Orient records from comparisons" button.

The page sends the request of its last Generate; the helper plans it with
``auto_orient`` through the planner the CLI's ``--auto_orient`` runs, so both
surfaces decide with one owner and draw the same orientation.
"""
from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw
from gbdraw.api import (
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    plan_linear_request,
    resolve_request,
)
from gbdraw.api.record_planning import record_is_reverse_complemented
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.auto_orientation import decide_auto_orientation

_EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
_QUERY = _EXAMPLES / "MjeNMV.gb"
_SUBJECT = _EXAMPLES / "MelaMJNV.gb"
_TABLE = _EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out"


def _reverse_hits_table(directory: Path) -> Path:
    """The tblastx table with the subject coordinates mirrored: MelaMJNV's hits are mostly reverse."""

    length = len(SeqIO.read(_SUBJECT, "genbank"))
    table = directory / "MjeNMV.MelaMJNV.mirrored.tsv"
    with _TABLE.open(encoding="utf-8") as source, table.open("w", encoding="utf-8", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        for row in csv.reader(source, delimiter="\t"):
            if row and not row[0].startswith("#"):
                row[8], row[9] = str(length + 1 - int(row[8])), str(length + 1 - int(row[9]))
                writer.writerow(row)
    return table


def _request(table: Path, *, subject_reversed: bool = False) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(_QUERY), record_key="card-query"),
            RecordInput(
                GenBankInputSource(_SUBJECT),
                record_key="card-subject",
                presentation=RecordPresentation(reverse_complement=subject_reversed),
            ),
        ),
        options=LinearDiagramOptions(blast_files=(str(table),)),
    )


def _committed(request: LinearDiagramRequest, directory: Path) -> tuple[dict[str, object], dict[str, str]]:
    """The request as the Web's last Generate commits it: planned, encoded, its resources staged."""

    encoded = encode_canonical_request(resolve_request(request))
    staged = directory / "resources"
    staged.mkdir(parents=True, exist_ok=True)
    paths = {}
    for index, resource in enumerate(encoded.resources, start=1):
        target = staged / f"{index:04d}.bin"
        target.write_bytes(
            resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes()
        )
        paths[resource.resource_id] = str(target)
    return encoded.payload, paths


def _planner_flags(request: LinearDiagramRequest) -> list[bool]:
    """The orientation the CLI planner draws with ``--auto_orient``."""

    plan = plan_linear_request(replace(request, options=replace(request.options, auto_orient=True)))
    return [record_is_reverse_complemented(record) for record in plan.records]


@pytest.mark.linear
def test_the_helper_decides_what_the_cli_planner_draws(tmp_path: Path) -> None:
    request = _request(_reverse_hits_table(tmp_path))
    payload, paths = _committed(request, tmp_path)
    reply = decide_auto_orientation(payload, paths, str(tmp_path / "out"))
    assert reply == {"records": [
        {"recordKey": "card-query", "reverseComplement": False},
        {"recordKey": "card-subject", "reverseComplement": True},
    ]}
    assert [record["reverseComplement"] for record in reply["records"]] == _planner_flags(request)


@pytest.mark.linear
def test_the_helper_decides_every_record_from_the_comparisons(tmp_path: Path) -> None:
    forward = _request(_TABLE)
    payload, paths = _committed(forward, tmp_path / "forward")
    assert [record["reverseComplement"] for record in decide_auto_orientation(
        payload, paths, str(tmp_path / "out"))["records"]] == [False, False] == _planner_flags(forward)
    # A checked Reverse complement the comparisons contradict is decided again:
    # the helper's answer is the planner's on the request with the flags cleared.
    explicit = _request(_TABLE, subject_reversed=True)
    payload, paths = _committed(explicit, tmp_path / "explicit")
    assert [record["reverseComplement"] for record in decide_auto_orientation(
        payload, paths, str(tmp_path / "out"))["records"]] == [False, False] == _planner_flags(forward)
    # On the command line an explicit reverse still wins.
    assert _planner_flags(explicit) == [False, True]
    # A record reversed by the rule and drawn so stays reversed while its hits are mostly reverse.
    mirrored = _reverse_hits_table(tmp_path)
    payload, paths = _committed(_request(mirrored, subject_reversed=True), tmp_path / "drawn-reversed")
    assert [record["reverseComplement"] for record in decide_auto_orientation(
        payload, paths, str(tmp_path / "out"))["records"]] == [False, True] == _planner_flags(_request(mirrored))


@pytest.mark.linear
def test_the_worker_json_boundary_names_its_operation(tmp_path: Path) -> None:
    from gbdraw.web_support.auto_orientation import decide_auto_orientation_json

    payload, paths = _committed(_request(_reverse_hits_table(tmp_path)), tmp_path)
    reply = json.loads(decide_auto_orientation_json(json.dumps(payload), json.dumps(paths), str(tmp_path / "out")))
    assert [record["reverseComplement"] for record in reply["records"]] == [False, True]
    helpers_js = Path(gbdraw.__file__).parent / "web" / "js" / "app" / "python-helpers.js"
    namespace: dict[str, object] = {}
    exec(helpers_js.read_text(encoding="utf-8").split("`", 1)[1].rsplit("`", 1)[0], namespace)
    call = namespace["call_web_json_helper"]
    assert callable(call)
    # The Worker passes the staged identities after the helper's arguments; a failure names the operation.
    error = json.loads(call(
        "decide_auto_orientation_json", json.dumps(payload), json.dumps(paths), str(tmp_path / "out"), "{}"
    ))["error"]
    assert error["operation"] == "decideAutoOrientation"


def test_a_circular_request_is_refused_with_a_diagnostic(tmp_path: Path) -> None:
    from gbdraw.api import CircularDiagramRequest

    encoded = encode_canonical_request(resolve_request(CircularDiagramRequest(
        records=(RecordInput(GenBankInputSource(_QUERY)),),
    )))
    staged = {}
    for index, resource in enumerate(encoded.resources, start=1):
        target = tmp_path / f"{index:04d}.bin"
        target.write_bytes(
            resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes()
        )
        staged[resource.resource_id] = str(target)
    with pytest.raises(ValidationError, match="requires a Linear request") as raised:
        decide_auto_orientation(encoded.payload, staged, str(tmp_path / "out"))
    assert raised.value.diagnostic is not None
