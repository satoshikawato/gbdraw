"""The Web's Align at feature asks Python for the starts and translations ``--anchor_feature`` draws with.

The Worker helper runs the planner and drawing of ``gbdraw linear
--anchor_feature``, so the starts and translations it returns are the ones a
CLI Session stores, and a request with them draws the CLI figure. Each
record's status and the reference are the planner's.
"""
from __future__ import annotations

import json
import re
from dataclasses import replace
from pathlib import Path

import pytest

from gbdraw.api import (
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    plan_linear_request,
    resolve_request,
)
from gbdraw.api.options import LinearMultiRecordOptions, LinearRecordTranslation
from gbdraw.api.request_render import build_request_diagram
from gbdraw.api.requests import RecordDisplayOptions, RecordPresentation
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import encode_canonical_request, encode_record_translations
from gbdraw.web_support.error_adapter import serialize_web_error
from gbdraw.web_support.feature_anchor import resolve_feature_anchor, resolve_feature_anchor_json

_DATA = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data"
_HUMAN = _DATA / "human-mitochondrion" / "HmmtDNA.gbk"
_ZEBRAFISH = _DATA / "metazoan-mitochondria-four" / "NC_002333.2.gb"


def _committed(*, linear: tuple[int, ...] = (), starts: tuple[int | None, ...] = (None, None),
               **options: object) -> LinearDiagramRequest:
    """A request as the Web commits it: planned, with record keys. ``linear`` records are drawn linear."""

    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(GenBankInputSource(path), display=(
                RecordDisplayOptions(is_circular=False) if index in linear
                else RecordDisplayOptions(start_coordinate=start)
            ))
            for index, (path, start) in enumerate(zip((_HUMAN, _ZEBRAFISH), starts))
        ),
        options=LinearDiagramOptions(**options),
    )
    return plan_linear_request(request).request


def _cli_resolved(request: LinearDiagramRequest, selector: str) -> LinearDiagramRequest:
    resolved = resolve_request(replace(request, options=replace(request.options, anchor_feature=selector)))
    assert isinstance(resolved, LinearDiagramRequest) and resolved.layout is not None
    return resolved


def _cli_translations(request: LinearDiagramRequest, selector: str) -> list[dict[str, object]]:
    layout = _cli_resolved(request, selector).layout
    assert layout is not None
    return list(encode_record_translations(layout.record_translations))


def _written(request: LinearDiagramRequest, reply: dict) -> LinearDiagramRequest:
    """``request`` with the reply's starts and translations written, as the Web's Apply writes them."""

    return replace(
        request,
        records=tuple(
            replace(record, display=replace(record.display, start_coordinate=item["startCoordinate"]))
            if item["status"] == "rotated" else record
            for record, item in zip(request.records, reply["records"], strict=True)
        ),
        layout=LinearMultiRecordOptions(record_translations=tuple(
            LinearRecordTranslation(record_key=item["recordKey"], x=item["x"], y=item["y"])
            for item in reply["recordTranslations"]
        )),
    )


@pytest.mark.linear
@pytest.mark.parametrize("config_overrides", [None, {"canvas.linear.align_center": True}])
@pytest.mark.parametrize("linear", [(), (0,), (0, 1)])
def test_the_helper_returns_the_starts_and_translations_the_cli_stores(
    config_overrides: dict[str, object] | None, linear: tuple[int, ...],
) -> None:
    request = _committed(config_overrides=config_overrides, linear=linear)
    reply = resolve_feature_anchor(request, "gene=COX1")
    assert reply["recordTranslations"] == _cli_translations(request, "gene=COX1")
    expected = {
        (): [("NC_012920.1", "rotated", 5904), ("NC_002333.2", "rotated", 6425)],
        (0,): [("NC_012920.1", "aligned", None), ("NC_002333.2", "rotated", 6425)],
        (0, 1): [("NC_012920.1", "reference", None), ("NC_002333.2", "aligned", None)],
    }[linear]
    assert [(record["recordId"], record["status"], record["startCoordinate"]) for record in reply["records"]] == expected
    assert [record["startCoordinate"] for record in reply["records"]] == [
        record.display.start_coordinate for record in _cli_resolved(request, "gene=COX1").records
    ]
    human, zebrafish = reply["records"]
    assert human["beforeX"] == zebrafish["beforeX"] == 0.0
    if linear == (0, 1):
        assert human["afterX"] == 0.0 and zebrafish["afterX"] < 0.0
    elif linear == (0,):
        assert human["afterX"] < 0.0 and zebrafish["afterX"] == 0.0
    elif config_overrides is None:
        assert human["afterX"] == zebrafish["afterX"] == 0.0
    # The written starts and translations, without the option, draw the anchored figure.
    anchored = replace(request, options=replace(request.options, anchor_feature="gene=COX1"))
    assert (
        build_request_diagram(_written(request, reply)).drawing.tostring()
        == build_request_diagram(anchored).drawing.tostring()
    )


@pytest.mark.linear
def test_the_helper_starts_from_the_translations_and_starts_the_request_holds() -> None:
    request = _committed(linear=(0, 1))
    keys = [str(record.record_key) for record in request.records]
    held = replace(request, layout=LinearMultiRecordOptions(record_translations=(
        LinearRecordTranslation(record_key=keys[0], x=30.0, y=4.0),
        LinearRecordTranslation(record_key=keys[1], x=-500.0, y=6.0),
    )))
    reply = resolve_feature_anchor(held, "gene=COX1")
    assert reply["recordTranslations"] == _cli_translations(held, "gene=COX1")
    assert [record["beforeX"] for record in reply["records"]] == [30.0, -500.0]
    assert reply["recordTranslations"][0] == {"recordKey": keys[0], "x": 30.0, "y": 4.0}
    # A rotated record's earlier start and sideways place are replaced: Apply again rotates it.
    rotated = _committed(starts=(100, 7000))
    rotated = replace(rotated, layout=LinearMultiRecordOptions(record_translations=(
        LinearRecordTranslation(record_key=keys[0], x=30.0, y=4.0),
        LinearRecordTranslation(record_key=keys[1], x=-500.0, y=6.0),
    )))
    again = resolve_feature_anchor(rotated, "gene=COX1")
    assert [(record["startCoordinate"], record["afterX"]) for record in again["records"]] == [(5904, 0.0), (6425, 0.0)]
    assert again["recordTranslations"] == _cli_translations(rotated, "gene=COX1")


@pytest.mark.linear
def test_a_record_without_the_feature_is_listed_as_unmatched() -> None:
    request = _committed(starts=(None, 100))
    reply = resolve_feature_anchor(request, "gene=RNR1")
    # Only the human record names its 12S rRNA gene RNR1 (648..1601); the zebrafish record keeps its start.
    assert [(record["status"], record["startCoordinate"]) for record in reply["records"]] == [
        ("rotated", 648), ("unmatched", 100),
    ]
    assert [record["afterX"] for record in reply["records"]] == [0.0, 0.0]


def _web_error(error: ValidationError) -> dict[str, object]:
    return serialize_web_error(error, operation="resolveFeatureAnchor", stage="helper")


@pytest.mark.linear
@pytest.mark.parametrize(
    ("selector", "change", "context"),
    [
        ("gene=nonexistent", None, {"field": "anchor_feature", "reason": "SELECTOR_NO_MATCH"}),
        ("type=tRNA", None, {"field": "anchor_feature", "reason": "SELECTOR_AMBIGUOUS"}),
        ("gene=COX1", "shared-row", {"field": "anchor_feature", "reason": "ANCHOR_SHARED_ROW"}),
        ("gene=COX1", "wrap", {"field": "wrap_length", "reason": "WRAP_ANCHOR_FEATURE"}),
        ("", None, {"field": "anchor_feature", "reason": "REQUIRED"}),
    ],
)
def test_the_web_receives_each_refusal_as_a_field_and_a_reason(
    selector: str, change: str | None, context: dict[str, object],
) -> None:
    # Under a wrap, the cut human record drawn linear would move to the rotated zebrafish record's left end.
    request = (
        _committed(config_overrides={"canvas.linear.wrap_length": 5000}, linear=(0,))
        if change == "wrap" else _committed()
    )
    if change == "shared-row":
        request = replace(
            request,
            records=tuple(replace(record, presentation=RecordPresentation(grid_row=1)) for record in request.records),
            layout=LinearMultiRecordOptions(),
        )
    with pytest.raises(ValidationError) as caught:
        resolve_feature_anchor(request, selector)
    payload = _web_error(caught.value)
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"] == context
    assert payload["operation"] == "resolveFeatureAnchor"


@pytest.mark.linear
def test_the_json_boundary_reads_the_staged_request(tmp_path: Path) -> None:
    request = _committed()
    encoded = encode_canonical_request(request)
    paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.resource_id
        if resource.source_path is not None:
            target.write_bytes(Path(resource.source_path).read_bytes())
        else:
            assert resource.content is not None
            target.write_bytes(resource.content)
        paths[resource.resource_id] = str(target)
    reply = json.loads(resolve_feature_anchor_json(
        json.dumps(encoded.payload), json.dumps("gene=COX1"), json.dumps(paths), str(tmp_path / "output"),
    ))
    assert reply["recordTranslations"] == _cli_translations(request, "gene=COX1")
    assert reply["selector"] == "gene=COX1"
    assert re.fullmatch(r"record-\d+", reply["records"][1]["recordKey"])
    assert [(record["status"], record["startCoordinate"]) for record in reply["records"]] == [
        ("rotated", 5904), ("rotated", 6425),
    ]


@pytest.mark.linear
def test_records_the_anchor_only_rotates_wrap_in_the_helper() -> None:
    request = _committed(config_overrides={"canvas.linear.wrap_length": 5000})
    reply = resolve_feature_anchor(request, "gene=COX1")
    assert [record["afterX"] for record in reply["records"]] == [0.0, 0.0]
    assert "data-gbdraw-wrap-row" in build_request_diagram(_written(request, reply)).drawing.tostring()
