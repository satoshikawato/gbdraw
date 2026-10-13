"""The label fit report the Web reads, and the rows it computes on demand (OD-F14)."""

from __future__ import annotations

from pathlib import Path
import warnings

from Bio import SeqIO
import pytest

from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.exceptions import LabelsDoNotFitError, ValidationError
from gbdraw.labels.fit import FONT_SIZE_PATHS, ROW_KINDS, Y_OFFSET_PATH
from gbdraw.labels.fit import LabelFit, LabelFitRecord, LabelFitReport, LabelFitRow
from gbdraw.session import build_session_document
from gbdraw.web_support.label_fit import (
    LabelFitRowSearches,
    label_fit_error_payload,
    label_fit_report_payload,
)
from gbdraw.web_support.request_render import render_canonical_web_request

TEST_INPUTS = Path(__file__).parent / "test_inputs"


def _report(*, drawn: bool) -> LabelFitReport:
    return LabelFitReport(
        records=(LabelFitRecord("rec-1", "circular_outer", 111, 14.0, LabelFit(20, 21)),),
        rows=(
            LabelFitRow("outward_vertical", "found", {Y_OFFSET_PATH: 1.15}, (1708.0, 1130.0)),
            LabelFitRow("smaller_font", "found", {path: 12.0 for path in FONT_SIZE_PATHS}, (1615.0, 988.0)),
            LabelFitRow("tighter_rows", "not_found", reason="Does not fit down to 0.1 px label spacing"),
            LabelFitRow("combination", "not_computed"),
        ),
        drawn=drawn,
    )


def test_report_payload_carries_counts_rows_and_python_wording() -> None:
    payload = label_fit_report_payload(_report(drawn=False))

    assert payload["schema"] == 1 and payload["drawn"] is False
    assert (payload["overlapPairs"], payload["crossingPairs"]) == (20, 21)
    assert payload["records"] == [{
        "recordId": "rec-1", "labelKind": "circular_outer", "labelKindTitle": "outer labels",
        "labelCount": 111, "fontSize": 14.0, "overlapPairs": 20, "crossingPairs": 21,
    }]
    vertical, font, tighter, combination = payload["rows"]
    assert vertical == {
        "kind": "outward_vertical", "title": "outward, vertical only", "status": "found",
        "values": {Y_OFFSET_PATH: 1.15},
        "settings": [{"name": "Outer Y Offset", "unit": "", "value": 1.15}],
        "canvasSize": {"width": 1708.0, "height": 1130.0}, "reason": None,
    }
    # One control sets both font size classes: one setting, both config paths.
    assert font["values"] == {path: 12.0 for path in FONT_SIZE_PATHS}
    assert font["settings"] == [{"name": "Label Font Size", "unit": "pt", "value": 12.0}]
    assert tighter["status"] == "not_found" and tighter["values"] is None and tighter["settings"] == []
    assert tighter["reason"] == "Does not fit down to 0.1 px label spacing"
    assert combination["status"] == "not_computed" and combination["canvasSize"] is None


def test_a_row_with_an_unknown_config_path_is_an_error() -> None:
    report = LabelFitReport(
        records=(),
        rows=(LabelFitRow("outward_vertical", "found", {"labels.unknown": 1.0}, (1.0, 1.0)),),
        drawn=True,
    )
    with pytest.raises(ValidationError, match="unknown config path"):
        label_fit_report_payload(report)


def test_the_error_payload_is_found_through_the_exception_chain() -> None:
    report = _report(drawn=False)
    try:
        try:
            raise LabelsDoNotFitError("labels overlap", report=report)
        except LabelsDoNotFitError as error:
            raise ValidationError("wrapped") from error
    except ValidationError as wrapped:
        assert label_fit_error_payload(wrapped) == label_fit_report_payload(report)
    assert label_fit_error_payload(ValueError("other")) is None


def _crossing_request() -> dict:
    record = SeqIO.read(str(TEST_INPUTS / "SARS-CoV-2.gbk"), "genbank")
    return build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(config_overrides={
                "labels.circular.scope": "outer",
                "labels.font_size.short": 14.0,
                "labels.font_size.long": 14.0,
            }),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    ).to_dict()


def test_web_render_draws_crossing_leaders_with_counts_only_and_computes_rows_later(tmp_path) -> None:
    document = _crossing_request()
    resources = tmp_path / "resources"
    resources.mkdir()
    from gbdraw.session_io import materialize_embedded_file

    resource_paths = {
        resource_id: materialize_embedded_file(entry, temp_dir=resources, role=resource_id, prefix_role=False)
        for resource_id, entry in document["resources"].items()
    }
    searches = LabelFitRowSearches()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        response = render_canonical_web_request(
            document["renderRequest"],
            resource_paths=resource_paths,
            output_directory=tmp_path / "out",
            _label_fit_searches=searches,
        )

    (entry,) = response["metadata"]["labelFitReports"]
    assert entry["resultIndex"] == 0 and entry["resultName"] == "fit.svg"
    report = entry["report"]
    assert report["drawn"] is True and report["overlapPairs"] == 0 and report["crossingPairs"] > 0
    # OD-F14: Generate searched no row.
    assert [row["kind"] for row in report["rows"]] == list(ROW_KINDS)
    assert {row["status"] for row in report["rows"]} == {"not_computed"}
    token = entry["rowToken"]
    assert isinstance(token, str)

    row = searches.row(token, "outward_vertical")
    assert row["kind"] == "outward_vertical" and row["status"] in {"found", "not_found"}
    if row["status"] == "found":
        assert set(row["values"] or {}) == {Y_OFFSET_PATH} and row["canvasSize"] is not None

    # A newer render of the mode releases the figure: its rows answer "not computed".
    searches.retain("circular", [])
    stale = searches.row(token, "outward_vertical")
    assert stale["status"] == "not_computed" and stale["reason"]
    with pytest.raises(ValidationError):
        searches.row(token, "sideways")


def test_a_render_without_searches_keeps_no_token(tmp_path) -> None:
    document = _crossing_request()
    from gbdraw.web_support.request_render import render_embedded_canonical_web_request

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        response = render_embedded_canonical_web_request(
            document["renderRequest"], resources=document["resources"], workspace=tmp_path / "ws"
        )
    (entry,) = response["metadata"]["labelFitReports"]
    assert entry["rowToken"] is None
