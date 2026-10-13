"""The label fit report the Web reads, and the rows it computes on demand (OD-F14)."""

from __future__ import annotations

from pathlib import Path
from typing import Any
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
from gbdraw.labels.circular_fit import ROW_KINDS
from gbdraw.labels.fit import (
    FONT_SIZE_PATHS,
    PLACEMENT_PATH,
    SCOPE_PATH,
    SEPARATE_FEATURES_PATH,
    SEPARATE_FEATURES_VALUES,
    Y_OFFSET_PATH,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
)
from gbdraw.session import build_session_document
from gbdraw.web_support.label_fit import (
    LabelFitRowSearches,
    label_fit_failure,
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
    payload = label_fit_report_payload(_report(drawn=False), "circular")

    assert payload["schema"] == 1 and payload["drawn"] is False
    assert (payload["overlapPairs"], payload["crossingPairs"]) == (20, 21)
    assert payload["records"] == [{
        "recordId": "rec-1", "labelKind": "circular_outer", "labelKindTitle": "outer labels",
        "labelCount": 111, "fontSize": 14.0, "overlapPairs": 20, "crossingPairs": 21,
        "definitionOverlaps": 0, "refused": False,
        "counts": "20 label pairs overlap, 21 leader pairs cross",
    }]
    vertical, font, tighter, combination = payload["rows"]
    assert vertical == {
        "kind": "outward_vertical", "title": "outward, vertical only", "status": "found",
        "values": {Y_OFFSET_PATH: 1.15},
        "settings": [{"name": "Outer Y Offset", "unit": "", "value": 1.15}],
        "canvasSize": {"width": 1708.0, "height": 1130.0}, "reason": None, "note": None,
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
        label_fit_report_payload(report, "circular")


def test_word_rows_show_the_web_option_and_a_refused_record_says_why() -> None:
    report = LabelFitReport(
        records=(LabelFitRecord("rec-2", "circular_inner", 145, 8.0, LabelFit(0, 0), refused=True),),
        rows=(
            LabelFitRow("labels_out", "found", {SCOPE_PATH: "outer"}, (900.0, 900.0)),
            LabelFitRow("radial", "found", {PLACEMENT_PATH: "radial"}, (950.0, 950.0)),
        ),
        drawn=False,
    )
    payload = label_fit_report_payload(report, "circular")
    (record,) = payload["records"]
    # OV-699: 0 overlaps and 0 crossings, and the reason instead of the counts.
    assert record["refused"] is True and record["labelKindTitle"] == "radial inner labels"
    assert record["counts"] == "they need more room than the circular tracks leave"
    outside, radial = payload["rows"]
    assert outside["values"] == {SCOPE_PATH: "outer"}
    assert outside["settings"] == [{"name": "Label Mode", "unit": "", "value": "Out"}]
    assert radial["settings"] == [{"name": "Circular Label Placement", "unit": "", "value": "Radial"}]


@pytest.mark.parametrize(
    ("mode", "kind", "name"),
    [("circular", "circular_embedded", "Resolve Overlaps"), ("linear", "linear", "Resolve Feature Overlaps")],
)
def test_the_separate_features_row_names_the_resolve_overlaps_control_of_the_mode(mode, kind, name) -> None:
    report = LabelFitReport(
        records=(LabelFitRecord("rec-3", kind, 3, 14.0, LabelFit(1, 0)),),
        rows=(LabelFitRow("separate_features", "found", SEPARATE_FEATURES_VALUES, (700.0, 700.0)),),
        drawn=False,
    )
    (row,) = label_fit_report_payload(report, mode)["rows"]
    assert row["title"] == "separate overlapping features"
    # The switch stays a boolean for the request mapping; the dialog shows the option word.
    assert row["values"] == {SEPARATE_FEATURES_PATH: True}
    assert row["settings"] == [{"name": name, "unit": "", "value": "On"}]
    off = LabelFitReport(
        records=(), rows=(LabelFitRow("separate_features", "found", {SEPARATE_FEATURES_PATH: False}, (1.0, 1.0)),),
        drawn=True,
    )
    with pytest.raises(ValidationError, match="unknown config path or value"):
        label_fit_report_payload(off, mode)


def test_every_control_a_row_names_has_that_accessible_name() -> None:
    """The dialog names a row's settings by their Web controls; each control carries that name (R12, OV-721)."""

    import re

    from gbdraw.web_support import label_fit

    html = (Path(__file__).parents[1] / "gbdraw" / "web" / "index.html").read_text(encoding="utf-8")
    names = {name for name, _ in label_fit._SETTING_NAMES.values()}
    names |= {name for name, _ in label_fit._MODE_SETTING_NAMES.values()}
    for name in sorted(names):
        text = re.escape(name)
        named = (
            re.search(rf'<(?:input|select)\b[^>]*\baria-label="{text}"', html)
            or re.search(rf'<label\b[^>]*\bfor="[^"]+"[^>]*>{text}</label>', html)
            or re.search(rf"<label\b[^>]*>\s*<input\b[^>]*>\s*<span>{text}</span>\s*</label>", html)
        )
        assert named, f"no control is named {name!r}"


class _Rows:
    def __init__(self) -> None:
        self.asked: list[str] = []

    def row(self, kind: LabelFitRowKind) -> LabelFitRow:
        self.asked.append(kind)
        return LabelFitRow(kind, "found", {path: 11.0 for path in FONT_SIZE_PATHS}, (800.0, 800.0))


def test_a_failure_with_rows_to_compute_is_kept_and_found_through_the_exception_chain() -> None:
    report = _report(drawn=False)
    rows = _Rows()
    error = LabelsDoNotFitError("labels overlap", report=report, row_search=rows)
    try:
        try:
            raise error
        except LabelsDoNotFitError as raised:
            raise ValidationError("wrapped") from raised
    except ValidationError as wrapped:
        assert label_fit_failure(wrapped) is error
    assert label_fit_failure(ValueError("other")) is None

    searches = LabelFitRowSearches()
    assert searches.failure_payload(error, "circular")["rowToken"] is None
    token = searches.retain_failure("circular", error)
    payload = searches.failure_payload(error, "circular")
    assert payload == {"report": label_fit_report_payload(report, "circular"), "rowToken": token}
    assert searches.row(str(token), "combination")["status"] == "found" and rows.asked == ["combination"]
    # A newer render of the mode releases the failure.
    searches.retain("circular", [])
    assert searches.row(str(token), "combination")["status"] == "not_computed"
    # A failure whose rows are all computed keeps no token.
    done = LabelsDoNotFitError("labels overlap", report=LabelFitReport(records=(), rows=(), drawn=False))
    assert searches.retain_failure("circular", done) is None


def _sars_request(overrides: dict[str, object]) -> dict:
    record = SeqIO.read(str(TEST_INPUTS / "SARS-CoV-2.gbk"), "genbank")
    return build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(config_overrides={"labels.circular.scope": "outer", **overrides}),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    ).to_dict()


def _crossing_request() -> dict:
    # External labels only: crossing leaders, no overlap (the embedded labels would overlap).
    return _sars_request({
        "labels.rendering": "external_only",
        "labels.font_size.short": 14.0,
        "labels.font_size.long": 14.0,
    })


def _web_render(document: dict, tmp_path: Path, searches: LabelFitRowSearches | None = None) -> dict:
    from gbdraw.session_io import materialize_embedded_file

    resources = tmp_path / "resources"
    resources.mkdir(parents=True, exist_ok=True)
    resource_paths = {
        resource_id: materialize_embedded_file(entry, temp_dir=resources, role=resource_id, prefix_role=False)
        for resource_id, entry in document["resources"].items()
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return render_canonical_web_request(
            document["renderRequest"],
            resource_paths=resource_paths,
            output_directory=tmp_path / "out",
            _label_fit_searches=searches,
        )


def test_overprinted_embedded_labels_offer_only_the_resolve_overlaps_control(tmp_path) -> None:
    """LP-D52 (E-ROW): SARS-CoV-2 with labels out at defaults stops; Resolve Overlaps fits it."""

    searches = LabelFitRowSearches()
    with pytest.raises(LabelsDoNotFitError) as raised:
        _web_render(_sars_request({}), tmp_path, searches)
    payload = searches.failure_payload(raised.value, "circular")
    found = [row for row in payload["report"]["rows"] if row["status"] == "found"]
    assert [row["kind"] for row in found] == ["separate_features"]
    assert found[0]["values"] == {SEPARATE_FEATURES_PATH: True}
    assert found[0]["settings"] == [{"name": "Resolve Overlaps", "unit": "", "value": "On"}]

    fitted = _web_render(_sars_request({SEPARATE_FEATURES_PATH: True}), tmp_path / "fitted")
    assert fitted["metadata"]["labelFitReports"] == []


def test_with_separate_strands_the_resolve_overlaps_row_leaves_crossings_and_draws_with_the_notice(tmp_path) -> None:
    """LP-D57 (OV-724): with the Web's Separate Strands on, Resolve Overlaps removes the
    overlap but leaders still cross; the row is offered, says so, and draws with the notice."""

    strands = {"canvas.strandedness": True}
    searches = LabelFitRowSearches()
    with pytest.raises(LabelsDoNotFitError) as raised:
        _web_render(_sars_request(strands), tmp_path, searches)
    rows = searches.failure_payload(raised.value, "circular")["report"]["rows"]
    offered = [row for row in rows if row["values"] is not None]
    assert [(row["kind"], row["status"]) for row in offered] == [("separate_features", "found_with_notice")]
    (separate,) = offered
    assert separate["values"] == {SEPARATE_FEATURES_PATH: True}
    assert separate["settings"] == [{"name": "Resolve Overlaps", "unit": "", "value": "On"}]
    assert separate["canvasSize"] is not None
    assert separate["note"] == "14 leader pairs still cross"
    assert {row["note"] for row in rows if row is not separate} == {None}
    assert "  --resolve_overlaps   (separate overlapping features; 14 leader pairs still cross)\n" in str(raised.value)

    drawn = _web_render(_sars_request({**strands, SEPARATE_FEATURES_PATH: True}), tmp_path / "drawn")
    (entry,) = drawn["metadata"]["labelFitReports"]
    assert entry["report"]["drawn"] is True
    assert (entry["report"]["overlapPairs"], entry["report"]["crossingPairs"]) == (0, 14)


def test_web_render_draws_crossing_leaders_with_counts_only_and_computes_rows_later(tmp_path) -> None:
    searches = LabelFitRowSearches()
    response = _web_render(_crossing_request(), tmp_path, searches)

    (entry,) = response["metadata"]["labelFitReports"]
    assert entry["resultIndex"] == 0 and entry["resultName"] == "fit.svg"
    report = entry["report"]
    assert report["drawn"] is True and report["overlapPairs"] == 0 and report["crossingPairs"] > 0
    # OD-F14: Generate searched no row.
    assert [row["kind"] for row in report["rows"]] == list(ROW_KINDS[:6])
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


def test_a_row_computed_after_the_conservation_file_is_gone_reads_no_file(tmp_path, monkeypatch, caplog) -> None:
    """OV-720: the kept figure reassembles from the conservation rows read once."""

    import logging

    import gbdraw.analysis.conservation as conservation
    from gbdraw.api.diagram import build_circular_diagram, label_fit_report, label_fit_row

    record = SeqIO.read(str(TEST_INPUTS / "SARS-CoV-2.gbk"), "genbank")
    ring = tmp_path / "ring.tsv"
    ring.write_text(
        "".join(
            f"{record.id}\tother\t95.0\t1000\t50\t0\t{start}\t{start + 999}\t{start}\t{start + 999}\t1e-50\t1500\n"
            for start in range(1, len(record.seq) - 1000, 3000)
        ),
        encoding="utf-8",
    )
    options = CircularDiagramOptions(
        config_overrides={
            "labels.circular.scope": "outer",
            "labels.rendering": "external_only",
            "labels.font_size.short": 14.0,
            "labels.font_size.long": 14.0,
        },
        conservation_blast_files=(str(ring),),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(record, options=options, _label_fit_warn_search_seconds=0.0)
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert {row.status for row in report.rows} == {"not_computed"}

    ring.unlink()
    reads: list[object] = []
    read = conservation.read_comparison_table

    def counted(*args: Any, **kwargs: Any) -> Any:
        reads.append(args)
        return read(*args, **kwargs)

    monkeypatch.setattr(conservation, "read_comparison_table", counted)
    with caplog.at_level(logging.WARNING), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        row = label_fit_row(drawing, "outward_vertical")
    assert row is not None and row.status in {"found", "not_found"}
    assert reads == []
    assert "Skipping conservation source" not in caplog.text
