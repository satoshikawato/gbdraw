"""The label fit report the Web reads, and the rows it computes on demand (OD-F14)."""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any
import warnings

from Bio import SeqIO
import pytest

from tests.utils.web_markup import read_english_index

from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularRequestTrackOptions,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.exceptions import LabelsDoNotFitError, LabelsMayNotFitWarning, ValidationError
from gbdraw.labels.circular_fit import ROW_KINDS
from gbdraw.layout.figure_size import LargerFigure
from gbdraw.labels.fit import (
    DEFINITION_FONT_PATH,
    FONT_SIZE_PATHS,
    PLACEMENT_PATH,
    SCOPE_PATH,
    SEPARATE_FEATURES_PATH,
    SEPARATE_FEATURES_VALUES,
    Y_OFFSET_PATH,
    LabelCapacityRecord,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelKind,
    LimitReason,
    PlainReason,
    UnfitReason,
    format_label_fit_report,
)
from gbdraw.session import build_session_document
from tests.utils.session_files import embedded_session
from tests.utils.session_drawing_view import drawing_view
from gbdraw.web_support.label_fit import (
    WEB_ROW_SEARCH_SECONDS,
    LabelFitRowSearches,
    label_capacity_notices,
    label_fit_failure,
    label_fit_report_payload,
)
from gbdraw.web_support.error_adapter import serialize_web_error
from gbdraw.web_support.request_render import render_canonical_web_request

TEST_INPUTS = Path(__file__).parent / "test_inputs"


def _embedded_document(document):
    return embedded_session(document.to_dict(), document.resources.sources)


def _report(*, drawn: bool) -> LabelFitReport:
    return LabelFitReport(
        records=(LabelFitRecord("rec-1", "circular_outer", 111, 14.0, LabelFit(20, 21)),),
        rows=(
            LabelFitRow("outward_vertical", "found", {Y_OFFSET_PATH: 1.15}, (1708.0, 1130.0)),
            LabelFitRow("smaller_font", "found", {path: 12.0 for path in FONT_SIZE_PATHS}, (1615.0, 988.0)),
            LabelFitRow("tighter_rows", "not_found", reason=LimitReason("spacing_limit", 0.1)),
            LabelFitRow("combination", "not_computed"),
        ),
        drawn=drawn,
    )


def test_report_payload_carries_counts_rows_and_python_wording() -> None:
    payload = label_fit_report_payload(_report(drawn=False), "circular")

    assert payload["schema"] == 1 and payload["drawn"] is False
    assert (payload["overlapPairs"], payload["crossingPairs"]) == (20, 21)
    assert payload["records"] == [{
        "recordId": "rec-1", "recordIndex": None, "labelKind": "circular_outer", "labelKindTitle": "outer labels",
        "labelCount": 111, "fontSize": 14.0, "overlapPairs": 20, "crossingPairs": 21,
        "definitionOverlaps": 0, "regionOverlaps": 0, "refused": False, "capacityRatio": None,
        "counts": "20 label pairs overlap, 21 leader pairs cross",
    }]
    assert payload["remedies"] == []
    vertical, font, tighter, combination = payload["rows"]
    assert vertical == {
        "kind": "outward_vertical", "title": "outward, vertical only", "status": "found",
        "values": {Y_OFFSET_PATH: 1.15}, "trackOptions": {},
        "settings": [{"name": "Outer Y Offset", "unit": "", "value": 1.15}],
        "canvasSize": {"width": 1708.0, "height": 1130.0}, "reason": None, "remaining": [],
    }
    # One control sets both font size classes: one setting, both config paths.
    assert font["values"] == {path: 12.0 for path in FONT_SIZE_PATHS}
    assert font["settings"] == [{"name": "Label Font Size", "unit": "pt", "value": 12.0}]
    assert tighter["status"] == "not_found" and tighter["values"] is None and tighter["settings"] == []
    assert tighter["reason"] == {"kind": "spacing_limit", "limit": 0.1}
    assert combination["status"] == "not_computed" and combination["canvasSize"] is None


def test_row_reasons_and_notes_reach_the_web_as_data_not_english() -> None:
    """OV-3109: why a row has no value and what a row that draws with the notice leaves are kinds, numbers and counts
    by field name, which the Web words in the shown language; the payload carries no English sentence for them."""

    from gbdraw.labels.fit import CombinationLimitReason

    notice = LabelFitRow(
        "outward_vertical", "found_with_notice", {Y_OFFSET_PATH: 1.3}, (1.0, 1.0),
        remaining=LabelFit(1, 2, 3, region_overlaps=1, region_cuts=2),
    )
    definition = LabelFitRow(
        "narrow_tracks", "found_with_notice", {DEFINITION_FONT_PATH: 12.0}, (1.0, 1.0),
        remaining=LabelFit(0, 4, definition_track_overlaps=2),
    )
    rows = (
        LabelFitRow("outward_vertical", "not_found", reason=LimitReason("y_offset_limit", 2.0)),
        LabelFitRow("combination", "not_found", reason=CombinationLimitReason(13.0, 2.0)),
        LabelFitRow("labels_out", "not_found", reason=UnfitReason(LabelFit(0, 1))),
        LabelFitRow("radial", "not_found", reason=UnfitReason(None)),
        LabelFitRow("outward_horizontal", "not_found", reason=PlainReason("outer_labels_only")),
        notice,
        definition,
    )
    report = LabelFitReport(records=(LabelFitRecord("rec-1", "circular_outer", 3, 14.0, LabelFit(1, 1)),), rows=rows,
                            drawn=False)
    payload_rows = label_fit_report_payload(report, "circular")["rows"]
    assert [row["reason"] for row in payload_rows] == [
        {"kind": "y_offset_limit", "limit": 2.0},
        {"kind": "combination_limit", "fontSize": 13.0, "limit": 2.0},
        {"kind": "labels_do_not_fit", "counts": [
            {"field": "overlap_pairs", "count": 0}, {"field": "crossing_pairs", "count": 1},
        ]},
        {"kind": "labels_do_not_fit", "counts": []},
        {"kind": "outer_labels_only"},
        None,
        None,
    ]
    # A label row names crossings, labels over the definition text and region label warnings (not its overlapping
    # pairs); a definition row names only its definition lines over a track.
    assert payload_rows[5]["remaining"] == [
        {"field": "crossing_pairs", "count": 2}, {"field": "definition_overlaps", "count": 3},
        {"field": "region_overlaps", "count": 1}, {"field": "region_cuts", "count": 2},
    ]
    assert payload_rows[6]["remaining"] == [{"field": "definition_track_overlaps", "count": 2}]
    assert all("note" not in row for row in payload_rows)
    assert not any(word in json.dumps(payload_rows) for word in ("Does not", "still", "Labels", "drawn with"))


def test_labels_over_capacity_carry_their_ratio_the_remedies_and_their_own_reason() -> None:
    """CAPACITY-STOP: the dialog lists Python's sentence per record and the remedies; the banner has its own reason."""

    capacity = LabelCapacityRecord(562, 6.0, 2.49, ("filter", "label_text", "linear"))
    report = LabelFitReport(
        records=(LabelFitRecord("AP027078.1", "circular_outer", 562, 6.0, LabelFit(0, 0), capacity=capacity),),
        rows=(),
        drawn=False,
    )
    payload = label_fit_report_payload(report, "circular")

    (record,) = payload["records"]
    # OV-2179: the ratio as shown (one decimal place, Python's rounding); the Web words the record from its values.
    assert record["capacityRatio"] == 2.5
    assert (record["labelKind"], record["labelCount"], record["fontSize"]) == ("circular_outer", 562, 6.0)
    assert payload["rows"] == [] and payload["remedies"] == ["filter", "label_text", "linear"]
    # At a requested size the Web gets the remedy and Python's estimate of each requested side, as its messages
    # write it; an Auto report has no estimate key, so its payload is unchanged.
    assert "largerFigure" not in payload
    sized = label_fit_report_payload(replace(report, larger_figure=LargerFigure("240.3 mm", None)), "circular")
    assert sized["remedies"] == ["filter", "label_text", "figure_size", "linear"]
    assert sized["largerFigure"] == {"width": "240.3 mm", "height": None}
    error = serialize_web_error(LabelsDoNotFitError("PRIVATE", report=report), operation="generate", stage="render")
    assert error["context"] == {"reason": "LABEL_CAPACITY", "labelCapacityOuter": 2.5}
    # Both sides: the largest room of each.
    inner = LabelCapacityRecord(40, 6.0, 1.84, ("filter", "label_text", "linear"))
    both = LabelFitReport(
        records=(
            *report.records,
            LabelFitRecord("AP027078.1", "circular_inner", 40, 6.0, LabelFit(0, 0), capacity=inner),
        ),
        rows=(),
        drawn=False,
    )
    error = serialize_web_error(LabelsDoNotFitError("PRIVATE", report=both), operation="generate", stage="render")
    assert error["context"] == {"reason": "LABEL_CAPACITY", "labelCapacityInner": 1.8, "labelCapacityOuter": 2.5}
    # Nothing to compute: the failure keeps no row token.
    searches = LabelFitRowSearches()
    failure = LabelsDoNotFitError("PRIVATE", report=report)
    assert searches.retain_failure("circular", failure) is None


def test_a_report_with_no_fitting_row_lists_the_ways_to_make_room() -> None:
    """NO-ROW-WAYS (OV-2180, PD-OI-097): when every row is known not to fit, the report names the ways to make room the
    capacity stop names (the same owner, so a smaller font is left out at the 6 pt floor), and the CLI text lists them."""

    not_found = (
        LabelFitRow("outward_both", "not_found", reason=LimitReason("xy_offset_limit", 2.0)),
        LabelFitRow("smaller_font", "not_found", reason=LimitReason("font_limit", 6.0)),
    )

    def report(font_size: float, rows: tuple[LabelFitRow, ...] = not_found, kind: LabelKind = "circular_outer") -> LabelFitReport:
        record = LabelFitRecord("AP027078.1", kind, 562, font_size, LabelFit(30, 2))
        return LabelFitReport(records=(record,), rows=rows, drawn=False)

    at_floor = report(6.0)
    assert at_floor.room_remedies == ("filter", "label_text", "linear")
    assert label_fit_report_payload(at_floor, "circular")["remedies"] == ["filter", "label_text", "linear"]
    assert report(8.0).room_remedies == ("filter", "label_text", "font", "linear")
    assert format_label_fit_report(at_floor).splitlines()[1:] == [
        "No single change of the label settings fits this map.",
        "To make room:",
        "  - label fewer features with --label_blacklist or --label_whitelist;",
        "  - shorten the label text with --qualifier_priority (for example gene names only);",
        "  - or draw a Linear map (gbdraw linear).",
    ]
    # A row that fits: no ways.
    found = report(6.0, (*not_found, LabelFitRow("outward_vertical", "found", {Y_OFFSET_PATH: 1.5}, (1.0, 1.0))))
    assert found.room_remedies == () and "To make room:" not in format_label_fit_report(found)
    assert label_fit_report_payload(found, "circular")["remedies"] == []
    # Rows not computed yet may fit: the report carries the ways for the dialog, which shows them once every row is
    # known not to fit; the CLI text names the rows not searched and lists no ways.
    pending = report(6.0, (not_found[0], LabelFitRow("combination", "not_computed")))
    assert pending.room_remedies == ("filter", "label_text", "linear")
    assert "To make room:" not in format_label_fit_report(pending)
    # Only the labels the capacity check reads (outer and inner Circular labels) have these ways.
    assert report(6.0, kind="linear").room_remedies == ()
    # At a requested size a larger figure is also a way here, named without an estimate (PD-OI-114 r3): the Web gets
    # the requested sides, and an Auto report keeps its payload and words (above).
    sized = replace(at_floor, larger_figure=LargerFigure(True, None))
    assert sized.room_remedies == ("filter", "label_text", "figure_size", "linear")
    payload = label_fit_report_payload(sized, "circular")
    assert payload["remedies"] == ["filter", "label_text", "figure_size", "linear"]
    assert payload["largerFigure"] == {"width": True, "height": None}
    assert "largerFigure" not in label_fit_report_payload(at_floor, "circular")
    assert "  - make the figure larger: increase --figure_width;" in format_label_fit_report(sized).splitlines()
    # No ways without a horizontal label kind, so none of them is a larger figure either.
    assert replace(report(6.0, kind="linear"), larger_figure=LargerFigure(True, None)).room_remedies == ()


def test_capacity_warnings_reach_the_notice_while_other_warnings_pass_on() -> None:
    notices: list[str] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with label_capacity_notices(notices.append):
            capacity = LabelCapacityRecord(195, 14.0, 1.36, ())
            warnings.warn(
                LabelsMayNotFitWarning("SENTENCE", record_id="NC_005859.1", label_kind="circular_outer", capacity=capacity)
            )
            warnings.warn(LabelsMayNotFitWarning("AGAIN", record_id="r", label_kind="circular_inner", capacity=capacity))
            warnings.warn(UserWarning("other"))
    # OV-2179: the values, not Python's English sentence, so the Web words them in the UI language.
    assert [json.loads(notice) for notice in notices] == [
        {"recordId": "NC_005859.1", "labelKind": "circular_outer", "labelCount": 195, "fontSize": 14.0,
         "capacityRatio": 1.4},
        {"recordId": "r", "labelKind": "circular_inner", "labelCount": 195, "fontSize": 14.0, "capacityRatio": 1.4},
    ]
    assert [str(item.message) for item in caught] == ["other"]


def test_each_fit_pass_sends_its_capacity_notice_while_the_caller_gets_the_kept_pass_once() -> None:
    """D-17: the Web overlay gets each figure-size fit pass's CAPACITY-WARN when it is warned (progress); the caller's
    warnings stay those of the kept pass, once (D-16)."""

    from gbdraw.labels.circular_capacity import CircularLabelCapacity, check_label_capacity
    from gbdraw.layout.figure_size import MeasuredBase, fit_base_length

    notices: list[str] = []
    events: list[str] = []

    def measure(base: int) -> MeasuredBase[int]:
        events.append(f"pass {base} with {len(notices)} notices")
        ratio = round(base / 1000, 2)
        check_label_capacity("r", [CircularLabelCapacity(
            "circular_outer", LabelCapacityRecord(195, 14.0, ratio, ()), 1.0, "near")])
        return MeasuredBase(base + 100.0 - 500.0, base, (base + 100.0,))

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with label_capacity_notices(notices.append):
            fit = fit_base_length(measure, first_base=500, slope=1.0)
    assert fit.base == 398
    assert events == ["pass 500 with 0 notices", "pass 398 with 1 notices"]
    assert [json.loads(notice)["capacityRatio"] for notice in notices] == [0.5, 0.4]
    assert caught == []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit_base_length(measure, first_base=500, slope=1.0)
    assert [item.message.capacity.ratio for item in caught] == [0.4]


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

    html = read_english_index(Path(__file__).parents[1] / "gbdraw" / "web" / "index.html")
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
        self.deadlines: list[float | None] = []

    def row(self, kind: LabelFitRowKind, *, deadline: float | None = None) -> LabelFitRow:
        self.asked.append(kind)
        self.deadlines.append(deadline)
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
    return drawing_view(_embedded_document(build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(config_overrides={"labels.circular.scope": "outer", **overrides}),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))


# Lambda with outer labels at 17 pt and no blacklist: 10 leader pairs cross and no labels overlap, so
# the figure is drawn with the notice. (16 pt did while labels were laid out in point-mode boxes, OV-165.
# SARS-CoV-2 with external labels only at 14 pt crossed only because labels next to the origin could wrap
# past their neighbours, OV-1366; it fits now.)
CROSSING_OVERRIDES: dict[str, object] = {
    "labels.circular.scope": "outer",
    "labels.filtering.blacklist_keywords": [],
    "labels.font_size.short": 17.0,
    "labels.font_size.long": 17.0,
}


def _crossing_record():
    return SeqIO.read(str(TEST_INPUTS / "NC_001416.gb"), "genbank")


def _crossing_options() -> CircularDiagramOptions:
    return CircularDiagramOptions(config_overrides=dict(CROSSING_OVERRIDES))


def _crossing_request() -> dict:
    return drawing_view(_embedded_document(build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(_crossing_record())),),
            options=_crossing_options(),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))


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
    # The labels fit. The definition, never wrapped, is drawn over the tracks with its notice (DF-D02, DF-D25).
    (entry,) = fitted["metadata"]["labelFitReports"]
    assert {record["labelKind"] for record in entry["report"]["records"]} == {"circular_definition"}


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
    # What the row leaves is data the Web words in its language (OV-3109); the CLI words it in English.
    assert separate["remaining"] == [{"field": "crossing_pairs", "count": 15}] and "note" not in separate
    assert [row["remaining"] for row in rows if row is not separate] == [[]] * (len(rows) - 1)
    assert "  --resolve_overlaps   (separate overlapping features; 15 leader pairs still cross)\n" in str(raised.value)

    drawn = _web_render(_sars_request({**strands, SEPARATE_FEATURES_PATH: True}), tmp_path / "drawn")
    (entry,) = drawn["metadata"]["labelFitReports"]
    assert entry["report"]["drawn"] is True
    assert (entry["report"]["overlapPairs"], entry["report"]["crossingPairs"]) == (0, 15)


def _mers_request(
    center_reserved_radius: float | None = None, *, track_narrowing: bool = True, definition_font: float | None = None
) -> dict:
    record = SeqIO.read(str(TEST_INPUTS / "MERS-CoV.gbk"), "genbank")
    overrides: dict[str, object] = {}
    if not track_narrowing:
        overrides["canvas.circular.track_narrowing"] = False
    if definition_font is not None:
        overrides[DEFINITION_FONT_PATH] = definition_font
    return drawing_view(_embedded_document(build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                tracks=CircularRequestTrackOptions(center_reserved_radius=center_reserved_radius),
                config_overrides=overrides or None,
            ),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))


def test_a_definition_over_the_tracks_offers_the_center_reserved_radius_that_narrows_them(tmp_path) -> None:
    """PD-OI-097 (N-a): MERS-CoV's definition is drawn on one line over a track. At 18 pt its band (the farthest
    corner of the bold line's ink, OV-1854) is cleared neither by the floor nor by any Center Reserved Radius, so the
    narrow row finds nothing and the font row offers 17.5 pt. At 17.5 pt with track narrowing off, the font row offers
    the size at which that line clears the tracks at their widths, and the narrow row writes the radius the
    definition reserves as Center Reserved Radius; drawn with it, the figure has no notice."""

    (default,) = _web_render(_mers_request(), tmp_path / "default")["metadata"]["labelFitReports"]
    assert default["report"]["drawn"] is True
    assert [row["kind"] for row in default["report"]["rows"]] == ["smaller_definition_font", "narrow_tracks"]

    searches = LabelFitRowSearches()
    (entry,) = _web_render(_mers_request(track_narrowing=False, definition_font=17.5), tmp_path, searches)[
        "metadata"
    ]["labelFitReports"]
    assert entry["report"]["drawn"] is True
    # The record is named by its drawn index too: "Edit Species / Strain…" finds its definition by it.
    assert [(record["recordId"], record["recordIndex"], record["labelKind"]) for record in entry["report"]["records"]] == [
        ("NC_019843.3", 0, "circular_definition")
    ]
    assert [row["kind"] for row in entry["report"]["rows"]] == ["smaller_definition_font", "narrow_tracks"]
    font = searches.row(entry["rowToken"], "smaller_definition_font")
    # 14.5 pt since the inside tracks moved outward with the tick labels inside the scale (OV-3300; 14 pt before, and
    # 13 pt before the GC skew windows were centred, OV-2870): the skew ring the line meets moved.
    assert (font["status"], font["values"]) == ("found", {DEFINITION_FONT_PATH: 14.5})
    narrow = searches.row(entry["rowToken"], "narrow_tracks")
    # The radius is a track option of the request, not a config path: the Web writes it by the request's name.
    assert (narrow["status"], narrow["values"], narrow["trackOptions"]) == (
        "found", {}, {"centerReservedRadius": 247.31}
    )
    assert font["trackOptions"] == {}
    assert narrow["settings"] == [{"name": "Center Reserved Radius", "unit": "px", "value": 247.31}]
    assert narrow["remaining"] == []

    narrowed = _web_render(
        _mers_request(247.31, track_narrowing=False, definition_font=17.5), tmp_path / "narrowed"
    )["metadata"]["labelFitReports"]
    assert narrowed == []


@pytest.mark.parametrize("radius", [200.0, 100.0])
def test_an_explicit_center_reserved_radius_leaves_out_the_narrow_row(tmp_path, radius: float) -> None:
    """PD-OI-097 (N-a): the narrow row writes ``center_reserved_radius`` only when none is set; with one set that
    leaves the definition over a track, the notice offers the smaller font alone."""

    (entry,) = _web_render(_mers_request(radius), tmp_path)["metadata"]["labelFitReports"]
    assert entry["report"]["drawn"] is True
    assert [row["kind"] for row in entry["report"]["rows"]] == ["smaller_definition_font"]


def test_web_render_draws_crossing_leaders_with_counts_only_and_computes_rows_later(tmp_path) -> None:
    searches = LabelFitRowSearches()
    response = _web_render(_crossing_request(), tmp_path, searches)

    (entry,) = response["metadata"]["labelFitReports"]
    assert entry["resultIndex"] == 0 and entry["resultName"] == "fit.svg"
    report = entry["report"]
    assert report["drawn"] is True and report["overlapPairs"] == 0 and report["crossingPairs"] > 0
    # Python owns the notice bar's words (the Owner's Q8 wording).
    assert report["notice"] == [f"Leader lines cross ({report['crossingPairs']} pairs)."]
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
    assert stale["status"] == "not_computed" and stale["reason"] == {"kind": "stale"}
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

    record = _crossing_record()
    ring = tmp_path / "ring.tsv"
    ring.write_text(
        "".join(
            f"{record.id}\tother\t95.0\t1000\t50\t0\t{start}\t{start + 999}\t{start}\t{start + 999}\t1e-50\t1500\n"
            for start in range(1, len(record.seq) - 1000, 3000)
        ),
        encoding="utf-8",
    )
    options = CircularDiagramOptions(config_overrides=dict(CROSSING_OVERRIDES), conservation_blast_files=(str(ring),))
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


# The OV-732 figures' band and its label, which stays clear of the feature labels and leaders. A region label that
# warns (over them, OV-3015) is one no outer label row clears, so each such row would verify every grid value up to
# its limit: about 120 full assemblies per figure in place of 5 (OV-3104). The tests pin the files, not that, and
# check that no region label warns.
OV732_BAND = "12000\t12500\tregion"  # start, end, label


def _write_input_files(directory: Path, record) -> dict[str, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    texts = {
        "depth": "".join(f"{record.id}\t{i}\t{10 + i % 50}\n" for i in range(1, len(record.seq) + 1)),
        "colors": "CDS\tproduct\t.*\t#D55E00\tCDS\n",
        "default_colors": "CDS\t#D0D0D0\n",
        "visibility": "*\tCDS\tproduct\t.*\tshow\n",
        "annotations": f"set_id\tid\tmark\tstart\tend\tlabel\ns\ta\tband\t{OV732_BAND}\n",
    }
    paths = {name: directory / f"{name}.tsv" for name in texts}
    for name, text in texts.items():
        paths[name].write_text(text, encoding="utf-8")
    return paths


def _crossing_options_with_files(
    paths: dict[str, Path], *, extra_overrides: dict[str, object] | None = None
) -> CircularDiagramOptions:
    from gbdraw.annotations.models import AnnotationOptions
    from gbdraw.api.options import ColorOptions

    return CircularDiagramOptions(
        config_overrides={**CROSSING_OVERRIDES, "canvas.show_depth": True, **(extra_overrides or {})},
        colors=ColorOptions(
            color_table_file=str(paths["colors"]), default_colors_file=str(paths["default_colors"])
        ),
        feature_visibility_table_file=str(paths["visibility"]),
        annotations=AnnotationOptions(table_file=str(paths["annotations"])),
        depth_file=str(paths["depth"]),
    )


def test_rows_computed_after_every_input_file_is_gone_match_the_rows_with_the_files(tmp_path) -> None:
    """OV-732: the builder reads each file input once; a row reassembles the figure without them."""

    from gbdraw.api.diagram import build_circular_diagram, label_fit_report, label_fit_row

    record = _crossing_record()
    drawings = {}
    for name in ("kept", "gone"):
        paths = _write_input_files(tmp_path / name, record)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            drawings[name] = build_circular_diagram(
                record, options=_crossing_options_with_files(paths), _label_fit_warn_search_seconds=0.0
            )
    for path in (tmp_path / "gone").iterdir():
        path.unlink()

    report = label_fit_report(drawings["gone"])
    assert report is not None and report.drawn and report.rows
    assert "circular_region" not in {record.label_kind for record in report.records}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for row in report.rows:
            assert label_fit_row(drawings["gone"], row.kind) == label_fit_row(drawings["kept"], row.kind)


# The OV-732 tests pin that rows compute without the input files, not how fast, so their searches have no
# bound: a bounded row stops once the next verification, as long as the last, would end after the bound
# (CircularLabelFitSearch), and the time of the same verifications varies widely from run to run.
# The bound itself is pinned by test_a_web_row_stops_at_its_bound_and_says_no_value_was_found.
NO_ROW_BOUND = math.inf


def test_web_rows_computed_after_the_staged_depth_file_is_gone_match_the_rows_with_it(tmp_path) -> None:
    """OV-732: the Web deletes the staged files after Generate; "Choose a fix" rows still compute."""

    import shutil

    record = _crossing_record()
    paths = _write_input_files(tmp_path / "inputs", record)
    document = drawing_view(_embedded_document(build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=_crossing_options_with_files(paths),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))
    searches = {name: LabelFitRowSearches(row_search_seconds=NO_ROW_BOUND) for name in ("kept", "gone")}
    entries = {}
    for name, search in searches.items():
        (entries[name],) = _web_render(document, tmp_path / name, search)["metadata"]["labelFitReports"]
    shutil.rmtree(tmp_path / "gone" / "resources")
    shutil.rmtree(tmp_path / "inputs")

    kinds = [row["kind"] for row in entries["gone"]["report"]["rows"]]
    assert kinds and "circular_region" not in {record["labelKind"] for record in entries["gone"]["report"]["records"]}
    for kind in kinds:
        gone = searches["gone"].row(entries["gone"]["rowToken"], kind)
        assert gone == searches["kept"].row(entries["kept"]["rowToken"], kind)
        assert gone["status"] != "not_computed"


# Two lambda records drawn together: a 12 pt definition stays clear of their tracks, so the rows are the label rows
# only (the definition rows redraw both records per font step and would triple the run, DF-D02).
TWO_RECORD_OVERRIDES: dict[str, object] = {"objects.definition.circular.font_size": 12}


def _two_records_with_files(directory: Path) -> tuple[list[Any], dict[str, Path]]:
    import copy

    base = _crossing_record()
    records = []
    for record_id in ("LAMBDA_A", "LAMBDA_B"):
        record = copy.deepcopy(base)
        record.id = record.name = record_id
        records.append(record)
    paths = _write_input_files(directory, records[0])
    paths["depth"].write_text(
        "".join(
            f"{record.id}\t{i}\t{10 + i % 50}\n" for record in records for i in range(1, len(record.seq) + 1)
        ),
        encoding="utf-8",
    )
    paths["annotations"].write_text(
        f"set_id\tid\tmark\trecord\tstart\tend\tlabel\ns\ta\tband\tLAMBDA_B\t{OV732_BAND}\n", encoding="utf-8"
    )
    return records, paths


def test_rows_of_two_circular_records_computed_after_every_input_file_is_gone(tmp_path) -> None:
    """OV-732: the multi-record builder reads the Depth, color, visibility, and annotation files once."""

    from gbdraw.api.diagram import build_circular_multi_diagram, label_fit_report, label_fit_row

    drawings = {}
    for name in ("kept", "gone"):
        records, paths = _two_records_with_files(tmp_path / name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            drawings[name] = build_circular_multi_diagram(
                records,
                options=_crossing_options_with_files(paths, extra_overrides=TWO_RECORD_OVERRIDES),
                _label_fit_warn_search_seconds=0.0,
            )
    for path in (tmp_path / "gone").iterdir():
        path.unlink()

    report = label_fit_report(drawings["gone"])
    assert report is not None and report.drawn and report.rows
    assert "circular_region" not in {record.label_kind for record in report.records}
    assert {record.record_id for record in report.records} == {"LAMBDA_A", "LAMBDA_B"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for row in report.rows:
            assert label_fit_row(drawings["gone"], row.kind) == label_fit_row(drawings["kept"], row.kind)


def test_web_rows_of_two_circular_records_compute_after_their_input_files_are_gone(tmp_path) -> None:
    """OV-732: on the Web, a multi-record figure's rows compute after the staged files are deleted."""

    import shutil

    records, paths = _two_records_with_files(tmp_path / "inputs")
    document = drawing_view(_embedded_document(build_session_document(
        CircularDiagramRequest(
            records=tuple(RecordInput(source=InMemoryRecordSource(record)) for record in records),
            options=_crossing_options_with_files(paths, extra_overrides=TWO_RECORD_OVERRIDES),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))
    searches = {name: LabelFitRowSearches(row_search_seconds=NO_ROW_BOUND) for name in ("kept", "gone")}
    entries = {}
    for name, search in searches.items():
        (entries[name],) = _web_render(document, tmp_path / name, search)["metadata"]["labelFitReports"]
    shutil.rmtree(tmp_path / "gone" / "resources")
    shutil.rmtree(tmp_path / "inputs")

    assert {record["recordId"] for record in entries["gone"]["report"]["records"]} == {"LAMBDA_A", "LAMBDA_B"}
    kinds = [row["kind"] for row in entries["gone"]["report"]["rows"]]
    assert kinds and "circular_region" not in {record["labelKind"] for record in entries["gone"]["report"]["records"]}
    for kind in kinds:
        gone = searches["gone"].row(entries["gone"]["rowToken"], kind)
        assert gone == searches["kept"].row(entries["kept"]["rowToken"], kind)
        assert gone["status"] != "not_computed"


def _default_config():
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.toml import load_config_toml

    return GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))


def test_a_row_assembly_that_fails_for_a_reason_other_than_labels_raises_it() -> None:
    """OV-733: only a label setting the engine cannot draw counts as "does not fit"."""

    from gbdraw.api.diagram import _checked_circular_drawing, build_circular_diagram, label_fit_row
    from gbdraw.exceptions import LabelLayoutError

    record = _crossing_record()
    options = _crossing_options()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(record, options=options, _label_fit_warn_search_seconds=0.0)

    def failing(error: Exception):
        def assemble(_cfg, overrides):
            if not overrides:
                return drawing
            raise error

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return _checked_circular_drawing(assemble, _default_config(), warn_search_seconds=0.0, map_item="map")

    with pytest.raises(ValidationError, match="Depth file does not exist"):
        label_fit_row(failing(ValidationError("Depth file does not exist: gone.tsv")), "smaller_font")
    row = label_fit_row(failing(LabelLayoutError("radial label packing did not converge")), "smaller_font")
    assert row is not None and row.status == "not_found"


def test_a_row_whose_tracks_do_not_fit_does_not_fit() -> None:
    """OV-1266: a row's settings that leave the tracks no room (TRACK_LAYOUT) do not fit; Generate goes on."""

    from gbdraw.api.diagram import (
        _checked_circular_drawing,
        _unfit_row_verification,
        build_circular_diagram,
        label_fit_row,
    )

    record = _crossing_record()
    options = _crossing_options()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(record, options=options, _label_fit_warn_search_seconds=0.0)
    error = ValidationError(
        "Circular track slot 'features' cannot fit inside",
        diagnostic={"code": "TRACK_LAYOUT", "reason": "CENTER_RESERVED"},
    )

    def assemble(_cfg, overrides):
        if not overrides:
            return drawing
        raise error

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        checked = _checked_circular_drawing(assemble, _default_config(), warn_search_seconds=0.0, map_item="map")
    row = label_fit_row(checked, "smaller_font")
    assert row is not None and row.status == "not_found"
    unfit = _unfit_row_verification(error)
    assert unfit is not None and unfit.outcome == "stop"
    assert unfit.reason == PlainReason("tracks_do_not_fit")


def test_a_failed_render_keeps_the_shown_results_rows_computable() -> None:
    """OV-737: a failure of a mode keeps its last drawn figure; a newer drawn render releases both."""

    from gbdraw.api.diagram import build_circular_diagram

    record = _crossing_record()
    options = _crossing_options()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(record, options=options, _label_fit_warn_search_seconds=0.0)
    searches = LabelFitRowSearches()
    (drawn,) = searches.retain("circular", [drawing])
    failure = LabelsDoNotFitError("labels overlap", report=_report(drawn=False), row_search=_Rows())
    failed = searches.retain_failure("circular", failure)
    assert drawn is not None and failed is not None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert searches.row(drawn, "outward_vertical")["status"] != "not_computed"
    assert searches.row(failed, "combination")["status"] == "found"
    # A failure without rows to compute keeps the drawn figure too.
    done = LabelsDoNotFitError("labels overlap", report=LabelFitReport(records=(), rows=(), drawn=False))
    assert searches.retain_failure("circular", done) is None
    assert searches.row(failed, "combination")["status"] == "not_computed"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert searches.row(drawn, "outward_vertical")["status"] != "not_computed"
        # A newer drawn render of the mode releases the older figure and any failure.
        searches.retain_failure("circular", failure)
        searches.retain("circular", [])
    assert searches.row(drawn, "outward_vertical")["status"] == "not_computed"
    assert searches.failure_payload(failure, "circular")["rowToken"] is None


def test_the_notice_bar_words_are_pythons_with_one_plural_rule() -> None:
    def notice(crossing: int, covering: int) -> list[str]:
        record = LabelFitRecord("NC_1", "circular_outer", 80, 14.0, LabelFit(0, crossing, covering))
        return label_fit_report_payload(LabelFitReport((record,), (), drawn=True), "circular")["notice"]

    assert notice(1, 0) == ["Leader lines cross (1 pair)."]
    assert notice(2, 3) == ["Leader lines cross (2 pairs).", "Labels cover the definition text (3)."]
    assert notice(0, 1) == ["Labels cover the definition text (1)."]


def test_a_web_row_stops_at_its_bound_and_says_no_value_was_found() -> None:
    """LP-D60: the Web searches each row for at most 30 s; a row that finds nothing says so."""

    import time

    from gbdraw.api.diagram import build_circular_diagram

    assert WEB_ROW_SEARCH_SECONDS == 30.0
    record = _crossing_record()
    options = _crossing_options()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(record, options=options, _label_fit_warn_search_seconds=0.0)
    # A bound of 0 s stops the search before its first probe, as 30 s stops a long one.
    searches = LabelFitRowSearches(row_search_seconds=0.0)
    (token,) = searches.retain("circular", [drawing])
    row = searches.row(str(token), "outward_vertical")
    assert row["status"] == "not_computed" and row["values"] is None
    assert row["reason"] == {"kind": "search_bound", "limit": 0.0}

    # The rows a stopped Generate left are searched with the same bound.
    rows = _Rows()
    failure = LabelsDoNotFitError("labels overlap", report=_report(drawn=False), row_search=rows)
    bounded = LabelFitRowSearches()
    failed = bounded.retain_failure("circular", failure)
    before = time.monotonic()
    assert bounded.row(str(failed), "combination")["status"] == "found"
    (deadline,) = rows.deadlines
    assert deadline is not None and before + WEB_ROW_SEARCH_SECONDS <= deadline <= time.monotonic() + WEB_ROW_SEARCH_SECONDS


def test_a_web_row_of_a_figure_whose_first_assembly_took_the_whole_bound_still_verifies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OV-2527: the first full assembly does not estimate a verification. A figure that took 30 s to draw (a slow
    machine, the browser) gets its rows; before, every row stopped at once with "No value found within 30 seconds."
    """

    import time

    import gbdraw.api.diagram as diagram_module
    from gbdraw.api.diagram import build_circular_diagram, label_fit_report

    clock = [0.0]
    monkeypatch.setattr(time, "monotonic", lambda: clock[0])
    real_assemble = diagram_module.assemble_circular_diagram_from_record
    assemblies: list[float] = []

    def timed_assemble(*args: object, **kwargs: object):
        # The first assembly takes the whole Web bound; each verification half a second.
        clock[0] += WEB_ROW_SEARCH_SECONDS if not assemblies else 0.5
        assemblies.append(clock[0])
        return real_assemble(*args, **kwargs)

    monkeypatch.setattr(diagram_module, "assemble_circular_diagram_from_record", timed_assemble)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = build_circular_diagram(
            _crossing_record(), options=_crossing_options(), _label_fit_warn_search_seconds=0.0
        )
    report = label_fit_report(drawing)
    assert report is not None and report.drawn and report.rows
    searches = LabelFitRowSearches()
    (token,) = searches.retain("circular", [drawing])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        row = searches.row(str(token), report.rows[0].kind)

    assert len(assemblies) > 1
    assert row["status"] != "not_computed", row["reason"]


@pytest.mark.parametrize("digit", ("١", "²", "３"))
def test_label_fit_row_token_with_non_ascii_digits_is_stale(digit: str) -> None:
    rows = _Rows()
    error = LabelsDoNotFitError("labels overlap", report=_report(drawn=False), row_search=rows)
    searches = LabelFitRowSearches()
    token = str(searches.retain_failure("circular", error))
    mode, generation, index = token.split(":")
    assert searches.row(token, "combination")["status"] == "found"
    for stale in (f"{mode}:{digit}:{index}", f"{mode}:{generation}:{digit}"):
        row = searches.row(stale, "combination")
        assert row["status"] == "not_computed" and row["reason"]


def test_a_linear_figure_with_region_labels_over_others_opens_the_notice_with_no_rows(tmp_path) -> None:
    """LP-D68, LP-D82: the Linear Result carries a report for the notice bar, with no fix rows."""

    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api import (
        AnnotationOptions,
        AnnotationSet,
        CoordinateSpan,
        LinearDiagramOptions,
        LinearDiagramRequest,
        RegionAnnotation,
    )

    record = SeqRecord(Seq("ACGT" * 1500), id="adjacent", name="adjacent", annotations={"topology": "linear"})
    record.features = [
        SeqFeature(FeatureLocation(2999, 3300, strand=1), type="CDS", qualifiers={"product": ["photosystem II protein D1"]}),
    ]
    regions = (
        RegionAnnotation("r0", CoordinateSpan(None, 1000, 1100), label="First region with a long label", mark="highlight"),
        RegionAnnotation("r1", CoordinateSpan(None, 1120, 1220), label="Second region with a long label", mark="highlight"),
    )
    document = drawing_view(_embedded_document(build_session_document(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(
                config_overrides={"labels.linear.scope": "all", "canvas.linear.track_layout": "above"},
                annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            ),
            output=RenderOutputRequest(output_prefix="fit", formats=("svg",)),
        )
    )))
    searches = LabelFitRowSearches()
    response = _web_render(document, tmp_path, searches)

    (entry,) = response["metadata"]["labelFitReports"]
    report = entry["report"]
    assert report["drawn"] is True and report["regionOverlaps"] == 1
    assert report["notice"] == ["Region labels overlap other labels or lines (1)."]
    (record_payload,) = report["records"]
    assert (record_payload["labelKind"], record_payload["labelKindTitle"], record_payload["counts"]) == (
        "linear_region", "region labels", "1 region label overlap"
    )
    assert report["rows"] == []
