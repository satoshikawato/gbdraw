"""What the Web reads of a wrapped figure: each refusal's reason, the export limit's facts, and the rows drawn."""

from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO
import pytest

from gbdraw.api.options import LinearDiagramOptions, LinearRecordTranslation
from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput, RenderOutputRequest
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ExportError, GbdrawError
from gbdraw.io.regions import parse_region_spec
from gbdraw.layout.linear_wrap import LinearWrapSummary, refuse_feature_track_on_the_axis, refuse_undrawable_wrap
from gbdraw.render.export_limits import check_export_limit
from gbdraw.web_support.error_adapter import serialize_web_error
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from gbdraw.session import build_session_document

WRAP = "canvas.linear.wrap_length"
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
SORANGIUM = LinearWrapSummary(50_000, 50_000, (13_033_779,))
SORANGIUM_SIZE = (2_000.0, 45_967.0)


def _config(**overrides: object) -> GbdrawConfig:
    return GbdrawConfig.from_dict(modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), overrides))


def _caught(call) -> GbdrawError:
    with pytest.raises(GbdrawError) as caught:
        call()
    return caught.value


def _context(error: BaseException, operation: str = "generate") -> dict:
    payload = serialize_web_error(error, operation=operation, stage="render")
    return {"code": payload["code"], **payload["context"]}


# OV-2465: the adapter kept only `WRAP`; every other wrap refusal reached the Web without its reason.
@pytest.mark.parametrize(
    ("refuse", "code", "reason"),
    [
        (lambda: refuse_undrawable_wrap(_config(), [0], ("blast_files",), row_length=10_000), "COMPARISON_INPUT", "WRAP"),
        (lambda: refuse_undrawable_wrap(_config(**{"canvas.linear.track_layout": "middle"}), [0], (), row_length=10_000),
         "INPUT_INVALID", "WRAP_TRACK_LAYOUT"),
        (lambda: refuse_undrawable_wrap(_config(**{"canvas.linear.track_layout": "above", "objects.scale.style": "bar"}), [0], (), row_length=10_000),
         "INPUT_INVALID", "WRAP_SCALE_STYLE"),
        (lambda: refuse_undrawable_wrap(
            _config(**{"canvas.linear.track_layout": "above", "objects.scale.style": "ruler", "canvas.linear.normalize_length": True}),
            [0], (), row_length=10_000), "INPUT_INVALID", "WRAP_NORMALIZE_LENGTH"),
        (lambda: refuse_undrawable_wrap(_config(**{"canvas.linear.track_layout": "above", "objects.scale.style": "ruler"}), [0, 0], (), row_length=10_000),
         "INPUT_INVALID", "WRAP_GRID_ROW"),
        (lambda: refuse_undrawable_wrap(
            _config(**{"canvas.linear.track_layout": "above", "objects.scale.style": "ruler"}), [0], (), row_length=10_000, anchor_feature=True),
         "INPUT_INVALID", "WRAP_ANCHOR_FEATURE"),
        (lambda: refuse_undrawable_wrap(
            _config(**{"canvas.linear.track_layout": "above", "objects.scale.style": "ruler", WRAP: 10_000}), [0], (),
            row_length=10_000, records=[SeqIO.read(LAMBDA, "genbank")],
            translations=[LinearRecordTranslation(record_key="record-1", x=5.0)]),
         "INPUT_INVALID", "WRAP_RECORD_POSITION"),
        # Q5 A: a Similarity alignment whose comparison inputs were removed moves records; the Web names Reset record
        # positions for it, not the comparison inputs it no longer has.
        (lambda: refuse_undrawable_wrap(_config(), [0], ("similarity_alignment",), row_length=10_000),
         "INPUT_INVALID", "WRAP_RECORD_POSITION"),
        (lambda: refuse_undrawable_wrap(_config(), [0], ("blast_files", "similarity_alignment"), row_length=10_000),
         "COMPARISON_INPUT", "WRAP"),
        (lambda: refuse_feature_track_on_the_axis("middle"), "INPUT_INVALID", "WRAP_TRACK_LAYOUT"),
    ],
)
def test_each_wrap_refusal_keeps_its_reason_on_the_web(refuse, code: str, reason: str) -> None:
    assert _context(_caught(refuse)) == {"code": code, "field": "wrap_length", "reason": reason}


def _export_refusal(fmt: str, size: tuple[float, float], wrap: LinearWrapSummary | None, dpi: float | None = None):
    error = _caught(lambda: check_export_limit(fmt, *size, wrap=wrap, dpi=dpi))
    assert isinstance(error, ExportError)
    return error


def test_a_png_refusal_names_its_limit_size_row_length_and_dpi() -> None:
    error = _export_refusal("png", SORANGIUM_SIZE, SORANGIUM, dpi=300)
    # The CLI text is unchanged (tests/test_export_limits.py pins it).
    assert "DPI of 68 or lower" in str(error) and "--wrap_length 500kb" in str(error)
    assert error.diagnostic == {
        "code": "EXPORT_LIMIT", "reason": "PNG_SIZE", "limit": 32_767, "width": 6_250, "height": 143_647,
        "rowLength": 500_000, "dpi": 68,
    }
    assert _context(error, "export-png") == {
        "code": "EXPORT_LIMIT", "reason": "PNG_SIZE", "limit": 32_767, "width": 6_250, "height": 143_647,
        "rowLength": 500_000, "dpi": 68,
    }


def test_a_wrapped_pdf_refusal_names_a_row_length_and_no_dpi() -> None:
    error = _export_refusal("pdf", SORANGIUM_SIZE, SORANGIUM, dpi=300)
    assert _context(error, "export-pdf") == {
        "code": "EXPORT_LIMIT", "reason": "PDF_SIZE", "limit": 14_400, "width": 1_500, "height": 34_475,
        "rowLength": 200_000,
    }


def test_a_refusal_with_no_advice_names_only_the_limit_and_size() -> None:
    many_records = LinearWrapSummary(50_000, 50_000, (60_000,) * 200)
    assert _context(_export_refusal("png", (2_000.0, 80_000.0), many_records), "export-png") == {
        "code": "EXPORT_LIMIT", "reason": "PNG_SIZE", "limit": 32_767, "width": 2_000, "height": 80_000,
    }
    assert _context(_export_refusal("png", (33_000.0, 100.0), None), "export-png") == {
        "code": "EXPORT_LIMIT", "reason": "PNG_SIZE", "limit": 32_767, "width": 33_000, "height": 100,
    }


@pytest.fixture(scope="module")
def web_helpers() -> dict:
    source = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/app/python-helpers.js").read_text()
    namespace: dict = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


def _web_export_check(helpers: dict, fmt: str, size: tuple[float, float], wrap: LinearWrapSummary | None,
                      dpi: float | None) -> dict:
    # What the Web sends (D-55, D-70): the committed request's length and each record's source length and region.
    request = {"format": fmt, "width": size[0], "height": size[1], "dpi": dpi,
               "linearWrap": None if wrap is None
               else {"requested": wrap.requested, "records": _whole(wrap.record_lengths)}}
    return json.loads(helpers["call_web_json_helper"]("check_web_export_limit_json", json.dumps(request)))


def test_the_web_export_check_reads_the_wrap_and_refuses_as_the_cli_does(web_helpers: dict) -> None:
    # The Result's wrap length and records decide the rows and the advice; the PNG scales by its DPI.
    png = _web_export_check(web_helpers, "png", SORANGIUM_SIZE, SORANGIUM, 300)["error"]
    assert (png["code"], png["operation"]) == ("EXPORT_LIMIT", "checkExportLimit")
    assert png["context"] == {"reason": "PNG_SIZE", "limit": 32_767, "width": 6_250, "height": 143_647,
                              "rowLength": 500_000, "dpi": 68}
    pdf = _web_export_check(web_helpers, "pdf", SORANGIUM_SIZE, SORANGIUM, None)["error"]
    assert pdf["context"] == {"reason": "PDF_SIZE", "limit": 14_400, "width": 1_500, "height": 34_475,
                              "rowLength": 200_000}
    # An unwrapped PDF stays as it was, also when a wrap length cut no record; a PNG that fits at its DPI passes.
    assert _web_export_check(web_helpers, "pdf", SORANGIUM_SIZE, None, None) == {"ok": True}
    uncut = LinearWrapSummary(20_000_000, 20_000_000, (13_033_779,))
    assert _web_export_check(web_helpers, "pdf", SORANGIUM_SIZE, uncut, None) == {"ok": True}
    assert _web_export_check(web_helpers, "png", SORANGIUM_SIZE, SORANGIUM, 68) == {"ok": True}
    # A PNG names its DPI; an unknown format is not checked as another.
    assert "error" in _web_export_check(web_helpers, "png", SORANGIUM_SIZE, None, None)
    assert "error" in _web_export_check(web_helpers, "eps", SORANGIUM_SIZE, None, None)


def _whole(record_lengths: tuple[int, ...] | list[int]) -> list[dict]:
    return [{"sourceLength": length, "region": None} for length in record_lengths]


def _web_export_leaf(helpers: dict, leaf: object, record_lengths: list[int]) -> dict:
    request = {"format": "pdf", "width": SORANGIUM_SIZE[0], "height": SORANGIUM_SIZE[1], "dpi": None,
               "linearWrap": {"requested": leaf, "records": _whole(record_lengths)}}
    return json.loads(helpers["call_web_json_helper"]("check_web_export_limit_json", json.dumps(request)))


def test_the_web_export_check_reads_the_stored_wrap_leaf_with_the_config_validator(web_helpers: dict) -> None:
    # Review 3 F1 (D-52): a loaded Result's committed request goes as stored, so a CLI or Python API Session's
    # "5kb" or "Auto" is checked as the CLI checks it; "none" and null are off.
    for leaf in ("50kb", "Auto", " 50000 ", 50_000):
        refusal = _web_export_leaf(web_helpers, leaf, [13_033_779])["error"]
        assert refusal["context"]["reason"] == "PDF_SIZE", leaf
    for leaf in ("none", None):
        assert _web_export_leaf(web_helpers, leaf, [13_033_779]) == {"ok": True}, leaf
    # A leaf the validator refuses is refused with its diagnostic, as a Generate of it would be.
    refused = _web_export_leaf(web_helpers, "50 parsecs", [13_033_779])["error"]
    assert refused["code"] == "INPUT_INVALID"
    assert refused["context"]["configPath"] == WRAP


def _web_export_records(helpers: dict, records: list[dict], requested: object = 50_000) -> dict:
    request = {"format": "pdf", "width": SORANGIUM_SIZE[0], "height": SORANGIUM_SIZE[1], "dpi": None,
               "linearWrap": {"requested": requested, "records": records}}
    return json.loads(helpers["call_web_json_helper"]("check_web_export_limit_json", json.dumps(request)))


def _cli_refusal(record_lengths: tuple[int, ...]) -> dict:
    with pytest.raises(ExportError) as caught:
        check_export_limit("pdf", *SORANGIUM_SIZE, wrap=LinearWrapSummary(50_000, 50_000, record_lengths))
    return serialize_web_error(caught.value, operation="checkExportLimit", stage="render")["context"]


def test_the_web_export_check_crops_each_record_by_the_drawings_rule(web_helpers: dict) -> None:
    # D-59 residual, D-70: a loaded Result's cropped record is known by its source length and region; Python crops
    # it with the drawing's own function, so the refusal and its advice are the CLI's for the drawn bases.
    region = {"sourceLength": 13_033_779, "region": {"start": 1_000_001, "end": 6_000_000}}
    refusal = _web_export_records(web_helpers, [region])["error"]
    assert refusal["context"] == _cli_refusal((5_000_000,))
    # An end past the record is drawn to the record end, as the drawing clamps it.
    past_end = {"sourceLength": 6_000_000, "region": {"start": 1_000_001, "end": 9_000_000}}
    assert _web_export_records(web_helpers, [past_end])["error"]["context"] == _cli_refusal((5_000_000,))
    # A crop that cuts nothing leaves the PDF as an unwrapped one, unchecked; its source record would be cut.
    short = {"sourceLength": 13_033_779, "region": {"start": 1, "end": 40_000}}
    assert _web_export_records(web_helpers, [short]) == {"ok": True}
    # A region the drawing cannot crop stops the export with its region diagnostic, as a Generate of it would.
    beyond = _web_export_records(web_helpers, [{"sourceLength": 60, "region": {"start": 1_000, "end": 2_000}}])["error"]
    assert (beyond["code"], beyond["context"]) == ("REGION_INVALID", {"field": "region", "reason": "RECORD_BOUNDS"})


def _summaries(tmp_path: Path, **overrides: object) -> list[dict]:
    document = build_session_document(
        LinearDiagramRequest(
            records=(RecordInput(GenBankInputSource(str(LAMBDA))),),
            options=LinearDiagramOptions(config_overrides=overrides),
            output=RenderOutputRequest(output_prefix="web", formats=("svg",)),
        )
    ).to_dict()
    metadata = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "workspace"
    )["metadata"]
    # D-31: the rows ride each Result's input summary; the reader-less list is gone.
    assert "linearWrap" not in metadata
    return metadata["inputSummaries"]


def test_each_input_summary_names_the_rows_its_result_drew(tmp_path: Path) -> None:
    [wrapped] = _summaries(tmp_path, **{WRAP: "auto"})
    assert wrapped["linearWrap"] == {
        "rowLength": 10_000,
        "rowCount": 5,
        "requested": "auto",
        "recordRows": [5],
        "drawnAxisLength": 10_000,
    }
    assert wrapped["records"] == [{"sourceLength": 48_502, "region": None}]
    [plain] = _summaries(tmp_path, **{WRAP: 100_000})
    assert plain["linearWrap"] is None
    # Nothing cut: the records still name what the Result drew (review 2 F2).
    assert plain["records"] == [{"sourceLength": 48_502, "region": None}]


def test_a_cropped_record_names_its_source_length_and_the_region_asked_for(tmp_path: Path) -> None:
    # Review 3 F2 (D-53): the Web tells a cropped record apart by its source length and region, as its card
    # states them; the region is the one the request asked for, also when its end runs past the record.
    def summary(start: int, end: int) -> dict:
        document = build_session_document(
            LinearDiagramRequest(
                records=(RecordInput(GenBankInputSource(str(LAMBDA)), region=parse_region_spec("NC_001416.1:5001-35000")),),
                options=LinearDiagramOptions(config_overrides={WRAP: "auto"}),
                output=RenderOutputRequest(output_prefix="web", formats=("svg",)),
            )
        ).to_dict()
        request = document["drawings"][0]["renderRequest"]
        # The region as a card states it: the Web sends the fields as typed (a Session writer clamps its own).
        request["records"][0]["region"].update(start=start, end=end)
        [entry] = render_embedded_canonical_web_request(
            request, resources=document["resources"], workspace=tmp_path / f"{start}-{end}",
        )["metadata"]["inputSummaries"]
        return entry

    cropped = summary(5_001, 35_000)
    assert cropped["records"] == [{"sourceLength": 48_502, "region": {"start": 5_001, "end": 35_000}}]
    assert cropped["linearWrap"]["recordRows"] == [len(range(0, 30_000, cropped["linearWrap"]["rowLength"]))]
    past_end = summary(1, 60_000)
    assert past_end["records"] == [{"sourceLength": 48_502, "region": {"start": 1, "end": 60_000}}]
