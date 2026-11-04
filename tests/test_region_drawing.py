from __future__ import annotations

import json
from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from gbdraw.annotations.models import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureIdentitySpan,
    RegionAnnotation,
)
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    DepthTrackInput,
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    MaterializedSession,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
    RegionSelection,
    SessionDrawing,
    SessionDrawingSelectionError,
    SessionFormatError,
    build_session_document,
    derive_region_drawing,
    load_session_document,
    materialize_session,
    render_request,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.region_drawing import derive_region_request
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import DiagramRequest
from gbdraw.auto_sizes import DrawnExtent, adapt_explicit_settings, auto_setting_values
from gbdraw.exceptions import ValidationError
from gbdraw.features.overrides import FeatureOverride
from gbdraw.io.record_select import RecordSelector
from tests.utils.session_drawing_view import drawing_view
from tests.utils.svg_compare import assert_nested_approx
from tests.utils.two_mode_session import two_mode_session

REPO = Path(__file__).resolve().parents[1]
VECTORS = json.loads((REPO / "tests" / "fixtures" / "auto_size_vectors.json").read_text(encoding="utf-8"))
MITO = REPO / "tests" / "test_inputs" / "HmmtDNA.gbk"
TWO_RECORDS = REPO / "tests" / "fixtures" / "web_batch_two_records.gb"


def _extent(payload: dict[str, Any]) -> DrawnExtent:
    return DrawnExtent(payload["mode"], tuple(payload["recordLengths"]))


@pytest.fixture(scope="module")
def cfg():
    return apply_config_overrides(None, None)


@pytest.mark.parametrize("case", VECTORS["autoValues"], ids=lambda case: case["name"])
def test_auto_value_vectors(cfg, case: dict[str, Any]) -> None:
    values = auto_setting_values(_extent(case), cfg)
    for setting, expected in case["expected"].items():
        assert values[setting] == pytest.approx(expected, abs=1e-9), setting


@pytest.mark.parametrize("case", VECTORS["adaptations"], ids=lambda case: case["name"])
def test_size_rule_vectors(cfg, case: dict[str, Any]) -> None:
    adaptation = adapt_explicit_settings(
        case["explicit"], source=_extent(case["source"]), target=_extent(case["target"]), cfg=cfg
    )
    assert dict(adaptation.kept) == case["expected"]["kept"]
    reset = [
        {
            "setting": item.setting,
            "explicit": item.explicit,
            "sourceAuto": item.source_auto,
            "targetAuto": item.target_auto,
        }
        for item in adaptation.reset
    ]
    assert_nested_approx(reset, case["expected"]["reset"], abs=1e-9)


@pytest.fixture
def materialize(tmp_path: Path) -> Iterator[Any]:
    """Build a Session from a request and materialize it for the test's lifetime."""
    contexts = []

    def open_session(request: DiagramRequest) -> MaterializedSession:
        context = materialize_session(build_session_document(request), output_directory=tmp_path / "out")
        contexts.append(context)
        return context.__enter__()

    yield open_session
    for context in reversed(contexts):
        context.__exit__(None, None, None)


def _mito_circular(**options: Any) -> CircularDiagramRequest:
    return CircularDiagramRequest(
        records=[RecordInput(GenBankInputSource(MITO), record_key="mito", **options.pop("record", {}))],
        options=CircularDiagramOptions(**options),
    )


@pytest.mark.parametrize("case", VECTORS["regions"], ids=lambda case: case["name"])
def test_region_vectors(materialize, case: dict[str, Any]) -> None:
    materialized = materialize(_mito_circular())
    selection = RegionSelection("mito", case["start"], case["end"])
    expected = case["expected"]
    if "error" in expected:
        message = "cross the origin" if expected["error"] == "crossesOrigin" else "is outside 1..16569"
        with pytest.raises(ValidationError, match=message):
            derive_region_drawing(materialized, [selection], margin=case["margin"])
        return
    region = derive_region_drawing(materialized, [selection], margin=case["margin"])
    [record] = region.drawing.request.records
    if expected["whole"]:
        assert record.region is None
        assert region.drawing.mode == "circular"
    else:
        assert (record.region.start, record.region.end) == (expected["start"], expected["end"])
        assert region.drawing.mode == "linear"


def test_region_of_a_circular_drawing_is_linear_with_adapted_sizes(materialize, tmp_path: Path) -> None:
    source = _mito_circular(
        config_overrides={
            "labels.font_size.short": 16,
            "labels.font_size.long": 16,
            "objects.legends.font_size.short": 18,
            "objects.legends.font_size.long": 18,
            "canvas.circular.track_type": "middle",
            "objects.depth.tick_font_size": 6,
        },
        window=500,
        step=50,
        evalue=1e-9,
    )
    materialized = materialize(source)
    region = derive_region_drawing(
        materialized, [RegionSelection("mito", 12337, 14148)], margin=600, name="ND5 region"
    )

    assert (region.drawing.mode, region.drawing.name, region.drawing.state) == ("linear", "ND5 region", None)
    request = region.drawing.request
    assert isinstance(request, LinearDiagramRequest)
    [record] = request.records
    assert (record.record_key, record.region.start, record.region.end) == ("mito", 11737, 14748)
    assert record.region.reverse_complement is False
    assert {item.setting for item in region.adaptation.reset} == {"labels.font_size", "objects.depth.tick_font_size"}
    assert region.adaptation.kept == {"objects.legends.font_size": 18, "window": 500, "step": 50}
    overrides = request.options.config_overrides
    assert overrides["objects.legends.font_size.short"] == 18
    assert not any(path.startswith(("labels.font_size.s", "labels.font_size.l", "canvas.circular")) for path in overrides)
    assert "objects.depth.tick_font_size" not in overrides
    assert region.dropped == ("Circular-only settings: canvas.circular.track_type.",)
    assert (request.options.window, request.options.step) == (500, 50)
    # Comparison thresholds start from the Linear defaults.
    assert request.options.evalue == pytest.approx(1e-2)

    rendered = render_request(request)
    assert all(path.is_file() for path in rendered.output_paths)
    document = build_session_document(request)
    assert drawing_view(document.to_dict())["renderRequest"]["records"][0]["region"] == {
        "selector": None,
        "start": 11737,
        "end": 14748,
        "reverseComplement": False,
    }


def test_a_region_drawing_is_written_through_the_drawings_api(materialize) -> None:
    materialized = materialize(_mito_circular())
    region = derive_region_drawing(materialized, [RegionSelection("mito", 100, 900)])

    assert (region.drawing.id, region.drawing.name, region.drawing.state) == (None, None, None)
    document = build_session_document(drawings=[region.drawing])
    assert document.drawings == (SessionDrawing("linear", "Linear", "linear", True),)
    # The current Session version names each drawing by its mode.
    named = derive_region_drawing(materialized, [RegionSelection("mito", 100, 900)], name="ND5 region")
    with pytest.raises(SessionFormatError, match="names each drawing by its mode"):
        build_session_document(drawings=[named.drawing])


def test_the_source_drawing_is_selected_by_id_or_name(tmp_path: Path) -> None:
    # Both drawings name their record "record-1": HmmtDNA (16,569 bp) in the
    # Circular drawing and lambda (48,502 bp) in the Linear drawing.
    document = load_session_document(two_mode_session())
    region = RegionSelection("record-1", 20000, 30000)
    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        with pytest.raises(SessionDrawingSelectionError, match="2 drawings"):
            derive_region_drawing(materialized, [region])
        with pytest.raises(ValidationError, match=r"outside 1\.\.16569"):
            derive_region_drawing(materialized, [region], drawing="circular")
        derived = derive_region_drawing(materialized, [region], drawing="Linear")
        [record] = derived.drawing.request.records
        assert record.source.path.name.endswith("NC_001416.gb")
    assert (record.region.start, record.region.end) == (20000, 30000)


def test_adapt_sizes_off_keeps_every_size(materialize) -> None:
    materialized = materialize(
        _mito_circular(config_overrides={"objects.legends.font_size.short": 30, "objects.legends.font_size.long": 30})
    )
    region = derive_region_drawing(
        materialized, [RegionSelection("mito", 1, 40)], mode="circular", adapt_sizes=False
    )
    assert region.adaptation.reset == ()
    assert region.adaptation.kept == {"objects.legends.font_size": 30}
    assert region.drawing.request.options.config_overrides["objects.legends.font_size.long"] == 30


def test_whole_record_keeps_the_mode_orientation_and_display(materialize) -> None:
    materialized = materialize(_mito_circular(record={"presentation": RecordPresentation(reverse_complement=True)}))
    region = derive_region_drawing(materialized, [RegionSelection("mito", 1, 16569)])
    [record] = region.drawing.request.records
    assert region.drawing.mode == "circular"
    assert isinstance(region.drawing.request, CircularDiagramRequest)
    assert record.region is None and record.presentation.reverse_complement is True
    assert region.adaptation.reset == ()


def test_a_crop_keeps_the_orientation_in_source_coordinates(materialize) -> None:
    materialized = materialize(_mito_circular(record={"presentation": RecordPresentation(reverse_complement=True)}))
    [record] = derive_region_drawing(materialized, [RegionSelection("mito", 100, 900)]).drawing.request.records
    assert (record.region.start, record.region.end, record.region.reverse_complement) == (100, 900, True)
    assert record.presentation.reverse_complement is False
    [forward] = derive_region_drawing(
        materialized, [RegionSelection("mito", 100, 900)], reverse_complement=False
    ).drawing.request.records
    assert forward.region.reverse_complement is False


def test_regions_of_several_records_follow_source_order(materialize, tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"TESTB\t{position}\t{position % 7}\n" for position in range(1, 4001)), encoding="utf-8")
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)],
        options=LinearDiagramOptions(depth_tracks=(DepthTrackInput(source=(None, str(depth))),)),
    )
    materialized = materialize(source)
    collection = resolve_request_records(source)
    keys = [item.record_key for item in collection.provenance]
    assert keys == ["pair:1", "pair:2"]
    lengths = [len(record) for record in collection.records]
    region = derive_region_drawing(
        materialized,
        [RegionSelection("pair:2", 1, min(50, lengths[1])), RegionSelection("pair:1", 1, min(50, lengths[0]))],
    )
    records = region.drawing.request.records
    assert [record.record_key for record in records] == ["pair:1", "pair:2"]
    # The source rows' own selectors (the Session stores each record by its ID),
    # each in its crop, as the Web writes a crop.
    assert [record.selector for record in records] == [None, None]
    assert [record.region.record_id for record in records] == ["TESTA", "TESTB"]
    assert all(record.cardinality is RecordCardinality.EXACTLY_ONE for record in records)
    [track] = region.drawing.request.options.depth_tracks
    assert track.source[0] is None and str(track.source[1]).endswith("depth.tsv")

    assert derive_region_drawing(materialized, [RegionSelection("pair:1", 1, 20)]).drawing.request.options.depth_tracks is None


def test_a_row_of_every_record_of_a_file_is_pinned_by_index(tmp_path: Path) -> None:
    document = build_session_document(
        LinearDiagramRequest(
            records=[RecordInput(GenBankInputSource(TWO_RECORDS), selector=RecordSelector("#1", None, 0), record_key="pair")]
        )
    ).to_dict()
    document["drawings"][0]["renderRequest"]["records"][0].update(cardinality="all", selector=None)
    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        region = derive_region_drawing(materialized, [RegionSelection("pair:2", 1, 20)])
    [record] = region.drawing.request.records
    assert (record.record_key, record.selector, record.region.record_index) == ("pair:2", None, 1)


def test_a_crop_of_a_record_of_a_file_is_stated_as_the_web_writes_it() -> None:
    # A card of a file of two records, with a record gap: each crop names its
    # record in the region, and the drawing, whose records take a row each,
    # states no layout (OV-2464).
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)],
        layout=LinearMultiRecordOptions(record_gap_px=40),
    )
    request = derive_region_request(source, [RegionSelection("pair:1", 2, 30), RegionSelection("pair:2", 5, 40)]).request
    assert request.layout is None
    assert [(record.selector, record.region.record_index) for record in request.records] == [(None, 0), (None, 1)]
    stated = build_session_document(request).to_dict()["drawings"][0]["renderRequest"]
    assert stated["layout"] == {}
    assert [(record["selector"], record["region"]["selector"] is not None) for record in stated["records"]] == [
        (None, True),
        (None, True),
    ]


def test_feature_edits_and_annotations_inside_the_region_are_carried(materialize) -> None:
    plain = _mito_circular()
    catalog = resolve_request_records(plain).provenance[0].source_feature_catalog
    assert catalog
    inside = next(feature for feature in catalog if feature.location_parts[0][0] >= 1000 and feature.location_parts[-1][1] <= 2000)
    outside = next(feature for feature in catalog if feature.location_parts[0][0] >= 10000)
    annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                id="marks",
                annotations=(
                    RegionAnnotation("in-feature", FeatureIdentitySpan("mito", inside.biological_feature_id)),
                    RegionAnnotation("out-feature", FeatureIdentitySpan("mito", outside.biological_feature_id)),
                    RegionAnnotation("in-span", CoordinateSpan(None, 1200, 1300)),
                    RegionAnnotation("edge-error", CoordinateSpan(None, 500, 1500, out_of_bounds="error")),
                    RegionAnnotation("wrap", CoordinateSpan(None, 16000, 100, wraps_origin=True)),
                    RegionAnnotation("by-id", CoordinateSpan(RecordSelector("NC_012920.1", "NC_012920.1", None), 1100, 1150)),
                ),
            ),
        )
    )
    source = _mito_circular(
        annotations=annotations,
        feature_overrides=(
            FeatureOverride("mito", inside.biological_feature_id, label_text="kept"),
            FeatureOverride("mito", outside.biological_feature_id, label_text="dropped"),
        ),
    )
    region = derive_region_drawing(materialize(source), [RegionSelection("mito", 1000, 2000)])
    options = region.drawing.request.options
    assert [row.label_text for row in options.feature_overrides] == ["kept"]
    [marks] = options.annotations.sets
    assert [item.id for item in marks.annotations] == ["in-feature", "in-span", "by-id"]
    assert marks.annotations[2].target.record.record_id == "NC_012920.1"
    assert "Annotation 'wrap' (crosses the origin)." in region.dropped
    assert "2 annotation(s) outside the regions." in region.dropped


def test_comparison_tables_are_not_carried(materialize) -> None:
    examples = REPO / "examples"
    source = LinearDiagramRequest(
        records=[
            RecordInput(GenBankInputSource(examples / "MjeNMV.gb"), record_key="a"),
            RecordInput(GenBankInputSource(examples / "MelaMJNV.gb"), record_key="b"),
        ],
        options=LinearDiagramOptions(blast_files=[str(examples / "MjeNMV.MelaMJNV.tblastx.out")]),
    )
    region = derive_region_drawing(
        materialize(source), [RegionSelection("a", 1, 30000), RegionSelection("b", 1, 30000)]
    )
    assert region.drawing.request.options.blast_files is None
    assert region.drawing.request.options.linear_comparisons is None
    assert "The comparison tables (they use the coordinates of the whole records)." in region.dropped


@pytest.mark.parametrize(
    ("selections", "options", "message"),
    [
        ([], {}, "at least one region"),
        ([RegionSelection("other", 1, 10)], {}, "no record 'other'"),
        ([RegionSelection("mito", 1, 10), RegionSelection("mito", 20, 30)], {}, "more than one region"),
        ([RegionSelection("mito", 1, 10)], {"margin": -1}, "margin"),
    ],
)
def test_invalid_regions_are_refused(materialize, selections, options, message) -> None:
    with pytest.raises(ValidationError, match=message):
        derive_region_drawing(materialize(_mito_circular()), selections, **options)


def test_a_circular_region_drawing_takes_one_record(materialize) -> None:
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)]
    )
    with pytest.raises(ValidationError, match="exactly one record"):
        derive_region_drawing(
            materialize(source),
            [RegionSelection("pair:1", 1, 10), RegionSelection("pair:2", 1, 10)],
            mode="circular",
        )


def test_circular_track_widths_follow_the_rule(materialize) -> None:
    from gbdraw.api import CircularRequestTrackOptions, CircularTrackSlot, ScalarSpec

    slots = (
        CircularTrackSlot("features", "features", width=ScalarSpec(30.0, "px")),
        CircularTrackSlot("gc_content", "gc_content", width=ScalarSpec(60.0, "px")),
    )
    source = _mito_circular(tracks=CircularRequestTrackOptions(circular_track_slots=slots))
    region = derive_region_drawing(materialize(source), [RegionSelection("mito", 1, 16000)], mode="circular")
    assert region.adaptation.reset == ()
    kept_slots = region.drawing.request.options.tracks.circular_track_slots
    # The request encoder writes a renderer by the planner's name, so the alias
    # gc_content reads back as dinucleotide_content (OV-653); the widths stay.
    assert [replace(slot) for slot in kept_slots] == [
        slots[0],
        replace(slots[1], renderer="dinucleotide_content"),
    ]


# The request-level core and the Worker helper of New drawing from a region.


def test_the_drawing_derivation_is_the_request_core(materialize) -> None:
    from gbdraw.session import session_to_request

    source = _mito_circular(config_overrides={"canvas.circular.track_type": "middle"})
    materialized = materialize(source)
    selection = [RegionSelection("mito", 100, 900)]
    drawing = derive_region_drawing(materialized, selection)
    core = derive_region_request(session_to_request(materialized), selection)
    assert (core.request, core.mode, core.adaptation) == (drawing.drawing.request, "linear", drawing.adaptation)
    assert [(item.kind, item.items) for item in core.not_carried] == [
        ("OTHER_MODE_SETTINGS", ("canvas.circular.track_type",))
    ]
    assert drawing.dropped == ("Circular-only settings: canvas.circular.track_type.",)


def test_an_origin_crossing_region_has_its_own_reason(materialize) -> None:
    from gbdraw.api.region_drawing import StartAfterEndRegion

    with pytest.raises(StartAfterEndRegion) as refusal:
        derive_region_drawing(materialize(_mito_circular()), [RegionSelection("mito", 16000, 200)])
    assert refusal.value.diagnostic == {"code": "REGION_INVALID", "field": "region", "reason": "ORIGIN_CROSSING"}
    assert (refusal.value.record_key, refusal.value.start, refusal.value.end) == ("mito", 16000, 200)


def test_a_start_after_the_end_of_a_linear_record_is_worded_without_the_origin(tmp_path: Path) -> None:
    """A linear record has no origin to cross: its refusal names the order only."""
    from gbdraw.api import RecordDisplayOptions
    from gbdraw.api.region_drawing import StartAfterEndRegion

    linear = _mito_circular(record={"display": RecordDisplayOptions(is_circular=False)})
    with pytest.raises(StartAfterEndRegion) as refusal:
        derive_region_request(linear, [RegionSelection("mito", 16000, 200)])
    assert refusal.value.diagnostic == {"code": "REGION_INVALID", "field": "region", "reason": "START_AFTER_END"}
    assert "origin" not in str(refusal.value)
    reply = _web_region(linear, [{"recordKey": "mito", "start": 16000, "end": 200}], tmp_path)[2]
    assert reply == {"status": "refused", "reason": "START_AFTER_END", "recordKey": "mito", "start": 16000, "end": 200}


def _web_source(request: DiagramRequest) -> tuple[dict[str, Any], dict[str, str]]:
    from gbdraw.session_request_codec import encode_canonical_request

    encoded = encode_canonical_request(request)
    return encoded.payload, {resource.resource_id: str(resource.source_path) for resource in encoded.resources}


def _web_region(request: DiagramRequest, regions: list[dict[str, Any]], tmp_path: Path, *, adapt_sizes: bool = True):
    from gbdraw.web_support.region_drawing import derive_web_region_drawing_json

    payload, paths = _web_source(request)
    reply = json.loads(derive_web_region_drawing_json(
        json.dumps(payload), json.dumps(paths), str(tmp_path / "out"), json.dumps(regions), json.dumps(adapt_sizes)
    ))
    return payload, paths, reply


def test_the_web_helper_derives_one_record_of_a_circular_drawing(tmp_path: Path) -> None:
    from gbdraw.session_request_codec import decode_canonical_request

    source = _mito_circular(config_overrides={"labels.font_size.short": 16, "labels.font_size.long": 16})
    payload, paths, reply = _web_region(source, [{"recordKey": "mito", "start": 100, "end": 900}], tmp_path)
    assert (reply["status"], reply["mode"]) == ("derived", "linear")
    request = reply["renderRequest"]
    # The new drawing names the source's resource: the same file, no new upload.
    assert request["records"][0]["source"] == payload["records"][0]["source"]
    assert request["records"][0]["region"]["start"] == 100
    assert [item["setting"] for item in reply["returnedToAuto"]] == ["labels.font_size"]
    assert reply["returnedToAuto"][0]["explicit"] == 16
    decoded = decode_canonical_request(request, resource_paths=paths, output_directory=tmp_path / "out")
    assert isinstance(decoded, LinearDiagramRequest)
    kept = _web_region(
        source, [{"recordKey": "mito", "start": 100, "end": 900}], tmp_path, adapt_sizes=False
    )[2]
    assert kept["returnedToAuto"] == []


def test_the_web_helper_keeps_a_whole_circular_record_circular(tmp_path: Path) -> None:
    reply = _web_region(_mito_circular(), [{"recordKey": "mito", "start": 1, "end": 16569}], tmp_path)[2]
    assert (reply["status"], reply["mode"], reply["renderRequest"]["mode"]) == ("derived", "circular", "circular")


def test_the_web_helper_draws_two_records_in_one_linear_drawing(tmp_path: Path) -> None:
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)],
    )
    payload, _paths, reply = _web_region(
        source, [{"recordKey": "pair:1", "start": 1, "end": 40}, {"recordKey": "pair:2", "start": 1, "end": 40}], tmp_path
    )
    assert (reply["status"], reply["mode"]) == ("derived", "linear")
    records = reply["renderRequest"]["records"]
    assert [record["recordKey"] for record in records] == ["pair:1", "pair:2"]
    # Each input keeps a resource of its own (the codec's rule); each reads the source's file.
    source_id = payload["records"][0]["source"]["resourceId"]
    assert {reply["resourceSources"][record["source"]["resourceId"]] for record in records} == {source_id}


def test_the_web_helper_refuses_an_origin_crossing_region(tmp_path: Path) -> None:
    reply = _web_region(_mito_circular(), [{"recordKey": "mito", "start": 16000, "end": 200}], tmp_path)[2]
    assert reply == {"status": "refused", "reason": "ORIGIN_CROSSING", "recordKey": "mito", "start": 16000, "end": 200}


def _origin_record(tmp_path: Path) -> Path:
    """A 1000 bp circular record with a CDS at 101..400 and a CDS across the origin."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ACGT" * 250), id="orig", name="orig", description="origin test")
    record.annotations.update({"molecule_type": "DNA", "topology": "circular"})
    record.features = [
        SeqFeature(FeatureLocation(100, 400, strand=1), type="CDS", qualifiers={"locus_tag": ["inside"]}),
        SeqFeature(
            CompoundLocation([FeatureLocation(900, 1000, strand=1), FeatureLocation(0, 50, strand=1)]),
            type="CDS", qualifiers={"locus_tag": ["across"]},
        ),
    ]
    path = tmp_path / "origin.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _feature_id(request: DiagramRequest, locus_tag: str) -> str:
    catalog = resolve_request_records(request).provenance[0].source_feature_catalog
    return next(
        feature.biological_feature_id for feature in catalog
        if dict(feature.qualifiers).get("locus_tag") == (locus_tag,)
    )


def test_a_feature_region_is_the_span_of_its_location(tmp_path: Path) -> None:
    from gbdraw.api.region_drawing import FeatureRegionSelection

    source = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_origin_record(tmp_path)), record_key="o")])
    by_feature = derive_region_request(source, [FeatureRegionSelection("o", _feature_id(source, "inside"))])
    by_coordinates = derive_region_request(source, [RegionSelection("o", 101, 400)])
    assert by_feature == by_coordinates
    assert by_feature.mode == "linear"


def test_a_feature_across_the_origin_is_refused(tmp_path: Path) -> None:
    from gbdraw.api.region_drawing import FeatureRegionSelection, StartAfterEndRegion

    source = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_origin_record(tmp_path)), record_key="o")])
    with pytest.raises(StartAfterEndRegion) as refusal:
        derive_region_request(source, [FeatureRegionSelection("o", _feature_id(source, "across"))])
    assert (refusal.value.record_key, refusal.value.start, refusal.value.end) == ("o", 901, 50)


def test_an_unknown_feature_is_refused(tmp_path: Path) -> None:
    from gbdraw.api.region_drawing import FeatureRegionSelection

    source = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_origin_record(tmp_path)), record_key="o")])
    with pytest.raises(ValidationError, match="no feature 'missing'"):
        derive_region_request(source, [FeatureRegionSelection("o", "missing")])


def test_the_web_helper_takes_a_feature_region(tmp_path: Path) -> None:
    source = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_origin_record(tmp_path)), record_key="o")])
    inside, across = _feature_id(source, "inside"), _feature_id(source, "across")
    reply = _web_region(source, [{"recordKey": "o", "featureId": inside}], tmp_path)[2]
    assert reply["status"] == "derived"
    assert reply["renderRequest"]["records"][0]["region"]["start"] == 101
    refused = _web_region(source, [{"recordKey": "o", "featureId": across}], tmp_path)[2]
    assert refused == {"status": "refused", "reason": "ORIGIN_CROSSING", "recordKey": "o", "start": 901, "end": 50}
