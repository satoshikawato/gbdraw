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
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import DiagramRequest
from gbdraw.auto_sizes import DrawnExtent, adapt_explicit_settings, auto_setting_values
from gbdraw.exceptions import ValidationError
from gbdraw.features.overrides import FeatureOverride
from gbdraw.io.record_select import RecordSelector
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
    assert reset == pytest.approx(case["expected"]["reset"], abs=1e-9)


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
    assert document.to_dict()["renderRequest"]["records"][0]["region"] == {
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
    # The source rows' own selectors (the Session stores each record by its ID).
    assert [record.selector.record_id for record in records] == ["TESTA", "TESTB"]
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
    document["renderRequest"]["records"][0].update(cardinality="all", selector=None)
    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        region = derive_region_drawing(materialized, [RegionSelection("pair:2", 1, 20)])
    [record] = region.drawing.request.records
    assert (record.record_key, record.selector.record_index) == ("pair:2", 1)


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
    assert [replace(slot) for slot in kept_slots] == list(slots)
