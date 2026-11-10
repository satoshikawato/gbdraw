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
from gbdraw.auto_sizes import DrawnExtent, adapt_explicit_settings, auto_setting_values, record_auto_values
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
    if "recordValues" in case:
        # A setting is the drawing's Auto value only when every record takes it.
        for setting in case["recordValues"][0]:
            shared = {record[setting] for record in case["recordValues"]}
            assert values.get(setting) == (shared.pop() if len(shared) == 1 else None), setting
        records = record_auto_values(_extent(case), cfg)
        assert [dict(record) for record in records] == case["recordValues"]


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
        with pytest.raises(ValidationError, match="is outside 1..16569"):
            derive_region_drawing(materialized, [selection], margin=case["margin"])
        return
    region = derive_region_drawing(materialized, [selection], margin=case["margin"])
    [record] = region.drawing.request.records
    if expected["whole"]:
        assert record.windows == ()
        assert region.drawing.mode == "circular"
    else:
        assert (record.windows[0].start, record.windows[0].end) == (expected["start"], expected["end"])
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
    assert (record.record_key, record.windows[0].start, record.windows[0].end) == ("mito", 11737, 14748)
    assert record.windows[0].reverse_complement is False
    # The Auto GC window of the 3,012 bp region is a tenth of its length, not the mitogenome's 1 kb: the hand-set
    # window and step return to Auto, as every size whose Auto differs does.
    assert {item.setting for item in region.adaptation.reset} == {
        "labels.font_size", "objects.depth.tick_font_size", "window", "step"
    }
    assert region.adaptation.kept == {"objects.legends.font_size": 18}
    overrides = request.options.config_overrides
    assert overrides["objects.legends.font_size.short"] == 18
    assert not any(path.startswith(("labels.font_size.s", "labels.font_size.l", "canvas.circular")) for path in overrides)
    assert "objects.depth.tick_font_size" not in overrides
    assert region.dropped == ("Circular-only settings: canvas.circular.track_type.",)
    assert (request.options.window, request.options.step) == (None, None)
    # Comparison thresholds start from the Linear defaults.
    assert request.options.evalue == pytest.approx(1e-2)

    rendered = render_request(request)
    assert all(path.is_file() for path in rendered.output_paths)
    document = build_session_document(request)
    assert drawing_view(document.to_dict())["renderRequest"]["records"][0]["windows"] == [{
        "start": 11737,
        "end": 14748,
        "reverseComplement": False,
        "column": None,
    }]


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
    assert (record.windows[0].start, record.windows[0].end) == (20000, 30000)


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
    assert record.windows == () and record.presentation.reverse_complement is True
    assert region.adaptation.reset == ()


def test_a_crop_keeps_the_orientation_in_source_coordinates(materialize) -> None:
    materialized = materialize(_mito_circular(record={"presentation": RecordPresentation(reverse_complement=True)}))
    [record] = derive_region_drawing(materialized, [RegionSelection("mito", 100, 900)]).drawing.request.records
    assert (record.windows[0].start, record.windows[0].end, record.windows[0].reverse_complement) == (100, 900, True)
    assert record.presentation.reverse_complement is False
    [forward] = derive_region_drawing(
        materialized, [RegionSelection("mito", 100, 900)], reverse_complement=False
    ).drawing.request.records
    assert forward.windows[0].reverse_complement is False


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


def test_a_region_carries_compatibility_depth_inputs_in_their_canonical_form(tmp_path: Path) -> None:
    # OV-2934: per-record depth_track_files followed the request unchanged, a row per source record, so the new
    # drawing of one record was refused at Generate. The region carries them as depth_tracks, as a Session stores them.
    first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"
    first.write_text("".join(f"TESTA\t{position}\t{position % 5}\n" for position in range(1, 3001)))
    second.write_text("".join(f"TESTB\t{position}\t{position % 7}\n" for position in range(1, 3001)))
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)],
        options=LinearDiagramOptions(depth_track_files=[[str(first)], [str(second)]], depth_track_heights=[12]),
    )
    region = derive_region_request(source, [RegionSelection("pair:2", 1, 50)])
    options = region.request.options
    assert options.depth_track_files is None and options.depth_track_heights is None
    [track] = options.depth_tracks
    assert (track.source, track.height) == ((str(second),), 12)
    output = replace(region.request.output, output_directory=tmp_path / "out")
    assert all(path.is_file() for path in render_request(replace(region.request, output=output)).output_paths)


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
    assert (record.record_key, record.selector.record_index) == ("pair:2", 1)


def test_a_crop_of_a_record_of_a_file_is_stated_as_the_web_writes_it() -> None:
    # A card of a file of two records, with a record gap: each crop names its
    # record in its selector, and the drawing, whose records take a row each,
    # states no layout (OV-2464).
    source = LinearDiagramRequest(
        records=[RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair", cardinality=RecordCardinality.ALL)],
        layout=LinearMultiRecordOptions(record_gap_px=40),
    )
    request = derive_region_request(source, [RegionSelection("pair:1", 2, 30), RegionSelection("pair:2", 5, 40)]).request
    assert request.layout is None
    assert [record.selector.record_index for record in request.records] == [0, 1]
    stated = build_session_document(request).to_dict()["drawings"][0]["renderRequest"]
    assert stated["layout"] == {}
    assert [(record["selector"] is not None, len(record["windows"])) for record in stated["records"]] == [
        (True, 1),
        (True, 1),
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


def test_a_losat_search_pairs_the_windows_of_the_drawn_records() -> None:
    from gbdraw.api import LosatSearchOptions

    examples = REPO / "examples"
    source = LinearDiagramRequest(
        records=[
            RecordInput(GenBankInputSource(examples / "MjeNMV.gb"), record_key="a"),
            RecordInput(GenBankInputSource(examples / "MelaMJNV.gb"), record_key="b"),
        ],
        options=LinearDiagramOptions(losat_search=LosatSearchOptions(program="tlosatx", pairs=[(0, 1)])),
    )
    request = derive_region_request(
        source, [RegionSelection("a", 1, 3000), RegionSelection("a", 9001, 12000), RegionSelection("b", 1, 3000)]
    ).request
    # Record a is drawn in two windows (units 0 and 1); each is searched against record b (unit 2).
    assert list(request.options.losat_search.pairs) == [(0, 2), (1, 2)]


@pytest.mark.parametrize(
    ("selections", "options", "message"),
    [
        ([], {}, "at least one region"),
        ([RegionSelection("other", 1, 10)], {}, "no record 'other'"),
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


def _circular_records(tmp_path: Path, big: int, small: int) -> list[RecordInput]:
    """Two circular records, BIG and SMALL (record keys ``pair:1`` and ``pair:2``), in one GenBank file."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    path = tmp_path / f"records_{big}_{small}.gb"
    SeqIO.write(
        [
            SeqRecord(
                Seq("ACGTTGCA" * (length // 8)), id=name, name=name, description=name,
                annotations={"molecule_type": "DNA", "topology": "circular"},
            )
            for name, length in (("BIG", big), ("SMALL", small))
        ],
        path,
        "genbank",
    )
    return [RecordInput(GenBankInputSource(path), record_key="pair", cardinality=RecordCardinality.ALL)]


def _circular_canvas_of_50_and_5_kb(tmp_path: Path, **options: Any) -> CircularDiagramRequest:
    """A Circular drawing of two circular records, 50 kb and 5 kb: each draws its GC content and ticks at its own
    length (Auto window 1,000 and 500 bp, scale interval 5,000 and 1,000 bp)."""
    return CircularDiagramRequest(
        records=_circular_records(tmp_path, 50_000, 5_000), options=CircularDiagramOptions(**options)
    )


def test_a_circular_record_compares_its_own_auto_window_and_scale_interval(tmp_path: Path) -> None:
    # Each record of a Circular drawing draws at its own length: the 5 kb record with Auto window 500 bp and scale
    # interval 1,000 bp, not the 50 kb record's 1,000 bp and 5,000 bp. A drawing of the whole 5 kb record has the same
    # Auto values, so its hand-set values stay (OV-3231: the scale interval follows the same rule).
    source = _circular_canvas_of_50_and_5_kb(
        tmp_path, window=300, step=30, config_overrides={"objects.scale.interval": 2_000}
    )
    region = derive_region_request(source, [RegionSelection("pair:2", 1, 5_000)])
    assert region.mode == "circular"
    assert region.adaptation.reset == ()
    assert {key: region.adaptation.kept[key] for key in ("window", "step", "objects.scale.interval")} == {
        "window": 300, "step": 30, "objects.scale.interval": 2_000,
    }
    assert (region.request.options.window, region.request.options.step) == (300, 30)

    # The 50 kb record's own Auto values change at 4 kb: they return to Auto, the 4 kb record's values shown.
    region = derive_region_request(source, [RegionSelection("pair:1", 1, 4_000)], mode="circular")
    assert {item.setting: (item.source_auto, item.target_auto) for item in region.adaptation.reset} == {
        "window": (1_000, 400), "step": (100, 40), "objects.scale.interval": (5_000, 1_000),
    }


def _gene_density() -> tuple[Any, ...]:
    from gbdraw.analysis.feature_density import FeatureDensityTrack
    from gbdraw.features.query import FeatureQuery

    return (FeatureDensityTrack("d", FeatureQuery("gene")),)


def test_a_circular_region_compares_the_feature_density_window_of_the_canvas(tmp_path: Path) -> None:
    # The records of one Circular canvas share one feature density window, the longest record's (1,000 bp); the
    # whole 5 kb record drawn alone takes 500 bp, so a hand-set window and step the density track reads return to
    # Auto, and the reset names the density track's values.
    source = _circular_canvas_of_50_and_5_kb(tmp_path, window=300, step=30, feature_density_tracks=_gene_density())
    region = derive_region_request(source, [RegionSelection("pair:2", 1, 5_000)])
    assert region.mode == "circular"
    assert {item.setting: (item.source_auto, item.target_auto) for item in region.adaptation.reset} == {
        "window": (1_000, 500), "step": (100, 50),
    }
    # A track with its own window and step reads neither: the hand-set values stay.
    from gbdraw.analysis.feature_density import FeatureDensityTrack
    from gbdraw.features.query import FeatureQuery

    own = (FeatureDensityTrack("d", FeatureQuery("gene"), window=2_000, step=200),)
    source = _circular_canvas_of_50_and_5_kb(tmp_path, window=300, step=30, feature_density_tracks=own)
    assert derive_region_request(source, [RegionSelection("pair:2", 1, 5_000)]).adaptation.reset == ()


@pytest.mark.parametrize(("enabled", "reset"), [(True, {"window": (1_000, 500), "step": (100, 50)}), (False, {})])
def test_a_circular_region_compares_the_density_window_only_of_drawn_density_tracks(
    tmp_path: Path, enabled: bool, reset: dict[str, tuple[int, int]]
) -> None:
    # A stated stack draws only the feature density tracks its enabled slots name: with the slot disabled, no track
    # reads the canvas's shared density window, so the hand-set window and step stay (the GC content window of the
    # 5 kb record is 500 bp in both drawings).
    from gbdraw.api import CircularRequestTrackOptions, CircularTrackSlot

    slots = (
        CircularTrackSlot("features", "features"),
        CircularTrackSlot("density_d", "feature_density", enabled=enabled, params={"track": "d"}),
        CircularTrackSlot("gc_content", "gc_content"),
    )
    source = _circular_canvas_of_50_and_5_kb(
        tmp_path, window=300, step=30, feature_density_tracks=_gene_density(),
        tracks=CircularRequestTrackOptions(circular_track_slots=slots),
    )
    region = derive_region_request(source, [RegionSelection("pair:2", 1, 5_000)])
    assert {item.setting: (item.source_auto, item.target_auto) for item in region.adaptation.reset} == reset


def test_a_region_of_a_circular_batch_compares_each_diagrams_own_size_class(tmp_path: Path) -> None:
    # OV-3253: each diagram of a batch takes the size class of its own record. The 5 kb diagram's Auto label font is
    # the short one, as is that of a drawing of the whole 5 kb record, so the hand-set font stays; the 300 kb
    # diagram's is the long one, so a drawing of 5 kb of it returns the font to Auto.
    from gbdraw.api.requests import CircularBatchOutputPolicy, CircularBatchRequest

    batch = CircularBatchRequest(
        records=_circular_records(tmp_path, 300_000, 5_000),
        options=CircularDiagramOptions(config_overrides={"labels.font_size.short": 11, "labels.font_size.long": 11}),
        output_policy=CircularBatchOutputPolicy(output_prefix=str(tmp_path / "o")),
    )
    region = derive_region_request(batch, [RegionSelection("pair:2", 1, 5_000)])
    assert region.adaptation.reset == ()
    overrides = region.request.options.config_overrides
    assert {key: overrides[key] for key in ("labels.font_size.short", "labels.font_size.long")} == {
        "labels.font_size.short": 11, "labels.font_size.long": 11,
    }
    region = derive_region_request(batch, [RegionSelection("pair:1", 1, 5_000)], mode="circular")
    assert [(item.setting, item.source_auto, item.target_auto) for item in region.adaptation.reset] == [
        ("labels.font_size", 8.0, 14.0)
    ]


def test_each_diagram_of_a_circular_batch_draws_its_density_track_at_its_own_window(cfg, tmp_path: Path) -> None:
    # The renderer's side of the batch comparison above: the 50 kb diagram samples per 1 kb, the 5 kb one per 500 bp
    # (the window record_auto_values gives each record of a batch); a canvas of both samples per 1 kb.
    import re

    from gbdraw.api import render_request
    from gbdraw.api.requests import CircularBatchOutputPolicy, CircularBatchRequest

    batch = CircularBatchRequest(
        records=_circular_records(tmp_path, 50_000, 5_000),
        options=CircularDiagramOptions(feature_density_tracks=_gene_density()),
        output_policy=CircularBatchOutputPolicy(output_prefix=str(tmp_path / "batch")),
    )
    result = render_request(batch)
    windows = [set(re.findall(r"gene per ([0-9.]+ k?bp)", Path(path).read_text())) for path in result.output_paths]
    extent = DrawnExtent("circular", (50_000, 5_000), density_window_settings=("window",), separate_records=True)
    assert [values["feature_density.window"] for values in record_auto_values(extent, cfg)] == [1_000, 500]
    assert windows == [{"1 kbp"}, {"500 bp"}]


def test_a_region_of_a_circular_batch_compares_each_diagrams_own_density_window(tmp_path: Path) -> None:
    # In separate diagrams each record draws its feature density track at its own length (500 bp for 5 kb): the
    # hand-set window and step stay.
    from gbdraw.api.requests import CircularBatchOutputPolicy, CircularBatchRequest

    canvas = _circular_canvas_of_50_and_5_kb(tmp_path, window=300, step=30, feature_density_tracks=_gene_density())
    batch = CircularBatchRequest(
        records=canvas.records, options=canvas.options, output_policy=CircularBatchOutputPolicy(output_prefix=str(tmp_path / "o"))
    )
    region = derive_region_request(batch, [RegionSelection("pair:2", 1, 5_000)])
    assert region.adaptation.reset == ()
    assert (region.request.options.window, region.request.options.step) == (300, 30)


def test_a_linear_region_of_circular_records_compares_each_records_own_auto_window(tmp_path: Path) -> None:
    # The records of a Linear drawing share the longest record's Auto window (1,000 bp at 20 kb); the 5 kb record
    # drew its GC content with 500 bp in the Circular source, so a hand-set window returns to Auto.
    source = _circular_canvas_of_50_and_5_kb(tmp_path, window=300, step=30)
    region = derive_region_request(source, [RegionSelection("pair:1", 1, 20_000), RegionSelection("pair:2", 1, 5_000)])
    assert region.mode == "linear"
    reset = {item.setting: (item.source_auto, item.target_auto) for item in region.adaptation.reset}
    # The two source records differ: no one source value to name.
    assert reset == {"window": (None, 1_000), "step": (None, 100)}
    # From the 50 kb record alone the Auto window is the same.
    region = derive_region_request(source, [RegionSelection("pair:1", 1, 20_000)])
    assert region.adaptation.reset == ()


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


def _mito_depth_source(tmp_path: Path, *depth_slots: Any, **depth_track: Any) -> CircularDiagramRequest:
    """The human mtDNA drawn in Circular with one Depth series and a Percent GC axis whose tick font is set by hand,
    as the Web and the CLI state the Percent mode: a config override."""
    from gbdraw.api import CircularRequestTrackOptions, CircularTrackSlot

    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"NC_012920.1\t{position}\t{position % 9}\n" for position in range(1, 16570)))
    slots = (
        CircularTrackSlot("features", "features"),
        *depth_slots,
        CircularTrackSlot("gc_content", "gc_content"),
    )
    return _mito_circular(
        depth_tracks=(DepthTrackInput(source=str(depth), **depth_track),),
        config_overrides={"objects.gc_content.mode": "percent", "objects.gc_content.tick_font_size": 7},
        tracks=CircularRequestTrackOptions(circular_track_slots=slots) if depth_slots else None,
    )


def _gc_font_adaptation(region: Any) -> tuple[Any, ...]:
    kept = region.adaptation.kept.get("objects.gc_content.tick_font_size")
    reset = [
        (item.explicit, item.source_auto, item.target_auto)
        for item in region.adaptation.reset
        if item.setting == "objects.gc_content.tick_font_size"
    ]
    return ("kept", kept) if kept is not None else ("reset", *reset)


def test_a_region_compares_the_gc_axis_font_the_source_draws_beside_a_narrow_depth_row(tmp_path: Path) -> None:
    # OV-2877: a Percent GC axis without its own font is drawn in the font of the first Depth row
    # (CircularScalarAxes.gc_tick_font_size); a 20 px row draws 5 px, the Linear GC axis 8 px, so the hand-set 7 px
    # returns to Auto. The Percent mode and the row's width come from the source's own settings.
    from gbdraw.api import CircularTrackSlot, ScalarSpec

    source = _mito_depth_source(tmp_path, CircularTrackSlot("depth", "depth", width=ScalarSpec(20.0, "px")))
    region = derive_region_request(source, [RegionSelection("mito", 12337, 14148)])
    assert _gc_font_adaptation(region) == ("reset", (7, 5.0, 8.0))


def test_a_region_compares_the_gc_axis_font_of_the_depth_series_own_font(tmp_path: Path) -> None:
    # A Depth series' own tick font is the font the Percent GC axis takes, at any row width.
    source = _mito_depth_source(tmp_path, tick_font_size=12)
    linear = derive_region_request(source, [RegionSelection("mito", 12337, 14148)])
    assert _gc_font_adaptation(linear) == ("reset", (7, 12.0, 8.0))
    # A whole-record Circular drawing draws the same Depth series beside its GC axis: the size is kept.
    whole = derive_region_request(source, [RegionSelection("mito", 1, 16569)])
    assert whole.mode == "circular"
    assert _gc_font_adaptation(whole) == ("kept", 7)


def test_a_region_compares_the_gc_axis_font_of_the_first_depth_row_of_the_stack(tmp_path: Path) -> None:
    # The first Depth row of the stack draws series 1 at a share of the radius (0.05 of 390 px: 19.5 px draws 5 px);
    # the Depth row of series 0 after it (60 px) would draw 8 px, the Linear GC axis's size.
    from gbdraw.api import CircularTrackSlot, ScalarSpec

    source = _mito_depth_source(
        tmp_path,
        CircularTrackSlot("depth_2", "depth", width=ScalarSpec(0.05, "factor"), params={"track_index": 1}),
        CircularTrackSlot("depth", "depth", width=ScalarSpec(60.0, "px")),
    )
    [series] = source.options.depth_tracks
    source = replace(source, options=replace(source.options, depth_tracks=(series, series)))
    region = derive_region_request(source, [RegionSelection("mito", 12337, 14148)])
    assert _gc_font_adaptation(region) == ("reset", (7, 5.0, 8.0))


def _adaptation_of(region: Any) -> tuple[Any, ...]:
    reset = [(item.setting, item.explicit, item.source_auto, item.target_auto) for item in region.adaptation.reset]
    return (dict(region.adaptation.kept), reset)


def test_circular_track_slots_given_as_strings_adapt_like_typed_slots(tmp_path: Path) -> None:
    # OV-2935: the CLI spelling of a slot ("id:renderer@width=20px") and the typed slot are one drawing; the Depth row's
    # hand-set width is reset in a Linear drawing, and the widths are kept in a whole-record Circular one, in both.
    from gbdraw.api import CircularRequestTrackOptions, CircularTrackSlot, ScalarSpec

    typed = (
        CircularTrackSlot("features", "features", width=ScalarSpec(30.0, "px")),
        CircularTrackSlot("depth", "depth", width=ScalarSpec(20.0, "px")),
        CircularTrackSlot("gc_content", "gc_content", width=ScalarSpec(60.0, "px")),
    )
    strings = ("features:features@width=30px", "depth:depth@width=20px", "gc_content:gc_content@width=60px")
    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"NC_012920.1\t{position}\t{position % 9}\n" for position in range(1, 16570)))

    def derive(slots: Any, selection: RegionSelection, mode: str) -> Any:
        source = _mito_circular(
            depth_tracks=(DepthTrackInput(source=str(depth)),),
            tracks=CircularRequestTrackOptions(circular_track_slots=slots),
        )
        return derive_region_request(source, [selection], mode=mode)

    for selection, mode in ((RegionSelection("mito", 12337, 14148), "linear"), (RegionSelection("mito", 1, 16000), "circular")):
        from_typed = derive(typed, selection, mode)
        from_strings = derive(strings, selection, mode)
        assert _adaptation_of(from_strings) == _adaptation_of(from_typed)
        assert from_typed.adaptation.kept or from_typed.adaptation.reset
    linear = derive(strings, *(RegionSelection("mito", 12337, 14148), "linear"))
    assert {item.setting for item in linear.adaptation.reset} >= {"circularSlot.depth.width", "circularSlot.features.width"}
    whole = derive(strings, RegionSelection("mito", 1, 16000), "circular")
    assert whole.adaptation.kept["circularSlot.features.width"] == 30.0
    slots = whole.request.options.tracks.circular_track_slots
    assert [(slot.id, slot.width) for slot in slots] == [(slot.id, slot.width) for slot in typed]


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


@pytest.mark.parametrize(("start", "end"), [(16000, 16570), (0, 200)])
def test_a_region_outside_the_record_names_the_record_bounds(materialize, start: int, end: int) -> None:
    from gbdraw.web_support.error_adapter import serialize_web_error

    with pytest.raises(ValidationError, match=r"is outside 1\.\.16569") as refusal:
        derive_region_drawing(materialize(_mito_circular()), [RegionSelection("mito", start, end)])
    assert refusal.value.diagnostic == {"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"}
    shown = serialize_web_error(refusal.value, operation="deriveRegionDrawing", stage="helper")
    assert (shown["code"], shown["context"]) == ("REGION_INVALID", {"field": "region", "reason": "RECORD_BOUNDS"})


def test_a_region_across_the_origin_of_a_circular_record_is_one_window(materialize) -> None:
    # D-34: a region whose start is after its end runs across the origin of a circular record; it is drawn as one window.
    from gbdraw.io.regions import RecordWindow

    [record] = derive_region_request(_mito_circular(), [RegionSelection("mito", 16000, 200)]).request.records
    assert record.windows == (RecordWindow(16000, 200),)
    assert record.display.is_circular is True
    derived = derive_region_drawing(materialize(_mito_circular()), [RegionSelection("mito", 16000, 200)])
    assert derived.drawing is not None


@pytest.mark.parametrize(
    ("regions", "margin", "windows"),
    [
        ([(16000, 200)], 100, [(15900, 300)]),
        ([(16000, 200), (150, 400)], 0, [(16000, 400)]),
        ([(16000, 200), (16400, 16569), (1, 100)], 0, [(16000, 200)]),
        ([(16000, 200), (5000, 6000)], 0, [(16000, 200), (5000, 6000)]),
        ([(16000, 200)], 8500, []),
        ([(1, 100), (16400, 16569)], 0, [(1, 100), (16400, 16569)]),
    ],
    ids=["margin", "joined-after", "inside", "two-windows", "whole-turn", "no-crossing-stays-apart"],
)
def test_regions_across_the_origin_join_round_the_circle(regions, margin: int, windows) -> None:
    from gbdraw.io.regions import RecordWindow

    [record] = derive_region_request(
        _mito_circular(), [RegionSelection("mito", start, end) for start, end in regions], margin=margin
    ).request.records
    assert [(window.start, window.end) for window in record.windows] == windows
    assert all(isinstance(window, RecordWindow) for window in record.windows)


ORIGIN_CROSSING = json.loads(
    (REPO / "tests" / "fixtures" / "region_origin_crossing_vectors.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize("case", ORIGIN_CROSSING, ids=lambda case: case["name"])
def test_python_draws_across_the_origin_exactly_the_spans_the_web_says_cross_it(case: dict[str, Any]) -> None:
    """The origin-crossing twin (tools/twins.json): tests/web/region-axis.test.mjs reads the same vectors."""
    from gbdraw.api import RecordDisplayOptions
    from gbdraw.api.request_render import plan_request

    axis = case["axisSource"]
    request = _mito_circular(record={
        "display": RecordDisplayOptions(start_coordinate=axis["anchorBase"]),
        "presentation": RecordPresentation(reverse_complement=axis["step"] == -1),
    })
    # The vector's axis is the one Python writes for this record and display.
    [drawn] = plan_request(request).build()._gbdraw_track_slot_geometry["records"]
    assert drawn["axis"]["source"] == axis
    span = case["sourceSpan"]
    selection = [RegionSelection("mito", span["start"], span["end"])]
    [record] = derive_region_request(request, selection).request.records
    if case["pythonOutcome"] == "acrossOrigin":
        # D-34: Python draws the span the Web says crosses the origin as one window across it.
        assert (record.windows[0].start, record.windows[0].end) == (span["start"], span["end"])
        assert record.windows[0].crosses_origin
        return
    assert case["pythonOutcome"] == "withinRecord"
    whole = (span["start"], span["end"]) == (1, axis["sourceLength"])
    assert (record.windows == ()) if whole else ((record.windows[0].start, record.windows[0].end) == (span["start"], span["end"]))


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
    assert request["records"][0]["windows"][0]["start"] == 100
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


def test_the_web_helper_names_only_the_resources_the_derived_request_reads(tmp_path: Path) -> None:
    examples = REPO / "examples"
    source = LinearDiagramRequest(
        records=[
            RecordInput(GenBankInputSource(examples / "MjeNMV.gb"), record_key="a"),
            RecordInput(GenBankInputSource(examples / "MelaMJNV.gb"), record_key="b"),
        ],
        options=LinearDiagramOptions(blast_files=[str(examples / "MjeNMV.MelaMJNV.tblastx.out")]),
    )
    payload, paths, reply = _web_region(
        source, [{"recordKey": "a", "start": 1, "end": 30000}, {"recordKey": "b", "start": 1, "end": 30000}], tmp_path
    )
    assert reply["status"] == "derived"
    # The comparison table is not carried, so the one-drawing document holds no copy of it.
    [table_id] = [resource_id for resource_id, path in paths.items() if path.endswith(".tblastx.out")]
    named = {record["source"]["resourceId"] for record in reply["renderRequest"]["records"]}
    assert table_id not in reply["resourceSources"].values()
    assert set(reply["resourceSources"]) == named


def test_the_web_helper_takes_an_origin_crossing_region(tmp_path: Path) -> None:
    reply = _web_region(_mito_circular(), [{"recordKey": "mito", "start": 16000, "end": 200}], tmp_path)[2]
    assert reply["status"] == "derived"
    [record] = reply["renderRequest"]["records"]
    assert record["windows"] == [{"start": 16000, "end": 200, "reverseComplement": False, "column": None}]


def test_the_web_helper_refuses_a_start_after_the_end_of_a_linear_record(tmp_path: Path) -> None:
    from gbdraw.api import RecordDisplayOptions

    linear = _mito_circular(record={"display": RecordDisplayOptions(is_circular=False)})
    reply = _web_region(linear, [{"recordKey": "mito", "start": 16000, "end": 200}], tmp_path)[2]
    assert reply == {"status": "refused", "reason": "START_AFTER_END", "recordKey": "mito", "start": 16000, "end": 200}


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


def test_a_feature_across_the_origin_is_one_window_across_it(tmp_path: Path) -> None:
    from gbdraw.api.region_drawing import FeatureRegionSelection
    from gbdraw.io.regions import RecordWindow

    source = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_origin_record(tmp_path)), record_key="o")])
    [record] = derive_region_request(
        source, [FeatureRegionSelection("o", _feature_id(source, "across"))]
    ).request.records
    assert record.windows == (RecordWindow(901, 50),)


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
    assert reply["renderRequest"]["records"][0]["windows"][0]["start"] == 101
    across_reply = _web_region(source, [{"recordKey": "o", "featureId": across}], tmp_path)[2]
    assert across_reply["status"] == "derived"
    assert [(w["start"], w["end"]) for w in across_reply["renderRequest"]["records"][0]["windows"]] == [(901, 50)]
