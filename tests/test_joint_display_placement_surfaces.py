"""Public C/D adapters and their shared persisted intent."""

import importlib
import copy
import gzip
import json
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

import gbdraw
import gbdraw.api as api
from gbdraw.api.request_render import build_request_diagram, plan_request
from gbdraw.api.requests import InMemoryRecordSource, RecordInput
from gbdraw.features.placement import FeaturePlacementOverride, FeaturePlacementTarget
from gbdraw.session_request_codec import CanonicalRequestDecodingError, encode_canonical_request, decode_canonical_request


def source_record():
    record = SeqRecord(Seq("ACGT" * 90), id="shared", annotations={
        "molecule_type": "DNA", "topology": "circular",
    })
    record.features = [SeqFeature(SimpleLocation(20, 105, strand=1), type="CDS",
                                 qualifiers={"locus_tag": ["chosen"]})]
    return record


def test_public_placement_types_are_the_existing_typed_classes():
    for namespace in (gbdraw, api):
        for cls in (FeaturePlacementTarget, FeaturePlacementOverride):
            assert getattr(namespace, cls.__name__) is cls
            assert cls.__name__ in namespace.__all__


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("surface", ["dataframe", "path", "exact"])
def test_joint_root_typed_cli_intent_and_svg(mode, surface, tmp_path):
    source_path = tmp_path / "source.gbk"
    SeqIO.write(source_record(), source_path, "genbank")
    source = SeqIO.read(source_path, "genbank")
    table = DataFrame([{"record": "#1", "feature_selector": "locus_tag=chosen",
                        "placement": "main", "level": ""}])
    table_path = tmp_path / "placement.tsv"
    table.to_csv(table_path, sep="\t", index=False)
    cls = getattr(api, f"{mode.title()}DiagramRequest")
    options_cls = getattr(api, f"{mode.title()}DiagramOptions")
    display = api.RecordDisplayOptions(start_coordinate=71)
    request = cls(records=(RecordInput(InMemoryRecordSource(source), display=display),),
                  options=options_cls(feature_placement_table=table,
                      config_overrides={"canvas.feature_overlap_tolerance_bp": 1}))
    plan = plan_request(request)
    exact = plan.request.options.feature_placements
    assert len(exact) == 1
    placements = {"dataframe": table, "path": table_path, "exact": exact}[surface]
    public_options = getattr(gbdraw, f"{mode.title()}Options")(
        features=gbdraw.FeatureOptions(placements=placements),
        config_overrides={"canvas.feature_overlap_tolerance_bp": 1})
    public = getattr(gbdraw, f"draw_{mode}")(
        source, options=public_options, record_displays=[display])
    # Use the facade's existing mode defaults for a whole-SVG comparison.
    interface = importlib.import_module("gbdraw.interface")
    typed_options = getattr(interface, f"_{mode}_options")(public_options, record_count=1)
    typed = build_request_diagram(replace(request, options=typed_options))
    assert public._drawing.tostring() == typed.drawing.tostring()
    cli = importlib.import_module(f"gbdraw.{mode}")
    result = getattr(cli, f"run_{mode}_from_namespace")(cli._get_args([
        "--gbk", str(source_path), "--display_start_coordinate", "71",
        "--feature_placement_table", str(table_path),
        "--feature_overlap_tolerance_bp", "1", "-o", str(tmp_path / mode),
    ]))
    cli_plan = plan_request(result.canonical_request)
    assert cli_plan.request.options.feature_placements == exact
    assert cli_plan.transforms[0].source_base_to_display_index(71) == 0
    assert result.outputs[0].svg_path.stat().st_size > 1000
    # The CLI dispatch and typed renderer consume the identical materialized request.
    from tests.utils.svg_compare import compare_svgs
    expected_svg = tmp_path / "typed-cli.svg"
    built_cli = build_request_diagram(result.canonical_request)
    drawing = built_cli.items[0].drawing if mode == "circular" else built_cli.drawing
    expected_svg.write_text(drawing.tostring())
    assert compare_svgs(result.outputs[0].svg_path, expected_svg).equal


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_joint_canonical_round_trip(mode, tmp_path):
    request = getattr(api, f"{mode.title()}DiagramRequest")(
        records=(RecordInput(InMemoryRecordSource(source_record()),
                             display=api.RecordDisplayOptions(None, 1)),),
        options=getattr(api, f"{mode.title()}DiagramOptions")(
            feature_placement_table=DataFrame([{
                "feature_selector": "locus_tag=chosen", "placement": "main",
            }]), config_overrides={"canvas.feature_overlap_tolerance_bp": 1}),
    )
    materialized = plan_request(request).request
    encoded = encode_canonical_request(materialized)
    assert encoded.payload["schema"] == 7
    assert encoded.payload["records"][0]["display"] == {
        "isCircular": None, "startCoordinate": 1,
    }
    rows = encoded.payload["diagramOptions"]["featurePlacements"]
    assert rows[0]["placement"] == {"kind": "main"}
    paths = {}
    for resource in encoded.resources:
        path = tmp_path / resource.name
        path.write_bytes(resource.content if resource.content is not None
                         else resource.source_path.read_bytes())
        paths[resource.resource_id] = path
    decoded = decode_canonical_request(encoded.payload, resource_paths=paths, output_directory=tmp_path)
    assert decoded.records[0].display == request.records[0].display
    assert decoded.options.feature_placements == materialized.options.feature_placements
    assert plan_request(decoded).transforms[0].source_base_to_display_index(1) == 0


def test_historical_v40_schema6_typed_promotion_without_render(tmp_path):
    fixture = Path(__file__).parent / "fixtures/sessions/test_linear_cli_sidecar_reuses0.v40-schema6.json.gz"
    historical = json.loads(gzip.decompress(fixture.read_bytes()))
    assert historical["version"] == 40
    assert historical["renderRequest"]["schema"] == 6
    document = api.load_session_document(fixture)
    with api.materialize_session(document, output_directory=tmp_path) as materialized:
        request = api.session_to_request(materialized)
        assert all(record.display == api.RecordDisplayOptions() for record in request.records)
        assert not request.options.feature_placements
        promoted = api.build_session_document(request).to_dict()
    assert promoted["version"] == 41
    assert promoted["renderRequest"]["schema"] == 7
    assert [r["cardinality"] for r in promoted["renderRequest"]["records"]] == [
        r["cardinality"] for r in historical["renderRequest"]["records"]
    ]


@pytest.mark.parametrize("mutation", ["display-field", "linear-start", "placement-field", "auto", "level", "side"])
def test_current_codec_rejects_non_requested_or_invalid_intent(mutation, tmp_path):
    request = api.LinearDiagramRequest(records=(RecordInput(InMemoryRecordSource(source_record())),))
    encoded = encode_canonical_request(request)
    payload = copy.deepcopy(encoded.payload)
    row = {"recordKey": payload["records"][0]["recordKey"], "biologicalFeatureId": "known",
           "placement": {"kind": "main"}}
    if mutation == "display-field":
        payload["records"][0]["display"]["pixelOffset"] = 1
    elif mutation == "linear-start":
        payload["records"][0]["display"] = {"isCircular": False, "startCoordinate": 1}
    else:
        if mutation == "placement-field":
            row["feature_track_id"] = 1
        elif mutation == "auto":
            row["placement"] = {"kind": "auto"}
        else:
            row["placement"] = {"kind": "lane", "side": "outward" if mutation == "side" else "above",
                                "level": 2 if mutation == "level" else 1}
        payload["diagramOptions"]["featurePlacements"] = [row]
    paths = {}
    for resource in encoded.resources:
        path = tmp_path / resource.name
        path.write_bytes(resource.content)
        paths[resource.resource_id] = path
    with pytest.raises(CanonicalRequestDecodingError):
        decode_canonical_request(payload, resource_paths=paths, output_directory=tmp_path)
