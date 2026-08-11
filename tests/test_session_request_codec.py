from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularTrackOptions,
    ColorOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearOutputOptions,
    LinearTrackOptions,
)
from gbdraw.api.config import load_default_config
from gbdraw.api.request_render import render_request
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import parse_region_spec
from gbdraw.config.models import GbdrawConfig
from gbdraw.exceptions import ValidationError
from gbdraw.features.shapes import resolve_feature_rendering
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
    UNKNOWN_FIELD_POLICY,
    CanonicalRequestDecodingError,
    CanonicalRequestEncodingError,
    EncodedCanonicalRequest,
    _read_typed_json_resource,
    decode_canonical_request,
    encode_canonical_request,
    encode_canonical_typed_resource,
)
from gbdraw.tracks import (
    CircularTrackSlot,
    ScalarSpec,
    normalize_circular_track_slots,
)


def _materialize_resources(
    encoded: EncodedCanonicalRequest,
    directory: Path,
) -> dict[str, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}
    for resource in encoded.resources:
        if resource.source_path is not None:
            result[resource.resource_id] = resource.source_path
            continue
        target = directory / resource.name
        target.write_bytes(resource.content or b"")
        result[resource.resource_id] = target
    return result


def _source_file(path: Path, content: str = "source\n") -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def test_typed_collinearity_anchor_schema2_and_schema1_compatibility(
    tmp_path: Path,
) -> None:
    anchor = CollinearityAnchor(
        query_protein_id="query-protein",
        subject_protein_id="subject-protein",
        query_record_index=0,
        subject_record_index=1,
        query_order=0,
        subject_order=0,
        query_start=10,
        query_end=20,
        subject_start=30,
        subject_end=40,
        identity=95.0,
        evalue=1e-10,
        bitscore=100.0,
        alignment_length=10,
        query_feature_svg_id="query-stable",
        subject_feature_svg_id="subject-stable",
        source="test",
        query_unit_id="query-unit",
        subject_unit_id="subject-unit",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name="Query",
        subject_display_name="Subject",
        query_view_feature_svg_id="query-view",
        subject_view_feature_svg_id="subject-view",
        query_feature_index=3,
        subject_feature_index=7,
    )
    current = json.loads(
        encode_canonical_typed_resource("anchor", anchor).decode("utf-8")
    )
    assert current["schema"] == 2

    def decode(payload: dict[str, Any], name: str) -> CollinearityAnchor:
        path = tmp_path / name
        path.write_text(json.dumps(payload), encoding="utf-8")
        return _read_typed_json_resource(
            "typed-anchor",
            value_kind="anchor",
            expected=CollinearityAnchor,
            path="test.anchor",
            resource_paths={"typed-anchor": path},
        )

    assert decode(current, "current.json") == anchor

    legacy = copy.deepcopy(current)
    legacy["schema"] = 1
    for field_name in (
        "queryViewFeatureSvgId",
        "subjectViewFeatureSvgId",
        "queryFeatureIndex",
        "subjectFeatureIndex",
    ):
        legacy["value"]["fields"].pop(field_name)
    decoded_legacy = decode(legacy, "legacy.json")
    assert decoded_legacy.query_view_feature_svg_id == ""
    assert decoded_legacy.subject_view_feature_svg_id == ""
    assert decoded_legacy.query_feature_index is None
    assert decoded_legacy.subject_feature_index is None

    missing_current = copy.deepcopy(current)
    missing_current["value"]["fields"].pop("queryFeatureIndex")
    with pytest.raises(CanonicalRequestDecodingError, match="Missing required"):
        decode(missing_current, "missing-current.json")

    missing_legacy_required = copy.deepcopy(legacy)
    missing_legacy_required["value"]["fields"].pop("queryProteinId")
    with pytest.raises(CanonicalRequestDecodingError, match="Missing required"):
        decode(missing_legacy_required, "missing-legacy.json")

    unknown_legacy = copy.deepcopy(legacy)
    unknown_legacy["value"]["fields"]["unknownIdentity"] = "bad"
    with pytest.raises(CanonicalRequestDecodingError, match="Unknown field"):
        decode(unknown_legacy, "unknown-legacy.json")


def _payload_for_schema(
    encoded: EncodedCanonicalRequest,
    schema: int,
) -> dict[str, Any]:
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    if schema < 5:
        payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey", None)
    nested_output = payload["diagramOptions"].get("output")
    if schema in {1, 2} and nested_output is not None:
        nested_output["outputPrefix"] = "ignored-legacy-prefix"
    return payload


def _standard_collinearity_payload(
    *,
    min_anchors: int = 1,
    max_gene_gap: int = 0,
    max_diagonal_drift: int = 0,
    max_conflicts: int = 0,
    max_paralog_links: int = 9,
) -> dict[str, Any]:
    return {
        "kind": "standard",
        "parameters": {
            "minAnchors": min_anchors,
            "maxGeneGap": max_gene_gap,
            "blockMergeGap": 50,
            "singletonMergeGap": 25,
            "maxDiagonalDrift": max_diagonal_drift,
            "maxConflictsInMergeGap": max_conflicts,
            "maxParalogLinksPerOrthogroup": max_paralog_links,
            "gapPenalty": 1.0,
            "nearbyDuplicateWindow": 0,
            "scoreMode": "constant",
            "constantAnchorScore": 50.0,
            "minBlockScore": None,
            "blockEvalue": None,
        },
    }


def _table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_type": ["CDS", "tRNA"],
            "caption": ["coding sequence", "transfer RNA"],
        }
    )


def test_schema5_round_trips_circular_batch_grouping_and_outputs(
    tmp_path: Path,
) -> None:
    sources = (
        _source_file(tmp_path / "a.gbk"),
        _source_file(tmp_path / "b.gbk"),
    )
    encoded = encode_canonical_request(
        CircularBatchRequest(
            records=tuple(
                RecordInput(source=GenBankInputSource(source))
                for source in sources
            ),
            outputs=(
                RenderOutputRequest(output_prefix="same-id"),
                RenderOutputRequest(output_prefix="same-id_2"),
            ),
        )
    )

    assert encoded.payload["schema"] == 5
    assert encoded.payload["grouping"] == "batch"
    assert [item["prefix"] for item in encoded.payload["output"]] == [
        "same-id",
        "same-id_2",
    ]

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "resources"),
        output_directory=tmp_path / "output",
    )

    assert isinstance(decoded, CircularBatchRequest)
    assert decoded.grouping == "batch"
    assert [item.output_prefix for item in decoded.outputs] == [
        "same-id",
        "same-id_2",
    ]


@pytest.mark.parametrize("schema", (3, 4))
def test_branch_only_request_schemas_are_rejected(
    tmp_path: Path,
    schema: int,
) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema

    with pytest.raises(
        CanonicalRequestDecodingError,
        match=f"Unsupported canonical request schema: {schema}",
    ):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / "resources",
            ),
            output_directory=tmp_path / "output",
        )


def test_codec_module_does_not_import_cli_or_legacy_session_owners() -> None:
    source_path = Path(__file__).parents[1] / "gbdraw" / "session_request_codec.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imported_modules = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    imported_modules.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )

    forbidden_prefixes = (
        "gbdraw.circular",
        "gbdraw.linear",
        "gbdraw.session_io",
        "gbdraw.cli",
    )
    assert not any(
        module == prefix or module.startswith(f"{prefix}.")
        for module in imported_modules
        for prefix in forbidden_prefixes
    )


def test_circular_request_payload_round_trip_uses_stable_resources(tmp_path: Path) -> None:
    gbk = _source_file(tmp_path / "source-a.gbk")
    gff = _source_file(tmp_path / "source-b.gff3")
    fasta = _source_file(tmp_path / "source-b.fna")
    table = _table()
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(gbk),
                selector=parse_record_selector("record-a"),
                presentation=RecordPresentation(
                    label="Record A", grid_row=1, grid_column=1
                ),
            ),
            RecordInput(
                source=GffFastaInputSource(gff, fasta),
                region=parse_region_spec("1-100:rc"),
                presentation=RecordPresentation(
                    subtitle="annotation", grid_row=1, grid_column=2
                ),
            ),
        ),
        options=CircularDiagramOptions(
            colors=ColorOptions(
                color_table=table,
                default_colors_palette="soft",
            ),
            tracks=CircularTrackOptions(
                circular_track_slots=(
                    CircularTrackSlot(
                        id="features",
                        renderer="features",
                        radius=ScalarSpec(0.8),
                    ),
                ),
                circular_track_axis_index=0,
            ),
            selected_features_set=("CDS", "tRNA"),
            feature_visibility_table=table,
            depth_track_tables=((table, None), (None, table)),
            conservation_dataframes=(table,),
            conservation_reference="query",
            species="Example species",
            evalue=1e-10,
            output=CircularOutputOptions(
                legend="left",
                plot_title_position="top",
            ),
        ),
        layout=CircularMultiRecordOptions(
            multi_record_size_mode="auto",
            multi_record_min_radius_ratio=0.6,
            multi_record_column_gap_ratio=0.2,
            multi_record_row_gap_ratio=0.1,
        ),
        output=RenderOutputRequest(
            output_prefix="canonical-circular",
            output_directory=tmp_path / "untrusted-output",
            formats=("svg", "interactive_svg"),
            overwrite=True,
            interactive_metadata_policy="required",
        ),
    )

    encoded = encode_canonical_request(request)

    assert encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA
    assert encoded.payload["mode"] == "circular"
    assert encoded.payload["comparisons"] == []
    assert encoded.payload["records"][0]["source"] == {
        "kind": "genbank",
        "resourceId": "record-1-genbank",
    }
    serialized = json.dumps(encoded.payload)
    assert str(tmp_path) not in serialized
    assert "outputDirectory" not in encoded.payload["output"]
    assert encoded.payload["output"]["prefix"] == "canonical-circular"
    assert encoded.payload["output"]["overwrite"] is False
    assert encoded.payload["diagramOptions"]["output"] == {
        "legend": "left",
        "plotTitlePosition": "top",
    }

    resource_paths = _materialize_resources(encoded, tmp_path / "materialized")
    replay_directory = tmp_path / "replay-output"
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=resource_paths,
        output_directory=replay_directory,
    )

    assert isinstance(decoded, CircularDiagramRequest)
    assert decoded.output.output_directory == replay_directory
    assert isinstance(decoded.options.feature_visibility_table, pd.DataFrame)
    assert isinstance(decoded.options.tracks.circular_track_slots[0], CircularTrackSlot)
    assert decoded.options.tracks.circular_track_slots[0].radius == ScalarSpec(0.8)
    assert decoded.options.output == CircularOutputOptions(
        legend="left",
        plot_title_position="top",
    )
    assert decoded.output.output_prefix == "canonical-circular"
    assert decoded.output.overwrite is False
    assert decoded.records[1].region == parse_region_spec("1-100:rc")
    assert encode_canonical_request(decoded).payload == encoded.payload

    legacy_unsafe_payload = copy.deepcopy(encoded.payload)
    legacy_unsafe_payload["output"]["overwrite"] = True
    legacy_unsafe = decode_canonical_request(
        legacy_unsafe_payload,
        resource_paths=resource_paths,
        output_directory=replay_directory,
    )
    assert legacy_unsafe.output.overwrite is False


@pytest.mark.parametrize("schema", sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS))
@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_supported_schemas_restore_mode_specific_nested_option_types(
    tmp_path: Path,
    schema: int,
    mode: str,
) -> None:
    source = _source_file(tmp_path / f"{mode}-{schema}.gbk")
    record = RecordInput(source=GenBankInputSource(source))
    if mode == "circular":
        request = CircularDiagramRequest(
            records=(record,),
            options=CircularDiagramOptions(
                tracks=CircularTrackOptions(center_reserved_radius=12),
                output=CircularOutputOptions(
                    legend="left",
                    plot_title_position="top",
                ),
            ),
            layout=CircularMultiRecordOptions(
                multi_record_positions=("#1@1",),
            ),
        )
    else:
        request = LinearDiagramRequest(
            records=(record,),
            options=LinearDiagramOptions(
                tracks=LinearTrackOptions(),
                output=LinearOutputOptions(
                    legend="left",
                    plot_title_position="top",
                ),
            ),
        )
    encoded = encode_canonical_request(request)

    decoded = decode_canonical_request(
        _payload_for_schema(encoded, schema),
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"materialized-{mode}-{schema}",
        ),
        output_directory=tmp_path / f"output-{mode}-{schema}",
    )

    if mode == "circular":
        assert isinstance(decoded, CircularDiagramRequest)
        assert isinstance(decoded.options.tracks, CircularTrackOptions)
        assert isinstance(decoded.options.output, CircularOutputOptions)
        assert decoded.layout is not None
        assert decoded.layout.multi_record_positions == ("#1@1",)
    else:
        assert isinstance(decoded, LinearDiagramRequest)
        assert isinstance(decoded.options.tracks, LinearTrackOptions)
        assert isinstance(decoded.options.output, LinearOutputOptions)


@pytest.mark.parametrize("schema", sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS))
@pytest.mark.parametrize(
    ("mode", "field_name", "default_value"),
    [
        ("circular", "depthTrackHeights", None),
        ("linear", "conservationReference", "auto"),
    ],
)
def test_supported_schemas_ignore_wrong_mode_shared_defaults(
    tmp_path: Path,
    schema: int,
    mode: str,
    field_name: str,
    default_value: object,
) -> None:
    source = _source_file(tmp_path / f"{mode}-{schema}.gbk")
    record = RecordInput(source=GenBankInputSource(source))
    request = (
        CircularDiagramRequest(records=(record,))
        if mode == "circular"
        else LinearDiagramRequest(records=(record,))
    )
    encoded = encode_canonical_request(request)
    payload = _payload_for_schema(encoded, schema)
    payload["diagramOptions"][field_name] = default_value

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"materialized-{mode}-{schema}",
        ),
        output_directory=tmp_path / f"output-{mode}-{schema}",
    )

    assert isinstance(
        decoded,
        CircularDiagramRequest if mode == "circular" else LinearDiagramRequest,
    )


@pytest.mark.parametrize(
    ("mode", "field_name", "value"),
    [
        ("circular", "depthTrackHeights", [20]),
        ("linear", "conservationReference", "query"),
    ],
)
def test_current_schema_rejects_nondefault_wrong_mode_shared_options(
    tmp_path: Path,
    mode: str,
    field_name: str,
    value: object,
) -> None:
    source = _source_file(tmp_path / f"{mode}.gbk")
    record = RecordInput(source=GenBankInputSource(source))
    request = (
        CircularDiagramRequest(records=(record,))
        if mode == "circular"
        else LinearDiagramRequest(records=(record,))
    )
    encoded = encode_canonical_request(request)
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"][field_name] = value

    with pytest.raises(CanonicalRequestDecodingError, match=field_name):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / f"materialized-{mode}",
            ),
            output_directory=tmp_path / f"output-{mode}",
        )


@pytest.mark.parametrize(
    ("mode", "field_name", "value", "message"),
    [
        ("circular", "linearTrackAxisIndex", 0, "Linear track values"),
        ("linear", "centerReservedRadius", 10, "Circular track values"),
    ],
)
def test_current_schema_rejects_populated_wrong_mode_track_fields(
    tmp_path: Path,
    mode: str,
    field_name: str,
    value: object,
    message: str,
) -> None:
    source = _source_file(tmp_path / f"{mode}.gbk")
    record = RecordInput(source=GenBankInputSource(source))
    request = (
        CircularDiagramRequest(
            records=(record,),
            options=CircularDiagramOptions(tracks=CircularTrackOptions()),
        )
        if mode == "circular"
        else LinearDiagramRequest(
            records=(record,),
            options=LinearDiagramOptions(tracks=LinearTrackOptions()),
        )
    )
    encoded = encode_canonical_request(request)
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["tracks"][field_name] = value

    with pytest.raises(CanonicalRequestDecodingError, match=message):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / f"materialized-{mode}",
            ),
            output_directory=tmp_path / f"output-{mode}",
        )


def test_linear_comparison_kinds_and_payload_round_trip(tmp_path: Path) -> None:
    gbk_a = _source_file(tmp_path / "a.gbk")
    gbk_b = _source_file(tmp_path / "b.gbk")
    nucleotide = _source_file(tmp_path / "nucleotide.tsv", "a\tb\n")
    protein_table = _table()
    member = OrthogroupMember(
        orthogroup_id="OG1",
        protein_id="protein-a",
        record_index=0,
        feature_index=1,
        record_id="record-a",
        label="Protein A",
        start=10,
        end=40,
        strand=1,
        feature_svg_id="feature-a",
        source_protein_id="source-a",
    )
    orthogroups = OrthogroupResult(
        orthogroups={"OG1": [member]},
        member_by_protein_id={"protein-a": member},
        names_by_orthogroup_id={"OG1": "Example group"},
    )
    anchor = CollinearityAnchor(
        query_protein_id="protein-a",
        subject_protein_id="protein-b",
        query_record_index=0,
        subject_record_index=1,
        query_order=0,
        subject_order=1,
        query_start=10,
        query_end=40,
        subject_start=20,
        subject_end=50,
        identity=91.5,
        evalue=1e-20,
        bitscore=100.0,
        alignment_length=30,
        query_feature_svg_id="feature-a",
        subject_feature_svg_id="feature-b",
        source="precomputed",
        query_unit_id="unit-a",
        subject_unit_id="unit-b",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name="Protein A",
        subject_display_name="Protein B",
    )
    collinearity = CollinearityResult(
        blocks=(
            CollinearityBlock(
                block_id="block-1",
                query_record_index=0,
                subject_record_index=1,
                orientation="plus",
                score=100.0,
                anchors=(anchor,),
            ),
        ),
        orthogroups=orthogroups,
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(source=GenBankInputSource(gbk_a)),
            RecordInput(source=GenBankInputSource(gbk_b)),
        ),
        options=LinearDiagramOptions(
            blast_files=(str(nucleotide),),
            protein_comparisons=(protein_table,),
            orthogroups=orthogroups,
            protein_blastp_mode="collinear",
            pairwise_match_style="curve",
            collinearity_blocks=collinearity,
            collinearity_params=LosslessCollinearityParameters(
                min_anchors=2,
                max_unit_gap=3,
            ),
            collinearity_unit_mode="locus",
            collinearity_search_scope="all",
            losatp_threads=2,
            output=LinearOutputOptions(
                legend="bottom",
                plot_title_position="bottom",
            ),
        ),
        output=RenderOutputRequest(output_prefix="canonical-linear"),
    )

    encoded = encode_canonical_request(request)

    assert [item["kind"] for item in encoded.payload["comparisons"]] == [
        "nucleotideBlast",
        "precomputedProteinComparison",
        "orthogroupResult",
        "collinearityResult",
        "generatedProteinComparison",
    ]
    assert "blastFiles" not in encoded.payload["diagramOptions"]
    assert "proteinBlastpMode" not in encoded.payload["diagramOptions"]
    assert encoded.payload["output"]["prefix"] == "canonical-linear"
    assert encoded.payload["diagramOptions"]["output"] == {
        "legend": "bottom",
        "plotTitlePosition": "bottom",
    }
    generated_pipeline = next(
        item
        for item in encoded.payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    assert generated_pipeline["settings"]["collinearityParams"]["kind"] == "lossless"

    resource_paths = _materialize_resources(encoded, tmp_path / "materialized")
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=resource_paths,
        output_directory=tmp_path / "replay",
    )

    assert isinstance(decoded, LinearDiagramRequest)
    assert decoded.options.blast_files == (str(nucleotide),)
    assert isinstance(decoded.options.protein_comparisons[0], pd.DataFrame)
    assert decoded.options.orthogroups == orthogroups
    assert decoded.options.collinearity_blocks == collinearity
    assert decoded.options.collinearity_params == request.options.collinearity_params
    assert decoded.options.output == request.options.output
    assert decoded.output.output_prefix == "canonical-linear"
    assert encode_canonical_request(decoded).payload == encoded.payload


@pytest.mark.parametrize("schema", (1, 2))
def test_supported_schemas_privately_migrate_standard_collinearity_parameters(
    tmp_path: Path,
    schema: int,
) -> None:
    source = _source_file(tmp_path / f"standard-{schema}.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=LinearDiagramOptions(protein_blastp_mode="collinear"),
        )
    )
    payload = _payload_for_schema(encoded, schema)
    pipeline = next(
        item
        for item in payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    if schema == 1:
        pipeline.pop("pairs")
    pipeline["settings"]["collinearityParams"] = _standard_collinearity_payload(
        min_anchors=3,
        max_gene_gap=7,
        max_diagonal_drift=4,
        max_conflicts=2,
    )

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"standard-{schema}-resources",
        ),
        output_directory=tmp_path / f"standard-{schema}-output",
    )

    assert decoded.options.collinearity_params == LosslessCollinearityParameters(
        min_anchors=3,
        max_unit_gap=7,
        max_diagonal_drift=4,
        max_conflicts=2,
        merge_orientation="either",
    )
    assert decoded.options.collinear_max_paralog_links_per_orthogroup == 9
    latest_pipeline = next(
        item
        for item in encode_canonical_request(decoded).payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    assert latest_pipeline["settings"]["collinearityParams"]["kind"] == "lossless"


def test_standard_collinearity_embedded_max_paralog_does_not_override_explicit_setting(
    tmp_path: Path,
) -> None:
    source = _source_file(tmp_path / "standard-explicit.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=LinearDiagramOptions(
                protein_blastp_mode="collinear",
                collinear_max_paralog_links_per_orthogroup=5,
            ),
        )
    )
    payload = _payload_for_schema(encoded, 2)
    pipeline = next(
        item
        for item in payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    pipeline["settings"]["collinearityParams"] = _standard_collinearity_payload()

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "standard-explicit-resources",
        ),
        output_directory=tmp_path / "standard-explicit-output",
    )

    assert decoded.options.collinear_max_paralog_links_per_orthogroup == 5


def test_standard_collinearity_embedded_max_paralog_is_inert_in_orthogroup_mode(
    tmp_path: Path,
) -> None:
    source = _source_file(tmp_path / "standard-orthogroup.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=LinearDiagramOptions(protein_blastp_mode="orthogroup"),
        )
    )
    payload = _payload_for_schema(encoded, 2)
    pipeline = next(
        item
        for item in payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    pipeline["settings"]["collinearityParams"] = _standard_collinearity_payload()

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "standard-orthogroup-resources",
        ),
        output_directory=tmp_path / "standard-orthogroup-output",
    )

    assert decoded.options.collinear_max_paralog_links_per_orthogroup == 2


def test_current_schema_rejects_standard_collinearity_parameters(
    tmp_path: Path,
) -> None:
    source = _source_file(tmp_path / "current-standard.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=LinearDiagramOptions(protein_blastp_mode="collinear"),
        )
    )
    pipeline = next(
        item
        for item in encoded.payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    pipeline["settings"]["collinearityParams"] = _standard_collinearity_payload()

    with pytest.raises(
        CanonicalRequestDecodingError,
        match="Unsupported collinearity parameter kind",
    ):
        decode_canonical_request(
            encoded.payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / "current-standard-resources",
            ),
            output_directory=tmp_path / "current-standard-output",
        )


def test_schema2_nested_output_prefix_is_required_but_ignored(
    tmp_path: Path,
) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=LinearDiagramOptions(
                output=LinearOutputOptions(
                    legend="left",
                    plot_title_position="top",
                ),
            ),
            output=RenderOutputRequest(output_prefix="outer-authority"),
        )
    )
    legacy_payload = copy.deepcopy(encoded.payload)
    legacy_payload["schema"] = 2
    legacy_payload.pop("grouping", None)
    legacy_payload["diagramOptions"]["output"]["outputPrefix"] = "ignored-nested"

    decoded = decode_canonical_request(
        legacy_payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "legacy"),
        output_directory=tmp_path / "output",
    )

    assert decoded.output.output_prefix == "outer-authority"
    assert decoded.options.output == LinearOutputOptions(
        legend="left",
        plot_title_position="top",
    )
    latest = encode_canonical_request(decoded).payload
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert latest["output"]["prefix"] == "outer-authority"
    assert "outputPrefix" not in latest["diagramOptions"]["output"]

    missing_legacy_prefix = copy.deepcopy(legacy_payload)
    missing_legacy_prefix["diagramOptions"]["output"].pop("outputPrefix")
    with pytest.raises(CanonicalRequestDecodingError, match="outputPrefix"):
        decode_canonical_request(
            missing_legacy_prefix,
            resource_paths=_materialize_resources(encoded, tmp_path / "missing"),
            output_directory=tmp_path / "output-missing",
        )


def test_current_schema_rejects_legacy_nested_output_prefix(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=CircularDiagramOptions(
                output=CircularOutputOptions(legend="none")
            ),
            output=RenderOutputRequest(output_prefix="outer-authority"),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["output"]["outputPrefix"] = "stale"

    with pytest.raises(CanonicalRequestDecodingError, match="outputPrefix"):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(encoded, tmp_path / "materialized"),
            output_directory=tmp_path / "output",
        )


@pytest.mark.parametrize("schema", [1, 2])
def test_legacy_circular_schemas_migrate_removed_layout_values(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                tracks=CircularTrackOptions(
                    circular_track_slots=(
                        CircularTrackSlot(
                            id="gc_content",
                            renderer="dinucleotide_content",
                        ),
                        "gc_skew:dinucleotide_skew@nt=GC",
                    ),
                ),
            ),
            layout=CircularMultiRecordOptions(multi_record_size_mode="auto"),
        )
    )
    legacy_payload = copy.deepcopy(encoded.payload)
    legacy_payload["schema"] = schema
    legacy_payload.pop("grouping", None)
    if schema == 1:
        legacy_payload["records"][0].pop("recordKey")
    legacy_payload["layout"]["multiRecordSizeMode"] = "sqrt"
    legacy_slot = legacy_payload["diagramOptions"]["tracks"]["circularTrackSlots"][0]
    legacy_slot["spacing"] = {"value": 4.0, "unit": "px"}
    legacy_payload["diagramOptions"]["tracks"]["circularTrackSlots"][1] = (
        "gc_skew:dinucleotide_skew@nt=GC,spacing=5px,"
        "strict=true,compress=true,reserve=true"
    )

    decoded = decode_canonical_request(
        legacy_payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"legacy-circular-{schema}",
        ),
        output_directory=tmp_path / "output",
    )

    assert isinstance(decoded, CircularDiagramRequest)
    assert decoded.layout is not None
    assert decoded.layout.multi_record_size_mode == "auto"
    assert decoded.options.tracks is not None
    decoded_slot = decoded.options.tracks.circular_track_slots[0]
    assert isinstance(decoded_slot, CircularTrackSlot)
    decoded_string_slot = decoded.options.tracks.circular_track_slots[1]
    assert isinstance(decoded_string_slot, CircularTrackSlot)
    normalized_slots = normalize_circular_track_slots(
        [decoded_slot, decoded_string_slot]
    )
    assert normalized_slots[0].legacy_spacing == ScalarSpec(4.0, "px")
    assert normalized_slots[1].legacy_spacing == ScalarSpec(5.0, "px")
    assert "__gbdraw_legacy_spacing" not in decoded_slot.params
    assert "__gbdraw_legacy_spacing" not in decoded_string_slot.params

    latest = encode_canonical_request(decoded).payload
    latest_slot = latest["diagramOptions"]["tracks"]["circularTrackSlots"][0]
    latest_string_slot = latest["diagramOptions"]["tracks"]["circularTrackSlots"][1]
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert latest["layout"]["multiRecordSizeMode"] == "auto"
    assert latest_slot["innerGapPx"] == pytest.approx(4.0)
    assert latest_slot["outerGapPx"] == pytest.approx(4.0)
    assert "__gbdraw_legacy_spacing" not in latest_slot["params"]
    assert {"spacing", "strict", "compress", "reserve"}.isdisjoint(latest_slot)
    assert latest_string_slot["innerGapPx"] == pytest.approx(5.0)
    assert latest_string_slot["outerGapPx"] == pytest.approx(5.0)
    assert "__gbdraw_legacy_spacing" not in latest_string_slot["params"]
    assert {"spacing", "strict", "compress", "reserve"}.isdisjoint(
        latest_string_slot
    )


def test_current_circular_writer_uses_only_canonical_layout_fields() -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                tracks=CircularTrackOptions(
                    circular_track_slots=(
                        CircularTrackSlot(
                            id="gc_content",
                            renderer="dinucleotide_content",
                            inner_gap_px=4.0,
                            outer_gap_px=6.0,
                        ),
                    ),
                ),
            ),
            layout=CircularMultiRecordOptions(multi_record_size_mode="auto"),
        )
    )

    slot = encoded.payload["diagramOptions"]["tracks"]["circularTrackSlots"][0]
    assert encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA
    assert encoded.payload["layout"]["multiRecordSizeMode"] == "auto"
    assert slot["innerGapPx"] == pytest.approx(4.0)
    assert slot["outerGapPx"] == pytest.approx(6.0)
    assert {"spacing", "strict", "compress", "reserve"}.isdisjoint(slot)


def test_current_schema_rejects_private_circular_track_params(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                tracks=CircularTrackOptions(
                    circular_track_slots=(
                        CircularTrackSlot(
                            id="gc_content",
                            renderer="dinucleotide_content",
                        ),
                    ),
                ),
            ),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["tracks"]["circularTrackSlots"][0]["params"][
        "__gbdraw_legacy_spacing"
    ] = "4px"

    with pytest.raises(CanonicalRequestDecodingError, match="private"):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / "current-schema-private-slot",
            ),
            output_directory=tmp_path / "output",
        )


def test_legacy_factor_spacing_replays_but_cannot_leak_into_current_schema(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = 2
    payload.pop("grouping", None)
    payload["diagramOptions"]["tracks"] = {
        "circularTrackSlots": [
            "gc_content:dinucleotide_content@spacing=0.1"
        ],
        "circularTrackAxisIndex": None,
        "linearTrackSlots": None,
        "linearTrackAxisIndex": None,
        "centerReservedRadius": None,
    }

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "legacy-factor-slot",
        ),
        output_directory=tmp_path / "output",
    )
    assert decoded.options.tracks is not None
    normalized = normalize_circular_track_slots(
        decoded.options.tracks.circular_track_slots
    )
    assert normalized[0].legacy_spacing == ScalarSpec(0.1, "factor")

    with pytest.raises(
        CanonicalRequestEncodingError,
        match="factor-based Circular spacing",
    ):
        encode_canonical_request(decoded)


def test_schema2_migrates_removed_feature_table_field(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    table = pd.DataFrame(
        [{"record_id": "record", "feature_type": "CDS", "visible": True}]
    )
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(feature_visibility_table=table),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = 2
    payload.pop("grouping", None)
    diagram_options = payload["diagramOptions"]
    diagram_options["featureTable"] = diagram_options.pop(
        "featureVisibilityTable"
    )

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "legacy-feature-table",
        ),
        output_directory=tmp_path / "output",
    )

    pd.testing.assert_frame_equal(
        decoded.options.feature_visibility_table,
        table,
    )
    latest = encode_canonical_request(decoded).payload["diagramOptions"]
    assert "featureTable" not in latest
    assert "featureVisibilityTable" in latest


def test_current_schema_rejects_removed_feature_table_field(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["featureTable"] = None

    with pytest.raises(CanonicalRequestDecodingError, match="featureTable"):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / "current-feature-table",
            ),
            output_directory=tmp_path / "output",
        )


def test_collinear_pipeline_ignores_legacy_derived_comparison_pairs(
    tmp_path: Path,
) -> None:
    gbk_a = _source_file(tmp_path / "a.gbk")
    gbk_b = _source_file(tmp_path / "b.gbk")
    request = LinearDiagramRequest(
        records=(
            RecordInput(source=GenBankInputSource(gbk_a)),
            RecordInput(source=GenBankInputSource(gbk_b)),
        ),
        options=LinearDiagramOptions(
            protein_blastp_mode="collinear",
            collinearity_search_scope="all",
        ),
    )
    encoded = encode_canonical_request(request)
    pipeline = next(
        item
        for item in encoded.payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    pipeline["pairs"] = [{"queryRecordIndex": 0, "subjectRecordIndex": 1}]

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "materialized-legacy"),
        output_directory=tmp_path / "replay-legacy",
    )

    assert decoded.options.protein_blastp_mode == "collinear"
    assert decoded.options.protein_comparison_pairs is None
    reencoded_pipeline = next(
        item
        for item in encode_canonical_request(decoded).payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    )
    assert reencoded_pipeline["pairs"] == []


def test_file_backed_options_and_typed_config_round_trip(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    table_file = _source_file(tmp_path / "rules.tsv", "key\tvalue\nCDS\tkeep\n")
    depth_file = _source_file(tmp_path / "depth.tsv", "record\t1\t5\n")
    blast_file = _source_file(tmp_path / "conservation.tsv", "a\tb\n")
    fasta_file = _source_file(tmp_path / "comparison.fna", ">comparison\nACGT\n")
    config = GbdrawConfig.from_dict(load_default_config())
    request = CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=CircularDiagramOptions(
            config=config,
            colors=ColorOptions(default_colors_file=str(table_file)),
            feature_visibility_table_file=str(table_file),
            label_whitelist_file=str(table_file),
            qualifier_priority_file=str(table_file),
            label_override_file=str(table_file),
            depth_track_files=((str(depth_file),),),
            conservation_blast_files=(str(blast_file),),
            conservation_fasta_files=(None, str(fasta_file)),
        ),
    )

    encoded = encode_canonical_request(request)
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "materialized-files"),
        output_directory=tmp_path / "output",
    )

    assert isinstance(decoded.options.config, GbdrawConfig)
    assert decoded.options.feature_visibility_table_file == str(table_file)
    assert set(encoded.payload["diagramOptions"]).isdisjoint(
        {
            "depthFile",
            "depthFiles",
            "depthTrackFiles",
            "depthTrackLabels",
        }
    )
    assert decoded.options.depth_tracks is not None
    assert decoded.options.depth_tracks[0].source == str(depth_file)
    assert decoded.options.depth_track_files is None
    assert decoded.options.conservation_blast_files == (str(blast_file),)
    assert decoded.options.conservation_fasta_files == (None, str(fasta_file))
    assert encode_canonical_request(decoded).payload == encoded.payload


@pytest.mark.parametrize(
    ("request_type", "options_type"),
    [
        (CircularDiagramRequest, CircularDiagramOptions),
        (LinearDiagramRequest, LinearDiagramOptions),
    ],
)
def test_schema5_round_trips_scale_visibility_in_full_config_and_override(
    tmp_path: Path,
    request_type,
    options_type,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    config_dict = load_default_config()
    config_dict["objects"]["scale"]["show"] = False
    request = request_type(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=options_type(
            config=GbdrawConfig.from_dict(config_dict),
            config_overrides={"objects.scale.show": False},
        ),
    )

    encoded = encode_canonical_request(request)
    diagram_options = encoded.payload["diagramOptions"]
    assert diagram_options["config"]["objects"]["scale"]["show"] is False
    assert (
        diagram_options["configOverrides"]["objects.scale.show"] is False
    )
    assert "show_scale" not in diagram_options
    assert "show_scale" not in diagram_options["configOverrides"]

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / request_type.__name__,
        ),
        output_directory=tmp_path / "output",
    )
    assert decoded.options.config.objects.scale.show is False
    assert decoded.options.config_overrides["objects.scale.show"] is False
    assert encode_canonical_request(decoded).payload == encoded.payload


@pytest.mark.parametrize(
    ("request_type", "options_type"),
    [
        (CircularDiagramRequest, CircularDiagramOptions),
        (LinearDiagramRequest, LinearDiagramOptions),
    ],
)
def test_schema5_round_trips_arrow_shape_and_geometry_overrides(
    tmp_path: Path,
    request_type,
    options_type,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    request = request_type(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=options_type(
            feature_shapes={"CDS": "arrow"},
            config_overrides={
                "objects.features.arrow_geometry.head_length_ratio": 1.25,
                "objects.features.arrow_geometry.shaft_width_ratio": 0.25,
            },
        ),
    )

    encoded = encode_canonical_request(request)
    diagram_options = encoded.payload["diagramOptions"]
    assert diagram_options["featureShapes"] == {"CDS": "arrow"}
    assert (
        diagram_options["configOverrides"][
            "objects.features.arrow_geometry.head_length_ratio"
        ]
        == 1.25
    )
    assert (
        diagram_options["configOverrides"][
            "objects.features.arrow_geometry.shaft_width_ratio"
        ]
        == 0.25
    )

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / request_type.__name__,
        ),
        output_directory=tmp_path / "output",
    )
    assert decoded.options.feature_shapes == {"CDS": "arrow"}
    assert decoded.options.config_overrides == request.options.config_overrides
    assert encode_canonical_request(decoded).payload == encoded.payload


@pytest.mark.parametrize(
    ("request_type", "options_type"),
    [
        (CircularDiagramRequest, CircularDiagramOptions),
        (LinearDiagramRequest, LinearDiagramOptions),
    ],
)
def test_schema5_full_config_missing_scale_visibility_defaults_visible(
    tmp_path: Path,
    request_type,
    options_type,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    encoded = encode_canonical_request(
        request_type(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=options_type(
                config=GbdrawConfig.from_dict(load_default_config()),
            ),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    del payload["diagramOptions"]["config"]["objects"]["scale"]["show"]

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / request_type.__name__,
        ),
        output_directory=tmp_path / "output",
    )

    assert decoded.options.config.objects.scale.show is True


def test_canonical_depth_tracks_round_trip_mixed_sources(tmp_path: Path) -> None:
    gbk_a = _source_file(tmp_path / "a.gbk")
    gbk_b = _source_file(tmp_path / "b.gbk")
    depth_file = _source_file(
        tmp_path / "a.depth.tsv",
        "a\t1\t5\n",
    )
    depth_table = pd.DataFrame(
        {
            "reference_name": ["b"],
            "position": [1],
            "depth": [8],
        }
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(source=GenBankInputSource(gbk_a)),
            RecordInput(source=GenBankInputSource(gbk_b)),
        ),
        options=LinearDiagramOptions(
            depth_tracks=(
                DepthTrackInput(
                    source=(depth_file, depth_table),
                    label="Coverage",
                    color="#345678",
                    height=24,
                    large_tick_interval=10,
                    small_tick_interval=5,
                    tick_font_size=9,
                ),
            ),
        ),
    )

    encoded = encode_canonical_request(request)
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "materialized-depth-tracks",
        ),
        output_directory=tmp_path / "output",
    )

    assert "depthTracks" in encoded.payload["diagramOptions"]
    assert decoded.options.depth_tracks is not None
    decoded_track = decoded.options.depth_tracks[0]
    assert isinstance(decoded_track, DepthTrackInput)
    assert isinstance(decoded_track.source, tuple)
    assert decoded_track.source[0] == str(depth_file)
    assert isinstance(decoded_track.source[1], pd.DataFrame)
    assert decoded_track.source[1].equals(depth_table)
    assert decoded_track.label == "Coverage"
    assert decoded_track.color == "#345678"
    assert decoded_track.height == 24
    assert encode_canonical_request(decoded).payload == encoded.payload


def test_canonical_depth_track_shape_and_text_are_strict(tmp_path: Path) -> None:
    gbk = _source_file(tmp_path / "record.gbk")
    depth = _source_file(tmp_path / "record.depth.tsv", "record\t1\t5\n")
    request = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(gbk)),),
        options=LinearDiagramOptions(
            depth_tracks=(DepthTrackInput(source=depth, label="Coverage"),),
        ),
    )
    encoded = encode_canonical_request(request)
    resource_paths = _materialize_resources(
        encoded,
        tmp_path / "materialized-strict-depth",
    )

    missing = copy.deepcopy(encoded.payload)
    del missing["diagramOptions"]["depthTracks"][0]["label"]
    with pytest.raises(CanonicalRequestDecodingError, match="Missing required.*label"):
        decode_canonical_request(
            missing,
            resource_paths=resource_paths,
            output_directory=tmp_path / "missing",
        )

    blank = copy.deepcopy(encoded.payload)
    blank["diagramOptions"]["depthTracks"][0]["label"] = " "
    with pytest.raises(CanonicalRequestDecodingError, match="non-empty string"):
        decode_canonical_request(
            blank,
            resource_paths=resource_paths,
            output_directory=tmp_path / "blank",
        )

    non_string = copy.deepcopy(encoded.payload)
    non_string["diagramOptions"]["depthTracks"][0]["color"] = 123
    with pytest.raises(CanonicalRequestDecodingError, match="string or null"):
        decode_canonical_request(
            non_string,
            resource_paths=resource_paths,
            output_directory=tmp_path / "non-string",
        )

    unknown = copy.deepcopy(encoded.payload)
    unknown["diagramOptions"]["depthTracks"][0]["legacyLabel"] = "old"
    with pytest.raises(CanonicalRequestDecodingError, match="Unknown field.*legacyLabel"):
        decode_canonical_request(
            unknown,
            resource_paths=resource_paths,
            output_directory=tmp_path / "unknown",
        )


def test_legacy_depth_matrix_is_written_as_canonical_tracks(tmp_path: Path) -> None:
    gbk_a = _source_file(tmp_path / "a.gbk")
    gbk_b = _source_file(tmp_path / "b.gbk")
    depth_a = _source_file(tmp_path / "a.depth.tsv", "a\t1\t5\n")
    depth_b = _source_file(tmp_path / "b.depth.tsv", "b\t1\t8\n")
    request = LinearDiagramRequest(
        records=(
            RecordInput(source=GenBankInputSource(gbk_a)),
            RecordInput(source=GenBankInputSource(gbk_b)),
        ),
        options=LinearDiagramOptions(
            depth_track_files=((str(depth_a), None), (None, str(depth_b))),
            depth_track_labels=("First", "Second"),
            depth_track_colors=("#112233", "#445566"),
            depth_track_heights=(18, 24),
            depth_track_large_tick_intervals=(10,),
            depth_track_small_tick_intervals=(2, 4),
            depth_track_tick_font_sizes=(8, 9),
        ),
    )

    encoded = encode_canonical_request(request)
    diagram_options = encoded.payload["diagramOptions"]

    assert set(diagram_options).isdisjoint(
        {
            "depthTrackFiles",
            "depthTrackLabels",
            "depthTrackColors",
            "depthTrackHeights",
            "depthTrackLargeTickIntervals",
            "depthTrackSmallTickIntervals",
            "depthTrackTickFontSizes",
        }
    )
    assert len(diagram_options["depthTracks"]) == 2
    assert diagram_options["depthTracks"][0]["source"][1] is None
    assert diagram_options["depthTracks"][1]["source"][0] is None

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / "materialized-legacy-depth",
        ),
        output_directory=tmp_path / "output",
    )
    assert decoded.options.depth_tracks is not None
    assert [track.label for track in decoded.options.depth_tracks] == [
        "First",
        "Second",
    ]
    assert [track.height for track in decoded.options.depth_tracks] == [18, 24]
    assert decoded.options.depth_track_files is None
    assert encode_canonical_request(decoded).payload == encoded.payload


def test_in_memory_record_is_encoded_as_a_genbank_resource(tmp_path: Path) -> None:
    record = SeqRecord(
        Seq("ATGCGC"),
        id="memory-record",
        annotations={"molecule_type": "DNA"},
    )
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),)
    )

    encoded = encode_canonical_request(request)

    assert encoded.payload["records"][0]["source"]["kind"] == "genbank"
    resource = encoded.resources[0]
    assert resource.content is not None
    assert b"LOCUS" in resource.content

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path),
        output_directory=tmp_path / "output",
    )
    assert isinstance(decoded.records[0].source, GenBankInputSource)


def test_annotation_targets_and_styles_round_trip(tmp_path: Path) -> None:
    record = SeqRecord(
        Seq("ATGCGC" * 20),
        id="memory-record",
        annotations={"molecule_type": "DNA"},
    )
    annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "review",
                (
                    RegionAnnotation(
                        "coordinates",
                        CoordinateSpan(
                            parse_record_selector("#1"),
                            2,
                            20,
                            coordinate_space="local",
                        ),
                        label="Review",
                        mark="band",
                        style=RegionAnnotationStyle(
                            fill="#88aaff",
                            hatch=HatchStyle(angle=30, cross=True),
                        ),
                        legend_label="Reviewed",
                        metadata={"owner": "curation"},
                    ),
                    RegionAnnotation(
                        "features",
                        FeatureSpan(
                            parse_record_selector("memory-record"),
                            (FeatureSelector("dnaA", "gene"),),
                            envelope="segments",
                            circular_path="forward",
                        ),
                    ),
                ),
                legend_label="Review regions",
            ),
        )
    )
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(annotations=annotations),
    )

    encoded = encode_canonical_request(request)
    payload = encoded.payload["diagramOptions"]["annotations"]
    assert payload["sets"][0]["annotations"][0]["target"]["kind"] == "coordinateSpan"
    assert payload["sets"][0]["annotations"][1]["target"]["kind"] == "featureSpan"

    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "annotations"),
        output_directory=tmp_path / "output",
    )
    assert decoded.options.annotations == annotations
    assert encode_canonical_request(decoded).payload == encoded.payload


def test_decode_requires_caller_owned_output_directory(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(records=(RecordInput(source=GenBankInputSource(source)),))
    )

    with pytest.raises(CanonicalRequestDecodingError, match="supplied by the caller"):
        decode_canonical_request(
            encoded.payload,
            resource_paths=_materialize_resources(encoded, tmp_path),
            output_directory=None,
        )


@pytest.mark.parametrize(
    ("request_type", "expected_options"),
    [
        (
            CircularDiagramRequest,
            {
                "evalue": 1e-5,
                "bitscore": 50.0,
                "identity": 70.0,
                "alignmentLength": 0,
                "configOverrides": {
                    "canvas.show_gc": True,
                    "canvas.show_skew": True,
                },
            },
        ),
        (
            LinearDiagramRequest,
            {
                "evalue": 1e-2,
                "bitscore": 50.0,
                "identity": 0.0,
                "alignmentLength": 0,
                "configOverrides": {
                    "canvas.show_gc": False,
                    "canvas.show_skew": False,
                    "objects.axis.linear.stroke_color": "lightgray",
                },
            },
        ),
    ],
)
def test_current_writer_serializes_resolved_mode_defaults(
    request_type,
    expected_options: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})

    encoded = encode_canonical_request(
        request_type(records=(RecordInput(source=InMemoryRecordSource(record)),))
    )

    for name, expected in expected_options.items():
        assert encoded.payload["diagramOptions"][name] == expected


@pytest.mark.parametrize("schema", [1, 2])
@pytest.mark.parametrize(
    ("request_type", "expected_overrides"),
    [
        (
            CircularDiagramRequest,
            {
                "canvas.show_gc": True,
                "canvas.show_skew": True,
            },
        ),
        (
            LinearDiagramRequest,
            {
                "canvas.show_gc": True,
                "canvas.show_skew": True,
                "objects.axis.linear.stroke_color": "gray",
            },
        ),
    ],
)
def test_sparse_supported_schemas_retain_historical_defaults(
    tmp_path: Path,
    schema: int,
    request_type,
    expected_overrides: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        request_type(records=(RecordInput(source=InMemoryRecordSource(record)),))
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey")
    for name in (
        "configOverrides",
        "selectedFeaturesSet",
        "evalue",
        "bitscore",
        "identity",
        "alignmentLength",
    ):
        payload["diagramOptions"].pop(name)

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"{request_type.__name__}-{schema}",
        ),
        output_directory=tmp_path / "output",
    )

    assert decoded.options.evalue == 1e-5
    assert decoded.options.bitscore == 50.0
    assert decoded.options.identity == 70.0
    assert decoded.options.alignment_length == 0
    assert decoded.options.selected_features_set == (
        "CDS",
        "rRNA",
        "tRNA",
        "tmRNA",
        "ncRNA",
        "misc_RNA",
        "repeat_region",
    )
    assert decoded.options.config_overrides == expected_overrides

    payload_with_null_overrides = copy.deepcopy(payload)
    payload_with_null_overrides["diagramOptions"]["configOverrides"] = {
        "show_gc": None,
        "show_skew": None,
        "linear_axis_stroke_color": None,
    }
    decoded_null_overrides = decode_canonical_request(
        payload_with_null_overrides,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"{request_type.__name__}-{schema}-null-overrides",
        ),
        output_directory=tmp_path / "output-null-overrides",
    )
    assert decoded_null_overrides.options.config_overrides == expected_overrides


@pytest.mark.parametrize(
    "request_type",
    [CircularDiagramRequest, LinearDiagramRequest],
)
def test_current_schema_sparse_options_use_current_mode_defaults(
    tmp_path: Path,
    request_type,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    request = request_type(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
    )
    encoded = encode_canonical_request(request)
    payload = copy.deepcopy(encoded.payload)
    for name in (
        "configOverrides",
        "selectedFeaturesSet",
        "evalue",
        "bitscore",
        "identity",
        "alignmentLength",
    ):
        payload["diagramOptions"].pop(name)

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / request_type.__name__,
        ),
        output_directory=tmp_path / "output",
    )

    assert decoded.options.evalue == request.options.evalue
    assert decoded.options.bitscore == request.options.bitscore
    assert decoded.options.identity == request.options.identity
    assert decoded.options.alignment_length == request.options.alignment_length
    assert decoded.options.selected_features_set == request.options.selected_features_set
    assert decoded.options.config_overrides == request.options.config_overrides


@pytest.mark.parametrize("schema", sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS))
def test_supported_schemas_preserve_explicit_serialized_defaults(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),)
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    if schema < 5:
        payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey")
    payload["diagramOptions"].update(
        {
            "evalue": 0.5,
            "bitscore": 12.0,
            "identity": 85.0,
            "alignmentLength": 123,
            "configOverrides": {
                "show_gc": False,
                "show_skew": False,
                "linear_axis_stroke_color": "hotpink",
            },
        }
    )

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"explicit-{schema}",
        ),
        output_directory=tmp_path / "output",
    )

    assert decoded.options.evalue == 0.5
    assert decoded.options.bitscore == 12.0
    assert decoded.options.identity == 85.0
    assert decoded.options.alignment_length == 123
    assert decoded.options.config_overrides == {
        "canvas.show_gc": False,
        "canvas.show_skew": False,
        "objects.axis.linear.stroke_color": "hotpink",
    }


@pytest.mark.parametrize("schema", sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS))
@pytest.mark.parametrize(
    ("request_type", "legacy_overrides", "scope_path", "expected_scope"),
    [
        (
            CircularDiagramRequest,
            {"show_labels": True, "allow_inner_labels": True},
            "labels.circular.scope",
            "both",
        ),
        (
            LinearDiagramRequest,
            {"show_labels": "orthogroup_top"},
            "labels.linear.scope",
            "orthogroup_top",
        ),
    ],
)
def test_supported_schemas_migrate_legacy_flat_label_settings(
    tmp_path: Path,
    schema: int,
    request_type,
    legacy_overrides: dict[str, object],
    scope_path: str,
    expected_scope: str,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    encoded = encode_canonical_request(
        request_type(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    if schema < 5:
        payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey")
    payload["diagramOptions"]["configOverrides"] = legacy_overrides

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"legacy-label-overrides-{request_type.__name__}-{schema}",
        ),
        output_directory=tmp_path / "output",
    )

    overrides = decoded.options.config_overrides
    assert overrides is not None
    assert overrides[scope_path] == expected_scope
    assert not set(legacy_overrides) & set(overrides)

    latest = encode_canonical_request(decoded).payload
    latest_overrides = latest["diagramOptions"]["configOverrides"]
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert latest_overrides[scope_path] == expected_scope
    assert not set(legacy_overrides) & set(latest_overrides)


@pytest.mark.parametrize("schema", sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS))
def test_supported_schemas_migrate_legacy_full_config_label_settings(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    config = GbdrawConfig.from_dict(load_default_config())
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(config=config),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    if schema < 5:
        payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey")
    legacy_config = payload["diagramOptions"]["config"]
    legacy_config["canvas"]["show_labels"] = True
    legacy_config["canvas"]["circular"]["show_labels"] = True
    legacy_config["canvas"]["linear"]["show_labels"] = True
    legacy_config["canvas"]["circular"]["allow_inner_labels"] = True
    legacy_config["canvas"]["linear"]["track_layout"] = "spreadout"
    legacy_config["labels"]["circular"].pop("scope")
    legacy_config["labels"]["linear"].pop("scope")
    legacy_config["labels"]["linear"]["placement"] = "on_feature"

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"legacy-label-config-{schema}",
        ),
        output_directory=tmp_path / "output",
    )

    decoded_config = decoded.options.config
    assert isinstance(decoded_config, GbdrawConfig)
    assert decoded_config.labels.circular.scope == "both"
    assert decoded_config.labels.linear.scope == "all"
    assert decoded_config.canvas.linear.track_layout == "above"
    assert decoded_config.labels.linear.placement == "above_feature"

    latest = encode_canonical_request(decoded).payload
    latest_config = latest["diagramOptions"]["config"]
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert latest_config["labels"]["circular"]["scope"] == "both"
    assert latest_config["labels"]["linear"]["scope"] == "all"
    assert latest_config["canvas"]["linear"]["track_layout"] == "above"
    assert latest_config["labels"]["linear"]["placement"] == "above_feature"
    assert "show_labels" not in latest_config["canvas"]
    assert "show_labels" not in latest_config["canvas"]["circular"]
    assert "show_labels" not in latest_config["canvas"]["linear"]
    assert "allow_inner_labels" not in latest_config["canvas"]["circular"]


@pytest.mark.parametrize("schema", [1, 2])
@pytest.mark.parametrize(
    ("legacy_layout", "current_layout"),
    [
        ("spreadout", "above"),
        ("tuckin", "below"),
    ],
)
def test_legacy_schemas_migrate_retired_config_overrides(
    tmp_path: Path,
    schema: int,
    legacy_layout: str,
    current_layout: str,
) -> None:
    record = SeqRecord(
        Seq("ATGC" * 100),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            output=RenderOutputRequest(
                output_prefix=f"legacy-{schema}-{legacy_layout}",
                overwrite=True,
            ),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    payload.pop("grouping", None)
    if schema == 1:
        payload["records"][0].pop("recordKey")
    payload["diagramOptions"]["configOverrides"] = {
        "depth_tick_interval": 12,
        "label_placement": "on_feature",
        "linear_track_layout": legacy_layout,
    }

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(
            encoded,
            tmp_path / f"legacy-overrides-{schema}-{legacy_layout}",
        ),
        output_directory=tmp_path / "output",
    )

    overrides = decoded.options.config_overrides
    assert overrides is not None
    assert "depth_tick_interval" not in overrides
    assert overrides["objects.depth.large_tick_interval"] == 12
    assert overrides["labels.linear.placement"] == "above_feature"
    assert overrides["canvas.linear.track_layout"] == current_layout

    rendered = render_request(decoded)
    assert rendered.output_paths == (
        tmp_path / "output" / f"legacy-{schema}-{legacy_layout}.svg",
    )

    latest = encode_canonical_request(decoded).payload
    latest_overrides = latest["diagramOptions"]["configOverrides"]
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert "depth_tick_interval" not in latest_overrides
    assert latest_overrides["objects.depth.large_tick_interval"] == 12
    assert latest_overrides["labels.linear.placement"] == "above_feature"
    assert latest_overrides["canvas.linear.track_layout"] == current_layout


@pytest.mark.parametrize(
    ("overrides", "retired"),
    [
        ({"depth_tick_interval": 12}, "depth_tick_interval"),
        ({"label_placement": "on_feature"}, "label_placement"),
        ({"linear_track_layout": "spreadout"}, "linear_track_layout"),
        ({"linear_track_layout": "tuckin"}, "linear_track_layout"),
    ],
)
def test_current_schema_rejects_retired_config_overrides(
    tmp_path: Path,
    overrides: dict[str, object],
    retired: str,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["configOverrides"] = overrides

    with pytest.raises(CanonicalRequestDecodingError, match=retired):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(
                encoded,
                tmp_path / "current-schema-retired-overrides",
            ),
            output_directory=tmp_path / "output",
        )


@pytest.mark.parametrize(
    ("overrides", "retired"),
    [
        ({"depth_tick_interval": 12}, "depth_tick_interval"),
        ({"label_placement": "on_feature"}, "label_placement"),
        ({"linear_track_layout": "spreadout"}, "linear_track_layout"),
        ({"linear_track_layout": "tuckin"}, "linear_track_layout"),
    ],
)
def test_current_schema_writer_rejects_retired_config_overrides(
    overrides: dict[str, object],
    retired: str,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    with pytest.raises(ValidationError, match=retired):
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(config_overrides=overrides),
        )


@pytest.mark.parametrize("schema", [1, 2])
def test_legacy_canonical_schema_preserves_repeat_rectangle_default(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(
                selected_features_set=("repeat_region",)
            ),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
    payload.pop("grouping", None)
    if schema == 1:
        for record_payload in payload["records"]:
            record_payload.pop("recordKey", None)
    payload["diagramOptions"].get("featureShapes", {}).pop("repeat_region", None)

    decoded = decode_canonical_request(
        payload,
        resource_paths=_materialize_resources(encoded, tmp_path / f"schema-{schema}"),
        output_directory=tmp_path / "output",
    )

    assert decoded.options.feature_shapes["repeat_region"] == "rectangle"


def test_current_canonical_schema_uses_underlay_default_and_round_trips_override(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    default_encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(
                selected_features_set=("repeat_region",)
            ),
        )
    )
    default_decoded = decode_canonical_request(
        default_encoded.payload,
        resource_paths=_materialize_resources(default_encoded, tmp_path / "default"),
        output_directory=tmp_path / "output-default",
    )
    assert resolve_feature_rendering(
        "repeat_region", default_decoded.options.feature_shapes
    ) == "underlay"

    explicit_encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=LinearDiagramOptions(
                selected_features_set=("repeat_region",),
                feature_shapes={"repeat_region": "underlay"},
            ),
        )
    )
    explicit_decoded = decode_canonical_request(
        explicit_encoded.payload,
        resource_paths=_materialize_resources(explicit_encoded, tmp_path / "explicit"),
        output_directory=tmp_path / "output-explicit",
    )
    assert explicit_decoded.options.feature_shapes == {"repeat_region": "underlay"}
    assert encode_canonical_request(explicit_decoded).payload == explicit_encoded.payload


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (lambda payload: payload.update(schema=6), "Unsupported canonical request schema"),
        (lambda payload: payload.update(mode="radial"), "Unsupported canonical request mode"),
        (lambda payload: payload.pop("output"), "Missing required field"),
        (lambda payload: payload.update(futureField=True), "Unknown field"),
        (
            lambda payload: payload["output"].update(outputDirectory="unsafe"),
            "Unknown field",
        ),
        (
            lambda payload: payload["records"][0]["source"].update(kind="remoteUrl"),
            "Unsupported record source kind",
        ),
    ],
)
def test_current_schema_rejects_unknown_or_invalid_structure(
    tmp_path: Path,
    mutator,
    message: str,
) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(records=(RecordInput(source=GenBankInputSource(source)),))
    )
    payload = copy.deepcopy(encoded.payload)
    mutator(payload)

    with pytest.raises(CanonicalRequestDecodingError, match=message):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(encoded, tmp_path),
            output_directory=tmp_path / "output",
        )


def test_decode_rejects_unknown_comparison_kind(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(records=(RecordInput(source=GenBankInputSource(source)),))
    )
    payload = copy.deepcopy(encoded.payload)
    payload["comparisons"].append({"kind": "futureComparison"})

    with pytest.raises(CanonicalRequestDecodingError, match="Unsupported comparison kind"):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(encoded, tmp_path),
            output_directory=tmp_path / "output",
        )


def test_decode_rejects_missing_resource_reference(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(records=(RecordInput(source=GenBankInputSource(source)),))
    )

    with pytest.raises(CanonicalRequestDecodingError, match="is missing"):
        decode_canonical_request(
            encoded.payload,
            resource_paths={},
            output_directory=tmp_path / "output",
        )


def test_encode_rejects_noncanonical_option_values(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    invalid_config = load_default_config()
    invalid_config["labels"]["filtering"]["extension"] = object()
    invalid_type = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=LinearDiagramOptions(window="10"),
    )
    invalid_json = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=LinearDiagramOptions(config=invalid_config),
    )
    empty_table = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=LinearDiagramOptions(feature_visibility_table=pd.DataFrame()),
    )

    with pytest.raises(CanonicalRequestEncodingError, match="typed contract"):
        encode_canonical_request(invalid_type)
    with pytest.raises(CanonicalRequestEncodingError, match="unsupported value type"):
        encode_canonical_request(invalid_json)
    with pytest.raises(CanonicalRequestEncodingError, match="at least one column"):
        encode_canonical_request(empty_table)


def test_schema_one_unknown_field_policy_is_explicit() -> None:
    assert UNKNOWN_FIELD_POLICY == "reject"
