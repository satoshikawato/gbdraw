from __future__ import annotations

import ast
import copy
import json
from pathlib import Path

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
    CircularMultiRecordOptions,
    ColorOptions,
    DiagramOptions,
    OutputOptions,
    TrackOptions,
)
from gbdraw.api.config import load_default_config
from gbdraw.api.request_render import render_request
from gbdraw.api.requests import (
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
from gbdraw.features.shapes import resolve_feature_rendering
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
    UNKNOWN_FIELD_POLICY,
    CanonicalRequestDecodingError,
    CanonicalRequestEncodingError,
    EncodedCanonicalRequest,
    decode_canonical_request,
    encode_canonical_request,
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


def _table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_type": ["CDS", "tRNA"],
            "caption": ["coding sequence", "transfer RNA"],
        }
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
        options=DiagramOptions(
            colors=ColorOptions(
                color_table=table,
                default_colors_palette="soft",
            ),
            tracks=TrackOptions(
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
            output=OutputOptions(legend="left", plot_title_position="top"),
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
    assert decoded.options.output == OutputOptions(
        legend="left",
        plot_title_position="top",
    )
    assert decoded.output.output_prefix == "canonical-circular"
    assert decoded.records[1].region == parse_region_spec("1-100:rc")
    assert encode_canonical_request(decoded).payload == encoded.payload


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
        options=DiagramOptions(
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
            output=OutputOptions(legend="bottom", plot_title_position="bottom"),
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


def test_schema3_nested_output_prefix_is_required_but_ignored(
    tmp_path: Path,
) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=DiagramOptions(
                output=OutputOptions(legend="left", plot_title_position="top"),
            ),
            output=RenderOutputRequest(output_prefix="outer-authority"),
        )
    )
    legacy_payload = copy.deepcopy(encoded.payload)
    legacy_payload["schema"] = 3
    legacy_payload["diagramOptions"]["output"]["outputPrefix"] = "ignored-nested"

    decoded = decode_canonical_request(
        legacy_payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "legacy"),
        output_directory=tmp_path / "output",
    )

    assert decoded.output.output_prefix == "outer-authority"
    assert decoded.options.output == OutputOptions(
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


def test_schema4_rejects_legacy_nested_output_prefix(tmp_path: Path) -> None:
    source = _source_file(tmp_path / "record.gbk")
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source)),),
            options=DiagramOptions(output=OutputOptions(legend="none")),
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


@pytest.mark.parametrize("schema", [1, 2, 3])
def test_legacy_circular_schemas_migrate_removed_layout_values(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=DiagramOptions(
                tracks=TrackOptions(
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
            options=DiagramOptions(
                tracks=TrackOptions(
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


def test_schema4_rejects_private_circular_track_params(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=DiagramOptions(
                tracks=TrackOptions(
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
                tmp_path / "schema4-private-slot",
            ),
            output_directory=tmp_path / "output",
        )


def test_legacy_factor_spacing_replays_but_cannot_leak_into_schema4(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = 3
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


def test_schema3_migrates_removed_feature_table_field(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    table = pd.DataFrame(
        [{"record_id": "record", "feature_type": "CDS", "visible": True}]
    )
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=DiagramOptions(feature_visibility_table=table),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = 3
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


def test_schema4_rejects_removed_feature_table_field(
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
        options=DiagramOptions(
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
    config = GbdrawConfig.from_dict(load_default_config())
    request = CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=DiagramOptions(
            config=config,
            colors=ColorOptions(default_colors_file=str(table_file)),
            feature_visibility_table_file=str(table_file),
            label_whitelist_file=str(table_file),
            qualifier_priority_file=str(table_file),
            label_override_file=str(table_file),
            depth_track_files=((str(depth_file), None),),
            conservation_blast_files=(str(blast_file),),
        ),
    )

    encoded = encode_canonical_request(request)
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "materialized-files"),
        output_directory=tmp_path / "output",
    )

    assert isinstance(decoded.options.config, dict)
    assert decoded.options.feature_visibility_table_file == str(table_file)
    assert decoded.options.depth_track_files == ((str(depth_file), None),)
    assert decoded.options.conservation_blast_files == (str(blast_file),)
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
        options=DiagramOptions(annotations=annotations),
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
                    "show_gc": True,
                    "show_skew": True,
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
                    "show_gc": False,
                    "show_skew": False,
                    "linear_axis_stroke_color": "lightgray",
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


@pytest.mark.parametrize("schema", [1, 2, 3])
@pytest.mark.parametrize(
    ("request_type", "expected_overrides"),
    [
        (
            CircularDiagramRequest,
            {
                "show_gc": True,
                "show_skew": True,
            },
        ),
        (
            LinearDiagramRequest,
            {
                "show_gc": True,
                "show_skew": True,
                "linear_axis_stroke_color": "gray",
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
def test_schema4_sparse_options_use_current_mode_defaults(
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
        "show_gc": False,
        "show_skew": False,
        "linear_axis_stroke_color": "hotpink",
    }


@pytest.mark.parametrize("schema", [1, 2, 3])
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
    assert overrides["depth_large_tick_interval"] == 12
    assert overrides["label_placement"] == "above_feature"
    assert overrides["linear_track_layout"] == current_layout

    rendered = render_request(decoded)
    assert rendered.output_paths == (
        tmp_path / "output" / f"legacy-{schema}-{legacy_layout}.svg",
    )

    latest = encode_canonical_request(decoded).payload
    latest_overrides = latest["diagramOptions"]["configOverrides"]
    assert latest["schema"] == CANONICAL_REQUEST_SCHEMA
    assert "depth_tick_interval" not in latest_overrides
    assert latest_overrides["depth_large_tick_interval"] == 12
    assert latest_overrides["label_placement"] == "above_feature"
    assert latest_overrides["linear_track_layout"] == current_layout


@pytest.mark.parametrize(
    ("overrides", "retired"),
    [
        ({"depth_tick_interval": 12}, "depth_tick_interval"),
        ({"label_placement": "on_feature"}, "label_placement=on_feature"),
        ({"linear_track_layout": "spreadout"}, "linear_track_layout=spreadout"),
        ({"linear_track_layout": "tuckin"}, "linear_track_layout=tuckin"),
    ],
)
def test_schema4_rejects_retired_config_overrides(
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
                tmp_path / "schema4-retired-overrides",
            ),
            output_directory=tmp_path / "output",
        )


@pytest.mark.parametrize(
    ("overrides", "retired"),
    [
        ({"depth_tick_interval": 12}, "depth_tick_interval"),
        ({"label_placement": "on_feature"}, "label_placement=on_feature"),
        ({"linear_track_layout": "spreadout"}, "linear_track_layout=spreadout"),
        ({"linear_track_layout": "tuckin"}, "linear_track_layout=tuckin"),
    ],
)
def test_schema4_writer_rejects_retired_config_overrides(
    overrides: dict[str, object],
    retired: str,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=DiagramOptions(config_overrides=overrides),
    )

    with pytest.raises(CanonicalRequestEncodingError, match=retired):
        encode_canonical_request(request)


@pytest.mark.parametrize("schema", [1, 2])
def test_legacy_canonical_schema_preserves_repeat_rectangle_default(
    tmp_path: Path,
    schema: int,
) -> None:
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    encoded = encode_canonical_request(
        LinearDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=DiagramOptions(selected_features_set=("repeat_region",)),
        )
    )
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = schema
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
            options=DiagramOptions(selected_features_set=("repeat_region",)),
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
            options=DiagramOptions(
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
        (lambda payload: payload.update(schema=5), "Unsupported canonical request schema"),
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
    invalid_type = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=DiagramOptions(window="10"),
    )
    invalid_json = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=DiagramOptions(config_overrides={"bad": object()}),
    )
    empty_table = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        options=DiagramOptions(feature_visibility_table=pd.DataFrame()),
    )

    with pytest.raises(CanonicalRequestEncodingError, match="typed contract"):
        encode_canonical_request(invalid_type)
    with pytest.raises(CanonicalRequestEncodingError, match="unsupported value type"):
        encode_canonical_request(invalid_json)
    with pytest.raises(CanonicalRequestEncodingError, match="at least one column"):
        encode_canonical_request(empty_table)


def test_schema_one_unknown_field_policy_is_explicit() -> None:
    assert UNKNOWN_FIELD_POLICY == "reject"
