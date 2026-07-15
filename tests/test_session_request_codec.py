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
from gbdraw.api.options import (
    CircularMultiRecordOptions,
    ColorOptions,
    DiagramOptions,
    TrackOptions,
)
from gbdraw.api.config import load_default_config
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
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    UNKNOWN_FIELD_POLICY,
    CanonicalRequestDecodingError,
    CanonicalRequestEncodingError,
    EncodedCanonicalRequest,
    decode_canonical_request,
    encode_canonical_request,
)
from gbdraw.tracks import CircularTrackSlot, ScalarSpec


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
        ),
        layout=CircularMultiRecordOptions(
            multi_record_size_mode="sqrt",
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
    assert encode_canonical_request(decoded).payload == encoded.payload


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
    ("mutator", "message"),
    [
        (lambda payload: payload.update(schema=2), "Unsupported canonical request schema"),
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
def test_schema_one_rejects_unknown_or_invalid_structure(
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
