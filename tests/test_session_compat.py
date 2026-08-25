from __future__ import annotations

import base64
import json
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
)
from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.linear import linear_main
import gbdraw.cli_utils.session as cli_session_module
from gbdraw.api.request_render import (
    CurrentRequestArtifacts,
    PreparedDiagramRequest,
    RequestRenderResult,
)
from gbdraw.api.requests import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.api.session_compat import (
    SessionCompatibleRequestRenderResult,
    adapt_session_request,
    canonical_payload_for_session_decode,
    render_session_compatible_request,
    rewrite_protein_artifact_references,
)
from gbdraw.exceptions import ValidationError
from gbdraw.session import (
    build_session_document,
    load_session_document,
    materialize_session,
    render_session,
    session_to_request,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA

_RELEASED_SCHEMA_V2_SESSION = (
    Path(__file__).parent
    / "fixtures"
    / "sessions"
    / "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz"
)
_VERSION_39_SESSION = (
    Path(__file__).parent
    / "fixtures"
    / "sessions"
    / "BGC0000708-BGC0000713.v39.gbdraw-session.json.gz"
)
_WSSV_CURRENT_SESSION = (
    Path(__file__).parent.parent
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "WSSV_genome_comparison.gbdraw-session.json"
)


def test_version_39_multiline_conservation_labels_expand_at_compat_boundary() -> None:
    payload = {
        "diagramOptions": {
            "conservationBlastFiles": [
                {"resourceId": "comparison-1", "representation": "file"},
                {"resourceId": "comparison-2", "representation": "file"},
            ],
            "conservationLabels": ["First comparison\nSecond comparison"],
        }
    }

    migrated = canonical_payload_for_session_decode(39, payload)

    assert migrated["diagramOptions"]["conservationLabels"] == [
        "First comparison",
        "Second comparison",
    ]
    assert payload["diagramOptions"]["conservationLabels"] == [
        "First comparison\nSecond comparison"
    ]


@pytest.mark.parametrize("version", (33, CURRENT_SESSION_VERSION))
def test_multiline_conservation_label_migration_is_version_39_only(
    version: int,
) -> None:
    labels = ["First comparison\nSecond comparison"]
    payload = {
        "diagramOptions": {
            "conservationBlastFiles": [
                {"resourceId": "comparison-1", "representation": "file"},
                {"resourceId": "comparison-2", "representation": "file"},
            ],
            "conservationLabels": labels,
        }
    }

    migrated = canonical_payload_for_session_decode(version, payload)

    assert migrated["diagramOptions"]["conservationLabels"] == labels


def test_version_39_multiline_conservation_label_count_mismatch_stays_strict() -> None:
    labels = ["First comparison\nSecond comparison"]
    payload = {
        "diagramOptions": {
            "conservationBlastFiles": [
                {"resourceId": "comparison-1", "representation": "file"},
                {"resourceId": "comparison-2", "representation": "file"},
                {"resourceId": "comparison-3", "representation": "file"},
            ],
            "conservationLabels": labels,
        }
    }

    migrated = canonical_payload_for_session_decode(39, payload)

    assert migrated["diagramOptions"]["conservationLabels"] == labels


def _linear_request(
    tmp_path: Path,
    *,
    output_prefix: str = "session-compatible",
) -> LinearDiagramRequest:
    record = SeqRecord(
        Seq("ATGCGCAT"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    return LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key="record-1",
            ),
        ),
        output=RenderOutputRequest(
            output_prefix=output_prefix,
            output_directory=tmp_path,
        ),
    )


def test_session_adapter_passes_plain_prepared_request_to_current_renderer(
    tmp_path: Path,
    monkeypatch,
) -> None:
    request = _linear_request(tmp_path)
    session_artifacts = build_session_document(request).to_dict()
    captured: list[type[PreparedDiagramRequest]] = []

    def fake_render(
        prepared: PreparedDiagramRequest,
        **_kwargs,
    ) -> RequestRenderResult:
        captured.append(type(prepared))
        return RequestRenderResult(
            mode=prepared.mode,
            request=prepared.request,
            records=prepared.records,
            drawing=prepared.drawing,
            output_paths=(),
            linear_metadata=prepared.linear_metadata,
            losat_cache_entries=prepared.losat_cache_entries,
            losat_derived_cache_entries=prepared.losat_derived_cache_entries,
            protein_identity_manifest=prepared.protein_identity_manifest,
        )

    monkeypatch.setattr(
        "gbdraw.api.session_compat.render_prepared_request",
        fake_render,
    )

    result = render_session_compatible_request(request, session_artifacts)

    assert captured == [PreparedDiagramRequest]
    assert isinstance(result, SessionCompatibleRequestRenderResult)


def _released_canonical_session(
    tmp_path: Path,
    *,
    version: int,
    request_schema: int,
) -> dict[str, object]:
    data = build_session_document(
        _linear_request(tmp_path, output_prefix=f"released-v{version}")
    ).to_dict()
    data["version"] = version
    payload = data["renderRequest"]
    assert isinstance(payload, dict)
    payload["schema"] = request_schema
    payload.pop("grouping", None)
    records = payload["records"]
    assert isinstance(records, list)
    for record in records:
        assert isinstance(record, dict)
        record.pop("cardinality", None)
    if request_schema == 1:
        records[0].pop("recordKey", None)
    nested_output = payload["diagramOptions"].get("output")
    if isinstance(nested_output, dict):
        nested_output["outputPrefix"] = "ignored-legacy-prefix"
    return data


def _released_cli_session(
    tmp_path: Path,
    *,
    version: int,
) -> dict[str, object]:
    source_path = tmp_path / f"released-v{version}.gb"
    record = SeqRecord(
        Seq("ATGCGCAT"),
        id=f"released-v{version}",
        annotations={"molecule_type": "DNA"},
    )
    SeqIO.write(record, source_path, "genbank")
    source_bytes = source_path.read_bytes()
    embedded = {
        "name": source_path.name,
        "type": "application/octet-stream",
        "size": len(source_bytes),
        "lastModified": 0,
        "data": base64.b64encode(source_bytes).decode("ascii"),
    }
    return {
        "format": "gbdraw-session",
        "version": version,
        "createdAt": "2026-07-30T00:00:00Z",
        "config": {
            "form": {"prefix": f"released-v{version}"},
            "adv": {},
        },
        "ui": {"mode": "linear", "lInputType": "gb"},
        "files": {"linearSeqs": [{"gb": embedded}]},
        "cliInvocation": {
            "schema": 1,
            "mode": "linear",
            "args": [
                "--gbk",
                source_path.name,
                "--format",
                "svg",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ],
            "renderFormats": ["svg"],
            "fileBindings": [
                {
                    "argIndex": 1,
                    "slot": "files.linearSeqs[0].gb",
                    "name": source_path.name,
                }
            ],
            "generatedBy": "gbdraw",
        },
    }


def _orthogroup_result_with_reference(reference: str) -> OrthogroupResult:
    member = OrthogroupMember(
        orthogroup_id="OG1",
        protein_id=reference,
        record_index=0,
        feature_index=0,
        record_id="record",
        label="Protein",
        start=0,
        end=8,
        strand=1,
        feature_svg_id="feature-1",
        source_protein_id=reference,
    )
    return OrthogroupResult(
        orthogroups={"OG1": [member]},
        member_by_protein_id={reference: member},
    )


def _collinearity_result_with_reference(reference: str) -> CollinearityResult:
    anchor = CollinearityAnchor(
        query_protein_id=reference,
        subject_protein_id="current-protein",
        query_record_index=0,
        subject_record_index=0,
        query_order=0,
        subject_order=1,
        query_start=0,
        query_end=3,
        subject_start=4,
        subject_end=7,
        identity=95.0,
        evalue=1e-20,
        bitscore=100.0,
        alignment_length=3,
        query_feature_svg_id="feature-1",
        subject_feature_svg_id="feature-2",
        source="precomputed",
        query_unit_id="unit-1",
        subject_unit_id="unit-2",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name="Query",
        subject_display_name="Subject",
    )
    return CollinearityResult(
        blocks=(
            CollinearityBlock(
                block_id="block-1",
                query_record_index=0,
                subject_record_index=0,
                orientation="plus",
                score=100.0,
                anchors=(anchor,),
            ),
        ),
    )


def test_unresolved_session_batch_preflights_all_resolved_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = tuple(
        SeqRecord(
            Seq("ATGCGCAT"),
            id=record_id,
            annotations={"molecule_type": "DNA"},
        )
        for record_id in ("first", "second")
    )
    output_prefix = tmp_path / "diagram"
    request = CircularBatchRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=f"record-{index}",
            )
            for index, record in enumerate(records, start=1)
        ),
        output_policy=CircularBatchOutputPolicy(output_prefix=output_prefix),
    )
    session_artifacts = build_session_document(request).to_dict()
    second_output = tmp_path / "diagram_2.svg"
    second_output.write_text("occupied", encoding="utf-8")
    monkeypatch.setattr(
        "gbdraw.api.session_compat.build_request_plan_diagram",
        lambda _plan, **_kwargs: pytest.fail(
            "resolved session batch outputs must be preflighted before building"
        ),
    )

    with pytest.raises(ValidationError, match="already exist"):
        render_session_compatible_request(request, session_artifacts)

    assert not (tmp_path / "diagram_1.svg").exists()
    assert second_output.read_text(encoding="utf-8") == "occupied"


def test_rewrite_protein_artifact_references_rejects_key_collisions() -> None:
    with pytest.raises(ValidationError, match="duplicate key"):
        rewrite_protein_artifact_references(
            {"legacy-a": 1, "legacy-b": 2},
            {"legacy-a": "h_same", "legacy-b": "h_same"},
        )


def test_rewrite_protein_artifact_references_updates_compound_ids() -> None:
    query_id = "p_r_old_0_9_1_deadbeefdead"
    subject_id = "p_r_other_10_19_1_cafebabecafe"

    rewritten = rewrite_protein_artifact_references(
        {
            "query_protein_id": f"{query_id};{subject_id}",
            "supporting_edge": f"{query_id}->{subject_id}:rbh",
        },
        {
            query_id: "h_aaaaaaaaaaaaaaaaaaaaaaaaaa",
            subject_id: "h_bbbbbbbbbbbbbbbbbbbbbbbbbb",
        },
    )

    assert rewritten == {
        "query_protein_id": (
            "h_aaaaaaaaaaaaaaaaaaaaaaaaaa;h_bbbbbbbbbbbbbbbbbbbbbbbbbb"
        ),
        "supporting_edge": (
            "h_aaaaaaaaaaaaaaaaaaaaaaaaaa->h_bbbbbbbbbbbbbbbbbbbbbbbbbb:rbh"
        ),
    }


@pytest.mark.parametrize("compound", [False, True], ids=["exact", "compound"])
@pytest.mark.parametrize("owner", ["comparison", "alignment-target"])
def test_feature_analysis_ids_fail_closed_across_protein_request_artifacts(
    owner: str,
    compound: bool,
    tmp_path: Path,
) -> None:
    base_request = _linear_request(tmp_path)
    session = build_session_document(base_request).to_dict()
    session["version"] = 33
    feature_analysis_id = f"f_{'a' * 64}"
    value = (
        f"record@instance|alias~{feature_analysis_id}"
        if compound
        else feature_analysis_id
    )
    if owner == "comparison":
        options = replace(
            base_request.options,
            protein_comparisons=(
                DataFrame({"query_protein_id": [value]}),
            ),
        )
    else:
        options = replace(
            base_request.options,
            align_orthogroup_feature=value,
        )
    request = replace(base_request, options=options)

    with pytest.raises(
        ValidationError,
        match="no verified session artifact resolved",
    ):
        adapt_session_request(request, session)


@pytest.mark.parametrize("compound", [False, True], ids=["exact", "compound"])
@pytest.mark.parametrize(
    "owner",
    ["orthogroup-result", "collinearity-result"],
)
def test_typed_protein_results_fail_closed_for_unresolved_analysis_ids(
    owner: str,
    compound: bool,
    tmp_path: Path,
) -> None:
    base_request = _linear_request(tmp_path)
    session = _released_canonical_session(
        tmp_path,
        version=33,
        request_schema=2,
    )
    feature_analysis_id = f"f_{'b' * 64}"
    reference = (
        f"record@instance|alias~{feature_analysis_id}"
        if compound
        else feature_analysis_id
    )
    if owner == "orthogroup-result":
        options = replace(
            base_request.options,
            orthogroups=_orthogroup_result_with_reference(reference),
        )
    else:
        options = replace(
            base_request.options,
            collinearity_blocks=_collinearity_result_with_reference(reference),
        )
    request = replace(base_request, options=options)

    with pytest.raises(
        ValidationError,
        match="no verified session artifact resolved",
    ):
        adapt_session_request(request, session)


@pytest.mark.parametrize(
    ("version", "request_schema"),
    (
        (27, None),
        (28, None),
        (29, None),
        (30, None),
        (31, 1),
        (32, 2),
        (33, 2),
    ),
    ids=lambda value: f"v{value}" if isinstance(value, int) else "cli",
)
def test_released_session_versions_write_current_reloadable_sidecars(
    tmp_path: Path,
    version: int,
    request_schema: int | None,
) -> None:
    source_session = (
        _released_cli_session(tmp_path, version=version)
        if request_schema is None
        else _released_canonical_session(
            tmp_path,
            version=version,
            request_schema=request_schema,
        )
    )
    session_path = tmp_path / f"released-v{version}.gbdraw-session.json"
    session_path.write_text(json.dumps(source_session), encoding="utf-8")
    output_prefix = tmp_path / f"rendered-v{version}"
    sidecar_path = tmp_path / f"rendered-v{version}.gbdraw-session.json.gz"

    linear_main(
        [
            "--session",
            str(session_path),
            "--output",
            str(output_prefix),
            "--format",
            "svg",
            "--session_output",
            str(sidecar_path),
        ]
    )

    assert output_prefix.with_suffix(".svg").is_file()
    saved = load_session_document(sidecar_path)
    assert saved.version == CURRENT_SESSION_VERSION
    assert saved.to_dict()["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA

    replay_directory = tmp_path / "replay"
    with materialize_session(saved, output_directory=replay_directory) as materialized:
        replayed = render_session(materialized)

    assert replayed.output_paths == (
        replay_directory / f"rendered-v{version}.svg",
    )
    assert replayed.output_paths[0].is_file()


def test_supported_v32_session_replays_through_compatibility_adapter(
    tmp_path: Path,
) -> None:
    legacy_raw = {
        "schema": 2,
        "kind": "raw-losat",
        "key": "legacy-protein-key",
        "text": "p_r_old\tp_r_other\n",
        "program": "blastp",
        "outfmt": "6",
        "args": [],
        "queryCanonicalHash": "old-query",
        "subjectCanonicalHash": "old-subject",
    }
    legacy_derived = {
        "schema": 1,
        "kind": "derived-losatp-payload",
        "key": "legacy-derived-key",
        "mode": "orthogroup",
        "payload": {"groups": []},
    }
    data = build_session_document(_linear_request(tmp_path)).to_dict()
    data["version"] = 32
    data["losatCache"] = {"entries": [legacy_raw]}
    data["losatDerivedCache"] = {"entries": [legacy_derived]}
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        result = render_session(materialized)

    assert isinstance(result, SessionCompatibleRequestRenderResult)
    assert result.output_paths == (tmp_path / "session-compatible.svg",)
    assert result.output_paths[0].is_file()
    assert result.legacy_protein_raw_candidates
    assert result.legacy_protein_raw_candidates[0]["state"] == "rejected"
    assert result.legacy_protein_derived_evidence == (legacy_derived,)
    assert result.protein_id_map == {}


def test_released_schema_v2_fixture_promotes_to_current_typed_artifacts(
    tmp_path: Path,
) -> None:
    document = load_session_document(_RELEASED_SCHEMA_V2_SESSION)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        adapted = adapt_session_request(
            session_to_request(materialized),
            document.to_dict(),
        )

    protein_entries = tuple(
        entry
        for entry in adapted.artifacts.losat_cache_entries
        if entry.get("identityKind") == "protein"
    )
    assert document.version == 33
    assert len(protein_entries) == 25
    assert {entry["schema"] for entry in protein_entries} == {4}
    assert adapted.migration_report.protein_id_map
    assert adapted.migration_report.protein_raw_candidates == ()


def test_released_schema_v2_alignment_target_promotes_with_protein_artifacts(
    tmp_path: Path,
) -> None:
    document = load_session_document(_RELEASED_SCHEMA_V2_SESSION)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        comparisons = request.options.protein_comparisons
        orthogroups = request.options.orthogroups
        assert comparisons
        assert orthogroups is not None
        legacy_target = str(comparisons[0].iloc[0]["query_protein_id"])
        request = replace(
            request,
            options=replace(
                request.options,
                align_orthogroup_feature=legacy_target,
            ),
        )

        adapted = adapt_session_request(request, document.to_dict())

    current_target = adapted.migration_report.protein_id_map[legacy_target]
    adapted_options = adapted.request.options
    assert adapted_options.align_orthogroup_feature == current_target
    assert adapted_options.orthogroups is not None
    assert current_target in adapted_options.orthogroups.member_by_protein_id
    assert adapted_options.protein_comparisons
    assert (
        adapted_options.protein_comparisons[0].iloc[0]["query_protein_id"]
        == current_target
    )


def test_released_schema_v2_fixture_cli_sidecar_is_current_and_rerenders(
    tmp_path: Path,
) -> None:
    output_prefix = tmp_path / "released-replay"
    sidecar_path = tmp_path / "released-replay.gbdraw-session.json.gz"

    linear_main(
        [
            "--session",
            str(_RELEASED_SCHEMA_V2_SESSION),
            "--output",
            str(output_prefix),
            "--format",
            "svg",
            "--session_output",
            str(sidecar_path),
        ]
    )

    assert output_prefix.with_suffix(".svg").is_file()
    saved = load_session_document(sidecar_path)
    payload = saved.to_dict()
    assert saved.version == CURRENT_SESSION_VERSION
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert "depth_tick_interval" not in payload["config"]["adv"]
    assert all(
        "tick_interval" not in track
        for track in payload["config"]["adv"]["depth_tracks"]
    )
    assert (
        "collinearMaxGeneGap"
        not in payload["config"]["losat"]["blastp"]
    )

    rerendered_prefix = tmp_path / "released-rerender"
    linear_main(
        [
            "--session",
            str(sidecar_path),
            "--output",
            str(rerendered_prefix),
            "--format",
            "svg",
        ]
    )
    assert rerendered_prefix.with_suffix(".svg").is_file()


def test_version_39_typed_replay_retains_dormant_comparison_resource(
    tmp_path: Path,
) -> None:
    source_document = load_session_document(_VERSION_39_SESSION)
    source_payload = source_document.to_dict()
    record_uids = [
        str(row["uid"])
        for row in source_payload["config"]["linearRecordLayout"]["rows"]
    ]
    sequence_files = []
    for uid, record in zip(
        record_uids,
        source_payload["renderRequest"]["records"],
        strict=True,
    ):
        resource_id = str(record["source"]["resourceId"])
        entry = dict(source_payload["resources"][resource_id])
        entry["name"] = f"{uid}.gb"
        sequence_files.append({"uid": uid, "gb": entry})
    dormant_bytes = b"query\tsubject\t97\n"
    dormant_file = {
        "name": "retained-non-adjacent.tsv",
        "type": "text/tab-separated-values",
        "size": len(dormant_bytes),
        "lastModified": 0,
        "encoding": "base64",
        "data": base64.b64encode(dormant_bytes).decode("ascii"),
    }
    migration_source = {
        "config": source_payload["config"],
        "ui": source_payload["ui"],
        "webFiles": source_payload["webFiles"],
        "resources": source_payload["resources"],
        "files": {
            "linearSeqs": sequence_files,
            "linearComparisons": [
                {
                    "id": "dormant-v39-upload",
                    "queryUid": record_uids[0],
                    "subjectUid": record_uids[2],
                    "source": "upload",
                    "file": dormant_file,
                }
            ],
        },
    }
    adjunct, web_file_inventory = (
        cli_session_module._project_session_adjunct_for_current_write(
            migration_source,
            source_version=39,
        )
    )

    with materialize_session(source_document, output_directory=tmp_path) as materialized:
        rewritten = build_session_document(
            session_to_request(materialized),
            adjunct=adjunct,
            web_file_inventory=web_file_inventory,
        ).to_dict()

    assert rewritten["config"]["linearComparisonPlan"]["edges"] == [
        {
            "id": "dormant-v39-upload",
            "queryUid": record_uids[0],
            "subjectUid": record_uids[2],
            "included": False,
            "fileActive": False,
            "losatFilenameActive": False,
            "source": "upload",
            "losatFilename": "",
        }
    ]
    comparison_bindings = rewritten["webFiles"]["bindings"]["linearComparisons"]
    assert [set(binding) for binding in comparison_bindings] == [{"id", "file"}]
    retained_resource_id = comparison_bindings[0]["file"]["resourceId"]
    assert base64.b64decode(
        rewritten["resources"][retained_resource_id]["data"]
    ) == dormant_bytes
    assert sum(
        base64.b64decode(resource["data"]) == dormant_bytes
        for resource in rewritten["resources"].values()
    ) == 1


def test_current_typed_replay_retains_web_only_losat_fastas(
    tmp_path: Path,
) -> None:
    source_document = load_session_document(_WSSV_CURRENT_SESSION)
    source_payload = source_document.to_dict()
    source_web_files = source_payload["webFiles"]
    source_resources = source_payload["resources"]
    source_ids = source_web_files["conservationLosatFastaSources"]
    expected_data = [source_resources[resource_id]["data"] for resource_id in source_ids]
    expected_names = [
        source_web_files["resourceOriginalNames"][resource_id]
        for resource_id in source_ids
    ]

    adjunct, web_file_inventory = (
        cli_session_module._project_session_adjunct_for_current_write(
            source_payload,
            source_version=source_document.version,
        )
    )
    assert web_file_inventory is not None

    with materialize_session(source_document, output_directory=tmp_path) as materialized:
        rewritten = build_session_document(
            session_to_request(materialized),
            adjunct=adjunct,
            web_file_inventory=web_file_inventory,
        ).to_dict()

    rewritten_web_files = rewritten["webFiles"]
    rewritten_resources = rewritten["resources"]
    rewritten_ids = rewritten_web_files["conservationLosatFastaSources"]
    fasta_bindings = rewritten_web_files["bindings"]["c_conservation_fastas"]
    assert rewritten_ids == [binding["resourceId"] for binding in fasta_bindings]
    assert [rewritten_resources[resource_id]["data"] for resource_id in rewritten_ids] == (
        expected_data
    )
    assert [
        rewritten_web_files["resourceOriginalNames"][resource_id]
        for resource_id in rewritten_ids
    ] == expected_names


def test_released_schema_v2_fixture_sidecar_collision_is_atomic(
    tmp_path: Path,
) -> None:
    output_prefix = tmp_path / "released-replay"
    sidecar_path = tmp_path / "occupied.gbdraw-session.json"
    sidecar_path.write_text("keep this session", encoding="utf-8")

    with pytest.raises(ValidationError, match="Session output already exists"):
        linear_main(
            [
                "--session",
                str(_RELEASED_SCHEMA_V2_SESSION),
                "--output",
                str(output_prefix),
                "--format",
                "svg",
                "--session_output",
                str(sidecar_path),
            ]
        )

    assert sidecar_path.read_text(encoding="utf-8") == "keep this session"
    assert not output_prefix.with_suffix(".svg").exists()


def test_session_adapter_rejects_unsupported_session_schema(tmp_path: Path) -> None:
    request = _linear_request(tmp_path)
    data = build_session_document(request).to_dict()
    data["version"] = 38

    with pytest.raises(ValidationError, match="Unsupported session version"):
        adapt_session_request(request, data)


def test_current_artifact_type_rejects_legacy_cache_schema() -> None:
    with pytest.raises(ValidationError, match="Unsupported current LOSAT artifact"):
        CurrentRequestArtifacts(
            losat_cache_entries=(
                {
                    "schema": 2,
                    "kind": "raw-losat",
                    "key": "legacy-protein-key",
                    "text": "",
                    "program": "blastp",
                },
            )
        )


def test_released_noncanonical_linear_cli_replay_promotes_current_sidecar(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "legacy-source.gb"
    record = SeqRecord(
        Seq("ATGCGCAT"),
        id="legacy-record",
        annotations={"molecule_type": "DNA"},
    )
    SeqIO.write(record, source_path, "genbank")
    source_bytes = source_path.read_bytes()
    embedded = {
        "name": source_path.name,
        "type": "application/octet-stream",
        "size": len(source_bytes),
        "lastModified": 0,
        "data": base64.b64encode(source_bytes).decode("ascii"),
    }
    legacy_raw = {
        "schema": 2,
        "kind": "raw-losat",
        "key": "legacy-protein-key",
        "text": "p_r_old\tp_r_other\n",
        "program": "blastp",
        "outfmt": "6",
        "args": [],
        "queryCanonicalHash": "old-query",
        "subjectCanonicalHash": "old-subject",
    }
    session = {
        "format": "gbdraw-session",
        "version": 30,
        "createdAt": "2026-07-30T00:00:00Z",
        "config": {
            "form": {"prefix": "legacy"},
            "adv": {
                "depth_tick_interval": 10,
                "depth_tracks": [{"tick_interval": 5}],
            },
            "losat": {"blastp": {"collinearMaxGeneGap": 2}},
        },
        "ui": {"mode": "linear", "lInputType": "gb"},
        "files": {"linearSeqs": [{"gb": embedded}]},
        "losatCache": {"entries": [legacy_raw]},
        "cliInvocation": {
            "schema": 1,
            "mode": "linear",
            "args": [
                "--gbk",
                source_path.name,
                "--format",
                "svg",
                "--no-gc",
                "--no-skew",
                "--legend",
                "none",
            ],
            "renderFormats": ["svg"],
            "fileBindings": [
                {
                    "argIndex": 1,
                    "slot": "files.linearSeqs[0].gb",
                    "name": source_path.name,
                }
            ],
            "generatedBy": "gbdraw",
        },
    }
    session_path = tmp_path / "legacy-v30.gbdraw-session.json"
    session_path.write_text(json.dumps(session), encoding="utf-8")
    output_prefix = tmp_path / "legacy-replay"
    sidecar_path = tmp_path / "legacy-replay.gbdraw-session.json.gz"

    linear_main(
        [
            "--session",
            str(session_path),
            "--output",
            str(output_prefix),
            "--format",
            "svg",
            "--session_output",
            str(sidecar_path),
        ]
    )

    assert output_prefix.with_suffix(".svg").is_file()
    saved = load_session_document(sidecar_path)
    payload = saved.to_dict()
    assert saved.version == CURRENT_SESSION_VERSION
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert "depth_tick_interval" not in payload["config"]["adv"]
    assert payload["config"]["adv"]["depth_large_tick_interval"] == 10
    assert payload["config"]["adv"]["depth_tracks"] == [
        {"large_tick_interval": 5}
    ]
    assert (
        payload["config"]["losat"]["blastp"]["collinearMaxUnitGap"]
        == 2
    )
    assert (
        "collinearMaxGeneGap"
        not in payload["config"]["losat"]["blastp"]
    )

    rerendered_prefix = tmp_path / "legacy-rerender"
    linear_main(
        [
            "--session",
            str(sidecar_path),
            "--output",
            str(rerendered_prefix),
            "--format",
            "svg",
        ]
    )
    assert rerendered_prefix.with_suffix(".svg").is_file()
