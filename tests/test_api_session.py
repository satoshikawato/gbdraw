from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
)
from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.api import (
    CircularDiagramRequest,
    DepthTrackInput,
    InMemoryRecordSource,
    LinearDiagramRequest,
    LinearDiagramOptions,
    LinearTrackOptions,
    LinearTrackSlot,
    RecordInput,
    RenderOutputRequest,
    ScalarSpec,
    SessionFormatError,
    SessionRenderError,
    SessionResourceError,
    SessionVersionError,
    build_session_document,
    load_session_document,
    materialize_session,
    render_session,
    save_session_document,
    session_to_request,
    normalize_request_records,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION


def _record(record_id: str = "record") -> RecordInput:
    seqrecord = SeqRecord(
        Seq("ATGCGCAT"),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )
    return RecordInput(source=InMemoryRecordSource(seqrecord))


@pytest.mark.parametrize(
    ("request_type", "mode"),
    ((CircularDiagramRequest, "circular"), (LinearDiagramRequest, "linear")),
)
def test_session_document_round_trip_owns_resource_lifetime(
    tmp_path: Path,
    request_type,
    mode: str,
) -> None:
    request = request_type(
        records=(_record(),),
        output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
    )
    document = build_session_document(request, title="round-trip")

    assert document.version == CURRENT_SESSION_VERSION
    assert document.mode == mode
    assert document.has_canonical_request is True
    assert document.to_dict()["resources"]["record-1-genbank"]["encoding"] == "base64"
    assert document.to_dict()["renderRequest"]["output"]["overwrite"] is False

    with materialize_session(
        document,
        output_directory=tmp_path,
        temporary_directory=tmp_path / "materialized",
    ) as materialized:
        decoded = session_to_request(materialized)
        assert decoded.output.overwrite is False
        resource_path = decoded.records[0].source.path
        assert resource_path.is_file()
        assert materialized.active is True
        rebuilt = build_session_document(decoded)
        assert [entry["name"] for entry in rebuilt.to_dict()["resources"].values()] == [
            entry["name"] for entry in document.to_dict()["resources"].values()
        ]

    assert materialized.active is False
    assert not resource_path.exists()
    with pytest.raises(SessionResourceError, match="no longer active"):
        session_to_request(materialized)


def test_session_document_save_load_and_render(tmp_path: Path) -> None:
    request = CircularDiagramRequest(
        records=(_record(),),
        output=RenderOutputRequest(
            output_prefix="canonical",
            output_directory=tmp_path,
            overwrite=True,
        ),
    )
    session_path = tmp_path / "canonical.gbdraw-session.json"
    save_session_document(session_path, request)
    assert session_path.read_text(encoding="utf-8").startswith(
        f'{{"format":"gbdraw-session","version":{CURRENT_SESSION_VERSION},'
    )
    document = load_session_document(session_path)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        result = render_session(materialized)

    assert result.mode == "circular"
    assert result.output_paths == (tmp_path / "canonical.svg",)
    assert result.output_paths[0].is_file()


def test_session_document_gzip_save_load(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw-session.json.gz"

    save_session_document(session_path, request)
    document = load_session_document(session_path)

    assert session_path.read_bytes().startswith(b"\x1f\x8b")
    assert document.version == CURRENT_SESSION_VERSION
    assert document.mode == "circular"


def test_session_document_requires_explicit_overwrite(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw-session.json"

    save_session_document(session_path, request)
    original = session_path.read_bytes()

    with pytest.raises(SessionFormatError, match="overwrite=True"):
        save_session_document(session_path, request)
    assert session_path.read_bytes() == original

    save_session_document(session_path, request, overwrite=True)
    assert load_session_document(session_path).mode == "circular"


def test_direct_session_render_ignores_stored_overwrite_permission(
    tmp_path: Path,
) -> None:
    request = CircularDiagramRequest(
        records=(_record(),),
        output=RenderOutputRequest(
            output_prefix="protected",
            output_directory=tmp_path,
            overwrite=True,
        ),
    )
    document = build_session_document(request)
    document._data["renderRequest"]["output"]["overwrite"] = True
    protected = tmp_path / "protected.svg"
    protected.write_text("keep this diagram", encoding="utf-8")

    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionRenderError, match="already exist"):
            render_session(materialized)

    assert protected.read_text(encoding="utf-8") == "keep this diagram"


def test_current_document_quarantines_legacy_protein_cache_on_save(
    tmp_path: Path,
) -> None:
    legacy_entry = {
        "schema": 2,
        "kind": "raw-losat",
        "key": "legacy-key",
        "text": "p_r_old\tp_r_other\n",
        "program": "blastp",
    }
    session_path = tmp_path / "legacy-round-trip.gbdraw-session.json"

    saved = save_session_document(
        session_path,
        LinearDiagramRequest(records=(_record(),)),
        adjunct={"losatCache": {"entries": [legacy_entry]}},
    )
    reloaded = load_session_document(session_path)

    assert saved.version == reloaded.version == CURRENT_SESSION_VERSION
    assert reloaded.to_dict()["losatCache"]["entries"] == []
    assert reloaded.to_dict()["legacyArtifacts"]["proteinRawCandidates"][
        "entries"
    ] == [
        {
            "state": "pending",
            "originalEntry": legacy_entry,
            "rejectionReason": None,
        }
    ]


def test_legacy_session_has_no_public_typed_conversion(tmp_path: Path) -> None:
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": 30,
            "files": {},
            "ui": {"mode": "circular"},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionVersionError, match="internal CLI replay only"):
            session_to_request(materialized)


def test_duplicate_json_resource_id_is_rejected(tmp_path: Path) -> None:
    session_path = tmp_path / "duplicate.gbdraw-session.json"
    session_path.write_text(
        '{"format":"gbdraw-session","version":31,'
        '"renderRequest":{},"resources":{"same":{},"same":{}}}',
        encoding="utf-8",
    )

    with pytest.raises(SessionFormatError, match="duplicate object key"):
        load_session_document(session_path)


def test_duplicate_sanitized_resource_name_is_rejected(tmp_path: Path) -> None:
    data = build_session_document(
        LinearDiagramRequest(records=(_record("a"), _record("b")))
    ).to_dict()
    resource_ids = list(data["resources"])
    data["resources"][resource_ids[1]]["name"] = data["resources"][resource_ids[0]]["name"]

    with pytest.raises(SessionResourceError, match="Duplicate canonical resource filename"):
        load_session_document(data)


def test_partial_materialization_failure_cleans_owned_directory(tmp_path: Path) -> None:
    data = build_session_document(
        LinearDiagramRequest(records=(_record("a"), _record("b")))
    ).to_dict()
    data["resources"]["record-2-genbank"]["data"] = "not-base64!"
    root = tmp_path / "materialized"

    with pytest.raises(SessionResourceError, match="record-2-genbank"):
        with materialize_session(
            data,
            output_directory=tmp_path,
            temporary_directory=root,
        ):
            pytest.fail("materialization unexpectedly succeeded")

    assert root.is_dir()
    assert list(root.iterdir()) == []


def test_session_document_returns_detached_payload() -> None:
    document = build_session_document(CircularDiagramRequest(records=(_record(),)))
    detached = document.to_dict()
    detached["renderRequest"]["mode"] = "linear"

    assert document.mode == "circular"
    json.dumps(document.to_dict())


def test_v32_web_slot_specs_drop_only_legacy_feature_geometry(
    tmp_path: Path,
) -> None:
    canonical_slots = (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3,legend_label=Genes",
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(
                linear_track_slots=canonical_slots,
                linear_track_axis_index=1,
            )
        ),
    )
    data = build_session_document(request).to_dict()
    data["version"] = 32
    encoded_slots = data["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ]
    encoded_slots[1] = (
        "features:features@side=below,h=48px,spacing=9px,"
        "z=3,legend_label=Genes"
    )
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.tracks.linear_track_slots == (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3,legend_label=Genes",
    )
    assert decoded.options.tracks.linear_track_axis_index == 1
    assert document.to_dict()["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ] == data["renderRequest"]["diagramOptions"]["tracks"]["linearTrackSlots"]


def test_v32_structured_slots_preserve_non_feature_geometry_and_fields(
    tmp_path: Path,
) -> None:
    feature_slot = LinearTrackSlot(
        id="features",
        renderer="features",
        enabled=False,
        side="below",
        z=3,
        params={"legend_label": "Genes"},
    )
    numeric_slot = LinearTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        side="above",
        height=ScalarSpec(27, "px"),
        spacing=ScalarSpec(4, "px"),
        params={"nt": "GC"},
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(
                linear_track_slots=(numeric_slot, feature_slot),
                linear_track_axis_index=1,
            )
        ),
    )
    data = build_session_document(request).to_dict()
    data["version"] = 32
    encoded_feature = data["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ][1]
    encoded_feature["height"] = {"value": 48, "unit": "px"}
    encoded_feature["spacing"] = {"value": 9, "unit": "px"}

    with materialize_session(data, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    preserved_numeric, migrated_feature = decoded.options.tracks.linear_track_slots
    assert migrated_feature == LinearTrackSlot(
        id="features",
        renderer="features",
        enabled=False,
        side="below",
        height=None,
        spacing=None,
        z=3,
        params={"legend_label": "Genes"},
    )
    assert preserved_numeric == numeric_slot
    assert decoded.options.tracks.linear_track_axis_index == 1


def test_current_session_preserves_v2_feature_geometry(tmp_path: Path) -> None:
    feature_slot = LinearTrackSlot(
        id="features",
        renderer="features",
        side="overlay",
        height=ScalarSpec(48, "px"),
        spacing=ScalarSpec(9, "px"),
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(linear_track_slots=(feature_slot,))
        ),
    )
    document = build_session_document(request)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.tracks.linear_track_slots == (feature_slot,)


def test_web_writer_payload_decodes_with_python_codec(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [node, "tests/web/session-request.test.mjs", "--print"],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": 31,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        records = normalize_request_records(request)

    assert request.output.output_prefix == "web-session"
    assert [record.id for record in records] == ["WEBTEST"]


def test_web_depth_writer_payload_decodes_with_python_codec(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [node, "tests/web/session-request.test.mjs", "--print-depth"],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)

    assert request.options.depth_tracks is not None
    assert [track.label for track in request.options.depth_tracks] == [
        "Sample A",
        "Sample B",
    ]
    assert request.options.depth_track_files is None
    assert request.options.depth_track_labels is None


def test_web_resolved_protein_writer_preserves_alignment_settings(
    tmp_path: Path,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--print-resolved-protein",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)

    assert isinstance(request, LinearDiagramRequest)
    assert request.options.protein_blastp_mode == "none"
    assert request.options.linear_comparisons is not None
    assert len(request.options.linear_comparisons) == 1
    assert request.options.align_orthogroup_feature == "resolved-feature-anchor"


@pytest.mark.parametrize("collinearity_value_kind", ("result", "blocks"))
def test_python_typed_protein_results_round_trip_through_web_projection(
    tmp_path: Path,
    collinearity_value_kind: str,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
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
    block = CollinearityBlock(
        block_id="block-1",
        query_record_index=0,
        subject_record_index=1,
        orientation="plus",
        score=100.0,
        anchors=(anchor,),
    )
    collinearity = (
        CollinearityResult(blocks=(block,), orthogroups=orthogroups)
        if collinearity_value_kind == "result"
        else (block,)
    )
    request = LinearDiagramRequest(
        records=(_record("record-a"), _record("record-b")),
        options=LinearDiagramOptions(
            orthogroups=orthogroups,
            collinearity_blocks=collinearity,
        ),
    )
    session_path = tmp_path / f"typed-{collinearity_value_kind}.json"
    save_session_document(session_path, request)

    completed = subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--round-trip-session",
            str(session_path),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    comparisons = canonical["renderRequest"]["comparisons"]
    orthogroup_entry = next(
        item for item in comparisons if item["kind"] == "orthogroupResult"
    )
    collinearity_entry = next(
        item for item in comparisons if item["kind"] == "collinearityResult"
    )
    assert orthogroup_entry["encoding"] == "canonicalJson"
    assert collinearity_entry["encoding"] == "canonicalJson"
    assert collinearity_entry["valueKind"] == collinearity_value_kind

    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
        }
    )
    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert isinstance(decoded, LinearDiagramRequest)
    assert decoded.options.orthogroups == orthogroups
    assert decoded.options.collinearity_blocks == collinearity


def test_python_depth_request_projects_in_web(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    shared_depth = pd.DataFrame(
        {
            "reference_name": ["first", "second"],
            "position": [1, 1],
            "depth": [5, 8],
        }
    )
    sparse_depth = pd.DataFrame(
        {
            "reference_name": ["first"],
            "position": [1],
            "depth": [13],
        }
    )
    request = LinearDiagramRequest(
        records=(_record("first"), _record("second")),
        options=LinearDiagramOptions(
            depth_tracks=(
                DepthTrackInput(
                    source=shared_depth,
                    label="Shared",
                    color="#112233",
                    height=18,
                    large_tick_interval=10,
                ),
                DepthTrackInput(
                    source=(sparse_depth, None),
                    label="Sparse",
                    color="#445566",
                    height=24,
                    small_tick_interval=5,
                    tick_font_size=9,
                ),
            ),
        ),
    )
    session_path = tmp_path / "python-canonical-depth.gbdraw-session.json"
    save_session_document(session_path, request)
    session = json.loads(session_path.read_text(encoding="utf-8"))

    assert len(session["renderRequest"]["diagramOptions"]["depthTracks"]) == 2
    subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--project-session",
            str(session_path),
        ],
        check=True,
        cwd=Path(__file__).parents[1],
    )
