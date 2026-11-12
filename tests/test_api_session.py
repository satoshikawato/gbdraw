from __future__ import annotations

import base64
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
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
    SessionDrawingSpec,
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
from gbdraw.api.options import losatp_analysis_mode
from gbdraw.exceptions import ValidationError
from gbdraw.session import SessionDocument
from gbdraw.session_container import BytesResourceSource
from gbdraw.session_io import CURRENT_SESSION_VERSION, load_session, read_source_bytes, split_resource_bytes
from tests.utils.session_files import embedded_session, load_session_mapping, read_session_payload, write_session_file
from tests.utils.session_drawing_view import build_mode_scoped_session_document, drawing_view, upgrade_to_mode_scoped
from tests.utils.web_render_result import web_result_svg, web_session_results


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
    # A new resource's bytes stay in their source until the Session file is written.
    assert "encoding" not in document.to_dict()["resources"]["record-1-genbank"]
    written = drawing_view(document.to_dict())
    assert written["renderRequest"]["output"]["overwrite"] is False
    assert written["results"] == []
    assert written["editorState"]["featureCatalog"] is None

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
    session_path = tmp_path / "canonical.gbdraw"
    save_session_document(session_path, request)
    with zipfile.ZipFile(session_path) as archive:
        names = archive.namelist()
        assert json.loads(archive.read("gbdraw-session.json"))["version"] == CURRENT_SESSION_VERSION
    assert names == ["gbdraw-session.json", "resources/record-1-genbank"]
    document = load_session_document(session_path)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        result = render_session(materialized)

    assert result.mode == "circular"
    assert result.output_paths == (tmp_path / "canonical.svg",)
    assert result.output_paths[0].is_file()


@pytest.mark.parametrize("name", ["canonical.json", "canonical.json.gz", "CANONICAL.JSON"])
def test_session_document_refuses_a_json_output_name(tmp_path: Path, name: str) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / name

    with pytest.raises(SessionFormatError) as refused:
        save_session_document(session_path, request)

    assert str(refused.value) == (
        "A Session is saved as a .gbdraw file. Use a name that ends with .gbdraw, "
        f"not .json or .json.gz: {session_path}"
    )
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name", ["canonical.zip", "canonical"])
def test_session_document_keeps_any_other_output_name(tmp_path: Path, name: str) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / name

    save_session_document(session_path, request)

    assert list(tmp_path.iterdir()) == [session_path]
    assert load_session_document(session_path).mode == "circular"


def test_session_file_stores_a_gzip_input_as_it_is(tmp_path: Path) -> None:
    import gzip
    import io

    from Bio import SeqIO

    from gbdraw.api import GenBankInputSource

    text = io.StringIO()
    SeqIO.write(_record().source.record, text, "genbank")
    genbank = tmp_path / "input.gb.gz"
    genbank.write_bytes(gzip.compress(text.getvalue().encode(), mtime=0))
    request = CircularDiagramRequest(records=(RecordInput(source=GenBankInputSource(genbank)),))
    session_path = tmp_path / "gzip.gbdraw"

    save_session_document(session_path, request)

    with zipfile.ZipFile(session_path) as archive:
        (entry,) = [info for info in archive.infolist() if info.filename.startswith("resources/")]
        assert entry.compress_type == zipfile.ZIP_STORED
        assert archive.read(entry) == genbank.read_bytes()
    with materialize_session(load_session_document(session_path), output_directory=tmp_path) as materialized:
        (path,) = materialized.resource_paths.values()
        assert path.read_bytes() == genbank.read_bytes()


def test_a_json_session_51_is_refused(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw"
    save_session_document(session_path, request)
    json_path = tmp_path / "canonical.gbdraw"
    with zipfile.ZipFile(session_path) as archive:
        json_path.write_bytes(archive.read("gbdraw-session.json"))

    with pytest.raises(SessionFormatError) as refused:
        load_session_document(json_path)

    assert str(refused.value) == (
        f"A Session version {CURRENT_SESSION_VERSION} file must be a .gbdraw file, not JSON: {json_path}"
    )


def test_session_document_requires_explicit_overwrite(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw"

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
    document._data["drawings"][0]["renderRequest"]["output"]["overwrite"] = True
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
    session_path = tmp_path / "legacy-round-trip.gbdraw"

    saved = save_session_document(
        session_path,
        drawings=[
            SessionDrawingSpec(
                LinearDiagramRequest(records=(_record(),)),
                state={"losatCache": {"entries": [legacy_entry]}},
            )
        ],
    )
    reloaded = load_session_document(session_path)

    assert saved.version == reloaded.version == CURRENT_SESSION_VERSION
    assert drawing_view(reloaded.to_dict())["losatCache"]["entries"] == []
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
    session_path = tmp_path / "duplicate.gbdraw"
    session_path.write_text(
        '{"format":"gbdraw-session","version":31,'
        '"renderRequest":{},"resources":{"same":{},"same":{}}}',
        encoding="utf-8",
    )

    with pytest.raises(SessionFormatError, match="duplicate object key"):
        load_session_document(session_path)


def test_duplicate_sanitized_resource_name_is_rejected(tmp_path: Path) -> None:
    built = build_session_document(LinearDiagramRequest(records=(_record("a"), _record("b"))))
    data = embedded_session(built.to_dict(), built.resources.sources)
    resource_ids = list(data["resources"])
    data["resources"][resource_ids[1]]["name"] = data["resources"][resource_ids[0]]["name"]

    with pytest.raises(SessionResourceError, match="Duplicate canonical resource filename"):
        load_session_mapping(data)


def _session_with_resource_field(tmp_path: Path, field: str, value) -> tuple[dict, Path]:
    """A saved JSON Session 48 whose GenBank resource also carries ``field``.

    ``value`` maps the SHA-256 hex digest of the resource bytes to the value.
    A Session file's resources carry no checksum (its entries' CRC-32 replaces it).
    """

    data = build_mode_scoped_session_document(CircularDiagramRequest(records=(_record(),)))
    resource = data["resources"]["record-1-genbank"]
    resource[field] = value(hashlib.sha256(base64.b64decode(resource["data"])).hexdigest())
    path = tmp_path / "resource-field.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return data, path


@pytest.mark.parametrize(
    "checksum",
    (
        lambda digest: f"sha256:{digest}",
        lambda digest: digest,
        lambda digest: f"SHA256:{digest.upper()}",
    ),
    ids=("prefixed", "bare-hex", "upper-case"),
)
def test_resource_checksum_loads_through_both_session_loaders(
    tmp_path: Path,
    checksum,
) -> None:
    data, path = _session_with_resource_field(tmp_path, "checksum", checksum)

    document = load_session_document(path)
    # The checksum is checked where the bytes are read, then dropped with them.
    assert "checksum" not in document.to_dict()["resources"]["record-1-genbank"]
    assert load_session_document(data).to_dict() == document.to_dict()
    assert load_session(path)["resources"] == document.to_dict()["resources"]
    with materialize_session(document, output_directory=tmp_path) as materialized:
        records = normalize_request_records(session_to_request(materialized))
    assert [record.id for record in records] == ["record"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("checksum", lambda digest: "sha256:" + "0" * 64, "checksum does not match"),
        ("checksum", lambda digest: "md5:" + digest[:32], "checksum does not match"),
        ("checksum", lambda digest: 7, "checksum does not match"),
        ("sha256", lambda digest: digest, r"unknown field\(s\): sha256"),
    ),
    ids=("wrong-digest", "other-algorithm", "not-a-string", "unknown-field"),
)
def test_both_session_loaders_reject_the_same_resource_descriptor(
    tmp_path: Path,
    field: str,
    value,
    message: str,
) -> None:
    data, path = _session_with_resource_field(tmp_path, field, value)

    with pytest.raises(SessionResourceError, match=message):
        load_session_document(data)
    with pytest.raises(SessionResourceError, match=message):
        load_session_document(path)
    with pytest.raises(ValidationError, match=message):
        load_session(path)


@pytest.mark.browser
def test_web_saved_resource_checksum_loads_through_both_session_loaders(
    tmp_path: Path,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    data, path = _session_with_resource_field(
        tmp_path, "checksum", lambda digest: f"sha256:{digest}"
    )
    completed = subprocess.run(
        [node, "tests/web/session-request.test.mjs", "--round-trip-session", str(path)],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    # The Web writer keeps the checksum of a resource it adopted from a Session.
    assert (
        canonical["resources"]["record-1-genbank"]["checksum"]
        == data["resources"]["record-1-genbank"]["checksum"]
    )
    # A JSON Session 48 (the mode-scoped layout ``data`` already has) holding what the Web wrote.
    saved_path = tmp_path / "web-saved.json"
    saved_path.write_text(
        json.dumps(
            {
                **data,
                "version": 48,
                "renderRequest": canonical["renderRequest"],
                "resources": canonical["resources"],
            }
        ),
        encoding="utf-8",
    )

    document = load_session_document(saved_path)
    assert load_session(saved_path)["resources"] == document.to_dict()["resources"]
    with materialize_session(document, output_directory=tmp_path) as materialized:
        records = normalize_request_records(session_to_request(materialized))
    assert [record.id for record in records] == ["record"]


def test_partial_materialization_failure_cleans_owned_directory(tmp_path: Path) -> None:
    built = build_session_document(LinearDiagramRequest(records=(_record("a"), _record("b"))))
    # The second resource's bytes are shorter than its descriptor says.
    document = SessionDocument._built(
        built.to_dict(), {**built.resources.sources, "record-2-genbank": BytesResourceSource(b"short")}
    )
    root = tmp_path / "materialized"

    with pytest.raises(SessionResourceError, match="record-2-genbank"):
        with materialize_session(
            document,
            output_directory=tmp_path,
            temporary_directory=root,
        ):
            pytest.fail("materialization unexpectedly succeeded")

    assert root.is_dir()
    assert list(root.iterdir()) == []


def test_session_document_returns_detached_payload() -> None:
    document = build_session_document(CircularDiagramRequest(records=(_record(),)))
    detached = document.to_dict()
    detached["drawings"][0]["renderRequest"]["mode"] = "linear"

    assert document.mode == "circular"
    json.dumps(document.to_dict())


def test_v32_web_slot_specs_drop_only_legacy_feature_geometry(
    tmp_path: Path,
) -> None:
    canonical_slots = (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3",
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
    data = build_mode_scoped_session_document(request)
    data["version"] = 32
    encoded_slots = data["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ]
    encoded_slots[1] = (
        "features:features@side=below,h=48px,spacing=9px,"
        "z=3"
    )
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.tracks.linear_track_slots == (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3",
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
        params={},
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
    data = build_mode_scoped_session_document(request)
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
        params={},
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


@pytest.mark.browser
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


@pytest.mark.browser
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
            "version": 48,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
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


@pytest.mark.browser
def test_web_resolved_protein_writer_preserves_typed_alignment_layout(
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
            "version": 48,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None, "alignmentResetReceipt": None},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)

    assert isinstance(request, LinearDiagramRequest)
    assert losatp_analysis_mode(request.options.losat_search) == "none"
    assert request.options.linear_comparisons is not None
    assert len(request.options.linear_comparisons) == 1
    assert not hasattr(request.options, "align_orthogroup_feature")
    assert request.similarity_alignment is not None
    assert request.similarity_alignment.group_id == "og-resolved"
    assert request.similarity_alignment.schema == 2
    assert request.similarity_alignment.records[2].status.value == "aligned"
    assert request.records[2].presentation.reverse_complement is False
    assert request.layout is not None
    assert [item.record_key for item in request.layout.record_translations] == [
        "first",
        "second",
        "third",
    ]


@pytest.mark.parametrize("collinearity_value_kind", ("result", "blocks"))
@pytest.mark.browser
def test_python_typed_protein_round_trip_with_explicit_adjacent_comparison(
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
    session_path = tmp_path / f"typed-{collinearity_value_kind}.gbdraw"
    save_session_document(session_path, request)

    completed = subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--round-trip-session",
            str(session_path),
            "--activate-adjacent-comparison",
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
            "version": 48,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
        }
    )
    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert isinstance(decoded, LinearDiagramRequest)
    from dataclasses import replace
    from gbdraw.analysis.protein_colinearity import materialize_ortholog_paths
    assert materialize_ortholog_paths(decoded.options.orthogroups) == orthogroups
    decoded_collinearity = decoded.options.collinearity_blocks
    if isinstance(decoded_collinearity, CollinearityResult):
        decoded_collinearity = replace(decoded_collinearity, orthogroups=materialize_ortholog_paths(decoded_collinearity.orthogroups))
    assert decoded_collinearity == collinearity


def test_two_depth_tracks_naming_one_resource_render_as_two_copies(tmp_path: Path) -> None:
    """The Web Save Session writes one Depth file used twice once (D-220, D-228); Python reads that Session."""
    depth = pd.DataFrame({"reference_name": ["record"] * 2, "position": [1, 5], "depth": [3, 9]})
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(depth_tracks=(
            DepthTrackInput(source=depth, label="A"), DepthTrackInput(source=depth, label="B"),
        )),
        output=RenderOutputRequest(output_prefix="depth", overwrite=True),
    )
    built = build_session_document(request)
    copies = embedded_session(built.to_dict(), built.resources.sources)
    (drawing,) = copies["drawings"]
    tracks = drawing["renderRequest"]["diagramOptions"]["depthTracks"]
    first, second = (track["source"]["resourceId"] for track in tracks)
    assert first != second
    assert copies["resources"][first]["data"] == copies["resources"][second]["data"]
    once = json.loads(json.dumps(copies))
    once["drawings"][0]["renderRequest"]["diagramOptions"]["depthTracks"][1]["source"]["resourceId"] = first
    del once["resources"][second]
    svgs = []
    for name, document in (("copies", copies), ("once", once)):
        session_path = tmp_path / f"{name}.gbdraw"
        write_session_file(session_path, document)
        output = tmp_path / name
        output.mkdir()
        with materialize_session(load_session_document(session_path), output_directory=output) as materialized:
            result = render_session(materialized)
        svgs.append(result.output_paths[0].read_text(encoding="utf-8"))
    assert svgs[0] == svgs[1]


@pytest.mark.browser
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
    session_path = tmp_path / "python-canonical-depth.gbdraw"
    save_session_document(session_path, request)
    session = drawing_view(read_session_payload(session_path))

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


def _record_local_collinear_session(tmp_path: Path, search_scope: str = "adjacent") -> Path:
    """Render the Issue #460 evidence through typed, Web, and session boundaries."""
    import xml.etree.ElementTree as ET
    from dataclasses import replace

    from tests.test_collinearity import _record_local_mixed_fixture
    from gbdraw.analysis.collinearity import build_orthogroup_collinearity_blocks_from_hits
    from gbdraw.api import render_request
    from gbdraw.web_support.feature_catalog import build_feature_catalog_item
    from gbdraw.web_support.request_render import render_canonical_web_request

    records, extraction, hits = _record_local_mixed_fixture()
    for record in records:
        record.name = record.id
        record.annotations["organism"] = "synthetic construct"
    groups = build_orthogroup_collinearity_blocks_from_hits(
        hits, extraction, records=records, search_scope=search_scope,
        infer_orthogroups=True,
    )
    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=InMemoryRecordSource(record), record_key=record.id) for record in records),
        options=LinearDiagramOptions(
            orthogroups=groups.orthogroups,
            collinearity_blocks=groups, collinearity_search_scope=search_scope,
        ),
        output=RenderOutputRequest(
            output_directory=tmp_path / "typed", output_prefix="mixed",
            formats=("svg", "interactive_svg"),
        ),
    )
    typed = render_request(request, include_feature_catalog=True)
    svg = (tmp_path / "typed" / "mixed.svg").read_text()
    catalog_item = build_feature_catalog_item(
        svg, typed.interactive_context, result_index=0, result_name="mixed.svg",
    )
    document = build_session_document(request)
    with materialize_session(document, output_directory=tmp_path / "decoded") as materialized:
        decoded = session_to_request(materialized)
        assert decoded.options.collinearity_blocks == groups
        assert decoded.options.orthogroups == groups.orthogroups
        web = render_canonical_web_request(
            drawing_view(document.to_dict())["renderRequest"], resource_paths=materialized.resource_paths,
            output_directory=tmp_path / "web",
        )
    assert web["featureCatalog"]["items"] == [catalog_item]
    semantic = {group["id"]: group for group in catalog_item["orthogroups"]}
    assert set(semantic) == {"og_1", "og_2"}
    assert semantic["og_1"]["scope"] == "cross_record"
    assert semantic["og_2"]["scope"] == "record_local"
    assert semantic["og_2"]["source_record_index"] == 0
    biological = {
        (feature["recordKey"], feature["biologicalFeatureId"]): feature
        for feature in catalog_item["biologicalFeatures"]
    }
    for group_id, expected in (
        ("og_1", {"a0": "anchor", "a1": "inparalog", "b0": "anchor"}),
        ("og_2", {"a3": "local_paralog", "a4": "local_paralog"}),
    ):
        members = semantic[group_id]["members"]
        assert {
            biological[(member["recordKey"], member["biologicalFeatureId"])]["qualifiers"]["protein_id"][0]: member["role"]
            for member in members
        } == expected
        if group_id == "og_2":
            assert {biological[(m["recordKey"], m["biologicalFeatureId"])]["record_id"] for m in members} == {"record_a"}

    def comparison_paths(source):
        return [element for element in ET.fromstring(source).iter()
                if element.tag.endswith("}path") and element.get("data-match-kind") == "collinear"]

    paths = comparison_paths(svg)
    assert len(groups.blocks) == len(paths) == 1
    assert paths[0].get("data-orthogroup-id") == "og_1"
    assert all(block.query_record_index != block.subject_record_index for block in groups.blocks)
    geometry_only = build_orthogroup_collinearity_blocks_from_hits(
        {pair: table for pair, table in hits.items() if pair[0] != pair[1]},
        extraction, records=records, search_scope=search_scope,
        infer_orthogroups=True,
    )
    baseline = render_request(replace(
        request,
        options=replace(request.options, orthogroups=geometry_only.orthogroups, collinearity_blocks=geometry_only),
        output=replace(request.output, output_directory=tmp_path / "geometry", formats=("svg",)),
    ))
    assert [path.get("d") for path in paths] == [
        path.get("d") for path in comparison_paths(baseline.output_paths[0].read_text())
    ]
    interactive_source = next(path for path in typed.output_paths if "interactive" in path.name).read_text()
    interactive_root = ET.fromstring(interactive_source)
    metadata = next(element for element in interactive_root.iter()
                    if element.get("id") == "gbdraw-interactive-feature-metadata")
    interactive = json.loads(metadata.text)
    assert interactive["items"][0]["orthogroups"] == catalog_item["orthogroups"]
    assert len(comparison_paths(interactive_source)) == 1
    assert [path.get("d") for path in comparison_paths(web_result_svg(web["results"][0]))] == [paths[0].get("d")]

    session_path = tmp_path / "mixed.gbdraw"
    saved = save_session_document(session_path, drawings=[SessionDrawingSpec(request, state={
        "results": web_session_results(web),
        "editorState": {"featureCatalog": web["featureCatalog"]},
    })])
    reloaded = load_session_document(session_path)
    assert reloaded.version == saved.version == CURRENT_SESSION_VERSION
    assert drawing_view(reloaded.to_dict())["editorState"]["featureCatalog"] == web["featureCatalog"]
    with materialize_session(reloaded, output_directory=tmp_path / "reload") as materialized:
        restored = session_to_request(materialized)
        assert restored.options.collinearity_blocks == groups
        rerendered = render_canonical_web_request(
            drawing_view(reloaded.to_dict())["renderRequest"], resource_paths=materialized.resource_paths,
            output_directory=tmp_path / "reload-web",
        )
    assert rerendered["featureCatalog"] == web["featureCatalog"]
    assert [path.get("d") for path in comparison_paths(web_result_svg(rerendered["results"][0]))] == [paths[0].get("d")]
    return session_path


@pytest.mark.parametrize("search_scope", ["adjacent", "all"])
def test_record_local_collinear_catalog_session_round_trip(tmp_path: Path, search_scope: str) -> None:
    _record_local_collinear_session(tmp_path, search_scope)


def test_composite_document_replays_only_committed_request(tmp_path):
    import base64
    payload = json.loads((Path(__file__).parent / 'fixtures/sessions/single.v41-bindings1.json').read_text())
    committed = load_session_document(payload)
    # Deliberately non-GenBank draft bytes: replay must never consume them.
    payload['resources']['draft'] = {
        'kind': 'web-file', 'name': 'draft.gb', 'type': 'text/plain',
        'encoding': 'base64', 'size': 5, 'data': base64.b64encode(b'draft').decode(),
    }
    leaf = {'resourceId': 'draft', 'name': 'draft.gb', 'type': '', 'lastModified': 0}
    payload['webFiles']['bindings'].update(schema=2, c_gb={
        'kind': 'composite', 'components': [leaf, {**leaf, 'name': 'again.gb'}],
        'name': 'logical.gb', 'type': 'text/plain', 'lastModified': 0,
    })
    document = load_session_document(payload)
    saved = tmp_path / 'composite.json'
    write_session_file(saved, document.to_dict(), document.resources.sources)
    restored = load_session_document(saved)
    assert restored.to_dict()['webFiles'] == payload['webFiles']
    with materialize_session(restored, output_directory=tmp_path / 'composite') as materialized:
        request = session_to_request(materialized)
        assert len(request.records) == 1
        result = render_session(materialized)
    with materialize_session(committed, output_directory=tmp_path / 'committed') as materialized:
        control = render_session(materialized)
    assert result.output_paths[0].read_bytes() == control.output_paths[0].read_bytes()


# OV-1145: Session 48 stores the protein results of a reverse-complemented
# Linear record in its search frame; an upgrade or re-save of a Session 46/47
# converts them.
_SESSION_FIXTURES = Path(__file__).parent / "fixtures" / "sessions"
_V47_BGC = _SESSION_FIXTURES / "BGC0000708-BGC0000713.v47.gbdraw-session.json.gz"
_PROTEIN_RESOURCES = ("comparison-canonical-orthogroups-1", "comparison-canonical-protein-4")


def _v47_bgc() -> dict:
    import gzip

    return json.loads(gzip.decompress(_V47_BGC.read_bytes()))


def _two_drawings(version: int) -> dict:
    from tests.utils.two_mode_session import two_mode_session

    from tests.test_session_contract import _with_current_request

    # The fixture stores the schema 9 request of its time; the Circular Gallery Session is at the current schema.
    return {**two_mode_session(_with_current_request(_v47_bgc())), "version": version}


@pytest.mark.parametrize("version", [46, 47])
def test_an_upgrade_stores_the_protein_results_of_a_reversed_record_in_the_search_frame(
    tmp_path: Path, version: int,
) -> None:
    from gbdraw.session import upgrade_session_document

    source = _two_drawings(version)
    upgrade = upgrade_session_document(source, temporary_directory=tmp_path)
    upgraded = upgrade.document.to_dict()
    source_document, source_sources = split_resource_bytes(source)

    assert upgrade.warnings == ()
    assert upgraded["version"] == CURRENT_SESSION_VERSION == 51
    assert [drawing.id for drawing in upgrade.document.drawings] == ["circular", "linear"]
    # Bytes change under the stored IDs, only for the protein results.
    assert upgraded["resources"].keys() == source["resources"].keys()
    changed = sorted(
        key for key, descriptor in upgraded["resources"].items()
        if descriptor != source_document["resources"][key]
        or upgrade.document.resources.read(key) != read_source_bytes(source_sources[key], descriptor["size"])
    )
    assert changed == sorted(f"linear-{key}" for key in _PROTEIN_RESOURCES)
    for key in changed:
        stored = upgrade.document.resources.read(key)
        assert upgraded["resources"][key]["size"] == len(stored)
    again = upgrade_session_document(upgrade.document, temporary_directory=tmp_path).document
    assert again.to_dict() == upgraded


def test_a_cli_resave_of_a_session_47_writes_the_search_frame_for_both_drawings(
    tmp_path: Path,
) -> None:
    from gbdraw.cli_utils.session import render_main
    from gbdraw.session import load_session_document, upgrade_session_document
    from tests.utils.session_files import write_session_file

    source = tmp_path / "two.gbdraw"
    write_session_file(source, _two_drawings(47))
    render_main(["--session", str(source), "-o", str(tmp_path / "fig"), "-f", "svg", "--save_session"])
    resaved = load_session_document(tmp_path / "fig.gbdraw")
    upgraded = upgrade_session_document(source, temporary_directory=tmp_path).document

    assert resaved.version == CURRENT_SESSION_VERSION
    assert [drawing.id for drawing in resaved.drawings] == ["circular", "linear"]

    def protein_results(document) -> tuple[object, ...]:
        # Both hold the same search-frame results.
        from gbdraw.session import materialize_session, session_to_request

        with materialize_session(document, output_directory=tmp_path) as materialized:
            options = session_to_request(materialized, drawing="linear").options
        return (
            [comparison.matches.to_dict("records") for comparison in options.linear_comparisons or ()],
            repr(options.orthogroups),
        )

    assert protein_results(resaved) == protein_results(upgraded)
    # The re-save keeps the IDs of the converted resources, as the upgrade does
    # (the base re-save wrote the table under a new ID).
    for key in (f"linear-{key}" for key in _PROTEIN_RESOURCES):
        assert resaved.resources.read(key) == upgraded.resources.read(key)


def test_python_load_reads_a_session_47_and_a_settings_only_session_47() -> None:
    from gbdraw.session import load_session_document
    from gbdraw.session_io import is_settings_only_session

    data = _v47_bgc()
    data["webFiles"]["bindings"] = {"schema": 2, "c_gb": None}
    assert load_session_document(data).version == 47
    # A settings-only Session 47: the current one, written by a 47 writer.
    current = upgrade_to_mode_scoped(_SESSION_FIXTURES / "settings-only.v42.json.gz")
    settings_only = {**current, "version": 47}
    assert is_settings_only_session(settings_only)
    assert load_session_document(settings_only).drawings == load_session_document(current).drawings


def test_session_save_warns_about_a_linear_record_once_and_refuses_an_existing_output_first(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1246: a save to an existing path built the whole document, and warned
    about the record's topology again, before it refused the output."""

    import logging

    import gbdraw.session as session_module

    seqrecord = SeqRecord(
        Seq("ATGCGCAT"),
        id="linear_record",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    request = CircularDiagramRequest(records=(RecordInput(source=InMemoryRecordSource(seqrecord)),))
    session_path = tmp_path / "linear.gbdraw"
    question = "Are you sure you want to visualize it as circular?"

    with caplog.at_level(logging.WARNING):
        save_session_document(session_path, request)
    assert caplog.text.count(question) == 1

    caplog.clear()
    builds: list[object] = []
    real_build = session_module.build_session_document
    monkeypatch.setattr(
        session_module,
        "build_session_document",
        lambda *args, **kwargs: builds.append(args) or real_build(*args, **kwargs),
    )
    with caplog.at_level(logging.WARNING), pytest.raises(SessionFormatError, match="already exists"):
        save_session_document(session_path, request)
    assert builds == []
    assert question not in caplog.text


def test_an_edited_document_of_a_session_file_loads_with_the_files_of_its_source(tmp_path: Path) -> None:
    """``to_dict()`` holds no bytes; the edited mapping loads with the
    resources of the document it came from, and a mapping without them is
    refused as it was before."""

    source = load_session_document(
        Path(__file__).resolve().parents[1] / "gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw"
    )
    edited = source.to_dict()
    edited["title"] = "edited"

    document = load_session_document(edited, resources=source.resources)

    assert document.to_dict()["title"] == "edited"
    assert {resource_id: document.resources.read(resource_id) for resource_id in document.to_dict()["resources"]} == {
        resource_id: source.resources.read(resource_id) for resource_id in source.to_dict()["resources"]
    }
    with pytest.raises(ValidationError, match="has no bytes"):
        load_session_document(edited)


RETIRED_JSON_51 = (
    Path(__file__).resolve().parent / "retired_sessions/refused/drawings-per-drawing-cards.v51.gbdraw-session.json.gz"
)


def test_a_json_session_51_that_dev_wrote_is_refused_by_every_reader(tmp_path: Path) -> None:
    """A JSON Session 51 never reached main; each reader refuses the one dev wrote (PD-OI-127)."""

    message = f"A Session version {CURRENT_SESSION_VERSION} file must be a .gbdraw file, not JSON: {RETIRED_JSON_51}"
    with pytest.raises(SessionFormatError) as refused:
        load_session_document(RETIRED_JSON_51)
    assert str(refused.value) == message
    with pytest.raises(ValidationError) as refused_payload:
        load_session(RETIRED_JSON_51)
    assert str(refused_payload.value) == message
    completed = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "render", "--session", str(RETIRED_JSON_51), "-o", str(tmp_path / "out")],
        capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert completed.returncode == 1
    assert message in completed.stderr
    assert not any(tmp_path.iterdir())


GZIP_AND_PLAIN = Path(__file__).resolve().parent / "fixtures/sessions/cli-circular-gzip-and-plain.v51.gbdraw"


def test_the_session_file_stores_a_gzip_input_as_it_is_and_deflates_a_plain_one(tmp_path: Path) -> None:
    """The positive fixture of the Session file container: its gzip input
    entry is stored, its plain table entry is deflated, and it replays the
    figure it holds."""

    with zipfile.ZipFile(GZIP_AND_PLAIN) as archive:
        entries = archive.infolist()
    assert entries[0].filename == "gbdraw-session.json"
    methods = {entry.filename: entry.compress_type for entry in entries[1:]}
    assert methods["resources/record-1-genbank"] == zipfile.ZIP_STORED
    assert methods["resources/colors-color-table"] == zipfile.ZIP_DEFLATED

    document = load_session_document(GZIP_AND_PLAIN)
    assert document.version == CURRENT_SESSION_VERSION
    assert document.resources.read("record-1-genbank")[:2] == b"\x1f\x8b"
    assert b"cytochrome b\t#e41a1c" in document.resources.read("colors-color-table")
    with materialize_session(document, output_directory=tmp_path) as materialized:
        result = render_session(materialized)
    (stored,) = drawing_view(document.to_dict())["results"]
    assert result.output_paths[0].read_text(encoding="utf-8") == stored["content"]


GALLERY_BASIC = Path(__file__).resolve().parents[1] / "gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw"


def _with_damaged_deflate(source: Path, entry_name: str, target: Path) -> Path:
    """``source`` with the deflated data of ``entry_name`` made unreadable (an invalid block type)."""

    data = bytearray(source.read_bytes())
    info = zipfile.ZipFile(source).getinfo(entry_name)
    assert info.compress_type == zipfile.ZIP_DEFLATED
    name_length, extra_length = int.from_bytes(data[info.header_offset + 26: info.header_offset + 28], "little"), int.from_bytes(
        data[info.header_offset + 28: info.header_offset + 30], "little")
    start = info.header_offset + 30 + name_length + extra_length
    data[start: start + 2] = b"\xff\xff"
    target.write_bytes(bytes(data))
    return target


def test_a_session_file_whose_document_data_is_damaged_is_refused_as_damaged(tmp_path: Path) -> None:
    damaged = _with_damaged_deflate(GALLERY_BASIC, "gbdraw-session.json", tmp_path / "damaged.gbdraw")

    with pytest.raises(SessionFormatError) as refused:
        load_session_document(damaged)

    assert str(refused.value) == "This file is not a gbdraw Session: it is incomplete or damaged."


def test_a_session_file_whose_input_data_is_damaged_is_refused_as_damaged(tmp_path: Path) -> None:
    damaged = _with_damaged_deflate(GALLERY_BASIC, "resources/record-1-genbank", tmp_path / "damaged.gbdraw")

    document = load_session_document(damaged)
    with pytest.raises(SessionResourceError, match="This file is not a gbdraw Session: it is incomplete or damaged."):
        with materialize_session(document, output_directory=tmp_path):
            pass
    completed = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "render", "--session", str(damaged), "-o", str(tmp_path / "out")],
        capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert completed.returncode == 1
    assert "This file is not a gbdraw Session: it is incomplete or damaged." in completed.stderr
    assert "Traceback" not in completed.stderr


def test_a_session_file_holding_an_older_document_version_is_refused(tmp_path: Path) -> None:
    """Only a Session 51 was ever written as a .gbdraw file (PD-OI-127)."""

    path = tmp_path / "older.gbdraw"
    with zipfile.ZipFile(GALLERY_BASIC) as source, zipfile.ZipFile(path, "w") as target:
        for info in source.infolist():
            data = source.read(info)
            if info.filename == "gbdraw-session.json":
                data = json.dumps({**json.loads(data), "version": 50}).encode("utf-8")
            target.writestr(info, data)

    with pytest.raises(SessionFormatError) as refused:
        load_session_document(path)

    assert str(refused.value) == "Unsupported session version: 50."


def test_a_session_51_mapping_with_embedded_inputs_is_refused() -> None:
    """A JSON Session 51 never reached main, as a file or as a mapping (PD-OI-127)."""

    source = load_session_document(GALLERY_BASIC)
    mapping = embedded_session(source.to_dict(), source.resources.sources)

    for read in (load_session_document, SessionDocument):
        with pytest.raises(ValidationError) as refused:
            read(mapping)
        assert str(refused.value) == (
            f"A Session version {CURRENT_SESSION_VERSION} file must be a .gbdraw file, not JSON: <mapping>"
        )


def test_a_mapping_without_its_resources_names_how_to_pass_them() -> None:
    mapping = load_session_document(GALLERY_BASIC).to_dict()
    message = (
        "Canonical resource 'record-1-genbank' has no bytes: pass resources=document.resources with a "
        "document's to_dict(), or save and load the .gbdraw file."
    )

    for read in (load_session_document, SessionDocument):
        with pytest.raises(ValidationError) as refused:
            read(mapping)
        assert str(refused.value) == message
