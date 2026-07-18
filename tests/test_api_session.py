from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    SessionFormatError,
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

    assert document.version == 32
    assert document.mode == mode
    assert document.has_canonical_request is True
    assert document.to_dict()["resources"]["record-1-genbank"]["encoding"] == "base64"

    with materialize_session(
        document,
        output_directory=tmp_path,
        temporary_directory=tmp_path / "materialized",
    ) as materialized:
        decoded = session_to_request(materialized)
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
        '{"format":"gbdraw-session","version":32,'
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
    assert document.version == 32
    assert document.mode == "circular"


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
