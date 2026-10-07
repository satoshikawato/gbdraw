"""The Session resource table: one allocator for every Python Session writer."""

from __future__ import annotations

import base64
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pytest

from gbdraw.api.options import CircularDiagramOptions
from gbdraw.api.request_render import resolve_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.session import (
    SessionResourceError,
    _build_session_document_from_resolved_request,
    build_session_document,
    load_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.session_resources import SessionResourceTable, canonical_resource_ids

FIXTURES = Path(__file__).parent / "fixtures" / "sessions"
GENBANK = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
RING_ROW = "q\ts\t90\t100\t0\t0\t1\t100\t1\t100\t1e-10\t200\n"


def _descriptor(name: str, content: bytes, *, kind: str = "web-file") -> dict[str, Any]:
    return {
        "kind": kind,
        "name": name,
        "type": "text/plain",
        "size": len(content),
        "lastModified": 0,
        "encoding": "base64",
        "data": base64.b64encode(content).decode("ascii"),
    }


def _request(*ring_files: Path) -> CircularDiagramRequest:
    return CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(GENBANK)),),
        options=CircularDiagramOptions(
            conservation_blast_files=tuple(map(str, ring_files)) or None,
        ),
        output=RenderOutputRequest(output_prefix="stored"),
    )


def _ring(directory: Path, name: str = "ring.tsv", row: str = RING_ROW) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(row, encoding="utf-8")
    return path


def _rings(payload: dict[str, Any]) -> list[str]:
    return [ref["resourceId"] for ref in payload["diagramOptions"]["conservationBlastFiles"]]


def test_equal_bytes_and_name_take_the_held_resource() -> None:
    held = _descriptor("HmmtDNA.gbk", GENBANK.read_bytes(), kind="genbank")
    table = SessionResourceTable({"genbank-1": held})

    encoded = encode_canonical_request(_request(), table=table)

    assert encoded.payload["records"][0]["source"] == {
        "kind": "genbank",
        "resourceId": "genbank-1",
    }
    assert encoded.resources == ()
    assert table.descriptors() == {"genbank-1": held}


def test_a_held_id_with_other_bytes_takes_the_next_number() -> None:
    table = SessionResourceTable({
        "record-1-genbank": _descriptor("other.gbk", b"other bytes", kind="genbank"),
    })

    encoded = encode_canonical_request(_request(), table=table)

    assert encoded.payload["records"][0]["source"]["resourceId"] == "record-1-genbank-2"
    assert [(item.resource_id, item.name) for item in encoded.resources] == [
        ("record-1-genbank-2", "HmmtDNA.gbk"),
    ]
    assert list(table.descriptors()) == ["record-1-genbank", "record-1-genbank-2"]


def test_a_held_name_with_other_bytes_numbers_the_name_and_pins_ring_labels(
    tmp_path: Path,
) -> None:
    ring = _ring(tmp_path)
    alone = encode_canonical_request(_request(ring))
    assert "conservationLabels" not in alone.payload["diagramOptions"]
    table = SessionResourceTable({"ring-1": _descriptor("ring.tsv", b"other rows\n")})

    encoded = encode_canonical_request(_request(ring), table=table)

    assert _rings(encoded.payload) == ["conservation-blast-files-1"]
    assert table.descriptors()["conservation-blast-files-1"]["name"] == "ring.2.tsv"
    # A replay labels the ring with its file name, so the label is stored.
    assert encoded.payload["diagramOptions"]["conservationLabels"] == ["ring.tsv"]


def test_equal_bytes_under_another_name_pin_ring_labels(tmp_path: Path) -> None:
    ring = _ring(tmp_path)
    table = SessionResourceTable({"ring-1": _descriptor("first-upload.tsv", ring.read_bytes())})

    encoded = encode_canonical_request(_request(ring), table=table)

    assert _rings(encoded.payload) == ["ring-1"]
    assert encoded.payload["diagramOptions"]["conservationLabels"] == ["ring.tsv"]


def test_inputs_of_one_request_keep_resources_of_their_own(tmp_path: Path) -> None:
    rings = (_ring(tmp_path / "a"), _ring(tmp_path / "b"))
    held = _descriptor("ring.tsv", rings[0].read_bytes())
    table = SessionResourceTable({"ring-1": held})

    encoded = encode_canonical_request(_request(*rings), table=table)

    # Equal bytes are still two inputs: the first takes the held resource and
    # the second keeps its own, as when the request is encoded alone.
    assert _rings(encoded.payload) == ["ring-1", "conservation-blast-files-2"]
    assert table.descriptors()["conservation-blast-files-2"]["name"] == "ring.2.tsv"
    assert _rings(encode_canonical_request(_request(*rings)).payload) == [
        "conservation-blast-files-1",
        "conservation-blast-files-2",
    ]


def test_a_web_binding_shares_any_resource_with_equal_bytes() -> None:
    held = _descriptor("record-1-genbank-HmmtDNA.gbk", b"LOCUS same\n", kind="genbank")
    table = SessionResourceTable({"record-1-genbank": held})

    assert table.bind(_descriptor("HmmtDNA.gbk", b"LOCUS same\n")) == "record-1-genbank"
    assert table.bind(_descriptor("HmmtDNA.gbk", b"LOCUS same\n"), preferred_id="genbank-1") == (
        "record-1-genbank"
    )
    assert table.bind(_descriptor("new.gb", b"LOCUS new\n")) == "resource-0001"
    assert table.descriptors()["resource-0001"]["name"] == "resource-0001-new.gb"


def test_a_resave_drops_the_resources_nothing_names() -> None:
    source = build_session_document(_request()).to_dict()
    resources = {
        **source["resources"],
        "orphan": _descriptor("orphan.txt", b"orphan"),
    }

    document = _build_session_document_from_resolved_request(
        resolve_request(_request()),
        adjunct={"webFiles": {"resourceOriginalNames": {
            "record-1-genbank": "HmmtDNA.gbk",
            "orphan": "orphan.txt",
        }}},
        resources=resources,
    ).to_dict()

    assert document["resources"] == source["resources"]
    assert document["webFiles"]["resourceOriginalNames"] == {"record-1-genbank": "HmmtDNA.gbk"}


def test_reference_walker_reads_the_three_reference_keys_at_any_depth() -> None:
    payload = {
        "records": [
            {"source": {"kind": "genbank", "resourceId": "genbank-1"}},
            {"source": {
                "kind": "gffFasta",
                "gffResourceId": "gff-1",
                "fastaResourceId": "fasta-1",
            }},
        ],
        "diagramOptions": {
            "depthTracks": [[{"resourceId": "depth-1", "representation": "file"}, None]],
            "label": {"resource": "not-a-reference", "id": "not-a-reference"},
            "blank": {"resourceId": "  "},
        },
    }

    assert canonical_resource_ids(payload) == {"genbank-1", "gff-1", "fasta-1", "depth-1"}
    assert canonical_resource_ids(None) == set()


def test_a_request_reference_must_resolve() -> None:
    session = build_session_document(_request()).to_dict()
    session["resources"] = {}

    with pytest.raises(SessionResourceError, match="missing canonical resource.*record-1-genbank"):
        load_session_document(session)


# SHA-256 of build_session_document(request) for each fixture's decoded
# request, written by the per-request allocator before the Session resource
# table existed (3522c0bc): a table without seed resources must reproduce the
# single-request writer byte for byte. A codec change that changes a fixture's
# single-request output updates its digest here.
_SINGLE_REQUEST_DIGESTS = {
    "synthetic_conservation.gbdraw-session.json.gz":
        "5169d265c8f034382a97eb0ec59e10c6f88b0a6be8b547a29461a987d61bf275",
    "q-frame-main-web-upload.v42.gbdraw-session.json.gz":
        "fb04d983c55dc29a6a55fa189bf34fd2fabdc5bdfea51172c95682a341178bcd",
    "composite-circular-three-files.v44-schema8.gbdraw-session.json.gz":
        "9ecce263f8d25346d6f9a2e830f8ad8617d1645be0a8ae1aa1f9dc439a5a61ca",
    "feature-edits-crop-rc.v44.gbdraw-session.json.gz":
        "e82961ce4fecb9a430cdf6cf05ab90d84dad00b7b262cc5aec2e3f3c8a58cfc5",
    "selected-feature-annotations.v44.gbdraw-session.json.gz":
        "0931b4ee5ae98d5bd99ca7deb9400e0610449d68173c32eab5b85fdf1ca048ee",
}


@pytest.mark.parametrize("fixture", sorted(_SINGLE_REQUEST_DIGESTS))
def test_a_single_request_session_is_byte_identical(fixture: str, tmp_path: Path) -> None:
    with materialize_session(FIXTURES / fixture, output_directory=tmp_path) as materialized:
        for path in materialized.resource_paths.values():
            # A file resource stores its modification time.
            os.utime(path, (1_700_000_000, 1_700_000_000))
        request = session_to_request(materialized)
        document = build_session_document(
            request,
            created_at=datetime(2026, 1, 2, tzinfo=timezone.utc),
        ).to_dict()
    text = json.dumps(document, sort_keys=True, separators=(",", ":"))

    assert hashlib.sha256(text.encode("utf-8")).hexdigest() == _SINGLE_REQUEST_DIGESTS[fixture]
