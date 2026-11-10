"""The Session resource table: one allocator for every Python Session writer."""

from __future__ import annotations

import base64
import hashlib
import json
import os
import shutil
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
    _mode_scoped_document_of_request,
    build_session_document,
    load_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.session_request_codec import CanonicalRequestResource, encode_canonical_request
from gbdraw.session_resources import SessionResourceTable, _resource_descriptor, canonical_resource_ids

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
    # D-03: every precomputed ring stores its drawn label (the file stem).
    assert alone.payload["diagramOptions"]["conservationLabels"] == ["ring"]
    table = SessionResourceTable({"ring-1": _descriptor("ring.tsv", b"other rows\n")})

    encoded = encode_canonical_request(_request(ring), table=table)

    assert _rings(encoded.payload) == ["conservation-blast-files-1"]
    assert table.descriptors()["conservation-blast-files-1"]["name"] == "ring.2.tsv"
    # The stored label keeps the drawn one although the resource is renamed.
    assert encoded.payload["diagramOptions"]["conservationLabels"] == ["ring"]


def test_equal_bytes_under_another_name_pin_ring_labels(tmp_path: Path) -> None:
    ring = _ring(tmp_path)
    table = SessionResourceTable({"ring-1": _descriptor("first-upload.tsv", ring.read_bytes())})

    encoded = encode_canonical_request(_request(ring), table=table)

    assert _rings(encoded.payload) == ["ring-1"]
    assert encoded.payload["diagramOptions"]["conservationLabels"] == ["ring"]


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
# single-request output updates its digest here; so does a Session version
# bump (46, 47, 48, 49, 50: only the version field changed; this checked document
# is the mode-scoped Session the writer produced before drawings, so it stays at
# the last mode-scoped version, 50, when the writer moves on to 51) and a change of the packaged
# config a stored full config is completed with (the label style keys, LP-D98;
# the Linear box default, LP-D106; the Circular track narrowing leaves, PD-OI-078;
# the Circular definition show leaf, PD-OI-109; the ring colour scheme leaf, PD-OI-110;
# the italic gene label leaves; the identity colour scale leaves; the figure size leaves, #1427;
# the unused Linear leaf canvas.linear.canvas_padding, removed; the scalar axis label
# style leaves, which a stored full config or configOverrides without them reads as
# numbered; request schema 10 writing record windows and each similarity ring's comparison
# genome as a list of files, so that fixture's three rings are lists of one file each; the file
# metadata gbdraw writes, a media type by resource kind and lastModified 0, OV-2964).
_SINGLE_REQUEST_DIGESTS = {
    "synthetic_conservation.gbdraw-session.json.gz":
        "300759bb064f75d913f97bef3ee67f4ffcd5e43e2ed9d553bb9e2edd6fa98472",
    "q-frame-main-web-upload.v42.gbdraw-session.json.gz":
        "b359ca523356fa2de78495c140f8dfbc512641b78c391f1d3a9470592dd5f50d",
    "composite-circular-three-files.v44-schema8.gbdraw-session.json.gz":
        "49f804c719e50914ca66c84f41ee76e2b332495564a3bb12f3f3ed003da6fab8",
    "feature-edits-crop-rc.v44.gbdraw-session.json.gz":
        "31a911105b1c256269a62125407195658853db749ece05f1149b3fd93222e5e2",
    "selected-feature-annotations.v44.gbdraw-session.json.gz":
        "e6eaf5cfcec72f56b446ed904d207fd70467ebe5c5d8fcef22da9fc1305b5e2f",
}


@pytest.mark.parametrize("fixture", sorted(_SINGLE_REQUEST_DIGESTS))
def test_a_single_request_session_is_byte_identical(fixture: str, tmp_path: Path) -> None:
    with materialize_session(FIXTURES / fixture, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        # The checked mode-scoped Session 48, before the move to drawings.
        document = _mode_scoped_document_of_request(
            request,
            resolved=False,
            created_at=datetime(2026, 1, 2, tzinfo=timezone.utc),
        ).to_dict()
    text = json.dumps(document, sort_keys=True, separators=(",", ":"))

    assert hashlib.sha256(text.encode("utf-8")).hexdigest() == _SINGLE_REQUEST_DIGESTS[fixture]


class _FixedNow(datetime):
    @classmethod
    def now(cls, tz: Any = None) -> "_FixedNow":
        return cls(2026, 1, 2, tzinfo=timezone.utc)


def _file_leaves(value: object) -> list[dict[str, Any]]:
    """Every ordinary file binding (a leaf with a ``resourceId``) under ``value``."""

    if isinstance(value, list):
        return [leaf for item in value for leaf in _file_leaves(item)]
    if not isinstance(value, dict):
        return []
    if "resourceId" in value:
        return [value]
    return [leaf for item in value.values() for leaf in _file_leaves(item)]


def test_a_cli_session_does_not_depend_on_its_inputs_file_metadata(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """OV-2964: the CLI writes the same Session bytes whatever its inputs' modification times and the
    machine's media type guesses; each file binding says of its file what its resource says."""

    from gbdraw.circular import circular_main
    import gbdraw.cli_utils.session as cli_session

    monkeypatch.setattr(cli_session, "datetime", _FixedNow)
    genbank = tmp_path / "inputs" / "record.gb"
    genbank.parent.mkdir()
    genbank.write_bytes(GENBANK.read_bytes())
    written = []
    for index, modified in enumerate((1_700_000_000, 1_800_000_000)):
        os.utime(genbank, (modified, modified))
        shutil.rmtree(tmp_path / "out", ignore_errors=True)
        session = tmp_path / f"run-{index}.gbdraw-session.json"
        circular_main([
            "--gbk", str(genbank), "-o", str(tmp_path / "out" / "fig"),
            "--session_output", str(session),
        ])
        written.append(session.read_bytes())

    assert written[0] == written[1]
    document = json.loads(written[0])
    resources = document["resources"]
    assert {(resource["type"], resource["lastModified"]) for resource in resources.values()} == {
        ("application/octet-stream", 0)
    }
    leaves = [leaf for drawing in document["drawings"] for leaf in _file_leaves(drawing.get("webFiles", {}))]
    assert leaves
    for leaf in leaves:
        assert (leaf["type"], leaf["lastModified"]) == (
            resources[leaf["resourceId"]]["type"], resources[leaf["resourceId"]]["lastModified"]
        )


@pytest.mark.parametrize(
    ("kind", "media_type"),
    [
        ("genbank", "application/octet-stream"),
        ("nucleotide-blast", "application/octet-stream"),
        ("canonical-tsv", "text/tab-separated-values"),
        ("collinearity-result", "application/json"),
        ("orthogroup-result", "application/json"),
    ],
)
def test_a_resource_gbdraw_writes_has_the_media_type_of_its_kind(kind: str, media_type: str) -> None:
    """OV-2964: a resource of a kind gbdraw makes itself has the media type the Web writers write for it,
    any input file ``application/octet-stream``, and no modification time."""

    resource = CanonicalRequestResource("resource-1", kind, "resource-1.json", content=b"{}")

    descriptor = _resource_descriptor(resource)

    assert (descriptor["type"], descriptor["lastModified"]) == (media_type, 0)
