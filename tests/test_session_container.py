"""The Session file container: the shared vectors, the writer's layout and its staging."""

from __future__ import annotations

import base64
import hashlib
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.render.output_paths import OutputExistsError
from gbdraw.session_container import (
    SESSION_CONTAINER_DAMAGED_MESSAGE,
    SESSION_CONTAINER_NO_DOCUMENT_MESSAGE,
    SESSION_DOCUMENT_ENTRY,
    BytesResourceSource,
    FileResourceSource,
    SessionContainerError,
    is_session_container,
    read_session_container,
    write_session_container,
)

VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "session_container_vectors.json").read_text(encoding="utf-8")
)
_DATA_DESCRIPTOR_FLAG = 0x0008


def _outcome(path: Path) -> dict[str, object]:
    try:
        container = read_session_container(path)
        with zipfile.ZipFile(path) as archive:
            document = archive.read(SESSION_DOCUMENT_ENTRY)
        resources: dict[str, object] = {}
        for resource_id in container.resource_ids():
            with container.resource(resource_id).open() as stream:
                content = stream.read()
            resources[resource_id] = {"size": len(content), "sha256": hashlib.sha256(content).hexdigest()}
    except SessionContainerError as error:
        return {"ok": False, "error": error.code, "message": str(error)}
    return {
        "ok": True,
        "document": {"size": len(document), "sha256": hashlib.sha256(document).hexdigest()},
        "resources": resources,
    }


@pytest.mark.parametrize("case", VECTORS["read"], ids=lambda case: case["id"])
def test_reader_reaches_each_vector_outcome(case: dict, tmp_path: Path) -> None:
    path = tmp_path / "case.gbdraw"
    path.write_bytes(base64.b64decode(case["file"]))
    outcome = _outcome(path)
    if case["expect"]["ok"]:
        assert outcome == case["expect"]
    else:
        message = VECTORS["messages"][case["expect"]["error"]]["python"].format(path=path)
        assert outcome == {**case["expect"], "message": message}


def test_vector_messages_are_the_owner_texts() -> None:
    assert VECTORS["messages"]["damaged"] == {
        "python": SESSION_CONTAINER_DAMAGED_MESSAGE, "web": SESSION_CONTAINER_DAMAGED_MESSAGE,
    }
    assert VECTORS["messages"]["no-document"] == {
        "python": SESSION_CONTAINER_NO_DOCUMENT_MESSAGE, "web": SESSION_CONTAINER_NO_DOCUMENT_MESSAGE,
    }


def _write_case(case: dict, path: Path) -> dict[str, bytes]:
    resources = {key: base64.b64decode(value) for key, value in case["resources"].items()}
    write_session_container(
        case["document"], {key: BytesResourceSource(value) for key, value in resources.items()}, path,
    )
    return resources


@pytest.mark.parametrize("case", VECTORS["write"], ids=lambda case: case["id"])
def test_writer_writes_each_vector_layout(case: dict, tmp_path: Path) -> None:
    path = tmp_path / "written.gbdraw"
    resources = _write_case(case, path)
    with zipfile.ZipFile(path) as archive:
        infos = archive.infolist()
        assert [[info.filename, info.compress_type] for info in infos] == case["expect"]["entries"]
        # The Python writer writes to a seekable file: sizes in the headers, no data descriptors.
        assert not any(info.flag_bits & _DATA_DESCRIPTOR_FLAG for info in infos)
        assert json.loads(archive.read(SESSION_DOCUMENT_ENTRY)) == case["document"]
    container = read_session_container(path)
    assert container.document == case["document"]
    for resource_id, content in resources.items():
        with container.resource(resource_id).open() as stream:
            assert stream.read() == content
    # The first entry's name follows the 30-byte local header, so a tool can sniff the file.
    head = path.read_bytes()[:64]
    assert is_session_container(head)
    assert head[30:30 + len(SESSION_DOCUMENT_ENTRY)] == SESSION_DOCUMENT_ENTRY.encode()


@pytest.mark.skipif(shutil.which("unzip") is None, reason="unzip is not installed")
def test_unzip_tests_a_written_file(tmp_path: Path) -> None:
    path = tmp_path / "written.gbdraw"
    _write_case(VECTORS["write"][0], path)
    result = subprocess.run(["unzip", "-t", str(path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


def test_writer_adds_zip64_fields_when_an_entry_needs_them(tmp_path: Path, monkeypatch) -> None:
    case = VECTORS["write"][0]
    path = tmp_path / "zip64.gbdraw"
    monkeypatch.setattr(zipfile, "ZIP64_LIMIT", 64)
    resources = _write_case(case, path)
    monkeypatch.undo()
    with zipfile.ZipFile(path) as archive:
        versions = {info.filename: info.extract_version for info in archive.infolist()}
    assert versions[SESSION_DOCUMENT_ENTRY] >= zipfile.ZIP64_VERSION
    container = read_session_container(path)
    for resource_id, content in resources.items():
        with container.resource(resource_id).open() as stream:
            assert stream.read() == content


def test_writer_streams_files_and_keeps_gzip_bytes(tmp_path: Path) -> None:
    case = VECTORS["write"][0]
    sources = {}
    for key, value in case["resources"].items():
        source_path = tmp_path / f"{key}.input"
        source_path.write_bytes(base64.b64decode(value))
        sources[key] = FileResourceSource(source_path)
    path = tmp_path / "files.gbdraw"
    write_session_container(case["document"], sources, path)
    with zipfile.ZipFile(path) as archive:
        assert archive.getinfo("resources/reads").compress_type == zipfile.ZIP_STORED
        assert archive.read("resources/reads") == base64.b64decode(case["resources"]["reads"])


def test_writer_refuses_bytes_of_another_size_and_writes_nothing(tmp_path: Path) -> None:
    case = VECTORS["write"][0]
    path = tmp_path / "changed.gbdraw"
    sources = {key: BytesResourceSource(base64.b64decode(value)) for key, value in case["resources"].items()}
    sources["genome"] = BytesResourceSource(b"changed since it was declared")
    with pytest.raises(ValidationError, match="Invalid embedded resource bytes or byte size."):
        write_session_container(case["document"], sources, path)
    assert list(tmp_path.iterdir()) == []


def test_writer_refuses_resources_the_document_does_not_declare(tmp_path: Path) -> None:
    case = VECTORS["write"][0]
    with pytest.raises(ValidationError, match="do not match"):
        write_session_container(case["document"], {"genome": BytesResourceSource(b"")}, tmp_path / "x.gbdraw")


def test_writer_keeps_an_existing_file_unless_told_to_replace_it(tmp_path: Path) -> None:
    path = tmp_path / "kept.gbdraw"
    path.write_bytes(b"kept")
    with pytest.raises(OutputExistsError):
        _write_case_no_overwrite(VECTORS["write"][0], path)
    assert path.read_bytes() == b"kept"
    assert list(tmp_path.iterdir()) == [path]


def _write_case_no_overwrite(case: dict, path: Path) -> None:
    resources = {key: BytesResourceSource(base64.b64decode(value)) for key, value in case["resources"].items()}
    write_session_container(case["document"], resources, path, overwrite=False)


def test_entries_are_found_by_resource_id_only(tmp_path: Path) -> None:
    case = next(case for case in VECTORS["read"] if case["id"] == "dot-dot-entries")
    folder = tmp_path / "inside"
    folder.mkdir()
    path = folder / "case.gbdraw"
    path.write_bytes(base64.b64decode(case["file"]))
    container = read_session_container(path)
    assert container.resource_ids() == ("genome", "reads", "empty")
    with pytest.raises(SessionContainerError):
        container.resource("../x")
    for resource_id in container.resource_ids():
        with container.resource(resource_id).open() as stream:
            assert stream.read() != b"outside"
    assert sorted(item.name for item in tmp_path.rglob("*")) == ["case.gbdraw", "inside"]


def test_content_decides_whether_a_file_is_a_container() -> None:
    assert is_session_container(b"PK\x03\x04rest")
    assert is_session_container(b"PK\x05\x06")
    assert not is_session_container(b"\x1f\x8b\x08\x00")
    assert not is_session_container(b'{"format":')
