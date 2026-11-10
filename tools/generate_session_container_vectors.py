#!/usr/bin/env python3
"""Write tests/fixtures/session_container_vectors.json, the Session file twin's vectors.

The Python owner (gbdraw/session_container.py) and the Web twin
(gbdraw/web/js/services/session-container.js) both read every ``read`` case
and must reach its outcome; both write every ``write`` case and must produce
its entries. The files are stored whole, base64, so each case pins the bytes
it was made from. Deflated bytes depend on the zlib build, so the script has
no --check: rerun it only to add or change a case, and review the outcomes.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import io
import json
import struct
import sys
import tempfile
import warnings
import zipfile
import zlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "tests" / "fixtures" / "session_container_vectors.json"
sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_container import (  # noqa: E402
    SESSION_CONTAINER_DAMAGED_MESSAGE,
    SESSION_CONTAINER_NO_DOCUMENT_MESSAGE,
    SESSION_DOCUMENT_ENTRY,
    BytesResourceSource,
    resource_entry_name,
    unsupported_session_version_message,
    write_session_container,
)

GENOME = b"LOCUS       demo 40 bp DNA linear\nORIGIN\n        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n"
READS = gzip.compress(b"chr1\t1\t5\nchr1\t2\t7\nchr1\t3\t6\n", compresslevel=9, mtime=0)
DATE_TIME = (1980, 1, 1, 0, 0, 0)


def _document(
    resources: dict[str, bytes], sizes: dict[str, int] | None = None, *, embed: bool = False,
) -> dict[str, object]:
    """A Session 51 document declaring ``resources``; ``embed`` also puts each one's bytes in its
    descriptor as base64, as a JSON Session does."""

    return {
        "format": "gbdraw-session",
        "version": 51,
        "createdAt": "2026-11-08T00:00:00.000Z",
        "resources": {
            resource_id: {
                "kind": "web-file",
                "name": f"{resource_id}.dat",
                "type": "application/octet-stream",
                "size": (sizes or {}).get(resource_id, len(data)),
                "lastModified": 0,
                **({"encoding": "base64", "data": base64.b64encode(data).decode("ascii")} if embed else {}),
            }
            for resource_id, data in resources.items()
        },
        "ui": {},
        "drawings": [],
    }


def _json(document: dict[str, object]) -> bytes:
    return json.dumps(document, ensure_ascii=False, separators=(",", ":")).encode("utf-8")


def _python_written(resources: dict[str, bytes]) -> bytes:
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "case.gbdraw"
        write_session_container(
            _document(resources), {key: BytesResourceSource(value) for key, value in resources.items()}, path,
        )
        return path.read_bytes()


class _Unseekable(io.RawIOBase):
    """A write-only stream: zipfile then writes data descriptors, as a streaming writer does."""

    def __init__(self) -> None:
        super().__init__()
        self.buffer = io.BytesIO()

    def writable(self) -> bool:
        return True

    def write(self, data: object) -> int:  # type: ignore[override]
        return self.buffer.write(data)  # type: ignore[arg-type]


def _zip(entries: list[tuple[str, bytes, int]], *, force_zip64: bool = False, unseekable: bool = False) -> bytes:
    output: io.BytesIO | _Unseekable = _Unseekable() if unseekable else io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # duplicate names are a case
        with zipfile.ZipFile(output, "w", allowZip64=True) as archive:  # type: ignore[arg-type]
            for name, data, method in entries:
                info = zipfile.ZipInfo(name, date_time=DATE_TIME)
                info.compress_type = method
                info.external_attr = 0o100644 << 16
                with archive.open(info, "w", force_zip64=force_zip64) as entry:
                    entry.write(data)
    return output.buffer.getvalue() if isinstance(output, _Unseekable) else output.getvalue()


def _session_entries(resources: dict[str, bytes], document: bytes | None = None) -> list[tuple[str, bytes, int]]:
    entries = [(SESSION_DOCUMENT_ENTRY, document if document is not None else _json(_document(resources)),
                zipfile.ZIP_DEFLATED)]
    for resource_id, data in resources.items():
        method = zipfile.ZIP_STORED if data[:2] == b"\x1f\x8b" else zipfile.ZIP_DEFLATED
        entries.append((resource_entry_name(resource_id), data, method))
    return entries


def _patch_entry(data: bytes, name: str, *, flags: int = 0, file_size: int | None = None) -> bytes:
    """Set general purpose flags or the uncompressed size of ``name`` in its local and central headers."""

    patched = bytearray(data)
    info = zipfile.ZipFile(io.BytesIO(data)).getinfo(name)
    local = info.header_offset
    central = data.index(b"PK\x01\x02")
    while data[central + 46: central + 46 + len(name)] != name.encode():
        central = data.index(b"PK\x01\x02", central + 4)
    for offset, flag_at, size_at in ((local, 6, 22), (central, 8, 24)):
        if flags:
            (current,) = struct.unpack_from("<H", patched, offset + flag_at)
            struct.pack_into("<H", patched, offset + flag_at, current | flags)
        if file_size is not None:
            struct.pack_into("<I", patched, offset + size_at, file_size)
    return bytes(patched)


def _damaged_deflate(data: bytes, name: str) -> bytes:
    """``data`` with the deflated data of entry ``name`` starting with an invalid block type."""

    info = zipfile.ZipFile(io.BytesIO(data)).getinfo(name)
    assert info.compress_type == zipfile.ZIP_DEFLATED
    name_length, extra_length = struct.unpack_from("<HH", data, info.header_offset + 26)
    start = info.header_offset + 30 + name_length + extra_length
    return data[:start] + b"\xff\xff" + data[start + 2:]


def _renamed(data: bytes, name: str, raw_name: bytes) -> bytes:
    """``data`` with entry ``name`` renamed to ``raw_name``, of the same length, in its local and central headers.

    ``zipfile`` writes neither a NUL in a name nor, on Windows, a backslash, so the bytes are replaced.
    """

    old = name.encode()
    assert len(old) == len(raw_name) and data.count(old) == 2
    return data.replace(old, raw_name)


def _entry_count(data: bytes, count: int) -> bytes:
    """``data`` whose end record declares ``count`` entries."""

    patched = bytearray(data)
    struct.pack_into("<HH", patched, data.rindex(b"PK\x05\x06") + 8, count, count)
    return bytes(patched)


def _replace_last_data(data: bytes, name: str, compressed: bytes) -> bytes:
    """``data`` whose last entry, ``name``, holds ``compressed`` as its data, its sizes and the directory moved to fit."""

    info = zipfile.ZipFile(io.BytesIO(data)).getinfo(name)
    name_length, extra_length = struct.unpack_from("<HH", data, info.header_offset + 26)
    start = info.header_offset + 30 + name_length + extra_length
    assert data[start + info.compress_size: start + info.compress_size + 4] == b"PK\x01\x02"
    patched = bytearray(data[:start] + compressed + data[start + info.compress_size:])
    shift = len(compressed) - info.compress_size
    struct.pack_into("<I", patched, info.header_offset + 18, len(compressed))
    central = data.index(b"PK\x01\x02", start) + shift
    while patched[central + 46: central + 46 + len(name)] != name.encode():
        central = patched.index(b"PK\x01\x02", central + 4)
    struct.pack_into("<I", patched, central + 20, len(compressed))
    end = patched.rindex(b"PK\x05\x06")
    (offset,) = struct.unpack_from("<I", patched, end + 16)
    struct.pack_into("<I", patched, end + 16, offset + shift)
    return bytes(patched)


def _raw_deflate(data: bytes, *, finish: bool = True) -> bytes:
    deflater = zlib.compressobj(9, zlib.DEFLATED, -zlib.MAX_WBITS)
    return deflater.compress(data) + deflater.flush(zlib.Z_FINISH if finish else zlib.Z_SYNC_FLUSH)


def _read_outcome(data: bytes) -> dict[str, object]:
    """What a reader reaches: the document and every resource it declares, read in full."""

    from gbdraw.session_container import SessionContainerError, read_session_container
    from gbdraw.session_io import check_session_file_version

    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "case.gbdraw"
        path.write_bytes(data)
        try:
            container = read_session_container(path)
            check_session_file_version(container)
            document = zipfile.ZipFile(path).read(SESSION_DOCUMENT_ENTRY)
            resources: dict[str, object] = {}
            for resource_id in container.resource_ids():
                with container.resource(resource_id).open_bytes() as stream:
                    content = stream.read()
                resources[resource_id] = {"size": len(content), "sha256": hashlib.sha256(content).hexdigest()}
        except SessionContainerError as error:
            return {"ok": False, "error": error.code}
    return {
        "ok": True,
        "document": {"size": len(document), "sha256": hashlib.sha256(document).hexdigest()},
        "resources": resources,
    }


def build() -> dict[str, object]:
    resources = {"genome": GENOME, "reads": READS, "empty": b""}
    written = _python_written(resources)
    stored = zipfile.ZipFile(io.BytesIO(written)).getinfo("resources/reads")
    name_length, extra_length = struct.unpack_from("<HH", written, stored.header_offset + 26)
    flip_at = stored.header_offset + 30 + name_length + extra_length + 12
    crc_broken = written[:flip_at] + bytes([written[flip_at] ^ 0xFF]) + written[flip_at + 1:]
    end = written.rindex(b"PK\x05\x06")
    (central_size, central_offset) = struct.unpack_from("<II", written, end + 12)
    no_central = written[:central_offset] + b"\0" * central_size + written[central_offset + central_size:]
    size_document = _json(_document(resources, {"reads": len(READS) + 1}))
    size_broken = _patch_entry(
        _zip(_session_entries(resources, size_document)), "resources/reads", file_size=len(READS) + 1,
    )
    read_cases = [
        ("python-written", "Written by write_session_container: a plain, a gzip (stored) and an empty resource.",
         written),
        ("zip64", "Every entry with ZIP64 fields, as a writer makes them for entries over 4 GiB.",
         _zip(_session_entries(resources), force_zip64=True)),
        ("data-descriptors", "Sizes and CRC-32 in data descriptors after each entry, as a streaming writer puts them.",
         _zip(_session_entries(resources), unseekable=True)),
        ("hand-rezipped", "Zipped again by hand: resources first, a directory entry and a file the document does "
         "not name; such entries are not read.",
         _zip([("resources/", b"", zipfile.ZIP_STORED), *_session_entries(resources)[1:],
               ("notes.txt", b"not a resource", zipfile.ZIP_DEFLATED), _session_entries(resources)[0]])),
        ("dot-dot-entries", "Entries named ../x and resources/../genome: resources are found by ID, so they are "
         "not read and never become paths.",
         _zip([*_session_entries(resources), ("../x", b"outside", zipfile.ZIP_STORED),
               ("resources/../genome", b"outside", zipfile.ZIP_STORED)])),
        ("truncated", "The Python-written file without its last 200 bytes: no end record.", written[:-200]),
        ("central-directory-missing", "The central directory overwritten with zeros; the end record still points "
         "at it.", no_central),
        ("crc-mismatch", "One byte of the stored gzip resource flipped.", crc_broken),
        ("size-mismatch", "The stored resource declares one byte more, in its headers and in the document, than "
         "its entry holds.", size_broken),
        ("declared-size-mismatch", "The document declares another size than the resource's entry has.",
         _zip(_session_entries(resources, _json(_document(resources, {"genome": len(GENOME) - 1}))))),
        ("missing-resource-entry", "The document declares a resource that has no entry.",
         _zip(_session_entries(resources)[:-1])),
        ("duplicate-entry", "Two entries named resources/genome.",
         _zip([*_session_entries(resources), ("resources/genome", GENOME, zipfile.ZIP_DEFLATED)])),
        ("encrypted-entry", "The plain resource's entry is marked encrypted.",
         _patch_entry(written, "resources/genome", flags=0x0001)),
        ("unsupported-method", "The plain resource is compressed with bzip2 (method 12).",
         _zip([*_session_entries(resources)[:1], ("resources/genome", GENOME, zipfile.ZIP_BZIP2),
               *_session_entries(resources)[2:]])),
        ("end-record-count", "The end record declares one entry more than the directory holds.",
         _entry_count(written, 5)),
        ("nul-in-resource-name", "The plain resource's entry is named resources/genome, a NUL and junk: its whole "
         "name is not resources/genome, so the resource has no entry.",
         _renamed(_zip([*_session_entries(resources)[:1], ("resources/genomeXjunk", GENOME, zipfile.ZIP_DEFLATED),
                        *_session_entries(resources)[2:]]),
                  "resources/genomeXjunk", b"resources/genome\x00junk")),
        ("nul-name-beside-resource", "Beside resources/genome, an entry named resources/genome, a NUL and junk: two "
         "whole names, so not a duplicate, and the second is not read.",
         _renamed(_zip([*_session_entries(resources), ("resources/genomeXjunk", b"other", zipfile.ZIP_DEFLATED)]),
                  "resources/genomeXjunk", b"resources/genome\x00junk")),
        ("nul-in-document-name", "The document's entry is named gbdraw-session.json, a NUL and junk: the ZIP has no "
         "gbdraw-session.json.",
         _renamed(_zip([(SESSION_DOCUMENT_ENTRY + "Xjunk", _json(_document(resources)), zipfile.ZIP_DEFLATED),
                        *_session_entries(resources)[1:]]),
                  SESSION_DOCUMENT_ENTRY + "Xjunk", SESSION_DOCUMENT_ENTRY.encode() + b"\x00junk")),
        ("backslash-resource-name", "The plain resource's entry is named resources\\genome: a resource is found by "
         "the name resources/genome only, so it has no entry.",
         _zip([*_session_entries(resources)[:1], ("resources\\genome", GENOME, zipfile.ZIP_DEFLATED),
               *_session_entries(resources)[2:]])),
        ("deflate-trailing-bytes", "Bytes after the end of the plain resource's deflate stream, inside its "
         "compressed size.",
         _replace_last_data(_zip([*_session_entries(resources)[:1], *_session_entries(resources)[2:],
                                  _session_entries(resources)[1]]),
                            "resources/genome", _raw_deflate(GENOME) + b"JUNK")),
        ("deflate-unfinished", "The plain resource's deflate stream has no final block: its data ends before the "
         "stream does, though its size and CRC-32 match.",
         _replace_last_data(_zip([*_session_entries(resources)[:1], *_session_entries(resources)[2:],
                                  _session_entries(resources)[1]]),
                            "resources/genome", _raw_deflate(GENOME, finish=False))),
        ("embedded-resource-bytes", "Every descriptor of gbdraw-session.json also embeds its bytes as base64 "
         "(encoding and data), as a JSON Session does; the entries are kept.",
         _zip(_session_entries(resources, _json(_document(resources, embed=True))))),
        ("damaged-deflate-document", "The deflated data of gbdraw-session.json starts with an invalid block type: "
         "it cannot be inflated.", _damaged_deflate(written, SESSION_DOCUMENT_ENTRY)),
        ("damaged-deflate-resource", "The deflated data of the plain resource starts with an invalid block type: "
         "it cannot be inflated when the resource is read.", _damaged_deflate(written, "resources/genome")),
        ("older-document-version", "gbdraw-session.json is a Session 50: only a Session 51 was ever written as a "
         "container.", _zip(_session_entries(resources, _json({**_document(resources), "version": 50})))),
        ("no-document", "A ZIP with resources and no gbdraw-session.json.", _zip(_session_entries(resources)[1:])),
        ("empty-zip", "A ZIP with no entries: only an end record.", _zip([])),
        ("invalid-json", "gbdraw-session.json is not JSON.", _zip(_session_entries(resources, b'{"format":'))),
        ("duplicate-json-key", "gbdraw-session.json names the key format twice: refused as not valid JSON.",
         _zip(_session_entries(resources, b'{"format":"gbdraw-session",' + _json(_document(resources))[1:]))),
    ]
    return {
        "schema": 1,
        "description": (
            "Vectors of the Session file container (.gbdraw). read: a file (base64) and what a reader reaches when "
            "it reads the document and every resource the document declares, in full: the document entry's size "
            "and SHA-256 and each resource's, or the refusal code. write: a document and its resources (base64); "
            "a writer writes these entries in this order, each with this compression method (0 stored, 8 deflated), "
            "and reads back the same bytes. Generated by tools/generate_session_container_vectors.py."
        ),
        "messages": {
            "damaged": {"python": SESSION_CONTAINER_DAMAGED_MESSAGE, "web": SESSION_CONTAINER_DAMAGED_MESSAGE},
            "no-document": {
                "python": SESSION_CONTAINER_NO_DOCUMENT_MESSAGE, "web": SESSION_CONTAINER_NO_DOCUMENT_MESSAGE,
            },
            "invalid-json": {
                "python": "Not a valid JSON session file: {path}",
                "web": "This file is not a gbdraw Session: it is not valid JSON.",
            },
            "unsupported-version": {
                "python": unsupported_session_version_message("{version}"),
                # The Web's refusal of a JSON Session 51 too (PD-OI-127).
                "web": "An input value is invalid. This gbdraw does not read this Session version. "
                "Load a Session file that this version of gbdraw supports.",
            },
        },
        "read": [
            {"id": case_id, "description": description, "file": base64.b64encode(data).decode("ascii"),
             "expect": _read_outcome(data)}
            for case_id, description, data in read_cases
        ],
        "write": [
            {
                "id": "plain-gzip-empty",
                "description": "A plain resource is deflated, a gzip resource stored, an empty resource deflated.",
                "document": _document(resources),
                "resources": {key: base64.b64encode(value).decode("ascii") for key, value in resources.items()},
                "expect": {
                    "entries": [[name, method] for name, _, method in _session_entries(resources)],
                },
            },
        ],
    }


def main() -> None:
    TARGET.write_text(json.dumps(build(), indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
