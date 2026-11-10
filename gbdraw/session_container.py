"""The Session file: one ZIP container, written and read here.

A Session file holds the Session document as its first entry,
``gbdraw-session.json``, and one entry ``resources/<id>`` per resource of the
document's ``resources``, in that order. An input that is gzip already is
stored as it is; the document and every other input are deflated. The
document names no entry: a resource's entry is found by its ID only, and an
entry name never becomes a path. The ZIP's CRC-32 and sizes are the checks;
the container carries no ``checksum``.

This module is the one owner of that layout, of the checks a reader makes
and of the texts it refuses with. The Web app reads and writes the same file
with a registered twin (``services/session-container.js``);
``tests/fixtures/session_container_vectors.json`` pins both to each other.
"""

from __future__ import annotations

import io
import json
import os
import struct
import tempfile
import zipfile
import zlib
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, BinaryIO, Literal, Protocol, runtime_checkable

from .exceptions import ValidationError
from .render.output_paths import OutputExistsError, commit_staged_output_file

if TYPE_CHECKING:
    from _typeshed import WriteableBuffer

# The entry that holds the Session document.
SESSION_DOCUMENT_ENTRY = "gbdraw-session.json"
# A resource's entry is this prefix and the resource ID.
SESSION_RESOURCE_ENTRY_PREFIX = "resources/"
# The compression methods a Session entry may use: stored and deflated.
SESSION_ENTRY_STORED = zipfile.ZIP_STORED
SESSION_ENTRY_DEFLATED = zipfile.ZIP_DEFLATED
SESSION_ENTRY_METHODS = (SESSION_ENTRY_STORED, SESSION_ENTRY_DEFLATED)
# The first bytes of a gzip stream: an input that starts with them is stored.
GZIP_SIGNATURE = b"\x1f\x8b"
# The first bytes of a ZIP file: a local entry, or the end record of an empty ZIP.
SESSION_CONTAINER_SIGNATURES = (b"PK\x03\x04", b"PK\x05\x06")

SESSION_CONTAINER_DAMAGED_MESSAGE = "This file is not a gbdraw Session: it is incomplete or damaged."
SESSION_CONTAINER_NO_DOCUMENT_MESSAGE = "This file is not a gbdraw Session: it has no gbdraw-session.json."
# What the Web names each refusal of a Session file by: one it cannot read, its
# JSON, and resources that do not match their bytes.
SESSION_FILE_UNREADABLE = {"code": "INPUT_UNREADABLE"}
SESSION_JSON_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "JSON_FORMAT"}
SESSION_RESOURCES_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}

# Every entry gets the same time, system and permissions, so the bytes depend
# on the content and the zlib build only. The Web writer gives its entries the
# same time (year, month, day, hour, minute, second).
SESSION_ENTRY_DATE_TIME = (1980, 1, 1, 0, 0, 0)
_UNIX_SYSTEM = 3
_FILE_ATTRIBUTES = 0o100644 << 16
# General purpose flags: bit 0 encrypted, bit 6 strong encryption; an entry
# with either is refused as damaged.
SESSION_ENTRY_ENCRYPTED_FLAGS = 0x0041
# General purpose flag bit 11: the entry's name is UTF-8 (else code page 437).
_UTF8_NAME_FLAG = 0x0800
_CHUNK_BYTES = 1 << 20
_LOCAL_HEADER_SIGNATURE = b"PK\x03\x04"
_LOCAL_HEADER_BYTES = 30
_END_RECORD_SIGNATURE = b"PK\x05\x06"
_END_RECORD_BYTES = 22
_MAX_COMMENT_BYTES = 0xFFFF
_ZIP64_LOCATOR_SIGNATURE = b"PK\x06\x07"
_ZIP64_LOCATOR_BYTES = 20
_ZIP64_END_BYTES = 56
# What zipfile raises for a file it cannot read as a ZIP, besides OSError;
# zlib.error is a deflated entry whose data is damaged.
_ZIP_ERRORS = (
    zipfile.BadZipFile, zipfile.LargeZipFile, EOFError, ValueError, NotImplementedError, struct.error, zlib.error,
)

SessionContainerErrorCode = Literal["damaged", "no-document", "invalid-json", "unsupported-version"]


class SessionContainerError(ValidationError):
    """A Session file that cannot be read: a damaged container, one without a
    document, a document that is not JSON (container or JSON Session), or a
    container whose document is of a version older than the current one.

    ``code`` names the refusal; the vectors give the text of each per surface.
    """

    def __init__(self, code: SessionContainerErrorCode, path: Path, *, version: int | None = None) -> None:
        super().__init__(
            SESSION_CONTAINER_NO_DOCUMENT_MESSAGE if code == "no-document"
            else invalid_session_json_message(path) if code == "invalid-json"
            else unsupported_session_version_message(version) if code == "unsupported-version"
            else SESSION_CONTAINER_DAMAGED_MESSAGE
        )
        self.code: SessionContainerErrorCode = code


def unsupported_session_version_message(version: object) -> str:
    """The refusal of a Session of a version this gbdraw does not read."""

    return f"Unsupported session version: {version}."


def invalid_session_json_message(path: str | Path) -> str:
    """The refusal of a Session file whose JSON is not readable."""

    return f"Not a valid JSON session file: {path}"


def reject_duplicate_json_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON object keys before the decoder can discard them."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError(f"Session JSON contains a duplicate object key: {key!r}.", diagnostic=SESSION_JSON_INVALID)
        result[key] = value
    return result


def parse_session_json(text: str, path: str | Path) -> object:
    """The JSON value of a Session document; unreadable JSON is refused as such."""

    try:
        return json.loads(text, object_pairs_hook=reject_duplicate_json_keys)
    except json.JSONDecodeError as exc:
        raise SessionContainerError("invalid-json", Path(path)) from exc


@runtime_checkable
class SessionResourceSource(Protocol):
    """The bytes of one resource, read once from the start when the file is written."""

    def open_bytes(self) -> BinaryIO: ...


@dataclass(frozen=True)
class FileResourceSource:
    """A resource whose bytes are a local file, as they are on disk."""

    path: Path

    def open_bytes(self) -> BinaryIO:
        return open(self.path, "rb")


@dataclass(frozen=True)
class BytesResourceSource:
    """A resource whose bytes are in memory."""

    data: bytes

    def open_bytes(self) -> BinaryIO:
        return io.BytesIO(self.data)


def is_session_container(head: bytes) -> bool:
    """Whether a file that starts with ``head`` is read as a Session container.

    The content decides, never the name: anything else is read as a JSON
    Session, plain or gzip.
    """

    return head[:4] in SESSION_CONTAINER_SIGNATURES


def resource_entry_name(resource_id: str) -> str:
    """The entry of a resource: the prefix and its ID."""

    return SESSION_RESOURCE_ENTRY_PREFIX + resource_id


def write_session_container(
    document: Mapping[str, object],
    resources: Mapping[str, SessionResourceSource],
    path: str | Path,
    *,
    overwrite: bool = True,
) -> None:
    """Write ``document`` and its resources as the Session file ``path``.

    ``document`` is serialized as given; its ``resources`` name the resources
    in the order they are written, and ``resources`` gives the bytes of each.
    A resource's bytes must have the ``size`` its descriptor declares. The file
    is staged beside ``path`` and replaces it only when it is complete.
    """

    declared = document.get("resources")
    if not isinstance(declared, Mapping) or set(declared) != set(resources):
        raise ValidationError("Session resources do not match the resource bytes to write.", diagnostic=SESSION_RESOURCES_INVALID)
    document_bytes = json.dumps(document, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    with staged_session_file(path, overwrite=overwrite) as output:
        with zipfile.ZipFile(output, "w", allowZip64=True) as archive:
            _write_entry(archive, SESSION_DOCUMENT_ENTRY, io.BytesIO(document_bytes), len(document_bytes))
            for resource_id, descriptor in declared.items():
                size = descriptor.get("size") if isinstance(descriptor, Mapping) else None
                if not isinstance(size, int) or isinstance(size, bool) or size < 0:
                    raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID)
                with resources[resource_id].open_bytes() as stream:
                    _write_entry(archive, resource_entry_name(resource_id), stream, size)


def _write_entry(archive: zipfile.ZipFile, name: str, stream: BinaryIO, size: int) -> None:
    head = stream.read(len(GZIP_SIGNATURE))
    info = zipfile.ZipInfo(name, date_time=SESSION_ENTRY_DATE_TIME)
    info.create_system = _UNIX_SYSTEM
    info.external_attr = _FILE_ATTRIBUTES
    info.compress_type = SESSION_ENTRY_STORED if head == GZIP_SIGNATURE else SESSION_ENTRY_DEFLATED
    # The declared size decides whether the entry needs ZIP64 fields.
    info.file_size = size
    written = len(head)
    with zipfile.ZipFile.open(archive, info, "w") as entry:
        entry.write(head)
        while chunk := stream.read(_CHUNK_BYTES):
            entry.write(chunk)
            written += len(chunk)
    if written != size:
        raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID)


@contextmanager
def staged_session_file(path: str | Path, *, overwrite: bool) -> Iterator[BinaryIO]:
    """The one stream a Session file is written to.

    The file is staged in the same directory and becomes ``path`` only when
    the block ends without an error.
    """

    output_path = Path(path)
    temp_path: Path | None = None
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temp_fd, temp_name = tempfile.mkstemp(
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
        )
        temp_path = Path(temp_name)
        with os.fdopen(temp_fd, "w+b") as raw_file:
            yield raw_file
            raw_file.flush()
            os.fsync(raw_file.fileno())
        commit_staged_output_file(temp_path, output_path, overwrite=overwrite)
    except FileExistsError as exc:
        raise OutputExistsError((output_path,), subject="Session output") from exc
    except OSError as exc:
        raise ValidationError(f"Could not write session sidecar: {output_path}") from exc
    finally:
        try:
            if temp_path is not None:
                temp_path.unlink(missing_ok=True)
        except OSError:
            pass


@dataclass(frozen=True)
class _Entry:
    """Where an entry is, as the central directory gave it when the file was read.

    ``name`` is the whole name the directory holds: ``zipfile`` cuts a name at
    its first NUL (and on Windows turns ``\\`` into ``/``), the Web reader
    does not, so both compare the whole name. ``end`` is where the next entry
    or the directory starts; the entry's data must end before it.
    """

    name: str
    raw_name: bytes
    header_offset: int
    compress_type: int
    compress_size: int
    file_size: int
    crc: int
    end: int


@dataclass(frozen=True)
class SessionContainerEntry:
    """One resource of a read Session file; its bytes stay in the file until opened."""

    path: Path
    entry: _Entry

    @property
    def size(self) -> int:
        return self.entry.file_size

    def open_bytes(self) -> BinaryIO:
        """A stream of the resource's bytes, refused as damaged when its size or CRC-32 is wrong."""

        return _open_checked(self.path, self.entry)


@dataclass(frozen=True)
class SessionContainer:
    """A Session file read by :func:`read_session_container`."""

    path: Path
    # The parsed ``gbdraw-session.json``, not validated as a Session yet.
    document: object
    _resources: Mapping[str, _Entry]

    def resource(self, resource_id: str) -> SessionContainerEntry:
        """The resource ``resource_id`` of the document, found by its ID."""

        entry = self._resources.get(resource_id)
        if entry is None:
            raise SessionContainerError("damaged", self.path)
        return SessionContainerEntry(self.path, entry)

    def resource_ids(self) -> tuple[str, ...]:
        return tuple(self._resources)


def read_session_container(path: str | Path) -> SessionContainer:
    """Read the Session file ``path`` and check what can be checked without its inputs.

    Refused as damaged: an unreadable or cut central directory, or one that
    holds another number of entries than its end record declares, two entries
    of one whole name, an encrypted entry or one compressed otherwise than
    stored or deflated, a wrong size or CRC-32 of the document, a resource of
    the document whose entry is missing or has another size than it declares,
    and a resource whose descriptor embeds its bytes (``encoding`` or
    ``data``, a JSON Session's form: a Session file's bytes are its entries).
    A ZIP without ``gbdraw-session.json``, and a document that is not UTF-8
    JSON, are refused as such; so is a document whose ``version`` is an
    integer below the current one, since only the current version was ever
    written as a container. A resource's own
    bytes are checked when they are read (:meth:`SessionContainerEntry.open_bytes`).
    Entries the document does not name are not read.
    """

    session_path = Path(path)
    entries = _central_directory(session_path)
    document_entry = entries.get(SESSION_DOCUMENT_ENTRY)
    if document_entry is None:
        raise SessionContainerError("no-document", session_path)
    with _open_checked(session_path, document_entry) as stream:
        document_bytes = stream.read()
    try:
        text = document_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SessionContainerError("invalid-json", session_path) from exc
    try:
        document = parse_session_json(text, session_path)
    except SessionContainerError:
        raise
    except ValidationError as exc:
        # A duplicate key: the document is not valid JSON, as the Web reader says (OV-2692).
        raise SessionContainerError("invalid-json", session_path) from exc
    return SessionContainer(session_path, document, _declared_resources(session_path, document, entries))


def _central_directory(path: Path) -> dict[str, _Entry]:
    try:
        with open(path, "rb") as file:
            with zipfile.ZipFile(file) as archive:
                infos = archive.infolist()
                directory_start = archive.start_dir
            declared_count = _declared_entry_count(file)
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {path}", diagnostic=SESSION_FILE_UNREADABLE) from exc
    except _ZIP_ERRORS as exc:
        raise SessionContainerError("damaged", path) from exc
    if declared_count != len(infos):
        raise SessionContainerError("damaged", path)
    by_offset = sorted(infos, key=lambda info: info.header_offset)
    ends = {id(info): following.header_offset for info, following in zip(by_offset, by_offset[1:])}
    entries: dict[str, _Entry] = {}
    for info in infos:
        name = info.orig_filename
        if (name in entries or info.flag_bits & SESSION_ENTRY_ENCRYPTED_FLAGS
                or info.compress_type not in SESSION_ENTRY_METHODS):
            raise SessionContainerError("damaged", path)
        raw_name = name.encode("utf-8" if info.flag_bits & _UTF8_NAME_FLAG else "cp437")
        entries[name] = _Entry(
            name, raw_name, info.header_offset, info.compress_type, info.compress_size, info.file_size, info.CRC,
            ends.get(id(info), directory_start),
        )
    return entries


def _declared_entry_count(file: BinaryIO) -> int:
    """The number of entries the end record declares, or the ZIP64 end record
    when there is one, found where ``zipfile`` found them; ``zipfile`` reads
    the directory by its size and never compares the count."""

    size = file.seek(0, os.SEEK_END)
    tail_start = max(0, size - _END_RECORD_BYTES - _MAX_COMMENT_BYTES)
    file.seek(tail_start)
    tail = file.read()
    end_at = tail.rfind(_END_RECORD_SIGNATURE)
    (count,) = struct.unpack_from("<H", tail, end_at + 10)
    locator_at = tail_start + end_at - _ZIP64_LOCATOR_BYTES
    if locator_at >= 0:
        file.seek(locator_at)
        if file.read(len(_ZIP64_LOCATOR_SIGNATURE)) == _ZIP64_LOCATOR_SIGNATURE:
            file.seek(locator_at - _ZIP64_END_BYTES + 32)
            (count,) = struct.unpack("<Q", file.read(8))
    return int(count)


def _declared_resources(path: Path, document: object, entries: Mapping[str, _Entry]) -> dict[str, _Entry]:
    """The entries of the resources the document declares, each with its declared size.

    A document whose ``resources`` is not an object declares none here; the
    Session reader refuses it.
    """

    resources = document.get("resources") if isinstance(document, dict) else None
    if not isinstance(resources, dict):
        return {}
    declared: dict[str, _Entry] = {}
    for resource_id, descriptor in resources.items():
        entry = entries.get(resource_entry_name(resource_id))
        size = descriptor.get("size") if isinstance(descriptor, dict) else None
        if (entry is None or entry.file_size != size or isinstance(size, bool)
                or "encoding" in descriptor or "data" in descriptor):
            raise SessionContainerError("damaged", path)
        declared[resource_id] = entry
    return declared


def _open_checked(path: Path, entry: _Entry) -> BinaryIO:
    """Open one entry of ``path`` where the directory placed it, as the Web reader opens it."""

    try:
        file = open(path, "rb")
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {path}", diagnostic=SESSION_FILE_UNREADABLE) from exc
    try:
        file.seek(entry.header_offset)
        header = file.read(_LOCAL_HEADER_BYTES)
        if len(header) != _LOCAL_HEADER_BYTES or header[:4] != _LOCAL_HEADER_SIGNATURE:
            raise SessionContainerError("damaged", path)
        name_length, extra_length = struct.unpack_from("<HH", header, 26)
        start = entry.header_offset + _LOCAL_HEADER_BYTES + name_length + extra_length
        if file.read(name_length) != entry.raw_name or start + entry.compress_size > entry.end:
            raise SessionContainerError("damaged", path)
        file.seek(start)
    except SessionContainerError:
        file.close()
        raise
    except OSError as exc:
        file.close()
        raise ValidationError(f"Could not read session file: {path}", diagnostic=SESSION_FILE_UNREADABLE) from exc
    return io.BufferedReader(_CheckedEntry(path, file, entry), _CHUNK_BYTES)


class _CheckedEntry(io.RawIOBase):
    """An entry's content, refused as damaged unless it is exactly its size with
    a matching CRC-32, and a deflated entry's data is one whole deflate stream
    that fills its compressed size, as the Web's ``DecompressionStream`` reads it
    (``zipfile`` reads past neither an unfinished stream nor bytes after it).
    The checks run when the last byte of content is read.
    """

    def __init__(self, path: Path, file: BinaryIO, entry: _Entry) -> None:
        super().__init__()
        self._path = path
        self._file = file
        self._left = entry.file_size
        self._compressed_left = entry.compress_size
        self._crc = 0
        self._expected_crc = entry.crc
        # Python's zipfile and the Web read a deflated entry of no bytes as empty content.
        self._inflater = (
            zlib.decompressobj(-zlib.MAX_WBITS)
            if entry.compress_type == SESSION_ENTRY_DEFLATED and entry.compress_size else None
        )
        self._checked = False

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: WriteableBuffer) -> int:
        view = memoryview(buffer).cast("B")
        if not self._left or not len(view):
            self._finish()
            return 0
        chunk = self._read(min(len(view), self._left))
        if not chunk:
            raise SessionContainerError("damaged", self._path)
        self._left -= len(chunk)
        self._crc = zlib.crc32(chunk, self._crc)
        view[: len(chunk)] = chunk
        if not self._left:
            self._finish()
        return len(chunk)

    def _finish(self) -> None:
        if self._left or self._checked:
            return
        inflater = self._inflater
        if (self._read(1) or self._crc != self._expected_crc or self._compressed_left
                or (inflater is not None and (not inflater.eof or inflater.unused_data))):
            raise SessionContainerError("damaged", self._path)
        self._checked = True

    def _read(self, size: int) -> bytes:
        """Up to ``size`` bytes of content; empty at the end of the entry's data."""

        try:
            inflater = self._inflater
            if inflater is None:
                return self._compressed(size)
            while not inflater.eof:
                data = inflater.unconsumed_tail or self._compressed(_CHUNK_BYTES)
                if not data:
                    return b""
                content = inflater.decompress(data, size)
                if content:
                    return content
            return b""
        except OSError as exc:
            raise ValidationError(f"Could not read session file: {self._path}", diagnostic=SESSION_FILE_UNREADABLE) from exc
        except zlib.error as exc:
            raise SessionContainerError("damaged", self._path) from exc

    def _compressed(self, size: int) -> bytes:
        data = self._file.read(min(size, self._compressed_left))
        self._compressed_left -= len(data)
        return data

    def close(self) -> None:
        if not self.closed:
            self._file.close()
        super().close()


__all__ = [
    "BytesResourceSource",
    "FileResourceSource",
    "GZIP_SIGNATURE",
    "SESSION_CONTAINER_DAMAGED_MESSAGE",
    "SESSION_CONTAINER_NO_DOCUMENT_MESSAGE",
    "SESSION_CONTAINER_SIGNATURES",
    "SESSION_DOCUMENT_ENTRY",
    "SESSION_ENTRY_DATE_TIME",
    "SESSION_ENTRY_DEFLATED",
    "SESSION_ENTRY_ENCRYPTED_FLAGS",
    "SESSION_ENTRY_METHODS",
    "SESSION_ENTRY_STORED",
    "SESSION_FILE_UNREADABLE",
    "SESSION_JSON_INVALID",
    "SESSION_RESOURCES_INVALID",
    "SESSION_RESOURCE_ENTRY_PREFIX",
    "SessionContainer",
    "SessionContainerEntry",
    "SessionContainerError",
    "SessionResourceSource",
    "invalid_session_json_message",
    "is_session_container",
    "parse_session_json",
    "read_session_container",
    "reject_duplicate_json_keys",
    "resource_entry_name",
    "staged_session_file",
    "unsupported_session_version_message",
    "write_session_container",
]
