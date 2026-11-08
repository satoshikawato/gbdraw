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
import shutil
import struct
import tempfile
import zipfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO, TYPE_CHECKING, BinaryIO, Literal, Protocol

from .exceptions import ValidationError
from .render.output_paths import OutputExistsError, commit_staged_output_file

if TYPE_CHECKING:
    from _typeshed import WriteableBuffer

# The entry that holds the Session document.
SESSION_DOCUMENT_ENTRY = "gbdraw-session.json"
# A resource's entry is this prefix and the resource ID.
SESSION_RESOURCE_ENTRY_PREFIX = "resources/"
# The compression methods a Session entry may use: stored and deflated.
SESSION_ENTRY_METHODS = (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED)
# The first bytes of a gzip stream: an input that starts with them is stored.
GZIP_SIGNATURE = b"\x1f\x8b"
# The first bytes of a ZIP file: a local entry, or the end record of an empty ZIP.
SESSION_CONTAINER_SIGNATURES = (b"PK\x03\x04", b"PK\x05\x06")

SESSION_CONTAINER_DAMAGED_MESSAGE = "This file is not a gbdraw Session: it is incomplete or damaged."
SESSION_CONTAINER_NO_DOCUMENT_MESSAGE = "This file is not a gbdraw Session: it has no gbdraw-session.json."

# Every entry gets the same time, system and permissions, so the bytes depend
# on the content and the zlib build only.
_ENTRY_DATE_TIME = (1980, 1, 1, 0, 0, 0)
_UNIX_SYSTEM = 3
_FILE_ATTRIBUTES = 0o100644 << 16
# General purpose flags: bit 0 encrypted, bit 6 strong encryption.
_ENCRYPTED_FLAGS = 0x0041
_CHUNK_BYTES = 1 << 20
# What zipfile raises for a file it cannot read as a ZIP, besides OSError.
_ZIP_ERRORS = (
    zipfile.BadZipFile, zipfile.LargeZipFile, EOFError, ValueError, NotImplementedError, struct.error,
)

SessionContainerErrorCode = Literal["damaged", "no-document", "invalid-json"]


class SessionContainerError(ValidationError):
    """A Session file that cannot be read: a damaged container, one without a
    document, or a document that is not JSON (container or JSON Session).

    ``code`` names the refusal; the vectors give the text of each per surface.
    """

    def __init__(self, code: SessionContainerErrorCode, path: Path) -> None:
        super().__init__(
            SESSION_CONTAINER_NO_DOCUMENT_MESSAGE if code == "no-document"
            else invalid_session_json_message(path) if code == "invalid-json"
            else SESSION_CONTAINER_DAMAGED_MESSAGE
        )
        self.code: SessionContainerErrorCode = code


def invalid_session_json_message(path: str | Path) -> str:
    """The refusal of a Session file whose JSON is not readable."""

    return f"Not a valid JSON session file: {path}"


def reject_duplicate_json_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON object keys before the decoder can discard them."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError(f"Session JSON contains a duplicate object key: {key!r}.")
        result[key] = value
    return result


def parse_session_json(text: str, path: str | Path) -> object:
    """The JSON value of a Session document; unreadable JSON is refused as such."""

    try:
        return json.loads(text, object_pairs_hook=reject_duplicate_json_keys)
    except json.JSONDecodeError as exc:
        raise SessionContainerError("invalid-json", Path(path)) from exc


class SessionResourceSource(Protocol):
    """The bytes of one resource, read once from the start when the file is written."""

    def open(self) -> BinaryIO: ...


@dataclass(frozen=True)
class FileResourceSource:
    """A resource whose bytes are a local file, as they are on disk."""

    path: Path

    def open(self) -> BinaryIO:
        return open(self.path, "rb")


@dataclass(frozen=True)
class BytesResourceSource:
    """A resource whose bytes are in memory."""

    data: bytes

    def open(self) -> BinaryIO:
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
        raise ValidationError("Session resources do not match the resource bytes to write.")
    document_bytes = json.dumps(document, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    with staged_session_file(path, overwrite=overwrite) as output:
        with zipfile.ZipFile(output, "w", allowZip64=True) as archive:
            _write_entry(archive, SESSION_DOCUMENT_ENTRY, io.BytesIO(document_bytes), len(document_bytes))
            for resource_id, descriptor in declared.items():
                size = descriptor.get("size") if isinstance(descriptor, Mapping) else None
                if not isinstance(size, int) or isinstance(size, bool) or size < 0:
                    raise ValidationError("Invalid embedded resource bytes or byte size.")
                with resources[resource_id].open() as stream:
                    _write_entry(archive, resource_entry_name(resource_id), stream, size)


def _write_entry(archive: zipfile.ZipFile, name: str, stream: BinaryIO, size: int) -> None:
    head = stream.read(len(GZIP_SIGNATURE))
    info = zipfile.ZipInfo(name, date_time=_ENTRY_DATE_TIME)
    info.create_system = _UNIX_SYSTEM
    info.external_attr = _FILE_ATTRIBUTES
    info.compress_type = zipfile.ZIP_STORED if head == GZIP_SIGNATURE else zipfile.ZIP_DEFLATED
    # The declared size decides whether the entry needs ZIP64 fields.
    info.file_size = size
    written = len(head)
    with archive.open(info, "w") as entry:
        entry.write(head)
        while chunk := stream.read(_CHUNK_BYTES):
            entry.write(chunk)
            written += len(chunk)
    if written != size:
        raise ValidationError("Invalid embedded resource bytes or byte size.")


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
    """Where an entry is, as the central directory gave it when the file was read."""

    name: str
    header_offset: int
    compress_type: int
    compress_size: int
    file_size: int
    crc: int


@dataclass(frozen=True)
class SessionContainerResource:
    """One resource of a read Session file; its bytes stay in the file until opened."""

    path: Path
    entry: _Entry

    @property
    def size(self) -> int:
        return self.entry.file_size

    def open(self) -> BinaryIO:
        """A stream of the resource's bytes, refused as damaged when its size or CRC-32 is wrong."""

        return _open_checked(self.path, self.entry)

    def copy_to(self, output: BinaryIO) -> None:
        with self.open() as stream:
            shutil.copyfileobj(stream, output, _CHUNK_BYTES)


@dataclass(frozen=True)
class SessionContainer:
    """A Session file read by :func:`read_session_container`."""

    path: Path
    # The parsed ``gbdraw-session.json``, not validated as a Session yet.
    document: object
    _resources: Mapping[str, _Entry]

    def resource(self, resource_id: str) -> SessionContainerResource:
        """The resource ``resource_id`` of the document, found by its ID."""

        entry = self._resources.get(resource_id)
        if entry is None:
            raise SessionContainerError("damaged", self.path)
        return SessionContainerResource(self.path, entry)

    def resource_ids(self) -> tuple[str, ...]:
        return tuple(self._resources)


def read_session_container(path: str | Path) -> SessionContainer:
    """Read the Session file ``path`` and check what can be checked without its inputs.

    Refused as damaged: an unreadable or cut central directory, two entries of
    one name, an encrypted entry or one compressed otherwise than stored or
    deflated, a wrong size or CRC-32 of the document, and a resource of the
    document whose entry is missing or has another size than it declares.
    A ZIP without ``gbdraw-session.json``, and a document that is not UTF-8
    JSON, are refused as such. A resource's own
    bytes are checked when they are read (:meth:`SessionContainerResource.open`).
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
    document = parse_session_json(text, session_path)
    return SessionContainer(session_path, document, _declared_resources(session_path, document, entries))


def _central_directory(path: Path) -> dict[str, _Entry]:
    try:
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {path}") from exc
    except _ZIP_ERRORS as exc:
        raise SessionContainerError("damaged", path) from exc
    entries: dict[str, _Entry] = {}
    for info in infos:
        if (info.filename in entries or info.flag_bits & _ENCRYPTED_FLAGS
                or info.compress_type not in SESSION_ENTRY_METHODS):
            raise SessionContainerError("damaged", path)
        entries[info.filename] = _Entry(
            info.filename, info.header_offset, info.compress_type, info.compress_size, info.file_size, info.CRC,
        )
    return entries


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
        if entry is None or entry.file_size != size or isinstance(size, bool):
            raise SessionContainerError("damaged", path)
        declared[resource_id] = entry
    return declared


def _open_checked(path: Path, entry: _Entry) -> BinaryIO:
    """Open one entry of ``path`` as it was when the file was read."""

    try:
        archive = zipfile.ZipFile(path)
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {path}") from exc
    except _ZIP_ERRORS as exc:
        raise SessionContainerError("damaged", path) from exc
    try:
        info = archive.NameToInfo.get(entry.name)
        if info is None or (
            info.header_offset, info.compress_type, info.compress_size, info.file_size, info.CRC,
        ) != (entry.header_offset, entry.compress_type, entry.compress_size, entry.file_size, entry.crc):
            raise SessionContainerError("damaged", path)
        stream = archive.open(info)
    except SessionContainerError:
        archive.close()
        raise
    except OSError as exc:
        archive.close()
        raise ValidationError(f"Could not read session file: {path}") from exc
    except (*_ZIP_ERRORS, RuntimeError) as exc:
        archive.close()
        raise SessionContainerError("damaged", path) from exc
    return io.BufferedReader(_CheckedEntry(path, archive, stream, entry.file_size), _CHUNK_BYTES)


class _CheckedEntry(io.RawIOBase):
    """An entry's bytes: refused as damaged unless exactly its size with a matching CRC-32.

    ``zipfile`` checks the CRC-32 when the entry ends; this counts the bytes,
    which it does not check for a stored entry.
    """

    def __init__(self, path: Path, archive: zipfile.ZipFile, stream: IO[bytes], size: int) -> None:
        super().__init__()
        self._path = path
        self._archive = archive
        self._stream = stream
        self._left = size

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: WriteableBuffer) -> int:
        view = memoryview(buffer).cast("B")
        try:
            chunk = self._stream.read(min(len(view), self._left + 1))
        except OSError as exc:
            raise ValidationError(f"Could not read session file: {self._path}") from exc
        except (*_ZIP_ERRORS, RuntimeError) as exc:
            raise SessionContainerError("damaged", self._path) from exc
        if len(chunk) > self._left or (not chunk and self._left):
            raise SessionContainerError("damaged", self._path)
        self._left -= len(chunk)
        view[: len(chunk)] = chunk
        return len(chunk)

    def close(self) -> None:
        if not self.closed:
            self._stream.close()
            self._archive.close()
        super().close()


__all__ = [
    "BytesResourceSource",
    "FileResourceSource",
    "GZIP_SIGNATURE",
    "SESSION_CONTAINER_DAMAGED_MESSAGE",
    "SESSION_CONTAINER_NO_DOCUMENT_MESSAGE",
    "SESSION_CONTAINER_SIGNATURES",
    "SESSION_DOCUMENT_ENTRY",
    "SESSION_ENTRY_METHODS",
    "SESSION_RESOURCE_ENTRY_PREFIX",
    "SessionContainer",
    "SessionContainerError",
    "SessionContainerResource",
    "SessionResourceSource",
    "invalid_session_json_message",
    "is_session_container",
    "parse_session_json",
    "reject_duplicate_json_keys",
    "read_session_container",
    "resource_entry_name",
    "staged_session_file",
    "write_session_container",
]
