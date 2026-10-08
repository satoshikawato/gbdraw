"""The resources of one Session document, allocated by content.

A Session stores input bytes as resources keyed by ID; the canonical request
and the Web file bindings name them. :class:`SessionResourceTable` allocates
those IDs and file names for every Python writer: the request codec, the CLI
re-save of a Session, and the Web file bindings.
"""

from __future__ import annotations

import base64
import hashlib
import mimetypes
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Collection, Mapping

from gbdraw import session_io
from gbdraw.io.filenames import safe_embedded_filename, unique_filename
from gbdraw.session_request_codec import (
    _RESOURCE_ID_RE,
    CanonicalRequestEncodingError,
    CanonicalRequestResource,
)

# The keys under which a canonical request or a Web file binding names a
# resource, as the Web reads them (services/canonical-resource-references.js).
_RESOURCE_REFERENCE_FIELDS = frozenset({"resourceId", "gffResourceId", "fastaResourceId"})


def canonical_resource_ids(value: object) -> set[str]:
    """The resource IDs that a canonical request or Web file bindings name."""

    found: set[str] = set()
    pending = [value]
    while pending:
        item = pending.pop()
        if isinstance(item, list):
            pending.extend(item)
        elif isinstance(item, Mapping):
            for key, child in item.items():
                if key in _RESOURCE_REFERENCE_FIELDS and isinstance(child, str) and child.strip():
                    found.add(child)
                else:
                    pending.append(child)
    return found


@dataclass
class _HeldResource:
    """A table entry: a Session descriptor, or a request resource to serialize."""

    size: int
    name: str
    descriptor: Mapping[str, Any] | None = None
    resource: CanonicalRequestResource | None = None
    digest: str | None = None


class SessionResourceTable:
    """The resources of one Session, seeded with the resources it already has.

    Bytes equal to a held resource (equal size, then SHA-256) take that
    resource's ID and file name, so a re-saved Session keeps the IDs of its
    unchanged inputs. A new resource takes its preferred ID while it is free,
    else ``<id>-2``, ``<id>-3``, ... A Session materializes its resources side
    by side by sanitized name, so a new name used before (``a/X.fna`` and
    ``b/X.fna``) takes the next number, ``X.2.fna``, the ``--losat_output_dir``
    rule. A table without seed resources gives a request the IDs and names it
    gets alone.
    """

    def __init__(self, resources: Mapping[str, Mapping[str, Any]] | None = None) -> None:
        self._held: dict[str, _HeldResource] = {}
        self._names: set[str] = set()
        self._payloads: dict[tuple[Any, Any], str] = {}
        self._next_number = 1
        for resource_id, descriptor in (resources or {}).items():
            self._hold_descriptor(resource_id, descriptor)

    def request(self) -> RequestResources:
        """Allocate the resources of one canonical request."""

        return RequestResources(self)

    def bind(self, entry: Mapping[str, Any], *, preferred_id: str = "") -> str:
        """The resource of one Web file's bytes; a binding carries its own name.

        A held resource with equal bytes is shared. A new resource keeps
        ``preferred_id`` and the descriptor's name when both are free and safe;
        otherwise it is ``resource-NNNN``, named ``resource-NNNN-<name>``, as the
        Web writer names it (``services/session-resources.js``).
        """

        content: list[bytes] = []

        def read() -> bytes:
            if not content:
                content.append(session_io._embedded_resource_bytes(entry))
            return content[0]

        held = self._held.get(preferred_id)
        existing = held.descriptor if held is not None else None
        if entry.get("checksum") and entry is not existing and entry.get("checksum") != (existing or {}).get("checksum"):
            read()
        if existing is not None and all(
            entry.get(field) == existing.get(field) for field in ("encoding", "data", "size")
        ):
            return preferred_id
        size = entry.get("size")
        if entry.get("encoding") == "base64":
            shared = self._payloads.get((size, entry["data"]))
            if shared is not None:
                return shared
        digest: list[str] = []

        def entry_digest() -> str:
            if not digest:
                digest.append(hashlib.sha256(read()).hexdigest())
            return digest[0]

        shared = self._find_equal(size, entry_digest, preferred_id=preferred_id)
        if shared is not None:
            return shared
        data = read()
        safe_name = safe_embedded_filename(entry.get("name"), fallback="resource.dat")
        if (_RESOURCE_ID_RE.fullmatch(preferred_id) and preferred_id not in self._held
                and safe_name == entry.get("name") and safe_name not in self._names):
            resource_id, name = preferred_id, safe_name
        else:
            while True:
                resource_id = f"resource-{self._next_number:04d}"
                self._next_number += 1
                name = f"{resource_id}-{safe_name}"
                if resource_id not in self._held and name not in self._names:
                    break
        self._hold_descriptor(
            resource_id,
            {
                **entry, "kind": str(entry.get("kind") or "web-file"), "name": name,
                "size": len(data),
                "type": str(entry.get("type") or "application/octet-stream"),
                "encoding": "base64", "data": (
                    entry["data"] if entry.get("encoding") != session_io.DEPTH_FILE_ENCODING
                    else base64.b64encode(data).decode("ascii")
                ),
            },
            digest=digest[0] if digest else None,
        )
        return resource_id

    def descriptors(self) -> dict[str, Mapping[str, Any]]:
        """Every held resource as a Session resource descriptor, in table order."""

        return {
            resource_id: (
                held.descriptor if held.descriptor is not None
                else _resource_descriptor(_request_resource(held))
            )
            for resource_id, held in self._held.items()
        }

    def _find_equal(
        self,
        size: object,
        digest: Callable[[], str],
        *,
        preferred_id: str = "",
        name: str | None = None,
        exclude: Collection[str] = (),
    ) -> str | None:
        """A held resource with these bytes: ``preferred_id``, then ``name``, then the first."""

        candidates = [
            resource_id for resource_id, held in self._held.items()
            if held.size == size and resource_id not in exclude
        ]
        if not candidates:
            return None
        wanted = digest()
        candidates.sort(key=lambda resource_id: (
            resource_id != preferred_id, self._held[resource_id].name != name,
        ))
        return next(
            (resource_id for resource_id in candidates if self._digest(resource_id) == wanted),
            None,
        )

    def _digest(self, resource_id: str) -> str:
        held = self._held[resource_id]
        if held.digest is None:
            held.digest = hashlib.sha256(
                session_io._embedded_resource_bytes(held.descriptor) if held.descriptor is not None
                else _resource_bytes(_request_resource(held))
            ).hexdigest()
        return held.digest

    def _free_id(self, preferred_id: str) -> str:
        resource_id, ordinal = preferred_id, 1
        while resource_id in self._held:
            ordinal += 1
            resource_id = f"{preferred_id}-{ordinal}"
        return resource_id

    def _hold_descriptor(
        self,
        resource_id: str,
        descriptor: Mapping[str, Any],
        *,
        digest: str | None = None,
    ) -> None:
        name = safe_embedded_filename(descriptor.get("name"))
        self._held[resource_id] = _HeldResource(
            size=descriptor["size"], name=name, descriptor=descriptor, digest=digest,
        )
        self._names.add(name)
        if descriptor.get("encoding") == "base64":
            self._payloads.setdefault((descriptor["size"], descriptor["data"]), resource_id)

    def _hold_resource(self, resource: CanonicalRequestResource, *, size: int) -> None:
        name = safe_embedded_filename(resource.name)
        self._held[resource.resource_id] = _HeldResource(size=size, name=name, resource=resource)
        self._names.add(name)


class RequestResources:
    """The resources one canonical request adds to a :class:`SessionResourceTable`.

    Each input keeps a resource of its own, as when the request is encoded
    alone: the planner reads each source file as one genome
    (``losat_source_ids``), so two inputs never become one file. An input
    takes a held resource with equal bytes that no other input of the request
    took (its preferred ID first, then its file name), else a new resource.
    """

    def __init__(self, table: SessionResourceTable) -> None:
        self._table = table
        self._requested: set[str] = set()
        self._taken: set[str] = set()
        self._added: list[CanonicalRequestResource] = []

    def add_path(self, resource_id: str, *, kind: str, value: object) -> str:
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise CanonicalRequestEncodingError(
                f"Resource {resource_id!r} must identify a materialized file."
            )
        path = Path(str(value))
        if not path.is_file():
            raise CanonicalRequestEncodingError(
                f"Canonical request resource is not a file: {path}."
            )
        return self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=path.name,
                source_path=path,
            )
        )

    def add_bytes(
        self,
        resource_id: str,
        *,
        kind: str,
        name: str,
        content: bytes,
    ) -> str:
        return self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=name,
                content=content,
            )
        )

    def added(self) -> tuple[CanonicalRequestResource, ...]:
        """The new resources of this request, in the order it added them."""

        return tuple(self._added)

    def _add(self, resource: CanonicalRequestResource) -> str:
        if resource.resource_id in self._requested:
            raise CanonicalRequestEncodingError(
                f"Duplicate canonical resource ID: {resource.resource_id}."
            )
        self._requested.add(resource.resource_id)
        table = self._table
        name = safe_embedded_filename(resource.name)
        size = _resource_size(resource)
        resource_id = table._find_equal(
            size,
            lambda: hashlib.sha256(_resource_bytes(resource)).hexdigest(),
            preferred_id=resource.resource_id,
            name=name,
            exclude=self._taken,
        )
        if resource_id is None:
            resource_id = table._free_id(resource.resource_id)
            unique_name = unique_filename(resource.name, table._names, key=safe_embedded_filename)
            if (resource_id, unique_name) != (resource.resource_id, resource.name):
                resource = replace(resource, resource_id=resource_id, name=unique_name)
            table._hold_resource(resource, size=size)
            self._added.append(resource)
        self._taken.add(resource_id)
        return resource_id


def _request_resource(held: _HeldResource) -> CanonicalRequestResource:
    assert held.resource is not None
    return held.resource


def _resource_size(resource: CanonicalRequestResource) -> int:
    if resource.content is not None:
        return len(resource.content)
    assert resource.source_path is not None
    try:
        return resource.source_path.stat().st_size
    except OSError as exc:
        raise CanonicalRequestEncodingError(
            f"Could not read canonical resource: {resource.source_path}"
        ) from exc


def _resource_bytes(resource: CanonicalRequestResource) -> bytes:
    if resource.content is not None:
        return resource.content
    assert resource.source_path is not None
    try:
        return resource.source_path.read_bytes()
    except OSError as exc:
        raise CanonicalRequestEncodingError(
            f"Could not read canonical resource: {resource.source_path}"
        ) from exc


def _resource_descriptor(resource: CanonicalRequestResource) -> dict[str, Any]:
    content = _resource_bytes(resource)
    last_modified = 0
    if resource.source_path is not None:
        try:
            last_modified = int(resource.source_path.stat().st_mtime * 1000)
        except OSError as exc:
            raise CanonicalRequestEncodingError(
                f"Could not read canonical resource: {resource.source_path}"
            ) from exc
    media_type = mimetypes.guess_type(resource.name)[0] or "application/octet-stream"
    return {
        "kind": resource.kind,
        "name": safe_embedded_filename(resource.name),
        "type": media_type,
        "size": len(content),
        "lastModified": last_modified,
        "encoding": "base64",
        "data": base64.b64encode(content).decode("ascii"),
    }


__all__ = ["RequestResources", "SessionResourceTable", "canonical_resource_ids"]
