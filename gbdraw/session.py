"""Public canonical session document lifecycle.

Canonical sessions carry a CLI-independent ``renderRequest`` and a mapping of
embedded resources. Resource paths exposed by :class:`MaterializedSession` are
temporary and are valid only while the materialization context is active.
"""

from __future__ import annotations

import base64
import copy
import json
import mimetypes
import re
import tempfile
from contextlib import AbstractContextManager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping

from gbdraw.exceptions import ValidationError
from gbdraw.render.output_paths import preflight_output_paths
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    CANONICAL_SESSION_MIN_VERSION,
    DEPTH_FILE_ENCODING,
    SESSION_FORMAT,
    _read_session_text,
    _attach_current_web_file_bindings,
    _reject_duplicate_json_keys,
    expand_session_feature_catalog,
    materialize_embedded_file,
    normalize_current_session_artifacts,
    safe_embedded_filename,
    validate_session,
    write_session_json,
)

if TYPE_CHECKING:
    from gbdraw.api.request_render import (
        CircularBatchRenderResult,
        RequestRenderResult,
    )
    from gbdraw.api.requests import DiagramRequest
    from gbdraw.session_request_codec import CanonicalRequestResource


_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
_RESOURCE_REQUIRED_FIELDS = frozenset(
    {"kind", "name", "type", "size", "encoding", "data"}
)
_RESOURCE_OPTIONAL_FIELDS = frozenset({"lastModified"})
class SessionError(ValidationError):
    """Base error for the public session bridge."""


class SessionFormatError(SessionError):
    """Raised when a session document or canonical envelope is malformed."""


class SessionVersionError(SessionError):
    """Raised when a session version has no requested conversion path."""


class SessionResourceError(SessionError):
    """Raised when an embedded canonical resource cannot be materialized."""


class SessionConversionError(SessionError):
    """Raised when a canonical payload cannot be converted to a typed request."""


class SessionRenderError(SessionError):
    """Raised when canonical session rendering fails."""


@dataclass(frozen=True)
class SessionDocument:
    """Validated session envelope detached from caller-owned mutable mappings."""

    _data: Mapping[str, Any]
    source_path: Path | None = None

    def __post_init__(self) -> None:
        cloned = expand_session_feature_catalog(copy.deepcopy(dict(self._data)))
        _validate_document(cloned)
        object.__setattr__(self, "_data", cloned)
        if self.source_path is not None:
            object.__setattr__(self, "source_path", Path(self.source_path))

    @property
    def version(self) -> int:
        """Session envelope version."""

        return int(self._data["version"])

    @property
    def mode(self) -> Literal["circular", "linear"] | None:
        """Canonical mode, when the document contains a canonical request."""

        request = self._data.get("renderRequest")
        if isinstance(request, Mapping) and request.get("mode") in {
            "circular",
            "linear",
        }:
            return request["mode"]  # type: ignore[return-value]
        return None

    @property
    def has_canonical_request(self) -> bool:
        """Whether this document can enter the canonical typed bridge."""

        return self.version >= CANONICAL_SESSION_MIN_VERSION and isinstance(
            self._data.get("renderRequest"), Mapping
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a detached JSON-compatible copy of the document."""

        return copy.deepcopy(dict(self._data))


@dataclass
class _MaterializationLifetime:
    active: bool = True


@dataclass(frozen=True)
class MaterializedSession:
    """Materialized canonical resources with context-bounded path lifetime.

    ``resource_paths`` and any paths stored in a decoded request become invalid
    as soon as the owning ``materialize_session(...)`` context exits.
    """

    document: SessionDocument
    temp_directory: Path
    resource_paths: Mapping[str, Path]
    output_directory: Path
    _lifetime: _MaterializationLifetime

    @property
    def active(self) -> bool:
        """Whether temporary resources are still owned by an active context."""

        return self._lifetime.active


class _SessionMaterializationContext(AbstractContextManager[MaterializedSession]):
    def __init__(
        self,
        document: SessionDocument,
        *,
        output_directory: str | Path,
        temporary_directory: str | Path | None,
    ) -> None:
        raw_output_directory = str(output_directory).strip()
        if not raw_output_directory:
            raise SessionResourceError("A replay output directory is required.")
        self._document = document
        self._output_directory = Path(raw_output_directory)
        self._temporary_directory = (
            Path(temporary_directory) if temporary_directory is not None else None
        )
        self._owner: tempfile.TemporaryDirectory[str] | None = None
        self._materialized: MaterializedSession | None = None

    def __enter__(self) -> MaterializedSession:
        if self._materialized is not None:
            raise SessionResourceError("A session materialization context is single-use.")
        try:
            if self._temporary_directory is not None:
                self._temporary_directory.mkdir(parents=True, exist_ok=True)
            self._owner = tempfile.TemporaryDirectory(
                prefix="gbdraw-session-v31-",
                dir=(
                    str(self._temporary_directory)
                    if self._temporary_directory is not None
                    else None
                ),
            )
            temp_directory = Path(self._owner.name)
            resource_paths = _materialize_resources(
                self._document,
                temp_directory=temp_directory,
            )
            self._materialized = MaterializedSession(
                document=self._document,
                temp_directory=temp_directory,
                resource_paths=resource_paths,
                output_directory=self._output_directory,
                _lifetime=_MaterializationLifetime(),
            )
            return self._materialized
        except SessionError:
            self._cleanup_after_failed_enter()
            raise
        except (OSError, ValidationError) as exc:
            self._cleanup_after_failed_enter()
            raise SessionResourceError(
                f"Canonical session resources could not be materialized: {exc}"
            ) from exc

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if self._materialized is not None:
            self._materialized._lifetime.active = False
        if self._owner is None:
            return False
        try:
            self._owner.cleanup()
        except OSError as cleanup_error:
            if exc_value is not None:
                if hasattr(exc_value, "add_note"):
                    exc_value.add_note(
                        f"Temporary session cleanup also failed: {cleanup_error}"
                    )
                return False
            raise SessionResourceError(
                "Temporary canonical session resources could not be cleaned up."
            ) from cleanup_error
        return False

    def _cleanup_after_failed_enter(self) -> None:
        if self._owner is None:
            return
        try:
            self._owner.cleanup()
        except OSError as cleanup_error:
            raise SessionResourceError(
                "Session materialization failed and temporary resources could not be cleaned up."
            ) from cleanup_error


def load_session_document(
    source: str | Path | Mapping[str, Any] | SessionDocument,
) -> SessionDocument:
    """Load and validate a session document without materializing resources."""

    if isinstance(source, SessionDocument):
        return source
    if isinstance(source, Mapping):
        try:
            return SessionDocument(source)
        except SessionError:
            raise
        except ValidationError as exc:
            raise _classify_validation_error(exc) from exc

    path = Path(source)
    try:
        payload = json.loads(
            _read_session_text(path),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except json.JSONDecodeError as exc:
        raise SessionFormatError(f"Not a valid JSON session file: {path}") from exc
    except OSError as exc:
        raise SessionFormatError(f"Could not read session file: {path}") from exc
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    if not isinstance(payload, Mapping):
        raise SessionFormatError("Session JSON must be an object.")
    try:
        return SessionDocument(payload, source_path=path)
    except SessionError:
        raise
    except ValidationError as exc:
        raise _classify_validation_error(exc) from exc


def materialize_session(
    document: SessionDocument | Mapping[str, Any] | str | Path,
    *,
    output_directory: str | Path,
    temporary_directory: str | Path | None = None,
) -> AbstractContextManager[MaterializedSession]:
    """Create the sole owner of temporary canonical session resources."""

    return _SessionMaterializationContext(
        load_session_document(document),
        output_directory=output_directory,
        temporary_directory=temporary_directory,
    )


def session_to_request(materialized: MaterializedSession) -> DiagramRequest:
    """Decode a canonical request within its resource lifetime."""

    if not isinstance(materialized, MaterializedSession):
        raise SessionConversionError("A MaterializedSession is required.")
    if not materialized.active:
        raise SessionResourceError(
            "Materialized session resources are no longer active; decode inside the context."
        )
    document = materialized.document
    if not document.has_canonical_request:
        raise SessionVersionError(
            "Sessions version 27 through 30 support internal CLI replay only and do not "
            "have a public typed-request conversion."
        )
    payload = document._data.get("renderRequest")
    assert isinstance(payload, Mapping)
    from gbdraw.session_request_codec import (
        CanonicalRequestCodecError,
        decode_canonical_request,
    )
    from gbdraw.api.session_compat import canonical_payload_for_session_decode

    try:
        return decode_canonical_request(
            canonical_payload_for_session_decode(document.version, payload),
            resource_paths=materialized.resource_paths,
            output_directory=materialized.output_directory,
        )
    except CanonicalRequestCodecError as exc:
        raise SessionConversionError(str(exc)) from exc


def render_session(
    materialized: MaterializedSession,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Decode and render a canonical session while its resources are active."""

    from gbdraw.api.session_compat import render_session_compatible_request

    try:
        return render_session_compatible_request(
            session_to_request(materialized),
            materialized.document.to_dict(),
        )
    except SessionError:
        raise
    except Exception as exc:
        raise SessionRenderError(f"Canonical session rendering failed: {exc}") from exc


def build_session_document(
    request: DiagramRequest,
    *,
    title: str | None = None,
    created_at: datetime | None = None,
    adjunct: Mapping[str, Any] | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
) -> SessionDocument:
    """Build a current-version document from one typed request.

    ``adjunct`` may contain Web/editor artifacts such as ``ui`` or ``results``;
    it cannot replace canonical envelope fields.
    """

    from gbdraw.session_request_codec import (
        CanonicalRequestCodecError,
        encode_canonical_request,
    )
    from gbdraw.api.request_render import resolve_request

    try:
        encoded = encode_canonical_request(resolve_request(request))
        resources = {
            resource.resource_id: _serialize_canonical_resource(resource)
            for resource in encoded.resources
        }
    except (CanonicalRequestCodecError, ValidationError) as exc:
        raise SessionConversionError(str(exc)) from exc

    adjunct_data = copy.deepcopy(dict(adjunct or {}))
    reserved = {"format", "version", "createdAt", "renderRequest", "resources"}
    conflicting = reserved & set(adjunct_data)
    if conflicting:
        raise SessionFormatError(
            "Session adjunct cannot replace canonical field(s): "
            + ", ".join(sorted(conflicting))
            + "."
        )
    timestamp = created_at or datetime.now(timezone.utc)
    data = {
        "format": SESSION_FORMAT,
        "version": CURRENT_SESSION_VERSION,
        "createdAt": timestamp.isoformat(),
        "renderRequest": encoded.payload,
        "resources": resources,
        "results": [],
        "editorState": {"featureCatalog": None},
    }
    data.update(adjunct_data)
    if web_file_inventory is not None:
        _attach_current_web_file_bindings(data, web_file_inventory)
    editor_state = data.get("editorState")
    if isinstance(editor_state, Mapping):
        normalized_editor_state = copy.deepcopy(dict(editor_state))
        normalized_editor_state.setdefault("featureCatalog", None)
        data["editorState"] = normalized_editor_state
    if title is not None:
        data["title"] = str(title)
    normalize_current_session_artifacts(data)
    return SessionDocument(data)


def save_session_document(
    path: str | Path,
    request: DiagramRequest,
    *,
    title: str | None = None,
    created_at: datetime | None = None,
    adjunct: Mapping[str, Any] | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> SessionDocument:
    """Build and write a current canonical session through a staged commit."""

    document = build_session_document(
        request,
        title=title,
        created_at=created_at,
        adjunct=adjunct,
        web_file_inventory=web_file_inventory,
    )
    try:
        output_path = Path(path)
        preflight_output_paths((output_path,), overwrite=True)
        if (output_path.exists() or output_path.is_symlink()) and not overwrite:
            raise ValidationError(
                f"Session output already exists: {output_path}. "
                "Pass overwrite=True to replace it."
            )
        write_session_json(path, document._data, overwrite=overwrite)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    return document


def with_request_output(
    request: DiagramRequest,
    *,
    output_prefix: str | None = None,
    output_directory: str | Path | None = None,
    formats: str | tuple[str, ...] | None = None,
    overwrite: bool | None = None,
) -> DiagramRequest:
    """Return a request with caller-owned replay output overrides."""

    from gbdraw.api.requests import CircularBatchRequest, RenderOutputRequest

    if isinstance(request, CircularBatchRequest):
        item_count = len(request.outputs)
        if output_prefix is None:
            prefixes = tuple(output.output_prefix for output in request.outputs)
        elif item_count == 1:
            prefixes = (output_prefix,)
        else:
            prefixes = tuple(
                f"{output_prefix}_{index}"
                for index in range(1, item_count + 1)
            )
        outputs = tuple(
            RenderOutputRequest(
                output_prefix=prefix,
                output_directory=(
                    output_directory
                    if output_directory is not None
                    else current.output_directory
                ),
                formats=formats if formats is not None else current.formats,
                overwrite=current.overwrite if overwrite is None else overwrite,
                interactive_metadata_policy=current.interactive_metadata_policy,
            )
            for prefix, current in zip(
                prefixes,
                request.outputs,
                strict=True,
            )
        )
        return replace(request, outputs=outputs)

    current = request.output
    updated = RenderOutputRequest(
        output_prefix=output_prefix or current.output_prefix,
        output_directory=(
            output_directory
            if output_directory is not None
            else current.output_directory
        ),
        formats=formats if formats is not None else current.formats,
        overwrite=current.overwrite if overwrite is None else overwrite,
        interactive_metadata_policy=current.interactive_metadata_policy,
    )
    return replace(request, output=updated)


def _validate_document(data: Mapping[str, Any]) -> None:
    try:
        validate_session(data)
    except ValidationError as exc:
        raise _classify_validation_error(exc) from exc
    if int(data.get("version", 0)) < CANONICAL_SESSION_MIN_VERSION:
        return
    resources = data.get("resources")
    assert isinstance(resources, Mapping)
    sanitized_names: set[str] = set()
    for resource_id, entry in resources.items():
        if not isinstance(resource_id, str) or not _RESOURCE_ID_RE.fullmatch(resource_id):
            raise SessionResourceError(
                f"Invalid canonical resource ID: {resource_id!r}."
            )
        if not isinstance(entry, Mapping):
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} must be an object."
            )
        fields = set(entry)
        missing = _RESOURCE_REQUIRED_FIELDS - fields
        unknown = fields - _RESOURCE_REQUIRED_FIELDS - _RESOURCE_OPTIONAL_FIELDS
        if missing:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} is missing field(s): "
                + ", ".join(sorted(missing))
                + "."
            )
        if unknown:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has unknown field(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        if not isinstance(entry.get("kind"), str) or not str(entry["kind"]).strip():
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires a kind."
            )
        name = safe_embedded_filename(entry.get("name"), fallback="")
        if not name:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires a safe basename."
            )
        if name in sanitized_names:
            raise SessionResourceError(
                f"Duplicate canonical resource filename after sanitization: {name!r}."
            )
        sanitized_names.add(name)
        encoding = entry.get("encoding")
        if encoding not in {"base64", DEPTH_FILE_ENCODING}:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has unsupported encoding {encoding!r}."
            )
        if encoding == "base64" and not isinstance(entry.get("data"), str):
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires base64 string data."
            )
        if encoding == DEPTH_FILE_ENCODING and not isinstance(entry.get("data"), Mapping):
            raise SessionResourceError(
                f"Canonical depth resource {resource_id!r} requires an object payload."
            )
        size = entry.get("size")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has invalid size metadata."
            )


def _materialize_resources(
    document: SessionDocument,
    *,
    temp_directory: Path,
) -> dict[str, Path]:
    if document.version < CANONICAL_SESSION_MIN_VERSION:
        return {}
    raw_resources = document._data.get("resources")
    assert isinstance(raw_resources, Mapping)
    result: dict[str, Path] = {}
    for resource_id, entry in raw_resources.items():
        assert isinstance(resource_id, str)
        assert isinstance(entry, Mapping)
        try:
            result[resource_id] = materialize_embedded_file(
                entry,
                temp_dir=temp_directory,
                role=resource_id,
                prefix_role=False,
            )
        except ValidationError as exc:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} could not be materialized: {exc}"
            ) from exc
    return result


def _serialize_canonical_resource(
    resource: CanonicalRequestResource,
) -> dict[str, Any]:
    if resource.content is not None:
        content = resource.content
        last_modified = 0
    else:
        assert resource.source_path is not None
        try:
            content = resource.source_path.read_bytes()
            last_modified = int(resource.source_path.stat().st_mtime * 1000)
        except OSError as exc:
            raise SessionResourceError(
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


def _classify_validation_error(exc: ValidationError) -> SessionError:
    message = str(exc)
    if "version" in message.lower():
        return SessionVersionError(message)
    return SessionFormatError(message)


__all__ = [
    "MaterializedSession",
    "SessionConversionError",
    "SessionDocument",
    "SessionError",
    "SessionFormatError",
    "SessionRenderError",
    "SessionResourceError",
    "SessionVersionError",
    "build_session_document",
    "load_session_document",
    "materialize_session",
    "render_session",
    "save_session_document",
    "session_to_request",
    "with_request_output",
]
