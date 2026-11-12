"""Write the Session files a test builds.

The current version is written as the Session file gbdraw writes (a
``.gbdraw`` container); an older version as the JSON it was saved in, plain
or gzip by the name, so a test can build an old Session to read.
"""

from __future__ import annotations

import base64
import gzip
import json
import zipfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from gbdraw.session_container import (
    SessionResourceSource,
    is_session_container,
    read_session_container,
    write_session_container,
)
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    compact_session_feature_catalog,
    expand_session_feature_catalog,
    read_source_bytes,
    split_resource_bytes,
    validate_session,
)


def write_session_file(
    path: str | Path,
    payload: Mapping[str, Any],
    sources: Mapping[str, SessionResourceSource] | None = None,
) -> None:
    """Validate ``payload`` and write it as ``path``.

    A current Session's resources may be embedded (base64, as in a JSON
    Session) or be descriptors whose bytes ``sources`` holds.
    """

    expanded = expand_session_feature_catalog(payload)
    if expanded.get("version") == CURRENT_SESSION_VERSION:
        document, split = split_resource_bytes(expanded, sources)
        validate_session(document, resource_sources=split)
        write_session_container(compact_session_feature_catalog(document), split, path)
        return
    if sources:
        expanded = embedded_session(expanded, sources)
    validate_session(expanded)
    text = json.dumps(compact_session_feature_catalog(expanded), ensure_ascii=False, separators=(",", ":"))
    data = text.encode("utf-8")
    Path(path).write_bytes(gzip.compress(data, mtime=0) if str(path).endswith(".gz") else data)


def write_unchecked_session_file(path: str | Path, payload: Mapping[str, Any]) -> None:
    """Write ``payload`` as ``path`` without validating it, so a test can give a reader a bad Session.

    A current Session becomes a Session file whose resources' bytes are the
    embedded ones; an older version is written as JSON.
    """

    if payload.get("version") == CURRENT_SESSION_VERSION:
        document, split = split_resource_bytes(payload)
        write_session_container(document, split, path)
        return
    Path(path).write_text(json.dumps(payload), encoding="utf-8")


def load_session_mapping(payload: Mapping[str, Any]) -> Any:
    """Load a test payload as gbdraw loads it from a file.

    A Session 51 whose resources embed their bytes (as tests and the Web's
    legacy upgrade build it) is read as its .gbdraw file is: its bytes move
    beside the document. An older JSON Session is read as it is.
    """

    from gbdraw.session import load_session_document
    from gbdraw.session_io import SessionResources

    if payload.get("version") != CURRENT_SESSION_VERSION:
        return load_session_document(payload)
    document, sources = split_resource_bytes(payload)
    return load_session_document(document, resources=SessionResources(document.get("resources", {}), sources))


def embedded_session(
    document: Mapping[str, Any],
    sources: Mapping[str, SessionResourceSource],
) -> dict[str, Any]:
    """``document`` with each resource's bytes embedded as base64, as a JSON Session holds them.

    Tests edit such a payload and write it back with :func:`write_session_file`;
    gbdraw itself never embeds bytes in a current Session.
    """

    resources = document.get("resources")
    if not isinstance(resources, Mapping):
        return dict(document)
    embedded = {}
    for resource_id, descriptor in resources.items():
        with sources[resource_id].open_bytes() as stream:
            data = stream.read()
        embedded[resource_id] = {**descriptor, "encoding": "base64", "data": base64.b64encode(data).decode("ascii")}
    return {**document, "resources": embedded}


def _is_embedded(payload: Mapping[str, Any]) -> bool:
    resources = payload.get("resources")
    return isinstance(resources, Mapping) and any(
        isinstance(entry, Mapping) and ("encoding" in entry or "data" in entry) for entry in resources.values()
    )


def json_document(payload: Mapping[str, Any]) -> dict[str, Any]:
    """``payload`` as a loaded document holds it: resource descriptors without their bytes."""

    return split_resource_bytes(payload)[0]


def read_session_payload(path: str | Path) -> dict[str, Any]:
    """Any Session file (a container, or an older JSON Session) as a JSON Session payload, not validated."""

    from gbdraw.session_io import _read_session_json

    session_path = Path(path)
    head = session_path.read_bytes()[:4]
    if is_session_container(head):
        container = read_session_container(session_path)
        document = container.document
        assert isinstance(document, dict)
        return embedded_session(
            document, {resource_id: container.resource(resource_id) for resource_id in container.resource_ids()}
        )
    payload, _sources = _read_session_json(session_path)
    assert isinstance(payload, dict)
    return payload


def read_session_text(path: str | Path) -> str:
    """The JSON text of a Session file's document (``gbdraw-session.json``)."""

    from gbdraw.session_container import SESSION_DOCUMENT_ENTRY

    with zipfile.ZipFile(path) as archive:
        return archive.read(SESSION_DOCUMENT_ENTRY).decode("utf-8")


def build_session_payload(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """``build_session_json`` as one payload with its resources' bytes embedded.

    A context's ``source_session`` may be embedded too; its bytes become its
    ``source_resources``.
    """

    import dataclasses

    from gbdraw.session_io import SessionBuildContext, build_session_json

    context = args[0] if args else kwargs.get("context")
    source = getattr(context, "source_session", None)
    if (
        isinstance(context, SessionBuildContext)
        and isinstance(source, Mapping)
        and context.source_resources is None
        and _is_embedded(source)
    ):
        document, split = split_resource_bytes(source)
        context = dataclasses.replace(context, source_session=document, source_resources=split)
        args = (context, *args[1:]) if args else args
        if not args:
            kwargs["context"] = context
    document, sources = build_session_json(*args, **kwargs)
    return embedded_session(document, sources)


def legacy_file_tree(value: Any) -> Any:
    """``value`` with every CLI resource file (``SessionResourceFile``) as the
    embedded entry an older JSON Session stores: its descriptor's metadata and
    its bytes as base64."""

    if isinstance(value, Mapping) and {"resourceId", "descriptor", "source"} <= value.keys():
        descriptor = value["descriptor"]
        data = read_source_bytes(value["source"], descriptor["size"])
        return {
            "name": value["name"],
            "type": value["type"],
            "size": len(data),
            "lastModified": value["lastModified"],
            "data": base64.b64encode(data).decode("ascii"),
        }
    if isinstance(value, Mapping):
        return {key: legacy_file_tree(item) for key, item in value.items()}
    if isinstance(value, list):
        return [legacy_file_tree(item) for item in value]
    return value
