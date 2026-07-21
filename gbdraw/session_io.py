#!/usr/bin/env python
# coding: utf-8

"""GUI session JSON loading, validation, materialization, and sidecar building."""

from __future__ import annotations

import base64
import binascii
import copy
import csv
import gzip
import io
import json
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePath, PureWindowsPath
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence

from .definition_line_styles import DEFINITION_LINE_KINDS, parse_definition_line_style_overrides
from .exceptions import ValidationError
from .io.regions import RegionSpec, parse_region_specs
from .render.formats import normalize_format_token

if TYPE_CHECKING:
    from .api.requests import DiagramRequest

SESSION_FORMAT = "gbdraw-session"
CURRENT_SESSION_VERSION = 35
CANONICAL_SESSION_MIN_VERSION = 31
SUPPORTED_SESSION_VERSIONS = frozenset(
    {27, 28, 29, 30, 31, 32, 33, 34, CURRENT_SESSION_VERSION}
)
PROTEIN_LOSAT_CACHE_SCHEMA = 3
NUCLEOTIDE_LOSAT_CACHE_SCHEMA = 2
LOSAT_DERIVED_CACHE_SCHEMA = 2
LEGACY_LOSAT_DERIVED_CACHE_SCHEMA = 1
PROTEIN_IDENTITY_MANIFEST_SCHEMA = 1
LEGACY_PROTEIN_CANDIDATE_SCHEMA = 1
DEPTH_FILE_ENCODING = "gbdraw-depth-table-v1"
DEPTH_FILE_SCHEMA = 1
JS_MAX_SAFE_INTEGER = 9_007_199_254_740_991

_DEPTH_COLUMNS = ("reference_name", "position", "depth")
_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
_SLOT_PART_RE = re.compile(r"([^\[\]]+)|\[(\d+)\]")


@dataclass(frozen=True)
class SessionFileBinding:
    argIndex: int
    slot: str
    name: str


@dataclass(frozen=True)
class SessionRunSpec:
    mode: Literal["circular", "linear"]
    args: tuple[str, ...]
    source_session: Mapping[str, Any]
    warnings: tuple[str, ...] = ()
    cli_invocation_args: tuple[str, ...] = ()
    file_bindings: tuple[SessionFileBinding, ...] = ()


@dataclass(frozen=True)
class SessionBuildContext:
    mode: Literal["circular", "linear"]
    output_prefix: str | None
    render_formats: tuple[str, ...]
    source_session: Mapping[str, Any] | None = None
    cli_invocation_args: tuple[str, ...] = ()
    file_bindings: tuple[SessionFileBinding | Mapping[str, Any], ...] = ()
    linear_record_metadata: tuple[Mapping[str, Any], ...] = ()


def load_session(path: str | Path) -> dict[str, Any]:
    """Load and validate a plain or gzip-compressed gbdraw GUI session JSON file."""

    session_path = Path(path)
    try:
        payload = json.loads(
            _read_session_text(session_path),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except json.JSONDecodeError as exc:
        raise ValidationError(f"Not a valid JSON session file: {session_path}") from exc
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {session_path}") from exc
    if not isinstance(payload, dict):
        raise ValidationError("Session JSON must be an object.")
    validate_session(payload)
    return payload


def _read_session_text(path: str | Path) -> str:
    """Read UTF-8 session JSON, detecting gzip by its file signature."""

    session_path = Path(path)
    with session_path.open("rb") as session_file:
        is_gzip = session_file.read(2) == b"\x1f\x8b"
    if is_gzip:
        with gzip.open(session_path, mode="rt", encoding="utf-8") as session_file:
            return session_file.read()
    return session_path.read_text(encoding="utf-8")


def validate_session(session: Mapping[str, Any]) -> None:
    """Validate the conservative session envelope used by the CLI."""

    if not isinstance(session, Mapping):
        raise ValidationError("Session JSON must be an object.")
    if session.get("format") != SESSION_FORMAT:
        raise ValidationError("Not a gbdraw-session JSON file.")
    version = session.get("version")
    if not isinstance(version, int):
        raise ValidationError("Session version is required and must be an integer.")
    if version > CURRENT_SESSION_VERSION:
        raise ValidationError(
            f"Session version {version} is newer than this gbdraw supports "
            f"({CURRENT_SESSION_VERSION})."
        )
    if version not in SUPPORTED_SESSION_VERSIONS:
        raise ValidationError(f"Unsupported session version: {version}.")
    if version >= CANONICAL_SESSION_MIN_VERSION:
        render_request = session.get("renderRequest")
        resources = session.get("resources")
        if not isinstance(render_request, Mapping):
            raise ValidationError(
                f"Session version {version} requires a canonical renderRequest object."
            )
        if not isinstance(resources, Mapping):
            raise ValidationError(
                f"Session version {version} requires a canonical resources object."
            )
        files = session.get("files")
        if files is not None and not isinstance(files, Mapping):
            raise ValidationError("Session files must be an object when present.")
    else:
        files = session.get("files")
        if files is None or not isinstance(files, Mapping):
            raise ValidationError("Session files are required for CLI regeneration.")
    if version == CURRENT_SESSION_VERSION:
        validate_current_session_artifacts(session)


def empty_protein_identity_manifest() -> dict[str, Any]:
    """Return an empty, valid protein identity manifest."""

    return {
        "schema": PROTEIN_IDENTITY_MANIFEST_SCHEMA,
        "proteinSets": {},
        "recordAnalyses": {},
        "recordInstances": {},
    }


def classify_raw_losat_cache_entry(entry: object) -> str:
    """Classify a raw LOSAT entry without guessing its schema owner."""

    if (
        not isinstance(entry, Mapping)
        or entry.get("kind") != "raw-losat"
        or not isinstance(entry.get("text"), str)
    ):
        return "invalid"
    schema = entry.get("schema")
    program = str(entry.get("program") or "").lower()
    identity_kind = entry.get("identityKind")
    from .analysis.protein_colinearity import (
        is_legacy_protein_losat_cache_entry,
        is_protein_losat_cache_entry,
    )

    if is_protein_losat_cache_entry(entry):
        return "protein-current"
    if is_legacy_protein_losat_cache_entry(entry):
        return "protein-legacy"
    if (
        schema == NUCLEOTIDE_LOSAT_CACHE_SCHEMA
        and program != "blastp"
        and identity_kind in {None, "nucleotide"}
    ):
        return "nucleotide-current"
    return "invalid"


def validate_current_session_artifacts(session: Mapping[str, Any]) -> None:
    """Validate version-35 cache, manifest, and legacy artifact boundaries."""

    cache_entries = _artifact_entries(session, "losatCache")
    protein_entries: list[Mapping[str, Any]] = []
    seen_cache_keys: set[str] = set()
    for index, entry in enumerate(cache_entries):
        classification = classify_raw_losat_cache_entry(entry)
        if classification == "protein-legacy":
            raise ValidationError(
                "Session version 35 cannot store schema-2 protein entries in losatCache; "
                "use legacyArtifacts.proteinRawCandidates."
            )
        if classification == "invalid":
            raise ValidationError(
                f"Invalid current LOSAT cache entry at losatCache.entries[{index}]."
            )
        assert isinstance(entry, Mapping)
        key = entry.get("key")
        if not isinstance(key, str) or not key:
            raise ValidationError(
                f"LOSAT cache entry at losatCache.entries[{index}] requires a key."
            )
        if key in seen_cache_keys:
            raise ValidationError(f"Duplicate LOSAT cache key: {key!r}.")
        seen_cache_keys.add(key)
        if classification == "protein-current":
            protein_entries.append(entry)

    derived_entries = _artifact_entries(session, "losatDerivedCache")
    seen_derived_keys: set[str] = set()
    for index, entry in enumerate(derived_entries):
        if not _is_derived_cache_entry(entry, schema=LOSAT_DERIVED_CACHE_SCHEMA):
            raise ValidationError(
                "Invalid current derived LOSATP cache entry at "
                f"losatDerivedCache.entries[{index}]."
            )
        assert isinstance(entry, Mapping)
        key = str(entry["key"])
        if key in seen_derived_keys:
            raise ValidationError(f"Duplicate derived LOSATP cache key: {key!r}.")
        seen_derived_keys.add(key)

    manifest = session.get("proteinIdentityManifest")
    if manifest is not None and not _is_valid_protein_identity_manifest(manifest):
        raise ValidationError("Invalid proteinIdentityManifest schema-1 artifact.")
    if protein_entries and manifest is None:
        raise ValidationError(
            "Current protein LOSATP cache entries require proteinIdentityManifest."
        )
    if protein_entries:
        assert isinstance(manifest, Mapping)
        for index, entry in enumerate(protein_entries):
            if not _protein_raw_entry_matches_manifest(entry, manifest):
                raise ValidationError(
                    "Protein LOSATP cache entry does not resolve through the manifest: "
                    f"losatCache.entries[{index}]."
                )

    legacy_artifacts = session.get("legacyArtifacts")
    if legacy_artifacts is None:
        return
    if not isinstance(legacy_artifacts, Mapping):
        raise ValidationError("Session legacyArtifacts must be an object when present.")
    candidates = legacy_artifacts.get("proteinRawCandidates")
    if candidates is not None:
        _validate_legacy_protein_candidate_envelope(candidates)
    derived_evidence = legacy_artifacts.get("proteinDerivedEvidence")
    if derived_evidence is not None:
        _validate_legacy_derived_evidence(derived_evidence)


def normalize_current_session_artifacts(
    session: dict[str, Any],
    *,
    losat_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    losat_derived_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    protein_identity_manifest: Mapping[str, Any] | None = None,
    legacy_protein_raw_candidates: Sequence[Mapping[str, Any]] | None = None,
    legacy_protein_derived_evidence: Sequence[Mapping[str, Any]] | None = None,
) -> None:
    """Normalize artifacts in-place for a current session writer.

    Legacy protein raw entries and derived schema-1 evidence are kept outside
    the current cache maps so a save-before-generate round trip is lossless.
    """

    source_raw_entries = (
        list(losat_cache_entries)
        if losat_cache_entries is not None
        else list(_artifact_entries(session, "losatCache"))
    )
    current_raw_entries: list[dict[str, Any]] = []
    imported_legacy_entries: list[dict[str, Any]] = []
    for index, entry in enumerate(source_raw_entries):
        classification = classify_raw_losat_cache_entry(entry)
        if classification in {"protein-current", "nucleotide-current"}:
            current_raw_entries.append(_json_clone(entry))
        elif classification == "protein-legacy":
            imported_legacy_entries.append(_json_clone(entry))
        else:
            raise ValidationError(
                f"Cannot write invalid LOSAT cache entry at index {index}."
            )
    session["losatCache"] = {"entries": current_raw_entries}

    source_derived_entries = (
        list(losat_derived_cache_entries)
        if losat_derived_cache_entries is not None
        else list(_artifact_entries(session, "losatDerivedCache"))
    )
    current_derived_entries: list[dict[str, Any]] = []
    imported_derived_evidence: list[dict[str, Any]] = []
    for index, entry in enumerate(source_derived_entries):
        if _is_derived_cache_entry(entry, schema=LOSAT_DERIVED_CACHE_SCHEMA):
            current_derived_entries.append(_json_clone(entry))
        elif _is_derived_cache_entry(
            entry, schema=LEGACY_LOSAT_DERIVED_CACHE_SCHEMA
        ):
            imported_derived_evidence.append(_json_clone(entry))
        else:
            raise ValidationError(
                f"Cannot write invalid derived LOSATP cache entry at index {index}."
            )
    session["losatDerivedCache"] = {"entries": current_derived_entries}

    source_manifest = (
        protein_identity_manifest
        if protein_identity_manifest is not None
        else session.get("proteinIdentityManifest")
    )
    if source_manifest is None:
        session["proteinIdentityManifest"] = empty_protein_identity_manifest()
    elif _is_valid_protein_identity_manifest(source_manifest):
        session["proteinIdentityManifest"] = _json_clone(source_manifest)
    else:
        raise ValidationError("Cannot write an invalid proteinIdentityManifest.")

    existing_legacy = session.get("legacyArtifacts")
    normalized_legacy = (
        _json_clone(existing_legacy) if isinstance(existing_legacy, Mapping) else {}
    )
    existing_candidates = normalized_legacy.get("proteinRawCandidates")
    candidate_entries = (
        list(legacy_protein_raw_candidates)
        if legacy_protein_raw_candidates is not None
        else _legacy_candidate_entries(existing_candidates)
    )
    candidate_entries.extend(
        {
            "state": "pending",
            "originalEntry": entry,
            "rejectionReason": None,
        }
        for entry in imported_legacy_entries
    )
    serializable_candidates = _normalize_legacy_candidate_entries(candidate_entries)
    if serializable_candidates:
        normalized_legacy["proteinRawCandidates"] = {
            "schema": LEGACY_PROTEIN_CANDIDATE_SCHEMA,
            "entries": serializable_candidates,
        }
    else:
        normalized_legacy.pop("proteinRawCandidates", None)

    existing_evidence = normalized_legacy.get("proteinDerivedEvidence")
    evidence_entries = (
        list(legacy_protein_derived_evidence)
        if legacy_protein_derived_evidence is not None
        else _legacy_derived_entries(existing_evidence)
    )
    evidence_entries.extend(imported_derived_evidence)
    normalized_evidence = _normalize_legacy_derived_entries(evidence_entries)
    if normalized_evidence:
        normalized_legacy["proteinDerivedEvidence"] = {
            "schema": LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
            "entries": normalized_evidence,
        }
    else:
        normalized_legacy.pop("proteinDerivedEvidence", None)

    if normalized_legacy:
        session["legacyArtifacts"] = normalized_legacy
    else:
        session.pop("legacyArtifacts", None)
    validate_current_session_artifacts(session)


def _artifact_entries(session: Mapping[str, Any], field: str) -> list[Any]:
    container = session.get(field)
    if container is None:
        return []
    if not isinstance(container, Mapping):
        raise ValidationError(f"Session {field} must be an object when present.")
    entries = container.get("entries", [])
    if not isinstance(entries, list):
        raise ValidationError(f"Session {field}.entries must be an array.")
    return entries


def _is_derived_cache_entry(entry: object, *, schema: int) -> bool:
    return (
        isinstance(entry, Mapping)
        and entry.get("schema") == schema
        and entry.get("kind") == "derived-losatp-payload"
        and isinstance(entry.get("key"), str)
        and bool(entry.get("key"))
        and isinstance(entry.get("payload"), Mapping)
    )


def _is_valid_protein_identity_manifest(manifest: object) -> bool:
    if not isinstance(manifest, Mapping):
        return False
    if manifest == empty_protein_identity_manifest():
        return True
    try:
        from .analysis.protein_colinearity import (
            validate_protein_identity_manifest,
        )

        validate_protein_identity_manifest(manifest)
    except (ImportError, ValidationError, TypeError, ValueError):
        return False
    return True


def _protein_raw_entry_matches_manifest(
    entry: Mapping[str, Any], manifest: Mapping[str, Any]
) -> bool:
    from .analysis.protein_colinearity import (
        validate_protein_raw_entry_references,
    )

    return validate_protein_raw_entry_references(entry, manifest)


def _validate_legacy_protein_candidate_envelope(envelope: object) -> None:
    if not isinstance(envelope, Mapping) or envelope.get(
        "schema"
    ) != LEGACY_PROTEIN_CANDIDATE_SCHEMA:
        raise ValidationError("Invalid legacy protein raw candidate envelope.")
    entries = envelope.get("entries")
    if not isinstance(entries, list):
        raise ValidationError("Legacy protein raw candidate entries must be an array.")
    from .analysis.protein_colinearity import (
        validate_legacy_protein_raw_candidate_envelope,
    )

    validate_legacy_protein_raw_candidate_envelope(envelope)
    for index, candidate in enumerate(entries):
        if (
            not isinstance(candidate, Mapping)
            or candidate.get("state") not in {"pending", "promoted", "rejected"}
            or classify_raw_losat_cache_entry(candidate.get("originalEntry"))
            != "protein-legacy"
            or (
                candidate.get("rejectionReason") is not None
                and not isinstance(candidate.get("rejectionReason"), str)
            )
        ):
            raise ValidationError(
                f"Invalid legacy protein raw candidate at entries[{index}]."
            )


def _validate_legacy_derived_evidence(envelope: object) -> None:
    if not isinstance(envelope, Mapping) or envelope.get(
        "schema"
    ) != LEGACY_LOSAT_DERIVED_CACHE_SCHEMA:
        raise ValidationError("Invalid legacy protein derived evidence envelope.")
    entries = envelope.get("entries")
    if not isinstance(entries, list) or not all(
        _is_derived_cache_entry(entry, schema=LEGACY_LOSAT_DERIVED_CACHE_SCHEMA)
        for entry in entries
    ):
        raise ValidationError("Invalid legacy protein derived evidence entries.")


def _legacy_candidate_entries(envelope: object) -> list[Mapping[str, Any]]:
    if envelope is None:
        return []
    _validate_legacy_protein_candidate_envelope(envelope)
    assert isinstance(envelope, Mapping)
    return list(envelope["entries"])


def _normalize_legacy_candidate_entries(
    entries: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    for candidate in entries:
        if not isinstance(candidate, Mapping) or candidate.get("state") == "promoted":
            continue
        probe = {
            "schema": LEGACY_PROTEIN_CANDIDATE_SCHEMA,
            "entries": [candidate],
        }
        _validate_legacy_protein_candidate_envelope(probe)
        clone = _json_clone(candidate)
        fingerprint = json.dumps(
            clone, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        normalized.append(clone)
    return normalized


def _legacy_derived_entries(envelope: object) -> list[Mapping[str, Any]]:
    if envelope is None:
        return []
    _validate_legacy_derived_evidence(envelope)
    assert isinstance(envelope, Mapping)
    return list(envelope["entries"])


def _normalize_legacy_derived_entries(
    entries: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, entry in enumerate(entries):
        if not _is_derived_cache_entry(
            entry, schema=LEGACY_LOSAT_DERIVED_CACHE_SCHEMA
        ):
            raise ValidationError(
                f"Invalid legacy derived LOSATP evidence at index {index}."
            )
        clone = _json_clone(entry)
        fingerprint = json.dumps(
            clone, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        normalized.append(clone)
    return normalized


def session_mode(session: Mapping[str, Any]) -> str | None:
    """Return the declared session mode when available."""

    render_request = session.get("renderRequest")
    if isinstance(render_request, Mapping):
        mode = render_request.get("mode")
        if mode in {"circular", "linear"}:
            return str(mode)

    cli_invocation = session.get("cliInvocation")
    if isinstance(cli_invocation, Mapping):
        mode = cli_invocation.get("mode")
        if mode in {"circular", "linear"}:
            return str(mode)
    ui = session.get("ui")
    if isinstance(ui, Mapping):
        mode = ui.get("mode")
        if mode in {"circular", "linear"}:
            return str(mode)
    return None


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Reject duplicate JSON object keys before the decoder can discard them."""

    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError(f"Session JSON contains a duplicate object key: {key!r}.")
        result[key] = value
    return result


def safe_embedded_filename(name: object, *, fallback: str = "embedded-file") -> str:
    """Return a basename-only filename safe for materializing embedded content."""

    raw_name = str(name or "").replace("\\", "/").split("/")[-1].strip()
    cleaned = _SAFE_FILENAME_RE.sub("_", raw_name).strip("._")
    return cleaned or fallback


def encode_depth_text(text: str) -> dict[str, Any] | None:
    """Encode samtools depth text using the browser session depth codec schema."""

    if not isinstance(text, str) or not text:
        return None
    crlf_count = text.count("\r\n")
    lf_count = text.count("\n")
    if crlf_count > 0 and crlf_count != lf_count:
        return None
    line_ending = "\r\n" if "\r\n" in text else "\n"
    normalized = text.replace("\r\n", "\n")
    if "\r" in normalized:
        return None
    final_newline = normalized.endswith("\n")
    lines = normalized.split("\n")
    if final_newline:
        lines.pop()
    if not lines or any(line == "" for line in lines):
        return None

    first_fields = lines[0].split("\t")
    if len(first_fields) != len(_DEPTH_COLUMNS):
        return None
    line_index = 0
    header: list[str] | None = None
    if _has_depth_header(first_fields):
        header = first_fields
        line_index = 1
    if line_index >= len(lines):
        return None

    records: list[dict[str, Any]] = []
    row_count = 0
    for line in lines[line_index:]:
        fields = line.split("\t")
        if len(fields) != len(_DEPTH_COLUMNS):
            return None
        position = _parse_positive_safe_integer(fields[1])
        depth_value = str(fields[2] or "").strip()
        if position is None or not _is_depth_text(depth_value):
            return None
        _append_depth_row(records, str(fields[0] or ""), position, depth_value)
        row_count += 1
    if row_count == 0:
        return None

    return {
        "schema": DEPTH_FILE_SCHEMA,
        "columns": list(_DEPTH_COLUMNS),
        "lineEnding": line_ending,
        "finalNewline": final_newline,
        "rowCount": row_count,
        "header": header,
        "records": records,
    }


def decode_depth_payload(payload: Mapping[str, Any]) -> str:
    """Decode a browser depth-file codec payload into TSV text."""

    if not isinstance(payload, Mapping) or payload.get("schema") != DEPTH_FILE_SCHEMA:
        raise ValidationError("Invalid embedded depth file.")
    records = payload.get("records")
    if not isinstance(records, list):
        raise ValidationError("Invalid embedded depth file records.")
    line_ending = "\r\n" if payload.get("lineEnding") == "\r\n" else "\n"
    lines: list[str] = []
    header = _decode_depth_header(payload.get("header"))
    if header is not None:
        lines.append(header)
    decoded_rows = 0
    for record in records:
        if not isinstance(record, Mapping) or not isinstance(record.get("runs"), list):
            raise ValidationError("Invalid embedded depth record.")
        reference_name = str(record.get("id") or "")
        for run in record["runs"]:
            decoded_rows += _decode_depth_run(reference_name, run, lines)

    declared_rows = payload.get("rowCount")
    if declared_rows is not None and declared_rows != decoded_rows:
        raise ValidationError("Embedded depth file row count does not match payload.")
    if not lines:
        return ""
    body = line_ending.join(lines)
    return body if payload.get("finalNewline") is False else f"{body}{line_ending}"


def serialize_file_entry(path: str | Path, *, depth: bool = False) -> dict[str, Any]:
    """Serialize a local file into the GUI session embedded-file shape."""

    file_path = Path(path)
    try:
        data = file_path.read_bytes()
    except OSError as exc:
        raise ValidationError(f"Could not read file for session embedding: {file_path}") from exc
    entry: dict[str, Any] = {
        "name": file_path.name or "file",
        "type": _guess_file_type(file_path),
        "size": len(data),
        "lastModified": int(file_path.stat().st_mtime * 1000),
    }
    if depth:
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            text = ""
        encoded_depth = encode_depth_text(text)
        if encoded_depth is not None:
            entry["encoding"] = DEPTH_FILE_ENCODING
            entry["data"] = encoded_depth
            return entry
    entry["data"] = base64.b64encode(data).decode("ascii")
    return entry


def materialize_embedded_file(
    entry: Mapping[str, Any],
    *,
    temp_dir: Path,
    role: str,
    prefix_role: bool = True,
) -> Path:
    """Decode one embedded session file into temp_dir and return its path."""

    temp_dir.mkdir(parents=True, exist_ok=True)
    if not isinstance(entry, Mapping):
        raise ValidationError(f"Embedded file for {role} is missing or invalid.")
    filename = safe_embedded_filename(entry.get("name"), fallback=f"{role}.dat")
    output_name = (
        f"{safe_embedded_filename(role)}-{filename}" if prefix_role else filename
    )
    output_path = temp_dir / output_name
    output_path = _assert_under_directory(output_path, temp_dir)

    if entry.get("encoding") == DEPTH_FILE_ENCODING:
        data = entry.get("data")
        if not isinstance(data, Mapping):
            raise ValidationError(f"Embedded depth payload for {role} is malformed.")
        text = decode_depth_payload(data)
        payload_bytes = text.encode("utf-8")
    else:
        data = entry.get("data")
        if not isinstance(data, str):
            raise ValidationError(f"Embedded file for {role} has no base64 data.")
        try:
            payload_bytes = base64.b64decode(data, validate=True)
        except (binascii.Error, ValueError) as exc:
            raise ValidationError(f"Embedded file for {role} has invalid base64 data.") from exc

    declared_size = entry.get("size")
    if declared_size is not None:
        try:
            expected_size = int(declared_size)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"Embedded file for {role} has invalid size metadata.") from exc
        if expected_size >= 0 and expected_size != len(payload_bytes):
            raise ValidationError(
                f"Embedded file size mismatch for {role}: expected {expected_size}, "
                f"decoded {len(payload_bytes)}."
            )
    try:
        output_path.write_bytes(payload_bytes)
    except OSError as exc:
        raise ValidationError(f"Could not materialize embedded file for {role}.") from exc
    return output_path


def session_to_cli_args(
    session: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    """Convert a GUI/CLI session into normal CLI arguments and temp files."""

    validate_session(session)
    if int(session.get("version", 0)) >= CANONICAL_SESSION_MIN_VERSION:
        raise ValidationError(
            "Canonical renderRequest sessions cannot be replayed through "
            "legacy CLI arguments."
        )
    if format_override is not None:
        format_override = ",".join(
            normalize_format_token(value) for value in format_override.split(",")
        )
    declared_mode = session_mode(session)
    if declared_mode and declared_mode != mode:
        raise ValidationError(
            f"Session mode is {declared_mode!r}; it cannot be used with the {mode} command."
        )

    cli_invocation = session.get("cliInvocation")
    if isinstance(cli_invocation, Mapping) and cli_invocation:
        return _session_cli_invocation_to_args(
            session,
            cli_invocation=cli_invocation,
            mode=mode,
            temp_dir=temp_dir,
            output_override=output_override,
            format_override=format_override,
        )
    return _gui_session_to_cli_args(
        session,
        mode=mode,
        temp_dir=temp_dir,
        output_override=output_override,
        format_override=format_override,
    )


def build_session_json(
    context: SessionBuildContext,
    *,
    svg_results: Sequence[tuple[str, str]],
    embedded_files: Mapping[str, Any],
    generated_at: datetime,
    losat_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    losat_derived_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    protein_identity_manifest: Mapping[str, Any] | None = None,
    legacy_protein_raw_candidates: Sequence[Mapping[str, Any]] | None = None,
    legacy_protein_derived_evidence: Sequence[Mapping[str, Any]] | None = None,
    canonical_request: DiagramRequest | None = None,
) -> dict[str, Any]:
    """Build a GUI-loadable session JSON payload from a CLI run."""

    if context.source_session is not None:
        payload: dict[str, Any] = _json_clone(context.source_session)
    else:
        payload = {}

    payload["format"] = SESSION_FORMAT
    payload["version"] = CURRENT_SESSION_VERSION
    payload["createdAt"] = generated_at.isoformat()
    if context.output_prefix:
        payload["title"] = Path(str(context.output_prefix)).name
    else:
        payload.setdefault("title", "gbdraw")

    config = payload.get("config")
    if not isinstance(config, dict):
        config = _minimal_config_from_cli_args(context)
        payload["config"] = config
    _update_config_prefix(config, context.output_prefix)

    ui = payload.get("ui")
    if not isinstance(ui, dict):
        ui = {}
        payload["ui"] = ui
    ui["mode"] = context.mode
    ui.setdefault("zoom", 1)
    ui.setdefault("selectedResultIndex", 0)
    ui.setdefault("canvasPan", {"x": 0, "y": 0})
    ui.setdefault("canvasPadding", {"top": 0, "right": 0, "bottom": 0, "left": 0})
    if context.mode == "circular":
        ui.setdefault("cInputType", _input_type_from_args(context.cli_invocation_args))
    else:
        ui.setdefault("lInputType", _input_type_from_args(context.cli_invocation_args))

    payload["files"] = _json_clone(embedded_files)
    payload["results"] = [
        {"name": name or f"Result {index + 1}", "content": content}
        for index, (name, content) in enumerate(svg_results)
    ]
    payload.setdefault("features", {})
    payload.setdefault("orthogroupState", {})
    if context.mode == "linear":
        _populate_linear_session_fields_from_cli_context(payload, context)
    payload["cliInvocation"] = {
        "schema": 1,
        "mode": context.mode,
        "args": [str(arg) for arg in context.cli_invocation_args],
        "renderFormats": [normalize_format_token(fmt) for fmt in context.render_formats],
        "fileBindings": [_binding_to_json(binding) for binding in context.file_bindings],
        "generatedBy": "gbdraw",
    }
    if canonical_request is not None:
        from .session import build_session_document

        canonical = build_session_document(
            canonical_request,
            created_at=generated_at,
        ).to_dict()
        payload["renderRequest"] = canonical["renderRequest"]
        payload["resources"] = canonical["resources"]
    elif context.source_session is not None and isinstance(
        context.source_session.get("renderRequest"), Mapping
    ) and isinstance(context.source_session.get("resources"), Mapping):
        payload["renderRequest"] = _json_clone(context.source_session["renderRequest"])
        payload["resources"] = _json_clone(context.source_session["resources"])
    else:
        raise ValidationError(
            f"A canonical typed request is required to write a version {CURRENT_SESSION_VERSION} session."
        )
    normalize_current_session_artifacts(
        payload,
        losat_cache_entries=losat_cache_entries,
        losat_derived_cache_entries=losat_derived_cache_entries,
        protein_identity_manifest=protein_identity_manifest,
        legacy_protein_raw_candidates=legacy_protein_raw_candidates,
        legacy_protein_derived_evidence=legacy_protein_derived_evidence,
    )
    validate_session(payload)
    return payload


def write_session_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    """Write plain or ``.gz`` session JSON with an atomic replacement."""

    if payload.get("version") == CURRENT_SESSION_VERSION:
        validate_session(payload)

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    try:
        if output_path.suffix.lower() == ".gz":
            with temp_path.open("wb") as raw_file:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw_file,
                    compresslevel=6,
                    mtime=0,
                ) as compressed_file:
                    with io.TextIOWrapper(compressed_file, encoding="utf-8") as text_file:
                        json.dump(
                            payload,
                            text_file,
                            ensure_ascii=False,
                            separators=(",", ":"),
                        )
        else:
            temp_path.write_text(
                json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
                encoding="utf-8",
            )
        temp_path.replace(output_path)
    except OSError as exc:
        raise ValidationError(f"Could not write session sidecar: {output_path}") from exc
    finally:
        try:
            if temp_path.exists():
                temp_path.unlink()
        except OSError:
            pass


def get_session_slot(session: Mapping[str, Any], slot: str) -> Any:
    """Resolve a slot path such as files.linearSeqs[0].gb inside a session."""

    current: Any = session
    for part in _parse_slot(slot):
        if isinstance(part, int):
            if not isinstance(current, Sequence) or isinstance(current, (str, bytes, bytearray)):
                raise ValidationError(f"Session file binding slot is not a list: {slot}")
            if part < 0 or part >= len(current):
                raise ValidationError(f"Session file binding slot index is out of range: {slot}")
            current = current[part]
        else:
            if not isinstance(current, Mapping) or part not in current:
                raise ValidationError(f"Session file binding slot is missing: {slot}")
            current = current[part]
    return current


def _session_cli_invocation_to_args(
    session: Mapping[str, Any],
    *,
    cli_invocation: Mapping[str, Any],
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    if cli_invocation.get("schema") != 1:
        raise ValidationError("Unsupported cliInvocation schema.")
    invocation_mode = cli_invocation.get("mode")
    if invocation_mode != mode:
        raise ValidationError(
            f"Session cliInvocation mode is {invocation_mode!r}; expected {mode!r}."
        )
    raw_args = cli_invocation.get("args")
    if not isinstance(raw_args, list) or not all(isinstance(arg, str) for arg in raw_args):
        raise ValidationError("Session cliInvocation args must be a string array.")

    invocation_args = [str(arg) for arg in raw_args]
    run_args = list(invocation_args)
    file_bindings = _normalize_file_bindings(cli_invocation.get("fileBindings"))
    for binding in file_bindings:
        if binding.argIndex < 0 or binding.argIndex >= len(run_args):
            raise ValidationError(
                f"cliInvocation.fileBindings argIndex {binding.argIndex} is out of range."
            )
        entry = get_session_slot(session, binding.slot)
        materialized = materialize_embedded_file(
            entry,
            temp_dir=temp_dir,
            role=f"arg{binding.argIndex}",
        )
        run_args[binding.argIndex] = str(materialized)

    _restore_cli_table_paths(session, run_args, temp_dir=temp_dir)

    run_args = _apply_option_override(run_args, "-o", "--output", output_override)
    run_args = _apply_option_override(run_args, "-f", "--format", format_override)
    invocation_args = _apply_option_override(invocation_args, "-o", "--output", output_override)
    invocation_args = _apply_option_override(invocation_args, "-f", "--format", format_override)
    session_version = int(session.get("version", 0))
    run_args = migrate_legacy_repeat_feature_shape_args(
        run_args,
        session_version=session_version,
    )
    invocation_args = migrate_legacy_repeat_feature_shape_args(
        invocation_args,
        session_version=session_version,
    )

    return SessionRunSpec(
        mode=mode,
        args=tuple(run_args),
        source_session=session,
        cli_invocation_args=tuple(invocation_args),
        file_bindings=tuple(file_bindings),
    )


def _gui_session_to_cli_args(
    session: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    config = session.get("config")
    if not isinstance(config, Mapping):
        raise ValidationError("GUI session config is required when cliInvocation is absent.")
    files = session.get("files")
    if not isinstance(files, Mapping):
        raise ValidationError("GUI session files are required.")
    ui = session.get("ui")
    if not isinstance(ui, Mapping):
        ui = {}
    form = config.get("form") if isinstance(config.get("form"), Mapping) else {}
    adv = dict(config.get("adv")) if isinstance(config.get("adv"), Mapping) else {}
    if int(session.get("version", 0)) <= 30:
        effective_features = adv.get("features")
        if not isinstance(effective_features, list):
            effective_features = [
                "CDS",
                "rRNA",
                "tRNA",
                "tmRNA",
                "ncRNA",
                "misc_RNA",
                "repeat_region",
            ]
        if "repeat_region" in effective_features:
            feature_shapes = dict(adv.get("feature_shapes") or {})
            feature_shapes.setdefault("repeat_region", "rectangle")
            adv["feature_shapes"] = feature_shapes

    run_args: list[str] = []
    invocation_args: list[str] = []
    bindings: list[SessionFileBinding] = []

    output_prefix = output_override or _string_or_none(form.get("prefix"))
    if output_prefix:
        _append_pair(run_args, invocation_args, "-o", output_prefix)
    _append_pair(run_args, invocation_args, "-f", format_override or "svg")
    _append_common_gui_args(run_args, invocation_args, config=config, form=form, adv=adv)

    if mode == "circular":
        _append_circular_gui_args(
            run_args,
            invocation_args,
            bindings,
            session=session,
            files=files,
            ui=ui,
            form=form,
            adv=adv,
            temp_dir=temp_dir,
        )
    else:
        _append_linear_gui_args(
            run_args,
            invocation_args,
            bindings,
            session=session,
            files=files,
            ui=ui,
            config=config,
            form=form,
            adv=adv,
            temp_dir=temp_dir,
        )

    return SessionRunSpec(
        mode=mode,
        args=tuple(run_args),
        source_session=session,
        cli_invocation_args=tuple(invocation_args),
        file_bindings=tuple(bindings),
    )


def _restore_cli_table_paths(
    session: Mapping[str, Any],
    run_args: list[str],
    *,
    temp_dir: Path,
) -> None:
    files = session.get("files")
    if not isinstance(files, Mapping):
        return
    cli_tables = files.get("cliTables")
    if not isinstance(cli_tables, list):
        return

    for table_entry in cli_tables:
        if not isinstance(table_entry, Mapping):
            continue
        try:
            arg_index = int(table_entry.get("argIndex"))
        except (TypeError, ValueError) as exc:
            raise ValidationError("files.cliTables argIndex must be an integer.") from exc
        if arg_index < 0 or arg_index >= len(run_args):
            raise ValidationError("files.cliTables argIndex is out of range.")

        table_path = Path(str(run_args[arg_index]))
        if not table_path.is_file():
            table_slot = str(table_entry.get("slot") or "").strip()
            if not table_slot:
                raise ValidationError("files.cliTables slot is required.")
            table_file_entry = get_session_slot(session, table_slot)
            table_path = materialize_embedded_file(
                table_file_entry,
                temp_dir=temp_dir,
                role=f"arg{arg_index}",
            )
            run_args[arg_index] = str(table_path)

        dependencies = table_entry.get("dependencies")
        if not isinstance(dependencies, list) or not dependencies:
            continue
        replacements: dict[tuple[int, str], str] = {}
        for dependency in dependencies:
            if not isinstance(dependency, Mapping):
                continue
            try:
                row_index = int(dependency.get("rowIndex"))
            except (TypeError, ValueError) as exc:
                raise ValidationError("files.cliTables dependencies rowIndex must be an integer.") from exc
            column = str(dependency.get("column") or "").strip()
            slot = str(dependency.get("slot") or "").strip()
            if not column or not slot:
                raise ValidationError("files.cliTables dependency entries require column and slot.")
            dependency_entry = get_session_slot(session, slot)
            materialized = materialize_embedded_file(
                dependency_entry,
                temp_dir=temp_dir,
                role=slot.replace(".", "_").replace("[", "_").replace("]", ""),
            )
            replacements[(row_index, column)] = _relative_path_for_table(
                materialized,
                table_path.parent,
            )
        if replacements:
            _rewrite_tsv_path_cells(table_path, replacements)


def _relative_path_for_table(path: Path, table_dir: Path) -> str:
    try:
        return os.path.relpath(str(path), str(table_dir))
    except ValueError:
        return str(path)


def _rewrite_tsv_path_cells(
    table_path: Path,
    replacements: Mapping[tuple[int, str], str],
) -> None:
    try:
        with table_path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
    except OSError as exc:
        raise ValidationError(f"Could not read restored TSV table: {table_path}") from exc

    header: list[str] | None = None
    output_rows: list[list[str]] = []
    data_row_index = 0
    for cells in rows:
        if not cells or all(str(cell).strip() == "" for cell in cells):
            continue
        if header is None:
            header = [str(cell).strip() for cell in cells]
            output_rows.append(header)
            continue
        values = [str(cell) for cell in cells]
        while len(values) < len(header):
            values.append("")
        for (target_row_index, column), replacement in replacements.items():
            if target_row_index != data_row_index or column not in header:
                continue
            values[header.index(column)] = replacement
        output_rows.append(values[: len(header)])
        data_row_index += 1

    if header is None:
        raise ValidationError(f"Restored TSV table has no header row: {table_path}")
    try:
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(output_rows)
    except OSError as exc:
        raise ValidationError(f"Could not rewrite restored TSV table: {table_path}") from exc


def _append_common_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    config: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    if _string_or_none(form.get("species")):
        _append_pair(run_args, invocation_args, "--species", str(form.get("species")))
    if _string_or_none(form.get("strain")):
        _append_pair(run_args, invocation_args, "--strain", str(form.get("strain")))
    if form.get("separate_strands") is True:
        _append_flag(run_args, invocation_args, "--separate_strands")
    features = adv.get("features")
    if isinstance(features, list) and features:
        _append_pair(run_args, invocation_args, "-k", ",".join(str(item) for item in features if item))
    feature_shapes = adv.get("feature_shapes")
    if isinstance(feature_shapes, Mapping):
        for feature_type, rendering in feature_shapes.items():
            _append_pair(
                run_args,
                invocation_args,
                "--feature_shape",
                f"{str(feature_type).strip()}={str(rendering).strip().lower()}",
            )
    for key, option in (
        ("window_size", "--window"),
        ("step_size", "--step"),
        ("nt", "--nt"),
        ("def_font_size", "--definition_font_size"),
        ("label_font_size", "--label_font_size"),
        ("block_stroke_width", "--block_stroke_width"),
        ("block_stroke_color", "--block_stroke_color"),
        ("line_stroke_width", "--line_stroke_width"),
        ("line_stroke_color", "--line_stroke_color"),
        ("axis_stroke_width", "--axis_stroke_width"),
        ("axis_stroke_color", "--axis_stroke_color"),
        ("legend_box_size", "--legend_box_size"),
        ("legend_font_size", "--legend_font_size"),
        ("scale_interval", "--scale_interval"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))
    if adv.get("resolve_overlaps") is True:
        _append_flag(run_args, invocation_args, "--resolve_overlaps")
    if adv.get("gc_content_mode") == "percent":
        _append_pair(run_args, invocation_args, "--gc_content_mode", "percent")
        for key, option in (
            ("gc_content_min_percent", "--gc_content_min_percent"),
            ("gc_content_max_percent", "--gc_content_max_percent"),
            ("gc_content_tick_interval", "--gc_content_large_tick_interval"),
            ("gc_content_small_tick_interval", "--gc_content_small_tick_interval"),
            ("gc_content_tick_font_size", "--gc_content_tick_font_size"),
        ):
            value = adv.get(key)
            if value not in (None, "", False):
                _append_pair(run_args, invocation_args, option, str(value))
        if adv.get("gc_content_show_axis") is False:
            _append_flag(run_args, invocation_args, "--hide_gc_content_axis")
        if adv.get("gc_content_show_ticks") is False:
            _append_flag(run_args, invocation_args, "--hide_gc_content_ticks")


def _append_circular_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    files: Mapping[str, Any],
    ui: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    if _string_or_none(form.get("track_type")):
        _append_pair(run_args, invocation_args, "--track_type", str(form.get("track_type")))
    if _string_or_none(form.get("legend")):
        _append_pair(run_args, invocation_args, "-l", str(form.get("legend")))
    plot_title = _string_or_none(form.get("plot_title"))
    if plot_title:
        _append_pair(run_args, invocation_args, "--plot_title", plot_title)
    if _string_or_none(adv.get("plot_title_position")):
        _append_pair(run_args, invocation_args, "--plot_title_position", str(adv.get("plot_title_position")))
    if adv.get("plot_title_font_size") not in (None, "", False):
        _append_pair(run_args, invocation_args, "--plot_title_font_size", str(adv.get("plot_title_font_size")))
    if adv.get("keep_full_definition_with_plot_title") is True:
        _append_flag(run_args, invocation_args, "--keep_full_definition_with_plot_title")
    if adv.get("center_reserved_radius") not in (None, "", False):
        _append_pair(run_args, invocation_args, "--center_reserved_radius", str(adv.get("center_reserved_radius")))
    labels_mode = str(form.get("labels_mode") or "none")
    if labels_mode == "out":
        _append_flag(run_args, invocation_args, "--labels")
    elif labels_mode == "both":
        _append_pair(run_args, invocation_args, "--labels", "both")
    if form.get("suppress_gc") is True:
        _append_flag(run_args, invocation_args, "--suppress_gc")
    if form.get("suppress_skew") is True:
        _append_flag(run_args, invocation_args, "--suppress_skew")
    if form.get("multi_record_canvas") is True:
        _append_flag(run_args, invocation_args, "--multi_record_canvas")
    for key, option in (
        ("multi_record_size_mode", "--multi_record_size_mode"),
        ("multi_record_min_radius_ratio", "--multi_record_min_radius_ratio"),
        ("multi_record_column_gap_ratio", "--multi_record_column_gap_ratio"),
        ("multi_record_row_gap_ratio", "--multi_record_row_gap_ratio"),
        ("feature_width_circular", "--feature_width"),
        ("depth_width_circular", "--depth_width"),
        ("gc_content_width_circular", "--gc_content_width"),
        ("gc_content_radius_circular", "--gc_content_radius"),
        ("gc_skew_width_circular", "--gc_skew_width"),
        ("gc_skew_radius_circular", "--gc_skew_radius"),
        ("tick_label_font_size", "--tick_label_font_size"),
        ("circular_label_spacing", "--circular_label_spacing"),
        ("circular_label_placement", "--label_placement"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))

    input_type = str(ui.get("cInputType") or "gb")
    if input_type == "gb":
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_gb",
            option="--gbk",
            temp_dir=temp_dir,
        )
    else:
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_gff",
            option="--gff",
            temp_dir=temp_dir,
        )
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_fasta",
            option="--fasta",
            temp_dir=temp_dir,
        )

    _append_depth_gui_options(
        run_args,
        invocation_args,
        bindings,
        session=session,
        slot_prefix="files.c_depth",
        option="--depth_track",
        temp_dir=temp_dir,
        show_depth=bool(form.get("show_depth")),
    )

    conservation_blasts = files.get("c_conservation_blasts")
    if isinstance(conservation_blasts, list) and conservation_blasts:
        _append_flag(run_args, invocation_args, "--conservation_blast")
        for index, entry in enumerate(conservation_blasts):
            if entry:
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.c_conservation_blasts[{index}]",
                    temp_dir=temp_dir,
                )


def _append_linear_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    files: Mapping[str, Any],
    ui: Mapping[str, Any],
    config: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    if _string_or_none(form.get("scale_style")):
        _append_pair(run_args, invocation_args, "--scale_style", str(form.get("scale_style")))
    if form.get("align_center") is True:
        _append_flag(run_args, invocation_args, "--align_center")
    if form.get("show_gc") is True:
        _append_flag(run_args, invocation_args, "--show_gc")
    if form.get("show_skew") is True:
        _append_flag(run_args, invocation_args, "--show_skew")
    if form.get("normalize_length") is True:
        _append_flag(run_args, invocation_args, "--normalize_length")
    if _string_or_none(form.get("legend")) and form.get("legend") != "right":
        _append_pair(run_args, invocation_args, "-l", str(form.get("legend")))
    labels_mode = str(form.get("show_labels_linear") or "none")
    if labels_mode == "all":
        _append_flag(run_args, invocation_args, "--show_labels")
    elif labels_mode == "first":
        _append_pair(run_args, invocation_args, "--show_labels", "first")
    elif labels_mode == "orthogroup_top":
        _append_pair(run_args, invocation_args, "--show_labels", "orthogroup_top")
    for key, option in (
        ("feature_height", "--feature_height"),
        ("gc_height", "--gc_height"),
        ("comparison_height", "--comparison_height"),
        ("scale_font_size", "--scale_font_size"),
        ("scale_stroke_width", "--scale_stroke_width"),
        ("scale_stroke_color", "--scale_stroke_color"),
        ("ruler_label_color", "--ruler_label_color"),
        ("pairwise_match_style", "--pairwise_match_style"),
        ("track_axis_gap", "--track_axis_gap"),
        ("label_placement", "--label_placement"),
        ("label_rendering", "--label_rendering"),
        ("label_rotation", "--label_rotation"),
        ("linear_label_spacing", "--linear_label_spacing"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))
    if _string_or_none(form.get("linear_track_layout")):
        _append_pair(run_args, invocation_args, "--track_layout", str(form.get("linear_track_layout")))
    if form.get("linear_ruler_on_axis") is True:
        _append_flag(run_args, invocation_args, "--ruler_on_axis")
    _append_linear_definition_line_style_args(run_args, invocation_args, adv)

    linear_seqs = files.get("linearSeqs")
    if not isinstance(linear_seqs, list) or not linear_seqs:
        raise ValidationError("Linear GUI session has no embedded sequence files.")
    input_type = str(ui.get("lInputType") or "gb")
    if input_type == "gb":
        _append_flag(run_args, invocation_args, "--gbk")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("gb"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].gb",
                    temp_dir=temp_dir,
                )
    else:
        _append_flag(run_args, invocation_args, "--gff")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("gff"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].gff",
                    temp_dir=temp_dir,
                )
        _append_flag(run_args, invocation_args, "--fasta")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("fasta"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].fasta",
                    temp_dir=temp_dir,
                )
    blast_slots = [
        index for index, seq in enumerate(linear_seqs)
        if isinstance(seq, Mapping) and seq.get("blast")
    ]
    if blast_slots:
        _append_flag(run_args, invocation_args, "-b")
        for index in blast_slots:
            _append_materialized_value(
                run_args,
                invocation_args,
                bindings,
                session=session,
                slot=f"files.linearSeqs[{index}].blast",
                temp_dir=temp_dir,
            )
    elif _gui_linear_losat_program(config, adv) == "blastp":
        _append_linear_gui_blastp_args(
            run_args,
            invocation_args,
            session=session,
            config=config,
            adv=adv,
        )
    _append_linear_gui_sequence_options(
        run_args,
        invocation_args,
        linear_seqs=linear_seqs,
    )
    if form.get("show_depth") is True:
        depth_rows = [
            _as_list(seq.get("depth") if isinstance(seq, Mapping) else None)
            for seq in linear_seqs
        ]
        track_count = max((len(row) for row in depth_rows), default=0)
        for track_index in range(track_count):
            _append_flag(run_args, invocation_args, "--depth_track")
            for record_index, row in enumerate(depth_rows):
                if track_index >= len(row) or not row[track_index]:
                    _append_value(run_args, invocation_args, "none")
                    continue
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{record_index}].depth"
                    + (f"[{track_index}]" if isinstance(linear_seqs[record_index].get("depth"), list) else ""),
                    temp_dir=temp_dir,
                )
        if track_count:
            _append_flag(run_args, invocation_args, "--show_depth")


def _append_linear_definition_line_style_args(
    run_args: list[str],
    invocation_args: list[str],
    adv: Mapping[str, Any],
) -> None:
    styles = adv.get("linear_definition_line_styles")
    if not isinstance(styles, Mapping):
        return
    for line_kind in DEFINITION_LINE_KINDS:
        raw_style = styles.get(line_kind)
        if not isinstance(raw_style, Mapping):
            continue
        parts: list[str] = []
        font_size = raw_style.get("font_size")
        if font_size not in (None, "", False):
            parts.append(f"size={font_size}")
        font_weight = str(raw_style.get("font_weight") or "").strip()
        if font_weight.lower() in {"auto", "none", "null", "default", "normal"}:
            font_weight = ""
        if font_weight:
            parts.append(f"weight={font_weight}")
        fill = str(raw_style.get("fill") or "").strip()
        if fill:
            parts.append(f"color={fill}")
        if parts:
            _append_pair(run_args, invocation_args, "--definition_line_style", f"{line_kind}:{','.join(parts)}")


def _append_linear_gui_sequence_options(
    run_args: list[str],
    invocation_args: list[str],
    *,
    linear_seqs: Sequence[Any],
) -> None:
    labels = [
        str(seq.get("definition") or "") if isinstance(seq, Mapping) else ""
        for seq in linear_seqs
    ]
    if any(label.strip() for label in labels):
        for label in labels:
            _append_pair(run_args, invocation_args, "--record_label", label)

    subtitles = [
        str(seq.get("record_subtitle") or "") if isinstance(seq, Mapping) else ""
        for seq in linear_seqs
    ]
    if any(subtitle.strip() for subtitle in subtitles):
        for subtitle in subtitles:
            _append_pair(run_args, invocation_args, "--record_subtitle", subtitle)

    record_selectors: list[str] = []
    reverse_flags: list[bool] = []
    region_specs: list[str] = []
    for index, seq in enumerate(linear_seqs):
        if not isinstance(seq, Mapping):
            record_selectors.append("")
            reverse_flags.append(False)
            continue
        record_selector = str(seq.get("region_record_id") or "").strip()
        record_selectors.append(record_selector)
        start = seq.get("region_start")
        end = seq.get("region_end")
        has_start = start not in (None, "")
        has_end = end not in (None, "")
        if has_start != has_end:
            raise ValidationError(
                f"Linear sequence #{index + 1} has an incomplete region start/end."
            )
        wants_reverse = bool(seq.get("region_reverse"))
        if has_start and has_end:
            try:
                start_int = int(start)
                end_int = int(end)
            except (TypeError, ValueError) as exc:
                raise ValidationError(
                    f"Linear sequence #{index + 1} has invalid region coordinates."
                ) from exc
            if start_int < 1 or end_int < 1:
                raise ValidationError(
                    f"Linear sequence #{index + 1} region coordinates must be >= 1."
                )
            suffix = ":rc" if wants_reverse else ""
            region_specs.append(f"#{index + 1}:{start_int}-{end_int}{suffix}")
            reverse_flags.append(False)
        else:
            reverse_flags.append(wants_reverse)

    if any(selector for selector in record_selectors):
        for selector in record_selectors:
            _append_pair(run_args, invocation_args, "--record_id", selector)
    if any(reverse_flags):
        for flag in reverse_flags:
            _append_pair(run_args, invocation_args, "--reverse_complement", "1" if flag else "0")
    for spec in region_specs:
        _append_pair(run_args, invocation_args, "--region", spec)


def _gui_linear_losat_program(config: Mapping[str, Any], adv: Mapping[str, Any]) -> str:
    blast_source = str(config.get("blastSource") or adv.get("blastSource") or "").strip().lower()
    losat_program = str(config.get("losatProgram") or adv.get("losatProgram") or "").strip().lower()
    if blast_source != "losat":
        return ""
    return losat_program


def _append_linear_gui_blastp_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    session: Mapping[str, Any],
    config: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    losat_cfg = config.get("losat")
    if not isinstance(losat_cfg, Mapping):
        return
    blastp_cfg = losat_cfg.get("blastp")
    if not isinstance(blastp_cfg, Mapping):
        return
    mode = str(blastp_cfg.get("mode") or "none").strip().lower()
    if mode not in {"pairwise", "orthogroup", "collinear"}:
        return
    _append_pair(run_args, invocation_args, "--protein_blastp_mode", mode)
    threads_per_job = str(losat_cfg.get("threadsPerJob") or "auto").strip().lower()
    if threads_per_job != "auto":
        try:
            parsed_threads = int(threads_per_job)
        except ValueError:
            parsed_threads = 0
        if parsed_threads >= 1:
            _append_pair(run_args, invocation_args, "--losatp_threads", str(parsed_threads))

    max_hits = blastp_cfg.get("maxHits")
    if max_hits not in (None, "", False):
        _append_pair(run_args, invocation_args, "--protein_blastp_max_hits", str(max_hits))
        if mode == "pairwise":
            _append_pair(run_args, invocation_args, "--protein_blastp_candidate_limit", str(max_hits))
    candidate_limit = blastp_cfg.get("candidateLimit")
    if mode != "pairwise" and candidate_limit not in (None, "", False):
        _append_pair(run_args, invocation_args, "--protein_blastp_candidate_limit", str(candidate_limit))

    for key, option in (
        ("min_bitscore", "--bitscore"),
        ("evalue", "--evalue"),
        ("identity", "--identity"),
        ("alignment_length", "--alignment_length"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))

    if mode == "orthogroup":
        orthogroup_state = session.get("orthogroupState")
        selected_target = (
            str(orthogroup_state.get("selectedOrthogroupAlignmentFeature") or "").strip()
            if isinstance(orthogroup_state, Mapping)
            else ""
        )
        if selected_target:
            _append_pair(run_args, invocation_args, "--align_orthogroup_feature", selected_target)

    if mode != "collinear":
        return
    for key, option in (
        ("collinearMinAnchors", "--collinear_min_anchors"),
        ("collinearMaxGeneGap", "--collinear_max_gene_gap"),
        ("collinearMaxDiagonalDrift", "--collinear_max_diagonal_drift"),
        ("collinearMaxConflictsInMergeGap", "--collinear_max_conflicts_in_merge_gap"),
        ("collinearUnitMode", "--collinear_unit_mode"),
        ("collinearSearchScope", "--collinear_search_scope"),
        ("collinearColorMode", "--collinear_color_mode"),
        ("collinearMaxParalogLinksPerOrthogroup", "--collinear_max_paralog_links_per_orthogroup"),
    ):
        value = blastp_cfg.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))


def _append_depth_gui_options(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot_prefix: str,
    option: str,
    temp_dir: Path,
    show_depth: bool,
) -> None:
    if not show_depth:
        return
    entries = _as_list(get_session_slot(session, slot_prefix))
    for index, entry in enumerate(entries):
        if not entry:
            continue
        _append_flag(run_args, invocation_args, option)
        slot = slot_prefix if len(entries) == 1 else f"{slot_prefix}[{index}]"
        _append_materialized_value(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot=slot,
            temp_dir=temp_dir,
        )
    if entries:
        _append_flag(run_args, invocation_args, "--show_depth")


def _append_materialized_file_option(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot: str,
    option: str,
    temp_dir: Path,
) -> None:
    _append_flag(run_args, invocation_args, option)
    _append_materialized_value(
        run_args,
        invocation_args,
        bindings,
        session=session,
        slot=slot,
        temp_dir=temp_dir,
    )


def _append_materialized_value(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot: str,
    temp_dir: Path,
) -> None:
    entry = get_session_slot(session, slot)
    path = materialize_embedded_file(
        entry,
        temp_dir=temp_dir,
        role=slot.replace(".", "_").replace("[", "_").replace("]", ""),
    )
    arg_index = len(run_args)
    name = safe_embedded_filename(entry.get("name") if isinstance(entry, Mapping) else "")
    run_args.append(str(path))
    invocation_args.append(name)
    bindings.append(SessionFileBinding(argIndex=arg_index, slot=slot, name=name))


def _append_flag(run_args: list[str], invocation_args: list[str], option: str) -> None:
    run_args.append(str(option))
    invocation_args.append(str(option))


def _append_pair(
    run_args: list[str],
    invocation_args: list[str],
    option: str,
    value: object,
) -> None:
    run_args.extend([str(option), str(value)])
    invocation_args.extend([str(option), str(value)])


def _append_value(run_args: list[str], invocation_args: list[str], value: object) -> None:
    run_args.append(str(value))
    invocation_args.append(str(value))


def _apply_option_override(
    args: list[str],
    short_option: str,
    long_option: str,
    value: str | None,
) -> list[str]:
    if value is None:
        return list(args)
    result: list[str] = []
    replaced = False
    replace_index: int | None = None
    for index, token in enumerate(args[:-1]):
        if token in {short_option, long_option}:
            replace_index = index + 1
    for index, token in enumerate(args):
        if replace_index is not None and index == replace_index:
            result.append(str(value))
            replaced = True
        else:
            result.append(token)
    if not replaced:
        result.extend([short_option, str(value)])
    return result


def _normalize_file_bindings(value: Any) -> list[SessionFileBinding]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValidationError("cliInvocation.fileBindings must be an array.")
    bindings: list[SessionFileBinding] = []
    for item in value:
        if not isinstance(item, Mapping):
            raise ValidationError("cliInvocation.fileBindings entries must be objects.")
        try:
            arg_index = int(item.get("argIndex"))
        except (TypeError, ValueError) as exc:
            raise ValidationError("cliInvocation.fileBindings argIndex must be an integer.") from exc
        slot = str(item.get("slot") or "").strip()
        if not slot:
            raise ValidationError("cliInvocation.fileBindings slot is required.")
        name = safe_embedded_filename(item.get("name"), fallback="file")
        bindings.append(SessionFileBinding(argIndex=arg_index, slot=slot, name=name))
    return bindings


def _binding_to_json(binding: SessionFileBinding | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(binding, SessionFileBinding):
        return {
            "argIndex": binding.argIndex,
            "slot": binding.slot,
            "name": binding.name,
        }
    return {
        "argIndex": int(binding.get("argIndex", 0)),
        "slot": str(binding.get("slot", "")),
        "name": safe_embedded_filename(binding.get("name"), fallback="file"),
    }


def _parse_slot(slot: str) -> list[str | int]:
    normalized = str(slot or "").strip()
    if not normalized:
        raise ValidationError("Session slot cannot be empty.")
    parts: list[str | int] = []
    for raw_part in normalized.split("."):
        if not raw_part:
            raise ValidationError(f"Invalid session slot: {slot}")
        position = 0
        for match in _SLOT_PART_RE.finditer(raw_part):
            if match.start() != position:
                raise ValidationError(f"Invalid session slot: {slot}")
            position = match.end()
            key, index = match.groups()
            if key is not None:
                parts.append(key)
            elif index is not None:
                parts.append(int(index))
        if position != len(raw_part):
            raise ValidationError(f"Invalid session slot: {slot}")
    return parts


def _assert_under_directory(path: Path, directory: Path) -> Path:
    resolved_directory = directory.resolve()
    resolved_path = path.resolve()
    try:
        resolved_path.relative_to(resolved_directory)
    except ValueError as exc:
        raise ValidationError("Embedded filename cannot be safely materialized.") from exc
    return resolved_path


def _decode_depth_header(header: Any) -> str | None:
    if header is None:
        return None
    if (
        not isinstance(header, list)
        or len(header) != len(_DEPTH_COLUMNS)
    ):
        raise ValidationError("Invalid embedded depth file header.")
    return "\t".join(str(value if value is not None else "") for value in header)


def _decode_depth_run(reference_name: str, run: Any, lines: list[str]) -> int:
    if not isinstance(run, list) or len(run) != 4 or not isinstance(run[3], list):
        raise ValidationError("Invalid embedded depth run.")
    start, step, count, depths = run
    for value in (start, step, count):
        if not isinstance(value, int) or value <= 0 or value > JS_MAX_SAFE_INTEGER:
            raise ValidationError("Invalid embedded depth coordinates.")
    if len(depths) != count:
        raise ValidationError("Invalid embedded depth coordinates.")
    for index, depth_value in enumerate(depths):
        position = start + step * index
        if position > JS_MAX_SAFE_INTEGER:
            raise ValidationError("Invalid embedded depth coordinate overflow.")
        lines.append(f"{reference_name}\t{position}\t{'' if depth_value is None else depth_value}")
    return count


def _parse_positive_safe_integer(value: object) -> int | None:
    text = str(value or "").strip()
    if not re.fullmatch(r"[+-]?\d+", text):
        return None
    parsed = int(text)
    if parsed <= 0 or parsed > JS_MAX_SAFE_INTEGER:
        return None
    return parsed


def _is_depth_text(value: object) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    try:
        parsed = float(text)
    except ValueError:
        return False
    return math.isfinite(parsed) and parsed >= 0


def _has_depth_header(fields: Sequence[str]) -> bool:
    return (
        len(fields) >= 3
        and (_parse_positive_safe_integer(fields[1]) is None or not _is_depth_text(fields[2]))
    )


def _append_depth_row(
    records: list[dict[str, Any]],
    reference_name: str,
    position: int,
    depth_value: str,
) -> None:
    if not records or records[-1]["id"] != reference_name:
        records.append({"id": reference_name, "runs": []})
    runs = records[-1]["runs"]
    if not runs:
        runs.append([position, 1, 1, [depth_value]])
        return
    run = runs[-1]
    start, step, count, depths = run
    if count == 1:
        next_step = position - start
        if next_step > 0:
            run[1] = next_step
            run[2] = 2
            depths.append(depth_value)
            return
        runs.append([position, 1, 1, [depth_value]])
        return
    if position == start + step * count:
        run[2] = count + 1
        depths.append(depth_value)
        return
    runs.append([position, 1, 1, [depth_value]])


def _guess_file_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return "text/tab-separated-values"
    if suffix in {".txt", ".gff", ".gff3", ".fa", ".fasta", ".fna", ".gb", ".gbk", ".gbff"}:
        return "text/plain"
    return "application/octet-stream"


def _json_clone(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value))
    except (TypeError, ValueError):
        return copy.deepcopy(value)


def _minimal_config_from_cli_args(context: SessionBuildContext) -> dict[str, Any]:
    args = list(context.cli_invocation_args)
    form: dict[str, Any] = {
        "prefix": context.output_prefix or _option_value(args, "-o", "--output") or "",
        "species": _option_value(args, "--species") or "",
        "strain": _option_value(args, "--strain") or "",
        "plot_title": _option_value(args, "--plot_title") or "",
        "separate_strands": "--separate_strands" in args,
    }
    adv: dict[str, Any] = {}
    losat: dict[str, Any] = {
        "outfmt": "6",
        "parallelWorkers": None,
        "executionMode": "auto",
        "totalThreadBudget": "safe",
        "threadsPerJob": "auto",
        "blastn": {"task": "megablast"},
        "blastp": {
            "mode": "orthogroup",
            "maxHits": 5,
            "candidateLimit": None,
            "orthogroupMembershipMode": "anchor_core_v1",
            "orthogroupMemberMaxHits": 5,
            "collinearMinAnchors": 1,
            "collinearMaxGeneGap": 0,
            "collinearMaxDiagonalDrift": 0,
            "collinearMaxConflictsInMergeGap": 1,
            "collinearMaxParalogLinksPerOrthogroup": 2,
            "collinearColorMode": "orientation",
            "collinearUnitMode": "auto",
            "collinearAnchorMode": "rbh",
            "collinearSearchScope": "adjacent",
        },
    }
    circular_conservation: dict[str, Any] | None = None
    if context.mode == "circular":
        form.update(
            {
                "track_type": _option_value(args, "--track_type") or "tuckin",
                "legend": _option_value(args, "-l", "--legend") or "left",
                "labels_mode": _optional_choice_option_value(
                    args,
                    "--labels",
                    choices=("none", "out", "both"),
                    default_when_present="out",
                    default_when_absent="none",
                ),
                "multi_record_canvas": "--multi_record_canvas" in args,
                "suppress_gc": "--suppress_gc" in args,
                "suppress_skew": "--suppress_skew" in args,
                "show_depth": "--show_depth" in args or "--depth" in args or "--depth_track" in args,
            }
        )
        adv["plot_title_position"] = _option_value(args, "--plot_title_position") or "none"
        circular_conservation = _populate_circular_cli_config(args, form, adv)
    else:
        form.update(
            {
                "legend": _option_value(args, "-l", "--legend") or "bottom",
                "scale_style": _option_value(args, "--scale_style") or "bar",
                "linear_track_layout": _option_value(args, "--track_layout") or "middle",
                "linear_ruler_on_axis": "--ruler_on_axis" in args,
                "align_center": "--align_center" in args,
                "keep_definition_left_aligned": "--keep_definition_left_aligned" in args,
                "show_gc": "--show_gc" in args,
                "show_skew": "--show_skew" in args,
                "show_depth": "--show_depth" in args or "--depth" in args or "--depth_track" in args,
                "normalize_length": "--normalize_length" in args,
                "show_labels_linear": _optional_choice_option_value(
                    args,
                    "--show_labels",
                    choices=("all", "first", "orthogroup_top", "none"),
                    default_when_present="all",
                    default_when_absent="none",
                ),
            }
        )
        adv["plot_title_position"] = _option_value(args, "--plot_title_position") or "bottom"
        _populate_linear_cli_config(args, form, adv, losat)
    features = _option_value(args, "-k", "--features")
    if features:
        adv["features"] = [item for item in features.split(",") if item]
    feature_shapes = _feature_shapes_from_cli_args(args)
    if feature_shapes:
        adv["feature_shapes"] = feature_shapes
    definition_line_styles = _definition_line_styles_from_cli_args(args)
    if definition_line_styles:
        adv["linear_definition_line_styles"] = definition_line_styles
    for key, option_names in {
        "nt": ("-n", "--nt"),
        "window_size": ("-w", "--window"),
        "step_size": ("-s", "--step"),
        "def_font_size": ("--definition_font_size",),
        "label_font_size": ("--label_font_size",),
        "plot_title_font_size": ("--plot_title_font_size",),
    }.items():
        value = _option_value(args, *option_names)
        if value is not None:
            adv[key] = value
    _populate_shared_cli_config(args, adv)
    filter_mode, blacklist_text = _label_filter_config_from_cli_args(args)
    blast_source = str(
        adv.get("blastSource")
        or ("upload" if _has_option(args, "-b", "--blast") else "losat")
    )
    losat_program = str(adv.get("losatProgram") or "blastn")
    return {
        "form": form,
        "adv": adv,
        "losat": losat,
        "cliOptions": _cli_options_payload(args),
        "colors": {},
        "palette": _option_value(args, "-p", "--palette") or "default",
        "rules": [],
        "qualifierPriorityRules": [],
        "filterMode": filter_mode,
        "whitelist": [],
        "blacklistText": blacklist_text,
        "blastSource": blast_source,
        "losatProgram": losat_program,
        "circularConservation": circular_conservation
        or {
            "enabled": False,
            "source": "losat",
            "losat_program": "blastn",
            "subject_gencode": 1,
            "reference": "auto",
            "labels": "",
            "series": [],
            "ring_width": None,
            "ring_gap": None,
        },
    }


def _populate_linear_session_fields_from_cli_context(
    payload: dict[str, Any],
    context: SessionBuildContext,
) -> None:
    args = tuple(str(arg) for arg in context.cli_invocation_args)
    _populate_linear_orthogroup_state_from_cli_args(payload, args)

    files = payload.get("files")
    if not isinstance(files, dict):
        return
    linear_seqs = files.get("linearSeqs")
    if not isinstance(linear_seqs, list) or not linear_seqs:
        return

    file_count = len(linear_seqs)
    selectors = _linear_input_selectors_from_cli_args(args, file_count)
    for index, selector in enumerate(selectors["record_ids"]):
        if not selector or index >= len(linear_seqs) or not isinstance(linear_seqs[index], dict):
            continue
        linear_seqs[index]["region_record_id"] = selector

    metadata = _normalize_linear_record_metadata(context.linear_record_metadata)
    labels = _linear_record_labels_from_cli_args(args, metadata)
    for source_index, label in labels.items():
        if 0 <= source_index < len(linear_seqs) and isinstance(linear_seqs[source_index], dict):
            linear_seqs[source_index]["definition"] = label

    subtitles = _linear_record_subtitles_from_cli_args(args, metadata)
    for source_index, subtitle in subtitles.items():
        if 0 <= source_index < len(linear_seqs) and isinstance(linear_seqs[source_index], dict):
            linear_seqs[source_index]["record_subtitle"] = subtitle

    region_fields = _linear_region_metadata_from_cli_args(args, metadata)
    region_sources = set(region_fields)
    for source_index, fields in region_fields.items():
        if 0 <= source_index < len(linear_seqs) and isinstance(linear_seqs[source_index], dict):
            linear_seqs[source_index].update(fields)

    for index, reverse_flag in enumerate(selectors["reverse_flags"]):
        if (
            reverse_flag
            and index not in region_sources
            and index < len(linear_seqs)
            and isinstance(linear_seqs[index], dict)
        ):
            linear_seqs[index]["region_reverse"] = True


def _populate_linear_orthogroup_state_from_cli_args(
    payload: dict[str, Any],
    args: Sequence[str],
) -> None:
    selected = _option_value(
        args,
        "--align_orthogroup_feature",
        "--align-orthogroup-feature",
    )
    if not selected:
        return
    protein_mode = _option_value(args, "--protein_blastp_mode", "--protein-blastp-mode")
    if protein_mode != "orthogroup":
        return
    orthogroup_state = payload.get("orthogroupState")
    if not isinstance(orthogroup_state, dict):
        orthogroup_state = {}
        payload["orthogroupState"] = orthogroup_state
    orthogroup_state["selectedOrthogroupAlignmentFeature"] = str(selected).strip()


def _normalize_linear_record_metadata(
    value: Sequence[Mapping[str, Any]] | None,
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    if not value:
        return normalized
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            continue
        try:
            loaded_index = int(item.get("loaded_index", index))
        except (TypeError, ValueError):
            loaded_index = index
        try:
            source_index = int(item.get("source_index"))
        except (TypeError, ValueError):
            source_index = -1
        try:
            source_loaded_count = int(item.get("source_loaded_count", 0))
        except (TypeError, ValueError):
            source_loaded_count = 0
        try:
            source_loaded_index = int(item.get("source_loaded_index", 0))
        except (TypeError, ValueError):
            source_loaded_index = 0
        normalized.append(
            {
                "loaded_index": loaded_index,
                "source_index": source_index,
                "source_loaded_index": source_loaded_index,
                "source_loaded_count": source_loaded_count,
                "record_id": str(item.get("record_id") or ""),
                "source_file": str(item.get("source_file") or ""),
                "source_basename": str(item.get("source_basename") or ""),
            }
        )
    return sorted(normalized, key=lambda item: int(item["loaded_index"]))


def _linear_record_labels_from_cli_args(
    args: Sequence[str],
    record_metadata: Sequence[Mapping[str, Any]],
) -> dict[int, str]:
    return _linear_record_text_values_from_cli_args(
        args,
        record_metadata,
        "--record_label",
        "--record-label",
    )


def _linear_record_subtitles_from_cli_args(
    args: Sequence[str],
    record_metadata: Sequence[Mapping[str, Any]],
) -> dict[int, str]:
    return _linear_record_text_values_from_cli_args(
        args,
        record_metadata,
        "--record_subtitle",
        "--record-subtitle",
    )


def _linear_record_text_values_from_cli_args(
    args: Sequence[str],
    record_metadata: Sequence[Mapping[str, Any]],
    *options: str,
) -> dict[int, str]:
    labels = _option_all_values(args, *options)
    mapped: dict[int, str] = {}
    if not labels or not record_metadata:
        return mapped
    for loaded_index, raw_label in enumerate(labels):
        if loaded_index >= len(record_metadata):
            break
        label = str(raw_label or "").strip()
        if not label:
            continue
        record = record_metadata[loaded_index]
        if int(record.get("source_loaded_count", 0) or 0) != 1:
            continue
        try:
            source_index = int(record.get("source_index", -1))
        except (TypeError, ValueError):
            source_index = -1
        if source_index >= 0:
            mapped[source_index] = label
    return mapped


def _linear_region_metadata_from_cli_args(
    args: Sequence[str],
    record_metadata: Sequence[Mapping[str, Any]],
) -> dict[int, dict[str, Any]]:
    raw_specs = _option_all_values(args, "--region")
    if not raw_specs or not record_metadata:
        return {}
    try:
        specs = parse_region_specs(raw_specs)
    except ValueError:
        return {}
    assignments = _assign_linear_region_specs_to_metadata(specs, record_metadata)
    reverse_flags = _linear_input_selectors_from_cli_args(
        args,
        _linear_file_count_from_metadata(record_metadata),
    )["reverse_flags"]
    mapped: dict[int, dict[str, Any]] = {}
    for loaded_index, spec in assignments.items():
        if loaded_index < 0 or loaded_index >= len(record_metadata):
            continue
        record = record_metadata[loaded_index]
        if int(record.get("source_loaded_count", 0) or 0) != 1:
            continue
        try:
            source_index = int(record.get("source_index", -1))
        except (TypeError, ValueError):
            source_index = -1
        if source_index < 0:
            continue
        if source_index < len(reverse_flags) and reverse_flags[source_index]:
            continue
        mapped[source_index] = {
            "region_start": int(spec.start),
            "region_end": int(spec.end),
            "region_reverse": bool(spec.reverse_complement),
        }
    return mapped


def _assign_linear_region_specs_to_metadata(
    specs: Sequence[RegionSpec],
    record_metadata: Sequence[Mapping[str, Any]],
) -> dict[int, RegionSpec]:
    total = len(record_metadata)
    if total == 0:
        return {}
    selectorless = [
        spec
        for spec in specs
        if spec.record_id is None and spec.record_index is None and spec.file_selector is None
    ]
    selectorful = [spec for spec in specs if spec not in selectorless]
    assignments: dict[int, RegionSpec] = {}
    if selectorless:
        if selectorful:
            return {}
        if total == 1 and len(specs) == 1:
            return {0: selectorless[0]}
        if len(specs) == total:
            return {index: spec for index, spec in enumerate(specs)}
        return {}

    id_to_indices: dict[str, list[int]] = {}
    file_aliases: list[set[str]] = []
    for index, record in enumerate(record_metadata):
        record_id = str(record.get("record_id") or "")
        id_to_indices.setdefault(record_id, []).append(index)
        aliases = _linear_metadata_file_aliases(record)
        file_aliases.append(aliases)

    def _find_file_indices(file_selector: str) -> list[int]:
        return [
            index
            for index, aliases in enumerate(file_aliases)
            if file_selector in aliases
        ]

    for spec in selectorful:
        file_indices: list[int] | None = None
        if spec.file_selector:
            file_indices = _find_file_indices(spec.file_selector)
            if not file_indices:
                return {}
        if spec.record_index is not None:
            if file_indices is None:
                target_index = spec.record_index
                if target_index < 0 or target_index >= total:
                    return {}
            else:
                if spec.record_index < 0 or spec.record_index >= len(file_indices):
                    return {}
                target_index = file_indices[spec.record_index]
            if target_index in assignments:
                return {}
            assignments[target_index] = spec
            continue

        record_id = spec.record_id or ""
        matches = (
            id_to_indices.get(record_id, [])
            if file_indices is None
            else [index for index in file_indices if str(record_metadata[index].get("record_id") or "") == record_id]
        )
        if len(matches) != 1:
            return {}
        target_index = matches[0]
        if target_index in assignments:
            return {}
        assignments[target_index] = spec
    return assignments


def _linear_metadata_file_aliases(record: Mapping[str, Any]) -> set[str]:
    aliases = {
        str(record.get("source_file") or ""),
        str(record.get("source_basename") or ""),
    }
    expanded: set[str] = set()
    for alias in aliases:
        if not alias:
            continue
        expanded.add(alias)
        try:
            expanded.add(PurePath(alias).name)
        except Exception:
            pass
        try:
            expanded.add(PureWindowsPath(alias).name)
        except Exception:
            pass
    return {alias for alias in expanded if alias}


def _linear_file_count_from_metadata(record_metadata: Sequence[Mapping[str, Any]]) -> int:
    max_source = -1
    for record in record_metadata:
        try:
            max_source = max(max_source, int(record.get("source_index", -1)))
        except (TypeError, ValueError):
            continue
    return max_source + 1


def _linear_input_selectors_from_cli_args(
    args: Sequence[str],
    file_count: int,
) -> dict[str, list[Any]]:
    record_ids = _option_all_values(args, "--record_id", "--record-id")[:file_count]
    reverse_values = _option_all_values(args, "--reverse_complement", "--reverse-complement")[:file_count]
    while len(record_ids) < file_count:
        record_ids.append("")
    while len(reverse_values) < file_count:
        reverse_values.append("")
    return {
        "record_ids": [_normalize_optional_cli_text(value) for value in record_ids],
        "reverse_flags": [_parse_cli_bool(value) for value in reverse_values],
    }


def _normalize_optional_cli_text(value: object) -> str:
    text = str(value or "").strip()
    if text.lower() in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
        return ""
    return text


def _parse_cli_bool(value: object) -> bool:
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"", "0", "false", "no", "n", "off", "none", "null", "-"}:
        return False
    return False


def _populate_shared_cli_config(args: Sequence[str], adv: dict[str, Any]) -> None:
    for key, option_names in {
        "block_stroke_width": ("--block_stroke_width", "--block-stroke-width"),
        "block_stroke_color": ("--block_stroke_color", "--block-stroke-color"),
        "line_stroke_width": ("--line_stroke_width", "--line-stroke-width"),
        "line_stroke_color": ("--line_stroke_color", "--line-stroke-color"),
        "axis_stroke_width": ("--axis_stroke_width", "--axis-stroke-width"),
        "axis_stroke_color": ("--axis_stroke_color", "--axis-stroke-color"),
        "legend_box_size": ("--legend_box_size", "--legend-box-size"),
        "legend_font_size": ("--legend_font_size", "--legend-font-size"),
        "scale_interval": ("--scale_interval", "--scale-interval"),
        "label_rendering": ("--label_rendering", "--label-rendering"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)
    if _has_option(args, "--resolve_overlaps", "--resolve-overlaps"):
        adv["resolve_overlaps"] = True
    _populate_gc_percent_cli_config(args, adv)
    _populate_depth_cli_config(args, adv)


def _populate_gc_percent_cli_config(args: Sequence[str], adv: dict[str, Any]) -> None:
    mode = _option_value(args, "--gc_content_mode", "--gc-content-mode")
    if mode is not None:
        adv["gc_content_mode"] = mode
    for key, option_names in {
        "gc_content_min_percent": ("--gc_content_min_percent", "--gc-content-min-percent"),
        "gc_content_max_percent": ("--gc_content_max_percent", "--gc-content-max-percent"),
        "gc_content_tick_interval": (
            "--gc_content_large_tick_interval",
            "--gc_content_tick_interval",
            "--gc-content-large-tick-interval",
            "--gc-content-tick-interval",
        ),
        "gc_content_small_tick_interval": (
            "--gc_content_small_tick_interval",
            "--gc-content-small-tick-interval",
        ),
        "gc_content_tick_font_size": ("--gc_content_tick_font_size", "--gc-content-tick-font-size"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)
    if _has_option(args, "--show_gc_content_axis", "--show-gc-content-axis"):
        adv["gc_content_show_axis"] = True
    if _has_option(args, "--hide_gc_content_axis", "--hide-gc-content-axis"):
        adv["gc_content_show_axis"] = False
    if _has_option(args, "--show_gc_content_ticks", "--show-gc-content-ticks"):
        adv["gc_content_show_ticks"] = True
    if _has_option(args, "--hide_gc_content_ticks", "--hide-gc-content-ticks"):
        adv["gc_content_show_ticks"] = False


def _populate_depth_cli_config(args: Sequence[str], adv: dict[str, Any]) -> None:
    for key, option_names in {
        "depth_color": ("--depth_color", "--depth-color"),
        "depth_height": ("--depth_height", "--depth-height"),
        "depth_width_circular": ("--depth_width", "--depth-width"),
        "depth_window_size": ("--depth_window", "--depth-window"),
        "depth_step_size": ("--depth_step", "--depth-step"),
        "depth_min": ("--depth_min", "--depth-min"),
        "depth_max": ("--depth_max", "--depth-max"),
        "depth_tick_interval": (
            "--depth_large_tick_interval",
            "--depth_tick_interval",
            "--depth-large-tick-interval",
            "--depth-tick-interval",
        ),
        "depth_small_tick_interval": ("--depth_small_tick_interval", "--depth-small-tick-interval"),
        "depth_tick_font_size": ("--depth_tick_font_size", "--depth-tick-font-size"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)
    if _has_option(args, "--share_depth_axis", "--share-depth-axis"):
        adv["depth_share_axis"] = True
    if _has_option(args, "--depth_log_scale", "--depth-log-scale"):
        adv["depth_normalize"] = True
    if _has_option(args, "--no_depth_log_scale", "--no-depth-log-scale"):
        adv["depth_normalize"] = False
    if _has_option(args, "--show_depth_axis", "--show-depth-axis"):
        adv["depth_show_axis"] = True
    if _has_option(args, "--hide_depth_axis", "--hide-depth-axis"):
        adv["depth_show_axis"] = False
    if _has_option(args, "--show_depth_ticks", "--show-depth-ticks"):
        adv["depth_show_ticks"] = True
    if _has_option(args, "--hide_depth_ticks", "--hide-depth-ticks"):
        adv["depth_show_ticks"] = False

    labels = _option_values(args, "--depth_track_label", "--depth-track-label")
    colors = _option_values(args, "--depth_track_color", "--depth-track-color")
    heights = _option_values(args, "--depth_track_height", "--depth-track-height")
    large_ticks = _option_values(
        args,
        "--depth_track_large_tick_interval",
        "--depth-track-large-tick-interval",
    )
    small_ticks = _option_values(
        args,
        "--depth_track_small_tick_interval",
        "--depth-track-small-tick-interval",
    )
    tick_fonts = _option_values(args, "--depth_track_tick_font_size", "--depth-track-tick-font-size")
    track_count = max(
        len(labels),
        len(colors),
        len(heights),
        len(large_ticks),
        len(small_ticks),
        len(tick_fonts),
        0,
    )
    if track_count:
        tracks: list[dict[str, Any]] = []
        for index in range(track_count):
            entry: dict[str, Any] = {}
            if index < len(labels):
                entry["label"] = labels[index]
            if index < len(colors):
                entry["color"] = colors[index]
            if index < len(heights):
                entry["height"] = heights[index]
            if index < len(large_ticks):
                entry["large_tick_interval"] = large_ticks[index]
            if index < len(small_ticks):
                entry["small_tick_interval"] = small_ticks[index]
            if index < len(tick_fonts):
                entry["tick_font_size"] = tick_fonts[index]
            tracks.append(entry)
        adv["depth_tracks"] = tracks


def _populate_circular_cli_config(
    args: Sequence[str],
    form: dict[str, Any],
    adv: dict[str, Any],
) -> dict[str, Any] | None:
    circular_conservation: dict[str, Any] | None = None
    if _has_option(args, "--conservation_blast", "--conservation-blast"):
        labels = _option_values(args, "--conservation_labels", "--conservation-labels")
        colors = _option_values(args, "--conservation_colors", "--conservation-colors")
        series: list[dict[str, Any]] = []
        for index in range(max(len(labels), len(colors))):
            entry: dict[str, Any] = {"sourceIndex": index}
            if index < len(labels):
                entry["label"] = labels[index]
            if index < len(colors):
                entry["color"] = colors[index]
            series.append(entry)
        adv["min_bitscore"] = _option_value(args, "--bitscore") or 50
        adv["evalue"] = _option_value(args, "--evalue") or "1e-5"
        adv["identity"] = _option_value(args, "--identity") or 70
        adv["alignment_length"] = _option_value(args, "--alignment_length", "--alignment-length") or 0
        circular_conservation = {
            "enabled": True,
            "source": "upload",
            "losat_program": "blastn",
            "subject_gencode": 1,
            "reference": _option_value(args, "--conservation_reference", "--conservation-reference") or "auto",
            "labels": "\n".join(labels),
            "series": series,
            "ring_width": _option_value(args, "--conservation_ring_width", "--conservation-ring-width"),
            "ring_gap": _option_value(args, "--conservation_ring_gap", "--conservation-ring-gap"),
        }
    for key, option_names in {
        "multi_record_size_mode": ("--multi_record_size_mode", "--multi-record-size-mode"),
        "multi_record_min_radius_ratio": ("--multi_record_min_radius_ratio", "--multi-record-min-radius-ratio"),
        "multi_record_column_gap_ratio": (
            "--multi_record_column_gap_ratio",
            "--multi-record-column-gap-ratio",
        ),
        "multi_record_row_gap_ratio": ("--multi_record_row_gap_ratio", "--multi-record-row-gap-ratio"),
        "tick_label_font_size": ("--tick_label_font_size", "--tick-label-font-size"),
        "circular_label_spacing": ("--circular_label_spacing", "--circular-label-spacing"),
        "circular_label_placement": ("--label_placement", "--label-placement"),
        "feature_width_circular": ("--feature_width", "--feature-width"),
        "gc_content_width_circular": ("--gc_content_width", "--gc-content-width"),
        "gc_content_radius_circular": ("--gc_content_radius", "--gc-content-radius"),
        "gc_skew_width_circular": ("--gc_skew_width", "--gc-skew-width"),
        "gc_skew_radius_circular": ("--gc_skew_radius", "--gc-skew-radius"),
        "center_reserved_radius": ("--center_reserved_radius", "--center-reserved-radius"),
        "outer_label_x_offset": ("--outer_label_x_radius_offset", "--outer-label-x-radius-offset"),
        "outer_label_y_offset": ("--outer_label_y_radius_offset", "--outer-label-y-radius-offset"),
        "inner_label_x_offset": ("--inner_label_x_radius_offset", "--inner-label-x-radius-offset"),
        "inner_label_y_offset": ("--inner_label_y_radius_offset", "--inner-label-y-radius-offset"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)
    if _has_option(args, "--keep_full_definition_with_plot_title", "--keep-full-definition-with-plot-title"):
        adv["keep_full_definition_with_plot_title"] = True
    circular_slots = _option_all_values(args, "--circular_track_slot", "--circular-track-slot")
    if circular_slots or _has_option(args, "--circular_track_order", "--circular-track-order"):
        adv["circular_track_slots_enabled"] = True
        adv["cli_circular_track_order"] = _option_value(args, "--circular_track_order", "--circular-track-order") or ""
        adv["cli_circular_track_slots"] = circular_slots
        _copy_option_value(args, adv, "circular_track_slots_axis_index", "--circular_track_axis_index", "--circular-track-axis-index")
    return circular_conservation


def _populate_linear_cli_config(
    args: Sequence[str],
    form: dict[str, Any],
    adv: dict[str, Any],
    losat: dict[str, Any],
) -> None:
    for key, option_names in {
        "feature_height": ("--feature_height", "--feature-height"),
        "gc_height": ("--gc_height", "--gc-height"),
        "comparison_height": ("--comparison_height", "--comparison-height"),
        "scale_font_size": ("--scale_font_size", "--scale-font-size"),
        "scale_stroke_width": ("--scale_stroke_width", "--scale-stroke-width"),
        "scale_stroke_color": ("--scale_stroke_color", "--scale-stroke-color"),
        "ruler_label_color": ("--ruler_label_color", "--ruler-label-color"),
        "pairwise_match_style": ("--pairwise_match_style", "--pairwise-match-style"),
        "track_axis_gap": ("--track_axis_gap", "--track-axis-gap"),
        "label_placement": ("--label_placement", "--label-placement"),
        "label_rotation": ("--label_rotation", "--label-rotation"),
        "linear_label_spacing": ("--linear_label_spacing", "--linear-label-spacing"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)
    if _has_option(args, "--show_replicon", "--show-replicon"):
        adv["linear_show_replicon"] = True
    if _has_option(args, "--hide_accession", "--hide-accession"):
        adv["linear_show_accession"] = False
    if _has_option(args, "--hide_length", "--hide-length"):
        adv["linear_show_length"] = False
    for key, option_names in {
        "min_bitscore": ("--bitscore",),
        "evalue": ("--evalue",),
        "identity": ("--identity",),
        "alignment_length": ("--alignment_length", "--alignment-length"),
    }.items():
        _copy_option_value(args, adv, key, *option_names)

    protein_mode = _option_value(args, "--protein_blastp_mode", "--protein-blastp-mode") or "none"
    if protein_mode != "none":
        adv["blastSource"] = "losat"
        adv["losatProgram"] = "blastp"
        losat["threadsPerJob"] = _option_value(args, "--losatp_threads", "--losatp-threads") or "auto"
        blastp = losat.setdefault("blastp", {})
        blastp["mode"] = protein_mode
        blastp["maxHits"] = _option_value(args, "--protein_blastp_max_hits", "--protein-blastp-max-hits") or 5
        candidate_limit = _option_value(
            args,
            "--protein_blastp_candidate_limit",
            "--protein-blastp-candidate-limit",
        )
        blastp["candidateLimit"] = candidate_limit if candidate_limit not in (None, "none") else None
        blastp["collinearMinAnchors"] = _option_value(args, "--collinear_min_anchors", "--collinear-min-anchors") or 1
        blastp["collinearMaxGeneGap"] = (
            _option_value(
                args,
                "--collinear_max_unit_gap",
                "--collinear_max_gene_gap",
                "--collinear-max-unit-gap",
                "--collinear-max-gene-gap",
            )
            or 0
        )
        blastp["collinearMaxDiagonalDrift"] = _option_value(
            args,
            "--collinear_max_diagonal_drift",
            "--collinear-max-diagonal-drift",
        ) or 0
        blastp["collinearMaxConflictsInMergeGap"] = _option_value(
            args,
            "--collinear_max_conflicts_in_merge_gap",
            "--collinear-max-conflicts-in-merge-gap",
        ) or 1
        blastp["collinearMaxParalogLinksPerOrthogroup"] = _option_value(
            args,
            "--collinear_max_paralog_links_per_orthogroup",
            "--collinear-max-paralog-links-per-orthogroup",
        ) or 2
        blastp["collinearColorMode"] = _option_value(args, "--collinear_color_mode", "--collinear-color-mode") or "orientation"
        blastp["collinearUnitMode"] = _option_value(args, "--collinear_unit_mode", "--collinear-unit-mode") or "auto"
        blastp["collinearSearchScope"] = _option_value(args, "--collinear_search_scope", "--collinear-search-scope") or "adjacent"
    elif _has_option(args, "-b", "--blast"):
        adv["blastSource"] = "upload"
    else:
        adv["blastSource"] = "losat"
        adv["losatProgram"] = "blastn"

    linear_slots = _option_all_values(args, "--linear_track_slot", "--linear-track-slot")
    if linear_slots or _has_option(args, "--linear_track_order", "--linear-track-order"):
        adv["linear_track_slots_enabled"] = True
        adv["cli_linear_track_order"] = _option_value(args, "--linear_track_order", "--linear-track-order") or ""
        adv["cli_linear_track_slots"] = linear_slots
        _copy_option_value(args, adv, "linear_track_slots_axis_index", "--linear_track_axis_index", "--linear-track-axis-index")


def _copy_option_value(
    args: Sequence[str],
    target: dict[str, Any],
    key: str,
    *names: str,
) -> None:
    value = _option_value(args, *names)
    if value is not None:
        target[key] = value


def _cli_options_payload(args: Sequence[str]) -> dict[str, Any]:
    raw_args = [str(arg) for arg in args]
    entries: list[dict[str, Any]] = []
    by_key: dict[str, list[Any]] = {}
    positionals: list[str] = []

    index = 0
    while index < len(raw_args):
        token = raw_args[index]
        if not _looks_like_cli_option(token):
            positionals.append(token)
            index += 1
            continue

        option = token
        values: list[str] = []
        if token.startswith("--") and "=" in token:
            option, value = token.split("=", 1)
            values.append(value)
            index += 1
        else:
            index += 1
            while index < len(raw_args) and not _looks_like_cli_option(raw_args[index]):
                values.append(raw_args[index])
                index += 1

        key = _cli_option_key(option)
        entry = {
            "option": option,
            "key": key,
            "values": values,
        }
        if not values:
            entry["flag"] = True
        entries.append(entry)
        by_key.setdefault(key, []).append(True if not values else (values[0] if len(values) == 1 else values))

    return {
        "schema": 1,
        "mode": "argv",
        "rawArgs": raw_args,
        "options": entries,
        "byKey": by_key,
        "positionals": positionals,
    }


def _looks_like_cli_option(token: object) -> bool:
    text = str(token)
    if text == "-":
        return False
    if not text.startswith("-"):
        return False
    if re.match(r"^-\d", text):
        return False
    return True


def _cli_option_key(option: str) -> str:
    aliases = {
        "-b": "blast",
        "-d": "default_colors",
        "-f": "format",
        "-k": "features",
        "-l": "legend",
        "-n": "nt",
        "-o": "output",
        "-p": "palette",
        "-s": "step",
        "-t": "table",
        "-w": "window",
    }
    if option in aliases:
        return aliases[option]
    return option.lstrip("-").replace("-", "_")


def migrate_legacy_repeat_feature_shape_args(
    args: Sequence[str],
    *,
    session_version: int,
) -> list[str]:
    """Preserve the old repeat rectangle for non-canonical v27-30 replay."""

    migrated = [str(arg) for arg in args]
    if int(session_version) > 30:
        return migrated
    features_raw = _option_value(migrated, "-k", "--features")
    effective_features = (
        {item.strip() for item in features_raw.split(",") if item.strip()}
        if features_raw is not None
        else {
            "CDS",
            "rRNA",
            "tRNA",
            "tmRNA",
            "ncRNA",
            "misc_RNA",
            "repeat_region",
        }
    )
    if (
        "repeat_region" in effective_features
        and "repeat_region" not in _feature_shapes_from_cli_args(migrated)
    ):
        insertion_index = next(
            (
                index
                for index, token in enumerate(migrated)
                if token in {"-f", "--format"} or token.startswith("--format=")
            ),
            len(migrated),
        )
        migrated[insertion_index:insertion_index] = [
            "--feature_shape",
            "repeat_region=rectangle",
        ]
    return migrated


def _feature_shapes_from_cli_args(args: Sequence[str]) -> dict[str, str]:
    shapes: dict[str, str] = {}
    for assignment in _option_all_values(args, "--feature_shape", "--feature-shape"):
        feature_type, separator, shape = str(assignment).partition("=")
        feature_type = feature_type.strip()
        shape = shape.strip().lower()
        if separator and feature_type and shape in {"arrow", "rectangle", "underlay"}:
            shapes[feature_type] = shape
    return shapes


def _definition_line_styles_from_cli_args(args: Sequence[str]) -> dict[str, dict[str, object]]:
    assignments = _option_all_values(args, "--definition_line_style", "--definition-line-style")
    if not assignments:
        return {}
    try:
        return parse_definition_line_style_overrides(assignments)
    except ValueError:
        return {}


def _label_filter_config_from_cli_args(args: Sequence[str]) -> tuple[str, str]:
    blacklist = _option_value(args, "--label_blacklist", "--label-blacklist")
    if blacklist is not None:
        try:
            is_file = Path(blacklist).is_file()
        except OSError:
            is_file = False
        return "Blacklist", "" if is_file else blacklist
    if _has_option(args, "--label_whitelist", "--label-whitelist"):
        return "Whitelist", ""
    return "None", ""


def _has_option(args: Sequence[str], *names: str) -> bool:
    for token in args:
        text = str(token)
        if any(text == name or text.startswith(f"{name}=") for name in names):
            return True
    return False


def _optional_choice_option_value(
    args: Sequence[str],
    name: str,
    *,
    choices: Sequence[str],
    default_when_present: str,
    default_when_absent: str,
) -> str:
    if not _has_option(args, name):
        return default_when_absent
    value = _option_value(args, name)
    return value if value in choices else default_when_present


def _option_all_values(args: Sequence[str], *names: str) -> list[str]:
    values: list[str] = []
    for index, token in enumerate(args):
        text = str(token)
        for name in names:
            if text == name:
                if index + 1 < len(args):
                    values.append(str(args[index + 1]))
                break
            prefix = f"{name}="
            if text.startswith(prefix):
                values.append(text[len(prefix):])
                break
    return values


def _option_values(args: Sequence[str], *names: str) -> list[str]:
    values: list[str] = []
    index = 0
    while index < len(args):
        text = str(args[index])
        matched_name = next((name for name in names if text == name or text.startswith(f"{name}=")), None)
        if matched_name is None:
            index += 1
            continue
        prefix = f"{matched_name}="
        if text.startswith(prefix):
            values.append(text[len(prefix):])
            index += 1
            continue
        index += 1
        while index < len(args) and not str(args[index]).startswith("-"):
            values.append(str(args[index]))
            index += 1
    return values


def _update_config_prefix(config: dict[str, Any], output_prefix: str | None) -> None:
    if not output_prefix:
        return
    form = config.setdefault("form", {})
    if isinstance(form, dict):
        form["prefix"] = output_prefix


def _input_type_from_args(args: Sequence[str]) -> str:
    return "gff" if "--gff" in args else "gb"


def _option_value(args: Sequence[str], *names: str) -> str | None:
    for index, token in enumerate(args):
        text = str(token)
        for name in names:
            prefix = f"{name}="
            if text.startswith(prefix):
                return text[len(prefix):]
            if text == name and index + 1 < len(args):
                return str(args[index + 1])
    return None


def _string_or_none(value: object) -> str | None:
    text = str(value or "").strip()
    return text or None


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


__all__ = [
    "CURRENT_SESSION_VERSION",
    "CANONICAL_SESSION_MIN_VERSION",
    "DEPTH_FILE_ENCODING",
    "DEPTH_FILE_SCHEMA",
    "LEGACY_LOSAT_DERIVED_CACHE_SCHEMA",
    "LEGACY_PROTEIN_CANDIDATE_SCHEMA",
    "LOSAT_DERIVED_CACHE_SCHEMA",
    "NUCLEOTIDE_LOSAT_CACHE_SCHEMA",
    "PROTEIN_IDENTITY_MANIFEST_SCHEMA",
    "PROTEIN_LOSAT_CACHE_SCHEMA",
    "SESSION_FORMAT",
    "SUPPORTED_SESSION_VERSIONS",
    "SessionBuildContext",
    "SessionFileBinding",
    "SessionRunSpec",
    "build_session_json",
    "classify_raw_losat_cache_entry",
    "decode_depth_payload",
    "encode_depth_text",
    "empty_protein_identity_manifest",
    "get_session_slot",
    "load_session",
    "materialize_embedded_file",
    "migrate_legacy_repeat_feature_shape_args",
    "normalize_current_session_artifacts",
    "safe_embedded_filename",
    "serialize_file_entry",
    "session_mode",
    "session_to_cli_args",
    "validate_session",
    "validate_current_session_artifacts",
    "write_session_json",
]
