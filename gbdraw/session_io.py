#!/usr/bin/env python
# coding: utf-8

"""GUI session JSON loading, validation, materialization, and sidecar building."""

from __future__ import annotations

import base64
import binascii
import copy
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence

from .analysis.protein_artifacts import (
    CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
    is_current_derived_protein_artifact,
    validate_current_derived_protein_artifacts,
)
from .definition_line_styles import DEFINITION_LINE_KINDS
from .exceptions import GbdrawError, ValidationError
from .render.formats import normalize_format_token
from .render.output_paths import commit_staged_output_file

if TYPE_CHECKING:
    from .analysis.protein_colinearity import ProteinIdentityManifest
    from .api.requests import DiagramRequest

SESSION_FORMAT = "gbdraw-session"
CURRENT_SESSION_VERSION = 41
CANONICAL_SESSION_MIN_VERSION = 31
SUPPORTED_SESSION_VERSIONS = frozenset(
    {27, 28, 29, 30, 31, 32, 33, 39, 40, CURRENT_SESSION_VERSION}
)
CANONICAL_REQUEST_SCHEMAS_BY_SESSION_VERSION = {
    31: frozenset({1}),
    32: frozenset({1, 2}),
    33: frozenset({2}),
    39: frozenset({5}),
    40: frozenset({5}),
    CURRENT_SESSION_VERSION: frozenset({6}),
}
CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39
PROTEIN_LOSAT_CACHE_SCHEMA = 4
NUCLEOTIDE_LOSAT_CACHE_SCHEMA = 2
LOSAT_DERIVED_CACHE_SCHEMA = CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA
LEGACY_LOSAT_DERIVED_CACHE_SCHEMA = 1
PROTEIN_IDENTITY_MANIFEST_SCHEMA = 2
LEGACY_PROTEIN_CANDIDATE_SCHEMA = 1
FEATURE_CATALOG_SCHEMA = 1
FEATURE_CATALOG_ENCODING = "biological-authority-v1"
CURRENT_FEATURE_CATALOG_SCHEMA = 3
CURRENT_SESSION_TOP_LEVEL_FIELDS = frozenset(
    {
        "format",
        "version",
        "createdAt",
        "title",
        "renderRequest",
        "resources",
        "webFiles",
        "config",
        "ui",
        "files",
        "results",
        "features",
        "editorState",
        "orthogroupState",
        "losatCache",
        "losatDerivedCache",
        "proteinIdentityManifest",
        "legacyArtifacts",
        "runMetadata",
        "cliInvocation",
    }
)
CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS = frozenset(
    {
        "extractedFeatures",
        "biologicalFeatures",
        "featureSelectorSafetyScope",
        "featureRecordIds",
        "featureCatalog",
    }
)
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


def _feature_catalog_key(feature: Mapping[str, Any]) -> tuple[int, str] | None:
    record_index = feature.get("record_idx")
    stable_id = (
        feature.get("stable_svg_id")
        or feature.get("stable_feature_id")
        or feature.get("svg_id")
    )
    if (
        not isinstance(record_index, int)
        or isinstance(record_index, bool)
        or not isinstance(stable_id, str)
        or not stable_id
    ):
        return None
    return record_index, stable_id


def _json_values_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, Mapping):
        return (
            left.keys() == right.keys()
            and all(_json_values_equal(left[key], right[key]) for key in left)
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _json_values_equal(left_value, right_value)
            for left_value, right_value in zip(left, right, strict=True)
        )
    return bool(left == right)


def _first_feature_qualifier(
    qualifiers: Mapping[str, Any],
    key: str,
) -> str:
    values = qualifiers.get(key)
    if not isinstance(values, list):
        return ""
    return next((value for value in values if value != ""), "")


def _feature_qualifiers_are_strings(value: Any) -> bool:
    return isinstance(value, Mapping) and all(
        isinstance(key, str)
        and isinstance(items, list)
        and all(isinstance(item, str) for item in items)
        for key, items in value.items()
    )


def _expand_compact_biological_feature(
    feature: Mapping[str, Any],
    *,
    index: int,
    profile: str,
) -> dict[str, Any]:
    expanded = copy.deepcopy(dict(feature))
    qualifiers = expanded.get("qualifiers")
    selector = expanded.get("selector")
    if not isinstance(qualifiers, Mapping) or not isinstance(selector, Mapping):
        raise ValidationError(
            "Compact session biological features require qualifier and selector objects."
        )
    selector_qualifiers = selector.get("qualifiers")
    if not _feature_qualifiers_are_strings(qualifiers) or (
        "qualifiers" in selector
        and not _feature_qualifiers_are_strings(selector_qualifiers)
    ):
        raise ValidationError(
            "Compact session feature qualifiers must contain string arrays."
        )
    normalized_qualifiers = dict(qualifiers)
    normalized_selector = dict(selector)
    svg_id = expanded.get("svg_id")
    if not isinstance(svg_id, str) or not svg_id:
        raise ValidationError("Compact session biological features require svg_id.")

    expanded.setdefault("id", f"f{index}")
    expanded.setdefault("stable_svg_id", svg_id)
    expanded.setdefault("stable_feature_id", svg_id)
    for field in (
        "protein_id",
        "locus_tag",
        "gene_id",
        "old_locus_tag",
        "gene",
        "product",
    ):
        expanded.setdefault(
            field,
            _first_feature_qualifier(normalized_qualifiers, field),
        )
    expanded.setdefault("source_protein_id", expanded.get("protein_id", ""))
    expanded.setdefault(
        "note",
        _first_feature_qualifier(normalized_qualifiers, "note")[:50],
    )
    expanded.setdefault("sequence_warnings", [])
    normalized_selector.setdefault(
        "qualifiers",
        copy.deepcopy(normalized_qualifiers),
    )
    normalized_selector.setdefault("hash", svg_id)
    expanded["selector"] = normalized_selector

    if profile == "rich-v1":
        translation = _first_feature_qualifier(
            normalized_qualifiers,
            "translation",
        )
        expanded.setdefault("amino_acid_sequence", translation)
        if not isinstance(expanded.get("nucleotide_sequence"), str):
            raise ValidationError(
                "Compact rich biological features require nucleotide sequences."
            )
    elif profile != "sanitized-v1":
        raise ValidationError(
            f"Unsupported compact feature catalog profile: {profile!r}."
        )
    return expanded


def _compact_biological_feature(
    feature: Mapping[str, Any],
    *,
    index: int,
    profile: str,
) -> dict[str, Any] | None:
    compact = copy.deepcopy(dict(feature))
    qualifiers = compact.get("qualifiers")
    selector = compact.get("selector")
    if not isinstance(qualifiers, Mapping) or not isinstance(selector, Mapping):
        return None
    selector_qualifiers = selector.get("qualifiers")
    if not _feature_qualifiers_are_strings(qualifiers) or (
        "qualifiers" in selector
        and not _feature_qualifiers_are_strings(selector_qualifiers)
    ):
        return None
    svg_id = compact.get("svg_id")
    if not isinstance(svg_id, str) or not svg_id:
        return None

    if compact.get("id") == f"f{index}":
        compact.pop("id")
    if compact.get("stable_svg_id") == svg_id:
        compact.pop("stable_svg_id")
    if compact.get("stable_feature_id") == svg_id:
        compact.pop("stable_feature_id")
    if (
        "source_protein_id" in compact
        and compact.get("source_protein_id") == compact.get("protein_id")
    ):
        compact.pop("source_protein_id")
    if compact.get("sequence_warnings") == []:
        compact.pop("sequence_warnings")

    compact_selector = dict(selector)
    if _json_values_equal(compact_selector.get("qualifiers"), qualifiers):
        compact_selector.pop("qualifiers")
    if compact_selector.get("hash") == svg_id:
        compact_selector.pop("hash")
    compact["selector"] = compact_selector

    for field in (
        "protein_id",
        "locus_tag",
        "gene_id",
        "old_locus_tag",
        "gene",
        "product",
    ):
        if compact.get(field) == _first_feature_qualifier(qualifiers, field):
            compact.pop(field)
    if compact.get("note") == _first_feature_qualifier(qualifiers, "note")[:50]:
        compact.pop("note")
    if (
        profile == "rich-v1"
        and compact.get("amino_acid_sequence")
        == _first_feature_qualifier(qualifiers, "translation")
    ):
        compact.pop("amino_acid_sequence")

    try:
        expanded = _expand_compact_biological_feature(
            compact,
            index=index,
            profile=profile,
        )
    except ValidationError:
        return None
    return compact if _json_values_equal(expanded, dict(feature)) else None


def _compact_feature_catalog(features: Mapping[str, Any]) -> Mapping[str, Any]:
    if "featureCatalog" in features:
        return features
    extracted = features.get("extractedFeatures")
    biological = features.get("biologicalFeatures")
    if (
        not isinstance(extracted, list)
        or not extracted
        or not isinstance(biological, list)
        or not biological
        or not all(isinstance(feature, Mapping) for feature in (*extracted, *biological))
    ):
        return features

    has_rich_sequences = [
        "nucleotide_sequence" in feature and "amino_acid_sequence" in feature
        for feature in biological
    ]
    if all(has_rich_sequences):
        profile = "rich-v1"
    elif not any(
        "nucleotide_sequence" in feature or "amino_acid_sequence" in feature
        for feature in biological
    ):
        profile = "sanitized-v1"
    else:
        return features

    biological_indexes: dict[tuple[int, str], int] = {}
    for index, feature in enumerate(biological):
        key = _feature_catalog_key(feature)
        if key is None or key in biological_indexes:
            return features
        biological_indexes[key] = index

    compact_biological: list[dict[str, Any]] = []
    for index, feature in enumerate(biological):
        compact = _compact_biological_feature(
            feature,
            index=index,
            profile=profile,
        )
        if compact is None:
            return features
        compact_biological.append(compact)

    references: list[list[Any]] = []
    referenced_biological_indexes: set[int] = set()
    for feature in extracted:
        key = _feature_catalog_key(feature)
        biological_index = biological_indexes.get(key) if key is not None else None
        feature_id = feature.get("id")
        if biological_index is None or not isinstance(feature_id, str):
            return features
        if biological_index in referenced_biological_indexes:
            return features
        referenced_biological_indexes.add(biological_index)
        projected = copy.deepcopy(dict(biological[biological_index]))
        projected.pop("feature_index", None)
        projected["id"] = feature_id
        rendered_id = feature.get("rendered_feature_svg_id")
        reference: list[Any] = [biological_index, feature_id]
        if isinstance(rendered_id, str):
            projected["rendered_feature_svg_id"] = rendered_id
            reference.append(rendered_id)
        if not _json_values_equal(projected, dict(feature)):
            return features
        references.append(reference)

    compact_features = dict(features)
    compact_features.pop("extractedFeatures", None)
    compact_features["biologicalFeatures"] = compact_biological
    compact_features["featureCatalog"] = {
        "schema": FEATURE_CATALOG_SCHEMA,
        "encoding": FEATURE_CATALOG_ENCODING,
        "profile": profile,
        "extracted": references,
    }
    return compact_features


def compact_session_feature_catalog(
    session: Mapping[str, Any],
) -> Mapping[str, Any]:
    """Compact the released v39 feature representation for compatibility."""

    if session.get("version") != 39:
        return session
    features = session.get("features")
    if not isinstance(features, Mapping):
        return session
    compact_features = _compact_feature_catalog(features)
    if compact_features is features:
        return session
    compact_session = dict(session)
    compact_session["features"] = compact_features
    return compact_session


def expand_session_feature_catalog(
    session: Mapping[str, Any],
) -> dict[str, Any]:
    """Expand the released v39 compact feature representation."""

    expanded_session = dict(session)
    if expanded_session.get("version") == CURRENT_SESSION_VERSION:
        return expanded_session
    features = expanded_session.get("features")
    if not isinstance(features, Mapping):
        return expanded_session
    if "featureCatalog" not in features:
        return expanded_session
    catalog = features.get("featureCatalog")
    if (
        not isinstance(catalog, Mapping)
        or type(catalog.get("schema")) is not int
        or catalog.get("schema") != FEATURE_CATALOG_SCHEMA
        or catalog.get("encoding") != FEATURE_CATALOG_ENCODING
        or catalog.get("profile") not in {"rich-v1", "sanitized-v1"}
        or not isinstance(catalog.get("extracted"), list)
        or "extractedFeatures" in features
    ):
        raise ValidationError("Invalid compact session feature catalog.")
    biological = features.get("biologicalFeatures")
    if not isinstance(biological, list) or not all(
        isinstance(feature, Mapping) for feature in biological
    ):
        raise ValidationError(
            "Compact session feature catalog requires biologicalFeatures."
        )

    references = catalog["extracted"]
    if len(references) > len(biological):
        raise ValidationError("Invalid compact extracted-feature reference.")
    validated_references: list[tuple[int, str, str | None]] = []
    referenced_biological_indexes: set[int] = set()
    for reference in references:
        if (
            not isinstance(reference, list)
            or len(reference) not in {2, 3}
            or not isinstance(reference[0], int)
            or isinstance(reference[0], bool)
            or not 0 <= reference[0] < len(biological)
            or not isinstance(reference[1], str)
            or (len(reference) == 3 and not isinstance(reference[2], str))
        ):
            raise ValidationError("Invalid compact extracted-feature reference.")
        biological_index = reference[0]
        if biological_index in referenced_biological_indexes:
            raise ValidationError("Invalid compact extracted-feature reference.")
        referenced_biological_indexes.add(biological_index)
        validated_references.append(
            (
                biological_index,
                reference[1],
                reference[2] if len(reference) == 3 else None,
            )
        )

    profile = str(catalog["profile"])
    expanded_biological = [
        _expand_compact_biological_feature(
            feature,
            index=index,
            profile=profile,
        )
        for index, feature in enumerate(biological)
    ]
    expanded_extracted: list[dict[str, Any]] = []
    for biological_index, feature_id, rendered_id in validated_references:
        projected = copy.deepcopy(expanded_biological[biological_index])
        projected.pop("feature_index", None)
        projected["id"] = feature_id
        if rendered_id is not None:
            projected["rendered_feature_svg_id"] = rendered_id
        expanded_extracted.append(projected)

    expanded_features = dict(features)
    expanded_features.pop("featureCatalog", None)
    expanded_features["biologicalFeatures"] = expanded_biological
    expanded_features["extractedFeatures"] = expanded_extracted
    expanded_session["features"] = expanded_features
    return expanded_session


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
    payload = expand_session_feature_catalog(payload)
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
        request_schema = render_request.get("schema")
        if not isinstance(request_schema, int) or isinstance(request_schema, bool):
            raise ValidationError("renderRequest.schema must be an integer.")
        from .session_request_codec import SUPPORTED_CANONICAL_REQUEST_SCHEMAS

        if request_schema not in SUPPORTED_CANONICAL_REQUEST_SCHEMAS:
            raise ValidationError(
                f"Unsupported canonical renderRequest schema: {request_schema}."
            )
        compatible_schemas = CANONICAL_REQUEST_SCHEMAS_BY_SESSION_VERSION.get(version)
        if compatible_schemas is None or request_schema not in compatible_schemas:
            expected = ", ".join(str(schema) for schema in sorted(compatible_schemas or ()))
            raise ValidationError(
                f"Session version {version} is incompatible with canonical "
                f"renderRequest schema {request_schema}; expected schema {expected}."
            )
        if not isinstance(resources, Mapping):
            raise ValidationError(
                f"Session version {version} requires a canonical resources object."
            )
        files = session.get("files")
        if version == CURRENT_SESSION_VERSION and "files" in session:
            raise ValidationError(
                f"Session version {version} cannot contain legacy files; "
                "use resources and webFiles."
            )
        if files is not None and not isinstance(files, Mapping):
            raise ValidationError("Session files must be an object when present.")
    else:
        if "renderRequest" in session:
            raise ValidationError(
                f"Session version {version} predates canonical renderRequest schemas."
            )
        files = session.get("files")
        if files is None or not isinstance(files, Mapping):
            raise ValidationError("Session files are required for CLI regeneration.")
    if version >= CURRENT_ARTIFACT_SESSION_MIN_VERSION:
        validate_current_session_artifacts(session)
    if version == CURRENT_SESSION_VERSION:
        _validate_current_comparison_authority(session)
        _validate_current_feature_catalog_authority(session)


def _validate_current_comparison_authority(
    session: Mapping[str, Any],
) -> None:
    """Reject retired comparison fields and validate an optional Web draft."""

    config_value = session.get("config")
    config = config_value if isinstance(config_value, Mapping) else {}
    ui_value = session.get("ui")
    ui = ui_value if isinstance(ui_value, Mapping) else {}
    web_files_value = session.get("webFiles")
    web_files = web_files_value if isinstance(web_files_value, Mapping) else {}
    has_bindings = "bindings" in web_files
    bindings_value = web_files.get("bindings")
    if has_bindings and not isinstance(bindings_value, Mapping):
        raise ValidationError("Session webFiles.bindings must be an object.")
    bindings = bindings_value if isinstance(bindings_value, Mapping) else {}
    if has_bindings and bindings.get("schema") != 1:
        raise ValidationError("Unsupported Web file binding schema.")
    adv_value = config.get("adv")
    adv = adv_value if isinstance(adv_value, Mapping) else {}
    layout_value = config.get("linearRecordLayout")
    layout = layout_value if isinstance(layout_value, Mapping) else None

    if "blastSource" in config or "blastSource" in adv:
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} cannot contain retired "
            "blastSource state."
        )
    if "blastSource" in ui:
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} cannot contain "
            "ui.blastSource."
        )
    if layout is not None and "comparisons" in layout:
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} cannot contain "
            "config.linearRecordLayout.comparisons."
        )

    linear_sequences = bindings.get("linearSeqs")
    if isinstance(linear_sequences, list):
        for sequence in linear_sequences:
            if not isinstance(sequence, Mapping):
                continue
            if "blast" in sequence:
                raise ValidationError(
                    f"Session version {CURRENT_SESSION_VERSION} cannot contain "
                    "per-record BLAST bindings."
                )
            if "losat_filename" in sequence:
                raise ValidationError(
                    f"Session version {CURRENT_SESSION_VERSION} cannot contain "
                    "per-record LOSAT filenames."
                )
    linear_metadata = web_files.get("linearRecordMetadata")
    if isinstance(linear_metadata, list):
        for metadata in linear_metadata:
            if not isinstance(metadata, Mapping):
                continue
            if "losatFilename" in metadata or "losat_filename" in metadata:
                raise ValidationError(
                    f"Session version {CURRENT_SESSION_VERSION} cannot contain "
                    "per-record LOSAT filenames."
                )
    if "linearCanonicalComparisons" in bindings:
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} cannot bind comparison "
            "artifacts outside "
            "the committed request."
        )

    comparison_bindings_value = bindings.get("linearComparisons")
    if (
        "linearComparisons" in bindings
        and not isinstance(comparison_bindings_value, list)
    ):
        raise ValidationError("Current comparison file bindings must be an array.")
    comparison_bindings = (
        comparison_bindings_value
        if isinstance(comparison_bindings_value, list)
        else []
    )
    has_web_draft = (
        "linearRecordLayout" in config or "linearComparisonPlan" in config
    )
    if not has_web_draft and not comparison_bindings:
        return

    plan = config.get("linearComparisonPlan")
    if not isinstance(plan, Mapping):
        raise ValidationError(
            "Current Web comparison draft requires config.linearComparisonPlan."
        )
    if plan.get("mode") not in {"none", "adjacent", "selected"}:
        raise ValidationError("config.linearComparisonPlan.mode is invalid.")
    if plan.get("defaultSource") not in {"losat", "upload"}:
        raise ValidationError(
            "config.linearComparisonPlan.defaultSource is invalid."
        )
    edges = plan.get("edges")
    if not isinstance(edges, list):
        raise ValidationError("config.linearComparisonPlan.edges must be an array.")
    allowed_edge_fields = {
        "id",
        "queryUid",
        "subjectUid",
        "included",
        "fileActive",
        "losatFilenameActive",
        "source",
        "losatFilename",
    }
    edge_ids: set[str] = set()
    for edge in edges:
        if not isinstance(edge, Mapping):
            raise ValidationError(
                "Each config.linearComparisonPlan edge must be an object."
            )
        unknown = sorted(str(field) for field in edge if field not in allowed_edge_fields)
        if unknown:
            raise ValidationError(
                "config.linearComparisonPlan edge contains retired or unknown "
                f"field(s): {', '.join(unknown)}."
            )
        edge_id = str(edge.get("id") or "").strip()
        query_uid = str(edge.get("queryUid") or "").strip()
        subject_uid = str(edge.get("subjectUid") or "").strip()
        if not edge_id or not query_uid or not subject_uid:
            raise ValidationError(
                "Current comparison-plan edges require stable IDs and endpoint UIDs."
            )
        if edge_id in edge_ids:
            raise ValidationError(
                f"Current comparison-plan edge ID is duplicated: {edge_id}."
            )
        edge_ids.add(edge_id)
        if (
            type(edge.get("included")) is not bool
            or type(edge.get("fileActive")) is not bool
            or type(edge.get("losatFilenameActive")) is not bool
            or edge.get("source") not in {"losat", "upload"}
            or not isinstance(edge.get("losatFilename"), str)
        ):
            raise ValidationError(
                "Current comparison-plan edge metadata is invalid."
            )

    bound_ids: set[str] = set()
    resources_value = session.get("resources")
    resources = resources_value if isinstance(resources_value, Mapping) else {}
    for binding in comparison_bindings:
        if not isinstance(binding, Mapping):
            raise ValidationError(
                "Each current comparison file binding must be an object."
            )
        unknown = sorted(
            str(field) for field in binding if field not in {"id", "file"}
        )
        if unknown:
            raise ValidationError(
                "Current comparison file binding duplicates plan metadata: "
                + ", ".join(unknown)
                + "."
            )
        edge_id = str(binding.get("id") or "").strip()
        if not edge_id or edge_id not in edge_ids:
            raise ValidationError(
                "Current comparison file binding must reference a plan edge ID."
            )
        if edge_id in bound_ids:
            raise ValidationError(
                f"Current comparison file binding is duplicated: {edge_id}."
            )
        file_binding = binding.get("file")
        if not isinstance(file_binding, Mapping):
            raise ValidationError(
                "Each current comparison file binding requires a file resource binding."
            )
        resource_id = str(file_binding.get("resourceId") or "").strip()
        if not resource_id:
            raise ValidationError(
                "Each current comparison file binding requires a resourceId."
            )
        if resource_id not in resources:
            raise ValidationError(
                "Current comparison file binding references a missing resource: "
                f"{resource_id}."
            )
        bound_ids.add(edge_id)
    for edge in edges:
        if edge.get("fileActive") and str(edge.get("id") or "") not in bound_ids:
            raise ValidationError(
                "Active comparison file is missing its Web file binding: "
                f"{edge.get('id')}."
            )


def _validate_current_feature_catalog_authority(
    session: Mapping[str, Any],
) -> None:
    """Require the current schema-3 catalog and reject duplicated derived payloads."""

    unknown_fields = sorted(
        str(field)
        for field in session
        if field not in CURRENT_SESSION_TOP_LEVEL_FIELDS
    )
    if unknown_fields:
        raise ValidationError(
            "Session contains unclassified top-level field(s): "
            + ", ".join(unknown_fields)
            + "."
        )

    features = session.get("features")
    if isinstance(features, Mapping):
        duplicated = CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS & set(features)
        if duplicated:
            raise ValidationError(
                f"Session version {CURRENT_SESSION_VERSION} cannot contain derived feature payloads in "
                f"features: {', '.join(sorted(duplicated))}."
            )
    elif "features" in session:
        raise ValidationError("Session features must be an object when present.")

    orthogroup_state = session.get("orthogroupState")
    if isinstance(orthogroup_state, Mapping):
        if "groups" in orthogroup_state:
            raise ValidationError(
                f"Session version {CURRENT_SESSION_VERSION} cannot contain derived orthogroupState.groups."
            )
    elif "orthogroupState" in session:
        raise ValidationError(
            "Session orthogroupState must be an object when present."
        )

    results = session.get("results")
    if "results" not in session or not isinstance(results, list):
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} requires a results array."
        )
    editor_state = session.get("editorState")
    if (
        not isinstance(editor_state, Mapping)
        or "featureCatalog" not in editor_state
    ):
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} "
            + (
                "Results require editorState.featureCatalog."
                if results
                else "requires editorState.featureCatalog."
            )
        )
    catalog = (
        editor_state.get("featureCatalog")
        if isinstance(editor_state, Mapping)
        else None
    )
    if not results:
        if catalog is not None and (
            not isinstance(catalog, Mapping)
            or catalog.get("schema") != CURRENT_FEATURE_CATALOG_SCHEMA
            or catalog.get("items") != []
        ):
            raise ValidationError(
                "An empty Result set requires an empty schema-3 feature catalog."
            )
        return
    if not isinstance(catalog, Mapping):
        raise ValidationError(
            f"Session version {CURRENT_SESSION_VERSION} Results require editorState.featureCatalog."
        )
    items = catalog.get("items")
    if (
        catalog.get("schema") != CURRENT_FEATURE_CATALOG_SCHEMA
        or not isinstance(items, list)
        or len(items) != len(results)
    ):
        raise ValidationError(
            "Session feature catalog must contain one schema-3 item per Result."
        )

    from .web_support.feature_catalog import select_feature_catalog_item

    for result_index, result in enumerate(results):
        if not isinstance(result, Mapping):
            raise ValidationError("Session results must contain objects.")
        raw_result_name = result.get("name")
        result_name = (
            raw_result_name.strip()
            if isinstance(raw_result_name, str)
            else ""
        )
        content = result.get("content")
        if (
            not result_name
            or result_name.lower().endswith(".interactive.svg")
            or not isinstance(content, str)
            or "<svg" not in content
        ):
            raise ValidationError(
                "Each current Session Result must be a named plain SVG."
            )
        if (
            "gbdraw-interactive-feature-metadata" in content
            or "gbdraw-interactive-feature-script" in content
        ):
            raise ValidationError(
                "Current Session Results must contain plain SVG only."
            )
        try:
            catalog_item = select_feature_catalog_item(
                catalog,
                result_index=result_index,
                result_name=result_name,
            )
        except GbdrawError as exc:
            raise ValidationError(str(exc)) from exc
        biological_features = catalog_item.get("biologicalFeatures")
        if not isinstance(biological_features, list) or any(
            not isinstance(feature, Mapping)
            or not str(feature.get("instanceHash") or "").strip()
            for feature in biological_features
        ):
            raise ValidationError(
                f"Session version {CURRENT_SESSION_VERSION} Results require "
                "featureCatalog biologicalFeatures[].instanceHash."
            )


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
    """Validate current cache, manifest, and legacy artifact boundaries."""

    validate_current_web_state_field_names(session.get("config"))
    session_version = session.get("version")
    cache_entries = _artifact_entries(session, "losatCache")
    protein_entries: list[Mapping[str, Any]] = []
    seen_cache_keys: set[str] = set()
    for index, entry in enumerate(cache_entries):
        classification = classify_raw_losat_cache_entry(entry)
        if classification == "protein-legacy":
            raise ValidationError(
                f"Session version {session_version} cannot store legacy protein "
                "entries in losatCache; "
                "use the matching legacyArtifacts candidate envelope."
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
        if not _is_current_derived_cache_entry(entry):
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
    validated_manifest = (
        _validated_protein_identity_manifest(manifest)
        if manifest is not None
        else None
    )
    if manifest is not None and validated_manifest is None:
        raise ValidationError("Invalid proteinIdentityManifest schema-2 artifact.")
    if protein_entries and manifest is None:
        raise ValidationError(
            "Current protein LOSATP cache entries require proteinIdentityManifest."
        )
    if protein_entries:
        assert validated_manifest is not None
        for index, entry in enumerate(protein_entries):
            if not _protein_raw_entry_matches_manifest(entry, validated_manifest):
                raise ValidationError(
                    "Protein LOSATP cache entry does not resolve through the manifest: "
                    f"losatCache.entries[{index}]."
                )
    if derived_entries:
        try:
            validate_current_derived_protein_artifacts(
                derived_entries,
                manifest,
            )
        except ValidationError as exc:
            raise ValidationError(
                "Current derived LOSATP cache contains unresolved protein references."
            ) from exc

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


def migrate_persisted_web_state_field_names(config: object) -> object:
    """Project released Web config field names without mutating persisted data."""

    if not isinstance(config, Mapping):
        return config

    migrated = dict(config)
    adv = config.get("adv")
    if isinstance(adv, Mapping):
        migrated_adv = dict(adv)
        if "depth_tick_interval" in migrated_adv:
            migrated_adv.setdefault(
                "depth_large_tick_interval",
                migrated_adv["depth_tick_interval"],
            )
            migrated_adv.pop("depth_tick_interval")
        depth_tracks = migrated_adv.get("depth_tracks")
        if isinstance(depth_tracks, list):
            migrated_tracks: list[Any] = []
            for track in depth_tracks:
                if not isinstance(track, Mapping) or "tick_interval" not in track:
                    migrated_tracks.append(track)
                    continue
                migrated_track = dict(track)
                migrated_track.setdefault(
                    "large_tick_interval",
                    migrated_track["tick_interval"],
                )
                migrated_track.pop("tick_interval")
                migrated_tracks.append(migrated_track)
            migrated_adv["depth_tracks"] = migrated_tracks
        migrated["adv"] = migrated_adv

    losat = config.get("losat")
    if isinstance(losat, Mapping):
        blastp = losat.get("blastp")
        if isinstance(blastp, Mapping) and "collinearMaxGeneGap" in blastp:
            migrated_blastp = dict(blastp)
            migrated_blastp.setdefault(
                "collinearMaxUnitGap",
                migrated_blastp["collinearMaxGeneGap"],
            )
            migrated_blastp.pop("collinearMaxGeneGap")
            migrated_losat = dict(losat)
            migrated_losat["blastp"] = migrated_blastp
            migrated["losat"] = migrated_losat
    return migrated


def validate_current_web_state_field_names(config: object) -> None:
    """Reject obsolete Web config names at current session write boundaries."""

    if not isinstance(config, Mapping):
        return
    adv = config.get("adv")
    if isinstance(adv, Mapping):
        if "depth_tick_interval" in adv:
            raise ValidationError(
                "Web state field adv.depth_tick_interval is obsolete; "
                "use adv.depth_large_tick_interval."
            )
        depth_tracks = adv.get("depth_tracks")
        if isinstance(depth_tracks, list):
            for index, track in enumerate(depth_tracks):
                if isinstance(track, Mapping) and "tick_interval" in track:
                    raise ValidationError(
                        f"Web state field adv.depth_tracks[{index}].tick_interval "
                        "is obsolete; use large_tick_interval."
                    )
    losat = config.get("losat")
    blastp = losat.get("blastp") if isinstance(losat, Mapping) else None
    if isinstance(blastp, Mapping) and "collinearMaxGeneGap" in blastp:
        raise ValidationError(
            "Web state field losat.blastp.collinearMaxGeneGap is obsolete; "
            "use losat.blastp.collinearMaxUnitGap."
        )


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

    Legacy protein artifacts are kept outside the current cache maps so a
    save-before-generate round trip is lossless.
    """

    source_manifest = (
        protein_identity_manifest
        if protein_identity_manifest is not None
        else session.get("proteinIdentityManifest")
    )
    if source_manifest is None:
        current_manifest = empty_protein_identity_manifest()
    elif _is_valid_protein_identity_manifest(source_manifest):
        current_manifest = _json_clone(source_manifest)
    else:
        raise ValidationError("Cannot write an invalid proteinIdentityManifest.")

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
    source_derived_entries = (
        list(losat_derived_cache_entries)
        if losat_derived_cache_entries is not None
        else list(_artifact_entries(session, "losatDerivedCache"))
    )
    current_derived_entries: list[dict[str, Any]] = []
    imported_derived_evidence: list[dict[str, Any]] = []
    for index, entry in enumerate(source_derived_entries):
        if _is_current_derived_cache_entry(entry):
            current_derived_entries.append(_json_clone(entry))
        elif _is_derived_cache_entry(
            entry, schema=LEGACY_LOSAT_DERIVED_CACHE_SCHEMA
        ):
            imported_derived_evidence.append(_json_clone(entry))
        else:
            raise ValidationError(
                f"Cannot write invalid derived LOSATP cache entry at index {index}."
            )
    session["losatCache"] = {"entries": current_raw_entries}
    session["losatDerivedCache"] = {"entries": current_derived_entries}
    session["proteinIdentityManifest"] = current_manifest

    existing_legacy = session.get("legacyArtifacts")
    normalized_legacy: dict[str, Any] = {}
    existing_candidates = (
        existing_legacy.get("proteinRawCandidates")
        if isinstance(existing_legacy, Mapping)
        else None
    )
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

    existing_evidence = (
        existing_legacy.get("proteinDerivedEvidence")
        if isinstance(existing_legacy, Mapping)
        else None
    )
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


def _is_current_derived_cache_entry(entry: object) -> bool:
    return is_current_derived_protein_artifact(entry)


def _validated_protein_identity_manifest(
    manifest: object,
) -> ProteinIdentityManifest | None:
    if not isinstance(manifest, Mapping):
        return None
    try:
        from .analysis.protein_colinearity import (
            validate_protein_identity_manifest,
        )

        return validate_protein_identity_manifest(manifest)
    except (ImportError, ValidationError, TypeError, ValueError):
        return None


def _is_valid_protein_identity_manifest(manifest: object) -> bool:
    return _validated_protein_identity_manifest(manifest) is not None


def _protein_raw_entry_matches_manifest(
    entry: Mapping[str, Any],
    manifest: ProteinIdentityManifest | Mapping[str, Any],
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


def _legacy_comparison_source(value: Any, fallback: str = "losat") -> str:
    source = str(value or "").strip().lower()
    if source in {"upload", "files", "file"}:
        return "upload"
    if source == "losat":
        return "losat"
    return fallback


def _stable_migrated_comparison_id(
    prefix: str,
    index: int,
    query_uid: str,
    subject_uid: str,
) -> str:
    def safe(value: str) -> str:
        token = _SAFE_FILENAME_RE.sub("-", value).strip("-")
        return token or "record"

    return (
        f"linear-comparison-migrated-{prefix}-{index + 1}-"
        f"{safe(query_uid)}-{safe(subject_uid)}"
    )


def _migrate_legacy_linear_comparison_draft(
    config: Mapping[str, Any],
    files: Mapping[str, Any],
    *,
    force_web_draft: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    migrated_config = _json_clone(dict(config))
    migrated_files = _json_clone(dict(files))
    linear_sequences_value = migrated_files.get("linearSeqs")
    linear_sequences = (
        linear_sequences_value if isinstance(linear_sequences_value, list) else []
    )
    legacy_rows: list[dict[str, Any]] = []
    sanitized_sequences: list[Any] = []
    for sequence in linear_sequences:
        if not isinstance(sequence, Mapping):
            sanitized_sequences.append(sequence)
            legacy_rows.append({"uid": "", "blast": None, "losatFilename": ""})
            continue
        row = dict(sequence)
        legacy_rows.append(
            {
                "uid": str(row.get("uid") or ""),
                "blast": row.pop("blast", None),
                "losatFilename": str(row.pop("losat_filename", "") or ""),
            }
        )
        sanitized_sequences.append(row)
    migrated_files["linearSeqs"] = sanitized_sequences

    comparisons_value = migrated_files.get("linearComparisons")
    file_comparisons = (
        [dict(item) for item in comparisons_value if isinstance(item, Mapping)]
        if isinstance(comparisons_value, list)
        else []
    )
    adv_value = migrated_config.get("adv")
    adv = dict(adv_value) if isinstance(adv_value, Mapping) else {}
    raw_source = migrated_config.get("blastSource", adv.get("blastSource"))
    global_source = _legacy_comparison_source(raw_source)
    legacy_none = str(raw_source or "").strip().lower() == "none"
    migrated_config.pop("blastSource", None)
    adv.pop("blastSource", None)
    migrated_config["adv"] = adv

    layout_value = migrated_config.get("linearRecordLayout")
    layout = dict(layout_value) if isinstance(layout_value, Mapping) else None
    explicit_value = layout.get("comparisons") if layout is not None else None
    explicit = (
        [dict(item) for item in explicit_value if isinstance(item, Mapping)]
        if isinstance(explicit_value, list)
        else None
    )
    if layout is not None:
        layout.pop("comparisons", None)
        migrated_config["linearRecordLayout"] = layout

    existing_plan = migrated_config.get("linearComparisonPlan")
    if isinstance(existing_plan, Mapping):
        plan = _json_clone(dict(existing_plan))
        edges_value = plan.get("edges")
        edges = edges_value if isinstance(edges_value, list) else []
        file_by_id = {
            str(item.get("id") or ""): item.get("file")
            for item in file_comparisons
            if str(item.get("id") or "") and item.get("file")
        }
        bindings = []
        sanitized_edges = []
        for edge in edges:
            if not isinstance(edge, Mapping):
                continue
            metadata = dict(edge)
            file_entry = metadata.pop("file", None) or file_by_id.get(
                str(metadata.get("id") or "")
            )
            sanitized_edges.append(metadata)
            if file_entry:
                bindings.append(
                    {"id": str(metadata.get("id") or ""), "file": file_entry}
                )
        plan["edges"] = sanitized_edges
        migrated_config["linearComparisonPlan"] = plan
        migrated_files["linearComparisons"] = bindings
        return migrated_config, migrated_files

    if not force_web_draft:
        migrated_config.pop("linearRecordLayout", None)
        migrated_config.pop("linearComparisonPlan", None)
        migrated_files["linearComparisons"] = []
        return migrated_config, migrated_files

    uid_index = {
        str(row.get("uid") or ""): index
        for index, row in enumerate(legacy_rows)
        if str(row.get("uid") or "")
    }
    legacy_binding_entries = [
        {
            "index": index,
            "comparison": comparison,
            "id": str(comparison.get("id") or ""),
            "queryUid": str(comparison.get("queryUid") or ""),
            "subjectUid": str(comparison.get("subjectUid") or ""),
            "file": comparison.get("file"),
        }
        for index, comparison in enumerate(file_comparisons)
    ]
    file_by_id: dict[str, Mapping[str, Any]] = {}
    file_by_pair: dict[tuple[str, str], Mapping[str, Any]] = {}
    for entry in legacy_binding_entries:
        file_entry = entry["file"]
        if not file_entry:
            continue
        comparison_id = str(entry["id"])
        query_uid = str(entry["queryUid"])
        subject_uid = str(entry["subjectUid"])
        if comparison_id:
            file_by_id.setdefault(comparison_id, entry)
        if query_uid and subject_uid:
            file_by_pair.setdefault((query_uid, subject_uid), entry)

    consumed_binding_indexes: set[int] = set()

    def consume_binding(entry: Mapping[str, Any] | None) -> Any:
        if not entry or not entry.get("file"):
            return None
        consumed_binding_indexes.add(int(entry["index"]))
        return entry["file"]

    def positional_file(index: int) -> Any:
        if not 0 <= index < len(legacy_rows) - 1:
            return None
        row = legacy_rows[index]
        next_row = legacy_rows[index + 1]
        endpoint_binding = file_by_pair.get(
            (str(row.get("uid") or ""), str(next_row.get("uid") or ""))
        )
        if row.get("blast"):
            # Canonical projection mirrored adjacent uploads into both legacy shapes.
            consume_binding(endpoint_binding)
            return row["blast"]
        return consume_binding(endpoint_binding)

    used_ids: set[str] = set()
    edges_with_files: list[dict[str, Any]] = []

    def add_edge(
        *,
        edge_id: Any,
        query_uid: str,
        subject_uid: str,
        source: str,
        included: bool,
        file_entry: Any,
        file_active: bool,
        losat_filename: str,
        losat_filename_active: bool,
        prefix: str,
        index: int,
    ) -> None:
        if not query_uid or not subject_uid:
            return
        base_id = str(edge_id or "").strip() or _stable_migrated_comparison_id(
            prefix, index, query_uid, subject_uid
        )
        unique_id = base_id
        suffix = 2
        while unique_id in used_ids:
            unique_id = f"{base_id}-{suffix}"
            suffix += 1
        used_ids.add(unique_id)
        edges_with_files.append(
            {
                "id": unique_id,
                "queryUid": query_uid,
                "subjectUid": subject_uid,
                "included": bool(included),
                "fileActive": bool(file_active),
                "losatFilenameActive": bool(losat_filename_active),
                "source": source,
                "losatFilename": str(losat_filename or ""),
                "file": file_entry,
            }
        )

    authoritative_explicit = bool(layout and layout.get("enabled")) and explicit is not None
    mode = "none" if legacy_none else "adjacent"
    used_payload_gaps: set[int] = set()
    if authoritative_explicit:
        mode = "selected" if explicit else "none"
        for index, comparison in enumerate(explicit):
            query_index_value = comparison.get("queryIndex")
            subject_index_value = comparison.get("subjectIndex")
            query_uid = str(comparison.get("queryUid") or "")
            subject_uid = str(comparison.get("subjectUid") or "")
            query_index = uid_index.get(query_uid)
            subject_index = uid_index.get(subject_uid)
            if query_index is None and isinstance(query_index_value, int):
                query_index = query_index_value
                if 0 <= query_index < len(legacy_rows):
                    query_uid = str(legacy_rows[query_index].get("uid") or "")
            if subject_index is None and isinstance(subject_index_value, int):
                subject_index = subject_index_value
                if 0 <= subject_index < len(legacy_rows):
                    subject_uid = str(legacy_rows[subject_index].get("uid") or "")
            adjacent_gap = (
                query_index
                if query_index is not None and subject_index == query_index + 1
                else None
            )
            file_entry = (
                consume_binding(file_by_id.get(str(comparison.get("id") or "")))
                or consume_binding(file_by_pair.get((query_uid, subject_uid)))
                or (positional_file(adjacent_gap) if adjacent_gap is not None else None)
            )
            losat_filename = (
                str(legacy_rows[adjacent_gap].get("losatFilename") or "")
                if adjacent_gap is not None
                else ""
            )
            if adjacent_gap is not None and (file_entry or losat_filename):
                used_payload_gaps.add(adjacent_gap)
            add_edge(
                edge_id=comparison.get("id"),
                query_uid=query_uid,
                subject_uid=subject_uid,
                source=_legacy_comparison_source(
                    comparison.get("source"), global_source
                ),
                included=True,
                file_entry=file_entry,
                file_active=bool(file_entry),
                losat_filename=losat_filename,
                losat_filename_active=bool(losat_filename),
                prefix="selected",
                index=index,
            )

    for index, row in enumerate(legacy_rows[:-1]):
        file_entry = positional_file(index)
        losat_filename = str(row.get("losatFilename") or "")
        if (not file_entry and not losat_filename) or index in used_payload_gaps:
            continue
        payload_active = mode == "adjacent"
        file_active = payload_active and global_source == "upload" and bool(file_entry)
        filename_active = (
            payload_active and global_source == "losat" and bool(losat_filename)
        )
        add_edge(
            edge_id=None,
            query_uid=str(row.get("uid") or ""),
            subject_uid=str(legacy_rows[index + 1].get("uid") or ""),
            source=global_source,
            included=file_active or filename_active,
            file_entry=file_entry,
            file_active=file_active,
            losat_filename=losat_filename,
            losat_filename_active=filename_active,
            prefix="adjacent",
            index=index,
        )

    for entry in legacy_binding_entries:
        entry_index = int(entry["index"])
        file_entry = entry.get("file")
        if not file_entry or entry_index in consumed_binding_indexes:
            continue
        comparison = entry["comparison"]
        query_uid = str(entry["queryUid"])
        subject_uid = str(entry["subjectUid"])
        query_index_value = comparison.get("queryIndex")
        subject_index_value = comparison.get("subjectIndex")
        if not query_uid and isinstance(query_index_value, int):
            if 0 <= query_index_value < len(legacy_rows):
                query_uid = str(legacy_rows[query_index_value].get("uid") or "")
        if not subject_uid and isinstance(subject_index_value, int):
            if 0 <= subject_index_value < len(legacy_rows):
                subject_uid = str(legacy_rows[subject_index_value].get("uid") or "")
        add_edge(
            edge_id=entry["id"],
            query_uid=query_uid,
            subject_uid=subject_uid,
            source=_legacy_comparison_source(
                comparison.get("source"), global_source
            ),
            included=False,
            file_entry=file_entry,
            file_active=False,
            losat_filename="",
            losat_filename_active=False,
            prefix="retained",
            index=entry_index,
        )

    migrated_config["linearComparisonPlan"] = {
        "mode": mode,
        "defaultSource": global_source,
        "edges": [
            {key: value for key, value in edge.items() if key != "file"}
            for edge in edges_with_files
        ],
    }
    migrated_files["linearComparisons"] = [
        {"id": edge["id"], "file": edge["file"]}
        for edge in edges_with_files
        if edge.get("file")
    ]
    return migrated_config, migrated_files


def migrate_legacy_linear_comparison_draft_for_current_writer(
    config: Mapping[str, Any],
    files: Mapping[str, Any],
    *,
    force_web_draft: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Project a supported pre-v40 comparison draft into final v40 ownership."""

    return _migrate_legacy_linear_comparison_draft(
        config,
        files,
        force_web_draft=force_web_draft,
    )


def _embedded_entry_bytes(entry: Mapping[str, Any]) -> bytes | None:
    encoding = entry.get("encoding")
    data = entry.get("data")
    if encoding != DEPTH_FILE_ENCODING and isinstance(data, str):
        try:
            return base64.b64decode(data, validate=True)
        except (binascii.Error, ValueError):
            return None
    if encoding == DEPTH_FILE_ENCODING and isinstance(data, Mapping):
        try:
            return decode_depth_payload(data).encode("utf-8")
        except ValidationError:
            return None
    return None


def _attach_current_web_file_bindings(
    payload: dict[str, Any],
    files: Mapping[str, Any],
) -> None:
    resources_value = payload.get("resources")
    if not isinstance(resources_value, dict):
        raise ValidationError("Current session resources must be an object.")
    resources = resources_value
    identities: dict[tuple[int, str], str] = {}
    for resource_id, resource in resources.items():
        if not isinstance(resource, Mapping):
            continue
        resource_bytes = _embedded_entry_bytes(resource)
        if resource_bytes is None:
            continue
        identity = (len(resource_bytes), hashlib.sha256(resource_bytes).hexdigest())
        identities.setdefault(identity, str(resource_id))

    next_number = 1

    def add_file(entry: Any) -> dict[str, Any] | None:
        nonlocal next_number
        if not isinstance(entry, Mapping):
            return None
        file_bytes = _embedded_entry_bytes(entry)
        if file_bytes is None:
            return None
        identity = (len(file_bytes), hashlib.sha256(file_bytes).hexdigest())
        resource_id = identities.get(identity)
        if resource_id is None:
            while True:
                candidate = f"resource-{next_number:04d}"
                next_number += 1
                if candidate not in resources:
                    resource_id = candidate
                    break
            name = safe_embedded_filename(
                entry.get("name"), fallback="resource.dat"
            )
            resources[resource_id] = {
                "kind": "web-file",
                "name": f"{resource_id}-{name}",
                "type": str(entry.get("type") or "application/octet-stream"),
                "size": len(file_bytes),
                "lastModified": int(entry.get("lastModified") or 0),
                "encoding": "base64",
                "data": base64.b64encode(file_bytes).decode("ascii"),
            }
            identities[identity] = resource_id
        return {
            "resourceId": resource_id,
            "name": str(entry.get("name") or "file"),
            "type": str(entry.get("type") or ""),
            "lastModified": int(entry.get("lastModified") or 0),
        }

    def add_value(value: Any) -> Any:
        if isinstance(value, list):
            return [add_value(item) for item in value]
        return add_file(value)

    linear_sequences_value = files.get("linearSeqs")
    linear_sequences = (
        linear_sequences_value if isinstance(linear_sequences_value, list) else []
    )
    sequence_bindings = []
    for sequence in linear_sequences:
        if not isinstance(sequence, Mapping):
            continue
        sequence_bindings.append(
            {
                "uid": str(sequence.get("uid") or ""),
                "gb": add_file(sequence.get("gb")),
                "gff": add_file(sequence.get("gff")),
                "fasta": add_file(sequence.get("fasta")),
                "depth": add_value(sequence.get("depth")),
                "losat_gencode": sequence.get("losat_gencode", 1),
                "definition": str(sequence.get("definition") or ""),
                "record_subtitle": str(sequence.get("record_subtitle") or ""),
                "region_record_id": str(sequence.get("region_record_id") or ""),
                "region_start": sequence.get("region_start"),
                "region_end": sequence.get("region_end"),
                "region_reverse": bool(sequence.get("region_reverse")),
            }
        )
    comparison_bindings = []
    comparisons_value = files.get("linearComparisons")
    comparisons = comparisons_value if isinstance(comparisons_value, list) else []
    for comparison in comparisons:
        if not isinstance(comparison, Mapping):
            continue
        binding = add_file(comparison.get("file"))
        comparison_id = str(comparison.get("id") or "")
        if comparison_id and binding is not None:
            comparison_bindings.append({"id": comparison_id, "file": binding})

    bindings = {
        "schema": 1,
        "c_gb": add_file(files.get("c_gb")),
        "c_gff": add_file(files.get("c_gff")),
        "c_fasta": add_file(files.get("c_fasta")),
        "c_depth": add_value(files.get("c_depth")),
        "c_conservation_blasts": add_value(
            files.get("c_conservation_blasts")
        ),
        "c_conservation_blasts_source": (
            "losat-cache"
            if files.get("c_conservation_blasts_source") == "losat-cache"
            else None
        ),
        "c_conservation_fastas": add_value(files.get("c_conservation_fastas")),
        "c_conservation_sequence_sources": add_value(
            files.get("c_conservation_sequence_sources")
        ),
        "d_color": add_file(files.get("d_color")),
        "t_color": add_file(files.get("t_color")),
        "blacklist": add_file(files.get("blacklist")),
        "whitelist": add_file(files.get("whitelist")),
        "qualifier_priority": add_file(files.get("qualifier_priority")),
        "linearSeqs": sequence_bindings,
        "linearComparisons": comparison_bindings,
    }
    web_files_value = payload.get("webFiles")
    web_files = (
        _json_clone(dict(web_files_value))
        if isinstance(web_files_value, Mapping)
        else {}
    )
    web_files.pop("bindings", None)
    metadata_value = web_files.get("linearRecordMetadata")
    if isinstance(metadata_value, list):
        web_files["linearRecordMetadata"] = [
            {
                key: value
                for key, value in metadata.items()
                if key not in {"losatFilename", "losat_filename"}
            }
            if isinstance(metadata, Mapping)
            else metadata
            for metadata in metadata_value
        ]
    resource_aliases: dict[str, str] = {}
    for source_field, binding_field in (
        ("conservationLosatFastaSources", "c_conservation_fastas"),
        ("conservationSequenceSources", "c_conservation_sequence_sources"),
    ):
        source_ids = web_files.get(source_field)
        rebound_files = bindings[binding_field]
        if not isinstance(source_ids, list) or not isinstance(rebound_files, list):
            continue
        rewritten_ids: list[str | None] = []
        for index, source_id_value in enumerate(source_ids):
            rebound = rebound_files[index] if index < len(rebound_files) else None
            rebound_id = (
                str(rebound.get("resourceId") or "").strip()
                if isinstance(rebound, Mapping)
                else ""
            )
            rewritten_ids.append(rebound_id or None)
            source_id = str(source_id_value or "").strip()
            if source_id and rebound_id:
                resource_aliases[source_id] = rebound_id
        web_files[source_field] = rewritten_ids
    original_names = web_files.get("resourceOriginalNames")
    if isinstance(original_names, Mapping):
        web_files["resourceOriginalNames"] = {
            resource_aliases.get(str(resource_id), str(resource_id)): name
            for resource_id, name in original_names.items()
            if resource_aliases.get(str(resource_id), str(resource_id)) in resources
        }
    web_files["bindings"] = bindings
    payload["webFiles"] = web_files


def build_session_json(
    context: SessionBuildContext,
    *,
    svg_results: Sequence[tuple[str, str]],
    embedded_files: Mapping[str, Any],
    generated_at: datetime,
    feature_catalog: Mapping[str, Any] | None = None,
    losat_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    losat_derived_cache_entries: Sequence[Mapping[str, Any]] | None = None,
    protein_identity_manifest: Mapping[str, Any] | None = None,
    legacy_protein_raw_candidates: Sequence[Mapping[str, Any]] | None = None,
    legacy_protein_derived_evidence: Sequence[Mapping[str, Any]] | None = None,
    canonical_request: DiagramRequest | None = None,
    _canonical_request_is_resolved: bool = False,
) -> dict[str, Any]:
    """Build a GUI-loadable session JSON payload from a CLI run."""

    source_version: int | None = None
    if context.source_session is not None:
        validate_session(context.source_session)
        source_version = int(context.source_session["version"])
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
        config = {"adv": {}}
        payload["config"] = config
    elif source_version is not None and source_version < CURRENT_SESSION_VERSION:
        migrated_config = migrate_persisted_web_state_field_names(config)
        assert isinstance(migrated_config, dict)
        config = migrated_config
        payload["config"] = config

    ui = payload.get("ui")
    if not isinstance(ui, dict):
        ui = {}
        payload["ui"] = ui
    ui["mode"] = context.mode
    ui.pop("blastSource", None)
    ui.setdefault("zoom", 1)
    ui.setdefault("selectedResultIndex", 0)
    ui.setdefault("canvasPan", {"x": 0, "y": 0})
    ui.setdefault("canvasPadding", {"top": 0, "right": 0, "bottom": 0, "left": 0})

    payload["files"] = _json_clone(embedded_files)
    payload["results"] = [
        {"name": name or f"Result {index + 1}", "content": content}
        for index, (name, content) in enumerate(svg_results)
    ]
    editor_state = payload.get("editorState")
    editor_state = (
        dict(editor_state) if isinstance(editor_state, Mapping) else {}
    )
    editor_state["featureCatalog"] = _json_clone(
        feature_catalog
        if feature_catalog is not None
        else {
            "schema": CURRENT_FEATURE_CATALOG_SCHEMA,
            "items": [
                {
                    "resultIndex": index,
                    "resultName": result["name"],
                    "recordKeys": [],
                    "features": [],
                    "biologicalFeatures": [],
                    "orthogroups": [],
                    "annotations": [],
                    "comparisonMatches": [],
                }
                for index, result in enumerate(payload["results"])
            ],
        }
    )
    payload["editorState"] = editor_state

    features = payload.get("features")
    features = dict(features) if isinstance(features, Mapping) else {}
    for key in CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS:
        features.pop(key, None)
    payload["features"] = features

    orthogroup_state = payload.get("orthogroupState")
    orthogroup_state = (
        dict(orthogroup_state)
        if isinstance(orthogroup_state, Mapping)
        else {}
    )
    orthogroup_state.pop("groups", None)
    payload["orthogroupState"] = orthogroup_state
    payload["cliInvocation"] = {
        "schema": 1,
        "mode": context.mode,
        "args": [str(arg) for arg in context.cli_invocation_args],
        "renderFormats": [normalize_format_token(fmt) for fmt in context.render_formats],
        "fileBindings": [_binding_to_json(binding) for binding in context.file_bindings],
        "generatedBy": "gbdraw",
    }
    if canonical_request is not None:
        from .session import (
            _build_session_document_from_resolved_request,
            build_session_document,
        )

        build_document = (
            _build_session_document_from_resolved_request
            if _canonical_request_is_resolved
            else build_session_document
        )
        canonical_document = build_document(
            canonical_request,
            created_at=generated_at,
        )
        canonical = canonical_document.to_dict()
        payload["renderRequest"] = canonical["renderRequest"]
        payload["resources"] = canonical["resources"]
    else:
        raise ValidationError(
            f"A canonical typed request is required to write a version {CURRENT_SESSION_VERSION} session."
        )
    files_value = payload.get("files")
    files_for_web = files_value if isinstance(files_value, Mapping) else {}
    if source_version is not None and source_version < CURRENT_SESSION_VERSION:
        force_web_comparison_draft = (
            isinstance(config.get("linearRecordLayout"), Mapping)
            or isinstance(config.get("linearComparisonPlan"), Mapping)
            or not isinstance(config.get("cliOptions"), Mapping)
        )
        config, files_for_web = _migrate_legacy_linear_comparison_draft(
            config,
            files_for_web,
            force_web_draft=force_web_comparison_draft,
        )
        payload["config"] = config
    _attach_current_web_file_bindings(payload, files_for_web)
    payload.pop("files", None)
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


def write_session_json(
    path: str | Path,
    payload: Mapping[str, Any],
    *,
    overwrite: bool = True,
) -> None:
    """Write plain or ``.gz`` session JSON through a same-directory stage."""

    expanded_payload = expand_session_feature_catalog(payload)
    validate_session(expanded_payload)
    serialized_payload = compact_session_feature_catalog(expanded_payload)

    output_path = Path(path)
    temp_path: Path | None = None
    temp_fd: int | None = None
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temp_fd, temp_name = tempfile.mkstemp(
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
        )
        temp_path = Path(temp_name)
        if output_path.suffix.lower() == ".gz":
            raw_file = os.fdopen(temp_fd, "wb")
            temp_fd = None
            with raw_file:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw_file,
                    compresslevel=9,
                    mtime=0,
                ) as compressed_file:
                    with io.TextIOWrapper(compressed_file, encoding="utf-8") as text_file:
                        json.dump(
                            serialized_payload,
                            text_file,
                            ensure_ascii=False,
                            separators=(",", ":"),
                        )
                raw_file.flush()
                os.fsync(raw_file.fileno())
        else:
            text_file = os.fdopen(temp_fd, "w", encoding="utf-8")
            temp_fd = None
            with text_file:
                json.dump(
                    serialized_payload,
                    text_file,
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
                text_file.flush()
                os.fsync(text_file.fileno())
        commit_staged_output_file(
            temp_path,
            output_path,
            overwrite=overwrite,
        )
    except FileExistsError as exc:
        raise ValidationError(
            f"Session output already exists: {output_path}. "
            "Pass overwrite=True to replace it."
        ) from exc
    except OSError as exc:
        raise ValidationError(f"Could not write session sidecar: {output_path}") from exc
    finally:
        if temp_fd is not None:
            try:
                os.close(temp_fd)
            except OSError:
                pass
        try:
            if temp_path is not None:
                temp_path.unlink(missing_ok=True)
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


def _canonicalize_legacy_session_cli_args(
    args: Sequence[str],
    *,
    mode: Literal["circular", "linear"],
) -> tuple[list[str], dict[int, int]]:
    replacements = {
        "--depth": "--depth_track",
        "--depth_tick_interval": "--depth_large_tick_interval",
        "--gc_content_tick_interval": "--gc_content_large_tick_interval",
        "--feature_table": "--feature_visibility_table",
        "--annotation-table": "--annotation_table",
        "--losatp-bin": "--losatp_bin",
        "--ncbi-blastp-bin": "--ncbi_blastp_bin",
        "--losatp-threads": "--losatp_threads",
        "--protein-blastp-mode": "--protein_blastp_mode",
        "--protein-blastp-max-hits": "--protein_blastp_max_hits",
        "--protein-blastp-candidate-limit": "--protein_blastp_candidate_limit",
        "--align-orthogroup-feature": "--align_orthogroup_feature",
        "--collinear-unit-mode": "--collinear_unit_mode",
        "--collinear-search-scope": "--collinear_search_scope",
        "--collinear-min-anchors": "--collinear_min_anchors",
        "--collinear_max_gene_gap": "--collinear_max_unit_gap",
        "--collinear-max-unit-gap": "--collinear_max_unit_gap",
        "--collinear-max-gene-gap": "--collinear_max_unit_gap",
        "--collinear-max-diagonal-drift": "--collinear_max_diagonal_drift",
        "--collinear-max-conflicts-in-merge-gap": "--collinear_max_conflicts_in_merge_gap",
        "--collinear-max-paralog-links-per-orthogroup": (
            "--collinear_max_paralog_links_per_orthogroup"
        ),
        "--collinear-color-mode": "--collinear_color_mode",
        "--keep-definition-left-aligned": "--keep_definition_left_aligned",
        "--pairwise-match-style": "--pairwise_match_style",
        "--definition-line-style": "--definition_line_style",
        "--record-subtitle": "--record_subtitle",
        "--circular-track-slot": "--circular_track_slot",
    }
    if mode == "circular":
        replacements.update(
            {
                "--suppress_gc": "--no-gc",
                "--suppress_skew": "--no-skew",
            }
        )
    else:
        replacements.update(
            {
                "--show_gc": "--gc",
                "--show_skew": "--skew",
            }
        )

    canonical_args: list[str] = []
    source_to_canonical_index: dict[int, int] = {}
    for source_index, raw_token in enumerate(args):
        token = str(raw_token)
        if token == "--show_depth":
            continue
        option, separator, inline_value = (
            token.partition("=")
            if token.startswith("--")
            else (token, "", "")
        )
        option = replacements.get(option, option)
        if separator:
            if mode == "circular" and option == "--circular_track_slot":
                inline_value = _migrate_legacy_circular_slot_cli_value(inline_value)
            if option == "--multi_record_size_mode" and inline_value == "sqrt":
                inline_value = "auto"
            elif mode == "linear" and option == "--label_placement":
                inline_value = {
                    "on_feature": "above_feature",
                }.get(inline_value, inline_value)
            elif mode == "linear" and option == "--track_layout":
                inline_value = {
                    "spreadout": "above",
                    "tuckin": "below",
                }.get(inline_value, inline_value)
            token = f"{option}={inline_value}"
        else:
            token = option
            if canonical_args:
                previous_option = canonical_args[-1].partition("=")[0]
                if previous_option == "--multi_record_size_mode" and token == "sqrt":
                    token = "auto"
                elif (
                    mode == "linear"
                    and previous_option == "--label_placement"
                    and token == "on_feature"
                ):
                    token = "above_feature"
                elif mode == "linear" and previous_option == "--track_layout":
                    token = {"spreadout": "above", "tuckin": "below"}.get(token, token)
                elif mode == "circular" and previous_option == "--circular_track_slot":
                    token = _migrate_legacy_circular_slot_cli_value(token)
        source_to_canonical_index[source_index] = len(canonical_args)
        canonical_args.append(token)
    return canonical_args, source_to_canonical_index


def _migrate_legacy_circular_slot_cli_value(value: str) -> str:
    """Move retired slot fields into the private persisted-data transport."""

    head, separator, raw_options = str(value).partition("@")
    if not separator:
        return str(value)
    migrated: list[str] = []
    for raw_part in raw_options.split(","):
        part = raw_part.strip()
        if not part or "=" not in part:
            migrated.append(part)
            continue
        raw_key, raw_value = part.split("=", 1)
        key = raw_key.strip().lower()
        if key in {"strict", "compress", "reserve"}:
            continue
        if key == "spacing":
            migrated.append(f"__gbdraw_legacy_spacing={raw_value.strip()}")
        else:
            migrated.append(part)
    return head if not migrated else f"{head}@{','.join(migrated)}"


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

    session_version = int(session.get("version", 0))
    migrate_legacy_cli = session_version < CURRENT_SESSION_VERSION
    _restore_cli_table_paths(
        session,
        run_args,
        temp_dir=temp_dir,
        migrate_legacy_cli=migrate_legacy_cli,
    )

    if migrate_legacy_cli:
        run_args, _ = _canonicalize_legacy_session_cli_args(run_args, mode=mode)
        invocation_args, invocation_index_map = _canonicalize_legacy_session_cli_args(
            invocation_args,
            mode=mode,
        )
        remapped_bindings: list[SessionFileBinding] = []
        for binding in file_bindings:
            if binding.argIndex not in invocation_index_map:
                raise ValidationError(
                    "cliInvocation.fileBindings cannot reference a removed legacy CLI flag."
                )
            remapped_bindings.append(
                SessionFileBinding(
                    argIndex=invocation_index_map[binding.argIndex],
                    slot=binding.slot,
                    name=binding.name,
                )
            )
        file_bindings = remapped_bindings

    run_args = _apply_option_override(run_args, "-o", "--output", output_override)
    run_args = _apply_option_override(run_args, "-f", "--format", format_override)
    invocation_args = _apply_option_override(invocation_args, "-o", "--output", output_override)
    invocation_args = _apply_option_override(invocation_args, "-f", "--format", format_override)
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
    _append_common_gui_args(run_args, invocation_args, form=form, adv=adv)

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
    migrate_legacy_cli: bool = False,
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

        table_kind = str(table_entry.get("kind") or "").strip()
        preceding_option = (
            str(run_args[arg_index - 1]) if arg_index > 0 else ""
        )
        if (
            migrate_legacy_cli
            and (
                table_kind == "circular_track"
                or preceding_option == "--circular_track_table"
            )
        ):
            _migrate_legacy_circular_track_table(table_path)

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


def _migrate_legacy_circular_track_table(table_path: Path) -> None:
    """Move a persisted Circular table's spacing column to reader-only params."""

    try:
        with table_path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
    except OSError as exc:
        raise ValidationError(
            f"Could not read restored Circular track table: {table_path}"
        ) from exc
    if not rows:
        return

    header_index = next(
        (
            index
            for index, cells in enumerate(rows)
            if cells and any(str(cell).strip() for cell in cells)
        ),
        None,
    )
    if header_index is None:
        return
    header = [str(cell).strip() for cell in rows[header_index]]
    if "spacing" not in header:
        return

    spacing_index = header.index("spacing")
    if "params" not in header:
        header.append("params")
    params_index = header.index("params")
    migrated_rows: list[list[str]] = []
    for index, cells in enumerate(rows):
        if index < header_index:
            migrated_rows.append([str(cell) for cell in cells])
            continue
        if index == header_index:
            migrated_rows.append(
                [column for column in header if column != "spacing"]
            )
            continue
        values = [str(cell) for cell in cells]
        while len(values) < len(header):
            values.append("")
        values = values[: len(header)]
        spacing = values[spacing_index].strip()
        if spacing:
            legacy_param = f"__gbdraw_legacy_spacing={spacing}"
            existing_params = values[params_index].strip()
            values[params_index] = (
                f"{existing_params},{legacy_param}"
                if existing_params
                else legacy_param
            )
        migrated_rows.append(
            [
                value
                for column, value in zip(header, values, strict=True)
                if column != "spacing"
            ]
        )
    try:
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(migrated_rows)
    except OSError as exc:
        raise ValidationError(
            f"Could not migrate restored Circular track table: {table_path}"
        ) from exc


def _append_common_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    if _string_or_none(form.get("species")):
        _append_pair(run_args, invocation_args, "--species", str(form.get("species")))
    if _string_or_none(form.get("strain")):
        _append_pair(run_args, invocation_args, "--strain", str(form.get("strain")))
    if form.get("separate_strands") is True:
        _append_flag(run_args, invocation_args, "--separate_strands")
    if form.get("show_scale") is False:
        _append_flag(run_args, invocation_args, "--hide_scale")
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
        ("arrow_head_length_ratio", "--arrow_head_length_ratio"),
        ("arrow_shaft_width_ratio", "--arrow_shaft_width_ratio"),
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
        _append_flag(run_args, invocation_args, "--no-gc")
    if form.get("suppress_skew") is True:
        _append_flag(run_args, invocation_args, "--no-skew")
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
            if (
                key == "multi_record_size_mode"
                and str(value).strip().lower() == "sqrt"
            ):
                value = "auto"
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
        _append_flag(run_args, invocation_args, "--gc")
    if form.get("show_skew") is True:
        _append_flag(run_args, invocation_args, "--skew")
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
            if key == "label_placement" and str(value).strip().lower() == "on_feature":
                value = "above_feature"
            _append_pair(run_args, invocation_args, option, str(value))
    if _string_or_none(form.get("linear_track_layout")):
        layout = str(form.get("linear_track_layout")).strip().lower()
        layout = {"spreadout": "above", "tuckin": "below"}.get(layout, layout)
        _append_pair(run_args, invocation_args, "--track_layout", layout)
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
        ("collinearMaxUnitGap", "--collinear_max_unit_gap"),
        ("collinearMaxDiagonalDrift", "--collinear_max_diagonal_drift"),
        ("collinearMaxConflictsInMergeGap", "--collinear_max_conflicts_in_merge_gap"),
        ("collinearUnitMode", "--collinear_unit_mode"),
        ("collinearSearchScope", "--collinear_search_scope"),
        ("collinearColorMode", "--collinear_color_mode"),
        ("collinearMaxParalogLinksPerOrthogroup", "--collinear_max_paralog_links_per_orthogroup"),
    ):
        value = blastp_cfg.get(key)
        if (
            key == "collinearMaxUnitGap"
            and value in (None, "")
            and int(session.get("version", 0)) < CURRENT_SESSION_VERSION
        ):
            value = blastp_cfg.get("collinearMaxGeneGap")
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
        if separator and feature_type and shape in {
            "arrow",
            "rectangle",
            "underlay",
        }:
            shapes[feature_type] = shape
    return shapes


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
    "CANONICAL_REQUEST_SCHEMAS_BY_SESSION_VERSION",
    "CURRENT_SESSION_VERSION",
    "CANONICAL_SESSION_MIN_VERSION",
    "DEPTH_FILE_ENCODING",
    "DEPTH_FILE_SCHEMA",
    "FEATURE_CATALOG_ENCODING",
    "FEATURE_CATALOG_SCHEMA",
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
    "compact_session_feature_catalog",
    "decode_depth_payload",
    "encode_depth_text",
    "empty_protein_identity_manifest",
    "expand_session_feature_catalog",
    "get_session_slot",
    "load_session",
    "materialize_embedded_file",
    "migrate_legacy_linear_comparison_draft_for_current_writer",
    "migrate_persisted_web_state_field_names",
    "migrate_legacy_repeat_feature_shape_args",
    "normalize_current_session_artifacts",
    "safe_embedded_filename",
    "serialize_file_entry",
    "session_mode",
    "session_to_cli_args",
    "validate_session",
    "validate_current_session_artifacts",
    "validate_current_web_state_field_names",
    "write_session_json",
]
