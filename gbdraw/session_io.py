#!/usr/bin/env python
# coding: utf-8

"""GUI session JSON loading, validation, materialization, and sidecar building."""

from __future__ import annotations

import base64
import binascii
import copy
import csv
import functools
import gzip
import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, NoReturn, Sequence, TextIO, cast

from .analysis.protein_artifacts import (
    CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
    is_current_derived_protein_artifact,
    validate_current_derived_protein_artifacts,
)
from .definition_line_styles import DEFINITION_LINE_KINDS
from .exceptions import GbdrawError, ValidationError
from .io.colors import named_color_hex
from .io.filenames import _SAFE_FILENAME_RE, safe_embedded_filename
from .render.formats import normalize_format_token
from .render.output_paths import commit_staged_output_file
from .web_support.mode_scoped_settings import (
    DIAGRAM_MODES,
    LOSAT_EXECUTION_FIELDS,
    MODE_SCOPED_SETTINGS,
    SLICE_CONTAINERS,
    DiagramMode,
    ModeScopedSetting,
    unmanaged_config_override_modes,
)

if TYPE_CHECKING:
    from .analysis.protein_colinearity import ProteinIdentityManifest
    from .api.requests import DiagramRequest

SESSION_FORMAT = "gbdraw-session"
CURRENT_SESSION_VERSION = 46
# Version 44 is the first with the current active-config and record-display
# draft shapes; version 46 keeps the Web draft of each diagram mode in
# ``modes`` (PD-OI-086) and keys per-feature edits by source identity.
TYPED_DRAFT_SESSION_MIN_VERSION = 44
MODE_SCOPED_SESSION_MIN_VERSION = 46
CURRENT_AUTHORITY_SESSION_MIN_VERSION = 40
CANONICAL_SESSION_MIN_VERSION = 31
SUPPORTED_SESSION_VERSIONS = frozenset(
    {27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 44, CURRENT_SESSION_VERSION}
)
CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39
PROTEIN_LOSAT_CACHE_SCHEMA = 4
NUCLEOTIDE_LOSAT_CACHE_SCHEMA = 2
LOSAT_DERIVED_CACHE_SCHEMA = CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA
LEGACY_LOSAT_DERIVED_CACHE_SCHEMA = 1
PROTEIN_IDENTITY_MANIFEST_SCHEMA = 2
LEGACY_PROTEIN_CANDIDATE_SCHEMA = 1
FEATURE_CATALOG_SCHEMA = 1
FEATURE_CATALOG_ENCODING = "biological-authority-v1"
CURRENT_FEATURE_CATALOG_SCHEMA = 5
# The catalog schema each current-authority Session version writes.
FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION = {44: 4, CURRENT_SESSION_VERSION: CURRENT_FEATURE_CATALOG_SCHEMA}
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
        "modes",
        # Session 46: the CLI options a Session 27-44 draft kept in config.
        "cliOptions",
        # E1 (Q0): the other diagram mode's committed Result set.
        "otherModeResult",
    }
)
# The flat Web draft of Sessions 44 and older; Session 46 keeps it per mode in
# ``modes`` (MODE_SCOPED_SETTINGS).
FLAT_DRAFT_TOP_LEVEL_FIELDS = frozenset({"config", "features"})
# Session 46 keeps the other diagram mode's Result set in ``otherModeResult``:
# the fields of one committed set, named as at the top level.
OTHER_MODE_RESULT_FIELDS = frozenset(
    {"renderRequest", "results", "editorState", "ui", "runMetadata", "cliInvocation"}
)
# The per-set part of the shared ``ui`` and ``editorState`` objects.
OTHER_MODE_RESULT_UI_FIELDS = frozenset(
    {
        "selectedResultIndex",
        "generatedLegendPosition",
        "generatedMultiRecordCanvas",
        "generatedCircularPlotTitlePosition",
        "appliedPaletteName",
        "appliedPaletteColors",
    }
)
OTHER_MODE_RESULT_EDITOR_FIELDS = frozenset(
    {"featureCatalog", "alignmentResetReceipt", "legend", "originalSvgStroke"}
)
OTHER_MODE_RESULT_LEGEND_FIELDS = frozenset({"originalOrder", "originalColors"})
# The Web reader names an unusable ``otherModeResult`` as a Session field error.
_OTHER_MODE_RESULT_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}
CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS = frozenset(
    {
        "extractedFeatures",
        "biologicalFeatures",
        "featureSelectorSafetyScope",
        "featureRecordIds",
        "featureCatalog",
    }
)
# Version 46 replaces the rendered-ID-keyed per-feature edit maps with
# features.featureOverrides (design Q4); its writer never stores them.
RETIRED_RENDERED_ID_FEATURE_FIELDS = frozenset(
    {
        "featureVisibilityOverrides",
        "labelVisibilityOverrides",
        "labelTextFeatureOverrides",
        "labelTextFeatureOverrideSources",
    }
)
# While a Session 41-44 draft is migrated, a draft row names the mode of its
# record key (``scope``), since both modes can use the same record key for the
# same feature; the split then moves each row into that mode's slice, where the
# row has no ``scope`` (Session 46).
DRAFT_SCOPES = ("circular", "linear")
# The fields of a per-feature edit draft row in a Session 46 mode slice.
FEATURE_OVERRIDE_DRAFT_FIELDS = frozenset(
    {
        "recordKey",
        "biologicalFeatureId",
        "featureVisibility",
        "labelVisibility",
        "labelText",
        "labelSourceText",
    }
)
# A Session 46 field that names an invalid shape is a Session field error.
_SESSION_FIELDS_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}
DEPTH_FILE_ENCODING = "gbdraw-depth-table-v1"
DEPTH_FILE_SCHEMA = 1
JS_MAX_SAFE_INTEGER = 9_007_199_254_740_991

_DEPTH_COLUMNS = ("reference_name", "position", "depth")
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
    version = expanded_session.get("version")
    # validate_session reports a missing or non-integer version.
    if not isinstance(version, int) or version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION:
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
    """Load and validate a plain or gzip-compressed gbdraw GUI session JSON file.

    It validates the document as ``gbdraw.session.load_session_document`` does,
    canonical resource descriptors included.
    """

    from .session import _validate_document

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
    _validate_document(payload)
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
    if "otherModeResult" in session and version < CURRENT_SESSION_VERSION:
        raise ValidationError(
            f"Session version {version} cannot contain otherModeResult.", diagnostic=_OTHER_MODE_RESULT_INVALID
        )
    if version >= CANONICAL_SESSION_MIN_VERSION:
        render_request = session.get("renderRequest")
        resources = session.get("resources")
        settings_only = is_settings_only_session(session)
        if not settings_only:
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
        if not isinstance(resources, Mapping):
            raise ValidationError(
                f"Session version {version} requires a canonical resources object."
            )
        files = session.get("files")
        if version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION and "files" in session:
            raise ValidationError(
                f"Session version {version} cannot contain legacy files; "
                "use resources and webFiles."
            )
        if files is not None and not isinstance(files, Mapping):
            raise ValidationError("Session files must be an object when present.")
    else:
        files = session.get("files")
        if files is None or not isinstance(files, Mapping):
            raise ValidationError("Session files are required for CLI regeneration.")
    _validate_web_file_bindings(session)
    if version >= CURRENT_ARTIFACT_SESSION_MIN_VERSION:
        validate_current_session_artifacts(session)
    if version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION:
        _validate_current_top_level_fields(session)
        _validate_mode_scoped_fields(session, version)
        _validate_current_retired_active_config_paths(session)
        _validate_current_comparison_authority(session)
        _validate_current_feature_catalog_authority(session, version)
        _validate_alignment_reset_receipt(session)
        _validate_other_mode_result(session, version)
    if version >= MODE_SCOPED_SESSION_MIN_VERSION:
        _validate_mode_slices(session)
    elif version >= 41:
        _validate_display_placement_drafts(session)
    if is_settings_only_session(session):
        _validate_settings_only_session(session)


def is_settings_only_session(session: Mapping[str, Any]) -> bool:
    """Recognize the explicit document variant, never a missing-resource error."""
    return session.get("version") in (42, 44, CURRENT_SESSION_VERSION) and "renderRequest" in session and session["renderRequest"] is None


def _validate_settings_only_session(session: Mapping[str, Any]) -> None:
    web_files = session.get("webFiles", {})
    bindings = web_files.get("bindings") if isinstance(web_files, Mapping) else None
    if (not isinstance(bindings, Mapping) or not isinstance(bindings.get("linearSeqs"), list)
            or not {"c_gb", "c_gff", "c_fasta"} <= bindings.keys()):
        raise ValidationError("Settings-only Session requires an explicit Web input inventory.")

    def has_input(value: Any) -> bool:
        return any(map(has_input, value)) if isinstance(value, list) else value is not None

    if (set(web_files) != {"bindings"}
            or any(has_input(bindings.get(key)) for key in (
                "c_gb", "c_gff", "c_fasta", "c_conservation_fastas", "c_conservation_sequence_sources"))
            or any(not isinstance(row, Mapping) or any(has_input(row.get(key))
                for key in ("gb", "gff", "fasta")) for row in bindings["linearSeqs"])):
        raise ValidationError("Settings-only Session cannot contain biological sources.")
    manifest = session.get("proteinIdentityManifest") or {}
    if (session.get("results") != [] or session.get("editorState", {}).get("featureCatalog") is not None
            or session.get("cliInvocation") is not None or session.get("runMetadata")
            or session.get("legacyArtifacts")
            or any((session.get(key) or {}).get("entries") for key in ("losatCache", "losatDerivedCache"))
            or any(manifest.get(key) for key in ("proteinSets", "recordAnalyses", "recordInstances"))
            or bindings.get("c_conservation_blasts_source") == "losat-cache"):
        raise ValidationError("Settings-only Session cannot contain committed render artifacts.")
    ui = session.get("ui")
    mode = ui.get("mode") if isinstance(ui, Mapping) else None
    # Session 46 keeps the shown mode's configuration in its mode slice.
    config = (
        _mode_slice_config(session, mode)
        if session.get("version", 0) >= MODE_SCOPED_SESSION_MIN_VERSION
        else session.get("config")
    )
    if (not isinstance(config, Mapping) or not isinstance(config.get("form"), Mapping)
            or not isinstance(config.get("adv"), Mapping) or mode not in DIAGRAM_MODES):
        raise ValidationError("Settings-only Session requires an active Web configuration and mode.")
    from .session_resources import canonical_resource_ids

    if set(session["resources"]) - canonical_resource_ids(bindings):
        raise ValidationError("Settings-only Session contains an unbound resource.")


def _validate_web_file_bindings(session: Mapping[str, Any]) -> None:
    """Admit the Web draft inventory independently of committed replay."""
    web_files = session.get("webFiles")
    if not isinstance(web_files, Mapping) or "bindings" not in web_files:
        return
    bindings = web_files["bindings"]
    if not isinstance(bindings, Mapping):
        raise ValidationError("Session webFiles.bindings must be an object.")
    schema = bindings.get("schema")
    if isinstance(schema, bool) or schema not in (1, 2):
        raise ValidationError("Unsupported Web file binding schema.")
    current = schema == 2
    if current and (session.get("version") not in (41, 42, 44, CURRENT_SESSION_VERSION) or "c_gb" not in bindings):
        raise ValidationError(
            f"Web binding schema 2 requires session 41, 42, 44, or {CURRENT_SESSION_VERSION} and c_gb."
        )
    resources = session.get("resources", {})

    def metadata(value: Mapping[str, Any]) -> None:
        modified = value.get("lastModified")
        if (not isinstance(value.get("name"), str)
                or not isinstance(value.get("type"), str)
                or isinstance(modified, bool)
                or not isinstance(modified, (int, float))
                or not math.isfinite(modified) or modified < 0):
            raise ValidationError("Invalid Web file binding metadata.")

    def leaf(value: Any) -> None:
        if not isinstance(value, Mapping) or "kind" in value or "components" in value:
            raise ValidationError("A component must be an ordinary Web file binding.")
        resource_id = value.get("resourceId")
        if not isinstance(resource_id, str) or not resource_id.strip():
            raise ValidationError("A Web file binding requires a resourceId.")
        if current and resource_id != resource_id.strip():
            raise ValidationError("A Web file binding requires a canonical resourceId.")
        resource_id = resource_id.strip()
        if resource_id not in resources:
            raise ValidationError(f"Web file binding references a missing resource: {resource_id}.")
        if current:
            if set(value) != {"resourceId", "name", "type", "lastModified"}:
                raise ValidationError("Invalid ordinary Web file binding fields.")
            metadata(value)
            descriptor = resources[resource_id]
            if (not isinstance(descriptor, Mapping)
                    or descriptor.get("encoding") not in ("base64", DEPTH_FILE_ENCODING)
                    or (descriptor.get("encoding") == "base64" and not isinstance(descriptor.get("data"), str))
                    or (descriptor.get("encoding") == DEPTH_FILE_ENCODING and not isinstance(descriptor.get("data"), Mapping))
                    or isinstance(descriptor.get("size"), bool)
                    or not isinstance(descriptor.get("size"), int) or descriptor["size"] < 0):
                raise ValidationError("Invalid Web file binding resource payload.")

    def value(binding: Any, composite_allowed: bool = False) -> None:
        if binding is None:
            return
        if isinstance(binding, list):
            for item in binding:
                value(item)
        elif isinstance(binding, Mapping) and "kind" in binding:
            if not current or not composite_allowed or binding["kind"] != "composite":
                raise ValidationError("Unsupported Web composite file binding.")
            if set(binding) != {"kind", "components", "name", "type", "lastModified"}:
                raise ValidationError("Invalid composite binding fields.")
            metadata(binding)
            components = binding["components"]
            if not isinstance(components, list) or len(components) < 2:
                raise ValidationError("A composite binding requires at least two components.")
            for component in components:
                leaf(component)
                if resources[component["resourceId"]]["encoding"] != "base64":
                    raise ValidationError("Composite components require base64 resources.")
        else:
            leaf(binding)

    slots = {
        "c_gb", "c_gff", "c_fasta", "c_depth", "c_conservation_blasts",
        "c_conservation_fastas", "c_conservation_sequence_sources", "d_color",
        "t_color", "blacklist", "whitelist", "qualifier_priority",
    }
    if current and set(bindings) - slots - {
        "schema", "c_conservation_blasts_source", "linearSeqs", "linearComparisons",
    }:
        raise ValidationError("Unknown Web binding inventory field.")
    for slot in slots:
        value(bindings.get(slot), slot == "c_gb")
    for sequence in bindings.get("linearSeqs", []):
        if isinstance(sequence, Mapping):
            for slot in ("gb", "gff", "fasta", "depth", "blast"):
                value(sequence.get(slot))
    for slot in ("linearComparisons", "linearCanonicalComparisons"):
        entries = bindings.get(slot)
        if isinstance(entries, list):
            for row in entries:
                if isinstance(row, Mapping):
                    value(row.get("file"))


def _validate_display_placement_drafts(session: Mapping[str, Any]) -> None:
    """Validate a flat Session 41-44 draft's editable intent apart from the committed request."""

    config = session.get("config", {})
    if not isinstance(config, Mapping):
        return
    _validate_display_placement_draft_config(
        config,
        mode=(session.get("renderRequest") or {}).get("mode") or session.get("ui", {}).get("mode"),
        scoped=True,
        current=session.get("version", 0) >= TYPED_DRAFT_SESSION_MIN_VERSION,
    )


def _validate_display_placement_draft_config(
    config: Mapping[str, Any],
    *,
    mode: str,
    scoped: bool,
    current: bool,
) -> None:
    """Validate the record display and feature placement drafts of one Web draft.

    A flat Session 41-44 draft names each record display row's mode
    (``scoped``); a Session 46 mode slice holds the rows of its own mode only.
    Placement rows are keyed by ``[recordKey, biologicalFeatureId]`` in both.
    """
    from .api.requests import RecordDisplayOptions
    from .features.placement import FeaturePlacementOverride, normalize_feature_placements

    drafts = config.get("recordDisplayDrafts", [])
    if not isinstance(drafts, list):
        raise ValidationError("config.recordDisplayDrafts must be an array.")
    expected_fields = {
        "sourceUid", "selector", "recordId", "topologyOverride", "startCoordinate",
    }
    if scoped:
        expected_fields.add("scope")
    if current:
        expected_fields |= {"reverseComplementOverride", "anchorIntent"}
    keys = set()
    for row in drafts:
        if not isinstance(row, Mapping) or set(row) != expected_fields:
            raise ValidationError("Invalid record display draft fields.")
        if (scoped and row["scope"] not in {"circular", "linear"}) or any(
            not isinstance(row[name], str) or "\0" in row[name]
            for name in ("sourceUid", "selector", "recordId")
        ) or not row["sourceUid"] or not re.fullmatch(r"#[1-9]\d*", row["selector"]):
            raise ValidationError("Record display drafts require a source UID and exact selector.")
        RecordDisplayOptions(row["topologyOverride"], None)
        RecordDisplayOptions(None, row["startCoordinate"])
        if current:
            reverse = row["reverseComplementOverride"]
            if reverse is not None and not isinstance(reverse, bool):
                raise ValidationError(
                    "Record display reverse override must be a boolean or null."
                )
            intent = row["anchorIntent"]
            if intent is not None:
                intent_fields = {
                    "schema", "recordKey", "biologicalFeatureId", "placement",
                    "anchor", "offsetBp", "orientForward",
                }
                offset = intent.get("offsetBp") if isinstance(intent, Mapping) else None
                if (
                    not isinstance(intent, Mapping)
                    or set(intent) != intent_fields
                    or intent.get("schema") != 1
                    or isinstance(intent.get("schema"), bool)
                    or not isinstance(intent.get("recordKey"), str)
                    or not intent["recordKey"]
                    or "\0" in intent["recordKey"]
                    or not isinstance(intent.get("biologicalFeatureId"), str)
                    or not intent["biologicalFeatureId"]
                    or "\0" in intent["biologicalFeatureId"]
                    or intent.get("placement") not in {"anchor", "feature-end"}
                    or (
                        intent.get("anchor") not in {"five-prime", "midpoint", "three-prime"}
                        if intent.get("placement") == "anchor"
                        else intent.get("anchor") is not None
                    )
                    or not isinstance(offset, int)
                    or isinstance(offset, bool)
                    or abs(offset) > 9_007_199_254_740_991
                    or not isinstance(intent.get("orientForward"), bool)
                ):
                    raise ValidationError("Invalid record display anchor intent.")
        key = (row.get("scope"), row["sourceUid"], row["selector"])
        if key in keys:
            raise ValidationError("Duplicate record display draft identity.")
        keys.add(key)
    placements = config.get("featurePlacementOverrides", {})
    if not isinstance(placements, Mapping):
        raise ValidationError("config.featurePlacementOverrides must be an object.")
    rows = normalize_feature_placements(
        tuple(FeaturePlacementOverride.from_mapping(row) for row in placements.values())
    )
    if set(placements) != {_draft_pair_key(row.record_key, row.biological_feature_id) for row in rows}:
        raise ValidationError("Feature placement draft keys must encode their exact identity as a JSON pair.")
    for row in rows:
        row.target.validate_mode(mode)
    adv = config.get("adv", {})
    if not isinstance(adv, Mapping):
        raise ValidationError("config.adv must be an object.")
    tolerance = adv.get("feature_overlap_tolerance_bp", 0)
    if isinstance(tolerance, bool) or not isinstance(tolerance, int) or tolerance < 0:
        raise ValidationError("Feature overlap tolerance must be a non-negative integer.")


def _draft_identity_key(scope: str, record_key: str, biological_feature_id: str) -> str:
    return json.dumps([scope, record_key, biological_feature_id], ensure_ascii=False, separators=(",", ":"))


def _draft_pair_key(record_key: str, biological_feature_id: str) -> str:
    return json.dumps([record_key, biological_feature_id], ensure_ascii=False, separators=(",", ":"))


def _validate_feature_override_drafts(features: Mapping[str, Any]) -> None:
    """Validate a Session 46 mode slice's per-feature edit drafts keyed by source identity.

    A draft row is a request ``featureOverrides`` row plus the Web-only
    ``labelSourceText``; a row may hold only that source text (a bulk label
    edit's target). The key encodes the identity as a JSON pair, as placement
    drafts do.
    """
    from .features.overrides import FeatureOverride

    invalid = _SESSION_FIELDS_INVALID
    drafts = features.get("featureOverrides", {})
    if not isinstance(drafts, Mapping):
        raise ValidationError("features.featureOverrides must be an object.", diagnostic=invalid)
    for key, row in drafts.items():
        if not isinstance(row, Mapping) or set(row) != FEATURE_OVERRIDE_DRAFT_FIELDS:
            raise ValidationError("Invalid feature override draft fields.", diagnostic=invalid)
        source_text = row["labelSourceText"]
        if source_text is not None and (
            not isinstance(source_text, str) or not source_text or "\0" in source_text
        ):
            raise ValidationError(
                "Feature override labelSourceText must be text or null.", diagnostic=invalid
            )
        edits = {name: row[name] for name in ("featureVisibility", "labelVisibility", "labelText")}
        if any(value is not None for value in edits.values()):
            override = FeatureOverride.from_mapping(
                {"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"], **edits}
            )
            identity = (override.record_key, override.biological_feature_id)
        elif source_text is None:
            raise ValidationError(
                "A feature override draft must hold an edit or a label source text.",
                diagnostic=invalid,
            )
        else:
            identity = (row["recordKey"], row["biologicalFeatureId"])
            if any(not isinstance(value, str) or not value or "\0" in value for value in identity):
                raise ValidationError(
                    "Feature override drafts require a record key and feature ID.",
                    diagnostic=invalid,
                )
        if key != _draft_pair_key(*identity):
            raise ValidationError(
                "Feature override draft keys must encode their identity as a JSON pair.",
                diagnostic=invalid,
            )


def _mode_slice(session: Mapping[str, Any], mode: object) -> Mapping[str, Any] | None:
    modes = session.get("modes")
    mode_slice = modes.get(mode) if isinstance(modes, Mapping) and isinstance(mode, str) else None
    return mode_slice if isinstance(mode_slice, Mapping) else None


def _mode_slice_config(session: Mapping[str, Any], mode: object) -> Mapping[str, Any] | None:
    mode_slice = _mode_slice(session, mode)
    config = mode_slice.get("config") if mode_slice is not None else None
    return config if isinstance(config, Mapping) else None


def _session_draft_configs(
    session: Mapping[str, Any],
) -> list[tuple[DiagramMode | None, Mapping[str, Any]]]:
    """Each Web draft config of a Session: its mode slices' (46), else the flat one."""

    if session.get("version", 0) >= MODE_SCOPED_SESSION_MIN_VERSION:
        return [
            (mode, config)
            for mode in DIAGRAM_MODES
            if (config := _mode_slice_config(session, mode)) is not None
        ]
    config = session.get("config")
    return [(None, config)] if isinstance(config, Mapping) else []


# The top-level homes that Session 46 moved into ``modes``: the flat draft,
# and the Legend edit, per-feature stroke, and per-mode ``ui`` keys.
_RETIRED_MODE_SCOPED_FIELDS: dict[str, frozenset[str]] = {
    "": FLAT_DRAFT_TOP_LEVEL_FIELDS,
    "editorState": frozenset({"featureStrokes"}),
    **{
        domain: frozenset(row.path for row in MODE_SCOPED_SETTINGS if row.domain == domain)
        for domain in ("editorState.legend", "ui")
    },
}


def _validate_mode_scoped_fields(session: Mapping[str, Any], version: int) -> None:
    """Admit ``modes`` in Session 46 only, and its draft nowhere else."""

    if version < MODE_SCOPED_SESSION_MIN_VERSION:
        for field in ("modes", "cliOptions"):
            if field in session:
                raise ValidationError(
                    f"Session version {version} cannot contain {field}.", diagnostic=_SESSION_FIELDS_INVALID
                )
        return
    if "cliOptions" in session and not isinstance(session["cliOptions"], Mapping):
        raise ValidationError("Session cliOptions must be an object.", diagnostic=_SESSION_FIELDS_INVALID)
    retired = sorted(
        f"{domain}.{field}" if domain else field
        for domain, fields in _RETIRED_MODE_SCOPED_FIELDS.items()
        if isinstance(container := _container_at(session, domain), Mapping)
        for field in fields & set(container)
    )
    if retired:
        raise ValidationError(
            f"Session version {version} keeps the Web draft of each mode in modes; "
            f"it cannot contain {', '.join(retired)}.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )
    ui = session.get("ui")
    execution = ui.get("losatExecution") if isinstance(ui, Mapping) else None
    if execution is not None and (
        not isinstance(execution, Mapping) or set(execution) - set(LOSAT_EXECUTION_FIELDS)
    ):
        raise ValidationError(
            "ui.losatExecution holds only the LOSAT execution settings.", diagnostic=_SESSION_FIELDS_INVALID
        )


def _container_at(source: Mapping[str, Any], domain: str) -> object:
    current: object = source
    for part in domain.split(".") if domain else ():
        current = current.get(part) if isinstance(current, Mapping) else None
    return current


def _validate_mode_slices(session: Mapping[str, Any]) -> None:
    """Validate each Session 46 mode slice against its own mode.

    A slice holds only registry fields (MODE_SCOPED_SETTINGS); a missing field
    is that mode's default, and a missing slice is all defaults. Each slice
    runs the draft validators with its own mode, so a lane placement or a
    GUI-unmanaged override of one mode never reaches the other (OV-106).
    """
    from .web_support.config_overrides import validate_and_project_web_config_overrides

    modes = session.get("modes")
    if modes is None:
        return
    if not isinstance(modes, Mapping) or set(modes) - set(DIAGRAM_MODES):
        raise ValidationError(
            "Session modes holds a circular and a linear slice only.", diagnostic=_SESSION_FIELDS_INVALID
        )
    for mode, mode_slice in modes.items():
        _validate_mode_slice_fields(mode_slice, mode)
        config = mode_slice.get("config", {})
        _validate_display_placement_draft_config(config, mode=mode, scoped=False, current=True)
        _validate_feature_override_drafts(mode_slice.get("features", {}))
        overrides = config.get("unmanagedConfigOverrides")
        if overrides:
            validate_and_project_web_config_overrides(mode=mode, overrides=overrides)


def _validate_mode_slice_fields(value: object, mode: str, container: str = "") -> None:
    """Reject a slice field that the registry does not name."""

    label = f"modes.{mode}" + (f".{container}" if container else "")
    if not isinstance(value, Mapping):
        raise ValidationError(f"Session {label} must be an object.", diagnostic=_SESSION_FIELDS_INVALID)
    unknown = sorted(str(field) for field in set(value) - SLICE_CONTAINERS[container])
    if unknown:
        raise ValidationError(
            f"Session {label} cannot contain {', '.join(unknown)}.", diagnostic=_SESSION_FIELDS_INVALID
        )
    for field, child in value.items():
        child_container = f"{container}.{field}" if container else field
        if child_container in SLICE_CONTAINERS:
            _validate_mode_slice_fields(child, mode, child_container)


def _validate_current_retired_active_config_paths(
    session: Mapping[str, Any],
) -> None:
    """Reject the two retired v40 circular-track draft paths."""

    config = session.get("config")
    if not isinstance(config, Mapping):
        return
    advanced = config.get("adv")
    if not isinstance(advanced, Mapping):
        return
    for field in ("cli_circular_track_order", "cli_circular_track_slots"):
        if field in advanced:
            raise ValidationError(
                f"Session version 40 cannot contain config.adv.{field}."
            )


def _validate_current_comparison_authority(
    session: Mapping[str, Any],
) -> None:
    """Reject retired v40 comparison fields and validate an optional Web draft."""

    config_value = (
        # Session 46 keeps the comparison draft in the Linear slice.
        _mode_slice_config(session, "linear")
        if session.get("version", 0) >= MODE_SCOPED_SESSION_MIN_VERSION
        else session.get("config")
    )
    config = config_value if isinstance(config_value, Mapping) else {}
    ui_value = session.get("ui")
    ui = ui_value if isinstance(ui_value, Mapping) else {}
    web_files_value = session.get("webFiles")
    web_files = web_files_value if isinstance(web_files_value, Mapping) else {}
    bindings_value = web_files.get("bindings")
    bindings = bindings_value if isinstance(bindings_value, Mapping) else {}
    adv_value = config.get("adv")
    adv = adv_value if isinstance(adv_value, Mapping) else {}
    layout_value = config.get("linearRecordLayout")
    layout = layout_value if isinstance(layout_value, Mapping) else None

    if "blastSource" in config or "blastSource" in adv:
        raise ValidationError(
            "Session version 40 cannot contain retired blastSource state."
        )
    if "blastSource" in ui:
        raise ValidationError(
            "Session version 40 cannot contain ui.blastSource."
        )
    if layout is not None and "comparisons" in layout:
        raise ValidationError(
            "Session version 40 cannot contain "
            "config.linearRecordLayout.comparisons."
        )

    linear_sequences = bindings.get("linearSeqs")
    if isinstance(linear_sequences, list):
        for sequence in linear_sequences:
            if not isinstance(sequence, Mapping):
                continue
            if "blast" in sequence:
                raise ValidationError(
                    "Session version 40 cannot contain per-record BLAST bindings."
                )
            if "losat_filename" in sequence:
                raise ValidationError(
                    "Session version 40 cannot contain per-record LOSAT filenames."
                )
    linear_metadata = web_files.get("linearRecordMetadata")
    if isinstance(linear_metadata, list):
        for metadata in linear_metadata:
            if not isinstance(metadata, Mapping):
                continue
            if "losatFilename" in metadata or "losat_filename" in metadata:
                raise ValidationError(
                    "Session version 40 cannot contain per-record LOSAT filenames."
                )
    if "linearCanonicalComparisons" in bindings:
        raise ValidationError(
            "Session version 40 cannot bind comparison artifacts outside "
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
        bound_ids.add(edge_id)
    for edge in edges:
        if edge.get("fileActive") and str(edge.get("id") or "") not in bound_ids:
            raise ValidationError(
                "Active comparison file is missing its Web file binding: "
                f"{edge.get('id')}."
            )


def _validate_alignment_reset_receipt(session: Mapping[str, Any]) -> None:
    """Admit optional artifact restoration history without interpreting render policy."""
    editor = session.get("editorState")
    request = session.get("renderRequest")
    plan = request.get("layout", {}).get("similarityAlignment") if isinstance(request, Mapping) else None
    schema = request.get("schema") if isinstance(request, Mapping) else None
    # One rule with services/session-authority.js: request schema 8 and later.
    if (isinstance(schema, int) and not isinstance(schema, bool) and schema >= 8 and plan
            and isinstance(editor, Mapping) and "alignmentResetReceipt" not in editor):
        raise ValidationError("Current alignment Session requires editorState.alignmentResetReceipt.")
    receipt = editor.get("alignmentResetReceipt") if isinstance(editor, Mapping) else None
    if receipt is None:
        return
    def invalid() -> NoReturn:
        raise ValidationError("Alignment reset receipt is malformed or stale.")

    if (not isinstance(receipt, Mapping)
            or set(receipt) != {"binding", "directions", "referenceDeltaX"}
            or not isinstance(receipt.get("binding"), str)
            or not re.fullmatch(r"[0-9a-f]{64}", receipt["binding"])
            or not isinstance(receipt.get("directions"), list)
            or not isinstance(plan, Mapping) or not isinstance(request, Mapping)
            or request.get("mode") != "linear"):
        invalid()
    eligible = {item["recordKey"] for item in plan["records"] if item["status"] != "skipped"}
    keys: set[str] = set()
    for delta in receipt["directions"]:
        if (not isinstance(delta, Mapping) or set(delta) != {"recordKey", "before", "after"}
                or delta["recordKey"] not in eligible or delta["recordKey"] in keys
                or not isinstance(delta["before"], bool) or not isinstance(delta["after"], bool)
                or delta["before"] == delta["after"]):
            invalid()
        keys.add(delta["recordKey"])
    delta = receipt["referenceDeltaX"]
    if delta is not None and (not isinstance(delta, Mapping)
            or set(delta) != {"recordKey", "deltaX"}
            or delta["recordKey"] != plan["reference"]["recordKey"]
            or isinstance(delta["deltaX"], bool) or not isinstance(delta["deltaX"], (int, float))
            or not math.isfinite(delta["deltaX"]) or delta["deltaX"] == 0):
        invalid()
    fingerprints: dict[str, str] = {}

    def source_identity(source: Mapping[str, Any]) -> dict[str, str]:
        identity = {}
        for key, value in source.items():
            if key == "kind":
                identity[key] = value
            else:
                if value not in fingerprints:
                    resource = session.get("resources", {}).get(value)
                    if not isinstance(resource, Mapping) or resource.get("encoding") != "base64":
                        invalid()
                    try:
                        content = base64.b64decode(resource["data"], validate=True)
                    except (ValueError, TypeError, KeyError) as exc:
                        raise ValidationError("Alignment reset receipt resource is invalid.") from exc
                    if len(content) != resource.get("size"):
                        invalid()
                    fingerprints[value] = hashlib.sha256(content).hexdigest()
                identity[key] = fingerprints[value]
        return identity

    records = [{
        "recordKey": record["recordKey"], "source": source_identity(record["source"]),
        "selector": record["selector"],
        "region": ({"start": record["region"]["start"], "end": record["region"]["end"]}
                   if record.get("region") else None), "display": record["display"],
    } for record in request["records"]]
    binding = {"plan": {**plan, "records": sorted(plan["records"], key=lambda item: item["recordKey"].encode("utf-16be"))},
               "records": sorted(records, key=lambda item: item["recordKey"].encode("utf-16be"))}
    digest = hashlib.sha256(json.dumps(binding, ensure_ascii=False, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()
    if receipt["binding"] != digest:
        raise ValidationError("Alignment reset receipt source or plan binding changed.")


def _validate_current_top_level_fields(session: Mapping[str, Any]) -> None:
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


def other_mode_result_view(session: Mapping[str, Any]) -> dict[str, Any]:
    """The Session as its ``otherModeResult`` set would be at the top level."""
    other = session.get("otherModeResult")
    if not isinstance(other, Mapping):
        raise ValidationError("Session has no otherModeResult.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    return {**session, **other}


def _validate_other_mode_result(session: Mapping[str, Any], version: int) -> None:
    """Admit the other diagram mode's committed Result set (Session 46)."""
    if "otherModeResult" not in session:
        return
    other = session["otherModeResult"]
    if not isinstance(other, Mapping) or set(other) - OTHER_MODE_RESULT_FIELDS:
        raise ValidationError("Session otherModeResult must contain only a committed Result set.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    top_request = session.get("renderRequest")
    request = other.get("renderRequest")
    if (
        not isinstance(top_request, Mapping)
        or not isinstance(request, Mapping)
        or request.get("mode") not in ("circular", "linear")
        or request.get("mode") == top_request.get("mode")
        or request.get("schema") != top_request.get("schema")
    ):
        raise ValidationError("Session otherModeResult requires a committed request of the other mode.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    view = other_mode_result_view(session)
    _validate_current_feature_catalog_authority(view, version)
    if not other.get("results"):
        raise ValidationError("Session otherModeResult requires a Result.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    _validate_alignment_reset_receipt(view)
    from .session_resources import canonical_resource_ids

    missing = sorted(canonical_resource_ids(request) - set(session.get("resources") or {}))
    if missing:
        raise ValidationError(
            "Session otherModeResult names missing resource(s): " + ", ".join(missing) + ".",
            diagnostic=_OTHER_MODE_RESULT_INVALID,
        )
    editor = other.get("editorState")
    legend = editor.get("legend", {}) if isinstance(editor, Mapping) else None
    ui = other.get("ui", {})
    _validate_other_mode_run_metadata(other.get("runMetadata", {}), other["results"])
    if (
        not isinstance(editor, Mapping)
        or set(editor) - OTHER_MODE_RESULT_EDITOR_FIELDS
        or not isinstance(legend, Mapping)
        or set(legend) - OTHER_MODE_RESULT_LEGEND_FIELDS
        or not isinstance(ui, Mapping)
        or set(ui) - OTHER_MODE_RESULT_UI_FIELDS
    ):
        raise ValidationError("Session otherModeResult editorState and ui hold only that Result set's fields.", diagnostic=_OTHER_MODE_RESULT_INVALID)


_ANNOTATION_WARNING_FIELDS = frozenset(
    {"code", "setId", "annotationId", "recordId", "recordIndex", "missingCount", "message", "resultIndex", "resultName"}
)
_ANNOTATION_WARNING_CODES = frozenset(
    {"feature_selector_unmatched", "empty_span", "out_of_bounds_skipped", "out_of_bounds_clipped"}
)
_COMPARISON_WARNING_FIELDS = frozenset(
    {
        "code", "queryRecordIndex", "subjectRecordIndex", "queryRecordId", "subjectRecordId",
        "rowCount", "exampleIds", "message", "resultIndex", "resultName",
    }
)
_FEATURE_IDENTITY_NOTICE_FIELDS = frozenset({"biologicalFeatureId", "kinds", "recordKey", "resultIndex", "status"})
_FEATURE_IDENTITY_NOTICE_STATUSES = frozenset({"crop_excluded", "absent", "unresolved"})
_FEATURE_IDENTITY_NOTICE_KINDS = frozenset({"placement", "feature_visibility", "label_visibility", "label_text"})
OTHER_MODE_RESULT_RUN_METADATA_FIELDS = frozenset(
    {"trackSlotGeometry", "annotationWarnings", "featureIdentityNotices", "comparisonWarnings"}
)


def _index(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _validate_other_mode_run_metadata(run_metadata: Any, results: list[Any]) -> None:
    """The other set's ``runMetadata``, checked as the Web reader checks it."""
    if not isinstance(run_metadata, Mapping) or set(run_metadata) - OTHER_MODE_RESULT_RUN_METADATA_FIELDS:
        raise ValidationError("Session otherModeResult.runMetadata holds only that Result set's metadata.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    for key in ("annotationWarnings", "comparisonWarnings", "featureIdentityNotices"):
        if key in run_metadata and not isinstance(run_metadata[key], list):
            raise ValidationError(f"Session otherModeResult.runMetadata.{key} must be an array.", diagnostic=_OTHER_MODE_RESULT_INVALID)

    def result_named(index: Any, name: Any) -> bool:
        return (
            _index(index) and index < len(results) and isinstance(results[index], Mapping)
            and results[index].get("name") == name
        )

    for warning in run_metadata.get("annotationWarnings", []):
        if (
            not isinstance(warning, Mapping) or set(warning) != _ANNOTATION_WARNING_FIELDS
            or warning["code"] not in _ANNOTATION_WARNING_CODES
            or not all(isinstance(warning[key], str) for key in ("setId", "annotationId", "recordId", "message", "resultName"))
            or not warning["setId"] or not warning["annotationId"]
            or not all(_index(warning[key]) for key in ("recordIndex", "missingCount", "resultIndex"))
            or (warning["code"] == "feature_selector_unmatched") != (warning["missingCount"] > 0)
            or not result_named(warning["resultIndex"], warning["resultName"])
        ):
            raise ValidationError("Annotation warnings do not match the successful Result metadata schema.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    for warning in run_metadata.get("comparisonWarnings", []):
        if (
            not isinstance(warning, Mapping) or set(warning) != _COMPARISON_WARNING_FIELDS
            or warning["code"] != "comparison_record_id_unmatched"
            or not all(isinstance(warning[key], str) for key in ("queryRecordId", "subjectRecordId", "message", "resultName"))
            or not warning["message"]
            or not isinstance(warning["exampleIds"], list)
            or not all(isinstance(item, str) for item in warning["exampleIds"])
            or not all(_index(warning[key]) for key in ("queryRecordIndex", "subjectRecordIndex", "resultIndex"))
            or not _index(warning["rowCount"]) or warning["rowCount"] < 1
            or not result_named(warning["resultIndex"], warning["resultName"])
        ):
            raise ValidationError("Comparison warnings do not match the successful Result metadata schema.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    for notice in run_metadata.get("featureIdentityNotices", []):
        if (
            not isinstance(notice, Mapping) or set(notice) != _FEATURE_IDENTITY_NOTICE_FIELDS
            or not isinstance(notice["recordKey"], str) or not notice["recordKey"]
            or not isinstance(notice["biologicalFeatureId"], str) or not notice["biologicalFeatureId"]
            or notice["status"] not in _FEATURE_IDENTITY_NOTICE_STATUSES
            or not isinstance(notice["kinds"], list) or not notice["kinds"]
            or not all(kind in _FEATURE_IDENTITY_NOTICE_KINDS for kind in notice["kinds"])
            or not _index(notice["resultIndex"]) or notice["resultIndex"] >= len(results)
        ):
            raise ValidationError("runMetadata.featureIdentityNotices contains an invalid notice.", diagnostic=_OTHER_MODE_RESULT_INVALID)


def _validate_current_feature_catalog_authority(
    session: Mapping[str, Any],
    version: int,
) -> None:
    """Require the version-owned catalog and reject duplicated payloads."""

    catalog_schema = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION.get(version, 3)

    features = session.get("features")
    if isinstance(features, Mapping):
        duplicated = CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS & set(features)
        if duplicated:
            raise ValidationError(
                "Session version 40 cannot contain derived feature payloads in "
                f"features: {', '.join(sorted(duplicated))}."
            )
    elif "features" in session:
        raise ValidationError("Session features must be an object when present.")

    orthogroup_state = session.get("orthogroupState")
    if isinstance(orthogroup_state, Mapping):
        if "groups" in orthogroup_state:
            raise ValidationError(
                "Session version 40 cannot contain derived orthogroupState.groups."
            )
    elif "orthogroupState" in session:
        raise ValidationError(
            "Session orthogroupState must be an object when present."
        )

    results = session.get("results")
    if "results" not in session or not isinstance(results, list):
        raise ValidationError("Session version 40 requires a results array.")
    editor_state = session.get("editorState")
    if (
        not isinstance(editor_state, Mapping)
        or "featureCatalog" not in editor_state
    ):
        raise ValidationError(
            "Session version 40 "
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
            or catalog.get("schema") != catalog_schema
            or catalog.get("items") != []
        ):
            raise ValidationError(
                "An empty Result set requires an empty feature catalog."
            )
        return
    if not isinstance(catalog, Mapping):
        raise ValidationError(
            "Session version 40 Results require editorState.featureCatalog."
        )
    items = catalog.get("items")
    if (
        catalog.get("schema") != catalog_schema
        or not isinstance(items, list)
        or len(items) != len(results)
    ):
        raise ValidationError(
            "Session feature catalog must contain one version-compatible item per Result."
        )

    from .web_support.feature_catalog import (
        promote_legacy_feature_catalog,
        select_feature_catalog_item,
    )

    result_names: list[str] = []
    for result in results:
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
        result_names.append(result_name)
    try:
        # A schema 3 catalog is read as promoted, which may infer what the
        # current rules require (OV-269); its promoted form is the one validated.
        if catalog_schema == 3:
            catalog = promote_legacy_feature_catalog(catalog)
            catalog_schema = CURRENT_FEATURE_CATALOG_SCHEMA
        for result_index, result_name in enumerate(result_names):
            select_feature_catalog_item(
                catalog,
                result_index=result_index,
                result_name=result_name,
                expected_schema=catalog_schema,
            )
    except GbdrawError as exc:
        raise ValidationError(str(exc)) from exc


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

    session_version = session.get("version")
    for _, config in _session_draft_configs(session):
        validate_current_web_state_field_names(
            config,
            include_linear_label_visibility=(
                isinstance(session_version, int)
                and session_version >= TYPED_DRAFT_SESSION_MIN_VERSION
            ),
        )
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
                validated_manifest,
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


_PLACEMENT_LANE_MODES = {
    "outward": "circular",
    "inward": "circular",
    "above": "linear",
    "below": "linear",
}


def _migrate_session_feature_placements(placements: object) -> object:
    """Scope the Session 41-44 placement drafts keyed by [recordKey, featureId].

    Such a row reached every request with its record key, so a Main row is kept
    for both modes and a lane row for the mode of its side. A row whose key does
    not encode its identity is kept as is, and the draft check rejects it. This
    is the twin of ``migrateSessionFeaturePlacements`` in the Web
    ``feature-edit-migration.js``; ``tests/fixtures/feature-placement-migration.json``
    pins both.
    """

    if not isinstance(placements, Mapping):
        return placements
    migrated: dict[str, Any] = {}
    for key, row in placements.items():
        fields = row if isinstance(row, Mapping) else {}
        target = fields.get("placement")
        target = target if isinstance(target, Mapping) else {}
        encoded = json.dumps(
            [fields.get("recordKey"), fields.get("biologicalFeatureId")],
            ensure_ascii=False,
            separators=(",", ":"),
        )
        modes: list[str] = []
        if key == encoded:
            if target.get("kind") == "main":
                modes = ["circular", "linear"]
            elif isinstance(target.get("side"), str) and target["side"] in _PLACEMENT_LANE_MODES:
                modes = [_PLACEMENT_LANE_MODES[target["side"]]]
        if not modes:
            migrated[key] = row
        for scope in modes:
            scoped = {"scope": scope, **fields}
            migrated[_feature_identity_key_of(scoped) or key] = scoped
    return migrated


# The Web readers' ``String(value ?? '')`` and ``.trim()`` (ECMAScript white
# space and line terminators), so the twins below read text as the Web does.
_JS_WHITESPACE = (
    "\t\n\v\f\r \u00a0\u1680\u2000\u2001\u2002\u2003\u2004\u2005\u2006"
    "\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000\ufeff"
)


def _js_string(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return str(int(value)) if value.is_integer() and abs(value) < 1e21 else repr(value)
    if isinstance(value, list):
        return ",".join(_js_string(item) for item in value)
    return "[object Object]"


def _js_text(value: object) -> str:
    return _js_string(value).strip(_JS_WHITESPACE)


def _feature_identity_key(scope: object, record_key: object, feature_id: object) -> str:
    """The draft key of one feature in one mode, or ``""`` for an invalid identity.

    The twin of ``featureIdentityKey`` in the Web ``feature-placement.js``.
    """

    if (
        scope in DRAFT_SCOPES
        and isinstance(record_key, str)
        and isinstance(feature_id, str)
        and record_key.strip(_JS_WHITESPACE)
        and feature_id.strip(_JS_WHITESPACE)
        and "\0" not in record_key + feature_id
    ):
        return _draft_identity_key(str(scope), record_key, feature_id)
    return ""


def _feature_identity_key_of(row: Mapping[str, Any]) -> str:
    record_key = row.get("record_key")
    feature_id = row.get("biological_feature_id")
    return _feature_identity_key(
        row.get("scope"),
        row.get("recordKey") if record_key is None else record_key,
        row.get("biologicalFeatureId") if feature_id is None else feature_id,
    )


# Session 44 and older kept per-feature edits in four maps keyed by rendered SVG
# ID (RETIRED_RENDERED_ID_FEATURE_FIELDS). A rendered ID is
# ``<hash>[_record_<n>][__instance_<s>_<digest>]``, optionally with a part suffix.
_RENDERED_PART_SUFFIX = re.compile(r"__(?:part|line)[0-9]+\Z")
_RENDERED_INSTANCE_SUFFIX = re.compile(r"__instance_([A-Za-z0-9_.-]+?)_[0-9a-f]{16}\Z")
_RENDERED_RECORD_INSTANCE = re.compile(r"record_([1-9][0-9]*)")
_RENDERED_SOURCE_INSTANCE = re.compile(r"0|[1-9][0-9]*")
_RENDERED_LINEAR_RECORD = re.compile(r"([^\n\r\u2028\u2029]*)_record_([1-9][0-9]*)\Z")
_FEATURE_VISIBILITY_EDIT_VALUES = {
    "on": "on",
    "off": "off",
    "exclude_matching": "exclude_matching",
    "suppress": "exclude_matching",
}
_FEATURE_EDIT_DRAFT_FIELDS = ("featureVisibility", "labelVisibility", "labelText", "labelSourceText")
_JS_MAX_SAFE_INTEGER = 2**53 - 1


def _rendered_id_without_part(rendered_id: object) -> str:
    return _RENDERED_PART_SUFFIX.sub("", _js_text(rendered_id), count=1)


def _parse_rendered_id(rendered_id: object) -> tuple[str, int | None, int | None]:
    """The stable hash, one-based record position, and source index of a rendered ID."""

    rest = _rendered_id_without_part(rendered_id)
    record_ordinal: int | None = None
    source_index: int | None = None
    match = _RENDERED_INSTANCE_SUFFIX.search(rest)
    while match:
        instance = match.group(1)
        record = _RENDERED_RECORD_INSTANCE.fullmatch(instance)
        if record:
            if record_ordinal is None:
                record_ordinal = int(record.group(1))
        elif _RENDERED_SOURCE_INSTANCE.fullmatch(instance) and source_index is None:
            source_index = int(instance)
        rest = rest[: match.start()]
        match = _RENDERED_INSTANCE_SUFFIX.search(rest)
    linear_record = _RENDERED_LINEAR_RECORD.match(rest)
    if linear_record:
        rest = linear_record.group(1)
        if record_ordinal is None:
            record_ordinal = int(linear_record.group(2))
    return rest, record_ordinal, source_index


def _safe_integer(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    if isinstance(value, int) and abs(value) <= _JS_MAX_SAFE_INTEGER:
        return value
    return None


def _mapping_list(value: object) -> list[Mapping[str, Any]]:
    """The entries of a JSON array, each read as an object (a non-object reads as empty)."""

    return [item if isinstance(item, Mapping) else {} for item in value] if isinstance(value, list) else []


@dataclass(frozen=True)
class FeatureEditMigration:
    """The Session 46 ``features`` of an older Session and what Load would report."""

    features: dict[str, Any]
    dropped_count: int
    narrowed_visibility_count: int


# An identity of the feature index: the draft key, the hash a rendered ID
# names, the source feature index, and the one-based record position.
_IndexedFeature = tuple[str, str, int | None, int]
_LINEAR_RECORD_SUFFIX = re.compile(r"_record_([1-9][0-9]*)(?=__|\Z)")


def _catalog_feature_index(
    catalog: object, mode: object
) -> tuple[dict[str, dict[str, None]], list[_IndexedFeature]]:
    """The rendered IDs and source features of a saved feature catalog (schema 3-5)."""

    rendered_by_id: dict[str, dict[str, None]] = {}
    biological: list[_IndexedFeature] = []
    items = catalog.get("items") if isinstance(catalog, Mapping) else None
    for item in _mapping_list(items):
        record_keys_value = item.get("recordKeys")
        record_keys = (
            [_js_text(record_key) for record_key in record_keys_value]
            if isinstance(record_keys_value, list)
            else []
        )
        for feature in _mapping_list(item.get("features")):
            key = _feature_identity_key(
                mode, _js_text(feature.get("recordKey")), _js_text(feature.get("biologicalFeatureId"))
            )
            svg_id = _rendered_id_without_part(feature.get("svgId"))
            if key and svg_id:
                rendered_by_id.setdefault(svg_id, {})[key] = None
        for feature in _mapping_list(item.get("biologicalFeatures")):
            record_key = _js_text(feature.get("recordKey"))
            feature_id = _js_text(feature.get("biologicalFeatureId"))
            key = _feature_identity_key(mode, record_key, feature_id)
            if key:
                biological.append(
                    (
                        key,
                        _js_text(feature.get("stableFeatureId")) or feature_id,
                        _safe_integer(feature.get("sourceFeatureIndex")),
                        record_keys.index(record_key) + 1 if record_key in record_keys else 0,
                    )
                )
    return rendered_by_id, biological


def _nonnegative_integer(value: object) -> int | None:
    """``Number(value)`` when it is a non-negative safe integer (``null`` and ``''`` are not)."""

    if value is None or value == "":
        return None
    number = _js_number(value)
    if math.isfinite(number) and number.is_integer() and 0 <= number <= _JS_MAX_SAFE_INTEGER:
        return int(number)
    return None


def _first(feature: Mapping[str, Any], *fields: str) -> object:
    """The first of ``fields`` that is set (JavaScript ``a ?? b``)."""

    for field in fields:
        value = feature.get(field)
        if value is not None:
            return value
    return None


def _promoted_request_records(request: Mapping[str, Any]) -> object:
    """The request records with the cardinality Web Load promotes schema 1-5 records to (OV-149)."""

    from .session_request_codec import legacy_record_cardinality

    records = request.get("records")
    schema = request.get("schema")
    if not isinstance(records, list) or not isinstance(schema, int) or schema >= 6:
        return records
    return [
        {**record, "cardinality": legacy_record_cardinality(request.get("mode"), record).value}
        if isinstance(record, Mapping)
        else record
        for record in records
    ]


def _legacy_feature_index(legacy: object, mode: object) -> list[_IndexedFeature]:
    """The source features of a Session without a feature catalog (31-33, 39).

    Such a Session keyed an edit by the rendered ID ``<drawn hash>[_record_<n>]
    [__instance_<s>_<digest>]``: the feature's hash in the drawn (cropped,
    reverse-complemented) record and the record's position in the Result.
    ``legacy`` holds the request ``records`` and ``features``: the Session's
    source features read again with its crops and orientations (each with its
    drawn hash beside its source hash) or else its saved feature metadata,
    whose drawn hash serves only records drawn untransformed;
    ``biologicalFeatures`` are every source feature read again. A feature's
    input is its request record (Linear: ``fileIdx``; Circular: one file) and
    ``record_idx`` its record in that input. Identities are the renderer's:
    the record key, ``<recordKey>:<n>`` for each record of an ALL input with
    several records, and the source hash, ``~<source index>`` when the record
    has it twice. The twin of ``legacyIndex`` in the Web
    ``feature-edit-migration.js``.
    """

    fields = legacy if isinstance(legacy, Mapping) else {}
    linear = mode == "linear"
    records_value = fields.get("records")
    records = (
        [record if isinstance(record, Mapping) else {} for record in records_value]
        if isinstance(records_value, list)
        else []
    )

    def describe(feature: Mapping[str, Any]) -> tuple[int | None, int | None, int | None, str, str]:
        ordinal = _LINEAR_RECORD_SUFFIX.search(_js_text(_first(feature, "svg_id", "svgId")))
        svg_ordinal = int(ordinal.group(1)) if ordinal else None
        file_index = _nonnegative_integer(feature.get("fileIdx"))
        if linear:
            input_index = file_index if file_index is not None else (
                svg_ordinal - 1 if svg_ordinal else None
            )
        else:
            input_index = 0
        record_index = (
            0
            if linear and file_index is None
            else _nonnegative_integer(_first(feature, "record_idx", "recordIndex"))
        )
        source_index = next(
            (
                index
                for index in (
                    _nonnegative_integer(feature.get(field))
                    for field in ("source_feature_index", "sourceFeatureIndex", "feature_index")
                )
                if index is not None
            ),
            None,
        )
        drawn_hashes = (
            selector.get("hash") if isinstance(selector, Mapping) else None
            for selector in (feature.get("drawn_selector"), feature.get("drawnSelector"))
        )
        return (
            input_index,
            record_index,
            source_index,
            _js_text(_first(feature, "stable_feature_id", "stableFeatureId", "stable_svg_id")),
            _js_text(next((value for value in drawn_hashes if value is not None), None)),
        )

    def described(value: object) -> list[tuple[int, int, int | None, str, str]]:
        entries: list[tuple[int, int, int | None, str, str]] = []
        for input_index, record_index, source_index, source_hash, drawn_hash in map(
            describe, _mapping_list(value)
        ):
            if input_index is not None and record_index is not None and source_hash:
                entries.append((input_index, record_index, source_index, source_hash, drawn_hash))
        return entries

    listed = described(fields.get("features"))
    biological_value = fields.get("biologicalFeatures")
    sources = (
        described(biological_value)
        if isinstance(biological_value, list) and biological_value
        else listed
    )
    record_counts: dict[int, int] = {}
    for input_index, record_index, *_ in (*listed, *sources):
        record_counts[input_index] = max(record_counts.get(input_index) or 1, record_index + 1)

    def record_of(input_index: int, record_index: int) -> Mapping[str, Any] | None:
        position = input_index if linear else (record_index if len(records) > 1 else 0)
        return records[position] if position < len(records) else None

    def record_key_of(input_index: int, record_index: int) -> str:
        record = record_of(input_index, record_index)
        record_key = _js_text(record.get("recordKey")) if record is not None else ""
        if not record_key:
            return ""
        assert record is not None
        if (
            record.get("cardinality") == "all"
            and (linear or len(records) == 1)
            and record_counts.get(input_index, 0) > 1
        ):
            return f"{record_key}:{record_index + 1}"
        return record_key

    offsets: dict[int, int] = {}
    offset = 0
    for input_index in range(len(records)):
        offsets[input_index] = offset
        offset += record_counts.get(input_index) or 1
    hash_counts: dict[tuple[str, str], int] = {}
    for input_index, record_index, _, source_hash, _ in sources:
        identity = (record_key_of(input_index, record_index), source_hash)
        hash_counts[identity] = hash_counts.get(identity, 0) + 1
    biological: list[_IndexedFeature] = []
    seen: set[tuple[str, str, int]] = set()
    for input_index, record_index, source_index, source_hash, drawn_hash in listed:
        record = record_of(input_index, record_index)
        transformed = record is not None and _drawn_transformed(record)
        drawn_hash = drawn_hash or ("" if transformed else source_hash)
        record_key = record_key_of(input_index, record_index)
        if not drawn_hash or not record_key:
            continue
        duplicated = hash_counts.get((record_key, source_hash), 0) > 1
        if duplicated and source_index is None:
            continue
        key = _feature_identity_key(
            mode, record_key, f"{source_hash}~{source_index}" if duplicated else source_hash
        )
        record_ordinal = (
            offsets.get(input_index, 0) + record_index + 1 if linear else record_index + 1
        )
        once = (key, drawn_hash, record_ordinal)
        if not key or once in seen:
            continue
        seen.add(once)
        biological.append((key, drawn_hash, source_index, record_ordinal))
    return biological


def migrate_session_feature_edits(
    features: object, *, mode: object, catalog: object, legacy: object = None
) -> FeatureEditMigration:
    """Move the rendered-ID edit maps of a Session 44 or older to draft rows.

    Each edit becomes a ``features.featureOverrides`` row keyed by
    ``[mode, recordKey, biologicalFeatureId]`` through the Session's saved
    feature catalog (schema 3, 4, or 5), or, for a Session without one, through
    ``legacy`` (its request records and source features; see
    :func:`_legacy_feature_index`). Rule 1: a rendered ID of the catalog
    names every identity drawn with it. Rule 2: otherwise the hash, record
    position, and source index of its suffixes must name exactly one indexed
    feature: the source hash of a catalog feature, the drawn hash of a feature
    of a Session without a catalog. Any other edit is dropped and counted (a
    label source text is not an edit of its own). A Feature visibility edit
    that now names fewer features than its hash did is counted as narrowed.
    When a label map is migrated the saved label table (``labelOverrideRows``)
    is cleared, as it was built from those maps.

    This is the twin of ``migrateSessionFeatureEdits`` in the Web
    ``feature-edit-migration.js``;
    ``tests/fixtures/feature-edit-migration-vectors.json`` pins both. Without
    a catalog or ``legacy`` every edit is dropped.
    """

    source = features if isinstance(features, Mapping) else {}
    rows: dict[str, dict[str, Any]] = {}
    dropped_count = 0
    narrowed_visibility_count = 0

    def non_empty_map(field: str) -> bool:
        value = source.get(field)
        return isinstance(value, Mapping) and len(value) > 0

    migrated_label_edits = non_empty_map("labelVisibilityOverrides") or non_empty_map(
        "labelTextFeatureOverrides"
    )
    if any(non_empty_map(field) for field in RETIRED_RENDERED_ID_FEATURE_FIELDS):
        rendered_by_id, biological = (
            _catalog_feature_index(catalog, mode)
            if _js_truthy(catalog)
            else ({}, _legacy_feature_index(legacy, mode))
        )
        by_hash: dict[str, dict[str, None]] | None = None

        def resolve(old_key: str) -> list[str]:
            drawn = rendered_by_id.get(_rendered_id_without_part(old_key))
            if drawn:
                return list(drawn)
            stable_id, record_ordinal, source_index = _parse_rendered_id(old_key)
            if not stable_id:
                return []
            candidates = [
                key
                for key, feature_stable_id, feature_source_index, feature_ordinal in biological
                if feature_stable_id == stable_id
                and (record_ordinal is None or feature_ordinal == record_ordinal)
                and (source_index is None or feature_source_index == source_index)
            ]
            return candidates if len(candidates) == 1 else []

        def identities_with_hash(stable_id: str) -> dict[str, None]:
            # Every identity the `hash` row of a Session before 46 reached.
            nonlocal by_hash
            if by_hash is None:
                by_hash = {}
                for svg_id, keys in rendered_by_id.items():
                    by_hash.setdefault(_parse_rendered_id(svg_id)[0], {}).update(keys)
                for key, feature_stable_id, _, _ in biological:
                    by_hash.setdefault(feature_stable_id, {})[key] = None
            return by_hash.get(stable_id, {})

        def row_for(key: str) -> dict[str, Any]:
            if key not in rows:
                scope, record_key, feature_id = json.loads(key)
                rows[key] = {
                    "scope": scope,
                    "recordKey": record_key,
                    "biologicalFeatureId": feature_id,
                    **dict.fromkeys(_FEATURE_EDIT_DRAFT_FIELDS),
                }
            return rows[key]

        def assign_feature_visibility(row: dict[str, Any], value: object) -> bool:
            visibility = _FEATURE_VISIBILITY_EDIT_VALUES.get(_js_text(value).lower())
            if visibility is None:
                return False
            if row["featureVisibility"] is None:
                row["featureVisibility"] = visibility
            return True

        def assign_label_visibility(row: dict[str, Any], value: object) -> bool:
            visibility = _js_text(value).lower()
            if visibility not in ("on", "off"):
                return False
            if row["labelVisibility"] is None:
                row["labelVisibility"] = visibility
            return True

        def assign_label_text(row: dict[str, Any], value: object) -> bool:
            label_text = re.sub(r"[\t\r\n\0]+", " ", _js_string(value))
            if label_text.strip(_JS_WHITESPACE):
                if row["labelText"] is None:
                    row["labelText"] = label_text
            elif row["labelVisibility"] != "on":
                # A blank text hid the label (its table row drew an empty label).
                row["labelVisibility"] = "off"
            return True

        def assign_label_source_text(row: dict[str, Any], value: object) -> bool:
            source_text = _js_string(value)
            if not source_text or "\0" in source_text:
                return False
            if row["labelSourceText"] is None:
                row["labelSourceText"] = source_text
            return True

        for field, assign in (
            ("featureVisibilityOverrides", assign_feature_visibility),
            ("labelVisibilityOverrides", assign_label_visibility),
            ("labelTextFeatureOverrides", assign_label_text),
            ("labelTextFeatureOverrideSources", assign_label_source_text),
        ):
            edits = source.get(field)
            for old_key, value in (edits.items() if isinstance(edits, Mapping) else ()):
                keys = resolve(old_key)
                if not keys or not all(assign(row_for(key), value) for key in keys):
                    if field != "labelTextFeatureOverrideSources":
                        dropped_count += 1
                elif field == "featureVisibilityOverrides" and any(
                    key not in keys for key in identities_with_hash(_parse_rendered_id(old_key)[0])
                ):
                    narrowed_visibility_count += 1
        # A source text alone is kept only for a feature a bulk label edit can reach.
        rows = {
            key: row
            for key, row in rows.items()
            if any(row[field] is not None for field in _FEATURE_EDIT_DRAFT_FIELDS)
        }

    migrated = {
        key: value for key, value in source.items() if key not in RETIRED_RENDERED_ID_FEATURE_FIELDS
    }
    migrated["featureOverrides"] = rows
    if migrated_label_edits:
        migrated["labelOverrideRows"] = []
    return FeatureEditMigration(migrated, dropped_count, narrowed_visibility_count)


# A record selector's position is read as a JavaScript array index: the
# canonical decimal text of the value.
_JS_ARRAY_INDEX = re.compile(r"0|[1-9][0-9]*")
_EXPANDED_RECORD_ORDINAL = re.compile(r"[1-9][0-9]*")
_SOURCE_INDEX_SUFFIX = re.compile(r"~[0-9]+\Z")


def _js_truthy(value: object) -> bool:
    if value is None or isinstance(value, bool):
        return bool(value)
    if isinstance(value, (int, float)):
        return value != 0 and value == value
    if isinstance(value, str):
        return value != ""
    # Every JavaScript object and array is true, an empty one too.
    return True


def _drawn_transformed(record: Mapping[str, Any]) -> bool:
    """Whether a request record is drawn cropped, reverse-complemented, or rotated.

    The twin of ``drawnTransformed`` in the Web ``feature-edit-migration.js``.
    """

    presentation = record.get("presentation")
    display = record.get("display")
    return (
        _js_truthy(record.get("region"))
        or (isinstance(presentation, Mapping) and _js_truthy(presentation.get("reverseComplement")))
        or (isinstance(display, Mapping) and _safe_integer(display.get("startCoordinate")) is not None)
    )


def _record_key_belongs_to_record(record_key: str, record: Mapping[str, Any]) -> bool:
    """Whether a record key names a request record (an ALL record owns ``<key>:<n>``).

    The twin of ``recordKeyBelongsToRequest`` in the Web ``feature-placement.js``
    for a decoded request, whose record keys are strings.
    """

    request_key = record.get("recordKey")
    if not isinstance(request_key, str):
        return False
    return record_key == request_key or (
        record.get("cardinality") == "all"
        and record_key.startswith(f"{request_key}:")
        and _EXPANDED_RECORD_ORDINAL.fullmatch(record_key[len(request_key) + 1 :]) is not None
    )


@dataclass(frozen=True)
class AnnotationTargetMigration:
    """The annotation sets of an older Session and how many targets moved."""

    annotation_sets: Any
    migrated_count: int


def migrate_session_annotation_targets(
    annotation_sets: object, *, mode: object, catalog: object, records: object
) -> AnnotationTargetMigration:
    """Move the certain ``hash=`` annotation targets of a Session 44 or older.

    A Session before 46 named a selected feature in an annotation by
    ``hash=<hash>`` (a featureSpan target with one hash selector), which the
    renderer matches in the drawn record. Such a target becomes a
    featureIdentity target in ``mode`` only when the figure cannot change: the
    record it binds (the saved catalog's records in order, as the renderer
    binds them) is drawn without a crop, reverse complement, or rotation by its
    request record in ``records``, and the hash names exactly one feature of
    the saved ``catalog``, in that record. A moved target keeps the saved
    ``envelope`` and ``circularPath`` it has. Every other target stays as
    saved; when none moves, ``annotation_sets`` is returned as is.

    This is the twin of ``migrateSessionAnnotationTargets`` in the Web
    ``feature-edit-migration.js`` (R-7);
    ``tests/fixtures/annotation-target-migration-vectors.json`` pins both.
    Without a catalog nothing moves, as in the Web app, which reads no sources
    again for these targets.
    """

    record_ids: dict[str, set[str]] = {}
    features_by_hash: dict[str, dict[tuple[str, str], None]] = {}
    items = catalog.get("items") if isinstance(catalog, Mapping) else None
    for item in _mapping_list(items):
        item_record_keys = item.get("recordKeys")
        for record_key in item_record_keys if isinstance(item_record_keys, list) else ():
            record_ids.setdefault(_js_text(record_key), set())
        for feature in _mapping_list(item.get("biologicalFeatures")):
            record_key = _js_text(feature.get("recordKey"))
            feature_id = _js_text(feature.get("biologicalFeatureId"))
            if not record_key or not feature_id:
                continue
            # A record listed only by a later item has no ID from this feature.
            if record_key in record_ids:
                record_id = feature.get("record_id")
                record_ids[record_key].add(
                    _js_text(feature.get("recordId") if record_id is None else record_id)
                )
            source_hash = _js_text(feature.get("stableFeatureId")) or _SOURCE_INDEX_SUFFIX.sub(
                "", feature_id, count=1
            )
            features_by_hash.setdefault(source_hash, {})[(record_key, feature_id)] = None
    catalog_record_keys = list(record_ids)
    request_records = (
        [record for record in records if isinstance(record, Mapping)] if isinstance(records, list) else []
    )

    def bound_record_key(selector: object) -> str:
        if selector is None:
            return catalog_record_keys[0] if len(catalog_record_keys) == 1 else ""
        if not isinstance(selector, Mapping):
            return ""
        if selector.get("kind") == "recordIndex":
            index = _js_string(selector.get("index"))
            position = int(index) if _JS_ARRAY_INDEX.fullmatch(index) else len(catalog_record_keys)
            return catalog_record_keys[position] if position < len(catalog_record_keys) else ""
        if selector.get("kind") != "recordId":
            return ""
        # A record without catalog features has no known ID, so the binding is not certain.
        if any(len(record_ids[record_key]) != 1 for record_key in catalog_record_keys):
            return ""
        record_id = _js_text(selector.get("value"))
        matches = [record_key for record_key in catalog_record_keys if record_id in record_ids[record_key]]
        return matches[0] if len(matches) == 1 else ""

    def identity_target(target: object) -> dict[str, Any] | None:
        if not isinstance(target, Mapping) or target.get("kind") != "featureSpan":
            return None
        selectors = target.get("selectors")
        selector = selectors[0] if isinstance(selectors, list) and len(selectors) == 1 else None
        if not isinstance(selector, Mapping) or selector.get("key") != "hash":
            return None
        record_key = bound_record_key(target.get("record"))
        matches = list(features_by_hash.get(_js_text(selector.get("value")), {}))
        if not record_key or len(matches) != 1 or matches[0][0] != record_key:
            return None
        request = next(
            (record for record in request_records if _record_key_belongs_to_record(record_key, record)),
            None,
        )
        if request is None or _drawn_transformed(request):
            return None
        migrated: dict[str, Any] = {
            "kind": "featureIdentity",
            "scope": mode,
            "recordKey": matches[0][0],
            "biologicalFeatureId": matches[0][1],
        }
        for field in ("envelope", "circularPath"):
            if field in target:
                migrated[field] = target[field]
        return migrated if _feature_identity_key_of(migrated) else None

    migrated_count = 0
    migrated_sets: list[Any] = []
    for annotation_set in annotation_sets if isinstance(annotation_sets, list) else []:
        annotations = annotation_set.get("annotations") if isinstance(annotation_set, Mapping) else None
        if not isinstance(annotations, list):
            migrated_sets.append(annotation_set)
            continue
        migrated_annotations: list[Any] = []
        for annotation in annotations:
            target = identity_target(annotation.get("target") if isinstance(annotation, Mapping) else None)
            if target is None:
                migrated_annotations.append(annotation)
                continue
            migrated_count += 1
            migrated_annotations.append({**annotation, "target": target})
        migrated_sets.append({**annotation_set, "annotations": migrated_annotations})
    return AnnotationTargetMigration(
        migrated_sets if migrated_count else annotation_sets, migrated_count
    )


def migrate_persisted_web_state_field_names(config: object) -> object:
    """Project released Web config into the current shape without mutation."""

    if not isinstance(config, Mapping):
        return config

    migrated = dict(config)
    adv = config.get("adv")
    if isinstance(adv, Mapping):
        migrated_adv = dict(adv)
        for current, legacy, label in (
            (
                "linear_accession_visibility",
                "linear_show_accession",
                "Linear Accession visibility",
            ),
            (
                "linear_length_visibility",
                "linear_show_length",
                "Linear Length / Coordinates visibility",
            ),
        ):
            if current in migrated_adv:
                value = str(migrated_adv[current]).strip().lower()
                if value not in {"auto", "show", "hide"}:
                    raise ValidationError(
                        f"{label} must be one of: auto, show, hide."
                    )
                migrated_adv[current] = value
            elif legacy in migrated_adv:
                value = migrated_adv[legacy]
                if not isinstance(value, bool):
                    raise ValidationError(f"{label} legacy value must be a boolean.")
                migrated_adv[current] = "show" if value else "hide"
            else:
                migrated_adv[current] = "show"
            migrated_adv.pop(legacy, None)
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
    if "featurePlacementOverrides" in config:
        migrated["featurePlacementOverrides"] = _migrate_session_feature_placements(
            config["featurePlacementOverrides"]
        )
    drafts = config.get("recordDisplayDrafts")
    if isinstance(drafts, list):
        migrated["recordDisplayDrafts"] = [
            (
                row
                if not isinstance(row, Mapping)
                or "reverseComplementOverride" in row
                or "anchorIntent" in row
                else {
                    **row,
                    "reverseComplementOverride": None,
                    "anchorIntent": None,
                }
            )
            for row in drafts
        ]
    return migrated


# --- Session 46: the split of a Session 27-44 draft into mode slices ---------

# A Session before the pairwise match style existed drew ribbons; the flat
# draft of such a Session holds that value for the active mode (the twin of
# ``withHistoricalPairwiseMatchStyleFallback`` in the Web ``services/config.js``).
_HISTORICAL_FLAT_PROFILE_VALUES = {"pairwise_match_style": "ribbon"}
# A registry value that a Session 44 or older saved elsewhere.
_FLAT_DRAFT_SOURCES = {("ui", "selectedFeatureRecordIdx"): ("features", "selectedFeatureRecordIdx")}
_ANNOTATION_RECORD_BINDING_KEY = "_gbdraw_web_target_record_key"
_JSON_DECODER = json.JSONDecoder()
_ABSENT = object()


def _diagram_mode(value: object) -> DiagramMode | None:
    return cast("DiagramMode", value) if value in DIAGRAM_MODES else None


def _depth_slots(value: object) -> list[Any]:
    # depthFileSlotsFromValue: a list is the series slots; a value is one slot.
    if isinstance(value, list):
        return list(value)
    return [value] if _js_truthy(value) else []


def session_depth_source_widths(bindings: object) -> dict[DiagramMode, int]:
    """Each mode's Depth series count in the Web file bindings, 0 when none is bound.

    Circular reads ``c_depth`` (one row of series slots per record, or one
    row); Linear reads each ``linearSeqs[].depth``. The width is the widest row
    when any slot holds a file, as ``reconcileDepthTrackStateAfterSessionFiles``
    counts it.
    """

    source = bindings if isinstance(bindings, Mapping) else {}
    circular = source.get("c_depth")
    circular_rows = (
        [_depth_slots(row) for row in circular]
        if isinstance(circular, list) and all(isinstance(row, list) for row in circular)
        else [_depth_slots(circular)]
    )
    sequences = source.get("linearSeqs")
    linear_rows = [
        _depth_slots(sequence.get("depth"))
        for sequence in (sequences if isinstance(sequences, list) else [])
        if isinstance(sequence, Mapping)
    ]

    def width(rows: list[list[Any]]) -> int:
        if not any(_js_truthy(slot) for row in rows for slot in row):
            return 0
        return max(len(row) for row in rows)

    return {"circular": width(circular_rows), "linear": width(linear_rows)}


def _mode_profile_values(mode_profiles: object, mode: DiagramMode) -> Mapping[str, Any]:
    profiles = mode_profiles.get("profiles") if isinstance(mode_profiles, Mapping) else None
    profile = profiles.get(mode) if isinstance(profiles, Mapping) else None
    values = profile.get("values") if isinstance(profile, Mapping) else None
    return values if isinstance(values, Mapping) else {}


# The comparison colors every Web palette holds (``DEFAULT_COMPARISON_COLORS``
# in the Web ``utils/color-utils.js``).
_DEFAULT_COMPARISON_COLORS = {
    "pairwise_match": "#d3d3d3",
    "pairwise_match_min": "#FFE7E7",
    "pairwise_match_max": "#FF7272",
    "collinear_block_plus_min": "#f0f1f5",
    "collinear_block_plus": "#8b9cc1",
    "collinear_block_minus_min": "#FFE7E7",
    "collinear_block_minus": "#E15759",
}
_INDEXED_COLLINEAR_BLOCK_COLOR = re.compile(r"collinear_block_[0-9]+\Z")


def _normalize_palette_colors(colors: Mapping[str, Any]) -> dict[str, Any]:
    """The twin of ``normalizePaletteColors`` in the Web ``utils/color-utils.js``."""

    normalized = {key: value for key, value in colors.items() if not _INDEXED_COLLINEAR_BLOCK_COLOR.match(key)}
    if _js_truthy(normalized.get("collinear_block_plus_max")) and not _js_truthy(
        normalized.get("collinear_block_plus")
    ):
        normalized["collinear_block_plus"] = normalized["collinear_block_plus_max"]
    for key, value in _DEFAULT_COMPARISON_COLORS.items():
        if not _js_truthy(normalized.get(key)):
            normalized[key] = value
    return normalized


@functools.cache
def _color_palettes() -> dict[str, Any]:
    if sys.version_info >= (3, 11):
        import tomllib
    else:
        import tomli as tomllib
    from importlib import resources

    with resources.files("gbdraw.data").joinpath("color_palettes.toml").open("rb") as handle:
        return tomllib.load(handle)


def _palette_colors(name: object) -> dict[str, Any]:
    """A palette's normalized colors, or none (``paletteColorsFromDefinitions``).

    The Web reads the palettes from ``gallery/palettes/palettes.json``, which
    ``tools/generate_palette_explorer_assets.py`` writes from the same
    ``gbdraw/data/color_palettes.toml``.
    """

    palette = _js_text(name)
    colors = _color_palettes().get(palette) if palette != "title" else None
    if not isinstance(colors, Mapping) or not colors:
        return {}
    return _normalize_palette_colors({str(key): str(value) for key, value in colors.items()})


def _merges_override_colors(config: object) -> bool:
    colors = config.get("colors") if isinstance(config, Mapping) else None
    return (
        isinstance(config, Mapping)
        and _js_truthy(config.get("colorsAreOverrides"))
        and isinstance(colors, Mapping)
        and len(colors) > 0
    )


def mode_split_palette_colors(config: object) -> dict[str, Any] | None:
    """The palette colors a draft's override colors merge into, else ``None``.

    The split's ``palette_colors`` context: the draft palette's colors (an
    empty name reads as the ``default`` a fresh Load selects) when the draft
    keeps ``colorsAreOverrides`` with colors after the older migrations.
    """

    if not _merges_override_colors(config):
        return None
    assert isinstance(config, Mapping)
    return _palette_colors(_js_text(config.get("palette")) or "default")


def _with_resolved_override_colors(
    config: Mapping[str, Any], palette_colors: Mapping[str, Any] | None
) -> dict[str, Any]:
    """The draft ``config`` with ``colorsAreOverrides`` resolved and dropped.

    With the flag and colors, the stored colors override ``palette_colors``,
    as Web Load reads them (``services/config.js`` applyConfigData); otherwise
    the stored colors are complete. Named colors stay as saved: Web Load
    resolves them afterwards, as it does for complete colors.
    """

    resolved = {key: value for key, value in config.items() if key != "colorsAreOverrides"}
    if _merges_override_colors(config):
        # ``normalizeColorMap``: each value trimmed (``resolveColorToHex`` keeps
        # a name it cannot resolve without a browser).
        colors = config["colors"]
        overrides = {key: _js_text(value) if _js_truthy(value) else "" for key, value in colors.items()}
        resolved["colors"] = _normalize_palette_colors({**(palette_colors or {}), **overrides})
    return resolved


def _profile_active_mode(draft: Mapping[str, Any], mode_profiles: object) -> DiagramMode | None:
    """The mode whose values a flat draft holds for the mode-profile fields.

    ``ui.mode``, else ``renderRequest.mode``, else ``modeProfiles.activeMode``:
    the twin of the rule of ``withHistoricalPairwiseMatchStyleFallback`` in the
    Web ``services/config.js``.
    """

    ui = draft.get("ui")
    request = draft.get("renderRequest")
    return (
        _diagram_mode(ui.get("mode") if isinstance(ui, Mapping) else None)
        or _diagram_mode(request.get("mode") if isinstance(request, Mapping) else None)
        or _diagram_mode(mode_profiles.get("activeMode") if isinstance(mode_profiles, Mapping) else None)
    )


def _unscoped_rows(rows: object, mode: DiagramMode) -> object:
    """The ``scope``d draft rows of ``mode``, keyed by identity pair and without ``scope``."""

    if not isinstance(rows, Mapping):
        return rows
    unscoped: dict[str, Any] = {}
    for key, row in rows.items():
        if isinstance(row, Mapping) and "scope" in row:
            if row["scope"] != mode:
                continue
            fields = {field: value for field, value in row.items() if field != "scope"}
            record_key, feature_id = fields.get("recordKey"), fields.get("biologicalFeatureId")
            if isinstance(record_key, str) and isinstance(feature_id, str):
                key = _draft_pair_key(record_key, feature_id)
            unscoped[key] = fields
        else:
            unscoped[key] = row
    return unscoped


def _annotation_binding_mode(annotation: object) -> DiagramMode | None:
    """The mode whose record an annotation's target names, if it names one."""

    if not isinstance(annotation, Mapping):
        return None
    target = annotation.get("target")
    if isinstance(target, Mapping) and target.get("kind") == "featureIdentity":
        return _diagram_mode(target.get("scope"))
    metadata = annotation.get("metadata")
    binding = metadata.get(_ANNOTATION_RECORD_BINDING_KEY) if isinstance(metadata, Mapping) else None
    if not isinstance(binding, str):
        return None
    # A record key starts with its source key, a JSON array whose first item
    # is the mode (``annotationSourceKey``).
    try:
        source_key, _ = _JSON_DECODER.raw_decode(binding.strip(_JS_WHITESPACE))
    except ValueError:
        return None
    return _diagram_mode(source_key[0]) if isinstance(source_key, list) and source_key else None


def _annotation_sets_of_mode(sets: list[Any], mode: DiagramMode) -> list[Any]:
    result: list[Any] = []
    for annotation_set in sets:
        annotations = annotation_set.get("annotations") if isinstance(annotation_set, Mapping) else None
        if not isinstance(annotations, list):
            result.append(_json_clone(annotation_set))
            continue
        kept = []
        for annotation in annotations:
            bound = _annotation_binding_mode(annotation)
            if bound is not None and bound != mode:
                continue
            annotation = _json_clone(annotation)
            target = annotation.get("target") if isinstance(annotation, dict) else None
            if isinstance(target, dict) and target.get("kind") == "featureIdentity":
                target.pop("scope", None)
            kept.append(annotation)
        result.append({**_json_clone(annotation_set), "annotations": kept})
    return result


def _split_setting(
    row: ModeScopedSetting,
    value: Any,
    *,
    committed: DiagramMode,
    widths: Mapping[DiagramMode, int],
) -> dict[DiagramMode, Any]:
    """The value each slice takes for one saved registry value."""

    if row.migrate == "own":
        return {cast("DiagramMode", row.modes): _json_clone(value)}
    if row.migrate == "result-mode":
        return {committed: _unscoped_rows(_json_clone(value), committed)}
    if row.migrate == "layout":
        slots = value if isinstance(value, Mapping) else {}
        return {mode: _json_clone(slots[mode]) for mode in DIAGRAM_MODES if mode in slots}
    split: dict[DiagramMode, Any] = {}
    for mode in DIAGRAM_MODES:
        if row.migrate == "show-if-source":
            # ``show_depth && hasSource``
            split[mode] = (widths[mode] > 0) if _js_truthy(value) else value
        elif row.migrate == "depth" and isinstance(value, list):
            split[mode] = _json_clone(value[: max(1, widths[mode])])
        elif row.migrate == "by-scope" and isinstance(value, list):
            split[mode] = [
                {field: item for field, item in _json_clone(draft).items() if field != "scope"}
                for draft in value
                if isinstance(draft, Mapping) and draft.get("scope") == mode
            ]
        elif row.migrate == "by-side":
            split[mode] = _unscoped_rows(_json_clone(value), mode)
        elif row.migrate == "by-leaf" and isinstance(value, Mapping):
            split[mode] = {
                path: _json_clone(leaf)
                for path, leaf in value.items()
                if mode in unmanaged_config_override_modes(str(path))
            }
        elif row.migrate == "by-binding" and isinstance(value, list):
            split[mode] = _annotation_sets_of_mode(value, mode)
        else:
            split[mode] = _json_clone(value)
    return split


_SHORT_HEX_COLOR = re.compile(r"#([0-9a-fA-F]{3})")
_LONG_HEX_COLOR = re.compile(r"#[0-9a-fA-F]{6}")


def _optional_hex_color(value: object) -> str | None:
    """``normalizeOptionalHexColor`` of the Web ``services/config.js``.

    A 3- or 6-digit hex code, or a color name through the shared named-color
    table (``gbdraw.io.colors``), becomes lowercase 6-digit hex; anything else
    is ``None`` (OV-160).
    """

    if value is None or value == "":
        return None
    text = _js_text(value)
    if text and not text.startswith("#"):
        text = named_color_hex(text) or text
    short = _SHORT_HEX_COLOR.fullmatch(text)
    if short:
        return "#" + "".join(char * 2 for char in short.group(1)).lower()
    return text.lower() if _LONG_HEX_COLOR.fullmatch(text) else None


def _stroke_width(value: object) -> int | float | None:
    """``normalizeStrokeWidth`` of the Web ``services/config.js``."""

    if value is None or value == "":
        return None
    number = _js_number(value)
    return _js_integral(number) if math.isfinite(number) and number >= 0 else None


def _stroke_override_map(source: object) -> dict[str, Any]:
    """``normalizeStrokeOverrideMap(source, { requireOverride: true })`` of the Web."""

    normalized: dict[str, Any] = {}
    if not isinstance(source, Mapping):
        return normalized
    for key, value in source.items():
        name = _js_text(key)
        if not name or not isinstance(value, Mapping):
            continue
        override: dict[str, Any] = {}
        color = _optional_hex_color(value.get("strokeColor"))
        width = _stroke_width(value.get("strokeWidth"))
        if color is not None:
            override["strokeColor"] = color
        if width is not None:
            override["strokeWidth"] = width
        if not override:
            continue
        if "originalStrokeColor" in value:
            override["originalStrokeColor"] = _optional_hex_color(value["originalStrokeColor"])
        if "originalStrokeWidth" in value:
            override["originalStrokeWidth"] = _stroke_width(value["originalStrokeWidth"])
        normalized[name] = override
    return normalized


def _legend_color_overrides(source: object) -> dict[str, str]:
    """``normalizeLegendColorOverrides`` of the Web ``services/config.js``."""

    normalized: dict[str, str] = {}
    if not isinstance(source, Mapping):
        return normalized
    for key, value in source.items():
        caption = _js_text(key)
        color = _optional_hex_color(value)
        if caption and color:
            normalized[caption] = color
    return normalized


def _without_stroke_disclosure(entries: object) -> object:
    """Legend rows without ``showStroke``: the Stroke options disclosure that
    earlier Sessions saved is view state (OV-157), which Web Load drops
    (``normalizeSessionLegendEntries``); other values stay as they are."""

    if not isinstance(entries, list):
        return entries
    return [
        {key: value for key, value in entry.items() if key != "showStroke"}
        if isinstance(entry, Mapping) else entry
        for entry in entries
    ]


def _with_normalized_editor_colors(draft: Mapping[str, Any]) -> Mapping[str, Any]:
    """The draft with the Legend and feature stroke and color edits and the
    original SVG stroke as Web Load reads them (``normalizeEditorStateData``),
    and its Legend rows without the Stroke options disclosure; absent fields
    stay absent."""

    editor = draft.get("editorState")
    if not isinstance(editor, Mapping):
        return draft
    changed = dict(editor)
    legend = editor.get("legend")
    if isinstance(legend, Mapping):
        changed_legend = dict(legend)
        if "strokeOverrides" in legend:
            changed_legend["strokeOverrides"] = _stroke_override_map(legend["strokeOverrides"])
        if "colorOverrides" in legend:
            changed_legend["colorOverrides"] = _legend_color_overrides(legend["colorOverrides"])
        for rows in ("entries", "deletedEntries", "dormantEntries"):
            if rows in legend:
                changed_legend[rows] = _without_stroke_disclosure(legend[rows])
        changed["legend"] = changed_legend
    strokes = editor.get("featureStrokes")
    if isinstance(strokes, Mapping) and "overrides" in strokes:
        changed["featureStrokes"] = {**strokes, "overrides": _stroke_override_map(strokes["overrides"])}
    svg_stroke = editor.get("originalSvgStroke")
    if isinstance(svg_stroke, Mapping):
        changed_stroke = dict(svg_stroke)
        if "color" in svg_stroke:
            changed_stroke["color"] = _optional_hex_color(svg_stroke["color"])
        if "width" in svg_stroke:
            changed_stroke["width"] = _stroke_width(svg_stroke["width"])
        changed["originalSvgStroke"] = changed_stroke
    return {**draft, "editorState": changed}


def split_draft_into_modes(
    draft: Mapping[str, Any],
    *,
    committed_mode: object,
    mode_profiles: object = None,
    depth_sources: Mapping[DiagramMode, int] | Mapping[str, object] | None = None,
    palette_colors: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Split the flat Web draft of a Session 27-44 into Session 46 mode slices.

    ``draft`` holds the Session's ``config``, ``features``, ``editorState``,
    and ``ui`` after the older normalizers (field names and placement rows,
    per-feature edits, annotation targets); other fields are kept as they are.
    The result has ``modes`` instead of the flat draft when the draft has a
    ``config`` (a Session written without one, by the CLI or the Python API,
    gets none, and its slice values are dropped): each registry row
    (MODE_SCOPED_SETTINGS) fills the slices by its ``migrate`` token, a value
    that a slice does not take is absent there (that mode's default), and a
    field that no row names is dropped. App-level settings move to the top
    level: ``config.losat``'s execution settings to ``ui.losatExecution``,
    ``config.adv.rich_feature_popup`` to ``ui.richFeaturePopup``, a boolean
    ``config.paletteInstantPreviewEnabled`` to ``ui`` (over a saved ``ui``
    value, as Load applies it last), and ``config.cliOptions`` to
    ``cliOptions``. ``config.colorsAreOverrides`` is resolved into ``colors``
    over ``palette_colors`` (``mode_split_palette_colors``) and dropped. The
    Legend and feature stroke and color edits are read as Web Load reads them
    (hex colors, color names through the shared table, else ``None``), and the
    Legend rows lose ``showStroke`` (OV-157).

    ``committed_mode`` is the saved Result's mode (``renderRequest.mode``, else
    ``ui.mode``). ``mode_profiles`` defaults to ``config.modeProfiles``, and
    ``depth_sources`` (each mode's Depth series count) to the counts in
    ``webFiles.bindings``, else in the legacy ``files``. The draft's own mode
    (``_profile_active_mode``) holds the flat values of the mode-profile fields.

    This is the twin of ``splitDraftIntoModes`` in the Web
    ``services/mode-scoped-migration.js``;
    ``tests/fixtures/sessions/mode-split-vectors.json`` pins both.
    """

    committed = _diagram_mode(committed_mode)
    if committed is None:
        raise ValidationError(
            "Splitting a Session draft by mode requires the saved Result's mode.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )
    draft = _with_normalized_editor_colors(draft)
    config_value = draft.get("config")
    if isinstance(config_value, Mapping) and "colorsAreOverrides" in config_value:
        config_value = _with_resolved_override_colors(config_value, palette_colors)
        draft = {**draft, "config": config_value}
    config: Mapping[str, Any] = config_value if isinstance(config_value, Mapping) else {}
    profiles = config.get("modeProfiles") if mode_profiles is None else mode_profiles
    active = _profile_active_mode(draft, profiles) or committed
    if depth_sources is None:
        web_files = draft.get("webFiles")
        bindings = web_files.get("bindings") if isinstance(web_files, Mapping) else draft.get("files")
        widths = session_depth_source_widths(bindings)
    else:
        widths = {mode: _safe_integer(depth_sources.get(mode)) or 0 for mode in DIAGRAM_MODES}
    adv_value = config.get("adv")
    adv: Mapping[str, Any] = adv_value if isinstance(adv_value, Mapping) else {}

    slices: dict[DiagramMode, dict[str, Any]] = {mode: {} for mode in DIAGRAM_MODES}
    for row in MODE_SCOPED_SETTINGS:
        source_domain, source_path = _FLAT_DRAFT_SOURCES.get((row.domain, row.path), (row.domain, row.path))
        container = _container_at(draft, source_domain)
        if not isinstance(container, Mapping):
            continue
        split: dict[DiagramMode, Any]
        if row.migrate == "profile":
            split = {}
            flat = container.get(source_path, _HISTORICAL_FLAT_PROFILE_VALUES.get(source_path, _ABSENT))
            if flat is not _ABSENT:
                split[active] = _json_clone(flat)
            other: DiagramMode = "linear" if active == "circular" else "circular"
            saved = _mode_profile_values(profiles, other)
            if source_path in saved:
                split[other] = _json_clone(saved[source_path])
        elif source_path in container:
            split = _split_setting(row, container[source_path], committed=committed, widths=widths)
        else:
            continue
        for mode, value in split.items():
            target = slices[mode]
            for part in row.domain.split("."):
                target = target.setdefault(part, {})
            target[row.path] = value

    result = {key: value for key, value in draft.items() if key not in FLAT_DRAFT_TOP_LEVEL_FIELDS}
    for domain, fields in _RETIRED_MODE_SCOPED_FIELDS.items():
        if not domain:
            continue
        head, _, rest = domain.partition(".")
        container = _container_at(result, domain)
        if not isinstance(container, Mapping) or not fields & set(container):
            continue
        kept = {key: value for key, value in container.items() if key not in fields}
        if rest:
            parent = dict(result[head])
            if kept:
                parent[rest] = kept
            else:
                parent.pop(rest, None)
            result[head] = parent
        else:
            result[head] = kept
    losat = config.get("losat")
    execution = {
        field: _json_clone(losat[field])
        for field in LOSAT_EXECUTION_FIELDS
        if isinstance(losat, Mapping) and field in losat
    }
    # App-level settings leave the draft: LOSAT execution, the rich feature
    # popup (a missing value reads true), and Instant Preview.
    app_ui: dict[str, Any] = {"losatExecution": execution} if execution else {}
    if "rich_feature_popup" in adv:
        app_ui["richFeaturePopup"] = _json_clone(adv["rich_feature_popup"])
    # Web Load applies ``ui``, then a boolean draft value, which wins.
    preview = config.get("paletteInstantPreviewEnabled")
    if isinstance(preview, bool):
        app_ui["paletteInstantPreviewEnabled"] = preview
    if app_ui:
        result["ui"] = {**(result.get("ui") or {}), **app_ui}
    # Session provenance, written only when the draft has it.
    if "cliOptions" in config:
        result["cliOptions"] = _json_clone(config["cliOptions"])
    if isinstance(config_value, Mapping):
        result["modes"] = {mode: slices[mode] for mode in DIAGRAM_MODES}
    return result


# --- Sessions 27-44: the Web Load value migrations of the draft config ------
#
# Twins of the config steps of ``migrateSessionDataToCurrent`` (Sessions before
# 40) and of the current-writer restore (40-44) in the Web ``services/config.js``.
# ``tools/generate_draft_value_migration_vectors.mjs`` runs the Web functions and
# writes the vectors these twins read.

_CURRENT_CIRCULAR_TRACK_SLOT_SCHEMA = 4
_LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA = 3
_CURRENT_LINEAR_TRACK_SLOT_SCHEMA = 2
_LEGACY_LINEAR_TRACK_SLOT_SCHEMA = 1
# A Session up to this version saved Linear slots with schema-1 meaning,
# whatever schema version it stored (LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION).
_LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32
# The Web writers of Sessions 27-33 saved slot and feature-shape forms that
# today's readers fill in (``withoutLegacyNullCircularSlotSpacing`` and
# ``migrateLegacyFeatureRenderingConfig``).
_LEGACY_SLOT_SHAPE_SESSION_VERSION = 33
_LINEAR_TRACK_RENDERERS = (
    "features", "dinucleotide_content", "dinucleotide_skew", "depth", "annotations", "spacer",
)
_LINEAR_TRACK_RENDERER_ALIASES = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
}


def _js_integer(value: object) -> bool:
    return (isinstance(value, int) and not isinstance(value, bool)) or (
        isinstance(value, float) and value.is_integer()
    )


def _current_option_value(value: object, fallback: str, supported: tuple[str, ...], label: str) -> str:
    """``requireCurrentValue`` in the Web ``services/current-option-values.js``."""

    normalized = ("" if value is None else _js_string(value)).strip(_JS_WHITESPACE).lower() or fallback
    if normalized not in supported:
        raise ValidationError(
            f"{label} must be one of: {', '.join(supported)}.", diagnostic=_SESSION_FIELDS_INVALID
        )
    return normalized


def _migrated_linear_track_layout(value: object) -> str:
    normalized = ("" if value is None else _js_string(value)).strip(_JS_WHITESPACE).lower() or "middle"
    normalized = {"spreadout": "above", "tuckin": "below"}.get(normalized, normalized)
    return _current_option_value(normalized, "middle", ("above", "middle", "below"), "Linear track layout")


def _migrated_linear_label_placement(value: object) -> str:
    normalized = ("" if value is None else _js_string(value)).strip(_JS_WHITESPACE).lower() or "auto"
    return _current_option_value(
        "above_feature" if normalized == "on_feature" else normalized,
        "auto",
        ("auto", "above_feature"),
        "Linear label placement",
    )


def _migrated_circular_multi_record_size_mode(value: object) -> str:
    normalized = ("" if value is None else _js_string(value)).strip(_JS_WHITESPACE).lower() or "auto"
    return _current_option_value(
        "auto" if normalized == "sqrt" else normalized,
        "auto",
        ("auto", "linear", "equal"),
        "Circular multi-record size mode",
    )


def migrate_persisted_web_option_values(config: object) -> object:
    """Move a Session 27-39 draft's retired option values to today's.

    Linear track layout ``spreadout``/``tuckin`` becomes ``above``/``below``,
    label placement ``on_feature`` becomes ``above_feature``, and multi-record
    size ``sqrt`` becomes ``auto``; any other value must be a current one. The
    twin of the option-value steps of ``migratePersistedWebOptionValues``.
    """

    if not isinstance(config, Mapping):
        return config
    migrated = dict(config)
    form = config.get("form")
    if isinstance(form, Mapping) and "linear_track_layout" in form:
        migrated["form"] = {**form, "linear_track_layout": _migrated_linear_track_layout(form["linear_track_layout"])}
    adv = config.get("adv")
    if isinstance(adv, Mapping):
        adv = dict(adv)
        if "label_placement" in adv:
            adv["label_placement"] = _migrated_linear_label_placement(adv["label_placement"])
        if "multi_record_size_mode" in adv:
            adv["multi_record_size_mode"] = _migrated_circular_multi_record_size_mode(adv["multi_record_size_mode"])
        migrated["adv"] = adv
    return migrated


def _require_current_circular_track_slots(config: Mapping[str, Any]) -> None:
    """Refuse Circular slots that the Web would still have to migrate.

    The Web migrates schema-3 Circular slots (``migrateImportedCircularTrackSlots``),
    but no Session 27-44 written by ``main`` or a release has them, so Python
    has no twin of that step and refuses such a draft instead of writing it.
    """

    from .session import SessionFormatError

    adv = config.get("adv")
    if not isinstance(adv, Mapping) or "circular_track_slots" not in adv:
        return
    if not isinstance(adv["circular_track_slots"], list):
        raise ValidationError("Custom Track Slots must be an array.", diagnostic=_SESSION_FIELDS_INVALID)
    schema = adv.get("circular_track_slots_schema_version")
    if _js_integer(schema) and schema == _CURRENT_CIRCULAR_TRACK_SLOT_SCHEMA:
        return
    if _js_integer(schema) and schema == _LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA:
        raise SessionFormatError(
            "Custom Track Slots of Circular slot schema 3 cannot be migrated by Python; "
            "open and save the Session in the Web app first.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )
    raise ValidationError(
        "Custom Track Slots use an obsolete schema. Recreate the slots with schema version "
        f"{_CURRENT_CIRCULAR_TRACK_SLOT_SCHEMA}.",
        diagnostic=_SESSION_FIELDS_INVALID,
    )


def _without_legacy_null_circular_slot_spacing(config: Mapping[str, Any]) -> Mapping[str, Any]:
    """Drop the ``spacing: null`` that Session 27-33 writers saved in each Circular slot.

    Only while Custom Track Slots are off, where the null is lossless. The twin
    of ``withoutLegacyNullCircularSlotSpacing``.
    """

    adv = config.get("adv")
    if (
        not isinstance(adv, Mapping)
        or _js_truthy(adv.get("circular_track_slots_enabled"))
        or not (
            _js_integer(adv.get("circular_track_slots_schema_version"))
            and adv["circular_track_slots_schema_version"] == _CURRENT_CIRCULAR_TRACK_SLOT_SCHEMA
        )
        or not isinstance(adv.get("circular_track_slots"), list)
    ):
        return config
    return {
        **config,
        "adv": {
            **adv,
            "circular_track_slots": [
                {key: value for key, value in slot.items() if key != "spacing"}
                if isinstance(slot, Mapping) and "spacing" in slot and slot["spacing"] is None
                else slot
                for slot in adv["circular_track_slots"]
            ],
        },
    }


def _linear_track_renderer(value: object) -> str:
    text = (_js_string(value) if _js_truthy(value) else "features").strip(_JS_WHITESPACE).lower()
    renderer = _LINEAR_TRACK_RENDERER_ALIASES.get(text) or text
    return renderer if renderer in _LINEAR_TRACK_RENDERERS else "features"


def migrate_imported_linear_track_slots(config: object, source_version: object) -> object:
    """Bring a draft's Linear slots to schema 2.

    A Session up to 32 saved them with schema-1 meaning whatever it stored;
    later ones state their schema. A schema-1 features slot loses its height and
    spacing, and every slot gets its own ``params`` object. The twin of
    ``migrateImportedLinearTrackSlots`` and ``migrateLinearTrackSlotsToCurrentSchema``.
    """

    adv = config.get("adv") if isinstance(config, Mapping) else None
    if not isinstance(config, Mapping) or not isinstance(adv, Mapping) or "linear_track_slots" not in adv:
        return config
    slots = adv["linear_track_slots"]
    if not isinstance(slots, list):
        raise ValidationError("Custom Track Slots must be an array.", diagnostic=_SESSION_FIELDS_INVALID)
    obsolete = ValidationError(
        "Custom Track Slots use an obsolete schema. Recreate the slots with schema version "
        f"{_CURRENT_LINEAR_TRACK_SLOT_SCHEMA}.",
        diagnostic=_SESSION_FIELDS_INVALID,
    )
    supported = (_LEGACY_LINEAR_TRACK_SLOT_SCHEMA, _CURRENT_LINEAR_TRACK_SLOT_SCHEMA)
    stored = adv.get("linear_track_slots_schema_version", _ABSENT)
    if stored is not _ABSENT and not (_js_integer(stored) and stored in supported):
        raise obsolete
    if _js_integer(source_version) and cast(int, source_version) <= _LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION:
        schema: object = _LEGACY_LINEAR_TRACK_SLOT_SCHEMA
    else:
        schema = stored
    if schema is _ABSENT or schema not in supported:
        raise obsolete

    def migrated(slot: object) -> object:
        if not isinstance(slot, Mapping):
            return slot
        params = slot.get("params")
        result = {**slot, "params": dict(params) if isinstance(params, Mapping) else {}}
        if schema == _LEGACY_LINEAR_TRACK_SLOT_SCHEMA and _linear_track_renderer(slot.get("renderer")) == "features":
            result.pop("height", None)
            result.pop("spacing", None)
        return result

    return {
        **config,
        "adv": {
            **adv,
            "linear_track_slots_schema_version": _CURRENT_LINEAR_TRACK_SLOT_SCHEMA,
            "linear_track_slots": [migrated(slot) for slot in slots],
        },
    }


def _with_legacy_repeat_region_shape(config: Mapping[str, Any]) -> Mapping[str, Any]:
    """Give a Session 27-33 draft's repeat regions the rectangle they were drawn as.

    The twin of ``migrateLegacyFeatureRenderingConfig`` for a legacy Session.
    """

    adv = config.get("adv")
    if not isinstance(adv, Mapping):
        return config
    features = adv.get("features")
    if isinstance(features, list) and "repeat_region" not in features:
        return config
    shapes = adv.get("feature_shapes")
    shapes = shapes if isinstance(shapes, Mapping) else {}
    if "repeat_region" in shapes:
        return config
    return {**config, "adv": {**adv, "feature_shapes": {**shapes, "repeat_region": "rectangle"}}}


def migrate_session_draft_values(config: object, source_version: int) -> object:
    """The Web Load value migrations of a Session 27-44 draft config, in its order.

    Before 40 (``migrateSessionDataToCurrent``): the option values, the
    Circular slot check, the Session 27-33 slot spacing, the Linear slots, and
    the Session 27-33 repeat-region shape. From 40 (the current-writer
    restore): the Circular slot check and the Linear slots, and a stored
    ``colors`` drops ``colorsAreOverrides`` (``restoreCurrentWriterActiveConfig``),
    so such colors load as saved.
    """

    if not isinstance(config, Mapping):
        return config
    if source_version < CURRENT_AUTHORITY_SESSION_MIN_VERSION:
        migrated = migrate_persisted_web_option_values(config)
        assert isinstance(migrated, Mapping)
        config = migrated
    elif "colors" in config:
        config = {key: value for key, value in config.items() if key != "colorsAreOverrides"}
    _require_current_circular_track_slots(config)
    legacy_shapes = source_version <= _LEGACY_SLOT_SHAPE_SESSION_VERSION
    if legacy_shapes:
        config = _without_legacy_null_circular_slot_spacing(config)
    migrated = migrate_imported_linear_track_slots(config, source_version)
    assert isinstance(migrated, Mapping)
    return _with_legacy_repeat_region_shape(migrated) if legacy_shapes else dict(migrated)


_DEFAULT_ANNOTATION_STYLE: dict[str, Any] = {
    "stroke": "#404040",
    "strokeWidth": 1.5,
    "strokeDasharray": [],
    "lineCap": "tick",
    "fill": "#94a3b8",
    "fillOpacity": 0.2,
    "hatch": None,
    "labelColor": "#202020",
    "labelFontSize": None,
    "labelOrientation": "auto",
    "labelPosition": "center",
    "labelOffset": 4,
}
_ANNOTATION_MARKS = ("line", "bracket", "band", "highlight")


def _js_number(value: object) -> float:
    """``Number(value)`` for JSON values (NaN for what JavaScript cannot read)."""

    if value is None:
        return 0.0
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip(_JS_WHITESPACE)
        if not text:
            return 0.0
        try:
            return float(text) if not re.search(r"[^0-9eE.+-]", text) else math.nan
        except ValueError:
            return math.nan
    if isinstance(value, list) and len(value) <= 1:
        return _js_number(value[0]) if value else 0.0
    return math.nan


def _js_integral(value: float) -> int | float:
    return int(value) if value.is_integer() else value


def _clean_id(value: object, fallback: str) -> str:
    return (_js_string(value) if _js_truthy(value) else "").strip(_JS_WHITESPACE) or fallback


def _draft_annotation_target(target: object) -> dict[str, Any]:
    source = _json_clone(target) if isinstance(target, Mapping) else {}
    envelope = "segments" if source.get("envelope") == "segments" else "outer_bounds"
    circular_path = source.get("circularPath") if source.get("circularPath") in ("forward", "reverse") else "shortest"
    if source.get("kind") == "featureIdentity":
        if not _feature_identity_key_of(source):
            raise ValidationError("Invalid selected-feature annotation target.", diagnostic=_SESSION_FIELDS_INVALID)
        return {
            "kind": "featureIdentity",
            "scope": source.get("scope"),
            "recordKey": source.get("recordKey"),
            "biologicalFeatureId": source.get("biologicalFeatureId"),
            "envelope": envelope,
            "circularPath": circular_path,
        }
    if source.get("kind") == "featureSpan":
        selectors = source.get("selectors")
        return {
            "kind": "featureSpan",
            "record": source.get("record"),
            "selectors": [
                {
                    "key": None
                    if not isinstance(selector, Mapping) or selector.get("key") in (None, "")
                    else _js_string(selector.get("key")),
                    "value": _js_string(selector.get("value"))
                    if isinstance(selector, Mapping) and _js_truthy(selector.get("value"))
                    else "",
                }
                for selector in selectors
            ]
            if isinstance(selectors, list)
            else [],
            "envelope": envelope,
            "circularPath": circular_path,
        }
    start = max(1.0, _js_number(source.get("start")) or 1.0)
    end = max(1.0, _js_number(source.get("end")) or 1.0)
    return {
        "kind": "coordinateSpan",
        "record": source.get("record"),
        "start": _js_integral(start),
        "end": _js_integral(end),
        "coordinateSpace": "local" if source.get("coordinateSpace") == "local" else "source",
        "wrapsOrigin": start > end,
        "outOfBounds": source.get("outOfBounds") if source.get("outOfBounds") in ("skip", "error") else "clip",
    }


def _draft_annotation_style(style: object) -> dict[str, Any]:
    return {**_DEFAULT_ANNOTATION_STYLE, **(_json_clone(style) if isinstance(style, Mapping) else {})}


def _draft_annotation_sets_of_request(sets: object, mode: DiagramMode) -> list[dict[str, Any]]:
    """The Web draft annotation sets of a request's sets, targets in ``mode``.

    The twin of ``draftAnnotationSetsOfRequest`` (``normalizeAnnotationSets``)
    in the Web ``services/annotation-state.js``.
    """

    used_set_ids: set[str] = set()
    result: list[dict[str, Any]] = []
    for set_index, raw_set in enumerate(sets if isinstance(sets, list) else []):
        source = raw_set if isinstance(raw_set, Mapping) else {}
        set_id = _clean_id(_clean_id(source.get("id"), "annotations"), f"annotations_{set_index + 1}")
        while set_id in used_set_ids:
            set_id = f"{set_id}_{set_index + 1}"
        used_set_ids.add(set_id)
        legend_label = source.get("legendLabel")
        used_item_ids: set[str] = set()
        annotations: list[dict[str, Any]] = []
        raw_items = source.get("annotations")
        for item_index, raw_item in enumerate(raw_items if isinstance(raw_items, list) else []):
            item = _json_clone(raw_item) if isinstance(raw_item, Mapping) else {}
            target = item.get("target")
            if isinstance(target, Mapping) and target.get("kind") == "featureIdentity":
                target = {"scope": mode, **target}
            item_id = _clean_id(item.get("id"), f"region_{item_index + 1}")
            while item_id in used_item_ids:
                item_id = f"{item_id}_{item_index + 1}"
            used_item_ids.add(item_id)
            lane = item.get("lane")
            metadata = item.get("metadata")
            annotations.append(
                {
                    "id": item_id,
                    "target": _draft_annotation_target(target),
                    "label": _js_string(item.get("label")) if _js_truthy(item.get("label")) else "",
                    "mark": item.get("mark") if item.get("mark") in _ANNOTATION_MARKS else "bracket",
                    "lane": None
                    if lane is None or lane == ""
                    else _js_integral(max(0.0, _js_number(lane) or 0.0)),
                    "style": None if item.get("style") is None else _draft_annotation_style(item.get("style")),
                    "legendLabel": None if item.get("legendLabel") is None else _js_string(item.get("legendLabel")),
                    "metadata": _json_clone(metadata) if isinstance(metadata, Mapping) else {},
                }
            )
        result.append(
            {
                "id": set_id,
                "annotations": annotations,
                "defaultStyle": _draft_annotation_style(source.get("defaultStyle")),
                "legendLabel": None if legend_label is None else _js_string(legend_label),
            }
        )
    return result


# The ``config`` keys the CLI writer of Sessions 40 and 41 derived from its
# options (``_cli_web_config`` on main 8228ffab..4e8c9380; later ones wrote
# ``{"adv": {}}``). The Web writer of those Sessions always added
# ``annotationSets``, ``linearComparisonPlan`` and ``webEdits``, so its draft
# never holds only these keys. Mirrored by ``CLI_WRITER_CONFIG_DOMAINS`` in
# gbdraw/web/js/services/session-active-config-contract.js.
_CLI_WRITER_CONFIG_DOMAINS = frozenset(
    {
        "form", "adv", "losat", "cliOptions", "colors", "palette", "rules",
        "qualifierPriorityRules", "filterMode", "whitelist", "blacklistText",
        "losatProgram", "circularConservation",
    }
)


def _holds_cli_writer_config(session: Mapping[str, Any]) -> bool:
    """A Session 40 or 41 the CLI wrote: its ``config`` is no Web draft (OV-269)."""

    config = session.get("config")
    invocation = session.get("cliInvocation")
    return (
        session.get("version") in (40, 41)
        and isinstance(invocation, Mapping)
        and invocation.get("generatedBy") == "gbdraw"
        and isinstance(config, Mapping)
        and set(config) <= _CLI_WRITER_CONFIG_DOMAINS
    )


@dataclass(frozen=True)
class SessionDraftMigration:
    """A Session 27-44 Web draft after the migrations that Web Load runs.

    ``session`` still has the flat shape (``config``, ``features``); the split
    into mode slices (``split_draft_into_modes``) follows. A config-less
    Session 40-44 whose request names features by ``hash=`` gets the request's
    annotation sets as its draft (OV-135).
    """

    session: dict[str, Any]
    dropped_feature_edit_count: int = 0
    narrowed_visibility_count: int = 0
    migrated_annotation_count: int = 0


def migrate_session_flat_draft(
    session: Mapping[str, Any], *, source_features: Mapping[str, Any] | None = None
) -> SessionDraftMigration:
    """Run the Session 27-44 draft migrations in Web Load's order.

    A Session 40 or 41 the CLI wrote first drops its option-derived ``config``
    (``_holds_cli_writer_config``): it holds no Web draft. Then field names
    and placement rows (``migrate_persisted_web_state_field_names``), the
    draft's option values, slots and shapes (``migrate_session_draft_values``),
    then per-feature edits (``migrate_session_feature_edits``) through the
    saved catalog or, without one, through the request records and the first
    non-empty of ``source_features["extractedFeatures"]`` (the sources read
    again, see :func:`gbdraw.session_migration.read_legacy_source_features`)
    and the saved feature metadata, then the ``hash=`` annotation targets of a
    Session 40-44 (``migrate_session_annotation_targets``). A Session 40-44
    without a draft takes the annotation sets of its request first, as Web
    Load builds its draft from the request, and keeps them as its draft only
    when a target moved.
    """

    migrated: dict[str, Any] = dict(session)
    if _holds_cli_writer_config(session):
        del migrated["config"]
    version = session.get("version")
    version = version if isinstance(version, int) else 0
    request_value = session.get("renderRequest")
    request: Mapping[str, Any] = request_value if isinstance(request_value, Mapping) else {}
    editor_state = session.get("editorState")
    catalog = editor_state.get("featureCatalog") if isinstance(editor_state, Mapping) else None
    config = migrated.get("config")
    if isinstance(config, Mapping):
        config = migrate_session_draft_values(migrate_persisted_web_state_field_names(config), version)
        migrated["config"] = config
    features = session.get("features")
    dropped = narrowed = 0
    if isinstance(features, Mapping):
        has_catalog = isinstance(catalog, Mapping)
        read_again = source_features if isinstance(source_features, Mapping) else {}
        legacy = None if has_catalog else {
            "records": _promoted_request_records(request),
            "features": next(
                (
                    candidates
                    for candidates in (
                        read_again.get("extractedFeatures"),
                        features.get("biologicalFeatures"),
                        features.get("extractedFeatures"),
                    )
                    if isinstance(candidates, list) and candidates
                ),
                [],
            ),
            "biologicalFeatures": read_again.get("biologicalFeatures") or [],
        }
        edits = migrate_session_feature_edits(
            features, mode=request.get("mode"), catalog=catalog if has_catalog else None, legacy=legacy
        )
        migrated_features = edits.features
        # The Web split keeps no empty draft that the Session did not save.
        for key in ("labelOverrideRows",) if has_catalog else ("labelOverrideRows", "featureOverrides"):
            if key not in features and migrated_features.get(key) in ([], {}):
                migrated_features.pop(key)
        migrated["features"] = migrated_features
        dropped, narrowed = edits.dropped_count, edits.narrowed_visibility_count
    moved = 0
    if CURRENT_AUTHORITY_SESSION_MIN_VERSION <= version < MODE_SCOPED_SESSION_MIN_VERSION:
        mode = _diagram_mode(request.get("mode"))
        options = request.get("diagramOptions")
        annotations = options.get("annotations") if isinstance(options, Mapping) else None
        draft_sets: object = (
            config.get("annotationSets")
            if isinstance(config, Mapping)
            else _draft_annotation_sets_of_request(annotations.get("sets"), mode)
            if mode is not None and isinstance(annotations, Mapping)
            else None
        )
        targets = migrate_session_annotation_targets(
            draft_sets, mode=request.get("mode"), catalog=catalog, records=request.get("records")
        )
        moved = targets.migrated_count
        if moved:
            migrated["config"] = {
                **(config if isinstance(config, Mapping) else {}),
                "annotationSets": targets.annotation_sets,
            }
    return SessionDraftMigration(migrated, dropped, narrowed, moved)


def validate_current_web_state_field_names(
    config: object,
    *,
    include_linear_label_visibility: bool = True,
) -> None:
    """Reject obsolete Web config names at current session write boundaries."""

    if not isinstance(config, Mapping):
        return
    adv = config.get("adv")
    if isinstance(adv, Mapping):
        if include_linear_label_visibility:
            for field in ("linear_show_accession", "linear_show_length"):
                if field in adv:
                    raise ValidationError(
                        f"Web state field adv.{field} is obsolete; "
                        "use the selected visibility mode."
                    )
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

    request = session.get("renderRequest")
    editor = session.get("editorState")
    if (isinstance(request, Mapping)
            and request.get("layout", {}).get("similarityAlignment") is not None
            and isinstance(editor, dict)):
        editor.setdefault("alignmentResetReceipt", None)

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
        file_by_id: dict[str, Any] = {
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
    legacy_binding_entries: list[dict[str, Any]] = [
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
    file_by_id = {}
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
        assert explicit is not None
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


def _embedded_resource_bytes(entry: Mapping[str, Any]) -> bytes:
    """The bytes of a resource descriptor, checked against its size and checksum.

    The one check of a declared ``checksum``: the SHA-256 hex digest of the
    bytes, bare or as ``sha256:<hex>``, in any case. Null or empty declares
    none, as in the Web reader (``services/session-resource-backing.js``).
    """

    data = _embedded_entry_bytes(entry)
    if data is None or len(data) != entry.get("size"):
        raise ValidationError("Invalid embedded resource bytes or byte size.")
    checksum = entry.get("checksum")
    if checksum not in (None, "") and (
        not isinstance(checksum, str)
        or checksum.strip().lower().removeprefix("sha256:") != hashlib.sha256(data).hexdigest()
    ):
        raise ValidationError("Embedded resource checksum does not match.")
    return data


def _project_web_file_binding(
    resources: Mapping[str, Any],
    binding: Any,
    *,
    schema: int | None,
) -> Any:
    """Transport an admitted binding without choosing destination identities.

    Schema None is the existing direct-source metadata-default context.
    Source descriptors remain alive, unmodified and encoded until assembly.
    """
    if binding is None:
        return None
    if isinstance(binding, list):
        return [_project_web_file_binding(resources, item, schema=schema) for item in binding]
    if binding.get("kind") == "composite":
        return {
            **binding,
            "components": [
                _project_web_file_binding(resources, part, schema=schema)
                for part in binding["components"]
            ],
        }
    resource_id = binding["resourceId"]
    resource = resources.get(resource_id)
    if not isinstance(resource, Mapping):
        raise ValidationError(f"Web file binding references a missing resource: {resource_id}.")
    metadata = (
        {key: binding[key] for key in ("name", "type", "lastModified")}
        if schema == 2 else {
            "name": str(binding.get("name") or resource.get("name") or "file"),
            "type": str(binding.get("type") or resource.get("type") or ""),
            "lastModified": int(binding.get("lastModified") or resource.get("lastModified") or 0),
        }
    )
    return {"resourceId": resource_id, "descriptor": resource, **metadata}


def _attach_current_web_file_bindings(
    payload: dict[str, Any],
    files: Mapping[str, Any],
) -> None:
    from .session_resources import SessionResourceTable

    resources_value = payload.get("resources")
    if not isinstance(resources_value, dict):
        raise ValidationError("Current session resources must be an object.")
    resources = resources_value
    for resource_id in resources:
        if not isinstance(resource_id, str) or resource_id != resource_id.strip() or not resource_id:
            raise ValidationError("Canonical resource IDs must be unique non-empty strings.")
    table = SessionResourceTable(resources)

    def allocate_file(entry: Mapping[str, Any], preferred_id: str, metadata: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "resourceId": table.bind(entry, preferred_id=preferred_id),
            "name": str(metadata.get("name", "file")),
            "type": str(metadata.get("type") or ""),
            "lastModified": metadata.get("lastModified", 0),
        }

    def add_file(value: Any) -> dict[str, Any] | None:
        if value is None:
            return None
        if not isinstance(value, Mapping) or "components" in value or value.get("kind") == "composite":
            raise ValidationError("Expected an ordinary Web inventory file.")
        if "descriptor" in value:
            if set(value) != {"descriptor", "resourceId", "name", "type", "lastModified"} or not isinstance(value["descriptor"], Mapping):
                raise ValidationError("Invalid resource-backed Web inventory file.")
            return allocate_file(value["descriptor"], value["resourceId"], value)
        if "resourceId" in value:
            raise ValidationError("Web inventory source references require a descriptor.")
        return allocate_file(value, "", value)

    def add_value(value: Any, *, composite_allowed: bool = False) -> Any:
        if isinstance(value, list):
            return [add_value(item) for item in value]
        if isinstance(value, Mapping) and value.get("kind") == "composite":
            if (not composite_allowed
                    or set(value) != {"kind", "components", "name", "type", "lastModified"}
                    or not isinstance(value["components"], list)
                    or len(value["components"]) < 2
                    or any(part is None for part in value["components"])):
                raise ValidationError("Invalid composite Web inventory file.")
            return {**value, "components": [add_file(part) for part in value["components"]]}
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
        "schema": 2,
        "c_gb": add_value(files.get("c_gb"), composite_allowed=True),
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
    resources.update(table.descriptors())
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
    explicit_names: dict[str, str] = {}

    def reference_value(value: Any, source: Any) -> Any:
        if isinstance(value, list):
            return [reference_value(item, source[index]) for index, item in enumerate(value)]
        if value is None:
            return None
        if isinstance(source, Mapping) and "resourceId" in source:
            resource_aliases[source["resourceId"]] = value["resourceId"]
        explicit_names.setdefault(value["resourceId"], value["name"])
        return value["resourceId"]

    for source_field, binding_field in (
        ("conservationLosatFastaSources", "c_conservation_fastas"),
        ("conservationSequenceSources", "c_conservation_sequence_sources"),
    ):
        if source_field not in web_files:
            continue
        rebound = bindings[binding_field]
        source = files.get(binding_field)
        web_files[source_field] = reference_value(
            rebound if isinstance(rebound, list) else [rebound] if rebound is not None else [],
            source if isinstance(source, list) else [source] if source is not None else [],
        )
    original_names = web_files.get("resourceOriginalNames")
    if isinstance(original_names, Mapping):
        web_files["resourceOriginalNames"] = {
            **{
                resource_aliases.get(str(resource_id), str(resource_id)): name
                for resource_id, name in original_names.items()
                if resource_aliases.get(str(resource_id), str(resource_id)) in resources
            },
            **explicit_names,
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
    source_composite = None
    if context.source_session is not None:
        validate_session(context.source_session)
        source_version = int(context.source_session["version"])
        source_web_files = context.source_session.get("webFiles")
        source_bindings = source_web_files.get("bindings") if isinstance(source_web_files, Mapping) else None
        explicit = source_bindings.get("c_gb") if isinstance(source_bindings, Mapping) else None
        if (
            isinstance(source_bindings, Mapping)
            and isinstance(explicit, Mapping)
            and explicit.get("kind") == "composite"
        ):
            from .session import _validate_document

            _validate_document(context.source_session)
            source_composite = _project_web_file_binding(
                context.source_session["resources"], explicit, schema=source_bindings["schema"],
            )
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

    if source_version is not None and source_version < MODE_SCOPED_SESSION_MIN_VERSION:
        payload = migrate_session_flat_draft(payload).session
    config = payload.get("config")
    if not isinstance(config, dict):
        config = dict(config) if isinstance(config, Mapping) else {}

    ui = payload.get("ui")
    ui = dict(ui) if isinstance(ui, Mapping) else {}
    payload["ui"] = ui
    ui["mode"] = context.mode
    ui.pop("blastSource", None)
    ui.setdefault("zoom", 1)
    ui.setdefault("selectedResultIndex", 0)
    ui.setdefault("canvasPan", {"x": 0, "y": 0})

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

    orthogroup_state = payload.get("orthogroupState")
    orthogroup_state = (
        dict(orthogroup_state)
        if isinstance(orthogroup_state, Mapping)
        else {}
    )
    orthogroup_state.pop("groups", None)
    orthogroup_state.pop("selectedOrthogroupAlignmentFeature", None)
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
    if source_composite is not None:
        files_for_web = {**files_for_web, "c_gb": source_composite}
    if source_version is not None and source_version < CURRENT_AUTHORITY_SESSION_MIN_VERSION:
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
    if source_version is not None and source_version < MODE_SCOPED_SESSION_MIN_VERSION:
        # The flat draft of an older source goes to the mode slices (Session 46).
        payload = split_draft_into_modes(
            payload,
            committed_mode=context.mode,
            depth_sources=session_depth_source_widths(files_for_web),
            palette_colors=mode_split_palette_colors(payload.get("config")),
        )
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
    _write_validated_session_json(path, expanded_payload, overwrite=overwrite)


def _write_validated_session_json(
    path: str | Path,
    expanded_payload: Mapping[str, Any],
    *,
    overwrite: bool,
) -> None:
    """Write a session that ``validate_session`` already accepted.

    Only writers whose payload cannot have changed since its validation call
    this: ``write_session_json`` after validating, a ``SessionDocument`` (validated
    when it was built), and the CLI sidecar, which writes what
    ``build_session_json`` returned.
    """

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
        text_file: TextIO
        if output_path.suffix.lower() == ".gz":
            raw_file = os.fdopen(temp_fd, "wb")
            temp_fd = None
            with raw_file:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw_file,
                    compresslevel=6,
                    mtime=0,
                ) as compressed_file:
                    # Streamed: the gzip bytes depend on the write chunks.
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
                # json.dumps encodes in C at once; json.dump encodes in Python.
                text_file.write(
                    json.dumps(
                        serialized_payload,
                        ensure_ascii=False,
                        separators=(",", ":"),
                    )
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


@dataclass(frozen=True)
class RetiredCliOption:
    """One retired CLI flag (design D4).

    Fresh runs reject the flag and name ``replacement``. Legacy session argv
    is rewritten before replay: ``renamed_to`` keeps the value, and
    ``value_rewrites`` replaces flag and value with zero or more tokens. A flag
    with neither has no legacy rewrite.
    """

    option: str
    modes: tuple[Literal["circular", "linear"], ...]
    replacement: str
    renamed_to: str | None = None
    value_rewrites: Mapping[str, tuple[str, ...]] | None = None
    # argparse nargs of the retired flag, so a fresh run reports the flag
    # itself rather than its extra values. A rename keeps every value.
    nargs: str | None = None

    def rewrite(self, value: str) -> tuple[str, ...] | None:
        if self.renamed_to is not None:
            return (self.renamed_to, value)
        if self.value_rewrites is not None:
            return self.value_rewrites.get(str(value).strip().lower())
        return None

    def message(self, value: object = None) -> str:
        replacement = self.replacement
        if self.value_rewrites is not None and isinstance(value, str):
            tokens = self.value_rewrites.get(value.strip().lower())
            if tokens is not None:
                replacement = (
                    "use " + " ".join(tokens)
                    if tokens
                    else "omit it (no protein comparison)"
                )
        return f"{self.option} was retired; {replacement}."


def _losatp_mode_rewrites() -> dict[str, tuple[str, ...]]:
    from gbdraw.api.options import LOSATP_MODE_WIRE

    rewrites: dict[str, tuple[str, ...]] = {"none": ()}
    for typed_mode, wire_mode in LOSATP_MODE_WIRE.items():
        if typed_mode != "none":
            rewrites[wire_mode] = ("--losat", "losatp", "--losatp_mode", typed_mode)
    return rewrites


# The one old -> new table for retired CLI flags (design 3.5). Both the
# fresh-run rejection and the legacy session argv rewrite read it.
RETIRED_CLI_OPTIONS: Mapping[str, RetiredCliOption] = {
    item.option: item
    for item in (
        RetiredCliOption(
            "--protein_blastp_mode",
            ("linear",),
            "use --losat losatp --losatp_mode {similarity_groups,collinear,pairwise}",
            value_rewrites=_losatp_mode_rewrites(),
        ),
        RetiredCliOption("--losatp_bin", ("linear",), "use --losat_bin", "--losat_bin"),
        RetiredCliOption(
            "--ncbi_blastp_bin", ("linear",), "use --ncbi_blast_bin", "--ncbi_blast_bin"
        ),
        RetiredCliOption(
            "--losatp_threads", ("linear",), "use --losat_threads", "--losat_threads"
        ),
        RetiredCliOption(
            "--protein_blastp_max_hits",
            ("linear",),
            "use --losatp_max_hits",
            "--losatp_max_hits",
        ),
        RetiredCliOption(
            "--protein_blastp_candidate_limit",
            ("linear",),
            "use --losatp_max_target_seqs",
            "--losatp_max_target_seqs",
        ),
        RetiredCliOption(
            "--align_orthogroup_feature",
            ("linear",),
            "use --similarity_alignment_feature",
            "--similarity_alignment_feature",
        ),
        RetiredCliOption(
            "--protein_blastp_output",
            ("linear",),
            "use --losat_output_dir DIR, which writes DIR/losatp.raw.tsv",
        ),
        RetiredCliOption(
            "--conservation_fasta",
            ("circular",),
            "use --conservation_sequence (FASTA, GenBank, or DDBJ)",
            "--conservation_sequence",
            nargs="+",
        ),
    )
}


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

    def extend(option_index: int, value_index: int, tokens: Sequence[str]) -> None:
        if not tokens:
            return
        source_to_canonical_index[option_index] = len(canonical_args)
        source_to_canonical_index[value_index] = len(canonical_args) + len(tokens) - 1
        canonical_args.extend(tokens)

    pending_retired: tuple[RetiredCliOption, int] | None = None
    for source_index, raw_token in enumerate(args):
        token = str(raw_token)
        if pending_retired is not None:
            retired, option_index = pending_retired
            pending_retired = None
            rewritten = retired.rewrite(token)
            extend(
                option_index,
                source_index,
                (retired.option, token) if rewritten is None else rewritten,
            )
            continue
        if token == "--show_depth":
            continue
        option, separator, inline_value = (
            token.partition("=")
            if token.startswith("--")
            else (token, "", "")
        )
        option = replacements.get(option, option)
        retired_option = RETIRED_CLI_OPTIONS.get(option)
        if retired_option is not None and mode in retired_option.modes:
            if not separator:
                pending_retired = (retired_option, source_index)
                continue
            rewritten = retired_option.rewrite(inline_value)
            if rewritten is None:
                rewritten = (f"{option}={inline_value}",)
            elif retired_option.renamed_to is not None:
                rewritten = (f"{retired_option.renamed_to}={inline_value}",)
            extend(source_index, source_index, rewritten)
            continue
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
    if pending_retired is not None:
        retired, option_index = pending_retired
        source_to_canonical_index[option_index] = len(canonical_args)
        canonical_args.append(retired.option)
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


def canonicalize_cli_invocation(
    args: Sequence[str],
    file_bindings: Sequence[SessionFileBinding],
    *,
    mode: Literal["circular", "linear"],
) -> tuple[list[str], list[SessionFileBinding]]:
    """Rewrite retired CLI flags in a recorded invocation and remap its bindings.

    Legacy session replay and the Gallery session refresh both use this, so a
    recorded ``cliInvocation`` reaches the current flag names the same way.
    """

    canonical_args, index_map = _canonicalize_legacy_session_cli_args(args, mode=mode)
    remapped: list[SessionFileBinding] = []
    for binding in file_bindings:
        if binding.argIndex not in index_map:
            raise ValidationError(
                "cliInvocation.fileBindings cannot reference a removed legacy CLI flag."
            )
        remapped.append(
            SessionFileBinding(
                argIndex=index_map[binding.argIndex],
                slot=binding.slot,
                name=binding.name,
            )
        )
    return canonical_args, remapped


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
    migrate_legacy_cli = session_version < CURRENT_AUTHORITY_SESSION_MIN_VERSION
    _restore_cli_table_paths(
        session,
        run_args,
        temp_dir=temp_dir,
        migrate_legacy_cli=migrate_legacy_cli,
    )

    if migrate_legacy_cli:
        run_args, _ = _canonicalize_legacy_session_cli_args(run_args, mode=mode)
        invocation_args, file_bindings = canonicalize_cli_invocation(
            invocation_args,
            file_bindings,
            mode=mode,
        )

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
    form_value = config.get("form")
    form = form_value if isinstance(form_value, Mapping) else {}
    adv_value = config.get("adv")
    adv = dict(adv_value) if isinstance(adv_value, Mapping) else {}
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
            arg_index = int(cast(Any, table_entry.get("argIndex")))
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
                row_index = int(cast(Any, dependency.get("rowIndex")))
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
                start_int = int(cast(Any, start))
                end_int = int(cast(Any, end))
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
    _append_pair(run_args, invocation_args, "--losat", "losatp")
    _append_pair(
        run_args,
        invocation_args,
        "--losatp_mode",
        _losatp_mode_rewrites()[mode][-1],
    )
    threads_per_job = str(losat_cfg.get("threadsPerJob") or "auto").strip().lower()
    if threads_per_job != "auto":
        try:
            parsed_threads = int(threads_per_job)
        except ValueError:
            parsed_threads = 0
        if parsed_threads >= 1:
            _append_pair(run_args, invocation_args, "--losat_threads", str(parsed_threads))

    max_hits = blastp_cfg.get("maxHits")
    if max_hits not in (None, "", False):
        _append_pair(run_args, invocation_args, "--losatp_max_hits", str(max_hits))
        if mode == "pairwise":
            _append_pair(run_args, invocation_args, "--losatp_max_target_seqs", str(max_hits))
    candidate_limit = blastp_cfg.get("candidateLimit")
    if mode != "pairwise" and candidate_limit not in (None, "", False):
        _append_pair(run_args, invocation_args, "--losatp_max_target_seqs", str(candidate_limit))

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
            _append_pair(run_args, invocation_args, "--similarity_alignment_feature", selected_target)

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
            and int(session.get("version", 0)) < CURRENT_AUTHORITY_SESSION_MIN_VERSION
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
            arg_index = int(cast(Any, item.get("argIndex")))
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
    "AnnotationTargetMigration",
    "CURRENT_SESSION_VERSION",
    "RETIRED_RENDERED_ID_FEATURE_FIELDS",
    "CANONICAL_SESSION_MIN_VERSION",
    "DEPTH_FILE_ENCODING",
    "DEPTH_FILE_SCHEMA",
    "FEATURE_CATALOG_ENCODING",
    "FEATURE_CATALOG_SCHEMA",
    "FeatureEditMigration",
    "LEGACY_LOSAT_DERIVED_CACHE_SCHEMA",
    "LEGACY_PROTEIN_CANDIDATE_SCHEMA",
    "LOSAT_DERIVED_CACHE_SCHEMA",
    "NUCLEOTIDE_LOSAT_CACHE_SCHEMA",
    "PROTEIN_IDENTITY_MANIFEST_SCHEMA",
    "PROTEIN_LOSAT_CACHE_SCHEMA",
    "MODE_SCOPED_SESSION_MIN_VERSION",
    "SESSION_FORMAT",
    "SUPPORTED_SESSION_VERSIONS",
    "SessionBuildContext",
    "SessionFileBinding",
    "SessionRunSpec",
    "build_session_json",
    "canonicalize_cli_invocation",
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
    "migrate_persisted_web_option_values",
    "migrate_persisted_web_state_field_names",
    "migrate_imported_linear_track_slots",
    "migrate_session_draft_values",
    "migrate_session_annotation_targets",
    "migrate_session_feature_edits",
    "migrate_legacy_repeat_feature_shape_args",
    "normalize_current_session_artifacts",
    "safe_embedded_filename",
    "serialize_file_entry",
    "session_depth_source_widths",
    "session_mode",
    "session_to_cli_args",
    "split_draft_into_modes",
    "mode_split_palette_colors",
    "migrate_session_flat_draft",
    "SessionDraftMigration",
    "validate_session",
    "validate_current_session_artifacts",
    "validate_current_web_state_field_names",
    "write_session_json",
]
