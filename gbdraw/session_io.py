#!/usr/bin/env python
# coding: utf-8

"""GUI session JSON loading, validation, materialization, and sidecar building."""

from __future__ import annotations

import base64
import binascii
import copy
import hashlib
import json
import math
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Mapping,
    MutableMapping,
    NamedTuple,
    NoReturn,
    Sequence,
    TypedDict,
    TypeGuard,
    cast,
)

from .analysis.protein_artifacts import (
    CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
    is_current_derived_protein_artifact,
    validate_current_derived_protein_artifacts,
)
from .analysis.protein_colinearity import (
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    is_json_exact_protein_identity_manifest,
)
from .comparisons.losat_jobs import (
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    losat_pair_filename,
)
from .core.numeric import (
    ASCII_NONNEGATIVE_INTEGER_PATTERN,
    JS_MAX_SAFE_INTEGER,
    JS_WHITESPACE,
    parse_ascii_nonnegative_integer,
)
from .core.js_values import (
    JS_ARRAY_INDEX,
    js_integer,
    js_integral,
    js_key_order,
    js_number,
    js_string,
    js_text,
    js_truthy,
    json_clone,
    mapping_list,
    nonnegative_integer,
    safe_integer,
)
from .io.input_text import GZIP_ERRORS, open_input_text
from .exceptions import GbdrawError, SessionFormatError, ValidationError
from .io.colors import is_user_color, named_color_hex, packaged_color_palettes
from .io.filenames import _SAFE_FILENAME_RE, safe_embedded_filename
from .session_container import (
    SESSION_FILE_UNREADABLE,
    SESSION_JSON_INVALID,
    SESSION_RESOURCES_INVALID,
    BytesResourceSource,
    FileResourceSource,
    SessionContainer,
    SessionContainerError,
    SessionResourceSource,
    invalid_session_json_message,
    is_session_container,
    parse_session_json,
    read_session_container,
    write_session_container,
)
from .session_legacy_values import (
    LEGACY_FEATURE_SHAPES,
    LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES,
    _option_value,
)
from .session_contract import (
    CURRENT_SESSION_VERSION,
    FeatureOverrideDraftRow,
    SessionContainerResource,
    FlatDraftFields,
    GbdrawSession,
    GbdrawSession51,
    LegacySessionFields,
    SessionOtherModeEditorState,
    SessionOtherModeLegend,
    SessionOtherModeResult,
    SessionAnnotationWarning,
    SessionComparisonWarning,
    SessionFeatureIdentityNotice,
    SessionOtherModeUi,
    SessionRunMetadata,
    keys_of,
    literal_values,
)
from .web_support.mode_scoped_settings import (
    DIAGRAM_MODES,
    LOSAT_EXECUTION_FIELDS,
    MODE_SCOPED_FIRST_SESSION_VERSION,
    MODE_SCOPED_SETTINGS,
    DiagramMode,
    ModeScopedSetting,
    slice_containers,
    unmanaged_config_override_modes,
)

if TYPE_CHECKING:
    from .analysis.protein_colinearity import ProteinIdentityManifest
    from .session_drawings import ModeScopedSourceFacts

SESSION_FORMAT = "gbdraw-session"
# Version 44 is the first with the current active-config and record-display
# draft shapes; version 46 keeps the Web draft of each diagram mode in
# ``modes`` (PD-OI-090) and keys per-feature edits by source identity;
# version 47 keys a slice's Legend color and stroke edits by the row key
# (``originalCaption``) instead of the shown caption (PD-OI-061 revision 5)
# and adds each slice's ``adv.font_family`` (PD-OI-098); version 48 stores
# the protein results of a reverse-complemented Linear record in its search
# frame (OV-1145) and what each file card's file put into a slice, for Remove
# on the card (OV-1280); version 49 adds each Circular slice's
# ``adv.show_definition`` (PD-OI-109); version 50 lists the drawings of a
# project (below) and adds each Linear drawing's ``form.linear_wrap_length``
# and the config leaf ``canvas.linear.wrap_length`` (a record drawn over
# several rows).
TYPED_DRAFT_SESSION_MIN_VERSION = 44
MODE_SCOPED_SESSION_MIN_VERSION = MODE_SCOPED_FIRST_SESSION_VERSION
ROW_KEYED_LEGEND_SESSION_MIN_VERSION = 47
# From version 48 a Session without a Web draft (CLI, Python API) holds each
# table only as its file card's file: the Web reads every row as that file's.
FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION = 48
# From version 50 a Session lists its drawings in ``drawings``
# (gbdraw.session_drawings owns the layout).
DRAWINGS_SESSION_MIN_VERSION = 51
# The first Session whose Linear comparison-plan edges may hold `windowLosatFilenames`, the names typed on the
# window pairs of a record pair (OV-2742). The Web drops the field from an older file on Load.
WINDOW_LOSAT_FILENAMES_SESSION_VERSION = 51
# The newest mode-scoped (46-50) layout: a Session 46-50 moves from it to
# drawings, and Python's writers build a new document in it before that move
# (gbdraw.session_drawings.drawings_document_from_mode_scoped).
MODE_SCOPED_SESSION_MAX_VERSION = DRAWINGS_SESSION_MIN_VERSION - 1
# Every mode-scoped (``modes``) version, read without Python.
MODE_SCOPED_SESSION_VERSIONS = frozenset(range(MODE_SCOPED_SESSION_MIN_VERSION, MODE_SCOPED_SESSION_MAX_VERSION + 1))
# From version 41 the record display and feature placement drafts of a flat draft are checked.
DISPLAY_PLACEMENT_DRAFT_SESSION_MIN_VERSION = 41
CURRENT_AUTHORITY_SESSION_MIN_VERSION = 40
CANONICAL_SESSION_MIN_VERSION = 31
SUPPORTED_SESSION_VERSIONS = frozenset(
    {27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 44, DRAWINGS_SESSION_MIN_VERSION} | MODE_SCOPED_SESSION_VERSIONS
)
CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39
LOSAT_DERIVED_CACHE_SCHEMA = CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA
LEGACY_LOSAT_DERIVED_CACHE_SCHEMA = 1
LEGACY_PROTEIN_CANDIDATE_SCHEMA = 1
FEATURE_CATALOG_SCHEMA = 1
FEATURE_CATALOG_ENCODING = "biological-authority-v1"
# Schema 6 stores no sequence: Python reads them when a reader needs them (D5). Schema 5, written by Sessions 46-50,
# stored them; it is promoted when read, and it stays the form an exported Interactive SVG embeds.
CURRENT_FEATURE_CATALOG_SCHEMA = 6
SEQUENCE_FEATURE_CATALOG_SCHEMA = 5
# The catalog schema each current-authority Session version writes.
FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION = {
    44: 4,
    **{version: SEQUENCE_FEATURE_CATALOG_SCHEMA for version in MODE_SCOPED_SESSION_VERSIONS},
    DRAWINGS_SESSION_MIN_VERSION: CURRENT_FEATURE_CATALOG_SCHEMA,
}
# The fields of a record display draft row (Session 41+), the two that Session
# 44 added, and the fields of its anchor intent.
RECORD_DISPLAY_DRAFT_FIELDS = frozenset({"sourceUid", "selector", "recordId", "topologyOverride", "startCoordinate"})
TYPED_RECORD_DISPLAY_DRAFT_FIELDS = frozenset({"reverseComplementOverride", "anchorIntent"})
RECORD_ANCHOR_INTENT_FIELDS = frozenset(
    {"schema", "recordKey", "biologicalFeatureId", "placement", "anchor", "offsetBp", "orientForward"}
)
# The field sets below are the keys of the Session document's declaration
# (gbdraw.session_contract).
# Every top-level key a Session of a supported version may hold.
CURRENT_SESSION_TOP_LEVEL_FIELDS = keys_of(GbdrawSession) | keys_of(LegacySessionFields)
# The flat Web draft of Sessions 44 and older; Session 46 keeps it per mode in
# ``modes`` (MODE_SCOPED_SETTINGS).
FLAT_DRAFT_TOP_LEVEL_FIELDS = keys_of(FlatDraftFields)
# Session 46 keeps the other diagram mode's Result set in ``otherModeResult``:
# the fields of one committed set, named as at the top level.
OTHER_MODE_RESULT_FIELDS = keys_of(SessionOtherModeResult)
# The per-set part of the shared ``ui`` and ``editorState`` objects.
OTHER_MODE_RESULT_UI_FIELDS = keys_of(SessionOtherModeUi)
OTHER_MODE_RESULT_EDITOR_FIELDS = keys_of(SessionOtherModeEditorState)
OTHER_MODE_RESULT_LEGEND_FIELDS = keys_of(SessionOtherModeLegend)
OTHER_MODE_RESULT_RUN_METADATA_FIELDS = keys_of(SessionRunMetadata)
# The rows of the run notice lists in ``runMetadata``: their keys and the values of their fixed fields, which
# ``_validate_other_mode_run_metadata`` and the Web's validators take.
ANNOTATION_WARNING_FIELDS = keys_of(SessionAnnotationWarning)
ANNOTATION_WARNING_CODES = literal_values(SessionAnnotationWarning, "code")
COMPARISON_WARNING_FIELDS = keys_of(SessionComparisonWarning)
COMPARISON_WARNING_CODES = literal_values(SessionComparisonWarning, "code")
FEATURE_IDENTITY_NOTICE_FIELDS = keys_of(SessionFeatureIdentityNotice)
FEATURE_IDENTITY_NOTICE_STATUSES = literal_values(SessionFeatureIdentityNotice, "status")
FEATURE_IDENTITY_NOTICE_KINDS = literal_values(SessionFeatureIdentityNotice, "kinds")
# Web file binding slots (``webFiles.bindings``, schema 2): the five file
# cards, which every drawing binds, and the slots of each diagram mode's inputs
# in a Session 51 drawing. A Circular drawing's input files are its File cards,
# ``circularSources`` (``{uid, gb, gff, fasta}`` rows, like ``linearSeqs``).
WEB_BINDING_CARD_SLOTS = frozenset({"d_color", "t_color", "blacklist", "whitelist", "qualifier_priority"})
WEB_BINDING_MODE_SLOTS: Mapping[str, frozenset[str]] = {
    "circular": frozenset({
        "circularSources", "c_depth", "c_conservation_blasts",
        "c_conservation_blasts_source", "c_conservation_fastas", "c_conservation_sequence_sources",
    }),
    "linear": frozenset({"linearSeqs", "linearComparisons"}),
}
# The slot that holds a mode's main input; schema 2 bindings always hold it.
WEB_BINDING_INPUT_SLOT: Mapping[str, str] = {"circular": "circularSources", "linear": "linearSeqs"}
# The one Circular input of Sessions 41-50, in place of ``circularSources``:
# ``c_gb`` (one GenBank file, or a ``composite`` of several) or ``c_gff`` and
# ``c_fasta``. The move to drawings makes File cards of them
# (:func:`gbdraw.session_drawings.circular_sources_of_mode_scoped`).
MODE_SCOPED_CIRCULAR_INPUT_SLOTS = frozenset({"c_gb", "c_gff", "c_fasta"})


class CircularSourceCard(TypedDict):
    """One Circular File card, a ``circularSources`` row: its uid and the
    file of each slot (a binding, or in a Web inventory a file entry) or None."""

    uid: str
    gb: Mapping[str, object] | None
    gff: Mapping[str, object] | None
    fasta: Mapping[str, object] | None


# The fields of one Circular File card (``circularSources`` row).
CIRCULAR_SOURCE_FIELDS = CircularSourceCard.__required_keys__
# The file slots of a Circular File card, every field but its uid, in the
# card's order; the Web reads them from the generated contracts.
CIRCULAR_SOURCE_FILE_SLOTS: tuple[str, ...] = tuple(field for field in CircularSourceCard.__annotations__ if field != "uid")
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
FEATURE_OVERRIDE_DRAFT_FIELDS = keys_of(FeatureOverrideDraftRow)
# A Session 46 field that names an invalid shape is a Session field error.
_SESSION_FIELDS_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}
DEPTH_FILE_ENCODING = "gbdraw-depth-table-v1"
DEPTH_FILE_SCHEMA = 1

_DEPTH_COLUMNS = ("reference_name", "position", "depth")
_SLOT_PART_RE = re.compile(rf"([^\[\]]+)|\[({ASCII_NONNEGATIVE_INTEGER_PATTERN})\]")
# A Circular source record ``#N`` (N from 1, ASCII digits as the Web reads it,
# OV-611): a record display draft's selector and a slice's OFF record key.
_SOURCE_RECORD_SELECTOR = re.compile(r"#[1-9][0-9]*")


@dataclass(frozen=True)
class SessionAdmission:
    """The admission checks of :func:`validate_session` that a reading runs.

    The CLI runs them all (:data:`CLI_SESSION_ADMISSION`). Web Load of a
    Session 27-44 reads the file with the Web's own readers first, which
    repair or drop the parts they refuse, and Python then converts what they
    kept (:data:`WEB_READ_SESSION_ADMISSION`, D-15-131): it does not check
    again that Results hold their feature catalog (a Result without one is
    set aside) or that protein LOSAT caches resolve through the manifest.
    """

    feature_catalog: bool
    losat_manifest: bool


CLI_SESSION_ADMISSION = SessionAdmission(feature_catalog=True, losat_manifest=True)
WEB_READ_SESSION_ADMISSION = SessionAdmission(feature_catalog=False, losat_manifest=False)


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
    # What the run reported beside its Results (``gbdraw.api.run_metadata``,
    # the producer of what a Web Save writes); the Session keeps it as it is,
    # in place of any ``runMetadata`` a source Session held.
    run_metadata: SessionRunMetadata | None = None
    # The bytes' sources of ``source_session``'s resources.
    source_resources: Mapping[str, SessionResourceSource] | None = None


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
    """Load and validate a Session file, or an older JSON Session; return its document.

    It validates the document as ``gbdraw.session.load_session_document`` does,
    canonical resource descriptors included. The document is JSON, without
    its resources' bytes; ``load_session_document`` keeps them.
    """

    from .session import load_session_document

    return load_session_document(path).to_dict()


def _read_session_json(path: str | Path) -> tuple[object, dict[str, SessionResourceSource]]:
    """Parse a Session file and give the sources of its resources' bytes; every
    unreadable file is a ``ValidationError``.

    The content decides, never the name: a ZIP is a Session file
    (``session_container``), whose resources' bytes stay in their entries;
    anything else is a JSON Session, plain or gzip, of a version before the
    current one, whose embedded resources give their bytes where the
    document is read (:func:`split_resource_bytes`). A truncated or corrupt
    gzip stream and text that is not UTF-8 are "not a valid JSON session
    file", as malformed JSON is (OV-1251); a duplicate key raises its own
    ``ValidationError``.
    """

    session_path = Path(path)
    try:
        with open(session_path, "rb") as session_file:
            head = session_file.read(4)
        if is_session_container(head):
            container = read_session_container(session_path)
            check_session_file_version(container)
            return container.document, {
                resource_id: container.resource(resource_id) for resource_id in container.resource_ids()
            }
        payload = parse_session_json(_read_session_text(session_path), session_path)
    except (UnicodeDecodeError, *GZIP_ERRORS) as exc:
        raise ValidationError(invalid_session_json_message(session_path), diagnostic=SESSION_JSON_INVALID) from exc
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {session_path}", diagnostic=SESSION_FILE_UNREADABLE) from exc
    if isinstance(payload, dict) and payload.get("version") == CURRENT_SESSION_VERSION:
        # The current version is written as a Session file only; its JSON form never reached main.
        raise ValidationError(retired_json_session_message(session_path), diagnostic=SESSION_RESOURCES_INVALID)
    return payload, {}


def check_session_file_version(container: SessionContainer) -> None:
    """Refuse a Session file whose document is older than the current version.

    Only the current version was ever written as a Session file, so an older
    document in one is refused as an unsupported version; a newer one is left
    to the Session admission, which refuses it as too new. The Web's twin is
    ``assertSessionFileVersion`` in services/session-authority.js.
    """

    document = container.document
    version = document.get("version") if isinstance(document, dict) else None
    if isinstance(version, int) and not isinstance(version, bool) and version < CURRENT_SESSION_VERSION:
        raise SessionContainerError("unsupported-version", container.path, version=version)


def _read_session_text(path: str | Path) -> str:
    """Read UTF-8 session JSON, plain or gzip (told by its signature, as every input file is), with no size cap."""

    with open_input_text(path, encoding="utf-8", capped=False) as session_file:
        return session_file.read()


def validate_session(
    session: Mapping[str, object],
    admission: SessionAdmission = CLI_SESSION_ADMISSION,
    *,
    resource_sources: Mapping[str, SessionResourceSource] | None = None,
) -> None:
    """Validate the conservative session envelope used by the CLI, with the
    checks of ``admission``.

    ``resource_sources`` hold the bytes of the document's resources; a JSON
    Session 27-50 as it was read (embedded resources) gives its own
    (:func:`split_resource_bytes`).
    """

    embedded = validate_session_resources(session, admission, resource_sources)
    if embedded is not None and embedded.refused is not None:
        raise embedded.refused


class EmbeddedResources(NamedTuple):
    """The resources of a checked JSON Session 27-50, decoded once.

    ``document`` holds JSON descriptors and ``sources`` their bytes
    (:func:`split_resource_bytes`); ``refused`` is the refusal of bytes that do
    not match their descriptor, which the caller raises after its own checks
    of the stored fields. ``document`` is a Session 27-50 as it was stored,
    unvalidated beyond its envelope: a JSON object of that version.
    """

    document: dict[str, object]
    sources: dict[str, SessionResourceSource]
    refused: ValidationError | None


def validate_session_resources(
    session: Mapping[str, object],
    admission: SessionAdmission,
    resource_sources: Mapping[str, SessionResourceSource] | None,
) -> EmbeddedResources | None:
    """:func:`validate_session`, without raising the refusal of embedded bytes;
    for a JSON Session 27-50 whose resources embed their bytes, also those
    resources (:class:`EmbeddedResources`)."""

    if not isinstance(session, Mapping):
        raise ValidationError("Session JSON must be an object.")
    if session.get("format") != SESSION_FORMAT:
        raise ValidationError("Not a gbdraw-session JSON file.")
    version = session.get("version")
    if not isinstance(version, int):
        raise ValidationError("Session version is required and must be an integer.")
    if version > max(SUPPORTED_SESSION_VERSIONS):
        raise ValidationError(
            f"Session version {version} is newer than this gbdraw supports "
            f"({CURRENT_SESSION_VERSION})."
        )
    if version not in SUPPORTED_SESSION_VERSIONS:
        raise ValidationError(f"Unsupported session version: {version}.")
    refused_bytes: ValidationError | None = None
    split: tuple[dict[str, object], dict[str, SessionResourceSource]] = (dict(session), {})
    embeds = False
    declared = session.get("resources")
    if resource_sources is None and isinstance(declared, Mapping) and any(
        isinstance(entry, Mapping) and ("encoding" in entry or "data" in entry)
        for entry in declared.values()
    ):
        embeds = True
        # A JSON Session 27-50 is checked as it was written, its bytes embedded;
        # bytes that do not match their descriptor are refused after the
        # document's own checks, which name the field a reader can fix.
        try:
            split = split_resource_bytes(session)
            resource_sources = split[1]
        except ValidationError as exc:
            refused_bytes = exc
            resource_sources = _decodable_resource_sources(declared)
    sources = resource_sources or {}
    if version >= DRAWINGS_SESSION_MIN_VERSION:
        from .session_drawings import checked_drawing_views

        for view in checked_drawing_views(session):
            _validate_document_body(view, version, admission, sources)
    else:
        _validate_document_body(session, version, admission, sources)
    return EmbeddedResources(split[0], split[1], refused_bytes) if embeds else None


def _validate_document_body(
    session: Mapping[str, object],
    version: int,
    admission: SessionAdmission,
    sources: Mapping[str, SessionResourceSource],
) -> None:
    """Validate a Session 27-50, or one Session 51 drawing seen as one."""

    if "otherModeResult" in session and version < MODE_SCOPED_SESSION_MIN_VERSION:
        raise ValidationError(
            f"Session version {version} cannot contain otherModeResult.", diagnostic=_OTHER_MODE_RESULT_INVALID
        )
    if version >= CANONICAL_SESSION_MIN_VERSION:
        render_request = session.get("renderRequest")
        resources = session.get("resources")
        settings_only = is_settings_only_session(session)
        # A Session 51 drawing not generated yet has no request and no Results.
        undrawn = version >= DRAWINGS_SESSION_MIN_VERSION and render_request is None
        if not settings_only and not undrawn:
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
        validate_current_session_artifacts(session, admission)
    if version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION:
        _validate_current_top_level_fields(session)
        _validate_mode_scoped_fields(session, version)
        _validate_current_retired_active_config_paths(session, version)
        _validate_current_comparison_authority(session, version)
        _validate_current_feature_payloads(session, version)
        if admission.feature_catalog:
            _validate_current_feature_catalog_authority(session, version)
        _validate_alignment_reset_receipt(session, sources)
        _validate_other_mode_result(session, version, sources)
    if version >= MODE_SCOPED_SESSION_MIN_VERSION:
        _validate_mode_slices(session, version)
    elif version >= DISPLAY_PLACEMENT_DRAFT_SESSION_MIN_VERSION:
        _validate_display_placement_drafts(session)
    if is_settings_only_session(session):
        _validate_settings_only_session(session)


def _is_web_draft_version(version: object, *older: int) -> bool:
    """Whether ``version`` is one of ``older`` or a mode-scoped (``modes``) version."""
    return not isinstance(version, bool) and (version in older or version in MODE_SCOPED_SESSION_VERSIONS)


def is_settings_only_session(session: Mapping[str, Any]) -> bool:
    """Recognize the explicit document variant, never a missing-resource error."""
    return (
        _is_web_draft_version(session.get("version"), 42, 44)
        and "renderRequest" in session
        and session["renderRequest"] is None
    )


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


def _validate_circular_source_rows(rows: object, leaf: Callable[[object], None]) -> None:
    """A Circular drawing's File cards: a list of ``{uid, gb, gff, fasta}``
    rows with distinct non-empty uids, each file a binding (``leaf``) or null.

    Load admits a card without the file its drawing's input type reads;
    Generate asks for it.
    """

    if not isinstance(rows, list):
        raise ValidationError("Session webFiles.bindings.circularSources must be a list.", diagnostic=_SESSION_FIELDS_INVALID)
    uids: set[str] = set()
    for row in rows:
        if not isinstance(row, Mapping) or set(row) != CIRCULAR_SOURCE_FIELDS:
            raise ValidationError("A Circular File card holds uid, gb, gff and fasta.", diagnostic=_SESSION_FIELDS_INVALID)
        uid = row["uid"]
        if not isinstance(uid, str) or not uid or uid != uid.strip() or uid in uids:
            raise ValidationError("Each Circular File card needs a distinct uid.", diagnostic=_SESSION_FIELDS_INVALID)
        uids.add(uid)
        for key in CIRCULAR_SOURCE_FILE_SLOTS:
            if row[key] is not None:
                leaf(row[key])


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
    version = session.get("version")
    drawing = version == DRAWINGS_SESSION_MIN_VERSION
    if current and not (drawing or _is_web_draft_version(version, 41, 42, 44)):
        raise ValidationError(
            f"Web binding schema 2 requires session 41, 42, 44, or {MODE_SCOPED_SESSION_MIN_VERSION} and later.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )
    # A drawings-layout drawing (its one-drawing view names its mode in
    # ``ui.mode``) binds its own mode's inputs; Sessions 41 to the last
    # mode-scoped version bind both modes', the Circular one in c_gb, c_gff and c_fasta.
    mode = session["ui"]["mode"] if drawing else "circular"
    mode_slots = WEB_BINDING_MODE_SLOTS if drawing else {
        **WEB_BINDING_MODE_SLOTS,
        "circular": WEB_BINDING_MODE_SLOTS["circular"] - {"circularSources"} | MODE_SCOPED_CIRCULAR_INPUT_SLOTS,
    }
    input_slot = WEB_BINDING_INPUT_SLOT[mode] if drawing else "c_gb"
    other_slots = frozenset().union(*(slots for key, slots in mode_slots.items() if key != mode))
    if current and (input_slot not in bindings or (drawing and other_slots & set(bindings))):
        raise ValidationError(
            f"Web binding schema 2 requires {input_slot}"
            + (" and no other mode's input." if drawing else ".")
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
        else:
            leaf(binding)

    # Every slot that holds one file binding or a list of them.
    slots = (mode_slots["circular"] | WEB_BINDING_CARD_SLOTS) - {"c_conservation_blasts_source", "circularSources"}
    if current and set(bindings) - {"schema"} - WEB_BINDING_CARD_SLOTS - frozenset().union(*mode_slots.values()):
        raise ValidationError("Unknown Web binding inventory field.")
    for slot in slots:
        # A composite is the c_gb of a Session 41-50 only; the move to
        # drawings makes one File card of each of its components.
        value(bindings.get(slot), slot == "c_gb")
    if "circularSources" in bindings:
        _validate_circular_source_rows(bindings["circularSources"], leaf)
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
    ui = session.get("ui")
    _validate_display_placement_draft_config(
        config,
        mode=(session.get("renderRequest") or {}).get("mode") or (ui.get("mode") if isinstance(ui, Mapping) else None),
        scoped=True,
        current=session.get("version", 0) >= TYPED_DRAFT_SESSION_MIN_VERSION,
    )


def _validate_display_placement_draft_config(
    config: Mapping[str, Any],
    *,
    mode: object,
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
    expected_fields = set(RECORD_DISPLAY_DRAFT_FIELDS)
    if scoped:
        expected_fields.add("scope")
    if current:
        expected_fields |= TYPED_RECORD_DISPLAY_DRAFT_FIELDS
    keys = set()
    for row in drafts:
        if not isinstance(row, Mapping) or set(row) != expected_fields:
            raise ValidationError("Invalid record display draft fields.")
        if (scoped and row["scope"] not in {"circular", "linear"}) or any(
            not isinstance(row[name], str) or "\0" in row[name]
            for name in ("sourceUid", "selector", "recordId")
        ) or not row["sourceUid"] or not _SOURCE_RECORD_SELECTOR.fullmatch(row["selector"]):
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
                offset = intent.get("offsetBp") if isinstance(intent, Mapping) else None
                if (
                    not isinstance(intent, Mapping)
                    or set(intent) != RECORD_ANCHOR_INTENT_FIELDS
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
                    or abs(offset) > JS_MAX_SAFE_INTEGER
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


def _validate_mode_slices(session: Mapping[str, Any], version: int) -> None:
    """Validate each Session 46 mode slice against its own mode.

    A slice holds only the registry fields of its version (``slice_containers``); a missing field
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
    web_files = session.get("webFiles")
    bindings = web_files.get("bindings") if isinstance(web_files, Mapping) else None
    def bound_uids(slot: str) -> set[object]:
        cards = bindings.get(slot) if isinstance(bindings, Mapping) else None
        return {card.get("uid") for card in cards if isinstance(card, Mapping)} if isinstance(cards, list) else set()

    containers = slice_containers(version)
    # A Session 51 Circular drawing names its records by identity, a Session
    # 46-50 by source selector.
    identities = version >= DRAWINGS_SESSION_MIN_VERSION
    for mode, mode_slice in modes.items():
        _validate_mode_slice_fields(mode_slice, mode, containers)
        config = mode_slice.get("config", {})
        if "recordsOff" in config:
            _validate_records_off(
                config["recordsOff"], mode,
                bound_uids("linearSeqs" if mode == "linear" else "circularSources"),
                identities=identities,
            )
        _validate_display_placement_draft_config(config, mode=mode, scoped=False, current=True)
        _validate_feature_override_drafts(mode_slice.get("features", {}))
        overrides = config.get("unmanagedConfigOverrides")
        if overrides:
            validate_and_project_web_config_overrides(mode=mode, overrides=overrides)


def _is_record_identity(key: str, uids: set[object]) -> bool:
    """Whether ``key`` names a record of a bound card: the card uid (a file of
    one record) or ``<uid>:<k>``, ``k`` its 1-based position in the file."""

    base, _sep, position = key.rpartition(":")
    return key in uids or (base in uids and re.fullmatch(r"[1-9][0-9]*", position) is not None)


def _validate_records_off(records_off: object, mode: str, bound_uids: set[object], *, identities: bool) -> None:
    """Validate a slice's OFF records (record selection D-05), as the Web reader does.

    Each key appears once: a Linear card uid bound in ``webFiles``; a
    Circular source selector ``#N``, and in a Session 51 a Circular record
    identity of a bound File card (:func:`_is_record_identity`). A Session 51
    keeps a ``#N`` that the move of a Session 46-50 could not place
    (``circular_sources_of_mode_scoped``, OV-2869) until the Web reads the
    cards' records.
    """

    def names_record(key: str) -> bool:
        if mode == "linear":
            return key in bound_uids
        return (identities and _is_record_identity(key, bound_uids)) or _SOURCE_RECORD_SELECTOR.fullmatch(key) is not None

    if (
        not isinstance(records_off, list)
        or len(set(map(str, records_off))) != len(records_off)
        or not all(
            isinstance(key, str)
            and key
            and key == key.strip()
            and names_record(key)
            for key in records_off
        )
    ):
        raise ValidationError(
            f"Session modes.{mode}.config.recordsOff must list distinct record keys of its mode.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )


def mode_slice_of_version(value: Mapping[str, object], version: int, container: str = "") -> dict[str, object]:
    """``value`` (a slice, or its ``container``) with only the keys a slice of a Session of ``version`` may hold
    (``slice_containers``): what a writer of that version writes from a current slice."""

    containers = slice_containers(version)
    kept: dict[str, object] = {}
    for key, child in value.items():
        if key not in containers[container]:
            continue
        path = f"{container}.{key}" if container else key
        kept[key] = mode_slice_of_version(child, version, path) if path in containers and isinstance(child, Mapping) else child
    return kept


def _validate_mode_slice_fields(
    value: object, mode: str, containers: Mapping[str, frozenset[str]], container: str = ""
) -> None:
    """Reject a slice field that the registry of the Session's version (``containers``) does not name."""

    label = f"modes.{mode}" + (f".{container}" if container else "")
    if not isinstance(value, Mapping):
        raise ValidationError(f"Session {label} must be an object.", diagnostic=_SESSION_FIELDS_INVALID)
    unknown = sorted(str(field) for field in set(value) - containers[container])
    if unknown:
        raise ValidationError(
            f"Session {label} cannot contain {', '.join(unknown)}.", diagnostic=_SESSION_FIELDS_INVALID
        )
    for field, child in value.items():
        child_container = f"{container}.{field}" if container else field
        if child_container in containers:
            _validate_mode_slice_fields(child, mode, containers, child_container)


def _validate_current_retired_active_config_paths(
    session: Mapping[str, Any],
    version: int,
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
                f"Session version {version} cannot contain config.adv.{field}."
            )


def _validate_comparison_plan_edge(
    edge: object,
    allowed_fields: set[str],
    edge_ids: set[str],
) -> None:
    """Check one edge of a current Web comparison plan and add its ID to ``edge_ids``."""

    if not isinstance(edge, Mapping):
        raise ValidationError(
            "Each config.linearComparisonPlan edge must be an object."
        )
    unknown = sorted(str(field) for field in edge if field not in allowed_fields)
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
    # The Web checks an entry with `isWindowLosatFilenameEntry` (services/linear-comparisons.js).
    window_names = edge.get("windowLosatFilenames", {})
    if not isinstance(window_names, Mapping) or not all(
        isinstance(key, str) and key.strip() and isinstance(name, str) and name.strip()
        for key, name in window_names.items()
    ):
        raise ValidationError(
            "config.linearComparisonPlan edge windowLosatFilenames must map "
            "unit edge keys to file names.",
            diagnostic=_SESSION_FIELDS_INVALID,
        )


def _validate_current_comparison_authority(
    session: Mapping[str, Any],
    version: int,
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
            f"Session version {version} cannot contain retired blastSource state."
        )
    if "blastSource" in ui:
        raise ValidationError(
            f"Session version {version} cannot contain ui.blastSource."
        )
    if layout is not None and "comparisons" in layout:
        raise ValidationError(
            f"Session version {version} cannot contain "
            "config.linearRecordLayout.comparisons."
        )

    linear_sequences = bindings.get("linearSeqs")
    if isinstance(linear_sequences, list):
        for sequence in linear_sequences:
            if not isinstance(sequence, Mapping):
                continue
            if "blast" in sequence:
                raise ValidationError(
                    f"Session version {version} cannot contain per-record BLAST bindings."
                )
            if "losat_filename" in sequence:
                raise ValidationError(
                    f"Session version {version} cannot contain per-record LOSAT filenames."
                )
    linear_metadata = web_files.get("linearRecordMetadata")
    if isinstance(linear_metadata, list):
        for metadata in linear_metadata:
            if not isinstance(metadata, Mapping):
                continue
            if "losatFilename" in metadata or "losat_filename" in metadata:
                raise ValidationError(
                    f"Session version {version} cannot contain per-record LOSAT filenames."
                )
    if "linearCanonicalComparisons" in bindings:
        raise ValidationError(
            f"Session version {version} cannot bind comparison artifacts outside "
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
        # The names typed on the window pairs of a record pair, by unit edge key (OV-2742).
        *(("windowLosatFilenames",) if version >= WINDOW_LOSAT_FILENAMES_SESSION_VERSION else ()),
    }
    edge_ids: set[str] = set()
    for edge in edges:
        _validate_comparison_plan_edge(edge, allowed_edge_fields, edge_ids)

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


def _validate_alignment_reset_receipt(
    session: Mapping[str, object],
    sources: Mapping[str, SessionResourceSource],
) -> None:
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
    from .api.record_planning import SourceFileContents

    def source_identity(source: Mapping[str, Any]) -> dict[str, str]:
        identity = {}
        for key, value in source.items():
            if key == "kind":
                identity[key] = value
            else:
                if value not in fingerprints:
                    descriptors = session.get("resources")
                    if not isinstance(descriptors, Mapping) or not isinstance(descriptors.get(value), Mapping):
                        invalid()
                    try:
                        content = read_resource_bytes(SessionResources(descriptors, sources), value)
                    except (ValidationError, TypeError, KeyError) as exc:
                        raise ValidationError("Alignment reset receipt resource is invalid.") from exc
                    fingerprints[value] = SourceFileContents.fingerprint((content,))
                identity[key] = fingerprints[value]
        return identity

    from .session_request_codec import stored_record_windows

    def region_binding(record: Mapping[str, object]) -> object:
        # One window binds as the region it was before windows (request schema 10), so stored receipts still match.
        spans = [{"start": window["start"], "end": window["end"]} for window in stored_record_windows(record)]
        return spans[0] if len(spans) == 1 else spans or None

    records = [{
        "recordKey": record["recordKey"], "source": source_identity(record["source"]),
        "selector": record["selector"],
        "region": region_binding(record), "display": record["display"],
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


def _validate_other_mode_result(
    session: Mapping[str, object],
    version: int,
    sources: Mapping[str, SessionResourceSource],
) -> None:
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
    _validate_current_feature_payloads(view, version)
    _validate_current_feature_catalog_authority(view, version)
    if not other.get("results"):
        raise ValidationError("Session otherModeResult requires a Result.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    _validate_alignment_reset_receipt(view, sources)
    from .session_resources import canonical_resource_ids

    declared = session.get("resources")
    missing = sorted(canonical_resource_ids(request) - set(declared if isinstance(declared, Mapping) else ()))
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
            not isinstance(warning, Mapping) or set(warning) != ANNOTATION_WARNING_FIELDS
            or warning["code"] not in ANNOTATION_WARNING_CODES
            or not all(isinstance(warning[key], str) for key in ("setId", "annotationId", "recordId", "message", "resultName"))
            or not warning["setId"] or not warning["annotationId"]
            or not all(_index(warning[key]) for key in ("recordIndex", "missingCount", "resultIndex"))
            or (warning["code"] == "feature_selector_unmatched") != (warning["missingCount"] > 0)
            or not result_named(warning["resultIndex"], warning["resultName"])
        ):
            raise ValidationError("Annotation warnings do not match the successful Result metadata schema.", diagnostic=_OTHER_MODE_RESULT_INVALID)
    for warning in run_metadata.get("comparisonWarnings", []):
        if (
            not isinstance(warning, Mapping) or set(warning) != COMPARISON_WARNING_FIELDS
            or warning["code"] not in COMPARISON_WARNING_CODES
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
            not isinstance(notice, Mapping) or set(notice) != FEATURE_IDENTITY_NOTICE_FIELDS
            or not isinstance(notice["recordKey"], str) or not notice["recordKey"]
            or not isinstance(notice["biologicalFeatureId"], str) or not notice["biologicalFeatureId"]
            or notice["status"] not in FEATURE_IDENTITY_NOTICE_STATUSES
            or not isinstance(notice["kinds"], list) or not notice["kinds"]
            or not all(kind in FEATURE_IDENTITY_NOTICE_KINDS for kind in notice["kinds"])
            or not _index(notice["resultIndex"]) or notice["resultIndex"] >= len(results)
        ):
            raise ValidationError("runMetadata.featureIdentityNotices contains an invalid notice.", diagnostic=_OTHER_MODE_RESULT_INVALID)


def _validate_current_feature_payloads(
    session: Mapping[str, Any],
    version: int,
) -> None:
    """Reject the derived feature payloads a current Session never stores, and
    require its results array."""

    features = session.get("features")
    if isinstance(features, Mapping):
        duplicated = CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS & set(features)
        if duplicated:
            raise ValidationError(
                f"Session version {version} cannot contain derived feature payloads in "
                f"features: {', '.join(sorted(duplicated))}."
            )
    elif "features" in session:
        raise ValidationError("Session features must be an object when present.")

    orthogroup_state = session.get("orthogroupState")
    if isinstance(orthogroup_state, Mapping):
        if "groups" in orthogroup_state:
            raise ValidationError(
                f"Session version {version} cannot contain derived orthogroupState.groups."
            )
    elif "orthogroupState" in session:
        raise ValidationError(
            "Session orthogroupState must be an object when present."
        )
    if not isinstance(session.get("results"), list):
        raise ValidationError(f"Session version {version} requires a results array.")


def _validate_current_feature_catalog_authority(
    session: Mapping[str, Any],
    version: int,
) -> None:
    """Require the version-owned catalog of the Results, each a named plain SVG."""

    catalog_schema = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION.get(version, 3)
    results = session.get("results")
    assert isinstance(results, list)
    editor_state = session.get("editorState")
    if (
        not isinstance(editor_state, Mapping)
        or "featureCatalog" not in editor_state
    ):
        raise ValidationError(
            f"Session version {version} "
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
            f"Session version {version} Results require editorState.featureCatalog."
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
            catalog = promote_legacy_feature_catalog(catalog, schema=SEQUENCE_FEATURE_CATALOG_SCHEMA)
            catalog_schema = SEQUENCE_FEATURE_CATALOG_SCHEMA
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


def validate_current_session_artifacts(
    session: Mapping[str, Any],
    admission: SessionAdmission = CLI_SESSION_ADMISSION,
) -> None:
    """Validate current cache, manifest, and legacy artifact boundaries; the
    manifest resolution of protein caches only with ``admission.losat_manifest``."""

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
    if admission.losat_manifest:
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


def _feature_identity_key(scope: object, record_key: object, feature_id: object) -> str:
    """The draft key of one feature in one mode, or ``""`` for an invalid identity.

    The twin of ``featureIdentityKey`` in the Web ``feature-placement.js``.
    """

    if (
        scope in DRAFT_SCOPES
        and isinstance(record_key, str)
        and isinstance(feature_id, str)
        and record_key.strip(JS_WHITESPACE)
        and feature_id.strip(JS_WHITESPACE)
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


def _rendered_id_without_part(rendered_id: object) -> str:
    return _RENDERED_PART_SUFFIX.sub("", js_text(rendered_id), count=1)


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
    for item in mapping_list(items):
        record_keys_value = item.get("recordKeys")
        record_keys = (
            [js_text(record_key) for record_key in record_keys_value]
            if isinstance(record_keys_value, list)
            else []
        )
        for feature in mapping_list(item.get("features")):
            key = _feature_identity_key(
                mode, js_text(feature.get("recordKey")), js_text(feature.get("biologicalFeatureId"))
            )
            svg_id = _rendered_id_without_part(feature.get("svgId"))
            if key and svg_id:
                rendered_by_id.setdefault(svg_id, {})[key] = None
        for feature in mapping_list(item.get("biologicalFeatures")):
            record_key = js_text(feature.get("recordKey"))
            feature_id = js_text(feature.get("biologicalFeatureId"))
            key = _feature_identity_key(mode, record_key, feature_id)
            if key:
                biological.append(
                    (
                        key,
                        js_text(feature.get("stableFeatureId")) or feature_id,
                        safe_integer(feature.get("sourceFeatureIndex")),
                        record_keys.index(record_key) + 1 if record_key in record_keys else 0,
                    )
                )
    return rendered_by_id, biological


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
    the record key, ``expanded_record_key`` for each record of an ALL input
    with several records, and the source hash, ``~<source index>`` when the record
    has it twice.
    """
    from .api.requests import expanded_record_key

    fields = legacy if isinstance(legacy, Mapping) else {}
    linear = mode == "linear"
    records_value = fields.get("records")
    records = (
        [record if isinstance(record, Mapping) else {} for record in records_value]
        if isinstance(records_value, list)
        else []
    )

    def describe(feature: Mapping[str, Any]) -> tuple[int | None, int | None, int | None, str, str]:
        ordinal = _LINEAR_RECORD_SUFFIX.search(js_text(_first(feature, "svg_id", "svgId")))
        svg_ordinal = int(ordinal.group(1)) if ordinal else None
        file_index = nonnegative_integer(feature.get("fileIdx"))
        if linear:
            input_index = file_index if file_index is not None else (
                svg_ordinal - 1 if svg_ordinal else None
            )
        else:
            input_index = 0
        record_index = (
            0
            if linear and file_index is None
            else nonnegative_integer(_first(feature, "record_idx", "recordIndex"))
        )
        source_index = next(
            (
                index
                for index in (
                    nonnegative_integer(feature.get(field))
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
            js_text(_first(feature, "stable_feature_id", "stableFeatureId", "stable_svg_id")),
            js_text(next((value for value in drawn_hashes if value is not None), None)),
        )

    def described(value: object) -> list[tuple[int, int, int | None, str, str]]:
        entries: list[tuple[int, int, int | None, str, str]] = []
        for input_index, record_index, source_index, source_hash, drawn_hash in map(
            describe, mapping_list(value)
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
        record_key = js_text(record.get("recordKey")) if record is not None else ""
        if not record_key:
            return ""
        assert record is not None
        if (
            record.get("cardinality") == "all"
            and (linear or len(records) == 1)
            and record_counts.get(input_index, 0) > 1
        ):
            return expanded_record_key(record_key, record_index + 1)
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


# The last Session version whose `hash` color rules a Web writer filled with
# the hash of the drawn feature (release 0.13.0: Session 30, first-parent main:
# Session 44); `hash=` names the source-record hash after it (OV-401, release
# D-39).
DRAWN_HASH_RULE_SESSION_MAX_VERSION = 44

# A feature a saved diagram drew: its source hash, the hash it was drawn with,
# its source index, and its one-based record position in its Result.
DrawnFeature = tuple[str, str, int | None, int]


def _drawn_features_of_index(index: Sequence[_IndexedFeature]) -> list[DrawnFeature]:
    """The features of a feature index (:func:`_legacy_feature_index`,
    :func:`_catalog_drawn_features`) with the source hash of their identity key."""

    return [
        (_SOURCE_INDEX_SUFFIX.sub("", json.loads(key)[2], count=1), drawn_hash, source_index, ordinal)
        for key, drawn_hash, source_index, ordinal in index
    ]


def _catalog_drawn_features(catalog: object, mode: object) -> list[_IndexedFeature]:
    """Each feature a saved catalog (schema 3 or 4) drew, by the hash it was drawn with.

    A catalog row's rendered ID (``svgId``) carries the feature's hash in its
    drawn (cropped, reverse-complemented) record. The record position is the
    rendered ID's or the row's place in its item's records.
    """

    drawn: list[_IndexedFeature] = []
    items = catalog.get("items") if isinstance(catalog, Mapping) else None
    for item in mapping_list(items):
        record_keys_value = item.get("recordKeys")
        record_keys = (
            [js_text(record_key) for record_key in record_keys_value]
            if isinstance(record_keys_value, list)
            else []
        )
        for feature in mapping_list(item.get("features")):
            record_key = js_text(feature.get("recordKey"))
            key = _feature_identity_key(mode, record_key, js_text(feature.get("biologicalFeatureId")))
            drawn_hash, record_ordinal, source_index = _parse_rendered_id(feature.get("svgId"))
            if not key or not drawn_hash:
                continue
            if record_ordinal is None:
                record_ordinal = record_keys.index(record_key) + 1 if record_key in record_keys else 0
            drawn.append((key, drawn_hash, source_index, record_ordinal))
    return drawn


_ANCHORED_SELECTOR_VALUE = re.compile(r"\^(.*)\$", re.DOTALL)


def _hash_selector_parts(value: object) -> tuple[bool, str, int | None, int | None]:
    """Whether a ``hash`` value is written ``^<hash>$``, and the hash, record
    position, and source index of its rendered ID (release 0.13.0 wrote the
    rendered ID, record and instance suffixes included)."""

    text = js_text(value)
    anchored = _ANCHORED_SELECTOR_VALUE.fullmatch(text)
    return (anchored is not None, *_parse_rendered_id(anchored.group(1) if anchored else text))


def source_hash_selector_value(
    value: object, drawn: Sequence[DrawnFeature]
) -> tuple[str | None, int]:
    """Read a ``hash`` value of a Session 44 or older in the source frame.

    Such a Session matched ``hash=`` against the hash of the drawn feature, so
    on a cropped or reverse-complemented record its value is not the source
    hash ``hash=`` names now (OV-401, release D-39). The value, written as
    ``<hash>`` or ``^<hash>$`` and maybe with the rendered ID's record and
    instance suffixes, names the features of ``drawn`` drawn with that hash,
    ignoring case as the matchers do, in the record and with the source index
    its suffixes give. Returns the value naming their one source hash, in the
    same form without the suffixes, and the number of their source hashes.
    The value is ``None`` unless they have one source hash, and when that is
    the value's own hash (a record drawn untransformed) and the value has no
    Linear record suffix ``_record_<n>``, which no matcher reads (OV-416): such
    a value means what it meant and stays as saved. The
    Session draft (:func:`migrate_session_hash_rules`) and the CLI replay of a
    request (``gbdraw.api.session_compat``) read with it.
    """

    anchored, stable_id, record_ordinal, source_index = _hash_selector_parts(value)
    wanted = stable_id.lower()
    sources = {
        source_hash
        for source_hash, drawn_hash, feature_source_index, feature_ordinal in drawn
        if wanted
        and drawn_hash.lower() == wanted
        and (record_ordinal is None or feature_ordinal == record_ordinal)
        and (source_index is None or feature_source_index == source_index)
    }
    if len(sources) != 1:
        return None, len(sources)
    source = sources.pop()
    if source.lower() == wanted and not _LINEAR_RECORD_SUFFIX.search(js_text(value)):
        return None, 1
    return (f"^{source}$" if anchored else source), 1


def unmapped_hash_selector_warning(count: int, version: object) -> str:
    """The warning for ``hash`` values of an older Session that could not be renamed."""

    return (
        f"{count} hash= value(s) from Session version {version} could not be renamed "
        "to the hash of a feature in its source record: the saved diagram drew a "
        "cropped or reverse-complemented record, and the value names features of "
        "several source hashes there, or the sources were not read again. hash= now "
        "names a feature by its hash in the source record, so each may now match "
        "another feature or none."
    )


def _legacy_linear_request_records(files: object) -> list[dict[str, Any]]:
    """The request records of the Linear cards of a Session before 31.

    Such a Session (release 0.13.0) has no request: each card gives its record
    key, crop, and orientation, as the hash readers read them. For the source
    read of
    :func:`gbdraw.session_migration.read_request_less_source_features`, which
    the Web does not run, a card also gives its GenBank ``source`` (named by
    its slot path) and its record ``selector``.
    """

    linear_seqs = files.get("linearSeqs") if isinstance(files, Mapping) else None
    records: list[dict[str, Any]] = []
    for index, seq in enumerate(mapping_list(linear_seqs)):
        start = nonnegative_integer(seq.get("region_start"))
        end = nonnegative_integer(seq.get("region_end"))
        cropped = start is not None or end is not None
        reverse = bool(seq.get("region_reverse"))
        record_id = js_text(seq.get("region_record_id")).strip()
        records.append(
            {
                "recordKey": js_text(seq.get("uid")),
                "cardinality": "exactly_one" if cropped or record_id else "all",
                "region": (
                    {"start": start, "end": end, "reverseComplement": reverse} if cropped else None
                ),
                "presentation": {"reverseComplement": not cropped and reverse},
                "source": (
                    {"kind": "genbank", "resourceId": f"files.linearSeqs[{index}].gb"}
                    if isinstance(seq.get("gb"), Mapping)
                    else None
                ),
                "selector": {"kind": "recordId", "value": record_id} if record_id else None,
            }
        )
    return records


def _hash_frame_changed(record: Mapping[str, Any]) -> bool:
    """Whether a request record is drawn cropped or reverse-complemented, so
    the drawn hashes of its features differ from their source hashes (a
    rotation or a region that only selects a record changes no hash).
    """

    from .session_request_codec import stored_record_windows

    presentation = record.get("presentation")
    return any(
        nonnegative_integer(window["start"]) is not None
        or nonnegative_integer(window["end"]) is not None
        or js_truthy(window["reverseComplement"])
        for window in stored_record_windows(record)
    ) or (isinstance(presentation, Mapping) and js_truthy(presentation.get("reverseComplement")))


def has_hash_color_rules(rules: object) -> bool:
    """Whether specific color rules (``qual`` and ``val``) hold a ``hash`` rule.
    """

    return isinstance(rules, list) and any(
        isinstance(rule, Mapping) and js_text(rule.get("qual")).lower() == "hash" for rule in rules
    )


def _drawn_hashes_unknown(records: object, features: object) -> bool:
    """Whether ``features`` cannot name what a ``hash`` value of a request
    record in ``records`` drew: a record is drawn cropped or
    reverse-complemented and no feature carries its drawn hash (they are not
    the sources read again)."""

    return any(_hash_frame_changed(record) for record in mapping_list(records)) and not any(
        isinstance(feature.get("drawn_selector") or feature.get("drawnSelector"), Mapping)
        for feature in mapping_list(features)
    )


def request_less_hash_rules_need_drawn_hashes(session: Mapping[str, Any]) -> bool:
    """Whether the ``hash`` color rules of a Session without a request and
    without a saved catalog need its sources read again: one of its Linear
    cards is drawn cropped or reverse-complemented and its saved feature
    metadata holds no drawn hashes, so :func:`migrate_session_hash_rules`
    would count a value it cannot name as unmapped."""

    editor_state = session.get("editorState")
    config = session.get("config")
    features = session.get("features")
    return (
        session.get("renderRequest") is None
        and not (isinstance(editor_state, Mapping) and isinstance(editor_state.get("featureCatalog"), Mapping))
        and isinstance(config, Mapping)
        and has_hash_color_rules(config.get("rules"))
        and _drawn_hashes_unknown(
            _legacy_linear_request_records(session.get("files")),
            _legacy_draft_features({}, features if isinstance(features, Mapping) else {}),
        )
    )


@dataclass(frozen=True)
class HashRuleMigration:
    """The ``hash`` color rules of an older Session and how many stay unrenamed."""

    rules: Any
    unmapped_count: int


def migrate_session_hash_rules(
    rules: object,
    *,
    mode: object,
    catalog: object,
    legacy: object = None,
    records: object = None,
) -> HashRuleMigration:
    """Name the features of the ``hash`` color rules of a Session 44 or older by source hash.

    ``rules`` are the Session's specific color rules (``config.rules``,
    ``qual`` and ``val``), whose ``hash`` values release 0.13.0 and main wrote
    with the drawn feature's hash (its rendered ID in 0.13.0). Each such value
    that names features drawn with that hash (:func:`source_hash_selector_value`,
    through the saved ``catalog`` or else ``legacy``, as
    :func:`migrate_session_feature_edits` reads them) gets their source hash,
    so the rule matches the features it matched (S6 of the one feature-hash
    builder); a value whose features have its own hash as source hash (an
    untransformed record) is kept. A value naming features of several source
    hashes is kept and counted as unmapped; so is a value naming no indexed feature when a
    request record in ``records`` is drawn cropped or reverse-complemented but
    ``legacy`` holds no drawn hashes (its features are not the sources read
    again). Every other value is kept as saved. A list without a change is
    returned as saved.
    """

    if not has_hash_color_rules(rules):
        return HashRuleMigration(rules, 0)
    assert isinstance(rules, list)
    if js_truthy(catalog):
        drawn = _drawn_features_of_index(_catalog_drawn_features(catalog, mode))
        drawn_hashes_unknown = False
    else:
        drawn = _drawn_features_of_index(_legacy_feature_index(legacy, mode))
        drawn_hashes_unknown = _drawn_hashes_unknown(
            records, legacy.get("features") if isinstance(legacy, Mapping) else None
        )
    unmapped_count = 0
    changed = False
    migrated: list[Any] = []
    for rule in rules:
        value = rule.get("val") if isinstance(rule, Mapping) else None
        if not isinstance(rule, Mapping) or js_text(rule.get("qual")).lower() != "hash" or not js_text(value):
            migrated.append(rule)
            continue
        resolved, source_count = source_hash_selector_value(value, drawn)
        if resolved is None:
            unmapped_count += source_count > 1 or (source_count == 0 and drawn_hashes_unknown)
            migrated.append(rule)
        else:
            changed = True
            migrated.append({**rule, "val": resolved})
    return HashRuleMigration(migrated if changed else rules, unmapped_count)


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

    ``tests/fixtures/feature-edit-migration-vectors.json`` pins it. Without
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
            if js_truthy(catalog)
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
            visibility = _FEATURE_VISIBILITY_EDIT_VALUES.get(js_text(value).lower())
            if visibility is None:
                return False
            if row["featureVisibility"] is None:
                row["featureVisibility"] = visibility
            return True

        def assign_label_visibility(row: dict[str, Any], value: object) -> bool:
            visibility = js_text(value).lower()
            if visibility not in ("on", "off"):
                return False
            if row["labelVisibility"] is None:
                row["labelVisibility"] = visibility
            return True

        def assign_label_text(row: dict[str, Any], value: object) -> bool:
            label_text = re.sub(r"[\t\r\n\0]+", " ", js_string(value))
            if label_text.strip(JS_WHITESPACE):
                if row["labelText"] is None:
                    row["labelText"] = label_text
            elif row["labelVisibility"] != "on":
                # A blank text hid the label (its table row drew an empty label).
                row["labelVisibility"] = "off"
            return True

        def assign_label_source_text(row: dict[str, Any], value: object) -> bool:
            source_text = js_string(value)
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


_EXPANDED_RECORD_ORDINAL = re.compile(r"[1-9][0-9]*")
_SOURCE_INDEX_SUFFIX = re.compile(r"~[0-9]+\Z")


def _drawn_transformed(record: Mapping[str, Any]) -> bool:
    """Whether a request record is drawn cropped, reverse-complemented, or rotated.
    """

    from .session_request_codec import stored_record_windows

    presentation = record.get("presentation")
    display = record.get("display")
    return (
        # A schema <= 9 region counts by JavaScript truthiness, as the Web reader reads it.
        js_truthy(record.get("region"))
        or bool(stored_record_windows(record))
        or (isinstance(presentation, Mapping) and js_truthy(presentation.get("reverseComplement")))
        or (isinstance(display, Mapping) and safe_integer(display.get("startCoordinate")) is not None)
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

    R-7; ``tests/fixtures/annotation-target-migration-vectors.json`` pins it.
    Without a catalog nothing moves, as in the Web app, which reads no sources
    again for these targets.
    """

    record_ids: dict[str, set[str]] = {}
    features_by_hash: dict[str, dict[tuple[str, str], None]] = {}
    items = catalog.get("items") if isinstance(catalog, Mapping) else None
    for item in mapping_list(items):
        item_record_keys = item.get("recordKeys")
        for record_key in item_record_keys if isinstance(item_record_keys, list) else ():
            record_ids.setdefault(js_text(record_key), set())
        for feature in mapping_list(item.get("biologicalFeatures")):
            record_key = js_text(feature.get("recordKey"))
            feature_id = js_text(feature.get("biologicalFeatureId"))
            if not record_key or not feature_id:
                continue
            # A record listed only by a later item has no ID from this feature.
            if record_key in record_ids:
                record_id = feature.get("record_id")
                record_ids[record_key].add(
                    js_text(feature.get("recordId") if record_id is None else record_id)
                )
            source_hash = js_text(feature.get("stableFeatureId")) or _SOURCE_INDEX_SUFFIX.sub(
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
            index = js_string(selector.get("index"))
            position = int(index) if JS_ARRAY_INDEX.fullmatch(index) else len(catalog_record_keys)
            return catalog_record_keys[position] if position < len(catalog_record_keys) else ""
        if selector.get("kind") != "recordId":
            return ""
        # A record without catalog features has no known ID, so the binding is not certain.
        if any(len(record_ids[record_key]) != 1 for record_key in catalog_record_keys):
            return ""
        record_id = js_text(selector.get("value"))
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
        matches = list(features_by_hash.get(js_text(selector.get("value")), {}))
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
    return [value] if js_truthy(value) else []


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
        if not any(js_truthy(slot) for row in rows for slot in row):
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
    if js_truthy(normalized.get("collinear_block_plus_max")) and not js_truthy(
        normalized.get("collinear_block_plus")
    ):
        normalized["collinear_block_plus"] = normalized["collinear_block_plus_max"]
    for key, value in _DEFAULT_COMPARISON_COLORS.items():
        if not js_truthy(normalized.get(key)):
            normalized[key] = value
    return normalized


def _palette_colors(name: object) -> dict[str, Any]:
    """A palette's normalized colors, or none (``paletteColorsFromDefinitions``).

    The Web reads the palettes from ``gallery/palettes/palettes.json``, which
    ``tools/generate_palette_explorer_assets.py`` writes from the same
    ``gbdraw/data/color_palettes.toml``.
    """

    palette = js_text(name)
    colors = packaged_color_palettes().get(palette)
    if not colors:
        return {}
    return _normalize_palette_colors(dict(colors))


def _merges_override_colors(config: object) -> bool:
    colors = config.get("colors") if isinstance(config, Mapping) else None
    return (
        isinstance(config, Mapping)
        and js_truthy(config.get("colorsAreOverrides"))
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
    return _palette_colors(js_text(config.get("palette")) or "default")


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
        overrides = {key: js_text(value) if js_truthy(value) else "" for key, value in colors.items()}
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
        source_key, _ = _JSON_DECODER.raw_decode(binding.strip(JS_WHITESPACE))
    except ValueError:
        return None
    return _diagram_mode(source_key[0]) if isinstance(source_key, list) and source_key else None


def _annotation_sets_of_mode(sets: list[Any], mode: DiagramMode) -> list[Any]:
    result: list[Any] = []
    for annotation_set in sets:
        annotations = annotation_set.get("annotations") if isinstance(annotation_set, Mapping) else None
        if not isinstance(annotations, list):
            result.append(json_clone(annotation_set))
            continue
        kept = []
        for annotation in annotations:
            bound = _annotation_binding_mode(annotation)
            if bound is not None and bound != mode:
                continue
            annotation = json_clone(annotation)
            target = annotation.get("target") if isinstance(annotation, dict) else None
            if isinstance(target, dict) and target.get("kind") == "featureIdentity":
                target.pop("scope", None)
            kept.append(annotation)
        result.append({**json_clone(annotation_set), "annotations": kept})
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
        return {cast("DiagramMode", row.modes): json_clone(value)}
    if row.migrate == "result-mode":
        return {committed: _unscoped_rows(json_clone(value), committed)}
    if row.migrate == "layout":
        slots = value if isinstance(value, Mapping) else {}
        return {mode: json_clone(slots[mode]) for mode in DIAGRAM_MODES if mode in slots}
    split: dict[DiagramMode, Any] = {}
    for mode in DIAGRAM_MODES:
        if row.migrate == "show-if-source":
            # ``show_depth && hasSource``
            split[mode] = (widths[mode] > 0) if js_truthy(value) else value
        elif row.migrate == "depth" and isinstance(value, list):
            split[mode] = json_clone(value[: max(1, widths[mode])])
        elif row.migrate == "by-scope" and isinstance(value, list):
            split[mode] = [
                {field: item for field, item in json_clone(draft).items() if field != "scope"}
                for draft in value
                if isinstance(draft, Mapping) and draft.get("scope") == mode
            ]
        elif row.migrate == "by-side":
            split[mode] = _unscoped_rows(json_clone(value), mode)
        elif row.migrate == "by-leaf" and isinstance(value, Mapping):
            split[mode] = {
                path: json_clone(leaf)
                for path, leaf in value.items()
                if mode in unmanaged_config_override_modes(str(path))
            }
        elif row.migrate == "by-binding" and isinstance(value, list):
            split[mode] = _annotation_sets_of_mode(value, mode)
        else:
            split[mode] = json_clone(value)
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
    text = js_text(value)
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
    number = js_number(value)
    return js_integral(number) if math.isfinite(number) and number >= 0 else None


def _stroke_override_map(source: object) -> dict[str, Any]:
    """``normalizeStrokeOverrideMap(source, { requireOverride: true })`` of the Web."""

    normalized: dict[str, Any] = {}
    if not isinstance(source, Mapping):
        return normalized
    for key in js_key_order(source):
        value = source[key]
        name = js_text(key)
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
    for key in js_key_order(source):
        value = source[key]
        caption = js_text(key)
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


# How a Session 46 Legend edit is read in Session 47 (``legend_edits_by_row_key``):
# ``moved`` to the key of the row whose caption it was stored under; kept under
# its key as ``unmatched`` (no row shows or holds the key) or as ``unshown``
# (the key of a renamed row that no row shows); or dropped as a ``collision``
# (a later edit reached the same row, ``rowKey``).
LegendEditRekeyOutcome = Literal["moved", "unmatched", "unshown", "collision"]
LegendEditStore = Literal["colorOverrides", "strokeOverrides"]


class _LegendEditRekeyRequired(TypedDict):
    store: LegendEditStore
    fromKey: str
    toKey: str | None
    outcome: LegendEditRekeyOutcome


class LegendEditRekey(_LegendEditRekeyRequired, total=False):
    """One Session 46 Legend edit that did not stay as it was; ``toKey`` is
    ``None`` when it was dropped (a ``collision`` that keeps no stroke field). ``rowKey``: the row a ``collision`` reached;
    ``slice``: the mode slice, or ``"draft"`` (``session_legend_edits_by_row_key``)."""

    rowKey: str
    slice: str


def _session_46_legend_rows(legend: Mapping[str, Any]) -> tuple[dict[str, str], set[str]]:
    """The rows a Session 46 Legend's edits reach, by the caption each shows,
    and the keys of the renamed rows; the twin of ``session46LegendRows`` in
    the Web ``services/session-drawings-upgrade.js``.

    As Session 46 Load read the rows (an entry with a caption and a usable
    color, the first of each caption in its list) and as its compile found
    them (``draftLegendRows`` at 9839de8e): a listed row before a dormant one
    (the later dormant entry of a caption wins), and a deleted row as its
    Restore shows it. A dormant entry that is not renamed, or whose key a
    listed or deleted row holds, is no row.
    """

    def loaded(rows: str) -> list[tuple[str, str, bool]]:
        entries = legend.get(rows)
        captions: set[str] = set()
        found: list[tuple[str, str, bool]] = []
        for entry in entries if isinstance(entries, list) else ():
            if not isinstance(entry, Mapping):
                continue
            caption = js_text(entry.get("caption") or "")
            if not caption or not is_user_color(js_text(entry.get("color"))) or caption in captions:
                continue
            captions.add(caption)
            key = js_text(entry.get("originalCaption") or entry.get("caption") or "")
            found.append((caption, key, entry.get("dormant") is True))
        return found

    entries = loaded("entries")
    listed = [row for row in entries if not row[2]]
    deleted = loaded("deletedEntries")
    listed_keys = {key for _caption, key, _dormant in listed}
    deleted_keys = {key for _caption, key, _dormant in deleted}
    dormant = [
        row for row in [*loaded("dormantEntries"), *(row for row in entries if row[2])]
        if row[0] != row[1] and row[1] not in listed_keys and row[1] not in deleted_keys
    ]
    row_of = {caption: key for caption, key, _dormant in [*dormant, *listed]}
    for caption, key, _dormant in deleted:
        row_of.setdefault(caption, key)
    renamed_keys = {key for caption, key, _dormant in [*listed, *dormant, *deleted] if caption != key}
    return row_of, renamed_keys


def legend_edits_by_row_key(legend: Mapping[str, Any]) -> tuple[dict[str, Any], list[LegendEditRekey]]:
    """A Session 46 Legend (a slice's ``editorState.legend``) read as Session 47,
    and what the reading did to each edit that did not stay as it was.

    Session 46 stored the Legend color and stroke edits (``colorOverrides``,
    ``strokeOverrides``) under the caption a row shows, Session 47 under the
    row key (``originalCaption``, else the caption). Each row key gets the
    value Session 46 applied to that row (the Web compile at 9839de8e,
    ``app/candidate-render.js``):

    - the stores are read as Session 46 Load read them
      (:func:`_legend_color_overrides`, :func:`_stroke_override_map`);
    - an edit reaches the row whose caption it is stored under
      (:func:`_session_46_legend_rows`; ``moved`` when the key differs), else
      the row of Python's key of that name: it stays under its key, as
      ``unshown`` when that is the key of a renamed row that no row shows,
      else as ``unmatched``;
    - the color edits apply in their order, then the strokes in the order of
      the color edits' keys and then their own; a later edit paints over an
      earlier one on the same row (``collision``, ``rowKey``), a stroke field
      by field (``strokeColor``, ``strokeWidth``).

    Each outcome is ``{store, fromKey, toKey, outcome}`` (``toKey`` is
    ``None`` when dropped). Other fields are kept as they are.

    The twin of ``legendEditsByRowKey`` in the Web
    ``services/session-drawings-upgrade.js``;
    ``tests/fixtures/sessions/legend-row-key-vectors.json`` pins both.
    """

    row_of, renamed_keys = _session_46_legend_rows(legend)
    colors = _legend_color_overrides(legend.get("colorOverrides"))
    strokes = _stroke_override_map(legend.get("strokeOverrides"))
    changed = dict(legend)
    outcomes: list[LegendEditRekey] = []

    def key_by_row(store: LegendEditStore, edits: Mapping[str, Any], order: list[str], strokes_paint: bool) -> None:
        if not isinstance(legend.get(store), Mapping):
            return
        last_of = {row_of.get(key) or key: key for key in order}

        def fields_of(key: str) -> list[str]:
            return list(edits[key]) if strokes_paint else ["color"]

        field_last_of = {(row_of.get(key) or key, field): key for key in order for field in fields_of(key)}
        keyed: dict[str, Any] = {}
        for key in order:
            row_key = row_of.get(key) or key
            keyed[row_key] = {**keyed.get(row_key, {}), **edits[key]} if strokes_paint else edits[key]
            if last_of[row_key] != key:
                kept = any(field_last_of[(row_key, field)] == key for field in fields_of(key))
                outcomes.append({"store": store, "fromKey": key, "toKey": row_key if kept else None,
                                 "outcome": "collision", "rowKey": row_key})
            elif row_key != key:
                outcomes.append({"store": store, "fromKey": key, "toKey": row_key, "outcome": "moved"})
            elif key not in row_of:
                outcomes.append({"store": store, "fromKey": key, "toKey": key,
                                 "outcome": "unshown" if key in renamed_keys else "unmatched"})
        changed[store] = keyed

    color_keys = js_key_order(colors)
    key_by_row("colorOverrides", colors, color_keys, False)
    stroke_keys = js_key_order(strokes)
    key_by_row("strokeOverrides", strokes,
               [key for key in color_keys if key in strokes] + [key for key in stroke_keys if key not in colors], True)
    return changed, outcomes


def session_legend_edits_by_row_key(
    session: Mapping[str, Any],
) -> tuple[dict[str, Any], list[LegendEditRekey]]:
    """The fields of a Session older than 47 with each Legend's edits keyed by
    row key (:func:`legend_edits_by_row_key`): every mode slice's, and the flat
    draft's of a Session before 46; and the outcomes, each naming its
    ``slice`` (the mode, or ``"draft"``)."""

    outcomes: list[LegendEditRekey] = []

    def read(holder: Mapping[str, Any], slice_name: str) -> Mapping[str, Any]:
        editor = holder.get("editorState")
        legend = editor.get("legend") if isinstance(editor, Mapping) else None
        if not isinstance(editor, Mapping) or not isinstance(legend, Mapping):
            return holder
        keyed, found = legend_edits_by_row_key(legend)
        for outcome in found:
            outcome["slice"] = slice_name
        outcomes.extend(found)
        return {**holder, "editorState": {**editor, "legend": keyed}}

    changed = dict(read(session, "draft"))
    modes = changed.get("modes")
    if isinstance(modes, Mapping):
        changed["modes"] = {
            mode: read(value, str(mode)) if isinstance(value, Mapping) else value
            for mode, value in modes.items()
        }
    return changed, outcomes


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

    The Web Load and Gallery publication run this split; Python alone owns it.
    ``tests/fixtures/sessions/mode-split-vectors.json`` pins it.
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
        widths = {mode: safe_integer(depth_sources.get(mode)) or 0 for mode in DIAGRAM_MODES}
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
                split[active] = json_clone(flat)
            other: DiagramMode = "linear" if active == "circular" else "circular"
            saved = _mode_profile_values(profiles, other)
            if source_path in saved:
                split[other] = json_clone(saved[source_path])
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
        field: json_clone(losat[field])
        for field in LOSAT_EXECUTION_FIELDS
        if isinstance(losat, Mapping) and field in losat
    }
    # App-level settings leave the draft: LOSAT execution, the rich feature
    # popup (a missing value reads true), and Instant Preview.
    app_ui: dict[str, Any] = {"losatExecution": execution} if execution else {}
    if "rich_feature_popup" in adv:
        app_ui["richFeaturePopup"] = json_clone(adv["rich_feature_popup"])
    # Web Load applies ``ui``, then a boolean draft value, which wins.
    preview = config.get("paletteInstantPreviewEnabled")
    if isinstance(preview, bool):
        app_ui["paletteInstantPreviewEnabled"] = preview
    if app_ui:
        result["ui"] = {**(result.get("ui") or {}), **app_ui}
    # Session provenance, written only when the draft has it.
    if "cliOptions" in config:
        result["cliOptions"] = json_clone(config["cliOptions"])
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


def _current_option_value(value: object, fallback: str, supported: tuple[str, ...], label: str) -> str:
    """``requireCurrentValue`` in the Web ``services/current-option-values.js``."""

    normalized = ("" if value is None else js_string(value)).strip(JS_WHITESPACE).lower() or fallback
    if normalized not in supported:
        raise ValidationError(
            f"{label} must be one of: {', '.join(supported)}.", diagnostic=_SESSION_FIELDS_INVALID
        )
    return normalized


def _migrated_linear_track_layout(value: object) -> str:
    normalized = ("" if value is None else js_string(value)).strip(JS_WHITESPACE).lower() or "middle"
    normalized = {"spreadout": "above", "tuckin": "below"}.get(normalized, normalized)
    return _current_option_value(normalized, "middle", ("above", "middle", "below"), "Linear track layout")


def _migrated_linear_label_placement(value: object) -> str:
    normalized = ("" if value is None else js_string(value)).strip(JS_WHITESPACE).lower() or "auto"
    return _current_option_value(
        "above_feature" if normalized == "on_feature" else normalized,
        "auto",
        ("auto", "above_feature"),
        "Linear label placement",
    )


def _migrated_circular_multi_record_size_mode(value: object) -> str:
    normalized = ("" if value is None else js_string(value)).strip(JS_WHITESPACE).lower() or "auto"
    return _current_option_value(
        LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES.get(normalized, normalized),
        "auto",
        ("auto", "linear", "equal"),
        "Circular multi-record size mode",
    )


def migrate_persisted_web_option_values(config: object) -> object:
    """Move a Session 27-39 draft's retired option values to today's.

    Linear track layout ``spreadout``/``tuckin`` becomes ``above``/``below``,
    label placement ``on_feature`` becomes ``above_feature``, and multi-record
    size ``sqrt`` becomes ``auto``; any other value must be a current one. The
    twin of the Web's strict readers (``migratePersisted*`` in
    ``current-option-values.js``); Web Load reads these values with its reader
    table, which reads an unknown one as the default with a row (E-OD02).
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

    adv = config.get("adv")
    if not isinstance(adv, Mapping) or "circular_track_slots" not in adv:
        return
    if not isinstance(adv["circular_track_slots"], list):
        raise ValidationError("Custom Track Slots must be an array.", diagnostic=_SESSION_FIELDS_INVALID)
    schema = adv.get("circular_track_slots_schema_version")
    if js_integer(schema) and schema == _CURRENT_CIRCULAR_TRACK_SLOT_SCHEMA:
        return
    if js_integer(schema) and schema == _LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA:
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
        or js_truthy(adv.get("circular_track_slots_enabled"))
        or not (
            js_integer(adv.get("circular_track_slots_schema_version"))
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
    text = (js_string(value) if js_truthy(value) else "features").strip(JS_WHITESPACE).lower()
    renderer = _LINEAR_TRACK_RENDERER_ALIASES.get(text) or text
    return renderer if renderer in _LINEAR_TRACK_RENDERERS else "features"


def migrate_imported_linear_track_slots(config: object, source_version: object) -> object:
    """Bring a draft's Linear slots to schema 2.

    A Session up to 32 saved them with schema-1 meaning whatever it stored;
    later ones state their schema. A schema-1 features slot loses its height and
    spacing, and every slot gets its own ``params`` object. The twin of
    ``migrateImportedLinearTrackSlots`` and ``migrateLinearTrackSlotsToCurrentSchema``;
    a stack it refuses, Web Load reads with ``readImportedTrackSlots`` (E-OD02).
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
    if stored is not _ABSENT and not (js_integer(stored) and stored in supported):
        raise obsolete
    if js_integer(source_version) and cast(int, source_version) <= _LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION:
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
    shapes = adv.get("feature_shapes")
    shapes = shapes if isinstance(shapes, Mapping) else {}
    missing = {
        feature_type: shape
        for feature_type, shape in LEGACY_FEATURE_SHAPES.items()
        if (not isinstance(features, list) or feature_type in features) and feature_type not in shapes
    }
    if not missing:
        return config
    return {**config, "adv": {**adv, "feature_shapes": {**shapes, **missing}}}


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


_ANNOTATION_MARKS = ("line", "bracket", "band", "highlight")


def _clean_id(value: object, fallback: str) -> str:
    return (js_string(value) if js_truthy(value) else "").strip(JS_WHITESPACE) or fallback


def _draft_annotation_target(target: object) -> dict[str, Any]:
    source = json_clone(target) if isinstance(target, Mapping) else {}
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
                    else js_string(selector.get("key")),
                    "value": js_string(selector.get("value"))
                    if isinstance(selector, Mapping) and js_truthy(selector.get("value"))
                    else "",
                }
                for selector in selectors
            ]
            if isinstance(selectors, list)
            else [],
            "envelope": envelope,
            "circularPath": circular_path,
        }
    start = max(1.0, js_number(source.get("start")) or 1.0)
    end = max(1.0, js_number(source.get("end")) or 1.0)
    return {
        "kind": "coordinateSpan",
        "record": source.get("record"),
        "start": js_integral(start),
        "end": js_integral(end),
        "coordinateSpace": "local" if source.get("coordinateSpace") == "local" else "source",
        "wrapsOrigin": start > end,
        "outOfBounds": source.get("outOfBounds") if source.get("outOfBounds") in ("skip", "error") else "clip",
    }


def _draft_annotation_style(style: object) -> dict[str, Any]:
    from .session_request_codec import describe_draft_annotation_style

    return {**describe_draft_annotation_style(), **(json_clone(style) if isinstance(style, Mapping) else {})}


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
            item = json_clone(raw_item) if isinstance(raw_item, Mapping) else {}
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
                    "label": js_string(item.get("label")) if js_truthy(item.get("label")) else "",
                    "mark": item.get("mark") if item.get("mark") in _ANNOTATION_MARKS else "bracket",
                    "lane": None
                    if lane is None or lane == ""
                    else js_integral(max(0.0, js_number(lane) or 0.0)),
                    "style": None if item.get("style") is None else _draft_annotation_style(item.get("style")),
                    "legendLabel": None if item.get("legendLabel") is None else js_string(item.get("legendLabel")),
                    "metadata": json_clone(metadata) if isinstance(metadata, Mapping) else {},
                    **({"strand": item["strand"]} if item.get("strand") in {"+", "-"} else {}),
                }
            )
        result.append(
            {
                "id": set_id,
                "annotations": annotations,
                "defaultStyle": _draft_annotation_style(source.get("defaultStyle")),
                "legendLabel": None if legend_label is None else js_string(legend_label),
                **({"keepHighlightOverlay": True} if source.get("keepHighlightOverlay") is True else {}),
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
    unmapped_hash_selector_count: int = 0


def _legacy_draft_features(read_again: Mapping[str, Any], saved_features: Mapping[str, Any]) -> list[Any]:
    """The features that name a catalog-less Session's edits and ``hash``
    rules: the first non-empty of the sources read again and the saved
    feature metadata."""

    return next(
        (
            candidates
            for candidates in (
                read_again.get("extractedFeatures"),
                saved_features.get("biologicalFeatures"),
                saved_features.get("extractedFeatures"),
            )
            if isinstance(candidates, list) and candidates
        ),
        [],
    )


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
    again, see :func:`gbdraw.session_migration.read_legacy_source_features`
    and, for a Session before 31,
    :func:`gbdraw.session_migration.read_request_less_source_features`) and
    the saved feature metadata, then the ``hash`` color rules of a Session
    44 or older through the same catalog or sources
    (``migrate_session_hash_rules``), then the ``hash=`` annotation targets of
    a Session 40-44 (``migrate_session_annotation_targets``). A Session 40-44
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
    has_catalog = isinstance(catalog, Mapping)
    read_again = source_features if isinstance(source_features, Mapping) else {}
    saved_features = features if isinstance(features, Mapping) else {}
    legacy = None if has_catalog else {
        "records": _promoted_request_records(request),
        "features": _legacy_draft_features(read_again, saved_features),
        "biologicalFeatures": read_again.get("biologicalFeatures") or [],
    }
    if isinstance(features, Mapping):
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
    unmapped = 0
    if version <= DRAWN_HASH_RULE_SESSION_MAX_VERSION:
        # The mode of the Session's figure: its request's, else (a Session
        # before 31) its last Generate's invocation's, else ``ui.mode``. A
        # Session before 31 has no request: its Linear cards give the records
        # when that figure is Linear, as the old Web Load read them (OV-2056).
        generated_mode = session_mode(session)
        rule_records = (
            _promoted_request_records(request)
            if request
            else _legacy_linear_request_records(session.get("files")) if generated_mode == "linear" else None
        )
        hash_rules = migrate_session_hash_rules(
            config.get("rules") if isinstance(config, Mapping) else None,
            mode=generated_mode,
            catalog=catalog if has_catalog else None,
            legacy=None if legacy is None else {**legacy, "records": rule_records},
            records=rule_records,
        )
        if isinstance(config, Mapping) and hash_rules.rules is not config.get("rules"):
            config = {**config, "rules": hash_rules.rules}
            migrated["config"] = config
        unmapped = hash_rules.unmapped_count
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
    return SessionDraftMigration(migrated, dropped, narrowed, moved, unmapped)


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
    else:
        proof = _validated_protein_identity_manifest(source_manifest)
        if proof is None:
            raise ValidationError("Cannot write an invalid proteinIdentityManifest.")
        # The validator returns a proof unchanged. A read-only proof of plain
        # JSON data is written as it is; anything else as a JSON clone.
        current_manifest = (
            proof
            if proof is source_manifest and is_json_exact_protein_identity_manifest(proof)
            else json_clone(source_manifest)
        )

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
            clone = json_clone(entry)
            # As the Web writes it: a shown entry keeps its text name or gets
            # the default `losat_pair_<n>.tsv` (OV-1007); a hidden entry is
            # cache only and has no name (OV-1013).
            clone["filename"] = (
                ""
                if entry.get("display") is False
                else losat_pair_filename(entry, index + 1)
            )
            current_raw_entries.append(clone)
        elif classification == "protein-legacy":
            imported_legacy_entries.append(json_clone(entry))
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
            current_derived_entries.append(json_clone(entry))
        elif _is_derived_cache_entry(
            entry, schema=LEGACY_LOSAT_DERIVED_CACHE_SCHEMA
        ):
            imported_derived_evidence.append(json_clone(entry))
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


def _protein_manifest_holders(
    session: MutableMapping[str, object],
) -> list[MutableMapping[str, object]]:
    """Return the mappings that hold a manifest: each drawing of a 51, else the payload."""

    drawings = session.get("drawings")
    if isinstance(drawings, list):
        return [drawing for drawing in drawings if isinstance(drawing, MutableMapping)]
    return [session]


def seal_protein_identity_manifest(session: MutableMapping[str, object]) -> None:
    """Replace the payload's ``proteinIdentityManifest`` by its checked proof.

    The caller owns ``session``. Every later check of the payload trusts the
    proof. Only a proof of plain JSON data replaces the value: a manifest that
    fails, or one that holds other Python values, stays as it is, so
    ``validate_session`` checks and reports it as it does for a raw payload.
    A Session 51 holds a manifest in each drawing: each is sealed.
    """

    from .analysis.protein_colinearity import validate_protein_identity_manifest

    for holder in _protein_manifest_holders(session):
        manifest = holder.get("proteinIdentityManifest")
        if not isinstance(manifest, Mapping):
            continue
        try:
            proof = validate_protein_identity_manifest(manifest)
        except Exception:
            # validate_session reports the failure itself.
            continue
        if is_json_exact_protein_identity_manifest(proof):
            holder["proteinIdentityManifest"] = proof


def thaw_protein_identity_manifest(session: MutableMapping[str, object]) -> None:
    """Replace each proof in the payload by a plain, mutable copy.

    A Session 51 holds a manifest in each drawing: each is thawed.
    """

    from .analysis.protein_colinearity import ProteinIdentityManifest

    for holder in _protein_manifest_holders(session):
        manifest = holder.get("proteinIdentityManifest")
        if isinstance(manifest, ProteinIdentityManifest):
            holder["proteinIdentityManifest"] = manifest.to_dict()


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
        clone = json_clone(candidate)
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
        clone = json_clone(entry)
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


def file_resource_entry(path: str | Path) -> SessionResourceFile:
    """A local input file as a Web inventory file whose bytes stay in the file.

    The Session file's entry is written from the file, as it is on disk, when
    the Session is saved; nothing is read into memory here. Its ``type`` and
    ``lastModified`` are an input file's (:func:`written_file_metadata`).
    """

    file_path = Path(path)
    try:
        stat = file_path.stat()
    except OSError as exc:
        raise ValidationError(f"Could not read file for the Session: {file_path}") from exc
    name = file_path.name or "file"
    metadata = written_file_metadata()
    descriptor: SessionContainerResource = {
        "kind": "web-file",
        "name": name,
        "type": metadata["type"],
        "size": stat.st_size,
        "lastModified": metadata["lastModified"],
    }
    return {
        "resourceId": "",
        "descriptor": descriptor,
        "source": FileResourceSource(file_path),
        "name": name,
        "type": descriptor["type"],
        "lastModified": descriptor["lastModified"],
    }


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


def materialize_resource(resources: SessionResources, resource_id: str, *, temp_dir: Path) -> Path:
    """Write one Session resource's bytes, streamed from its source, to its file name in ``temp_dir``.

    A Session file's entry is checked by its size and CRC-32 as it is read;
    every source by the size its descriptor declares.
    """

    descriptor = resources.descriptors[resource_id]
    source = resources.source(resource_id)
    temp_dir.mkdir(parents=True, exist_ok=True)
    filename = safe_embedded_filename(descriptor.get("name"), fallback=f"{resource_id}.dat")
    output_path = _assert_under_directory(temp_dir / filename, temp_dir)
    written = 0
    try:
        with source.open_bytes() as stream, open(output_path, "wb") as output:
            while chunk := stream.read(1 << 20):
                output.write(chunk)
                written += len(chunk)
    except OSError as exc:
        raise ValidationError(f"Could not materialize embedded file for {resource_id}.") from exc
    if written != descriptor["size"]:
        raise ValidationError(
            f"Embedded file size mismatch for {resource_id}: expected {descriptor['size']}, decoded {written}.",
            diagnostic=SESSION_RESOURCES_INVALID,
        )
    return output_path


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


def _legacy_cli_protein_mode(config: Mapping[str, Any]) -> str:
    """The last ``--protein_blastp_mode`` of a legacy draft's ``cliOptions.rawArgs``, else ``none``."""

    cli_options = config.get("cliOptions")
    raw_args = cli_options.get("rawArgs") if isinstance(cli_options, Mapping) else None
    args = [str(arg) for arg in raw_args] if isinstance(raw_args, list) else []
    names = ("--protein_blastp_mode", "--protein-blastp-mode")
    values = [
        value
        for index, token in enumerate(args)
        if token.split("=", 1)[0] in names
        and (value := _option_value(args[index:index + 2], *names)) is not None
    ]
    return (values[-1] if values else "none").strip().lower()


class _LegacyRecordRow(NamedTuple):
    """What one ``linearSeqs`` row of a legacy draft held for its comparisons.

    ``blast`` is the upload compared with the next row, ``losat_filename`` the
    LOSAT result file name of that pair; a row that is not an object holds
    neither and has no uid.
    """

    uid: str
    blast: object
    losat_filename: str


class _LegacyComparisonSetting(NamedTuple):
    """The draft-wide comparison source of a legacy draft and the plan mode it implies."""

    default_source: str
    mode: Literal["none", "adjacent"]


@dataclass(frozen=True)
class _LegacyComparisonFile:
    """A ``linearComparisons`` entry of a legacy draft that binds a file, at its list position."""

    index: int
    comparison: dict[str, object]
    id: str
    query_uid: str
    subject_uid: str
    file: object


@dataclass(frozen=True)
class _MigratedEdgePayload:
    """The source and the result inputs of one migrated comparison edge."""

    source: str
    included: bool
    file: object
    file_active: bool
    losat_filename: str
    losat_filename_active: bool


class _MigratedComparisonEdge(TypedDict):
    id: str
    queryUid: str
    subjectUid: str
    included: bool
    fileActive: bool
    losatFilenameActive: bool
    source: str
    losatFilename: str


class _MigratedComparisonPlan(TypedDict):
    mode: Literal["none", "adjacent", "selected"]
    defaultSource: str
    edges: list[_MigratedComparisonEdge]


class _ComparisonFileBinding(TypedDict):
    """One ``files.linearComparisons`` entry of the current writer: the file of plan edge ``id``."""

    id: str
    file: object


def _take_legacy_record_rows(files: dict[str, object]) -> list[_LegacyRecordRow]:
    """Remove ``blast`` and ``losat_filename`` from each ``linearSeqs`` row of ``files`` and return what they held."""

    sequences_value = files.get("linearSeqs")
    rows: list[_LegacyRecordRow] = []
    sanitized_sequences: list[object] = []
    for sequence in sequences_value if isinstance(sequences_value, list) else []:
        if not isinstance(sequence, Mapping):
            sanitized_sequences.append(sequence)
            rows.append(_LegacyRecordRow(uid="", blast=None, losat_filename=""))
            continue
        row = dict(sequence)
        rows.append(
            _LegacyRecordRow(
                uid=str(row.get("uid") or ""),
                blast=row.pop("blast", None),
                losat_filename=str(row.pop("losat_filename", "") or ""),
            )
        )
        sanitized_sequences.append(row)
    files["linearSeqs"] = sanitized_sequences
    return rows


def _take_legacy_comparison_setting(config: dict[str, object]) -> _LegacyComparisonSetting:
    """Remove ``blastSource`` and ``adv.blastSource`` from ``config`` and return the setting they held."""

    adv_value = config.get("adv")
    adv = dict(adv_value) if isinstance(adv_value, Mapping) else {}
    raw_source = config.get("blastSource", adv.get("blastSource"))
    setting = _LegacyComparisonSetting(
        default_source=_legacy_comparison_source(raw_source),
        mode="none" if str(raw_source or "").strip().lower() == "none" else "adjacent",
    )
    config.pop("blastSource", None)
    adv.pop("blastSource", None)
    config["adv"] = adv
    return setting


def _take_layout_comparisons(config: dict[str, object]) -> list[dict[str, object]] | None:
    """Remove ``linearRecordLayout.comparisons`` from ``config``.

    Return its object items when the layout is enabled and the value is a
    list (the selected comparisons, authoritative even when empty); else None.
    """

    layout_value = config.get("linearRecordLayout")
    if not isinstance(layout_value, Mapping):
        return None
    layout = dict(layout_value)
    explicit_value = layout.pop("comparisons", None)
    config["linearRecordLayout"] = layout
    if not layout.get("enabled") or not isinstance(explicit_value, list):
        return None
    return [dict(item) for item in explicit_value if isinstance(item, Mapping)]


def _carry_existing_comparison_plan(
    existing_plan: Mapping[str, object],
    file_comparisons: list[dict[str, object]],
) -> tuple[dict[str, object], list[_ComparisonFileBinding]]:
    """An existing plan with each edge's ``file`` (or the file bound to its id) moved to a binding."""

    plan = json_clone(dict(existing_plan))
    edges_value = plan.get("edges")
    file_by_id: dict[str, object] = {
        str(item.get("id") or ""): item.get("file")
        for item in file_comparisons
        if str(item.get("id") or "") and item.get("file")
    }
    bindings: list[_ComparisonFileBinding] = []
    sanitized_edges = []
    for edge in edges_value if isinstance(edges_value, list) else []:
        if not isinstance(edge, Mapping):
            continue
        metadata = dict(edge)
        file_entry = metadata.pop("file", None) or file_by_id.get(str(metadata.get("id") or ""))
        sanitized_edges.append(metadata)
        if file_entry:
            bindings.append({"id": str(metadata.get("id") or ""), "file": file_entry})
    plan["edges"] = sanitized_edges
    return plan, bindings


def _cli_only_comparison_plan(config: Mapping[str, object]) -> _MigratedComparisonPlan | None:
    """The plan of a CLI-only draft: none (No comparison), or the adjacent LOSATP plan.

    A CLI-only draft gains no synthetic Web comparison. A legacy CLI writer
    selected LOSATP exactly when the last ``--protein_blastp_mode`` was not
    ``none``, and compared adjacent records.
    """

    if _legacy_cli_protein_mode(config) == "none":
        return None
    return {"mode": "adjacent", "defaultSource": "losat", "edges": []}


class _LegacyComparisonPlanBuilder:
    """The comparison plan of a legacy Web draft and the files its edges bind.

    The selected comparisons of an enabled layout come first, then each
    adjacent pair whose rows held a file or a LOSAT file name, then each bound
    file no edge took (kept, not included). Each edge takes a draft-unique id.
    """

    def __init__(
        self,
        rows: list[_LegacyRecordRow],
        file_comparisons: list[dict[str, object]],
        setting: _LegacyComparisonSetting,
    ) -> None:
        self.rows = rows
        self.setting = setting
        self.uid_index = {row.uid: index for index, row in enumerate(rows) if row.uid}
        self.comparison_files = [
            _LegacyComparisonFile(
                index=index,
                comparison=comparison,
                id=str(comparison.get("id") or ""),
                query_uid=str(comparison.get("queryUid") or ""),
                subject_uid=str(comparison.get("subjectUid") or ""),
                file=comparison.get("file"),
            )
            for index, comparison in enumerate(file_comparisons)
            if comparison.get("file")
        ]
        self.file_by_id: dict[str, _LegacyComparisonFile] = {}
        self.file_by_pair: dict[tuple[str, str], _LegacyComparisonFile] = {}
        for comparison_file in self.comparison_files:
            if comparison_file.id:
                self.file_by_id.setdefault(comparison_file.id, comparison_file)
            if comparison_file.query_uid and comparison_file.subject_uid:
                pair = (comparison_file.query_uid, comparison_file.subject_uid)
                self.file_by_pair.setdefault(pair, comparison_file)
        self.consumed_indexes: set[int] = set()
        self.used_ids: set[str] = set()
        self.edges: list[_MigratedComparisonEdge] = []
        self.bindings: list[_ComparisonFileBinding] = []

    def build(
        self, explicit: list[dict[str, object]] | None
    ) -> tuple[_MigratedComparisonPlan, list[_ComparisonFileBinding]]:
        mode: Literal["none", "adjacent", "selected"] = self.setting.mode
        used_payload_gaps: set[int] = set()
        if explicit is not None:
            mode = "selected" if explicit else "none"
            used_payload_gaps = self._add_selected_edges(explicit)
        self._add_adjacent_edges(mode, used_payload_gaps)
        self._add_retained_edges()
        plan: _MigratedComparisonPlan = {
            "mode": mode,
            "defaultSource": self.setting.default_source,
            "edges": self.edges,
        }
        return plan, self.bindings

    def _consume(self, comparison_file: _LegacyComparisonFile | None) -> object:
        if comparison_file is None:
            return None
        self.consumed_indexes.add(comparison_file.index)
        return comparison_file.file

    def _positional_file(self, index: int) -> object:
        """The file of the adjacent pair at gap ``index`` (a row with a next row); it consumes the binding of that pair."""

        row = self.rows[index]
        endpoint_binding = self.file_by_pair.get((row.uid, self.rows[index + 1].uid))
        if row.blast:
            # Canonical projection mirrored adjacent uploads into both legacy shapes.
            self._consume(endpoint_binding)
            return row.blast
        return self._consume(endpoint_binding)

    def _selected_endpoint(
        self, comparison: Mapping[str, object], uid_key: str, index_key: str
    ) -> tuple[str, int | None]:
        """A selected comparison's endpoint uid and row index: by uid, else by its stored index."""

        uid = str(comparison.get(uid_key) or "")
        index = self.uid_index.get(uid)
        index_value = comparison.get(index_key)
        if index is None and isinstance(index_value, int):
            index = index_value
            if 0 <= index < len(self.rows):
                uid = self.rows[index].uid
        return uid, index

    def _add_selected_edges(self, explicit: list[dict[str, object]]) -> set[int]:
        """Add the selected comparisons; return the adjacent gaps whose payload they took."""

        used_payload_gaps: set[int] = set()
        for index, comparison in enumerate(explicit):
            query_uid, query_index = self._selected_endpoint(comparison, "queryUid", "queryIndex")
            subject_uid, subject_index = self._selected_endpoint(
                comparison, "subjectUid", "subjectIndex"
            )
            # An adjacent pair is a gap between two rows: its file and LOSAT file name are those of the row.
            adjacent_gap = (
                query_index
                if query_index is not None
                and subject_index == query_index + 1
                and 0 <= query_index < len(self.rows) - 1
                else None
            )
            file_entry = (
                self._consume(self.file_by_id.get(str(comparison.get("id") or "")))
                or self._consume(self.file_by_pair.get((query_uid, subject_uid)))
                or (self._positional_file(adjacent_gap) if adjacent_gap is not None else None)
            )
            losat_filename = (
                self.rows[adjacent_gap].losat_filename if adjacent_gap is not None else ""
            )
            if adjacent_gap is not None and (file_entry or losat_filename):
                used_payload_gaps.add(adjacent_gap)
            self._add_edge(
                comparison.get("id"),
                (query_uid, subject_uid),
                ("selected", index),
                _MigratedEdgePayload(
                    source=_legacy_comparison_source(
                        comparison.get("source"), self.setting.default_source
                    ),
                    included=True,
                    file=file_entry,
                    file_active=bool(file_entry),
                    losat_filename=losat_filename,
                    losat_filename_active=bool(losat_filename),
                ),
            )
        return used_payload_gaps

    def _add_adjacent_edges(self, mode: str, used_payload_gaps: set[int]) -> None:
        """Add each adjacent pair whose rows held a file or a LOSAT file name and no selected edge took."""

        source = self.setting.default_source
        for index, row in enumerate(self.rows[:-1]):
            file_entry = self._positional_file(index)
            if (not file_entry and not row.losat_filename) or index in used_payload_gaps:
                continue
            payload_active = mode == "adjacent"
            file_active = payload_active and source == "upload" and bool(file_entry)
            filename_active = payload_active and source == "losat" and bool(row.losat_filename)
            self._add_edge(
                None,
                (row.uid, self.rows[index + 1].uid),
                ("adjacent", index),
                _MigratedEdgePayload(
                    source=source,
                    included=file_active or filename_active,
                    file=file_entry,
                    file_active=file_active,
                    losat_filename=row.losat_filename,
                    losat_filename_active=filename_active,
                ),
            )

    def _uid_at(self, index_value: object) -> str:
        if isinstance(index_value, int) and 0 <= index_value < len(self.rows):
            return self.rows[index_value].uid
        return ""

    def _add_retained_edges(self) -> None:
        """Keep each bound file no edge took as an edge that is not included."""

        for comparison_file in self.comparison_files:
            if comparison_file.index in self.consumed_indexes:
                continue
            comparison = comparison_file.comparison
            self._add_edge(
                comparison_file.id,
                (
                    comparison_file.query_uid or self._uid_at(comparison.get("queryIndex")),
                    comparison_file.subject_uid or self._uid_at(comparison.get("subjectIndex")),
                ),
                ("retained", comparison_file.index),
                _MigratedEdgePayload(
                    source=_legacy_comparison_source(
                        comparison.get("source"), self.setting.default_source
                    ),
                    included=False,
                    file=comparison_file.file,
                    file_active=False,
                    losat_filename="",
                    losat_filename_active=False,
                ),
            )

    def _add_edge(
        self,
        edge_id: object,
        pair: tuple[str, str],
        fallback: tuple[str, int],
        payload: _MigratedEdgePayload,
    ) -> None:
        """Add an edge between ``pair`` under ``edge_id``, or the stable id of ``fallback``, made unique."""

        query_uid, subject_uid = pair
        if not query_uid or not subject_uid:
            return
        base_id = str(edge_id or "").strip() or _stable_migrated_comparison_id(
            *fallback, query_uid, subject_uid
        )
        unique_id = base_id
        suffix = 2
        while unique_id in self.used_ids:
            unique_id = f"{base_id}-{suffix}"
            suffix += 1
        self.used_ids.add(unique_id)
        self.edges.append(
            {
                "id": unique_id,
                "queryUid": query_uid,
                "subjectUid": subject_uid,
                "included": bool(payload.included),
                "fileActive": bool(payload.file_active),
                "losatFilenameActive": bool(payload.losat_filename_active),
                "source": payload.source,
                "losatFilename": payload.losat_filename,
            }
        )
        if payload.file:
            self.bindings.append({"id": unique_id, "file": payload.file})


def migrate_legacy_linear_comparison_draft_for_current_writer(
    config: Mapping[str, Any],
    files: Mapping[str, Any],
    *,
    force_web_draft: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Project a supported pre-v40 comparison draft into final v40 ownership.

    The draft of a Session 27-39 held its comparisons in ``linearSeqs[].blast`` and
    ``losat_filename``, ``blastSource`` (or ``adv.blastSource``),
    ``linearRecordLayout.comparisons`` and ``linearComparisons``; the result
    holds them as ``config.linearComparisonPlan`` and the files its edges bind
    as ``files.linearComparisons``. An existing plan is carried with its edge
    files moved to bindings. A CLI-only draft (not ``force_web_draft``) gains
    no Web comparison.
    """

    migrated_config = json_clone(dict(config))
    migrated_files = json_clone(dict(files))
    rows = _take_legacy_record_rows(migrated_files)
    comparisons_value = migrated_files.get("linearComparisons")
    file_comparisons = (
        [dict(item) for item in comparisons_value if isinstance(item, Mapping)]
        if isinstance(comparisons_value, list)
        else []
    )
    setting = _take_legacy_comparison_setting(migrated_config)
    explicit = _take_layout_comparisons(migrated_config)
    existing_plan = migrated_config.get("linearComparisonPlan")
    plan: Mapping[str, object] | None
    bindings: list[_ComparisonFileBinding]
    if isinstance(existing_plan, Mapping):
        plan, bindings = _carry_existing_comparison_plan(existing_plan, file_comparisons)
    elif not force_web_draft:
        migrated_config.pop("linearRecordLayout", None)
        migrated_config.pop("linearComparisonPlan", None)
        plan, bindings = _cli_only_comparison_plan(migrated_config), []
    else:
        plan, bindings = _LegacyComparisonPlanBuilder(rows, file_comparisons, setting).build(explicit)
    if plan is not None:
        migrated_config["linearComparisonPlan"] = plan
    migrated_files["linearComparisons"] = bindings
    return migrated_config, migrated_files


class _SessionResourceDescriptorRequired(TypedDict):
    size: int
    # Base64 text, or the encoded depth payload under ``DEPTH_FILE_ENCODING``.
    data: str | Mapping[str, object]


class SessionResourceDescriptor(_SessionResourceDescriptorRequired, total=False):
    """An embedded file: a resource of a JSON Session 27-50, or an inline Web file.

    Its bytes are in the JSON. They leave it once, where the document is read
    (:func:`split_resource_bytes`); a Session document in memory and in a
    Session file holds :class:`SessionContainerResource` descriptors only.
    """

    kind: str
    name: str
    type: str
    lastModified: int
    encoding: str
    checksum: str | None


class WrittenFileMetadata(TypedDict):
    """What a file descriptor gbdraw writes says of its file besides its name and bytes."""

    type: str
    lastModified: int


# The media type of each resource kind gbdraw makes itself: canonical TSV
# tables and analysis results. The Web writers read both constants from the
# generated contracts (tools/generate_web_contracts.py).
WRITTEN_RESOURCE_MEDIA_TYPES: Mapping[str, str] = {
    "canonical-tsv": "text/tab-separated-values",
    "collinearity-result": "application/json",
    "orthogroup-result": "application/json",
}
INPUT_FILE_MEDIA_TYPE = "application/octet-stream"


def written_file_metadata(kind: str | None = None) -> WrittenFileMetadata:
    """The ``type`` and ``lastModified`` a Python Session writer stores for a file of resource ``kind``.

    The one owner of that metadata (OV-2964): a kind gbdraw makes itself has
    its media type, any input file (``kind`` None or an input kind) has
    ``application/octet-stream``, and ``lastModified`` is 0, so a Session
    depends only on its inputs' names and bytes and its options, never on the
    machine's media type guess or an input's modification time. A Web upload
    keeps the browser's metadata; its writer is the Web's.
    """

    return {"type": WRITTEN_RESOURCE_MEDIA_TYPES.get(kind or "", INPUT_FILE_MEDIA_TYPE), "lastModified": 0}


def is_session_resource_descriptor(value: Mapping[str, object]) -> TypeGuard[SessionResourceDescriptor]:
    """Whether ``value`` has the required fields of an embedded file.

    Its bytes are checked against them where they are read
    (:func:`_embedded_resource_bytes`).
    """

    return isinstance(value.get("size"), int) and isinstance(value.get("data"), (str, Mapping))


def missing_resource_bytes_message(resource_id: str) -> str:
    """The refusal of a resource whose bytes were not given: what a mapping from
    ``to_dict()`` needs, since the document holds no input bytes."""

    return (
        f"Canonical resource {resource_id!r} has no bytes: pass resources=document.resources with a "
        "document's to_dict(), or save and load the .gbdraw file."
    )


def retired_json_session_message(source: object) -> str:
    """The refusal of a JSON Session of the current version, whose JSON form never reached main;
    ``source`` is the file's path, or ``<mapping>`` for a mapping."""

    return f"A Session version {CURRENT_SESSION_VERSION} file must be a .gbdraw file, not JSON: {source}"


@dataclass(frozen=True)
class SessionResources:
    """The resources of one Session: the document's descriptors and their bytes' sources, by ID.

    The descriptors are JSON (``SessionContainerResource``); the sources stay
    beside the document, never in it. :func:`read_resource_bytes` is the one
    reader of a resource's bytes.
    """

    descriptors: Mapping[str, SessionContainerResource]
    sources: Mapping[str, SessionResourceSource]

    def source(self, resource_id: str) -> SessionResourceSource:
        source = self.sources.get(resource_id)
        if source is None or resource_id not in self.descriptors:
            raise ValidationError(missing_resource_bytes_message(resource_id), diagnostic=SESSION_RESOURCES_INVALID)
        return source

    def read(self, resource_id: str) -> bytes:
        """The bytes of resource ``resource_id`` (:func:`read_resource_bytes`)."""

        return read_resource_bytes(self, resource_id)


class SessionResourceFile(TypedDict):
    """A Web inventory file whose bytes stay at their source until the Session file is written.

    A binding of a read document (:func:`_project_web_file_binding`) or a CLI
    input file (:func:`file_resource_entry`); :func:`_attach_current_web_file_bindings`
    binds it by its descriptor and source. ``resourceId`` is the preferred ID
    (empty for a new file).
    """

    resourceId: str
    descriptor: SessionContainerResource
    source: SessionResourceSource
    name: str
    type: str
    lastModified: int | float


class BoundWebFile(TypedDict):
    """A Web inventory file bound to a resource of the document it is written into
    (``webFiles.bindings``): the resource ID and the file's metadata as the inventory gave it."""

    resourceId: str
    name: str
    type: str
    # As the inventory holds it; the Session reader checks it.
    lastModified: object


def read_resource_bytes(resources: SessionResources, resource_id: str) -> bytes:
    """The bytes of one Session resource, found by its ID and checked against its size."""

    return read_source_bytes(resources.source(resource_id), resources.descriptors[resource_id]["size"])


def read_source_bytes(source: SessionResourceSource, size: int) -> bytes:
    """The bytes of ``source``, refused unless there are exactly ``size`` of them.

    A Session file's entry also checks its CRC-32 (``session_container``).
    """

    try:
        with source.open_bytes() as stream:
            data = stream.read(size + 1)
    except OSError as exc:
        raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID) from exc
    if len(data) != size:
        raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID)
    return data


def split_resource_bytes(
    session: Mapping[str, object],
    sources: Mapping[str, SessionResourceSource] | None = None,
) -> tuple[dict[str, object], dict[str, SessionResourceSource]]:
    """``session`` with JSON descriptors only, and the source of each resource's bytes.

    The one place a document's resources take their in-memory form. A
    descriptor of a JSON Session 27-50 (``encoding`` and ``data``) gives its
    bytes, checked against its size and any ``checksum``, and keeps its file
    metadata; any other descriptor must be a ``SessionContainerResource``
    whose bytes ``sources`` holds (a Session file's entries, or the sources
    of the document it was built from). The descriptors are not checked
    here: the Session reader checks each as a ``SessionContainerResource``
    afterwards, so they stay unchecked JSON values.
    """

    resources = session.get("resources")
    if not isinstance(resources, Mapping):
        return dict(session), {}
    descriptors: dict[str, object] = {}
    split: dict[str, SessionResourceSource] = {}
    for resource_id, entry in resources.items():
        if isinstance(entry, Mapping) and ("encoding" in entry or "data" in entry):
            try:
                content = _embedded_resource_bytes(entry)
            except ValidationError as exc:
                raise ValidationError(f"Canonical resource {resource_id!r} is invalid: {exc}", diagnostic=SESSION_RESOURCES_INVALID) from exc
            modified = entry.get("lastModified", 0)
            descriptors[resource_id] = {
                "kind": entry.get("kind"),
                "name": entry.get("name"),
                "type": entry.get("type"),
                "size": entry.get("size"),
                # A JSON Session may store a fraction of a millisecond; a time
                # it cannot hold is read as 0, as the Web recomputes it.
                "lastModified": (int(modified) if isinstance(modified, (int, float))
                                 and not isinstance(modified, bool) and math.isfinite(modified)
                                 and modified >= 0 else 0),
            }
            split[resource_id] = BytesResourceSource(content)
            continue
        descriptors[resource_id] = entry
        source = (sources or {}).get(resource_id)
        if source is None:
            raise ValidationError(missing_resource_bytes_message(resource_id), diagnostic=SESSION_RESOURCES_INVALID)
        split[resource_id] = source
    return {**session, "resources": descriptors}, split


def _decodable_resource_sources(resources: Mapping[str, object]) -> dict[str, SessionResourceSource]:
    """The sources of the embedded resources whose bytes decode, unchecked against their descriptors."""

    sources: dict[str, SessionResourceSource] = {}
    for resource_id, entry in resources.items():
        content = _embedded_entry_bytes(entry) if isinstance(entry, Mapping) else None
        if content is not None:
            sources[resource_id] = BytesResourceSource(content)
    return sources


def _embedded_entry_bytes(entry: Mapping[str, object]) -> bytes | None:
    encoding = entry.get("encoding")
    data = entry.get("data")
    if encoding in (None, "base64") and isinstance(data, str):
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


def _embedded_resource_bytes(entry: Mapping[str, object]) -> bytes:
    """The bytes of a resource descriptor, checked against its size and checksum.

    The one check of a declared ``checksum``: the SHA-256 hex digest of the
    bytes, bare or as ``sha256:<hex>``, in any case. Null or empty declares
    none, as in the Web reader (``services/session-resource-backing.js``).
    """

    data = _embedded_entry_bytes(entry)
    if data is None or len(data) != entry.get("size"):
        raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID)
    checksum = entry.get("checksum")
    if checksum not in (None, "") and (
        not isinstance(checksum, str)
        or checksum.strip().lower().removeprefix("sha256:") != hashlib.sha256(data).hexdigest()
    ):
        raise ValidationError("Embedded resource checksum does not match.")
    return data


def _project_web_file_binding(
    resources: SessionResources,
    binding: Any,
    *,
    schema: int | None,
) -> Any:
    """Transport an admitted binding without choosing destination identities.

    Schema None is the existing direct-source metadata-default context.
    Each leaf carries its resource's descriptor and bytes' source, unmodified,
    until assembly (:func:`_attach_current_web_file_bindings`).
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
    resource = resources.descriptors.get(resource_id)
    if not isinstance(resource, Mapping):
        raise ValidationError(f"Web file binding references a missing resource: {resource_id}.")
    leaf: SessionResourceFile = {
        "resourceId": resource_id,
        "descriptor": resource,
        "source": resources.source(resource_id),
        "name": binding["name"] if schema == 2 else str(binding.get("name") or resource.get("name") or "file"),
        "type": binding["type"] if schema == 2 else str(binding.get("type") or resource.get("type") or ""),
        "lastModified": (
            binding["lastModified"] if schema == 2
            else int(binding.get("lastModified") or resource.get("lastModified") or 0)
        ),
    }
    return leaf


def _web_card_binding(
    sequence: Mapping[str, object],
    add_file: Callable[[object], object],
    add_value: Callable[[object], object],
) -> dict[str, object]:
    """One Linear card's bindings, the card shape the Web writes (``linearBindings`` in services/session-resources.js)."""

    return {
        "uid": str(sequence.get("uid") or ""),
        "gb": add_file(sequence.get("gb")),
        "gff": add_file(sequence.get("gff")),
        "fasta": add_file(sequence.get("fasta")),
        "depth": add_value(sequence.get("depth")),
        "losat_gencode": sequence.get("losat_gencode", 1),
        "definition": str(sequence.get("definition") or ""),
        "record_subtitle": str(sequence.get("record_subtitle") or ""),
        "file_definition": str(sequence.get("file_definition") or ""),
        "file_subtitle": str(sequence.get("file_subtitle") or ""),
        "inferred_definition": str(sequence.get("inferred_definition") or ""),
        "region_record_id": str(sequence.get("region_record_id") or ""),
        "region_start": sequence.get("region_start"),
        "region_end": sequence.get("region_end"),
        "region_reverse": bool(sequence.get("region_reverse")),
        # The card's window list (several windows, or one across the origin), as the Web writes it.
        "windows": [dict(window) for window in mapping_list(sequence.get("windows"))],
    }


def _circular_input_bindings(
    files: Mapping[str, object],
    add_file: Callable[[object], Mapping[str, object] | None],
    add_value: Callable[..., object],
) -> dict[str, object]:
    """The Circular input of a Web inventory as bindings: File cards (Session 51,
    and the CLI writer's own mode-scoped payload), or the one input of a Session
    41-50 source, which the move to drawings makes cards of."""

    circular_cards = files.get("circularSources")
    if circular_cards is not None and not isinstance(circular_cards, list):
        raise ValidationError(
            "Web inventory circularSources must be a list.", diagnostic=_SESSION_FIELDS_INVALID
        )
    circular_input: dict[str, object] = (
        {
            "circularSources": [
                CircularSourceCard(
                    uid=str(card.get("uid") or ""),
                    gb=add_file(card.get("gb")),
                    gff=add_file(card.get("gff")),
                    fasta=add_file(card.get("fasta")),
                )
                for card in circular_cards
                if isinstance(card, Mapping)
            ]
        }
        if circular_cards is not None
        else {
            "c_gb": add_value(files.get("c_gb"), composite_allowed=True),
            "c_gff": add_file(files.get("c_gff")),
            "c_fasta": add_file(files.get("c_fasta")),
        }
    )
    if circular_cards is not None and MODE_SCOPED_CIRCULAR_INPUT_SLOTS & set(files):
        raise ValidationError(
            "A Web inventory holds File cards or c_gb, c_gff and c_fasta, not both.", diagnostic=_SESSION_FIELDS_INVALID
        )
    return circular_input


def _attach_current_web_file_bindings(
    payload: dict[str, object],
    files: Mapping[str, object],
    sources: dict[str, SessionResourceSource],
) -> None:
    """Bind the Web files ``files`` to resources of ``payload``.

    ``sources`` holds the bytes' source of each resource of ``payload``; the
    sources of the resources the bindings add go to it.
    """

    from .session_resources import SessionResourceTable

    resources_value = payload.get("resources")
    if not isinstance(resources_value, dict):
        raise ValidationError("Current session resources must be an object.")
    resources = resources_value
    for resource_id in resources:
        if not isinstance(resource_id, str) or resource_id != resource_id.strip() or not resource_id:
            raise ValidationError("Canonical resource IDs must be unique non-empty strings.")
    table = SessionResourceTable(SessionResources(resources, sources))

    def metadata_of(resource_id: str, metadata: Mapping[str, object]) -> BoundWebFile:
        return {
            "resourceId": resource_id,
            "name": str(metadata.get("name", "file")),
            "type": str(metadata.get("type") or ""),
            "lastModified": metadata.get("lastModified", 0),
        }

    def allocate_file(entry: Mapping[str, object], preferred_id: str, metadata: Mapping[str, object]) -> BoundWebFile:
        if not is_session_resource_descriptor(entry):
            raise ValidationError("Invalid embedded resource bytes or byte size.", diagnostic=SESSION_RESOURCES_INVALID)
        return metadata_of(table.bind(entry, preferred_id=preferred_id), metadata)

    def add_file(value: object) -> BoundWebFile | None:
        if value is None:
            return None
        if not isinstance(value, Mapping) or "components" in value or value.get("kind") == "composite":
            raise ValidationError("Expected an ordinary Web inventory file.")
        if "descriptor" in value:
            if (set(value) != {"descriptor", "source", "resourceId", "name", "type", "lastModified"}
                    or not isinstance(value["descriptor"], Mapping)):
                raise ValidationError("Invalid resource-backed Web inventory file.")
            # The descriptor of a checked document (:func:`_project_web_file_binding`).
            descriptor = cast(SessionContainerResource, value["descriptor"])
            return metadata_of(
                table.bind_resource(descriptor, value["source"], preferred_id=value["resourceId"]),
                value,
            )
        if "resourceId" in value:
            raise ValidationError("Web inventory source references require a descriptor.")
        return allocate_file(value, "", value)

    def add_value(value: object, *, composite_allowed: bool = False) -> object:
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
    sequence_bindings = [
        _web_card_binding(sequence, add_file, add_value)
        for sequence in linear_sequences
        if isinstance(sequence, Mapping)
    ]
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

    circular_input = _circular_input_bindings(files, add_file, add_value)
    bindings = {
        "schema": 2,
        **circular_input,
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
    bound = table.resources()
    resources.update(bound.descriptors)
    sources.update(bound.sources)
    payload["webFiles"] = _rebound_web_files(payload.get("webFiles"), bindings, files, resources)


def _rebound_web_files(
    web_files_value: object,
    bindings: dict[str, object],
    files: Mapping[str, object],
    resources: Mapping[str, object],
) -> dict[str, object]:
    """The ``webFiles`` of a payload whose Web files ``files`` were bound to its
    resources as ``bindings``: its source references and original names follow
    the resource IDs the binding gave."""

    web_files = (
        json_clone(dict(web_files_value))
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

    def reference_value(value: object, source: object) -> object:
        if isinstance(value, list):
            parts = source if isinstance(source, list) else []
            return [reference_value(item, parts[index]) for index, item in enumerate(value)]
        if value is None:
            return None
        # A binding that add_file made (BoundWebFile).
        assert isinstance(value, Mapping)
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
    return web_files


def drawings_document_of_mode_scoped(
    data: Mapping[str, object],
    *,
    stored_result_modes: frozenset[str] = frozenset(),
) -> tuple[GbdrawSession51, ModeScopedSourceFacts]:
    """A Session 46-50 as the current Session: its move to drawings
    (:func:`gbdraw.session_drawings.drawings_document_from_mode_scoped`, which only restructures the document, as its
    Web twin does), then each catalog promoted from schema 5, which stored the sequences, to the current schema, which
    holds none (D5)."""

    from .session_drawings import drawings_document_from_mode_scoped
    from .web_support.feature_catalog import promote_legacy_feature_catalog

    document, facts = drawings_document_from_mode_scoped(data, stored_result_modes=stored_result_modes)
    drawings = document.get("drawings")
    for drawing in drawings if isinstance(drawings, list) else ():
        editor = drawing.get("editorState") if isinstance(drawing, dict) else None
        if not isinstance(editor, dict):
            continue
        catalog = editor.get("featureCatalog")
        if isinstance(catalog, Mapping) and catalog.get("schema") == SEQUENCE_FEATURE_CATALOG_SCHEMA:
            editor["featureCatalog"] = promote_legacy_feature_catalog(catalog)
    return document, facts


JSON_SESSION_OUTPUT_SUFFIXES = (".json", ".json.gz")


def refuse_json_session_output(path: str | Path) -> None:
    """Refuse a Session output name that says JSON: a Session is written as a Session file only."""

    if str(path).lower().endswith(JSON_SESSION_OUTPUT_SUFFIXES):
        raise ValidationError(
            "A Session is saved as a .gbdraw file. Use a name that ends with .gbdraw, "
            f"not .json or .json.gz: {path}"
        )


def _write_validated_session_file(
    path: str | Path,
    expanded_payload: Mapping[str, object],
    sources: Mapping[str, SessionResourceSource],
    *,
    overwrite: bool,
) -> None:
    """Write a Session 51 that ``validate_session`` already accepted as the Session file ``path``.

    Only writers whose payload cannot have changed since its validation call
    this: a ``SessionDocument`` (validated when it was built) and the CLI
    sidecar, which writes what ``build_session_json`` returned. Each
    resource's bytes go to its entry, read from its source in ``sources``.
    """

    refuse_json_session_output(path)
    serialized_payload = compact_session_feature_catalog(expanded_payload)
    if serialized_payload.get("version") != CURRENT_SESSION_VERSION:
        raise ValidationError(f"A Session file holds a Session version {CURRENT_SESSION_VERSION}.", diagnostic=SESSION_RESOURCES_INVALID)
    resources = serialized_payload.get("resources", {})
    write_session_container(
        serialized_payload,
        {resource_id: SessionResources(resources, sources).source(resource_id) for resource_id in resources},
        path,
        overwrite=overwrite,
    )


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
                number = parse_ascii_nonnegative_integer(index)
                if number is None:
                    raise ValidationError(f"Invalid session slot: {slot}")
                parts.append(number)
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


__all__ = [
    "AnnotationTargetMigration",
    "HashRuleMigration",
    "CURRENT_SESSION_VERSION",
    "check_session_file_version",
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
    "FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION",
    "MODE_SCOPED_SESSION_MAX_VERSION",
    "MODE_SCOPED_SESSION_MIN_VERSION",
    "ROW_KEYED_LEGEND_SESSION_MIN_VERSION",
    "SESSION_FORMAT",
    "SUPPORTED_SESSION_VERSIONS",
    "SessionBuildContext",
    "SessionFileBinding",
    "SessionRunSpec",
    "classify_raw_losat_cache_entry",
    "compact_session_feature_catalog",
    "decode_depth_payload",
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
    "migrate_session_hash_rules",
    "request_less_hash_rules_need_drawn_hashes",
    "migrate_session_feature_edits",
    "normalize_current_session_artifacts",
    "legend_edits_by_row_key",
    "safe_embedded_filename",
    "session_legend_edits_by_row_key",
    "file_resource_entry",
    "session_depth_source_widths",
    "session_mode",
    "split_draft_into_modes",
    "mode_split_palette_colors",
    "migrate_session_flat_draft",
    "SessionDraftMigration",
    "validate_session",
    "validate_current_session_artifacts",
    "validate_current_web_state_field_names",
]
