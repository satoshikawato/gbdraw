#!/usr/bin/env python3
"""Generate the browser's constants and typedefs of persisted contracts, and
their JSON Schema and model diagrams for the documentation."""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Mapping, get_args, get_type_hints

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "gbdraw" / "web" / "js" / "contracts" / "schema.generated.js"
REQUEST_TARGET = TARGET.with_name("request.generated.js")
SESSION_TARGET = TARGET.with_name("session.generated.js")
REQUEST_SHAPE_TARGET = TARGET.with_name("request-shape.generated.js")
SESSION_SHAPE_TARGET = TARGET.with_name("session-shape.generated.js")
RULE_SYNTHESIS_TARGET = TARGET.with_name("rule-synthesis.generated.js")
CARD_WINDOWS_TARGET = TARGET.with_name("card-windows.generated.js")
LOCATE_TARGET = TARGET.with_name("locate.generated.js")
DERIVED_RESULT_SHAPE_TARGET = TARGET.with_name("derived-result-shape.generated.js")
# The built-in color tables are data, not a persisted contract: a Session saves
# the rows a table adds, so the module stays outside `contracts/`.
BUILTIN_COLOR_TABLES_TARGET = TARGET.parent.parent / "builtin-color-tables.generated.js"
# The colours of the similarity rings without a colour of their own: data, like the
# built-in color tables, so the module stays outside `contracts/` too.
CONSERVATION_COLORS_TARGET = TARGET.parent.parent / "utils" / "conservation-colors.generated.js"
# The default hit colours of Sequence search: data, like the conservation colours.
SEQUENCE_SEARCH_COLORS_TARGET = TARGET.parent.parent / "utils" / "sequence-search-colors.generated.js"
# The programs of Sequence search: data owned by the same Python module.
SEQUENCE_SEARCH_PROGRAMS_TARGET = TARGET.parent.parent / "utils" / "sequence-search-programs.generated.js"
# The names of the records a Result keeps of what Python drew; Python writes them on a Legend it draws with
# row edits, and the Result executor writes them, so the module stays outside `contracts/`.
RESULT_RECORDS_TARGET = TARGET.parent.parent / "services" / "result-records.generated.js"
# The attribute names Python writes on drawn groups that the Web reads in a Result, outside `contracts/` too.
DRAWN_ATTRIBUTES_TARGET = TARGET.parent.parent / "services" / "drawn-attributes.generated.js"
LABEL_FIT_WORDS_TARGET = REPO_ROOT / "tools" / "i18n" / "label-fit-words.generated.mjs"
COLLINEARITY_WORDS_TARGET = REPO_ROOT / "tools" / "i18n" / "collinearity-words.generated.mjs"
# Standard JSON Schema of the stored Session and the canonical request, for
# readers outside gbdraw: documentation, not a Web module.
SESSION_JSON_SCHEMA_TARGET = REPO_ROOT / "docs" / "REFERENCE" / "schemas" / "gbdraw-session.schema.json"
REQUEST_JSON_SCHEMA_TARGET = SESSION_JSON_SCHEMA_TARGET.with_name("gbdraw-render-request.schema.json")
# The page of the schemas, with model diagrams of the same descriptions.
SCHEMA_README_TARGET = SESSION_JSON_SCHEMA_TARGET.with_name("README.md")
sys.path.insert(0, str(REPO_ROOT))

# (Web name, Python module that defines the value, Python name). The Web keeps
# its own names; the Python module is the only place a value is written.
_CONSTANTS: tuple[tuple[str, str, str], ...] = (
    ("SESSION_VERSION", "gbdraw.session_io", "CURRENT_SESSION_VERSION"),
    ("SUPPORTED_SESSION_VERSIONS", "gbdraw.session_io", "SUPPORTED_SESSION_VERSIONS"),
    ("TYPED_DRAFT_SESSION_VERSION", "gbdraw.session_io", "TYPED_DRAFT_SESSION_MIN_VERSION"),
    ("MODE_SCOPED_SESSION_VERSION", "gbdraw.session_io", "MODE_SCOPED_SESSION_MIN_VERSION"),
    ("ROW_KEYED_LEGEND_SESSION_VERSION", "gbdraw.session_io", "ROW_KEYED_LEGEND_SESSION_MIN_VERSION"),
    ("FILE_IMPORT_ORIGIN_SESSION_VERSION", "gbdraw.session_io", "FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION"),
    ("DISPLAY_PLACEMENT_DRAFT_SESSION_MIN_VERSION", "gbdraw.session_io", "DISPLAY_PLACEMENT_DRAFT_SESSION_MIN_VERSION"),
    ("CURRENT_AUTHORITY_SESSION_MIN_VERSION", "gbdraw.session_io", "CURRENT_AUTHORITY_SESSION_MIN_VERSION"),
    ("CURRENT_ARTIFACT_SESSION_MIN_VERSION", "gbdraw.session_io", "CURRENT_ARTIFACT_SESSION_MIN_VERSION"),
    ("CANONICAL_SESSION_MIN_VERSION", "gbdraw.session_io", "CANONICAL_SESSION_MIN_VERSION"),
    ("DRAWN_HASH_RULE_SESSION_MAX_VERSION", "gbdraw.session_io", "DRAWN_HASH_RULE_SESSION_MAX_VERSION"),
    ("MAIN_DISPLAY_FRAME_SESSION_VERSION", "gbdraw.api.session_compat", "MAIN_DISPLAY_FRAME_SESSION_VERSION"),
    ("DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION", "gbdraw.api.session_compat",
     "DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION"),
    ("CANONICAL_REQUEST_SCHEMA", "gbdraw.session_request_codec", "CANONICAL_REQUEST_SCHEMA"),
    ("DISPLAY_PLACEMENT_SCHEMA", "gbdraw.session_request_codec", "DISPLAY_PLACEMENT_SCHEMA"),
    ("FEATURE_OVERRIDE_SCHEMA", "gbdraw.session_request_codec", "FEATURE_OVERRIDE_SCHEMA"),
    ("STORED_REGION_FIELDS", "gbdraw.session_request_codec", "STORED_REGION_FIELDS"),
    ("SUPPORTED_CANONICAL_REQUEST_SCHEMAS", "gbdraw.session_request_codec", "SUPPORTED_CANONICAL_REQUEST_SCHEMAS"),
    ("SIMILARITY_ALIGNMENT_PLAN_SCHEMA", "gbdraw.layout.similarity_alignment", "SIMILARITY_ALIGNMENT_PLAN_SCHEMA"),
    ("FEATURE_CATALOG_SCHEMA", "gbdraw.session_io", "CURRENT_FEATURE_CATALOG_SCHEMA"),
    ("ANCHOR_FEATURE_CATALOG_SCHEMA", "gbdraw.web_support.feature_catalog", "ANCHOR_FEATURE_CATALOG_SCHEMA"),
    ("LEGACY_FEATURE_CATALOG_SCHEMA", "gbdraw.web_support.feature_catalog", "LEGACY_FEATURE_CATALOG_SCHEMA"),
    # The schema 5 catalog of Sessions 46-50 stored the sequences; its promotion to FEATURE_CATALOG_SCHEMA drops them.
    ("SEQUENCE_FEATURE_CATALOG_SCHEMA", "gbdraw.session_io", "SEQUENCE_FEATURE_CATALOG_SCHEMA"),
    # The catalog schema an exported Interactive SVG embeds, with every sequence (the export profile, D5-02).
    ("INTERACTIVE_EMBED_SCHEMA", "gbdraw.web_support.feature_catalog", "INTERACTIVE_EMBED_SCHEMA"),
    # What a catalog of schema 5 or older stored of the sequences: the promotion drops it, schema 6 refuses it.
    ("STORED_SEQUENCE_KEYS", "gbdraw.web_support.feature_catalog", "STORED_SEQUENCE_KEYS"),
    # Why a match span's bases cannot be read (the `readMatchSequences` reply's reasons).
    ("MATCH_SPAN_REASONS", "gbdraw.web_support.feature_catalog", "MATCH_SPAN_REASONS"),
    ("FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION", "gbdraw.session_io", "FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION"),
    ("PROTEIN_LOSAT_CACHE_SCHEMA", "gbdraw.analysis.protein_colinearity", "PROTEIN_LOSAT_CACHE_SCHEMA"),
    ("LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA", "gbdraw.analysis.protein_colinearity", "LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA"),
    ("PROTEIN_IDENTITY_MANIFEST_SCHEMA", "gbdraw.analysis.protein_colinearity", "PROTEIN_IDENTITY_MANIFEST_SCHEMA"),
    ("ORTHOGROUP_INFERENCE_VERSION", "gbdraw.analysis.protein_colinearity", "ORTHOGROUP_INFERENCE_VERSION"),
    ("NUCLEOTIDE_LOSAT_CACHE_SCHEMA", "gbdraw.comparisons.losat_jobs", "NUCLEOTIDE_LOSAT_CACHE_SCHEMA"),
    ("CIRCULAR_CONSERVATION_FLOW", "gbdraw.comparisons.circular_losat", "CIRCULAR_CONSERVATION_FLOW"),
    ("LOSAT_DERIVED_CACHE_SCHEMA", "gbdraw.analysis.protein_artifacts", "CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA"),
    ("LEGACY_LOSAT_DERIVED_CACHE_SCHEMA", "gbdraw.session_io", "LEGACY_LOSAT_DERIVED_CACHE_SCHEMA"),
    ("LEGACY_PROTEIN_CANDIDATE_SCHEMA", "gbdraw.session_io", "LEGACY_PROTEIN_CANDIDATE_SCHEMA"),
    ("DEPTH_FILE_SCHEMA", "gbdraw.session_io", "DEPTH_FILE_SCHEMA"),
    ("DEPTH_FILE_ENCODING", "gbdraw.session_io", "DEPTH_FILE_ENCODING"),
    ("GALLERY_SESSION_DRAWN_DRAWINGS", "gbdraw.session_contract", "GALLERY_SESSION_DRAWN_DRAWINGS"),
    # The Session file's entries and its refusals (the twin's vectors pin the rest).
    ("SESSION_DOCUMENT_ENTRY", "gbdraw.session_container", "SESSION_DOCUMENT_ENTRY"),
    ("SESSION_RESOURCE_ENTRY_PREFIX", "gbdraw.session_container", "SESSION_RESOURCE_ENTRY_PREFIX"),
    ("SESSION_CONTAINER_DAMAGED_MESSAGE", "gbdraw.session_container", "SESSION_CONTAINER_DAMAGED_MESSAGE"),
    ("SESSION_CONTAINER_NO_DOCUMENT_MESSAGE", "gbdraw.session_container", "SESSION_CONTAINER_NO_DOCUMENT_MESSAGE"),
    ("SESSION_ENTRY_STORED", "gbdraw.session_container", "SESSION_ENTRY_STORED"),
    ("SESSION_ENTRY_DEFLATED", "gbdraw.session_container", "SESSION_ENTRY_DEFLATED"),
    ("SESSION_ENTRY_ENCRYPTED_FLAGS", "gbdraw.session_container", "SESSION_ENTRY_ENCRYPTED_FLAGS"),
    ("SESSION_ENTRY_DATE_TIME", "gbdraw.session_container", "SESSION_ENTRY_DATE_TIME"),
    ("LEGEND_POSITIONS", "gbdraw.api.options", "LEGEND_POSITIONS"),
    # A "This feature only" color rule names a feature instance; the Source
    # recipe names the CLI record of its record key.
    ("FEATURE_INSTANCE_SEPARATOR", "gbdraw.features.ids", "FEATURE_INSTANCE_SEPARATOR"),
    ("CLI_RECORD_KEY_PREFIX", "gbdraw.api.record_planning", "CLI_RECORD_KEY_PREFIX"),
    ("EXPANDED_RECORD_KEY_SEPARATOR", "gbdraw.api.requests", "EXPANDED_RECORD_KEY_SEPARATOR"),
    ("CONFIG_UNSET_NUMBER_SPELLINGS", "gbdraw.core.numeric", "UNSET_NUMBER_SPELLINGS"),
    # The most bytes a gzip input file expands to; the Web stops reading past it.
    ("MAX_EXPANDED_INPUT_BYTES", "gbdraw.io.input_text", "MAX_EXPANDED_INPUT_BYTES"),
    # The most characters of an input file's name a refusal names; the page cuts at the same count.
    ("INPUT_NAME_CHARACTERS", "gbdraw.web_support.error_adapter", "INPUT_NAME_CHARACTERS"),
    # Session field sets: the keys of the classes of gbdraw/session_contract.py.
    ("FLAT_DRAFT_TOP_LEVEL_FIELDS", "gbdraw.session_io", "FLAT_DRAFT_TOP_LEVEL_FIELDS"),
    ("OTHER_MODE_RESULT_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_FIELDS"),
    ("OTHER_MODE_RESULT_UI_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_UI_FIELDS"),
    ("OTHER_MODE_RESULT_EDITOR_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_EDITOR_FIELDS"),
    ("OTHER_MODE_RESULT_LEGEND_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_LEGEND_FIELDS"),
    ("OTHER_MODE_RESULT_RUN_METADATA_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_RUN_METADATA_FIELDS"),
    # The rows of the run notice lists in `runMetadata`: their keys and the values of their fixed fields.
    ("ANNOTATION_WARNING_FIELDS", "gbdraw.session_io", "ANNOTATION_WARNING_FIELDS"),
    ("ANNOTATION_WARNING_CODES", "gbdraw.session_io", "ANNOTATION_WARNING_CODES"),
    ("COMPARISON_WARNING_FIELDS", "gbdraw.session_io", "COMPARISON_WARNING_FIELDS"),
    ("COMPARISON_WARNING_CODES", "gbdraw.session_io", "COMPARISON_WARNING_CODES"),
    ("FEATURE_IDENTITY_NOTICE_FIELDS", "gbdraw.session_io", "FEATURE_IDENTITY_NOTICE_FIELDS"),
    ("FEATURE_IDENTITY_NOTICE_STATUSES", "gbdraw.session_io", "FEATURE_IDENTITY_NOTICE_STATUSES"),
    ("FEATURE_IDENTITY_NOTICE_KINDS", "gbdraw.session_io", "FEATURE_IDENTITY_NOTICE_KINDS"),
    ("CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS", "gbdraw.session_io", "CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS"),
    # The file slots of a Circular File card (a `circularSources` row) besides its uid, in the card's order.
    ("CIRCULAR_SOURCE_FILE_SLOTS", "gbdraw.session_io", "CIRCULAR_SOURCE_FILE_SLOTS"),
    ("LEGACY_SVG_ATTRIBUTE_RENAMES", "gbdraw.render.interactive_svg", "LEGACY_SVG_ATTRIBUTE_RENAMES"),
    # A former palette name reads as the current one (Session Load).
    ("PALETTE_NAME_ALIASES", "gbdraw.io.colors", "PALETTE_NAME_ALIASES"),
    # The PNG DPI choices a Session's ui.downloadDpi names; the Web keeps their labels.
    ("SESSION_PNG_DPI_CHOICES", "gbdraw.session", "SESSION_PNG_DPI_CHOICES"),
    # The DPI of a PNG without one: the Web's DPI menu starts there, and the CLI and Python draw at it.
    ("DEFAULT_PNG_DPI", "gbdraw.render.export_limits", "DEFAULT_PNG_DPI"),
    # What an absent axis label style reads as: every writer before the setting drew numbered ticks.
    ("LEGACY_SCALAR_AXIS_LABELS", "gbdraw.session_legacy_values", "LEGACY_SCALAR_AXIS_LABELS"),
    # The params keys each track renderer takes without reading them; the Source recipe leaves them out.
    ("CIRCULAR_SLOT_UNREAD_PARAM_KEYS", "gbdraw.tracks.circular", "CIRCULAR_SLOT_UNREAD_PARAM_KEYS"),
    # The label fit row statuses Apply can choose, the preferred first: the dialog selects by this order.
    ("LABEL_FIT_OFFERED_ROW_STATUSES", "gbdraw.labels.fit", "OFFERED_ROW_STATUSES"),
    ("LINEAR_SLOT_UNREAD_PARAM_KEYS", "gbdraw.tracks.linear", "LINEAR_SLOT_UNREAD_PARAM_KEYS"),
    # The media type a Session writer stores for a resource kind gbdraw makes, and for an input file (OV-2964).
    ("WRITTEN_RESOURCE_MEDIA_TYPES", "gbdraw.session_io", "WRITTEN_RESOURCE_MEDIA_TYPES"),
    ("INPUT_FILE_MEDIA_TYPE", "gbdraw.session_io", "INPUT_FILE_MEDIA_TYPE"),
    # The first Session whose comparison-plan edges hold window-pair names; Load drops them from an older file.
    ("WINDOW_LOSAT_FILENAMES_SESSION_VERSION", "gbdraw.session_io", "WINDOW_LOSAT_FILENAMES_SESSION_VERSION"),
    ("AUTOMATIC_ANNOTATION_LANE_SIDES", "gbdraw.annotations.planning", "AUTOMATIC_ANNOTATION_LANE_SIDES"),
    # Input summary warnings of an input drawn in full (OV-2954): no "not drawn" heading.
    ("INPUT_WARNINGS_OF_DRAWN_INPUT", "gbdraw.api.input_summary", "INPUT_WARNINGS_OF_DRAWN_INPUT"),
    # The genetic codes LOSAT reads, and how a refusal names them: Generate refuses any other before a search (OV-3326).
    ("LOSAT_GENETIC_CODES", "gbdraw.comparisons.losat_runtime", "LOSAT_GENETIC_CODES"),
    ("LOSAT_GENETIC_CODE_TEXT", "gbdraw.comparisons.losat_runtime", "LOSAT_GENETIC_CODE_TEXT"),
)
# (Web name, Python module, function): values a Python function derives from
# its owner, written in this order after the constants.
_DERIVED: tuple[tuple[str, str, str], ...] = (
    ("EXPECTED_WEB_RUNTIME_CAPABILITIES", "gbdraw.web_support.capabilities", "get_web_runtime_capabilities"),
    ("CONFIG_UNSET_NUMBER_PATHS", "gbdraw.config.modify", "config_paths_reading_unset_spellings"),
    ("LEGACY_REQUEST_TABLES", "gbdraw.session_request_codec", "describe_legacy_request_tables"),
    ("SESSION_DRAWINGS_LAYOUT", "gbdraw.session_drawings", "describe_session_drawings_layout"),
    ("FRESH_PROTEIN_COMPARISON_DEFAULTS", "gbdraw.session_request_codec", "describe_fresh_protein_comparison"),
    # The sides of a request without `diagramOptions.output`, as Python draws it (OV-2491, R17-38).
    ("ABSENT_CIRCULAR_OUTPUT", "gbdraw.session_request_codec", "describe_absent_circular_output"),
    ("ABSENT_LINEAR_OUTPUT", "gbdraw.session_request_codec", "describe_absent_linear_output"),
    # The files of each record source kind (resource kind, source field), in the order the source names them.
    ("RECORD_SOURCE_FILES", "gbdraw.session_request_codec", "describe_record_source_files"),
    # The words of --feature_density_track and where a request track holds each value: the Source recipe's spelling.
    ("FEATURE_DENSITY_TRACK_WORDS", "gbdraw.session_request_codec", "describe_feature_density_track_words"),
    # The Auto organism name rule's data; the Web fast path runs its steps (R4).
    ("ORGANISM_NAME_RULE", "gbdraw.core.record_metadata", "organism_name_rule_data"),
    # The raw-cache args keep one spelling; the LOSAT Wasm reads the runtime spelling.
    ("LOSAT_SEARCH_FLAGS", "gbdraw.comparisons.losat_runtime", "losat_search_flags"),
    # The style of a new annotation set: Python's style defaults with the highlight fill (OV-2813).
    ("DRAFT_ANNOTATION_STYLE", "gbdraw.session_request_codec", "describe_draft_annotation_style"),
)
# (Web name, Python module, dataclass, field): the values of a field typed as a
# string Literal (typing.get_args, in the owner's order) and the field's
# default, written after the derived values. The Web shows these values in the
# UI language (services/pairwise-match-popup.js).
_LITERAL_FIELDS: tuple[tuple[str, str, str, str], ...] = (
    ("ORTHOGROUP_MEMBER_ROLE", "gbdraw.analysis.protein_colinearity", "OrthogroupMember", "role"),
    ("ORTHOGROUP_MEMBER_CONFIDENCE", "gbdraw.analysis.protein_colinearity", "OrthogroupMember", "confidence"),
    # How a depth or GC Percent axis is labelled (Web: the Numbered Ticks checkboxes).
    ("SCALAR_AXIS_LABELS", "gbdraw.config.models.objects", "ObjectsDepthConfig", "axis_labels"),
)
# (Web name, Python module, name): the regular expression a Python reader reads
# a text with, as an ECMA-262 pattern that matches the same texts
# (``_portable_pattern``), written last. The Web matches it against the whole
# text, as Python's ``fullmatch`` does.
_PATTERNS: tuple[tuple[str, str, str], ...] = (
    # A `#N` record selector, as the renderer reads it (OV-2960).
    ("RECORD_INDEX_SELECTOR_PATTERN", "gbdraw.io.record_select", "RECORD_INDEX_SELECTOR_PATTERN"),
)
_IDENTIFIER = re.compile(r"[A-Za-z_$][A-Za-z0-9_$]*|0|[1-9][0-9]*")


def _scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and math.isfinite(value):
        return repr(value)
    if not isinstance(value, (int, str)):
        raise TypeError(f"Unsupported Web contract value: {value!r}")
    if isinstance(value, int):
        return str(value)
    return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"


def _frozen(value: object, indent: str = "") -> str:
    """Render JSON-compatible data as nested ``Object.freeze`` literals."""
    inner = indent + "  "
    if isinstance(value, dict):
        fields = []
        for key, item in value.items():
            # A key that is no identifier or index is quoted (a dotted config path).
            name = str(key) if _IDENTIFIER.fullmatch(str(key)) else json.dumps(key)
            fields.append(f"{inner}{name}: {_frozen(item, inner)}")
        return "Object.freeze({\n" + ",\n".join(fields) + f"\n{indent}}})"
    if isinstance(value, list):
        if any(isinstance(item, (dict, list)) for item in value):
            # One nested item per line (the rows of a table).
            items = ",\n".join(f"{inner}{_frozen(item, inner)}" for item in value)
            return "Object.freeze([\n" + items + f"\n{indent}])"
        return "Object.freeze([" + ", ".join(_scalar(item) for item in value) + "])"
    return _scalar(value)


def _declaration(web_name: str, value: object) -> str:
    if isinstance(value, tuple) and value and all(isinstance(item, str) for item in value):
        # Strings in the owner's order, each typed as itself (a typedef may name them).
        return (
            f"/** @type {{readonly [{', '.join(_scalar(item) for item in value)}]}} */\n"
            f"export const {web_name} = {_frozen(list(value))};\n"
        )
    if isinstance(value, tuple) and value and all(type(item) is int for item in value):
        # A fixed tuple of numbers, such as a date and time.
        return (
            f"/** @type {{readonly [{', '.join('number' for _ in value)}]}} */\n"
            f"export const {web_name} = {_frozen(list(value))};\n"
        )
    if isinstance(value, (frozenset, set)):
        items = "string" if all(isinstance(item, str) for item in value) else "number"
        return (
            f"/** @type {{readonly {items}[]}} */\n"
            f"export const {web_name} = {_frozen(sorted(value))};\n"
        )
    if isinstance(value, Mapping) and all(isinstance(item, tuple) for item in value.values()):
        # Strings per key, in the owner's order.
        keys = " | ".join(_scalar(key) for key in value)
        return (
            f"/** @type {{Readonly<Record<{keys}, readonly string[]>>}} */\n"
            f"export const {web_name} = {_frozen({key: list(item) for key, item in value.items()})};\n"
        )
    if isinstance(value, Mapping) and all(isinstance(item, str) for item in value.values()):
        # One string per key, in the owner's order.
        keys = " | ".join(_scalar(key) for key in value)
        return (
            f"/** @type {{Readonly<Record<{keys}, string>>}} */\n"
            f"export const {web_name} = {_frozen(dict(value))};\n"
        )
    if isinstance(value, dict):
        return (
            "/** @type {Readonly<Record<number, number>>} */\n"
            f"export const {web_name} = {_frozen(dict(sorted(value.items())))};\n"
        )
    return f"export const {web_name} = {_scalar(value)};\n"


def _owner(module: str, name: str) -> str:
    return f"// {module.replace('.', '/')}.py: {name}\n"


def _literal_field_declaration(web_name: str, owner: type, field_name: str) -> str:
    """The values and the default of a dataclass field typed as a string Literal."""
    values = get_args(get_type_hints(owner)[field_name])
    default = next(field.default for field in dataclasses.fields(owner) if field.name == field_name)
    if not values or not all(isinstance(value, str) for value in values) or default not in values:
        raise TypeError(f"{owner.__name__}.{field_name} is not a string Literal field with a listed default")
    union = f"readonly ({' | '.join(_scalar(value) for value in values)})[]"
    return (
        f"/** @type {{Readonly<{{ values: {union}, default: {_scalar(default)} }}>}} */\n"
        f"export const {web_name} = Object.freeze({{\n"
        f"  values: /** @type {{{union}}} */ ({_frozen(list(values))}),\n"
        f"  default: {_scalar(default)}\n"
        "});\n"
    )


def render_web_contracts_module() -> str:
    parts = [
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from the Python owner named\n"
        "// above each constant; do not edit manually.\n"
    ]
    for web_name, module, name in _CONSTANTS:
        value = getattr(importlib.import_module(module), name)
        parts.append(_owner(module, name) + _declaration(web_name, value))
    for web_name, module, name in _DERIVED:
        value = getattr(importlib.import_module(module), name)()
        declaration = (
            _declaration(web_name, value)
            if isinstance(value, (frozenset, set))
            else f"export const {web_name} = {_frozen(value)};\n"
        )
        parts.append(_owner(module, f"{name}()") + declaration)
    for web_name, module, class_name, field_name in _LITERAL_FIELDS:
        owner = getattr(importlib.import_module(module), class_name)
        parts.append(
            _owner(module, f"{class_name}.{field_name}") + _literal_field_declaration(web_name, owner, field_name)
        )
    for web_name, module, name in _PATTERNS:
        pattern = getattr(importlib.import_module(module), name)
        parts.append(_owner(module, name) + f"export const {web_name} = {_scalar(_portable_pattern(pattern))};\n")
    return "\n".join(parts)


_JSDOC_SCALARS = {
    "string": "string",
    "integer": "number",
    "number": "number",
    "boolean": "boolean",
    "null": "null",
    "unknown": "unknown",
}


def _jsdoc_literal(value: object) -> str:
    if value is None or isinstance(value, bool):
        return {None: "null", True: "true", False: "false"}[value]
    return _scalar(value)


def _jsdoc_type(shape: dict[str, Any]) -> str:
    """Render one shape of ``gbdraw.contract_shapes`` as a JSDoc type."""
    kind = shape["kind"]
    if kind in _JSDOC_SCALARS:
        return _JSDOC_SCALARS[kind]
    if kind == "ref":
        return shape["name"]
    if kind == "literal":
        # A literal with no value admits nothing: an optional field that must be absent.
        return " | ".join(_jsdoc_literal(value) for value in shape["values"]) or "never"
    if kind == "union":
        return " | ".join(_jsdoc_type(option) for option in shape["options"])
    if kind == "array":
        items = _jsdoc_type(shape["items"])
        several = shape["items"]["kind"] == "union" or len(shape["items"].get("values", ())) > 1
        return f"({items})[]" if several else f"{items}[]"
    if kind == "tuple":
        return "[" + ", ".join(_jsdoc_type(item) for item in shape["items"]) + "]"
    if kind == "record":
        return f"Record<string, {_jsdoc_type(shape['values'])}>"
    if kind != "object":
        raise TypeError(f"Unsupported contract shape: {kind!r}")
    if not shape["fields"] and not shape["open"]:
        return "Record<string, never>"
    return "{ " + ", ".join(_members(shape)) + " }"


def _members(shape: dict[str, Any]) -> list[str]:
    """The members of an object type; an open object ends with an index signature."""
    members = [
        f"{_property_name(name)}{'?' if name in shape['optional'] else ''}: {_jsdoc_type(item)}"
        for name, item in shape["fields"].items()
    ]
    return [*members, "[key: string]: unknown"] if shape["open"] else members


def _property_name(name: str) -> str:
    if not _IDENTIFIER.fullmatch(name):
        raise TypeError(f"Unsupported Web contract key: {name!r}")
    return name


def _typedefs(types: dict[str, dict[str, Any]]) -> str:
    """Render named shapes as JSDoc typedefs; an object lists its properties.

    An open object is written as a type literal: ``@property`` has no index signature.
    """
    parts = []
    for name, shape in types.items():
        if shape["kind"] == "object" and shape["open"]:
            members = ",\n".join(f" *   {member}" for member in _members(shape))
            parts.append(f"/**\n * @typedef {{{{\n{members}\n * }}}} {name}\n */\n")
            continue
        if shape["kind"] != "object" or not shape["fields"]:
            parts.append(f"/** @typedef {{{_jsdoc_type(shape)}}} {name} */\n")
            continue
        lines = [f" * @typedef {{object}} {name}"]
        for field, item in shape["fields"].items():
            key = f"[{_property_name(field)}]" if field in shape["optional"] else _property_name(field)
            lines.append(f" * @property {{{_jsdoc_type(item)}}} {key}")
        parts.append("/**\n" + "\n".join(lines) + "\n */\n")
    return "".join(parts)


def render_request_typedefs_module() -> str:
    from gbdraw.session_request_codec import describe_canonical_request

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_canonical_request()\n"
        "// in gbdraw/session_request_codec.py; do not edit manually.\n"
        "// The canonical render request of the current schema, as Python writes and reads it.\n"
        "\n"
        + _typedefs(describe_canonical_request()["types"])
        + "export {};\n"
    )


def render_rule_synthesis_typedefs_module() -> str:
    from gbdraw.web_support.rule_matching import describe_rule_helpers

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_rule_helpers()\n"
        "// in gbdraw/web_support/rule_matching.py; do not edit manually.\n"
        "// The diagram helpers `synthesizeColorRules` and `evaluateRules`: what they read\n"
        "// and answer (not persisted).\n"
        "\n"
        + _typedefs(describe_rule_helpers()["types"])
        + "export {};\n"
    )


def render_locate_typedefs_module() -> str:
    from gbdraw.web_support.locate import describe_locate

    description = describe_locate()
    published = {**description, "types": ecma_types(description["types"])}
    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_locate()\n"
        "// in gbdraw/web_support/locate.py; do not edit manually.\n"
        "// The diagram helpers of Find in other files: their inputs and replies (not persisted), as typedefs and,\n"
        "// for gbdraw/web/js/utils/contract-shape.js, the replies' description as data.\n"
        "\n"
        "/** @import { RootedContractDescription } from '../utils/contract-shape.js' */\n"
        "\n"
        + _typedefs(description["types"])
        + "\n/** @type {RootedContractDescription} */\n"
        + f"export const LOCATE_DESCRIPTION = {_frozen_json(json.loads(json.dumps(published)))};\n"
    )


def render_card_windows_typedefs_module() -> str:
    from gbdraw.comparisons.losat_jobs import PER_FILE_WARNING_TEXT
    from gbdraw.web_support.linear_units import describe_card_window_plan

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_card_window_plan()\n"
        "// in gbdraw/web_support/linear_units.py and PER_FILE_WARNING_TEXT in gbdraw/comparisons/losat_jobs.py;\n"
        "// do not edit manually.\n"
        "// The diagram helpers `planCardWindows` and `planRecordPlaces`: their requests and answers, and the answer of\n"
        "// `planGenomeScales` (not persisted).\n"
        "\n"
        + _typedefs(describe_card_window_plan()["types"])
        + "\n/** The text of a `planGenomeScales` warning (`GenomeScaleWarning.text`), shown with `tm()`: a catalog key. */\n"
        + f"export const GENOME_SCALE_WARNING_TEXT = {json.dumps(PER_FILE_WARNING_TEXT)};\n"
    )


def _references(shape: object) -> set[str]:
    """The names a shape refers to."""
    if isinstance(shape, dict):
        own = {shape["name"]} if shape.get("kind") == "ref" else set()
        return own.union(*(_references(item) for item in shape.values()))
    return set().union(*(_references(item) for item in shape)) if isinstance(shape, list) else set()


def render_session_typedefs_module() -> str:
    from gbdraw.session_contract import (
        describe_session_admission,
        describe_session_container_resource,
        describe_session_document,
    )
    from gbdraw.session_request_codec import describe_canonical_request

    description = describe_session_document()
    types = description["types"]
    admission = describe_session_admission()
    shared = set(types) & set(admission["types"])
    if shared:
        raise RuntimeError(f"The admission description reuses writer type name(s): {sorted(shared)}")
    bands = " | ".join(band["root"] for band in admission["bands"])
    # A name the Session does not define is the request's; a leaf may import a leaf (R13).
    imported = sorted(_references(types) - set(types))
    unknown = set(imported) - set(describe_canonical_request()["types"])
    if unknown:
        raise RuntimeError(f"The Session description refers to undescribed type(s): {sorted(unknown)}")
    root = types[description["root"]]["fields"]
    # A closed band holds no other SessionInput key (the reader drops it): the
    # typedef says so, as the open bands below 46 say it of `otherModeResult`.
    band_types = dict(admission["types"])
    for band in admission["bands"]:
        shape = band_types[band["root"]]
        if not shape["open"]:
            absent = [key for key in types["SessionInput"]["fields"] if key not in shape["fields"]]
            band_types[band["root"]] = {
                **shape,
                "fields": {**shape["fields"], **{key: {"kind": "literal", "values": []} for key in absent}},
                "optional": [*shape["optional"], *absent],
            }
    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_session_document()\n"
        "// in gbdraw/session_contract.py; do not edit manually.\n"
        f"// {description['root']} is a Session of the current version as its writers write it.\n"
        "// A reader takes a SessionInput: the keys a Session of any supported version may\n"
        "// hold, with unvalidated values; a SessionEnvelope is one whose format and version\n"
        "// were read. SessionChecked is a stored Session of any supported version as\n"
        "// validate_session proves it (describe_session_admission()): one type per band.\n"
        "\n"
        f"/** @import {{ {', '.join(imported)} }} from './{REQUEST_TARGET.name}' */\n"
        + _typedefs(types)
        + _typedefs(describe_session_container_resource()["types"])
        + "/** @typedef {SessionInput & { format: "
        + f"{_jsdoc_type(root['format'])}, version: {_jsdoc_type(root['version'])}"
        + " }} SessionEnvelope */\n"
        + _typedefs(band_types)
        + f"/** @typedef {{{bands}}} SessionChecked */\n"
        "export {};\n"
    )


_JSON_KEY = re.compile(r"[A-Za-z_$][A-Za-z0-9_$]*")


def _frozen_json(value: object, indent: str = "") -> str:
    """Render parsed JSON as deeply frozen literals, a container on one line when it fits.

    A key that is not an identifier is quoted; ``__proto__`` is refused, since
    a literal would read it as the prototype.
    """
    if value is None or isinstance(value, (bool, int, float, str)):
        return json.dumps(value, allow_nan=False)
    if isinstance(value, dict):
        if "__proto__" in value:
            raise TypeError("A __proto__ key cannot be rendered as an object literal key.")
        items = [(f"{key if _JSON_KEY.fullmatch(key) else json.dumps(key)}: ", item) for key, item in value.items()]
        brackets = "{", "}"
    elif isinstance(value, list):
        items = [("", item) for item in value]
        brackets = "[", "]"
    else:
        raise TypeError(f"Unsupported Web contract value: {value!r}")
    inner = indent + "  "
    line = ", ".join(prefix + _frozen_json(item, inner) for prefix, item in items)
    if "\n" not in line and len(indent) + len(line) <= 100:
        return f"Object.freeze({brackets[0]}{line}{brackets[1]})"
    rows = ",\n".join(f"{inner}{prefix}{_frozen_json(item, inner)}" for prefix, item in items)
    return f"Object.freeze({brackets[0]}\n{rows}\n{indent}{brackets[1]})"


def render_description_module(
    description: Mapping[str, Any], *, export_name: str, owner: str, description_type: str = "ContractDescription"
) -> str:
    """Render a description of ``gbdraw.contract_shapes`` as data that
    ``shapeViolations`` (gbdraw/web/js/utils/contract-shape.js) checks a parsed
    value against. ``owner`` names the Python function it comes from;
    ``description_type`` is the checker's typedef of it: ``ContractDescription``
    (``{"root", "types"}``), ``BandedContractDescription`` (``{"bands", "types"}``)
    or ``RootedContractDescription`` (``{"roots", "types"}``). Its patterns are
    spelled for ECMA-262 (:func:`ecma_types`)."""
    published = {**description, "types": ecma_types(description["types"])}
    return (
        "// @ts-check\n"
        f"// Generated by tools/generate_web_contracts.py from {owner}; do not edit manually.\n"
        "// The description as data, for gbdraw/web/js/utils/contract-shape.js.\n"
        "\n"
        f"/** @import {{ {description_type} }} from '../utils/contract-shape.js' */\n"
        "\n"
        f"/** @type {{{description_type}}} */\n"
        f"export const {export_name} = {_frozen_json(json.loads(json.dumps(published)))};\n"
    )


def render_request_shape_module() -> str:
    from gbdraw.session_request_codec import describe_canonical_request

    return render_description_module(
        describe_canonical_request(),
        export_name="CANONICAL_REQUEST_DESCRIPTION",
        owner="describe_canonical_request() in gbdraw/session_request_codec.py",
    )


def render_session_shape_module() -> str:
    from gbdraw.session_contract import describe_session_admission

    return render_description_module(
        describe_session_admission(),
        export_name="SESSION_ADMISSION_DESCRIPTION",
        owner="describe_session_admission() in gbdraw/session_contract.py",
        description_type="BandedContractDescription",
    )


def render_derived_result_shape_module() -> str:
    from gbdraw.analysis.protein_artifacts import describe_derived_analysis_results

    return render_description_module(
        describe_derived_analysis_results(),
        export_name="DERIVED_RESULT_DESCRIPTION",
        owner="describe_derived_analysis_results() in gbdraw/analysis/protein_artifacts.py",
        description_type="RootedContractDescription",
    )


JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
_SHAPE_BOUNDS = ("minimum", "exclusiveMinimum", "maximum", "minItems", "pattern")


def _runs(values: list[int]) -> list[list[int]]:
    """Sorted integers as runs of consecutive values."""
    runs: list[list[int]] = []
    for value in values:
        if runs and value == runs[-1][-1] + 1:
            runs[-1].append(value)
        else:
            runs.append([value])
    return runs


def _portable_pattern(pattern: str) -> str:
    """``pattern`` as an ECMA-262 pattern that matches what Python's ``re.search`` matches.

    ECMA-262, the dialect of JSON Schema, reads ``\\s`` as other characters
    than Python (it adds U+FEFF and leaves out U+001C-U+001F and U+0085), so
    ``\\s`` and ``\\S`` are spelled as the characters Python matches; no
    character above U+FFFF is one. A pattern whose ``\\s`` sits in a class, or
    whose ``$`` does not end ``\\s*`` (Python's ``$`` also matches before a
    final newline), has no such spelling here and raises.
    """
    if re.search(r"\[[^\]]*\\[sS]|\\\\[sS]|(?<!\\s\*)\$", pattern):
        raise ValueError(f"No portable spelling of the pattern {pattern!r}.")
    spaces = [code for code in range(0x10000) if re.fullmatch(r"\s", chr(code))]
    body = "".join(
        f"\\u{run[0]:04x}" + (f"-\\u{run[-1]:04x}" if len(run) > 1 else "") for run in _runs(spaces)
    )
    return pattern.replace("\\S", f"[^{body}]").replace("\\s", f"[{body}]")


def _ecma_shape(shape: Mapping[str, Any]) -> Any:
    """``shape`` with each pattern spelled for ECMA-262 (:func:`_portable_pattern`)."""
    kind = shape["kind"]
    if kind == "string" and "pattern" in shape:
        return {**shape, "pattern": _portable_pattern(shape["pattern"])}
    if kind in {"union", "tuple"}:
        key = "options" if kind == "union" else "items"
        return {**shape, key: [_ecma_shape(item) for item in shape[key]]}
    if kind in {"array", "record"}:
        key = "items" if kind == "array" else "values"
        return {**shape, key: _ecma_shape(shape[key])}
    if kind == "object":
        return {**shape, "fields": {name: _ecma_shape(item) for name, item in shape["fields"].items()}}
    return shape


def ecma_types(types: Mapping[str, Any]) -> dict[str, Any]:
    """The named shapes of a description as the files this tool writes publish them.

    Their readers, the Web checker ``shapeViolations``
    (gbdraw/web/js/utils/contract-shape.js) and JSON Schema, read a pattern as
    ECMA-262 does, so each pattern is spelled to match what Python's
    ``shape_violations`` matches.
    """
    return {name: _ecma_shape(shape) for name, shape in types.items()}


def json_schema(shape: Mapping[str, Any]) -> Any:
    """One shape of ``gbdraw.contract_shapes`` as a JSON Schema 2020-12 subschema.

    A ``ref`` points into ``$defs``; a bound keeps its name, which JSON Schema
    shares, so a pattern must already be spelled for ECMA-262 (:func:`ecma_types`).
    Python's checker refuses a number with a fraction (``1.0``) where an
    integer or a literal is declared; JSON Schema reads ``1.0`` as ``1``.
    """
    kind = shape["kind"]
    bounds = {name: shape[name] for name in _SHAPE_BOUNDS if name in shape}
    if kind == "unknown":
        return True
    if kind == "ref":
        return {"$ref": f"#/$defs/{shape['name']}"}
    if kind in {"string", "integer", "number", "boolean", "null"}:
        return {"type": kind, **bounds}
    if kind == "literal":
        # A literal with no value admits nothing: an optional field that must be absent.
        values = shape["values"]
        return {"const": values[0]} if len(values) == 1 else {"enum": values} if values else False
    if kind == "union":
        return {"anyOf": [json_schema(option) for option in shape["options"]]}
    if kind == "array":
        return {"type": "array", "items": json_schema(shape["items"]), **bounds}
    if kind == "tuple":
        items = [json_schema(item) for item in shape["items"]]
        prefix = {"prefixItems": items} if items else {}
        return {"type": "array", **prefix, "minItems": len(items), "maxItems": len(items)}
    if kind == "record":
        return {"type": "object", "additionalProperties": json_schema(shape["values"])}
    if kind != "object":
        raise TypeError(f"Unsupported contract shape: {kind!r}")
    defaults = shape.get("defaults", {})
    properties = {}
    for name, item in shape["fields"].items():
        schema = json_schema(item)
        if name in defaults:
            schema = {**({} if schema is True else schema), "default": defaults[name]}
        properties[name] = schema
    required = [name for name in shape["fields"] if name not in shape["optional"]]
    return {
        "type": "object",
        **({"properties": properties} if properties else {}),
        **({"required": required} if required else {}),
        **({} if shape["open"] else {"additionalProperties": False}),
    }


def _json_schema_document(
    root: Mapping[str, Any], types: Mapping[str, Any], *, title: str, owner: str, about: str
) -> str:
    document = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$comment": f"Generated by tools/generate_web_contracts.py from {owner}; do not edit manually.",
        "title": title,
        "description": about,
        **root,
        "$defs": {name: json_schema(shape) for name, shape in ecma_types(types).items()},
    }
    return json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def _versions_text(versions: list[int]) -> str:
    """``[27, 28, 29, 30, 40, 42]`` as ``27–30, 40, 42``."""
    return ", ".join(f"{run[0]}–{run[-1]}" if len(run) > 2 else ", ".join(map(str, run)) for run in _runs(versions))


def render_session_json_schema() -> str:
    from gbdraw.session_contract import describe_session_admission
    from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA

    admission = describe_session_admission()
    versions = [version for band in admission["bands"] for version in band["versions"]]
    return _json_schema_document(
        {"anyOf": [{"$ref": f"#/$defs/{band['root']}"} for band in admission["bands"]]},
        admission["types"],
        title=f"gbdraw Session file, versions {_versions_text(versions)}",
        owner="describe_session_admission() in gbdraw/session_contract.py",
        about=(
            "What gbdraw checks of a stored Session before it reads it, with one definition per band of "
            "versions. The gbdraw command line and Python API refuse a file outside this schema; they also "
            "check values, references and rules across fields that the schema does not state. A Session 51 "
            f"drawing's renderRequest of schema {CANONICAL_REQUEST_SCHEMA} follows gbdraw-render-request.schema.json."
        ),
    )


def render_request_json_schema() -> str:
    from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA, describe_canonical_request

    description = describe_canonical_request()
    return _json_schema_document(
        {"$ref": f"#/$defs/{description['root']}"},
        description["types"],
        title=f"gbdraw canonical render request, schema {CANONICAL_REQUEST_SCHEMA}",
        owner="describe_canonical_request() in gbdraw/session_request_codec.py",
        about=(
            f"A canonical render request of schema {CANONICAL_REQUEST_SCHEMA} as gbdraw writes it: "
            f"the renderRequest of a Session drawing. The gbdraw command line and Python API refuse a schema "
            f"{CANONICAL_REQUEST_SCHEMA} request outside this schema, except a retired value of "
            "collinearityAnchorMode, which they read as its default with a warning; they also check values and "
            "references that the schema does not state. A request of an older schema is converted, not described."
        ),
    )


# A class with more fields than this lists none: the schema holds them.
_DIAGRAM_MEMBERS = 20


def _diagram_type(shape: Mapping[str, Any]) -> str:
    """A short type text for a diagram member: no character Mermaid reads as syntax."""
    kind = shape["kind"]
    if kind == "ref":
        return shape["name"]
    if kind == "literal":
        values = shape["values"]
        if not values:
            return "absent"
        return " or ".join(json.dumps(value).replace('"', "'") for value in values) if len(values) <= 3 else (
            f"one of {len(values)} values"
        )
    if kind == "union":
        return " or ".join(_diagram_type(option) for option in shape["options"])
    if kind in {"array", "record"}:
        items = _diagram_type(shape["items" if kind == "array" else "values"])
        return f"list of {items}" if kind == "array" else f"map of {items}"
    if kind == "tuple":
        return f"{len(shape['items'])} items"
    return "any" if kind == "unknown" else kind


def _class_diagram(
    types: Mapping[str, Any], roots: list[str], *, external: Mapping[str, Any] | None = None, members: bool = True
) -> str:
    """A Mermaid class diagram of the named types reachable from ``roots``.

    A field that refers to a named type is an association labelled with the
    field; a named union is a class its options specialize. A type of
    ``external`` is drawn without members (another diagram shows it), and so
    is every type without ``members``.
    """
    lines = ["classDiagram", "  direction LR"]
    seen: list[str] = []
    pending = list(roots)
    while pending:
        name = pending.pop(0)
        if name in seen:
            continue
        seen.append(name)
        shape = types[name] if external is None or name not in external else {"kind": "ref"}
        if shape["kind"] == "ref":
            lines.append(f"  class {name}")
            continue
        if shape["kind"] == "union":
            lines.append(f"  class {name}")
            lines.append(f"  <<union>> {name}")
            options = [option["name"] for option in shape["options"] if option["kind"] == "ref"]
            lines.extend(f"  {name} <|-- {option}" for option in options)
            pending.extend(options)
            continue
        fields = shape.get("fields", {})
        listed = [] if not members else [f"    {len(fields)} fields"] if len(fields) > _DIAGRAM_MEMBERS else [
            f"    +{field}{'?' if field in shape['optional'] else ''} {_diagram_type(item)}"
            for field, item in fields.items()
        ]
        lines.extend([f"  class {name} {{", *listed, "  }"] if listed else [f"  class {name}"])
        # The fields that refer to each type: one association, named by its fields or their number.
        joined: dict[str, list[str]] = {}
        for field, item in fields.items():
            for target in sorted(_references(item)):
                joined.setdefault(target, []).append(field)
        for target, names in joined.items():
            lines.append(f"  {name} --> {target} : {', '.join(names) if len(names) <= 2 else f'{len(names)} fields'}")
            pending.append(target)
    return "\n".join(lines)


def _session_bands_diagram(admission: Mapping[str, Any]) -> str:
    """A flowchart of the version bands: each band's versions, its top level and the keys it requires."""
    lines = ["flowchart LR"]
    previous = None
    for band in admission["bands"]:
        root = admission["types"][band["root"]]
        required = ", ".join(sorted(field for field in root["fields"] if field not in root["optional"]))
        top_level = "open top level" if root["open"] else f"closed top level of {len(root['fields'])} keys"
        node = f"v{band['minVersion']}"
        label = f"Session {_versions_text(band['versions'])}<br/>{top_level}<br/>requires {required}"
        lines.append(f'  {node}["{label}"]')
        if previous is not None:
            lines.append(f"  {previous} --> {node}")
        previous = node
    return "\n".join(lines)


def render_schema_readme() -> str:
    """The page of the schemas: the files and the model diagrams of their descriptions."""
    from gbdraw import session_io as sio
    from gbdraw.session_contract import describe_session_admission, describe_session_document
    from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA, describe_canonical_request

    admission = describe_session_admission()
    writer = describe_session_document()
    request = describe_canonical_request()
    project = admission["bands"][-1]
    versions = sorted(sio.SUPPORTED_SESSION_VERSIONS)
    mode_scoped = _versions_text([
        version for version in versions
        if sio.MODE_SCOPED_SESSION_MIN_VERSION <= version <= sio.MODE_SCOPED_SESSION_MAX_VERSION
    ])
    figures = [
        (
            "Version bands of the Session reader",
            "Each box is a band of `gbdraw-session.schema.json`: the versions it covers and the top-level "
            "keys it requires (`describe_session_admission()`).",
            _session_bands_diagram(admission),
        ),
        (
            f"A Session {_versions_text(project['versions'])} file",
            "The project and its drawings as the reader checks them; it checks each drawing's parts on the "
            "drawing's one-drawing view (`describe_session_admission()`).",
            _class_diagram(admission["types"], [project["root"]]),
        ),
        (
            f"A Session {mode_scoped} document",
            f"The layout that gbdraw writes for a Session {mode_scoped} (`describe_session_document()`); the "
            f"one-drawing view of a Session {_versions_text(project['versions'])} drawing has this layout. The "
            "reader checks less than it holds: the Session schema states what. A class with more than "
            f"{_DIAGRAM_MEMBERS} fields shows only their number.",
            _class_diagram(writer["types"], [writer["root"]], external=request["types"]),
        ),
        (
            f"Canonical render request, schema {CANONICAL_REQUEST_SCHEMA}",
            "The types of `gbdraw-render-request.schema.json` and the fields that join them "
            "(`describe_canonical_request()`); the schema lists every field.",
            _class_diagram(request["types"], [request["root"]], members=False),
        ),
    ]
    parts = [
        "[Documentation home](../../DOCS.md) | [Tutorials](../../TUTORIALS/README.md) | "
        "[Technical documentation](../README.md) | [FAQ](../../FAQ.md)",
        "",
        "# Session and request JSON Schema",
        "",
        "<!-- Generated by tools/generate_web_contracts.py from gbdraw/session_contract.py and "
        "gbdraw/session_request_codec.py; do not edit manually. -->",
        "",
        f"- [`{SESSION_JSON_SCHEMA_TARGET.name}`]({SESSION_JSON_SCHEMA_TARGET.name}): a stored Session of "
        f"any supported version ({_versions_text(versions)}).",
        f"- [`{REQUEST_JSON_SCHEMA_TARGET.name}`]({REQUEST_JSON_SCHEMA_TARGET.name}): a canonical render "
        f"request of schema {CANONICAL_REQUEST_SCHEMA}.",
        "",
        "Both are JSON Schema draft 2020-12. What gbdraw refuses outside them, and what it reads anyway, is in "
        "[Session and request structure](../session-and-request-compatibility.md#session-and-request-structure).",
        "",
    ]
    for title, caption, diagram in figures:
        parts += [f"## {title}", "", caption, "", "```mermaid", diagram, "```", ""]
    return "\n".join(parts)


def render_builtin_color_tables_module() -> str:
    """The built-in color tables, by name: each one's title and rows as Python reads them."""
    from gbdraw.io.colors import builtin_color_table_data

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/data/color_tables/\n"
        "// (gbdraw/io/colors.py: builtin_color_table_data()); do not edit manually.\n"
        "\n"
        "/**\n"
        " * A built-in color table: its title and its rows, each\n"
        " * [feature_type, qualifier_key, value, color, caption].\n"
        " * @typedef {{ title: string, rows: ReadonlyArray<readonly string[]> }} BuiltinColorTable\n"
        " */\n"
        "/** @type {Readonly<Record<string, BuiltinColorTable>>} */\n"
        f"export const BUILTIN_COLOR_TABLES = {_frozen(builtin_color_table_data())};\n"
    )


def render_conservation_colors_module() -> str:
    """The similarity ring colour schemes, the default scheme, the colour of a ring without one, and their config paths."""
    from gbdraw.analysis.conservation import CONSERVATION_COLOR_SCHEMES
    from gbdraw.api.config import apply_config_overrides
    from gbdraw.config.models.objects import CONSERVATION_COLOR_SCHEME_PATH
    from gbdraw.config.modify import canonical_config_override_paths

    scheme_ids = " | ".join(_scalar(scheme) for scheme in CONSERVATION_COLOR_SCHEMES)
    packaged = apply_config_overrides(None, None).objects.conservation
    max_color_path = "objects.conservation.max_color"
    if max_color_path not in canonical_config_override_paths():
        raise SystemExit(f"{max_color_path} is no config override path")
    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/analysis/conservation.py\n"
        "// (CONSERVATION_COLOR_SCHEMES), the packaged config (objects.conservation.max_color) and\n"
        "// gbdraw/config/models/objects.py (CONSERVATION_COLOR_SCHEME_PATH); do not edit manually.\n"
        "\n"
        f"/** @typedef {{{scheme_ids}}} ConservationColorScheme A similarity ring colour scheme (Ring colors). */\n"
        "/**\n"
        " * The colours of each scheme, by drawn position, repeating: a ring without a colour of its own\n"
        " * takes `cycle[schemeIndex % cycle.length]`; `same` has none (the ring draws the configured colour).\n"
        " */\n"
        "/** @type {Readonly<Record<ConservationColorScheme, readonly string[]>>} */\n"
        f"export const CONSERVATION_COLOR_SCHEMES = {_frozen({key: list(value) for key, value in CONSERVATION_COLOR_SCHEMES.items()})};\n"
        "/** The scheme a drawing takes when it names none: the packaged default (`objects.conservation.color_scheme`). */\n"
        "/** @type {ConservationColorScheme} */\n"
        f"export const CONSERVATION_DEFAULT_COLOR_SCHEME = {_scalar(packaged.color_scheme)};\n"
        "/** The colour a ring without one draws under `same`: the packaged default (`objects.conservation.max_color`). */\n"
        f"export const CONSERVATION_SAME_COLOR = {_scalar(packaged.max_color)};\n"
        "/** The config path of the scheme a CLI or Python Session stores. */\n"
        f"export const CONSERVATION_COLOR_SCHEME_PATH = {_scalar(CONSERVATION_COLOR_SCHEME_PATH)};\n"
        "/** The config path of the colour a ring without one draws under `same`, which a drawing may keep preserved. */\n"
        f"export const CONSERVATION_MAX_COLOR_PATH = {_scalar(max_color_path)};\n"
    )


def render_sequence_search_programs_module() -> str:
    """The programs of Sequence search, in the order a search runs them."""
    from gbdraw.annotations import search_hits

    programs = " | ".join(_scalar(program) for program in search_hits.SEARCH_PROGRAMS)
    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/annotations/search_hits.py; do not edit manually.\n"
        "\n"
        f"/** @typedef {{{programs}}} SearchProgram The program that searched a query. */\n"
        "/** The programs in the order a search runs them: the order of their rows in the raw table and its name. */\n"
        "/** @type {readonly SearchProgram[]} */\n"
        f"export const SEARCH_PROGRAMS = {_frozen(list(search_hits.SEARCH_PROGRAMS))};\n"
    )


def render_sequence_search_colors_module() -> str:
    """The default hit colours of Sequence search and the Legend caption of its gradient."""
    from gbdraw.annotations import search_hits

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/annotations/search_hits.py; do not edit manually.\n"
        "\n"
        "/** The default hit colour of each query, in query order (repeated past the last). */\n"
        "/** @type {readonly string[]} */\n"
        f"export const SEARCH_QUERY_COLORS = {_frozen(list(search_hits.SEARCH_QUERY_COLORS))};\n"
        "/** The gradient by query position: its colour at the query's start and at its end. */\n"
        f"export const SEARCH_GRADIENT_START_COLOR = {_scalar(search_hits.SEARCH_GRADIENT_START_COLOR)};\n"
        f"export const SEARCH_GRADIENT_END_COLOR = {_scalar(search_hits.SEARCH_GRADIENT_END_COLOR)};\n"
        "/** The Legend caption of every gradient hit. */\n"
        f"export const QUERY_POSITION_LEGEND_LABEL = {_scalar(search_hits.QUERY_POSITION_LEGEND_LABEL)};\n"
    )


def render_result_records_module() -> str:
    """The attribute names of the records a Result keeps of what Python drew."""
    from gbdraw.svg import result_records

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/svg/result_records.py; do not edit manually.\n"
        "\n"
        "/** The prefix of the record of an attribute's drawn value: `data-gbdraw-base-<attribute>`. */\n"
        f"export const RESULT_BASE_ATTRIBUTE_PREFIX = {_scalar(result_records.RESULT_BASE_ATTRIBUTE_PREFIX)};\n"
        "/** The attributes whose drawn value a Result records. */\n"
        "/** @type {readonly string[]} */\n"
        f"export const RESULT_BASE_ATTRIBUTES = {_frozen(list(result_records.RESULT_BASE_ATTRIBUTES))};\n"
        "/** The record of Python's row order on a feature Legend group (its keys, as JSON). */\n"
        f"export const LEGEND_ORDER_RECORD = {_scalar(result_records.LEGEND_ORDER_RECORD)};\n"
        "/** The mark of a hidden Legend row that an earlier row of its name and swatch stands for: that row's key. */\n"
        f"export const LEGEND_MERGED_ROW = {_scalar(result_records.LEGEND_MERGED_ROW)};\n"
    )


def render_drawn_attributes_module() -> str:
    """The attribute names Python writes on drawn groups that the Web reads in a Result."""
    from gbdraw.legend.metrics import GRADIENT_SCALE_NAME_ATTRIBUTES
    from gbdraw.render.drawers.annotation_marks import ANNOTATION_ID_ATTRIBUTE, ANNOTATION_SET_ID_ATTRIBUTE
    from gbdraw.render.drawers.circular.conservation import (
        HIT_DISPLAY_SPAN_ATTRIBUTES,
        RING_LEGEND_NAME_ATTRIBUTE,
        RING_RADIUS_ATTRIBUTES,
    )
    from gbdraw.render.drawers.linear.wrap_rows import WRAP_ROW_ATTRIBUTE
    from gbdraw.render.interactive_svg import (
        MATCH_FRAGMENT_OWN_ATTRIBUTES,
        STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT,
        STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT,
    )
    from gbdraw.web_support.feature_catalog import MATCH_SPAN_REASON_TEXTS

    minimum, maximum = GRADIENT_SCALE_NAME_ATTRIBUTES
    span_start, span_end = HIT_DISPLAY_SPAN_ATTRIBUTES
    inner_radius, outer_radius = RING_RADIUS_ATTRIBUTES
    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/render/drawers/linear/wrap_rows.py\n"
        "// (WRAP_ROW_ATTRIBUTE), gbdraw/legend/metrics.py (GRADIENT_SCALE_NAME_ATTRIBUTES) and\n"
        "// gbdraw/render/drawers/annotation_marks.py (the annotation attributes) and\n"
        "// gbdraw/render/drawers/circular/conservation.py (the similarity-ring attributes) and\n"
        "// gbdraw/render/interactive_svg.py (the attributes of a match fragment and the slots of the standalone\n"
        "// script) and gbdraw/web_support/feature_catalog.py (MATCH_SPAN_REASON_TEXTS); do not edit manually.\n"
        "\n"
        "/** A drawn region annotation's group names its annotation and its set, in both modes. */\n"
        f"export const ANNOTATION_ID_ATTRIBUTE = {_scalar(ANNOTATION_ID_ATTRIBUTE)};\n"
        f"export const ANNOTATION_SET_ID_ATTRIBUTE = {_scalar(ANNOTATION_SET_ID_ATTRIBUTE)};\n"
        "/** A Linear row of a wrapped record is a group of its own; an annotation's marks on that row are its children. */\n"
        f"export const WRAP_ROW_ATTRIBUTE = {_scalar(WRAP_ROW_ATTRIBUTE)};\n"
        "/** The names of the ends of a Legend gradient row's scale, both or neither (an identity row has none). */\n"
        f"export const GRADIENT_SCALE_NAME_ATTRIBUTES = Object.freeze({{ min: {_scalar(minimum)}, max: {_scalar(maximum)} }});\n"
        "/** A similarity-ring hit's drawn span, in drawn bases: the display index of its first base and of the base after its last. */\n"
        f"export const HIT_DISPLAY_SPAN_ATTRIBUTES = Object.freeze({{ start: {_scalar(span_start)}, end: {_scalar(span_end)} }});\n"
        "/** A similarity ring's drawn radii (px, in its record's frame). */\n"
        f"export const RING_RADIUS_ATTRIBUTES = Object.freeze({{ inner: {_scalar(inner_radius)}, outer: {_scalar(outer_radius)} }});\n"
        "/** A similarity ring's name as its Legend row draws it, written also when no Legend is drawn. */\n"
        f"export const RING_LEGEND_NAME_ATTRIBUTE = {_scalar(RING_LEGEND_NAME_ATTRIBUTE)};\n"
        "/** The attributes a fragment of one comparison match holds on its own; its other data- attributes agree. */\n"
        "/** @type {readonly string[]} */\n"
        f"export const MATCH_FRAGMENT_OWN_ATTRIBUTES = {_frozen(list(MATCH_FRAGMENT_OWN_ATTRIBUTES))};\n"
        "/** Where the standalone interactive script's asset takes MATCH_FRAGMENT_OWN_ATTRIBUTES. */\n"
        f"export const STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT = {_scalar(STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT)};\n"
        "/** Where the standalone interactive script's asset takes MATCH_SPAN_REASON_TEXTS. */\n"
        f"export const STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT = {_scalar(STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT)};\n"
        "/**\n"
        " * The English sentence of each reason a match span cannot be read: the match popup shows it translated (a\n"
        " * catalog key) and an exported Interactive SVG shows it as it is. Not a persisted contract: a wording fix\n"
        " * changes no saved file.\n"
        " */\n"
        + _declaration("MATCH_SPAN_REASON_TEXTS", MATCH_SPAN_REASON_TEXTS)
    )


def render_label_fit_words_module() -> str:
    """The English words of a label fit row that the Web UI shows translated: catalog keys."""
    from gbdraw.labels.fit import ROW_TITLES
    from gbdraw.web_support.label_fit import label_fit_setting_words

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/labels/fit.py (ROW_TITLES)\n"
        "// and gbdraw/web_support/label_fit.py (label_fit_setting_words); do not edit manually.\n"
        "\n"
        "/** The words Python sends for a label fit row, which the dialog shows with `tm()`: catalog keys. */\n"
        "/** @type {readonly string[]} */\n"
        f"export const LABEL_FIT_WORDS = {_frozen([*ROW_TITLES.values(), *label_fit_setting_words()])};\n"
    )


def render_collinearity_words_module() -> str:
    """The block kinds and color modes Python writes on a match, which the popup shows in the UI language."""
    from gbdraw.analysis.collinearity import CollinearityBlockKind, CollinearityColorMode

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from gbdraw/analysis/collinearity.py\n"
        "// (CollinearityBlockKind, CollinearityColorMode); do not edit manually.\n"
        "\n"
        "/** The block kinds Python writes on a match, which the popup shows with `tmx()`: catalog keys. */\n"
        "/** @type {readonly string[]} */\n"
        f"export const COLLINEARITY_BLOCK_KINDS = {_frozen(list(get_args(CollinearityBlockKind)))};\n"
        "\n"
        "/** The color modes Python writes on a match; a test holds the Web's label table to them. */\n"
        "/** @type {readonly string[]} */\n"
        f"export const COLLINEARITY_COLOR_MODES = {_frozen(list(get_args(CollinearityColorMode)))};\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the generated contract files are out of sync.",
    )
    args = parser.parse_args()

    modules = {
        TARGET: render_web_contracts_module(),
        REQUEST_TARGET: render_request_typedefs_module(),
        SESSION_TARGET: render_session_typedefs_module(),
        REQUEST_SHAPE_TARGET: render_request_shape_module(),
        SESSION_SHAPE_TARGET: render_session_shape_module(),
        RULE_SYNTHESIS_TARGET: render_rule_synthesis_typedefs_module(),
        CARD_WINDOWS_TARGET: render_card_windows_typedefs_module(),
        LOCATE_TARGET: render_locate_typedefs_module(),
        DERIVED_RESULT_SHAPE_TARGET: render_derived_result_shape_module(),
        BUILTIN_COLOR_TABLES_TARGET: render_builtin_color_tables_module(),
        CONSERVATION_COLORS_TARGET: render_conservation_colors_module(),
        SEQUENCE_SEARCH_COLORS_TARGET: render_sequence_search_colors_module(),
        SEQUENCE_SEARCH_PROGRAMS_TARGET: render_sequence_search_programs_module(),
        RESULT_RECORDS_TARGET: render_result_records_module(),
        DRAWN_ATTRIBUTES_TARGET: render_drawn_attributes_module(),
        LABEL_FIT_WORDS_TARGET: render_label_fit_words_module(),
        COLLINEARITY_WORDS_TARGET: render_collinearity_words_module(),
        SESSION_JSON_SCHEMA_TARGET: render_session_json_schema(),
        REQUEST_JSON_SCHEMA_TARGET: render_request_json_schema(),
        SCHEMA_README_TARGET: render_schema_readme(),
    }
    if args.check:
        stale = [
            target.name
            for target, expected in modules.items()
            if not target.is_file() or target.read_text(encoding="utf-8") != expected
        ]
        if stale:
            print(
                f"Generated contract files are stale ({', '.join(stale)}). "
                "Run: python tools/generate_web_contracts.py"
            )
            return 1
        return 0

    for target, expected in modules.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(expected, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
