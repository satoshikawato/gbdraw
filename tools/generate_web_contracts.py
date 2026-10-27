#!/usr/bin/env python3
"""Generate the browser's constants and typedefs of persisted contracts."""

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
    ("SUPPORTED_CANONICAL_REQUEST_SCHEMAS", "gbdraw.session_request_codec", "SUPPORTED_CANONICAL_REQUEST_SCHEMAS"),
    ("SIMILARITY_ALIGNMENT_PLAN_SCHEMA", "gbdraw.layout.similarity_alignment", "SIMILARITY_ALIGNMENT_PLAN_SCHEMA"),
    ("FEATURE_CATALOG_SCHEMA", "gbdraw.session_io", "CURRENT_FEATURE_CATALOG_SCHEMA"),
    ("ANCHOR_FEATURE_CATALOG_SCHEMA", "gbdraw.web_support.feature_catalog", "ANCHOR_FEATURE_CATALOG_SCHEMA"),
    ("LEGACY_FEATURE_CATALOG_SCHEMA", "gbdraw.web_support.feature_catalog", "LEGACY_FEATURE_CATALOG_SCHEMA"),
    ("FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION", "gbdraw.session_io", "FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION"),
    ("PROTEIN_LOSAT_CACHE_SCHEMA", "gbdraw.analysis.protein_colinearity", "PROTEIN_LOSAT_CACHE_SCHEMA"),
    ("LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA", "gbdraw.analysis.protein_colinearity", "LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA"),
    ("PROTEIN_IDENTITY_MANIFEST_SCHEMA", "gbdraw.analysis.protein_colinearity", "PROTEIN_IDENTITY_MANIFEST_SCHEMA"),
    ("ORTHOGROUP_INFERENCE_VERSION", "gbdraw.analysis.protein_colinearity", "ORTHOGROUP_INFERENCE_VERSION"),
    ("NUCLEOTIDE_LOSAT_CACHE_SCHEMA", "gbdraw.comparisons.losat_jobs", "NUCLEOTIDE_LOSAT_CACHE_SCHEMA"),
    ("LOSAT_DERIVED_CACHE_SCHEMA", "gbdraw.analysis.protein_artifacts", "CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA"),
    ("LEGACY_LOSAT_DERIVED_CACHE_SCHEMA", "gbdraw.session_io", "LEGACY_LOSAT_DERIVED_CACHE_SCHEMA"),
    ("LEGACY_PROTEIN_CANDIDATE_SCHEMA", "gbdraw.session_io", "LEGACY_PROTEIN_CANDIDATE_SCHEMA"),
    ("DEPTH_FILE_SCHEMA", "gbdraw.session_io", "DEPTH_FILE_SCHEMA"),
    ("DEPTH_FILE_ENCODING", "gbdraw.session_io", "DEPTH_FILE_ENCODING"),
    ("LEGEND_POSITIONS", "gbdraw.api.options", "LEGEND_POSITIONS"),
    # A "This feature only" color rule names a feature instance; the Source
    # recipe names the CLI record of its record key.
    ("FEATURE_INSTANCE_SEPARATOR", "gbdraw.features.ids", "FEATURE_INSTANCE_SEPARATOR"),
    ("CLI_RECORD_KEY_PREFIX", "gbdraw.api.record_planning", "CLI_RECORD_KEY_PREFIX"),
    ("EXPANDED_RECORD_KEY_SEPARATOR", "gbdraw.api.record_planning", "EXPANDED_RECORD_KEY_SEPARATOR"),
    ("CONFIG_UNSET_NUMBER_SPELLINGS", "gbdraw.config.models.objects", "UNSET_NUMBER_SPELLINGS"),
    # Session field sets: the keys of the classes of gbdraw/session_contract.py.
    ("FLAT_DRAFT_TOP_LEVEL_FIELDS", "gbdraw.session_io", "FLAT_DRAFT_TOP_LEVEL_FIELDS"),
    ("OTHER_MODE_RESULT_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_FIELDS"),
    ("OTHER_MODE_RESULT_UI_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_UI_FIELDS"),
    ("OTHER_MODE_RESULT_EDITOR_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_EDITOR_FIELDS"),
    ("OTHER_MODE_RESULT_LEGEND_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_LEGEND_FIELDS"),
    ("OTHER_MODE_RESULT_RUN_METADATA_FIELDS", "gbdraw.session_io", "OTHER_MODE_RESULT_RUN_METADATA_FIELDS"),
    ("CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS", "gbdraw.session_io", "CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS"),
    ("LEGACY_SVG_ATTRIBUTE_RENAMES", "gbdraw.render.interactive_svg", "LEGACY_SVG_ATTRIBUTE_RENAMES"),
    # A former palette name reads as the current one (Session Load).
    ("PALETTE_NAME_ALIASES", "gbdraw.io.colors", "PALETTE_NAME_ALIASES"),
)
# (Web name, Python module, function): values a Python function derives from
# its owner, written in this order after the constants.
_DERIVED: tuple[tuple[str, str, str], ...] = (
    ("EXPECTED_WEB_RUNTIME_CAPABILITIES", "gbdraw.web_support.capabilities", "get_web_runtime_capabilities"),
    ("CONFIG_UNSET_NUMBER_PATHS", "gbdraw.config.modify", "config_paths_reading_unset_spellings"),
    ("LEGACY_REQUEST_TABLES", "gbdraw.session_request_codec", "describe_legacy_request_tables"),
    # The Auto organism name rule's data; the Web fast path runs its steps (R4).
    ("ORGANISM_NAME_RULE", "gbdraw.core.record_metadata", "organism_name_rule_data"),
)
# (Web name, Python module, dataclass, field): the values of a field typed as a
# string Literal (typing.get_args, in the owner's order) and the field's
# default, written after the derived values. The Web shows these values in the
# UI language (app/pairwise-match-popup.js).
_LITERAL_FIELDS: tuple[tuple[str, str, str, str], ...] = (
    ("ORTHOGROUP_MEMBER_ROLE", "gbdraw.analysis.protein_colinearity", "OrthogroupMember", "role"),
    ("ORTHOGROUP_MEMBER_CONFIDENCE", "gbdraw.analysis.protein_colinearity", "OrthogroupMember", "confidence"),
)
_IDENTIFIER = re.compile(r"[A-Za-z_$][A-Za-z0-9_$]*|0|[1-9][0-9]*")


def _scalar(value: object) -> str:
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
        return "Object.freeze([" + ", ".join(_scalar(item) for item in value) + "])"
    return _scalar(value)


def _declaration(web_name: str, value: object) -> str:
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
    from gbdraw.web_support.rule_matching import describe_rule_synthesis

    return (
        "// @ts-check\n"
        "// Generated by tools/generate_web_contracts.py from describe_rule_synthesis()\n"
        "// in gbdraw/web_support/rule_matching.py; do not edit manually.\n"
        "// The diagram helper `synthesizeColorRules`: its request and its answer (not persisted).\n"
        "\n"
        + _typedefs(describe_rule_synthesis()["types"])
        + "export {};\n"
    )


def _references(shape: object) -> set[str]:
    """The names a shape refers to."""
    if isinstance(shape, dict):
        own = {shape["name"]} if shape.get("kind") == "ref" else set()
        return own.union(*(_references(item) for item in shape.values()))
    return set().union(*(_references(item) for item in shape)) if isinstance(shape, list) else set()


def render_session_typedefs_module() -> str:
    from gbdraw.session_contract import describe_session_admission, describe_session_document
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
    (``{"root", "types"}``) or ``BandedContractDescription`` (``{"bands", "types"}``)."""
    return (
        "// @ts-check\n"
        f"// Generated by tools/generate_web_contracts.py from {owner}; do not edit manually.\n"
        "// The description as data, for gbdraw/web/js/utils/contract-shape.js.\n"
        "\n"
        f"/** @import {{ {description_type} }} from '../utils/contract-shape.js' */\n"
        "\n"
        f"/** @type {{{description_type}}} */\n"
        f"export const {export_name} = {_frozen_json(json.loads(json.dumps(description)))};\n"
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the generated browser contract modules are out of sync.",
    )
    args = parser.parse_args()

    modules = {
        TARGET: render_web_contracts_module(),
        REQUEST_TARGET: render_request_typedefs_module(),
        SESSION_TARGET: render_session_typedefs_module(),
        REQUEST_SHAPE_TARGET: render_request_shape_module(),
        SESSION_SHAPE_TARGET: render_session_shape_module(),
        RULE_SYNTHESIS_TARGET: render_rule_synthesis_typedefs_module(),
    }
    if args.check:
        stale = [
            target.name
            for target, expected in modules.items()
            if not target.is_file() or target.read_text(encoding="utf-8") != expected
        ]
        if stale:
            print(
                f"Browser contract modules are stale ({', '.join(stale)}). "
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
