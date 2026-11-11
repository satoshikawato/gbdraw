"""The published JSON Schema files agree with the structural checker of gbdraw.contract_shapes.

`tools/generate_web_contracts.py` writes docs/REFERENCE/schemas/ (the two
schemas and their page of model diagrams) from the descriptions; a standard
validator (`jsonschema`, the `dev` extra) reads the schemas here.
"""

from __future__ import annotations

import functools
import json
from typing import Any, cast

import pytest

from gbdraw.contract_shapes import shape_violations
from gbdraw.session_contract import describe_session_admission
from gbdraw.session_io import validate_session
from gbdraw.session_request_codec import _as_decoded_settings
from gbdraw.exceptions import ValidationError
from tests.test_canonical_request_shape import _stored_requests
from tests.test_contract_shape_vectors import VECTORS as SHAPE_VECTORS
from tests.test_session_admission_vectors import VECTORS as ADMISSION_VECTORS, _edited
from tests.test_web_contracts_generated import _generator
from tests.utils.contract_shape_check import REPO_ROOT, load_json, stored_json_files

# The published documentation; documentation-only pull requests run these tests.
pytestmark = pytest.mark.recipe

ADMISSION = describe_session_admission()
BAND_ROOT = {version: band["root"] for band in ADMISSION["bands"] for version in band["versions"]}
# JSON Schema reads 1.0 as the integer 1, as JSON.parse does; Python's json reads a float.
READ_AS_THE_WEB_READS = {"integral-fraction-at-an-integer"}
# A JSON Schema number is any JSON number text, 1e400 too. Python's json and JSON.parse read that text as
# infinity, which the checkers refuse (no number field takes one; OV-3237); the published schema cannot say so.
SCHEMA_TAKES_ANY_NUMBER_TEXT = {"number-that-overflows-at-a-number"}
# Sessions that hold a retired value, which stored_json_files() leaves out.
RETIRED_SESSIONS = sorted((REPO_ROOT / "tests" / "retired_sessions").glob("*.gbdraw-session.json*"))


def _as_schema_10(value: object) -> object:
    """A stored schema 9 request as the decoder reads it: each record's ``region`` is one of its ``windows``."""

    if isinstance(value, dict):
        if value.get("schema") == 9 and "diagramOptions" in value and isinstance(value.get("records"), list):
            value = {**value, "schema": 10, "records": [
                {
                    **{key: item for key, item in record.items() if key != "region"},
                    "windows": [] if record.get("region") is None else [{
                        **{key: record["region"][key] for key in ("start", "end", "reverseComplement")},
                        "column": None,
                    }],
                }
                for record in value["records"]
            ]}
        return {key: _as_schema_10(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_as_schema_10(item) for item in value]
    return value


def _sessions() -> list[tuple[str, dict[str, Any]]]:
    found = []
    for path in [*stored_json_files(), *RETIRED_SESSIONS]:
        document = load_json(path)
        if isinstance(document, dict) and document.get("format") == "gbdraw-session":
            # A retired Session stores the schema 9 request of its time (its region selectors are null).
            found.append((
                str(path.relative_to(REPO_ROOT)),
                cast(dict[str, Any], _as_schema_10(document)) if path in RETIRED_SESSIONS else document,
            ))
    return found


def _validator(schema: dict[str, Any]) -> Any:
    import jsonschema

    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(schema)


def _published(target: str) -> dict[str, Any]:
    return json.loads(getattr(_generator(), target).read_text(encoding="utf-8"))


@functools.cache
def _published_validator(target: str) -> Any:
    return _validator(_published(target))


def _with_root(document: dict[str, Any], name: str) -> dict[str, Any]:
    """The schema document whose root is its definition ``name``."""
    defs = document["$defs"]
    return {"$schema": document["$schema"], "$ref": f"#/$defs/{name}", "$defs": defs}


def _shape_schema(shape: dict[str, Any]) -> dict[str, Any]:
    generator = _generator()
    published = generator.ecma_types({**SHAPE_VECTORS["types"], "value": shape})
    defs = {name: generator.json_schema(item) for name, item in published.items()}
    return {"$schema": generator.JSON_SCHEMA_DIALECT, "$ref": "#/$defs/value", "$defs": defs}


def _in_band(document: object) -> bool:
    """Whether the structural checker places ``document`` inside the band of its version."""
    if not isinstance(document, dict) or document.get("format") != "gbdraw-session":
        return False
    root = BAND_ROOT.get(document.get("version"))
    return root is not None and not shape_violations(document, ADMISSION["types"][root], ADMISSION["types"], "")


def test_the_published_schemas_and_diagrams_are_current() -> None:
    generator = _generator()
    published = {
        generator.SESSION_JSON_SCHEMA_TARGET: generator.render_session_json_schema(),
        generator.REQUEST_JSON_SCHEMA_TARGET: generator.render_request_json_schema(),
        generator.SCHEMA_README_TARGET: generator.render_schema_readme(),
    }
    stale = [path.name for path, text in published.items() if path.read_text(encoding="utf-8") != text]

    assert stale == [], "Run: python tools/generate_web_contracts.py"


@pytest.mark.parametrize("vector", SHAPE_VECTORS["vectors"], ids=lambda vector: vector["id"])
def test_a_shape_vector_fits_its_schema_when_the_checker_reports_nothing(vector: dict[str, Any]) -> None:
    assert _validator(_shape_schema(vector["shape"])).is_valid(vector["value"]) == (vector["violations"] == [])


@pytest.mark.parametrize("difference", SHAPE_VECTORS["jsonDifferences"], ids=lambda difference: difference["id"])
def test_the_schema_reads_a_json_difference_as_python_does(difference: dict[str, Any]) -> None:
    reading = "web" if difference["id"] in READ_AS_THE_WEB_READS else "python"
    fits = _validator(_shape_schema(difference["shape"])).is_valid(json.loads(difference["text"]))

    assert fits == (difference[reading] == [] or difference["id"] in SCHEMA_TAKES_ANY_NUMBER_TEXT)


@pytest.mark.parametrize("vector", SHAPE_VECTORS["requestTypes"], ids=lambda vector: vector["id"])
def test_a_request_type_vector_fits_the_request_schema_when_the_checker_reports_nothing(
    vector: dict[str, Any],
) -> None:
    schema = _with_root(_published("REQUEST_JSON_SCHEMA_TARGET"), vector["type"])

    assert _validator(schema).is_valid(vector["value"]) == (vector["violations"] == [])


def test_a_stored_session_fits_the_session_schema_when_it_is_inside_its_band() -> None:
    validator = _published_validator("SESSION_JSON_SCHEMA_TARGET")
    disagreements = [name for name, document in _sessions() if validator.is_valid(document) != _in_band(document)]

    assert disagreements == []


def test_every_stored_session_python_accepts_fits_the_session_schema() -> None:
    validator = _published_validator("SESSION_JSON_SCHEMA_TARGET")
    accepted: list[str] = []
    for name, document in _sessions():
        try:
            validate_session(document)
        except ValidationError:
            continue
        accepted.append(name)

    assert [name for name, document in _sessions() if name in accepted and not validator.is_valid(document)] == []
    assert len(accepted) >= 50 and any(name.startswith("tests/retired_sessions/") for name in accepted)


@pytest.mark.parametrize("vector", ADMISSION_VECTORS, ids=lambda vector: vector["id"])
def test_an_admission_vector_fits_the_session_schema_when_it_is_inside_its_band(vector: dict[str, Any]) -> None:
    document = _edited(vector)

    assert _published_validator("SESSION_JSON_SCHEMA_TARGET").is_valid(document) == _in_band(document)


def test_every_stored_request_fits_the_request_schema_once_its_retired_anchor_modes_are_read() -> None:
    """The decoder reads a retired collinearityAnchorMode before it checks the
    description; that is the one stored form outside the schema that it reads."""
    validator = _published_validator("REQUEST_JSON_SCHEMA_TARGET")
    refused: list[str] = []
    repaired: list[str] = []
    modes: list[str] = []
    for name, document in _sessions():
        for request in _stored_requests(document):
            modes.append(request["mode"])
            if not validator.is_valid(_as_decoded_settings(request)):
                refused.append(name)
            elif not validator.is_valid(request):
                repaired.append(name)

    assert refused == []
    assert repaired and all(name.startswith("tests/retired_sessions/") for name in repaired)
    assert modes.count("circular") >= 4 and modes.count("linear") >= 4


def _patterns(value: object) -> list[str]:
    if isinstance(value, dict):
        own = [value["pattern"]] if isinstance(value.get("pattern"), str) else []
        return own + [pattern for item in value.values() for pattern in _patterns(item)]
    return [pattern for item in value for pattern in _patterns(item)] if isinstance(value, list) else []


def test_a_published_pattern_names_the_characters_python_reads_as_spaces() -> None:
    """ECMA-262, the pattern dialect of JSON Schema, reads \\s as other characters than Python."""
    patterns = [pattern for target in ("SESSION_JSON_SCHEMA_TARGET", "REQUEST_JSON_SCHEMA_TARGET")
                for pattern in _patterns(_published(target))]

    assert patterns and not any("\\s" in pattern or "\\S" in pattern for pattern in patterns)
    for unportable in ("[\\s,]", "^x$"):
        with pytest.raises(ValueError):
            _generator()._portable_pattern(unportable)
