"""The published JSON Schema files agree with the structural checker of gbdraw.contract_shapes.

`tools/generate_web_contracts.py` writes docs/REFERENCE/schemas/ and the model
diagrams of docs/REFERENCE/session-and-request-compatibility.md from the same
descriptions; a standard validator (`jsonschema`, the `dev` extra) reads the
schemas here.
"""

from __future__ import annotations

import functools
import json
from typing import Any

import pytest

from gbdraw.contract_shapes import shape_violations
from gbdraw.session_contract import describe_session_admission
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
    defs = {name: generator.json_schema(item) for name, item in SHAPE_VECTORS["types"].items()}
    defs["value"] = generator.json_schema(shape)
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
        generator.CONTRACT_MODEL_PAGE: generator.render_contract_model_page(),
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

    assert fits == (difference[reading] == [])


@pytest.mark.parametrize("vector", SHAPE_VECTORS["requestTypes"], ids=lambda vector: vector["id"])
def test_a_request_type_vector_fits_the_request_schema_when_the_checker_reports_nothing(
    vector: dict[str, Any],
) -> None:
    schema = _with_root(_published("REQUEST_JSON_SCHEMA_TARGET"), vector["type"])

    assert _validator(schema).is_valid(vector["value"]) == (vector["violations"] == [])


def test_a_stored_session_fits_the_session_schema_when_it_is_inside_its_band() -> None:
    validator = _published_validator("SESSION_JSON_SCHEMA_TARGET")
    disagreements: list[str] = []
    fitting = 0
    for path in stored_json_files():
        document = load_json(path)
        if not isinstance(document, dict) or document.get("format") != "gbdraw-session":
            continue
        fits = validator.is_valid(document)
        fitting += fits
        if fits != _in_band(document):
            disagreements.append(str(path.relative_to(REPO_ROOT)))

    assert disagreements == []
    assert fitting >= 50


@pytest.mark.parametrize("vector", ADMISSION_VECTORS, ids=lambda vector: vector["id"])
def test_an_admission_vector_fits_the_session_schema_when_it_is_inside_its_band(vector: dict[str, Any]) -> None:
    document = _edited(vector)

    assert _published_validator("SESSION_JSON_SCHEMA_TARGET").is_valid(document) == _in_band(document)


def test_every_stored_request_of_the_current_schema_fits_the_request_schema() -> None:
    validator = _published_validator("REQUEST_JSON_SCHEMA_TARGET")
    refused: list[str] = []
    modes: list[str] = []
    for path in stored_json_files():
        for request in _stored_requests(load_json(path)):
            modes.append(request["mode"])
            if not validator.is_valid(request):
                refused.append(str(path.relative_to(REPO_ROOT)))

    assert refused == []
    assert modes.count("circular") >= 4 and modes.count("linear") >= 4
