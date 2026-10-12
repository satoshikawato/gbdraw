"""Check a JSON value against a description of ``gbdraw.contract_shapes``."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Mapping

REPO_ROOT = Path(__file__).resolve().parents[2]

_SCALAR_CHECKS = {
    "string": lambda value: isinstance(value, str),
    "boolean": lambda value: isinstance(value, bool),
    "integer": lambda value: isinstance(value, int) and not isinstance(value, bool),
    "number": lambda value: isinstance(value, (int, float)) and not isinstance(value, bool),
    "null": lambda value: value is None,
}


def violations(
    value: object, shape: Mapping[str, Any], types: Mapping[str, Mapping[str, Any]], path: str
) -> list[str]:
    """Where ``value`` leaves ``shape``: an unknown field, a missing one, or a wrong value.

    ``types`` holds the named shapes a ``ref`` points to.
    """

    kind = shape["kind"]
    if kind == "unknown":
        return []
    if kind == "ref":
        return violations(value, types[shape["name"]], types, path)
    if kind in _SCALAR_CHECKS:
        return [] if _SCALAR_CHECKS[kind](value) else [f"{path}: {value!r} is not {kind}"]
    if kind == "literal":
        allowed = any(value == item and type(value) is type(item) for item in shape["values"])
        return [] if allowed else [f"{path}: {value!r} is not one of {shape['values']!r}"]
    if kind == "union":
        found = [violations(value, option, types, path) for option in shape["options"]]
        return min(found, key=len) if all(found) else []
    if kind in {"array", "tuple"}:
        if not isinstance(value, list):
            return [f"{path}: {value!r} is not an array"]
        if kind == "tuple" and len(value) != len(shape["items"]):
            return [f"{path}: expected {len(shape['items'])} items"]
        items = shape["items"] if kind == "tuple" else [shape["items"]] * len(value)
        return [
            problem
            for index, (item, item_shape) in enumerate(zip(value, items))
            for problem in violations(item, item_shape, types, f"{path}[{index}]")
        ]
    if not isinstance(value, dict):
        return [f"{path}: {value!r} is not an object"]
    if kind == "record":
        return [
            problem
            for key, item in value.items()
            for problem in violations(item, shape["values"], types, f"{path}.{key}")
        ]
    assert kind == "object", kind
    fields = shape["fields"]
    problems = [] if shape["open"] else [f"{path}.{key}: unknown field" for key in value if key not in fields]
    problems += [
        f"{path}.{key}: required field is absent"
        for key in fields
        if key not in value and key not in shape["optional"]
    ]
    for key, item in value.items():
        if key in fields:
            problems += violations(item, fields[key], types, f"{path}.{key}")
    return problems


def stored_json_files() -> list[Path]:
    """The JSON files of the folders that hold Sessions and requests: fixtures,
    test inputs, documentation Sessions and the Gallery Sessions."""

    return sorted(
        path
        for folder in (
            "tests/fixtures", "tests/web/fixtures", "tests/test_inputs", "docs/images",
            "gbdraw/web/gallery/sessions",
        )
        for path in (REPO_ROOT / folder).rglob("*.json*")
        if path.suffix in {".json", ".gz"}
    )


def load_json(path: Path) -> object:
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
