"""JavaScript's reading of a JSON value: ``String()``, ``Number()``, truthiness, the integer test, and key order;
and the JSON readers on top of them: a copy of a value, an array read as objects, and a safe or non-negative integer.

Python readers of Web-written Session values use these so that they read a value as the Web does.
"""

from __future__ import annotations

import copy
import json
import math
import re
from typing import Any, Mapping

from .numeric import JS_MAX_SAFE_INTEGER, JS_WHITESPACE


# The Web readers' ``String(value ?? '')`` and ``.trim()`` (JS_WHITESPACE).
def js_string(value: object) -> str:
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
        return ",".join(js_string(item) for item in value)
    return "[object Object]"


def js_text(value: object) -> str:
    return js_string(value).strip(JS_WHITESPACE)


# A JavaScript array index: the canonical decimal text of the value.
JS_ARRAY_INDEX = re.compile(r"0|[1-9][0-9]*")


def js_truthy(value: object) -> bool:
    if value is None or isinstance(value, bool):
        return bool(value)
    if isinstance(value, (int, float)):
        return value != 0 and value == value
    if isinstance(value, str):
        return value != ""
    # Every JavaScript object and array is true, an empty one too.
    return True


def js_key_order(record: Mapping[str, object]) -> list[str]:
    """The keys of a JSON object in the order JavaScript lists them: array
    indexes ascending, then the other keys in insertion order."""

    def is_index(key: object) -> bool:
        return isinstance(key, str) and JS_ARRAY_INDEX.fullmatch(key) is not None and int(key) < 2**32 - 1

    indexes = sorted((key for key in record if is_index(key)), key=int)
    return indexes + [key for key in record if not is_index(key)]


def js_integer(value: object) -> bool:
    return (isinstance(value, int) and not isinstance(value, bool)) or (
        isinstance(value, float) and value.is_integer()
    )


def js_number(value: object) -> float:
    """``Number(value)`` for JSON values (NaN for what JavaScript cannot read)."""

    if value is None:
        return 0.0
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip(JS_WHITESPACE)
        if not text:
            return 0.0
        try:
            return float(text) if not re.search(r"[^0-9eE.+-]", text) else math.nan
        except ValueError:
            return math.nan
    if isinstance(value, list) and len(value) <= 1:
        return js_number(value[0]) if value else 0.0
    return math.nan


def js_integral(value: float) -> int | float:
    return int(value) if value.is_integer() else value


def safe_integer(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    if isinstance(value, int) and abs(value) <= JS_MAX_SAFE_INTEGER:
        return value
    return None


def mapping_list(value: object) -> list[Mapping[str, object]]:
    """The entries of a JSON array, each read as an object (a non-object reads as empty)."""

    return [item if isinstance(item, Mapping) else {} for item in value] if isinstance(value, list) else []


def nonnegative_integer(value: object) -> int | None:
    """``Number(value)`` when it is a non-negative safe integer (``null`` and ``''`` are not)."""

    if value is None or value == "":
        return None
    number = js_number(value)
    if math.isfinite(number) and number.is_integer() and 0 <= number <= JS_MAX_SAFE_INTEGER:
        return int(number)
    return None


def json_clone(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value))
    except (TypeError, ValueError):
        return copy.deepcopy(value)


__all__ = [
    "js_string",
    "js_text",
    "JS_ARRAY_INDEX",
    "js_truthy",
    "js_key_order",
    "js_integer",
    "js_number",
    "js_integral",
    "safe_integer",
    "mapping_list",
    "nonnegative_integer",
    "json_clone",
]
