"""Neutral descriptions of persisted JSON contracts.

A shape is plain JSON-compatible data, so a generator for another language
renders it without importing the Python types it came from:

- ``{"kind": "string" | "integer" | "number" | "boolean" | "null" | "unknown"}``
- ``{"kind": "literal", "values": [...]}``: one of the listed JSON scalars
- ``{"kind": "array", "items": shape}`` and ``{"kind": "tuple", "items": [shape, ...]}``
- ``{"kind": "record", "values": shape}``: an object with any string keys
- ``{"kind": "union", "options": [shape, ...]}``: a nullable value is a union with null
- ``{"kind": "object", "fields": {name: shape}, "optional": [name, ...]}``: no other field
- ``{"kind": "ref", "name": name}``: a shape named in the description's ``types``
"""

from __future__ import annotations

from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from enum import Enum
import types
from typing import Any, Callable, Literal, Mapping, Union, get_args, get_origin

Shape = dict[str, Any]

STRING: Shape = {"kind": "string"}
INTEGER: Shape = {"kind": "integer"}
NUMBER: Shape = {"kind": "number"}
BOOLEAN: Shape = {"kind": "boolean"}
NULL: Shape = {"kind": "null"}
UNKNOWN: Shape = {"kind": "unknown"}
_SCALARS: dict[object, Shape] = {
    str: STRING,
    int: INTEGER,
    float: NUMBER,
    bool: BOOLEAN,
    type(None): NULL,
}


def literal(*values: object) -> Shape:
    """One of these JSON scalars."""

    for value in values:
        if value is not None and not isinstance(value, (str, int, bool)):
            raise TypeError(f"A literal must be a string, an integer, a boolean or null: {value!r}")
    return {"kind": "literal", "values": list(values)}


def array(items: Shape) -> Shape:
    return {"kind": "array", "items": items}


def record(values: Shape) -> Shape:
    return {"kind": "record", "values": values}


def tuple_of(*items: Shape) -> Shape:
    """An array of exactly these items; no item is the empty array."""

    return {"kind": "tuple", "items": list(items)}


def ref(name: str) -> Shape:
    return {"kind": "ref", "name": name}


def union(*options: Shape) -> Shape:
    """One of ``options``; nested unions are flattened and repeats dropped."""

    flat: list[Shape] = []
    for option in options:
        for item in option["options"] if option["kind"] == "union" else (option,):
            if item not in flat:
                flat.append(item)
    return flat[0] if len(flat) == 1 else {"kind": "union", "options": flat}


def nullable(shape: Shape) -> Shape:
    return union(shape, NULL)


def obj(fields: Mapping[str, Shape], *, optional: SequenceABC[str] = ()) -> Shape:
    unknown = set(optional) - set(fields)
    if unknown:
        raise ValueError(f"Optional field(s) not in the object: {sorted(unknown)}")
    return {"kind": "object", "fields": dict(fields), "optional": sorted(optional)}


def shape_of_hint(
    hint: object,
    special: Callable[[object], Shape | None] = lambda _hint: None,
) -> Shape:
    """The JSON shape of a value of type ``hint``.

    ``special`` names the shape of the classes a contract writes its own way
    (a dataclass, a table, a path); a hint with no rule raises.
    """

    known = special(hint)
    if known is not None:
        return known
    if hint is Any or hint is object:
        return UNKNOWN
    if hint in _SCALARS:
        return _SCALARS[hint]
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        return union(*(shape_of_hint(branch, special) for branch in args))
    if origin is Literal:
        return literal(*args)
    if origin is tuple and args and Ellipsis not in args:
        return tuple_of(*(shape_of_hint(item, special) for item in args))
    if origin in {tuple, list, SequenceABC}:
        return array(shape_of_hint(args[0], special) if args else UNKNOWN)
    if origin in {dict, MappingABC} or hint is dict:
        return record(shape_of_hint(args[1], special) if len(args) == 2 else UNKNOWN)
    if isinstance(hint, type) and issubclass(hint, Enum):
        return literal(*(member.value for member in hint))
    raise TypeError(f"No JSON shape is defined for the type {hint!r}.")
