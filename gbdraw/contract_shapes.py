"""Neutral descriptions of persisted JSON contracts.

A shape is plain JSON-compatible data, so a generator for another language
renders it without importing the Python types it came from:

- ``{"kind": "string" | "integer" | "number" | "boolean" | "null" | "unknown"}``
- ``{"kind": "literal", "values": [...]}``: one of the listed JSON scalars
- ``{"kind": "array", "items": shape}`` and ``{"kind": "tuple", "items": [shape, ...]}``
- ``{"kind": "record", "values": shape}``: an object with any string keys
- ``{"kind": "union", "options": [shape, ...]}``: a nullable value is a union with null
- ``{"kind": "object", "fields": {name: shape}, "optional": [name, ...], "open": bool}``:
  no other field, unless it is open; ``"defaults": {name: value}``, present when
  not empty, gives the JSON value a reader uses for a field that is absent (an
  optional one) or that it cannot read (any one)
- ``{"kind": "ref", "name": name}``: a shape named in the description's ``types``

A value domain that Python checks beyond the kind is written on the shape
(:func:`bounded`): ``minimum``, ``exclusiveMinimum`` and ``maximum`` on a
number or an integer, ``minItems`` on an array, and ``pattern`` (a regular
expression that must match somewhere in the text) on a string. A
literal may name the other spellings its reader takes as one of its values
(``spellings``: ``aliases``, ``foldHyphens``, ``blank``; empty, it reads case
and outer spaces only); the check does not read them.
"""

from __future__ import annotations

from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from enum import Enum
import json
import re
import types
from typing import Any, Callable, Literal, Mapping, NamedTuple, TypedDict, Union, get_args, get_origin, get_type_hints

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


class Spellings(TypedDict, total=False):
    """The other spellings of a name that its reader takes (:func:`read_spelling`)."""

    aliases: Mapping[str, str]  # a spelling, after case, spaces and hyphens are read, to its name
    foldHyphens: bool  # "-" reads as "_"
    blank: str  # the name an empty or unset value reads as


def read_spelling(value: object, names: SequenceABC[str], spellings: Spellings) -> str | None:
    """The name of ``names`` that ``value`` spells, or None.

    The text of the value without case and outer spaces, with ``-`` read as
    ``_`` when ``foldHyphens``; an alias reads as its name; an unset value
    (falsy) reads as ``blank``. The Web reads a stored request the same way
    (``readSpelling`` in gbdraw/web/js/utils/contract-shape.js).
    """

    text = str(value or spellings.get("blank", "")).strip().lower()
    if spellings.get("foldHyphens", False):
        text = text.replace("-", "_")
    text = spellings.get("aliases", {}).get(text, text)
    return text if text in names else None


def literal(*values: object, spellings: Spellings | None = None) -> Shape:
    """One of these JSON scalars.

    ``spellings`` names the other spellings of a value that its reader takes
    (:func:`read_spelling`); a reader that takes one writes the value.
    """

    for value in values:
        if value is not None and not isinstance(value, (str, int, bool)):
            raise TypeError(f"A literal must be a string, an integer, a boolean or null: {value!r}")
    shape: Shape = {"kind": "literal", "values": list(values)}
    if spellings is not None:
        shape["spellings"] = json.loads(json.dumps(spellings, allow_nan=False, sort_keys=True))
    return shape


_BOUNDS: dict[str, frozenset[str]] = {
    "minimum": frozenset({"integer", "number"}),
    "exclusiveMinimum": frozenset({"integer", "number"}),
    "maximum": frozenset({"integer", "number"}),
    "minItems": frozenset({"array"}),
    "pattern": frozenset({"string"}),
}


def bounded(shape: Shape, **bounds: object) -> Shape:
    """``shape`` with the value domain its owner checks.

    Each bound goes on the options of a union it applies to (``minItems`` on
    the array, ``minimum`` on the number); a bound that applies to none raises.
    """

    unknown = set(bounds) - set(_BOUNDS)
    if unknown:
        raise ValueError(f"Unknown bound(s): {sorted(unknown)}")
    options = shape["options"] if shape["kind"] == "union" else [shape]
    used: set[str] = set()
    result = []
    for option in options:
        own = {name: value for name, value in bounds.items() if option["kind"] in _BOUNDS[name]}
        used |= set(own)
        result.append({**option, **own})
    if used != set(bounds):
        raise ValueError(f"Bound(s) that apply to no option: {sorted(set(bounds) - used)}")
    return {**shape, "options": result} if shape["kind"] == "union" else result[0]


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


def non_null(shape: Shape) -> Shape:
    """``shape`` without its null alternative."""

    return union(*(
        item
        for item in (shape["options"] if shape["kind"] == "union" else (shape,))
        if item != NULL
    ))


def obj(
    fields: Mapping[str, Shape],
    *,
    optional: SequenceABC[str] = (),
    open: bool = False,
    defaults: Mapping[str, object] | None = None,
) -> Shape:
    """An object with these fields; an ``open`` one may hold other keys, of any value.

    ``defaults`` maps a field to the JSON value its reader uses in its place:
    when an optional field is absent, and when a stored value is unreadable
    (a required field's default is the writer's empty form of its part, or
    the repair target of a retired value: ``collinearityAnchorMode`` takes
    ``rbh`` for every unreadable stored value, ``one_to_one`` among them). A
    field whose absence means "unset" has none.
    """

    unknown = set(optional) - set(fields)
    if unknown:
        raise ValueError(f"Optional field(s) not in the object: {sorted(unknown)}")
    undeclared = set(defaults or ()) - set(fields)
    if undeclared:
        raise ValueError(f"Default(s) of field(s) not in the object: {sorted(undeclared)}")
    shape: Shape = {"kind": "object", "fields": dict(fields), "optional": sorted(optional), "open": open}
    if defaults:
        shape["defaults"] = {name: json.loads(json.dumps(defaults[name], allow_nan=False)) for name in sorted(defaults)}
    return shape


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


def typed_dict_shape(
    cls: type,
    special: Callable[[object], Shape | None] = lambda _hint: None,
    *,
    open: bool = False,
) -> Shape:
    """The object a ``TypedDict`` class declares: its keys are the wire names."""

    optional: frozenset[str] = getattr(cls, "__optional_keys__")
    return obj(
        {name: shape_of_hint(hint, special) for name, hint in get_type_hints(cls).items()},
        optional=sorted(optional),
        open=open,
    )


_SCALAR_CHECKS: dict[str, Callable[[object], bool]] = {
    "string": lambda value: isinstance(value, str),
    "boolean": lambda value: isinstance(value, bool),
    "integer": lambda value: isinstance(value, int) and not isinstance(value, bool),
    "number": lambda value: isinstance(value, (int, float)) and not isinstance(value, bool),
    "null": lambda value: value is None,
}

# Why a value leaves its shape; the Web checker (gbdraw/web/js/utils/contract-shape.js)
# reports the same codes, and tests/fixtures/contract_shape_vectors.json lists them.
ViolationReason = Literal[
    "unknown-key",  # a closed object holds a key it does not declare
    "missing-key",  # an object lacks a field that is not optional
    "wrong-kind",  # another JSON kind than the shape's
    "null-not-declared",  # null where the shape does not allow it
    "non-integer",  # a number with a fraction where an integer is declared
    "not-in-set",  # a scalar that is none of a literal's values
    "tuple-length",  # an array with another length than its tuple
    "out-of-range",  # a number outside its minimum, exclusiveMinimum or maximum
    "too-short",  # an array with fewer items than minItems
    "pattern",  # a text its pattern does not match
]


class Violation(NamedTuple):
    """Where a value leaves its shape and why; never the value, which can be private text."""

    path: str
    reason: ViolationReason
    message: str


def _kind_violation(value: object, kind: str, path: str) -> Violation:
    """Why ``value`` is not of a scalar, array or object ``kind``."""

    article = {"array": "an array", "tuple": "an array", "object": "an object", "record": "an object"}
    message = f"expected {article.get(kind, kind)}"
    if value is None:
        return Violation(path, "null-not-declared", message)
    if kind == "integer" and isinstance(value, float):
        return Violation(path, "non-integer", message)
    return Violation(path, "wrong-kind", message)


def _domain_violations(value: Any, shape: Mapping[str, Any], path: str) -> list[Violation]:
    """Where a value of the shape's kind leaves the domain its bounds name."""

    if "minimum" in shape and value < shape["minimum"]:
        return [Violation(path, "out-of-range", f"expected at least {shape['minimum']}")]
    if "exclusiveMinimum" in shape and value <= shape["exclusiveMinimum"]:
        return [Violation(path, "out-of-range", f"expected more than {shape['exclusiveMinimum']}")]
    if "maximum" in shape and value > shape["maximum"]:
        return [Violation(path, "out-of-range", f"expected at most {shape['maximum']}")]
    if "minItems" in shape and len(value) < shape["minItems"]:
        return [Violation(path, "too-short", f"expected at least {shape['minItems']} items")]
    if "pattern" in shape and re.search(shape["pattern"], value) is None:
        return [Violation(path, "pattern", f"expected text matching {shape['pattern']}")]
    return []


def shape_violations(
    value: object, shape: Mapping[str, Any], types: Mapping[str, Mapping[str, Any]], path: str
) -> list[Violation]:
    """Where a parsed JSON ``value`` leaves ``shape``, in order.

    ``types`` holds the named shapes a ``ref`` points to. Of a union's
    options, the first one with the fewest violations is reported.
    """

    kind = shape["kind"]
    if kind == "unknown":
        return []
    if kind == "ref":
        return shape_violations(value, types[shape["name"]], types, path)
    if kind in _SCALAR_CHECKS:
        if not _SCALAR_CHECKS[kind](value):
            return [_kind_violation(value, kind, path)]
        return _domain_violations(value, shape, path)
    if kind == "literal":
        if any(value == item and type(value) is type(item) for item in shape["values"]):
            return []
        reason: ViolationReason = "null-not-declared" if value is None else "not-in-set"
        return [Violation(path, reason, f"expected one of {json.dumps(shape['values'], separators=(',', ':'))}")]
    if kind == "union":
        found = [shape_violations(value, option, types, path) for option in shape["options"]]
        return min(found, key=len) if all(found) else []
    if kind in {"array", "tuple"}:
        if not isinstance(value, list):
            return [_kind_violation(value, kind, path)]
        if kind == "tuple" and len(value) != len(shape["items"]):
            return [Violation(path, "tuple-length", f"expected {len(shape['items'])} items")]
        items = shape["items"] if kind == "tuple" else [shape["items"]] * len(value)
        return _domain_violations(value, shape, path) + [
            problem
            for index, (item, item_shape) in enumerate(zip(value, items))
            for problem in shape_violations(item, item_shape, types, f"{path}[{index}]")
        ]
    if kind not in {"record", "object"}:
        raise ValueError(f"Unknown shape kind: {kind!r}")
    if not isinstance(value, dict):
        return [_kind_violation(value, kind, path)]
    if kind == "record":
        return [
            problem
            for key, item in value.items()
            for problem in shape_violations(item, shape["values"], types, f"{path}.{key}")
        ]
    fields = shape["fields"]
    problems = [] if shape["open"] else [
        Violation(f"{path}.{key}", "unknown-key", "unknown field") for key in value if key not in fields
    ]
    problems += [
        Violation(f"{path}.{key}", "missing-key", "required field is absent")
        for key in fields
        if key not in value and key not in shape["optional"]
    ]
    for key, item in value.items():
        if key in fields:
            problems += shape_violations(item, fields[key], types, f"{path}.{key}")
    return problems


def violations(
    value: object, shape: Mapping[str, Any], types: Mapping[str, Mapping[str, Any]], path: str
) -> list[str]:
    """:func:`shape_violations` as ``"<path>: <message>"``."""

    return [f"{problem.path}: {problem.message}" for problem in shape_violations(value, shape, types, path)]
