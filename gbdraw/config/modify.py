#!/usr/bin/env python
# coding: utf-8

"""Schema-derived updates for the typed gbdraw configuration."""

from __future__ import annotations

from copy import deepcopy
from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from dataclasses import asdict, fields, is_dataclass
from functools import lru_cache
from numbers import Integral, Real
from types import UnionType
from typing import Any, Literal, Mapping, Union, get_args, get_origin, get_type_hints

from gbdraw.exceptions import ValidationError

from .models import GbdrawConfig  # type: ignore[reportMissingImports]


def _nested_dataclass_type(annotation: object) -> type | None:
    if isinstance(annotation, type) and is_dataclass(annotation):
        return annotation
    if get_origin(annotation) not in {Union, UnionType}:
        return None
    nested = {
        member
        for member in get_args(annotation)
        if isinstance(member, type) and is_dataclass(member)
    }
    non_dataclass = {
        member
        for member in get_args(annotation)
        if member is not type(None)
        and not (isinstance(member, type) and is_dataclass(member))
    }
    if len(nested) == 1 and not non_dataclass:
        return next(iter(nested))
    return None


@lru_cache(maxsize=1)
def _config_override_schema() -> tuple[dict[str, object], frozenset[str]]:
    """Return leaf annotations and non-leaf paths from ``GbdrawConfig``."""

    leaves: dict[str, object] = {}
    branches: set[str] = set()

    def visit(config_type: type, prefix: tuple[str, ...] = ()) -> None:
        annotations = get_type_hints(config_type)
        for config_field in fields(config_type):
            path_parts = (*prefix, config_field.name)
            path = ".".join(path_parts)
            annotation = annotations[config_field.name]
            nested_type = _nested_dataclass_type(annotation)
            if nested_type is None:
                leaves[path] = annotation
                continue
            branches.add(path)
            visit(nested_type, path_parts)

    visit(GbdrawConfig)
    return leaves, frozenset(branches)


def canonical_config_override_paths() -> frozenset[str]:
    """Return every canonical dotted leaf path accepted as an override."""

    leaves, _branches = _config_override_schema()
    return frozenset(leaves)


def config_to_raw_dict(config: GbdrawConfig) -> dict[str, Any]:
    """Return a lossless raw representation of a typed configuration."""

    if not isinstance(config, GbdrawConfig):
        raise ValidationError("config must be GbdrawConfig")
    raw = asdict(config)
    # ``LabelsFilteringConfig.raw`` owns extension values such as DataFrames.
    # ``asdict`` nests that mapping under ``raw`` and deep-copies its derived
    # fields separately, which is not the input shape expected by from_dict().
    raw["labels"]["filtering"] = deepcopy(config.labels.filtering.as_dict())
    return raw


def _literal_value_matches(value: object, candidate: object) -> bool:
    return type(value) is type(candidate) and value == candidate


def _matches_literal_domain(annotation: object, value: object) -> bool:
    origin = get_origin(annotation)
    if origin is Literal:
        return any(
            _literal_value_matches(value, candidate)
            for candidate in get_args(annotation)
        )
    if origin in {Union, UnionType}:
        for member in get_args(annotation):
            if get_origin(member) is Literal and _matches_literal_domain(member, value):
                return True
            if member is type(None) and value is None:
                return True
            if isinstance(member, type) and isinstance(value, member):
                return True
        return False
    return True


def _has_literal_domain(annotation: object) -> bool:
    origin = get_origin(annotation)
    if origin is Literal:
        return True
    if origin in {Union, UnionType}:
        return any(_has_literal_domain(member) for member in get_args(annotation))
    return False


def _literal_domain_values(annotation: object) -> tuple[object, ...]:
    origin = get_origin(annotation)
    if origin is Literal:
        return get_args(annotation)
    if origin in {Union, UnionType}:
        return tuple(
            value
            for member in get_args(annotation)
            for value in _literal_domain_values(member)
        )
    return ()


def _matches_annotation(annotation: object, value: object) -> bool:
    if annotation is Any:
        return True
    origin = get_origin(annotation)
    if origin is Literal:
        return _matches_literal_domain(annotation, value)
    if origin in {Union, UnionType}:
        return any(
            _matches_annotation(member, value)
            for member in get_args(annotation)
        )
    if annotation is bool:
        return type(value) is bool
    if annotation is int:
        return not isinstance(value, bool) and isinstance(value, Integral)
    if annotation is float:
        return not isinstance(value, bool) and isinstance(value, Real)
    if annotation is str:
        return isinstance(value, str)
    if origin in {dict, Mapping, MappingABC}:
        if not isinstance(value, MappingABC):
            return False
        key_annotation, value_annotation = get_args(annotation) or (Any, Any)
        return all(
            _matches_annotation(key_annotation, key)
            and _matches_annotation(value_annotation, item)
            for key, item in value.items()
        )
    if origin in {list, tuple, SequenceABC}:
        if isinstance(value, (str, bytes)) or not isinstance(value, SequenceABC):
            return False
        item_annotations = get_args(annotation)
        if not item_annotations:
            return True
        if origin is tuple and len(item_annotations) > 1 and item_annotations[-1] is not Ellipsis:
            return len(value) == len(item_annotations) and all(
                _matches_annotation(item_annotation, item)
                for item_annotation, item in zip(item_annotations, value)
            )
        item_annotation = item_annotations[0]
        return all(_matches_annotation(item_annotation, item) for item in value)
    if annotation is type(None):
        return value is None
    if isinstance(annotation, type):
        return isinstance(value, annotation)
    return True


def _validate_override_path(path: object) -> tuple[str, object]:
    leaves, branches = _config_override_schema()
    if not isinstance(path, str) or not path:
        raise ValidationError("Config override paths must be non-empty strings.")
    if path in branches:
        raise ValidationError(
            f"Config override path {path!r} is not a leaf field."
        )
    try:
        return path, leaves[path]
    except KeyError as exc:
        raise ValidationError(f"Unknown config override path: {path!r}.") from exc


def _set_raw_config_path(
    raw_config: dict[str, Any],
    path: str,
    value: object,
) -> None:
    # Filtering ``raw`` is the one typed field whose serialized representation
    # replaces its parent mapping rather than adding a same-named child.
    if path == "labels.filtering.raw":
        if not isinstance(value, Mapping):
            raise ValidationError(
                "Config override 'labels.filtering.raw' must be a mapping."
            )
        raw_config["labels"]["filtering"] = deepcopy(dict(value))
        return

    keys = path.split(".")
    target: dict[str, Any] = raw_config
    for key in keys[:-1]:
        child = target.setdefault(key, {})
        if not isinstance(child, dict):
            raise ValidationError(
                f"Config override path {path!r} crosses non-mapping field {key!r}."
            )
        target = child
    target[keys[-1]] = deepcopy(value)


def validate_config_overrides(
    overrides: Mapping[str, object] | None,
) -> tuple[tuple[str, object], ...]:
    """Validate a canonical override mapping without applying it."""

    if overrides is not None and not isinstance(overrides, Mapping):
        raise ValidationError("config overrides must be a mapping or None.")

    validated: list[tuple[str, object]] = []
    for raw_path, value in dict(overrides or {}).items():
        path, annotation = _validate_override_path(raw_path)
        if _has_literal_domain(annotation) and not _matches_literal_domain(
            annotation,
            value,
        ):
            allowed = ", ".join(
                repr(candidate)
                for candidate in _literal_domain_values(annotation)
            )
            raise ValidationError(
                f"Invalid value {value!r} for config override {path!r}; "
                f"allowed literal values: {allowed}."
            )
        if not _matches_annotation(annotation, value):
            raise ValidationError(
                f"Invalid value {value!r} for config override {path!r}; "
                f"expected {annotation!r}."
            )
        validated.append((path, value))
    return tuple(validated)


def _apply_validated_config_overrides(
    config_dict: Mapping[str, Any],
    overrides: Mapping[str, object] | None,
) -> dict[str, Any]:
    """Apply schema-valid overrides without reparsing the source or result."""

    updated = deepcopy(dict(config_dict))
    validated = list(validate_config_overrides(overrides))

    # Replacing filtering.raw first makes typed sibling paths deterministic,
    # independent of the input mapping's insertion order.
    validated.sort(key=lambda item: item[0] != "labels.filtering.raw")
    for path, value in validated:
        _set_raw_config_path(updated, path, value)
    return updated


def modify_config_dict(
    config_dict: Mapping[str, Any],
    overrides: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    """Apply canonical dotted leaf overrides to a raw configuration.

    Accepted paths are derived recursively from the ``GbdrawConfig`` dataclass
    hierarchy. Flat aliases are intentionally not accepted.
    """

    if not isinstance(config_dict, Mapping):
        raise ValidationError("config_dict must be a mapping.")

    # Validate the source before updating it, while preserving raw top-level
    # metadata that is outside the typed rendering configuration.
    try:
        GbdrawConfig.from_dict(config_dict)
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ValidationError(f"Invalid configuration: {exc}") from exc

    updated = _apply_validated_config_overrides(config_dict, overrides)

    try:
        GbdrawConfig.from_dict(updated)
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ValidationError(
            f"Invalid value for config override: {exc}"
        ) from exc
    return updated


__all__ = [
    "canonical_config_override_paths",
    "config_to_raw_dict",
    "modify_config_dict",
    "validate_config_overrides",
]
