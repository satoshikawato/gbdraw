from __future__ import annotations

from collections.abc import Mapping as MappingABC
from types import MappingProxyType
from typing import Any, Mapping


def freeze_layout_value(value: Any) -> Any:
    if isinstance(value, MappingABC):
        return MappingProxyType(
            {
                key: freeze_layout_value(item)
                for key, item in value.items()
            }
        )
    if isinstance(value, (list, tuple)):
        return tuple(freeze_layout_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(freeze_layout_value(item) for item in value)
    return value


def freeze_layout_mapping(value: Mapping[Any, Any]) -> Mapping[Any, Any]:
    frozen = freeze_layout_value(value)
    if not isinstance(frozen, MappingABC):  # pragma: no cover - type invariant
        raise TypeError("expected a mapping")
    return frozen


__all__ = ["freeze_layout_mapping", "freeze_layout_value"]
