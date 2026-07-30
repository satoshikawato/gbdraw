"""Deterministic helpers for SVG-internal identifiers."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any

_SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")
_VALID_ID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")


def _safe_fragment(value: object, fallback: str = "item") -> str:
    text = str(value or "").strip()
    safe = _SAFE_ID_RE.sub("_", text).strip("_")
    return safe or fallback


def _canonical_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_canonical_value(item) for item in value]
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def stable_svg_id(
    prefix: object,
    *semantic_parts: object,
    namespace: object | None = None,
) -> str:
    """Return a version-local deterministic ID from semantic inputs."""

    payload = json.dumps(
        _canonical_value(semantic_parts),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    identifier = f"{_safe_fragment(prefix)}_{digest}"
    if namespace is not None:
        identifier = f"{identifier}_{_safe_fragment(namespace, 'scope')}"
    return identifier


def instance_svg_id(base_id: object, instance_id: object) -> str:
    """Namespace a semantic ID for one deterministic rendered instance."""

    payload = json.dumps(
        _canonical_value(instance_id),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    return (
        f"{_safe_fragment(base_id)}"
        f"__instance_{_safe_fragment(instance_id)}_{digest}"
    )


def track_slot_svg_id(
    slot_id: object,
    *,
    renderer: object,
    slot_index: int,
) -> str:
    """Return a deterministic DOM ID for one resolved user track slot."""

    return stable_svg_id(
        "track_slot",
        str(renderer).strip().lower(),
        str(slot_id),
        int(slot_index),
        namespace=slot_id,
    )


def definition_group_svg_id(
    record_id: object,
    *,
    mode: object,
    record_index: int = 0,
    record_count: int = 1,
) -> str:
    """Return a deterministic, valid DOM ID for one record definition group.

    Already-valid record IDs retain the established readable spelling. Unsafe
    IDs include a digest of the raw value so inputs that normalize to the same
    fragment remain distinct. Multi-record Linear IDs retain their historical
    suffix order.
    """

    raw_record_id = str(record_id or "").strip()
    normalized_mode = str(mode or "").strip().lower() or "diagram"
    safe_record_id = _safe_fragment(raw_record_id, "record")
    if not _VALID_ID_RE.fullmatch(safe_record_id):
        safe_record_id = f"record_{safe_record_id}"

    if _VALID_ID_RE.fullmatch(raw_record_id):
        base_id = f"{raw_record_id}_definition"
    else:
        digest = hashlib.sha256(
            json.dumps(
                _canonical_value((normalized_mode, raw_record_id)),
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()[:12]
        base_id = f"{safe_record_id}_{digest}_definition"

    if int(record_count) <= 1:
        return base_id
    return f"{base_id}_record_{int(record_index) + 1}"


def record_group_svg_id(
    record_id: object,
    *,
    mode: object,
    record_index: int = 0,
    record_count: int = 1,
) -> str:
    """Return the deterministic renderer-owned ID for one record group."""

    return stable_svg_id(
        "record_group",
        str(mode).strip().lower(),
        str(record_id or "").strip(),
        int(record_index),
        int(record_count),
    )


__all__ = [
    "definition_group_svg_id",
    "instance_svg_id",
    "record_group_svg_id",
    "stable_svg_id",
    "track_slot_svg_id",
]
