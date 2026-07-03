"""Stable SVG feature identifier helpers."""

from __future__ import annotations

import hashlib
import re
from typing import Any, Iterable

_SAFE_SVG_ID_FRAGMENT_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def _hash_feature_key(
    feature_type: str,
    start: int,
    end: int,
    strand: object,
    *,
    record_id: str | None = None,
) -> str:
    if record_id is not None:
        key = f"{record_id}:{feature_type}:{start}:{end}:{strand}"
    else:
        key = f"{feature_type}:{start}:{end}:{strand}"
    return "f" + hashlib.md5(key.encode()).hexdigest()[:8]


def _hash_feature_parts_key(
    feature_type: str,
    parts: Iterable[tuple[int, int, object]],
    *,
    record_id: str | None = None,
) -> str:
    materialized_parts = list(parts)
    normalized_parts = [
        f"{int(start)}:{int(end)}:{strand}"
        for start, end, strand in materialized_parts
    ]
    if len(normalized_parts) == 1:
        start, end, strand = materialized_parts[0]
        return _hash_feature_key(
            feature_type,
            int(start),
            int(end),
            strand,
            record_id=record_id,
        )
    location_key = ";".join(normalized_parts)
    if record_id is not None:
        key = f"{record_id}:{feature_type}:{location_key}"
    else:
        key = f"{feature_type}:{location_key}"
    return "f" + hashlib.md5(key.encode()).hexdigest()[:8]


def _normal_hash_strand(strand: object) -> object:
    if strand in (None, "", "none", "None", "undefined"):
        return None
    if isinstance(strand, str):
        normalized = strand.strip().lower()
        if normalized in {"positive", "plus", "+", "forward", "1"}:
            return 1
        if normalized in {"negative", "minus", "-", "reverse", "-1"}:
            return -1
        if normalized in {"undefined", "none", ""}:
            return None
    if isinstance(strand, (int, float)):
        try:
            return int(strand)
        except Exception:
            return strand
    return strand


def _location_parts(location: Any) -> list[tuple[int, int, object]]:
    if hasattr(location, "parts") and location.parts:
        raw_parts = location.parts
    else:
        raw_parts = [location]
    parts: list[tuple[int, int, object]] = []
    for part in raw_parts:
        parts.append((
            int(part.start),
            int(part.end),
            _normal_hash_strand(part.strand),
        ))
    return parts


def _feature_object_coordinate_parts(feature_object: Any) -> list[tuple[int, int, object]]:
    coords = getattr(feature_object, "coordinates", None)
    parts: list[tuple[int, int, object]] = []
    if coords:
        for coord in coords:
            try:
                parts.append((
                    int(getattr(coord, "start")),
                    int(getattr(coord, "end")),
                    _normal_hash_strand(getattr(coord, "strand", None)),
                ))
            except Exception:
                continue
    if parts:
        return parts

    location = getattr(feature_object, "location", None)
    if not location:
        return []
    for part in location:
        kind = getattr(part, "kind", None)
        if kind is None and isinstance(part, tuple) and len(part) >= 1:
            kind = part[0]
        if kind not in {"block", None}:
            continue
        try:
            start = getattr(part, "start", part[3])
            end = getattr(part, "end", part[4])
            strand = getattr(part, "strand", part[2])
            parts.append((int(start), int(end), _normal_hash_strand(strand)))
        except Exception:
            continue
    return parts


def compute_feature_hash(feature: Any, record_id: str | None = None) -> str:
    """Compute the stable data-gbdraw-feature-id for a BioPython feature."""

    parts = _location_parts(feature.location)
    return compute_feature_hash_from_location_parts(
        str(getattr(feature, "type", "") or ""),
        parts,
        record_id=record_id,
    )


def compute_feature_object_hash(
    feature_object: Any,
    record_id: str | None = None,
) -> str | None:
    """Compute the stable data-gbdraw-feature-id for a rendered FeatureObject."""

    parts = _feature_object_coordinate_parts(feature_object)
    if not parts:
        return None
    feature_type = str(
        getattr(feature_object, "feature_type", None)
        or getattr(feature_object, "type", "")
        or ""
    )
    return compute_feature_hash_from_location_parts(
        feature_type,
        parts,
        record_id=record_id if record_id is not None else getattr(feature_object, "record_id", None),
    )


def compute_feature_hash_from_parts(
    feature_type: str,
    start: int,
    end: int,
    strand: object,
    *,
    record_id: str | None = None,
) -> str:
    """Compute the stable data-gbdraw-feature-id from explicit feature parts."""

    return _hash_feature_key(
        feature_type,
        int(start),
        int(end),
        strand,
        record_id=record_id,
    )


def compute_feature_hash_from_location_parts(
    feature_type: str,
    parts: Iterable[tuple[int, int, object]],
    *,
    record_id: str | None = None,
) -> str:
    """Compute a stable feature id from all rendered location parts."""

    return _hash_feature_parts_key(feature_type, parts, record_id=record_id)


def make_svg_safe_id_fragment(value: object, fallback: str = "item") -> str:
    """Return a compact SVG-id-safe fragment without changing already-safe ids."""

    text = str(value or "").strip()
    safe = _SAFE_SVG_ID_FRAGMENT_RE.sub("_", text).strip("_")
    return safe or fallback


def make_linear_dom_id(
    base_id: object,
    *,
    record_index: int,
    record_count: int,
    suffix: str | None = None,
) -> str:
    """Return a linear SVG DOM id, preserving single-record legacy ids when possible."""

    base = make_svg_safe_id_fragment(base_id, "record")
    if suffix:
        base = f"{base}_{make_svg_safe_id_fragment(suffix, 'group')}"
    if int(record_count) <= 1:
        return base
    return f"{base}_record_{int(record_index) + 1}"


def make_linear_rendered_feature_id(
    *,
    record_index: int,
    stable_feature_id: str | None,
    record_count: int,
) -> str | None:
    """Return the rendered linear feature id for one displayed record instance."""

    if not stable_feature_id:
        return None
    stable_id = make_svg_safe_id_fragment(stable_feature_id, "")
    if not stable_id:
        return None
    if int(record_count) <= 1:
        return stable_id
    return f"{stable_id}_record_{int(record_index) + 1}"


__all__ = [
    "compute_feature_hash",
    "compute_feature_hash_from_parts",
    "compute_feature_hash_from_location_parts",
    "compute_feature_object_hash",
    "make_linear_dom_id",
    "make_linear_rendered_feature_id",
    "make_svg_safe_id_fragment",
]
