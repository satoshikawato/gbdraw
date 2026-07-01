"""Stable SVG feature identifier helpers."""

from __future__ import annotations

import hashlib
import re
from typing import Any

_SAFE_SVG_ID_FRAGMENT_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def _first_location_coordinates(location: Any) -> tuple[int, int, object]:
    if hasattr(location, "parts") and location.parts:
        first_part = location.parts[0]
        return int(first_part.start), int(first_part.end), first_part.strand
    return int(location.start), int(location.end), location.strand


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


def compute_feature_hash(feature: Any, record_id: str | None = None) -> str:
    """Compute the stable data-gbdraw-feature-id for a BioPython feature."""

    start, end, strand = _first_location_coordinates(feature.location)
    return compute_feature_hash_from_parts(
        str(getattr(feature, "type", "") or ""),
        start,
        end,
        strand,
        record_id=record_id,
    )


def compute_feature_object_hash(feature_object: Any) -> str | None:
    """Compute the stable data-gbdraw-feature-id for a rendered FeatureObject."""

    coords = getattr(feature_object, "coordinates", None)
    if not coords:
        return None
    coord = coords[0]
    return compute_feature_hash_from_parts(
        str(
            getattr(feature_object, "feature_type", None)
            or getattr(feature_object, "type", "")
            or ""
        ),
        int(coord.start),
        int(coord.end),
        coord.strand,
        record_id=getattr(feature_object, "record_id", None),
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
    "compute_feature_object_hash",
    "make_linear_dom_id",
    "make_linear_rendered_feature_id",
    "make_svg_safe_id_fragment",
]
