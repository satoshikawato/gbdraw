"""Stable SVG feature identifier helpers."""

from __future__ import annotations

import hashlib
from typing import Any


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


__all__ = [
    "compute_feature_hash",
    "compute_feature_hash_from_parts",
    "compute_feature_object_hash",
]
