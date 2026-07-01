#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Any, Optional

from Bio.SeqFeature import SeqFeature

from .ids import compute_feature_hash, compute_feature_hash_from_parts


def normalize_qualifier_values(raw_values: Any) -> list[str]:
    if raw_values is None:
        return []
    if isinstance(raw_values, (list, tuple, set)):
        return [str(value) for value in raw_values if value is not None]
    return [str(raw_values)]


def get_feature_type(feature: Any) -> str:
    feature_type = getattr(feature, "type", None)
    if feature_type is None:
        feature_type = getattr(feature, "feature_type", "")
    return str(feature_type or "")


def get_feature_qualifiers(feature: Any) -> dict:
    qualifiers = getattr(feature, "qualifiers", {})
    if qualifiers is None:
        return {}
    if isinstance(qualifiers, dict):
        return qualifiers
    try:
        return dict(qualifiers)
    except Exception:
        return {}


def get_feature_record_id(feature: Any, record_id: Optional[str]) -> Optional[str]:
    if record_id is not None and str(record_id).strip() != "":
        return str(record_id)
    feature_record_id = getattr(feature, "record_id", None)
    if feature_record_id is None or str(feature_record_id).strip() == "":
        return None
    return str(feature_record_id)


def iter_feature_coordinates(feature: Any) -> list[Any]:
    coordinates = getattr(feature, "coordinates", None)
    if coordinates is None:
        return []
    try:
        return list(coordinates)
    except Exception:
        return []


def extract_first_coordinate_part(feature: Any) -> Optional[tuple[int, int, Any]]:
    for part in iter_feature_coordinates(feature):
        try:
            start = int(getattr(part, "start"))
            end = int(getattr(part, "end"))
            strand = getattr(part, "strand", None)
        except Exception:
            continue
        return start, end, strand
    return None


def extract_coordinate_bounds(feature: Any) -> Optional[tuple[int, int]]:
    starts: list[int] = []
    ends: list[int] = []
    for part in iter_feature_coordinates(feature):
        try:
            starts.append(int(getattr(part, "start")))
            ends.append(int(getattr(part, "end")))
        except Exception:
            continue
    if not starts or not ends:
        return None
    return min(starts), max(ends)


def extract_first_location_part(feature: Any) -> Optional[tuple[int, int, Any]]:
    if isinstance(feature, SeqFeature):
        loc = feature.location
        if hasattr(loc, "parts") and loc.parts:
            part = loc.parts[0]
            return int(part.start), int(part.end), part.strand
        return int(loc.start), int(loc.end), loc.strand

    coordinate_part = extract_first_coordinate_part(feature)
    if coordinate_part is not None:
        return coordinate_part

    location = getattr(feature, "location", None)
    if not location:
        return None

    selected = None
    for part in location:
        kind = getattr(part, "kind", None)
        if kind is None and isinstance(part, tuple) and len(part) >= 1:
            kind = part[0]
        if kind == "block":
            selected = part
            break
    if selected is None:
        selected = location[0]

    try:
        start = getattr(selected, "start", selected[3])
        end = getattr(selected, "end", selected[4])
        strand = getattr(selected, "strand", selected[2])
    except Exception:
        return None
    try:
        return int(start), int(end), strand
    except Exception:
        return None


def normalize_strand_for_hash(strand: Any) -> Any:
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
        try:
            return int(normalized)
        except ValueError:
            return strand
    if isinstance(strand, (int, float)):
        try:
            return int(strand)
        except Exception:
            return strand
    return strand


def normalize_strand_token(strand: Any) -> str:
    if strand in (None, "", "none", "None", "undefined"):
        return "undefined"
    if isinstance(strand, str):
        normalized = strand.strip().lower()
        if normalized in {"positive", "plus", "+", "forward", "1"}:
            return "+"
        if normalized in {"negative", "minus", "-", "reverse", "-1"}:
            return "-"
        if normalized in {"undefined", "none", ""}:
            return "undefined"
        return "undefined"
    if isinstance(strand, (int, float)):
        try:
            numeric = int(strand)
        except Exception:
            return "undefined"
        if numeric == 1:
            return "+"
        if numeric == -1:
            return "-"
        return "undefined"
    return "undefined"


def get_feature_hash(feature: Any, record_id: Optional[str]) -> Optional[str]:
    resolved_record_id = get_feature_record_id(feature, record_id)
    if isinstance(feature, SeqFeature):
        return compute_feature_hash(feature, record_id=resolved_record_id)

    first_part = extract_first_location_part(feature)
    feature_type = get_feature_type(feature)
    if not first_part or not feature_type:
        return None
    start, end, strand = first_part
    normalized_strand = normalize_strand_for_hash(strand)
    return compute_feature_hash_from_parts(
        feature_type,
        start,
        end,
        normalized_strand,
        record_id=resolved_record_id,
    )


def get_feature_location_str(feature: Any) -> Optional[str]:
    if isinstance(feature, SeqFeature):
        try:
            return f"{int(feature.location.start)}..{int(feature.location.end)}"
        except Exception:
            return None

    coordinate_bounds = extract_coordinate_bounds(feature)
    if coordinate_bounds is not None:
        return f"{coordinate_bounds[0]}..{coordinate_bounds[1]}"

    location = getattr(feature, "location", None)
    if not location:
        return None

    starts: list[int] = []
    ends: list[int] = []
    for part in location:
        kind = getattr(part, "kind", None)
        if kind is None and isinstance(part, tuple) and len(part) >= 1:
            kind = part[0]
        if kind not in {"block", None}:
            continue
        try:
            start = int(getattr(part, "start", part[3]))
            end = int(getattr(part, "end", part[4]))
        except Exception:
            continue
        starts.append(start)
        ends.append(end)
    if not starts or not ends:
        return None
    return f"{min(starts)}..{max(ends)}"


def get_feature_position_str(feature: Any) -> Optional[str]:
    location = get_feature_location_str(feature)
    if not location:
        return None

    strand_source = None
    first_part = extract_first_location_part(feature)
    if first_part is not None:
        strand_source = first_part[2]
    strand = normalize_strand_token(strand_source)
    return f"{location}:{strand}"


def get_feature_record_location_str(
    feature: Any,
    record_id: Optional[str],
) -> Optional[str]:
    resolved_record_id = get_feature_record_id(feature, record_id)
    if not resolved_record_id:
        return None

    position = get_feature_position_str(feature)
    if not position:
        return None

    return f"{resolved_record_id}:{position}"


def get_qualifier_values(qualifiers: dict, qualifier_key: str) -> list[str]:
    for key, values in qualifiers.items():
        if str(key).lower() != str(qualifier_key).lower():
            continue
        return normalize_qualifier_values(values)
    return []


def build_feature_selector_values(
    feature: Any,
    record_id: Optional[str] = None,
) -> dict[str, object]:
    qualifiers: dict[str, list[str]] = {}
    for key, values in get_feature_qualifiers(feature).items():
        normalized_key = str(key or "").strip().lower()
        if not normalized_key:
            continue
        normalized_values = normalize_qualifier_values(values)
        if not normalized_values:
            continue
        qualifiers.setdefault(normalized_key, []).extend(normalized_values)

    selector: dict[str, object] = {"qualifiers": qualifiers}
    feature_hash = get_feature_hash(feature, record_id)
    if feature_hash:
        selector["hash"] = feature_hash
    location = get_feature_location_str(feature)
    if location:
        selector["location"] = location
    record_location = get_feature_record_location_str(feature, record_id)
    if record_location:
        selector["record_location"] = record_location
    return selector


__all__ = [
    "build_feature_selector_values",
    "get_feature_hash",
    "get_feature_location_str",
    "get_feature_qualifiers",
    "get_feature_record_id",
    "get_feature_record_location_str",
    "get_feature_type",
    "get_qualifier_values",
    "normalize_qualifier_values",
    "normalize_strand_token",
]
