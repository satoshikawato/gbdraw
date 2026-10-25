"""Stable SVG feature identifier helpers."""

from __future__ import annotations

from collections import Counter, defaultdict

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from gbdraw.core.record_metadata import (
    _iter_source_features,
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_location_parts,
    _stamped_source_feature_hash,
)
from gbdraw.exceptions import ValidationError

import hashlib
import re
from typing import Any, Iterable, NewType

# Frame-tagged feature hashes. Both are plain strings at run time; the types
# keep the two frames from being mixed (OV-412). A SOURCE hash is the hash of
# the feature in its source record: the SVG ID base, the catalog's stable
# feature ID, and what a `hash` selector matches. A DRAWN hash is the hash of
# the processed (cropped, reversed) location, which only comparison rows
# (`*_view_feature_svg_id`) and Session <= 44 readers name a feature by.
SourceFeatureHash = NewType("SourceFeatureHash", str)
DrawnFeatureHash = NewType("DrawnFeatureHash", str)

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


def parse_source_feature_hash(raw: str) -> SourceFeatureHash | None:
    """Tag a stable feature ID read from an SVG attribute or a JSON payload.

    ``raw`` is already normalized by the reader; an empty string is no ID.
    """

    return SourceFeatureHash(raw) if raw else None


def compute_feature_hash(feature: Any, record_id: str | None = None) -> DrawnFeatureHash:
    """Compute the stable data-gbdraw-feature-id for a BioPython feature."""

    parts = _location_parts(feature.location)
    return DrawnFeatureHash(
        compute_feature_hash_from_location_parts(
            str(getattr(feature, "type", "") or ""),
            parts,
            record_id=record_id,
        )
    )


def source_feature_location_parts(
    feature: Any,
    record: Any = None,
) -> tuple[tuple[int, int, int | None], ...]:
    """The feature's location parts in its source record.

    A cropped or reverse-complemented record keeps each feature's source parts
    (``core/record_metadata.py``); any other record maps its parts through its
    coordinate map (none: the record is its own source).
    """

    coord_base, coord_step = (1, 1) if record is None else _read_coord_map(record)
    return _source_feature_location_parts(feature) or _mapped_feature_location_parts(
        feature,
        coord_base=coord_base,
        coord_step=coord_step,
    )


def computed_source_feature_hash_and_parts(
    feature: SeqFeature,
    record: SeqRecord | None = None,
    *,
    record_id: str | None = None,
) -> tuple[SourceFeatureHash | None, tuple[tuple[int, int, int | None], ...]]:
    """Compute the feature hash and the source location parts it hashes.

    Only loading a source (``capture_source_features``) computes; readers use
    ``source_feature_hash_and_parts``, which returns the hash computed then.
    The hash is None when the feature has no source parts. ``record_id``
    defaults to the record's ID.
    """

    parts = source_feature_location_parts(feature, record)
    if not parts:
        return None, parts
    return (
        SourceFeatureHash(
            compute_feature_hash_from_location_parts(
                str(getattr(feature, "type", "") or ""),
                parts,
                record_id=_hash_record_id(record, record_id),
            )
        ),
        parts,
    )


def _hash_record_id(record: object, record_id: str | None) -> str | None:
    return record_id if record_id is not None else getattr(record, "id", None)


def source_feature_hash_and_parts(
    feature: SeqFeature,
    record: SeqRecord | None = None,
    *,
    record_id: str | None = None,
) -> tuple[SourceFeatureHash | None, tuple[tuple[int, int, int | None], ...]]:
    """The feature hash and the source location parts it hashes.

    The hash computed when the feature's source was loaded, carried by every
    crop, reverse complement or type filter copy, is returned for the record ID
    it was computed for; any other record ID computes. The hash is None when the
    feature has no source parts. ``record_id`` defaults to the record's ID.
    """

    stamped = _stamped_source_feature_hash(feature, _hash_record_id(record, record_id))
    if stamped is None:
        return computed_source_feature_hash_and_parts(feature, record, record_id=record_id)
    return SourceFeatureHash(stamped), source_feature_location_parts(feature, record)


def compute_source_feature_hash(
    feature: Any,
    record: Any = None,
    *,
    record_id: str | None = None,
) -> SourceFeatureHash | None:
    """The feature hash: the catalog's stable feature ID, the base of every
    SVG ID drawn for the feature, and the value a `hash` selector matches,
    whatever crop or reverse complement its record has (OV-401, OV-412).

    ``record_id`` defaults to the record's ID. Like
    ``source_feature_hash_and_parts``, it returns the hash computed when the
    source was loaded.
    """

    stamped = _stamped_source_feature_hash(feature, _hash_record_id(record, record_id))
    if stamped is not None:
        return SourceFeatureHash(stamped)
    return computed_source_feature_hash_and_parts(feature, record, record_id=record_id)[0]


def shared_source_feature_hashes(record: SeqRecord) -> frozenset[SourceFeatureHash]:
    """The source hashes that more than one feature of ``record`` has.

    Identical features (the same type and location) share their hash; the
    rendered ID of each adds its source index as an instance suffix
    (``gbdraw.svg.ids.instance_svg_id``). The count reads every source feature
    of the record, drawn or not, so hiding one of identical features leaves the
    rendered ID of the other as it was (OV-1485).
    """

    counts = Counter(
        compute_source_feature_hash(feature, record)
        for feature in _iter_source_features(record.features)
    )
    return frozenset(
        feature_hash for feature_hash, count in counts.items() if feature_hash and count > 1
    )


def compute_feature_object_hash(
    feature_object: Any,
    record_id: str | None = None,
) -> DrawnFeatureHash | None:
    """The hash of a FeatureObject's drawn location.

    Only comparison rows name a drawn feature by it (`*_view_feature_svg_id`,
    a frame tag). SVG IDs and `hash` selectors use the source hash
    (`compute_source_feature_hash`, `FeatureObject.feature_hash`).
    """

    parts = _feature_object_coordinate_parts(feature_object)
    if not parts:
        return None
    feature_type = str(
        getattr(feature_object, "feature_type", None)
        or getattr(feature_object, "type", "")
        or ""
    )
    return DrawnFeatureHash(
        compute_feature_hash_from_location_parts(
            feature_type,
            parts,
            record_id=record_id if record_id is not None else getattr(feature_object, "record_id", None),
        )
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


def disambiguate_feature_ids(
    identities: Iterable[tuple[str, str, int | None]],
    *,
    description: str = "Source feature catalog",
) -> tuple[str, ...]:
    """Keep legacy public IDs, disambiguating equal hashes by source ordinal."""
    rows = tuple(identities)
    counts = Counter((record_key, base) for record_key, base, _ in rows)
    indexes: dict[tuple[str, str], list[int | None]] = defaultdict(list)
    for record_key, base, index in rows:
        indexes[(record_key, base)].append(index)
    for identity, values in indexes.items():
        if len(values) > 1 and (None in values or len(set(values)) != len(values)):
            raise ValidationError(
                f"{description} contains duplicate biological source identity {identity!r}."
            )
    resolved = tuple(
        f"{base}~{index}" if counts[(record_key, base)] > 1 else base
        for record_key, base, index in rows
    )
    if len({(row[0], value) for row, value in zip(rows, resolved)}) != len(rows):
        raise ValidationError(
            f"{description} generated a duplicate canonical biological feature identity."
        )
    return resolved


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
    "compute_source_feature_hash",
    "source_feature_location_parts",
    "make_linear_dom_id",
    "make_linear_rendered_feature_id",
    "make_svg_safe_id_fragment",
]
