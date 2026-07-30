from __future__ import annotations

from typing import NamedTuple

from Bio.SeqRecord import SeqRecord


_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"
_SOURCE_FEATURE_INDEX_ATTR = "_gbdraw_source_feature_index"
_SOURCE_FEATURE_PARTS_ATTR = "_gbdraw_source_feature_location_parts"


def _feature_source_index_map(features: object) -> dict[int, int]:
    """Map nested feature objects to their flattened source-order ordinals."""

    indexes: dict[int, int] = {}

    def walk(items: object) -> None:
        for feature in items or ():  # type: ignore[union-attr]
            indexes[id(feature)] = len(indexes)
            walk(getattr(feature, "sub_features", None))

    walk(features)
    return indexes


def _source_feature_index(feature: object) -> int | None:
    """Return a pre-transform source ordinal attached to a feature, if any."""

    try:
        index = int(getattr(feature, _SOURCE_FEATURE_INDEX_ATTR))
    except (AttributeError, TypeError, ValueError):
        return None
    return index if index >= 0 else None


def _source_feature_location_parts(
    feature: object,
) -> tuple[tuple[int, int, int | None], ...] | None:
    """Return pre-transform biological location parts, if attached."""

    value = getattr(feature, _SOURCE_FEATURE_PARTS_ATTR, None)
    if not isinstance(value, tuple):
        return None
    parts: list[tuple[int, int, int | None]] = []
    for part in value:
        if not isinstance(part, tuple) or len(part) != 3:
            return None
        try:
            start = int(part[0])
            end = int(part[1])
        except (TypeError, ValueError):
            return None
        strand = part[2]
        if strand not in {-1, 1}:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts) or None


def _mapped_feature_location_parts(
    feature: object,
    *,
    coord_base: int,
    coord_step: int,
) -> tuple[tuple[int, int, int | None], ...]:
    """Map a feature's current parts into its biological source coordinates."""

    location = getattr(feature, "location", None)
    raw_parts = list(getattr(location, "parts", None) or [location])
    parts: list[tuple[int, int, int | None]] = []
    for part in raw_parts:
        if part is None:
            continue
        try:
            start, end = _absolute_display_interval(
                int(part.start),
                int(part.end),
                coord_base,
                coord_step,
            )
        except (AttributeError, TypeError, ValueError):
            continue
        strand = getattr(part, "strand", None)
        if strand in {-1, 1}:
            strand = int(strand) * (1 if int(coord_step) >= 0 else -1)
        else:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts)


def _copy_source_feature_identity(
    source: object,
    target: object,
    *,
    fallback_index: int,
    coord_base: int,
    coord_step: int,
) -> None:
    """Carry source ordinal and biological parts across a transformation."""

    index = _source_feature_index(source)
    setattr(
        target,
        _SOURCE_FEATURE_INDEX_ATTR,
        int(fallback_index) if index is None else index,
    )
    parts = _source_feature_location_parts(source)
    if parts is None:
        parts = _mapped_feature_location_parts(
            source,
            coord_base=coord_base,
            coord_step=coord_step,
        )
    if parts:
        setattr(target, _SOURCE_FEATURE_PARTS_ATTR, parts)


def _read_coord_map(record: object) -> tuple[int, int]:
    annotations = getattr(record, "annotations", None) or {}
    try:
        base = int(annotations.get(_COORD_BASE_KEY, 1))
    except (TypeError, ValueError):
        base = 1
    try:
        step = int(annotations.get(_COORD_STEP_KEY, 1))
    except (TypeError, ValueError):
        step = 1
    if step == 0:
        step = 1
    return base, (1 if step > 0 else -1)


def _write_coord_map(record: object, *, base: int, step: int) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}  # type: ignore[attr-defined]
    record.annotations[_COORD_BASE_KEY] = int(base)  # type: ignore[attr-defined]
    record.annotations[_COORD_STEP_KEY] = 1 if int(step) >= 0 else -1  # type: ignore[attr-defined]


def _absolute_display_interval(
    start: int,
    end: int,
    coord_base: int,
    coord_step: int,
) -> tuple[int, int]:
    if end <= start:
        coord = coord_base + (coord_step * start)
        return coord - 1, coord
    first_coord = coord_base + (coord_step * start)
    last_coord = coord_base + (coord_step * (end - 1))
    return min(first_coord, last_coord) - 1, max(first_coord, last_coord)


class RecordSourceMetadata(NamedTuple):
    organism: str
    strain: str
    replicon: str | None
    organelle: str | None


def infer_record_source_metadata(record: SeqRecord) -> RecordSourceMetadata:
    """Extract source-feature metadata used in definition labels."""
    annotations = getattr(record, "annotations", None) or {}
    organism = str(annotations.get("organism", "") or "").strip()
    strain = ""
    replicon: str | None = None
    organelle: str | None = None

    for feature in getattr(record, "features", []):
        if getattr(feature, "type", None) != "source":
            continue

        qualifiers = getattr(feature, "qualifiers", {}) or {}
        if "organism" in qualifiers and qualifiers["organism"]:
            organism = str(qualifiers["organism"][0]).strip()
        if "isolate" in qualifiers and qualifiers["isolate"]:
            strain = str(qualifiers["isolate"][0]).strip()
        elif "strain" in qualifiers and qualifiers["strain"]:
            strain = str(qualifiers["strain"][0]).strip()

        if "chromosome" in qualifiers and qualifiers["chromosome"]:
            replicon = f"Chromosome {str(qualifiers['chromosome'][0]).strip()}"
        elif "plasmid" in qualifiers and qualifiers["plasmid"]:
            replicon = str(qualifiers["plasmid"][0]).strip()

        if "organelle" in qualifiers and qualifiers["organelle"]:
            organelle = str(qualifiers["organelle"][0]).strip()

    return RecordSourceMetadata(
        organism=organism,
        strain=strain,
        replicon=replicon,
        organelle=organelle,
    )


__all__ = ["RecordSourceMetadata", "infer_record_source_metadata"]
