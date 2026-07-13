from __future__ import annotations

from typing import NamedTuple

from Bio.SeqRecord import SeqRecord


_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"


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
