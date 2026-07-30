"""Surface-neutral serialization of resolved track-slot geometry."""

from __future__ import annotations

from typing import Any, Literal, Mapping, Sequence


def collect_track_slot_geometry_records(
    canvas: Any,
    *,
    result_index: int,
    result_name: str,
) -> list[dict[str, Any]]:
    """Copy track-slot geometry from a canvas and add result identity."""

    geometry = getattr(canvas, "_gbdraw_track_slot_geometry", None)
    if not isinstance(geometry, Mapping):
        return []
    records = geometry.get("records")
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        return []
    copied_records: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, Mapping):
            continue
        copied_record = dict(record)
        copied_record["resultIndex"] = int(result_index)
        copied_record["resultName"] = str(result_name)
        slots = copied_record.get("slots")
        if isinstance(slots, Sequence) and not isinstance(slots, (str, bytes)):
            copied_record["slots"] = [
                dict(slot) for slot in slots if isinstance(slot, Mapping)
            ]
        else:
            copied_record["slots"] = []
        copied_records.append(copied_record)
    return copied_records


def build_track_slot_geometry_run_metadata(
    *,
    mode: Literal["circular", "linear"],
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build optional run metadata for resolved track-slot geometry."""

    if not records:
        return {}
    return {
        "trackSlotGeometry": {
            "schema": 1,
            "mode": mode,
            "source": "resolved",
            "records": [dict(record) for record in records],
        }
    }


__all__ = [
    "build_track_slot_geometry_run_metadata",
    "collect_track_slot_geometry_records",
]
