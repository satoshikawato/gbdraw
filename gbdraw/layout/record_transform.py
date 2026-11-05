"""The transform a sequence record is drawn by, read from the record itself.

``record_coordinates`` stays free of record objects; this module reads a
record's coordinate map and annotations. ``detected_topology`` is the one
reading of a record's ``topology`` annotation (the record planner's and this
module's); ``unset_record_transform`` is the transform the planner binds to a
record given in memory with no display option, the one transform of a record
no plan bound one to.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.layout.record_coordinates import RecordDisplayTransform

DetectedTopology = Literal["circular", "linear", "unknown"]


def detected_topology(record: SeqRecord, source_kind: str) -> DetectedTopology:
    """The topology a record's file states; a GFF3/FASTA source states none."""

    value = str(record.annotations.get("topology", "")).strip().lower()
    if source_kind == "gff_fasta":
        return "unknown"
    if value == "circular":
        return "circular"
    if value == "linear":
        return "linear"
    return "unknown"


def unset_record_transform(record: SeqRecord) -> RecordDisplayTransform:
    """The transform of a record no plan bound one to: the one the record planner binds to the record given in
    memory with no display option (its coordinate map, the complete source length a crop recorded, and circular
    when its file says so)."""

    base, step = _read_coord_map(record)
    length = len(record.seq)
    is_circular = detected_topology(record, "memory") == "circular"
    source_length = record.annotations.get("gbdraw_source_length")
    if record.annotations.get("gbdraw_region_applied") and isinstance(source_length, int):
        return RecordDisplayTransform(source_length, base, 1 if step == 1 else -1, None, is_circular)
    if step == 1:
        return RecordDisplayTransform(base + length - 1, base, 1, None, is_circular)
    return RecordDisplayTransform(max(base, length), base, -1, None, is_circular)


def bind_record_transforms(
    records: Sequence[SeqRecord], transforms: Sequence[RecordDisplayTransform] | None
) -> Sequence[RecordDisplayTransform]:
    """The transforms records are drawn by: the bound ones, else each record's unset transform."""

    return transforms if transforms is not None else [unset_record_transform(record) for record in records]


def spans_may_cross_origin(record: SeqRecord, transform: RecordDisplayTransform) -> bool:
    """Whether a span of the record may cross its origin: the record is drawn as circular (its bound transform)
    and no region was applied to it (a crop has no origin, even one that keeps every base)."""

    return transform.is_circular and not bool(record.annotations.get("gbdraw_region_applied"))
