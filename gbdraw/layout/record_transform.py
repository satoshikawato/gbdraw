"""The transform a sequence record is drawn by, read from the record itself.

``record_coordinates`` stays free of record objects; this module reads a
record's coordinate map and annotations to give the transform of a record no
plan bound one to, and whether a span of the record may cross its origin.
"""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.layout.record_coordinates import RecordDisplayTransform


def unset_record_transform(record: SeqRecord) -> RecordDisplayTransform:
    """The transform of a record no plan bound one to, from its coordinate map."""

    base, step = _read_coord_map(record)
    length = len(record.seq)
    if step == 1:
        return RecordDisplayTransform(base + length - 1, base, 1, None, False)
    return RecordDisplayTransform(max(base, length), base, -1, None, False)


def bound_record_transform(record: SeqRecord, transform: RecordDisplayTransform | None) -> RecordDisplayTransform:
    """The transform a record is drawn by: its bound one, else the one its coordinate map gives."""

    return transform if transform is not None else unset_record_transform(record)


def spans_may_cross_origin(record: SeqRecord, transform: RecordDisplayTransform | None) -> bool:
    """Whether a span of the record may cross its origin: the record is drawn as circular (its transform, else its
    file's topology) and no region was applied to it (a crop has no origin, even one that keeps every base)."""

    is_circular = (
        transform.is_circular if transform is not None
        else str(record.annotations.get("topology", "")).strip().lower() == "circular"
    )
    return is_circular and not bool(record.annotations.get("gbdraw_region_applied"))
