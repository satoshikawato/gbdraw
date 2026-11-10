#!/usr/bin/env python
# coding: utf-8

"""Record selection and reverse-complement helpers for linear mode."""

from __future__ import annotations

import copy
from dataclasses import dataclass
import logging
import re
from typing import Sequence

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from ..analysis.windows import WINDOW_SOURCE_ANNOTATION, WindowSource, reversed_window_source, window_source
from ..core.numeric import parse_ascii_nonnegative_integer
from ..core.record_metadata import (
    _copy_source_feature_identity,
    _feature_source_index_map,
    _holder_of,
    _read_coord_map,
    _source_period,
    _write_coord_map,
    record_identity_key,
)
from .record_sequences import sequence_view

logger = logging.getLogger(__name__)

# A ``#N`` record selector, as ``parse_record_selector`` reads it: spaces (``str.strip``'s, which are ``\s``'s)
# around it and after ``#``, an optional ``+`` (Linear Record Order read ``#+2`` as record 2, OV-2988), and N in
# ASCII digits, leading zeros allowed (``#01`` is record 1). The one reader of the form: Record Order in both modes,
# ``--record_id``, ``--region``, and the Web, which reads this pattern from the generated contracts
# (``RECORD_INDEX_SELECTOR_PATTERN``).
RECORD_INDEX_SELECTOR_PATTERN = r"\s*#\s*\+?([0-9]+)\s*"
_RECORD_INDEX_SELECTOR = re.compile(RECORD_INDEX_SELECTOR_PATTERN)


def record_index_number(text: str) -> int | None:
    """The N of a ``#N`` record selector text, 0 included (a reader refuses it); None for any other text."""

    match = _RECORD_INDEX_SELECTOR.fullmatch(text)
    return parse_ascii_nonnegative_integer(match.group(1)) if match else None


def record_selector_number(text: object) -> int | None:
    """The record, from 1, that a ``#N`` record selector names as the renderer reads it
    (:func:`parse_record_selector`): ``#01`` and ``# 1`` are record 1. None for a record ID
    or any other value, and for ``#0``, which the renderer refuses."""

    number = record_index_number(text) if isinstance(text, str) else None
    return number if number else None


@dataclass(frozen=True)
class RecordSelector:
    raw: str
    record_id: str | None
    record_index: int | None

    def label(self) -> str:
        if self.record_index is not None:
            return f"#{self.record_index + 1}"
        if self.record_id:
            return self.record_id
        return "(unspecified)"


def parse_record_selector(text: str | None) -> RecordSelector | None:
    if text is None:
        return None
    raw = str(text).strip()
    if not raw or raw.lower() in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
        return None
    if raw.startswith("#"):
        number = record_index_number(raw)
        if number is None:
            raise ValueError(f"Invalid record selector '{raw}'. Use #<number> or record_id.")
        idx = number - 1
        if idx < 0:
            raise ValueError(f"Record index must be >= 1 in selector '{raw}'.")
        return RecordSelector(raw=raw, record_id=None, record_index=idx)
    return RecordSelector(raw=raw, record_id=raw, record_index=None)


def parse_record_selectors(specs: Sequence[str] | None) -> list[RecordSelector]:
    if not specs:
        return []
    parsed: list[RecordSelector] = []
    for spec in specs:
        selector = parse_record_selector(spec)
        if selector is None:
            continue
        parsed.append(selector)
    return parsed


def select_record(
    records: Sequence[SeqRecord],
    selector: RecordSelector | None,
    *,
    log: logging.Logger | None = None,
) -> list[SeqRecord]:
    if selector is None:
        return list(records)
    log = log or logger
    total = len(records)
    if selector.record_index is not None:
        idx = selector.record_index
        if idx < 0 or idx >= total:
            raise ValueError(
                f"Record selector {selector.label()} is out of range (loaded {total} record(s))."
            )
        return [records[idx]]
    record_id = selector.record_id or ""
    matches = [rec for rec in records if rec.id == record_id]
    if not matches:
        raise ValueError(f"Record selector '{record_id}' did not match any record ID.")
    if len(matches) > 1:
        raise ValueError(
            f"Record selector '{record_id}' matched multiple records. Use #index to disambiguate."
        )
    return [matches[0]]


def drawn_record_units(records: Sequence[SeqRecord]) -> tuple[tuple[int, ...], ...]:
    """The drawn records grouped by the record they draw, in order: the units of a record with windows share one
    group (``record_identity_key``); any other drawn record is a group of its own."""

    groups: dict[object, list[int]] = {}
    for index, record in enumerate(records):
        key = record_identity_key(record)
        groups.setdefault(index if key is None else key, []).append(index)
    return tuple(tuple(indexes) for indexes in groups.values())


def select_record_units(records: Sequence[SeqRecord], selector: RecordSelector | None) -> tuple[int, ...]:
    """The drawn units of the one record ``selector`` names among drawn ``records`` (D-01).

    ``#index`` counts records, not units, and a record ID names the record
    whose units carry it; the units of a record with windows are selected
    together. ``None`` selects the only record. Raises ``ValueError`` when the
    selector names no record or more than one.
    """

    groups = drawn_record_units(records)
    if selector is None:
        if len(groups) != 1:
            raise ValueError("A record selector is required when more than one record is drawn.")
        return groups[0]
    if selector.record_index is not None:
        if not 0 <= selector.record_index < len(groups):
            raise ValueError(
                f"Record selector {selector.label()} is out of range ({len(groups)} record(s))."
            )
        return groups[selector.record_index]
    record_id = selector.record_id or ""
    matches = [group for group in groups if any(records[index].id == record_id for index in group)]
    if not matches:
        raise ValueError(f"Record selector '{record_id}' did not match any record ID.")
    if len(matches) > 1:
        raise ValueError(
            f"Record selector '{record_id}' matched multiple records. Use #index to disambiguate."
        )
    return matches[0]


def _feature_start(feature: SeqFeature) -> float:
    """The start a reverse complement orders features by; one with no location or no numeric start sorts last."""
    try:
        return int(feature.location.start)
    except (AttributeError, TypeError):
        return float("inf")


def edited_crop_draw_order(crop: SeqRecord) -> list[SeqFeature]:
    """The features of ``crop``, a crop edited since its cut, in the order its Session replay draws them.

    The replay cuts again the record the crop was cut from, with the edits on it
    (``session_request_codec._with_crop_edits``), where each feature takes its source feature index
    (``_feature_source_index_map``) as its place: a forward cut lists the features by that index, a
    reverse-complemented one by start, ties by that index, as ``reverse_records`` lists them (OV-3229). A copy of a
    cut feature, and a feature that carries an index the crop did not cut (``WindowSource.cut_indexes``), is a
    feature added at its place (``_holder_of``), as the replay reads it (OV-3247, OV-3256).

    Which of two holders of one carried index is the cut feature is decided on ``crop.features``, the user's list;
    the added features then take the indexes after the cut ones in that list's order. A reader of the drawn crop
    reads the indexes from the returned list, where added features take theirs in the order drawn, as on replay.
    """

    cut = window_source(crop).cut_indexes
    carried = _feature_source_index_map(crop.features)
    holders = [_holder_of(feature, carried[id(feature)], cut) for feature in crop.features]
    indexes = _feature_source_index_map(holders)
    if _read_coord_map(crop)[1] < 0:
        return sorted(holders, key=lambda feature: (_feature_start(feature), indexes[id(feature)]))
    return sorted(holders, key=lambda feature: indexes[id(feature)])


def in_draw_order(record: SeqRecord) -> SeqRecord:
    """``record`` with its features in the order a figure draws them and a Session writes them.

    A crop edited since its cut is a shallow copy that lists its features in ``edited_crop_draw_order``, the order
    its replay draws them in; a crop saved as its own bytes is written in that order too. Any other record, a crop
    as cut too, is ``record`` itself, in its own order. Readers of the drawn crop call ``_feature_source_index_map``
    with the returned record's list, not the caller's: an added feature takes its index in the order drawn.
    """

    stored = record.annotations.get(WINDOW_SOURCE_ANNOTATION)
    if not (
        record.annotations.get("gbdraw_region_applied")
        and isinstance(stored, WindowSource)
        and not stored.holds_for(record)
    ):
        return record
    ordered = copy.copy(record)
    ordered.features = edited_crop_draw_order(record)
    return ordered


def complemented_as_rna(record: SeqRecord) -> bool:
    """Whether ``record``'s reverse complement complements it as RNA (A to U), as Biopython's
    ``SeqRecord.reverse_complement`` does: when its ``molecule_type`` says RNA. One that says protein has none."""

    molecule_type = str(record.annotations.get("molecule_type", ""))
    if "protein" in molecule_type:
        raise ValueError("Proteins do not have complements!")
    return "RNA" in molecule_type


def _reverse_complement_record(record: SeqRecord, features: list[SeqFeature]) -> SeqRecord:
    """``record.reverse_complement(id=True, name=True, description=True, features=features, annotations=True,
    letter_annotations=True, dbxrefs=True)``, its sequence a reverse-strand view of ``record``'s
    (:func:`complemented_as_rna`)."""

    reversed_record = SeqRecord(
        sequence_view(record.seq, strand=-1, rna=complemented_as_rna(record)),
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs[:],
        features=features,
        annotations=record.annotations.copy(),
    )
    for key, value in record.letter_annotations.items():
        reversed_record.letter_annotations[key] = value[::-1]
    return reversed_record


def reverse_records(
    records: Sequence[SeqRecord],
    reverse_flag: bool,
    *,
    log: logging.Logger | None = None,
) -> list[SeqRecord]:
    if not reverse_flag:
        return list(records)
    log = log or logger
    reversed_records: list[SeqRecord] = []
    for rec in records:
        base, step = _read_coord_map(rec)
        record_length = len(rec.seq)
        rc_base = base + (step * max(0, record_length - 1))
        rc_step = -step
        source_indexes = _feature_source_index_map(rec.features)
        period = _source_period(rec)

        def reverse_feature(feature):
            reversed_feature = feature._flip(record_length)
            _copy_source_feature_identity(
                feature,
                reversed_feature,
                index=source_indexes[id(feature)],
                coord_base=base,
                coord_step=step,
                period=period,
            )
            sub_features = getattr(feature, "sub_features", None)
            if sub_features is not None:
                reversed_sub_features = [
                    reverse_feature(sub_feature)
                    for sub_feature in sub_features
                ]
                reversed_sub_features.sort(key=_feature_start)
                reversed_feature.sub_features = reversed_sub_features
            return reversed_feature

        reversed_features = [
            reverse_feature(feature)
            for feature in rec.features
        ]
        if rc_step > 0:
            # Reversed back to its source's orientation, a record lists its features in source order, as its source
            # and a forward cut of it do (OV-3271).
            order = {id(reversed_feature): source_indexes[id(feature)]
                     for feature, reversed_feature in zip(rec.features, reversed_features, strict=True)}
            reversed_features.sort(key=lambda feature: order[id(feature)])
        else:
            reversed_features.sort(key=_feature_start)
        reversed_record = _reverse_complement_record(rec, reversed_features)
        _write_coord_map(reversed_record, base=rc_base, step=rc_step)
        reversed_window_source(reversed_record, rec)
        reversed_records.append(reversed_record)
    if reverse_flag and reversed_records:
        log.info("Reverse complemented %d record(s).", len(reversed_records))
    return reversed_records


__all__ = [
    "RecordSelector",
    "parse_record_selector",
    "parse_record_selectors",
    "drawn_record_units",
    "select_record",
    "select_record_units",
    "reverse_records",
]
