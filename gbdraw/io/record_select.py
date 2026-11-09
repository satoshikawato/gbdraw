#!/usr/bin/env python
# coding: utf-8

"""Record selection and reverse-complement helpers for linear mode."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
from typing import Sequence

from Bio.SeqRecord import SeqRecord

from ..analysis.windows import reversed_window_source
from ..core.numeric import parse_ascii_nonnegative_integer
from ..core.record_metadata import (
    _copy_source_feature_identity,
    _feature_source_index_map,
    _read_coord_map,
    _write_coord_map,
    record_identity_key,
)

logger = logging.getLogger(__name__)

# A ``#N`` record selector, as ``parse_record_selector`` reads it: spaces (``str.strip``'s, which are ``\s``'s)
# around it and after ``#``, and N in ASCII digits, leading zeros allowed (``#01`` is record 1). The Web reads
# the same pattern from the generated contracts (``RECORD_INDEX_SELECTOR_PATTERN``).
RECORD_INDEX_SELECTOR_PATTERN = r"\s*#\s*([0-9]+)\s*"
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

        def feature_start(feature):
            try:
                return int(feature.location.start)
            except TypeError:
                return None

        def reverse_feature(feature):
            reversed_feature = feature._flip(record_length)
            _copy_source_feature_identity(
                feature,
                reversed_feature,
                fallback_index=source_indexes[id(feature)],
                coord_base=base,
                coord_step=step,
            )
            sub_features = getattr(feature, "sub_features", None)
            if sub_features is not None:
                reversed_sub_features = [
                    reverse_feature(sub_feature)
                    for sub_feature in sub_features
                ]
                reversed_sub_features.sort(key=feature_start)
                reversed_feature.sub_features = reversed_sub_features
            return reversed_feature

        reversed_features = [
            reverse_feature(feature)
            for feature in rec.features
        ]
        reversed_features.sort(key=feature_start)
        reversed_record = rec.reverse_complement(
            id=True,
            name=True,
            description=True,
            features=reversed_features,
            annotations=True,
            letter_annotations=True,
            dbxrefs=True,
        )
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
