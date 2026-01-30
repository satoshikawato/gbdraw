#!/usr/bin/env python
# coding: utf-8

"""Record selection and reverse-complement helpers for linear mode."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

logger = logging.getLogger(__name__)


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
        idx_text = raw[1:].strip()
        if not idx_text.isdigit():
            raise ValueError(f"Invalid record selector '{raw}'. Use #<number> or record_id.")
        idx = int(idx_text) - 1
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
        reversed_records.append(
            rec.reverse_complement(
                id=True,
                name=True,
                description=True,
                features=True,
                annotations=True,
                letter_annotations=True,
                dbxrefs=True,
            )
        )
    if reverse_flag and reversed_records:
        log.info("INFO: Reverse complemented %d record(s).", len(reversed_records))
    return reversed_records


__all__ = [
    "RecordSelector",
    "parse_record_selector",
    "parse_record_selectors",
    "select_record",
    "reverse_records",
]
