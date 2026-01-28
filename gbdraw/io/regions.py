#!/usr/bin/env python
# coding: utf-8

"""Region parsing and cropping helpers for linear mode."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from ..crop_genbank import check_start_end_coords, crop_and_shift_features

logger = logging.getLogger(__name__)

_REGION_RE = re.compile(
    r"^(?:(?P<selector>[^:]+):)?(?P<start>\d+)(?:\.\.|-)(?P<end>\d+)(?::(?P<strand>rc|rev|reverse|minus|-))?$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class RegionSpec:
    raw: str
    record_id: str | None
    record_index: int | None
    start: int
    end: int
    reverse_complement: bool

    def selector_label(self) -> str:
        if self.record_index is not None:
            return f"#{self.record_index + 1}"
        if self.record_id:
            return self.record_id
        return "(by-order)"


def _extract_parenthetical(selector: str) -> str | None:
    match = re.search(r"\(([^()]+)\)\s*$", selector)
    if not match:
        return None
    inner = match.group(1).strip()
    return inner or None


def parse_region_spec(spec: str) -> RegionSpec:
    if spec is None:
        raise ValueError("Region spec is missing.")
    text = str(spec).strip()
    if not text:
        raise ValueError("Region spec is empty.")

    match = _REGION_RE.match(text)
    if not match:
        raise ValueError(
            "Invalid region spec: '{0}'. Expected format: record_id:start-end[:rc] or #index:start-end[:rc].".format(
                text
            )
        )

    selector_raw = match.group("selector")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start < 1 or end < 1:
        raise ValueError(f"Region coordinates must be >= 1: '{text}'.")

    reverse = False
    strand_raw = match.group("strand")
    if strand_raw:
        reverse = True

    if start > end:
        reverse = True
        start, end = end, start

    record_id = None
    record_index = None

    if selector_raw:
        selector_raw = selector_raw.strip()
        if selector_raw.startswith("#"):
            index_str = selector_raw[1:]
            if not index_str.isdigit():
                raise ValueError(
                    f"Invalid record index in region spec '{text}'. Use #<number>:start-end."
                )
            record_index = int(index_str) - 1
            if record_index < 0:
                raise ValueError(
                    f"Record index must be >= 1 in region spec '{text}'."
                )
        else:
            record_id = _extract_parenthetical(selector_raw) or selector_raw

    return RegionSpec(
        raw=text,
        record_id=record_id,
        record_index=record_index,
        start=start,
        end=end,
        reverse_complement=reverse,
    )


def parse_region_specs(specs: Sequence[str] | None) -> list[RegionSpec]:
    if not specs:
        return []
    parsed: list[RegionSpec] = []
    for spec in specs:
        if spec is None:
            continue
        parsed.append(parse_region_spec(spec))
    return parsed


def _crop_record_to_region(
    record: SeqRecord,
    spec: RegionSpec,
    *,
    log: logging.Logger,
) -> SeqRecord:
    start_0, end_0 = check_start_end_coords(record, spec.start, spec.end)

    if start_0 != spec.start - 1 or end_0 != spec.end:
        log.warning(
            "WARNING: Region %s for record %s was clamped to %s..%s (record length %s).",
            spec.raw,
            record.id,
            start_0 + 1,
            end_0,
            len(record.seq),
        )

    new_seq = record.seq[start_0:end_0]
    new_record = SeqRecord(
        new_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=list(getattr(record, "dbxrefs", []) or []),
        annotations=dict(getattr(record, "annotations", {}) or {}),
    )

    if getattr(record, "letter_annotations", None):
        sliced = {}
        for key, values in record.letter_annotations.items():
            try:
                sliced[key] = values[start_0:end_0]
            except Exception:
                continue
        if sliced:
            new_record.letter_annotations = sliced

    new_record.features = crop_and_shift_features(record.features, start_0, end_0)

    if spec.reverse_complement:
        new_record = new_record.reverse_complement(
            id=True,
            name=True,
            description=True,
            features=True,
            annotations=True,
            letter_annotations=True,
            dbxrefs=True,
        )

    return new_record


def apply_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
    *,
    log: logging.Logger | None = None,
) -> list[SeqRecord]:
    if not specs:
        return list(records)
    log = log or logger

    total = len(records)
    selectorless = [spec for spec in specs if spec.record_id is None and spec.record_index is None]
    selectorful = [spec for spec in specs if spec.record_id is not None or spec.record_index is not None]

    assignments: dict[int, RegionSpec] = {}

    if selectorless:
        if selectorful:
            raise ValueError(
                "Region specs without a record selector cannot be mixed with selector-based specs."
            )
        if total == 1 and len(specs) == 1:
            assignments[0] = selectorless[0]
        elif len(specs) == total:
            for idx, spec in enumerate(specs):
                assignments[idx] = spec
        else:
            raise ValueError(
                "Region specs without a record selector require a single record or one spec per record in input order."
            )
    else:
        id_to_indices: dict[str, list[int]] = {}
        for idx, rec in enumerate(records):
            id_to_indices.setdefault(rec.id, []).append(idx)

        for spec in selectorful:
            if spec.record_index is not None:
                idx = spec.record_index
                if idx < 0 or idx >= total:
                    raise ValueError(
                        f"Region spec '{spec.raw}' targets record index #{idx + 1}, but only {total} record(s) were loaded."
                    )
                if idx in assignments:
                    raise ValueError(
                        f"Multiple region specs target record #{idx + 1}."
                    )
                assignments[idx] = spec
                continue

            record_id = spec.record_id or ""
            matches = id_to_indices.get(record_id, [])
            if not matches:
                raise ValueError(
                    f"Region spec '{spec.raw}' did not match any record ID."
                )
            if len(matches) > 1:
                raise ValueError(
                    f"Region spec '{spec.raw}' matched multiple records with ID '{record_id}'. Use #index to disambiguate."
                )
            idx = matches[0]
            if idx in assignments:
                raise ValueError(
                    f"Multiple region specs target record '{record_id}'."
                )
            assignments[idx] = spec

    new_records: list[SeqRecord] = []
    for idx, record in enumerate(records):
        spec = assignments.get(idx)
        if not spec:
            new_records.append(record)
            continue
        log.info(
            "INFO: Cropping record %s to %s..%s%s.",
            record.id,
            spec.start,
            spec.end,
            " (reverse complement)" if spec.reverse_complement else "",
        )
        new_records.append(_crop_record_to_region(record, spec, log=log))

    return new_records


__all__ = [
    "RegionSpec",
    "apply_region_specs",
    "parse_region_spec",
    "parse_region_specs",
]
