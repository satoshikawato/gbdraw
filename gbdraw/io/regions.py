#!/usr/bin/env python
# coding: utf-8

"""Region parsing and cropping helpers for linear mode."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
from pathlib import PurePath, PureWindowsPath
from typing import Sequence

from Bio.SeqRecord import SeqRecord

from ..core.numeric import ASCII_NONNEGATIVE_INTEGER_PATTERN, parse_ascii_nonnegative_integer
from ..core.record_metadata import _read_coord_map, _write_coord_map
from ..crop_genbank import check_start_end_coords, crop_and_shift_features
from .record_select import reverse_records

logger = logging.getLogger(__name__)

_REGION_COORD_RE = re.compile(
    rf"(?P<start>{ASCII_NONNEGATIVE_INTEGER_PATTERN})(?:\.\.|-)(?P<end>{ASCII_NONNEGATIVE_INTEGER_PATTERN})"
    r"(?::(?P<strand>rc|rev|reverse|minus|-))?$",
    re.IGNORECASE,
)

_REGION_APPLIED_KEY = "gbdraw_region_applied"






@dataclass(frozen=True)
class RegionSpec:
    raw: str
    file_selector: str | None
    record_id: str | None
    record_index: int | None
    start: int
    end: int
    reverse_complement: bool

    def selector_label(self) -> str:
        parts: list[str] = []
        if self.file_selector:
            parts.append(self.file_selector)
        if self.record_index is not None:
            parts.append(f"#{self.record_index + 1}")
        elif self.record_id:
            parts.append(self.record_id)
        if parts:
            return ":".join(parts)
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

    match = _REGION_COORD_RE.search(text)
    start = end = None
    if match and match.end() == len(text):
        start = parse_ascii_nonnegative_integer(match.group("start"))
        end = parse_ascii_nonnegative_integer(match.group("end"))
    if match is None or start is None or end is None:
        raise ValueError(
            "Invalid region spec: '{0}'. Expected format: record_id:start-end[:rc], "
            "#index:start-end[:rc], or file:record_selector:start-end[:rc].".format(text)
        )

    if start < 1 or end < 1:
        raise ValueError(f"Region coordinates must be >= 1: '{text}'.")

    reverse = False
    strand_raw = match.group("strand")
    if strand_raw:
        reverse = True

    if start > end:
        reverse = True
        start, end = end, start

    selector_text = text[: match.start()].rstrip(":")
    file_selector = None
    record_selector = None
    if selector_text:
        if ":" in selector_text:
            file_part, record_part = selector_text.rsplit(":", 1)
            file_part = file_part.strip()
            record_part = record_part.strip()
            if not record_part:
                raise ValueError(
                    f"File-scoped region spec '{text}' is missing a record selector."
                )
            file_selector = file_part or None
            record_selector = record_part
        else:
            record_selector = selector_text.strip()

    record_id = None
    record_index = None
    if record_selector:
        if record_selector.startswith("#"):
            index_str = record_selector[1:].strip()
            number = parse_ascii_nonnegative_integer(index_str)
            if number is None:
                raise ValueError(
                    f"Invalid record index in region spec '{text}'. Use #<number>:start-end."
                )
            record_index = number - 1
            if record_index < 0:
                raise ValueError(
                    f"Record index must be >= 1 in region spec '{text}'."
                )
        else:
            record_id = _extract_parenthetical(record_selector) or record_selector

    return RegionSpec(
        raw=text,
        file_selector=file_selector,
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
    coord_base, coord_step = _read_coord_map(record)

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
    if not record.annotations.get(_REGION_APPLIED_KEY):
        # Capture the complete source length before the first crop; later crops
        # retain this provenance through the copied annotations.
        new_record.annotations["gbdraw_source_length"] = len(record)

    if getattr(record, "letter_annotations", None):
        sliced = {}
        for key, values in record.letter_annotations.items():
            try:
                sliced[key] = values[start_0:end_0]
            except Exception:
                continue
        if sliced:
            new_record.letter_annotations = sliced

    new_record.features = crop_and_shift_features(
        record.features,
        start_0,
        end_0,
        coord_base=coord_base,
        coord_step=coord_step,
    )
    cropped_base = coord_base + (coord_step * start_0)
    cropped_step = coord_step
    _write_coord_map(new_record, base=cropped_base, step=cropped_step)

    if spec.reverse_complement:
        new_record = reverse_records((new_record,), True, log=log)[0]
        cropped_length = len(new_record.seq)
        rc_base = cropped_base + (cropped_step * max(0, cropped_length - 1))
        rc_step = -cropped_step
        _write_coord_map(new_record, base=rc_base, step=rc_step)

    if getattr(new_record, "annotations", None) is None:
        new_record.annotations = {}
    new_record.annotations[_REGION_APPLIED_KEY] = True

    return new_record


def assign_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
) -> dict[int, RegionSpec]:
    """Return the region spec that applies to each record index."""

    total = len(records)
    selectorless = [
        spec
        for spec in specs
        if spec.record_id is None and spec.record_index is None and spec.file_selector is None
    ]
    selectorful = [spec for spec in specs if spec not in selectorless]

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
        record_file_aliases: list[set[str]] = []
        for idx, rec in enumerate(records):
            id_to_indices.setdefault(rec.id, []).append(idx)
            aliases: set[str] = set()
            if getattr(rec, "annotations", None):
                raw = rec.annotations.get("gbdraw_source_file")
                if raw:
                    aliases.add(str(raw))
                raw_base = rec.annotations.get("gbdraw_source_basename")
                if raw_base:
                    aliases.add(str(raw_base))
            expanded: set[str] = set()
            for alias in aliases:
                expanded.add(alias)
                try:
                    expanded.add(PurePath(alias).name)
                except Exception:
                    pass
                try:
                    expanded.add(PureWindowsPath(alias).name)
                except Exception:
                    pass
            record_file_aliases.append({a for a in expanded if a})

        def _find_file_indices(file_selector: str) -> list[int]:
            matches: list[int] = []
            for idx, aliases in enumerate(record_file_aliases):
                if file_selector in aliases:
                    matches.append(idx)
            return matches

        for spec in selectorful:
            file_indices: list[int] | None = None
            if spec.file_selector:
                file_indices = _find_file_indices(spec.file_selector)
                if not file_indices:
                    available = sorted(
                        {alias for aliases in record_file_aliases for alias in aliases}
                    )
                    raise ValueError(
                        f"Region spec '{spec.raw}' did not match any input file. "
                        f"Available files: {', '.join(available) if available else '(unknown)'}."
                    )

            if spec.record_index is not None:
                if file_indices is None:
                    idx = spec.record_index
                    if idx < 0 or idx >= total:
                        raise ValueError(
                            f"Region spec '{spec.raw}' targets record index #{idx + 1}, "
                            f"but only {total} record(s) were loaded."
                        )
                else:
                    if spec.record_index < 0 or spec.record_index >= len(file_indices):
                        raise ValueError(
                            f"Region spec '{spec.raw}' targets record index #{spec.record_index + 1} "
                            f"in file '{spec.file_selector}', but only {len(file_indices)} record(s) were loaded."
                        )
                    idx = file_indices[spec.record_index]
                if idx in assignments:
                    raise ValueError(
                        f"Multiple region specs target record {spec.selector_label()}."
                    )
                assignments[idx] = spec
                continue

            record_id = spec.record_id or ""
            if file_indices is None:
                matches = id_to_indices.get(record_id, [])
            else:
                matches = [idx for idx in file_indices if records[idx].id == record_id]
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

    return assignments


def crop_record_to_region(
    record: SeqRecord,
    spec: RegionSpec,
    *,
    log: logging.Logger | None = None,
) -> SeqRecord:
    """Crop one record to ``spec``, then reverse-complement the crop when it says so."""

    log = log or logger
    log.info(
        "INFO: Cropping record %s to %s..%s%s.",
        record.id,
        spec.start,
        spec.end,
        " (reverse complement)" if spec.reverse_complement else "",
    )
    return _crop_record_to_region(record, spec, log=log)


def apply_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
    *,
    log: logging.Logger | None = None,
) -> list[SeqRecord]:
    if not specs:
        return list(records)
    assignments = assign_region_specs(records, specs)
    return [
        crop_record_to_region(record, assignments[idx], log=log)
        if idx in assignments
        else record
        for idx, record in enumerate(records)
    ]


__all__ = [
    "RegionSpec",
    "apply_region_specs",
    "assign_region_specs",
    "crop_record_to_region",
    "parse_region_spec",
    "parse_region_specs",
]
