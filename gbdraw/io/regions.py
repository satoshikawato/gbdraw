#!/usr/bin/env python
# coding: utf-8

"""Region parsing and cropping helpers for linear mode."""

from __future__ import annotations

from dataclasses import dataclass, replace
import logging
import re
from pathlib import PurePath, PureWindowsPath
from typing import Sequence

from Bio.SeqFeature import CompoundLocation, ExactPosition, Location, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from ..analysis.windows import WINDOW_SOURCE_ANNOTATION, keep_window_source
from ..core.numeric import ASCII_NONNEGATIVE_INTEGER_PATTERN, parse_ascii_nonnegative_integer
from ..core.record_metadata import (
    _copy_source_feature_identity,
    _feature_source_index_map,
    _read_coord_map,
    _source_intervals,
    _source_period,
    _wrap_source_base,
    _write_coord_map,
)
from ..exceptions import ValidationError
from ..crop_genbank import check_start_end_coords, crop_and_shift_features
from ..layout.record_transform import detected_topology, drawn_circular
from .record_select import RecordSelector, record_index_number, reverse_records
from .record_sequences import sequence_view, two_turns

logger = logging.getLogger(__name__)

_REGION_COORD_RE = re.compile(
    rf"(?P<start>{ASCII_NONNEGATIVE_INTEGER_PATTERN})(?:\.\.|-)(?P<end>{ASCII_NONNEGATIVE_INTEGER_PATTERN})"
    r"(?P<origin>:origin)?(?::(?P<strand>rc|rev|reverse|minus|-))?$",
    re.IGNORECASE,
)

_REGION_APPLIED_KEY = "gbdraw_region_applied"


def region_text(start: int, end: int, reverse_complement: bool) -> str:
    """The region string of a span without a record selector: the one writer of the ``START-END[:origin][:rc]`` form.

    A span whose ``start`` is after its ``end`` runs across the origin and is
    written with ``:origin`` (PD-OI-122); without it the larger number first
    reads as the reversed span.
    """

    return f"{start}-{end}" + (":origin" if start > end else "") + (":rc" if reverse_complement else "")


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

    def record_selector(self) -> RecordSelector | None:
        """The record this spec names (by ID or index), or ``None``; the file selector is not a record's."""

        if self.record_id is not None:
            return RecordSelector(raw=self.record_id, record_id=self.record_id, record_index=None)
        if self.record_index is not None:
            return RecordSelector(raw=f"#{self.record_index + 1}", record_id=None, record_index=self.record_index)
        return None


@dataclass(frozen=True)
class RecordWindow:
    """One span of a record, drawn as its own unit of a Linear row (D-01).

    ``start`` and ``end`` are 1-based source coordinates of the record;
    ``reverse_complement`` draws the span reversed; ``column`` orders the
    record's units with the other units of its row. Windows of one record
    leave out the interval between them. A window whose ``start`` is after its
    ``end`` runs from ``start`` across the origin to ``end``; only a circular
    record has one (D-28), and the crop refuses it on a linear record.
    """

    start: int
    end: int
    reverse_complement: bool = False
    column: int | None = None

    def __post_init__(self) -> None:
        for name in ("start", "end"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"Window {name} must be a positive integer.")
        if not isinstance(self.reverse_complement, bool):
            raise ValueError("Window reverse_complement must be a boolean.")
        if self.column is not None and (
            not isinstance(self.column, int) or isinstance(self.column, bool) or self.column < 1
        ):
            raise ValueError("Window column must be a positive integer or None.")

    @property
    def crosses_origin(self) -> bool:
        """Whether the window runs across the origin of its (circular) record."""
        return self.start > self.end

    @classmethod
    def from_region(cls, spec: RegionSpec, *, column: int | None = None) -> RecordWindow:
        """The window ``spec`` crops; its record selector stays with the spec's reader."""

        return cls(spec.start, spec.end, spec.reverse_complement, column)

    def region(self) -> RegionSpec:
        """The crop this window draws, as an unqualified region spec.

        The spec of a window across the origin keeps ``start`` after ``end``, and
        its ``raw`` text says ``:origin`` (:func:`region_text`).
        """

        return RegionSpec(
            raw=region_text(self.start, self.end, self.reverse_complement),
            file_selector=None,
            record_id=None,
            record_index=None,
            start=self.start,
            end=self.end,
            reverse_complement=self.reverse_complement,
        )


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
            "Invalid region spec: '{0}'. Expected format: record_id:start-end[:origin][:rc], "
            "#index:start-end[:origin][:rc], or file:record_selector:start-end[:origin][:rc]; "
            ":origin marks a span across the origin of a circular record.".format(text)
        )

    if start < 1 or end < 1:
        raise ValueError(f"Region coordinates must be >= 1: '{text}'.")

    reverse = False
    strand_raw = match.group("strand")
    if strand_raw:
        reverse = True

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
            number = record_index_number(record_selector)
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

    if match.group("origin"):
        # A span across the origin keeps its start after its end (D-28); the crop checks the record is circular.
        if start <= end:
            raise ValueError(
                f"Region '{text}' is marked :origin, but its start ({start}) is not greater than its end ({end}); "
                "a span across the origin starts after it ends."
            )
    elif start > end:
        # Without :origin the larger number first is the reversed span, as in v0.14.5.
        reverse = True
        start, end = end, start

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


def region_slice(record_length: int, start: int, end: int, record_id: str) -> tuple[int, int]:
    """The 0-based slice a region of a record draws, its end clamped as the crop clamps it.

    A region whose start is after its end runs across the origin (D-28): its
    slice is on the record read twice round, ``start - 1`` to ``length + end``.
    It reads only the record's length, so the crop and the Web export check of
    a drawn record's bases share it (D-70); whether the record is circular is
    the crop's to check.
    """

    if start <= end:
        return check_start_end_coords(record_length, start, end, record_id)
    if not 1 <= end < start <= record_length:
        raise ValidationError(
            f"Region {start}..{end} across the origin lies outside record {record_id} ({record_length} bp). "
            f"Choose a start and an end within 1..{record_length}.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
        )
    return start - 1, record_length + end


def _origin_crossing_view(record: SeqRecord, spec: RegionSpec, *, circular: bool) -> tuple[SeqRecord, int, int]:
    """The record read twice round, and the 0-based bounds of ``spec`` on it, for a span across the origin (D-28).

    Each feature appears once on each turn; a feature across the origin is one
    continuous location on the turn it starts on, so a window keeps it whole.
    The features keep the identity of the source feature they copy.
    """

    length = len(record.seq)
    if not circular:
        raise ValidationError(
            f"Region {spec.start}..{spec.end} of record {record.id} crosses the origin of a linear record; "
            "only a circular record has a span across its origin.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "START_AFTER_END"},
        )
    if record.annotations.get(_REGION_APPLIED_KEY) or _read_coord_map(record) != (1, 1):
        raise ValidationError(
            f"Region {spec.start}..{spec.end} of record {record.id} crosses the origin of a record that is "
            "already cropped or reversed.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "START_AFTER_END"},
        )
    start_0, end_0 = region_slice(length, spec.start, spec.end, str(record.id))
    source_indexes = _feature_source_index_map(record.features)
    turns: list[SeqFeature] = []
    for feature in record.features:
        index = source_indexes[id(feature)]
        shown = _two_turns(_unwrapped_location(feature.location, length), length, (start_0, end_0))
        if shown is None:
            continue
        copy = SeqFeature(
            shown,
            type=feature.type,
            qualifiers=feature.qualifiers,
            id=feature.id,
        )
        _copy_source_feature_identity(feature, copy, index=index, coord_base=1, coord_step=1, period=None)
        source_indexes[id(copy)] = index
        turns.append(copy)
    view = SeqRecord(
        two_turns(record.seq),
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=list(getattr(record, "dbxrefs", []) or []),
        annotations=dict(getattr(record, "annotations", {}) or {}),
    )
    view.features = turns
    if getattr(record, "letter_annotations", None):
        try:
            view.letter_annotations = {key: values + values for key, values in record.letter_annotations.items()}
        except (TypeError, ValueError):
            pass
    return view, start_0, end_0


def _unwrapped_location(location: Location, length: int) -> Location:
    """``location`` read as one pass round a circle of ``length`` bases: the parts past the origin move one turn on.

    Parts are in biological order. On the plus strand a part that starts before
    the one read just before it lies past the origin, and it and the parts after
    it move by ``length``; on the minus strand, read downwards, the parts read
    before such a jump move instead.
    """

    parts = list(getattr(location, "parts", None) or [location])
    minus = all(part.strand == -1 for part in parts)
    crossed = next(
        (
            index
            for index in range(1, len(parts))
            if (parts[index].start > parts[index - 1].start if minus else parts[index].start < parts[index - 1].start)
        ),
        None,
    )
    if crossed is None:
        return location
    before, after = parts[:crossed], parts[crossed:]
    moved = [part + length for part in before] + after if minus else before + [part + length for part in after]
    # The two parts that meet at the origin are one run of bases: join them, so the feature draws as one part there.
    last, first = moved[crossed - 1], moved[crossed]
    low, high = (first, last) if minus else (last, first)
    if low.end == high.start and low.strand == high.strand:
        moved[crossed - 1 : crossed + 1] = [SimpleLocation(low.start, high.end, strand=low.strand)]
    if len(moved) == 1:
        return moved[0]
    return CompoundLocation(moved, operator=getattr(location, "operator", "join"))


def _wrapped_location(location: Location, length: int) -> Location:
    """``location``, read on a circle of ``length`` bases from its origin onwards for up to two turns, on one turn:
    the inverse of ``_unwrapped_location``. A part past the origin moves one turn back, and a part across it is split
    there (``_source_intervals``) into a part that ends at the origin and one that starts after it, in biological
    order. The cut points at the origin are exact; the part's own ends keep their positions (a partial end stays
    partial)."""

    parts: list[Location] = []
    for part in getattr(location, "parts", None) or [location]:
        start = int(part.start)
        intervals = _source_intervals(start, int(part.end), 1, 1, length)
        if len(intervals) == 1:
            parts.append(part if intervals[0][0] == start else part._shift(intervals[0][0] - start))
            continue
        # A position plus an offset keeps its kind (minus gives a plain int).
        pieces = [
            SimpleLocation(part.start, ExactPosition(length), strand=part.strand),
            SimpleLocation(ExactPosition(0), part.end + (-length), strand=part.strand),
        ]
        parts.extend(reversed(pieces) if part.strand == -1 else pieces)
    if len(parts) == 1:
        return parts[0]
    return CompoundLocation(parts, operator=getattr(location, "operator", None) or "join")


def _two_turns(location: Location, length: int, window: tuple[int, int]) -> Location | None:
    """``location`` on the turns of the record read twice round that ``window`` shows, in biological order.

    A window across the origin shows at most one turn of each base, so a feature
    it shows on both turns (one that runs over the left-out arc) is one feature
    whose parts lie at both ends of the window; parts that meet are joined. A
    turn the window does not reach is left out whole, so it marks no cut end.
    """

    parts = list(getattr(location, "parts", None) or [location])
    low, high = window
    shown = [
        turn
        for turn in ([*parts], [part + length for part in parts])
        if any(part.start < high and part.end > low for part in turn)
    ]
    if all(part.strand == -1 for part in parts):
        shown.reverse()
    turns = [part for turn in shown for part in turn]
    if not turns:
        return None
    joined: list[Location] = []
    for part in turns:
        previous = joined[-1] if joined else None
        if previous is not None and previous.strand == part.strand and (
            previous.end == part.start or part.end == previous.start
        ):
            joined[-1] = SimpleLocation(min(previous.start, part.start), max(previous.end, part.end), strand=part.strand)
        else:
            joined.append(part)
    return joined[0] if len(joined) == 1 else CompoundLocation(joined, operator=getattr(location, "operator", "join"))


def _crop_record_to_region(
    record: SeqRecord,
    spec: RegionSpec,
    *,
    log: logging.Logger,
    circular: bool = False,
) -> SeqRecord:
    source_length = len(record.seq)
    crosses_origin = spec.start > spec.end
    source = record
    if crosses_origin:
        source, start_0, end_0 = _origin_crossing_view(record, spec, circular=circular)
    else:
        start_0, end_0 = region_slice(len(record.seq), spec.start, spec.end, str(record.id))
    coord_base, coord_step = _read_coord_map(record)

    if not crosses_origin and (start_0 != spec.start - 1 or end_0 != spec.end):
        log.warning(
            "Region %s for record %s was clamped to %s..%s (record length %s).",
            spec.raw,
            record.id,
            start_0 + 1,
            end_0,
            len(record.seq),
        )

    # The crop reads its bases from the record it is cut from (a view, not a copy).
    new_record = SeqRecord(
        sequence_view(source.seq, start_0, end_0),
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
    # Set once the crop is complete (below).
    new_record.annotations.pop(WINDOW_SOURCE_ANNOTATION, None)

    if getattr(source, "letter_annotations", None):
        sliced = {}
        for key, values in source.letter_annotations.items():
            try:
                sliced[key] = values[start_0:end_0]
            except Exception:
                continue
        if sliced:
            new_record.letter_annotations = sliced

    # A crop across the origin, or a crop of one, reads source bases round the circle (``_wrap_source_base``): the
    # base of its map is a source base, also where it lies past the source end (OV-3270).
    period = source_length if crosses_origin else _source_period(record)
    new_record.features = crop_and_shift_features(
        source.features,
        start_0,
        end_0,
        coord_base=coord_base,
        coord_step=coord_step,
        period=_source_period(record),
    )
    cropped_base = _wrap_source_base(coord_base + (coord_step * start_0), period)
    cropped_step = coord_step
    _write_coord_map(new_record, base=cropped_base, step=cropped_step)

    if spec.reverse_complement:
        new_record = reverse_records((new_record,), True, log=log)[0]
        cropped_length = len(new_record.seq)
        rc_base = _wrap_source_base(cropped_base + (cropped_step * max(0, cropped_length - 1)), period)
        rc_step = -cropped_step
        _write_coord_map(new_record, base=rc_base, step=rc_step)

    if getattr(new_record, "annotations", None) is None:
        new_record.annotations = {}
    new_record.annotations[_REGION_APPLIED_KEY] = True
    # The crop's windows read their flanking bases from the record it is cut from (OV-2890).
    keep_window_source(new_record, record)

    return new_record


def assign_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
) -> dict[int, tuple[RegionSpec, ...]]:
    """Return the region specs that apply to each record index, in the order given.

    Several specs that select one record are its windows (D-01).
    """

    total = len(records)
    selectorless = [
        spec
        for spec in specs
        if spec.record_id is None and spec.record_index is None and spec.file_selector is None
    ]
    selectorful = [spec for spec in specs if spec not in selectorless]

    assignments: dict[int, tuple[RegionSpec, ...]] = {}

    if selectorless:
        if selectorful:
            raise ValueError(
                "Region specs without a record selector cannot be mixed with selector-based specs."
            )
        if total == 1:
            assignments[0] = tuple(selectorless)
        elif len(specs) == total:
            for idx, spec in enumerate(specs):
                assignments[idx] = (spec,)
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
                assignments[idx] = (*assignments.get(idx, ()), spec)
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
            assignments[idx] = (*assignments.get(idx, ()), spec)

    return assignments


def crop_record_to_region(
    record: SeqRecord,
    spec: RegionSpec,
    *,
    log: logging.Logger | None = None,
    circular: bool = False,
) -> SeqRecord:
    """Crop one record to ``spec``, then reverse-complement the crop when it says so.

    A spec whose start is after its end runs across the origin; ``circular``
    says whether the record is drawn as circular, which such a spec requires.
    """

    log = log or logger
    log.info(
        "Cropping record %s to %s..%s%s.",
        record.id,
        spec.start,
        spec.end,
        " (reverse complement)" if spec.reverse_complement else "",
    )
    return _crop_record_to_region(record, spec, log=log, circular=circular)


def apply_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
    *,
    log: logging.Logger | None = None,
) -> list[SeqRecord]:
    """Crop each record to the specs that select it; a record with several specs gives one record per window."""

    if not specs:
        return list(records)
    windows, _cropped = crop_then_reverse_records(
        records, specs, [False] * len(records), is_circular=None, source_kind="memory", log=log
    )
    return [unit for units in windows for unit in units]


def window_drawn_reverse(window_reverse: bool, record_reverse: bool) -> bool:
    """Whether a window is drawn reverse-complemented: by its own flag or its record's, both asking for one reverse
    (the rule :func:`crop_record_windows` draws by)."""

    return window_reverse or record_reverse


def crop_record_windows(
    record: SeqRecord,
    specs: Sequence[RegionSpec],
    reverse: bool,
    *,
    log: logging.Logger | None = None,
    circular: bool = False,
) -> tuple[SeqRecord, ...]:
    """Crop ``record`` to each of its windows, then reverse-complement each (OV-1279).

    The one owner of that order: a record's ``reverse`` flag applies after its
    crop, as the spec's own ``:rc`` does, so both spellings draw one window and
    asking for both reverses once. A record without a window is only reversed.
    Windows of one record may overlap (D-18); a record with several windows
    orients each window by its own flag (D-01). A window across the origin
    (start after end) needs ``circular``, the record drawn as circular (D-28);
    it is cropped from the record read round the origin, so its features stay
    whole.
    """

    log = log or logger
    if not specs:
        return (reverse_records((record,), reverse, log=log)[0],)
    if len(specs) > 1:
        if reverse:
            raise ValueError(
                f"Record {record.id} has several windows; reverse-complement each window, not the record."
            )
    return tuple(
        crop_record_to_region(
            record, replace(spec, reverse_complement=window_drawn_reverse(spec.reverse_complement, reverse)), log=log,
            circular=circular,
        )
        for spec in specs
    )


def crop_then_reverse_records(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
    reverse: Sequence[bool],
    *,
    is_circular: bool | None,
    source_kind: str,
    log: logging.Logger | None = None,
) -> tuple[list[tuple[SeqRecord, ...]], frozenset[int]]:
    """Crop each record to the windows the specs select, then orient them (:func:`crop_record_windows`).

    Returns each record's drawn windows (one record when no spec selects it)
    and the indexes of the cropped records. A window across the origin
    (``:origin``) needs the record drawn as circular, decided as the drawing
    decides it (:func:`drawn_circular`): by its display option ``is_circular``,
    else by the topology its file of ``source_kind`` states.
    """

    assignments = assign_region_specs(records, specs) if specs else {}
    oriented = [
        crop_record_windows(
            record,
            assignments.get(index, ()),
            flag,
            log=log,
            circular=drawn_circular(is_circular, detected_topology(record, source_kind)),
        )
        for index, (record, flag) in enumerate(zip(records, reverse, strict=True))
    ]
    return oriented, frozenset(assignments)


__all__ = [
    "RecordWindow",
    "RegionSpec",
    "apply_region_specs",
    "assign_region_specs",
    "crop_record_to_region",
    "crop_record_windows",
    "crop_then_reverse_records",
    "parse_region_spec",
    "parse_region_specs",
    "region_text",
    "window_drawn_reverse",
]
