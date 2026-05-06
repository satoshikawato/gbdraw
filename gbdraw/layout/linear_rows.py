#!/usr/bin/env python
# coding: utf-8

"""Row-aware layout primitives for linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass, replace
import copy
import re
from pathlib import Path
from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from ..canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..exceptions import ValidationError


@dataclass(frozen=True)
class LinearInputRow:
    row_index: int
    source_id: str
    source_path: str
    label: str
    records: tuple[SeqRecord, ...]


@dataclass(frozen=True)
class LinearRecordRef:
    row_index: int
    record_index: int
    source_id: str
    record_id: str
    unique_key: str
    label: str


@dataclass(frozen=True)
class LinearRecordLayout:
    ref: LinearRecordRef
    x: float
    width: float
    scale: float
    record_length: int


@dataclass(frozen=True)
class LinearRowLayout:
    row_index: int
    y_axis: float
    x_origin: float
    width: float
    records: tuple[LinearRecordLayout, ...]


@dataclass(frozen=True)
class LinearInterval:
    record: LinearRecordLayout
    start: float
    end: float
    start_x: float
    end_x: float


_SOURCE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def make_source_id(source_path: str, existing: set[str] | None = None, fallback_index: int = 0) -> str:
    """Create a deterministic, comparison-friendly source identifier."""

    existing = existing if existing is not None else set()
    stem = Path(str(source_path)).stem.strip() or f"source{fallback_index + 1}"
    base = _SOURCE_ID_RE.sub("_", stem).strip("._-") or f"source{fallback_index + 1}"
    source_id = base
    suffix = 2
    while source_id in existing:
        source_id = f"{base}_{suffix}"
        suffix += 1
    existing.add(source_id)
    return source_id


def make_linear_input_row(
    *,
    row_index: int,
    source_id: str,
    source_path: str,
    label: str | None,
    records: Sequence[SeqRecord],
) -> LinearInputRow:
    row_records = tuple(records)
    if not row_records:
        raise ValidationError(f"Linear input row {row_index + 1} contains no records.")
    row_label = str(label or "").strip() or Path(str(source_path)).name or source_id
    return LinearInputRow(
        row_index=int(row_index),
        source_id=str(source_id),
        source_path=str(source_path),
        label=row_label,
        records=row_records,
    )


def flatten_linear_rows(rows: Sequence[LinearInputRow]) -> list[SeqRecord]:
    return [record for row in rows for record in row.records]


def linear_row_record_index_groups(
    rows: Sequence[LinearInputRow],
) -> tuple[list[SeqRecord], list[list[int]]]:
    """Flatten row records and return the flat record indices contained in each row."""

    flat_records: list[SeqRecord] = []
    index_groups: list[list[int]] = []
    for row in rows:
        row_indices: list[int] = []
        for record in row.records:
            row_indices.append(len(flat_records))
            flat_records.append(record)
        index_groups.append(row_indices)
    return flat_records, index_groups


def record_index_to_linear_row_index(rows: Sequence[LinearInputRow]) -> dict[int, int]:
    """Return a flat-record-index to source-row-index lookup."""

    _flat_records, index_groups = linear_row_record_index_groups(rows)
    lookup: dict[int, int] = {}
    for row_pair_index, record_indices in enumerate(index_groups):
        for record_index in record_indices:
            lookup[int(record_index)] = int(row_pair_index)
    return lookup


def adjacent_record_pair_to_linear_row_pair(
    rows: Sequence[LinearInputRow],
) -> dict[tuple[int, int], int]:
    """Map flat record pairs between adjacent source rows to their row-pair index."""

    _flat_records, index_groups = linear_row_record_index_groups(rows)
    pair_index_by_record_pair: dict[tuple[int, int], int] = {}
    for row_pair_index in range(max(0, len(index_groups) - 1)):
        query_indices = index_groups[row_pair_index]
        subject_indices = index_groups[row_pair_index + 1]
        for query_index in query_indices:
            for subject_index in subject_indices:
                pair_index_by_record_pair[(int(query_index), int(subject_index))] = int(row_pair_index)
    return pair_index_by_record_pair


def _record_label(record: SeqRecord) -> str:
    annotations = getattr(record, "annotations", None) or {}
    for key in (
        "gbdraw_record_label",
        "gbdraw_segment_label",
        "gbdraw_accession_label",
        "gbdraw_original_record_id",
    ):
        label = annotations.get(key)
        if label is not None and str(label).strip():
            return str(label).strip()
    return str(record.id)


def clone_rows_for_rendering(rows: Sequence[LinearInputRow]) -> tuple[list[LinearInputRow], dict[str, str]]:
    """Return rows whose record IDs are globally unique, plus unique->original ID mapping."""

    render_rows: list[LinearInputRow] = []
    id_map: dict[str, str] = {}
    for row in rows:
        render_records: list[SeqRecord] = []
        for record_index, record in enumerate(row.records):
            original_id = str(record.id)
            unique_key = f"{row.source_id}:{record_index + 1}:{original_id}"
            render_record = SeqRecord(
                record.seq,
                id=unique_key,
                name=str(getattr(record, "name", "") or original_id),
                description=str(getattr(record, "description", "") or original_id),
                dbxrefs=list(getattr(record, "dbxrefs", []) or []),
                annotations=dict(getattr(record, "annotations", {}) or {}),
            )
            render_record.features = list(getattr(record, "features", []) or [])
            if getattr(record, "letter_annotations", None):
                render_record.letter_annotations = dict(record.letter_annotations)
            render_record.annotations["gbdraw_original_record_id"] = original_id
            render_record.annotations["gbdraw_source_id"] = row.source_id
            render_record.annotations["gbdraw_record_unique_key"] = unique_key
            render_record.annotations.setdefault("gbdraw_accession_label", original_id)
            render_record.annotations.setdefault("gbdraw_segment_label", _record_label(record))
            id_map[unique_key] = original_id
            render_records.append(render_record)
        render_rows.append(
            LinearInputRow(
                row_index=row.row_index,
                source_id=row.source_id,
                source_path=row.source_path,
                label=row.label,
                records=tuple(render_records),
            )
        )
    return render_rows, id_map


def build_linear_layout_index(
    rows: Sequence[LinearInputRow],
    *,
    alignment_width: float,
    record_gap: float,
    normalize_length: bool,
    align_center: bool,
    y_axes: Sequence[float] | None = None,
) -> "LinearLayoutIndex":
    if not rows:
        raise ValidationError("rows is empty")

    normalized_gap = max(0.0, float(record_gap))
    row_lengths = [sum(len(record.seq) for record in row.records) for row in rows]
    if any(length <= 0 for length in row_lengths):
        raise ValidationError("Linear rows must contain records with positive sequence lengths.")

    max_gap_total = max((max(0, len(row.records) - 1) * normalized_gap) for row in rows)
    if normalize_length:
        global_scale = None
    else:
        available_bp_width = max(1.0, float(alignment_width) - max_gap_total)
        global_scale = available_bp_width / float(max(row_lengths))

    layouts: list[LinearRowLayout] = []
    for row_pos, row in enumerate(rows):
        gap_total = max(0, len(row.records) - 1) * normalized_gap
        row_total_length = row_lengths[row_pos]
        if normalize_length:
            bp_width = max(1.0, float(alignment_width) - gap_total)
            scale = bp_width / float(row_total_length)
        else:
            scale = float(global_scale or 0.0)

        record_layouts: list[LinearRecordLayout] = []
        row_width = sum(len(record.seq) * scale for record in row.records) + gap_total
        x_origin = max(0.0, (float(alignment_width) - row_width) / 2.0) if align_center else 0.0
        current_x = x_origin
        for record_index, record in enumerate(row.records):
            record_length = len(record.seq)
            record_id = str(
                (getattr(record, "annotations", None) or {}).get(
                    "gbdraw_original_record_id",
                    record.id,
                )
            )
            unique_key = str(
                (getattr(record, "annotations", None) or {}).get(
                    "gbdraw_record_unique_key",
                    f"{row.source_id}:{record_index + 1}:{record_id}",
                )
            )
            ref = LinearRecordRef(
                row_index=row.row_index,
                record_index=record_index,
                source_id=row.source_id,
                record_id=record_id,
                unique_key=unique_key,
                label=_record_label(record),
            )
            width = record_length * scale
            record_layouts.append(
                LinearRecordLayout(
                    ref=ref,
                    x=current_x,
                    width=width,
                    scale=scale,
                    record_length=record_length,
                )
            )
            current_x += width + normalized_gap
        layouts.append(
            LinearRowLayout(
                row_index=row.row_index,
                y_axis=float(y_axes[row_pos]) if y_axes is not None else 0.0,
                x_origin=x_origin,
                width=row_width,
                records=tuple(record_layouts),
            )
        )
    return LinearLayoutIndex(tuple(layouts))


class LinearLayoutIndex:
    """Coordinate resolver for row-aware linear rendering."""

    def __init__(self, rows: tuple[LinearRowLayout, ...]) -> None:
        self.rows = rows
        self._by_row_index = {row.row_index: row for row in rows}

    def with_y_axes(self, y_axes: Sequence[float]) -> "LinearLayoutIndex":
        if len(y_axes) != len(self.rows):
            raise ValidationError(
                f"Expected {len(self.rows)} row y positions; got {len(y_axes)}."
            )
        return LinearLayoutIndex(
            tuple(
                LinearRowLayout(
                    row_index=row.row_index,
                    y_axis=float(y_axes[index]),
                    x_origin=row.x_origin,
                    width=row.width,
                    records=row.records,
                )
                for index, row in enumerate(self.rows)
            )
        )

    def resolve_row(self, row_index: int) -> LinearRowLayout:
        try:
            return self._by_row_index[int(row_index)]
        except KeyError as exc:
            raise ValidationError(f"Unknown linear row index {row_index}.") from exc

    def resolve_record(self, row_index: int, record_id_or_key: object) -> LinearRecordLayout:
        row = self.resolve_row(row_index)
        query = str(record_id_or_key or "").strip()
        if not query:
            raise ValidationError(f"Missing comparison record ID for row {row_index + 1}.")

        exact = [record for record in row.records if record.ref.unique_key == query]
        if len(exact) == 1:
            return exact[0]

        qualified = [record for record in row.records if f"{record.ref.source_id}:{record.ref.record_id}" == query]
        if len(qualified) == 1:
            return qualified[0]
        if len(qualified) > 1:
            raise ValidationError(
                f"Comparison record ID '{query}' is ambiguous in row {row_index + 1}; use a unique record key."
            )

        raw = [record for record in row.records if record.ref.record_id == query]
        if len(raw) == 1:
            return raw[0]
        if len(raw) > 1:
            raise ValidationError(
                f"Comparison record ID '{query}' is ambiguous in row {row_index + 1}; use a unique record key."
            )

        available = ", ".join(record.ref.record_id for record in row.records)
        raise ValidationError(
            f"Comparison record ID '{query}' was not found in row {row_index + 1}. "
            f"Available IDs: {available or '(none)'}."
        )

    def resolve_position(self, row_index: int, record_id_or_key: object, position: object) -> tuple[float, float]:
        record = self.resolve_record(row_index, record_id_or_key)
        try:
            numeric_position = float(position)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"Invalid comparison coordinate '{position}'.") from exc
        if numeric_position < 1 or numeric_position > record.record_length:
            raise ValidationError(
                f"Comparison coordinate {numeric_position:g} is outside record "
                f"'{record.ref.record_id}' length {record.record_length}."
            )
        return record.x + (numeric_position * record.scale), self.resolve_row(row_index).y_axis

    def resolve_interval(
        self,
        row_index: int,
        record_id_or_key: object,
        start: object,
        end: object,
    ) -> LinearInterval:
        record = self.resolve_record(row_index, record_id_or_key)
        try:
            start_value = float(start)
            end_value = float(end)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"Invalid comparison interval '{start}'..'{end}'.") from exc
        for value in (start_value, end_value):
            if value < 1 or value > record.record_length:
                raise ValidationError(
                    f"Comparison coordinate {value:g} is outside record "
                    f"'{record.ref.record_id}' length {record.record_length}."
                )
        return LinearInterval(
            record=record,
            start=start_value,
            end=end_value,
            start_x=record.x + (start_value * record.scale),
            end_x=record.x + (end_value * record.scale),
        )


def segment_canvas_config(
    canvas_config: LinearCanvasConfigurator,
    segment: LinearRecordLayout,
) -> LinearCanvasConfigurator:
    """Return a shallow canvas copy sized to one record segment."""

    segment_config = copy.copy(canvas_config)
    segment_config.alignment_width = float(segment.width)
    segment_config.longest_genome = int(segment.record_length)
    segment_config.normalize_length = True
    return segment_config


def segment_cfg(cfg: GbdrawConfig) -> GbdrawConfig:
    """Return a config copy that makes single-segment renderers use full segment width."""

    return replace(
        cfg,
        canvas=replace(
            cfg.canvas,
            linear=replace(cfg.canvas.linear, normalize_length=True),
        ),
    )


__all__ = [
    "LinearInputRow",
    "LinearInterval",
    "LinearLayoutIndex",
    "LinearRecordLayout",
    "LinearRecordRef",
    "LinearRowLayout",
    "adjacent_record_pair_to_linear_row_pair",
    "build_linear_layout_index",
    "clone_rows_for_rendering",
    "flatten_linear_rows",
    "linear_row_record_index_groups",
    "make_linear_input_row",
    "make_source_id",
    "record_index_to_linear_row_index",
    "segment_canvas_config",
    "segment_cfg",
]
