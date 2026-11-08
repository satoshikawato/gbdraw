#!/usr/bin/env python
# coding: utf-8

"""Focused TSV readers for row-coupled CLI inputs, and the writers of the records and comparisons tables."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Sequence

from ..core.numeric import parse_ascii_integer
from ..exceptions import ValidationError
from ..io.input_text import GZIP_ERRORS, open_input_text
from ..io.regions import parse_region_spec
from ..tracks import (
    normalize_listed_circular_track_slots,
    parse_circular_track_slots,
)
from ..tracks.circular import _LEGACY_SPACING_PARAM
from ..tracks.parsing import split_kv_list


InputKind = Literal["gbk", "gff_fasta", "fasta"]
# How a records table names each input kind in its errors.
_INPUT_KIND_ROWS: dict[InputKind, str] = {"gbk": "GenBank rows", "gff_fasta": "GFF3/FASTA rows", "fasta": "FASTA rows"}


@dataclass(frozen=True)
class TablePathDependency:
    row_index: int
    row_number: int
    column: str
    path: str


@dataclass(frozen=True)
class ConservationTableRow:
    row_index: int
    row_number: int
    blast: str
    label: str
    color: str
    # The comparison genome files of the ring, in the listed order (a ``;``
    # separated cell); empty without one.
    comparison_sequences: tuple[str, ...] = ()
    losat_gencode: int | None = None


@dataclass(frozen=True)
class ConservationTable:
    table_path: str
    rows: tuple[ConservationTableRow, ...]
    has_label_column: bool
    has_color_column: bool
    path_dependencies: tuple[TablePathDependency, ...]
    has_comparison_sequence_column: bool = False

    @property
    def conservation_blast_files(self) -> list[str]:
        return [row.blast for row in self.rows]

    @property
    def labels(self) -> list[str] | None:
        if not self.has_label_column:
            return None
        return [row.label for row in self.rows]

    @property
    def colors(self) -> list[str] | None:
        if not self.has_color_column:
            return None
        return [row.color for row in self.rows]

    @property
    def losat_gencodes(self) -> list[int] | None:
        """Per-ring TLOSATX tables when the column sets any (empty cells use 1)."""

        if not any(row.losat_gencode is not None for row in self.rows):
            return None
        return [row.losat_gencode or 1 for row in self.rows]

    @property
    def comparison_sequence_files(self) -> list[tuple[str, ...]] | None:
        if not self.has_comparison_sequence_column:
            return None
        return [row.comparison_sequences for row in self.rows]

    @property
    def dependency_paths(self) -> tuple[str, ...]:
        return tuple(dep.path for dep in self.path_dependencies)


@dataclass(frozen=True)
class ComparisonTableRow:
    row_index: int
    row_number: int
    # Empty for a ``source=losat`` row.
    blast: str
    query: str
    subject: str
    source: str = "table"


@dataclass(frozen=True)
class ComparisonTable:
    table_path: str
    rows: tuple[ComparisonTableRow, ...]
    path_dependencies: tuple[TablePathDependency, ...]

    @property
    def dependency_paths(self) -> tuple[str, ...]:
        return tuple(dep.path for dep in self.path_dependencies)


@dataclass(frozen=True)
class CircularTrackTable:
    table_path: str
    slot_specs: tuple[str, ...]
    axis_index: int
    path_dependencies: tuple[TablePathDependency, ...] = ()


@dataclass(frozen=True)
class RecordsTableRow:
    row_index: int
    row_number: int
    gbk: str
    gff: str
    fasta: str
    record_label: str
    record_subtitle: str
    record_id: str
    region: str
    reverse_complement: bool
    order: int | None
    row: int | None
    column: int | None
    topology: str | None = None
    display_start: int | None = None
    losat_gencode: int | None = None
    # The value of the ``record_group`` column; rows sharing a non-empty value are the windows of one record (D-32).
    record_group: str = ""


@dataclass(frozen=True)
class RecordsTableRecord:
    """The rows of one record a records table draws, in table order.

    A row without a ``record_group`` value is a record of its own; rows sharing a
    value are the windows of one record (each row's region, reverse_complement and
    column are its window's). The record's label, subtitle, topology and LOSAT
    genetic code are the non-empty value its rows give; the reader refuses rows of
    one record that give different ones.
    """

    rows: tuple[RecordsTableRow, ...]

    @property
    def first(self) -> RecordsTableRow:
        return self.rows[0]

    @property
    def is_windowed(self) -> bool:
        """Whether the record is drawn through windows: rows that share a ``record_group`` value."""
        return bool(self.first.record_group)

    @property
    def record_label(self) -> str:
        return next((row.record_label for row in self.rows if row.record_label), "")

    @property
    def record_subtitle(self) -> str:
        return next((row.record_subtitle for row in self.rows if row.record_subtitle), "")

    @property
    def topology(self) -> str | None:
        return next((row.topology for row in self.rows if row.topology is not None), None)

    @property
    def losat_gencode(self) -> int | None:
        return next((row.losat_gencode for row in self.rows if row.losat_gencode is not None), None)


# What the rows of one ``record_group`` must name alike, and what they may give once for the record (D-32).
_RECORD_GROUP_SOURCE_FIELDS = ("gbk", "gff", "fasta", "record_id", "row")
_RECORD_GROUP_RECORD_FIELDS = ("record_label", "record_subtitle", "topology", "losat_gencode")


def _group_records_table_rows(
    table_path: Path, rows: Sequence[RecordsTableRow]
) -> tuple[RecordsTableRecord, ...]:
    """The table's records in table order: rows sharing a ``record_group`` value are one record's windows."""

    groups: list[list[RecordsTableRow]] = []
    grouped: dict[str, list[RecordsTableRow]] = {}
    for row in rows:
        if not row.record_group:
            groups.append([row])
            continue
        if not row.region:
            raise ValidationError(
                _cell_error(
                    table_path, row.row_number, "region",
                    f"value is required in a row with {RECORD_GROUP_COLUMN} {row.record_group!r}",
                ),
                diagnostic={"code": "TABLE_INVALID", "field": "region", "row": row.row_number},
            )
        if row.record_group in grouped:
            grouped[row.record_group].append(row)
            continue
        grouped[row.record_group] = [row]
        groups.append(grouped[row.record_group])
    for group in grouped.values():
        first = group[0]
        for name in _RECORD_GROUP_SOURCE_FIELDS:
            values = {getattr(row, name) for row in group}
            if len(values) > 1:
                raise ValidationError(
                    f"{table_path}: rows {', '.join(str(row.row_number) for row in group)} share "
                    f"{RECORD_GROUP_COLUMN} {first.record_group!r} and must name the same {name}.",
                    diagnostic={"code": "TABLE_INVALID", "field": name, "row": group[1].row_number},
                )
        for name in _RECORD_GROUP_RECORD_FIELDS:
            values = {getattr(row, name) for row in group} - {"", None}
            if len(values) > 1:
                raise ValidationError(
                    f"{table_path}: rows {', '.join(str(row.row_number) for row in group)} share "
                    f"{RECORD_GROUP_COLUMN} {first.record_group!r}, one record, and give different {name} values.",
                    diagnostic={"code": "TABLE_INVALID", "field": name, "row": group[1].row_number},
                )
    return tuple(RecordsTableRecord(tuple(group)) for group in groups)


@dataclass(frozen=True)
class RecordsTable:
    table_path: str
    input_kind: InputKind
    rows: tuple[RecordsTableRow, ...]
    path_dependencies: tuple[TablePathDependency, ...]
    # The records the rows draw, in table order (one per row unless rows share a ``record_group``).
    records: tuple[RecordsTableRecord, ...] = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "records", _group_records_table_rows(Path(self.table_path), self.rows))

    @property
    def gbk_files(self) -> list[str]:
        return [row.gbk for row in self.rows]

    @property
    def gff_files(self) -> list[str]:
        return [row.gff for row in self.rows]

    @property
    def fasta_files(self) -> list[str]:
        return [row.fasta for row in self.rows]

    @property
    def record_labels(self) -> list[str]:
        return [row.record_label for row in self.rows]

    @property
    def record_subtitles(self) -> list[str]:
        return [row.record_subtitle for row in self.rows]

    @property
    def record_ids(self) -> list[str]:
        return [row.record_id for row in self.rows]

    @property
    def reverse_flags(self) -> list[bool]:
        return [row.reverse_complement for row in self.rows]

    @property
    def has_multi_record_placement(self) -> bool:
        return any(row.row is not None or row.column is not None for row in self.rows)

    @property
    def dependency_paths(self) -> tuple[str, ...]:
        return tuple(dep.path for dep in self.path_dependencies)

    def row_scoped_region_specs(self) -> list[str]:
        return [
            f"#{loaded_index}:{row.region}"
            for loaded_index, row in enumerate(self.rows, start=1)
            if row.region
        ]


@dataclass(frozen=True)
class _RawRow:
    row_index: int
    row_number: int
    values: dict[str, str]


@dataclass(frozen=True)
class _TrackRow:
    row_index: int
    row_number: int
    slot_id: str
    renderer: str
    side: str | None
    params: tuple[tuple[str, str], ...]
    values: dict[str, str]


_CONSERVATION_COLUMNS = frozenset({"blast", "label", "color", "comparison_sequence", "losat_gencode"})
# Retired conservation-table columns (design D18, OD-2) -> replacement.
_RETIRED_CONSERVATION_COLUMNS = {"comparison_fasta": "comparison_sequence"}
_COMPARISON_COLUMNS = frozenset({"blast", "query", "subject", "source"})
_COMPARISON_SOURCES = frozenset({"table", "losat"})
_CIRCULAR_TRACK_COLUMNS = frozenset(
    {
        "id",
        "renderer",
        "side",
        "r",
        "w",
        "inner_gap_px",
        "outer_gap_px",
        "z",
        "params",
    }
)
# The records-table column whose shared value makes rows the windows of one record (D-32).
RECORD_GROUP_COLUMN = "record_group"
_RECORDS_COLUMNS = frozenset(
    {
        RECORD_GROUP_COLUMN,
        "topology",
        "display_start",
        "gbk",
        "gff",
        "fasta",
        "record_label",
        "record_subtitle",
        "record_id",
        "region",
        "reverse_complement",
        "order",
        "row",
        "column",
        "losat_gencode",
    }
)
_BOOLEAN_TRUE = {"1", "true", "yes", "y", "on"}
_BOOLEAN_FALSE = {"", "0", "false", "no", "n", "off", "none", "null", "-"}
_TRACK_TABLE_SIDES = {"outside", "axis", "inside"}
_CIRCULAR_TRACK_STRUCTURAL_PARAM_KEYS = frozenset(
    {
        "id",
        "renderer",
        "type",
        "side",
        "r",
        "radius",
        "w",
        "width",
        "spacing",
        "inner_gap_px",
        "outer_gap_px",
        "z",
        "z_index",
        "zindex",
        "enabled",
        "show",
        "visible",
        "strict",
        "compress",
        "reserve",
    }
)
_FEATURE_LANE_PARAM_KEYS = frozenset({"lane_direction", "lanes"})


# Separates the files of one comparison genome in a ``comparison_sequence`` cell.
COMPARISON_SEQUENCE_SEPARATOR = ";"


def read_conservation_table(path: str, *, losat: bool = False) -> ConservationTable:
    """Read a ring manifest.

    A ``comparison_sequence`` cell lists the files of the ring's comparison
    genome, separated by ``;``. Precomputed rings need ``blast``. With ``losat`` (``--losat``) the search
    builds each ring from ``comparison_sequence``, so ``blast`` is rejected; the
    optional ``losat_gencode`` column sets each TLOSATX ring's table.
    """

    table_path, header, rows = _read_tsv_table(
        path,
        allowed_columns=_CONSERVATION_COLUMNS | frozenset(_RETIRED_CONSERVATION_COLUMNS),
        table_name="conservation table",
    )
    for retired, replacement in _RETIRED_CONSERVATION_COLUMNS.items():
        if retired in header:
            raise ValidationError(
                f"{table_path}: conservation table column {retired} was retired; "
                f"use {replacement} (FASTA, GenBank, or DDBJ).",
                diagnostic={"code": "TABLE_INVALID", "field": retired},
            )
    if losat:
        if "blast" in header:
            raise ValidationError(
                f"{table_path}: conservation table column blast cannot be combined with "
                "--losat; the search builds each ring from comparison_sequence.",
                diagnostic={
                    "code": "COMPARISON_INPUT",
                    "reason": "RING_LOSAT_INPUT",
                    "field": "blast",
                },
            )
        _require_columns(table_path, header, ("comparison_sequence",))
    else:
        _require_columns(table_path, header, ("blast",))
    if not rows:
        raise ValidationError(f"{table_path}: conservation table has no data rows.")

    parsed_rows: list[ConservationTableRow] = []
    dependencies: list[TablePathDependency] = []
    for row in rows:
        required = "comparison_sequence" if losat else "blast"
        if not row.values.get(required, "").strip():
            raise ValidationError(
                _cell_error(table_path, row.row_number, required, "value is required"),
                diagnostic={"code": "TABLE_INVALID", "reason": "REQUIRED", "field": required, "row": row.row_number},
            )
        blast_raw = row.values.get("blast", "").strip()
        blast_path = _resolve_table_path(table_path, blast_raw) if blast_raw else ""
        gencode_raw = row.values.get("losat_gencode", "").strip()
        losat_gencode: int | None = None
        if gencode_raw:
            losat_gencode = parse_ascii_integer(gencode_raw) or 0
            if losat_gencode <= 0:
                raise ValidationError(
                    _cell_error(
                        table_path, row.row_number, "losat_gencode",
                        "must be a positive integer",
                    ),
                    diagnostic={
                        "code": "TABLE_INVALID",
                        "reason": "POSITIVE_INTEGER",
                        "field": "losat_gencode",
                        "row": row.row_number,
                    },
                )
        comparison_sequence_raw = row.values.get("comparison_sequence", "").strip()
        comparison_sequence_parts = (
            [part.strip() for part in comparison_sequence_raw.split(COMPARISON_SEQUENCE_SEPARATOR)]
            if comparison_sequence_raw
            else []
        )
        if not all(comparison_sequence_parts):
            raise ValidationError(
                _cell_error(
                    table_path, row.row_number, "comparison_sequence",
                    f"lists an empty path; separate the files of one genome with "
                    f"{COMPARISON_SEQUENCE_SEPARATOR!r}",
                ),
                diagnostic={
                    "code": "TABLE_INVALID",
                    "reason": "EMPTY_PATH",
                    "field": "comparison_sequence",
                    "row": row.row_number,
                },
            )
        comparison_sequence_paths = tuple(
            _resolve_table_path(table_path, part) for part in comparison_sequence_parts
        )
        parsed_rows.append(
            ConservationTableRow(
                row_index=row.row_index,
                row_number=row.row_number,
                blast=blast_path,
                label=row.values.get("label", "").strip(),
                color=row.values.get("color", "").strip(),
                comparison_sequences=comparison_sequence_paths,
                losat_gencode=losat_gencode,
            )
        )
        if blast_path:
            dependencies.append(
                TablePathDependency(
                    row_index=row.row_index,
                    row_number=row.row_number,
                    column="blast",
                    path=blast_path,
                )
            )
        dependencies.extend(
            TablePathDependency(
                row_index=row.row_index,
                row_number=row.row_number,
                column="comparison_sequence",
                path=path,
            )
            for path in comparison_sequence_paths
        )

    return ConservationTable(
        table_path=str(table_path),
        rows=tuple(parsed_rows),
        has_label_column="label" in header,
        has_color_column="color" in header,
        has_comparison_sequence_column="comparison_sequence" in header,
        path_dependencies=tuple(dependencies),
    )


def read_circular_track_table(
    path: str,
    *,
    _allow_legacy_transport: bool = False,
) -> CircularTrackTable:
    table_path, header, rows = _read_tsv_table(
        path,
        allowed_columns=_CIRCULAR_TRACK_COLUMNS,
        table_name="circular track table",
    )
    _require_columns(table_path, header, ("id", "renderer"))
    if not rows:
        raise ValidationError(f"{table_path}: circular track table has no data rows.")

    track_rows: list[_TrackRow] = []
    for row in rows:
        slot_id = row.values.get("id", "").strip()
        renderer = row.values.get("renderer", "").strip()
        if not slot_id:
            raise ValidationError(
                _cell_error(table_path, row.row_number, "id", "value is required")
            )
        if not renderer:
            raise ValidationError(
                _cell_error(table_path, row.row_number, "renderer", "value is required")
            )
        side = _normalize_track_table_side(
            table_path,
            row.row_number,
            row.values.get("side", ""),
        )
        params = _parse_circular_track_table_params(
            table_path,
            row.row_number,
            renderer,
            row.values.get("params", ""),
            _allow_legacy_transport=_allow_legacy_transport,
        )
        track_rows.append(
            _TrackRow(
                row_index=row.row_index,
                row_number=row.row_number,
                slot_id=slot_id,
                renderer=renderer,
                side=side,
                params=params,
                values=row.values,
            )
        )

    explicit_axis_rows = [row for row in track_rows if row.side == "axis"]
    if len(explicit_axis_rows) > 1:
        row_numbers = ", ".join(str(row.row_number) for row in explicit_axis_rows[:2])
        raise ValidationError(
            f"{table_path}: circular track table allows only one side=axis row; found rows {row_numbers}."
        )
    axis_row_index: int | None = None
    if explicit_axis_rows:
        axis_row_index = explicit_axis_rows[0].row_index
    else:
        for track_row in track_rows:
            if track_row.side is None and track_row.renderer.strip().lower() == "features":
                axis_row_index = track_row.row_index
                break

    normalized_rows: list[_TrackRow] = []
    for track_row in track_rows:
        side = track_row.side
        if axis_row_index is not None and track_row.row_index == axis_row_index:
            side = "axis"
        elif side is None:
            side = "inside"
        normalized_rows.append(
            _TrackRow(
                row_index=track_row.row_index,
                row_number=track_row.row_number,
                slot_id=track_row.slot_id,
                renderer=track_row.renderer,
                side=side,
                params=track_row.params,
                values=track_row.values,
            )
        )

    outside_rows = [row for row in normalized_rows if row.side == "outside"]
    axis_rows = [row for row in normalized_rows if row.side == "axis"]
    inside_rows = [row for row in normalized_rows if row.side == "inside"]
    ordered_rows = outside_rows + axis_rows + inside_rows
    axis_index = len(outside_rows)

    specs: list[str] = []
    for track_row in ordered_rows:
        specs.append(_circular_track_row_to_spec(table_path, track_row))
    try:
        slots = parse_circular_track_slots(
            specs,
            _allow_legacy_transport=_allow_legacy_transport,
        )
        normalize_listed_circular_track_slots(slots, axis_index)
    except Exception as exc:
        raise ValidationError(f"{table_path}: invalid circular track table: {exc}") from exc
    return CircularTrackTable(
        table_path=str(table_path),
        slot_specs=tuple(specs),
        axis_index=axis_index,
    )


def read_records_table(path: str) -> RecordsTable:
    table_path, _header, rows = _read_tsv_table(
        path,
        allowed_columns=_RECORDS_COLUMNS,
        table_name="records table",
    )
    if not rows:
        raise ValidationError(f"{table_path}: records table has no data rows.")

    parsed_rows: list[RecordsTableRow] = []
    dependencies: list[TablePathDependency] = []
    row_kinds: dict[InputKind, None] = {}
    seen_grid_cells: set[tuple[int, int]] = set()
    for row in rows:
        values = row.values
        gbk_raw = values.get("gbk", "").strip()
        gff_raw = values.get("gff", "").strip()
        fasta_raw = values.get("fasta", "").strip()
        if gbk_raw and (gff_raw or fasta_raw):
            raise ValidationError(
                f"{table_path}: row {row.row_number} cannot mix gbk with gff/fasta."
            )
        if not (gbk_raw or gff_raw or fasta_raw):
            raise ValidationError(
                f"{table_path}: row {row.row_number} must provide gbk, fasta, or both gff and fasta."
            )
        if gff_raw and not fasta_raw:
            raise ValidationError(
                _cell_error(table_path, row.row_number, "fasta", "value is required for GFF3/FASTA rows")
            )
        row_kinds.setdefault("gbk" if gbk_raw else "gff_fasta" if gff_raw else "fasta")

        gbk_path = _resolve_table_path(table_path, gbk_raw) if gbk_raw else ""
        gff_path = _resolve_table_path(table_path, gff_raw) if gff_raw else ""
        fasta_path = _resolve_table_path(table_path, fasta_raw) if fasta_raw else ""
        for dependency_column, resolved_path in (
            ("gbk", gbk_path),
            ("gff", gff_path),
            ("fasta", fasta_path),
        ):
            if resolved_path:
                dependencies.append(
                    TablePathDependency(
                        row_index=row.row_index,
                        row_number=row.row_number,
                        column=dependency_column,
                        path=resolved_path,
                    )
                )

        order = _parse_optional_positive_int(table_path, row, "order")
        grid_row = _parse_optional_positive_int(table_path, row, "row")
        column = _parse_optional_positive_int(table_path, row, "column")
        if grid_row is not None and column is not None:
            key = (grid_row, column)
            if key in seen_grid_cells:
                raise ValidationError(
                    f"{table_path}: duplicate records table placement row={grid_row}, column={column}."
                )
            seen_grid_cells.add(key)
        region = values.get("region", "").strip()
        _validate_records_table_region(table_path, row.row_number, region, values.get("record_id", "").strip())
        topology = values.get("topology", "").strip().lower()
        if topology not in {"", "auto", "circular", "linear"}:
            raise ValidationError(_cell_error(
                table_path, row.row_number, "topology", "expected auto, circular, or linear"
            ))
        display_start = _parse_optional_positive_int(table_path, row, "display_start")
        if display_start is not None and (topology == "linear" or region):
            raise ValidationError(_cell_error(
                table_path, row.row_number, "display_start",
                "cannot be combined with linear topology or a region crop",
            ))
        parsed_rows.append(
            RecordsTableRow(
                row_index=row.row_index,
                row_number=row.row_number,
                gbk=gbk_path,
                gff=gff_path,
                fasta=fasta_path,
                record_label=values.get("record_label", "").strip(),
                record_subtitle=values.get("record_subtitle", "").strip(),
                record_id=values.get("record_id", "").strip(),
                region=region,
                reverse_complement=_parse_table_bool(table_path, row, "reverse_complement"),
                order=order,
                row=grid_row,
                column=column,
                topology=topology if topology in {"circular", "linear"} else None,
                display_start=display_start,
                losat_gencode=_parse_optional_positive_int(table_path, row, "losat_gencode"),
                record_group=values.get(RECORD_GROUP_COLUMN, "").strip(),
            )
        )

    if len(row_kinds) > 1:
        raise ValidationError(
            f"{table_path}: records table cannot mix "
            + " with ".join(rows for kind, rows in _INPUT_KIND_ROWS.items() if kind in row_kinds) + "."
        )
    _validate_records_table_placement(table_path, parsed_rows)
    ordered_rows = tuple(
        sorted(
            parsed_rows,
            key=lambda item: (
                item.order is None,
                item.order if item.order is not None else item.row_index,
                item.row_index,
            ),
        )
    )
    return RecordsTable(
        table_path=str(table_path),
        input_kind=next(iter(row_kinds)),
        rows=ordered_rows,
        path_dependencies=tuple(dependencies),
    )


def read_comparisons_table(path: str) -> ComparisonTable:
    """Read a Linear comparison manifest with explicit record selectors."""

    table_path, header, rows = _read_tsv_table(
        path,
        allowed_columns=_COMPARISON_COLUMNS,
        table_name="comparisons table",
    )
    required = {"blast", "query", "subject"}
    missing = required - set(header)
    if missing:
        raise ValidationError(
            f"{table_path}: comparisons table is missing required column(s): "
            f"{', '.join(sorted(missing))}."
        )
    if not rows:
        raise ValidationError(f"{table_path}: comparisons table has no data rows.")

    parsed: list[ComparisonTableRow] = []
    dependencies: list[TablePathDependency] = []
    for row in rows:
        values = row.values
        blast_raw = values.get("blast", "").strip()
        query = values.get("query", "").strip()
        subject = values.get("subject", "").strip()
        source = values.get("source", "").strip().lower() or "table"
        if source not in _COMPARISON_SOURCES:
            raise ValidationError(
                _cell_error(table_path, row.row_number, "source", "expected table or losat"),
                diagnostic={"code": "TABLE_INVALID", "field": "source", "row": row.row_number},
            )
        required_cells: tuple[tuple[str, str], ...] = (("query", query), ("subject", subject))
        if source == "table":
            required_cells = (("blast", blast_raw), *required_cells)
        elif blast_raw:
            raise ValidationError(
                _cell_error(
                    table_path, row.row_number, "blast",
                    "must be empty in a source=losat row; gbdraw runs the search",
                ),
                diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN", "row": row.row_number},
            )
        for column, value in required_cells:
            if not value:
                raise ValidationError(
                    _cell_error(table_path, row.row_number, column, "value is required")
                )
        blast = _resolve_table_path(table_path, blast_raw) if blast_raw else ""
        if blast and not Path(blast).is_file():
            raise ValidationError(
                _cell_error(table_path, row.row_number, "blast", f"file does not exist: {blast}")
            )
        if query == subject:
            raise ValidationError(
                _cell_error(
                    table_path,
                    row.row_number,
                    "subject",
                    "query and subject must identify different records",
                )
            )
        if blast:
            dependencies.append(
                TablePathDependency(
                    row_index=row.row_index,
                    row_number=row.row_number,
                    column="blast",
                    path=blast,
                )
            )
        parsed.append(
            ComparisonTableRow(
                row_index=row.row_index,
                row_number=row.row_number,
                blast=blast,
                query=query,
                subject=subject,
                source=source,
            )
        )
    return ComparisonTable(
        table_path=str(table_path),
        rows=tuple(parsed),
        path_dependencies=tuple(dependencies),
    )


def _read_tsv_table(
    path: str,
    *,
    allowed_columns: frozenset[str],
    table_name: str,
) -> tuple[Path, tuple[str, ...], list[_RawRow]]:
    table_path = Path(str(path))
    try:
        with open_input_text(table_path, encoding="utf-8-sig", newline="") as handle:
            rows = list(enumerate(csv.reader(handle, delimiter="\t"), start=1))
    except (OSError, *GZIP_ERRORS) as exc:  # Also a gzip table cut off or damaged.
        raise ValidationError(f"Could not read {table_name}: {table_path}") from exc

    header: tuple[str, ...] | None = None
    raw_rows: list[_RawRow] = []
    for line_number, cells in rows:
        if _cells_are_blank(cells):
            continue
        if header is None:
            header = tuple(cell.strip() for cell in cells)
            _validate_header(table_path, header, allowed_columns)
            continue
        if len(cells) > len(header) and any(cell.strip() for cell in cells[len(header):]):
            raise ValidationError(
                f"{table_path}: row {line_number} has more columns than the header."
            )
        values = {
            column: (cells[index].strip() if index < len(cells) else "")
            for index, column in enumerate(header)
        }
        if all(value == "" for value in values.values()):
            continue
        raw_rows.append(
            _RawRow(
                row_index=len(raw_rows),
                row_number=line_number,
                values=values,
            )
        )
    if header is None:
        raise ValidationError(f"{table_path}: {table_name} is empty or has no header row.")
    return table_path, header, raw_rows


def _validate_header(
    table_path: Path,
    header: Sequence[str],
    allowed_columns: frozenset[str],
) -> None:
    if not header or all(not column for column in header):
        raise ValidationError(f"{table_path}: table header is empty.")
    seen: set[str] = set()
    for column in header:
        if not column:
            raise ValidationError(f"{table_path}: table header contains an empty column name.")
        if column in seen:
            raise ValidationError(f"{table_path}: duplicate table column {column!r}.")
        seen.add(column)
        if column not in allowed_columns:
            raise ValidationError(f"{table_path}: unknown table column {column!r}.")


def _require_columns(table_path: Path, header: Sequence[str], required: Sequence[str]) -> None:
    missing = [column for column in required if column not in header]
    if missing:
        joined = ", ".join(missing)
        raise ValidationError(f"{table_path}: missing required column(s): {joined}.")


def _cells_are_blank(cells: Sequence[str]) -> bool:
    return not cells or all(str(cell).strip() == "" for cell in cells)


def _resolve_table_path(table_path: Path, raw_value: str) -> str:
    candidate = Path(str(raw_value))
    if candidate.is_absolute():
        return str(candidate)
    return str((table_path.parent / candidate).resolve(strict=False))


def _parse_optional_positive_int(
    table_path: Path,
    row: _RawRow,
    column: str,
) -> int | None:
    raw = row.values.get(column, "").strip()
    if not raw:
        return None
    parsed = parse_ascii_integer(raw)
    if parsed is None or parsed <= 0:
        raise ValidationError(
            _cell_error(table_path, row.row_number, column, f"expected a positive integer, got {raw!r}")
        )
    return parsed


def _parse_table_bool(table_path: Path, row: _RawRow, column: str) -> bool:
    raw = row.values.get(column, "").strip()
    normalized = raw.lower()
    if normalized in _BOOLEAN_TRUE:
        return True
    if normalized in _BOOLEAN_FALSE:
        return False
    raise ValidationError(
        _cell_error(
            table_path,
            row.row_number,
            column,
            f"expected a boolean value, got {raw!r}",
        )
    )


def _normalize_track_table_side(
    table_path: Path,
    row_number: int,
    raw_value: str,
) -> str | None:
    raw = str(raw_value or "").strip().lower()
    if not raw:
        return None
    if raw not in _TRACK_TABLE_SIDES:
        raise ValidationError(
            _cell_error(
                table_path,
                row_number,
                "side",
                "expected one of outside, axis, inside",
            )
        )
    return raw


def _circular_track_row_to_spec(table_path: Path, row: _TrackRow) -> str:
    if row.side == "axis":
        if row.renderer.strip().lower() != "features":
            raise ValidationError(
                _cell_error(
                    table_path,
                    row.row_number,
                    "side",
                    "side=axis is only supported for renderer=features",
                )
            )
    options: list[str] = []
    for column, option in (
        ("r", "r"),
        ("w", "w"),
        ("inner_gap_px", "inner_gap_px"),
        ("outer_gap_px", "outer_gap_px"),
        ("z", "z"),
    ):
        value = row.values.get(column, "").strip()
        if value:
            options.append(f"{option}={value}")
    if row.side == "axis":
        options.extend(["side=overlay", "lane_direction=split"])
    elif row.side:
        options.append(f"side={row.side}")
    options.extend(f"{key}={value}" for key, value in row.params)
    suffix = f"@{','.join(options)}" if options else ""
    return f"{row.slot_id}:{row.renderer}{suffix}"


def _parse_circular_track_table_params(
    table_path: Path,
    row_number: int,
    renderer: str,
    raw_params: str,
    *,
    _allow_legacy_transport: bool = False,
) -> tuple[tuple[str, str], ...]:
    params = str(raw_params or "").strip()
    if not params:
        return ()
    try:
        parsed_params = split_kv_list(params)
    except ValueError as exc:
        raise ValidationError(
            _cell_error(table_path, row_number, "params", str(exc))
        ) from exc
    reserved_keys = set(_CIRCULAR_TRACK_STRUCTURAL_PARAM_KEYS)
    if renderer.strip().lower() == "features":
        reserved_keys.update(_FEATURE_LANE_PARAM_KEYS)
    for key, _value in parsed_params:
        normalized_key = key.strip().lower()
        if normalized_key.startswith("_") and not (
            _allow_legacy_transport
            and normalized_key == _LEGACY_SPACING_PARAM
        ):
            raise ValidationError(
                _cell_error(
                    table_path,
                    row_number,
                    "params",
                    f"private key {key!r} is not allowed",
                )
            )
        if normalized_key in reserved_keys:
            raise ValidationError(
                _cell_error(
                    table_path,
                    row_number,
                    "params",
                    f"structural key {key!r} is not allowed; use the dedicated table columns",
                )
            )
    return tuple(parsed_params)


def _validate_records_table_region(
    table_path: Path,
    row_number: int,
    raw_region: str,
    record_id: str,
) -> None:
    if not raw_region:
        return
    try:
        parsed = parse_region_spec(raw_region)
    except ValueError as exc:
        record = f" of record {record_id}" if record_id else ""
        raise ValidationError(
            _cell_error(
                table_path,
                row_number,
                "region",
                f"invalid region value {raw_region!r}{record}: {exc}",
            )
        ) from exc
    if (
        parsed.record_id is not None
        or parsed.record_index is not None
        or parsed.file_selector is not None
    ):
        raise ValidationError(
            _cell_error(
                table_path,
                row_number,
                "region",
                "value must not include a record or file selector",
            )
        )


def _validate_records_table_placement(
    table_path: Path,
    rows: Sequence[RecordsTableRow],
) -> None:
    has_row = any(row.row is not None for row in rows)
    has_column = any(row.column is not None for row in rows)
    if has_row and any(row.row is None for row in rows):
        raise ValidationError(
            f"{table_path}: if any records table row value is present, every row must provide row."
        )
    if has_column and not has_row:
        raise ValidationError(f"{table_path}: column requires row.")
    if has_column and any(row.column is None for row in rows if row.row is not None):
        raise ValidationError(
            f"{table_path}: if any records table column value is present, every placed row must provide column."
        )


def _cell_error(table_path: Path | str, row_number: int, column: str, message: str) -> str:
    return f"{table_path}: row {row_number}, column {column!r}: {message}."


# Writers. A table written here reads back through its reader unchanged.


@dataclass(frozen=True)
class RecordsTableEntry:
    """One row of a records table to write (:func:`records_table_text`); the values are the cells."""

    gbk: str = ""
    gff: str = ""
    fasta: str = ""
    record_label: str = ""
    record_id: str = ""
    region: str = ""
    reverse_complement: bool = False
    row: int | None = None
    column: int | None = None
    record_group: str = ""


@dataclass(frozen=True)
class ComparisonsTableEntry:
    """One row of a comparisons table to write (:func:`comparisons_table_text`)."""

    query: str
    subject: str
    source: Literal["table", "losat"] = "losat"
    blast: str = ""


@dataclass(frozen=True)
class ConservationTableEntry:
    """One row of a conservation table to write (:func:`conservation_table_text`)."""

    blast: str
    comparison_sequences: tuple[str, ...] = ()
    label: str = ""
    color: str = ""


def _table_text(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buffer.getvalue()


def _cell(value: int | None) -> str:
    return "" if value is None else str(value)


def records_table_text(entries: Sequence[RecordsTableEntry]) -> str:
    """A records table of ``entries`` in order: the input columns they use, then label, record, region,
    orientation, placement and record group."""

    inputs = [name for name in ("gbk", "gff", "fasta") if any(getattr(entry, name) for entry in entries)]
    header = [*inputs, "record_label", "record_id", "region", "reverse_complement", "row", "column", RECORD_GROUP_COLUMN]
    return _table_text(
        header,
        [
            [
                *(getattr(entry, name) for name in inputs),
                entry.record_label,
                entry.record_id,
                entry.region,
                "1" if entry.reverse_complement else "0",
                _cell(entry.row),
                _cell(entry.column),
                entry.record_group,
            ]
            for entry in entries
        ],
    )


def conservation_table_text(entries: Sequence[ConservationTableEntry]) -> str:
    """A conservation table of ``entries`` in order; a ring's genome files share one ``;`` cell."""

    for entry in entries:
        for path in entry.comparison_sequences:
            if COMPARISON_SEQUENCE_SEPARATOR in path:
                raise ValidationError(
                    f"The conservation table cannot name {path}: a comparison_sequence cell separates "
                    f"files with {COMPARISON_SEQUENCE_SEPARATOR!r}. Rename the file.",
                    diagnostic={"code": "INPUT_INVALID", "reason": "TABLE_PATH", "field": "comparison_sequence"},
                )
    return _table_text(
        ("blast", "comparison_sequence", "label", "color"),
        [
            [entry.blast, COMPARISON_SEQUENCE_SEPARATOR.join(entry.comparison_sequences), entry.label, entry.color]
            for entry in entries
        ],
    )


def comparisons_table_text(entries: Sequence[ComparisonsTableEntry]) -> str:
    """A comparisons table of ``entries`` in order."""

    return _table_text(
        ("blast", "query", "subject", "source"),
        [[entry.blast, entry.query, entry.subject, entry.source] for entry in entries],
    )


__all__ = [
    "CircularTrackTable",
    "COMPARISON_SEQUENCE_SEPARATOR",
    "ComparisonsTableEntry",
    "ConservationTable",
    "ConservationTableEntry",
    "ConservationTableRow",
    "ComparisonTable",
    "ComparisonTableRow",
    "RECORD_GROUP_COLUMN",
    "RecordsTable",
    "RecordsTableEntry",
    "RecordsTableRecord",
    "RecordsTableRow",
    "TablePathDependency",
    "comparisons_table_text",
    "conservation_table_text",
    "read_circular_track_table",
    "read_conservation_table",
    "read_comparisons_table",
    "read_records_table",
    "records_table_text",
]
