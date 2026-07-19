#!/usr/bin/env python
# coding: utf-8

"""Focused TSV readers for row-coupled CLI inputs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

from ..exceptions import ValidationError
from ..io.regions import parse_region_spec
from ..tracks import (
    normalize_circular_track_slots_with_axis,
    parse_circular_track_slots,
)
from ..tracks.parsing import split_kv_list


InputKind = Literal["gbk", "gff_fasta"]


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
    comparison_fasta: str = ""


@dataclass(frozen=True)
class ConservationTable:
    table_path: str
    rows: tuple[ConservationTableRow, ...]
    has_label_column: bool
    has_color_column: bool
    path_dependencies: tuple[TablePathDependency, ...]
    has_comparison_fasta_column: bool = False

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
    def comparison_fasta_files(self) -> list[str | None] | None:
        if not self.has_comparison_fasta_column:
            return None
        return [row.comparison_fasta or None for row in self.rows]

    @property
    def dependency_paths(self) -> tuple[str, ...]:
        return tuple(dep.path for dep in self.path_dependencies)


@dataclass(frozen=True)
class ComparisonTableRow:
    row_index: int
    row_number: int
    blast: str
    query: str
    subject: str


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


@dataclass(frozen=True)
class RecordsTable:
    table_path: str
    input_kind: InputKind
    rows: tuple[RecordsTableRow, ...]
    path_dependencies: tuple[TablePathDependency, ...]

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

    def multi_record_positions(self) -> list[str]:
        if not self.has_multi_record_placement:
            return []
        positioned = [
            (loaded_index, row)
            for loaded_index, row in enumerate(self.rows, start=1)
            if row.row is not None
        ]
        positioned.sort(
            key=lambda item: (
                int(item[1].row or 0),
                int(item[1].column) if item[1].column is not None else item[0],
                item[0],
            )
        )
        return [f"#{loaded_index}@{row.row}" for loaded_index, row in positioned]


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


_CONSERVATION_COLUMNS = frozenset({"blast", "label", "color", "comparison_fasta"})
_COMPARISON_COLUMNS = frozenset({"blast", "query", "subject"})
_CIRCULAR_TRACK_COLUMNS = frozenset(
    {
        "id",
        "renderer",
        "side",
        "r",
        "w",
        "spacing",
        "inner_gap_px",
        "outer_gap_px",
        "z",
        "params",
    }
)
_RECORDS_COLUMNS = frozenset(
    {
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


def read_conservation_table(path: str) -> ConservationTable:
    table_path, header, rows = _read_tsv_table(
        path,
        allowed_columns=_CONSERVATION_COLUMNS,
        table_name="conservation table",
    )
    _require_columns(table_path, header, ("blast",))
    if not rows:
        raise ValidationError(f"{table_path}: conservation table has no data rows.")

    parsed_rows: list[ConservationTableRow] = []
    dependencies: list[TablePathDependency] = []
    for row in rows:
        blast_raw = row.values.get("blast", "").strip()
        if not blast_raw:
            raise ValidationError(
                _cell_error(table_path, row.row_number, "blast", "value is required")
            )
        blast_path = _resolve_table_path(table_path, blast_raw)
        comparison_fasta_raw = row.values.get("comparison_fasta", "").strip()
        comparison_fasta_path = (
            _resolve_table_path(table_path, comparison_fasta_raw)
            if comparison_fasta_raw
            else ""
        )
        parsed_rows.append(
            ConservationTableRow(
                row_index=row.row_index,
                row_number=row.row_number,
                blast=blast_path,
                label=row.values.get("label", "").strip(),
                color=row.values.get("color", "").strip(),
                comparison_fasta=comparison_fasta_path,
            )
        )
        dependencies.append(
            TablePathDependency(
                row_index=row.row_index,
                row_number=row.row_number,
                column="blast",
                path=blast_path,
            )
        )
        if comparison_fasta_path:
            dependencies.append(
                TablePathDependency(
                    row_index=row.row_index,
                    row_number=row.row_number,
                    column="comparison_fasta",
                    path=comparison_fasta_path,
                )
            )

    return ConservationTable(
        table_path=str(table_path),
        rows=tuple(parsed_rows),
        has_label_column="label" in header,
        has_color_column="color" in header,
        has_comparison_fasta_column="comparison_fasta" in header,
        path_dependencies=tuple(dependencies),
    )


def read_circular_track_table(path: str) -> CircularTrackTable:
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
        for row in track_rows:
            if row.side is None and row.renderer.strip().lower() == "features":
                axis_row_index = row.row_index
                break

    normalized_rows: list[_TrackRow] = []
    for row in track_rows:
        side = row.side
        if axis_row_index is not None and row.row_index == axis_row_index:
            side = "axis"
        elif side is None:
            side = "inside"
        normalized_rows.append(
            _TrackRow(
                row_index=row.row_index,
                row_number=row.row_number,
                slot_id=row.slot_id,
                renderer=row.renderer,
                side=side,
                params=row.params,
                values=row.values,
            )
        )

    outside_rows = [row for row in normalized_rows if row.side == "outside"]
    axis_rows = [row for row in normalized_rows if row.side == "axis"]
    inside_rows = [row for row in normalized_rows if row.side == "inside"]
    ordered_rows = outside_rows + axis_rows + inside_rows
    axis_index = len(outside_rows)

    specs: list[str] = []
    for row in ordered_rows:
        specs.append(_circular_track_row_to_spec(table_path, row))
    try:
        slots = parse_circular_track_slots(specs)
        normalize_circular_track_slots_with_axis(slots, axis_index)
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
    has_gbk_rows = False
    has_gff_rows = False
    seen_grid_cells: set[tuple[int, int]] = set()
    for row in rows:
        values = row.values
        gbk_raw = values.get("gbk", "").strip()
        gff_raw = values.get("gff", "").strip()
        fasta_raw = values.get("fasta", "").strip()
        row_has_gbk = bool(gbk_raw)
        row_has_gff = bool(gff_raw or fasta_raw)
        if row_has_gbk and row_has_gff:
            raise ValidationError(
                f"{table_path}: row {row.row_number} cannot mix gbk with gff/fasta."
            )
        if not row_has_gbk and not row_has_gff:
            raise ValidationError(
                f"{table_path}: row {row.row_number} must provide either gbk or both gff and fasta."
            )
        if row_has_gff and not (gff_raw and fasta_raw):
            missing = "gff" if not gff_raw else "fasta"
            raise ValidationError(
                _cell_error(table_path, row.row_number, missing, "value is required for GFF3/FASTA rows")
            )

        gbk_path = _resolve_table_path(table_path, gbk_raw) if gbk_raw else ""
        gff_path = _resolve_table_path(table_path, gff_raw) if gff_raw else ""
        fasta_path = _resolve_table_path(table_path, fasta_raw) if fasta_raw else ""
        has_gbk_rows = has_gbk_rows or bool(gbk_path)
        has_gff_rows = has_gff_rows or bool(gff_path)
        for column, resolved_path in (
            ("gbk", gbk_path),
            ("gff", gff_path),
            ("fasta", fasta_path),
        ):
            if resolved_path:
                dependencies.append(
                    TablePathDependency(
                        row_index=row.row_index,
                        row_number=row.row_number,
                        column=column,
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
        _validate_records_table_region(table_path, row.row_number, region)
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
            )
        )

    if has_gbk_rows and has_gff_rows:
        raise ValidationError(
            f"{table_path}: records table cannot mix GenBank rows with GFF3/FASTA rows."
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
        input_kind="gbk" if has_gbk_rows else "gff_fasta",
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
        for column, value in (("blast", blast_raw), ("query", query), ("subject", subject)):
            if not value:
                raise ValidationError(
                    _cell_error(table_path, row.row_number, column, "value is required")
                )
        blast = _resolve_table_path(table_path, blast_raw)
        if not Path(blast).is_file():
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
        handle = table_path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise ValidationError(f"Could not read {table_name}: {table_path}") from exc

    with handle:
        reader = csv.reader(handle, delimiter="\t")
        header: tuple[str, ...] | None = None
        raw_rows: list[_RawRow] = []
        for line_number, cells in enumerate(reader, start=1):
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
    try:
        parsed = int(raw)
    except ValueError as exc:
        raise ValidationError(
            _cell_error(table_path, row.row_number, column, f"expected a positive integer, got {raw!r}")
        ) from exc
    if parsed <= 0:
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
        ("spacing", "spacing"),
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
) -> None:
    if not raw_region:
        return
    try:
        parsed = parse_region_spec(raw_region)
    except ValueError as exc:
        raise ValidationError(
            _cell_error(
                table_path,
                row_number,
                "region",
                f"invalid region value {raw_region!r}: {exc}",
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


__all__ = [
    "CircularTrackTable",
    "ConservationTable",
    "ConservationTableRow",
    "ComparisonTable",
    "ComparisonTableRow",
    "RecordsTable",
    "RecordsTableRow",
    "TablePathDependency",
    "read_circular_track_table",
    "read_conservation_table",
    "read_comparisons_table",
    "read_records_table",
]
