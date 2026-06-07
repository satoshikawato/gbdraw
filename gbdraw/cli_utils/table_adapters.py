"""CLI-boundary adapters for table-based arguments."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any, Literal

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.analysis.depth import read_depth_tsv  # type: ignore[reportMissingImports]
from gbdraw.analysis.depth_tracks import DepthTrackSpec  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import reverse_records
from gbdraw.io.regions import apply_region_specs, parse_region_spec
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot, ScalarSpec
from gbdraw.tracks.circular import SUPPORTED_CIRCULAR_TRACK_RENDERERS
from gbdraw.tracks.linear import SUPPORTED_LINEAR_TRACK_RENDERERS

from .tables import (
    DisplayRecordContext,
    HeaderedTableRow,
    SourceRecordContext,
    parse_optional_bool_cell,
    parse_optional_float_cell,
    parse_optional_int_cell,
    read_headered_tsv_table,
    resolve_table_path,
    table_error,
)


Mode = Literal["linear", "circular"]

_INPUT_TABLE_REQUIRED = ("input_id", "input_type")
_INPUT_TABLE_OPTIONAL = (
    "gbk",
    "gff",
    "fasta",
    "record_id",
    "region",
    "reverse_complement",
    "label",
    "order",
    "expand_records",
)
_DEPTH_TABLE_REQUIRED = ("record_id", "track_id", "file")
_DEPTH_TABLE_OPTIONAL = (
    "track_label",
    "track_color",
    "track_height",
    "track_width",
    "track_large_tick_interval",
    "track_small_tick_interval",
    "track_tick_font_size",
    "order",
)
_TRACK_TABLE_REQUIRED = ("slot_id", "renderer")
_TRACK_TABLE_OPTIONAL = (
    "order",
    "side",
    "track_id",
    "track_index",
    "height",
    "radius",
    "width",
    "spacing",
    "z",
    "enabled",
    "nt",
    "dinucleotide",
    "source_index",
    "lane_direction",
    "tick_label_layout",
)
_LINEAR_RENDERER_ALIASES = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
}
_CIRCULAR_RENDERER_ALIASES = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
    "conservation": "sequence_conservation",
    "circular_comparison": "sequence_conservation",
}


LoadGbkRecords = Callable[..., list[SeqRecord]]
LoadGffRecords = Callable[..., list[SeqRecord]]


@dataclass(frozen=True)
class InputTableRow:
    row_number: int
    input_id: str
    input_type: str
    gbk: str
    gff: str
    fasta: str
    record_id: str
    region: str
    reverse_complement: bool | None
    label: str
    order: int | None
    expand_records: bool


@dataclass(frozen=True)
class DepthTrackMetadata:
    label: str | None = None
    fill_color: str | None = None
    height: float | None = None
    width: float | None = None
    large_tick_interval: float | None = None
    small_tick_interval: float | None = None
    tick_font_size: float | None = None


@dataclass(frozen=True)
class DepthTrackTableResult:
    track_ids: tuple[str, ...]
    metadata_by_track_id: Mapping[str, DepthTrackMetadata]
    record_depth_tracks: list[list[DepthTrackSpec]]


@dataclass(frozen=True)
class TrackTableResult:
    slots: list[LinearTrackSlot] | list[CircularTrackSlot]
    axis_index: int | None = None


def _require_nonempty(row: HeaderedTableRow, column: str, table_name: str) -> str:
    value = row.cell(column)
    if not value:
        raise table_error(table_name, "value cannot be empty", row_number=row.row_number, column=column)
    return value


def _read_input_table(path: str) -> list[InputTableRow]:
    rows = read_headered_tsv_table(
        path,
        required=_INPUT_TABLE_REQUIRED,
        optional=_INPUT_TABLE_OPTIONAL,
        table_name="input_table",
    )
    parsed: list[InputTableRow] = []
    seen_exact_ids: set[str] = set()
    for row in rows:
        input_id = _require_nonempty(row, "input_id", "input_table")
        input_type = _require_nonempty(row, "input_type", "input_table").lower()
        if input_type not in {"gbk", "gff"}:
            raise table_error(
                "input_table",
                "input_type must be one of gbk or gff",
                row_number=row.row_number,
                column="input_type",
            )
        gbk = resolve_table_path(path, row.cell("gbk"))
        gff = resolve_table_path(path, row.cell("gff"))
        fasta = resolve_table_path(path, row.cell("fasta"))
        if input_type == "gbk":
            if not gbk:
                raise table_error("input_table", "gbk is required when input_type=gbk", row_number=row.row_number, column="gbk")
            if gff or fasta:
                raise table_error("input_table", "gff/fasta must be empty when input_type=gbk", row_number=row.row_number)
        if input_type == "gff":
            if not gff:
                raise table_error("input_table", "gff is required when input_type=gff", row_number=row.row_number, column="gff")
            if not fasta:
                raise table_error("input_table", "fasta is required when input_type=gff", row_number=row.row_number, column="fasta")
            if gbk:
                raise table_error("input_table", "gbk must be empty when input_type=gff", row_number=row.row_number, column="gbk")
        reverse_complement = parse_optional_bool_cell(
            row.cell("reverse_complement"),
            "reverse_complement",
            table_name="input_table",
            row_number=row.row_number,
            column="reverse_complement",
        )
        expand_records = bool(
            parse_optional_bool_cell(
                row.cell("expand_records"),
                "expand_records",
                table_name="input_table",
                row_number=row.row_number,
                column="expand_records",
            )
            or False
        )
        order = parse_optional_int_cell(
            row.cell("order"),
            "order",
            table_name="input_table",
            row_number=row.row_number,
            column="order",
        )
        if not expand_records:
            if input_id in seen_exact_ids:
                raise table_error(
                    "input_table",
                    f"duplicate input_id '{input_id}'",
                    row_number=row.row_number,
                    column="input_id",
                )
            seen_exact_ids.add(input_id)
        parsed.append(
            InputTableRow(
                row_number=row.row_number,
                input_id=input_id,
                input_type=input_type,
                gbk=gbk,
                gff=gff,
                fasta=fasta,
                record_id=row.cell("record_id"),
                region=row.cell("region"),
                reverse_complement=reverse_complement,
                label=row.cell("label"),
                order=order,
                expand_records=expand_records,
            )
        )
    return parsed


def _attach_input_id(record: SeqRecord, input_id: str) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations["gbdraw_input_id"] = input_id


def _attach_record_label(record: SeqRecord, label: str) -> None:
    text = str(label or "").strip()
    if not text:
        return
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations["gbdraw_record_label"] = text


def _derive_expanded_input_ids(prefix: str, records: Sequence[SeqRecord]) -> list[str]:
    record_ids = [str(getattr(record, "id", "") or "").strip() for record in records]
    if all(record_ids) and len(set(record_ids)) == len(record_ids):
        return [f"{prefix}:{record_id}" for record_id in record_ids]
    return [f"{prefix}:#{index + 1}" for index in range(len(records))]


def _load_source_records(
    row: InputTableRow,
    *,
    mode: Mode,
    load_gbk_records: LoadGbkRecords,
    load_gff_records: LoadGffRecords,
    selected_features_set: Sequence[str] | None,
    keep_all_features: bool,
) -> list[SeqRecord]:
    if row.input_type == "gbk":
        return load_gbk_records([row.gbk], mode, False)
    return load_gff_records(
        [row.gff],
        [row.fasta],
        mode,
        selected_features_set,
        keep_all_features=keep_all_features,
        load_comparison=False,
    )


def _apply_linear_row_region(record: SeqRecord, row: InputTableRow) -> SeqRecord:
    if not row.region:
        return record
    try:
        spec = parse_region_spec(row.region)
    except ValueError as exc:
        raise table_error(
            "input_table",
            str(exc),
            row_number=row.row_number,
            column="region",
        ) from exc
    if spec.file_selector is not None or spec.record_id is not None or spec.record_index is not None:
        raise table_error(
            "input_table",
            "region must be row-local, for example 1000-50000 or 1000-50000:rc",
            row_number=row.row_number,
            column="region",
        )
    if row.reverse_complement is True and spec.reverse_complement:
        raise table_error(
            "input_table",
            "reverse_complement=true cannot be combined with a region ending in :rc in the first implementation",
            row_number=row.row_number,
            column="region",
        )
    try:
        return apply_region_specs([record], [spec])[0]
    except ValueError as exc:
        raise table_error(
            "input_table",
            str(exc),
            row_number=row.row_number,
            column="region",
        ) from exc


def load_input_table_records(
    path: str,
    *,
    mode: Mode,
    load_gbk_records: LoadGbkRecords,
    load_gff_records: LoadGffRecords,
    selected_features_set: Sequence[str] | None = None,
    keep_all_features: bool = False,
) -> list[SeqRecord]:
    """Load displayed records from an --input_table file."""

    rows = sorted(_read_input_table(path), key=lambda item: (item.order if item.order is not None else item.row_number, item.row_number))
    out: list[SeqRecord] = []
    displayed_ids: set[str] = set()

    for row in rows:
        if mode == "circular":
            if row.region:
                raise table_error(
                    "input_table",
                    "region is not supported in circular --input_table in the first implementation",
                    row_number=row.row_number,
                    column="region",
                )
            if row.reverse_complement is not None:
                raise table_error(
                    "input_table",
                    "reverse_complement is not supported in circular --input_table in the first implementation",
                    row_number=row.row_number,
                    column="reverse_complement",
                )
        if row.expand_records and (row.record_id or row.region or row.label):
            raise table_error(
                "input_table",
                "expand_records=true cannot be combined with record_id, region, or label in the first implementation",
                row_number=row.row_number,
            )

        source_records = _load_source_records(
            row,
            mode=mode,
            load_gbk_records=load_gbk_records,
            load_gff_records=load_gff_records,
            selected_features_set=selected_features_set,
            keep_all_features=keep_all_features,
        )
        if not source_records:
            raise table_error("input_table", "source row loaded no records", row_number=row.row_number)

        if row.expand_records:
            records = reverse_records(source_records, bool(row.reverse_complement)) if mode == "linear" else list(source_records)
            derived_ids = _derive_expanded_input_ids(row.input_id, records)
            for record, derived_id in zip(records, derived_ids):
                if derived_id in displayed_ids:
                    raise table_error(
                        "input_table",
                        f"derived displayed input_id '{derived_id}' is duplicated",
                        row_number=row.row_number,
                        column="input_id",
                    )
                _attach_input_id(record, derived_id)
                displayed_ids.add(derived_id)
                out.append(record)
            continue

        if row.record_id:
            selected_index = SourceRecordContext(source_records).resolve_source_record_selector(
                row.record_id,
                row_number=row.row_number,
                column="record_id",
            )
            selected_record = source_records[selected_index]
        else:
            if len(source_records) != 1:
                raise table_error(
                    "input_table",
                    "row loaded multiple records; add record_id or set expand_records=true",
                    row_number=row.row_number,
                    column="record_id",
                )
            selected_record = source_records[0]

        if mode == "linear":
            selected_record = _apply_linear_row_region(selected_record, row)
            if row.reverse_complement is True:
                selected_record = reverse_records([selected_record], True)[0]
        _attach_input_id(selected_record, row.input_id)
        _attach_record_label(selected_record, row.label)
        if row.input_id in displayed_ids:
            raise table_error(
                "input_table",
                f"displayed input_id '{row.input_id}' is duplicated",
                row_number=row.row_number,
                column="input_id",
            )
        displayed_ids.add(row.input_id)
        out.append(selected_record)

    if not out:
        raise ValidationError("input_table: no records were loaded.")
    return out


def _default_depth_track_label(track_id: str, track_count: int, track_index: int) -> str:
    if track_count == 1:
        return "Depth"
    if track_id:
        return str(track_id)
    return f"Depth {track_index + 1}"


def _metadata_value(
    row: HeaderedTableRow,
    column: str,
    *,
    mode: Mode,
) -> str | float | int | None:
    raw = row.cell(column)
    if not raw:
        return None
    if column == "track_height":
        if mode == "circular":
            raise table_error(
                "depth_track_table",
                "track_height is only valid in linear mode",
                row_number=row.row_number,
                column=column,
            )
        return parse_optional_float_cell(raw, column, positive=True, table_name="depth_track_table", row_number=row.row_number, column=column)
    if column == "track_width":
        if mode == "linear":
            raise table_error(
                "depth_track_table",
                "track_width is only valid in circular mode",
                row_number=row.row_number,
                column=column,
            )
        return parse_optional_float_cell(raw, column, positive=True, table_name="depth_track_table", row_number=row.row_number, column=column)
    if column in {"track_large_tick_interval", "track_small_tick_interval", "track_tick_font_size"}:
        return parse_optional_float_cell(raw, column, positive=True, table_name="depth_track_table", row_number=row.row_number, column=column)
    return raw


def _set_track_metadata(
    metadata: dict[str, dict[str, str | float | int | None]],
    *,
    track_id: str,
    field: str,
    value: str | float | int | None,
    row: HeaderedTableRow,
) -> None:
    if value is None:
        return
    per_track = metadata.setdefault(track_id, {})
    previous = per_track.get(field)
    if previous is not None and previous != value:
        raise table_error(
            "depth_track_table",
            (
                f"conflicting {field} metadata for track_id '{track_id}'. "
                "Depth table metadata is track-scoped; use a distinct track_id for per-record styling."
            ),
            row_number=row.row_number,
            column=field,
        )
    per_track[field] = value


def load_depth_track_table(
    path: str,
    *,
    mode: Mode,
    records: Sequence[SeqRecord],
) -> DepthTrackTableResult:
    """Load --depth_track_table into named record-major DepthTrackSpec rows."""

    rows = read_headered_tsv_table(
        path,
        required=_DEPTH_TABLE_REQUIRED,
        optional=_DEPTH_TABLE_OPTIONAL,
        table_name="depth_track_table",
    )
    if not rows:
        raise ValidationError("depth_track_table: table contains no data rows.")

    display_context = DisplayRecordContext(records, table_name="depth_track_table")
    first_seen: dict[str, int] = {}
    orders: dict[str, int] = {}
    metadata: dict[str, dict[str, str | float | int | None]] = {}
    wildcard_files: dict[str, str] = {}
    specific_files: dict[tuple[int, str], str] = {}

    for appearance, row in enumerate(rows):
        record_selector = _require_nonempty(row, "record_id", "depth_track_table")
        track_id = _require_nonempty(row, "track_id", "depth_track_table")
        file_path = resolve_table_path(path, _require_nonempty(row, "file", "depth_track_table"))
        first_seen.setdefault(track_id, appearance)
        order = parse_optional_int_cell(
            row.cell("order"),
            "order",
            table_name="depth_track_table",
            row_number=row.row_number,
            column="order",
        )
        if order is not None:
            if track_id in orders and orders[track_id] != order:
                raise table_error(
                    "depth_track_table",
                    f"conflicting order metadata for track_id '{track_id}'",
                    row_number=row.row_number,
                    column="order",
                )
            orders[track_id] = order
        for column in (
            "track_label",
            "track_color",
            "track_height",
            "track_width",
            "track_large_tick_interval",
            "track_small_tick_interval",
            "track_tick_font_size",
        ):
            _set_track_metadata(
                metadata,
                track_id=track_id,
                field=column,
                value=_metadata_value(row, column, mode=mode),
                row=row,
            )

        if record_selector == "*":
            if track_id in wildcard_files:
                raise table_error(
                    "depth_track_table",
                    f"duplicate wildcard row for track_id '{track_id}'",
                    row_number=row.row_number,
                    column="record_id",
                )
            wildcard_files[track_id] = file_path
            continue

        record_index = display_context.resolve_display_record_selector(
            record_selector,
            row_number=row.row_number,
            column="record_id",
        )
        key = (record_index, track_id)
        if key in specific_files:
            raise table_error(
                "depth_track_table",
                f"duplicate row for record selector '{record_selector}' and track_id '{track_id}'",
                row_number=row.row_number,
            )
        specific_files[key] = file_path

    track_ids = tuple(
        sorted(
            first_seen,
            key=lambda track_id: (orders.get(track_id, first_seen[track_id]), first_seen[track_id]),
        )
    )
    metadata_by_track_id: dict[str, DepthTrackMetadata] = {}
    for track_id in track_ids:
        values = metadata.get(track_id, {})
        metadata_by_track_id[track_id] = DepthTrackMetadata(
            label=values.get("track_label") if isinstance(values.get("track_label"), str) else None,
            fill_color=values.get("track_color") if isinstance(values.get("track_color"), str) else None,
            height=float(values["track_height"]) if values.get("track_height") is not None else None,
            width=float(values["track_width"]) if values.get("track_width") is not None else None,
            large_tick_interval=float(values["track_large_tick_interval"]) if values.get("track_large_tick_interval") is not None else None,
            small_tick_interval=float(values["track_small_tick_interval"]) if values.get("track_small_tick_interval") is not None else None,
            tick_font_size=float(values["track_tick_font_size"]) if values.get("track_tick_font_size") is not None else None,
        )

    record_depth_tracks: list[list[DepthTrackSpec]] = []
    for record_index, _record in enumerate(records):
        row_specs: list[DepthTrackSpec] = []
        for track_index, track_id in enumerate(track_ids):
            file_path = specific_files.get((record_index, track_id), wildcard_files.get(track_id))
            if not file_path:
                continue
            track_meta = metadata_by_track_id[track_id]
            row_specs.append(
                DepthTrackSpec(
                    id=track_id,
                    label=track_meta.label or _default_depth_track_label(track_id, len(track_ids), track_index),
                    table=read_depth_tsv(file_path),
                    fill_color=track_meta.fill_color,
                    height=track_meta.height,
                    large_tick_interval=track_meta.large_tick_interval,
                    small_tick_interval=track_meta.small_tick_interval,
                    tick_font_size=track_meta.tick_font_size,
                )
            )
        record_depth_tracks.append(row_specs)

    return DepthTrackTableResult(
        track_ids=track_ids,
        metadata_by_track_id=metadata_by_track_id,
        record_depth_tracks=record_depth_tracks,
    )


def _normalize_renderer(renderer: str, *, mode: Mode) -> str:
    raw = str(renderer or "").strip().lower()
    if mode == "linear":
        return _LINEAR_RENDERER_ALIASES.get(raw, raw)
    return _CIRCULAR_RENDERER_ALIASES.get(raw, raw)


def _parse_linear_px_scalar(
    value: str,
    *,
    field_name: str,
    table_name: str,
    row_number: int,
    allow_zero: bool,
) -> ScalarSpec | None:
    text = str(value or "").strip()
    if not text:
        return None
    if text.endswith("%"):
        raise table_error(
            table_name,
            f"{field_name} only accepts px or unitless px values",
            row_number=row_number,
            column=field_name,
        )
    raw_value = text[:-2] if text.endswith("px") else text
    parsed = parse_optional_float_cell(
        raw_value,
        field_name,
        positive=not allow_zero,
        nonnegative=allow_zero,
        table_name=table_name,
        row_number=row_number,
        column=field_name,
    )
    return ScalarSpec(float(parsed), "px") if parsed is not None else None


def _parse_scalar_spec(
    value: str,
    *,
    field_name: str,
    table_name: str,
    row_number: int,
) -> ScalarSpec | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return ScalarSpec.parse(text)
    except Exception as exc:
        raise table_error(
            table_name,
            f"{field_name} must use the existing scalar syntax (px, %, or unitless factor)",
            row_number=row_number,
            column=field_name,
        ) from exc


def _depth_track_index_for_id(track_id: str, depth_track_ids: Sequence[str]) -> int:
    try:
        return list(depth_track_ids).index(track_id)
    except ValueError as exc:
        raise ValidationError(f"depth track_id '{track_id}' is not available.") from exc


def _slot_params_from_track_table_row(
    row: HeaderedTableRow,
    *,
    renderer: str,
    mode: Mode,
    depth_track_ids: Sequence[str],
) -> dict[str, Any]:
    params: dict[str, Any] = {}
    track_id = row.cell("track_id")
    raw_track_index = row.cell("track_index")
    if track_id:
        if renderer != "depth":
            raise table_error(
                "track_table",
                "track_id is only valid for depth renderer rows",
                row_number=row.row_number,
                column="track_id",
            )
        if not depth_track_ids:
            raise table_error(
                "track_table",
                "track_id requires depth input",
                row_number=row.row_number,
                column="track_id",
            )
        if track_id not in depth_track_ids:
            raise table_error(
                "track_table",
                f"unknown depth track_id '{track_id}'",
                row_number=row.row_number,
                column="track_id",
            )
        params["track_id"] = track_id
    if raw_track_index:
        parsed_index = parse_optional_int_cell(
            raw_track_index,
            "track_index",
            nonnegative=True,
            table_name="track_table",
            row_number=row.row_number,
            column="track_index",
        )
        if parsed_index is not None:
            params["track_index"] = parsed_index
    if track_id:
        derived_index = _depth_track_index_for_id(track_id, depth_track_ids)
        if "track_index" in params and int(params["track_index"]) != derived_index:
            raise table_error(
                "track_table",
                f"track_index={params['track_index']} does not point at track_id '{track_id}'",
                row_number=row.row_number,
                column="track_index",
            )
        params["track_index"] = derived_index

    nt = row.cell("nt")
    dinucleotide = row.cell("dinucleotide")
    if nt and dinucleotide and nt.upper() != dinucleotide.upper():
        raise table_error(
            "track_table",
            "nt and dinucleotide specify different values",
            row_number=row.row_number,
            column="dinucleotide",
        )
    if nt or dinucleotide:
        params["nt"] = (nt or dinucleotide).upper()
    source_index = parse_optional_int_cell(
        row.cell("source_index"),
        "source_index",
        nonnegative=True,
        table_name="track_table",
        row_number=row.row_number,
        column="source_index",
    )
    if source_index is not None:
        params["source_index"] = source_index
    for column in ("lane_direction", "tick_label_layout"):
        value = row.cell(column)
        if value:
            params[column] = value
    return params


def _read_track_table_rows(path: str) -> list[HeaderedTableRow]:
    rows = read_headered_tsv_table(
        path,
        required=_TRACK_TABLE_REQUIRED,
        optional=_TRACK_TABLE_OPTIONAL,
        table_name="track_table",
    )
    return sorted(
        rows,
        key=lambda row: (
            parse_optional_int_cell(
                row.cell("order"),
                "order",
                table_name="track_table",
                row_number=row.row_number,
                column="order",
            )
            if row.cell("order")
            else row.row_number,
            row.row_number,
        ),
    )


def _circular_width_from_depth_metadata(
    *,
    slot: CircularTrackSlot,
    metadata_by_track_id: Mapping[str, DepthTrackMetadata] | None,
) -> CircularTrackSlot:
    if slot.width is not None or str(slot.renderer) != "depth":
        return slot
    track_id = str((slot.params or {}).get("track_id", "")).strip()
    if not track_id or not metadata_by_track_id:
        return slot
    width = metadata_by_track_id.get(track_id, DepthTrackMetadata()).width
    if width is None:
        return slot
    return replace(slot, width=ScalarSpec(float(width), "px"))


def load_track_table_slots(
    path: str,
    *,
    mode: Mode,
    axis_before: str | None = None,
    depth_track_ids: Sequence[str] = (),
    depth_metadata_by_track_id: Mapping[str, DepthTrackMetadata] | None = None,
) -> TrackTableResult:
    """Convert --track_table rows directly into track slot dataclasses."""

    rows = _read_track_table_rows(path)
    seen_slot_ids: set[str] = set()
    slots: list[LinearTrackSlot] | list[CircularTrackSlot]
    linear_slots: list[LinearTrackSlot] = []
    circular_slots: list[CircularTrackSlot] = []
    supported = SUPPORTED_LINEAR_TRACK_RENDERERS if mode == "linear" else SUPPORTED_CIRCULAR_TRACK_RENDERERS

    for row in rows:
        slot_id = _require_nonempty(row, "slot_id", "track_table")
        if slot_id in seen_slot_ids:
            raise table_error("track_table", f"duplicate slot_id '{slot_id}'", row_number=row.row_number, column="slot_id")
        seen_slot_ids.add(slot_id)
        renderer = _normalize_renderer(_require_nonempty(row, "renderer", "track_table"), mode=mode)
        if renderer not in supported:
            raise table_error(
                "track_table",
                f"unknown {mode} track renderer '{renderer}'",
                row_number=row.row_number,
                column="renderer",
            )
        enabled = parse_optional_bool_cell(
            row.cell("enabled"),
            "enabled",
            table_name="track_table",
            row_number=row.row_number,
            column="enabled",
        )
        z = parse_optional_int_cell(
            row.cell("z"),
            "z",
            table_name="track_table",
            row_number=row.row_number,
            column="z",
        )
        params = _slot_params_from_track_table_row(
            row,
            renderer=renderer,
            mode=mode,
            depth_track_ids=depth_track_ids,
        )
        side = row.cell("side") or None

        if mode == "linear":
            if row.cell("radius"):
                raise table_error("track_table", "radius is only valid in circular mode", row_number=row.row_number, column="radius")
            if row.cell("width"):
                raise table_error("track_table", "width is only valid in circular mode", row_number=row.row_number, column="width")
            height = _parse_linear_px_scalar(
                row.cell("height"),
                field_name="height",
                table_name="track_table",
                row_number=row.row_number,
                allow_zero=False,
            )
            spacing = _parse_linear_px_scalar(
                row.cell("spacing"),
                field_name="spacing",
                table_name="track_table",
                row_number=row.row_number,
                allow_zero=True,
            )
            linear_slots.append(
                LinearTrackSlot(
                    id=slot_id,
                    renderer=renderer,
                    enabled=True if enabled is None else bool(enabled),
                    side=side,
                    height=height,
                    spacing=spacing,
                    z=0 if z is None else int(z),
                    params=params,
                )
            )
        else:
            if row.cell("height"):
                raise table_error("track_table", "height is only valid in linear mode", row_number=row.row_number, column="height")
            slot = CircularTrackSlot(
                id=slot_id,
                renderer=renderer,
                enabled=True if enabled is None else bool(enabled),
                side=side,
                radius=_parse_scalar_spec(row.cell("radius"), field_name="radius", table_name="track_table", row_number=row.row_number),
                width=_parse_scalar_spec(row.cell("width"), field_name="width", table_name="track_table", row_number=row.row_number),
                spacing=_parse_scalar_spec(row.cell("spacing"), field_name="spacing", table_name="track_table", row_number=row.row_number),
                z=0 if z is None else int(z),
                params=params,
            )
            circular_slots.append(
                _circular_width_from_depth_metadata(
                    slot=slot,
                    metadata_by_track_id=depth_metadata_by_track_id,
                )
            )

    if mode == "linear":
        slots = [slot for slot in linear_slots if slot.enabled]
    else:
        slots = [slot for slot in circular_slots if slot.enabled]
    if not slots:
        raise ValidationError("track_table: no enabled rows.")

    axis_index = None
    if axis_before:
        normalized_axis_before = str(axis_before).strip()
        for index, slot in enumerate(slots):
            if slot.id == normalized_axis_before:
                axis_index = index
                break
        if axis_index is None:
            raise ValidationError(
                f"track_table: --track_table_axis_before '{normalized_axis_before}' does not name an enabled slot_id."
            )
    return TrackTableResult(slots=slots, axis_index=axis_index)


def apply_depth_track_ids_to_slots(
    slots: Sequence[LinearTrackSlot] | Sequence[CircularTrackSlot] | None,
    *,
    depth_track_ids: Sequence[str],
    depth_metadata_by_track_id: Mapping[str, DepthTrackMetadata] | None = None,
) -> list[LinearTrackSlot] | list[CircularTrackSlot] | None:
    """Attach stable depth track IDs to depth slots generated from legacy order helpers."""

    if slots is None or not depth_track_ids:
        return list(slots) if slots is not None else None
    out: list[LinearTrackSlot | CircularTrackSlot] = []
    for slot in slots:
        if str(slot.renderer) != "depth":
            out.append(slot)
            continue
        params = dict(slot.params or {})
        raw_track_id = str(params.get("track_id", "") or "").strip()
        raw_index = params.get("track_index", 0)
        try:
            track_index = int(raw_index or 0)
        except (TypeError, ValueError):
            track_index = 0
        if not raw_track_id and 0 <= track_index < len(depth_track_ids):
            params["track_id"] = depth_track_ids[track_index]
        elif raw_track_id and raw_track_id in depth_track_ids:
            params["track_index"] = list(depth_track_ids).index(raw_track_id)
        updated = replace(slot, params=params)
        if isinstance(updated, CircularTrackSlot):
            updated = _circular_width_from_depth_metadata(
                slot=updated,
                metadata_by_track_id=depth_metadata_by_track_id,
            )
        out.append(updated)
    return out  # type: ignore[return-value]


__all__ = [
    "DepthTrackMetadata",
    "DepthTrackTableResult",
    "InputTableRow",
    "TrackTableResult",
    "apply_depth_track_ids_to_slots",
    "load_depth_track_table",
    "load_input_table_records",
    "load_track_table_slots",
]
