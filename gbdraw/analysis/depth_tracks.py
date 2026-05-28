"""Shared helpers for depth coverage track inputs and precomputation."""

from __future__ import annotations

import copy
from collections.abc import Sequence as SequenceABC
from dataclasses import dataclass
from typing import Callable, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.depth import depth_df, read_depth_tsv  # type: ignore[reportMissingImports]
from gbdraw.configurators import DepthConfigurator  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]


@dataclass(frozen=True)
class DepthTrackSpec:
    """One normalized depth table for one displayed record and logical track."""

    id: str
    label: str
    table: DataFrame
    fill_color: str | None = None


@dataclass(frozen=True)
class DepthTrackData:
    """One precomputed depth dataframe and its effective drawing config."""

    id: str
    label: str
    df: DataFrame
    config: DepthConfigurator


DepthDfBuilder = Callable[..., DataFrame]


def _has_any(values: Sequence[object | None]) -> bool:
    return any(value is not None for value in values)


def _is_nested_row(value: object) -> bool:
    return isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray))


def _expand_track_metadata(
    values: Sequence[str] | None,
    *,
    track_count: int,
    default_values: Sequence[str | None],
    field_name: str,
) -> list[str | None]:
    if track_count <= 0:
        return []
    if values is None:
        return list(default_values)
    items = [str(value).strip() for value in values]
    if len(items) == 1:
        return [items[0] or None for _ in range(track_count)]
    if len(items) != track_count:
        raise ValidationError(
            f"{field_name} count must be one or equal to the number of depth tracks ({track_count})."
        )
    return [item or None for item in items]


def _default_depth_track_ids(track_count: int) -> list[str]:
    if track_count <= 0:
        return []
    if track_count == 1:
        return ["depth"]
    return [f"depth_{index + 1}" for index in range(track_count)]


def _default_depth_track_labels(track_count: int) -> list[str | None]:
    if track_count <= 0:
        return []
    if track_count == 1:
        return ["Depth"]
    return [f"Depth {index + 1}" for index in range(track_count)]


def _record_major_sources(
    sources: Sequence[Sequence[object | None]],
    *,
    record_count: int,
    field_name: str,
) -> tuple[list[list[object | None]], bool]:
    raw_rows = list(sources)
    if not raw_rows:
        raise ValidationError(f"{field_name} must include at least one logical depth track.")
    if not all(_is_nested_row(row) for row in raw_rows):
        raise ValidationError(f"{field_name} must be a record-major nested sequence.")
    rows = [list(row) for row in raw_rows]
    if len(rows) == 1 and record_count > 1:
        return rows, True
    if len(rows) != record_count:
        raise ValidationError(
            f"{field_name} must contain one record row or one row per record ({record_count}); got {len(rows)}."
        )
    return rows, False


def _tracks_from_record_major_values(
    values: Sequence[Sequence[object | None]],
    *,
    records: Sequence[SeqRecord],
    labels: Sequence[str] | None,
    colors: Sequence[str] | None,
    field_name: str,
    values_are_files: bool,
) -> list[list[DepthTrackSpec]]:
    record_count = len(records)
    rows, single_shared_row = _record_major_sources(
        values,
        record_count=record_count,
        field_name=field_name,
    )
    track_count = max((len(row) for row in rows), default=0)
    if track_count <= 0:
        raise ValidationError(f"{field_name} must include at least one logical depth track.")

    track_ids = _default_depth_track_ids(track_count)
    track_labels = _expand_track_metadata(
        labels,
        track_count=track_count,
        default_values=_default_depth_track_labels(track_count),
        field_name="depth_track_labels",
    )
    track_colors = _expand_track_metadata(
        colors,
        track_count=track_count,
        default_values=[None for _ in range(track_count)],
        field_name="depth_track_colors",
    )
    record_tracks: list[list[DepthTrackSpec]] = [[] for _ in range(record_count)]

    for track_index in range(track_count):
        raw_track_values = [
            row[track_index] if track_index < len(row) else None
            for row in rows
        ]
        present = [(index, value) for index, value in enumerate(raw_track_values) if value is not None]
        if not present:
            raise ValidationError(f"Depth track {track_index + 1} is empty.")

        if single_shared_row or len(present) == 1:
            source_values = [present[0][1] for _ in range(record_count)]
        elif len(present) == record_count:
            source_values = list(raw_track_values)
        else:
            raise ValidationError(
                f"Depth track {track_index + 1} must have one table/file or one per record "
                f"({record_count}); got {len(present)}."
            )

        loaded_table: DataFrame | None = None
        if values_are_files and len(present) == 1:
            loaded_table = read_depth_tsv(str(present[0][1]))

        for record_index, source_value in enumerate(source_values):
            if source_value is None:
                continue
            table = (
                loaded_table
                if loaded_table is not None
                else read_depth_tsv(str(source_value))
                if values_are_files
                else source_value
            )
            if not isinstance(table, DataFrame):
                raise ValidationError(f"Depth track {track_index + 1} must contain pandas DataFrame values.")
            record_tracks[record_index].append(
                DepthTrackSpec(
                    id=track_ids[track_index],
                    label=str(track_labels[track_index] or _default_depth_track_labels(track_count)[track_index]),
                    table=table,
                    fill_color=track_colors[track_index],
                )
            )

    return record_tracks


def normalize_depth_tracks(
    records: Sequence[SeqRecord],
    *,
    depth_table: DataFrame | None = None,
    depth_file: str | None = None,
    depth_tables: Sequence[DataFrame] | None = None,
    depth_files: Sequence[str] | None = None,
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None,
    depth_track_files: Sequence[Sequence[str | None]] | None = None,
    depth_track_labels: Sequence[str] | None = None,
    depth_track_colors: Sequence[str] | None = None,
) -> list[list[DepthTrackSpec]] | None:
    """Normalize all public depth inputs into record-major depth track specs."""

    record_count = len(records)
    if record_count <= 0:
        raise ValidationError("records is empty")

    legacy_present = _has_any((depth_table, depth_file, depth_tables, depth_files))
    new_present = _has_any(
        (depth_track_tables, depth_track_files, depth_track_labels, depth_track_colors)
    )
    if legacy_present and new_present:
        raise ValidationError("Legacy depth inputs cannot be combined with depth_track_* inputs.")
    if not legacy_present and not new_present:
        return None

    if new_present:
        if depth_track_tables is not None and depth_track_files is not None:
            raise ValidationError("Pass either depth_track_tables or depth_track_files, not both.")
        if depth_track_tables is None and depth_track_files is None:
            raise ValidationError("depth_track_labels/depth_track_colors require depth_track_tables or depth_track_files.")
        if depth_track_tables is not None:
            return _tracks_from_record_major_values(
                depth_track_tables,
                records=records,
                labels=depth_track_labels,
                colors=depth_track_colors,
                field_name="depth_track_tables",
                values_are_files=False,
            )
        return _tracks_from_record_major_values(
            depth_track_files or (),
            records=records,
            labels=depth_track_labels,
            colors=depth_track_colors,
            field_name="depth_track_files",
            values_are_files=True,
        )

    if depth_table is not None and depth_file is not None:
        raise ValidationError("Pass either depth_table or depth_file, not both.")
    if (depth_table is not None or depth_file is not None) and (
        depth_tables is not None or depth_files is not None
    ):
        raise ValidationError("Use depth_table/depth_file or depth_tables/depth_files, not both.")
    if depth_tables is not None and depth_files is not None:
        raise ValidationError("Pass either depth_tables or depth_files, not both.")

    if depth_table is not None or depth_file is not None:
        table = depth_table if depth_table is not None else read_depth_tsv(str(depth_file))
        return _tracks_from_record_major_values(
            [[table]],
            records=records,
            labels=None,
            colors=None,
            field_name="depth_table",
            values_are_files=False,
        )

    if depth_tables is not None:
        return _tracks_from_record_major_values(
            [list(depth_tables)] if len(depth_tables) == 1 else [[table] for table in depth_tables],
            records=records,
            labels=None,
            colors=None,
            field_name="depth_tables",
            values_are_files=False,
        )

    return _tracks_from_record_major_values(
        [list(depth_files or ())] if len(depth_files or ()) == 1 else [[path] for path in (depth_files or ())],
        records=records,
        labels=None,
        colors=None,
        field_name="depth_files",
        values_are_files=True,
    )


def depth_track_count(record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None) -> int:
    if not record_depth_tracks:
        return 0
    return max((len(row) for row in record_depth_tracks), default=0)


def depth_track_data_count(record_depth_tracks: Sequence[Sequence[DepthTrackData]] | None) -> int:
    if not record_depth_tracks:
        return 0
    return max((len(row) for row in record_depth_tracks), default=0)


def clone_depth_config(
    base_config: DepthConfigurator,
    *,
    fill_color: str | None = None,
    window: int | None = None,
    step: int | None = None,
) -> DepthConfigurator:
    config = copy.copy(base_config)
    if fill_color:
        config.fill_color = str(fill_color)
    if window is not None:
        config.window = int(window)
    if step is not None:
        config.step = int(step)
    return config


def build_depth_track_dataframes(
    records: Sequence[SeqRecord],
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
    *,
    base_config: DepthConfigurator | None,
    depth_df_builder: DepthDfBuilder = depth_df,
    window_steps: Sequence[tuple[int, int]] | None = None,
) -> list[list[DepthTrackData]]:
    """Build depth dataframes once for each record/track pair."""

    if base_config is None or not record_depth_tracks:
        return [[] for _ in records]
    if len(record_depth_tracks) != len(records):
        raise ValidationError(
            f"Expected depth tracks for {len(records)} record(s); got {len(record_depth_tracks)}."
        )

    output: list[list[DepthTrackData]] = []
    for record_index, (record, track_specs) in enumerate(zip(records, record_depth_tracks)):
        if window_steps is not None and record_index < len(window_steps):
            window, step = window_steps[record_index]
        else:
            window, step = int(base_config.window), int(base_config.step)
        row: list[DepthTrackData] = []
        for spec in track_specs:
            config = clone_depth_config(
                base_config,
                fill_color=spec.fill_color,
                window=int(window),
                step=int(step),
            )
            row.append(
                DepthTrackData(
                    id=spec.id,
                    label=spec.label,
                    df=depth_df_builder(
                        record,
                        spec.table,
                        int(window),
                        int(step),
                        normalize=bool(config.normalize),
                        min_depth=config.min_depth,
                        max_depth=config.max_depth,
                    ),
                    config=config,
                )
            )
        output.append(row)

    if bool(getattr(base_config, "share_axis", False)) and base_config.max_depth is None:
        track_count = max((len(row) for row in output), default=0)
        for track_index in range(track_count):
            max_values = [
                float(row[track_index].df["depth"].max())
                for row in output
                if track_index < len(row)
                and not row[track_index].df.empty
                and "depth" in row[track_index].df.columns
            ]
            if not max_values:
                continue
            shared_max = max(max_values)
            for row in output:
                if track_index < len(row):
                    row[track_index].config.max_depth = shared_max

    return output


def sync_depth_track_legend_entries(
    legend_table: dict,
    depth_tracks: Sequence[DepthTrackData] | None,
) -> dict:
    """Replace the singleton depth legend entry with depth-track-aware entries."""

    if not depth_tracks:
        return legend_table
    out = dict(legend_table)
    out.pop("Depth", None)
    for index, track in enumerate(depth_tracks):
        label = str(track.label or ("Depth" if len(depth_tracks) == 1 else f"Depth {index + 1}"))
        unique_label = label
        suffix = 2
        while unique_label in out:
            unique_label = f"{label} ({suffix})"
            suffix += 1
        out[unique_label] = {
            "type": "solid",
            "fill": track.config.fill_color,
            "stroke": track.config.stroke_color,
            "width": track.config.stroke_width,
        }
    return out


__all__ = [
    "DepthTrackData",
    "DepthTrackSpec",
    "build_depth_track_dataframes",
    "clone_depth_config",
    "depth_track_count",
    "depth_track_data_count",
    "normalize_depth_tracks",
    "sync_depth_track_legend_entries",
]
