"""Shared helpers for depth coverage track inputs and precomputation."""

from __future__ import annotations

import copy
from collections.abc import Sequence as SequenceABC
from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Literal, Sequence, TypeVar

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.analysis.depth import (
    DepthAggregation,
    DepthCoverage,
    DepthFileFormat,
    depth_df,
    parse_depth_tsv,
    read_cached_depth_coverage,
)
from gbdraw.configurators import DepthConfigurator
from gbdraw.core.numeric import ascii_float
from gbdraw.exceptions import GbdrawError, ValidationError
from gbdraw.io.colors import check_user_color
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import bind_record_transforms
from gbdraw.layout.series_style import SeriesBaseline, SeriesStyle
from gbdraw.legend.sources import DepthLegendSource, LegendSlotSource
from gbdraw.render_stages import render_stage
from gbdraw.tracks.parsing import depth_slot_track_index, slot_legend_label

if TYPE_CHECKING:
    from gbdraw.api.options import DepthTrackInput
    from gbdraw.legend.table import LegendTable


# A bound of one series' value range: a number, ``auto`` (from its values, whatever the global setting), or None
# (the global ``min_depth``/``max_depth``).
DepthSeriesBound = float | Literal["auto"] | None


def series_bound(bound: DepthSeriesBound, global_bound: float | None) -> float | None:
    """A series' bound as its drawing config holds it: its own number, None for ``auto``, else the global one."""

    if bound is None:
        return global_bound
    return None if bound == "auto" else float(bound)


@dataclass(frozen=True)
class DepthTrackStyle:
    """The fields of a depth track spec its drawing config reads (not its id, label or source). The fields after
    the ticks default to how every depth track drew before they existed."""

    fill_color: str | None
    large_tick_interval: float | None
    small_tick_interval: float | None
    tick_font_size: float | None
    series_style: SeriesStyle = "area"
    negative_color: str | None = None
    baseline: SeriesBaseline = "auto"
    min_value: DepthSeriesBound = None
    max_value: DepthSeriesBound = None


@dataclass(frozen=True)
class DepthTrackSpec:
    """One depth source for one displayed record and logical track.

    ``table`` is a parsed file (its reader takes the file format), or an API table that ``depth_df`` converts. ``agg`` is the
    value of a window; ``style``, ``negative_color``, ``baseline`` and the value range ``min_value``/``max_value``
    (``DepthSeriesBound``) say how the series draws. Every default draws as before these fields existed.
    """

    id: str
    label: str
    table: DepthCoverage | DataFrame
    track_index: int = 0
    fill_color: str | None = None
    height: float | None = None
    large_tick_interval: float | None = None
    small_tick_interval: float | None = None
    tick_font_size: float | None = None
    agg: DepthAggregation = "mean"
    style: SeriesStyle = "area"
    negative_color: str | None = None
    baseline: SeriesBaseline = "auto"
    min_value: DepthSeriesBound = None
    max_value: DepthSeriesBound = None

    @property
    def drawing_style(self) -> DepthTrackStyle:
        return DepthTrackStyle(
            fill_color=self.fill_color,
            large_tick_interval=self.large_tick_interval,
            small_tick_interval=self.small_tick_interval,
            tick_font_size=self.tick_font_size,
            series_style=self.style,
            negative_color=self.negative_color,
            baseline=self.baseline,
            min_value=self.min_value,
            max_value=self.max_value,
        )


@dataclass(frozen=True)
class DepthTrackData:
    """One precomputed depth dataframe and its effective drawing config."""

    id: str
    label: str
    df: DataFrame
    config: DepthConfigurator
    track_index: int = 0
    height: float | None = None


DepthDfBuilder = Callable[..., DataFrame]
# Record-major normalized specs a builder read once (``freeze_depth_tracks``).
ResolvedDepthTracks = tuple[tuple[DepthTrackSpec, ...], ...]
DepthTrack = TypeVar("DepthTrack", DepthTrackSpec, DepthTrackData, DepthLegendSource)


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


def _expand_track_float_metadata(
    values: Sequence[float | str | None] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[float | None]:
    if track_count <= 0:
        return []
    if values is None:
        return [None for _ in range(track_count)]
    items: list[float | None] = []
    for value in values:
        if value is None:
            items.append(None)
            continue
        if isinstance(value, str) and value.strip().lower() in {"", "auto", "none", "null", "-"}:
            items.append(None)
            continue
        numeric = ascii_float(value)
        if numeric <= 0:
            raise ValidationError(f"{field_name} values must be > 0.")
        items.append(numeric)
    if len(items) == 1:
        return [items[0] for _ in range(track_count)]
    if len(items) != track_count:
        raise ValidationError(
            f"{field_name} count must be one or equal to the number of depth tracks ({track_count})."
        )
    return items


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


def _read_depth_table(path: str, file_format: DepthFileFormat = "auto") -> DepthCoverage:
    """Read one depth file in ``file_format``; in a Worker transaction, give its coverage its resource identity.

    There the decoded-resource layer owns the coverage, one per resource;
    elsewhere (CLI, Python API) the process cache does. Either way every record
    and render shares it read-only.
    """

    # Deferred: gbdraw.api imports this module.
    from gbdraw.api.prepared import get_or_build_decoded_resource, prepared_resource_identity

    if prepared_resource_identity(path) is None:
        return read_cached_depth_coverage(path, file_format)
    coverage: DepthCoverage = get_or_build_decoded_resource(
        path, ("depth-coverage-v1", file_format), lambda: parse_depth_tsv(path, file_format)
    )
    return coverage


def _name_series(exc: GbdrawError, series_index: int) -> None:
    """Let a failure of one logical series name the series, not its file."""

    if exc.diagnostic is not None:
        exc.diagnostic = {**exc.diagnostic, "seriesIndex": series_index}


def _read_depth_series(path: str, series_index: int, file_format: DepthFileFormat = "auto") -> DepthCoverage:
    """Read one logical series source; failures name the series, not the file."""

    try:
        return _read_depth_table(path, file_format)
    except GbdrawError as exc:
        _name_series(exc, series_index)
        raise


def _tracks_from_record_major_values(
    values: Sequence[Sequence[object | None]],
    *,
    records: Sequence[SeqRecord],
    labels: Sequence[str] | None,
    colors: Sequence[str] | None,
    heights: Sequence[float | str | None] | None,
    large_tick_intervals: Sequence[float | str | None] | None,
    small_tick_intervals: Sequence[float | str | None] | None,
    tick_font_sizes: Sequence[float | str | None] | None,
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
    track_heights = _expand_track_float_metadata(
        heights,
        track_count=track_count,
        field_name="depth_track_heights",
    )
    track_large_tick_intervals = _expand_track_float_metadata(
        large_tick_intervals,
        track_count=track_count,
        field_name="depth_track_large_tick_intervals",
    )
    track_small_tick_intervals = _expand_track_float_metadata(
        small_tick_intervals,
        track_count=track_count,
        field_name="depth_track_small_tick_intervals",
    )
    track_tick_font_sizes = _expand_track_float_metadata(
        tick_font_sizes,
        track_count=track_count,
        field_name="depth_track_tick_font_sizes",
    )
    record_tracks: list[list[DepthTrackSpec]] = [[] for _ in range(record_count)]
    loaded_files: dict[str, DepthCoverage] = {}

    for track_index in range(track_count):
        raw_track_values = [
            row[track_index] if track_index < len(row) else None
            for row in rows
        ]
        present = [(index, value) for index, value in enumerate(raw_track_values) if value is not None]
        if not present:
            raise ValidationError(f"Depth track {track_index + 1} is empty.")

        if single_shared_row:
            source_values = [present[0][1] for _ in range(record_count)]
        else:
            source_values = list(raw_track_values)

        for record_index, source_value in enumerate(source_values):
            if source_value is None:
                continue
            table: DepthCoverage | DataFrame
            if values_are_files:
                source_path = str(source_value)
                if source_path not in loaded_files:
                    loaded_files[source_path] = _read_depth_series(source_path, track_index)
                table = loaded_files[source_path]
            elif isinstance(source_value, (DataFrame, DepthCoverage)):
                # A coverage is the file of ``depth_file``, read before this call.
                table = source_value
            else:
                raise ValidationError(f"Depth track {track_index + 1} must contain pandas DataFrame values.")
            record_tracks[record_index].append(
                DepthTrackSpec(
                    id=track_ids[track_index],
                    label=str(track_labels[track_index] or _default_depth_track_labels(track_count)[track_index]),
                    table=table,
                    track_index=track_index,
                    fill_color=track_colors[track_index],
                    height=track_heights[track_index],
                    large_tick_interval=track_large_tick_intervals[track_index],
                    small_tick_interval=track_small_tick_intervals[track_index],
                    tick_font_size=track_tick_font_sizes[track_index],
                )
            )

    return record_tracks


def _tracks_from_logical_inputs(
    inputs: Sequence[DepthTrackInput],
    *,
    records: Sequence[SeqRecord],
) -> list[list[DepthTrackSpec]]:
    """Expand canonical logical tracks into sparse record-major specs."""

    from gbdraw.api.options import DepthTrackInput

    record_count = len(records)
    track_count = len(inputs)
    record_tracks: list[list[DepthTrackSpec]] = [[] for _ in range(record_count)]
    loaded_files: dict[str, DepthCoverage] = {}
    track_ids = _default_depth_track_ids(track_count)
    default_labels = _default_depth_track_labels(track_count)

    for track_index, track in enumerate(inputs):
        if not isinstance(track, DepthTrackInput):
            raise ValidationError(
                f"depth_tracks[{track_index}] must be a DepthTrackInput."
            )
        source = track.source
        if isinstance(source, SequenceABC) and not isinstance(
            source,
            (str, bytes, bytearray),
        ) and not isinstance(source, DataFrame):
            source_values = list(source)
            if len(source_values) != record_count:
                raise ValidationError(
                    f"depth_tracks[{track_index}].source must contain one source "
                    f"per displayed record ({record_count}); got {len(source_values)}."
                )
        else:
            source_values = [source for _ in range(record_count)]

        if not any(value is not None for value in source_values):
            raise ValidationError(f"Depth track {track_index + 1} is empty.")

        for record_index, source_value in enumerate(source_values):
            if source_value is None:
                continue
            table: DepthCoverage | DataFrame
            if isinstance(source_value, DataFrame):
                table = source_value
            elif isinstance(source_value, (str, Path)):
                source_path = str(source_value)
                if source_path not in loaded_files:
                    loaded_files[source_path] = _read_depth_series(source_path, track_index)
                table = loaded_files[source_path]
            else:
                raise ValidationError(
                    f"depth_tracks[{track_index}].source[{record_index}] must be "
                    "a path, DataFrame, or None."
                )
            record_tracks[record_index].append(
                DepthTrackSpec(
                    id=track_ids[track_index],
                    label=str(track.label or default_labels[track_index]),
                    table=table,
                    track_index=track_index,
                    fill_color=track.color,
                    height=track.height,
                    large_tick_interval=track.large_tick_interval,
                    small_tick_interval=track.small_tick_interval,
                    tick_font_size=track.tick_font_size,
                )
            )

    return record_tracks


def normalize_depth_tracks(
    records: Sequence[SeqRecord],
    *,
    depth_tracks: Sequence[DepthTrackInput] | None = None,
    depth_table: DataFrame | None = None,
    depth_file: str | None = None,
    depth_tables: Sequence[DataFrame] | None = None,
    depth_files: Sequence[str] | None = None,
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None,
    depth_track_files: Sequence[Sequence[str | None]] | None = None,
    depth_track_labels: Sequence[str] | None = None,
    depth_track_colors: Sequence[str] | None = None,
    depth_track_heights: Sequence[float | str | None] | None = None,
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None,
) -> list[list[DepthTrackSpec]] | None:
    """Normalize all public depth inputs into record-major depth track specs."""

    record_count = len(records)
    if record_count <= 0:
        raise ValidationError("records is empty")

    singular_plural_present = _has_any(
        (depth_table, depth_file, depth_tables, depth_files)
    )
    matrix_present = _has_any(
        (
            depth_track_tables,
            depth_track_files,
            depth_track_labels,
            depth_track_colors,
            depth_track_heights,
            depth_track_large_tick_intervals,
            depth_track_small_tick_intervals,
            depth_track_tick_font_sizes,
        )
    )
    compatibility_present = singular_plural_present or matrix_present
    if depth_tracks is not None and compatibility_present:
        raise ValidationError(
            "depth_tracks cannot be combined with compatibility depth inputs."
        )
    if singular_plural_present and matrix_present:
        raise ValidationError(
            "Singular/plural depth inputs cannot be combined with depth_track_* "
            "compatibility inputs."
        )
    if depth_tracks is None and not compatibility_present:
        return None

    if depth_tracks is not None:
        if isinstance(depth_tracks, (str, bytes)) or not isinstance(
            depth_tracks,
            SequenceABC,
        ):
            raise ValidationError(
                "depth_tracks must be a sequence of DepthTrackInput values."
            )
        if not depth_tracks:
            raise ValidationError(
                "depth_tracks must include at least one logical depth track."
            )
        return _tracks_from_logical_inputs(depth_tracks, records=records)

    if matrix_present:
        if depth_track_tables is not None and depth_track_files is not None:
            raise ValidationError("Pass either depth_track_tables or depth_track_files, not both.")
        if depth_track_tables is None and depth_track_files is None:
            raise ValidationError("depth_track metadata requires depth_track_tables or depth_track_files.")
        if depth_track_tables is not None:
            return _tracks_from_record_major_values(
                depth_track_tables,
                records=records,
                labels=depth_track_labels,
                colors=depth_track_colors,
                heights=depth_track_heights,
                large_tick_intervals=depth_track_large_tick_intervals,
                small_tick_intervals=depth_track_small_tick_intervals,
                tick_font_sizes=depth_track_tick_font_sizes,
                field_name="depth_track_tables",
                values_are_files=False,
            )
        return _tracks_from_record_major_values(
            depth_track_files or (),
            records=records,
            labels=depth_track_labels,
            colors=depth_track_colors,
            heights=depth_track_heights,
            large_tick_intervals=depth_track_large_tick_intervals,
            small_tick_intervals=depth_track_small_tick_intervals,
            tick_font_sizes=depth_track_tick_font_sizes,
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
        table = depth_table if depth_table is not None else _read_depth_table(str(depth_file))
        return _tracks_from_record_major_values(
            [[table]],
            records=records,
            labels=None,
            colors=None,
            heights=None,
            large_tick_intervals=None,
            small_tick_intervals=None,
            tick_font_sizes=None,
            field_name="depth_table",
            values_are_files=False,
        )

    if depth_tables is not None:
        return _tracks_from_record_major_values(
            [list(depth_tables)] if len(depth_tables) == 1 else [[table] for table in depth_tables],
            records=records,
            labels=None,
            colors=None,
            heights=None,
            large_tick_intervals=None,
            small_tick_intervals=None,
            tick_font_sizes=None,
            field_name="depth_tables",
            values_are_files=False,
        )

    return _tracks_from_record_major_values(
        [list(depth_files or ())] if len(depth_files or ()) == 1 else [[path] for path in (depth_files or ())],
        records=records,
        labels=None,
        colors=None,
        heights=None,
        large_tick_intervals=None,
        small_tick_intervals=None,
        tick_font_sizes=None,
        field_name="depth_files",
        values_are_files=True,
    )


def freeze_depth_tracks(
    record_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
) -> ResolvedDepthTracks | None:
    """Normalized specs, frozen for assembly inputs (``_resolved_depth_tracks``).

    A builder that assembles a figure more than once reads each depth file
    once; on the Web its staged files are gone after the first assembly.
    """

    if record_tracks is None:
        return None
    return tuple(tuple(row) for row in record_tracks)


def thaw_depth_tracks(resolved: ResolvedDepthTracks | None) -> list[list[DepthTrackSpec]] | None:
    """The builder's resolved specs as an assembly reads them; None when it resolved none.

    An assembly reads ``thaw_depth_tracks(inputs._resolved_depth_tracks) or
    normalize_depth_tracks(...)``: a resolved value holds a row per record, so
    it is never empty.
    """

    return None if resolved is None else [list(row) for row in resolved]


def index_depth_track_row(
    row: Sequence[DepthTrack],
    *,
    record_index: int | None = None,
) -> dict[int, DepthTrack]:
    """Index one sparse record row by logical depth track index."""

    indexed: dict[int, DepthTrack] = {}
    record_context = "" if record_index is None else f" for record {record_index}"
    for track in row:
        raw_index = track.track_index
        if isinstance(raw_index, bool) or not isinstance(raw_index, Integral):
            raise ValidationError(
                f"Depth track{record_context} has non-integer track_index={raw_index!r}."
            )
        track_index = int(raw_index)
        if track_index < 0:
            raise ValidationError(
                f"Depth track{record_context} has negative track_index={track_index}."
            )
        if track_index in indexed:
            raise ValidationError(
                f"Depth track row{record_context} contains duplicate track_index={track_index}."
            )
        indexed[track_index] = track
    return indexed


def _logical_depth_track_count(
    record_depth_tracks: Sequence[Sequence[DepthTrack]] | None,
) -> int:
    indices = depth_track_indices(record_depth_tracks)
    return (indices[-1] + 1) if indices else 0


def depth_track_indices(
    record_depth_tracks: Sequence[Sequence[DepthTrack]] | None,
) -> tuple[int, ...]:
    """Return all present logical depth track indices in ascending order."""

    indices: set[int] = set()
    for record_index, row in enumerate(record_depth_tracks or ()):
        indexed = index_depth_track_row(row, record_index=record_index)
        indices.update(indexed)
    return tuple(sorted(indices))


def depth_track_count(record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None) -> int:
    return _logical_depth_track_count(record_depth_tracks)


def depth_track_heights(
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
) -> list[float | None]:
    heights: list[float | None] = [None for _ in range(depth_track_count(record_depth_tracks))]
    for record_index, row in enumerate(record_depth_tracks or ()):
        for track_index, spec in index_depth_track_row(row, record_index=record_index).items():
            if heights[track_index] is not None:
                continue
            if spec.height is not None:
                heights[track_index] = float(spec.height)
    return heights


def depth_track_data_count(record_depth_tracks: Sequence[Sequence[DepthTrackData]] | None) -> int:
    return _logical_depth_track_count(record_depth_tracks)


def representative_depth_tracks(
    record_depth_tracks: Sequence[Sequence[DepthTrackData]] | None,
) -> list[DepthTrackData]:
    """Return the first present data item for every logical depth track."""

    representatives: dict[int, DepthTrackData] = {}
    for record_index, row in enumerate(record_depth_tracks or ()):
        for track_index, track in index_depth_track_row(row, record_index=record_index).items():
            representatives.setdefault(track_index, track)
    return [representatives[index] for index in sorted(representatives)]


def _checked_series_color(color: str, field: str) -> str:
    # Every depth track color (CLI, options, Session, slot) is checked here.
    check_user_color(
        color,
        where="for a depth track",
        diagnostic={"code": "INPUT_INVALID", "field": field, "reason": "COLOR"},
    )
    return str(color)


def clone_depth_config(
    base_config: DepthConfigurator,
    *,
    fill_color: str | None = None,
    window: int | None = None,
    step: int | None = None,
    large_tick_interval: float | None = None,
    small_tick_interval: float | None = None,
    tick_font_size: float | None = None,
) -> DepthConfigurator:
    config = copy.copy(base_config)
    if fill_color:
        config.fill_color = _checked_series_color(fill_color, "color")
    if window is not None:
        config.window = int(window)
    if step is not None:
        config.step = int(step)
    if large_tick_interval is not None:
        config.large_tick_interval = float(large_tick_interval)
    if small_tick_interval is not None:
        config.small_tick_interval = float(small_tick_interval)
    if tick_font_size is not None:
        config.tick_font_size = float(tick_font_size)
    return config


def depth_track_config(
    style: DepthTrackStyle,
    base_config: DepthConfigurator,
    window: int,
    step: int,
    shared_bounds: SharedBounds | None,
) -> DepthConfigurator:
    """Build one track's drawing config (render stage ``depth_track_config``).

    It reads the track's style, not its label, so a label change keeps the
    config and the depth content drawn with it (OV-1229): its colour, ticks,
    value range (``series_bound``), series style, negative colour and baseline.
    ``shared_bounds`` are the auto bounds over every record's frame of this
    logical track when the axis is shared (``_shared_bounds``); the config is
    complete when it is returned, so nothing writes into it afterwards (M5).
    """

    config = clone_depth_config(
        base_config,
        fill_color=style.fill_color,
        window=window,
        step=step,
        large_tick_interval=style.large_tick_interval,
        small_tick_interval=style.small_tick_interval,
        tick_font_size=style.tick_font_size,
    )
    config.min_depth = series_bound(style.min_value, base_config.min_depth)
    config.max_depth = series_bound(style.max_value, base_config.max_depth)
    config.series_style = style.series_style
    config.negative_color = (
        None if style.negative_color is None else _checked_series_color(style.negative_color, "negativeColor")
    )
    config.baseline = style.baseline
    if shared_bounds is not None:
        shared_min, shared_max = shared_bounds
        if config.min_depth is None and shared_min is not None:
            config.min_depth = shared_min
        if config.max_depth is None:
            config.max_depth = shared_max
    return config


SharedBounds = tuple[float | None, float]


def _shared_bounds(frames: Sequence[Sequence[tuple[DepthTrackSpec, DataFrame]]]) -> dict[int, SharedBounds]:
    """The auto bounds of each logical track over every record's frame: its largest value (0 when none is above
    0), and its smallest value when one is below 0 (else None: the axis starts at 0 as before)."""

    smallest: dict[int, float] = {}
    largest: dict[int, float] = {}
    for row in frames:
        for spec, frame in row:
            if frame.empty or "depth" not in frame.columns:
                continue
            index = spec.track_index
            smallest[index] = min(smallest.get(index, 0.0), float(frame["depth"].min()))
            largest[index] = max(largest.get(index, 0.0), float(frame["depth"].max()))
    return {index: (smallest[index] if smallest[index] < 0 else None, largest[index]) for index in largest}


def _series_frame(
    depth_df_builder: DepthDfBuilder,
    record: SeqRecord,
    spec: DepthTrackSpec,
    window: int,
    step: int,
    *,
    circular: bool,
    base_config: DepthConfigurator,
) -> DataFrame:
    """One series' binned frame of one record (render stage ``depth_frame``); a failure names the series."""

    try:
        frame: DataFrame = render_stage(
            "depth_frame",
            depth_df_builder,
            record,
            spec.table,
            window,
            step,
            circular=circular,
            min_depth=series_bound(spec.min_value, base_config.min_depth),
            max_depth=series_bound(spec.max_value, base_config.max_depth),
            agg=spec.agg,
            log_scale=base_config.normalize,
        )
    except GbdrawError as exc:
        _name_series(exc, spec.track_index)
        raise
    return frame


def build_depth_track_dataframes(
    records: Sequence[SeqRecord],
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
    record_transforms: Sequence[RecordDisplayTransform] | None,
    *,
    base_config: DepthConfigurator | None,
    depth_df_builder: DepthDfBuilder = depth_df,
    window_steps: Sequence[tuple[int, int]] | None = None,
) -> list[list[DepthTrackData]]:
    """Build depth dataframes once for each record/track pair.

    Two render stages per record and track: ``depth_frame`` (the binned
    dataframe: the series' window value ``agg``, clipped to its value range)
    and ``depth_track_config`` (its drawing config, taking the shared axis
    bounds as an argument). A record's windows wrap round its ends when it is
    drawn as circular (its transform, ``record_transforms``). A log-scaled
    track refuses a series with a negative value, naming the series.
    """

    if base_config is None or not record_depth_tracks:
        return [[] for _ in records]
    if len(record_depth_tracks) != len(records):
        raise ValidationError(
            f"Expected depth tracks for {len(records)} record(s); got {len(record_depth_tracks)}."
        )

    record_window_steps: list[tuple[int, int]] = []
    frames: list[list[tuple[DepthTrackSpec, DataFrame]]] = []
    transforms = bind_record_transforms(records, record_transforms)
    for record_index, (record, track_specs) in enumerate(zip(records, record_depth_tracks)):
        index_depth_track_row(track_specs, record_index=record_index)
        if window_steps is not None and record_index < len(window_steps):
            window, step = window_steps[record_index]
        else:
            window, step = int(base_config.window), int(base_config.step)
        record_window_steps.append((int(window), int(step)))
        frames.append(
            [
                (
                    spec,
                    _series_frame(
                        depth_df_builder,
                        record,
                        spec,
                        int(window),
                        int(step),
                        circular=transforms[record_index].is_circular,
                        base_config=base_config,
                    ),
                )
                for spec in track_specs
            ]
        )

    shared_bounds = _shared_bounds(frames) if base_config.share_axis else {}
    return [
        [
            DepthTrackData(
                id=spec.id,
                label=spec.label,
                df=frame,
                config=render_stage(
                    "depth_track_config",
                    depth_track_config,
                    spec.drawing_style,
                    base_config,
                    window,
                    step,
                    shared_bounds.get(spec.track_index),
                ),
                track_index=spec.track_index,
                height=spec.height,
            )
            for spec, frame in row
        ]
        for row, (window, step) in zip(frames, record_window_steps)
    ]


def depth_track_legend_name(track: DepthLegendSource, track_count: int) -> str:
    """The Legend row name of a Depth series without a slot ``legend_label``: its label, else ``Depth`` (one
    series) or ``Depth N``."""

    return str(track.label or ("Depth" if track_count == 1 else f"Depth {int(track.track_index) + 1}"))


def sync_depth_track_legend_entries(
    out: LegendTable,
    depth_tracks: Sequence[DepthLegendSource] | None,
    *,
    slots: Sequence[LegendSlotSource],
) -> LegendTable:
    """Replace the singleton Depth legend with logical-track-aware entries.

    ``slots`` is authoritative: only enabled Depth slots are represented, in
    slot order, and their legend label overrides are honored.
    """

    if not depth_tracks:
        return out
    out.discard_default_track("Depth")
    track_count = _logical_depth_track_count([depth_tracks])
    tracks_by_index = index_depth_track_row(depth_tracks)
    for slot in slots:
        if slot.renderer != "depth":
            continue
        track = tracks_by_index.get(depth_slot_track_index(slot.id, slot.params))
        if track is None:
            continue
        out.add_numeric_source(slot_legend_label(slot, depth_track_legend_name(track, track_count)), track.paint)
    return out


__all__ = [
    "DepthTrackData",
    "DepthTrackSpec",
    "DepthSeriesBound",
    "DepthTrackStyle",
    "ResolvedDepthTracks",
    "build_depth_track_dataframes",
    "clone_depth_config",
    "depth_track_config",
    "depth_track_count",
    "depth_track_data_count",
    "depth_track_heights",
    "depth_track_legend_name",
    "depth_track_indices",
    "freeze_depth_tracks",
    "index_depth_track_row",
    "normalize_depth_tracks",
    "representative_depth_tracks",
    "series_bound",
    "sync_depth_track_legend_entries",
    "thaw_depth_tracks",
]
