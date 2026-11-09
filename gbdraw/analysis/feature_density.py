"""Feature density tracks: how many features of a query lie in each window along a record, or what share of the
window's bases they cover.

A track is data, like a depth track: an id, a :class:`~gbdraw.features.query.FeatureQuery`, a metric, a drawing
style and its window. A ``feature_density`` slot names its track by id (``@track=<id>``), so the query is never
written inside the slot text.

Values come from the sorted ends of the matching features' parts at the window bounds of
:func:`gbdraw.analysis.windows.sliding_windows` (``numpy.searchsorted``), never from an array of every base:

- ``count``: the features whose parts overlap the window; a feature of several parts (a join) counts once.
- ``coverage``: the percentage of the window's bases inside the union of the matching features' parts.

A circular record's windows wrap round its origin; on a linear record or a crop the windows stop at the record's
ends, so a window there reads fewer bases (a crop does not hold the features beyond its ends).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, cast, get_args

import numpy as np
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray

from gbdraw.analysis.windows import sliding_windows
from gbdraw.core.color import interpolate_color, normalize_hex_color
from gbdraw.core.numeric import parse_ascii_integer
from gbdraw.exceptions import InputFileError, ValidationError
from gbdraw.features.query import FeatureQuery, format_bases
from gbdraw.io.colors import resolve_color_to_hex
from gbdraw.io.table_text import read_table_lines
from gbdraw.layout.record_coordinates import DisplayFragment, RecordDisplayTransform, SourceInterval

FeatureDensityMetric = Literal["count", "coverage"]
FeatureDensityStyle = Literal["heatmap", "histogram"]
FEATURE_DENSITY_METRICS: tuple[str, ...] = get_args(FeatureDensityMetric)
FEATURE_DENSITY_STYLES: tuple[str, ...] = get_args(FeatureDensityStyle)

# The fields of one track, as the CLI's ``key=value`` words and the track table's columns name them.
FEATURE_DENSITY_TRACK_FIELDS: tuple[str, ...] = (
    "id",
    "type",
    "qualifier",
    "value",
    "min_length",
    "max_length",
    "metric",
    "style",
    "window",
    "step",
    "color",
    "caption",
)
_FIELD_ALIASES = {"feature_type": "type"}
# The color steps of a heatmap band: fewer when a count track's largest value is smaller.
HEATMAP_LEVELS = 8
# The pale end of a heatmap band's colors and Legend gradient.
HEATMAP_LOW_COLOR = "#eef2f7"
_SLOT_ID_FORBIDDEN = frozenset(",@:=;\t\n ")


# The diagnostic of a refused feature density track: the Web words it as an invalid input.
_INVALID_TRACK = {"code": "INPUT_INVALID"}


@dataclass(frozen=True)
class FeatureDensityTrack:
    """One feature density track: which features (``query``), what per window (``metric``), how it is drawn
    (``style``), and its window and step in bases (``None``: the GC content window of the record)."""

    id: str
    query: FeatureQuery
    metric: FeatureDensityMetric = "count"
    style: FeatureDensityStyle = "heatmap"
    window: int | None = None
    step: int | None = None
    color: str | None = None
    caption: str | None = None

    def __post_init__(self) -> None:
        track_id = str(self.id or "").strip()
        if not track_id or _SLOT_ID_FORBIDDEN & set(track_id):
            raise ValidationError(
                f"A feature density track id must be a non-empty word without spaces or any of , @ : = ; (got {self.id!r}).",
                diagnostic=_INVALID_TRACK,
            )
        object.__setattr__(self, "id", track_id)
        if self.metric not in FEATURE_DENSITY_METRICS:
            raise ValidationError(
                f"Feature density track '{track_id}': metric must be one of {', '.join(FEATURE_DENSITY_METRICS)}.",
                diagnostic=_INVALID_TRACK,
            )
        if self.style not in FEATURE_DENSITY_STYLES:
            raise ValidationError(
                f"Feature density track '{track_id}': style must be one of {', '.join(FEATURE_DENSITY_STYLES)}.",
                diagnostic=_INVALID_TRACK,
            )
        for name in ("window", "step"):
            size = getattr(self, name)
            if size is not None and (isinstance(size, bool) or int(size) != size or int(size) <= 0):
                raise ValidationError(f"Feature density track '{track_id}': {name} must be a whole number > 0.", diagnostic=_INVALID_TRACK)
        if self.color is not None:
            # A heatmap's colors are steps from HEATMAP_LOW_COLOR to this hex: a hex of 3 or 6 digits or a color name.
            color = str(self.color).strip()
            try:
                normalize_hex_color(resolve_color_to_hex(color))
            except (TypeError, ValueError, ValidationError) as error:
                raise ValidationError(
                    f"Feature density track '{track_id}': color {color!r} is not an SVG color name or a hex of 3 or "
                    "6 digits.",
                    diagnostic={"code": "INPUT_INVALID", "field": "color"},
                ) from error
            object.__setattr__(self, "color", color)
        caption = str(self.caption).strip() if self.caption is not None else ""
        object.__setattr__(self, "caption", caption or None)


def _optional_int(value: str | None, *, field_name: str, track: str) -> int | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return parse_ascii_integer(text)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"Feature density track '{track}': {field_name} must be a whole number, not {text!r}.", diagnostic=_INVALID_TRACK) from error


def feature_density_track_from_fields(fields: Mapping[str, str], *, where: str) -> FeatureDensityTrack:
    """One track from its fields by name (:data:`FEATURE_DENSITY_TRACK_FIELDS`); ``where`` names the input."""

    normalized: dict[str, str] = {}
    for raw_key, raw_value in fields.items():
        key = _FIELD_ALIASES.get(str(raw_key).strip().lower(), str(raw_key).strip().lower())
        if key not in FEATURE_DENSITY_TRACK_FIELDS:
            raise ValidationError(
                f"{where}: unknown feature density field '{raw_key}'. Use {', '.join(FEATURE_DENSITY_TRACK_FIELDS)}.",
                diagnostic=_INVALID_TRACK,
            )
        normalized[key] = str(raw_value).strip()
    track_id = normalized.get("id", "")
    if not track_id:
        raise ValidationError(f"{where}: a feature density track needs an id (id=...).", diagnostic=_INVALID_TRACK)
    if not normalized.get("type"):
        raise ValidationError(f"{where}: feature density track '{track_id}' needs a feature type (type=..., or * for every type).", diagnostic=_INVALID_TRACK)
    query = FeatureQuery(
        normalized["type"],
        normalized.get("qualifier") or None,
        normalized.get("value") or None,
        _optional_int(normalized.get("min_length"), field_name="min_length", track=track_id),
        _optional_int(normalized.get("max_length"), field_name="max_length", track=track_id),
    )
    return FeatureDensityTrack(
        id=track_id,
        query=query,
        metric=cast(FeatureDensityMetric, normalized.get("metric") or "count"),
        style=cast(FeatureDensityStyle, normalized.get("style") or "heatmap"),
        window=_optional_int(normalized.get("window"), field_name="window", track=track_id),
        step=_optional_int(normalized.get("step"), field_name="step", track=track_id),
        color=normalized.get("color") or None,
        caption=normalized.get("caption") or None,
    )


def parse_feature_density_words(words: Sequence[str], *, where: str) -> FeatureDensityTrack:
    """One track from ``key=value`` words (``--feature_density_track id=genes type=gene``)."""

    fields: dict[str, str] = {}
    for word in words:
        key, separator, value = str(word).partition("=")
        if not separator:
            raise ValidationError(f"{where}: '{word}' is not a key=value word.", diagnostic=_INVALID_TRACK)
        if key.strip().lower() in fields:
            raise ValidationError(f"{where}: '{key}' is given twice.", diagnostic=_INVALID_TRACK)
        fields[key.strip().lower()] = value
    return feature_density_track_from_fields(fields, where=where)


def read_feature_density_table(path: str) -> tuple[FeatureDensityTrack, ...]:
    """The tracks of a tab-separated table whose first line names its columns
    (:data:`FEATURE_DENSITY_TRACK_FIELDS`; ``id`` and ``type`` are required)."""

    try:
        lines = read_table_lines(path)
    except FileNotFoundError as error:
        raise InputFileError(f"Feature density table not found: {error}") from error
    rows = [(number, line.rstrip("\r\n")) for number, line in enumerate(lines, start=1) if line.strip()]
    if not rows:
        raise ValidationError(f"Feature density table '{path}' has no header line.", diagnostic=_INVALID_TRACK)
    header = [cell.strip() for cell in rows[0][1].split("\t")]
    tracks = []
    for number, line in rows[1:]:
        cells = line.split("\t")
        if len(cells) > len(header):
            raise ValidationError(f"Feature density table '{path}' line {number}: more cells than columns.", diagnostic=_INVALID_TRACK)
        tracks.append(
            feature_density_track_from_fields(dict(zip(header, cells)), where=f"Feature density table '{path}' line {number}")
        )
    return unique_feature_density_tracks(tracks)


def unique_feature_density_tracks(tracks: Iterable[FeatureDensityTrack]) -> tuple[FeatureDensityTrack, ...]:
    """The tracks in order; two tracks of one id are refused."""

    result = tuple(tracks)
    seen: set[str] = set()
    for track in result:
        if track.id in seen:
            raise ValidationError(f"Two feature density tracks have the id '{track.id}'.", diagnostic=_INVALID_TRACK)
        seen.add(track.id)
    return result


@dataclass(frozen=True)
class FeatureDensityValues:
    """One track's values on one record: the value of each window (``samples`` as in the GC tracks) and the
    cells a heatmap band paints, ``bounds[k] .. bounds[k + 1]`` around sample ``k``."""

    samples: NDArray[np.int64]
    values: NDArray[np.float64]
    bounds: NDArray[np.int64]


def _feature_parts(
    record: SeqRecord, query: FeatureQuery, *, wrap: bool
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """The merged parts of each matching feature: starts, ends, and the index of the feature each part is of,
    sorted by feature and then by start. A wrapping record's parts appear once more a record length before
    and after (its windows reach past the origin)."""

    length = len(record.seq)
    record_id = str(record.id)
    starts: list[int] = []
    ends: list[int] = []
    owners: list[int] = []
    shifts = (-length, 0, length) if wrap else (0,)
    owner = 0
    for feature in record.features:
        if feature.type == "source" or not query.matches(feature, record_id=record_id):
            continue
        location = feature.location
        if location is None:
            continue
        parts = sorted(
            (max(0, int(part.start)) + shift, min(length, int(part.end)) + shift)
            for part in location.parts
            for shift in shifts
        )
        merged: list[list[int]] = []
        for start, end in parts:
            if end <= start:
                continue
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for start, end in merged:
            starts.append(start)
            ends.append(end)
            owners.append(owner)
        owner += 1
    return (
        np.asarray(starts, dtype=np.int64),
        np.asarray(ends, dtype=np.int64),
        np.asarray(owners, dtype=np.int64),
    )


def _window_counts(
    window_starts: NDArray[np.int64],
    window_ends: NDArray[np.int64],
    starts: NDArray[np.int64],
    ends: NDArray[np.int64],
    owners: NDArray[np.int64],
) -> NDArray[np.float64]:
    """How many features overlap each window, a feature of several parts once.

    The windows' starts and ends both rise, so the windows a part overlaps are one index range. The ranges of one
    feature's parts that contain a window are consecutive parts, so the window counts that feature once when each
    pair of consecutive parts both containing it takes one away.
    """

    count = len(window_starts)
    diff = np.zeros(count + 1, dtype=np.int64)
    if len(starts) == 0:
        return np.zeros(count, dtype=np.float64)
    first = np.searchsorted(window_ends, starts, side="right")
    stop = np.searchsorted(window_starts, ends, side="left")
    hit = first < stop
    np.add.at(diff, first[hit], 1)
    np.add.at(diff, stop[hit], -1)
    same = (owners[1:] == owners[:-1]) & hit[1:] & hit[:-1]
    pair_first, pair_stop = first[1:], stop[:-1]
    shared = same & (pair_first < pair_stop)
    np.add.at(diff, pair_first[shared], -1)
    np.add.at(diff, pair_stop[shared], 1)
    return np.cumsum(diff[:-1]).astype(np.float64)


def _union(starts: NDArray[np.int64], ends: NDArray[np.int64]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """The union of intervals as sorted disjoint intervals."""

    if len(starts) == 0:
        return starts, ends
    order = np.argsort(starts, kind="stable")
    starts, ends = starts[order], ends[order]
    reach = np.maximum.accumulate(ends)
    opens = np.ones(len(starts), dtype=bool)
    opens[1:] = starts[1:] > reach[:-1]
    # Each union interval ends where the reach stands at its last interval.
    lasts = np.flatnonzero(np.append(opens[1:], True))
    return starts[opens], reach[lasts]


def _covered_bases(
    window_starts: NDArray[np.int64],
    window_ends: NDArray[np.int64],
    starts: NDArray[np.int64],
    ends: NDArray[np.int64],
) -> NDArray[np.float64]:
    """The bases of each window inside the union of the intervals."""

    union_starts, union_ends = _union(starts, ends)
    if len(union_starts) == 0:
        return np.zeros(len(window_starts), dtype=np.float64)
    below = np.concatenate(([0], np.cumsum(union_ends - union_starts)))

    def covered_below(points: NDArray[np.int64]) -> NDArray[np.int64]:
        index = np.searchsorted(union_starts, points, side="right") - 1
        at = np.maximum(index, 0)
        inside = np.clip(points - union_starts[at], 0, union_ends[at] - union_starts[at])
        return np.where(index >= 0, below[at] + inside, 0)

    return (covered_below(window_ends) - covered_below(window_starts)).astype(np.float64)


def feature_density_values(
    record: SeqRecord,
    track: FeatureDensityTrack,
    *,
    window: int,
    step: int,
    wraps: bool,
) -> FeatureDensityValues:
    """The values of ``track`` on ``record`` with ``window`` and ``step`` (the track's own when it sets them).

    ``wraps`` says whether the windows wrap round the record's origin (``spans_may_cross_origin``: a whole record
    drawn as circular); otherwise they stop at the record's ends.
    """

    length = len(record.seq)
    window = int(track.window or window)
    step = int(track.step or step)
    wrap = wraps
    if wrap:
        # A window cannot read more than the whole circle.
        window = min(window, length)
    windows = sliding_windows(record, window, step, circular=wrap)
    window_starts, window_ends = windows.starts, windows.ends
    if not wrap:
        window_starts = np.clip(window_starts, 0, length)
        window_ends = np.clip(window_ends, 0, length)
    starts, ends, owners = _feature_parts(record, track.query, wrap=wrap)
    if track.metric == "count":
        values = _window_counts(window_starts, window_ends, starts, ends, owners)
    else:
        bases = np.maximum(window_ends - window_starts, 1)
        values = 100.0 * _covered_bases(window_starts, window_ends, starts, ends) / bases
    samples = windows.samples
    middles = (samples[1:] + samples[:-1]) // 2
    bounds = np.concatenate(([0], middles, [length])).astype(np.int64)
    return FeatureDensityValues(samples, values, bounds)


def feature_density_scale_max(values: Iterable[FeatureDensityValues]) -> float:
    """The top of a track's scale over the records it is drawn on: its largest value, at least 1 for a count and
    above 0 for a coverage, so an empty track draws nothing rather than dividing by 0."""

    top = max((float(np.max(item.values)) for item in values if len(item.values)), default=0.0)
    return top if top > 0 else 1.0


@dataclass(frozen=True)
class FeatureDensityLayer:
    """One track as drawn on one record: its values, the window and step they were read with, and the top of the
    scale the record's figure shares."""

    track: FeatureDensityTrack
    window: int
    step: int
    values: FeatureDensityValues
    scale_max: float

    @property
    def levels(self) -> int:
        """The color steps of the heatmap band: :data:`HEATMAP_LEVELS`, or a count track's largest value when that
        is smaller (each count its own color)."""

        if self.track.metric == "count":
            return max(1, min(HEATMAP_LEVELS, math.ceil(self.scale_max)))
        return HEATMAP_LEVELS

    def cell_levels(self) -> NDArray[np.int64]:
        """Each cell's color step, 0 (no fill) for a window without features up to :attr:`levels`."""

        fraction = np.clip(self.values.values / float(self.scale_max), 0.0, 1.0)
        return np.where(self.values.values > 0, np.maximum(1, np.ceil(fraction * self.levels)), 0).astype(np.int64)

    def runs(self) -> tuple[tuple[int, int, int], ...]:
        """The heatmap cells with equal adjacent cells merged: (start, end, level) in local positions, cells of
        level 0 left out."""

        levels = self.cell_levels()
        bounds = self.values.bounds
        if len(levels) == 0:
            return ()
        changes = np.flatnonzero(levels[1:] != levels[:-1]) + 1
        firsts = np.concatenate(([0], changes))
        lasts = np.concatenate((changes, [len(levels)]))
        return tuple(
            (int(bounds[first]), int(bounds[last]), int(levels[first]))
            for first, last in zip(firsts, lasts)
            if levels[first] > 0 and bounds[last] > bounds[first]
        )

    def scale_text(self, value: float) -> str:
        if self.track.metric == "coverage":
            return f"{value:.3g}%"
        return f"{value:g}"

    def caption(self) -> str:
        """The Legend name of the track: its own caption, else its query in words and what it counts per window;
        a histogram, which draws no axis, adds the top of its scale."""

        text = self.track.caption or (
            f"{self.track.query.describe()} per {format_bases(self.window)}"
            if self.track.metric == "count"
            else f"{self.track.query.describe()}, % of {format_bases(self.window)}"
        )
        if self.track.style == "histogram":
            return f"{text} (max {self.scale_text(self.scale_max)})"
        return text


def feature_density_layers(
    records: Sequence[SeqRecord],
    tracks: Sequence[FeatureDensityTrack],
    *,
    windows: Sequence[tuple[int, int]],
    wraps: Sequence[bool],
) -> list[dict[str, FeatureDensityLayer]]:
    """Every track on every record of one figure, by record and then track id; a track's scale is shared by the
    records drawn together. ``windows`` gives each record's (window, step) for tracks that set none, and ``wraps``
    whether its windows wrap round its origin."""

    per_record = [
        {
            track.id: (
                feature_density_values(record, track, window=window, step=step, wraps=record_wraps),
                int(track.window or window),
                int(track.step or step),
            )
            for track in tracks
        }
        for record, (window, step), record_wraps in zip(records, windows, wraps, strict=True)
    ]
    scales = {track.id: feature_density_scale_max(row[track.id][0] for row in per_record) for track in tracks}
    return [
        {
            track.id: FeatureDensityLayer(track, row[track.id][1], row[track.id][2], row[track.id][0], scales[track.id])
            for track in tracks
        }
        for row in per_record
    ]


@dataclass(frozen=True)
class DisplayCell:
    """One heatmap run on the display axis: its row (0 unwrapped), its display bounds and color step."""

    row: int
    start: int
    end: int
    level: int


def feature_density_display_cells(
    layer: FeatureDensityLayer, record_transform: RecordDisplayTransform | None
) -> tuple[DisplayCell, ...]:
    """The heatmap runs of ``layer`` where the record's transform draws them (through the same projection the
    feature parts take, ``RecordDisplayTransform.project_local_parts``)."""

    runs = layer.runs()
    if record_transform is None or not record_transform.projects:
        return tuple(DisplayCell(0, start, end, level) for start, end, level in runs)
    cells: list[DisplayCell] = []
    for start, end, level in runs:
        for fragment in record_transform.project_local_parts([SourceInterval(start, end)]):
            assert isinstance(fragment, DisplayFragment)
            cells.append(DisplayCell(fragment.row, fragment.display_start, fragment.display_end, level))
    return tuple(cells)


def heatmap_level_colors(color: str, levels: int) -> tuple[str, ...]:
    """The fill of each color step 1..``levels``: from a pale tint of ``color`` up to ``color`` itself."""

    top = resolve_color_to_hex(color)
    return tuple(interpolate_color(HEATMAP_LOW_COLOR, top, level / levels) for level in range(1, levels + 1))


__all__ = [
    "DisplayCell",
    "FEATURE_DENSITY_METRICS",
    "FEATURE_DENSITY_STYLES",
    "FEATURE_DENSITY_TRACK_FIELDS",
    "FeatureDensityLayer",
    "FeatureDensityMetric",
    "FeatureDensityStyle",
    "FeatureDensityTrack",
    "FeatureDensityValues",
    "HEATMAP_LEVELS",
    "HEATMAP_LOW_COLOR",
    "feature_density_display_cells",
    "feature_density_layers",
    "feature_density_scale_max",
    "feature_density_track_from_fields",
    "feature_density_values",
    "heatmap_level_colors",
    "parse_feature_density_words",
    "read_feature_density_table",
    "unique_feature_density_tracks",
]
