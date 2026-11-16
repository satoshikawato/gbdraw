"""The Linear drawer of a numeric series track (GC content, GC skew, depth): one drawer for every kind.

``draw_linear_series`` paints one row of a ``NumericSeriesPlot`` (the whole track of a record that does not wrap).
The axis faces left from the track's start; the Linear left room reads the same axis value (``ScalarAxis``).
"""

from __future__ import annotations

from gbdraw.svg.elements import ClipPath, Group, Line, Path, Rect, Text

from ....layout.linear_coords import normalize_position_to_linear_track
from ....layout.record_coordinates import RecordDisplayTransform
from ....layout.scalar_axis import (
    LINEAR_SCALAR_AXIS_LABEL_GAP_PX,
    LINEAR_SCALAR_AXIS_STROKE_PX,
    LINEAR_SCALAR_AXIS_TICK_PX,
    ScalarAxis,
    ScalarRow,
    linear_scalar_axis_tick_font_size_px,
)
from ....svg.ids import stable_svg_id
from ....svg.linear_tracks import LinearTrackFrame, linear_series_path_desc
from ....svg.path_data import path_point
from ...numeric_series import AXIS_STROKE, NumericSeriesPlot, SeriesBandFrame

SMALL_TICK_PX = 2.0


def _track_end_x(frame: LinearTrackFrame) -> float:
    return normalize_position_to_linear_track(frame.record_len, frame.record_len, frame.alignment_width)


def _band_rect(band: SeriesBandFrame, frame: LinearTrackFrame, group_id: str, *, border: bool) -> Rect:
    width = max(0.0, float(_track_end_x(frame)) - float(frame.start_x))
    insert, size = (frame.start_x, frame.start_y), (width, frame.track_height)
    if border:
        return Rect(id=f"{group_id}_percent_border", insert=insert, size=size, fill="none", stroke=band.border_color,
                    stroke_width=band.border_width)
    return Rect(id=f"{group_id}_percent_background", insert=insert, size=size, fill=band.background_color,
                fill_opacity=band.background_opacity, stroke="none")


def _axis_line(start: tuple[float, float], end: tuple[float, float]) -> Line:
    return Line(start=start, end=end, stroke=AXIS_STROKE, stroke_width=LINEAR_SCALAR_AXIS_STROKE_PX)


def _axis_text(text: str, x: float, y: float, font_size: float, font_family: str) -> Text:
    return Text(text, insert=(x, y), fill=AXIS_STROKE, font_size=font_size, font_family=font_family,
                text_anchor="end", dominant_baseline="central")


def _axis_group(axis: ScalarAxis, frame: LinearTrackFrame, axis_group_id: str) -> Group:
    """The axis: the baseline and the axis line up the track's start, its small and large ticks to the left, the
    numbers of the large ticks and the range label left of them."""

    start_x, start_y, _record_len, _width, track_height = frame
    baseline_y = float(start_y) + float(track_height)
    axis_group = Group(id=axis_group_id)
    axis_group.add(_axis_line((start_x, baseline_y), (_track_end_x(frame), baseline_y)))
    axis_group.add(_axis_line((start_x, start_y), (start_x, baseline_y)))
    font_size = linear_scalar_axis_tick_font_size_px(axis, track_height)
    labels = axis.labels(track_height, font_size)
    label_x = start_x - LINEAR_SCALAR_AXIS_TICK_PX - LINEAR_SCALAR_AXIS_LABEL_GAP_PX
    for value in axis.small_tick_values(labels):
        tick_y = baseline_y - (float(track_height) * axis.fraction(value))
        axis_group.add(_axis_line((start_x - SMALL_TICK_PX, tick_y), (start_x, tick_y)))
    for tick in labels.ticks:
        tick_y = baseline_y - tick.offset_px
        axis_group.add(_axis_line((start_x - LINEAR_SCALAR_AXIS_TICK_PX, tick_y), (start_x, tick_y)))
        if tick.label is not None:
            axis_group.add(_axis_text(tick.label, label_x, tick_y, font_size, axis.font_family))
    if labels.range_label is not None:
        axis_group.add(
            _axis_text(labels.range_label, label_x, float(start_y) + 0.5 * float(track_height), font_size, axis.font_family)
        )
    return axis_group


def _series_paths(plot: NumericSeriesPlot, desc: str, frame: LinearTrackFrame, group_id: str,
                  row: ScalarRow | None) -> list[ClipPath | Path]:
    """The series path, and for a series with a fill below its baseline its clip and the clipped copy."""

    paint = plot.paint
    style = {"stroke": paint.stroke, "fill_opacity": paint.fill_opacity, "fill_rule": "evenodd"}
    if paint.writes_stroke_width("linear"):
        style["stroke_width"] = paint.stroke_width
    main = Path(d=desc, fill=paint.fill, **style)
    if paint.below_fill is None:
        return [main]
    start_x, start_y, record_len, alignment_width, track_height = frame
    clip_id = stable_svg_id(
        f"{group_id}_clipper_line", "linear-skew-clip", group_id, plot.clip_key, record_len, alignment_width,
        track_height, start_x, start_y, desc,
        # Two rows may draw the same path; each later row refers to its own clip.
        *((row.row,) if row is not None and row.row else ()),
    )
    clip = ClipPath(id=clip_id)
    clip.add(
        Path(
            d=(
                f"M{path_point(start_x, start_y)} "
                f"L{path_point(alignment_width, start_y)} "
                f"L{path_point(alignment_width, start_y + track_height)} "
                f"L{path_point(start_x, start_y + track_height)} z"
            ),
            fill="white",
            stroke="none",
        )
    )
    below = Path(d=desc, fill=paint.below_fill, clip_path=f"url(#{clip_id})", clip_rule="nonzero", **style)
    return [clip, main, below]


def draw_linear_series(
    group: Group,
    plot: NumericSeriesPlot,
    frame: LinearTrackFrame,
    group_id: str,
    record_transform: RecordDisplayTransform | None = None,
    row: ScalarRow | None = None,
) -> Group:
    """Draw ``plot`` into ``group``: the band background, the series path (and its part below the baseline), the
    band border, and the axis."""

    if plot.samples is None:
        return group
    desc = linear_series_path_desc(plot.samples, frame, record_transform, row)
    if not desc:
        return group
    if plot.frame is not None:
        group.add(_band_rect(plot.frame, frame, group_id, border=False))
    for element in _series_paths(plot, desc, frame, group_id, row):
        group.add(element)
    if plot.frame is not None and plot.frame.border_width > 0:
        group.add(_band_rect(plot.frame, frame, group_id, border=True))
    if plot.axis is not None:
        group.add(_axis_group(plot.axis, frame, f"{group_id}_axis"))
    return group


__all__ = ["draw_linear_series"]
