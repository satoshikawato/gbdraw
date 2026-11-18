"""The Circular drawer of a numeric series track (GC content, GC skew, depth): one drawer for every kind.

``draw_circular_series`` paints a ``NumericSeriesPlot`` and ``circular_series_ink`` returns what that paints, both
from the same value, so what is drawn and what the label placement reserves for it cannot drift apart.
"""

from __future__ import annotations

import math

from gbdraw.svg.elements import Circle, ClipPath, Group, Path

from ....auto_sizes import scalar_axis_tick_font_size_px
from ....layout.circular import RadialBand
from ....layout.circular_depth_axis import (
    DEPTH_AXIS_SMALL_TICK_SIZE_PX,
    DEPTH_AXIS_STROKE_WIDTH_PX,
    DEPTH_AXIS_TICK_SIZE_PX,
    CircularRangeLabel,
    ScalarAxisTick,
    circular_range_label,
    scalar_axis_ink,
)
from ....layout.circular_ink import CircularInk, CircularInkProfile, CircularInkRing, stroke_pad_px
from ....layout.record_coordinates import CirclePlace
from ....layout.scalar_axis import ScalarAxis
from ....svg.circular_arcs import circle_arc
from ....svg.circular_tracks import (
    CircularTrackPlace,
    circle_record_turns,
    circle_unit_arc,
    circle_unit_turns,
    circular_series_baseline_radius,
    circular_series_line,
    circular_series_path_desc,
    generate_circle_path_desc,
)
from ....svg.ids import stable_svg_id
from ...numeric_series import AXIS_STROKE, NumericSeriesPlot, SeriesBandFrame
from .scalar_axis import circular_range_label_elements, circular_scalar_axis_group


def _band_radii(place: CircularTrackPlace) -> tuple[float, float]:
    """The baseline and outer radius of a band series: the band is centred on the track's anchor circle."""

    baseline_radius = max(0.0, (float(place.radius) * float(place.norm_factor)) - (0.5 * float(place.track_width)))
    return baseline_radius, baseline_radius + float(place.track_width)


def _axis_circle(axis: ScalarAxis, place: CircularTrackPlace) -> CirclePlace | None:
    """The shared circle the axis follows: a record of a circle several records share draws the axis once and its
    baseline as its own border arc, when the axis is drawn once on such a circle; else none."""

    return place.record_transform.circle if axis.once_on_shared_circle else None


def _axis_ticks(
    axis: ScalarAxis, baseline_radius: float, place: CircularTrackPlace
) -> tuple[list[ScalarAxisTick], CircularRangeLabel | None]:
    """The axis ticks as drawn, small ticks first, then large ticks with their numbers; and the range label this
    record draws: on a circle several records share, the label lies on the first record that holds it
    (``on_shared_circle``), never over a gap between records, and only that record draws it, after its own band,
    so no later record paints over it."""

    track_width = place.track_width
    font_size = scalar_axis_tick_font_size_px("circular", axis.style, track_width)
    labels = axis.labels(track_width, font_size)
    ticks = [
        ScalarAxisTick(baseline_radius + float(track_width) * axis.fraction(value), DEPTH_AXIS_SMALL_TICK_SIZE_PX, None)
        for value in axis.small_tick_values(labels)
    ]
    ticks.extend(ScalarAxisTick(baseline_radius + tick.offset_px, DEPTH_AXIS_TICK_SIZE_PX, tick.label) for tick in labels.ticks)
    range_label = circular_range_label(labels.range_label, baseline_radius + float(track_width), font_size=font_size)
    circle = _axis_circle(axis, place)
    if range_label is not None and circle is not None:
        record_turns = circle_record_turns(circle)
        range_label = range_label.on_shared_circle(record_turns, axis.font_family)
        own_start, own_end = record_turns[circle.position]
        if range_label is not None and not own_start <= range_label.start_rad / (2.0 * math.pi) < own_end:
            range_label = None
    return ticks, range_label


def _frame_ink(frame: SeriesBandFrame, place: CircularTrackPlace) -> list[CircularInk]:
    baseline_radius, outer_radius = _band_radii(place)
    # The arc of the ring a record that shares its circle paints; the whole ring else.
    turns = circle_unit_turns(place.record_transform)
    inks: list[CircularInk] = []
    if float(place.track_width) > 0 and frame.paints_background:
        inks.append(CircularInkRing(RadialBand(baseline_radius, outer_radius), *turns))
    if frame.paints_border:
        half = frame.border_overhang_px
        inks.extend(CircularInkRing(RadialBand(edge - half, edge + half), *turns) for edge in (baseline_radius, outer_radius))
    return inks


def _axis_ink(axis: ScalarAxis, place: CircularTrackPlace) -> list[CircularInk]:
    baseline_radius, outer_radius = _band_radii(place)
    ticks, range_label = _axis_ticks(axis, baseline_radius, place)
    circle = _axis_circle(axis, place)
    if circle is None or circle.draws_scalar_axes:
        return list(
            scalar_axis_ink(
                baseline_radius,
                outer_radius,
                ticks,
                font_size=scalar_axis_tick_font_size_px("circular", axis.style, place.track_width),
                font_family=axis.font_family,
                range_label=range_label,
                baseline_circle=circle is None,
            )
        )
    return [] if range_label is None else [range_label.ink(axis.font_family)]


def circular_series_ink(plot: NumericSeriesPlot, place: CircularTrackPlace) -> tuple[CircularInk, ...]:
    """What ``draw_circular_series`` paints with the same arguments: the area between the baseline and the line,
    with half its stroke (a ``line`` series: the line and half its stroke); a deviation series also its baseline
    (with a stroke); a band series its frame when it shows and its axis, the labels measured in the figure's font
    as they are written."""

    if plot.samples is None:
        return ()
    line = circular_series_line(plot.samples, place)
    if line is None:
        return ()
    pad = stroke_pad_px(plot.paint.stroke, plot.paint.stroke_width)
    record_len = place.record_len
    if plot.samples.deviation is not None:
        base = float(place.radius) * float(place.norm_factor)
        area = CircularInkProfile.between(line.positions, line.radii, record_length=record_len, base_px=base, pad_px=pad)
        # The baseline closes along the record's own arc (``circle_unit_turns``), the whole circle for a record alone.
        baseline = CircularInkRing(RadialBand(base - pad, base + pad), *circle_unit_turns(place.record_transform))
        return (area, baseline) if pad > 0 else (area,)
    inks: list[CircularInk] = [
        CircularInkProfile.along(line.positions, line.radii, record_length=record_len, pad_px=pad)
        if plot.samples.style == "line"
        else CircularInkProfile.between(line.positions, line.radii, record_length=record_len,
                                        base_px=line.baseline_radius, pad_px=pad)
    ]
    if plot.frame is not None:
        inks.extend(_frame_ink(plot.frame, place))
    if plot.axis is not None:
        inks.extend(_axis_ink(plot.axis, place))
    return tuple(inks)


def _add_frame_background(group: Group, frame: SeriesBandFrame, place: CircularTrackPlace, group_id: str) -> None:
    baseline_radius, _outer_radius = _band_radii(place)
    middle_radius = baseline_radius + (0.5 * float(place.track_width))
    unit_arc = circle_unit_arc(place.record_transform)
    style = {
        "id": f"{group_id}_percent_background",
        "fill": "none",
        "stroke": frame.background_color,
        "stroke_width": float(place.track_width),
        "stroke_opacity": frame.background_opacity,
    }
    if unit_arc is not None:
        # A record that shares its circle: the background is the record's own arc of the ring.
        group.add(Path(d=circle_arc(middle_radius, *unit_arc), **style))
    else:
        group.add(Circle(center=(0, 0), r=middle_radius, **style))


def _add_frame_border(group: Group, frame: SeriesBandFrame, place: CircularTrackPlace, group_id: str) -> None:
    if frame.border_width <= 0:
        return
    unit_arc = circle_unit_arc(place.record_transform)
    border_group = Group(id=f"{group_id}_percent_border")
    for border_radius in _band_radii(place):
        style = {"fill": "none", "stroke": frame.border_color, "stroke_width": frame.border_width}
        if unit_arc is not None:
            border_group.add(Path(d=circle_arc(border_radius, *unit_arc), **style))
        else:
            border_group.add(Circle(center=(0, 0), r=border_radius, **style))
    group.add(border_group)


def _add_axis(group: Group, axis: ScalarAxis, place: CircularTrackPlace, axis_group_id: str) -> None:
    baseline_radius, outer_radius = _band_radii(place)
    ticks, range_label = _axis_ticks(axis, baseline_radius, place)
    font_size = scalar_axis_tick_font_size_px("circular", axis.style, place.track_width)
    circle = _axis_circle(axis, place)
    if circle is None or circle.draws_scalar_axes:
        group.add(
            circular_scalar_axis_group(
                axis_group_id,
                baseline_radius,
                outer_radius,
                ticks,
                range_label,
                stroke=AXIS_STROKE,
                stroke_width=DEPTH_AXIS_STROKE_WIDTH_PX,
                font_size=font_size,
                font_family=axis.font_family,
                # A record of a shared circle draws its baseline as its own border arc.
                baseline_circle=circle is None,
            )
        )
        return
    if range_label is None:
        return
    # The range label of the circle's axis, drawn by a record other than the axis's when its arc holds it.
    label_group = Group(id=axis_group_id)
    for element in circular_range_label_elements(
        range_label, axis_group_id, stroke=AXIS_STROKE, font_size=font_size, font_family=axis.font_family
    ):
        label_group.add(element)
    group.add(label_group)


def _series_paths(plot: NumericSeriesPlot, desc: str, place: CircularTrackPlace, group_id: str,
                  record_id: str | None) -> list[ClipPath | Path]:
    """The series path, and for a series with a paint of its own below its baseline the clip of the disc inside the
    baseline and the clipped copy in that paint."""

    paint = plot.paint
    style = {"stroke": paint.stroke, "fill_opacity": paint.fill_opacity, "fill_rule": "evenodd"}
    main = Path(d=desc, fill=paint.fill, **style)
    elements: list[ClipPath | Path] = [main]
    if plot.samples is not None and paint.paints_below:
        baseline_radius = circular_series_baseline_radius(plot.samples, place)
        clip_id = stable_svg_id(
            "clipper_circle", "circular-skew-clip", group_id, record_id, place.record_len, place.track_width,
            place.norm_factor, plot.clip_key, desc,
            # A line's path does not hold its baseline: two baselines are two clips.
            *(() if plot.samples.deviation is not None else (baseline_radius,)),
        )
        clip = ClipPath(id=clip_id)
        clip.add(Path(d=generate_circle_path_desc(baseline_radius, 1.0), fill="white", stroke="none"))
        below = Path(d=desc, fill=paint.low_fill, clip_path=f"url(#{clip_id})", clip_rule="nonzero",
                     **{**style, "stroke": paint.low_stroke})
        elements = [clip, main, below]
    if paint.writes_stroke_width:
        for element in elements:
            if isinstance(element, Path):
                element["stroke-width"] = paint.stroke_width
    return elements


def draw_circular_series(
    group: Group, plot: NumericSeriesPlot, place: CircularTrackPlace, group_id: str, record_id: str | None
) -> Group:
    """Draw ``plot`` into ``group``: the band background, the series path (and its part below the baseline), the
    band border, and the axis, which a record of a shared circle draws once or each record draws (the axis's
    rule)."""

    if plot.samples is None:
        return group
    desc = circular_series_path_desc(plot.samples, place)
    if not desc:
        return group
    if plot.frame is not None:
        _add_frame_background(group, plot.frame, place, group_id)
    for element in _series_paths(plot, desc, place, group_id, record_id):
        group.add(element)
    if plot.frame is not None:
        _add_frame_border(group, plot.frame, place, group_id)
    if plot.axis is not None:
        _add_axis(group, plot.axis, place, f"{group_id}_axis")
    return group


__all__ = ["circular_series_ink", "draw_circular_series"]
