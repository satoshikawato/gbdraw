#!/usr/bin/env python
# coding: utf-8

import math
from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from gbdraw.svg.elements import Circle


from gbdraw.layout.record_coordinates import CirclePlace, RecordDisplayTransform
from ..layout.scalar_axis import SeriesSamples, held_to_end, project_scalar_line
from .circular_arcs import circle_arc, circle_arc_to, circle_point, position_angle, span_angle
from .path_data import path_point
from .track_paths import TRACK_PATH_TOLERANCE_PX, track_line_to, track_polyline


class CircularTrackPlace(NamedTuple):
    """Where a Circular numeric series track draws: the canvas radius, the record's drawn length, the track's
    width, its anchor as a share of the radius, and the record's transform."""

    radius: float
    record_len: int
    track_width: float
    norm_factor: float
    record_transform: RecordDisplayTransform


class CircularAreaLine(NamedTuple):
    """The line a Circular numeric series path draws: display positions, their radii, the baseline radius the
    area closes on, and where each drawn segment after the first starts."""

    positions: NDArray[np.float64]
    radii: NDArray[np.float64]
    baseline_radius: float
    segment_starts: tuple[int, ...] = ()

class CircularDrawnLine(NamedTuple):
    """The points a Circular track path draws for its line: x and y, and where each drawn segment after the
    first starts."""

    xs: NDArray[np.float64]
    ys: NDArray[np.float64]
    segment_starts: tuple[int, ...]


def circular_drawn_line(
    positions: NDArray[np.float64],
    radii: NDArray[np.float64],
    total_len: float,
    segment_starts: tuple[int, ...] = (),
) -> CircularDrawnLine:
    """The drawn points of a Circular track line: between two samples of one drawn segment, the angle and the
    radius change linearly, so the line follows the circle.

    Where two samples are so far apart on the circle that the straight piece between them would leave that
    line by more than ``TRACK_PATH_TOLERANCE_PX`` (a short record has few windows), points of the line are
    added between them; samples closer than that are joined directly, as before (OV-2863). The track path
    then keeps only the points its tolerance needs (:mod:`gbdraw.svg.track_paths`).
    """

    count = len(positions)
    pieces = np.ones(max(count - 1, 0), dtype=np.intp)
    if count > 1:
        span = np.abs(np.diff(positions)) * (2.0 * math.pi / float(total_len))
        reach = np.maximum(np.abs(radii[:-1]), np.abs(radii[1:]))
        # The longest angle whose chord stays within the tolerance of the circle of the larger radius.
        cos_limit = np.clip(1.0 - TRACK_PATH_TOLERANCE_PX / np.maximum(reach, TRACK_PATH_TOLERANCE_PX), -1.0, 1.0)
        max_angle = 2.0 * np.arccos(cos_limit)
        needed = np.ceil(np.divide(span, max_angle, out=np.ones_like(span), where=max_angle > 0))
        pieces = np.maximum(needed, 1.0).astype(np.intp)
        breaks = np.asarray(segment_starts, dtype=np.intp)
        breaks = breaks[(breaks > 0) & (breaks < count)]
        pieces[breaks - 1] = 1
    offsets = np.concatenate((np.zeros(1, dtype=np.intp), np.cumsum(pieces, dtype=np.intp)))
    gap = np.repeat(np.arange(len(pieces), dtype=np.intp), pieces)
    along = (np.arange(len(gap), dtype=np.intp) - offsets[gap]) / pieces[gap]
    drawn_positions = np.concatenate((positions[gap] + (positions[gap + 1] - positions[gap]) * along, positions[-1:]))
    drawn_radii = np.concatenate((radii[gap] + (radii[gap + 1] - radii[gap]) * along, radii[-1:]))
    angles_rad = np.radians(360.0 * (drawn_positions / float(total_len)) - 90.0)
    return CircularDrawnLine(
        drawn_radii * np.cos(angles_rad),
        drawn_radii * np.sin(angles_rad),
        tuple(int(offsets[start]) for start in segment_starts),
    )


def circle_unit_arc(record_transform: RecordDisplayTransform) -> tuple[float, float] | None:
    """Where a record that shares its circle with other records lies on it: its start angle and its span in
    degrees, clockwise. ``None`` when the record has the circle to itself."""

    place = record_transform.circle
    if place is None:
        return None
    return position_angle(place.offset, place.length), span_angle(record_transform.drawn_length, place.length)


def circle_unit_turns(record_transform: RecordDisplayTransform) -> tuple[float, float]:
    """The arc a record paints, in turns clockwise from 12 o'clock: its own arc on a circle it shares with other
    records, else the whole circle."""

    arc = circle_unit_arc(record_transform)
    if arc is None:
        return 0.0, 1.0
    start = (arc[0] + 90.0) / 360.0
    return start, start + arc[1] / 360.0


def circle_record_turns(place: CirclePlace) -> tuple[tuple[float, float], ...]:
    """The arcs of every record drawn on the circle of ``place``, in turns clockwise from 12 o'clock and in drawing
    order, so a record's own arc is at its ``position``."""

    return tuple((offset / place.length, (offset + length) / place.length) for offset, length in place.record_spans)


def within_record_arc(first_turn: float, last_turn: float, record_transform: RecordDisplayTransform | None) -> bool:
    """Whether text painted from ``first_turn`` to ``last_turn`` (turns clockwise from 12 o'clock) lies within its
    record's arc, so that it does not run over the gap onto a neighbour (OV-3113). A record alone holds any."""

    if record_transform is None or record_transform.circle is None:
        return True
    start, end = circle_unit_turns(record_transform)
    return start <= first_turn and last_turn <= end


def _circle_unit_samples(
    positions: ArrayLike, values: ArrayLike, record_transform: RecordDisplayTransform
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[int, ...]]:
    """A track's samples of a record that shares a circle, at their places on it: held on to the record's end
    there, as a Linear track holds to the end of its axis, since nothing joins one record's end to its start."""

    projected = project_scalar_line(positions, values, record_transform)
    if record_transform.start_coordinate is not None:
        return projected.positions, projected.values, projected.segment_starts
    place = record_transform.circle
    assert place is not None
    held_positions, held_values = held_to_end(
        projected.positions, projected.values, float(place.offset + record_transform.drawn_length)
    )
    return held_positions, held_values, projected.segment_starts


def generate_circle_path_desc(radius: float, norm_factor: float) -> str:
    """The closed circle of ``radius * norm_factor`` from 12 o'clock, clockwise (two arcs)."""
    return circle_arc(radius * norm_factor, -90.0, 360.0) + "z"


def _closes_round_its_origin(record_transform: RecordDisplayTransform, record_len: int) -> bool:
    """Whether a Circular line runs on from its last sample round to its first: the record drawn is a whole
    circular molecule. A crop of a circular record ends at two places on the molecule, so its line, like a linear
    record's, holds its last value to 12 o'clock and steps there (D-64); its transform's length is the source's."""

    return record_transform.is_circular and record_len == record_transform.length


def _closed_line_desc(line: CircularAreaLine, total_len: float) -> str:
    """The closed line of a GC content (deviation) or skew path: it starts at 12 o'clock on the radius of the
    first sample, runs through every sample, and closes back to 12 o'clock along the circle; the baseline
    circle follows as its own subpath. The line of a record that is not a whole circular molecule already holds
    its last value to 12 o'clock (``held_to_end``, ``_closes_round_its_origin``), so its close is a radial step
    there, not a ramp across bases (OV-2920)."""

    positions, radii, _baseline_radius, segment_starts = line
    drawn = circular_drawn_line(
        np.concatenate(([0.0], positions)),
        np.concatenate((radii[:1], radii)),
        total_len,
        tuple(start + 1 for start in segment_starts),
    )
    # The close runs from the last sample to 12 o'clock on the first radius; ``z`` draws its last piece, and
    # the points between are their own drawn segment so the last sample stays a kept end of the line.
    close = circular_drawn_line(np.array([positions[-1], total_len]), np.array([radii[-1], radii[0]]), total_len)
    xs = np.concatenate((drawn.xs, close.xs[1:-1]))
    ys = np.concatenate((drawn.ys, close.ys[1:-1]))
    closing_start = (len(drawn.xs),) if len(close.xs) > 2 else ()
    return track_polyline(xs, ys, (*drawn.segment_starts, *closing_start)) + "z"


def _placed_samples(
    positions: NDArray[np.float64], values: ArrayLike, place: CircularTrackPlace, closes_at_record_end: bool
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[int, ...]]:
    """A track's samples at their display positions: on the record's arc of a shared circle, projected by a
    display start, held to the end of a record that is not a whole circular molecule, else as they are, with the
    first sample repeated at the record's end when the track closes there."""

    transform = place.record_transform
    if transform.circle is not None:
        return _circle_unit_samples(positions, values, transform)
    if transform.start_coordinate is not None:
        return project_scalar_line(positions, values, transform)
    if not _closes_round_its_origin(transform, place.record_len):
        return (*held_to_end(positions, values, float(place.record_len)), ())
    value_array = np.asarray(values, dtype=float)
    if closes_at_record_end and len(positions) > 0 and np.isclose(positions[0], 0.0):
        positions = np.append(positions, float(place.record_len))
        value_array = np.append(value_array, value_array[0])
    return positions, value_array, ()


def _band_bottom_radius(place: CircularTrackPlace) -> float:
    """The bottom of a band series: half the track's width inside the circle at ``radius * norm_factor``."""

    return max(0.0, (float(place.radius) * float(place.norm_factor)) - (0.5 * float(place.track_width)))


def circular_series_baseline_radius(samples: SeriesSamples, place: CircularTrackPlace) -> float:
    """The radius of the baseline a numeric series path closes on, and below which (inside it) the part below the
    baseline is painted: the circle at ``radius * norm_factor`` for a deviation series, else the band's bottom
    raised by the band fraction ``samples.baseline``."""

    if samples.deviation is not None:
        return place.radius * place.norm_factor
    return _band_bottom_radius(place) + (float(place.track_width) * samples.baseline)


def circular_series_line(samples: SeriesSamples, place: CircularTrackPlace) -> CircularAreaLine | None:
    """The line a numeric series path draws: display positions, their radii and the baseline radius.

    A deviation series lies about the circle at ``radius * norm_factor``, its largest deviation half the track's
    width from it (a window without a value on the circle, OV-2799); a band series rises from the band's bottom,
    half the track's width inside that circle, by the fraction of the width each value fills, and closes on its
    baseline (``circular_series_baseline_radius``)."""

    scale = samples.deviation
    positions = samples.positions.to_numpy(dtype=float)
    if scale is not None:
        if place.record_len <= 0:
            return None
        norm_radius: float = circular_series_baseline_radius(samples, place)
        diff = (samples.values - scale.mean).fillna(0.0)
        if scale.max_diff == 0:
            radii = np.full(len(samples.values), norm_radius)
        else:
            radii = (norm_radius + (0.5 * place.track_width * (diff / scale.max_diff))).to_numpy(dtype=float)
        placed = _placed_samples(positions, radii, place, samples.closes_at_record_end)
        return CircularAreaLine(placed[0], placed[1], norm_radius, placed[2])
    if place.record_len <= 0 or place.track_width <= 0:
        return None
    bottom_radius = _band_bottom_radius(place)
    values = np.clip(samples.values.to_numpy(dtype=float), 0.0, 1.0)
    positions, values, segment_starts = _placed_samples(positions, values, place, samples.closes_at_record_end)
    return CircularAreaLine(
        positions,
        bottom_radius + (float(place.track_width) * values),
        circular_series_baseline_radius(samples, place),
        segment_starts,
    )


def _area_desc(line: CircularAreaLine, record_len: int) -> str:
    """The area between ``line`` and its baseline circle: from the baseline at the first sample along the line to
    the last, back down to the baseline and along it to the start (a Circular wrap kept)."""

    positions, radii, baseline_radius, segment_starts = line
    drawn = circular_drawn_line(positions, radii, float(record_len), segment_starts)
    first = position_angle(positions[0], record_len)
    last = position_angle(positions[-1], record_len)
    return (
        f"M {path_point(*circle_point(baseline_radius, first))}"
        + track_line_to(*drawn)
        + f"L {path_point(*circle_point(baseline_radius, last))} "
        + circle_arc_to(baseline_radius, last, first - last)
        + "z"
    )


def circular_series_path_desc(samples: SeriesSamples, place: CircularTrackPlace) -> str:
    """The path of a numeric series track: the area between its line and its baseline. A deviation series of a
    record alone closes its line round the circle and draws the baseline circle as its own subpath (filled even-odd);
    every other line closes down to its baseline (``_area_desc``). A ``line`` series is its line alone, open."""

    line = circular_series_line(samples, place)
    if line is None or len(line.positions) == 0:
        return ""
    if samples.style == "line":
        return track_polyline(*circular_drawn_line(line.positions, line.radii, float(place.record_len), line.segment_starts))
    if samples.deviation is not None and place.record_transform.circle is None:
        return _closed_line_desc(line, float(place.record_len)) + generate_circle_path_desc(place.radius, place.norm_factor)
    return _area_desc(line, place.record_len)


def draw_circle_path(radius: float, stroke_color: str, stroke_width: float) -> Circle:
    """
    Draws a circle path for the circular canvas.
    """
    circle_path = Circle(
        center=(0, 0),
        r=radius,
        stroke=stroke_color,
        stroke_width=stroke_width,
        fill="none",
    )
    return circle_path


__all__ = [
    "CircularAreaLine",
    "CircularTrackPlace",
    "circle_unit_arc",
    "circle_unit_turns",
    "circular_series_baseline_radius",
    "circular_series_line",
    "circular_series_path_desc",
    "draw_circle_path",
    "generate_circle_path_desc",
    "circle_record_turns",
    "within_record_arc",
]


