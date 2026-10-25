"""The one writer of arcs on a circle around the origin in SVG path data.

An SVG arc is drawn from its endpoints, so its centre comes from the rounded
endpoint numbers. Near 180 degrees (and near 360) a small endpoint error moves
the centre far, and the drawn arc bows off its circle (OV-1840). Every circle
arc gbdraw draws is therefore split here, so that no drawn arc is within 30
degrees of 180 or near 360:

- a span of at most ``ARC_SINGLE_MAX_DEG`` is one arc;
- a span below ``ARC_HALVES_MAX_DEG`` is two equal halves (each at most 150 degrees);
- a longer span is ``span / 2 - ARC_UNEQUAL_SHIFT_DEG`` then ``span / 2 + ARC_UNEQUAL_SHIFT_DEG``
  (a full circle is 120 + 240 degrees).

Angles are in degrees in the SVG frame: 0 at 3 o'clock, positive clockwise.
"""

from __future__ import annotations

import math

from .path_data import path_number, path_point

ARC_SINGLE_MAX_DEG = 20.0
ARC_HALVES_MAX_DEG = 300.0
ARC_UNEQUAL_SHIFT_DEG = 60.0


def circle_point(radius: float, angle_deg: float) -> tuple[float, float]:
    """The point at ``angle_deg`` on the circle of ``radius`` around the origin."""

    angle = math.radians(angle_deg)
    return float(radius) * math.cos(angle), float(radius) * math.sin(angle)


def arc_pieces(span_deg: float) -> tuple[float, ...]:
    """The signed spans of the arcs that draw ``span_deg`` (see the module rule)."""

    size = abs(float(span_deg))
    if size <= ARC_SINGLE_MAX_DEG:
        return (float(span_deg),)
    sign = 1.0 if span_deg > 0 else -1.0
    if size < ARC_HALVES_MAX_DEG:
        return (span_deg / 2.0, span_deg / 2.0)
    return (sign * (size / 2.0 - ARC_UNEQUAL_SHIFT_DEG), sign * (size / 2.0 + ARC_UNEQUAL_SHIFT_DEG))


def circle_arc_to(radius: float, start_deg: float, span_deg: float) -> str:
    """The ``A`` commands that continue a path standing at ``start_deg`` on the circle of ``radius`` along
    ``span_deg`` (positive clockwise), each with its own large-arc and sweep flags."""

    r = path_number(radius)
    commands = []
    angle = float(start_deg)
    for piece in arc_pieces(span_deg):
        angle += piece
        large = 1 if abs(piece) > 180.0 else 0
        sweep = 1 if piece > 0 else 0
        commands.append(f"A {r},{r} 0 {large} {sweep} {path_point(*circle_point(radius, angle))}")
    return " ".join(commands)


def circle_arc(radius: float, start_deg: float, span_deg: float) -> str:
    """``M`` at ``start_deg`` on the circle of ``radius``, then :func:`circle_arc_to`."""

    return f"M {path_point(*circle_point(radius, start_deg))} {circle_arc_to(radius, start_deg, span_deg)}"


def position_angle(position: float, total_length: float) -> float:
    """The angle of a sequence position on a Circular record (0 at 12 o'clock, clockwise)."""

    return 360.0 * (float(position) / float(total_length)) - 90.0


__all__ = [
    "ARC_HALVES_MAX_DEG",
    "ARC_SINGLE_MAX_DEG",
    "ARC_UNEQUAL_SHIFT_DEG",
    "arc_pieces",
    "circle_arc",
    "circle_arc_to",
    "circle_point",
    "position_angle",
]
