"""What a Circular track paints, as radial spans along the circle, and whether a box meets it.

A track's owner describes its paint: a ring between two radii over an arc (a feature block, the axis circle, a
conservation hit, an annotation band), the area between a circle and a line drawn through samples (GC content, GC
skew, depth), or a convex polygon (a tick mark, an axis tick label, an annotation label). ``box_meets_ink``
casts rays from the centre through a box and asks whether any ray meets the box and the paint at the same radius.
The centre definition text warns only when one of its line boxes meets a track's paint (D-18).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .circular import RadialBand
from .spatial import Aabb

# Rays through a box are this far apart at most, in turns of the circle (0.25 degrees).
_RAY_STEP_TURNS = 0.25 / 360.0


@dataclass(frozen=True, eq=False)
class CircularInkRing:
    """Paint between two radii over the arc from ``start_turn`` to ``end_turn`` (display turns from the top,
    clockwise; the whole circle by default)."""

    band: RadialBand
    start_turn: float = 0.0
    end_turn: float = 1.0

    def spans(self, turns: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        if self.end_turn - self.start_turn >= 1.0:
            inside = np.ones(turns.shape, dtype=bool)
        elif self.start_turn <= self.end_turn:
            inside = (turns >= self.start_turn) & (turns <= self.end_turn)
        else:
            inside = (turns >= self.start_turn) | (turns <= self.end_turn)
        inner = np.where(inside, float(self.band.inner_px), np.inf)
        outer = np.where(inside, float(self.band.outer_px), -np.inf)
        return inner, outer

    @property
    def nearest_px(self) -> float:
        return float(self.band.inner_px)

    def ray_turns(self) -> NDArray[np.float64]:
        if self.end_turn - self.start_turn >= 1.0:
            return np.empty(0)
        return np.array([self.start_turn, self.end_turn], dtype=float)


@dataclass(frozen=True, eq=False)
class CircularInkProfile:
    """The area between the circle at ``base_px`` and a line through samples: at ``turns[i]`` the paint spans
    ``inner_px[i]`` to ``outer_px[i]``, straight between samples."""

    turns: NDArray[np.float64]
    inner_px: NDArray[np.float64]
    outer_px: NDArray[np.float64]

    @classmethod
    def between(
        cls,
        positions: ArrayLike,
        radii: ArrayLike,
        *,
        record_length: int,
        base_px: float,
        pad_px: float = 0.0,
    ) -> CircularInkProfile:
        """The paint between the circle at ``base_px`` and the line through ``radii`` at display ``positions``,
        widened by ``pad_px`` (half a stroke)."""

        line = np.asarray(radii, dtype=float)
        turns = np.asarray(positions, dtype=float) / float(record_length)
        return cls(
            turns=np.mod(turns, 1.0),
            inner_px=np.minimum(line, base_px) - pad_px,
            outer_px=np.maximum(line, base_px) + pad_px,
        )

    def spans(self, turns: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        if self.turns.size == 0:
            return np.full(turns.shape, np.inf), np.full(turns.shape, -np.inf)
        inner = np.interp(turns, self.turns, self.inner_px, period=1.0)
        outer = np.interp(turns, self.turns, self.outer_px, period=1.0)
        return inner, outer

    @property
    def nearest_px(self) -> float:
        return float(self.inner_px.min()) if self.inner_px.size else math.inf

    def ray_turns(self) -> NDArray[np.float64]:
        return self.turns


@dataclass(frozen=True, eq=False)
class CircularInkPolygon:
    """Paint inside a convex polygon, in coordinates around the circle centre (y down, as drawn)."""

    points: tuple[tuple[float, float], ...]

    @classmethod
    def segment(cls, start: tuple[float, float], end: tuple[float, float], half_width: float) -> CircularInkPolygon:
        """A straight stroke from ``start`` to ``end`` with butt caps, ``half_width`` to each side."""

        (x1, y1), (x2, y2) = start, end
        length = math.hypot(x2 - x1, y2 - y1)
        if length <= 0.0:
            nx, ny = half_width, 0.0
        else:
            nx, ny = -(y2 - y1) / length * half_width, (x2 - x1) / length * half_width
        return cls(((x1 + nx, y1 + ny), (x2 + nx, y2 + ny), (x2 - nx, y2 - ny), (x1 - nx, y1 - ny)))

    @classmethod
    def box(cls, box: Aabb) -> CircularInkPolygon:
        return cls(((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y)))

    def spans(self, turns: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        angles = 2.0 * math.pi * turns - 0.5 * math.pi
        dx, dy = np.cos(angles), np.sin(angles)
        points = np.asarray(self.points, dtype=float)
        centroid = points.mean(axis=0)
        entry = np.zeros(turns.shape)
        exit_ = np.full(turns.shape, np.inf)
        for index in range(len(points)):
            p0, p1 = points[index], points[(index + 1) % len(points)]
            normal = np.array([p1[1] - p0[1], p0[0] - p1[0]])
            if float(np.dot(normal, centroid - p0)) > 0.0:
                normal = -normal
            # Inside: normal . (t * d - p0) <= 0, so t * (normal . d) <= normal . p0.
            along = normal[0] * dx + normal[1] * dy
            limit = float(np.dot(normal, p0))
            safe = np.where(np.abs(along) < 1e-12, 1.0, along)
            bound = limit / safe
            exit_ = np.where(along > 1e-12, np.minimum(exit_, bound), exit_)
            entry = np.where(along < -1e-12, np.maximum(entry, bound), entry)
            outside = (np.abs(along) <= 1e-12) & (limit < 0.0)
            exit_ = np.where(outside, -np.inf, exit_)
        missed = exit_ < entry
        return np.where(missed, np.inf, entry), np.where(missed, -np.inf, exit_)

    @property
    def nearest_px(self) -> float:
        edges = list(zip(self.points, (*self.points[1:], self.points[0])))
        sides = {math.copysign(1.0, x1 * y2 - x2 * y1) for (x1, y1), (x2, y2) in edges}
        return 0.0 if len(sides) == 1 else min(_segment_distance(p0, p1) for p0, p1 in edges)

    def ray_turns(self) -> NDArray[np.float64]:
        points = np.asarray(self.points, dtype=float)
        return _turns_of(np.arctan2(points[:, 1], points[:, 0]))


CircularInk: TypeAlias = CircularInkRing | CircularInkProfile | CircularInkPolygon


def _segment_distance(start: tuple[float, float], end: tuple[float, float]) -> float:
    """How near the centre the segment from ``start`` to ``end`` comes."""

    (x1, y1), (x2, y2) = start, end
    dx, dy = x2 - x1, y2 - y1
    length = dx * dx + dy * dy
    t = 0.0 if length <= 0.0 else min(1.0, max(0.0, -(x1 * dx + y1 * dy) / length))
    return math.hypot(x1 + t * dx, y1 + t * dy)


def paint_shows(color: str, opacity: float = 1.0) -> bool:
    """Whether a fill or stroke of ``color`` at ``opacity`` paints at all."""

    return str(color).strip().lower() != "none" and float(opacity) > 0


def stroke_pad_px(color: str, width: float) -> float:
    """How far a stroke paints past its line: half its width, when it paints at all."""

    return 0.5 * float(width) if float(width) > 0 and paint_shows(color) else 0.0


def _turns_of(angles: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.mod((angles + 0.5 * math.pi) / (2.0 * math.pi), 1.0)


def _box_turn_range(box: Aabb) -> tuple[float, float] | None:
    """The arc of turns a box covers seen from the centre, or ``None`` when the box holds the centre."""

    if box.min_x <= 0.0 <= box.max_x and box.min_y <= 0.0 <= box.max_y:
        return None
    corners = np.array(
        [
            math.atan2(y, x)
            for x, y in ((box.min_x, box.min_y), (box.min_x, box.max_y), (box.max_x, box.min_y), (box.max_x, box.max_y))
        ]
    )
    middle = math.atan2(0.5 * (box.min_y + box.max_y), 0.5 * (box.min_x + box.max_x))
    relative = np.mod(corners - middle + math.pi, 2.0 * math.pi) - math.pi
    start = float(_turns_of(np.array([middle + float(relative.min())]))[0])
    return start, start + float(relative.max() - relative.min()) / (2.0 * math.pi)


def _ray_turns(box: Aabb, inks: Sequence[CircularInk]) -> NDArray[np.float64]:
    covered = _box_turn_range(box)
    start, end = (0.0, 1.0) if covered is None else covered
    count = max(2, int(math.ceil((end - start) / _RAY_STEP_TURNS)) + 1)
    candidates = [np.linspace(start, end, count)]
    # Where each edge comes nearest the centre, and where each ink starts, ends, or turns.
    candidates.append(_turns_of(np.array([0.0, 0.5 * math.pi, math.pi, -0.5 * math.pi])))
    candidates.extend(ink.ray_turns() for ink in inks)
    turns = np.mod(np.concatenate(candidates), 1.0)
    if covered is None:
        return turns
    offset = np.mod(turns - start, 1.0)
    return turns[offset <= (end - start) + 1e-12]


def _ray_box(box: Aabb, turns: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Where each ray from the centre enters and leaves ``box`` (``inf``/``-inf`` when it misses)."""

    angles = 2.0 * math.pi * turns - 0.5 * math.pi
    entries: list[NDArray[np.float64]] = []
    exits: list[NDArray[np.float64]] = []
    for direction, low, high in ((np.cos(angles), box.min_x, box.max_x), (np.sin(angles), box.min_y, box.max_y)):
        parallel = np.abs(direction) < 1e-12
        safe = np.where(parallel, 1.0, direction)
        near = np.minimum(low / safe, high / safe)
        far = np.maximum(low / safe, high / safe)
        holds = low <= 0.0 <= high
        entries.append(np.where(parallel, -np.inf if holds else np.inf, near))
        exits.append(np.where(parallel, np.inf if holds else -np.inf, far))
    entry = np.maximum(np.maximum(entries[0], entries[1]), 0.0)
    exit_ = np.minimum(exits[0], exits[1])
    missed = exit_ < entry
    return np.where(missed, np.inf, entry), np.where(missed, -np.inf, exit_)


def box_meets_ink(box: Aabb, inks: Sequence[CircularInk], *, threshold_px: float) -> bool:
    """Whether ``box`` (centred on the circle) reaches more than ``threshold_px`` into any of ``inks``."""

    inner = Aabb(
        box.min_x + threshold_px, box.min_y + threshold_px, box.max_x - threshold_px, box.max_y - threshold_px
    ) if box.max_x - box.min_x > 2 * threshold_px and box.max_y - box.min_y > 2 * threshold_px else None
    if inner is None:
        return False
    # Only paint that comes nearer the centre than the box's farthest corner can meet it.
    reach = max(math.hypot(x, y) for x in (inner.min_x, inner.max_x) for y in (inner.min_y, inner.max_y))
    inks = [ink for ink in inks if ink.nearest_px <= reach]
    if not inks:
        return False
    turns = _ray_turns(inner, inks)
    entry, exit_ = _ray_box(inner, turns)
    for ink in inks:
        ink_inner, ink_outer = ink.spans(turns)
        if bool(np.any((exit_ >= ink_inner) & (entry <= ink_outer) & (ink_inner <= ink_outer))):
            return True
    return False


__all__ = [
    "CircularInk",
    "CircularInkPolygon",
    "CircularInkProfile",
    "CircularInkRing",
    "box_meets_ink",
    "paint_shows",
    "stroke_pad_px",
]
