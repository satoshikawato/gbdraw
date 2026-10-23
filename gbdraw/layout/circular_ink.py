"""What a Circular track paints, as radial spans along the circle, and whether a box meets it.

A track's owner describes its paint: a ring between two radii over an arc (a feature block, the axis circle, a
band), or the area between a circle and a line drawn through samples (GC content, GC skew). ``box_meets_ink``
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

    def ray_turns(self) -> NDArray[np.float64]:
        return self.turns


CircularInk: TypeAlias = CircularInkRing | CircularInkProfile


def stroke_pad_px(color: str, width: float) -> float:
    """How far a stroke paints past its line: half its width, when it paints at all."""

    return 0.5 * float(width) if float(width) > 0 and str(color).strip().lower() != "none" else 0.0


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
    if inner is None or not inks:
        return False
    turns = _ray_turns(inner, inks)
    entry, exit_ = _ray_box(inner, turns)
    for ink in inks:
        ink_inner, ink_outer = ink.spans(turns)
        if bool(np.any((exit_ >= ink_inner) & (entry <= ink_outer) & (ink_inner <= ink_outer))):
            return True
    return False


__all__ = ["CircularInk", "CircularInkProfile", "CircularInkRing", "box_meets_ink", "stroke_pad_px"]
