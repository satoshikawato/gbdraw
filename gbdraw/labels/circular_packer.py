"""An order-preserving packer of one side's horizontal Circular labels (CP-D04, CP-D05).

It is tried only for a record whose placement does not fit, and its placement is
drawn only when the fit check then finds neither overlaps nor crossings
(``prepare_label_list``). The labels keep their cyclic genome order on the side's
label ellipse (:class:`~gbdraw.labels.circular_capacity.LabelArc`), parameterised by

    psi: anchor = (center_x + x_radius sin psi, -y_radius cos psi),

0 at 12 o'clock, clockwise. Each label's target is where its feature's ray meets
the ellipse. Neighbouring boxes (``_label_y_offsets`` rows and the ``start_x > 0``
text side of ``_label_x_offsets``) stay a row gap apart vertically or 1 px apart
side by side. The packing minimises the squared distance of each label from its
target: pool-adjacent-violators on ``psi - D`` with the neighbour steps ``D``
measured again a few times, then bounded sweeps.

Two looks, named by what the reader sees, tried in this order (CP-D05):

- ``"half"`` (look 2): every label stays on its feature's genome half, as today
  (no label passes 12 or 6 o'clock).
- ``"pole"`` (look 1), only when no ``"half"`` placement fits: first whole runs at
  the ends of a genome half pass 12 or 6 o'clock to the other side (the splits
  of least displacement first); when no such split packs under the crossing
  limit, one cycle of labels, each moving only inside the cone of angles at
  which its straight leader from the elbow still runs outward.

A packing's crossing leaders are counted with the fit check's own test
(:func:`~gbdraw.labels.circular_leaders.leader_crossing_count`).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

from .circular import LEGACY_MIN_ROW_GAP_PX, _assign_leader_start_points, _label_x_offsets, _label_y_offsets
from .circular_capacity import LabelArc
from .circular_leaders import leader_crossing_count
from .circular_types import CircularLabel

PackerLook: TypeAlias = Literal["half", "pole"]
# CP-D05: look 2 (every label on its genome half) first, look 1 (runs pass a pole) only when look 2 does not fit.
PACKER_LOOKS: tuple[PackerLook, ...] = ("half", "pole")

# Side by side, two boxes stay this far apart (LP-D37).
SIDE_GAP_PX = 1.0
# "pole": the cone's angular margin, and its widenings tried in order.
CONE_MARGIN_DEG = 1.0
CONE_SCALES = (1.0, 1.25, 1.5, 2.0, 3.0, 5.0)
_CONE_STEP_RAD = math.radians(0.25)
# "pole": the end runs that may pass a pole hold at most this many labels.
MAX_RUN_SHIFT = 20
_STEP_START_RAD = 0.002
_STEP_BISECTIONS = 18
_PACK_ITERATIONS = 4
_CLOSING_SWEEPS = 3
_TWO_PI = 2.0 * math.pi

Box: TypeAlias = tuple[float, float, float, float]
Span: TypeAlias = tuple[float, float]


@dataclass(frozen=True)
class _PackItem:
    """One label with its box offsets from the anchor, its target and its elbow."""

    label: CircularLabel
    rows: Span  # (top, bottom) of the box below the anchor's y
    right_text: Span  # (left, right) of the box from the anchor's x, for an anchor right of the axis
    left_text: Span  # the same for an anchor left of the axis
    target: float  # psi where the feature's ray meets the ellipse, in (-pi, pi]
    elbow_radius: float
    elbow_angle: float


class _Ring:
    """The label ellipse of one side with the gaps two neighbouring boxes keep."""

    def __init__(self, arc: LabelArc, row_gap: float) -> None:
        self.arc = arc
        self.row_gap = row_gap

    def anchor(self, psi: float) -> tuple[float, float]:
        return self.arc.center_x + self.arc.x_radius * math.sin(psi), -self.arc.y_radius * math.cos(psi)

    def box(self, item: _PackItem, psi: float) -> Box:
        x, y = self.anchor(psi)
        left, right = item.right_text if x > 0 else item.left_text
        return x + left, y + item.rows[0], x + right, y + item.rows[1]

    def boxes_apart(self, a: Box, b: Box) -> bool:
        return (
            b[1] - a[3] >= self.row_gap
            or a[1] - b[3] >= self.row_gap
            or b[0] - a[2] >= SIDE_GAP_PX
            or a[0] - b[2] >= SIDE_GAP_PX
        )

    def apart(self, first: _PackItem, first_psi: float, second: _PackItem, second_psi: float) -> bool:
        return self.boxes_apart(self.box(first, first_psi), self.box(second, second_psi))

    def step(self, fixed: _PackItem, fixed_psi: float, moving: _PackItem, *, forward: bool = True) -> float:
        """The least angle ``moving`` keeps from ``fixed`` (after it when ``forward``, else before it)."""

        fixed_box = self.box(fixed, fixed_psi)

        def clear(delta: float) -> bool:
            if forward:
                return self.boxes_apart(fixed_box, self.box(moving, fixed_psi + delta))
            return self.boxes_apart(self.box(moving, fixed_psi - delta), fixed_box)

        if clear(0.0):
            return 0.0
        high = _STEP_START_RAD
        while not clear(high):
            high *= 2.0
            if high > math.pi:
                return math.pi
        low = 0.0
        for _ in range(_STEP_BISECTIONS):
            middle = 0.5 * (low + high)
            if clear(middle):
                high = middle
            else:
                low = middle
        return high

    def sweep_forward(self, items: Sequence[_PackItem], psis: Sequence[float], lower: Sequence[float]) -> list[float]:
        out = list(psis)
        for k in range(len(items)):
            out[k] = max(out[k], lower[k])
            if k:
                out[k] = max(out[k], out[k - 1] + self.step(items[k - 1], out[k - 1], items[k]))
        return out

    def sweep_backward(self, items: Sequence[_PackItem], psis: Sequence[float], upper: Sequence[float]) -> list[float]:
        out = list(psis)
        for k in range(len(items) - 1, -1, -1):
            out[k] = min(out[k], upper[k])
            if k < len(items) - 1:
                out[k] = min(out[k], out[k + 1] - self.step(items[k + 1], out[k + 1], items[k], forward=False))
        return out

    def pack(
        self,
        items: Sequence[_PackItem],
        targets: Sequence[float],
        lower: Sequence[float],
        upper: Sequence[float],
    ) -> list[float] | None:
        """Angles in order within ``lower``/``upper`` closest to ``targets``, or None when they do not fit."""

        count = len(items)
        upper = list(upper)
        psis = [min(max(target, lower[k]), upper[k]) for k, target in enumerate(targets)]
        for _ in range(_PACK_ITERATIONS):
            steps = [0.0]
            for k in range(1, count):
                steps.append(steps[-1] + self.step(items[k - 1], psis[k - 1], items[k]))
            floor: list[float] = []
            running = -math.inf
            for k in range(count):
                running = max(running, lower[k] - steps[k])
                floor.append(running)
            ceiling = [0.0] * count
            running = math.inf
            for k in range(count - 1, -1, -1):
                running = min(running, upper[k] - steps[k])
                ceiling[k] = running
            fitted = _isotonic([targets[k] - steps[k] for k in range(count)])
            psis = [max(min(fitted[k], ceiling[k]), floor[k]) + steps[k] for k in range(count)]
        for _ in range(_CLOSING_SWEEPS):
            # The cycle closes: the last label stays clear of the first one, one turn later.
            if count > 1:
                upper[-1] = min(
                    upper[-1],
                    psis[0] + _TWO_PI - self.step(items[0], psis[0] + _TWO_PI, items[-1], forward=False),
                )
            psis = self.sweep_forward(items, psis, lower)
            psis = self.sweep_backward(items, psis, upper)
        psis = self.sweep_forward(items, psis, lower)
        if any(psis[k] > upper[k] + 1e-9 for k in range(count)):
            return None
        if count > 1 and not self.apart(items[-1], psis[-1], items[0], psis[0] + _TWO_PI):
            return None
        return psis


def _isotonic(values: Sequence[float]) -> list[float]:
    """The least-squares non-decreasing fit of ``values`` (pool adjacent violators)."""

    sums: list[float] = []
    counts: list[int] = []
    for value in values:
        sums.append(value)
        counts.append(1)
        while len(sums) > 1 and sums[-2] / counts[-2] > sums[-1] / counts[-1]:
            total, count = sums.pop(), counts.pop()
            sums[-1] += total
            counts[-1] += count
    out: list[float] = []
    for total, count in zip(sums, counts):
        out.extend([total / count] * count)
    return out


def _item(label: CircularLabel, arc: LabelArc, total_length: int, is_outer: bool) -> _PackItem:
    top, _, bottom, _ = _label_y_offsets(float(label["middle"]), total_length, float(label["height_px"]), 0.0, is_outer)
    right_text, left_text = _label_x_offsets(float(label["width_px"]), 0.0, is_outer)
    # The feature's ray from the record's centre meets the ellipse at distance `reach`.
    feature_x, feature_y = float(label["feature_middle_x"]), float(label["feature_middle_y"])
    norm = math.hypot(feature_x, feature_y) or 1.0
    ux, uy = feature_x / norm, feature_y / norm
    a, b, cx = arc.x_radius, arc.y_radius, arc.center_x
    qa = ux * ux / (a * a) + uy * uy / (b * b)
    qb = -2.0 * ux * cx / (a * a)
    qc = cx * cx / (a * a) - 1.0
    reach = (-qb + math.sqrt(max(0.0, qb * qb - 4.0 * qa * qc))) / (2.0 * qa)
    elbow_x, elbow_y = float(label["middle_x"]), float(label["middle_y"])
    return _PackItem(
        label=label,
        rows=(top, bottom),
        right_text=(right_text[0], right_text[2]),
        left_text=(left_text[0], left_text[2]),
        target=math.atan2((reach * ux - cx) / a, -reach * uy / b),
        elbow_radius=math.hypot(elbow_x, elbow_y),
        elbow_angle=math.atan2(elbow_y, elbow_x),
    )


def _placed(ring: _Ring, items: Sequence[_PackItem], psis: Sequence[float], total_length: int) -> list[CircularLabel]:
    """Copies of the labels at these angles, in genome order, with their leader starts."""

    placed: list[CircularLabel] = []
    for item, psi in zip(items, psis):
        label = item.label.copy()
        label["start_x"], label["start_y"] = ring.anchor(psi)
        placed.append(label)
    placed.sort(key=lambda label: float(label["middle"]))
    return _assign_leader_start_points(placed, total_length)


def _side_limits(arc: LabelArc) -> tuple[Span, Span]:
    """The angles of the right and the left side, whose anchors keep the side gap from the vertical axis."""

    def angle(x: float) -> float:
        return math.asin(max(-1.0, min(1.0, (x - arc.center_x) / arc.x_radius)))

    right = angle(SIDE_GAP_PX)
    left = angle(-SIDE_GAP_PX)
    return (right, math.pi - right), (math.pi - left, _TWO_PI + left)


class _HalfLook:
    """Two sides split at the genome halves; end runs pass a pole only as wholes (shifts ``top``, ``bottom``).

    The split (0, 0) is look 2; any other split is look 1's first family.
    """

    def __init__(self, items: Sequence[_PackItem], arc: LabelArc, total_length: int) -> None:
        self.items = items
        self.count = len(items)
        self.half = next(
            (k for k, item in enumerate(items) if float(item.label["middle"]) >= total_length / 2), self.count
        )
        self.sides = _side_limits(arc)

    def valid(self, shift: int) -> bool:
        return -(self.count - self.half) < shift < self.half

    def split(self, top: int, bottom: int) -> tuple[list[int], list[int]]:
        """Item indices of the right side (top to bottom) and the left side (bottom to top).

        ``top`` > 0 moves the first ``top`` labels to the left side over 12 o'clock, < 0 the last ones to the
        right side; ``bottom`` > 0 moves the last ``bottom`` labels of the first half to the left side over
        6 o'clock, < 0 the first ones of the second half to the right side.
        """

        count, half = self.count, self.half
        right = (
            (list(range(count + top, count)) if top < 0 else [])
            + list(range(max(top, 0), half - max(bottom, 0)))
            + (list(range(half, half - bottom)) if bottom < 0 else [])
        )
        left = (
            (list(range(half - bottom, half)) if bottom > 0 else [])
            + list(range(half + max(-bottom, 0), count - max(-top, 0)))
            + (list(range(0, top)) if top > 0 else [])
        )
        return right, left

    def side(self, indices: Sequence[int], limits: Span) -> tuple[list[_PackItem], list[float]]:
        """The side's items and their targets within its limits, made non-decreasing."""

        low, high = limits
        middle = 0.5 * (low + high)
        targets: list[float] = []
        for k in indices:
            target = self.items[k].target % _TWO_PI
            target += _TWO_PI * round((middle - target) / _TWO_PI)
            target = min(max(target, low), high)
            targets.append(max(target, targets[-1]) if targets else target)
        return [self.items[k] for k in indices], targets

    def excess(self, ring: _Ring, top: int, bottom: int) -> float:
        """How far the fullest side overflows when its labels stack from its start (<= 0: both sides fit)."""

        worst = -1.0
        for indices, limits in zip(self.split(top, bottom), self.sides):
            if indices:
                seq, _ = self.side(indices, limits)
                start = [limits[0]] * len(seq)
                worst = max(worst, ring.sweep_forward(seq, start, start)[-1] - limits[1])
        return worst

    def greedy(self, ring: _Ring) -> tuple[int, int] | None:
        """A split that fits both sides, by single steps from (0, 0) that lower the excess; None when stuck."""

        state, current = (0, 0), self.excess(ring, 0, 0)
        while current > 0:
            moves = [
                (state[0] + d_top, state[1] + d_bottom)
                for d_top, d_bottom in ((1, 0), (-1, 0), (0, 1), (0, -1))
                if self.valid(state[0] + d_top) and self.valid(state[1] + d_bottom)
            ]
            best = min(((self.excess(ring, *move), move) for move in moves), default=None)
            if best is None or best[0] >= current:
                return None
            current, state = best
        return state

    def _prefix_excess(self, ring: _Ring, indices: Sequence[int], limits: Span) -> list[float]:
        """The overflow of every prefix of ``indices`` stacked from the side's start (``excess`` of each prefix)."""

        if not indices:
            return []
        seq, _ = self.side(indices, limits)
        start = [limits[0]] * len(seq)
        return [psi - limits[1] for psi in ring.sweep_forward(seq, start, start)]

    def states(self, ring: _Ring) -> list[tuple[int, int]]:
        """Per top shift up to ``MAX_RUN_SHIFT``, the smallest bottom shift that fits; fewest moved first.

        The values equal ``excess`` per state: for one top shift the right sides of all bottom shifts are
        prefixes of one sequence, and for one bottom shift the left sides of all top shifts are, so one
        sweep per shift replaces one per state.
        """

        shifts = [shift for shift in range(-MAX_RUN_SHIFT, MAX_RUN_SHIFT + 1) if self.valid(shift)]
        if not shifts:
            return []
        right_limits, left_limits = self.sides
        right = {top: self._prefix_excess(ring, self.split(top, shifts[0])[0], right_limits) for top in shifts}
        left = {bottom: self._prefix_excess(ring, self.split(shifts[-1], bottom)[1], left_limits) for bottom in shifts}

        def fits(top: int, bottom: int) -> bool:
            right_side, left_side = self.split(top, bottom)
            return (not right_side or right[top][len(right_side) - 1] <= 0) and (
                not left_side or left[bottom][len(left_side) - 1] <= 0
            )

        found: list[tuple[int, int]] = []
        for top in shifts:
            bottom = next((bottom for bottom in sorted(shifts, key=abs) if fits(top, bottom)), None)
            if bottom is not None:
                found.append((top, bottom))
        found.sort(key=lambda state: abs(state[0]) + abs(state[1]))
        return found

    def packed(self, ring: _Ring, state: tuple[int, int]) -> tuple[list[_PackItem], list[float], list[float]] | None:
        seq_all: list[_PackItem] = []
        psis_all: list[float] = []
        targets_all: list[float] = []
        for indices, limits in zip(self.split(*state), self.sides):
            if not indices:
                continue
            seq, targets = self.side(indices, limits)
            psis = ring.pack(seq, targets, [limits[0]] * len(seq), [limits[1]] * len(seq))
            if psis is None:
                return None
            seq_all += seq
            psis_all += psis
            targets_all += targets
        return seq_all, psis_all, targets_all


def _fewest_crossings(placements: Iterable[list[CircularLabel]], max_crossings: int) -> list[CircularLabel] | None:
    """The first placement with the fewest crossing leader pairs, at most ``max_crossings``; one without ends it."""

    best: list[CircularLabel] | None = None
    bound = max_crossings + 1
    for placed in placements:
        count = leader_crossing_count(placed, limit=bound)
        if count < bound:
            best, bound = placed, count
            if count == 0:
                break
    return best


_ON_HALVES = (0, 0)


def _half_placements(
    items: Sequence[_PackItem], arc: LabelArc, row_gaps: Sequence[float], total_length: int
) -> Iterator[list[CircularLabel]]:
    """Look 2: per row gap, the split at the genome halves (no label passes a pole), when it packs."""

    look = _HalfLook(items, arc, total_length)
    for row_gap in row_gaps:
        ring = _Ring(arc, row_gap)
        packing = look.packed(ring, _ON_HALVES)
        if packing is not None:
            yield _placed(ring, packing[0], packing[1], total_length)


def _run_placements(
    items: Sequence[_PackItem], arc: LabelArc, row_gaps: Sequence[float], total_length: int
) -> Iterator[list[CircularLabel]]:
    """Look 1, first family: per row gap, the splits whose end runs pass a pole, by least squared displacement."""

    look = _HalfLook(items, arc, total_length)
    for row_gap in row_gaps:
        ring = _Ring(arc, row_gap)
        fallback = look.greedy(ring)
        if fallback is None:
            continue
        states = [state for state in look.states(ring) if state != _ON_HALVES]
        if not states and fallback != _ON_HALVES:
            states = [fallback]
        packings = [packing for state in states if (packing := look.packed(ring, state)) is not None]
        packings.sort(key=lambda packing: sum((psi - t) ** 2 for psi, t in zip(packing[1], packing[2])))
        for seq, psis, _targets in packings:
            yield _placed(ring, seq, psis, total_length)


def _cone(arc: LabelArc, item: _PackItem, target: float, scale: float) -> Span:
    """The angles around ``target`` at which the label's straight leader from its elbow still runs outward."""

    def outward(psi: float) -> bool:
        x = arc.center_x + arc.x_radius * math.sin(psi)
        y = -arc.y_radius * math.cos(psi)
        radius = math.hypot(x, y)
        if radius <= item.elbow_radius:
            return False
        turn = (math.atan2(y, x) - item.elbow_angle + math.pi) % _TWO_PI - math.pi
        reach = math.acos(min(1.0, item.elbow_radius / radius)) - math.radians(CONE_MARGIN_DEG)
        return abs(turn) <= scale * max(0.0, reach)

    if not outward(target):
        return target, target
    low = high = target
    while outward(low - _CONE_STEP_RAD) and low - _CONE_STEP_RAD > target - math.pi:
        low -= _CONE_STEP_RAD
    while outward(high + _CONE_STEP_RAD) and high + _CONE_STEP_RAD < target + math.pi:
        high += _CONE_STEP_RAD
    return low, high


def _cycle_placements(
    items: Sequence[_PackItem], arc: LabelArc, row_gaps: Sequence[float], total_length: int
) -> Iterator[list[CircularLabel]]:
    """Look 1, second family: one cycle cut at the widest gap between targets, per cone widening and row gap."""

    count = len(items)
    raw = [item.target % _TWO_PI for item in items]
    spans = [(raw[(k + 1) % count] - raw[k]) % _TWO_PI for k in range(count)]
    cut = (max(range(count), key=spans.__getitem__) + 1) % count if count > 1 else 0
    order = [*items[cut:], *items[:cut]]
    targets = [raw[(cut + k) % count] for k in range(count)]
    for k in range(1, count):
        while targets[k] < targets[k - 1]:
            targets[k] += _TWO_PI
    for scale in CONE_SCALES:
        cones = [_cone(arc, item, target, scale) for item, target in zip(order, targets)]
        lower = [cone[0] for cone in cones]
        upper = [cone[1] for cone in cones]
        for row_gap in row_gaps:
            ring = _Ring(arc, row_gap)
            psis = ring.pack(order, targets, lower, upper)
            if psis is not None:
                yield _placed(ring, order, psis, total_length)


def pack_side(
    labels: Sequence[CircularLabel],
    arc: LabelArc,
    look: PackerLook,
    *,
    total_length: int,
    is_outer: bool,
    spacing_px: float,
    max_crossings: int,
) -> list[CircularLabel] | None:
    """One side's labels packed in genome order on ``arc`` in ``look``, or None when the look finds no placement.

    The row gaps tried are the circular label spacing, then the 0.1 px of the dense rows (LP-D45 order).
    The result is a family's first packing of fewest crossing leader pairs, at most ``max_crossings``
    (one without crossings ends the search); look 1 tries its cycle family only when its run family
    has none. Returns copies in genome order with ``start_x``/``start_y`` and leader starts; the input
    is not changed.
    """

    ordered = sorted(labels, key=lambda label: float(label["middle"]))
    if not ordered:
        return []
    items = [_item(label, arc, total_length, is_outer) for label in ordered]
    row_gaps = sorted({max(float(spacing_px), LEGACY_MIN_ROW_GAP_PX), LEGACY_MIN_ROW_GAP_PX}, reverse=True)
    families = (_half_placements,) if look == "half" else (_run_placements, _cycle_placements)
    for family in families:
        packed = _fewest_crossings(family(items, arc, row_gaps, total_length), max_crossings)
        if packed is not None:
            return packed
    return None


__all__ = ["PACKER_LOOKS", "PackerLook", "pack_side"]
