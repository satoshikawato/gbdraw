"""Local repair of slight label defects in a horizontal Circular layout (LP-D50, LP-D57).

Before the fit check counts a record, the final layout of its horizontal outer
and inner labels is repaired when it has 1 to ``MAX_DEFECTS`` defects:
overlapping labels, crossing leaders and labels over the centre definition
text, counted with the fit check's own predicates (:mod:`gbdraw.labels.fit`).
The repair moves labels only locally:

- elbow: a leader's bend point moves along its feature's radial line (one
  label, or a run as a ramp so that neighbouring leaders nest);
- toward the ring / away: one label or a run moves radially; toward the ring
  it stops at the label's own elbow circle;
- push: one label moves along its circle and the following labels of its side
  move just far enough to stay clear;
- arc: a run of neighbouring labels shifts along their circles;
- contact: the leader meets the LP-D23 point of the ink box that its arrival
  selects (LP-D52, LP-D55, :func:`~gbdraw.labels.circular_leaders.arrival_contact`).

Labels keep their side and the order of their features (LP-D50); a text moves
at most ``MAX_TEXT_SHIFT_PX`` and an elbow ``MAX_ELBOW_SHIFT_PX``. The repaired
layout is kept only when the record then has no defect at all; otherwise the
layout stays as it was. The search is bounded by ``MAX_EVALUATIONS`` candidate
evaluations, so its result does not depend on the machine, and it stops early
once the last move's rate (evaluations per cleared defect) cannot clear the
remaining defects within that budget.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

from ..layout.spatial import Aabb, touching_aabb_indices
from ..layout.text_geometry import Point
from .circular import (
    FITS,
    HEMISPHERE_AXIS_EPS_PX,
    _current_half_from_x,
    _label_bbox_max_radius,
    _label_bbox_min_radius,
    _label_bbox_rect,
    _label_fit_box,
    _leader_start_point,
    _segment_intersects_rect_interior,
    angle_from_middle_unwrapped,
    assign_leader_start_points,
    circular_label_fits,
    circular_label_shape,
    normalize_angle_near_reference,
)
from .circular_leaders import arrival_contact, contact_end, label_ink_box
from .fit import LabelShape, covers_definition_text, polylines_cross, shapes_overlap

# LP-D57: repair records with up to this many defects, within this many evaluations.
MAX_DEFECTS = 100
MAX_EVALUATIONS = 60_000
# LP-D57 Q3: how far a text and an elbow may move from where the engine put them.
MAX_TEXT_SHIFT_PX = 24.0
MAX_ELBOW_SHIFT_PX = 12.0
# Neighbours keep their order with at least this angular gap (or today's, if smaller).
MIN_ORDER_GAP_DEG = 0.05

ARC_STEPS_DEG = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
RADIAL_STEPS_PX = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0)
RUN_RADIAL_STEPS_PX = (2.0, 4.0, 6.0, 8.0, 10.0, 12.0)
ELBOW_STEPS_PX = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0)
ELBOW_RAMP_PX = (0.5, 1.0, 1.5, 2.0, 3.0)
# Runs reach this many labels beyond a defect cluster and hold at most MAX_RUN labels.
RUN_PAD = 4
MAX_RUN = 12
# A push moves at most this many following labels.
MAX_PUSH = 24
# A one-label fix that clears its cluster within this displacement ends the cluster's search.
CLEAN_STOP_PX = 6.0

_EPS = 1e-9

Side: TypeAlias = Literal["out", "in"]
DefectKind: TypeAlias = Literal["overlap", "crossing", "definition"]
Defect: TypeAlias = tuple[DefectKind, int, int]
Moved: TypeAlias = dict[int, dict]
# (min_x, min_y, max_x, max_y) of the broad phase.
Box: TypeAlias = tuple[float, float, float, float]


@dataclass(frozen=True)
class _Geometry:
    """One label as the counts see it, with the boxes the broad phase compares."""

    shape: LabelShape
    box: Box
    leader_box: Box | None
    reach: Box


def _points_box(points: Sequence[Point]) -> Box:
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    return (min(xs), min(ys), max(xs), max(ys))


def _apart(first: Box, second: Box) -> bool:
    return first[2] < second[0] or second[2] < first[0] or first[3] < second[1] or second[3] < first[1]


def _angle(x: float, y: float) -> float:
    return math.degrees(math.atan2(y, x))


def _clockwise(first: float, second: float) -> float:
    """The angle from ``first`` to ``second`` in the direction of increasing SVG angle, in [0, 360)."""

    return (second - first) % 360.0


class _Repair:
    """The search state of one record: copies of its labels, never the caller's."""

    def __init__(self, labels: Sequence[dict], total_length: int, definition_text: Sequence[Aabb]) -> None:
        self.total_length = int(total_length)
        self.definition_text = tuple(definition_text)
        self.labels = [dict(label) for label in labels]
        self.original = [dict(label) for label in labels]
        self.count = len(self.labels)
        self.movable = [not label.get("is_embedded") for label in self.labels]
        self.external = [index for index in range(self.count) if self.movable[index]]
        self.side: list[Side] = ["in" if label.get("is_inner") else "out" for label in self.labels]
        self.order: dict[Side, list[int]] = {
            side: sorted(
                (index for index in self.external if self.side[index] == side),
                key=lambda index: float(self.labels[index]["middle"]),
            )
            for side in ("out", "in")
        }
        self.position = {index: (side, at) for side, indices in self.order.items() for at, index in enumerate(indices)}
        self.half = {
            index: _current_half_from_x(float(self.labels[index]["start_x"]), axis_eps_px=HEMISPHERE_AXIS_EPS_PX)
            for index in self.external
        }
        self.elbow_radius = {index: self._elbow_radius(self.labels[index]) for index in self.external}
        self.feature_radius = {
            index: math.hypot(
                float(self.labels[index].get("feature_anchor_x", 0.0)),
                float(self.labels[index].get("feature_anchor_y", 0.0)),
            )
            for index in self.external
        }
        self.min_radius = {index: _label_bbox_min_radius(self.labels[index], self.total_length) for index in self.external}
        self.max_radius = {index: _label_bbox_max_radius(self.labels[index], self.total_length) for index in self.external}
        self.angle = {
            index: _angle(float(self.labels[index]["start_x"]), float(self.labels[index]["start_y"]))
            for index in self.external
        }
        self.evaluations = 0
        self._refresh()
        # Pairs that a move within the bounds can bring together.
        tallest = max((float(self.labels[index]["height_px"]) for index in self.external), default=0.0)
        expand = 2.0 * MAX_TEXT_SHIFT_PX + 2.0 * tallest + 4.0
        reaches = [Aabb(*geometry.reach) for geometry in self.geometry]
        self.peers: list[list[int]] = [
            [
                other
                for other in touching_aabb_indices(
                    reaches,
                    Aabb(reach.min_x - expand, reach.min_y - expand, reach.max_x + expand, reach.max_y + expand),
                    skip=index,
                )
                if self.movable[index] or self.movable[other]
            ]
            for index, reach in enumerate(reaches)
        ]

    # ---- geometry -------------------------------------------------------
    @staticmethod
    def _elbow_radius(label: dict) -> float:
        return math.hypot(float(label["middle_x"]), float(label["middle_y"]))

    def _geometry_of(self, label: dict) -> _Geometry:
        shape = circular_label_shape(label, self.total_length)
        box = _points_box([point for part in shape.parts for point in part]) if shape.parts else (0.0, 0.0, 0.0, 0.0)
        if not shape.leader:
            return _Geometry(shape, box, None, box)
        leader_box = _points_box(shape.leader)
        reach = (
            min(box[0], leader_box[0]),
            min(box[1], leader_box[1]),
            max(box[2], leader_box[2]),
            max(box[3], leader_box[3]),
        )
        return _Geometry(shape, box, leader_box, reach)

    def _refresh(self) -> None:
        """Recompute every engine leader end (they depend on the neighbours) and the geometry."""

        assign_leader_start_points([self.labels[index] for index in self.external], self.total_length)
        self.geometry = [self._geometry_of(label) for label in self.labels]
        self.text_rects = {index: _label_bbox_rect(self.labels[index], self.total_length, 0.0) for index in self.external}
        self._local_cache: dict[tuple[int, ...], frozenset[Defect]] = {}
        self._through_cache: dict[tuple[int, ...], int] = {}

    # ---- counts (the fit check's predicates) ------------------------------
    def _pair_defects(self, first: int, first_geometry: _Geometry, second: int, second_geometry: _Geometry) -> list[Defect]:
        if _apart(first_geometry.reach, second_geometry.reach):
            return []
        low, high = (first, second) if first < second else (second, first)
        found: list[Defect] = []
        if (
            first_geometry.shape.parts
            and second_geometry.shape.parts
            and not _apart(first_geometry.box, second_geometry.box)
            and shapes_overlap(first_geometry.shape.parts, second_geometry.shape.parts, touching=True)
        ):
            found.append(("overlap", low, high))
        if (
            first_geometry.leader_box is not None
            and second_geometry.leader_box is not None
            and not _apart(first_geometry.leader_box, second_geometry.leader_box)
            and polylines_cross(first_geometry.shape.leader, second_geometry.shape.leader)
        ):
            found.append(("crossing", low, high))
        return found

    def defects(self) -> set[Defect]:
        found: set[Defect] = set()
        for index in range(self.count):
            for other in self.peers[index]:
                if other > index:
                    found.update(self._pair_defects(index, self.geometry[index], other, self.geometry[other]))
            if covers_definition_text(self.geometry[index].shape, self.definition_text):
                found.add(("definition", index, index))
        return found

    def _local(self, geometry: dict[int, _Geometry]) -> frozenset[Defect]:
        """Defects that involve a label of ``geometry``, with those labels in their candidate places."""

        found: set[Defect] = set()
        for index, own in geometry.items():
            for other in self.peers[index]:
                if other in geometry and other < index:
                    continue
                found.update(self._pair_defects(index, own, other, geometry.get(other, self.geometry[other])))
            if covers_definition_text(own.shape, self.definition_text):
                found.add(("definition", index, index))
        return frozenset(found)

    def _current_local(self, indices: Sequence[int]) -> frozenset[Defect]:
        key = tuple(sorted(indices))
        cached = self._local_cache.get(key)
        if cached is None:
            cached = self._local({index: self.geometry[index] for index in key})
            self._local_cache[key] = cached
        return cached

    def _through(self, moved: Moved, geometry: dict[int, _Geometry]) -> int:
        """Leaders drawn through another label's text, for pairs that involve a moved label."""

        rects = {index: _label_bbox_rect(label, self.total_length, 0.0) for index, label in moved.items()}
        found: set[tuple[int, int]] = set()
        for index in geometry:
            for other in self.peers[index]:
                for leader_of, text_of in ((index, other), (other, index)):
                    if not self.movable[text_of]:
                        continue  # an embedded label's text sits on its feature
                    leader_geometry = geometry.get(leader_of, self.geometry[leader_of])
                    leader = leader_geometry.shape.leader
                    rect = rects.get(text_of, self.text_rects[text_of])
                    if leader_geometry.leader_box is None or _apart(leader_geometry.leader_box, rect):
                        continue
                    if any(
                        _segment_intersects_rect_interior(start, end, rect) for start, end in zip(leader, leader[1:])
                    ):
                        found.add((leader_of, text_of))
        return len(found)

    def _current_through(self, indices: Sequence[int]) -> int:
        key = tuple(sorted(indices))
        cached = self._through_cache.get(key)
        if cached is None:
            cached = self._through({index: self.labels[index] for index in key}, {index: self.geometry[index] for index in key})
            self._through_cache[key] = cached
        return cached

    # ---- rules ------------------------------------------------------------
    def _displacement(self, index: int, label: dict) -> tuple[float, float]:
        original = self.original[index]
        text = math.hypot(
            float(label["start_x"]) - float(original["start_x"]), float(label["start_y"]) - float(original["start_y"])
        )
        elbow = math.hypot(
            float(label["middle_x"]) - float(original["middle_x"]),
            float(label["middle_y"]) - float(original["middle_y"]),
        )
        return text, elbow

    def _order_kept(self, moved: Moved) -> bool:
        """No moved label passes a neighbour of its side; a pair out of order today stays as it is."""

        checked: set[tuple[int, int]] = set()
        for index in moved:
            side, at = self.position[index]
            indices = self.order[side]
            if len(indices) < 2:
                continue
            for first, second in (
                (indices[(at - 1) % len(indices)], index),
                (index, indices[(at + 1) % len(indices)]),
            ):
                if first == second or (first, second) in checked:
                    continue
                checked.add((first, second))
                first_label = moved.get(first, self.labels[first])
                second_label = moved.get(second, self.labels[second])
                gap = _clockwise(
                    _angle(float(first_label["start_x"]), float(first_label["start_y"])),
                    _angle(float(second_label["start_x"]), float(second_label["start_y"])),
                )
                today = _clockwise(self.angle[first], self.angle[second])
                if today > 180.0:
                    if abs(gap - today) > 1e-6:
                        return False
                    continue
                if gap > 180.0 or gap < min(MIN_ORDER_GAP_DEG, today) - _EPS:
                    return False
        return True

    def _within_bounds(self, index: int, label: dict) -> bool:
        half = self.half[index]
        if half != 0 and _current_half_from_x(float(label["start_x"]), axis_eps_px=HEMISPHERE_AXIS_EPS_PX) != half:
            return False
        text, elbow = self._displacement(index, label)
        if text > MAX_TEXT_SHIFT_PX + _EPS or elbow > MAX_ELBOW_SHIFT_PX + _EPS:
            return False
        elbow_radius = self._elbow_radius(label)
        elbow_moved = abs(elbow_radius - self.elbow_radius[index]) > _EPS
        if self.side[index] == "out":
            min_radius = _label_bbox_min_radius(label, self.total_length)
            # Q4: never closer to the ring than the label's own elbow circle (or than it already was).
            if min_radius < min(self.min_radius[index], self.elbow_radius[index]) - 1e-6:
                return False
            return not elbow_moved or self.feature_radius[index] + 1.0 <= elbow_radius <= min_radius
        max_radius = _label_bbox_max_radius(label, self.total_length)
        if max_radius > max(self.max_radius[index], self.elbow_radius[index]) + 1e-6:
            return False
        return not elbow_moved or max_radius <= elbow_radius <= self.feature_radius[index] - 1.0

    # ---- moves --------------------------------------------------------------
    def _finish(self, label: dict) -> dict:
        """Set the leader end of a moved label; a chosen contact follows the arrival rule."""

        box = label_ink_box(label, self.total_length) if "leader_contact" in label else None
        if box is not None:
            label["leader_contact"] = arrival_contact(label, box)
            end = contact_end(label, box, label["leader_contact"])
        else:
            label.pop("leader_contact", None)
            end = _leader_start_point(label, self.total_length)
        label["leader_start_x"], label["leader_start_y"] = float(end[0]), float(end[1])
        return label

    def _with_start(self, index: int, x: float, y: float) -> dict:
        label = dict(self.labels[index])
        label["start_x"] = float(x)
        label["start_y"] = float(y)
        target = float(
            label.get(
                "target_angle_unwrapped",
                angle_from_middle_unwrapped(float(label["middle"]), self.total_length),
            )
        )
        label["angle_unwrapped"] = normalize_angle_near_reference(_angle(x, y), target)
        return label

    def _arc(self, run: Sequence[int], deltas: Sequence[float]) -> Moved:
        moved: Moved = {}
        for index, delta in zip(run, deltas):
            label = self.labels[index]
            x, y = float(label["start_x"]), float(label["start_y"])
            radius = math.hypot(x, y)
            angle = math.radians(_angle(x, y) + delta)
            moved[index] = self._with_start(index, radius * math.cos(angle), radius * math.sin(angle))
        return moved

    def _radial(self, run: Sequence[int], delta: float) -> Moved:
        moved: Moved = {}
        for index in run:
            label = self.labels[index]
            x, y = float(label["start_x"]), float(label["start_y"])
            radius = math.hypot(x, y)
            if radius <= _EPS or radius + delta <= 1.0:
                return {}
            factor = (radius + delta) / radius
            moved[index] = self._with_start(index, x * factor, y * factor)
        return moved

    def _elbow(self, run: Sequence[int], deltas: Sequence[float]) -> Moved:
        moved: Moved = {}
        for index, delta in zip(run, deltas):
            if abs(delta) < 1e-12:
                continue
            label = dict(self.labels[index])
            x, y = float(label["middle_x"]), float(label["middle_y"])
            radius = math.hypot(x, y)
            if radius <= _EPS:
                return {}
            label["middle_x"] = x * (radius + delta) / radius
            label["middle_y"] = y * (radius + delta) / radius
            moved[index] = label
        return moved

    def _push(self, index: int, delta: float) -> Moved:
        """Move one label by ``delta`` degrees and push the next labels of its side just clear."""

        side, at = self.position[index]
        indices = self.order[side]
        moved = self._arc([index], [delta])
        previous = index
        step = 1 if delta > 0 else -1
        for offset in range(1, min(MAX_PUSH, len(indices) - 1) + 1):
            following = indices[(at + step * offset) % len(indices)]
            if following in moved:
                break
            previous_box = _points_box(_label_fit_box(moved[previous], self.total_length))

            def hits(turn: float, following: int = following, previous_box: Box = previous_box) -> bool:
                label = self._arc([following], [turn])[following] if turn else self.labels[following]
                return not _apart(_points_box(_label_fit_box(label, self.total_length)), previous_box)

            if not hits(0.0):
                break
            low, high = 0.0, delta
            while hits(high) and abs(high) < 3.0 * abs(delta):  # neighbours on a turning arc may need more
                low, high = high, high + delta
            if not hits(high):
                for _ in range(8):
                    middle = 0.5 * (low + high)
                    if hits(middle):
                        low = middle
                    else:
                        high = middle
            moved.update(self._arc([following], [high]))
            previous = following
        return moved

    def _contact(self, index: int) -> Moved:
        """The leader meets the point its arrival selects, unless the last piece then runs through its text."""

        label = dict(self.labels[index])
        box = label_ink_box(label, self.total_length)
        if box is None:
            return {}
        contact = arrival_contact(label, box)
        if label.get("leader_contact") == contact:
            return {}
        end = contact_end(label, box, contact)
        elbow = (float(label["middle_x"]), float(label["middle_y"]))
        if _segment_intersects_rect_interior(elbow, end, (box.min_x, box.min_y, box.max_x, box.max_y)):
            return {}
        label["leader_contact"] = contact
        return {index: label}

    def _runs(self, cluster: set[int]) -> Iterator[list[int]]:
        by_side: dict[Side, list[int]] = {}
        for index in cluster:
            if self.movable[index]:
                by_side.setdefault(self.side[index], []).append(index)
        for side, members in sorted(by_side.items()):
            order = self.order[side]
            size = len(order)
            positions = sorted(self.position[index][1] for index in members)
            # Start after the largest cyclic gap, so a cluster across the start of the order stays contiguous.
            gaps = [((positions[(at + 1) % len(positions)] - positions[at]) % size or size, at) for at in range(len(positions))]
            start = positions[(max(gaps)[1] + 1) % len(positions)]
            offsets = {(position - start) % size for position in positions}
            low, high = min(offsets), max(offsets)
            for first in range(low - RUN_PAD, high + 1):
                for last in range(max(first, low), min(high + RUN_PAD, first + MAX_RUN - 1) + 1):
                    if last - first + 1 >= size or not any(first <= offset <= last for offset in offsets):
                        continue
                    yield [order[(start + offset) % size] for offset in range(first, last + 1)]

    def _candidates(self, cluster: set[int]) -> Iterator[tuple[int, str, Moved]]:
        """(tier, kind, moved): tier 1 moves one label (and pushes), tier 2 moves a run."""

        runs = list(self._runs(cluster))
        window = sorted({index for index in cluster if self.movable[index]}.union(*runs))
        for index in window:
            toward_ring = -1.0 if self.side[index] == "out" else 1.0
            for step in ELBOW_STEPS_PX:
                for sign in (1.0, -1.0):
                    yield 1, "elbow", self._elbow([index], [sign * step])
            for step in RADIAL_STEPS_PX:
                for sign in (1.0, -1.0):
                    yield 1, ("ring" if sign == toward_ring else "away"), self._radial([index], sign * step)
        for index in window:
            yield 1, "contact", self._contact(index)
        for index in sorted(index for index in cluster if self.movable[index]):
            for step in ARC_STEPS_DEG:
                for sign in (1.0, -1.0):
                    yield 1, "arc", self._push(index, sign * step)
        for run in runs:
            for step in ARC_STEPS_DEG:
                for sign in (1.0, -1.0):
                    yield 2, "arc", self._arc(run, [sign * step] * len(run))
            if len(run) < 2:
                continue
            toward_ring = -1.0 if self.side[run[0]] == "out" else 1.0
            for step in RUN_RADIAL_STEPS_PX:
                for sign in (1.0, -1.0):
                    yield 2, ("ring" if sign == toward_ring else "away"), self._radial(run, sign * step)
            last = len(run) - 1
            for step in ELBOW_RAMP_PX:
                for sign in (1.0, -1.0):
                    yield 2, "elbow", self._elbow(run, [sign * step * (last - at) for at in range(len(run))])
                    yield 2, "elbow", self._elbow(run, [sign * step * at for at in range(len(run))])

    # ---- search ---------------------------------------------------------------
    def run(self) -> bool:
        """Repair cluster by cluster, smallest first; stop at the first cluster no move improves.

        The search also stops when, at the rate of its last move, the
        remaining defects would need more than ``MAX_EVALUATIONS``.
        """

        found = self.defects()
        last_evaluations = 0
        while found and self.evaluations < MAX_EVALUATIONS:
            cluster = min(_clusters(found), key=lambda members: (len(members), sorted(members)))
            if not any(self.movable[index] for index in cluster):
                return False
            best: tuple[tuple[int, int, float, int], Moved, bool] | None = None
            tier_done = 0
            for tier, _kind, moved in self._candidates(cluster):
                if tier != tier_done and best is not None and best[2] and best[0][2] <= CLEAN_STOP_PX:
                    break  # a small clean fix of one label: runs are not tried
                tier_done = tier
                if self.evaluations >= MAX_EVALUATIONS:
                    break
                if not moved or not all(self._within_bounds(index, label) for index, label in moved.items()):
                    continue
                if not self._order_kept(moved):
                    continue
                self.evaluations += 1
                for label in moved.values():
                    self._finish(label)
                geometry = {index: self._geometry_of(label) for index, label in moved.items()}
                before = self._current_local(list(moved))
                after = self._local(geometry)
                gain = len(before) - len(after)
                if gain <= 0:
                    continue
                through = max(0, self._through(moved, geometry) - self._current_through(list(moved)))
                shift = max(max(self._displacement(index, label)) for index, label in moved.items())
                key = (-gain, through, shift, len(moved))
                if best is None or key < best[0]:
                    best = (key, moved, not after)
            if best is None:
                return False
            for index, label in best[1].items():
                self.labels[index] = label
            self._refresh()
            cleared = -best[0][0]
            found = self.defects()
            spent = self.evaluations - last_evaluations
            last_evaluations = self.evaluations
            if found and (self.evaluations * cleared + len(found) * spent) > MAX_EVALUATIONS * cleared:
                return False
        return not found


def _clusters(defects: set[Defect]) -> list[set[int]]:
    parent: dict[int, int] = {}

    def find(index: int) -> int:
        while parent.setdefault(index, index) != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for _kind, first, second in sorted(defects):
        parent[find(first)] = find(second)
    groups: dict[int, set[int]] = {}
    for _kind, first, second in defects:
        groups.setdefault(find(first), set()).update((first, second))
    return list(groups.values())


def defect_count(labels: Sequence[dict], total_length: int, definition_text: Sequence[Aabb] = ()) -> int:
    """Overlapping pairs, crossing pairs and labels over the definition text of one record, as the fit check counts them."""

    fit = sum(circular_label_fits(labels, total_length, definition_text=definition_text).values(), FITS)
    return fit.overlap_pairs + fit.crossing_pairs + fit.definition_overlaps


def repair_circular_labels(
    labels: Sequence[dict],
    total_length: int,
    *,
    definition_text: Sequence[Aabb] = (),
) -> list[dict] | None:
    """The record's labels with their slight defects repaired, or None to keep them as they are.

    ``labels`` are one record's drawn labels (embedded and horizontal
    external labels); the result keeps their order. None when the record has
    no defect, more than ``MAX_DEFECTS``, radial labels, or defects the repair
    cannot clear completely.
    """

    if any(label.get("placement") == "radial" or label.get("oriented_corners") for label in labels):
        return None
    if all(label.get("is_embedded") for label in labels):
        return None  # nothing to move; one label can still cover the definition text
    count = defect_count(labels, total_length, definition_text)
    if not 0 < count <= MAX_DEFECTS:
        return None
    repair = _Repair(labels, total_length, definition_text)
    if not repair.run():
        return None
    repaired = repair.labels
    if defect_count(repaired, total_length, definition_text):
        return None
    return repaired


__all__ = ["MAX_DEFECTS", "MAX_EVALUATIONS", "defect_count", "repair_circular_labels"]
