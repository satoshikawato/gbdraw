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

After a repair, the labels around the repaired defects yield room (LP-D59):
the window is the labels of the original defects and ``YIELD_WINDOWS`` labels
on either side of them; it takes the leader ends of the arrival rule, its
yielders move up away from a crowd (upward yield, LP-D56), and a descent
lowers a least-squares cost of leader length, slant and displacement
(LP-D50 addendum 4). The yield keeps every rule of the repair and never adds a
defect or a leader through another label's text; if the record ends with
either, the repaired layout is kept. Records the repair does not change never
reach it.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

from ..layout.spatial import (
    Aabb,
    aabbs_touch,
    bounding_aabb,
    finite_aabb,
    touching_aabb_indices,
    union_aabbs,
)
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
    _segments_through_rect_interior,
    angle_from_middle_unwrapped,
    assign_leader_start_points,
    circular_label_fits,
    circular_label_shape,
    normalize_angle_near_reference,
)
from .circular_leaders import arrival_contact, contact_end, label_ink_box
from .circular_types import CircularLabel
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

# LP-D59 yield (study sections 6 and 7): cost weights per label, in px² (slant in degrees).
YIELD_LAMBDA = 16.0  # per deg² of leader slant
YIELD_MU = 16.0  # per px² of text or elbow displacement from the engine's place
YIELD_W_LEVEL = 4.0  # per px² of height change of a last piece that arrives at a side of its text
YIELD_W_RULE = 20_000.0  # a window label whose leader end is not the arrival rule's point
YIELD_W_THROUGH = 10_000.0  # a leader through another label's text (removing one is rewarded)
YIELD_MIN_GAIN = 0.5
# The window: the labels of the original defects and this many neighbours on their side, growing.
YIELD_WINDOWS = (4, 8, 12)
YIELD_EVALUATIONS = 20_000
YIELD_ARC_STEPS_DEG = (0.25, 0.5, 1.0, 2.0)
YIELD_RADIAL_STEPS_PX = (2.0, 4.0, 8.0)
YIELD_ELBOW_STEPS_PX = (1.0, 2.0, 4.0, 8.0)
# Upward yield (LP-D56): yielders next to a crowd move this far away from it; the crowd follows, tapered.
LIFT_STEPS_PX = (6.0, 12.0, 18.0, 24.0)
LIFT_YIELDERS = (2, 4, 6)
LIFT_CROWD = (4, 8)
LIFT_ROUNDS = 4
# Neighbouring rows stay this far apart (or as far as the engine put them, if closer and within the cap).
ROW_GAP_PX = 2.0
ROW_GAP_CAP_PX = 30.0

_EPS = 1e-9

Side: TypeAlias = Literal["out", "in"]
DefectKind: TypeAlias = Literal["overlap", "crossing", "definition"]
Defect: TypeAlias = tuple[DefectKind, int, int]
Moved: TypeAlias = dict[int, CircularLabel]
# (leader of, text of): a leader drawn through another label's text.
Through: TypeAlias = frozenset[tuple[int, int]]


@dataclass(frozen=True)
class _Geometry:
    """One label as the counts see it, with the boxes the broad phase compares.

    ``box`` is None without text parts, ``leader_box`` without a leader
    (:func:`~gbdraw.layout.spatial.aabbs_touch` then reads it as touching).
    """

    shape: LabelShape
    box: Aabb | None
    leader_box: Aabb | None
    reach: Aabb | None


def _angle(x: float, y: float) -> float:
    return math.degrees(math.atan2(y, x))


def _slant(anchor: tuple[float, float], end: tuple[float, float]) -> float:
    """Degrees between the chord from a leader's feature anchor to its end and the radial line at the anchor."""

    chord = math.hypot(end[0] - anchor[0], end[1] - anchor[1])
    radius = math.hypot(*anchor)
    if chord <= 1e-6 or radius <= 1e-6:
        return 0.0
    dot = ((end[0] - anchor[0]) * anchor[0] + (end[1] - anchor[1]) * anchor[1]) / (chord * radius)
    return math.degrees(math.acos(min(1.0, abs(dot))))


def _clockwise(first: float, second: float) -> float:
    """The angle from ``first`` to ``second`` in the direction of increasing SVG angle, in [0, 360)."""

    return (second - first) % 360.0


class _Repair:
    """The search state of one record: copies of its labels, never the caller's."""

    def __init__(self, labels: Sequence[CircularLabel], total_length: int, definition_text: Sequence[Aabb]) -> None:
        self.total_length = int(total_length)
        self.definition_text = tuple(definition_text)
        self.labels = [label.copy() for label in labels]
        self.original = [label.copy() for label in labels]
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
        self.yield_evaluations = 0
        # The labels of the original defects (set by run) and the yield's current window.
        self.seeds: frozenset[int] = frozenset()
        self.window: frozenset[int] = frozenset()
        # The engine's row gap of each pair of neighbours (key: moved label, neighbour).
        self._original_gaps: dict[tuple[int, int], float] = {}
        self._refresh()
        # Pairs that a move within the bounds can bring together.
        tallest = max((float(self.labels[index]["height_px"]) for index in self.external), default=0.0)
        expand = 2.0 * MAX_TEXT_SHIFT_PX + 2.0 * tallest + 4.0
        reaches = [geometry.reach for geometry in self.geometry]
        self.peers: list[list[int]] = [
            [
                other
                for other in touching_aabb_indices(
                    reaches, reach.expanded(expand) if reach is not None else None, skip=index
                )
                if self.movable[index] or self.movable[other]
            ]
            for index, reach in enumerate(reaches)
        ]

    # ---- geometry -------------------------------------------------------
    @staticmethod
    def _elbow_radius(label: CircularLabel) -> float:
        return math.hypot(float(label["middle_x"]), float(label["middle_y"]))

    def _geometry_of(self, label: CircularLabel) -> _Geometry:
        shape = circular_label_shape(label, self.total_length)
        text = [point for part in shape.parts for point in part]
        if not shape.leader:
            box = bounding_aabb(text)
            return _Geometry(shape, box, None, box)
        box = bounding_aabb(text)
        leader_box = bounding_aabb(shape.leader)
        reach = union_aabbs((box, leader_box)) if box is not None and leader_box is not None else None
        return _Geometry(shape, box, leader_box, reach)

    def _fit_box(self, label: CircularLabel) -> Aabb | None:
        """The label's fit box (`_label_fit_box`) as a box; None for a non-finite bound."""

        (min_x, min_y), _, (max_x, max_y), _ = _label_fit_box(label, self.total_length)
        return finite_aabb(min_x, min_y, max_x, max_y)

    def _text_rect(self, label: CircularLabel) -> tuple[tuple[float, float, float, float], Aabb]:
        rect = _label_bbox_rect(label, self.total_length, 0.0)
        return rect, Aabb(*rect)

    def _refresh(self) -> None:
        """Recompute every engine leader end (they depend on the neighbours) and the geometry."""

        assign_leader_start_points([self.labels[index] for index in self.external], self.total_length)
        self.geometry = [self._geometry_of(label) for label in self.labels]
        self.text_rects = {index: self._text_rect(self.labels[index]) for index in self.external}
        self._local_cache: dict[tuple[int, ...], frozenset[Defect]] = {}
        self._through_cache: dict[tuple[int, ...], Through] = {}

    # ---- counts (the fit check's predicates) ------------------------------
    def _pair_defects(self, first: int, first_geometry: _Geometry, second: int, second_geometry: _Geometry) -> list[Defect]:
        if not aabbs_touch(first_geometry.reach, second_geometry.reach):
            return []
        low, high = (first, second) if first < second else (second, first)
        found: list[Defect] = []
        if (
            first_geometry.shape.parts
            and second_geometry.shape.parts
            and aabbs_touch(first_geometry.box, second_geometry.box)
            and shapes_overlap(first_geometry.shape.parts, second_geometry.shape.parts, touching=True)
        ):
            found.append(("overlap", low, high))
        if (
            first_geometry.shape.leader
            and second_geometry.shape.leader
            and aabbs_touch(first_geometry.leader_box, second_geometry.leader_box)
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

    def _through(self, moved: Moved, geometry: dict[int, _Geometry]) -> Through:
        """Leaders drawn through another label's text, for pairs that involve a moved label."""

        rects = {index: self._text_rect(label) for index, label in moved.items()}
        found: set[tuple[int, int]] = set()
        for index in geometry:
            for other in self.peers[index]:
                for leader_of, text_of in ((index, other), (other, index)):
                    if not self.movable[text_of]:
                        continue  # an embedded label's text sits on its feature
                    leader_geometry = geometry.get(leader_of, self.geometry[leader_of])
                    leader = leader_geometry.shape.leader
                    rect, rect_box = rects.get(text_of, self.text_rects[text_of])
                    if not leader or not aabbs_touch(leader_geometry.leader_box, rect_box):
                        continue
                    if _segments_through_rect_interior(zip(leader, leader[1:]), rect):
                        found.add((leader_of, text_of))
        return frozenset(found)

    def _current_through(self, indices: Sequence[int]) -> Through:
        key = tuple(sorted(indices))
        cached = self._through_cache.get(key)
        if cached is None:
            cached = self._through({index: self.labels[index] for index in key}, {index: self.geometry[index] for index in key})
            self._through_cache[key] = cached
        return cached

    # ---- rules ------------------------------------------------------------
    def _displacement(self, index: int, label: CircularLabel) -> tuple[float, float]:
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

    def _within_bounds(self, index: int, label: CircularLabel) -> bool:
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
    def _finish(self, label: CircularLabel) -> CircularLabel:
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

    def _with_start(self, index: int, x: float, y: float) -> CircularLabel:
        label = self.labels[index].copy()
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
            label = self.labels[index].copy()
            x, y = float(label["middle_x"]), float(label["middle_y"])
            radius = math.hypot(x, y)
            if radius <= _EPS:
                return {}
            label["middle_x"] = x * (radius + delta) / radius
            label["middle_y"] = y * (radius + delta) / radius
            moved[index] = label
        return moved

    def _push(self, index: int, delta: float, within: frozenset[int] | None = None) -> Moved:
        """Move one label by ``delta`` degrees and push the next labels of its side just clear.

        Empty when the push would move a label outside ``within`` (the yield's window).
        """

        side, at = self.position[index]
        indices = self.order[side]
        moved = self._arc([index], [delta])
        previous = index
        step = 1 if delta > 0 else -1
        for offset in range(1, min(MAX_PUSH, len(indices) - 1) + 1):
            following = indices[(at + step * offset) % len(indices)]
            if following in moved:
                break
            previous_box = self._fit_box(moved[previous])
            # The following label turns on its circle as `_arc` turns it; only its fit box is tested,
            # on one copy (no label copy per test).
            turned = self.labels[following].copy()
            x, y = float(turned["start_x"]), float(turned["start_y"])
            radius, angle = math.hypot(x, y), _angle(x, y)

            def hits(
                turn: float,
                following: int = following,
                previous_box: Aabb | None = previous_box,
                turned: CircularLabel = turned,
                radius: float = radius,
                angle: float = angle,
            ) -> bool:
                if not turn:
                    return aabbs_touch(self.geometry[following].box, previous_box)
                radians = math.radians(angle + turn)
                turned["start_x"], turned["start_y"] = radius * math.cos(radians), radius * math.sin(radians)
                return aabbs_touch(self._fit_box(turned), previous_box)

            if not hits(0.0):
                break
            if within is not None and following not in within:
                return {}
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

    def _takes_rule_contact(self, label: CircularLabel) -> bool:
        """``label`` (a copy) takes the point its arrival selects, unless the last piece then runs through its text."""

        box = label_ink_box(label, self.total_length)
        if box is None:
            return False
        contact = arrival_contact(label, box)
        end = contact_end(label, box, contact)
        elbow = (float(label["middle_x"]), float(label["middle_y"]))
        if _segment_intersects_rect_interior(elbow, end, (box.min_x, box.min_y, box.max_x, box.max_y)):
            return False
        label["leader_contact"] = contact
        return True

    def _contact(self, index: int) -> Moved:
        """The leader meets the point its arrival selects, unless the last piece then runs through its text."""

        label = self.labels[index].copy()
        previous = label.get("leader_contact")
        if not self._takes_rule_contact(label) or label.get("leader_contact") == previous:
            return {}
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
        self.seeds = frozenset(index for _kind, *pair in found for index in pair if self.movable[index])
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
                through = max(0, len(self._through(moved, geometry)) - len(self._current_through(list(moved))))
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

    # ---- yield (LP-D59) -------------------------------------------------------
    def _yield_cost(self, index: int, label: CircularLabel, geometry: _Geometry) -> float:
        """One label's share of the yield's least-squares cost (study section 6)."""

        text, elbow_shift = self._displacement(index, label)
        cost = YIELD_MU * (text * text + elbow_shift * elbow_shift)
        if index in self.window and "leader_contact" not in label:
            cost += YIELD_W_RULE
        leader = geometry.shape.leader
        if len(leader) < 3:
            return cost
        anchor, elbow, end = leader[0], leader[1], leader[-1]
        length = math.dist(anchor, elbow) + math.dist(elbow, end)
        cost += length * length + YIELD_LAMBDA * _slant(anchor, end) ** 2
        box = geometry.box
        if box is not None and not box.min_x <= elbow[0] <= box.max_x:
            cost += YIELD_W_LEVEL * (elbow[1] - end[1]) ** 2  # the last piece arrives at a side
        return cost

    def _row_gap(self, first: CircularLabel, second: CircularLabel) -> float:
        a = self._fit_box(first)
        b = self._fit_box(second)
        if a is None or b is None:
            return math.inf
        return max(b.min_y - a.max_y, a.min_y - b.max_y, b.min_x - a.max_x, a.min_x - b.max_x, 0.0)

    def _row_gaps_kept(self, moved: Moved) -> bool:
        """Neighbouring rows of a moved label stay ``ROW_GAP_PX`` apart, or as far as the engine put them."""

        for index in moved:
            side, at = self.position[index]
            indices = self.order[side]
            for other in (indices[(at - 1) % len(indices)], indices[(at + 1) % len(indices)]):
                if other == index or (other in moved and other < index):
                    continue
                today = self._original_gaps.get((index, other))
                if today is None:
                    today = self._original_gaps[index, other] = self._row_gap(self.original[index], self.original[other])
                if today > ROW_GAP_CAP_PX:
                    continue
                gap = self._row_gap(moved[index], moved.get(other, self.labels[other]))
                if gap < min(ROW_GAP_PX, today) - 1e-6:
                    return False
        return True

    def _window_runs(self, width: int) -> list[list[int]]:
        """The seeds and ``width`` labels on either side of them, as runs of neighbours in each side's order."""

        runs: list[list[int]] = []
        for side, order in self.order.items():
            size = len(order)
            positions = {self.position[index][1] for index in self.seeds if self.side[index] == side}
            if not positions:
                continue
            taken = {(position + step) % size for position in positions for step in range(-width, width + 1)}
            if len(taken) >= size:
                runs.append(list(order))
                continue
            start = next(at for at in range(size) if at in taken and (at - 1) % size not in taken)
            run: list[int] = []
            for step in range(size):
                at = (start + step) % size
                if at in taken:
                    run.append(order[at])
                elif run:
                    runs.append(run)
                    run = []
            if run:
                runs.append(run)
        return runs

    def _try_yield(self, moved: Moved) -> float | None:
        """The cost change of a yield move, or None when it breaks a rule.

        The rules: window labels only, the repair's bounds and order, row
        gaps, no new defect, and no leader through another label's text that
        the moved labels' pairs do not have now (both directions, the
        repair's test). A moved label takes the arrival rule's point.
        """

        if not moved or not moved.keys() <= self.window:
            return None
        if not all(self._within_bounds(index, label) for index, label in moved.items()):
            return None
        if not self._order_kept(moved) or not self._row_gaps_kept(moved):
            return None
        if not all(self._takes_rule_contact(label) for label in moved.values()):
            return None
        self.yield_evaluations += 1
        for label in moved.values():
            self._finish(label)
        geometry = {index: self._geometry_of(label) for index, label in moved.items()}
        if self._local(geometry) - self._current_local(list(moved)):
            return None
        before = self._current_through(list(moved))
        after = self._through(moved, geometry)
        if after - before:
            return None
        change = sum(
            self._yield_cost(index, label, geometry[index]) - self._yield_cost(index, self.labels[index], self.geometry[index])
            for index, label in moved.items()
        )
        return change + YIELD_W_THROUGH * (len(after) - len(before))

    def _apply_yield(self, moved: Moved, allowed: Through) -> bool:
        """Apply a yield move; undo it when the leader ends the engine reassigns to other labels then add a defect
        or a leader through text outside ``allowed``."""

        previous = {index: self.labels[index] for index in moved}
        ends = [(label.get("leader_start_x"), label.get("leader_start_y")) for label in self.labels]
        for index, label in moved.items():
            self.labels[index] = label
        self._refresh()
        changed = [
            index
            for index, label in enumerate(self.labels)
            if index not in moved and (label.get("leader_start_x"), label.get("leader_start_y")) != ends[index]
        ]
        if not changed or (not self._current_local(changed) and self._current_through(changed) <= allowed):
            return True
        for index, label in previous.items():
            self.labels[index] = label
        self._refresh()
        return False

    def _rule_contacts(self, members: Sequence[int], allowed: Through) -> int:
        """Window labels take the arrival rule's point where that adds no defect and no leader through text."""

        applied = 0
        for index in members:
            if "leader_contact" in self.labels[index]:
                continue
            label = self.labels[index].copy()
            if not self._takes_rule_contact(label):
                continue
            moved = {index: self._finish(label)}
            geometry = {index: self._geometry_of(label)}
            self.yield_evaluations += 1
            if self._local(geometry) - self._current_local([index]):
                continue
            if self._through(moved, geometry) - self._current_through([index]):
                continue
            applied += self._apply_yield(moved, allowed)
        return applied

    def _lift(self, yielders: Sequence[int], crowd: Sequence[int], distance: float) -> Moved:
        """Yielders move ``distance`` px vertically away from the crowd; the crowd follows with a linear taper."""

        def centre(index: int) -> float:
            box = self.geometry[index].box
            return float(self.labels[index]["start_y"]) if box is None else 0.5 * (box.min_y + box.max_y)

        if not yielders or not crowd:
            return {}
        away = sum(map(centre, yielders)) / len(yielders) > sum(map(centre, crowd)) / len(crowd)
        sign = 1.0 if away else -1.0
        moved: Moved = {}
        shifts = [(index, distance) for index in yielders]
        shifts += [(index, distance * (len(crowd) - at) / (len(crowd) + 1)) for at, index in enumerate(crowd)]
        for index, shift in shifts:
            label = self.labels[index]
            moved[index] = self._with_start(index, float(label["start_x"]), float(label["start_y"]) + sign * shift)
        return moved

    def _lift_pass(self, runs: Sequence[Sequence[int]], allowed: Through) -> int:
        """Upward yield (LP-D56): best-improvement lifts over the split point, yielders, crowd and distance."""

        applied = 0
        for _ in range(LIFT_ROUNDS):
            best: tuple[float, Moved] | None = None
            for run in runs:
                for split in range(1, len(run)):
                    for yielder_count in LIFT_YIELDERS:
                        for crowd_count in LIFT_CROWD:
                            for yielders, crowd in (
                                (run[max(0, split - yielder_count) : split], run[split : split + crowd_count]),
                                (run[split : split + yielder_count], run[max(0, split - crowd_count) : split][::-1]),
                            ):
                                for distance in LIFT_STEPS_PX:
                                    if self.yield_evaluations >= YIELD_EVALUATIONS:
                                        break
                                    moved = self._lift(yielders, crowd, distance)
                                    change = self._try_yield(moved)
                                    if change is not None and change < -YIELD_MIN_GAIN and (best is None or change < best[0]):
                                        best = (change, moved)
            if best is None or not self._apply_yield(best[1], allowed):
                break
            applied += 1
        return applied

    def _yield_moves(self, index: int) -> Iterator[Moved]:
        for step in YIELD_ARC_STEPS_DEG:
            for sign in (1.0, -1.0):
                yield self._push(index, sign * step, self.window)
        for step in YIELD_RADIAL_STEPS_PX:
            for sign in (1.0, -1.0):
                yield self._radial([index], sign * step)
        for step in YIELD_ELBOW_STEPS_PX:
            for sign in (1.0, -1.0):
                yield self._elbow([index], [sign * step])

    def make_room(self) -> bool:
        """Yield room around the repaired defects (LP-D59); False when the record ends with a leader through
        another label's text that the repaired layout did not have (the caller then keeps the repaired layout)."""

        allowed = self._current_through(self.external)
        for width in YIELD_WINDOWS:
            runs = self._window_runs(width)
            members = [index for run in runs for index in run]
            self.window = frozenset(members)
            applied = self._rule_contacts(members, allowed) + self._lift_pass(runs, allowed)
            improved = True
            while improved and self.yield_evaluations < YIELD_EVALUATIONS:
                improved = False
                for index in members:
                    for moved in self._yield_moves(index):
                        if self.yield_evaluations >= YIELD_EVALUATIONS:
                            break
                        change = self._try_yield(moved)
                        if change is not None and change < -YIELD_MIN_GAIN and self._apply_yield(moved, allowed):
                            applied += 1
                            improved = True
                            break
            if not applied:
                break
        return not self.defects() and self._current_through(self.external) <= allowed


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


def defect_count(labels: Sequence[CircularLabel], total_length: int, definition_text: Sequence[Aabb] = ()) -> int:
    """Overlapping pairs, crossing pairs and labels over the definition text of one record, as the fit check counts them."""

    fit = sum(circular_label_fits(labels, total_length, definition_text=definition_text).values(), FITS)
    return fit.overlap_pairs + fit.crossing_pairs + fit.definition_overlaps


def repair_circular_labels(
    labels: Sequence[CircularLabel],
    total_length: int,
    *,
    definition_text: Sequence[Aabb] = (),
    defects: int | None = None,
) -> list[CircularLabel] | None:
    """The record's labels with their slight defects repaired, or None to keep them as they are.

    ``labels`` are one record's drawn labels (embedded and horizontal
    external labels); the result keeps their order. None when the record has
    no defect, more than ``MAX_DEFECTS``, radial labels, or defects the repair
    cannot clear completely. A repaired record then yields room around its
    repaired defects (LP-D59). ``defects``: the caller's :func:`defect_count`
    of ``labels``, when it has one.
    """

    if any(label.get("placement") == "radial" or label.get("oriented_corners") for label in labels):
        return None
    if all(label.get("is_embedded") for label in labels):
        return None  # nothing to move; one label can still cover the definition text
    count = defect_count(labels, total_length, definition_text) if defects is None else defects
    if not 0 < count <= MAX_DEFECTS:
        return None
    repair = _Repair(labels, total_length, definition_text)
    if not repair.run():
        return None
    repaired = [label.copy() for label in repair.labels]
    if defect_count(repaired, total_length, definition_text):
        return None
    # LP-D59: only a repaired record yields; all or nothing.
    if repair.make_room() and not defect_count(repair.labels, total_length, definition_text):
        return repair.labels
    return repaired


__all__ = ["MAX_DEFECTS", "MAX_EVALUATIONS", "defect_count", "repair_circular_labels"]
