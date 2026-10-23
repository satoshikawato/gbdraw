"""Leader lines of external circular labels, as the label drawer draws them.

A leader runs from its feature to an elbow and on to the label. The engine
picks the end on the label's font box (``circular._leader_start_point``); a
label whose ``leader_contact`` the local repair chose (LP-D52) ends instead
at one of the LP-D23 points of the box of its painted letters (the ink box,
LP-D25), 2.5 px short of it, and the point follows where the leader arrives
(LP-D55, :func:`arrival_contact`).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Literal, TypeAlias

from ..core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink
from ..layout.spatial import Aabb
from ..layout.text_geometry import Point
from .circular_types import CircularLabel
from .fit import polyline_crossing_count

LeaderSegment: TypeAlias = tuple[Point, Point]
TextAnchor: TypeAlias = Literal["middle", "start", "end"]
TextBaseline: TypeAlias = Literal["text-after-edge", "middle", "hanging"]
# A contact point as fractions (x, y) of the ink box, from its top-left corner.
LeaderContact: TypeAlias = tuple[float, float]

# LP-D25: a leader stops this far before its contact point.
CONTACT_GAP_PX = 2.5
# LP-D55: a leader arriving within this angle of horizontal meets the middle of the near side.
LEVEL_ARRIVAL_DEG = 25.0
# LP-D55: an elbow further than this share of the box height from its middle, on the
# feature's side, arrives from that side and meets the corner there.
ARRIVAL_BAND_SHARE = 0.25

_ANCHOR_SHARE: Mapping[TextAnchor, float] = {"start": 0.0, "middle": 0.5, "end": 1.0}


def horizontal_label_alignment(
    total_length: int,
    middle: float,
    start_x: float,
    *,
    is_inner: bool,
) -> tuple[TextAnchor, TextBaseline]:
    """The ``text-anchor`` and ``dominant-baseline`` a horizontal external label is drawn with."""

    angle = (360.0 * (middle / total_length)) % 360
    anchor: TextAnchor
    if start_x > 0:
        anchor = "end" if is_inner else "start"
    else:
        anchor = "start" if is_inner else "end"
    near_top = angle < 10 or angle >= 350
    if is_inner:
        return anchor, "hanging" if near_top else "middle"
    if near_top:
        return anchor, "text-after-edge"
    return anchor, "hanging" if 170 <= angle < 190 else "middle"


def label_ink_box(label: CircularLabel, total_length: int) -> Aabb | None:
    """The box of a horizontal label's painted letters, in record coordinates; None when it cannot be measured."""

    family = str(label.get("font_family", ""))
    ink = text_ink(str(label.get("label_text", "")), family)
    metrics = get_font_vertical_metrics(family)
    if ink is None or metrics is None or "font_size" not in label:
        return None
    start_x = float(label["start_x"])
    anchor, baseline = horizontal_label_alignment(
        total_length, float(label["middle"]), start_x, is_inner=bool(label.get("is_inner", False))
    )
    shift = dominant_baseline_shift_em(baseline, metrics)
    if shift is None:
        return None
    size = float(label["font_size"])
    origin_x = start_x - _ANCHOR_SHARE[anchor] * ink.advance * size
    baseline_y = float(label["start_y"]) + shift * size
    return Aabb(
        origin_x + ink.left * size,
        baseline_y - ink.top * size,
        origin_x + ink.right * size,
        baseline_y - ink.bottom * size,
    )


def contact_points(box: Aabb) -> tuple[LeaderContact, ...]:
    """The 12 LP-D23 points: 4 corners, the middles of the short sides, the middle and quarter points of the long sides."""

    corners: tuple[LeaderContact, ...] = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    if box.width >= box.height:
        return (
            *corners,
            (0.0, 0.5),
            (1.0, 0.5),
            *((share, edge) for edge in (0.0, 1.0) for share in (0.25, 0.5, 0.75)),
        )
    return (
        *corners,
        (0.5, 0.0),
        (0.5, 1.0),
        *((edge, share) for edge in (0.0, 1.0) for share in (0.25, 0.5, 0.75)),
    )


def _point(box: Aabb, contact: LeaderContact) -> Point:
    return (box.min_x + contact[0] * box.width, box.min_y + contact[1] * box.height)


def arrival_contact(label: CircularLabel, box: Aabb) -> LeaderContact:
    """The LP-D23 point a leader meets, from where it arrives (LP-D55, LP-D56).

    An elbow beside the box: the corner on the feature's side when the elbow
    lies on that side of the box's middle band (tRNA-Ser: from below, the
    bottom corner); else the middle of the near side when the last piece is
    within ``LEVEL_ARRIVAL_DEG`` of horizontal (tRNA-Asp, ribosomal protein
    S2); else the corner on the elbow's side. An elbow above or below the
    box: the nearest point of the facing long side.
    """

    elbow_x, elbow_y = float(label["middle_x"]), float(label["middle_y"])
    if box.width < box.height:
        return min(contact_points(box), key=lambda contact: math.dist(_point(box, contact), (elbow_x, elbow_y)))
    middle_y = 0.5 * (box.min_y + box.max_y)
    if box.min_x <= elbow_x <= box.max_x:
        edge = 0.0 if elbow_y < middle_y else 1.0
        return min(
            ((share, edge) for share in (0.25, 0.5, 0.75)),
            key=lambda contact: abs(box.min_x + contact[0] * box.width - elbow_x),
        )
    near = 0.0 if elbow_x < box.min_x else 1.0
    feature_below = float(label.get("feature_anchor_y", elbow_y)) > middle_y
    toward_feature = elbow_y - middle_y if feature_below else middle_y - elbow_y
    if toward_feature > ARRIVAL_BAND_SHARE * box.height:
        return (near, 1.0 if feature_below else 0.0)
    side_x = box.min_x + near * box.width
    if math.degrees(math.atan2(abs(elbow_y - middle_y), abs(elbow_x - side_x))) < LEVEL_ARRIVAL_DEG:
        return (near, 0.5)
    return (near, 1.0 if elbow_y > middle_y else 0.0)


def contact_end(label: CircularLabel, box: Aabb, contact: LeaderContact) -> Point:
    """Where a leader meeting ``contact`` of ``box`` ends: ``CONTACT_GAP_PX`` short of it, toward the elbow."""

    point_x, point_y = _point(box, contact)
    elbow_x, elbow_y = float(label["middle_x"]), float(label["middle_y"])
    distance = math.hypot(elbow_x - point_x, elbow_y - point_y)
    if distance <= CONTACT_GAP_PX:
        return point_x, point_y
    share = CONTACT_GAP_PX / distance
    return point_x + (elbow_x - point_x) * share, point_y + (elbow_y - point_y) * share


def chosen_contact_end(label: CircularLabel, total_length: int) -> Point | None:
    """The leader end of a label whose contact the repair chose, or None for the engine's own end."""

    contact = label.get("leader_contact")
    if contact is None:
        return None
    box = label_ink_box(label, total_length)
    if box is None:
        return None
    return contact_end(label, box, (float(contact[0]), float(contact[1])))


def leader_polyline(label: CircularLabel) -> tuple[Point, Point, Point]:
    """Return the drawn leader polyline: feature anchor -> elbow -> label contact."""

    feature = (
        float(label.get("feature_anchor_x", label.get("feature_middle_x", 0.0))),
        float(label.get("feature_anchor_y", label.get("feature_middle_y", 0.0))),
    )
    elbow = (float(label.get("middle_x", 0.0)), float(label.get("middle_y", 0.0)))
    contact = (
        float(label.get("leader_start_x", label.get("start_x", 0.0))),
        float(label.get("leader_start_y", label.get("start_y", 0.0))),
    )
    return (feature, elbow, contact)


def leader_segments(label: CircularLabel) -> tuple[LeaderSegment, LeaderSegment]:
    """Return the two segments of the drawn leader polyline."""

    feature, elbow, contact = leader_polyline(label)
    return ((feature, elbow), (elbow, contact))


def leader_crossing_count(labels: Sequence[CircularLabel]) -> int:
    """Count label pairs whose leader polylines cross (touching does not count)."""

    return polyline_crossing_count([leader_polyline(label) for label in labels])


def rim_text_placement(label: CircularLabel, record_length: int) -> tuple[Point, str, str, float]:
    """Where an outer or inner label's text is drawn: its point, anchor, baseline and turn (radial labels only)."""

    if label.get("placement") == "radial":
        return (
            (float(label.get("text_x", label["start_x"])), float(label.get("text_y", label["start_y"]))),
            str(label.get("text_anchor", "middle")),
            str(label.get("dominant_baseline", "middle")),
            float(label.get("rotation_deg", 0.0)),
        )
    anchor, baseline = horizontal_label_alignment(
        record_length, label["middle"], label["start_x"], is_inner=label.get("is_inner", False)
    )
    return (float(label["start_x"]), float(label["start_y"])), anchor, baseline, 0.0


__all__ = [
    "CONTACT_GAP_PX",
    "LeaderContact",
    "LeaderSegment",
    "TextAnchor",
    "TextBaseline",
    "arrival_contact",
    "chosen_contact_end",
    "contact_end",
    "contact_points",
    "horizontal_label_alignment",
    "label_ink_box",
    "leader_crossing_count",
    "leader_polyline",
    "leader_segments",
    "rim_text_placement",
]
