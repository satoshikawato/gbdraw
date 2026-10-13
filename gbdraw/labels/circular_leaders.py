"""Leader lines of external circular labels, as the label drawer draws them."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, TypeAlias

from ..layout.text_geometry import Point
from .fit import polyline_crossing_count

LeaderSegment: TypeAlias = tuple[Point, Point]


def leader_polyline(label: Mapping[str, Any]) -> tuple[Point, Point, Point]:
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


def leader_segments(label: Mapping[str, Any]) -> tuple[LeaderSegment, LeaderSegment]:
    """Return the two segments of the drawn leader polyline."""

    feature, elbow, contact = leader_polyline(label)
    return ((feature, elbow), (elbow, contact))


def leader_crossing_count(labels: Sequence[Mapping[str, Any]]) -> int:
    """Count label pairs whose leader polylines cross (touching does not count)."""

    return polyline_crossing_count([leader_polyline(label) for label in labels])


__all__ = ["LeaderSegment", "leader_crossing_count", "leader_polyline", "leader_segments"]
