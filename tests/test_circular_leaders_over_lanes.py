"""A Circular label's leader runs over the features of other lanes between its label and its own feature (OV-3625).

With predicted ORFs in their own lanes outside the circle (Owner D-131 H1=A), the leaders of the annotated features
inside cross those lanes. Painted under them, a leader vanished under a predicted ORF and came back as short pieces
between the lanes, so its label seemed to name a predicted ORF. The rule is general: the part of a leader from its
elbow to its own feature's edge is painted above the features of other lanes it crosses, and the rest stays under
the features, so a leader still ends under its own feature. Leaders that cross no feature of another lane are drawn
as before, byte for byte.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.labels.circular_leaders import leader_over_point, leaders_over_lanes
from gbdraw.labels.circular_types import CircularLabel
from tests.utils.cli_in_process import run_gbdraw

LAMBDA = Path(__file__).resolve().parent / "test_inputs" / "NC_001416.gb"
SVG = "{http://www.w3.org/2000/svg}"
TOTAL = 10_000


def _point(radius: float, position: float) -> tuple[float, float]:
    angle = math.radians(360.0 * position / TOTAL - 90.0)
    return radius * math.cos(angle), radius * math.sin(angle)


def _label(position: float, *, inner: bool = False) -> CircularLabel:
    """An external label of a feature in the lane 380-400 px, its elbow at 520 px (outer) or 250 px (inner)."""

    elbow = _point(250.0 if inner else 520.0, position)
    anchor = _point(390.0, position)
    label = {
        "is_embedded": False,
        "is_inner": inner,
        "middle_x": elbow[0],
        "middle_y": elbow[1],
        "feature_middle_x": anchor[0],
        "feature_middle_y": anchor[1],
        "feature_anchor_x": anchor[0],
        "feature_anchor_y": anchor[1],
        "feature_inner_radius_px": 380.0,
        "feature_outer_radius_px": 400.0,
    }
    return label  # type: ignore[return-value]


def test_a_leader_goes_over_the_features_of_other_lanes_between_its_label_and_its_feature() -> None:
    intervals = [
        (900.0, 1100.0, 380.0, 400.0),  # the labelled feature's own lane
        (950.0, 1300.0, 410.0, 430.0),  # a lane outside it, across the leader at 1000
        (2900.0, 3100.0, 410.0, 430.0),  # outside, but away from the leader at 2000
        (5900.0, 6100.0, 340.0, 360.0),  # a lane inside, across the inner label's leader at 6000
        (9950.0, 10_000.0, 410.0, 430.0),  # outside, across the leader at 9990 (wraps round the origin)
    ]
    labels = [_label(1000.0), _label(2000.0), _label(6000.0, inner=True), _label(6000.0), _label(9990.0)]

    marked = leaders_over_lanes(labels, intervals, TOTAL)

    assert [label.get("leader_over_to_px") for label in marked] == [400.0, None, 380.0, None, 400.0]
    # A copy: the placed labels are not changed.
    assert all("leader_over_to_px" not in label for label in labels)


def test_a_neighbour_in_the_own_lane_does_not_move_the_leader_over_it() -> None:
    marked = leaders_over_lanes([_label(1000.0)], [(990.0, 1010.0, 380.0, 400.0), (995.0, 1200.0, 380.0, 400.0)], TOTAL)

    assert "leader_over_to_px" not in marked[0]


def test_the_leader_splits_where_it_meets_its_own_feature_edge() -> None:
    (label,) = leaders_over_lanes([_label(1000.0)], [(950.0, 1300.0, 410.0, 430.0)], TOTAL)

    split = leader_over_point(label)

    assert split is not None
    assert math.hypot(*split) == pytest.approx(400.0)
    assert leader_over_point(_label(1000.0)) is None


def _path_points(d: str) -> list[tuple[str, str]]:
    """The points a feature path passes through: those of ``M`` and ``L`` and the end of each arc (``A`` gives two
    radii, three flags and its end point)."""

    points: list[tuple[str, str]] = []
    for command, numbers in re.findall(r"([MLA])([^MLAz]*)", d):
        values = re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", numbers)
        if command == "A":
            values = values[5:7]
        points.extend(zip(values[0::2], values[1::2]))
    return points


def _paint_order(svg: str) -> tuple[list[tuple[int, str, list[tuple[float, float]]]], list[tuple[int, float, float]]]:
    """The leader parts (paint index, group id, points) and the features outside the circle (paint index, radii)."""

    root = ET.fromstring(svg)
    leaders: list[tuple[int, str, list[tuple[float, float]]]] = []
    features: list[tuple[int, float, float]] = []
    axis = float(next(element for element in root.iter(f"{SVG}circle")).get("r", "0"))
    index = 0

    def walk(element: ET.Element, group: str) -> None:
        nonlocal index
        index += 1
        group = element.get("id", group) if element.tag == f"{SVG}g" else group
        if element.tag == f"{SVG}line" and element.get("data-label-feature-id"):
            x1, y1, x2, y2 = (float(element.get(name, "0")) for name in ("x1", "y1", "x2", "y2"))
            leaders.append((index, group, [(x1 + (x2 - x1) * k / 40, y1 + (y2 - y1) * k / 40) for k in range(41)]))
        elif element.tag == f"{SVG}path" and element.get("data-gbdraw-feature-id"):
            radii = [math.hypot(float(x), float(y)) for x, y in _path_points(element.get("d", ""))]
            features.append((index, min(radii), max(radii)))
        for child in element:
            walk(child, group)

    walk(root, "")
    return leaders, [feature for feature in features if feature[1] > axis]


@pytest.mark.circular
def test_leaders_run_over_the_predicted_lanes_and_under_the_annotated_features(tmp_path: Path) -> None:
    record = SeqIO.read(LAMBDA, "genbank")[0:9000]
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [feature for feature in record.features if feature.type == "CDS"]
    SeqIO.write(record, tmp_path / "crop.gb", "genbank")

    def draw(*args: str) -> str:
        result = run_gbdraw(
            ["circular", "--gbk", "crop.gb", "--labels", "out", *args, "-o", "out", "-f", "svg", "--overwrite"],
            cwd=tmp_path,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return (tmp_path / "out.svg").read_text()

    leaders, predicted = _paint_order(draw("--predict_orfs"))
    assert predicted
    lanes_from = min(inner for _index, inner, _outer in predicted)
    lanes_to = max(outer for _index, _inner, outer in predicted)
    last_lane_paint = max(index for index, _inner, _outer in predicted)
    first_lane_paint = min(index for index, _inner, _outer in predicted)
    over = [points for index, _group, points in leaders if index > last_lane_paint]
    under = [points for index, _group, points in leaders if index < first_lane_paint]
    # Every leader point in the predicted lanes is painted above them; every leader part painted under them stays
    # out of the lanes (inside the circle, under the annotated features, or out past the lanes with the labels).
    assert over
    for points in under:
        assert all(not lanes_from <= math.hypot(x, y) <= lanes_to for x, y in points)
    # Each part painted above stops at the circle: it never reaches into the annotated features.
    for points in over:
        assert min(math.hypot(x, y) for x, y in points) >= lanes_from - 30.0

    # Without lanes between a label and its feature, the leaders are painted under the features as before.
    plain, _ = _paint_order(draw())
    assert {group for _index, group, _points in plain} == {"label_leaders"}
