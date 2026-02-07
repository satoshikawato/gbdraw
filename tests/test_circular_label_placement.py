import math

from gbdraw.labels.circular import (
    angle_from_middle,
    improved_label_placement_fc,
    place_labels_on_arc_fc,
    x_overlap,
    y_overlap,
)


def _angle_of_label(label: dict) -> float:
    return math.degrees(math.atan2(label["start_y"], label["start_x"])) % 360.0


def _angle_diff(a: float, b: float) -> float:
    return abs(((a - b + 180.0) % 360.0) - 180.0)


def _make_label(middle: float, width_px: float = 80.0, height_px: float = 14.0) -> dict:
    return {
        "middle": middle,
        "width_px": width_px,
        "height_px": height_px,
        "is_inner": False,
        "feature_middle_x": 0.0,
        "feature_middle_y": 0.0,
    }


def _count_overlaps(labels: list[dict], total_length: int) -> int:
    overlaps = 0
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if y_overlap(labels[i], labels[j], total_length, 0.1) and x_overlap(labels[i], labels[j], 1):
                overlaps += 1
    return overlaps


def _target_delta_unwrapped(label: dict, total_length: int) -> float:
    angle = label.get("angle_unwrapped")
    if angle is None:
        angle = _angle_of_label(label)
    target = label.get("target_angle_unwrapped", 360.0 * (label["middle"] / total_length) - 90.0)
    return abs(angle - target)


def test_place_labels_stays_near_feature_angle_for_sparse_labels() -> None:
    total_length = 16569
    labels = [
        _make_label(100),
        _make_label(3000),
        _make_label(6000),
        _make_label(10000),
        _make_label(14000),
        _make_label(16000),
    ]

    placed = place_labels_on_arc_fc(
        labels,
        center_x=0.0,
        center_y=0.0,
        x_radius=430.0,
        y_radius=430.0,
        start_angle=0.0,
        end_angle=360.0,
        total_length=total_length,
    )

    for label in placed:
        expected = angle_from_middle(label["middle"], total_length)
        actual = _angle_of_label(label)
        assert _angle_diff(actual, expected) <= 5.0


def test_improved_label_placement_limits_angle_drift_under_density() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [_make_label(11800 + i * 110, width_px=92.0) for i in range(8)]

    for label in labels:
        angle = angle_from_middle(label["middle"], total_length)
        label["feature_middle_x"] = feature_radius * math.cos(math.radians(angle))
        label["feature_middle_y"] = feature_radius * math.sin(math.radians(angle))

    placed = place_labels_on_arc_fc(
        labels,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        start_angle=0.0,
        end_angle=360.0,
        total_length=total_length,
    )
    improved = improved_label_placement_fc(
        placed,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        feature_radius=feature_radius,
        total_length=total_length,
        start_angle=0.0,
        end_angle=360.0,
    )

    for label in improved:
        expected = angle_from_middle(label["middle"], total_length)
        actual = _angle_of_label(label)
        assert _angle_diff(actual, expected) <= 26.0


def test_improved_label_placement_resolves_dense_cluster_overlaps() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [_make_label(3600 + i * 85, width_px=82.0) for i in range(7)]

    for label in labels:
        angle = angle_from_middle(label["middle"], total_length)
        label["feature_middle_x"] = feature_radius * math.cos(math.radians(angle))
        label["feature_middle_y"] = feature_radius * math.sin(math.radians(angle))

    placed = place_labels_on_arc_fc(
        labels,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        start_angle=0.0,
        end_angle=360.0,
        total_length=total_length,
    )
    improved = improved_label_placement_fc(
        placed,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        feature_radius=feature_radius,
        total_length=total_length,
        start_angle=0.0,
        end_angle=360.0,
    )

    assert _count_overlaps(improved, total_length) == 0


def test_improved_label_placement_preserves_feature_order() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [_make_label(2200 + i * 70, width_px=84.0) for i in range(10)]

    for label in labels:
        angle = angle_from_middle(label["middle"], total_length)
        label["feature_middle_x"] = feature_radius * math.cos(math.radians(angle))
        label["feature_middle_y"] = feature_radius * math.sin(math.radians(angle))

    placed = place_labels_on_arc_fc(
        labels,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        start_angle=0.0,
        end_angle=360.0,
        total_length=total_length,
    )
    improved = improved_label_placement_fc(
        placed,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        feature_radius=feature_radius,
        total_length=total_length,
        start_angle=0.0,
        end_angle=360.0,
    )

    unwrapped_angles = [label["angle_unwrapped"] for label in improved]
    assert all(unwrapped_angles[i] < unwrapped_angles[i + 1] for i in range(len(unwrapped_angles) - 1))


def test_dense_triplet_keeps_center_label_closest_to_feature() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [
        _make_label(2450, width_px=92.0),
        _make_label(2520, width_px=92.0),
        _make_label(2590, width_px=92.0),
    ]

    for label in labels:
        angle = angle_from_middle(label["middle"], total_length)
        label["feature_middle_x"] = feature_radius * math.cos(math.radians(angle))
        label["feature_middle_y"] = feature_radius * math.sin(math.radians(angle))

    placed = place_labels_on_arc_fc(
        labels,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        start_angle=0.0,
        end_angle=360.0,
        total_length=total_length,
    )
    improved = improved_label_placement_fc(
        placed,
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        feature_radius=feature_radius,
        total_length=total_length,
        start_angle=0.0,
        end_angle=360.0,
    )

    left_delta = _target_delta_unwrapped(improved[0], total_length)
    middle_delta = _target_delta_unwrapped(improved[1], total_length)
    right_delta = _target_delta_unwrapped(improved[2], total_length)
    assert middle_delta <= left_delta + 0.05
    assert middle_delta <= right_delta + 0.05
