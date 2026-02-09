import math
from pathlib import Path

from Bio import SeqIO

from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.circular import (
    angle_from_middle,
    improved_label_placement_fc,
    minimum_bbox_gap_px,
    place_labels_on_arc_fc,
    prepare_label_list,
    x_overlap,
    y_overlap,
)
from gbdraw.labels.filtering import preprocess_label_filtering


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


def _count_overlaps_with_min_gap(labels: list[dict], total_length: int) -> int:
    overlaps = 0
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            min_gap_px = minimum_bbox_gap_px(labels[i], labels[j], base_margin_px=0.0)
            if y_overlap(labels[i], labels[j], total_length, min_gap_px) and x_overlap(labels[i], labels[j], min_gap_px):
                overlaps += 1
    return overlaps


def _target_delta_unwrapped(label: dict, total_length: int) -> float:
    angle = label.get("angle_unwrapped")
    if angle is None:
        angle = _angle_of_label(label)
    target = label.get("target_angle_unwrapped", 360.0 * (label["middle"] / total_length) - 90.0)
    return abs(angle - target)


def _load_mjenmv_external_labels_without_blacklist() -> tuple[list[dict], int]:
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=True,
        track_type="tuckin",
        resolve_overlaps=False,
        allow_inner_labels=False,
        label_blacklist="",
    )
    cfg = GbdrawConfig.from_dict(config_dict)

    default_colors = load_default_colors("", "default")
    color_table = None
    color_table, default_colors = preprocess_color_tables(color_table, default_colors)
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())

    selected_features = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]
    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        selected_features,
        default_colors,
        cfg.canvas.strandedness,
        cfg.canvas.resolve_overlaps,
        label_filtering,
    )

    labels = prepare_label_list(
        feature_dict,
        len(record.seq),
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        config_dict,
        cfg=cfg,
    )
    external_labels = [label for label in labels if not label.get("is_embedded")]
    return external_labels, len(record.seq)


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


def test_dense_mito_triplet_respects_minimum_bbox_gap() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [
        _make_label(12160, width_px=86.0),
        _make_label(12205, width_px=86.0),
        _make_label(12250, width_px=86.0),
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

    assert _count_overlaps_with_min_gap(improved, total_length) == 0


def test_dense_wide_labels_expand_search_beyond_default_shift() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [
        _make_label(7500, width_px=260.0, height_px=24.0),
        _make_label(7560, width_px=260.0, height_px=24.0),
        _make_label(7620, width_px=260.0, height_px=24.0),
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
        max_angle_shift_deg=25.0,
    )

    deltas = [_target_delta_unwrapped(label, total_length) for label in improved]
    assert max(deltas) > 0.5
    assert _count_overlaps_with_min_gap(improved, total_length) == 0


def test_dense_chain_cluster_moves_neighbor_labels_together() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    labels = [
        _make_label(1585, width_px=180.0, height_px=24.0),
        _make_label(1625, width_px=260.0, height_px=28.0),
        _make_label(1665, width_px=180.0, height_px=20.0),
        _make_label(1705, width_px=320.0, height_px=16.0),
        _make_label(1745, width_px=260.0, height_px=16.0),
        _make_label(1785, width_px=320.0, height_px=20.0),
        _make_label(1825, width_px=220.0, height_px=24.0),
        _make_label(1865, width_px=120.0, height_px=20.0),
        _make_label(1905, width_px=220.0, height_px=24.0),
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
        max_angle_shift_deg=25.0,
    )

    assert _count_overlaps_with_min_gap(improved, total_length) == 0


def test_wraparound_dense_cluster_resolves_after_target_reset_retry() -> None:
    total_length = 16569
    label_radius = 430.0
    feature_radius = 390.0
    middles = [146, 182, 320, 367, 446, 460, 16136, 16215, 16347, 16429, 16500, 16505]
    widths = [174.4, 192.1, 91.3, 159.4, 217.5, 188.7, 105.9, 112.3, 170.8, 199.2, 100.4, 197.7]
    heights = [26.7, 22.9, 28.7, 27.7, 26.5, 24.8, 29.4, 24.0, 25.9, 26.4, 20.8, 20.1]
    labels = [_make_label(m, width_px=w, height_px=h) for m, w, h in zip(middles, widths, heights)]

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
        max_angle_shift_deg=25.0,
    )

    assert _count_overlaps_with_min_gap(improved, total_length) == 0


def test_mjenmv_dense_labels_without_blacklist_have_no_outer_overlaps() -> None:
    external_labels, total_length = _load_mjenmv_external_labels_without_blacklist()
    assert len(external_labels) == 109
    assert _count_overlaps(external_labels, total_length) == 0
