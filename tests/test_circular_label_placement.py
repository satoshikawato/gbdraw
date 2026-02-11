import math
from pathlib import Path
from types import SimpleNamespace

from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.labels.circular as circular_labels_module
import gbdraw.render.groups.circular.labels as circular_labels_group_module
import gbdraw.render.groups.circular.seq_record as circular_seq_record_group_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
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


def _load_hmmtdna_external_labels(*, label_font_size: float = 22.0) -> tuple[list[dict], int]:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=True,
        track_type="tuckin",
        resolve_overlaps=False,
        allow_inner_labels=False,
        label_font_size=label_font_size,
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


def _make_legend_collision_fixture() -> tuple[list[dict], int, SimpleNamespace, SimpleNamespace]:
    total_length = 4000
    labels = [
        {
            "middle": 1000,
            "start_x": 250.0,
            "start_y": 0.0,
            "middle_x": 180.0,
            "middle_y": 0.0,
            "feature_middle_x": 120.0,
            "feature_middle_y": 0.0,
            "width_px": 80.0,
            "height_px": 20.0,
            "is_inner": False,
            "is_embedded": False,
        }
    ]
    canvas_config = SimpleNamespace(
        legend_position="right",
        legend_offset_x=760.0,
        legend_offset_y=460.0,
        total_width=1000.0,
        total_height=1000.0,
        offset_x=500.0,
        offset_y=500.0,
    )
    legend_config = SimpleNamespace(
        legend_width=120.0,
        legend_height=120.0,
        color_rect_size=20.0,
    )
    return labels, total_length, canvas_config, legend_config


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


def test_hmmtdna_font22_labels_remain_close_without_overlaps() -> None:
    y_overlap_calls = 0
    original_y_overlap = circular_labels_module.y_overlap

    def counting_y_overlap(*args, **kwargs):
        nonlocal y_overlap_calls
        y_overlap_calls += 1
        return original_y_overlap(*args, **kwargs)

    circular_labels_module.y_overlap = counting_y_overlap
    try:
        external_labels, total_length = _load_hmmtdna_external_labels(label_font_size=22.0)
    finally:
        circular_labels_module.y_overlap = original_y_overlap

    assert len(external_labels) > 0
    assert _count_overlaps(external_labels, total_length) == 0
    assert _count_overlaps_with_min_gap(external_labels, total_length) == 0
    assert y_overlap_calls < 500000

    max_target_delta = max(_target_delta_unwrapped(label, total_length) for label in external_labels)
    assert max_target_delta <= 30.0

    max_leader_length = max(
        math.hypot(label["start_x"] - label["middle_x"], label["start_y"] - label["middle_y"])
        for label in external_labels
    )
    assert max_leader_length <= 220.0


def test_circular_assembly_reuses_precalculated_labels_once() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
    selected_features = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=True,
        track_type="tuckin",
        resolve_overlaps=False,
        allow_inner_labels=False,
        label_font_size=22.0,
    )

    counts = {"assemble": 0, "labels_group": 0, "seq_record_group": 0}
    original_assemble_prepare = circular_assemble_module.prepare_label_list
    original_labels_prepare = circular_labels_group_module.prepare_label_list
    original_seq_prepare = circular_seq_record_group_module.prepare_label_list

    def wrap_counter(counter_key: str, original_func):
        def wrapped(*args, **kwargs):
            counts[counter_key] += 1
            return original_func(*args, **kwargs)

        return wrapped

    circular_assemble_module.prepare_label_list = wrap_counter("assemble", original_assemble_prepare)
    circular_labels_group_module.prepare_label_list = wrap_counter("labels_group", original_labels_prepare)
    circular_seq_record_group_module.prepare_label_list = wrap_counter("seq_record_group", original_seq_prepare)
    try:
        _ = assemble_circular_diagram_from_record(
            record,
            config_dict=config_dict,
            selected_features_set=selected_features,
            output_prefix="tmp",
            legend="left",
        )
    finally:
        circular_assemble_module.prepare_label_list = original_assemble_prepare
        circular_labels_group_module.prepare_label_list = original_labels_prepare
        circular_seq_record_group_module.prepare_label_list = original_seq_prepare

    assert counts["assemble"] == 1
    assert counts["labels_group"] == 0
    assert counts["seq_record_group"] == 0


def test_label_legend_collision_prefers_label_shift() -> None:
    labels, total_length, canvas_config, legend_config = _make_legend_collision_fixture()
    original_legend_x = canvas_config.legend_offset_x
    original_legend_y = canvas_config.legend_offset_y
    original_start_x = labels[0]["start_x"]

    circular_assemble_module._resolve_label_legend_collisions(labels, total_length, canvas_config, legend_config)

    assert not circular_assemble_module._labels_collide_with_legend(labels, total_length, canvas_config, legend_config)
    assert canvas_config.legend_offset_x == original_legend_x
    assert canvas_config.legend_offset_y == original_legend_y
    assert labels[0]["start_x"] < original_start_x


def test_label_legend_collision_expands_canvas_when_fallback_needed() -> None:
    labels, total_length, canvas_config, legend_config = _make_legend_collision_fixture()
    original_width = canvas_config.total_width
    original_height = canvas_config.total_height

    original_shift = circular_assemble_module._try_shift_labels_away_from_legend
    original_move = circular_assemble_module._try_move_legend_away_from_labels
    circular_assemble_module._try_shift_labels_away_from_legend = lambda *args, **kwargs: False
    circular_assemble_module._try_move_legend_away_from_labels = lambda *args, **kwargs: False
    try:
        circular_assemble_module._resolve_label_legend_collisions(labels, total_length, canvas_config, legend_config)
    finally:
        circular_assemble_module._try_shift_labels_away_from_legend = original_shift
        circular_assemble_module._try_move_legend_away_from_labels = original_move

    assert canvas_config.total_width > original_width or canvas_config.total_height > original_height
    assert not circular_assemble_module._labels_collide_with_legend(labels, total_length, canvas_config, legend_config)
