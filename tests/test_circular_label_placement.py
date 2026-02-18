import math
from pathlib import Path
from types import SimpleNamespace

from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.labels.circular as circular_labels_module
import gbdraw.render.groups.circular.labels as circular_labels_group_module
import gbdraw.render.groups.circular.seq_record as circular_seq_record_group_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.canvas.circular import CircularCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.sequence import determine_length_parameter
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.io.colors import load_default_colors
from gbdraw.layout.circular import calculate_feature_position_factors_circular
from gbdraw.layout.common import calculate_cds_ratio
from gbdraw.labels.circular import (
    angle_from_middle,
    improved_label_placement_fc,
    minimum_bbox_gap_px,
    place_labels_on_arc_fc,
    prepare_label_list,
    rearrange_labels_fc,
    x_overlap,
    y_overlap,
)
from gbdraw.labels.filtering import get_label_text, preprocess_label_filtering
from gbdraw.svg.arrows import calculate_circular_arrow_length


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


def _sum_leader_length(labels: list[dict]) -> float:
    return sum(
        math.hypot(
            float(label.get("leader_start_x", label["start_x"])) - float(label["middle_x"]),
            float(label.get("leader_start_y", label["start_y"])) - float(label["middle_y"]),
        )
        for label in labels
    )


def _label_bbox_min_radius(label: dict, total_length: int) -> float:
    min_x, max_x = circular_labels_module._label_x_bounds(label, minimum_margin=0.0)
    min_y, max_y = circular_labels_module._label_y_bounds(label, total_length, minimum_margin=0.0)
    closest_x = min(max(0.0, min_x), max_x)
    closest_y = min(max(0.0, min_y), max_y)
    return math.hypot(closest_x, closest_y)


def _label_bbox_boundary_min_clearance_against_features(
    label: dict,
    total_length: int,
    feature_intervals: list[tuple[float, float, float]],
    *,
    samples_per_edge: int = 25,
) -> float:
    min_x, max_x = circular_labels_module._label_x_bounds(label, minimum_margin=0.0)
    min_y, max_y = circular_labels_module._label_y_bounds(label, total_length, minimum_margin=0.0)
    if max_x <= min_x + 1e-9 or max_y <= min_y + 1e-9:
        return float("inf")

    min_clearance = float("inf")
    sample_count = max(3, int(samples_per_edge))
    for idx in range(sample_count):
        t = float(idx) / float(sample_count - 1)
        x = min_x + ((max_x - min_x) * t)
        y = min_y + ((max_y - min_y) * t)
        boundary_points = ((x, min_y), (x, max_y), (min_x, y), (max_x, y))

        for point_x, point_y in boundary_points:
            point_radius = math.hypot(point_x, point_y)
            point_position = ((math.degrees(math.atan2(point_y, point_x)) + 90.0) % 360.0) / 360.0 * float(total_length)
            required_feature_outer = circular_labels_module._max_outer_feature_radius_at_position(
                feature_intervals,
                point_position,
                total_length,
            )
            required_radius = required_feature_outer + circular_labels_module.MIN_OUTER_LABEL_TEXT_CLEARANCE_PX
            min_clearance = min(min_clearance, point_radius - required_radius)

    return min_clearance


def _count_half_plane_mismatches(labels: list[dict], total_length: int, axis_neutral_deg: float = 8.0) -> int:
    mismatch_count = 0
    for label in labels:
        feature_middle_x = float(label.get("feature_middle_x", 0.0))
        if abs(feature_middle_x) > 1e-6:
            preferred_half = 1 if feature_middle_x > 0 else -1
        else:
            target_angle = (360.0 * (label["middle"] / total_length) - 90.0) % 360.0
            vertical_axis_distance = min(_angle_diff(target_angle, 90.0), _angle_diff(target_angle, 270.0))
            if vertical_axis_distance <= axis_neutral_deg:
                continue
            preferred_half = 1 if math.cos(math.radians(target_angle)) >= 0 else -1

        start_x = float(label["start_x"])
        if abs(start_x) <= 1.0:
            continue
        actual_half = 1 if start_x > 0 else -1
        if actual_half != preferred_half:
            mismatch_count += 1
    return mismatch_count


def _target_delta_unwrapped(label: dict, total_length: int) -> float:
    angle = label.get("angle_unwrapped")
    if angle is None:
        angle = _angle_of_label(label)
    target = label.get("target_angle_unwrapped", 360.0 * (label["middle"] / total_length) - 90.0)
    return abs(angle - target)


def _label_unwrapped_angle_for_order_test(label: dict, total_length: int) -> float:
    angle = math.atan2(label["start_y"], label["start_x"])
    target = (2.0 * math.pi * (label["middle"] / total_length)) - (0.5 * math.pi)
    return angle + (2.0 * math.pi) * round((target - angle) / (2.0 * math.pi))


def _load_mjenmv_external_labels_with_config(
    *,
    strandedness: bool = True,
    resolve_overlaps: bool = False,
    track_type: str = "tuckin",
    label_blacklist: str = "",
) -> tuple[list[dict], int, GbdrawConfig]:
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=strandedness,
        track_type=track_type,
        resolve_overlaps=resolve_overlaps,
        allow_inner_labels=False,
        label_blacklist=label_blacklist,
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
    return external_labels, len(record.seq), cfg


def _load_mjenmv_external_labels_without_blacklist() -> tuple[list[dict], int]:
    external_labels, total_length, _ = _load_mjenmv_external_labels_with_config(
        strandedness=True,
        resolve_overlaps=False,
        track_type="tuckin",
        label_blacklist="",
    )
    return external_labels, total_length


def _load_hmmtdna_external_labels(*, label_font_size: float = 22.0) -> tuple[list[dict], int]:
    external_labels, total_length, _ = _load_hmmtdna_external_labels_with_config(
        label_font_size=label_font_size,
        strandedness=True,
        resolve_overlaps=False,
        track_type="tuckin",
    )
    return external_labels, total_length


def _load_hmmtdna_external_labels_with_config(
    *,
    label_font_size: float = 22.0,
    strandedness: bool = True,
    resolve_overlaps: bool = False,
    track_type: str = "tuckin",
) -> tuple[list[dict], int, GbdrawConfig]:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=strandedness,
        track_type=track_type,
        resolve_overlaps=resolve_overlaps,
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
    return external_labels, len(record.seq), cfg


def _load_hmmtdna_external_labels_with_feature_width(
    *,
    feature_width: float,
    label_font_size: float = 14.0,
    strandedness: bool = False,
    resolve_overlaps: bool = True,
    track_type: str = "middle",
) -> tuple[list[dict], int, GbdrawConfig, float]:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=strandedness,
        track_type=track_type,
        resolve_overlaps=resolve_overlaps,
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

    canvas_config = CircularCanvasConfigurator(
        output_prefix="tmp",
        config_dict=config_dict,
        legend="left",
        gb_record=record,
        cfg=cfg,
    )

    feature_track_ratio_factor_override = float(feature_width) / (
        float(canvas_config.radius) * float(canvas_config.track_ratio)
    )
    feature_radius_mapper, _, _ = circular_assemble_module._build_feature_radius_mapper(
        feature_dict,
        len(record.seq),
        canvas_config=canvas_config,
        cfg=cfg,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
    )
    default_anchor_px, default_arc_outer_px = circular_assemble_module._default_outer_label_arena(
        canvas_config=canvas_config,
        cfg=cfg,
    )
    if feature_radius_mapper is not None:
        default_anchor_px = float(feature_radius_mapper(default_anchor_px))
        default_arc_outer_px = float(feature_radius_mapper(default_arc_outer_px))
    if default_arc_outer_px < default_anchor_px:
        default_anchor_px, default_arc_outer_px = default_arc_outer_px, default_anchor_px
    outer_arena = (default_anchor_px, default_arc_outer_px)

    labels = prepare_label_list(
        feature_dict,
        len(record.seq),
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        config_dict,
        cfg=cfg,
        outer_arena=outer_arena,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
    )
    external_labels = [label for label in labels if not label.get("is_embedded")]
    return external_labels, len(record.seq), cfg, feature_track_ratio_factor_override


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


def _make_order_sensitive_legend_collision_fixture() -> tuple[list[dict], int, SimpleNamespace, SimpleNamespace]:
    total_length = 4000
    radius = 250.0
    labels = [
        {
            "middle": 1000,
            "start_x": radius * math.cos(0.0),
            "start_y": radius * math.sin(0.0),
            "middle_x": 180.0,
            "middle_y": 0.0,
            "feature_middle_x": 120.0,
            "feature_middle_y": 0.0,
            "width_px": 30.0,
            "height_px": 14.0,
            "is_inner": False,
            "is_embedded": False,
        },
        {
            "middle": 1051,
            "start_x": radius * math.cos(0.08),
            "start_y": radius * math.sin(0.08),
            "middle_x": 180.0,
            "middle_y": 0.0,
            "feature_middle_x": 120.0,
            "feature_middle_y": 0.0,
            "width_px": 30.0,
            "height_px": 14.0,
            "is_inner": False,
            "is_embedded": False,
        },
    ]
    canvas_config = SimpleNamespace(
        legend_position="right",
        legend_offset_x=730.0,
        legend_offset_y=520.0,
        total_width=1000.0,
        total_height=1000.0,
        offset_x=500.0,
        offset_y=500.0,
    )
    legend_config = SimpleNamespace(
        legend_width=120.0,
        legend_height=60.0,
        color_rect_size=20.0,
    )
    return labels, total_length, canvas_config, legend_config


def _make_label_leader_collision_fixture() -> tuple[list[dict], int]:
    """Return a minimal pair with a reproducible label-vs-leader collision."""
    total_length = 16569
    labels = [
        {
            "label_text": "tRNA-Gly",
            "middle": 10024.0,
            "start_x": -283.0202696366569,
            "start_y": 377.70662320104606,
            "middle_x": -281.28144819583736,
            "middle_y": 362.69543300936743,
            "feature_middle_x": -266.2517031780199,
            "feature_middle_y": 343.3154848747775,
            "feature_anchor_x": -266.2517031780199,
            "feature_anchor_y": 343.3154848747775,
            "width_px": 67.2,
            "height_px": 14.0,
            "is_inner": False,
            "is_embedded": False,
        },
        {
            "label_text": "NADH dehydrogenase subunit 3",
            "middle": 10231.0,
            "start_x": -355.20979826148726,
            "start_y": 374.414669818108,
            "middle_x": -338.77430088295495,
            "middle_y": 372.40950724473345,
            "feature_middle_x": -262.435772218122,
            "feature_middle_y": 288.4917077842585,
            "feature_anchor_x": -274.9014713984828,
            "feature_anchor_y": 302.1950639040108,
            "width_px": 235.2,
            "height_px": 14.0,
            "is_inner": False,
            "is_embedded": False,
        },
    ]
    return labels, total_length


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


def test_improved_label_placement_prefers_same_half_plane_when_overlap_tied() -> None:
    total_length = 20000
    label_radius = 430.0
    feature_radius = 390.0
    n = 12
    width_px = 120.0
    height_px = 24.0
    center = 18402
    spread = 60
    middles = [(center + int((idx - (n // 2)) * spread / n)) % total_length for idx in range(n)]
    labels = [_make_label(middle, width_px=width_px, height_px=height_px) for middle in middles]

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
    assert _count_half_plane_mismatches(improved, total_length) <= 2


def test_improved_label_placement_prefers_shorter_leaders_when_overlap_tied() -> None:
    total_length = 20000
    label_radius = 430.0
    feature_radius = 390.0
    middles = [2500, 2700, 2900]
    labels: list[dict] = []
    for middle in middles:
        target_angle = angle_from_middle(middle, total_length)
        shifted_angle = target_angle + 8.0
        labels.append(
            {
                "middle": middle,
                "width_px": 68.0,
                "height_px": 16.0,
                "is_inner": False,
                "middle_x": label_radius * math.cos(math.radians(target_angle)),
                "middle_y": label_radius * math.sin(math.radians(target_angle)),
                "start_x": label_radius * math.cos(math.radians(shifted_angle)),
                "start_y": label_radius * math.sin(math.radians(shifted_angle)),
                "feature_middle_x": feature_radius * math.cos(math.radians(target_angle)),
                "feature_middle_y": feature_radius * math.sin(math.radians(target_angle)),
            }
        )

    before_overlaps = _count_overlaps(labels, total_length)
    before_sum = _sum_leader_length(labels)
    improved = improved_label_placement_fc(
        [label.copy() for label in labels],
        center_x=0.0,
        center_y=0.0,
        x_radius=label_radius,
        y_radius=label_radius,
        feature_radius=feature_radius,
        total_length=total_length,
        start_angle=0.0,
        end_angle=360.0,
    )

    after_overlaps = _count_overlaps(improved, total_length)
    after_sum = _sum_leader_length(improved)
    unwrapped_angles = [label["angle_unwrapped"] for label in improved]

    assert after_overlaps <= before_overlaps
    assert all(unwrapped_angles[i] < unwrapped_angles[i + 1] for i in range(len(unwrapped_angles) - 1))
    assert after_sum < before_sum


def test_rearrange_labels_legacy_prefers_same_half_plane_when_overlap_tied() -> None:
    total_length = 200000
    feature_radius = 390.0
    n = 80
    width_px = 220.0
    height_px = 22.0
    center = 100000
    spread = 6000
    middles = [(center + int((idx - (n // 2)) * spread / n)) % total_length for idx in range(n)]
    labels: list[dict] = []
    for middle in middles:
        angle = angle_from_middle(middle, total_length)
        labels.append(
            {
                "middle": middle,
                "width_px": width_px,
                "height_px": height_px,
                "is_inner": False,
                "feature_middle_x": feature_radius * math.cos(math.radians(angle)),
                "feature_middle_y": feature_radius * math.sin(math.radians(angle)),
            }
        )

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    rearranged = rearrange_labels_fc(
        labels,
        feature_radius,
        total_length,
        "long",
        config_dict,
        strands="separate",
        is_outer=True,
        cfg=cfg,
    )

    assert _count_overlaps(rearranged, total_length) == 0
    assert _count_half_plane_mismatches(rearranged, total_length) <= 1


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


def test_effective_outer_middle_anchor_clearance_scales_with_feature_width_and_caps() -> None:
    base_clearance = float(circular_labels_module.MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX)

    zero_width_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=0.0,
    )
    assert math.isclose(zero_width_clearance, base_clearance + 2.0, rel_tol=1e-9, abs_tol=1e-9)

    small_width_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=10.0,
    )
    assert math.isclose(small_width_clearance, base_clearance + 3.0, rel_tol=1e-9, abs_tol=1e-9)

    capped_width_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(capped_width_clearance, base_clearance + 4.0, rel_tol=1e-9, abs_tol=1e-9)

    non_middle_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="tuckin",
        feature_band_width_px=80.0,
    )
    assert math.isclose(non_middle_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)

    stranded_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=True,
        strandedness=True,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(stranded_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)

    no_resolve_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=False,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(no_resolve_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)


def test_effective_outer_text_clearance_scales_with_feature_width_and_caps() -> None:
    base_clearance = float(circular_labels_module.MIN_OUTER_LABEL_TEXT_CLEARANCE_PX)

    zero_width_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=0.0,
    )
    assert math.isclose(zero_width_clearance, base_clearance + 1.0, rel_tol=1e-9, abs_tol=1e-9)

    small_width_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=10.0,
    )
    assert math.isclose(small_width_clearance, base_clearance + 1.5, rel_tol=1e-9, abs_tol=1e-9)

    capped_width_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(capped_width_clearance, base_clearance + 2.0, rel_tol=1e-9, abs_tol=1e-9)

    non_middle_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=True,
        strandedness=False,
        track_type="tuckin",
        feature_band_width_px=80.0,
    )
    assert math.isclose(non_middle_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)

    stranded_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=True,
        strandedness=True,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(stranded_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)

    no_resolve_clearance = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=False,
        strandedness=False,
        track_type="middle",
        feature_band_width_px=80.0,
    )
    assert math.isclose(no_resolve_clearance, base_clearance, rel_tol=1e-9, abs_tol=1e-9)


def test_mjenmv_resolve_overlaps_middle_has_no_outer_overlaps() -> None:
    external_labels, total_length, _ = _load_mjenmv_external_labels_with_config(
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
        label_blacklist="",
    )

    assert len(external_labels) == 109
    assert _count_overlaps(external_labels, total_length) == 0
    assert _count_overlaps_with_min_gap(external_labels, total_length) == 0


def test_mjenmv_resolve_overlaps_middle_keeps_wsv134_anchor_outside_wsv133_feature() -> None:
    external_labels, total_length, cfg = _load_mjenmv_external_labels_with_config(
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
        label_blacklist="",
    )
    assert external_labels

    wsv134_label = next((label for label in external_labels if label.get("label_text") == "wsv134-like protein"), None)
    wsv133_label = next((label for label in external_labels if label.get("label_text") == "wsv133-like protein"), None)
    assert wsv134_label is not None
    assert wsv133_label is not None

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    cds_ratio, _ = calculate_cds_ratio(cfg.canvas.circular.track_ratio, length_param, track_ratio_factor)
    feature_band_width_px = float(cfg.canvas.circular.radius) * float(cds_ratio)
    expected_anchor_clearance = circular_labels_module._effective_outer_middle_anchor_clearance_px(
        resolve_overlaps=bool(cfg.canvas.resolve_overlaps),
        strandedness=bool(cfg.canvas.strandedness),
        track_type=str(cfg.canvas.circular.track_type),
        feature_band_width_px=feature_band_width_px,
    )

    wsv134_middle_radius = math.hypot(float(wsv134_label["middle_x"]), float(wsv134_label["middle_y"]))
    wsv133_outer_radius = math.hypot(
        float(wsv133_label.get("feature_anchor_x", wsv133_label["feature_middle_x"])),
        float(wsv133_label.get("feature_anchor_y", wsv133_label["feature_middle_y"])),
    )
    radial_gap = wsv134_middle_radius - wsv133_outer_radius
    required_gap = expected_anchor_clearance + float(circular_labels_module.OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX)

    assert radial_gap >= required_gap - 1.0


def test_mjenmv_resolve_overlaps_middle_top_bottom_leader_anchor_selection() -> None:
    external_labels, total_length, _ = _load_mjenmv_external_labels_with_config(
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
        label_blacklist="",
    )
    labels = sorted((label.copy() for label in external_labels), key=lambda label: float(label["middle"]))
    assert labels

    global_collisions = circular_labels_module._count_label_leader_line_collisions(
        labels,
        total_length,
        margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
    )
    assert global_collisions == 0

    top_bottom_count = 0
    midpoint_preferred_count = 0
    fallback_count = 0

    for idx, label in enumerate(labels):
        meta = circular_labels_module._leader_start_meta(label, total_length)
        if meta is None:
            continue
        side, fixed_coord, lower, upper, _ = meta
        if side not in ("top", "bottom"):
            continue

        top_bottom_count += 1
        lower_bound = min(float(lower), float(upper))
        upper_bound = max(float(lower), float(upper))
        edge_y = float(fixed_coord)
        midpoint_x = 0.5 * (lower_bound + upper_bound)
        projected_x = min(
            max(float(label.get("middle_x", midpoint_x)), lower_bound),
            upper_bound,
        )
        chosen_x = min(
            max(float(label.get("leader_start_x", label["start_x"])), lower_bound),
            upper_bound,
        )
        chosen_y = float(label.get("leader_start_y", label["start_y"]))
        assert math.isclose(chosen_y, edge_y, abs_tol=1e-9)

        def local_collision(candidate_x: float) -> int:
            candidate = label.copy()
            candidate["leader_start_x"] = float(candidate_x)
            candidate["leader_start_y"] = edge_y
            return circular_labels_module._count_local_leader_line_collisions(
                labels,
                total_length,
                idx,
                margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
                candidate=candidate,
            )

        midpoint_collisions = local_collision(midpoint_x)
        projected_collisions = local_collision(projected_x)
        chosen_collisions = local_collision(chosen_x)
        best_collision = min(midpoint_collisions, projected_collisions, chosen_collisions)

        if midpoint_collisions == best_collision and midpoint_collisions <= projected_collisions:
            assert math.isclose(chosen_x, midpoint_x, abs_tol=1e-9)
            midpoint_preferred_count += 1

        if midpoint_collisions > min(projected_collisions, chosen_collisions):
            assert not math.isclose(chosen_x, midpoint_x, abs_tol=1e-9)
            assert chosen_collisions <= midpoint_collisions
            assert not math.isclose(chosen_x, lower_bound, abs_tol=1e-9)
            assert not math.isclose(chosen_x, upper_bound, abs_tol=1e-9)
            fallback_count += 1

    assert top_bottom_count > 0
    assert midpoint_preferred_count > 0
    assert fallback_count > 0


def test_hmmtdna_resolve_overlaps_middle_keeps_trna_lys_text_outside_feature_tracks() -> None:
    external_labels, total_length, cfg = _load_hmmtdna_external_labels_with_config(
        label_font_size=14.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )
    assert external_labels

    trna_lys = next((label for label in external_labels if label.get("label_text") == "tRNA-Lys"), None)
    assert trna_lys is not None

    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
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

    feature_intervals = circular_labels_module._build_outer_feature_radius_intervals(
        feature_dict,
        total_length,
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        cfg,
    )
    assert feature_intervals

    label_intervals = circular_labels_module._label_genome_intervals_for_clearance(trna_lys, total_length)
    assert label_intervals
    local_outer_radius = circular_labels_module._max_outer_feature_radius_for_intervals(
        feature_intervals,
        label_intervals,
        total_length,
    )

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = float(cfg.canvas.circular.track_ratio_factors[length_param][0])
    cds_ratio, _ = calculate_cds_ratio(cfg.canvas.circular.track_ratio, length_param, track_ratio_factor)
    feature_band_width_px = float(cfg.canvas.circular.radius) * float(cds_ratio)
    text_clearance_px = circular_labels_module._effective_outer_text_clearance_px(
        resolve_overlaps=bool(cfg.canvas.resolve_overlaps),
        strandedness=bool(cfg.canvas.strandedness),
        track_type=str(cfg.canvas.circular.track_type),
        feature_band_width_px=feature_band_width_px,
    )

    bbox_min_radius = circular_labels_module._label_bbox_min_radius(
        trna_lys,
        total_length,
        minimum_margin=0.0,
    )
    required_radius = (
        local_outer_radius
        + text_clearance_px
        + float(circular_labels_module.OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX)
    )
    assert bbox_min_radius >= required_radius - 1.0


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
    assert _count_overlaps_with_min_gap(external_labels, total_length) <= 1
    assert y_overlap_calls < 500000

    max_target_delta = max(_target_delta_unwrapped(label, total_length) for label in external_labels)
    assert max_target_delta <= 30.0

    max_leader_length = max(
        math.hypot(
            float(label.get("leader_start_x", label["start_x"])) - float(label["middle_x"]),
            float(label.get("leader_start_y", label["start_y"])) - float(label["middle_y"]),
        )
        for label in external_labels
    )
    assert max_leader_length <= 220.0
    assert _sum_leader_length(external_labels) <= 1764.0

    trna_lys_lengths = [
        math.hypot(
            float(label.get("leader_start_x", label["start_x"])) - float(label["middle_x"]),
            float(label.get("leader_start_y", label["start_y"])) - float(label["middle_y"]),
        )
        for label in external_labels
        if label.get("label_text") == "tRNA-Lys"
    ]
    assert trna_lys_lengths
    assert min(trna_lys_lengths) <= 145.0


def test_hmmtdna_resolve_overlaps_keeps_outer_labels_outside_feature_tracks() -> None:
    external_labels, total_length, cfg = _load_hmmtdna_external_labels_with_config(
        label_font_size=22.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )

    assert external_labels
    assert any(int(label.get("track_id", 0)) > 0 for label in external_labels)
    assert _count_overlaps(external_labels, total_length) == 0
    assert _count_overlaps_with_min_gap(external_labels, total_length) <= 1

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = cfg.canvas.circular.track_ratio_factors[length_param][0]
    cds_ratio, offset = calculate_cds_ratio(cfg.canvas.circular.track_ratio, length_param, track_ratio_factor)

    min_middle_clearance = float("inf")
    min_start_clearance = float("inf")
    for label in external_labels:
        track_id = int(label.get("track_id", 0))
        factors = calculate_feature_position_factors_circular(
            total_length,
            str(label.get("strand", "positive")),
            cfg.canvas.circular.track_ratio,
            cds_ratio,
            offset,
            cfg.canvas.circular.track_type,
            cfg.canvas.strandedness,
            track_id,
        )
        feature_outer_radius = float(cfg.canvas.circular.radius) * float(factors[2])
        middle_radius = math.hypot(float(label["middle_x"]), float(label["middle_y"]))
        start_radius = math.hypot(float(label["start_x"]), float(label["start_y"]))

        min_middle_clearance = min(min_middle_clearance, middle_radius - feature_outer_radius)
        min_start_clearance = min(min_start_clearance, start_radius - feature_outer_radius)

    assert min_middle_clearance >= circular_labels_module.MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX - 1.0
    assert min_start_clearance >= circular_labels_module.MIN_OUTER_LABEL_TEXT_CLEARANCE_PX - 1.0


def test_hmmtdna_resolve_overlaps_feature_width_keeps_middle_anchor_outside_feature_tracks() -> None:
    external_labels, total_length, cfg, feature_track_ratio_factor_override = (
        _load_hmmtdna_external_labels_with_feature_width(feature_width=75.0)
    )
    assert external_labels
    assert any(int(label.get("track_id", 0)) > 0 for label in external_labels)
    assert _count_overlaps(external_labels, total_length) == 0
    assert _count_overlaps_with_min_gap(external_labels, total_length) == 0

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    cds_ratio, offset = calculate_cds_ratio(
        cfg.canvas.circular.track_ratio,
        length_param,
        feature_track_ratio_factor_override,
    )

    min_middle_clearance = float("inf")
    for label in external_labels:
        track_id = int(label.get("track_id", 0))
        factors = calculate_feature_position_factors_circular(
            total_length,
            str(label.get("strand", "positive")),
            cfg.canvas.circular.track_ratio,
            cds_ratio,
            offset,
            cfg.canvas.circular.track_type,
            cfg.canvas.strandedness,
            track_id,
        )
        feature_outer_radius = float(cfg.canvas.circular.radius) * float(factors[2])
        middle_radius = math.hypot(float(label["middle_x"]), float(label["middle_y"]))
        min_middle_clearance = min(min_middle_clearance, middle_radius - feature_outer_radius)

    assert min_middle_clearance >= circular_labels_module.MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX - 1.0


def test_hmmtdna_resolve_overlaps_feature_width_keeps_middle_anchor_outside_local_feature_tracks() -> None:
    external_labels, total_length, cfg, feature_track_ratio_factor_override = (
        _load_hmmtdna_external_labels_with_feature_width(feature_width=75.0)
    )
    assert external_labels

    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
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

    feature_intervals = circular_labels_module._build_outer_feature_radius_intervals(
        feature_dict,
        total_length,
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        cfg,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
    )
    assert feature_intervals

    min_middle_clearance = float("inf")
    for label in external_labels:
        label_intervals = circular_labels_module._label_genome_intervals_for_clearance(label, total_length)
        assert label_intervals

        local_outer = circular_labels_module._max_outer_feature_radius_for_intervals(
            feature_intervals,
            label_intervals,
            total_length,
        )
        required_middle_radius = local_outer + circular_labels_module.MIN_OUTER_LABEL_ANCHOR_CLEARANCE_PX
        middle_radius = math.hypot(float(label["middle_x"]), float(label["middle_y"]))
        min_middle_clearance = min(min_middle_clearance, middle_radius - required_middle_radius)

    assert min_middle_clearance >= -1.0


def test_hmmtdna_resolve_overlaps_feature_width_font22_keeps_bbox_margin_from_local_feature_tracks() -> None:
    external_labels, total_length, cfg, feature_track_ratio_factor_override = (
        _load_hmmtdna_external_labels_with_feature_width(feature_width=75.0, label_font_size=22.0)
    )
    assert external_labels

    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
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

    feature_intervals = circular_labels_module._build_outer_feature_radius_intervals(
        feature_dict,
        total_length,
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        cfg,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
    )
    assert feature_intervals

    min_bbox_clearance = float("inf")
    for label in external_labels:
        bbox_clearance = _label_bbox_boundary_min_clearance_against_features(
            label,
            total_length,
            feature_intervals,
        )
        min_bbox_clearance = min(min_bbox_clearance, bbox_clearance)

    assert min_bbox_clearance >= circular_labels_module.OUTER_LABEL_FEATURE_CLEARANCE_SAFETY_PX - 1.5


def test_hmmtdna_resolve_overlaps_short_directional_features_use_center_anchor() -> None:
    external_labels, total_length, cfg = _load_hmmtdna_external_labels_with_config(
        label_font_size=14.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )
    assert external_labels

    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
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

    length_param = determine_length_parameter(total_length, cfg.labels.length_threshold.circular)
    track_ratio_factor = cfg.canvas.circular.track_ratio_factors[length_param][0]
    cds_ratio, offset = calculate_cds_ratio(cfg.canvas.circular.track_ratio, length_param, track_ratio_factor)
    circular_arrow_length_bp = calculate_circular_arrow_length(total_length)

    short_directional_feature_keys: set[tuple[float, float, float]] = set()
    non_short_feature_keys: set[tuple[float, float, float]] = set()
    for feature_object in feature_dict.values():
        if get_label_text(feature_object, label_filtering) == "":
            continue

        longest_segment_length = 0
        longest_segment_start = 0
        longest_segment_end = 0
        longest_segment_middle = 0.0
        coordinate_strand = "undefined"

        feature_location_count = 0
        for coordinate in feature_object.coordinates:
            if feature_object.location[feature_location_count].kind == "line":
                feature_location_count += 1
                continue

            coordinate_start = int(coordinate.start)
            coordinate_end = int(coordinate.end)
            coordinate_strand = str(feature_object.location[feature_location_count].strand)
            interval_length = abs(int(coordinate_end - coordinate_start) + 1)
            interval_middle = float(coordinate_end + coordinate_start) / 2.0
            feature_location_count += 1
            if interval_length > longest_segment_length:
                longest_segment_start = coordinate_start
                longest_segment_end = coordinate_end
                longest_segment_middle = interval_middle
                longest_segment_length = interval_length

        track_id = int(getattr(feature_object, "feature_track_id", 0))
        factors = calculate_feature_position_factors_circular(
            total_length,
            coordinate_strand,
            cfg.canvas.circular.track_ratio,
            cds_ratio,
            offset,
            cfg.canvas.circular.track_type,
            cfg.canvas.strandedness,
            track_id,
        )
        feature_middle_x = (float(cfg.canvas.circular.radius) * float(factors[1])) * math.cos(
            math.radians(360.0 * (longest_segment_middle / total_length) - 90)
        )
        feature_middle_y = (float(cfg.canvas.circular.radius) * float(factors[1])) * math.sin(
            math.radians(360.0 * (longest_segment_middle / total_length) - 90)
        )
        key = (round(feature_middle_x, 6), round(feature_middle_y, 6), round(longest_segment_middle, 6))
        longest_segment_bp = abs(int(longest_segment_end - longest_segment_start))
        is_short_directional = bool(getattr(feature_object, "is_directional", False)) and (
            longest_segment_bp < circular_arrow_length_bp
        )
        if is_short_directional:
            short_directional_feature_keys.add(key)
        else:
            non_short_feature_keys.add(key)

    assert short_directional_feature_keys

    short_anchors_checked = 0
    shifted_non_short_anchors = 0
    for label in external_labels:
        key = (
            round(float(label["feature_middle_x"]), 6),
            round(float(label["feature_middle_y"]), 6),
            round(float(label["middle"]), 6),
        )
        anchor_x = float(label.get("feature_anchor_x", label["feature_middle_x"]))
        anchor_y = float(label.get("feature_anchor_y", label["feature_middle_y"]))
        middle_x = float(label["feature_middle_x"])
        middle_y = float(label["feature_middle_y"])

        if key in short_directional_feature_keys:
            short_anchors_checked += 1
            assert math.isclose(anchor_x, middle_x, abs_tol=1e-6)
            assert math.isclose(anchor_y, middle_y, abs_tol=1e-6)
        elif key in non_short_feature_keys:
            if math.hypot(anchor_x - middle_x, anchor_y - middle_y) > 0.5:
                shifted_non_short_anchors += 1

    assert short_anchors_checked > 0
    assert shifted_non_short_anchors > 0


def test_hmmtdna_resolve_overlaps_default_font_keeps_label_bboxes_outside_local_feature_tracks() -> None:
    external_labels, total_length, cfg = _load_hmmtdna_external_labels_with_config(
        label_font_size=14.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )
    assert external_labels
    assert _count_overlaps(external_labels, total_length) == 0
    assert _count_overlaps_with_min_gap(external_labels, total_length) <= 1

    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
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

    feature_intervals = circular_labels_module._build_outer_feature_radius_intervals(
        feature_dict,
        total_length,
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        cfg,
    )
    assert feature_intervals

    min_bbox_clearance = float("inf")
    for label in external_labels:
        sampled_positions = circular_labels_module._sample_label_genome_positions(label, total_length)
        assert sampled_positions

        required_feature_outer = max(
            circular_labels_module._max_outer_feature_radius_at_position(feature_intervals, position, total_length)
            for position in sampled_positions
        )
        required_radius = required_feature_outer + circular_labels_module.MIN_OUTER_LABEL_TEXT_CLEARANCE_PX
        bbox_clearance = _label_bbox_min_radius(label, total_length) - required_radius
        min_bbox_clearance = min(min_bbox_clearance, bbox_clearance)

    assert min_bbox_clearance >= -1.0


def test_hmmtdna_resolve_overlaps_middle_has_no_label_leader_line_collisions() -> None:
    external_labels, total_length, _ = _load_hmmtdna_external_labels_with_config(
        label_font_size=14.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )
    assert external_labels

    collision_count = circular_labels_module._count_label_leader_line_collisions(
        external_labels,
        total_length,
        margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
    )
    assert collision_count == 0


def test_resolve_label_leader_line_collisions_reduces_fixture_collisions() -> None:
    labels, total_length = _make_label_leader_collision_fixture()
    labels = [label.copy() for label in labels]
    labels = circular_labels_module.assign_leader_start_points(labels, total_length)

    before_collisions = circular_labels_module._count_label_leader_line_collisions(
        labels,
        total_length,
        margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
    )
    before_overlaps = _count_overlaps_with_min_gap(labels, total_length)
    assert before_collisions > 0

    optimized = circular_labels_module._resolve_label_leader_line_collisions(
        labels,
        total_length,
        margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
        step_deg=circular_labels_module.LEADER_LABEL_SHIFT_STEP_DEG,
        max_shift_deg=circular_labels_module.LEADER_LABEL_MAX_SHIFT_DEG,
        max_passes=circular_labels_module.LEADER_LABEL_MAX_PASSES,
        min_order_gap_deg=circular_labels_module.LEADER_LABEL_MIN_ORDER_GAP_DEG,
    )

    after_collisions = circular_labels_module._count_label_leader_line_collisions(
        optimized,
        total_length,
        margin_px=circular_labels_module.LEADER_LABEL_COLLISION_MARGIN_PX,
    )
    after_overlaps = _count_overlaps_with_min_gap(optimized, total_length)
    unwrapped_angles = circular_labels_module._derive_monotonic_unwrapped_angles(optimized, total_length)

    assert after_collisions < before_collisions
    assert after_collisions == 0
    assert after_overlaps <= before_overlaps
    assert all(unwrapped_angles[i] < unwrapped_angles[i + 1] for i in range(len(unwrapped_angles) - 1))


def test_leader_start_lies_on_bbox_perimeter_without_corner_snap() -> None:
    external_labels, total_length = _load_hmmtdna_external_labels(label_font_size=22.0)
    assert external_labels

    for label in external_labels:
        assert "leader_start_x" in label
        assert "leader_start_y" in label

        leader_x = float(label["leader_start_x"])
        leader_y = float(label["leader_start_y"])
        min_x, max_x = circular_labels_module._label_x_bounds(label, minimum_margin=0.0)
        min_y, max_y = circular_labels_module._label_y_bounds(label, total_length, minimum_margin=0.0)

        assert min_x - 1e-9 <= leader_x <= max_x + 1e-9
        assert min_y - 1e-9 <= leader_y <= max_y + 1e-9

        on_vertical_edge = math.isclose(leader_x, min_x, abs_tol=1.0) or math.isclose(leader_x, max_x, abs_tol=1.0)
        on_horizontal_edge = math.isclose(leader_y, min_y, abs_tol=1.0) or math.isclose(leader_y, max_y, abs_tol=1.0)
        assert on_vertical_edge or on_horizontal_edge


def test_hmmtdna_resolve_overlaps_recomputes_leader_start_after_label_shifts() -> None:
    external_labels, total_length, _ = _load_hmmtdna_external_labels_with_config(
        label_font_size=14.0,
        strandedness=False,
        resolve_overlaps=True,
        track_type="middle",
    )
    assert external_labels

    for label in external_labels:
        assert "leader_start_x" in label
        assert "leader_start_y" in label

        leader_x = float(label["leader_start_x"])
        leader_y = float(label["leader_start_y"])
        min_x, max_x = circular_labels_module._label_x_bounds(label, minimum_margin=0.0)
        min_y, max_y = circular_labels_module._label_y_bounds(label, total_length, minimum_margin=0.0)

        assert min_x - 1e-9 <= leader_x <= max_x + 1e-9
        assert min_y - 1e-9 <= leader_y <= max_y + 1e-9

        on_vertical_edge = math.isclose(leader_x, min_x, abs_tol=1.0) or math.isclose(leader_x, max_x, abs_tol=1.0)
        on_horizontal_edge = math.isclose(leader_y, min_y, abs_tol=1.0) or math.isclose(leader_y, max_y, abs_tol=1.0)
        assert on_vertical_edge or on_horizontal_edge


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
    original_radius = math.hypot(labels[0]["start_x"], labels[0]["start_y"])

    circular_assemble_module._resolve_label_legend_collisions(labels, total_length, canvas_config, legend_config)

    assert not circular_assemble_module._labels_collide_with_legend(labels, total_length, canvas_config, legend_config)
    assert canvas_config.legend_offset_x == original_legend_x
    assert canvas_config.legend_offset_y == original_legend_y
    shifted_radius = math.hypot(labels[0]["start_x"], labels[0]["start_y"])
    assert math.isclose(shifted_radius, original_radius, abs_tol=0.5)
    assert abs(labels[0]["start_y"]) > 0.1


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


def test_label_legend_collision_keeps_feature_order_by_moving_neighbor_block() -> None:
    labels, total_length, canvas_config, legend_config = _make_order_sensitive_legend_collision_fixture()
    initial_low_y = labels[0]["start_y"]
    initial_high_y = labels[1]["start_y"]

    assert circular_assemble_module._legend_collision_indices(labels, total_length, canvas_config, legend_config) == [1]

    circular_assemble_module._resolve_label_legend_collisions(labels, total_length, canvas_config, legend_config)

    assert not circular_assemble_module._labels_collide_with_legend(labels, total_length, canvas_config, legend_config)

    low_unwrapped = _label_unwrapped_angle_for_order_test(labels[0], total_length)
    high_unwrapped = _label_unwrapped_angle_for_order_test(labels[1], total_length)
    assert low_unwrapped < high_unwrapped

    # Neighbor label should move together to preserve order during legend avoidance.
    assert abs(labels[0]["start_y"] - initial_low_y) > 0.1
    assert abs(labels[1]["start_y"] - initial_high_y) > 0.1


def test_expand_canvas_to_fit_external_labels_keeps_all_labels_inside() -> None:
    total_length = 4000
    labels = [
        {
            "middle": 1000,
            "start_x": 470.0,
            "start_y": 0.0,
            "width_px": 120.0,
            "height_px": 20.0,
            "is_inner": False,
            "is_embedded": False,
        }
    ]
    canvas_config = SimpleNamespace(
        total_width=1000.0,
        total_height=1000.0,
        offset_x=500.0,
        offset_y=500.0,
        legend_offset_x=0.0,
        legend_offset_y=0.0,
    )

    expanded = circular_assemble_module._expand_canvas_to_fit_external_labels(labels, total_length, canvas_config)
    assert expanded is True

    bounds = circular_assemble_module._external_label_bounds_on_canvas(labels, total_length, canvas_config)
    assert bounds is not None
    pad = circular_assemble_module.LABEL_CANVAS_PADDING_PX
    assert bounds[0] >= pad - 1e-6
    assert bounds[1] >= pad - 1e-6
    assert bounds[2] <= float(canvas_config.total_width) - pad + 1e-6
    assert bounds[3] <= float(canvas_config.total_height) - pad + 1e-6
