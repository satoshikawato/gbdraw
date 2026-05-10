from __future__ import annotations

import math
from pathlib import Path

import pytest
from Bio import SeqIO
from pandas import DataFrame
from svgwrite.container import Group

from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import FeatureDrawingConfigurator
from gbdraw.diagrams.linear.precalc import _precalculate_label_dimensions
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering
from gbdraw.labels.linear import (
    _find_lowest_available_track_indexed,
    _insert_external_label_index,
    calculate_label_bounds,
    find_lowest_available_track,
    prepare_label_list_linear,
)
from gbdraw.render.drawers.linear.labels import LabelDrawer


def _prepare_linear_labels(
    *,
    label_placement: str,
    label_rotation: float = 0.0,
    label_font_size: float = 14.0,
    linear_label_spacing: float | None = None,
    separate_strands: bool = False,
    track_layout: str = "middle",
    input_filename: str = "MjeNMV.gb",
    qualifier_priority: tuple[str, str] | None = None,
) -> list[dict]:
    input_path = Path(__file__).parent / "test_inputs" / input_filename
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels="all",
        strandedness=separate_strands,
        label_blacklist="",
        label_font_size=label_font_size,
        linear_label_spacing=linear_label_spacing,
        label_placement=label_placement,
        label_rotation=label_rotation,
        linear_track_layout=track_layout,
    )
    if qualifier_priority is not None:
        feature_type, priorities = qualifier_priority
        config_dict["labels"]["filtering"]["qualifier_priority_df"] = DataFrame(
            [{"feature_type": feature_type, "priorities": priorities}]
        )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )

    default_colors = load_default_colors("", "default")
    color_table = None
    color_table, default_colors = preprocess_color_tables(color_table, default_colors)
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        cfg.objects.features.features_drawn,
        default_colors,
        canvas_cfg.strandedness,
        canvas_cfg.resolve_overlaps,
        label_filtering,
    )

    return prepare_label_list_linear(
        feature_dict,
        len(record.seq),
        canvas_cfg.alignment_width,
        1.0,
        canvas_cfg.cds_height,
        canvas_cfg.strandedness,
        canvas_cfg.track_layout,
        canvas_cfg.track_axis_gap,
        config_dict,
        cfg=cfg,
    )


def _prepare_linear_labels_for_records(
    *,
    input_filenames: tuple[str, ...],
    show_labels: str = "first",
    label_placement: str,
    label_rotation: float = 0.0,
    label_font_size: float = 14.0,
    separate_strands: bool = False,
    track_layout: str = "middle",
    qualifier_priority: tuple[str, str] | None = None,
) -> list[dict]:
    input_dir = Path(__file__).parent / "test_inputs"
    records = [SeqIO.read(str(input_dir / filename), "genbank") for filename in input_filenames]

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=show_labels,
        strandedness=separate_strands,
        label_blacklist="",
        label_font_size=label_font_size,
        label_placement=label_placement,
        label_rotation=label_rotation,
        linear_track_layout=track_layout,
    )
    if qualifier_priority is not None:
        feature_type, priorities = qualifier_priority
        config_dict["labels"]["filtering"]["qualifier_priority_df"] = DataFrame(
            [{"feature_type": feature_type, "priorities": priorities}]
        )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=max(len(record.seq) for record in records),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        config_dict=config_dict,
        canvas_config=canvas_cfg,
        cfg=cfg,
    )

    _, all_labels_by_record, _ = _precalculate_label_dimensions(
        records,
        feature_cfg,
        canvas_cfg,
        config_dict,
        cfg=cfg,
    )
    return all_labels_by_record[records[0].id]


def _rotated_y_bounds_from_anchor(label: dict) -> tuple[float, float]:
    radians = math.radians(float(label["rotation_deg"]))
    width_px = float(label["width_px"])
    height_px = float(label["height_px"])
    text_anchor = str(label.get("text_anchor", "middle"))
    if text_anchor == "start":
        x_values = (0.0, width_px)
    elif text_anchor == "end":
        x_values = (-width_px, 0.0)
    else:
        x_values = (-width_px / 2.0, width_px / 2.0)
    half_height = height_px / 2.0
    y_values = (-half_height, half_height)
    y_offsets = [(x * math.sin(radians)) + (y * math.cos(radians)) for x in x_values for y in y_values]
    return min(y_offsets), max(y_offsets)


def _feature_adjacent_label_y(label: dict) -> float:
    label_vertical_gap = max(1.0, float(label["height_px"]) * 0.05)
    if bool(label.get("above_feature_place_above", True)):
        return float(label["feature_top_y"]) - label_vertical_gap - float(label["label_contact_y_offset"])
    return float(label["feature_bottom_y"]) + label_vertical_gap - float(label["label_contact_y_offset"])


def _label_bounds_overlap(label1: dict, label2: dict, min_gap_px: float = 0.0) -> bool:
    left1, right1, top1, bottom1 = calculate_label_bounds(label1)
    left2, right2, top2, bottom2 = calculate_label_bounds(label2)
    return not (
        right1 + min_gap_px <= left2
        or right2 + min_gap_px <= left1
        or bottom1 + min_gap_px <= top2
        or bottom2 + min_gap_px <= top1
    )


def _assert_no_label_bounds_overlap(labels: list[dict]) -> None:
    placed = []
    for label in sorted(labels, key=lambda item: float(item["feature_anchor_x"])):
        for other in placed:
            assert not _label_bounds_overlap(label, other)
        placed.append(label)


@pytest.mark.linear
def test_linear_above_feature_placement_embeds_all_labels() -> None:
    auto_labels = _prepare_linear_labels(label_placement="auto")
    above_feature_labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
    )

    assert auto_labels
    assert above_feature_labels
    assert any(not label["is_embedded"] for label in auto_labels)
    assert all(label["is_embedded"] for label in above_feature_labels)
    assert all(label["track_id"] == "track_0" for label in above_feature_labels)
    assert all(label["rotation_deg"] == pytest.approx(-45.0) for label in above_feature_labels)
    assert all(label["text_anchor"] == "start" for label in above_feature_labels)
    for label in above_feature_labels:
        _, y_max_offset = _rotated_y_bounds_from_anchor(label)
        rotated_bottom_y = float(label["middle_y"]) + y_max_offset
        assert rotated_bottom_y < float(label["feature_top_y"])
    assert {label["label_text"] for label in auto_labels} == {label["label_text"] for label in above_feature_labels}


@pytest.mark.linear
def test_linear_auto_ignores_label_rotation() -> None:
    labels = _prepare_linear_labels(
        label_placement="auto",
        label_rotation=45.0,
    )
    assert labels
    assert all(label["rotation_deg"] == pytest.approx(0.0) for label in labels)
    assert all(label["text_anchor"] == "middle" for label in labels)


@pytest.mark.linear
def test_linear_above_feature_anchor_is_feature_midpoint() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=30.0,
    )
    assert labels

    for label in labels:
        feature_mid_x = (float(label["feature_start_x"]) + float(label["feature_end_x"])) / 2.0
        label_contact_x = float(label["middle_x"]) + float(label["label_contact_x_offset"])
        assert label_contact_x == pytest.approx(feature_mid_x)
        assert label["text_anchor"] == "start"


@pytest.mark.linear
def test_linear_above_feature_negative_strand_is_below_and_rotates_opposite() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        separate_strands=True,
        input_filename="MG1655.gbk",
    )
    assert labels

    negative_labels = [label for label in labels if label["strand"] == "negative"]
    positive_labels = [label for label in labels if label["strand"] == "positive"]
    assert negative_labels
    assert positive_labels

    for label in positive_labels:
        assert label["rotation_deg"] == pytest.approx(-45.0)
        assert label["text_anchor"] == "start"
        _, y_max_offset = _rotated_y_bounds_from_anchor(label)
        rotated_bottom_y = float(label["middle_y"]) + y_max_offset
        assert rotated_bottom_y < float(label["feature_top_y"])

    for label in negative_labels:
        assert label["rotation_deg"] == pytest.approx(45.0)
        assert label["text_anchor"] == "start"
        y_min_offset, _ = _rotated_y_bounds_from_anchor(label)
        rotated_top_y = float(label["middle_y"]) + y_min_offset
        assert rotated_top_y > float(label["feature_bottom_y"])


@pytest.mark.linear
def test_linear_label_drawer_applies_rotation_transform() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        separate_strands=True,
        input_filename="MG1655.gbk",
    )
    assert labels

    # Verify both signs are carried through into SVG transforms.
    positive_label = next(label for label in labels if label["strand"] == "positive")
    negative_label = next(label for label in labels if label["strand"] == "negative")

    group = Group(id="test")
    drawer = LabelDrawer(config_dict={})
    group = drawer.draw(positive_label, group)
    group = drawer.draw(negative_label, group)

    positive_text = group.elements[-2]
    negative_text = group.elements[-1]
    positive_transform = str(positive_text.attribs.get("transform", ""))
    negative_transform = str(negative_text.attribs.get("transform", ""))
    positive_anchor = str(positive_text.attribs.get("text-anchor", ""))
    negative_anchor = str(negative_text.attribs.get("text-anchor", ""))

    assert "rotate(-45" in positive_transform
    assert "rotate(45" in negative_transform
    assert positive_anchor == "start"
    assert negative_anchor == "start"


@pytest.mark.linear
def test_linear_above_feature_rotated_labels_do_not_overlap() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
    )
    assert labels

    _assert_no_label_bounds_overlap(labels)


@pytest.mark.linear
def test_linear_above_feature_bgc0000708_gene_labels_do_not_overlap() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        input_filename="BGC0000708.gbk",
        qualifier_priority=("CDS", "gene"),
    )
    assert labels
    assert any(label.get("leader_line") for label in labels)
    assert any(not label.get("leader_line") for label in labels)
    _assert_no_label_bounds_overlap(labels)


@pytest.mark.linear
def test_linear_above_feature_bgc_multi_record_gene_labels_compact_without_overlap() -> None:
    labels = _prepare_linear_labels_for_records(
        input_filenames=(
            "BGC0000708.gbk",
            "BGC0000709.gbk",
            "BGC0000712.gbk",
            "BGC0000713.gbk",
        ),
        show_labels="first",
        label_placement="above_feature",
        label_rotation=45.0,
        qualifier_priority=("CDS", "gene"),
    )
    assert labels
    _assert_no_label_bounds_overlap(labels)

    collision_gap = max(1.0, max(float(label["height_px"]) for label in labels) * 0.05)
    rotated_heights = [
        calculate_label_bounds(label)[3] - calculate_label_bounds(label)[2]
        for label in labels
        if bool(label.get("above_feature_place_above", True))
    ]
    old_global_step = max(float(collision_gap), max(rotated_heights) + float(collision_gap), 4.0)

    shifted_livx = next(label for label in labels if label["label_text"] == "livX" and label.get("leader_line"))
    livx_original_y = _feature_adjacent_label_y(shifted_livx)
    livx_shift = livx_original_y - float(shifted_livx["middle_y"])

    assert livx_shift > 0.0
    assert livx_shift < old_global_step

    livy = next(label for label in labels if label["label_text"] == "livY")
    assert float(livy["middle_y"]) == pytest.approx(_feature_adjacent_label_y(livy))
    for neighbor in [label for label in labels if label["label_text"] in {"livW", "livZ"}]:
        assert not _label_bounds_overlap(livy, neighbor)


@pytest.mark.linear
def test_linear_above_feature_shifted_labels_get_leader_lines() -> None:
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
    )
    shifted = [label for label in labels if label.get("leader_line")]
    assert shifted
    for label in shifted:
        assert float(label["leader_start_x"]) == pytest.approx(float(label["feature_anchor_x"]))
        assert float(label["leader_end_x"]) == pytest.approx(float(label["feature_anchor_x"]))
        assert abs(float(label["leader_end_y"]) - float(label["leader_start_y"])) > 0.0


@pytest.mark.linear
def test_linear_precalc_includes_above_feature_rotated_embedded_labels() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gb"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels="all",
        strandedness=False,
        label_blacklist="",
        label_font_size=14.0,
        label_placement="above_feature",
        label_rotation=45.0,
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        config_dict=config_dict,
        canvas_config=canvas_cfg,
        cfg=cfg,
    )

    required_label_height, _, record_label_heights = _precalculate_label_dimensions(
        [record],
        feature_cfg,
        canvas_cfg,
        config_dict,
        cfg=cfg,
    )

    assert required_label_height > 0
    assert required_label_height > canvas_cfg.original_vertical_offset
    assert record_label_heights[record.id] == pytest.approx(required_label_height)


@pytest.mark.linear
def test_linear_canvas_alignment_width_starts_at_figure_width() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )

    assert canvas_cfg.alignment_width == pytest.approx(canvas_cfg.fig_width)


@pytest.mark.linear
def test_linear_precalculated_labels_match_final_alignment_width_positions() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels="all",
        strandedness=False,
        label_blacklist="",
        linear_track_layout="middle",
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        config_dict=config_dict,
        canvas_config=canvas_cfg,
        cfg=cfg,
    )

    _, all_labels, _ = _precalculate_label_dimensions(
        [record],
        feature_cfg,
        canvas_cfg,
        config_dict,
        cfg=cfg,
    )
    precalculated_tuples = sorted(
        (
            float(label["feature_start_x"]),
            float(label["feature_end_x"]),
            float(label["middle_x"]),
        )
        for label in all_labels.get(record.id, [])
        if label["label_text"] == "tRNA-Phe"
    )

    color_table, default_colors = preprocess_color_tables(
        feature_cfg.color_table,
        feature_cfg.default_colors,
    )
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        feature_cfg.selected_features_set,
        default_colors,
        canvas_cfg.strandedness,
        canvas_cfg.resolve_overlaps,
        label_filtering,
        directional_feature_types=feature_cfg.directional_feature_types,
        feature_visibility_rules=feature_cfg.feature_visibility_rules,
    )

    canvas_cfg.alignment_width = canvas_cfg.fig_width
    expected_tuples = sorted(
        (
            float(label["feature_start_x"]),
            float(label["feature_end_x"]),
            float(label["middle_x"]),
        )
        for label in prepare_label_list_linear(
            feature_dict,
            len(record.seq),
            canvas_cfg.alignment_width,
            1.0,
            canvas_cfg.cds_height,
            canvas_cfg.strandedness,
            canvas_cfg.track_layout,
            canvas_cfg.track_axis_gap,
            config_dict,
            cfg=cfg,
        )
        if label["label_text"] == "tRNA-Phe"
    )

    assert precalculated_tuples
    assert len(precalculated_tuples) == len(expected_tuples)
    assert precalculated_tuples == pytest.approx(expected_tuples)


@pytest.mark.linear
def test_linear_auto_labels_in_below_layout_are_placed_below_features() -> None:
    labels = _prepare_linear_labels(
        label_placement="auto",
        track_layout="below",
    )
    assert labels
    external_labels = [label for label in labels if not label["is_embedded"]]
    assert external_labels
    for label in external_labels:
        assert float(label["middle_y"]) > float(label["feature_bottom_y"])


@pytest.mark.linear
def test_label_spacing_config_defaults_and_overrides() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    assert cfg.labels.spacing.circular == pytest.approx(3.0)
    assert cfg.labels.spacing.linear == pytest.approx(3.0)

    modified = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        circular_label_spacing=7.5,
        linear_label_spacing=9.0,
    )
    modified_cfg = GbdrawConfig.from_dict(modified)
    assert modified_cfg.labels.spacing.circular == pytest.approx(7.5)
    assert modified_cfg.labels.spacing.linear == pytest.approx(9.0)


@pytest.mark.linear
def test_linear_external_label_spacing_override_increases_track_gap() -> None:
    default_labels = _prepare_linear_labels(
        label_placement="auto",
        track_layout="below",
        linear_label_spacing=3.0,
    )
    wider_labels = _prepare_linear_labels(
        label_placement="auto",
        track_layout="below",
        linear_label_spacing=12.0,
    )

    default_track_positions = sorted(
        {float(label["middle_y"]) for label in default_labels if not label["is_embedded"]}
    )
    wider_track_positions = sorted(
        {float(label["middle_y"]) for label in wider_labels if not label["is_embedded"]}
    )

    assert len(default_track_positions) >= 2
    assert len(default_track_positions) == len(wider_track_positions)
    assert (default_track_positions[1] - default_track_positions[0]) < (
        wider_track_positions[1] - wider_track_positions[0]
    )


@pytest.mark.linear
def test_linear_external_label_track_index_matches_legacy_scan() -> None:
    labels = [
        {"start": 0.0, "end": 50.0},
        {"start": 20.0, "end": 70.0},
        {"start": 80.0, "end": 110.0},
        {"start": 45.0, "end": 90.0},
        {"start": 112.0, "end": 140.0},
    ]
    legacy_tracks: dict[str, list[dict]] = {}
    indexed_tracks: dict[str, list[dict]] = {}
    indexed_track_indexes = {}
    indexed_label_by_id: dict[int, dict] = {}
    legacy_assignments: list[int] = []
    indexed_assignments: list[int] = []

    for label in labels:
        legacy_track = find_lowest_available_track(legacy_tracks, label)
        legacy_assignments.append(legacy_track)
        legacy_tracks.setdefault(f"track_{legacy_track}", []).append(label)

        indexed_track = _find_lowest_available_track_indexed(
            indexed_tracks,
            indexed_track_indexes,
            indexed_label_by_id,
            label,
            bucket_size=16.0,
        )
        indexed_assignments.append(indexed_track)
        indexed_tracks.setdefault(f"track_{indexed_track}", []).append(label)
        label_id = len(indexed_label_by_id)
        indexed_label_by_id[label_id] = label
        _insert_external_label_index(
            indexed_track_indexes,
            f"track_{indexed_track}",
            label_id,
            label,
            bucket_size=16.0,
        )

    assert indexed_assignments == legacy_assignments


@pytest.mark.linear
def test_linear_above_feature_spacing_override_does_not_move_embedded_labels() -> None:
    default_labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        linear_label_spacing=3.0,
    )
    wider_labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        linear_label_spacing=12.0,
    )

    default_positions = sorted((label["label_text"], float(label["middle_y"])) for label in default_labels)
    wider_positions = sorted((label["label_text"], float(label["middle_y"])) for label in wider_labels)

    assert [item[0] for item in default_positions] == [item[0] for item in wider_positions]
    assert [item[1] for item in default_positions] == pytest.approx([item[1] for item in wider_positions])
