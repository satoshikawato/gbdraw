from __future__ import annotations

import math
from pathlib import Path

import pytest
from Bio import SeqIO
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
from gbdraw.labels.linear import prepare_label_list_linear
from gbdraw.render.drawers.linear.labels import LabelDrawer


def _prepare_linear_labels(
    *,
    label_placement: str,
    label_rotation: float = 0.0,
    label_font_size: float = 14.0,
    separate_strands: bool = False,
    input_filename: str = "MjeNMV.gb",
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
        label_placement=label_placement,
        label_rotation=label_rotation,
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
        config_dict,
        cfg=cfg,
    )


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
        assert float(label["middle_x"]) == pytest.approx(feature_mid_x)
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
