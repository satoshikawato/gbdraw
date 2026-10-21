from __future__ import annotations

import math
import re
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Callable, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import ParamSpec, TypeVar

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite.container import Group

import gbdraw.labels.linear as linear_labels_module
import gbdraw.linear as linear_cli_module
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import FeatureDrawingConfigurator
from gbdraw.diagrams.linear.orthogroup_alignment import (
    OrthogroupLabelEligibility,
    build_orthogroup_label_eligibility,
)
from gbdraw.diagrams.linear.precalc import _precalculate_label_dimensions
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.labels.label_style import ROW_LABEL_GAP_EM
from gbdraw.features.ids import compute_feature_hash
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering
from gbdraw.layout.spatial import Aabb, AabbIndex, IntervalIndex, candidate_aabb_pairs
from gbdraw.labels.linear import (
    _find_lowest_available_track_indexed,
    _insert_external_label_index,
    _label_corners,
    _row_label_aabb,
    _rotated_label_rects_intersect,
    calculate_label_bounds,
    LinearLabelFrame,
    find_lowest_available_track,
    linear_label_fit,
    linear_record_labels,
    place_linear_label_row,
)
from gbdraw.core.text import TEXT_ANCHOR_SHARE, calculate_bbox_dimensions, calculate_svg_bbox_dimensions, text_ink_offsets
from gbdraw.exceptions import LabelLeadersCrossWarning
from gbdraw.render.drawers.linear.labels import LabelDrawer
from tests.utils.linear_render_context import alone_sequence_widths, make_linear_record_render_context

_P = ParamSpec("_P")
_R = TypeVar("_R")


def _synthetic_label_record(record_id: str, length: int, label: str, extra_label: str | None = None) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, name=record_id, description=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(100, 300, strand=1),
            type="CDS",
            qualifiers={"product": [label]},
        )
    ]
    if extra_label is not None:
        record.features.append(
            SeqFeature(
                FeatureLocation(500, 700, strand=1),
                type="CDS",
                qualifiers={"product": [extra_label]},
            )
        )
    return record


def _orthogroup_member(record_index: int, feature_svg_id: str | None) -> SimpleNamespace:
    return SimpleNamespace(
        record_index=record_index,
        feature_svg_id=feature_svg_id,
        start=0,
        end=1,
        protein_id=str(feature_svg_id or ""),
        source_protein_id=str(feature_svg_id or ""),
        representative=False,
    )


def _comparison_row(**metadata: object) -> dict[str, object]:
    row = dict.fromkeys(COMPARISON_COLUMNS, 0)
    row.update(
        {
            "query": "query",
            "subject": "subject",
            "identity": 100.0,
            "alignment_length": 100,
            # In the search frame on the first CDS of each synthetic record (PD-OI-073).
            "qstart": 101,
            "qend": 300,
            "sstart": 101,
            "send": 300,
            "evalue": 0.0,
            "bitscore": 100.0,
        }
    )
    row.update(metadata)
    return row


def _extract_label_font_sizes(svg_content: str, label_texts: set[str]) -> dict[str, float]:
    root = ET.fromstring(svg_content)
    sizes: dict[str, float] = {}
    for element in root.iter():
        tag = str(element.tag).rsplit("}", 1)[-1]
        if tag not in {"text", "textPath"}:
            continue
        text = "".join(element.itertext()).strip()
        if text not in label_texts:
            continue
        font_size = element.attrib.get("font-size")
        if font_size is None:
            continue
        sizes[text] = float(str(font_size).replace("px", "").replace("pt", ""))
    return sizes


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
    if not input_path.exists():
        input_path = Path(__file__).parent.parent / "examples" / input_filename
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        {
            "labels.linear.scope": 'all',
            "canvas.strandedness": separate_strands,
            "labels.filtering.blacklist_keywords": [],
            "labels.font_size.short": label_font_size,
            "labels.font_size.long": label_font_size,
            "labels.font_size.linear.short": label_font_size,
            "labels.font_size.linear.long": label_font_size,
            "labels.linear.placement": label_placement,
            "labels.linear.rotation": label_rotation,
            "canvas.linear.track_layout": track_layout,
            **({'labels.spacing.linear': linear_label_spacing} if linear_label_spacing is not None else {}),
        },
    )
    if qualifier_priority is not None:
        feature_type, priorities = qualifier_priority
        config_dict["labels"]["filtering"]["qualifier_priority_df"] = DataFrame(
            [{"feature_type": feature_type, "priorities": priorities}]
        )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        profile=profile,
        legend="none",
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

    return _labels_alone(feature_dict, len(record.seq), canvas_cfg, profile)


def _labels_alone(
    feature_dict: dict, record_length: int, canvas_cfg: LinearCanvasConfigurator, profile: LinearRenderProfile
) -> list[dict]:
    """One record's feature labels placed alone in its row."""

    record_labels = linear_record_labels(
        feature_dict,
        record_length,
        canvas_cfg.alignment_width,
        canvas_cfg.cds_height,
        canvas_cfg.strandedness,
        canvas_cfg.track_layout,
        canvas_cfg.track_axis_gap,
        profile,
    )
    return place_linear_label_row([(record_labels, 0.0)], profile=profile, track_layout=canvas_cfg.track_layout)[0]


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
    all_records: bool = False,
    label_frames: list[LinearLabelFrame] | None = None,
) -> list:
    input_dir = Path(__file__).parent / "test_inputs"
    records = [SeqIO.read(str(input_dir / filename), "genbank") for filename in input_filenames]

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        {
            "labels.linear.scope": show_labels,
            "canvas.strandedness": separate_strands,
            "labels.filtering.blacklist_keywords": [],
            "labels.font_size.short": label_font_size,
            "labels.font_size.long": label_font_size,
            "labels.font_size.linear.short": label_font_size,
            "labels.font_size.linear.long": label_font_size,
            "labels.linear.placement": label_placement,
            "labels.linear.rotation": label_rotation,
            "canvas.linear.track_layout": track_layout,
        },
    )
    if qualifier_priority is not None:
        feature_type, priorities = qualifier_priority
        config_dict["labels"]["filtering"]["qualifier_priority_df"] = DataFrame(
            [{"feature_type": feature_type, "priorities": priorities}]
        )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=max(len(record.seq) for record in records),
        profile=profile,
        legend="none",
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        profile=profile,
        canvas_config=canvas_cfg,
    )

    _, all_labels_by_record, _ = _precalculate_label_dimensions(
        records,
        feature_cfg,
        canvas_cfg,
        render_context=make_linear_record_render_context(profile),
        sequence_widths=alone_sequence_widths(records, canvas_cfg),
        label_frames=label_frames,
    )
    return all_labels_by_record if all_records else all_labels_by_record[0]


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


@pytest.mark.linear
@pytest.mark.parametrize(
    ("text_anchor", "expected"),
    (
        ("start", (100.0, 120.0, 45.0, 55.0)),
        ("middle", (90.0, 110.0, 45.0, 55.0)),
        ("end", (80.0, 100.0, 45.0, 55.0)),
    ),
)
def test_unrotated_label_bounds_skip_rotated_corner_geometry(
    monkeypatch: pytest.MonkeyPatch,
    text_anchor: str,
    expected: tuple[float, float, float, float],
) -> None:
    def unexpected_rotated_geometry(*_args, **_kwargs):
        raise AssertionError("unrotated labels should use the exact axis-aligned path")

    monkeypatch.setattr(
        linear_labels_module,
        "_rotated_bounds_from_anchor",
        unexpected_rotated_geometry,
    )

    bounds = calculate_label_bounds(
        {
            "width_px": 20.0,
            "height_px": 10.0,
            "rotation_deg": 0.0,
            "text_anchor": text_anchor,
            "middle_x": 100.0,
            "middle_y": 50.0,
        }
    )

    assert bounds == expected


@pytest.mark.linear
def test_rotated_label_bounds_retain_oriented_geometry() -> None:
    bounds = calculate_label_bounds(
        {
            "width_px": 20.0,
            "height_px": 10.0,
            "rotation_deg": 90.0,
            "text_anchor": "start",
            "middle_x": 100.0,
            "middle_y": 50.0,
        }
    )

    assert bounds == pytest.approx((95.0, 105.0, 50.0, 70.0))


def _assert_no_label_bounds_overlap(labels: list[dict], min_gap_px: float = 0.0) -> None:
    placed = []
    for label in sorted(labels, key=lambda item: float(item["feature_anchor_x"])):
        for other in placed:
            if min_gap_px:
                assert not _label_bounds_overlap(label, other, min_gap_px=min_gap_px)
            else:
                assert not _rotated_label_rects_intersect(label, other)
        placed.append(label)


@pytest.mark.linear
def test_orthogroup_label_eligibility_tracks_members_and_top_record() -> None:
    orthogroups = SimpleNamespace(
        orthogroups={
            "og_1": [
                _orthogroup_member(1, "feature_r1"),
                _orthogroup_member(0, "feature_r0"),
            ],
            "og_2": [
                _orthogroup_member(2, "feature_r2"),
                _orthogroup_member(1, "feature_r1_a"),
                _orthogroup_member(1, "feature_r1_b"),
            ],
            "og_3": [
                _orthogroup_member(0, None),
                _orthogroup_member(1, ""),
            ],
        }
    )

    assert build_orthogroup_label_eligibility(None) == OrthogroupLabelEligibility({}, {})
    eligibility = build_orthogroup_label_eligibility(orthogroups)
    assert eligibility.member_ids_by_record == {
        0: {"feature_r0"},
        1: {"feature_r1", "feature_r1_a", "feature_r1_b"},
        2: {"feature_r2"},
    }
    assert eligibility.top_member_ids_by_record == {
        0: {"feature_r0"},
        1: {"feature_r1_a", "feature_r1_b"},
    }
    comparisons = [
        DataFrame(
            [
                _comparison_row(
                    orthogroup_id="og_1",
                    query_record_index=1,
                    query_feature_svg_id="feature_r1",
                    subject_record_index=0,
                    subject_feature_svg_id="feature_r0",
                ),
                _comparison_row(
                    orthogroup_id="og_2",
                    query_record_index=2,
                    query_feature_svg_id="feature_r2",
                    subject_record_index=1,
                    subject_feature_svg_id="feature_r1_a",
                ),
                _comparison_row(
                    orthogroup_id="og_2",
                    query_record_index=1,
                    query_feature_svg_id="feature_r1_b",
                    subject_record_index=2,
                    subject_feature_svg_id="feature_r2_b",
                ),
            ]
        )
    ]

    comparison_eligibility = build_orthogroup_label_eligibility(comparisons=comparisons)
    assert comparison_eligibility.member_ids_by_record == {
        0: {"feature_r0"},
        1: {"feature_r1", "feature_r1_a", "feature_r1_b"},
        2: {"feature_r2", "feature_r2_b"},
    }
    assert comparison_eligibility.top_member_ids_by_record == {
        0: {"feature_r0"},
        1: {"feature_r1_a", "feature_r1_b"},
    }


@pytest.mark.linear
def test_linear_orthogroup_top_labels_suppress_only_lower_orthogroup_members() -> None:
    records = [
        _synthetic_label_record("record_0", 1000, "top_orthogroup_label", "top_non_orthogroup_label"),
        _synthetic_label_record("record_1", 1000, "lower_orthogroup_label", "lower_non_orthogroup_label"),
    ]
    top_feature_id = compute_feature_hash(records[0].features[0], record_id=records[0].id)
    lower_feature_id = compute_feature_hash(records[1].features[0], record_id=records[1].id)
    comparisons = [
        DataFrame(
            [
                _comparison_row(
                    orthogroup_id="og_1",
                    query_record_index=1,
                    query_feature_svg_id=lower_feature_id,
                    subject_record_index=0,
                    subject_feature_svg_id=top_feature_id,
                )
            ]
        )
    ]
    config_dict = modify_config_dict(
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "labels.linear.scope": 'orthogroup_top',
            "labels.filtering.blacklist_keywords": [],
        },
    )

    svg_content = assemble_linear_diagram_from_records(
        records,
        cfg=GbdrawConfig.from_dict(config_dict),
        protein_comparisons=comparisons,
        selected_features_set=["CDS"],
        legend="none",
    ).tostring()

    assert "top_orthogroup_label" in svg_content
    assert "lower_orthogroup_label" not in svg_content
    assert "top_non_orthogroup_label" in svg_content
    assert "lower_non_orthogroup_label" in svg_content


@pytest.mark.linear
def test_linear_orthogroup_top_labels_respect_per_feature_overrides() -> None:
    records = [
        _synthetic_label_record("record_0", 1000, "top_orthogroup_label", "top_non_orthogroup_label"),
        _synthetic_label_record("record_1", 1000, "lower_orthogroup_label", "lower_non_orthogroup_label"),
    ]
    top_feature_id = compute_feature_hash(records[0].features[0], record_id=records[0].id)
    lower_feature_id = compute_feature_hash(records[1].features[0], record_id=records[1].id)
    comparisons = [
        DataFrame(
            [
                _comparison_row(
                    orthogroup_id="og_1",
                    query_record_index=1,
                    query_feature_svg_id=lower_feature_id,
                    subject_record_index=0,
                    subject_feature_svg_id=top_feature_id,
                )
            ]
        )
    ]
    config_dict = modify_config_dict(
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "labels.linear.scope": 'orthogroup_top',
            "labels.filtering.blacklist_keywords": [],
        },
    )
    config_dict["labels"]["filtering"]["label_override_df"] = DataFrame(
        [
            ["record_1", "CDS", "hash", f"^{lower_feature_id}$", "forced_lower_label"],
            ["record_0", "CDS", "hash", f"^{top_feature_id}$", ""],
        ],
        columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
    )

    svg_content = assemble_linear_diagram_from_records(
        records,
        cfg=GbdrawConfig.from_dict(config_dict),
        protein_comparisons=comparisons,
        selected_features_set=["CDS"],
        legend="none",
    ).tostring()

    assert "forced_lower_label" in svg_content
    assert "lower_orthogroup_label" not in svg_content
    assert "top_orthogroup_label" not in svg_content
    assert "top_non_orthogroup_label" in svg_content
    assert "lower_non_orthogroup_label" in svg_content


@pytest.mark.linear
def test_linear_auto_label_font_size_resolves_once_per_diagram() -> None:
    config_dict = modify_config_dict(
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "labels.linear.scope": 'all',
            "labels.filtering.blacklist_keywords": [],
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    threshold = int(cfg.labels.length_threshold.linear)
    records = [
        _synthetic_label_record("short_record", threshold - 1, "linear_short_auto_label"),
        _synthetic_label_record("long_record", threshold + 1, "linear_long_auto_label"),
    ]

    canvas_config = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=max(len(record.seq) for record in records),
        profile=LinearRenderProfile(cfg),
        legend="none",
    )
    # The labels follow the size class of the longest record in the label scope, as the features do.
    expected_font_size = cfg.labels.font_size.linear.for_length_param("long")
    assert canvas_config.length_param == "long"
    svg_content = assemble_linear_diagram_from_records(
        records,
        cfg=cfg,
        selected_features_set=["CDS"],
        legend="none",
    ).tostring()

    font_sizes = _extract_label_font_sizes(
        svg_content,
        {"linear_short_auto_label", "linear_long_auto_label"},
    )
    assert set(font_sizes) == {"linear_short_auto_label", "linear_long_auto_label"}
    assert all(size == pytest.approx(expected_font_size) for size in font_sizes.values())


@pytest.mark.linear
def test_linear_labels_of_a_short_record_beside_a_long_one_use_the_long_size(
    gbdraw_runner, examples_dir: Path, test_inputs_dir: Path, tmp_path: Path
) -> None:
    """OV-482: 306,008 bp beside 50,466 bp drew short-class labels on long-class features."""
    returncode, output, svg_path = gbdraw_runner.run(
        "linear",
        [examples_dir / "MjeNMV.gb", test_inputs_dir / "BGC0000709.gbk"],
        "mixed_size_classes",
        tmp_path,
        extra_args=["--show_labels", "all"],
    )
    assert returncode == 0, output
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    label_font_sizes = {
        float(element.attrib["font-size"])
        for element in ET.parse(svg_path).getroot().iter()
        if element.tag.endswith("}text") and "data-label-feature-id" in element.attrib
    }
    assert label_font_sizes == {cfg.labels.font_size.linear.long}


@pytest.mark.linear
def test_linear_label_filter_does_not_narrow_the_label_scope(
    gbdraw_runner, examples_dir: Path, test_inputs_dir: Path, tmp_path: Path
) -> None:
    """The label scope is `--show_labels`: a whitelist that leaves labels only on the short record keeps the long size."""
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text("CDS\tproduct\ttransport protein\n", encoding="utf-8")
    returncode, output, svg_path = gbdraw_runner.run(
        "linear",
        [examples_dir / "MjeNMV.gb", test_inputs_dir / "BGC0000708.gbk"],
        "whitelist_on_short_record",
        tmp_path,
        extra_args=["--show_labels", "all", "--label_whitelist", str(whitelist)],
    )
    assert returncode == 0, output
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    labels = [
        element
        for element in ET.parse(svg_path).getroot().iter()
        if element.tag.endswith("}text") and "data-label-feature-id" in element.attrib
    ]
    assert {"".join(label.itertext()).strip() for label in labels} == {"transport protein"}
    assert {float(label.attrib["font-size"]) for label in labels} == {cfg.labels.font_size.linear.long}


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
    positive_label["feature_id"] = "linear-label-positive"
    negative_label["feature_id"] = "linear-label-negative"

    group = Group(id="test")
    drawer = LabelDrawer()
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
        label_font_size=18.0,
        input_filename="BGC0000708.gbk",
        qualifier_priority=("CDS", "gene"),
    )
    assert labels
    assert any(label.get("leader_line") for label in labels)
    assert any(not label.get("leader_line") for label in labels)
    _assert_no_label_bounds_overlap(labels)

    for label_text in ("livP", "livF"):
        label = next(label for label in labels if label["label_text"] == label_text)
        assert float(label["middle_y"]) == pytest.approx(_feature_adjacent_label_y(label))
        assert label["leader_line"] is False


@pytest.mark.linear
def test_linear_above_feature_bgc_multi_record_gene_labels_compact_without_overlap() -> None:
    labels = _prepare_linear_labels_for_records(
        input_filenames=(
            "BGC0000708.gbk",
            "BGC0000709.gbk",
            "BGC0000711.gbk",
            "BGC0000712.gbk",
            "BGC0000713.gbk",
        ),
        show_labels="first",
        label_placement="above_feature",
        label_rotation=45.0,
        label_font_size=18.0,
        qualifier_priority=("CDS", "gene"),
    )
    assert labels

    collision_gap = max(1.0, max(float(label["height_px"]) for label in labels) * 0.05)
    _assert_no_label_bounds_overlap(labels)

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

    livf = next(label for label in labels if label["label_text"] == "livF")
    assert float(livf["middle_y"]) == pytest.approx(_feature_adjacent_label_y(livf))
    assert livf["leader_line"] is False

    livy = next(label for label in labels if label["label_text"] == "livY")
    assert float(livy["middle_y"]) == pytest.approx(_feature_adjacent_label_y(livy))
    for neighbor in [label for label in labels if label["label_text"] in {"livW", "livZ"}]:
        assert not _rotated_label_rects_intersect(livy, neighbor)

    livp = next(label for label in labels if label["label_text"] == "livP")
    assert float(livp["middle_y"]) == pytest.approx(_feature_adjacent_label_y(livp))
    assert livp["leader_line"] is False



@pytest.mark.linear
@pytest.mark.parametrize("separate_strands", (False, True), ids=["one-side", "separate-strands"])
def test_linear_above_feature_horizontal_labels_take_the_auto_first_row_with_a_leader(separate_strands: bool) -> None:
    """LP-D120 (OV-1396): a horizontal above-feature label keeps the clearance off its feature that the auto
    placement keeps for its first row, with a leader from the feature edge to the label, as in auto; the
    lowest row is that first row (mirrored below the feature for a negative-strand label with separate
    strands). Before, the lowest row sat on the features without a leader."""

    spacing = 3.0
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        input_filename="NC_001416.gb",
        linear_label_spacing=spacing,
        separate_strands=separate_strands,
    )
    pitch = max(float(label["height_px"]) for label in labels) + spacing
    lowest = {True: 0, False: 0}
    for label in labels:
        above = bool(label["above_feature_place_above"])
        edge = float(label["feature_top_y"] if above else label["feature_bottom_y"])
        first_row_y = edge - pitch if above else edge + pitch
        middle_y = float(label["middle_y"])
        assert (middle_y <= first_row_y + 1e-6) if above else (middle_y >= first_row_y - 1e-6), label["label_text"]
        lowest[above] += middle_y == pytest.approx(first_row_y)
        assert label["leader_line"] is True, label["label_text"]
        assert (float(label["leader_start_x"]), float(label["leader_start_y"])) == pytest.approx(
            (float(label["feature_anchor_x"]), edge)
        )
        contact_y = middle_y + (0.5 if above else -0.5) * float(label["height_px"])
        assert (float(label["leader_end_x"]), float(label["leader_end_y"])) == pytest.approx(
            (float(label["feature_anchor_x"]), contact_y)
        )
    assert lowest[True] > 5
    assert (lowest[False] > 5) is separate_strands
    if not separate_strands:
        auto = _prepare_linear_labels(
            label_placement="auto", input_filename="NC_001416.gb", linear_label_spacing=spacing
        )
        auto_first_row = [float(label["middle_y"]) for label in auto if label["track_id"] == "track_1"]
        first_row_y = min(float(label["feature_top_y"]) for label in labels) - pitch
        assert auto_first_row and auto_first_row == pytest.approx([first_row_y] * len(auto_first_row))


@pytest.mark.linear
@pytest.mark.parametrize("separate_strands", (False, True), ids=["one-side", "separate-strands"])
def test_linear_above_feature_horizontal_labels_stack_at_the_auto_row_pitch(separate_strands: bool) -> None:
    """Study option P (row pitch above the first row): a horizontal above-feature label that moves off
    another label goes to a whole auto row off its own feature (row k at ``linear_row_y(edge, pitch, k)``),
    not the tallest label's height plus the packer's stacking gap; and the rows stay free of overlaps."""

    spacing = 3.0
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        input_filename="NC_001416.gb",
        linear_label_spacing=spacing,
        separate_strands=separate_strands,
    )
    pitch = max(float(label["height_px"]) for label in labels) + spacing
    rows = []
    for label in labels:
        above = bool(label["above_feature_place_above"])
        edge = float(label["feature_top_y"] if above else label["feature_bottom_y"])
        row = abs(float(label["middle_y"]) - edge) / pitch
        assert row == pytest.approx(round(row), abs=1e-6), label["label_text"]
        rows.append(round(row))
    assert max(rows) > 2
    _assert_no_label_bounds_overlap(labels)


@pytest.mark.linear
@pytest.mark.parametrize("separate_strands", (False, True), ids=["one-side", "separate-strands"])
def test_linear_above_feature_turned_labels_keep_their_place_on_the_feature(separate_strands: bool) -> None:
    """LP-D120: a turned above-feature label keeps its place on its feature (no leader, the small gap) unless
    it moves off another label; the first-row clearance is only for horizontal labels."""

    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=45.0,
        input_filename="NC_001416.gb",
        separate_strands=separate_strands,
    )
    on_feature = [label for label in labels if not label["leader_line"]]
    assert len(on_feature) > 5
    for label in on_feature:
        assert float(label["middle_y"]) == pytest.approx(_feature_adjacent_label_y(label))
    for label in labels:
        if label["leader_line"]:
            assert abs(float(label["middle_y"]) - float(label["feature_middle_y"])) > abs(
                _feature_adjacent_label_y(label) - float(label["feature_middle_y"])
            )


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
@pytest.mark.parametrize("label_rotation", (0.0, 45.0))
def test_linear_above_feature_resolves_overlaps_above_300_labels(
    label_rotation: float,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OV-695: a record with more than 300 above-feature labels gets the same overlap fix.

    OV-841: a label climbing a stack makes one broad-phase query and solves each
    neighbour's overlap range once, not once per round. Measured on this input:
    queries 11,073 (0 deg) and 7,967 (45 deg) before, one per label after; ranges at
    45 deg 73,930 before, 10,560 after (cap 25,000).
    """
    calls = {"query": 0, "range": 0}
    resolved_labels = 0
    inside_resolver = False

    def counting(name: str, original: Callable[_P, _R]) -> Callable[_P, _R]:
        def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _R:
            if inside_resolver:
                calls[name] += 1
            return original(*args, **kwargs)

        return wrapper

    resolve = linear_labels_module._resolve_above_feature_label_overlaps

    def counted_resolve(
        labels: list[dict], min_gap_px: float, keep_clear: Sequence[Aabb] = (), row_pitch: float | None = None
    ) -> None:
        nonlocal resolved_labels, inside_resolver
        resolved_labels += len(labels)
        inside_resolver = True
        try:
            resolve(labels, min_gap_px, keep_clear, row_pitch)
        finally:
            inside_resolver = False

    monkeypatch.setattr(linear_labels_module, "_resolve_above_feature_label_overlaps", counted_resolve)
    for index_class in (AabbIndex, IntervalIndex):
        monkeypatch.setattr(index_class, "query", counting("query", index_class.query))
    monkeypatch.setattr(
        linear_labels_module,
        "_vertical_overlap_shift_range",
        counting("range", linear_labels_module._vertical_overlap_shift_range),
    )
    labels = _prepare_linear_labels(
        label_placement="above_feature",
        label_rotation=label_rotation,
        separate_strands=True,
        input_filename="AP027078.gb",
    )
    monkeypatch.undo()
    assert resolved_labels > 300
    assert calls["query"] <= resolved_labels
    if label_rotation:
        assert calls["range"] <= 25_000

    overlapping = [
        (left, right)
        for left, right in candidate_aabb_pairs([_row_label_aabb(label) for label in labels], bucket_size=64.0)
        if _rotated_label_rects_intersect(labels[left], labels[right])
    ]
    assert overlapping == []
    assert any(label.get("leader_line") for label in labels)


@pytest.mark.linear
def test_linear_precalc_includes_above_feature_rotated_embedded_labels() -> None:
    input_path = Path(__file__).parent.parent / "examples" / "MjeNMV.gb"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        {
            "labels.linear.scope": 'all',
            "canvas.strandedness": False,
            "labels.filtering.blacklist_keywords": [],
            "labels.font_size.short": 14.0,
            "labels.font_size.long": 14.0,
            "labels.font_size.linear.short": 14.0,
            "labels.font_size.linear.long": 14.0,
            "labels.linear.placement": 'above_feature',
            "labels.linear.rotation": 45.0,
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        profile=profile,
        legend="none",
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        profile=profile,
        canvas_config=canvas_cfg,
    )

    required_label_height, _, record_label_heights = _precalculate_label_dimensions(
        [record],
        feature_cfg,
        canvas_cfg,
        render_context=make_linear_record_render_context(profile),
        sequence_widths=alone_sequence_widths([record], canvas_cfg),
    )

    assert required_label_height > 0
    assert required_label_height > canvas_cfg.original_vertical_offset
    assert record_label_heights[0] == pytest.approx(required_label_height)


@pytest.mark.linear
def test_linear_canvas_alignment_width_starts_at_figure_width() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        profile=profile,
        legend="none",
    )

    assert canvas_cfg.alignment_width == pytest.approx(canvas_cfg.fig_width)


@pytest.mark.linear
def test_linear_precalculated_labels_match_final_alignment_width_positions() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
    record = SeqIO.read(str(input_path), "genbank")

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        {
            "labels.linear.scope": 'all',
            "canvas.strandedness": False,
            "labels.filtering.blacklist_keywords": [],
            "canvas.linear.track_layout": 'middle',
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        profile=profile,
        legend="none",
    )
    feature_cfg = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        profile=profile,
        canvas_config=canvas_cfg,
    )

    _, all_labels, _ = _precalculate_label_dimensions(
        [record],
        feature_cfg,
        canvas_cfg,
        render_context=make_linear_record_render_context(profile),
        sequence_widths=alone_sequence_widths([record], canvas_cfg),
    )
    precalculated_tuples = sorted(
        (
            float(label["feature_start_x"]),
            float(label["feature_end_x"]),
            float(label["middle_x"]),
        )
        for label in all_labels[0]
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
        for label in _labels_alone(feature_dict, len(record.seq), canvas_cfg, profile)
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
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "labels.spacing.circular": 7.5,
            "labels.spacing.linear": 9.0,
        },
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
    # At 10 px each label keeps 1.25 px free on each side: the last one, 1 px after its neighbour,
    # takes another row.
    labels = [
        {"start": start, "end": end, "font_size": 10.0}
        for start, end in ((0.0, 50.0), (20.0, 70.0), (80.0, 110.0), (45.0, 90.0), (112.0, 140.0), (141.0, 160.0))
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
    assert indexed_assignments[-1] != indexed_assignments[-2]


@pytest.mark.linear
def test_linear_labels_in_one_row_keep_a_quarter_em_between_them() -> None:
    """Neighbouring labels in one label row keep at least ``ROW_LABEL_GAP_EM`` of their font size apart,
    so two names never read as one word (HmmtDNA at 24 pt drew "tRNA-Ser" 0.31 px from "tRNA-Lys")."""

    labels = _prepare_linear_labels(
        label_placement="auto", label_font_size=24.0, separate_strands=True, track_layout="above",
        input_filename="HmmtDNA.gbk",
    )
    rows: dict[float, list[dict]] = {}
    for label in labels:
        if not label["is_embedded"]:
            rows.setdefault(round(float(label["middle_y"]), 6), []).append(label)
    close = []
    for row in rows.values():
        row.sort(key=lambda label: float(label["start"]))
        for left, right in zip(row, row[1:]):
            if float(right["start"]) - float(left["end"]) < ROW_LABEL_GAP_EM * 24.0 - 1e-6:
                close.append((left["label_text"], right["label_text"]))

    assert len(rows) > 1 and close == []


def _box_along_text(label: dict, room_px: float) -> list[tuple[float, float]]:
    """A label's box with ``room_px`` added at both ends along its text (its own reading direction)."""

    from gbdraw.layout.text_geometry import anchor_x_bounds, rotate_point

    x_min, x_max = anchor_x_bounds(float(label["width_px"]), label["text_anchor"])
    half_height = 0.5 * float(label["height_px"])
    return [
        (float(label["middle_x"]) + px, float(label["middle_y"]) + py)
        for px, py in (
            rotate_point(point, float(label["rotation_deg"]))
            for point in (
                (x_min - room_px, -half_height),
                (x_max + room_px, -half_height),
                (x_max + room_px, half_height),
                (x_min - room_px, half_height),
            )
        )
    ]


@pytest.mark.linear
@pytest.mark.parametrize("label_rotation", (0.0, 45.0), ids=["horizontal", "turned"])
def test_linear_above_feature_labels_on_one_line_keep_a_quarter_em_along_their_text(label_rotation: float) -> None:
    """OV-1397 (horizontal) and OV-1363 (turned): above-feature labels that lie on one line keep the
    ``ROW_LABEL_GAP_EM`` row gap along their text, as the auto rows do (LP-D104). Lambda at 24 px drew
    "minor tail protein" 4.4 px from "Lom outer membrane protein" and turned "integrase" 3.1 px along
    the text from "ead/Ea22-like family protein"."""

    from gbdraw.layout.text_geometry import convex_polygons_intersect

    labels = _prepare_linear_labels(
        label_placement="above_feature", label_rotation=label_rotation, label_font_size=24.0,
        input_filename="NC_001416.gb",
    )
    half_gap = 0.5 * ROW_LABEL_GAP_EM * 24.0 - 1e-6
    boxes = [(label["label_text"], _box_along_text(label, half_gap)) for label in labels]
    close = [
        (left_text, right_text)
        for index, (left_text, left_box) in enumerate(boxes)
        for right_text, right_box in boxes[index + 1:]
        if convex_polygons_intersect(left_box, right_box, touching=False)
    ]
    assert len(labels) > 50 and close == []


@pytest.mark.linear
@pytest.mark.parametrize("label_rotation", (0.0, 45.0), ids=["horizontal", "turned"])
def test_linear_fit_counts_above_feature_labels_closer_than_the_row_gap(label_rotation: float) -> None:
    """The fit count reads the row gap of above-feature labels as the packer keeps it: two labels end to
    end on one line 0.1 em apart overlap, 0.3 em apart they fit."""

    from gbdraw.labels.linear import linear_row_label_fit

    labels = _prepare_linear_labels(
        label_placement="above_feature", label_rotation=label_rotation, label_font_size=24.0,
        input_filename="NC_001416.gb",
    )
    first = labels[0]
    radians = math.radians(float(first["rotation_deg"]))

    def next_on_line(gap_em: float) -> dict:
        step = float(first["width_px"]) + gap_em * 24.0
        return {
            **first,
            "middle_x": float(first["middle_x"]) + step * math.cos(radians),
            "middle_y": float(first["middle_y"]) + step * math.sin(radians),
            "leader_line": False,
        }

    alone = {**first, "leader_line": False}
    assert linear_row_label_fit([([alone, next_on_line(0.1)], 0.0)]).overlap_pairs == 1
    assert linear_row_label_fit([([alone, next_on_line(0.3)], 0.0)]).overlap_pairs == 0


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


_SVG_NS = "{http://www.w3.org/2000/svg}"
_TRANSLATE = re.compile(r"translate\(\s*([-+0-9.eE]+)(?:[\s,]+([-+0-9.eE]+))?\s*\)")


def _painted_feature_label_boxes(svg_path: Path) -> list[tuple[str, tuple[float, float, float, float]]]:
    """Each drawn feature label's painted box in canvas px: the advance the renderer anchors the text
    on, by the ink's top and bottom (left, right, top, bottom). Only translations move the labels here."""

    boxes: list[tuple[str, tuple[float, float, float, float]]] = []

    def walk(element: ET.Element, dx: float, dy: float) -> None:
        transform = element.get("transform", "").strip()
        if transform:
            assert _TRANSLATE.sub("", transform).strip() == "", transform
            for match in _TRANSLATE.finditer(transform):
                dx += float(match.group(1))
                dy += float(match.group(2) or 0.0)
        if element.tag == f"{_SVG_NS}text" and element.get("data-label-feature-id"):
            text = "".join(element.itertext())
            family = str(element.get("font-family"))
            size = float(str(element.get("font-size")))
            anchor = str(element.get("text-anchor", "start"))
            ink = text_ink_offsets(
                text, family, size, dominant_baseline=str(element.get("dominant-baseline")), text_anchor=anchor
            )
            assert ink is not None
            advance, _height = calculate_svg_bbox_dimensions(text, family, size, 96)
            left = dx + float(str(element.get("x"))) - TEXT_ANCHOR_SHARE[anchor] * advance
            y = dy + float(str(element.get("y")))
            boxes.append((text, (left, left + advance, y + ink[2], y + ink[3])))
        for child in element:
            walk(child, dx, dy)

    walk(ET.parse(svg_path).getroot(), 0.0, 0.0)
    return boxes


def _painted_overlaps(svg_path: Path) -> list[tuple[str, str]]:
    boxes = _painted_feature_label_boxes(svg_path)
    return [
        (first_text, second_text)
        for index, (first_text, first) in enumerate(boxes)
        for second_text, second in boxes[index + 1:]
        if min(first[1], second[1]) - max(first[0], second[0]) > 1e-6
        and min(first[3], second[3]) - max(first[2], second[2]) > 1e-6
    ]


def _draw_linear(tmp_path: Path, args: list[str]) -> Path:
    output = tmp_path / "out"
    linear_cli_module.linear_main([*args, "-o", str(output), "-f", "svg"])
    return output.with_suffix(".svg")


@pytest.mark.linear
def test_linear_labels_are_sized_by_their_advance_so_drawn_letters_keep_apart(tmp_path: Path) -> None:
    """OV-1294: lambda's 24 px labels were sized as 24 pt ink widths, so "CIII anti-termination" and
    "replication initiation O-like" shared ink in one row while the count reported none."""

    lambda_gbk = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        svg = _draw_linear(tmp_path, ["--gbk", str(lambda_gbk), "--show_labels"])

    assert _painted_overlaps(svg) == []


@pytest.mark.linear
def test_linear_labels_of_records_in_one_row_keep_apart(tmp_path: Path) -> None:
    """OV-1188: the labels of two records in one row were packed and counted per record, so labels
    at the junction overlapped (HmmtDNA twice at 40 pt: 3 pairs)."""

    hmm = str(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk")
    svg = _draw_linear(
        tmp_path,
        [
            "--gbk", hmm, hmm, "--show_labels", "all", "--label_font_size", "40",
            "--multi_record_position", "#1@1", "--multi_record_position", "#2@1", "--linear_record_gap", "5",
        ],
    )

    assert _painted_overlaps(svg) == []


def _row_label(start: float, width: float, *, middle_y: float = -30.0) -> dict:
    return {
        "label_text": "x",
        "middle": start + width / 2.0,
        "middle_x": start + width / 2.0,
        "middle_y": middle_y,
        "start": start,
        "end": start + width,
        "width_px": width,
        "height_px": 10.0,
        "font_size": 10.0,
        "is_embedded": False,
        "feature_middle_y": 0.0,
        "rotation_deg": 0.0,
        "text_anchor": "middle",
    }


@pytest.mark.linear
def test_linear_label_fit_counts_labels_of_records_in_one_row_together() -> None:
    """OV-1188: a label of the second record over one of the first, where the row draws them, counts
    for the second record; the same labels in two rows fit."""

    first = [_row_label(90.0, 20.0)]
    second = [_row_label(0.0, 20.0)]
    records = [("first", first), ("second", second)]

    one_row = linear_label_fit(
        records, 10.0, frames=[LinearLabelFrame(row=1, x=0.0), LinearLabelFrame(row=1, x=95.0)]
    )
    two_rows = linear_label_fit(
        records, 10.0, frames=[LinearLabelFrame(row=1, x=0.0), LinearLabelFrame(row=2, x=95.0)]
    )

    assert [(record.record_id, record.fit.overlap_pairs) for record in one_row.records] == [("second", 1)]
    assert two_rows.records == ()


@pytest.mark.linear
def test_linear_labels_of_records_in_one_row_are_packed_together() -> None:
    """OV-1188: the labels of the second record in a row take the rows the first record's labels leave
    free where the row draws them; the first record keeps the rows it has alone."""

    alone = _prepare_linear_labels_for_records(
        input_filenames=("HmmtDNA.gbk", "HmmtDNA.gbk"),
        show_labels="all",
        label_placement="auto",
        label_font_size=40.0,
        all_records=True,
    )
    # The second record starts 5 px after the first one ends, as ``--linear_record_gap 5`` draws it.
    width = max(float(label["feature_end_x"]) for label in alone[0])
    frames = [LinearLabelFrame(row=1, x=0.0), LinearLabelFrame(row=1, x=width + 5.0)]
    packed = _prepare_linear_labels_for_records(
        input_filenames=("HmmtDNA.gbk", "HmmtDNA.gbk"),
        show_labels="all",
        label_placement="auto",
        label_font_size=40.0,
        all_records=True,
        label_frames=frames,
    )

    assert packed[0] == alone[0]
    row_labels = [
        {**label, "middle_x": float(label["middle_x"]) + frame.x}
        for labels, frame in zip(packed, frames)
        for label in labels
    ]
    _assert_no_label_bounds_overlap(row_labels)
    assert any(
        moved["track_id"] != kept["track_id"] for moved, kept in zip(packed[1], alone[1], strict=True)
    )


def _svg_translation(element: ET.Element) -> tuple[float, float]:
    """The sum of the ``translate(...)`` steps in an element's ``transform``."""

    dx = dy = 0.0
    for tx, ty in re.findall(r"translate\(([-\d.e]+)[ ,]+([-\d.e]+)\)", element.attrib.get("transform", "")):
        dx += float(tx)
        dy += float(ty)
    return dx, dy


def _svg_text_boxes(svg: str, cfg: GbdrawConfig) -> tuple[list[list[tuple[float, float]]], list[list[tuple[float, float]]]]:
    """The drawn feature labels and definition lines of a one-record Linear SVG as polygons in canvas
    coordinates: each text's advance by its font box, anchored and turned as the SVG draws it."""

    font_family = cfg.objects.text.font_family
    dpi = cfg.canvas.dpi
    labels: list[list[tuple[float, float]]] = []
    definitions: list[list[tuple[float, float]]] = []

    def walk(element: ET.Element, offset: tuple[float, float], in_definition: bool) -> None:
        dx, dy = _svg_translation(element)
        offset = (offset[0] + dx, offset[1] + dy)
        in_definition = in_definition or element.attrib.get("data-gbdraw-role") == "record-definition"
        if element.tag.endswith("}text"):
            text = "".join(element.itertext()).strip()
            is_label = "data-label-feature-id" in element.attrib
            if text and (is_label or in_definition):
                size = float(element.attrib["font-size"])
                width, _ = calculate_svg_bbox_dimensions(text, font_family, size, dpi)
                _, height = calculate_bbox_dimensions(text, font_family, size, dpi)
                x, y = float(element.attrib.get("x", 0.0)), float(element.attrib.get("y", 0.0))
                left = {"start": 0.0, "end": -width}.get(element.attrib.get("text-anchor", "start"), -0.5 * width)
                corners = [(left, -0.5 * height), (left + width, -0.5 * height), (left + width, 0.5 * height), (left, 0.5 * height)]
                turn = re.search(r"rotate\(([-\d.e]+)", element.attrib.get("transform", ""))
                angle = math.radians(float(turn.group(1))) if turn else 0.0
                polygon = [
                    (offset[0] + x + cx * math.cos(angle) - cy * math.sin(angle), offset[1] + y + cx * math.sin(angle) + cy * math.cos(angle))
                    for cx, cy in corners
                ]
                (labels if is_label else definitions).append(polygon)
        for child in element:
            walk(child, offset, in_definition)

    walk(ET.fromstring(svg), (0.0, 0.0), False)
    return labels, definitions


@pytest.mark.linear
@pytest.mark.parametrize("label_rotation", (0.0, 45.0))
def test_linear_above_feature_labels_keep_clear_of_the_definition_text(label_rotation: float) -> None:
    """OV-1395: the lowest row of above-feature labels left the first label over the record's
    definition text (lambda: "terminase small subunit" over "NC_001416.1")."""

    from gbdraw.layout.text_geometry import convex_polygons_intersect

    record = SeqIO.read(str(Path(__file__).parent / "test_inputs" / "NC_001416.gb"), "genbank")
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {
            "labels.linear.scope": "all",
            "labels.linear.placement": "above_feature",
            "labels.linear.rotation": label_rotation,
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = assemble_linear_diagram_from_records([record], cfg=cfg, legend="none")
    labels, definitions = _svg_text_boxes(drawing.tostring(), cfg)

    assert labels and len(definitions) == 2
    covering = [
        label for label in labels if any(convex_polygons_intersect(label, line, touching=False) for line in definitions)
    ]
    assert covering == []
    assert sum(record.fit.definition_overlaps for record in drawing._gbdraw_linear_label_fit.records) == 0


@pytest.mark.linear
def test_linear_above_feature_label_moves_off_a_keep_clear_box() -> None:
    """A record's ``keep_clear`` box (its definition text) is placed before its labels: an
    above-feature label over it moves away from the axis with a leader, as off another label."""


    record = SeqRecord(Seq("A" * 1000), id="r", name="r", description="r")
    record.features = [SeqFeature(FeatureLocation(0, 40, strand=1), type="CDS", qualifiers={"product": ["a long first label"]})]
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {"labels.linear.scope": "all", "labels.linear.placement": "above_feature"},
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(num_of_entries=1, longest_genome=1000, profile=profile, legend="none")
    color_table, default_colors = preprocess_color_tables(None, load_default_colors("", "default"))
    feature_dict, _ = create_feature_dict(
        record, color_table, cfg.objects.features.features_drawn, default_colors,
        canvas_cfg.strandedness, canvas_cfg.resolve_overlaps,
        preprocess_label_filtering(cfg.labels.filtering.as_dict()),
    )
    record_labels = linear_record_labels(
        feature_dict, 1000, canvas_cfg.alignment_width, canvas_cfg.cds_height,
        canvas_cfg.strandedness, canvas_cfg.track_layout, canvas_cfg.track_axis_gap, profile,
    )
    (alone,) = place_linear_label_row([(record_labels, 0.0)], profile=profile, track_layout=canvas_cfg.track_layout)[0]
    left, _right, top, bottom = calculate_label_bounds(alone)
    assert left < 0.0

    header = Aabb(left - 50.0, top - 2.0, -5.0, bottom + 30.0)
    (placed,) = place_linear_label_row(
        [(record_labels._replace(keep_clear=(header,)), 0.0)], profile=profile, track_layout=canvas_cfg.track_layout
    )[0]
    _left, _right, _top, placed_bottom = calculate_label_bounds(placed)
    assert placed_bottom < header.min_y
    assert placed["leader_line"] is True


@pytest.mark.linear
@pytest.mark.parametrize("track_layout", ("above", "below"))
def test_linear_auto_label_over_a_keep_clear_box_takes_the_next_row(track_layout: str) -> None:
    """OV-1398: an ``auto`` label whose row would put it over its record's definition text (a
    ``keep_clear`` box) takes the next row, as off another label (render-stage L2 with
    ``--region``: "baculoviral IAP repeat-containing protein" over "LC738868.1"). A box that the
    row leaves clear moves nothing."""

    from gbdraw.layout.text_geometry import convex_polygons_intersect

    record = SeqRecord(Seq("A" * 1000), id="r", name="r", description="r")
    record.features = [SeqFeature(FeatureLocation(0, 40, strand=1), type="CDS", qualifiers={"product": ["a long first label"]})]
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {"labels.linear.scope": "all", "canvas.linear.track_layout": track_layout},
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_cfg = LinearCanvasConfigurator(num_of_entries=1, longest_genome=1000, profile=profile, legend="none")
    color_table, default_colors = preprocess_color_tables(None, load_default_colors("", "default"))
    feature_dict, _ = create_feature_dict(
        record, color_table, cfg.objects.features.features_drawn, default_colors,
        canvas_cfg.strandedness, canvas_cfg.resolve_overlaps,
        preprocess_label_filtering(cfg.labels.filtering.as_dict()),
    )
    record_labels = linear_record_labels(
        feature_dict, 1000, canvas_cfg.alignment_width, canvas_cfg.cds_height,
        canvas_cfg.strandedness, canvas_cfg.track_layout, canvas_cfg.track_axis_gap, profile,
    )
    (alone,) = place_linear_label_row([(record_labels, 0.0)], profile=profile, track_layout=canvas_cfg.track_layout)[0]
    assert alone["track_id"] == "track_1" and not alone["is_embedded"]
    left, _right, top, bottom = calculate_label_bounds(alone)
    assert left < 0.0

    def box_polygon(box: Aabb) -> tuple[tuple[float, float], ...]:
        return ((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y))

    # The definition text beside the features, reaching 1 px into the first row's letters.
    toward_axis = (bottom - 1.0, bottom + 40.0) if track_layout == "above" else (top - 40.0, top + 1.0)
    header = Aabb(left - 50.0, toward_axis[0], -5.0, toward_axis[1])
    (placed,) = place_linear_label_row(
        [(record_labels._replace(keep_clear=(header,)), 0.0)], profile=profile, track_layout=canvas_cfg.track_layout
    )[0]
    assert placed["track_id"] == "track_2"
    assert not convex_polygons_intersect(_label_corners(placed), box_polygon(header), touching=False)

    clear = Aabb(left - 50.0, bottom + 1.0, -5.0, bottom + 40.0) if track_layout == "above" else Aabb(left - 50.0, top - 40.0, -5.0, top - 1.0)
    (unmoved,) = place_linear_label_row(
        [(record_labels._replace(keep_clear=(clear,)), 0.0)], profile=profile, track_layout=canvas_cfg.track_layout
    )[0]
    assert unmoved["track_id"] == "track_1" and unmoved["middle_y"] == alone["middle_y"]


@pytest.mark.linear
def test_linear_label_fit_counts_a_label_over_the_definition_text() -> None:
    """The fit count sees a feature label the placement could not keep off the definition text: it
    warns ("labels cover the definition text"), it does not stop."""


    label = {
        "label_text": "x", "middle_x": 10.0, "middle_y": -20.0, "width_px": 40.0, "height_px": 10.0,
        "rotation_deg": 0.0, "text_anchor": "middle", "is_embedded": True, "font_size": 10.0,
        "feature_middle_y": 0.0,
    }
    clear = linear_label_fit([("r", [label])], 10.0, definition_boxes=[(Aabb(-50.0, -30.0, -20.0, -10.0),)])
    covered = linear_label_fit([("r", [label])], 10.0, definition_boxes=[(Aabb(-50.0, -30.0, -5.0, -10.0),)])

    assert clear.records == ()
    assert [record.fit.definition_overlaps for record in covered.records] == [1]
    assert covered.outcome == "warn"


@pytest.mark.linear
@pytest.mark.parametrize("shared_row", (False, True))
@pytest.mark.parametrize("translation_x", (-240.0, 160.0))
def test_linear_label_keep_clear_box_is_the_drawn_definition_box(
    monkeypatch: pytest.MonkeyPatch, translation_x: float, shared_row: bool
) -> None:
    """The definition box the labels keep clear of is the box the fit count reads: with
    ``keep_definition_left_aligned`` the text starts at the column left of the translated records,
    so a record translation moves the text relative to its record before the labels are placed
    (a record alone in its row: its own text; a shared row: the row's heading)."""

    import gbdraw.diagrams.linear.assemble as assemble_module
    from gbdraw.api.options import LinearMultiRecordOptions, LinearRecordTranslation
    from gbdraw.api.requests import RecordPresentation

    captured: dict[str, object] = {}

    def capturing(name: str, keyword: str, original: Callable[_P, _R]) -> Callable[_P, _R]:
        def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _R:
            if kwargs.get(keyword) is not None:
                captured[name] = list(kwargs[keyword])
            return original(*args, **kwargs)

        return wrapper

    monkeypatch.setattr(
        assemble_module,
        "_precalculate_label_dimensions",
        capturing("packer", "label_keep_clear", assemble_module._precalculate_label_dimensions),
    )
    monkeypatch.setattr(
        assemble_module, "linear_label_fit", capturing("drawn", "definition_boxes", assemble_module.linear_label_fit)
    )
    records = [
        _synthetic_label_record("first", 2000, "alpha"),
        _synthetic_label_record("second", 1500, "beta"),
        _synthetic_label_record("third", 1200, "gamma"),
    ]
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {"canvas.linear.keep_definition_left_aligned": True, "labels.linear.scope": "all"},
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assemble_linear_diagram_from_records(
            records,
            cfg=cfg,
            legend="none",
            layout=LinearMultiRecordOptions(
                multi_record_positions=("#1@1", "#2@2", "#3@2") if shared_row else None,
                record_translations=(LinearRecordTranslation("record-2", x=translation_x),),
            ),
            _record_presentations=[RecordPresentation(label="Same row") for _record in records],
        )

    packer, drawn = captured["packer"], captured["drawn"]
    assert isinstance(packer, list) and isinstance(drawn, list) and any(packer)
    if shared_row:
        # The fit count also reads each record's own header, above its labels' band.
        assert all(set(boxes) <= set(drawn_boxes) for boxes, drawn_boxes in zip(packer, drawn))
    else:
        assert packer == drawn


@pytest.mark.linear
@pytest.mark.parametrize("shared_row", (True, False))
def test_linear_horizontal_places_labels_and_vertical_plans_are_solved_once(
    monkeypatch: pytest.MonkeyPatch, shared_row: bool
) -> None:
    """Every record gets its x and width from one horizontal solve, in both row modes, which the
    labels, the keep-clear boxes, the collision bands and the drawing all read; the row stacking
    keeps them (OV-1400). The labels and the vertical plans are built once, at those widths."""

    import gbdraw.diagrams.linear.assemble as assemble_module
    from gbdraw.api.options import LinearMultiRecordOptions

    horizontal: list[object] = []
    stacked: list[object] = []
    label_widths: list[object] = []
    plans: list[object] = []

    def counting(calls: list[object], original, read=lambda result, kwargs: result):
        def solve(*args: object, **kwargs: object):
            result = original(*args, **kwargs)
            calls.append(read(result, kwargs))
            return result

        return solve

    monkeypatch.setattr(
        assemble_module,
        "solve_linear_horizontal_layout",
        counting(horizontal, assemble_module.solve_linear_horizontal_layout),
    )
    monkeypatch.setattr(
        assemble_module, "solve_linear_layout", counting(stacked, assemble_module.solve_linear_layout)
    )
    monkeypatch.setattr(
        assemble_module,
        "_precalculate_label_dimensions",
        counting(
            label_widths,
            assemble_module._precalculate_label_dimensions,
            lambda _result, kwargs: kwargs.get("sequence_widths"),
        ),
    )
    monkeypatch.setattr(
        assemble_module,
        "_build_linear_record_vertical_plans",
        counting(plans, assemble_module._build_linear_record_vertical_plans),
    )
    records = [
        _synthetic_label_record("first", 2000, "alpha"),
        _synthetic_label_record("second", 1500, "beta"),
        _synthetic_label_record("third", 1200, "gamma"),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assemble_linear_diagram_from_records(
            records,
            cfg=GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")),
            legend="none",
            layout=(
                LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@2", "#3@2"))
                if shared_row
                else None
            ),
        )

    assert len(horizontal) == 1
    assert len(plans) == 1
    places = horizontal[0].placements  # type: ignore[attr-defined]
    assert label_widths == [[p.sequence_width for p in sorted(places, key=lambda p: p.record_index)]]
    assert len(stacked) == (1 if shared_row else 0)
    if shared_row:
        drawn = stacked[0].placements  # type: ignore[attr-defined]
        assert [(p.record_index, p.x, p.sequence_width, p.px_per_bp) for p in places] == [
            (p.record_index, p.x, p.sequence_width, p.px_per_bp) for p in drawn
        ]


@pytest.mark.linear
@pytest.mark.parametrize("shared_row", (False, True))
@pytest.mark.parametrize("translation_x", (-240.0, 160.0))
def test_linear_row_gap_reads_the_drawn_collision_bands(
    monkeypatch: pytest.MonkeyPatch, translation_x: float, shared_row: bool
) -> None:
    """The bands the vertical gap between rows is solved from are the bands the drawing uses: at
    each record's final translation, with the definition column's left and text anchor of the
    drawing (OV-1400; a record alone in its row: the gap loop; a shared row: the row solve)."""

    import gbdraw.diagrams.linear.assemble as assemble_module
    import gbdraw.layout.linear_multi_record as multi_record_module
    from gbdraw.api.options import LinearMultiRecordOptions, LinearRecordTranslation
    from gbdraw.api.requests import RecordPresentation
    from gbdraw.layout.linear import CollisionBand, resolve_axis_gap

    read: list[CollisionBand] = []
    drawn: list[CollisionBand] = []

    def reading_gap(current: Sequence[CollisionBand], following: Sequence[CollisionBand], **kwargs: object):
        read.extend((*current, *following))
        return resolve_axis_gap(current, following, **kwargs)  # type: ignore[arg-type]

    def placing(*args: object, **kwargs: object):
        drawn.extend(band for bands in kwargs["collision_bands"] for band in bands)  # type: ignore[attr-defined]
        return original_place(*args, **kwargs)

    original_place = assemble_module._place_linear_annotation_labels_by_record
    monkeypatch.setattr(multi_record_module, "resolve_axis_gap", reading_gap)
    monkeypatch.setattr(assemble_module, "resolve_axis_gap", reading_gap)
    monkeypatch.setattr(assemble_module, "_place_linear_annotation_labels_by_record", placing)
    records = [
        _synthetic_label_record("first", 2000, "alpha"),
        _synthetic_label_record("second", 1500, "beta"),
        _synthetic_label_record("third", 1200, "gamma"),
    ]
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {
            "canvas.linear.keep_definition_left_aligned": True,
            "canvas.linear.align_center": True,
            "labels.linear.scope": "all",
        },
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assemble_linear_diagram_from_records(
            records,
            cfg=GbdrawConfig.from_dict(config_dict),
            legend="none",
            layout=LinearMultiRecordOptions(
                multi_record_positions=("#1@1", "#2@2", "#3@2") if shared_row else None,
                record_translations=(LinearRecordTranslation("record-2", x=translation_x),),
            ),
            _record_presentations=[RecordPresentation(label="Same row") for _record in records],
        )

    def key(band: CollisionBand) -> tuple[str, float, float, float, float]:
        return (band.kind, *(round(value, 6) for value in (band.x_start, band.x_end, band.top_y, band.bottom_y)))

    read_definitions = {key(band) for band in read if band.kind == "definition"}
    assert read_definitions
    assert read_definitions <= {key(band) for band in drawn if band.kind == "definition"}
    assert {key(band) for band in read} <= {key(band) for band in drawn}
