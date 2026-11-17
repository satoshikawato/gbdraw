from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import pytest

import gbdraw.labels.circular as circular_labels_module
from gbdraw.api.config import apply_config_overrides, load_default_config
from gbdraw.circular import _get_args, run_circular_from_namespace
from gbdraw.exceptions import ValidationError
from gbdraw.labels.circular_radial import place_radial_labels, radial_text_orientation


def _label(index: int, angle_deg: float, *, inner: bool = False, width: float = 80.0) -> dict:
    middle = (angle_deg % 360.0) / 360.0 * 10_000.0
    radius = 380.0 if not inner else 270.0
    return {
        "stable_id": f"feature-{index}",
        "input_order": index,
        "label_text": f"label {index}",
        "middle": middle,
        "preferred_middle": middle,
        "width_px": width,
        "height_px": 12.0,
        "middle_x": radius,
        "middle_y": 0.0,
        "feature_anchor_x": 340.0 if not inner else 300.0,
        "feature_anchor_y": 0.0,
        "feature_inner_radius_px": 290.0,
        "feature_outer_radius_px": 350.0,
        "is_inner": inner,
    }


@pytest.mark.parametrize("angle", [0.0, 90.0, 180.0, 270.0, 359.999999])
@pytest.mark.parametrize("inner", [False, True])
def test_radial_orientation_is_readable_and_stable(angle: float, inner: bool) -> None:
    rotation, anchor, baseline = radial_text_orientation(angle, is_inner=inner)
    assert -90.0 <= rotation <= 90.0
    assert anchor in {"start", "end"}
    assert baseline == "middle"


def test_radial_dense_seam_cluster_is_collision_free_and_deterministic() -> None:
    labels = [_label(index, angle) for index, angle in enumerate((358.0, 359.0, 0.0, 0.0, 1.0, 2.0))]
    first = place_radial_labels(
        labels,
        (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
    )
    second = place_radial_labels(
        deepcopy(labels),
        (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
    )

    assert first.collision_free
    assert len(first.labels) == len(labels)
    assert [label["placed_angle_deg"] for label in first.labels] == [
        label["placed_angle_deg"] for label in second.labels
    ]
    assert all(label["placement"] == "radial" for label in first.labels)
    assert all(len(label["oriented_corners"]) == 4 for label in first.labels)


def _label_on_feature(index: int, angle_deg: float, *, inner: bool) -> dict:
    """A label whose leader starts on its own feature, at the feature's angle (no shared anchor)."""

    label = _label(index, angle_deg, inner=inner)
    direction = math.radians(angle_deg) - 0.5 * math.pi
    anchor_radius, elbow_radius = (300.0, 280.0) if inner else (340.0, 380.0)
    label["feature_anchor_x"] = anchor_radius * math.cos(direction)
    label["feature_anchor_y"] = anchor_radius * math.sin(direction)
    label["middle_x"] = elbow_radius * math.cos(direction)
    label["middle_y"] = elbow_radius * math.sin(direction)
    return label


@pytest.mark.parametrize("inner", [False, True])
def test_radial_crowded_cluster_spreads_to_both_sides_of_its_features(inner: bool) -> None:
    """A cluster too crowded for its features' angles spreads around them, not one way only (OV-3648).

    Clockwise packing alone keeps the first label on its feature and pushes every other one clockwise, which slants
    the leaders of the whole cluster; whether it ran alone used to depend on whether its leaders crossed, so a small
    change of leader ends turned a cluster by several degrees.
    """

    angles = [90.0 + 0.25 * index for index in range(9)]
    labels = [_label_on_feature(index, angle, inner=inner) for index, angle in enumerate(angles)]
    result = place_radial_labels(
        () if inner else labels,
        labels if inner else (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
        inner_reserved_outer_radius_px=100.0,
    )
    assert result.collision_free
    offsets = {
        label["stable_id"]: ((label["placed_angle_deg"] - label["preferred_angle_deg"] + 180.0) % 360.0) - 180.0
        for label in result.labels
    }
    spread = max(offsets.values()) - min(offsets.values())
    assert spread > 1.0, offsets  # the cluster had to spread
    assert min(offsets.values()) < -0.25 * spread and max(offsets.values()) > 0.25 * spread, offsets
    assert abs(sum(offsets.values()) / len(offsets)) < 0.1 * spread, offsets


def test_long_inner_label_reports_required_radius_growth_without_truncation() -> None:
    inner = _label(0, 45.0, inner=True, width=500.0)
    result = place_radial_labels(
        (),
        (inner,),
        total_length=10_000,
        spacing_px=3.0,
        inner_reserved_outer_radius_px=100.0,
    )
    assert result.required_radius_growth_px > 0.0
    assert result.labels[0]["label_text"] == inner["label_text"]
    assert result.labels[0]["width_px"] == 500.0


def test_mode_label_config_defaults_and_circular_placement_validation() -> None:
    old_config = load_default_config()
    old_config["labels"].pop("circular")
    old_config["labels"].pop("linear")
    labels = apply_config_overrides(old_config, None).labels
    assert (labels.circular.scope, labels.circular.placement) == ("none", "horizontal")
    assert (labels.linear.scope, labels.linear.placement, labels.linear.rotation) == (
        "none",
        "auto",
        0.0,
    )
    assert apply_config_overrides(
        None,
        {"labels.circular.placement": "radial"},
    ).labels.circular.placement == "radial"
    with pytest.raises(ValidationError, match="horizontal.*radial"):
        apply_config_overrides(
            None,
            {"labels.circular.placement": "diagonal"},
        )


def test_circular_cli_placement_default_and_radial_parse() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    base_args = ["--gbk", str(input_path)]
    assert _get_args(base_args).label_placement == "horizontal"
    assert _get_args([*base_args, "--label_placement", "radial"]).label_placement == "radial"


def test_radial_solver_keeps_500_labels_without_collisions_or_cutoffs() -> None:
    label_count = 500
    total_length = 100_000
    labels = []
    for index in range(label_count):
        angle_deg = 360.0 * index / label_count
        angle_rad = math.radians(angle_deg - 90.0)
        labels.append(
            {
                "stable_id": f"feature-{index}",
                "input_order": index,
                "label_text": f"label {index}",
                "middle": angle_deg / 360.0 * total_length,
                "preferred_middle": angle_deg / 360.0 * total_length,
                "width_px": 22.0,
                "height_px": 8.0,
                "middle_x": 360.0 * math.cos(angle_rad),
                "middle_y": 360.0 * math.sin(angle_rad),
                "feature_anchor_x": 340.0 * math.cos(angle_rad),
                "feature_anchor_y": 340.0 * math.sin(angle_rad),
                "feature_inner_radius_px": 330.0,
                "feature_outer_radius_px": 350.0,
                "is_inner": False,
            }
        )

    result = place_radial_labels(
        labels,
        (),
        total_length=total_length,
        spacing_px=2.0,
        outer_reserved_radius_px=360.0,
    )

    assert len(result.labels) == label_count
    assert result.collision_free


def test_candidate_extraction_runs_once_across_inner_radius_preflight(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls = 0
    original = circular_labels_module.build_circular_label_candidates

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(circular_labels_module, "build_circular_label_candidates", counted)
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    args = _get_args(
        [
            "--gbk",
            str(input_path),
            "--labels",
            "both",
            "--label_placement",
            "radial",
            "--label_rendering",
            "external_only",
            "--legend",
            "none",
            # The full name is drawn over the GC rows on one line (DF-D25), and the definition rows of its warning
            # draw the figure again; the abbreviation keeps this to the one drawing whose preflight is counted.
            "--species",
            "MjeNMV",
            "--format",
            "svg",
            "--output",
            str(tmp_path / "radial-candidate-cache"),
        ]
    )

    run_circular_from_namespace(args)

    assert calls == 1


# The ribosomal protein cluster under ycf2 in the rotated tobacco plastome (docs recipe H-CLI-14, OV-3648): rpoA,
# rps11, rpl36, rps8, rpl14, rpl16, rps3, rpl22, rps19, rpl2, rpl23 and tRNA-Ile, by feature id (the inverted-repeat
# copies of rpl2, rpl23, rps19 and tRNA-Ile carry other ids).
_TOBACCO_SOUTH_FEATURE_IDS = {
    "fab2c2b1d", "f46c497f8", "f7cae9638", "fbaf469e4", "ff70975fb", "f2126f55f",
    "f7371cc8d", "f33c446c4", "f3d0f7434", "f2c0d3c71", "f4e9aae7d", "fc77059df",
}


def test_rotated_tobacco_plastome_keeps_its_south_inner_labels_over_their_features(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """H-CLI-14's inner labels below ycf2 spread around their features (OV-3648).

    The recipe places one feature outward, so inner leaders start on each feature's drawn outline; that removed a
    leader crossing from the clockwise packing, and the packing alone then turned the whole cluster up to 12 degrees
    one way, slanting every leader. The cluster spreads to both sides of its features instead.
    """

    results = []
    original = circular_labels_module.place_radial_labels

    def captured(*args, **kwargs):
        result = original(*args, **kwargs)
        results.append(result)
        return result

    monkeypatch.setattr(circular_labels_module, "place_radial_labels", captured)
    repository = Path(__file__).parents[1]
    placements = tmp_path / "placements.tsv"
    placements.write_text(
        "record\tfeature_selector\tplacement\tlevel\nNC_001879.2\tprotein_id=NP_054479.1\toutward\t1\n",
        encoding="utf-8",
    )
    priority = tmp_path / "qualifier_priority.tsv"
    priority.write_text("CDS\tgene,old_locus_tag\n", encoding="utf-8")
    args = _get_args(
        [
            "--gbk", str(repository / "tests/test_inputs/NC_001879.gbk"),
            "-t", "builtin:plastid",
            "-k", "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
            "--track_type", "tuckin",
            "--display_start_coordinate", "5500",
            "--feature_placement_table", str(placements),
            "--feature_overlap_tolerance_bp", "1",
            "--resolve_overlaps",
            "--gc",
            "--no-skew",
            "--labels", "both",
            "--label_placement", "radial",
            "--outer_label_x_radius_offset", "0.9",
            "--outer_label_y_radius_offset", "0.9",
            "--inner_label_x_radius_offset", "0.975",
            "--inner_label_y_radius_offset", "0.975",
            "--qualifier_priority", str(priority),
            "--annotation_table", str(repository / "examples/nicotiana-tabacum-regions.tsv"),
            "--circular_track_slot", "features:features@side=overlay,lane_direction=split",
            "--circular_track_slot",
            "plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,inner_gap_px=1,"
            "outer_gap_px=1,show_labels=true,padding_px=1,overflow=compress",
            "--circular_track_slot",
            "gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08,nt=GC,legend_label=GC content",
            "--legend", "none",
            "--format", "svg",
            "--output", str(tmp_path / "rotated_placed_chloroplast"),
        ]
    )

    run_circular_from_namespace(args)

    cluster = [
        label for label in results[-1].labels
        if label["feature_id"] in _TOBACCO_SOUTH_FEATURE_IDS
    ]
    assert len(cluster) == len(_TOBACCO_SOUTH_FEATURE_IDS)
    offsets = [
        ((float(label["placed_angle_deg"]) - float(label["preferred_angle_deg"]) + 180.0) % 360.0) - 180.0
        for label in cluster
    ]
    # Before the fix: every offset 0 to +12 degrees, mean +5.6. Spread around the features: -6 to +6, mean 0.
    assert min(offsets) < -2.0 and max(offsets) < 8.0, offsets
    assert abs(sum(offsets) / len(offsets)) < 2.0, offsets
