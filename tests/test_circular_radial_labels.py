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


def test_circular_placement_config_default_override_and_validation() -> None:
    old_config = load_default_config()
    old_config["labels"].pop("circular")
    assert apply_config_overrides(old_config, None).labels.circular.placement == "horizontal"
    assert apply_config_overrides(None, {"circular_label_placement": "radial"}).labels.circular.placement == "radial"
    with pytest.raises(ValidationError, match="horizontal.*radial"):
        apply_config_overrides(None, {"circular_label_placement": "diagonal"})


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
            "--format",
            "svg",
            "--output",
            str(tmp_path / "radial-candidate-cache"),
        ]
    )

    run_circular_from_namespace(args)

    assert calls == 1
