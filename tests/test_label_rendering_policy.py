from __future__ import annotations

import math

import pytest

import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.config.models import (
    CircularRenderProfile,
    GbdrawConfig,
    LinearRenderProfile,
)
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.text import calculate_svg_bbox_dimensions
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.labels.circular import prepare_label_list
from gbdraw.labels.linear import linear_record_labels, place_linear_label_row
from gbdraw.linear import _get_args as get_linear_args

circular_assemble_module


def _feature(feature_id: str, start: int, end: int, label: str, strand: str = "positive") -> FeatureObject:
    location = [FeatureLocationPart("block", "001", strand, start, end, True)]
    return FeatureObject(
        feature_id=feature_id,
        location=location,
        is_directional=True,
        color="#54bcf8",
        note="",
        label_text=label,
        coordinates=location,
        type="CDS",
        qualifiers={"product": [label]},
    )


def _config(label_rendering: str) -> GbdrawConfig:
    config_dict = modify_config_dict(
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "labels.rendering": label_rendering,
            "labels.filtering.blacklist_keywords": [],
            "labels.circular.scope": "outer",
            "labels.linear.scope": "all",
        },
    )
    return GbdrawConfig.from_dict(config_dict)


def _linear_policy_labels(label_rendering: str) -> list[dict]:
    cfg = _config(label_rendering)
    features = {
        "fits": _feature("fits", 1, 300, "A"),
        "external": _feature("external", 500, 510, "very long label that cannot fit inside a ten base feature"),
    }
    profile = LinearRenderProfile(cfg)
    record_labels = linear_record_labels(features, 1000, 1000, 1.0, 20, False, "middle", None, profile)
    return place_linear_label_row([(record_labels, 0.0)], profile=profile, track_layout="middle")[0]


def _circular_policy_labels(label_rendering: str) -> list[dict]:
    cfg = _config(label_rendering)
    features = {
        "fits": _feature("fits", 1, 400, "A"),
        "external": _feature("external", 500, 510, "very long label that cannot fit inside a ten base feature"),
    }
    return prepare_label_list(
        features,
        1000,
        390,
        0.19,
        CircularRenderProfile(cfg),
    )


@pytest.mark.circular
def test_circular_label_rendering_policy_filters_embedding_classes() -> None:
    auto_labels = _circular_policy_labels("auto")
    embedded_only_labels = _circular_policy_labels("embedded_only")
    external_only_labels = _circular_policy_labels("external_only")

    assert any(label["is_embedded"] for label in auto_labels)
    assert any(not label["is_embedded"] for label in auto_labels)
    assert embedded_only_labels
    assert all(label["is_embedded"] for label in embedded_only_labels)
    assert external_only_labels
    assert all(not label["is_embedded"] for label in external_only_labels)


@pytest.mark.linear
def test_linear_label_rendering_policy_filters_embedding_classes() -> None:
    auto_labels = _linear_policy_labels("auto")
    embedded_only_labels = _linear_policy_labels("embedded_only")
    external_only_labels = _linear_policy_labels("external_only")

    assert any(label["is_embedded"] for label in auto_labels)
    assert any(not label["is_embedded"] for label in auto_labels)
    assert embedded_only_labels
    assert all(label["is_embedded"] for label in embedded_only_labels)
    assert external_only_labels
    assert all(not label["is_embedded"] for label in external_only_labels)


@pytest.mark.linear
def test_linear_label_goes_inside_its_feature_only_with_a_quarter_em_to_spare() -> None:
    """A Linear label is drawn inside its feature only when its advance plus 0.25 em fits (Owner D-20)."""
    cfg = _config("auto")
    font_size = 16.0
    label = "dnaK"
    advance_px, _ = calculate_svg_bbox_dimensions(label, cfg.objects.text.font_family, font_size, cfg.canvas.dpi)
    margin_px = 0.25 * font_size
    # One px per bp here; a feature of n bp is n px long.
    tight_bp = math.floor(advance_px) + 2
    roomy_bp = math.ceil(advance_px + margin_px) + 1
    assert advance_px < tight_bp < advance_px + margin_px < roomy_bp
    features = {
        "tight": _feature("tight", 1, tight_bp, label),
        "roomy": _feature("roomy", 500, 500 + roomy_bp - 1, label),
    }
    profile = LinearRenderProfile(cfg)
    record_labels = linear_record_labels(
        features, 1000, 1000, 1.0, 20, False, "middle", None, profile, label_font_size=font_size
    )
    labels = place_linear_label_row([(record_labels, 0.0)], profile=profile, track_layout="middle")[0]

    embedded_by_start = {round(entry["feature_start_x"]): entry["is_embedded"] for entry in labels}
    assert embedded_by_start == {1: False, 500: True}


@pytest.mark.linear
def test_linear_label_rendering_policy_rejects_above_feature_conflict() -> None:
    with pytest.raises(SystemExit):
        get_linear_args(
            [
                "--gbk",
                "examples/MjeNMV.gb",
                "--label_placement",
                "above_feature",
                "--label_rendering",
                "external_only",
            ]
        )
