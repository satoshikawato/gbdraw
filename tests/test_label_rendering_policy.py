from __future__ import annotations

import pytest

import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.labels.circular import prepare_label_list
from gbdraw.labels.linear import prepare_label_list_linear
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


def _config(label_rendering: str) -> tuple[dict, GbdrawConfig]:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        label_rendering=label_rendering,
        label_blacklist="",
        show_labels=True,
    )
    return config_dict, GbdrawConfig.from_dict(config_dict)


def _linear_policy_labels(label_rendering: str) -> list[dict]:
    config_dict, cfg = _config(label_rendering)
    features = {
        "fits": _feature("fits", 1, 300, "A"),
        "external": _feature("external", 500, 510, "very long label that cannot fit inside a ten base feature"),
    }
    return prepare_label_list_linear(
        features,
        1000,
        1000,
        1.0,
        20,
        False,
        "middle",
        None,
        config_dict,
        cfg=cfg,
    )


def _circular_policy_labels(label_rendering: str) -> list[dict]:
    config_dict, cfg = _config(label_rendering)
    features = {
        "fits": _feature("fits", 1, 400, "A"),
        "external": _feature("external", 500, 510, "very long label that cannot fit inside a ten base feature"),
    }
    return prepare_label_list(
        features,
        1000,
        390,
        0.19,
        config_dict,
        cfg=cfg,
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
def test_linear_label_rendering_policy_rejects_above_feature_conflict() -> None:
    with pytest.raises(SystemExit):
        get_linear_args(
            [
                "--gbk",
                "tests/test_inputs/MjeNMV.gb",
                "--label_placement",
                "above_feature",
                "--label_rendering",
                "external_only",
            ]
        )
