from __future__ import annotations

import inspect

import pandas as pd
import pytest

from gbdraw.api.config import apply_config_overrides, load_default_config
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import (
    canonical_config_override_paths,
    modify_config_dict,
)
from gbdraw.exceptions import ValidationError


def test_override_paths_are_canonical_dataclass_leaves() -> None:
    paths = canonical_config_override_paths()

    assert "canvas.show_gc" in paths
    assert "objects.scale.show" in paths
    assert "labels.circular.scope" in paths
    assert "labels.linear.scope" in paths
    assert "objects.axis.linear.stroke_color" in paths
    assert "objects.features.arrow_geometry.head_length_ratio" in paths
    assert "objects.features.arrow_geometry.shaft_width_ratio" in paths
    assert "objects.definition.linear.line_styles.name.font_size" in paths
    assert "canvas" not in paths
    assert "objects.axis.linear" not in paths
    assert "show_gc" not in paths
    assert set(inspect.signature(modify_config_dict).parameters) == {
        "config_dict",
        "overrides",
    }


@pytest.mark.parametrize(
    "path",
    [
        "show_gc",
        "canvas.show_labels",
        "canvas.circular.show_labels",
        "canvas.linear.show_labels",
        "canvas.unknown",
    ],
)
def test_unknown_override_path_is_rejected(path: str) -> None:
    with pytest.raises(ValidationError, match="Unknown config override path"):
        modify_config_dict(load_default_config(), {path: False})


@pytest.mark.parametrize("path", ["canvas", "objects.axis.linear"])
def test_non_leaf_override_path_is_rejected(path: str) -> None:
    with pytest.raises(ValidationError, match="not a leaf field"):
        modify_config_dict(load_default_config(), {path: {}})


def test_canonical_override_preserves_unmodeled_raw_metadata() -> None:
    config = load_default_config()

    modified = modify_config_dict(
        config,
        {
            "canvas.show_gc": False,
            "objects.axis.linear.stroke_color": 'dimgray',
        },
    )

    assert modified["canvas"]["show_gc"] is False
    assert modified["objects"]["axis"]["linear"]["stroke_color"] == "dimgray"
    assert modified["title"] == config["title"]
    assert modified["png_output"] == config["png_output"]


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("canvas.linear.track_layout", "spreadout"),
        ("labels.linear.placement", "on_feature"),
        ("labels.circular.scope", "first"),
        ("labels.linear.scope", "outer"),
    ],
)
def test_literal_domains_reject_reader_only_values(path: str, value: str) -> None:
    with pytest.raises(ValidationError, match="Invalid value"):
        modify_config_dict(load_default_config(), {path: value})


def test_typed_override_preserves_filtering_extensions_and_dataframes() -> None:
    config = load_default_config()
    whitelist = pd.DataFrame(
        [["CDS", "product", "polymerase"]],
        columns=["feature_type", "qualifier", "keyword"],
    )
    filtering = config["labels"]["filtering"]
    filtering["whitelist_df"] = whitelist
    filtering["extension_key"] = {"keep": True}

    updated = apply_config_overrides(
        GbdrawConfig.from_dict(config),
        {"labels.circular.scope": "outer"},
    )

    updated_filtering = updated.labels.filtering.as_dict()
    pd.testing.assert_frame_equal(updated_filtering["whitelist_df"], whitelist)
    assert updated_filtering["extension_key"] == {"keep": True}
    assert updated.labels.circular.scope == "outer"


def test_filtering_raw_override_is_applied_before_typed_siblings() -> None:
    updated = apply_config_overrides(
        None,
        {
            "labels.filtering.blacklist_keywords": ["specific"],
            "labels.filtering.raw": {
                "blacklist_keywords": ["raw"],
                "qualifier_priority": {},
                "extension_key": {"keep": True},
            },
        },
    )

    filtering = updated.labels.filtering.as_dict()
    assert filtering["blacklist_keywords"] == ["specific"]
    assert filtering["extension_key"] == {"keep": True}


def test_default_override_resolution_parses_typed_config_once(monkeypatch) -> None:
    original_from_dict = GbdrawConfig.from_dict
    calls = 0

    def counting_from_dict(
        cls: type[GbdrawConfig],
        config_dict: dict,
    ) -> GbdrawConfig:
        nonlocal calls
        calls += 1
        return original_from_dict(config_dict)

    monkeypatch.setattr(
        GbdrawConfig,
        "from_dict",
        classmethod(counting_from_dict),
    )

    updated = apply_config_overrides(
        None,
        {"labels.circular.scope": "outer"},
    )

    assert updated.labels.circular.scope == "outer"
    assert calls == 1


def test_typed_options_reject_flat_override_aliases_eagerly() -> None:
    with pytest.raises(ValidationError, match="Unknown config override path"):
        CircularDiagramOptions(config_overrides={"show_gc": False})


def test_scale_visibility_defaults_visible_and_normalizes_raw_config_values() -> None:
    default_config = load_default_config()
    assert default_config["objects"]["scale"]["show"] is True
    assert GbdrawConfig.from_dict(default_config).objects.scale.show is True

    missing_config = load_default_config()
    del missing_config["objects"]["scale"]["show"]
    assert GbdrawConfig.from_dict(missing_config).objects.scale.show is True

    compatible_raw_config = load_default_config()
    compatible_raw_config["objects"]["scale"]["show"] = "off"
    assert GbdrawConfig.from_dict(compatible_raw_config).objects.scale.show is False


@pytest.mark.parametrize(
    "options_type",
    [CircularDiagramOptions, LinearDiagramOptions],
)
def test_scale_visibility_is_a_shared_strict_boolean_override(
    options_type: type[CircularDiagramOptions] | type[LinearDiagramOptions],
) -> None:
    options = options_type(config_overrides={"objects.scale.show": False})
    assert options.config_overrides["objects.scale.show"] is False

    with pytest.raises(ValidationError, match="objects.scale.show"):
        options_type(config_overrides={"objects.scale.show": 0})


def test_scale_style_rejects_unknown_raw_and_override_values() -> None:
    invalid_raw_config = load_default_config()
    invalid_raw_config["objects"]["scale"]["style"] = "unknown"
    with pytest.raises(ValidationError, match="scale_style"):
        GbdrawConfig.from_dict(invalid_raw_config)

    with pytest.raises(ValidationError, match="allowed literal values"):
        modify_config_dict(
            load_default_config(),
            {"objects.scale.style": "unknown"},
        )


def test_arrow_geometry_defaults_and_missing_block_compatibility() -> None:
    default_config = load_default_config()
    assert default_config["objects"]["features"]["arrow_geometry"] == {
        "head_length_ratio": "auto",
        "shaft_width_ratio": 1.0,
    }

    typed = GbdrawConfig.from_dict(default_config)
    assert typed.objects.features.arrow_geometry.head_length_ratio == "auto"
    assert typed.objects.features.arrow_geometry.shaft_width_ratio == 1.0

    del default_config["objects"]["features"]["arrow_geometry"]
    compatible = GbdrawConfig.from_dict(default_config)
    assert compatible.objects.features.arrow_geometry.head_length_ratio == "auto"
    assert compatible.objects.features.arrow_geometry.shaft_width_ratio == 1.0


@pytest.mark.parametrize(
    ("path", "value", "expected"),
    [
        ("objects.features.arrow_geometry.head_length_ratio", "auto", "auto"),
        ("objects.features.arrow_geometry.head_length_ratio", 1, 1.0),
        ("objects.features.arrow_geometry.head_length_ratio", 1.25, 1.25),
        ("objects.features.arrow_geometry.shaft_width_ratio", 0.25, 0.25),
        ("objects.features.arrow_geometry.shaft_width_ratio", 1, 1.0),
    ],
)
def test_arrow_geometry_dotted_overrides(path: str, value: object, expected: object) -> None:
    updated = modify_config_dict(load_default_config(), {path: value})
    geometry = GbdrawConfig.from_dict(updated).objects.features.arrow_geometry
    field_name = path.rsplit(".", 1)[-1]
    assert getattr(geometry, field_name) == expected


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("objects.features.arrow_geometry.head_length_ratio", True),
        ("objects.features.arrow_geometry.head_length_ratio", 0),
        ("objects.features.arrow_geometry.head_length_ratio", -0.5),
        ("objects.features.arrow_geometry.head_length_ratio", float("nan")),
        ("objects.features.arrow_geometry.head_length_ratio", float("inf")),
        ("objects.features.arrow_geometry.shaft_width_ratio", True),
        ("objects.features.arrow_geometry.shaft_width_ratio", 0),
        ("objects.features.arrow_geometry.shaft_width_ratio", -0.5),
        ("objects.features.arrow_geometry.shaft_width_ratio", 1.01),
        ("objects.features.arrow_geometry.shaft_width_ratio", float("nan")),
        ("objects.features.arrow_geometry.shaft_width_ratio", float("inf")),
    ],
)
def test_arrow_geometry_rejects_invalid_overrides(path: str, value: object) -> None:
    with pytest.raises(ValidationError, match="arrow|Invalid value"):
        modify_config_dict(load_default_config(), {path: value})
