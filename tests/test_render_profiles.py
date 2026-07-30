from dataclasses import FrozenInstanceError, replace

import pytest

from gbdraw.api.options import (
    CircularDiagramOptions,
    LinearDiagramOptions,
    resolve_circular_diagram_options,
)
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from gbdraw.config.models import (
    CircularRenderProfile,
    GbdrawConfig,
    LinearRenderProfile,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from tests.utils.linear_render_context import make_linear_record_render_context


def _config(
    *,
    circular_scope: object = "none",
    linear_scope: object = "none",
) -> GbdrawConfig:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["labels"]["circular"]["scope"] = circular_scope
    raw["labels"]["linear"]["scope"] = linear_scope
    return GbdrawConfig.from_dict(raw)


@pytest.mark.parametrize(
    ("circular", "linear"),
    [
        ("outer", "all"),
        ("none", "none"),
        ("both", "first"),
        ("outer", "orthogroup_top"),
    ],
)
def test_render_profiles_read_independent_mode_specific_label_settings(
    circular: str,
    linear: str,
) -> None:
    cfg = _config(
        circular_scope=circular,
        linear_scope=linear,
    )

    assert CircularRenderProfile(cfg).label_scope == circular
    assert LinearRenderProfile(cfg).label_scope == linear


@pytest.mark.parametrize("setting", ["first", "orthogroup_top"])
def test_config_rejects_linear_label_policies_in_circular_schema(
    setting: str,
) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["labels"]["circular"]["scope"] = setting

    with pytest.raises(ValidationError, match="labels.circular.scope"):
        GbdrawConfig.from_dict(raw)


@pytest.mark.parametrize("setting", ["first", "orthogroup_top", "sometimes"])
def test_circular_option_resolution_rejects_linear_label_policies(
    setting: str,
) -> None:
    with pytest.raises(ValidationError, match="labels.circular.scope"):
        resolve_circular_diagram_options(
            CircularDiagramOptions(
                config_overrides={"labels.circular.scope": setting},
            )
        )


def test_circular_options_reject_linear_label_scope_override() -> None:
    with pytest.raises(
        ValidationError,
        match="Circular config overrides cannot target Linear label settings",
    ):
        CircularDiagramOptions(
            config_overrides={"labels.linear.scope": "first"},
        )


def test_linear_options_reject_circular_label_scope_override() -> None:
    with pytest.raises(
        ValidationError,
        match="Linear config overrides cannot target Circular label settings",
    ):
        LinearDiagramOptions(
            config_overrides={"labels.circular.scope": "outer"},
        )


def test_circular_profile_never_interprets_linear_label_policy() -> None:
    cfg = _config(
        circular_scope="none",
        linear_scope="orthogroup_top",
    )

    assert CircularRenderProfile(cfg).label_scope == "none"
    assert LinearRenderProfile(cfg).label_scope == "orthogroup_top"


def test_render_profiles_capture_active_canvas_values() -> None:
    cfg = _config()
    cfg = replace(
        cfg,
        canvas=replace(
            cfg.canvas,
            show_gc=True,
            show_skew=False,
            show_depth=True,
            strandedness=False,
            resolve_overlaps=True,
        ),
    )

    profile = LinearRenderProfile(cfg)

    assert profile.config is cfg
    assert profile.show_gc is True
    assert profile.show_skew is False
    assert profile.show_depth is True
    assert profile.strandedness is False
    assert profile.resolve_overlaps is True


def test_linear_canvas_active_values_are_read_only_profile_projections() -> None:
    profile = LinearRenderProfile(_config(linear_scope="all"))
    canvas = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=100,
        profile=profile,
        legend="none",
    )

    assert canvas.profile is profile
    for field_name in (
        "show_gc",
        "show_skew",
        "show_depth",
        "strandedness",
        "resolve_overlaps",
        "track_layout",
    ):
        assert getattr(canvas, field_name) == getattr(profile, field_name)
        with pytest.raises(AttributeError):
            setattr(canvas, field_name, None)
    assert not hasattr(canvas, "show_labels")


def test_circular_canvas_active_values_are_read_only_profile_projections() -> None:
    profile = CircularRenderProfile(_config(circular_scope="outer"))
    canvas = CircularCanvasConfigurator(
        output_prefix="test",
        profile=profile,
        legend="none",
        gb_record=type("_Record", (), {"seq": "A" * 100})(),
    )

    assert canvas.profile is profile
    for field_name in (
        "show_gc",
        "show_skew",
        "show_depth",
        "strandedness",
        "resolve_overlaps",
    ):
        assert getattr(canvas, field_name) == getattr(profile, field_name)
        with pytest.raises(AttributeError):
            setattr(canvas, field_name, None)
    assert not hasattr(canvas, "show_labels")


def test_linear_record_render_context_owns_resolved_active_layout() -> None:
    cfg = _config(linear_scope="all")
    cfg = replace(cfg, canvas=replace(cfg.canvas, show_depth=True))
    profile = LinearRenderProfile(cfg)
    context = make_linear_record_render_context(profile, depth_available=True)

    assert context.profile is profile
    assert context.feature_track_layout == profile.track_layout
    assert context.depth_enabled is True
    assert isinstance(context.track_layout.slots, tuple)
    assert isinstance(context.track_layout.depth_track_offsets, tuple)
    with pytest.raises(FrozenInstanceError):
        context.depth_enabled = False  # type: ignore[misc]
    with pytest.raises(TypeError):
        context.track_layout.slots[0].params["mutable"] = True  # type: ignore[index]
    with pytest.raises(TypeError):
        context.track_layout.slots[0].params["nested"]["mutable"] = True  # type: ignore[index]
    assert context.track_layout.slots[0].params["nested"]["values"] == (1,)


def test_render_profiles_are_frozen() -> None:
    profile = CircularRenderProfile(_config())

    with pytest.raises(FrozenInstanceError):
        profile.show_gc = False  # type: ignore[misc]


@pytest.mark.parametrize(
    ("section", "key", "message"),
    [
        (None, "show_labels", "canvas.show_labels"),
        ("circular", "show_labels", "canvas.circular.show_labels"),
        ("circular", "allow_inner_labels", "canvas.circular.allow_inner_labels"),
        ("linear", "show_labels", "canvas.linear.show_labels"),
    ],
)
def test_config_rejects_retired_canvas_label_settings(
    section: str | None,
    key: str,
    message: str,
) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    target = raw["canvas"] if section is None else raw["canvas"][section]
    target[key] = True

    with pytest.raises(ValidationError, match=message):
        GbdrawConfig.from_dict(raw)


@pytest.mark.parametrize("setting", [True, False, "sometimes"])
def test_config_rejects_invalid_linear_label_scope(setting: object) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["labels"]["linear"]["scope"] = setting

    with pytest.raises(ValidationError, match="labels.linear.scope"):
        GbdrawConfig.from_dict(raw)


@pytest.mark.parametrize(
    ("section", "key", "value", "message"),
    [
        ("linear", "placement", "on_feature", "labels.linear.placement"),
        ("linear", "placement", "garbage", "labels.linear.placement"),
        ("circular", "placement", None, "circular label placement"),
    ],
)
def test_config_rejects_invalid_label_placement(
    section: str,
    key: str,
    value: object,
    message: str,
) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["labels"][section][key] = value

    with pytest.raises(ValidationError, match=message):
        GbdrawConfig.from_dict(raw)


@pytest.mark.parametrize("value", [None, "garbage", True])
def test_config_rejects_invalid_label_rendering(value: object) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["labels"]["rendering"] = value

    with pytest.raises(ValidationError, match="labels.rendering"):
        GbdrawConfig.from_dict(raw)
