from __future__ import annotations

import json

import pytest

from gbdraw.api.config import load_default_config
from gbdraw.exceptions import ValidationError
from gbdraw.web_support.config_overrides import (
    validate_and_project_web_config_overrides,
    validate_web_config_overrides_json,
)


UNMANAGED_OPACITY_PATH = "objects.gc_content.percent_background_opacity"


def test_known_typed_unmanaged_leaf_is_preserved_exactly() -> None:
    projected = validate_and_project_web_config_overrides(
        mode="circular",
        overrides={UNMANAGED_OPACITY_PATH: 0.42},
        managed_paths=["objects.scale.show"],
    )

    assert projected == {UNMANAGED_OPACITY_PATH: 0.42}


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"objects.gc_content.unknown": 1}, "Unknown config override path"),
        ({"objects.gc_content": {}}, "not a leaf field"),
        (
            {"objects.__proto__.percent_background_opacity": 0.42},
            "unsafe key",
        ),
        (
            {"labels.filtering.raw": {"constructor": {"leak": True}}},
            "unsafe key",
        ),
        ({UNMANAGED_OPACITY_PATH: "0.42"}, "Invalid value"),
        ({UNMANAGED_OPACITY_PATH: 1.1}, "finite number in.*0.0, 1.0"),
        ({"objects.scale.style": "unknown"}, "allowed literal values"),
    ],
)
def test_invalid_overlay_inputs_are_rejected_explicitly(
    overrides: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValidationError, match=message):
        validate_and_project_web_config_overrides(
            mode="circular",
            overrides=overrides,
            managed_paths=[],
        )


@pytest.mark.parametrize(
    ("mode", "path"),
    [
        ("circular", "objects.blast_match.curve_tension"),
        ("linear", "objects.axis.circular.stroke_color"),
    ],
)
def test_active_mode_incompatible_leaf_is_rejected(mode: str, path: str) -> None:
    with pytest.raises(ValidationError, match="cannot target"):
        validate_and_project_web_config_overrides(
            mode=mode,  # type: ignore[arg-type]
            overrides={path: 0.25 if path.endswith("curve_tension") else "black"},
            managed_paths=[],
        )


def test_gui_managed_paths_are_filtered_and_cannot_be_stored_as_unmanaged() -> None:
    projected = validate_and_project_web_config_overrides(
        mode="circular",
        overrides={
            "objects.scale.show": False,
            UNMANAGED_OPACITY_PATH: 0.42,
        },
        managed_paths=["objects.scale.show"],
    )
    assert projected == {UNMANAGED_OPACITY_PATH: 0.42}

    with pytest.raises(ValidationError, match="GUI-managed.*objects.scale.show"):
        validate_and_project_web_config_overrides(
            mode="circular",
            overrides={"objects.scale.show": False},
            managed_paths=["objects.scale.show"],
            require_unmanaged_only=True,
        )


def test_full_typed_config_projects_only_changed_unmanaged_leaves() -> None:
    config = load_default_config()
    config["objects"]["gc_content"]["percent_background_opacity"] = 0.42

    projected = validate_and_project_web_config_overrides(
        mode="circular",
        config=config,
        managed_paths=["objects.scale.show"],
    )

    assert projected == {UNMANAGED_OPACITY_PATH: 0.42}


def test_json_boundary_returns_validated_overlay() -> None:
    result = json.loads(validate_web_config_overrides_json(
        "linear",
        "null",
        json.dumps({"objects.blast_match.curve_tension": 0.25}),
        "[]",
    ))

    assert result == {
        "overrides": {"objects.blast_match.curve_tension": 0.25}
    }


def test_json_boundary_returns_actionable_safe_rejection() -> None:
    result = json.loads(validate_web_config_overrides_json(
        "circular",
        "null",
        json.dumps({"objects.gc_content.unknown": "biological secret"}),
        "[]",
    ))

    assert "objects.gc_content.unknown" in result["error"]
    assert "biological secret" not in result["error"]
