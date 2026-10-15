from __future__ import annotations

import json

import pytest

from gbdraw.api.config import load_default_config
from gbdraw.exceptions import ValidationError
from gbdraw.web_support.config_overrides import (
    validate_and_project_stored_request_config_overrides,
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


def test_full_typed_circular_config_ignores_inactive_linear_leaf() -> None:
    config = load_default_config()
    config["labels"]["linear"]["scope"] = "all"
    config["objects"]["gc_content"]["percent_background_opacity"] = 0.42

    projected = validate_and_project_web_config_overrides(
        mode="circular",
        config=config,
        managed_paths=["objects.scale.show"],
    )

    assert projected == {UNMANAGED_OPACITY_PATH: 0.42}


def test_explicit_circular_unmanaged_override_rejects_linear_label_scope() -> None:
    with pytest.raises(ValidationError, match="cannot target Linear label settings"):
        validate_and_project_web_config_overrides(
            mode="circular",
            overrides={"labels.linear.scope": "all"},
            managed_paths=[],
            require_unmanaged_only=True,
        )


def test_json_boundary_returns_validated_overlay() -> None:
    result = json.loads(validate_web_config_overrides_json(
        "linear",
        "null",
        json.dumps({"objects.blast_match.curve_tension": 0.25}),
        "[]",
        "storedRequest",
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
        "storedRequest",
    ))

    assert result["error"] == {
        "code": "INPUT_INVALID", "operation": "validateConfigOverrides", "stage": "helper",
        "context": {"field": "configOverrides", "reason": "UNKNOWN_CONFIG_PATH"}
    }
    assert "objects.gc_content.unknown" not in json.dumps(result)
    assert "biological secret" not in json.dumps(result)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_gui_filter_leaf_does_not_create_a_duplicate_raw_override(mode):
    projected = validate_and_project_web_config_overrides(
        mode=mode,
        overrides={"labels.filtering.blacklist_keywords": [], UNMANAGED_OPACITY_PATH: 0.42},
        managed_paths=["labels.filtering.blacklist_keywords"],
    )
    assert projected == {UNMANAGED_OPACITY_PATH: 0.42}


def test_explicit_raw_filtering_and_unmanaged_filter_changes_remain_preserved():
    raw = {"blacklist_keywords": [], "qualifier_priority": {"gene": ["gene"]}}
    projected = validate_and_project_web_config_overrides(
        mode="circular", overrides={"labels.filtering.raw": raw},
        managed_paths=["labels.filtering.blacklist_keywords"],
    )
    assert projected["labels.filtering.raw"] == raw
    config = load_default_config()
    config["labels"]["filtering"] = raw
    projected = validate_and_project_web_config_overrides(
        mode="circular", config=config, managed_paths=["labels.filtering.blacklist_keywords"],
    )
    assert projected["labels.filtering.raw"] == raw


def _legacy_label_config() -> dict:
    """A stored full config of a request schema 1 or 2 (Sessions 31-33, OV-790)."""
    config = load_default_config()
    config["canvas"]["show_labels"] = True
    config["canvas"]["circular"]["allow_inner_labels"] = True
    config["canvas"]["linear"]["track_layout"] = "spreadout"
    config["labels"]["circular"].pop("scope")
    config["labels"]["linear"].pop("scope")
    config["labels"]["linear"]["placement"] = "on_feature"
    return config


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_stored_full_config_is_read_as_the_request_decoder_reads_it(mode):
    # The Load check reads a stored config with the decoder's migration of
    # retired config keys (OV-790), as decode_canonical_request does.
    from gbdraw.session_request_codec import migrate_legacy_full_config

    legacy = _legacy_label_config()
    migrated = migrate_legacy_full_config(legacy, mode=mode, path="diagramOptions.config")
    assert "show_labels" not in migrated["canvas"]
    projected = validate_and_project_stored_request_config_overrides(mode=mode, config=legacy, managed_paths=[])
    assert projected == validate_and_project_web_config_overrides(
        mode=mode, config=migrated, managed_paths=[]
    )
    assert projected[f"labels.{mode}.scope"] == ("both" if mode == "circular" else "all")
    assert legacy["canvas"]["show_labels"] is True


def test_stored_full_config_with_conflicting_label_values_is_refused():
    legacy = _legacy_label_config()
    legacy["canvas"]["show_labels"] = False
    legacy["canvas"]["circular"]["show_labels"] = True
    result = json.loads(validate_web_config_overrides_json("circular", json.dumps(legacy), "null", "[]", "storedRequest"))
    # The decoder's label-conflict check refuses it; its error names no field.
    assert result["error"] == {
        "code": "VALIDATION_UNCLASSIFIED", "operation": "validateConfigOverrides", "stage": "helper", "context": {},
    }
    with pytest.raises(ValidationError, match="diagramOptions.config.canvas contains conflicting Circular label values"):
        validate_and_project_stored_request_config_overrides(mode="circular", config=legacy, managed_paths=[])


@pytest.mark.parametrize(
    ("mode", "flat", "canonical"),
    [
        ("circular", {"show_gc": False}, {"canvas.show_gc": False}),
        ("linear", {"linear_axis_stroke_color": "black"}, {"objects.axis.linear.stroke_color": "black"}),
        ("circular", {"gc_content_percent_border_width": 2}, {"objects.gc_content.percent_border_width": 2}),
    ],
)
def test_flat_override_names_are_read_as_the_request_decoder_reads_them(mode, flat, canonical):
    # The Load check reads persisted flat override names with the decoder's
    # migration (OV-799), as decode_canonical_request does at every schema.
    assert validate_and_project_stored_request_config_overrides(mode=mode, overrides=flat, managed_paths=[]) == (
        validate_and_project_web_config_overrides(mode=mode, overrides=canonical, managed_paths=[])
    )
    with pytest.raises(ValidationError, match="cannot contain both"):
        validate_and_project_stored_request_config_overrides(
            mode=mode, overrides={**flat, **{path: "other" for path in canonical}}, managed_paths=[],
        )


@pytest.mark.parametrize(
    ("flat", "managed"),
    [({"gc_stroke_color": "red"}, []), ({"label_blacklist": "foo,bar"}, ["labels.filtering.blacklist_keywords"])],
)
def test_a_draft_check_reads_no_flat_override_name(flat, managed):
    # Only a stored request's check reads flat names (OV-799, D-OD11); the
    # check of a draft's unmanaged overrides keeps its reading and its error.
    draft = json.loads(validate_web_config_overrides_json(
        "circular", "null", json.dumps(flat), json.dumps(managed), "draft",
    ))
    assert draft["error"] == {
        "code": "INPUT_INVALID", "operation": "validateConfigOverrides", "stage": "helper",
        "context": {"field": "configOverrides", "reason": "UNKNOWN_CONFIG_PATH"},
    }
    with pytest.raises(ValidationError, match="Unknown config override path"):
        validate_and_project_web_config_overrides(mode="circular", overrides=flat, managed_paths=managed)
    stored = json.loads(validate_web_config_overrides_json(
        "circular", "null", json.dumps({"gc_stroke_color": "red"}), "[]", "storedRequest"))
    assert stored == {"overrides": {"objects.gc_content.stroke_color": "red"}}


def test_json_boundary_requires_its_check():
    # `source` has no default at either end (D-OD16 N2): a call without it is a programming error.
    with pytest.raises(TypeError):
        validate_web_config_overrides_json("circular", "null", "{}", "[]")


@pytest.mark.parametrize("source", [True, False, "", "stored", None])
def test_json_boundary_names_its_check(source):
    # One field names the check (D-OD12): a value that names neither is refused.
    result = json.loads(validate_web_config_overrides_json("circular", "null", "{}", "[]", source))
    assert result["error"]["code"] == "HELPER_PROTOCOL"
    named = json.loads(validate_web_config_overrides_json(
        "circular", "null", json.dumps({"gc_stroke_color": "red"}), "[]", "storedRequest"))
    assert named == {"overrides": {"objects.gc_content.stroke_color": "red"}}
