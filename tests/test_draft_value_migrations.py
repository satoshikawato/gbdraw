"""The Python twins of the Web Load value migrations of a Session 27-44 draft.

``tools/generate_draft_value_migration_vectors.mjs`` runs the Web functions of
``services/config.js`` and writes the vectors; this file reads them.
"""

from __future__ import annotations

import copy
import gzip
import json
import subprocess
from pathlib import Path
from typing import Any, Callable

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session import SessionFormatError
from gbdraw.session_io import (
    _holds_cli_writer_config,
    _with_legacy_repeat_region_shape,
    _without_legacy_null_circular_slot_spacing,
    migrate_imported_linear_track_slots,
    migrate_persisted_web_option_values,
    migrate_session_draft_values,
    migrate_session_flat_draft,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures"
FIELDS_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


def _cases(name: str) -> list[dict[str, Any]]:
    vectors = json.loads((FIXTURES / name).read_text(encoding="utf-8"))
    assert vectors["schemaVersion"] == 1
    return vectors["cases"]


def _fixture_config(path: str) -> Any:
    data = (REPO_ROOT / path).read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)["config"]


def _check(case: dict[str, Any], step: Callable[[Any], Any]) -> None:
    source = _fixture_config(case["fixture"]) if "fixture" in case else copy.deepcopy(case["input"]["config"])
    before = copy.deepcopy(source)
    if "error" in case:
        with pytest.raises(ValidationError) as excinfo:
            step(source)
        assert str(excinfo.value) == case["error"]
        assert excinfo.value.diagnostic == FIELDS_INVALID
    else:
        assert json.dumps(step(source)) == json.dumps(case["expected"]["config"])
    assert source == before


def test_vectors_are_the_web_functions_output() -> None:
    result = subprocess.run(
        ["node", "tools/generate_draft_value_migration_vectors.mjs", "--check"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


_OPTION_CASES = _cases("draft-option-value-migration-vectors.json")
_SLOT_CASES = _cases("linear-track-slot-migration-vectors.json")
_SHAPE_CASES = _cases("session-33-draft-shape-migration-vectors.json")
_SHAPE_STEPS = {
    "circular-slot-spacing": _without_legacy_null_circular_slot_spacing,
    "repeat-region-shape": _with_legacy_repeat_region_shape,
}


@pytest.mark.parametrize("case", _OPTION_CASES, ids=[case["name"] for case in _OPTION_CASES])
def test_option_values_migrate_as_the_web_does(case: dict[str, Any]) -> None:
    _check(case, migrate_persisted_web_option_values)


@pytest.mark.parametrize("case", _SLOT_CASES, ids=[case["name"] for case in _SLOT_CASES])
def test_linear_track_slots_migrate_as_the_web_does(case: dict[str, Any]) -> None:
    _check(case, lambda config: migrate_imported_linear_track_slots(config, case["sessionVersion"]))


@pytest.mark.parametrize("case", _SHAPE_CASES, ids=[case["name"] for case in _SHAPE_CASES])
def test_session_33_shapes_migrate_as_the_web_does(case: dict[str, Any]) -> None:
    _check(case, _SHAPE_STEPS[case["step"]])


def test_draft_values_run_in_the_web_load_order() -> None:
    config = _fixture_config("tests/fixtures/sessions/BGC0000708-BGC0000713.v30.gbdraw-session.json.gz")
    config["form"]["linear_track_layout"] = "spreadout"

    expected = _with_legacy_repeat_region_shape(
        migrate_imported_linear_track_slots(
            _without_legacy_null_circular_slot_spacing(migrate_persisted_web_option_values(config)), 30
        )
    )
    assert migrate_session_draft_values(config, 30) == expected
    assert expected["form"]["linear_track_layout"] == "above"
    # From Session 40 the option values and the Session 27-33 shapes are current.
    assert migrate_session_draft_values({"form": {"linear_track_layout": "spreadout"}}, 40) == {
        "form": {"linear_track_layout": "spreadout"}
    }


def test_session_40_44_stored_colors_load_as_saved() -> None:
    # restoreCurrentWriterActiveConfig drops colorsAreOverrides when the
    # stored config has colors, so they are not merged with the palette.
    stored = {"palette": "default", "colors": {"CDS": "#123456"}, "colorsAreOverrides": True}
    assert migrate_session_draft_values(stored, 44) == {"palette": "default", "colors": {"CDS": "#123456"}}
    # Before 40 the flag reaches the split, which merges the colors.
    assert migrate_session_draft_values(stored, 39)["colorsAreOverrides"] is True


@pytest.mark.parametrize("version", [30, 44])
def test_circular_slot_schema_3_is_refused_not_written(version: int) -> None:
    # The Web migrates schema-3 Circular slots, which no Session 27-44 from
    # main or a release holds; Python refuses them rather than write a 46.
    session = {
        "format": "gbdraw-session",
        "version": version,
        "ui": {"mode": "circular"},
        "config": {"adv": {"circular_track_slots_schema_version": 3, "circular_track_slots": ["features"]}},
    }
    with pytest.raises(SessionFormatError, match="Circular slot schema 3") as excinfo:
        migrate_session_flat_draft(session)
    assert excinfo.value.diagnostic == FIELDS_INVALID


CLI_WRITER_CONFIG_VECTORS = json.loads(
    (FIXTURES / "cli-writer-config-vectors.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize(
    "case", CLI_WRITER_CONFIG_VECTORS, ids=[case["name"] for case in CLI_WRITER_CONFIG_VECTORS]
)
def test_cli_writer_config_vectors(case: dict[str, Any]) -> None:
    """OV-269: shared with tests/web/session-active-config-contract.test.mjs."""

    assert _holds_cli_writer_config(case["session"]) is case["holdsCliWriterConfig"]
