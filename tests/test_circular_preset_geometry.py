from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ORACLE_COMMIT = "0228e6fe896768c6dc58513945d935e96578f5d0"
ORACLE_PATH = REPO_ROOT / "tests" / "fixtures" / "circular_preset_oracle" / f"{ORACLE_COMMIT}.json"


def _load_capture_module():
    module_path = REPO_ROOT / "tools" / "capture_circular_preset_oracle.py"
    spec = importlib.util.spec_from_file_location("capture_circular_preset_oracle", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load oracle capture module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _case_key(case: dict[str, Any]) -> tuple[str, str, bool, str]:
    return (
        str(case["input"]),
        str(case["preset"]),
        bool(case["strandedness"]),
        str(case["visibility"]),
    )


def _assert_band_matches(actual: dict[str, Any] | None, expected: dict[str, Any] | None) -> None:
    if expected is None:
        assert actual is None
        return
    assert actual is not None
    for key in ("inner_px", "outer_px", "width_px", "center_px"):
        assert actual[key] == pytest.approx(expected[key], abs=1e-6)


def _assert_slot_matches(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    for key in ("slot_index", "id", "renderer", "side", "z", "explicit_anchor", "explicit_width", "compressed"):
        assert actual[key] == expected[key]
    for key in ("anchor_radius_px", "anchor_offset_px", "requested_width_px", "resolved_width_px"):
        if expected[key] is None:
            assert actual[key] is None
        else:
            assert actual[key] == pytest.approx(expected[key], abs=1e-6)
    _assert_band_matches(actual["packing_band_px"], expected["packing_band_px"])
    _assert_band_matches(actual["draw_band_px"], expected["draw_band_px"])
    _assert_band_matches(actual["reserved_band_px"], expected["reserved_band_px"])
    assert actual["params"] == expected["params"]


def test_circular_preset_geometry_matches_fixed_oracle() -> None:
    capture_module = _load_capture_module()
    expected = json.loads(ORACLE_PATH.read_text(encoding="utf-8"))
    actual = capture_module.capture(
        [
            Path("tests/test_inputs/MG1655.gbk"),
            Path("tests/test_inputs/AP027280.gb"),
        ]
    )

    assert expected["base_commit"] == ORACLE_COMMIT
    assert actual["schema"] == expected["schema"] == 1

    expected_cases = {_case_key(case): case for case in expected["cases"]}
    actual_cases = {_case_key(case): case for case in actual["cases"]}
    assert actual_cases.keys() == expected_cases.keys()

    for key, expected_case in expected_cases.items():
        actual_case = actual_cases[key]
        assert actual_case["axis_radius_px"] == pytest.approx(expected_case["axis_radius_px"], abs=1e-6)
        assert actual_case["outer_content_radius_px"] == pytest.approx(
            expected_case["outer_content_radius_px"],
            abs=1e-6,
        )
        assert len(actual_case["slots"]) == len(expected_case["slots"])
        for actual_slot, expected_slot in zip(actual_case["slots"], expected_case["slots"]):
            _assert_slot_matches(actual_slot, expected_slot)
