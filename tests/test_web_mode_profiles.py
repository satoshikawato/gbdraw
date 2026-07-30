from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"


def test_generated_web_mode_profiles_match_python_source() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "tools/generate_mode_profiles.py",
            "--check",
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_web_mode_profile_helpers() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    result = subprocess.run(
        [node, "tests/web/mode-profiles.test.mjs"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_web_mode_profile_consumers_use_mode_specific_defaults() -> None:
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    reset_source = (WEB_ROOT / "js" / "services" / "reset.js").read_text(
        encoding="utf-8"
    )
    setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(
        encoding="utf-8"
    )
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(
        encoding="utf-8"
    )
    request_source = (WEB_ROOT / "js" / "services" / "session-request.js").read_text(
        encoding="utf-8"
    )

    assert "createDefaultAdv = (profileMode = 'circular')" in state_source
    assert "...comparisonStateForMode(profileMode)" in state_source
    assert "features: [...MODE_DEFAULT_FEATURE_TYPES]" in state_source
    assert "trackDefaultsForMode('circular')" in state_source
    assert "trackDefaultsForMode('linear')" in state_source
    assert "managedAdvStateForMode(profileMode).axis_stroke_color" in state_source
    assert "'data-gbdraw-role'" in state_source
    assert "'data-gbdraw-orientation'" in state_source
    assert "createDefaultAdv(state.mode.value)" in reset_source
    assert "modeProfileStateManager?.reset?" in reset_source
    assert "modeProfileStateManager.invalidate(nextMode)" in setup_source
    assert (
        "modeProfileStateManager.transition(adv, previousMode, nextMode)"
        in setup_source
    )

    circular_normalization = run_source.split(
        "const normalizeConservationState = () => {", 1
    )[1].split("const runCircularLosatConservation", 1)[0]
    assert "DEFAULT_LINEAR_BLAST_FILTERS" not in circular_normalization
    assert "DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS" in circular_normalization
    assert "comparisonFiltersForMode('circular')" in run_source
    assert "comparisonFiltersForMode('linear')" in run_source
    assert not (WEB_ROOT / "js" / "app" / "cli-args.js").exists()
    assert "effectiveLinearAxisColor({" in request_source
    assert "state.modeProfileStateManager?.isManaged?." in request_source
    assert "import { PAIRWISE_LEGEND_SELECTOR }" in run_source
    assert "querySelectorAll(PAIRWISE_LEGEND_SELECTOR)" in run_source
    assert (
        """querySelectorAll(
      '[data-gbdraw-role="comparison-legend"]"""
        not in run_source
    )
