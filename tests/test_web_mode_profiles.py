from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

from gbdraw.legend.table import configure_pairwise_identity_legend_from_comparisons


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


def test_web_mode_profile_consumers_use_mode_specific_defaults() -> None:
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    contract_source = (
        WEB_ROOT / "js" / "services" / "session-active-config-contract.js"
    ).read_text(encoding="utf-8")
    reset_source = (WEB_ROOT / "js" / "services" / "reset.js").read_text(
        encoding="utf-8"
    )
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(
        encoding="utf-8"
    )
    candidate_source = (WEB_ROOT / "js" / "app" / "candidate-render.js").read_text(
        encoding="utf-8"
    )
    request_source = (WEB_ROOT / "js" / "services" / "session-request.js").read_text(
        encoding="utf-8"
    )
    sanitization_source = (
        WEB_ROOT / "js" / "services" / "svg-sanitization.js"
    ).read_text(encoding="utf-8")

    assert "createDefaultAdv = (mode = 'circular')" in contract_source
    assert "...comparisonStateForMode(mode)" in contract_source
    assert "features: [...MODE_DEFAULT_FEATURE_TYPES]" in contract_source
    assert "trackDefaultsForMode('circular')" in contract_source
    assert "trackDefaultsForMode('linear')" in contract_source
    assert "managedAdvStateForMode(mode).axis_stroke_color" in contract_source
    assert "'data-gbdraw-role'" in sanitization_source
    assert "'data-gbdraw-orientation'" in sanitization_source
    assert "createDefaultAdv(state.mode.value)" in reset_source
    assert "modeProfileStateManager?.reset?" in reset_source
    # The mode transition swaps the profiles; tests/web/session-active-mode
    # ("each mode keeps its own title and fonts") covers it in the browser.

    # Each mode resolves its thresholds with its own defaults (X-02: one
    # resolver on the generated domains; Generate keeps the draft).
    mode_profiles_source = (WEB_ROOT / "js" / "mode-profiles.js").read_text(encoding="utf-8")
    resolver = mode_profiles_source.split("export const resolveComparisonThresholds", 1)[1]
    assert "comparisonFiltersForMode(mode)" in resolver.split("export const", 1)[0]
    assert "resolveComparisonThresholds(drawing.adv, 'circular')" in run_source
    assert "resolveComparisonThresholds(drawing.adv, 'linear')" in run_source
    assert "normalizeBlastThreshold" not in run_source
    assert "resolveComparisonThresholds(drawing.adv, state.mode.value)" in request_source
    assert "comparisonFiltersForMode('linear')" in run_source
    assert not (WEB_ROOT / "js" / "app" / "cli-args.js").exists()
    assert "effectiveLinearAxisColor({" in request_source
    assert "drawing.modeProfileStateManager?.isManaged?." in request_source
    blast_config = SimpleNamespace()
    color_modes = configure_pairwise_identity_legend_from_comparisons(
        blast_config,
        None,
        additional_color_modes=("orientation",),
    )
    assert color_modes == {"orientation"}
    assert blast_config.hide_pairwise_identity_legend is True
    assert "shouldSuppressPairwiseIdentityLegend" not in run_source
    assert "suppressPairwiseIdentityLegend" not in candidate_source
    assert "PAIRWISE_LEGEND_SELECTOR" not in candidate_source
