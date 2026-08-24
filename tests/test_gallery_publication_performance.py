from __future__ import annotations

import json
from pathlib import Path

import pytest

from tools.measure_gallery_publication_performance import enforce_performance_gate


BASELINE = json.loads(
    (
        Path(__file__).parent
        / "performance"
        / "gallery-publication-baseline.json"
    ).read_text(encoding="utf-8")
)
WORKFLOW = (
    Path(__file__).parents[1] / ".github" / "workflows" / "gallery-publication.yml"
).read_text(encoding="utf-8")


def test_fixed_gallery_performance_baseline_enforces_125_percent_gate() -> None:
    for operation, values in BASELINE["operations"].items():
        trials = [
            {"wallSeconds": wall, "peakRssKiB": rss}
            for wall, rss in zip(
                values["wallSeconds"],
                values["peakRssKiB"],
                strict=True,
            )
        ]
        assert enforce_performance_gate(BASELINE, operation, trials)["passed"]

    failing_trials = [
        {
            "wallSeconds": BASELINE["operations"]["projection"][
                "wallMedianSeconds"
            ]
            * 1.26,
            "peakRssKiB": BASELINE["operations"]["projection"][
                "peakRssMaxKiB"
            ],
        }
    ] * 3
    with pytest.raises(RuntimeError, match="wall-time regression"):
        enforce_performance_gate(BASELINE, "projection", failing_trials)


def test_runtime_baseline_compares_trials_from_the_same_environment() -> None:
    baseline_trials = [
        {"wallSeconds": wall, "peakRssKiB": rss}
        for wall, rss in zip(
            (17.2, 17.4, 17.3),
            (3_100_000, 3_120_000, 3_110_000),
            strict=True,
        )
    ]
    passing_trials = [
        {"wallSeconds": wall, "peakRssKiB": rss}
        for wall, rss in zip(
            (20.0, 20.5, 20.2),
            (3_700_000, 3_750_000, 3_720_000),
            strict=True,
        )
    ]
    summary = enforce_performance_gate(
        BASELINE,
        "projection",
        passing_trials,
        baseline_trials=baseline_trials,
    )
    assert summary["passed"]
    assert summary["wallLimitSeconds"] == 21.62
    assert summary["peakRssLimitKiB"] == 3_900_000

    failing_trials = [
        {"wallSeconds": 21.8, "peakRssKiB": 3_750_000},
    ] * 3
    with pytest.raises(RuntimeError, match="wall-time regression"):
        enforce_performance_gate(
            BASELINE,
            "projection",
            failing_trials,
            baseline_trials=baseline_trials,
        )


def test_complete_refresh_performance_gate_is_manual_only() -> None:
    assert "complete_refresh:" in WORKFLOW
    assert (
        "if: ${{ github.event_name == 'workflow_dispatch' "
        "&& inputs.complete_refresh }}"
    ) in WORKFLOW
