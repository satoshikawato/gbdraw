from __future__ import annotations

import json
from pathlib import Path

import pytest

from tools.measure_gallery_publication_performance import enforce_performance_gate
from tools import measure_gallery_publication_performance as measurement


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


@pytest.mark.parametrize("metric", ["peakRssKiB", "wallSeconds"])
def test_projection_failure_retains_all_trials_and_measures_fixed_base(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, metric: str
) -> None:
    baseline_trials = [
        {"wallSeconds": 1.2, "peakRssKiB": 1_082_000},
        {"wallSeconds": 1.3, "peakRssKiB": 1_083_000},
        {"wallSeconds": 1.2, "peakRssKiB": 1_081_000},
    ]
    head_trials = [dict(trial) for trial in baseline_trials]
    if metric == "peakRssKiB":
        # A single high trial must fail even when the median is unchanged.
        head_trials[1][metric] = 1_400_000
    else:
        for trial in head_trials:
            trial[metric] = 1.6
    calls = []

    def run_trial(operation: str, revision: str, index: int) -> dict:
        calls.append((operation, revision, index))
        trials = baseline_trials if revision == BASELINE["baseSha"] else head_trials
        return trials[index - 1]

    monkeypatch.setattr(measurement, "_run_trial", run_trial)
    monkeypatch.setattr(measurement.subprocess, "check_output", lambda *a, **kw: "v20.20.2\n")
    output = tmp_path / "failed.json"
    # Omission of --baseline-revision must not compare the new profile against
    # historical measurements collected with a different Node/GC configuration.
    assert measurement.main([
        "projection", "--revision", "candidate", "--output", str(output)
    ]) == 1
    assert calls == [
        ("projection", revision, index)
        for revision in (BASELINE["baseSha"], "candidate")
        for index in (1, 2, 3)
    ]
    report = json.loads(output.read_text())
    assert report["baselineTrials"] == baseline_trials
    assert report["trials"] == head_trials
    assert report["summary"]["passed"] is False
    assert "regression" in report["summary"]["error"]
    assert report["nodeVersion"] == "v20.20.2"


@pytest.mark.parametrize("trials", [0, 1, 2, 4])
def test_gallery_gate_rejects_other_trial_counts(trials: int) -> None:
    with pytest.raises(SystemExit) as error:
        measurement.main(["projection", "--trials", str(trials)])
    assert error.value.code == 2
