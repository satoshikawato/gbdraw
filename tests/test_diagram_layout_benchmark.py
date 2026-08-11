from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from tools import benchmark_diagram_layout as benchmark


REPO_ROOT = Path(__file__).resolve().parents[1]


def _report(
    *,
    label: str,
    times_ms: dict[str, float] | None = None,
    noise_pct: dict[str, float] | None = None,
) -> dict:
    times_ms = times_ms or {name: 100.0 for name in benchmark.CASE_NAMES}
    noise_pct = noise_pct or {name: 2.0 for name in benchmark.CASE_NAMES}
    return {
        "benchmark": benchmark.BENCHMARK_NAME,
        "schemaVersion": benchmark.SCHEMA_VERSION,
        "label": label,
        "settings": {
            "profile": "standard",
            "warmupsPerCase": 3,
            "timingSamplesPerCase": 21,
            "memorySamplesPerCase": 3,
            "timingStatistic": "median",
            "noiseStatistic": "median absolute deviation / median",
            "memoryStatistic": "maximum tracemalloc peak",
        },
        "fixture": benchmark._record_specs("standard"),
        "cases": {
            name: {
                "medianAssemblyTimeMs": times_ms[name],
                "timeNoisePct": noise_pct[name],
                "peakMemoryBytes": 1_000,
            }
            for name in benchmark.CASE_NAMES
        },
    }


def test_case_summary_records_median_noise_and_peak_memory() -> None:
    summary = benchmark._summarize_case(
        [90.0, 95.0, 100.0, 105.0, 110.0],
        [1_000, 1_500, 1_250],
    )

    assert summary["medianAssemblyTimeMs"] == 100.0
    assert summary["timeNoisePct"] == 5.0
    assert summary["peakMemoryBytes"] == 1_500
    assert summary["timingSamplesMs"] == [90.0, 95.0, 100.0, 105.0, 110.0]
    assert summary["peakMemorySamplesBytes"] == [1_000, 1_500, 1_250]


def test_comparison_fails_only_for_regression_above_threshold_and_noise() -> None:
    baseline = _report(label="head")
    current = _report(
        label="working-tree",
        times_ms={
            "singleCircular": 104.999,
            "multiCircular": 108.0,
            "multiLinear": 108.0,
        },
        noise_pct={
            "singleCircular": 1.0,
            "multiCircular": 9.0,
            "multiLinear": 3.0,
        },
    )

    comparison, has_regression = benchmark.compare_reports(baseline, current)

    assert has_regression is True
    assert comparison["outcome"] == "fail"
    assert comparison["cases"]["singleCircular"]["timeDecision"] == "pass"
    noisy = comparison["cases"]["multiCircular"]
    assert noisy["timeDecision"] == "noise_exempt"
    assert noisy["measuredNoisePct"] == 9.0
    assert noisy["effectiveThresholdPct"] == 9.0
    assert comparison["cases"]["multiLinear"]["timeDecision"] == "regression"


def test_comparison_reports_pass_with_noise_exemption() -> None:
    baseline = _report(label="head")
    current = _report(
        label="working-tree",
        times_ms={name: 108.0 for name in benchmark.CASE_NAMES},
        noise_pct={name: 10.0 for name in benchmark.CASE_NAMES},
    )

    comparison, has_regression = benchmark.compare_reports(baseline, current)

    assert has_regression is False
    assert comparison["outcome"] == "pass_with_noise_exemption"
    assert {
        case["timeDecision"] for case in comparison["cases"].values()
    } == {"noise_exempt"}


def test_noise_equal_to_regression_does_not_exempt_regression() -> None:
    baseline = _report(label="head")
    current = _report(
        label="working-tree",
        times_ms={name: 108.0 for name in benchmark.CASE_NAMES},
        noise_pct={name: 8.0 for name in benchmark.CASE_NAMES},
    )

    comparison, has_regression = benchmark.compare_reports(baseline, current)

    assert has_regression is True
    assert comparison["outcome"] == "fail"
    assert {
        case["timeDecision"] for case in comparison["cases"].values()
    } == {"regression"}


def test_exact_threshold_is_not_a_regression() -> None:
    baseline = _report(label="head")
    current = _report(
        label="working-tree",
        times_ms={name: 105.0 for name in benchmark.CASE_NAMES},
    )

    comparison, has_regression = benchmark.compare_reports(baseline, current)

    assert has_regression is False
    assert comparison["outcome"] == "pass"
    assert {
        case["timeDecision"] for case in comparison["cases"].values()
    } == {"pass"}


def test_comparison_rejects_different_fixture_profiles() -> None:
    baseline = _report(label="head")
    current = _report(label="working-tree")
    current["fixture"] = benchmark._record_specs("quick")

    with pytest.raises(benchmark.BenchmarkDataError, match="fixtures differ"):
        benchmark.compare_reports(baseline, current)


@pytest.mark.parametrize("invalid", [None, "fast", float("nan"), float("inf")])
def test_comparison_rejects_invalid_metrics(invalid: object) -> None:
    baseline = _report(label="head")
    current = _report(label="working-tree")
    current["cases"]["multiLinear"]["medianAssemblyTimeMs"] = invalid

    with pytest.raises(benchmark.BenchmarkDataError, match="finite"):
        benchmark.compare_reports(baseline, current)


def test_json_writer_rejects_nonfinite_numbers(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Out of range float values"):
        benchmark._write_json({"invalid": float("nan")}, tmp_path / "result.json")


def test_source_root_rejects_package_loaded_from_another_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target_root = tmp_path / "target"
    api_dir = target_root / "gbdraw" / "api"
    api_dir.mkdir(parents=True)
    (api_dir / "diagram.py").write_text("# fixture\n", encoding="utf-8")
    monkeypatch.setitem(
        sys.modules,
        "gbdraw",
        SimpleNamespace(__file__=str(tmp_path / "other" / "gbdraw" / "__init__.py")),
    )

    with pytest.raises(benchmark.BenchmarkDataError, match="fresh process"):
        benchmark._load_api(target_root)


def test_json_writer_is_sorted_and_newline_terminated(tmp_path: Path) -> None:
    output = tmp_path / "result.json"

    benchmark._write_json({"z": 1, "a": 2}, output)

    assert output.read_text(encoding="utf-8") == '{\n  "a": 2,\n  "z": 1\n}\n'


def test_quick_profile_runs_all_three_real_assemblers(tmp_path: Path) -> None:
    output = tmp_path / "quick.json"

    completed = subprocess.run(
        [
            sys.executable,
            "tools/benchmark_diagram_layout.py",
            "run",
            "--source-root",
            str(REPO_ROOT),
            "--quick",
            "--label",
            "pytest",
            "--output",
            str(output),
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert completed.returncode == 0, completed.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert set(report["cases"]) == set(benchmark.CASE_NAMES)
    assert report["settings"] == {
        "profile": "quick",
        "warmupsPerCase": 0,
        "timingSamplesPerCase": 1,
        "memorySamplesPerCase": 1,
        "timingStatistic": "median",
        "noiseStatistic": "median absolute deviation / median",
        "memoryStatistic": "maximum tracemalloc peak",
    }
    for case in report["cases"].values():
        assert case["medianAssemblyTimeMs"] > 0.0
        assert case["peakMemoryBytes"] > 0
