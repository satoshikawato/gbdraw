"""The baseline runner must reject misleading comparisons and mixed imports."""
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("protein_benchmark", ROOT / "tools/benchmark_protein_comparison.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def report():
    return {"benchmark": benchmark.NAME, "schema": 1, "measurement": "timing",
            "settings": {"policy": benchmark.POLICY}, "dependencies": {"python": "test"},
            "source": {"head": "test"}, "cases": {"case": {"fixture": {"seed": 1},
            "stages": {"stage": {"semanticSha256": "same", "samples": [100.] * 7,
                                  "median": 100., "noisePct": 0.}}}}}


def test_comparison_rejects_scientific_differences_and_reports_noise():
    base = report()
    changed = copy.deepcopy(base)
    stage = changed["cases"]["case"]["stages"]["stage"]
    stage["semanticSha256"] = "different"
    assert benchmark.compare_reports(base, changed)[1] == 1
    stage["semanticSha256"] = "same"
    stage.update(samples=[120.] * 7, median=120.)
    assert benchmark.compare_reports(base, changed)[1] == 1
    stage["noisePct"] = 8.
    assert benchmark.compare_reports(base, changed)[1] == 2
    stage.update(samples=[100.], median=100., noisePct=0.)
    assert benchmark.compare_reports(base, changed)[1] == 2


def test_comparison_preserves_archived_sampling_policy():
    current = report()
    current["cases"]["case"]["stages"]["stage"]["samples"] = [100.] * 3
    assert benchmark.compare_reports(current, current)[1] == 0
    archived = copy.deepcopy(current)
    archived["settings"]["policy"]["samples"] = 7
    assert benchmark.compare_reports(archived, archived)[1] == 2


@pytest.mark.parametrize("key", ["fixture", "dependencies", "settings"])
def test_comparison_refuses_different_input_or_environment(key):
    base = report()
    changed = copy.deepcopy(base)
    owner = changed["cases"]["case"] if key == "fixture" else changed
    owner[key] = {"changed": True}
    with pytest.raises(ValueError, match="incomparable"):
        benchmark.compare_reports(base, changed)


def test_oracle_preserves_types_order_and_nonfinite_values():
    def encode(x):
        return benchmark.json_bytes(benchmark.canonical(x))
    df = pd.DataFrame({"a": [1., float("nan"), float("inf")]})
    assert encode(df) != encode(df.iloc[::-1])
    assert encode([1, 2]) != encode((1, 2))
    assert encode({"a": 1, "b": 2}) != encode({"b": 2, "a": 1})
    assert encode(pd.DataFrame({"a": [1]})) != encode(pd.DataFrame({"a": [1.]}))
    json.loads(encode(df))


def process_report(sample, sha="same"):
    current = report()
    current["settings"] = {"policy": benchmark.POLICY, "seed": 1, "warmups": 1, "samples": 1}
    current["processMaxRssKiB"] = 1
    current["cases"]["case"]["stages"]["stage"] = {"semanticSha256": sha, "samples": [sample],
                                                   **benchmark.spread([sample])}
    return current


def test_process_reports_merge_into_one_sample_per_process():
    merged = benchmark.merge_process_reports([process_report(x) for x in (100., 104., 98.)])
    stage = merged["cases"]["case"]["stages"]["stage"]
    assert stage["samples"] == [100., 104., 98.]
    assert (stage["median"], stage["mad"], stage["noisePct"]) == (100., 2., 2.)
    assert merged["settings"]["samples"] == merged["settings"]["processes"] == 3
    assert merged["settings"]["policy"]["maxNoisePct"] == 10.
    assert benchmark.compare_reports(merged, merged)[1] == 0
    # Spread between processes up to 10% still decides; above it the stage is inconclusive.
    noisy = benchmark.merge_process_reports([process_report(x) for x in (100., 108., 92.)])
    assert noisy["cases"]["case"]["stages"]["stage"]["noisePct"] == 8.
    assert benchmark.compare_reports(merged, noisy)[1] == 0
    noisier = benchmark.merge_process_reports([process_report(x) for x in (100., 112., 88.)])
    assert benchmark.compare_reports(merged, noisier)[1] == 2
    slower = benchmark.merge_process_reports([process_report(x) for x in (112., 111., 113.)])
    assert benchmark.compare_reports(merged, slower)[1] == 1
    # Samples from one process and samples from separate processes are not compared.
    with pytest.raises(ValueError, match="incomparable reports: settings"):
        benchmark.compare_reports(report(), merged)


@pytest.mark.parametrize("change", ["semantic", "samples"])
def test_merge_refuses_changed_results_and_multi_sample_reports(change):
    reports = [process_report(100.) for _ in range(3)]
    if change == "semantic":
        reports[2]["cases"]["case"]["stages"]["stage"]["semanticSha256"] = "other"
    else:
        for each in reports:
            each["settings"]["samples"] = 3
    with pytest.raises(ValueError, match="semantic results changed|one timing sample"):
        benchmark.merge_process_reports(reports)


def test_compare_sources_alternates_trees_in_fresh_processes(tmp_path, monkeypatch):
    calls = []

    def fake_run(command, check):
        root = command[command.index("--source-root") + 1]
        output = Path(command[command.index("--output") + 1])
        calls.append((root, output.name, command[command.index("--samples") + 1]))
        output.write_text(json.dumps(process_report(100. if root == "base" else 120.)))

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)
    assert benchmark.main(["compare-sources", "--baseline-root", "base", "--current-root", "head",
                           "--cases", "hsp-edges", "--output-dir", str(tmp_path)]) == 1
    assert calls == [(root, f"{side}-{index}.json", "1") for index in (1, 2, 3)
                     for root, side in (("base", "baseline"), ("head", "current"))]
    assert json.loads((tmp_path / "comparison.json").read_text())["exitCode"] == 1
    assert json.loads((tmp_path / "current.json").read_text())["cases"]["case"]["stages"]["stage"]["samples"] == [120.] * 3


def test_compare_sources_runs_real_processes(tmp_path):
    result = subprocess.run([sys.executable, str(ROOT / "tools/benchmark_protein_comparison.py"),
                             "compare-sources", "--baseline-root", str(ROOT), "--current-root", str(ROOT),
                             "--cases", "hsp-edges", "--warmups", "0", "--samples", "1",
                             "--output-dir", str(tmp_path)], capture_output=True, text=True)
    # One sample per tree is fewer than the policy's three, so the comparison is inconclusive.
    assert result.returncode == 2, result.stderr
    assert sorted(path.name for path in (tmp_path / "processes").iterdir()) == ["baseline-1.json", "current-1.json"]
    comparison = json.loads((tmp_path / "comparison.json").read_text())
    assert comparison["comparisons"] and all(row["semanticEqual"] for row in comparison["comparisons"].values())


def test_cli_compares_archived_reports(tmp_path):
    source = tmp_path / "baseline.json.gz"
    output = tmp_path / "comparison.json.gz"
    source.write_bytes(gzip.compress(json.dumps(report()).encode(), mtime=0))
    subprocess.run([sys.executable, str(ROOT / "tools/benchmark_protein_comparison.py"),
                    "compare", "--baseline", str(source), "--current", str(source),
                    "--output", str(output)], check=True)
    assert benchmark.read_report(output)["exitCode"] == 0


def test_foreign_import_is_rejected_before_loading_requested_source(tmp_path):
    code = """
import importlib.util, pathlib, sys, types
spec = importlib.util.spec_from_file_location('bench', sys.argv[1])
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
fake = types.ModuleType('gbdraw.analysis.foreign')
fake.__file__ = sys.argv[3]
sys.modules[fake.__name__] = fake
bench.load_source(pathlib.Path(sys.argv[2]))
"""
    result = subprocess.run([sys.executable, "-c", code, str(ROOT / "tools/benchmark_protein_comparison.py"),
                             str(ROOT), str(tmp_path / "foreign.py")], capture_output=True, text=True)
    assert result.returncode != 0
    assert "foreign gbdraw import" in result.stderr


def test_fresh_process_oracles_repeat_and_measure_existing_work(tmp_path):
    outputs = []
    for i in range(2):
        output = tmp_path / f"run-{i}.json"
        subprocess.run([sys.executable, str(ROOT / "tools/benchmark_protein_comparison.py"),
                        "run", "--source-root", str(ROOT), "--measure", "probe",
                        "--cases", "hsp-edges", "cache-64", "path-8", "merge-edges", "support-edges",
                        "--output", str(output)], check=True)
        outputs.append(json.loads(output.read_text()))
    assert benchmark.compare_reports(*outputs)[1] == 0
    cases = outputs[0]["cases"]
    cache = cases["cache-64"]["stages"]["cache_cycle"]["summary"]
    # A comparison larger than the payload cache's 64 entries keeps every table: a display-only change reuses all
    # of them, a threshold change filters again, and the next conversion replaces the hold.
    assert (cache["coldMisses"], cache["displayOnlyHits"], cache["displayOnlyMisses"]) == (64, 64, 0)
    assert (cache["thresholdMisses"], cache["heldTablesAfterSmaller"]) == (64, 1)
    assert cases["path-8"]["stages"]["selector"]["summary"]["paths"] == 64
    merge = cases["merge-edges"]["stages"]
    assert merge["max_conflicts_0"]["semanticSha256"] != merge["max_conflicts_1"]["semanticSha256"]
    assert cases["support-edges"]["stages"]["selector"]["summary"]["members"] == 6
