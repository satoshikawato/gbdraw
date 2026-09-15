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
    assert cases["cache-64"]["stages"]["cache_cycle"]["summary"]["warmParses"] == 64
    assert cases["path-8"]["stages"]["selector"]["summary"]["paths"] == 64
    merge = cases["merge-edges"]["stages"]
    assert merge["max_conflicts_0"]["semanticSha256"] != merge["max_conflicts_1"]["semanticSha256"]
    assert cases["support-edges"]["stages"]["selector"]["summary"]["members"] == 6
