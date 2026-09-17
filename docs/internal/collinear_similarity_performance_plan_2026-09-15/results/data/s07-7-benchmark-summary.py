"""Apply the unchanged benchmark policy to the six authorized source pairs."""
from collections import Counter
import hashlib
import gzip
import importlib.util
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
spec = importlib.util.spec_from_file_location("benchmark", ROOT / "tools/benchmark_protein_comparison.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
cases = ("gallery-collinear", "gallery-collinear-off", "gallery-collinear-vibrio",
         "gallery-collinear-vibrio-off", "gallery-orthogroup", "gallery-orthogroup-unbounded")
measured_runner_hash = hashlib.sha256(gzip.decompress((OUT / "s07-7-measured-runner.py.gz").read_bytes())).hexdigest()
rows = {}
reports = {}
for case in cases:
    row = {}
    case_reports = {}
    for measure in ("timing", "memory"):
        pair = []
        for label in ("baseline", "current"):
            path = OUT / f"s07-7-{label}-{case}-{measure}.json.gz"
            report = bench.read_report(path)
            assert report["source"]["runnerSha256"] == measured_runner_hash
            reports[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            pair.append(report)
        comparison, code = bench.compare_reports(*pair)
        b, c = [r["cases"][case]["stages"]["post_search"] for r in pair]
        case_reports[measure] = pair
        row[measure] = comparison["comparisons"][f"{case}/post_search"] | {
            "comparisonExitCode": code, "baselineMedian": b["median"],
            "currentMedian": c["median"], "baselineNoisePct": b["noisePct"],
            "currentNoisePct": c["noisePct"], "sampleCount": len(b["samples"]),
            "semanticSha256": b["semanticSha256"], "semanticBytes": b["semanticBytes"],
            "summary": b["summary"]}
        if measure == "timing":
            row[measure]["medianRegression"] = row[measure]["timeDeltaPct"] > bench.POLICY["regressionPct"]
        else:
            row[measure].update(peakDeltaPct=100 * (c["median"] / b["median"] - 1),
                                baselineRetainedBytes=b["retainedBytes"], currentRetainedBytes=c["retainedBytes"])
    for timing, memory in zip(case_reports["timing"], case_reports["memory"]):
        _, code = bench.compare_reports(timing, memory, semantic_only=True)
        assert code == 0, case
    rows[case] = row
summary = {"baseline": "S07.6 final dirty source; not HEAD alone",
           "candidate": "S07.7 plus local member/RBH/rank reductions; final source",
           "boundary": "prepared filtered raw hit tables -> complete ordered native post-search result",
           "excluded": ["raw search", "parse/filter", "fixture preparation", "result hashing", "render", "browser Generate"],
           "policy": case_reports["timing"][0]["settings"]["policy"],
           "measuredRunnerSha256": measured_runner_hash,
           "futureSampling": "Three samples by default; no automatic noise repeats; archived decisions use their recorded policy",
           "affinity": [3], "PYTHONHASHSEED": "0",
           "timingDecisions": dict(Counter(r["timing"]["decision"] for r in rows.values())),
           "allSemanticEqual": all(r[m]["semanticEqual"] for r in rows.values() for m in ("timing", "memory")),
           "cases": rows, "reportSha256": reports,
           "commandLedger": "s07-7-benchmark-commands.jsonl",
           "limitations": ["one independent tracemalloc observation per source/case; prepared inputs excluded, output retained",
                           "host snapshots before/after each command; no continuous contention or thermal monitoring",
                           "combined A+B+C effects; no per-change timing attribution",
                           "does not assess S07.6 unit speedup relative to S07 or Web end-to-end latency"]}
quick_paths = [OUT / f"s07-7-{label}-gallery-collinear-vibrio-timing-quick3.json.gz"
               for label in ("baseline", "current")]
if all(path.exists() for path in quick_paths):
    quick = [bench.read_report(path) for path in quick_paths]
    comparison, code = bench.compare_reports(*quick)
    stages = [report["cases"]["gallery-collinear-vibrio"]["stages"]["post_search"] for report in quick]
    assert all(len(stage["samples"]) == 3 for stage in stages)
    assert all(stage["semanticSha256"] == rows["gallery-collinear-vibrio"]["timing"]["semanticSha256"]
               for stage in stages)
    summary["quickComparison"] = {"case": "gallery-collinear-vibrio", "samplesPerSource": 3,
                                  "reason": "User requested three samples and stopping the 21-sample repeat",
                                  "comparison": comparison, "comparisonExitCode": code,
                                  "baseline": stages[0], "current": stages[1],
                                  "status": "Descriptive only; fewer than seven samples cannot pass the unchanged gate"}
    for path in quick_paths:
        reports[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary["limitations"].append("The initial noisy Vibrio ON pair remains inconclusive; the three-sample check does not replace it")
(OUT / "s07-7-benchmark-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({k: summary[k] for k in ("timingDecisions", "allSemanticEqual")}, indent=2))
for case, row in rows.items():
    t, m = row["timing"], row["memory"]
    print(f"{case}: {t['baselineMedian']:.3f} -> {t['currentMedian']:.3f} ms "
          f"({t['timeDeltaPct']:+.2f}%, {t['decision']}); "
          f"MAD/median {t['baselineNoisePct']:.2f}/{t['currentNoisePct']:.2f}%; "
          f"peak {m['baselineMedian']} -> {m['currentMedian']} bytes ({m['peakDeltaPct']:+.2f}%)")
