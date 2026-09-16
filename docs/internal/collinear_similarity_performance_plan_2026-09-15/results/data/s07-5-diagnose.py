"""S07.5 planning-only fit experiment and unit allocation inventory.

No production patch. Timing, inventory and tracemalloc are separate processes.
The existing benchmark owns Gallery loading, source identity and canonical hashes.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gc
import gzip
import importlib.metadata
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
import sys
import time
import tracemalloc
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT))
from tools import benchmark_protein_comparison as bench  # noqa: E402

pc, cc = bench.load_source(ROOT)
from gbdraw.analysis import collinearity_units as units  # noqa: E402


def narrow_fit_rows(hits, *, top_fraction):
    """Candidate: retain pandas ordering and scalar math, carry four columns."""
    if hits.empty:
        return []
    ordered = hits.loc[:, ["length_product", "query", "subject", "bitscore"]].sort_values(
        ["length_product", "query", "subject"],
        ascending=[True, True, True], kind="mergesort",
    )
    bin_size = max(4, int(math.ceil(len(ordered) / 8)))
    points = []
    for start in range(0, len(ordered), bin_size):
        block = ordered.iloc[start:start + bin_size]
        count = max(1, int(math.ceil(len(block) * float(top_fraction))))
        selected = block.sort_values(
            ["bitscore", "query", "subject"],
            ascending=[False, True, True], kind="mergesort",
        ).head(count)
        for row in selected.itertuples(index=False):
            length = pc._row_float(row, "length_product", 0.0)
            score = pc._row_float(row, "bitscore", 0.0)
            if length <= 0.0 or score <= 0.0:
                continue
            points.append((math.log10(length), math.log10(score)))
    return points


def semantic(value):
    return bench.digest(bench.json_bytes(bench.canonical(value)))


def self_check():
    """Independent small list oracle, including bin boundaries and full ties."""
    import pandas as pd
    rng = random.Random(20260916)
    checks = 0
    for size in [0, 1, 3, 4, 7, 8, 15, 16, 31, 32, 33, 65]:
        for fraction in [0.1, 0.25, 0.5, 1.0]:
            rows = [dict(query=f"q{rng.randrange(3)}", subject=f"s{rng.randrange(3)}",
                         bitscore=rng.choice([0.0, 10.0, 20.0]),
                         length_product=rng.choice([1.0, 100.0, 1000.0]),
                         ignored=i) for i in range(size)]
            frame = pd.DataFrame(rows, columns=["query", "subject", "bitscore", "length_product", "ignored"])
            ordered = sorted(rows, key=lambda r: (r["length_product"], r["query"], r["subject"]))
            width = max(4, math.ceil(size / 8))
            expected = []
            for start in range(0, size, width):
                block = ordered[start:start + width]
                selected = sorted(block, key=lambda r: (-r["bitscore"], r["query"], r["subject"]))
                for row in selected[:max(1, math.ceil(len(block) * fraction))]:
                    if row["bitscore"] > 0 and row["length_product"] > 0:
                        expected.append((math.log10(row["length_product"]), math.log10(row["bitscore"])))
            before = semantic(frame)
            assert pc._select_normalized_fit_rows(frame, top_fraction=fraction) == expected
            assert narrow_fit_rows(frame, top_fraction=fraction) == expected
            assert semantic(frame) == before
            checks += 1
    return checks


def prepare(name):
    gallery = bench.GALLERIES[{"hep": 0, "similarity": 1, "vibrio": 2}[name]]
    _, extraction, records, raw, inventory = bench.gallery_input(ROOT, gallery, pc)
    settings = dict(bitscore=50, evalue=0.01, identity=0, alignment_length=0)
    tables = {pair: pc.filter_protein_hits_by_thresholds(pc.parse_losatp_outfmt6(text), **settings)
              for pair, text in raw.items()}
    if name != "similarity":
        tables = cc._filter_hit_tables_by_search_scope(tables, record_count=len(records), scope="adjacent")
    members = {pair: pc._select_member_candidate_hits_per_query(frame, max_hits=5)
               for pair, frame in tables.items()}
    frames = []
    for frame in members.values():
        aggregate = pc._aggregate_hsps_by_protein_pair(pc._coerce_outfmt6_numeric_columns(frame), extraction.protein_map)
        frames.append(aggregate.loc[aggregate["min_coverage"].astype(float) >= 0.0].reset_index(drop=True))
    inventory.update(effectiveMode="orthogroup" if name == "similarity" else "collinear",
                     inference=True, scope="all" if name == "similarity" else "adjacent",
                     memberMaxHits=5, benchmarkThresholds=settings,
                     semanticTables=len(tables), semanticRows=sum(map(len, tables.values())),
                     memberRows=sum(map(len, members.values())), aggregatedRows=sum(map(len, frames)),
                     fitInputSha256=semantic(frames))
    return frames, extraction, records, inventory


def run_fit(frames):
    return [pc._select_normalized_fit_rows(frame, top_fraction=pc._ORTHOGROUP_SPLIT_TOP_FRACTION)
            for frame in frames]


def run_candidate(frames):
    return [narrow_fit_rows(frame, top_fraction=pc._ORTHOGROUP_SPLIT_TOP_FRACTION) for frame in frames]


def inventory_case(frames, extraction, records):
    counts = Counter()
    original_type = units.CollinearityUnit
    original_alias = units._add_alias
    def construct(*args, **kwargs):
        counts["unitConstructions"] += 1
        return original_type(*args, **kwargs)
    def alias(aliases, value):
        counts["aliasAttempts"] += 1
        counts["aliasListLengthAtAttempt"] += len(aliases)
        return original_alias(aliases, value)
    with patch.object(units, "CollinearityUnit", construct), patch.object(units, "_add_alias", alias):
        index = units.build_collinearity_unit_index(extraction, records=records)
    all_units = [unit for record in index.units_by_record for unit in record]
    counts.update(units=len(all_units), proteins=len(extraction.protein_map),
                  memberSortCalls=len(all_units), memberSortInputReferences=sum(len(u.cds_members) for u in all_units),
                  aliases=sum(len(u.aliases) for u in all_units),
                  aliasSetdefaultAttempts=sum(len(u.aliases) for u in all_units),
                  distinctRecordAliases=sum(len(a) + len(b) for a, b in zip(index.aliases_by_record, index.ambiguous_aliases_by_record)),
                  ambiguousRecordAliases=sum(map(len, index.ambiguous_aliases_by_record)),
                  collapsedUnits=sum(len(u.cds_members) > 1 for u in all_units))
    assert counts["unitConstructions"] == 2 * counts["units"]
    assert semantic(run_fit(frames)) == semantic(run_candidate(frames))
    models = [pc._fit_expected_bitscore_model(f, top_fraction=pc._ORTHOGROUP_SPLIT_TOP_FRACTION) for f in frames]
    with patch.object(pc, "_select_normalized_fit_rows", narrow_fit_rows):
        candidate_models = [pc._fit_expected_bitscore_model(f, top_fraction=pc._ORTHOGROUP_SPLIT_TOP_FRACTION) for f in frames]
    assert semantic(models) == semantic(candidate_models)
    return dict(unitCounts=dict(counts), unitSha256=semantic(index), fitPointsSha256=semantic(run_fit(frames)),
                exactModelParity=True, modelSha256=semantic(models),
                fullFrameDeepBytes=sum(int(f.memory_usage(deep=True).sum()) for f in frames),
                fourColumnDeepBytes=sum(int(f.loc[:, ["length_product", "query", "subject", "bitscore"]].memory_usage(deep=True).sum()) for f in frames),
                columns=sorted(set(len(f.columns) for f in frames)),
                bins=sum(math.ceil(len(f) / max(4, math.ceil(len(f) / 8))) for f in frames if len(f)),
                fitPoints=sum(map(len, run_fit(frames))))


def measure_fit(frames, measure):
    output = {}
    for label, operation in [("original", run_fit), ("narrow", run_candidate)]:
        if measure == "timing":
            operation(frames)
            samples, hashes = [], []
            for _ in range(7):
                gc.collect()
                start = time.perf_counter_ns()
                value = operation(frames)
                samples.append((time.perf_counter_ns() - start) / 1e6)
                hashes.append(semantic(value))
            median = statistics.median(samples)
            mad = statistics.median(abs(s - median) for s in samples)
            output[label] = dict(samplesMs=samples, medianMs=median, madMs=mad,
                                 noisePct=100 * mad / median, sampleSemanticSha256=hashes)
        else:
            gc.collect()
            tracemalloc.start()
            value = operation(frames)
            retained, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            output[label] = dict(retainedBytes=retained, peakBytes=peak, semanticSha256=semantic(value))
    if measure == "timing":
        assert len(set(output["original"]["sampleSemanticSha256"] + output["narrow"]["sampleSemanticSha256"])) == 1
        ratio = output["narrow"]["medianMs"] / output["original"]["medianMs"]
        output["changePct"] = 100 * (ratio - 1)
        output["medianRegression"] = ratio > 1.1
        output["decision"] = "inconclusive" if max(output[k]["noisePct"] for k in ["original", "narrow"]) > 5 else "regression" if ratio > 1.1 else "pass"
    else:
        assert output["original"]["semanticSha256"] == output["narrow"]["semanticSha256"]
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--measure", choices=["inventory", "timing", "memory"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = ["gbdraw/analysis/protein_colinearity.py", "gbdraw/analysis/collinearity.py",
             "gbdraw/analysis/collinearity_units.py", "tools/benchmark_protein_comparison.py"]
    report = dict(purpose="S07.5 planning-only; isolated fit-row selection, not full normalization/analysis/Web",
                  measure=args.measure, command=sys.argv, selfCheckCases=self_check(),
                  head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  sourceSha256={p: bench.digest((ROOT / p).read_bytes()) for p in paths},
                  scriptSha256=bench.digest(Path(__file__).read_bytes()),
                  dependencies={p: importlib.metadata.version(p) for p in ["pandas", "numpy", "biopython"]},
                  python=sys.version, boundary="prepared per-table aggregated/min-coverage-qualified frames -> ordered fit points; includes candidate projection",
                  policy=dict(warmups=1, samples=7, regressionPct=10, maxNoisePct=5), cases={})
    for name in ["hep", "similarity", "vibrio"]:
        frames, extraction, records, fixture = prepare(name)
        result = inventory_case(frames, extraction, records) if args.measure == "inventory" else measure_fit(frames, args.measure)
        report["cases"][name] = dict(fixture=fixture, result=result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({k: v["result"] for k, v in report["cases"].items()}, indent=2))


if __name__ == "__main__":
    main()
