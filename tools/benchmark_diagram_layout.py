#!/usr/bin/env python3
"""Benchmark diagram assembly time and memory for the layout performance gate.

The standard profile uses three warmups, 21 timing samples, and three separate
``tracemalloc`` samples per case.  ``--quick`` uses no warmup and one sample of
each kind with smaller fixtures; it exists for contract tests, not release
decisions.  Timing noise is the median absolute deviation divided by the median.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import statistics
import sys
import time
import tracemalloc
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any


BENCHMARK_NAME = "gbdraw-diagram-layout-assembly"
SCHEMA_VERSION = 1
DEFAULT_THRESHOLD_PCT = 5.0
CASE_NAMES = ("singleCircular", "multiCircular", "multiLinear")


@dataclass(frozen=True)
class BenchmarkSettings:
    profile: str
    warmups: int
    timing_samples: int
    memory_samples: int


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    assemble: Callable[[], object]


class BenchmarkDataError(ValueError):
    """Raised when benchmark arguments or comparison data are invalid."""


def _record_specs(profile: str) -> dict[str, list[tuple[int, int]]]:
    if profile == "quick":
        return {
            "singleCircular": [(1_200, 4)],
            "multiCircular": [(1_200, 4), (1_500, 5)],
            "multiLinear": [(1_300, 4), (1_600, 5)],
        }
    return {
        "singleCircular": [(24_000, 48)],
        "multiCircular": [
            (12_000, 28),
            (15_000, 34),
            (9_000, 22),
            (18_000, 40),
        ],
        "multiLinear": [
            (14_000, 32),
            (17_000, 38),
            (11_000, 26),
            (20_000, 44),
        ],
    }


def _load_api(source_root: Path) -> tuple[Any, Any, Any, Any]:
    root = source_root.resolve()
    if not (root / "gbdraw" / "api" / "diagram.py").is_file():
        raise BenchmarkDataError(f"source root does not contain gbdraw: {root}")
    loaded_package = sys.modules.get("gbdraw")
    if loaded_package is not None:
        loaded_file = getattr(loaded_package, "__file__", None)
        if loaded_file is None or not Path(loaded_file).resolve().is_relative_to(root):
            raise BenchmarkDataError(
                "gbdraw is already imported from a different source root; "
                "benchmark each tree in a fresh process"
            )
    sys.path.insert(0, str(root))
    # These three functions are the narrowest common assembly boundary in HEAD
    # and the working tree.  Package-root draw_* also measures request compilation
    # and Diagram wrapping, which are outside this layout assembly gate.
    from gbdraw.api.config import apply_config_overrides
    from gbdraw.api.diagram import (
        assemble_circular_diagram_from_record,
        assemble_circular_diagram_from_records,
        assemble_linear_diagram_from_records,
    )

    return (
        apply_config_overrides,
        assemble_circular_diagram_from_record,
        assemble_circular_diagram_from_records,
        assemble_linear_diagram_from_records,
    )


def _make_record(record_id: str, length: int, feature_count: int) -> Any:
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    motif = "ATGCGTAC"
    sequence = (motif * ((length + len(motif) - 1) // len(motif)))[:length]
    record = SeqRecord(
        Seq(sequence),
        id=record_id,
        name=record_id,
        description=f"Deterministic layout benchmark record {record_id}",
    )
    record.annotations["molecule_type"] = "DNA"
    feature_types = ("CDS", "tRNA", "CDS", "rRNA", "repeat_region")
    stride = max(1, length // feature_count)
    feature_length = max(30, min(360, (stride * 2) // 3))
    features = []
    for index in range(feature_count):
        start = min(length - 2, index * stride + max(1, stride // 8))
        end = min(length - 1, start + feature_length)
        if end <= start:
            continue
        feature_type = feature_types[index % len(feature_types)]
        features.append(
            SeqFeature(
                FeatureLocation(start, end, strand=1 if index % 2 == 0 else -1),
                type=feature_type,
                qualifiers={
                    "gene": [f"gene_{index + 1:03d}"],
                    "product": [f"benchmark {feature_type} product {index + 1:03d}"],
                },
            )
        )
    record.features = features
    return record


def _build_cases(source_root: Path, profile: str) -> list[BenchmarkCase]:
    from pandas import DataFrame

    (
        apply_config_overrides,
        assemble_single_circular,
        assemble_multi_circular,
        assemble_multi_linear,
    ) = _load_api(source_root)
    specs = _record_specs(profile)
    records = {
        case_name: [
            _make_record(f"{case_name}_{index + 1}", length, feature_count)
            for index, (length, feature_count) in enumerate(case_specs)
        ]
        for case_name, case_specs in specs.items()
    }
    config = apply_config_overrides(
        None,
        {
            "canvas.show_gc": False,
            "canvas.show_skew": False,
            "labels.circular.scope": "outer",
            "labels.linear.scope": "all",
        },
    )
    default_colors = DataFrame(
        [
            ("CDS", "#54bcf8"),
            ("rRNA", "#71ee7d"),
            ("tRNA", "#e8b441"),
            ("repeat_region", "#d3d3d3"),
            ("default", "#d3d3d3"),
            ("skew_high", "#6dded3"),
            ("skew_low", "#ad72e3"),
            ("gc_content", "#a1a1a1"),
            ("pairwise_match", "#d3d3d3"),
            ("pairwise_match_min", "#FFE7E7"),
            ("pairwise_match_max", "#FF7272"),
            ("collinear_block_plus_min", "#f0f1f5"),
            ("collinear_block_plus", "#8b9cc1"),
            ("collinear_block_minus_min", "#FFE7E7"),
            ("collinear_block_minus", "#E15759"),
        ],
        columns=["feature_type", "color"],
    )

    return [
        BenchmarkCase(
            "singleCircular",
            lambda: assemble_single_circular(
                records["singleCircular"][0],
                cfg=config,
                default_colors=default_colors,
                legend="right",
                plot_title="Single Circular layout benchmark",
                plot_title_position="top",
            ),
        ),
        BenchmarkCase(
            "multiCircular",
            lambda: assemble_multi_circular(
                records["multiCircular"],
                cfg=config,
                default_colors=default_colors,
                legend="right",
                plot_title="Multi-record Circular layout benchmark",
                plot_title_position="top",
            ),
        ),
        BenchmarkCase(
            "multiLinear",
            lambda: assemble_multi_linear(
                records["multiLinear"],
                cfg=config,
                default_colors=default_colors,
                legend="right",
                plot_title="Multi-record Linear layout benchmark",
                plot_title_position="top",
            ),
        ),
    ]


def _time_once(assemble: Callable[[], object]) -> float:
    gc.collect()
    started = time.perf_counter_ns()
    drawing = assemble()
    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000.0
    if drawing is None:
        raise BenchmarkDataError("assembly returned None")
    del drawing
    return elapsed_ms


def _peak_memory_once(assemble: Callable[[], object]) -> int:
    gc.collect()
    tracemalloc.start()
    try:
        drawing = assemble()
        if drawing is None:
            raise BenchmarkDataError("assembly returned None")
        _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        del drawing
        return peak_bytes
    finally:
        tracemalloc.stop()


def _relative_mad_pct(samples: Sequence[float]) -> float:
    median = statistics.median(samples)
    if median <= 0.0:
        return 0.0
    mad = statistics.median(abs(sample - median) for sample in samples)
    return 100.0 * mad / median


def _summarize_case(
    timing_samples_ms: Sequence[float],
    peak_memory_samples_bytes: Sequence[int],
) -> dict[str, Any]:
    if not timing_samples_ms or not peak_memory_samples_bytes:
        raise BenchmarkDataError("each case needs timing and memory samples")
    return {
        "medianAssemblyTimeMs": round(statistics.median(timing_samples_ms), 6),
        "timeNoisePct": round(_relative_mad_pct(timing_samples_ms), 6),
        "peakMemoryBytes": max(peak_memory_samples_bytes),
        "timingSamplesMs": [round(value, 6) for value in timing_samples_ms],
        "peakMemorySamplesBytes": list(peak_memory_samples_bytes),
    }


def run_benchmarks(
    *,
    source_root: Path,
    settings: BenchmarkSettings,
    label: str,
) -> dict[str, Any]:
    cases = _build_cases(source_root, settings.profile)
    results: dict[str, Any] = {}
    for case in cases:
        for _ in range(settings.warmups):
            _time_once(case.assemble)
        timing_samples = [
            _time_once(case.assemble) for _ in range(settings.timing_samples)
        ]
        memory_samples = [
            _peak_memory_once(case.assemble) for _ in range(settings.memory_samples)
        ]
        results[case.name] = _summarize_case(timing_samples, memory_samples)

    return {
        "benchmark": BENCHMARK_NAME,
        "schemaVersion": SCHEMA_VERSION,
        "label": label,
        "settings": {
            "profile": settings.profile,
            "warmupsPerCase": settings.warmups,
            "timingSamplesPerCase": settings.timing_samples,
            "memorySamplesPerCase": settings.memory_samples,
            "timingStatistic": "median",
            "noiseStatistic": "median absolute deviation / median",
            "memoryStatistic": "maximum tracemalloc peak",
        },
        "fixture": _record_specs(settings.profile),
        "cases": results,
    }


def _pct_change(current: float, baseline: float) -> float:
    if baseline <= 0.0:
        raise BenchmarkDataError("baseline measurements must be greater than zero")
    return 100.0 * (current / baseline - 1.0)


def _validate_comparable(
    baseline: Mapping[str, Any], current: Mapping[str, Any]
) -> None:
    for report_name, report in (("baseline", baseline), ("current", current)):
        if report.get("benchmark") != BENCHMARK_NAME:
            raise BenchmarkDataError(f"{report_name} benchmark name is invalid")
        if report.get("schemaVersion") != SCHEMA_VERSION:
            raise BenchmarkDataError(f"{report_name} schema version is invalid")
        cases = report.get("cases")
        if not isinstance(cases, Mapping) or set(cases) != set(CASE_NAMES):
            raise BenchmarkDataError(f"{report_name} cases are incomplete")
        for case_name in CASE_NAMES:
            case = cases[case_name]
            if not isinstance(case, Mapping):
                raise BenchmarkDataError(
                    f"{report_name} case {case_name} must be an object"
                )
            _finite_metric(
                case.get("medianAssemblyTimeMs"),
                name=f"{report_name}.{case_name}.medianAssemblyTimeMs",
                positive=True,
            )
            _finite_metric(
                case.get("timeNoisePct"),
                name=f"{report_name}.{case_name}.timeNoisePct",
                positive=False,
            )
            _positive_int_metric(
                case.get("peakMemoryBytes"),
                name=f"{report_name}.{case_name}.peakMemoryBytes",
            )
    if baseline.get("fixture") != current.get("fixture"):
        raise BenchmarkDataError("baseline and current fixtures differ")
    if baseline.get("settings") != current.get("settings"):
        raise BenchmarkDataError("baseline and current settings differ")


def _finite_metric(value: Any, *, name: str, positive: bool) -> float:
    if isinstance(value, bool):
        raise BenchmarkDataError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise BenchmarkDataError(f"{name} must be a finite number") from exc
    if not math.isfinite(number) or (number <= 0.0 if positive else number < 0.0):
        qualifier = "positive and " if positive else "non-negative and "
        raise BenchmarkDataError(f"{name} must be {qualifier}finite")
    return number


def _positive_int_metric(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise BenchmarkDataError(f"{name} must be a positive integer")
    return value


def compare_reports(
    baseline: Mapping[str, Any],
    current: Mapping[str, Any],
    *,
    threshold_pct: float = DEFAULT_THRESHOLD_PCT,
) -> tuple[dict[str, Any], bool]:
    threshold_pct = _finite_metric(
        threshold_pct,
        name="thresholdPct",
        positive=False,
    )
    _validate_comparable(baseline, current)
    compared_cases: dict[str, Any] = {}
    has_regression = False
    has_noise_exemption = False
    for case_name in CASE_NAMES:
        baseline_case = baseline["cases"][case_name]
        current_case = current["cases"][case_name]
        baseline_time = _finite_metric(
            baseline_case["medianAssemblyTimeMs"],
            name=f"baseline.{case_name}.medianAssemblyTimeMs",
            positive=True,
        )
        current_time = _finite_metric(
            current_case["medianAssemblyTimeMs"],
            name=f"current.{case_name}.medianAssemblyTimeMs",
            positive=True,
        )
        delta_pct = _pct_change(current_time, baseline_time)
        measured_noise_pct = max(
            float(baseline_case["timeNoisePct"]),
            float(current_case["timeNoisePct"]),
        )
        effective_threshold_pct = max(threshold_pct, measured_noise_pct)
        threshold_time = baseline_time * (1.0 + threshold_pct / 100.0)
        if current_time <= threshold_time:
            decision = "pass"
        elif measured_noise_pct > delta_pct:
            decision = "noise_exempt"
            has_noise_exemption = True
        else:
            decision = "regression"
            has_regression = True
        baseline_peak = _positive_int_metric(
            baseline_case["peakMemoryBytes"],
            name=f"baseline.{case_name}.peakMemoryBytes",
        )
        current_peak = _positive_int_metric(
            current_case["peakMemoryBytes"],
            name=f"current.{case_name}.peakMemoryBytes",
        )
        compared_cases[case_name] = {
            "baselineMedianAssemblyTimeMs": round(baseline_time, 6),
            "currentMedianAssemblyTimeMs": round(current_time, 6),
            "timeDeltaPct": round(delta_pct, 6),
            "measuredNoisePct": round(measured_noise_pct, 6),
            "effectiveThresholdPct": round(effective_threshold_pct, 6),
            "timeDecision": decision,
            "baselinePeakMemoryBytes": int(baseline_peak),
            "currentPeakMemoryBytes": int(current_peak),
            "peakMemoryDeltaPct": round(_pct_change(current_peak, baseline_peak), 6),
        }

    outcome = "fail" if has_regression else "pass"
    if has_noise_exemption and not has_regression:
        outcome = "pass_with_noise_exemption"
    comparison = {
        "benchmark": BENCHMARK_NAME,
        "schemaVersion": SCHEMA_VERSION,
        "baselineLabel": str(baseline.get("label", "unspecified")),
        "currentLabel": str(current.get("label", "unspecified")),
        "thresholdPct": round(threshold_pct, 6),
        "outcome": outcome,
        "cases": compared_cases,
    }
    return comparison, has_regression


def _write_json(payload: Mapping[str, Any], output: Path | None) -> None:
    rendered = json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n"
    if output is None:
        sys.stdout.write(rendered)
    else:
        output.write_text(rendered, encoding="utf-8")


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BenchmarkDataError(f"cannot read benchmark JSON {path}: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise BenchmarkDataError(f"benchmark JSON must be an object: {path}")
    return payload


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def _non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="run all three assembly cases")
    run_parser.add_argument("--source-root", type=Path, default=Path.cwd())
    run_parser.add_argument("--label", default="unspecified")
    run_parser.add_argument("--output", type=Path)
    run_parser.add_argument("--quick", action="store_true")
    run_parser.add_argument("--warmups", type=_non_negative_int)
    run_parser.add_argument("--samples", type=_positive_int)
    run_parser.add_argument("--memory-samples", type=_positive_int)

    compare_parser = subparsers.add_parser(
        "compare", help="compare two JSON reports and enforce the time gate"
    )
    compare_parser.add_argument("--baseline", type=Path, required=True)
    compare_parser.add_argument("--current", type=Path, required=True)
    compare_parser.add_argument("--threshold-pct", type=float, default=DEFAULT_THRESHOLD_PCT)
    compare_parser.add_argument("--output", type=Path)
    return parser


def _resolve_settings(args: argparse.Namespace) -> BenchmarkSettings:
    profile = "quick" if args.quick else "standard"
    return BenchmarkSettings(
        profile=profile,
        warmups=args.warmups if args.warmups is not None else (0 if args.quick else 3),
        timing_samples=args.samples if args.samples is not None else (1 if args.quick else 21),
        memory_samples=(
            args.memory_samples
            if args.memory_samples is not None
            else (1 if args.quick else 3)
        ),
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "run":
            report = run_benchmarks(
                source_root=args.source_root,
                settings=_resolve_settings(args),
                label=args.label,
            )
            _write_json(report, args.output)
            return 0
        comparison, has_regression = compare_reports(
            _read_json(args.baseline),
            _read_json(args.current),
            threshold_pct=args.threshold_pct,
        )
        _write_json(comparison, args.output)
        return 1 if has_regression else 0
    except BenchmarkDataError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
