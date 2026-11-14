#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, TypedDict

REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = REPO_ROOT / "tests" / "performance" / "gallery-publication-baseline.json"
TIME_COMMAND = Path("/usr/bin/time")
GALLERY_SESSIONS = Path("gbdraw/web/gallery/sessions")
PROJECTION_SESSION_ID = "vibrio-harveyi-group-collinear"
PEAK_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


class Trial(TypedDict):
    wallSeconds: float
    peakRssKiB: int


class MeasuredTrial(Trial):
    command: list[str]


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _projection_session(checkout: Path) -> str:
    """Return the Vibrio Gallery Session file as the checked-out revision stores it.

    The fixed base may store the Session in an older file format (a
    .gbdraw-session.json.gz file before the .gbdraw Session file), so each
    revision projects its own file of the same Session.
    """
    matches = sorted((checkout / GALLERY_SESSIONS).glob(f"{PROJECTION_SESSION_ID}.gbdraw*"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one {PROJECTION_SESSION_ID} Gallery Session file in "
            f"{GALLERY_SESSIONS} of {checkout}, found {len(matches)}."
        )
    return matches[0].relative_to(checkout).as_posix()


def _operation_command(operation: str, checkout: Path) -> list[str]:
    if operation == "refresh":
        return [sys.executable, "tools/refresh_gallery_sessions.py", "--no-assets"]
    if operation == "projection":
        return [
            shutil.which("node") or "node",
            # Avoid timing-dependent overlap of inflated Buffers and JSON.parse
            # allocations. Apply the same GC profile to the fixed base and head.
            "--no-incremental-marking",
            "tests/web/session-request.test.mjs",
            "--project-session",
            _projection_session(checkout),
        ]
    raise ValueError(f"Unsupported Gallery performance operation: {operation}")


def _run_trial(operation: str, revision: str, index: int) -> MeasuredTrial:
    with tempfile.TemporaryDirectory(prefix=f"gbdraw-{operation}-perf-") as tmpdir:
        trial_root = Path(tmpdir)
        worktree = trial_root / "worktree"
        subprocess.run(
            ["git", "worktree", "add", "--detach", str(worktree), revision],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        try:
            command = _operation_command(operation, worktree)
            started = time.perf_counter()
            result = subprocess.run(
                [str(TIME_COMMAND), "-v", *command],
                cwd=worktree,
                text=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=False,
            )
            wall_seconds = time.perf_counter() - started
            if result.returncode != 0:
                raise RuntimeError(
                    f"Gallery {operation} performance trial {index} failed:\n{result.stderr}"
                )
            match = PEAK_RSS_RE.search(result.stderr)
            if match is None:
                raise RuntimeError("Could not read peak RSS from /usr/bin/time output.")
            trial: MeasuredTrial = {
                "command": command,
                "wallSeconds": round(wall_seconds, 2),
                "peakRssKiB": int(match.group(1)),
            }
            print(json.dumps({"revision": revision, "trial": index, **trial}), flush=True)
            return trial
        finally:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(worktree)],
                cwd=REPO_ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
            )


def enforce_performance_gate(
    baseline: Mapping[str, Any],
    operation: str,
    trials: list[Trial],
    *,
    baseline_trials: list[Trial],
) -> dict[str, Any]:
    ratio = float(baseline["maximumRegressionRatio"])
    expected_wall_values = [float(trial["wallSeconds"]) for trial in baseline_trials]
    expected_rss_values = [int(trial["peakRssKiB"]) for trial in baseline_trials]
    wall_values = [float(trial["wallSeconds"]) for trial in trials]
    rss_values = [int(trial["peakRssKiB"]) for trial in trials]
    summary = {
        "wallSeconds": wall_values,
        "peakRssKiB": rss_values,
        "wallMedianSeconds": round(statistics.median(wall_values), 2),
        "peakRssMaxKiB": max(rss_values),
    }
    wall_limit = statistics.median(expected_wall_values) * ratio
    rss_limit = int(max(expected_rss_values) * ratio)
    if summary["wallMedianSeconds"] > wall_limit:
        raise RuntimeError(
            f"Gallery {operation} wall-time regression: "
            f"{summary['wallMedianSeconds']} > {wall_limit:.2f} seconds"
        )
    if summary["peakRssMaxKiB"] > rss_limit:
        raise RuntimeError(
            f"Gallery {operation} peak-RSS regression: "
            f"{summary['peakRssMaxKiB']} > {rss_limit} KiB"
        )
    return {
        **summary,
        "wallLimitSeconds": round(wall_limit, 2),
        "peakRssLimitKiB": rss_limit,
        "passed": True,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Measure the fixed Gallery performance gate.")
    parser.add_argument("operation", choices=("refresh", "projection"))
    parser.add_argument("--revision", default="HEAD")
    parser.add_argument("--baseline-revision")
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.trials != 3:
        parser.error("The Gallery performance gate requires exactly three trials.")
    baseline = _load_json(BASELINE_PATH)
    if args.baseline_revision and args.baseline_revision != baseline["baseSha"]:
        parser.error("--baseline-revision must match the checked-in fixed base SHA.")
    # Measure the fixed base in the same environment, executable and flags as
    # the head; no recorded measurement stands in for it.
    baseline_trials = [
        _run_trial(args.operation, baseline["baseSha"], index + 1)
        for index in range(args.trials)
    ]
    trials = [
        _run_trial(args.operation, args.revision, index + 1)
        for index in range(args.trials)
    ]
    report = {
        "schema": 1,
        "baseSha": baseline["baseSha"],
        "headRevision": args.revision,
        "operation": args.operation,
        "pythonVersion": platform.python_version(),
        "platform": platform.platform(),
        "baselineTrials": baseline_trials,
        "trials": trials,
    }
    if args.operation == "projection":
        report["nodeVersion"] = subprocess.check_output(
            [shutil.which("node") or "node", "--version"], text=True
        ).strip()
    try:
        report["summary"] = enforce_performance_gate(
            baseline,
            args.operation,
            trials,
            baseline_trials=baseline_trials,
        )
    except RuntimeError as error:
        report["summary"] = {"passed": False, "error": str(error)}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["summary"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
