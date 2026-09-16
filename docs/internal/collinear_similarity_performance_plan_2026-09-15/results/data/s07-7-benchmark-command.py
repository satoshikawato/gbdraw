"""Run only the authorized S07.6-final/current post-search comparison pairs.

The existing benchmark owns preparation, measurement, serialization and policy.
This recorder adds fresh processes, source guards and host/command provenance.
Run with taskset -c 3 env PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 python ...
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
BASE = ROOT.parent / "collinear-s07-6-20260916"
RUNNER = ROOT / "tools/benchmark_protein_comparison.py"
CASES = ("gallery-collinear", "gallery-collinear-off", "gallery-collinear-vibrio",
         "gallery-collinear-vibrio-off", "gallery-orthogroup", "gallery-orthogroup-unbounded")
LOGS = OUT / "s07-7-benchmark-logs"
LOGS.mkdir(exist_ok=True)


def host():
    return {"ps": subprocess.check_output(
        ["ps", "-eo", "pid,ppid,pcpu,etimes,comm", "--sort=-pcpu"], text=True),
        "loadavg": Path("/proc/loadavg").read_text()}


def verify_source(root, label):
    hashes = json.loads((OUT / "s07-7-start.json").read_text())["productionSha256"]
    if label == "current":
        hashes["gbdraw/analysis/protein_colinearity.py"] = (
            "608abd7f9e0ff266e0746c7a961712c74d358ae85c27653355aa8ccce4fbc17d")
    for path, expected in hashes.items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--measure", choices=("timing", "memory", "probe"))
    parser.add_argument("--sources", nargs="+", choices=("baseline", "current"), default=("baseline", "current"))
    parser.add_argument("--tag", default="")
    parser.add_argument("--samples", type=int)
    args = parser.parse_args()
    assert sorted(os.sched_getaffinity(0)) == [3]
    assert os.environ.get("PYTHONHASHSEED") == "0"
    assert os.environ.get("PYTHONDONTWRITEBYTECODE") == "1"
    for measure in ((args.measure,) if args.measure else ("timing", "memory")):
        for case in args.cases:
            for label, root in (("baseline", BASE), ("current", ROOT)):
                if label not in args.sources:
                    continue
                verify_source(root, label)
                name = f"{label}-{case}-{measure}" + (f"-{args.tag}" if args.tag else "")
                output = OUT / f"s07-7-{name}.json.gz"
                # Existing observations must never be silently overwritten.
                if output.exists():
                    raise FileExistsError(output)
                samples = args.samples if args.samples is not None else 3
                command = [sys.executable, "-B", str(RUNNER), "run", "--source-root", str(root),
                           "--cases", case, "--stages", "post_search", "--measure", measure,
                           "--warmups", "1", "--samples", str(samples), "--output", str(output)]
                entry = {"name": name, "argv": command, "cwd": str(ROOT),
                         "startUtc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                         "affinity": sorted(os.sched_getaffinity(0)), "PYTHONHASHSEED": "0",
                         "PYTHONDONTWRITEBYTECODE": "1", "hostBefore": host(),
                         "runnerSha256": hashlib.sha256(RUNNER.read_bytes()).hexdigest(),
                         "log": str((LOGS / f"{name}.log").relative_to(OUT))}
                print(f"START {name}", flush=True)
                tick = time.monotonic()
                with (LOGS / f"{name}.log").open("w") as log:
                    result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                entry.update(endUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                             wallSeconds=time.monotonic() - tick, exitCode=result.returncode,
                             hostAfter=host())
                with (OUT / "s07-7-benchmark-commands.jsonl").open("a") as stream:
                    stream.write(json.dumps(entry) + "\n")
                print(f"END {name} exit={result.returncode}", flush=True)
                if result.returncode:
                    return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
