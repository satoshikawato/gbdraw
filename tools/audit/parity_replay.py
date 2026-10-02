"""Replay the Source recipe commands saved by parity-sweep.audit.spec.js.

Each probe directory has a manifest.json, the GUI Result SVG for every probe,
and the helper files Run info offered. This script runs each recipe with the
CLI from this checkout in a fresh work directory and compares the CLI SVG with
the GUI SVG through tests/utils/svg_compare.compare_svgs, ignoring the same
binding attributes as the Gallery publication parity check plus the root
baseProfile attribute.

Usage:
    python tools/audit/parity_replay.py <probe-dir> <input> [<input> ...]
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tests.utils.svg_compare import compare_svgs  # noqa: E402

IGNORED_ATTRIBUTES = {
    # The CLI SVG (svgwrite) has a root baseProfile attribute; the GUI Result does not.
    "baseProfile",
    "data-label-feature-id",
    "data-gbdraw-label-binding-schema",
    "data-record-key",
    "data-record-translation-x",
    "data-record-translation-y",
}


def replay_probe(probe_dir: Path, entry: dict, inputs: list[Path]) -> dict:
    name = entry["name"]
    command = entry.get("command") or ""
    if entry.get("status") != "ok" or not command.startswith("gbdraw"):
        return {
            "name": name,
            "outcome": "NO_RECIPE",
            "gui_status": entry.get("status"),
            "error": entry.get("err"),
            "recipe": entry.get("recipe"),
        }
    work = probe_dir / f"cli_{name}"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    for path in inputs:
        shutil.copy(path, work)
    helper_dir = probe_dir / f"helpers_{name}"
    if helper_dir.is_dir():
        for path in helper_dir.iterdir():
            if path.name != "helpers.zip" and path.is_file():
                shutil.copy(path, work)
    argv = [sys.executable, "-m", "gbdraw.cli", *shlex.split(command)[1:]]
    environment = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    result = subprocess.run(
        argv, cwd=work, env=environment, capture_output=True, text=True, check=False
    )
    svgs = sorted(work.glob("*.svg"))
    if result.returncode != 0 or not svgs:
        tail = result.stderr.strip().splitlines()[-1:] if result.stderr else []
        return {"name": name, "outcome": "CLI_FAILED", "returncode": result.returncode, "stderr": tail}
    comparison = compare_svgs(
        probe_dir / f"{name}.svg",
        svgs[0],
        max_differences=8,
        ignored_attributes=IGNORED_ATTRIBUTES,
    )
    return {
        "name": name,
        "outcome": "MATCH" if comparison.equal else "DIFF",
        "same_as_baseline_in_gui": entry.get("sameAsBaseline"),
        "message": comparison.message,
        "differences": comparison.differences,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("probe_dir", type=Path)
    parser.add_argument("inputs", type=Path, nargs="+")
    args = parser.parse_args()
    probe_dir = args.probe_dir.resolve()
    manifest = json.loads((probe_dir / "manifest.json").read_text(encoding="utf-8"))
    inputs = [path.resolve() for path in args.inputs]
    rows = [replay_probe(probe_dir, entry, inputs) for entry in manifest if "name" in entry]
    for row in rows:
        print(f"{row['name']}: {row['outcome']}")
        for line in row.get("differences", [])[:4]:
            print(f"    {line}")
    (probe_dir / "replay-report.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
