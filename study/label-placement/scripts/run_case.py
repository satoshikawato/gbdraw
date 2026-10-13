#!/usr/bin/env python3
"""Render one circular case with today's placement or a study prototype, and write metrics.

Usage:
  run_case.py --out DIR --name NAME [--design A1|A2|A3|B|C|D|E] [--opts k=v,k=v] [--formats svg,png] -- <gbdraw circular args>
Run with: PYTHONPATH=/home/kawato/gbdraw-work/.worktrees/label-study /home/kawato/micromamba/bin/python3.13 -P run_case.py ...
Writes DIR/NAME.svg (+ .png), DIR/NAME.labels.json (the drawn label list), DIR/NAME.metrics.json.
Do not pass -o / -f in the gbdraw args.
"""
import argparse
import json
import os
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
os.environ.setdefault("TMPDIR", "/home/kawato/gbdraw-baselines/tmp")

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--name", required=True)
ap.add_argument("--design", default="")
ap.add_argument("--opts", default="")
ap.add_argument("--formats", default="svg,png")
ap.add_argument("rest", nargs=argparse.REMAINDER)
ns = ap.parse_args()
rest = ns.rest[1:] if ns.rest and ns.rest[0] == "--" else ns.rest
out = pathlib.Path(ns.out).resolve()
out.mkdir(parents=True, exist_ok=True)
os.environ["GBDRAW_LABEL_PROTO"] = ns.design
os.environ["GBDRAW_LABEL_PROTO_OPTS"] = ns.opts

import gbdraw.labels.circular as lc
import gbdraw.diagrams.circular.assemble as asm
import gbdraw.labels.proto as proto
import lp_metrics

calls = []
_orig = asm.prepare_label_list


def _capture(*a, **k):
    result = _orig(*a, **k)
    calls.append((int(a[1]), result))
    return result


asm.prepare_label_list = _capture
from gbdraw.cli import main

os.chdir(out)
sys.argv = ["gbdraw", *rest, "-o", ns.name, "-f", ns.formats, "--overwrite"]
t0 = time.time()
rc = main()
elapsed = time.time() - t0
if calls:
    total_length, labels = calls[-1]
    metrics = lp_metrics.compute(labels, total_length, lc, svg_path=str(out / f"{ns.name}.svg"))
    metrics.update(design=ns.design or "current", opts=ns.opts, seconds=round(elapsed, 1),
                   proto=dict(proto.LAST_RESULT), argv=rest)
    (out / f"{ns.name}.metrics.json").write_text(json.dumps(metrics, indent=1))
    keep = [{k: v for k, v in lab.items() if isinstance(v, (int, float, str, bool, list, tuple)) or v is None} for lab in labels]
    (out / f"{ns.name}.labels.json").write_text(json.dumps({"total_length": total_length, "labels": keep}))
    print(json.dumps({k: v for k, v in metrics.items() if k not in ("argv",)}))
sys.exit(rc or 0)
