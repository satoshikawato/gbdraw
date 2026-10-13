#!/usr/bin/env python3
"""Run several (case, design) renders in parallel processes.

Usage: run_matrix.py --out DIR [--jobs 4] [--opts k=v,...] CASE:DESIGN [CASE:DESIGN ...]
  DESIGN is 'cur' (today's placement), 'R' (today's --label_placement radial) or a prototype id (A1 A2 A3 B C D E).
  Output: DIR/<case>__<design>.{svg,png,metrics.json,labels.json,log}; prints one summary line per run.
"""
import argparse, json, pathlib, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = pathlib.Path(__file__).resolve().parent
WT = __import__('os').environ.get('LP_WT', '/home/kawato/gbdraw-work/.worktrees/label-study')   # LP_WT: another prototype worktree
ROOTS = {"TI": f"{WT}/tests/test_inputs", "EX": "/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs",
         "IN": str(HERE.parent / "work" / "inputs")}
ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--opts", default=""); ap.add_argument("--suffix", default=""); ap.add_argument("runs", nargs="+")
ns = ap.parse_args()
spec = json.loads((HERE / "cases.json").read_text())
out = pathlib.Path(ns.out).resolve(); out.mkdir(parents=True, exist_ok=True)


def one(run):
    case, design = run.split(":")
    c = spec["cases"][case]
    root, rel = c["gbk"].split("/", 1)
    args = ["circular", "--gbk", f"{ROOTS[root]}/{rel}", *spec["common"], "--label_font_size", str(c["font"])]
    name = f"{case}__{design}{ns.suffix}"
    cmd = ["/home/kawato/micromamba/bin/python3.13", "-P", str(HERE / "run_case.py"), "--out", str(out), "--name", name]
    if design == "R":
        args += ["--label_placement", "radial"]
    elif design != "cur":
        cmd += ["--design", design, "--opts", ns.opts]
    proc = subprocess.run(cmd + ["--", *args], capture_output=True, text=True, env={**__import__("os").environ, "PYTHONPATH": WT,
                          "TMPDIR": "/home/kawato/gbdraw-baselines/tmp"})
    (out / f"{name}.log").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    mpath = out / f"{name}.metrics.json"
    if proc.returncode != 0 or not mpath.exists():
        return f"{name}: FAILED rc={proc.returncode} (see {out / (name + '.log')}): " + (proc.stderr.strip().splitlines() or ["?"])[-1]
    m = json.loads(mpath.read_text())
    keys = ("labels_external", "font_size", "text_overprint_pairs", "text_overlap_pairs", "leader_crossing_pairs", "leader_through_text",
            "leaders_over_features", "leader_len_mean", "leader_len_max", "canvas_px", "seconds")
    return f"{name}: " + " ".join(f"{k}={m.get(k)}" for k in keys) + f" proto={m.get('proto')}"


with ThreadPoolExecutor(ns.jobs) as pool:
    for line in pool.map(one, ns.runs):
        print(line, flush=True)
