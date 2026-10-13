#!/usr/bin/env python3
"""Run `gbdraw circular ...` in-process and dump the label list that prepare_label_list returns.

Usage: PYTHONPATH=<worktree> python3 -P dump_labels.py <out.json> circular --gbk ... (normal CLI args)
Also records which candidate (legacy/current) rearrange_labels_fc picked and the packer angles.
"""
import json
import math
import sys

out_json = sys.argv[1]
argv = sys.argv[2:]

import gbdraw.labels.circular as lc
import gbdraw.diagrams.circular.assemble as asm

captured = {"calls": []}
_orig_prepare = lc.prepare_label_list
_orig_legacy_place = lc._legacy_place_labels_on_arc_fc
_orig_rearrange = lc.rearrange_labels_fc
_orig_score = lc._placement_score


def _legacy_place(labels, cx, cy, xr, yr, sa, ea, total_length):
    result = _orig_legacy_place(labels, cx, cy, xr, yr, sa, ea, total_length)
    angles = []
    for lab in result:
        ang = math.degrees(math.atan2((lab["start_y"] - cy) / yr, (lab["start_x"] - cx) / xr))
        angles.append(round(ang, 3))
    captured.setdefault("legacy_place", []).append(
        {"n": len(result), "x_radius": xr, "y_radius": yr, "angles_after_sequential_pack": angles,
         "texts": [lab["label_text"] for lab in result]}
    )
    return result


def _rearrange(labels, *a, **k):
    result = _orig_rearrange(labels, *a, **k)
    captured.setdefault("rearrange", []).append(
        {"n_in": len(labels), "is_outer": k.get("is_outer"), "score": list(_orig_score(result, a[1])) if result else None}
    )
    return result


def _prepare(*a, **k):
    result = _orig_prepare(*a, **k)
    keep = []
    for lab in result:
        keep.append({key: val for key, val in lab.items() if isinstance(val, (int, float, str, bool)) or val is None})
    captured["calls"].append({"total_length": a[1], "radius": a[2], "labels": keep})
    return result


lc._legacy_place_labels_on_arc_fc = _legacy_place
lc.rearrange_labels_fc = _rearrange
lc.prepare_label_list = _prepare
asm.prepare_label_list = _prepare

from gbdraw.cli import main

sys.argv = ["gbdraw"] + argv
try:
    rc = main()
finally:
    with open(out_json, "w") as fh:
        json.dump(captured, fh)
sys.exit(rc)
