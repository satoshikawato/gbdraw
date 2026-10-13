#!/usr/bin/env python3
"""Measurement only: run scope_run.py with the 300-label cap of Linear _resolve_above_feature_label_overlaps removed
(the rest of the function unchanged), and time the resolver. usage: same as scope_run.py"""
import runpy, sys, time, pathlib
import gbdraw.labels.linear as ll
src_timer = {'s': 0.0, 'n': 0}
_orig_place = ll._place_linear_label_without_overlap


def resolve_nocap(labels, min_gap_px):
    t0 = time.perf_counter()
    if not labels:
        return
    for place_above in (True, False):
        side_labels = [label for label in labels if bool(label.get("above_feature_place_above", True)) is place_above]
        if not side_labels:
            continue
        placed = []; placed_by_id = {}
        placed_index = ll.AabbIndex(bucket_size=ll._above_feature_label_bucket_size(side_labels, min_gap_px))
        direction = -1.0 if place_above else 1.0
        for label in sorted(side_labels, key=lambda item: float(item.get("feature_anchor_x", item["middle_x"]))):
            if ll._place_linear_label_without_overlap(label, placed, min_gap_px, direction, placed_index=placed_index, placed_by_id=placed_by_id):
                label["leader_line"] = True
                ll._update_linear_label_leader_end(label)
            pid = len(placed_by_id); placed_by_id[pid] = label
            placed_index.insert(pid, ll._linear_label_aabb(label)); placed.append(label)
    src_timer['s'] += time.perf_counter() - t0; src_timer['n'] += 1
    print('resolver_nocap', len(labels), 'labels', round(time.perf_counter() - t0, 2), 's', file=sys.stderr)


ll._resolve_above_feature_label_overlaps = resolve_nocap
sys.argv = [str(pathlib.Path(__file__).with_name('scope_run.py'))] + sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
