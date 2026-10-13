"""Hard criteria of the label placement study over rendered metrics files (orchestrator's check).

Usage: check_hard.py [--allow FIELD]... [--zero FIELD]... [--end] GLOB [GLOB ...]   (globs of *.metrics.json, quoted)
Prints one line per file that fails and a total. LP-D13: labels_back_over_height replaces labels_toward_axis.
--zero adds a field that must be 0; --end adds the three LP-D15 fields (where a leader ends).
"""
import argparse, glob, json, os

ZERO = ['leader_crossing_pairs', 'text_overlap_pairs', 'leader_through_text', 'labels_back_over_height',
        'labels_on_objects', 'labels_near_objects', 'leaders_on_objects', 'leader_acute_bends', 'leader_bends_lt110']
END = ['leaders_end_far_half', 'leaders_end_off_middle', 'leader_ends_unclear', 'leader_end_gap_ink_over5']
ap = argparse.ArgumentParser(); ap.add_argument('--allow', action='append', default=[]); ap.add_argument('--zero', action='append', default=[])
ap.add_argument('--end', action='store_true'); ap.add_argument('globs', nargs='+')
ns = ap.parse_args()
ZERO = ZERO + ns.zero + (END if ns.end else [])
files = sorted({f for g in ns.globs for f in glob.glob(g)})
bad = 0
for f in files:
    m = json.load(open(f)); name = os.path.basename(f)[:-len('.metrics.json')]
    why = [f'{k}={m.get(k)}' for k in ZERO if k not in ns.allow and m.get(k) != 0]
    dropped = (m.get('proto') or {}).get('dropped')
    if dropped:
        why.append(f'dropped={len(dropped) if hasattr(dropped, "__len__") else dropped}')
    gap = m.get('label_gap_min_px')
    if gap is not None and gap < 3.0:
        why.append(f'label_gap_min_px={gap}')
    if why:
        bad += 1; print(f'FAIL {name}: ' + ' '.join(why))
print(f'{len(files)} files, {bad} failing' + (f' (allowed: {" ".join(ns.allow)})' if ns.allow else ''))
