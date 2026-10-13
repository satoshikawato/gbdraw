#!/usr/bin/env python3
"""Run scope_run.py on cases of scope_cases.json (4 at a time); an existing NAME.scope.json is kept.

usage: scope_sweep.py [--out work/SCOPE] [--jobs 4] [--timeout 1500] [CASE ...]   (no CASE: every case)
CASE may carry a variant after ':' as '+'-joined parts: f<pt> (--label_font_size), y<k>/x<k> (outer offsets),
iy<k>/ix<k> (inner offsets), s<px> (--circular_label_spacing), ls<px> (--linear_label_spacing), fh<px>
(--feature_height), w<px> (--canvas_width, Linear), r<deg> (--label_rotation).
"""
import argparse, json, os, pathlib, subprocess, time
from concurrent.futures import ThreadPoolExecutor

HERE = pathlib.Path(__file__).resolve().parent
T = HERE.parent
WT = '/home/kawato/gbdraw-work/.worktrees/lp-scope'
ap = argparse.ArgumentParser()
ap.add_argument('--out', default=str(T / 'work' / 'SCOPE')); ap.add_argument('--jobs', type=int, default=4)
ap.add_argument('--timeout', type=int, default=1500); ap.add_argument('runs', nargs='*')
ns = ap.parse_args()
cases = json.loads((HERE / 'scope_cases.json').read_text())['cases']
out = pathlib.Path(ns.out); out.mkdir(parents=True, exist_ok=True)
OPT = {'f': '--label_font_size', 'y': '--outer_label_y_radius_offset', 'x': '--outer_label_x_radius_offset',
       'iy': '--inner_label_y_radius_offset', 'ix': '--inner_label_x_radius_offset', 's': '--circular_label_spacing',
       'ls': '--linear_label_spacing', 'fh': '--feature_height', 'w': '--canvas_width', 'r': '--label_rotation'}


def variant(args, var):
    a = list(args)
    for part in var.split('+'):
        if part == 'cur':
            continue
        if part in ('out', 'radial', 'mid'):            # --labels out | --label_placement radial | --track_type middle
            opt, val = {'out': ('--labels', 'out'), 'radial': ('--label_placement', 'radial'), 'mid': ('--track_type', 'middle')}[part]
            if opt in a:
                a[a.index(opt) + 1] = val
            else:
                a += [opt, val]
            continue
        if part.startswith('ixy'):                      # both inner offsets
            for opt in ('--inner_label_x_radius_offset', '--inner_label_y_radius_offset'):
                if opt in a:
                    a[a.index(opt) + 1] = part[3:]
                else:
                    a += [opt, part[3:]]
            continue
        key = next(k for k in sorted(OPT, key=len, reverse=True) if part.startswith(k))
        opt, val = OPT[key], part[len(key):]
        if opt in a:
            a[a.index(opt) + 1] = val
        else:
            a += [opt, val]
    return a


def one(run):
    case, _, var = run.partition(':')
    name = case + (f'__{var}' if var else '')
    if (out / f'{name}.scope.json').exists():
        return f'{name}: kept'
    args = variant(cases[case], var) if var else cases[case]
    cmd = ['/home/kawato/micromamba/bin/python3.13', '-P', str(HERE / 'scope_run.py'), '--out', str(out), '--name', name, '--', *args]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=ns.timeout,
                           env={**os.environ, 'PYTHONPATH': WT, 'TMPDIR': '/home/kawato/gbdraw-baselines/tmp'})
        rc, so, se = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        rc, so, se = 'timeout', str(e.stdout or ''), str(e.stderr or '')
    (out / f'{name}.log').write_text(f'{" ".join(cmd)}\n{so}\n--- stderr ---\n{se}')
    if not (out / f'{name}.scope.json').exists():
        tail = [l for l in (se or '').strip().splitlines() if l.strip()][-1:] or ['?']
        (out / f'{name}.failed.json').write_text(json.dumps({'rc': str(rc), 'error': tail[0][:400], 'seconds': round(time.time() - t0, 1), 'argv': args}))
        return f'{name}: FAILED rc={rc} {tail[0][:200]}'
    m = json.loads((out / f'{name}.scope.json').read_text())
    return f'{name}: kinds={m["labels_by_kind"]} overlaps={m["overlap_pairs"]} crossings={m["crossing_pairs"]} s={time.time() - t0:.0f}'


runs = ns.runs or list(cases)
print(time.strftime('%F %T'), 'start', len(runs), 'runs', flush=True)
with ThreadPoolExecutor(ns.jobs) as pool:
    for line in pool.map(one, runs):
        print(time.strftime('%T'), line, flush=True)
print(time.strftime('%F %T'), 'done', flush=True)
