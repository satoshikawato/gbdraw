#!/usr/bin/env python3
"""Fit study (LP-D41): render today's placement with one change at a time and record the metrics.

Usage: fit_sweep.py --out DIR [--jobs 6] [--wt WORKTREE] [--formats svg] CASE:VARIANT [CASE:VARIANT ...]
CASE: a key of page/extra_cases.json or scripts/cases.json.
VARIANT: '+'-joined parts, applied to the case's own arguments:
  cur        no change
  f<pt>      --label_font_size <pt> (replaces the case's value)
  x<k>       --outer_label_x_radius_offset <k>
  y<k>       --outer_label_y_radius_offset <k>
  xy<k>      --outer_label_x_radius_offset <k> and --outer_label_y_radius_offset <k>
  s<px>      --circular_label_spacing <px>
Example: p1_f14:f12.5  p1_f14:y1.1  p1_f14:f13+y1.05
Output: DIR/<case>__<variant>.{svg,metrics.json,labels.json,log}; an existing metrics file is kept (delete it to rerun).
"""
import argparse, json, os, pathlib, subprocess, time
from concurrent.futures import ThreadPoolExecutor

HERE = pathlib.Path(__file__).resolve().parent
TASK = HERE.parent
ROOTS = {'EX': '/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs',
         'IN': str(TASK / 'work' / 'inputs')}
ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True); ap.add_argument('--jobs', type=int, default=6)
ap.add_argument('--wt', default='/home/kawato/gbdraw-work/.worktrees/lp-fit')
ap.add_argument('--formats', default='svg'); ap.add_argument('runs', nargs='+')
ns = ap.parse_args()
ROOTS['TI'] = f'{ns.wt}/tests/test_inputs'
extra = json.loads((TASK / 'page' / 'extra_cases.json').read_text())['cases']
std = json.loads((HERE / 'cases.json').read_text())
out = pathlib.Path(ns.out).resolve(); out.mkdir(parents=True, exist_ok=True)


def case_args(case):
    if case in extra:
        return list(extra[case]['args'])
    c = std['cases'][case]
    root, rel = c['gbk'].split('/', 1)
    return ['circular', '--gbk', f'{ROOTS[root]}/{rel}', *std['common'], '--label_font_size', str(c['font'])]


def set_opt(args, opt, value):
    if opt in args:
        i = args.index(opt); args[i + 1] = value
    else:
        args += [opt, value]


def apply_variant(args, variant):
    for part in variant.split('+'):
        if part == 'cur':
            continue
        if part.startswith('xy'):
            set_opt(args, '--outer_label_x_radius_offset', part[2:]); set_opt(args, '--outer_label_y_radius_offset', part[2:])
        elif part[0] == 'x':
            set_opt(args, '--outer_label_x_radius_offset', part[1:])
        elif part[0] == 'y':
            set_opt(args, '--outer_label_y_radius_offset', part[1:])
        elif part[0] == 'f':
            set_opt(args, '--label_font_size', part[1:])
        elif part[0] == 's':
            set_opt(args, '--circular_label_spacing', part[1:])
        else:
            raise SystemExit(f'unknown variant part {part!r}')
    return args


def one(run):
    case, variant = run.split(':', 1)
    name = f'{case}__{variant}'
    mpath = out / f'{name}.metrics.json'
    if mpath.exists():
        return f'{name}: kept'
    args = apply_variant(case_args(case), variant)
    cmd = ['/home/kawato/micromamba/bin/python3.13', '-P', str(HERE / 'run_case_dev.py'), '--out', str(out), '--name', name,
           '--formats', ns.formats, '--', *args]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env={**os.environ, 'PYTHONPATH': ns.wt, 'TMPDIR': '/home/kawato/gbdraw-baselines/tmp'})
    (out / f'{name}.log').write_text(proc.stdout + '\n--- stderr ---\n' + proc.stderr)
    if proc.returncode != 0 or not mpath.exists():
        return f'{name}: FAILED rc={proc.returncode}: ' + (proc.stderr.strip().splitlines() or ['?'])[-1][:300]
    m = json.loads(mpath.read_text())
    keys = ('labels_external', 'text_overprint_pairs', 'text_overlap_pairs', 'leader_crossing_pairs', 'leader_through_text', 'canvas_px')
    return f'{name}: ' + ' '.join(f'{k}={m.get(k)}' for k in keys) + f' s={time.time() - t0:.0f}'


with ThreadPoolExecutor(ns.jobs) as pool:
    for line in pool.map(one, ns.runs):
        print(line, flush=True)
