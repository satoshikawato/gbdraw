#!/usr/bin/env python3
"""Render (case, design) pairs for the Owner page: extra cases (page/extra_cases.json) or standard ones (scripts/cases.json).

Usage: run_extra.py --out DIR [--jobs 4] [--opts k=v,...] [--suffix _x] CASE:DESIGN [CASE:DESIGN ...]
DESIGN: cur | R | A1 | A2 | A3 | B | C | D | E. Output as scripts/run_matrix.py: DIR/<case>__<design><suffix>.*
"""
import argparse, json, os, pathlib, subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = pathlib.Path(__file__).resolve().parent
S = HERE.parent / 'scripts'
WT = __import__('os').environ.get('LP_WT', '/home/kawato/gbdraw-work/.worktrees/label-study')   # LP_WT: another prototype worktree
ROOTS = {'TI': f'{WT}/tests/test_inputs', 'EX': '/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs',
         'IN': str(HERE.parent / 'work' / 'inputs')}
ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--jobs', type=int, default=4)
ap.add_argument('--opts', default=''); ap.add_argument('--suffix', default=''); ap.add_argument('runs', nargs='+')
ns = ap.parse_args()
extra = json.loads((HERE / 'extra_cases.json').read_text())['cases']
std = json.loads((S / 'cases.json').read_text())
out = pathlib.Path(ns.out).resolve(); out.mkdir(parents=True, exist_ok=True)


def case_args(case):
    if case in extra:
        return list(extra[case]['args'])
    c = std['cases'][case]
    root, rel = c['gbk'].split('/', 1)
    return ['circular', '--gbk', f'{ROOTS[root]}/{rel}', *std['common'], '--label_font_size', str(c['font'])]


def one(run):
    case, design = run.split(':')
    args = case_args(case)
    name = f'{case}__{design}{ns.suffix}'
    cmd = ['/home/kawato/micromamba/bin/python3.13', '-P', str(S / 'run_case.py'), '--out', str(out), '--name', name]
    if design == 'bare':      # the same map without labels: the other objects, for scripts/lp_objects.py
        kept, i = [], 0
        while i < len(args):
            if args[i] == '--labels':
                i += 2 if i + 1 < len(args) and not args[i + 1].startswith('-') else 1
                continue
            kept.append(args[i]); i += 1
        args = kept + ['--labels', 'none']
    elif design == 'R':
        args += ['--label_placement', 'radial']
    elif design not in ('cur', 'bare'):
        cmd += ['--design', design, '--opts', ns.opts]
    proc = subprocess.run(cmd + ['--', *args], capture_output=True, text=True,
                          env={**os.environ, 'PYTHONPATH': WT, 'TMPDIR': '/home/kawato/gbdraw-baselines/tmp'})
    (out / f'{name}.log').write_text(proc.stdout + '\n--- stderr ---\n' + proc.stderr)
    mpath = out / f'{name}.metrics.json'
    if proc.returncode != 0 or not mpath.exists():
        return f'{name}: FAILED rc={proc.returncode}: ' + (proc.stderr.strip().splitlines() or ['?'])[-1][:300]
    m = json.loads(mpath.read_text())
    keys = ('labels_external', 'labels_embedded', 'font_size', 'text_overprint_pairs', 'leader_crossing_pairs', 'leader_len_max', 'canvas_px', 'seconds')
    return f'{name}: ' + ' '.join(f'{k}={m.get(k)}' for k in keys)


with ThreadPoolExecutor(ns.jobs) as pool:
    for line in pool.map(one, ns.runs):
        print(line, flush=True)
