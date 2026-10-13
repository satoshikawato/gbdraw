#!/usr/bin/env python3
"""LP-D30: the same cases at larger label fonts, prototype Q round 3 with and without the font-relative row rule.

Each (case, font) is rendered twice: <case>__J_Q3_f<font> (port=1,arcfirst=3) and <case>__J_QA_f<font> (…,rowem=0.22).
The name starts with the case, so scripts/lp_objects.py finds the case's bare twin. usage: run_font_sweep.py OUT FONT[,FONT...] CASE ...
"""
import json, os, pathlib, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
HERE = pathlib.Path(__file__).resolve().parent; S = HERE.parent / 'scripts'
WT = os.environ.get('LP_WT', '/home/kawato/gbdraw-work/.worktrees/label-study')
ROOTS = {'TI': f'{WT}/tests/test_inputs', 'EX': '/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs',
         'IN': str(HERE.parent / 'work' / 'inputs')}
extra = json.loads((HERE / 'extra_cases.json').read_text())['cases']; std = json.loads((S / 'cases.json').read_text())
VARIANTS = (('Q3', 'port=1,arcfirst=3'), ('QA', 'port=1,arcfirst=3,rowem=0.22'), ('QB', 'port=1,arcfirst=3,rowem=0.22,rowmin=4'))
if os.environ.get('LP_VARIANTS'):   # e.g. LP_VARIANTS=QB
    VARIANTS = tuple(v for v in VARIANTS if v[0] in os.environ['LP_VARIANTS'].split(','))
out = pathlib.Path(sys.argv[1]).resolve(); out.mkdir(parents=True, exist_ok=True)
fonts = sys.argv[2].split(','); cases = sys.argv[3:]


def base(case):
    if case in extra:
        return list(extra[case]['args'])
    c = std['cases'][case]; root, rel = c['gbk'].split('/', 1)
    return ['circular', '--gbk', f'{ROOTS[root]}/{rel}', *std['common']]


def one(job):
    case, font, tag, opts = job
    name = f'{case}__J_{tag}_f{font}'
    cmd = ['/home/kawato/micromamba/bin/python3.13', '-P', str(S / 'run_case.py'), '--out', str(out), '--name', name,
           '--design', 'J', '--opts', opts, '--', *base(case), '--label_font_size', font]
    p = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': WT, 'TMPDIR': '/home/kawato/gbdraw-baselines/tmp'})
    (out / f'{name}.log').write_text(p.stdout + '\n--- stderr ---\n' + p.stderr)
    return f'{name}: rc={p.returncode}' + ('' if p.returncode == 0 else ' ' + (p.stderr.strip().splitlines() or ['?'])[-1][:200])


jobs = [(c, f, t, o) for c in cases for f in fonts for t, o in VARIANTS]
with ThreadPoolExecutor(4) as pool:
    for line in pool.map(one, jobs):
        print(line)
