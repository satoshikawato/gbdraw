#!/usr/bin/env python3
"""Render one study case (page/extra_cases.json) with extra gbdraw arguments appended (LP-D29: label gap sweep).

usage: run_case_plus.py --out DIR --name NAME [--design J --opts k=v,...] CASE -- <extra gbdraw args>
"""
import json, os, pathlib, subprocess, sys
HERE = pathlib.Path(__file__).resolve().parent
WT = os.environ.get('LP_WT', '/home/kawato/gbdraw-work/.worktrees/label-study')
a = sys.argv[1:]; cut = a.index('--'); head, extra_args = a[:cut], a[cut + 1:]
case = head[-1]; passed = head[:-1]
args = json.loads((HERE / 'extra_cases.json').read_text())['cases'][case]['args']
cmd = ['/home/kawato/micromamba/bin/python3.13', '-P', str(HERE.parent / 'scripts' / 'run_case.py'), *passed, '--', *args, *extra_args]
proc = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': WT, 'TMPDIR': '/home/kawato/gbdraw-baselines/tmp'})
out = pathlib.Path(passed[passed.index('--out') + 1]); name = passed[passed.index('--name') + 1]
(out / f'{name}.log').write_text(proc.stdout + '\n--- stderr ---\n' + proc.stderr)
print(name, 'rc', proc.returncode, (proc.stderr.strip().splitlines() or [''])[-1][:200] if proc.returncode else '')
sys.exit(proc.returncode)
