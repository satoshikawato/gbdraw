#!/usr/bin/env python3
"""LP-D28: render study cases (page/extra_cases.json) with an unmodified gbdraw tree (today's origin/dev).

usage: PYTHONPATH=<tree> python3.13 -P gq_dev_render.py OUT_DIR FORMAT[,FORMAT] CASE [CASE ...]
Writes OUT_DIR/<case>__dev.<ext>. Relative input names are looked up in <tree>/tests/test_inputs.
"""
import json, os, pathlib, sys, time
HERE = pathlib.Path(__file__).resolve().parent
os.environ.setdefault('TMPDIR', '/home/kawato/gbdraw-baselines/tmp')
out = pathlib.Path(sys.argv[1]).resolve(); out.mkdir(parents=True, exist_ok=True)
fmt, cases = sys.argv[2], sys.argv[3:]
extra = json.loads((HERE / 'extra_cases.json').read_text())['cases']
import gbdraw
from gbdraw.cli import main
tree = pathlib.Path(gbdraw.__file__).resolve().parent.parent
print('tree', tree, flush=True)
os.chdir(out)
for case in cases:
    args = [str(tree / 'tests' / 'test_inputs' / a) if a.endswith(('.gbff', '.gb', '.gbk')) and '/' not in a else a for a in extra[case]['args']]
    sys.argv = ['gbdraw', *args, '-o', f'{case}__dev', '-f', fmt, '--overwrite']
    t0 = time.time()
    try:
        rc = main()
    except SystemExit as e:
        rc = e.code
    print(case, 'rc', rc, 'seconds', round(time.time() - t0, 1), flush=True)
