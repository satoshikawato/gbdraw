#!/usr/bin/env python3
"""LP-D28: check that the study's Gallery cases still equal the current Gallery recipes on origin/dev.

For every Circular item of gbdraw/web/gallery/examples.json (origin/dev): compare its `command` token by token
with the args of the study case (page/extra_cases.json), paths reduced to basenames, and compare every input file
of the case byte for byte with the same path on origin/dev. usage: gq_check_recipes.py
"""
import json, pathlib, shlex, subprocess, sys
REPO = '/home/kawato/gbdraw-work'
T = pathlib.Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
DT = T / 'work' / 'dev-tree'
CASE = {'HmmtDNA_basic_circular': 'gal_hmmtdna_basic_circular', 'HmmtDNA_ATskew': 'gal_hmmtdna_atskew',
        'tobacco-chloroplast': 'gal_tobacco_chloroplast', 'Vnig_TUMSAT-TG-2018': 'gal_vnig_tumsat_tg_2018'}


def git(*a, binary=False):
    out = subprocess.run(['git', '-C', REPO, *a], capture_output=True, check=True).stdout
    return out if binary else out.decode()


print('origin/dev', git('rev-parse', '--short', 'origin/dev').strip(), '| dev-tree', (DT / 'COMMIT').read_text().strip()[:8])
items = json.loads(git('show', 'origin/dev:gbdraw/web/gallery/examples.json'))
extra = json.loads((T / 'page' / 'extra_cases.json').read_text())['cases']
bad = 0
for it in items:
    mode = 'circular' if it['command'].startswith('gbdraw circular') else 'linear'
    if mode != 'circular':
        print(f"{it['id']:32s} linear (not touched by Circular label placement)")
        continue
    toks = shlex.split(it['command'])[1:]
    want, i = [], 0
    while i < len(toks):
        if toks[i] in ('-o', '-f'):
            i += 2
            continue
        want.append(toks[i]); i += 1
    case = CASE[it['id']]
    have = extra[case]['args']
    norm = [pathlib.Path(a).name if a.startswith('/') else a for a in have]
    same = norm == want
    files = [a for a in have if a.startswith('/')]
    fres = []
    for f in files:
        p = pathlib.Path(f)
        try:
            rel = p.relative_to(DT)
            ok = p.read_bytes() == git('show', f'origin/dev:{rel}', binary=True)
        except Exception as e:  # noqa: BLE001
            ok = False; rel = f'{p} ({type(e).__name__})'
        fres.append((str(rel), ok))
    labels = [toks[k + 1] if k + 1 < len(toks) and not toks[k + 1].startswith('-') else '(flag only)' for k, t in enumerate(toks) if t == '--labels']
    print(f"{it['id']:32s} case {case:30s} args equal: {same} | inputs equal to origin/dev: {all(o for _, o in fres)} ({len(fres)}) | --labels {labels or 'absent'} | placement {'radial' if 'radial' in toks else 'default'}")
    if not same:
        bad += 1
        for a, b in zip(norm + ['<end>'] * 50, want + ['<end>'] * 50):
            if a != b:
                print('    first difference: case', a, '| gallery', b); break
    for rel, ok in fres:
        if not ok:
            bad += 1; print('    input differs or missing:', rel)
sys.exit(1 if bad else 0)
