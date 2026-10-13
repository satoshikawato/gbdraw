#!/usr/bin/env python3
"""Table J2 -> P per case (brief P1). usage: table_p.py [--md OUT.md] A=DIR:SUFFIX B=DIR:SUFFIX   (A = J2, B = P)
Columns: straight, leader length, lopsided, mixed, off a contact point, wrong direction, far half, ink gap over 5 px, canvas height,
label gap min, last piece of bent leaders (median, max), seconds. Sums over the cases of the first set."""
import glob, json, math, os, sys
T = '/home/kawato/gbdraw-baselines/label-placement-20261011'
sys.path.insert(0, f'{T}/scripts'); sys.path.insert(0, f'{T}/work/J')
import lop, fan_order, port_check, last_piece

args = sys.argv[1:]; out = None
if '--md' in args:
    out = args[args.index('--md') + 1]; del args[args.index('--md'):args.index('--md') + 2]
sets = [(a.split('=')[0], *a.split('=')[1].rsplit(':', 1)) for a in args]
cases = sorted(os.path.basename(f).split('__')[0] for f in glob.glob(f'{T}/{sets[0][1]}/*__{sets[0][2]}.metrics.json'))

def row(d, suf, c):
    base = f'{T}/{d}/{c}__{suf}'
    if not os.path.exists(base + '.metrics.json'):
        return None
    m = json.load(open(base + '.metrics.json')); lp = base + '.labels.json'
    cl = lop.clusters(lp); mixed = sum(fan_order.mixed(x['s']) for x in cl)
    pr = port_check.check(lp); off = sum(r[1] is None for r in pr); wrong = sum(r[2] is False for r in pr)
    ls = sorted(n for _, n, _ in last_piece.rows(lp))
    kink = sum(1 for l in json.load(open(lp))['labels'] if not l.get('is_embedded') and l.get('leader_points') and 0.5 <= abs(lop.turns(l['leader_points'])) < 8.0)
    lastp = f'{ls[len(ls) // 2]:.0f}/{ls[-1]:.0f}' if ls else '-'
    return dict(drawn=m['leaders_drawn'], straight=m['leaders_straight'], length=m['leader_len_sum'], lop=sum(x['lop'] for x in cl),
                mixed=mixed, clusters=len(cl), off=off, wrong=wrong, far=m['leaders_end_far_half'], over5=m['leader_end_gap_ink_over5'],
                height=m['canvas_px'][1], gap=m['label_gap_min_px'], last=lastp, nb=len(ls), lastmax=ls[-1] if ls else 0, sec=m['seconds'], kink=kink,
                hard=m.get('labels_back_over_height', 0))

lines = ['| case | ' + ' | '.join(f'{n}: straight / length / lopsided / mixed / off pt / wrong dir / far half / ink>5 / canvas h / gap min / last piece med/max / turns 0.5-8 deg / s' for n, _d, _s in sets) + ' |',
         '|---|' + '---|' * len(sets)]
tot = {n: {} for n, _d, _s in sets}
for c in cases:
    cells = []
    for n, d, suf in sets:
        r = row(d, suf, c)
        if r is None:
            cells.append('-'); continue
        for k, v in r.items():
            if isinstance(v, (int, float)):
                tot[n][k] = tot[n].get(k, 0) + v
        tot[n]['lastmax'] = max(tot[n].get('lastmax_', 0), r['lastmax']); tot[n]['lastmax_'] = tot[n]['lastmax']
        cells.append(f"{r['straight']}/{r['drawn']} / {r['length']:.0f} / {r['lop']} / {r['mixed']}/{r['clusters']} / {r['off']} / {r['wrong']} / {r['far']} / {r['over5']} / {r['height']} / {r['gap']} / {r['last']} / {r['kink']} / {r['sec']}")
    lines.append(f'| {c} | ' + ' | '.join(cells) + ' |')
cells = []
for n, _d, _s in sets:
    t = tot[n]
    cells.append(f"{t['straight']}/{t['drawn']} / {t['length']:.0f} / {t['lop']} / {t['mixed']}/{t['clusters']} / {t['off']} / {t['wrong']} / {t['far']} / {t['over5']} / {t['height']:.0f} / min {'-'} / bent {t['nb']} / {t['kink']} / {t['sec']:.0f}")
lines.append('| **SUM** | ' + ' | '.join(cells) + ' |')
text = '\n'.join(lines)
print(text)
if out:
    open(out, 'w').write(text + '\n')
