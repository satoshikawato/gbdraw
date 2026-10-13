"""LP-D15 table: where leaders end, per case and design. Usage: end_table.py NAME=DIR:SUFFIX [NAME=DIR:SUFFIX ...]
(for example H3=work/J:J J2=work/J:J_sym2). Prints one row per case found in the first set, then the sums."""
import glob, json, os, sys
T = '/home/kawato/gbdraw-baselines/label-placement-20261011'
sets = [(a.split('=')[0], *a.split('=')[1].split(':')) for a in sys.argv[1:]]
first = sets[0]
cases = sorted(os.path.basename(f).split('__')[0] for f in glob.glob(f'{T}/{first[1]}/*__{first[2]}.metrics.json'))
load = lambda d, suf, c: json.load(open(f'{T}/{d}/{c}__{suf}.metrics.json')) if os.path.exists(f'{T}/{d}/{c}__{suf}.metrics.json') else None
tot = {n: dict(far=0, off=0, unclear=0, over5=0, drawn=0, n=0, mx=0.0) for n, _d, _s in sets}
print('| case | leaders | ' + ' | '.join(f'{n}: far half / off middle / unclear / ink gap median, max / over 5 px' for n, _d, _s in sets) + ' |')
print('|---|---|' + '---|' * len(sets))
for c in cases:
    cells, drawn = [], None
    for n, d, suf in sets:
        m = load(d, suf, c)
        if m is None or 'leaders_end_far_half' not in m:
            cells.append('-'); continue
        drawn = m.get('leaders_drawn')
        t = tot[n]; t['far'] += m['leaders_end_far_half']; t['off'] += m.get('leaders_end_off_middle') or 0; t['unclear'] += m['leader_ends_unclear'] or 0
        t['over5'] += m['leader_end_gap_ink_over5'] or 0; t['drawn'] += m.get('leaders_drawn', 0); t['n'] += 1
        t['mx'] = max(t['mx'], m['leader_end_gap_ink_max'] or 0.0)
        cells.append(f"{m['leaders_end_far_half']} / {m.get('leaders_end_off_middle', '-')} / {m['leader_ends_unclear']} / {m['leader_end_gap_ink_median']}, {m['leader_end_gap_ink_max']} / {m['leader_end_gap_ink_over5']}")
    print(f'| {c} | {drawn} | ' + ' | '.join(cells) + ' |')
print('| SUM | ' + ' | '.join([''] + [f"{t['far']} / {t['off']} / {t['unclear']} / max {t['mx']} / {t['over5']} (of {t['drawn']} leaders, {t['n']} cases)" for t in tot.values()]) + ' |')
