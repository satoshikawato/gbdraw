#!/usr/bin/env python3
"""Stage check for inner labels (scope study): does an overlap right after the inner arc placement predict the final
drawn inner overlaps and crossings, and how long does each inner stage take?
usage: scope_stage.py CASE[:VARIANT] ...   (cases and variants of scope_sweep.py; prints one JSON line per run)
Counts use the scope definitions (boxes overlap > 0.5 px in x and y; proper crossing of leader polylines feature
anchor -> elbow -> contact). Stages: A = rearrange_labels_fc(inner) result, B = after _rebalance_inner_labels_strict_order,
C = the final label list returned by prepare_label_list (inner labels only), D = drawn (assemble may move labels later).
"""
import json, math, os, pathlib, sys, time
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault('TMPDIR', '/home/kawato/gbdraw-baselines/tmp')
import importlib.util
spec = importlib.util.spec_from_file_location('sweep_variant', HERE / 'scope_sweep.py')
src = (HERE / 'scope_sweep.py').read_text()
ns = {}
exec(src[src.index('OPT = {'):src.index('def one(')], ns)
variant = ns['variant']
cases = json.loads((HERE / 'scope_cases.json').read_text())['cases']
import gbdraw.labels.circular as lc


def boxes_and_leads(labels, T):
    bx, ld = [], []
    for l in labels:
        bx.append(lc._label_bbox_rect(l, T, 0.0))
        ld.append([(float(l.get('feature_anchor_x', l['feature_middle_x'])), float(l.get('feature_anchor_y', l['feature_middle_y']))),
                   (float(l['middle_x']), float(l['middle_y'])),
                   (float(l.get('leader_start_x', l['start_x'])), float(l.get('leader_start_y', l['start_y'])))])
    return bx, ld


def _cross(p, q, r, s):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2, d3, d4 = o(r, s, p), o(r, s, q), o(p, q, r), o(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-6 for d in (d1, d2, d3, d4))


def counts(labels, T):
    labels = [l for l in labels if l.get('is_inner') and not l.get('is_embedded')]
    bx, ld = boxes_and_leads(labels, T)
    n = len(labels); ov = cr = 0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = bx[i], bx[j]
            if min(a[2], b[2]) - max(a[0], b[0]) > 0.5 and min(a[3], b[3]) - max(a[1], b[1]) > 0.5:
                ov += 1
            if any(_cross(p, q, r, s) for p, q in zip(ld[i], ld[i][1:]) for r, s in zip(ld[j], ld[j][1:])):
                cr += 1
    return [ov, cr]


for run in sys.argv[1:]:
    case, _, var = run.partition(':')
    args = variant(cases[case], var) if var else list(cases[case])
    stage = {}
    T_holder = {}
    o_re = lc.rearrange_labels_fc
    o_rb = lc._rebalance_inner_labels_strict_order
    o_pl = lc.prepare_label_list

    def re_w(labels, *a, **k):
        t0 = time.perf_counter(); r = o_re(labels, *a, **k)
        if not k.get('is_outer', a[5] if len(a) > 5 else True):
            stage['A'] = counts([dict(x) for x in r], a[1]); stage['tA'] = round(time.perf_counter() - t0, 2)
        return r

    def rb_w(labels, T, *a, **k):
        t0 = time.perf_counter(); r = o_rb(labels, T, *a, **k)
        stage['B'] = counts([dict(x) for x in r], T); stage['tB'] = round(time.perf_counter() - t0, 2)
        return r

    def pl_w(fd, T, *a, **k):
        t0 = time.perf_counter(); r = o_pl(fd, T, *a, **k)
        stage['C'] = counts([dict(x) for x in r], T); stage['tAll'] = round(time.perf_counter() - t0, 2)
        T_holder['T'] = T; T_holder['labels'] = r
        return r
    lc.rearrange_labels_fc = re_w; lc._rebalance_inner_labels_strict_order = rb_w
    import gbdraw.diagrams.circular.assemble as asm
    asm.prepare_label_list = pl_w
    from gbdraw.cli import main
    out = pathlib.Path.cwd()
    os.chdir(out)
    sys.argv = ['gbdraw', *args, '-o', 'stage', '-f', 'svg', '--overwrite']
    try:
        main()
    except SystemExit:
        pass
    if 'labels' in T_holder:
        stage['D'] = counts(T_holder['labels'], T_holder['T'])
    lc.rearrange_labels_fc = o_re; lc._rebalance_inner_labels_strict_order = o_rb; asm.prepare_label_list = o_pl
    print(json.dumps({'run': run, **stage}), flush=True)
