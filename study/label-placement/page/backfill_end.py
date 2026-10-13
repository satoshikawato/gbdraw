"""Add the LP-D15 measures (where a leader ends: leaders_end_far_half, leader_end_gap_*, leader_ends_unclear) to
metrics.json files rendered before they existed.
Usage (PYTHONPATH=<study worktree>): backfill_end.py [--only SUBSTRING] DIR [DIR ...]"""
import json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
T = Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
sys.path.insert(0, str(T / 'scripts'))
import lp_metrics
from gbdraw.labels import circular as lc


def one(mj):
    mj = Path(mj)
    stem = str(mj)[:-len('.metrics.json')]
    lj, svg = Path(stem + '.labels.json'), Path(stem + '.svg')
    if not lj.exists():
        return f'skip {mj.name}: no labels.json'
    data = json.loads(lj.read_text()); met = json.loads(mj.read_text())
    ext = [l for l in data.get('labels', []) if not l.get('is_embedded')]
    try:
        boxes = [lp_metrics.label_box(l, data['total_length'], lc) for l in ext]
        leaders = [lp_metrics.label_leader(l) for l in ext]
    except Exception as error:
        return f'skip {mj.name}: {type(error).__name__} {error}'
    rotated = [abs(float(l.get('rotation_deg') or 0.0)) > 1e-6 and l.get('oriented_corners') is not None for l in ext]
    met.update(lp_metrics.end_fields(ext, boxes, leaders, rotated, str(svg) if svg.exists() else None))
    mj.write_text(json.dumps(met, indent=1))
    return None


if __name__ == '__main__':
    args = sys.argv[1:]
    only = None
    if args and args[0] == '--only':
        only, args = args[1], args[2:]
    files = [str(mj) for d in args for mj in sorted(Path(d).glob('*.metrics.json')) if only is None or only in mj.name]
    with ProcessPoolExecutor(4) as pool:
        notes = [n for n in pool.map(one, files) if n]
    print('\n'.join(notes)); print('backfilled', len(files) - len(notes), 'of', len(files))
