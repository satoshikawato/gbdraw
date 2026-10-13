"""Add the LP-D13 measure of LP-D09 (label_back_px_max, labels_back_over_height, label_back_list) to metrics.json
files rendered before it existed. Usage (in a study worktree, PYTHONPATH=$PWD): backfill_back.py DIR [DIR ...]"""
import json, sys
from pathlib import Path
T = Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
sys.path.insert(0, str(T / 'scripts'))
import lp_metrics
from gbdraw.labels import circular as lc
done = 0
for d in sys.argv[1:]:
    for mj in sorted(Path(d).glob('*.metrics.json')):
        lj = Path(str(mj)[:-len('.metrics.json')] + '.labels.json')
        if not lj.exists():
            continue
        data = json.loads(lj.read_text()); met = json.loads(mj.read_text())
        ext = [l for l in data.get('labels', []) if not l.get('is_embedded')]
        try:
            boxes = [lp_metrics.label_box(l, data['total_length'], lc) for l in ext]
            leaders = [lp_metrics.label_leader(l) for l in ext]
        except Exception as error:
            print('skip', mj.name, type(error).__name__, error); continue
        rotated = [abs(float(l.get('rotation_deg') or 0.0)) > 1e-6 and l.get('oriented_corners') is not None for l in ext]
        back = lp_metrics.back_px(ext, boxes, leaders, rotated)
        met['label_back_px_max'] = back[0][0] if back else 0.0
        met.pop('labels_back_over_half_height', None)
        met['labels_back_over_height'] = sum(1 for _b, over, _t in back if over)
        met['label_back_list'] = [f'{t}: {b}' for b, _o, t in back[:10] if b > 0]
        mj.write_text(json.dumps(met, indent=1)); done += 1
print('backfilled', done)
