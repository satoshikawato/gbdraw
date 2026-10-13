"""Add label_gap_min_px (smallest distance between two label boxes) to metrics.json files rendered before the metric existed."""
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
        except Exception as error:
            print('skip', mj.name, type(error).__name__, error); continue
        rotated = {i for i, l in enumerate(ext) if abs(float(l.get('rotation_deg') or 0.0)) > 1e-6}
        met['label_gap_min_px'] = lp_metrics.min_box_gap(boxes, rotated)
        mj.write_text(json.dumps(met, indent=1)); done += 1
print('backfilled', done)
