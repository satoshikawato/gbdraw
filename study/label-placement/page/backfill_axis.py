"""Add labels_toward_axis (LP-D09: labels that run back toward the ring) to metrics.json files rendered before the metric existed."""
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
        turned = lp_metrics.toward_axis(ext, boxes, leaders, rotated)
        met['labels_toward_axis'] = len(turned); met['labels_toward_axis_list'] = turned[:40]
        sharp, angle = lp_metrics.acute_bends(leaders)
        met['leader_acute_bends'] = len(sharp); met['leader_acute_list'] = [str(ext[i].get('label_text', '')) for i in sharp][:40]
        met['leader_min_angle_deg'] = angle
        tight, _ = lp_metrics.acute_bends(leaders, limit=110.0)
        met['leader_bends_lt110'] = len(tight); met['leader_bends_lt110_list'] = [str(ext[i].get('label_text', '')) for i in tight][:40]
        mj.write_text(json.dumps(met, indent=1)); done += 1
print('backfilled', done)
