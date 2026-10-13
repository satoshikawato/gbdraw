"""Add the leader turn fields (LP-D11: how far each leader turns) to metrics.json files rendered before the metric existed.
Usage (in the worktree, PYTHONPATH=$PWD): backfill_turn.py DIR [DIR ...]"""
import json, sys
from pathlib import Path
T = Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
sys.path.insert(0, str(T / 'scripts'))
import lp_metrics
done = 0
for d in sys.argv[1:]:
    for mj in sorted(Path(d).glob('*.metrics.json')):
        lj = Path(str(mj)[:-len('.metrics.json')] + '.labels.json')
        if not lj.exists():
            continue
        data = json.loads(lj.read_text()); met = json.loads(mj.read_text())
        ext = [l for l in data.get('labels', []) if not l.get('is_embedded')]
        try:
            leaders = [lp_metrics.label_leader(l) for l in ext]
        except Exception as error:
            print('skip', mj.name, type(error).__name__, error); continue
        met.update(lp_metrics.turn_fields(ext, leaders))
        mj.write_text(json.dumps(met, indent=1)); done += 1
print('backfilled', done)
