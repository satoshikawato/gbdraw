"""Validate S07.6 evidence and regenerate the bounded comparison report.

Missing timing is explicitly pending, never a pass. Run the documented timing
commands only after an actual-host contention check and all local work settles.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import statistics

DATA = Path(__file__).resolve().parent


def read(name):
    raw = (DATA/name).read_bytes()
    return json.loads(gzip.decompress(raw) if name.endswith('.gz') else raw)


def stage(report, case, name):
    return report['cases'][case]['stages'][name]


def compare(base, current, case, name):
    assert base['dependencies'] == current['dependencies']
    assert base['cases'][case]['fixture'] == current['cases'][case]['fixture'], case
    b, c = stage(base, case, name), stage(current, case, name)
    assert b['semanticSha256'] == c['semanticSha256'], (case, name)
    for entry in (b, c):
        assert set(entry['sampleSemanticSha256']) == {b['semanticSha256']}
    result = {'case': case, 'stage': name, 'semanticSha256': b['semanticSha256'],
              'semanticBytes': b['semanticBytes'], 'outputEqual': True}
    assert b['semanticBytes'] == c['semanticBytes']
    if base['measurement'] == current['measurement'] == 'timing':
        required = 21 if name == 'post_search' and '-vibrio' not in case else 7
        for report, entry in ((base, b), (current, c)):
            assert report['settings']['warmups'] == 1
            assert len(entry['samples']) == required
            median = statistics.median(entry['samples'])
            mad = statistics.median(abs(x-median) for x in entry['samples'])
            assert median == entry['median'] and mad == entry['mad']
        delta = 100*(c['median']/b['median']-1)
        noisy = max(b['noisePct'], c['noisePct']) > 5
        result.update(beforeMs=b['median'], afterMs=c['median'], changePct=delta,
                      beforeNoisePct=b['noisePct'], afterNoisePct=c['noisePct'],
                      medianRegression=delta > 10,
                      decision='inconclusive' if noisy else 'regression' if delta > 10 else 'pass')
    elif base['measurement'] == current['measurement'] == 'memory':
        result.update(beforePeak=b['median'], afterPeak=c['median'],
                      beforeRetained=b['retainedBytes'], afterRetained=c['retainedBytes'])
    elif base['measurement'] == current['measurement'] == 'probe':
        result.update(beforeOperations=b['operations'], afterOperations=c['operations'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-prefix', default='s07-6-current')
    args = parser.parse_args()
    memory, counters, timings, pending = [], [], [], []
    evidence_files = set()
    def report(name):
        name = name.replace('s07-6-current-', args.candidate_prefix+'-')
        evidence_files.add(name)
        return read(name)
    candidate = report('s07-6-current-memory.json.gz')
    current_probe = report('s07-6-current-probe.json.gz')
    for kind, cases in [('gallery', ['gallery-collinear', 'gallery-collinear-off']),
                        ('vibrio', ['gallery-collinear-vibrio', 'gallery-collinear-vibrio-off'])]:
        base_memory = report('s07-current-memory.json.gz' if kind == 'gallery' else 's07-vibrio-current-memory.json.gz')
        base_timing = report('s07-current-gallery-timing.json.gz' if kind == 'gallery' else 's07-vibrio-current-timing.json.gz')
        for case in cases:
            memory.append(compare(base_memory, candidate, case, 'post_search'))
            compare(base_timing, current_probe, case, 'post_search')
        new_timing = f'{args.candidate_prefix}-{kind}-timing.json.gz'
        if (DATA/new_timing).exists():
            changed = report(new_timing)
            timings.extend(compare(base_timing, changed, case, 'post_search') for case in cases)
        else:
            pending.append(new_timing)
    baseline_units = report('s07-6-baseline-unit-memory.json.gz')
    current_units = report('s07-6-current-unit-memory.json.gz')
    baseline_probe = report('s07-6-baseline-probe.json.gz')
    for case in ('gallery-collinear', 'gallery-collinear-vibrio', 'units-multicds'):
        memory.append(compare(baseline_units, current_units, case, 'unit_index'))
        counters.append(compare(baseline_probe, current_probe, case, 'unit_index'))
    unit_names = ['s07-6-baseline-unit-timing.json.gz', f'{args.candidate_prefix}-unit-timing.json.gz']
    if all((DATA/name).exists() for name in unit_names):
        before, after = map(report, unit_names)
        timings.extend(compare(before, after, case, 'unit_index') for case in baseline_units['cases'])
    else:
        pending.extend(name for name in unit_names if not (DATA/name).exists())
    owner = 'gbdraw/analysis/collinearity_units.py'
    measured_sha = candidate['source']['filesSha256'][owner]
    final_sha = hashlib.sha256((DATA.parents[4]/owner).read_bytes()).hexdigest()
    summary = {'measuredProductionSha256': measured_sha, 'workingProductionSha256': final_sha,
               'measuredSourceMatchesWorkingCandidate': measured_sha == final_sha,
               'measurementDisposition': ('Final candidate measurements' if measured_sha == final_sha else
                                         'Earlier candidate observations; final follow-up measurements deferred by user'),
               'scope': 'S07.6 only; Collinear post_search and prepared extraction -> unit_index',
               'policy': {'regressionPct': 10, 'maxNoisePct': 5, 'outputDifference': 'fail regardless of speed'},
               'memory': memory, 'counters': counters, 'timing': timings, 'pendingTimingReports': pending,
               'performanceAcceptance': 'NOT_ESTABLISHED' if pending else 'see numerical decisions and host evidence',
               'evidenceSha256': {name: hashlib.sha256((DATA/name).read_bytes()).hexdigest() for name in sorted(evidence_files)}}
    (DATA/('s07-6-comparison.json' if args.candidate_prefix == 's07-6-current' else args.candidate_prefix+'-comparison.json')).write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'memoryPairs': len(memory), 'counterPairs': len(counters), 'timingPairs': len(timings), 'pending': pending}))


if __name__ == '__main__':
    main()
