"""Summarize stored S08 observations; never execute or repeat a measurement."""
import gzip
import hashlib
import json
from pathlib import Path
import statistics

DATA = Path(__file__).resolve().parent
ROOT = DATA.parents[4]


def read(name):
    raw = (DATA / name).read_bytes()
    return json.loads(gzip.decompress(raw) if name.endswith('.gz') else raw)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def partition(row):
    """Disjoint wall-time boundaries, plus explicitly nested init observations."""
    messages = row['messages']
    render, = [m for m in messages if m['type'] == 'run'
               and 'diagram-generation-worker.js' in m['worker']]
    conversions = [m for m in messages if m.get('operation') == 'convertLosatpPairsToGenomicPayload']
    assert len(conversions) <= 1 and len(row['searches']) <= 1
    conversion = conversions[0] if conversions else None
    search = row['searches'][0] if row['searches'] else None
    assert not search or conversion
    first = search or conversion or render
    parts = {
        'runtimeInputIdentityPreparation': (first['start'] - row['start']) / 1000,
        'rawSearchBatch': (search['end'] - search['start']) / 1000 if search else 0,
        'rawPacking': (conversion['start'] - search['end']) / 1000 if search else 0,
        'conversionRoundTrip': (conversion['end'] - conversion['start']) / 1000 if conversion else 0,
        'typedToRenderRequest': (render['start'] - conversion['end']) / 1000 if conversion else 0,
        'renderRoundTrip': (render['end'] - render['start']) / 1000,
        'resultHistoryDomFrames': (row['end'] - render['end']) / 1000,
    }
    assert min(parts.values()) >= 0
    assert abs(sum(parts.values()) - row['seconds']) < 1e-8
    return {
        'seconds': parts,
        'total': row['seconds'],
        'nestedInitializations': [
            {'worker': m['worker'].rsplit('/', 1)[-1],
             'seconds': (m['end'] - m['start']) / 1000,
             'startAfterGenerateSeconds': (m['start'] - row['start']) / 1000}
            for m in messages if m['type'] == 'init'],
        'telemetry': row['telemetry'],
        'sourceJobs': sum(len(s['jobs']) for s in row['searches']),
        'completedJobs': sum(s['completed'] for s in row['searches']),
        'observedLosatRunMessages': sum(m['type'] == 'run' and 'losat-worker.js' in m['worker'] for m in messages),
        'contextWorkersConstructed': row['workerCount'],
    }


def main():
    evidence = read('s08-browser-timing.json.gz')
    artifacts = read('s08-artifacts.json')
    assert evidence['wheelSha256'] == artifacts['wheelSha256']
    assert evidence['runnerSha256'] == sha(ROOT / 'tests/run_protein_local_browser_acceptance.py')
    for path, expected in evidence['sourceSha256'].items():
        assert sha(ROOT / path) == expected, path
    for path, expected in artifacts['installedFilesMatched'].items():
        assert sha(ROOT / path) == expected, path
        assert sha(ROOT / '.venv/s08-cli/lib/python3.13/site-packages' / path) == expected, path
    assert evidence['warmups'] == 0
    assert len(evidence['runs']) == 16
    assert all(not r['external'] and not r['pageErrors'] and r['liveWorkersAfterDispose'] == 0
               for r in evidence['runs'])
    out = {
        'evidenceSha256': sha(DATA / 's08-browser-timing.json.gz'),
        'artifactLedgerSha256': sha(DATA / 's08-artifacts.json'),
        'sourceAndInstalledFilesReverified': True,
        'units': 'seconds',
        'boundary': 'Generate action through SVG DOM and two animation frames',
        'nestedIntervalsAreNotAdded': True,
        'browserBaseline': None,
        'caseDefinitions': {},
        'cases': {},
    }
    for case, values in evidence['cases'].items():
        first = next(r for r in evidence['runs'] if r['case'] == case)
        out['caseDefinitions'][case] = {
            'fixtureSha256': values['fixtureSha256'],
            'records': first['preparation']['records'],
            'settings': first['preparation']['settings'],
            'rawAndSavedGeometryEqual': values['raw-cold'][0]['geometrySha256'] == values['saved-raw'][0]['geometrySha256'],
        }
        out['cases'][case] = {}
        for state in ('raw-cold', 'saved-raw', 'derived-warm'):
            rows = values[state]
            assert len(rows) == (1 if state == 'raw-cold' else 3)
            samples = [r['seconds'] for r in rows]
            median = statistics.median(samples)
            assert len({r['geometrySha256'] for r in rows}) == 1
            index = min(range(len(rows)), key=lambda i: abs(samples[i] - median))
            out['cases'][case][state] = {
                'samples': samples,
                'median': median,
                'mad': statistics.median(abs(v - median) for v in samples) if len(rows) > 1 else None,
                'representativeSampleIndex': index,
                'representative': partition(rows[index]),
                'allPartitions': [partition(r) for r in rows],
            }
            print(case, state, json.dumps(out['cases'][case][state]['samples']), 'median', median)
    (DATA / 's08-browser-summary.json').write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
