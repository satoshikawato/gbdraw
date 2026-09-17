"""Regenerate S07 comparisons from raw archived and current measurement reports."""
import gzip
import hashlib
import json
from pathlib import Path
import statistics

DATA = Path(__file__).resolve().parent


def read(name):
    path = DATA / name
    return json.loads(gzip.decompress(path.read_bytes()))


def compare(base_name, current_name, names, stages=None):
    base, current = read(base_name), read(current_name)
    assert base['dependencies'] == current['dependencies']
    rows = []
    for name in names:
        before, after = base['cases'][name], current['cases'][name]
        assert before['fixture'] == after['fixture'], name
        keys = stages if stages is not None else before['stages']
        for stage in keys:
            b, a = before['stages'][stage], after['stages'][stage]
            assert b['semanticSha256'] == a['semanticSha256'], (name, stage)
            assert set(b['sampleSemanticSha256']) == set(a['sampleSemanticSha256'])
            row = {'case': name, 'stage': stage, 'baselineReport': base_name,
                   'currentReport': current_name, 'semanticEqual': True,
                   'semanticSha256': a['semanticSha256'], 'outputCanonicalBytes': a['semanticBytes'],
                   'before': b['median'], 'after': a['median'],
                   'deltaPct': 100 * (a['median'] / b['median'] - 1)}
            if base['measurement'] == 'timing':
                assert current['measurement'] == 'timing'
                assert len(b['samples']) == len(a['samples']) >= 7
                row['samplesPerSource'] = len(a['samples'])
                for label, values in [('baseline', b), ('current', a)]:
                    median = statistics.median(values['samples'])
                    mad = statistics.median(abs(x - median) for x in values['samples'])
                    assert median == values['median']
                    row[label + 'Mad'] = mad
                    row[label + 'NoisePct'] = 100 * mad / median
                row['medianRegression'] = row['deltaPct'] > 10
                row['decision'] = ('inconclusive' if max(row['baselineNoisePct'], row['currentNoisePct']) > 5
                                   else 'regression' if row['medianRegression'] else 'pass')
            else:
                assert base['measurement'] == current['measurement'] == 'memory'
                row['beforeRetainedBytes'] = b.get('retainedBytes')
                row['afterRetainedBytes'] = a.get('retainedBytes')
            rows.append(row)
    return rows


def main():
    merges = ['merge-300', 'merge-600', 'merge-1200', 'merge-edges']
    gallery = ['gallery-collinear', 'gallery-collinear-off']
    vibrio = ['gallery-collinear-vibrio', 'gallery-collinear-vibrio-off']
    timing = []
    for b, a, names, stages in [
        ('baseline-timing.json.gz', 's07-current-cluster-timing.json.gz', merges, None),
        ('s07-baseline-isolated-timing.json.gz', 's07-current-isolated-timing.json.gz', merges + gallery, None),
        ('s06-current-timing-final.json.gz', 's07-current-gallery-timing.json.gz', gallery, ['post_search']),
        ('s07-vibrio-baseline-timing.json.gz', 's07-vibrio-current-timing.json.gz', vibrio, None),
    ]:
        timing.extend(compare(b, a, names, stages))
    memory = []
    for b, a, names, stages in [
        ('memory.json.gz', 's07-current-memory.json.gz', merges, None),
        ('s07-baseline-isolated-memory.json.gz', 's07-current-memory.json.gz', merges + gallery, None),
        ('s06-current-memory-final.json.gz', 's07-current-memory.json.gz', gallery, ['post_search']),
        ('s07-vibrio-baseline-memory.json.gz', 's07-vibrio-current-memory.json.gz', vibrio, None),
    ]:
        memory.extend(compare(b, a, names, stages))
    b = read('s07-baseline-final-counters.json.gz')
    a = read('s07-current-final-counters.json.gz')
    counters = []
    for name, case in b['cases'].items():
        assert case['fixture'] == a['cases'][name]['fixture']
        for stage, before in case['stages'].items():
            after = a['cases'][name]['stages'][stage]
            assert before['semanticSha256'] == after['semanticSha256']
            counters.append({'case': name, 'stage': stage,
                'before': {k: v for k, v in before['operations'].items() if k.startswith('merge.')},
                'after': {k: v for k, v in after['operations'].items() if k.startswith('merge.')}})
    # Normal post-search processing must not return to exhaustive path traversal.
    for filename in ['s07-current-probe.json.gz', 's07-vibrio-current-probe.json.gz']:
        for case in read(filename)['cases'].values():
            for stage in case['stages'].values():
                assert all(v == 0 for k, v in stage['operations'].items() if k.startswith('paths.'))
    report = {'timing': timing, 'memory': memory, 'counters': counters,
              'performanceAcceptance': 'accepted-by-user',
              'acceptanceReceipt': 's07-acceptance.json',
              'independentRemeasurementPerformed': False,
              'limitations': 'Archived full-stage observations retain S01/S06 drift caveats. External Git, archive/compression and editor processes were observed. Instrumentation, memory and timing are separate runs; no S07 build/test/other benchmark overlaps S07 timing.',
              'normalPathEnumerationCalls': 0}
    filenames = {r[k] for rows in (timing, memory) for r in rows for k in ('baselineReport', 'currentReport')}
    filenames.update(['s07-acceptance.json', 's07-baseline-final-counters.json.gz', 's07-current-final-counters.json.gz',
                      's07-current-probe.json.gz', 's07-vibrio-current-probe.json.gz'])
    hotspots = []
    for filename in ['s07-current-probe.json.gz', 's07-vibrio-current-probe.json.gz']:
        for name, case in read(filename)['cases'].items():
            if not name.startswith('gallery-'):
                continue
            rows = {r['name']: r for r in case['stages']['post_search']['profile']}
            total = rows['analyze']['cumulativeSeconds']
            functions = (['_normalize_directional_hit_tables', '_build_anchor_core_orthogroups']
                         if not name.endswith('-off') else
                         ['build_collinearity_unit_index', '_select_member_candidate_hits_per_query',
                          '_select_orthogroup_edges'])
            hotspots.append({'case': name, 'sourceReport': filename,
                'instrumentedTotalSeconds': total,
                'functions': [{'name': f, 'cumulativeSeconds': rows[f]['cumulativeSeconds'],
                               'profilePct': 100 * rows[f]['cumulativeSeconds'] / total}
                              for f in functions]})
    report['profileHotspots'] = hotspots
    report['profileCaveat'] = 'Instrumented cumulative function shares, not uninstrumented timing percentages. Nested child costs must not be added to their parents.'
    report['reportSha256'] = {n: hashlib.sha256((DATA / n).read_bytes()).hexdigest() for n in sorted(filenames)}
    (DATA / 's07-comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    for row in timing:
        print(row['case'], row['stage'], round(row['before'], 4), round(row['after'], 4),
              round(row['deltaPct'], 2), row['decision'])


if __name__ == '__main__':
    main()
