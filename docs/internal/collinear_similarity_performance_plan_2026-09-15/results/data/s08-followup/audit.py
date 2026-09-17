"""Read-only source/artifact/preservation checks; write only this follow-up's report."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[5]
S08 = ROOT.parent / 'collinear-s08-20260916'
INSTALLED = ROOT / '.venv/raw-validation/lib/python3.13/site-packages'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args, cwd=ROOT):
    return subprocess.check_output(['git', *args], cwd=cwd)


def read_report(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def main():
    inherited = json.loads((OUT / 'inheritance.json').read_text())
    assert git('rev-parse', 'HEAD').decode().strip() == inherited['inheritedHead']
    assert git('rev-parse', 'origin/dev').decode().strip() == inherited['base']
    names = git('ls-files', '-z', 'gbdraw').decode().split('\0')[:-1]
    source_names = [n for n in names if n.endswith(('.py', '.js'))]
    source = {n: sha((ROOT / n).read_bytes()) for n in source_names}
    source_delta = [n for n in names if (ROOT / n).read_bytes() != (S08 / n).read_bytes()]
    assert source_delta == ['gbdraw/web/js/app/losat-cache.js', 'gbdraw/web/js/app/run-analysis.js'], source_delta
    py_names = [n for n in source_names if n.endswith('.py')]
    prior = read_report(OUT.parent / 's08-replacement-packaged.json.gz')
    assert {n: source[n] for n in py_names} == prior['sourceSha256']
    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as z:
        for n in py_names:
            assert z.read(n) == (ROOT / n).read_bytes(), n
    installed_names = [n for n in names if (INSTALLED / n).is_file()]
    for n in installed_names:
        assert (INSTALLED / n).read_bytes() == (ROOT / n).read_bytes(), n
    assert (INSTALLED / 'gbdraw/web' / wheel.name).read_bytes() == wheel.read_bytes()
    preservation = {}
    for name, before in inherited['preservationStart'].items():
        where = ROOT.parents[1] / ('' if name == 'shared-dev' else name)
        after = {'head': git('rev-parse', 'HEAD', cwd=where).decode().strip(),
                 'diffSha256': sha(git('diff', '--binary', 'HEAD', cwd=where)),
                 'status': git('status', '--porcelain', '--untracked-files=normal', cwd=where).decode()}
        assert after == before, name
        preservation[name] = after
    lifecycle = read_report(OUT / 'lifecycle.json.gz')
    assert 'inProgress' not in lifecycle
    steps = lifecycle['runs'][0]['steps']
    old = read_report(OUT.parent / 's08-lifecycle-packaged-observed.json.gz')['inProgress']['steps']
    equality = []
    for current in steps:
        # The archived prefix's clear-cache was a plain typed repeat, with no
        # raw searches or derived provenance. It is not this operation.
        previous = next((row for row in old if row['label'] == current['label']
                         and row['label'] != 'clear-cache'), None)
        if previous:
            assert current['geometrySha256'] == previous['geometrySha256'], current['label']
            assert current['svgSha256'] == previous['svgSha256'], current['label']
            assert current['provenance'] == previous['provenance'], current['label']
            equality.append({'label': current['label'], 'geometrySha256': current['geometrySha256'],
                             'svgSha256': current['svgSha256'], 'provenanceEqual': True,
                             'searchBatchCount': len(current['searches'])})
    cleared = next(row for row in steps if row['label'] == 'clear-cache')
    retry = next(row for row in steps if row['label'] == 'member-change-retry')
    assert cleared['telemetry']['uniqueJobs'] == 13 and cleared['telemetry']['cacheHits'] == 0
    assert cleared['searches']
    for key in ('svgSha256', 'geometrySha256', 'provenance'):
        assert cleared[key] == retry[key], key
    assert cleared['rawKeys'] == cleared['provenance']['upstreamRawKeys']
    assert set(cleared['rawKeys']).issubset(retry['rawKeys'])
    assert len(retry['rawKeys']) == 38 and len(cleared['rawKeys']) == 13
    saved = steps[0]
    assert saved['telemetry']['cacheHits'] == 13 and not saved['searches']
    regenerated = lifecycle['runs'][1]['regenerated']
    assert not regenerated['searches']
    assert regenerated['provenance'] == steps[-1]['provenance']
    packaged = read_report(OUT / 'offline-package.json.gz')
    for run in [*lifecycle['runs'], *packaged['runs']]:
        assert not run['external'] and not run['pageErrors']
        assert run['liveWorkersAfterDispose'] == 0
    for run in packaged['runs']:
        assert run['linear']['geometrySha256'] == saved['geometrySha256']
        assert run['linear']['provenance'] == saved['provenance']
        assert not run['linear']['searches']
    report = {
        'base': inherited['base'], 'inheritedHead': inherited['inheritedHead'],
        'sourceSha256': source, 'sourceDigest': sha(json.dumps(source, sort_keys=True).encode()),
        'newProductionFiles': source_delta, 'unchangedPythonFiles': len(py_names),
        'wheelSha256': sha(wheel.read_bytes()), 'installedMatchingFiles': len(installed_names),
        'preservedWorktrees': preservation, 's08GeometryAndProvenanceEquality': equality,
        'sourceReplacementFreshLoadProvenanceEqual': True,
        'clearCache': {'rawSearchJobs': 13, 'matchesCurrentRetrySvgGeometryProvenance': True,
                       'retainedRawEntriesBefore': 38, 'retainedRawEntriesAfter': 13,
                       'activeRawKeysMatchProvenance': True,
                       'excludedHistoricalComparison': 'S08 observed prefix clear-cache used committed typed resources; no raw search/provenance.'},
        'packagedSourceGeometryAndProvenanceEqual': True,
        'reviewedDiffSha256': sha(git('diff', '--binary', 'HEAD', '--', 'gbdraw', 'tests')),
        'changedTrackedFiles': git('diff', '--name-only', 'HEAD').decode().splitlines(),
        'timing': 'No new accepted timing samples or warmups; command durations are operational logs.'
    }
    (OUT / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('sourceSha256', 'preservedWorktrees')}, indent=2))


if __name__ == '__main__':
    main()
