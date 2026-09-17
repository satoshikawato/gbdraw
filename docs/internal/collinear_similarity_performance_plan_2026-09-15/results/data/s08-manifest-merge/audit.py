"""Read-only source/parity audit; never refresh historical evidence."""
from __future__ import annotations

import difflib
import gzip
import hashlib
import json
from pathlib import Path
from runpy import run_path
import subprocess
import zipfile

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[5]
DATA = OUT.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_report(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args])


def main():
    inheritance = json.loads((OUT / 'inheritance.json').read_text())
    source = Path(inheritance['source'])
    assert git(ROOT, 'rev-parse', 'HEAD').decode().strip() == inheritance['head']
    assert git(ROOT, 'rev-parse', 'origin/dev').decode().strip() == inheritance['base']
    branch = git(ROOT, 'branch', '--show-current').decode().strip()
    assert branch == 'perf/s08-manifest-merge-20260917'
    assert not git(ROOT, 'for-each-ref', '--format=%(upstream)', f'refs/heads/{branch}').strip()
    for name, expected in inheritance['files'].items():
        assert sha(source / name) == expected, name
    assert hashlib.sha256(git(source, 'diff', '--binary', 'HEAD')).hexdigest() == inheritance['inheritedTrackedDiffSha256']
    # Historical evidence and the RW-01/RW-02 result remain byte-identical.
    historical = [name for name in inheritance['files']
                  if '/s08-followup/' in name or name.endswith('S08_FOLLOWUP_RAW_VALIDATION.md')]
    for name in historical:
        assert sha(ROOT / name) == inheritance['files'][name], name
    preserved = json.loads((OUT / 'preserved-start.json').read_text())
    for name, before in preserved.items():
        path = Path(name)
        after = {'head': git(path, 'rev-parse', 'HEAD').decode().strip(),
                 'status': git(path, 'status', '--porcelain').decode(),
                 'trackedDiffSha256': hashlib.sha256(git(path, 'diff', 'HEAD', '--binary')).hexdigest()}
        assert before == after, name
    old = json.loads((DATA / 's08-followup/audit.json').read_text())
    hashes = {name: sha(ROOT / name) for name in old['sourceSha256']}
    changed = [name for name, value in hashes.items() if value != old['sourceSha256'][name]]
    assert changed == ['gbdraw/web/js/app/losat-cache.js', 'gbdraw/web/js/app/run-analysis.js']
    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as zf:
        python_files = [name for name in zf.namelist() if name.startswith('gbdraw/') and name.endswith('.py')]
        for name in python_files:
            assert zf.read(name) == (ROOT / name).read_bytes(), name
            assert hashes[name] == old['sourceSha256'][name], name
    installed = ROOT / '.venv/manifest-merge-package'
    tracked = git(ROOT, 'ls-files', 'gbdraw').decode().splitlines()
    support = run_path(str(ROOT / 'gbdraw/_build_support.py'))
    packaged_data = {str(path.relative_to(ROOT))
                     for pattern in support['get_package_data_patterns'](include_browser_wheel=True)
                     for path in (ROOT / 'gbdraw').glob(pattern) if path.is_file()}
    installed_files = [name for name in tracked if name.endswith('.py') or name in packaged_data]
    assert set(installed_files) == {name for name in tracked if (installed / name).is_file()}
    for name in installed_files:
        assert sha(ROOT / name) == sha(installed / name), name
    assert sha(wheel) == sha(installed / 'gbdraw/web' / wheel.name)
    baseline = read_report(DATA / 's08-followup/lifecycle.json.gz')['runs'][0]['steps'][0]
    rows = []
    for filename in ['browser-source.json.gz', 'browser-package.json.gz']:
        report = read_report(OUT / filename)
        assert report['runnerSha256'] == sha(ROOT / 'tests/run_protein_local_browser_acceptance.py')
        assert report['wheelSha256'] == sha(wheel)
        assert len(report['runs']) == 2
        for run in report['runs']:
            assert run['fixtureSha256'] == sha(ROOT / 'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz')
            assert run['workersAfterPreviewLoad'] == run['liveWorkersAfterDispose'] == 0
            assert not run['external'] and not run['pageErrors']
            for key in ['savedRaw', 'savedRawRetry']:
                value = run[key]
                for field in ['svgSha256', 'geometrySha256', 'provenance']:
                    assert value[field] == baseline[field], (filename, key, field)
                assert not value['searches'] and value['telemetry']['cacheHits'] == 13
            failed = run['invalidManifest']
            assert failed['preserved'] and not failed['searches']
            assert failed['result'] == {'status': 'error'}
            rows.append({'artifact': filename, 'viewport': run['viewport'],
                         'svgSha256': run['savedRaw']['svgSha256'],
                         'geometrySha256': run['savedRaw']['geometrySha256'],
                         'provenanceEqualToRW01RW02': True, 'savedRawHits': 13,
                         'searches': 0, 'failurePreservedResultHistoryCachesManifest': True,
                         'error': failed['error'], 'liveWorkersAfterDispose': 0})
    classifications = {'production': [], 'tests': [], 'documentation': []}
    current_names = git(ROOT, 'diff', 'HEAD', '--name-only').decode().splitlines()
    current_names.append('docs/internal/collinear_similarity_performance_plan_2026-09-15/results/S08_FOLLOWUP_MANIFEST_MERGE.md')
    for category in classifications:
        patches = []
        for name in current_names:
            target_category = 'production' if name.startswith('gbdraw/') else 'tests' if name.startswith('tests/') else 'documentation'
            if target_category != category:
                continue
            before = (source / name).read_text() if (source / name).exists() else ''
            after = (ROOT / name).read_text()
            if before == after:
                continue
            classifications[category].append(name)
            patches.extend(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                           fromfile='inherited/' + name, tofile='rw03/' + name))
        (OUT / f'new-{category}.patch').write_text(''.join(patches))
    result = {'base': inheritance['base'], 'inheritedHead': inheritance['head'], 'branch': branch,
              'upstream': None, 'sourceSha256': hashes,
              'newProductionHashes': {name: hashes[name] for name in changed},
              'unchangedSourceFiles': len(hashes) - len(changed), 'unchangedPythonFiles': len(python_files),
              'wheelSha256': sha(wheel), 'installedMatchingFiles': len(installed_files),
              'historicalFilesPreserved': len(historical), 'preservedWorktrees': list(preserved),
              'browserParity': rows, 'newDiffClassification': classifications,
              'testsSha256': {name: sha(ROOT / name) for name in classifications['tests']},
              'operationCounts': {'R': 3, 'beforeInput': 6, 'afterInput': 3, 'beforeMerged': 1, 'afterMerged': 1,
                                  'method': 'instrument real validator; execute actual run-analysis merge block'},
              'timing': 'No timing samples or warmups. Verification/build duration is not performance evidence.'}
    (OUT / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'sourceSha256'}, indent=2))


if __name__ == '__main__':
    main()
