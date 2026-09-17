"""Verify RW-04 source, inherited evidence, package parity and worktree isolation."""
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
PLAN = DATA.parent.parent


def digest(value):
    return hashlib.sha256(value).hexdigest()


def sha(path):
    return digest(path.read_bytes())


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args])


def report(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def category(name):
    if name.startswith('gbdraw/'):
        return 'production'
    if name.startswith(('tests/', 'tools/')):
        return 'tests-and-tools'
    if '/results/data/' in name:
        return 'generated-evidence'
    return 'documentation'


def main():
    inheritance = json.loads((OUT / 'inheritance.json').read_text())
    head, base = inheritance['inheritedCommit'], inheritance['base']
    assert git(ROOT, 'rev-parse', 'HEAD').decode().strip() == head
    assert git(ROOT, 'rev-parse', 'origin/dev').decode().strip() == base
    assert git(ROOT, 'branch', '--show-current').decode().strip() == inheritance['branch']
    assert not git(ROOT, 'for-each-ref', '--format=%(upstream)',
                   'refs/heads/' + inheritance['branch']).strip()
    preserved = json.loads((OUT / 'preserved-start.json').read_text())
    for name, before in preserved.items():
        after = {'head': git(name, 'rev-parse', 'HEAD').decode().strip(),
                 'status': git(name, 'status', '--porcelain').decode(),
                 'trackedDiffSha256': digest(git(name, 'diff', 'HEAD', '--binary'))}
        assert before == after, name

    previous_path = DATA / 's08-manifest-merge/audit.json'
    assert sha(previous_path) == inheritance['priorAuditSha256']
    previous = json.loads(previous_path.read_text())
    hashes = {name: sha(ROOT / name) for name in previous['sourceSha256']}
    changed = [name for name in hashes if hashes[name] != previous['sourceSha256'][name]]
    assert changed == ['gbdraw/web/js/app/losat-cache.js']
    assert git(ROOT, 'diff', '--name-only', head, '--', 'gbdraw').decode().splitlines() == changed
    historical = git(ROOT, 'ls-files', str(PLAN.relative_to(ROOT))).decode().splitlines()
    historical = [name for name in historical
                  if not name.endswith(('/MASTER_PLAN.md', '/S08_REDUNDANT_WORK.md'))]
    assert not git(ROOT, 'diff', head, '--', *historical).strip()

    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as archive:
        python_files = [name for name in archive.namelist()
                        if name.startswith('gbdraw/') and name.endswith('.py')]
        for name in python_files:
            assert archive.read(name) == (ROOT / name).read_bytes(), name
            assert hashes[name] == previous['sourceSha256'][name], name
    installed = ROOT / '.venv/alias-validation-package'
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

    parity = []
    for filename in ('browser-source.json.gz', 'browser-package.json.gz'):
        current = report(OUT / filename)
        old = report(DATA / 's08-manifest-merge' / filename)
        assert current['runnerSha256'] == sha(ROOT / 'tests/run_protein_local_browser_acceptance.py')
        assert current['wheelSha256'] == sha(wheel)
        assert len(current['runs']) == len(old['runs']) == 2
        for run, prior in zip(current['runs'], old['runs'], strict=True):
            assert run['viewport'] == prior['viewport']
            assert run['fixtureSha256'] == prior['fixtureSha256'] == sha(
                ROOT / 'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz')
            assert run['workersAfterPreviewLoad'] == run['liveWorkersAfterDispose'] == 0
            assert not run['external'] and not run['pageErrors']
            assert run['manifestFault'] == 'blank-alias'
            for step in ('savedRaw', 'savedRawRetry'):
                for field in ('svgSha256', 'geometrySha256', 'provenance'):
                    assert run[step][field] == prior[step][field], (filename, step, field)
                assert run[step]['telemetry']['cacheHits'] == 13 and not run[step]['searches']
            assert run['resolvedRepeat']['svgSha256'] == run['savedRaw']['svgSha256']
            assert not run['resolvedRepeat']['searches']
            failed = run['invalidManifest']
            assert failed['preserved'] and failed['result'] == {'status': 'error'}
            assert failed['error'] == prior['invalidManifest']['error']
            assert not failed['searches'] and failed['workerCount'] == failed['liveWorkers'] == 1
            parity.append({'artifact': filename, 'viewport': run['viewport'],
                           'svgSha256': run['savedRaw']['svgSha256'],
                           'geometrySha256': run['savedRaw']['geometrySha256'],
                           'provenanceEqualToRW03': True, 'savedRawHits': 13, 'searches': 0,
                           'blankAliasFailurePreservedResultHistoryCachesManifest': True,
                           'error': failed['error'], 'liveWorkersAfterDispose': 0})

    inherited = {key: [] for key in ('production', 'tests-and-tools', 'documentation', 'generated-evidence')}
    for name in git(ROOT, 'diff', '--name-only', base, head).decode().splitlines():
        inherited[category(name)].append(name)
    current = {key: [] for key in inherited}
    names = git(ROOT, 'diff', '--name-only', head).decode().splitlines()
    names += git(ROOT, 'ls-files', '--others', '--exclude-standard').decode().splitlines()
    for name in sorted(set(names)):
        current[category(name)].append(name)
    for kind in ('production', 'tests-and-tools', 'documentation'):
        patches = []
        for name in current[kind]:
            tracked_before = git(ROOT, 'ls-tree', '--name-only', head, '--', name).strip()
            before = git(ROOT, 'show', f'{head}:{name}').decode() if tracked_before else ''
            after = (ROOT / name).read_text()
            patches.extend(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                           fromfile='inherited/' + name, tofile='rw04/' + name))
        (OUT / f'new-{kind}.patch').write_text(''.join(patches))
    assert '2 valid feature aliases, 4 NFC/trim operations' in (OUT / 'focused-red.log').read_text()
    assert '2 valid feature aliases, 2 NFC/trim operations' in (OUT / 'focused-green.log').read_text()
    result = {'base': base, 'inheritedCommit': head, 'branch': inheritance['branch'], 'upstream': None,
              'sourceSha256': hashes, 'sourceAggregateSha256': digest(json.dumps(hashes, sort_keys=True).encode()),
              'changedSourceFiles': changed, 'unchangedSourceFiles': len(hashes) - len(changed),
              'unchangedPythonFiles': len(python_files), 'wheelSha256': sha(wheel),
              'installedMatchingFiles': len(installed_files), 'historicalFilesPreserved': len(historical),
              'preservedWorktrees': list(preserved), 'browserParity': parity,
              'inheritedDiffClassification': inherited, 'newDiffClassification': current,
              'testsSha256': {name: sha(ROOT / name) for name in current['tests-and-tools']},
              'operationCounts': {'scope': 'one manifest validation; two valid feature aliases',
                                  'beforeNfcTrim': 4, 'afterNfcTrim': 2,
                                  'manifestValidationCountChanged': False},
              'timing': 'Zero timing samples and warmups. Verification/build duration is not performance evidence.'}
    (OUT / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {'sourceSha256', 'inheritedDiffClassification', 'newDiffClassification'}}, indent=2))


if __name__ == '__main__':
    main()
