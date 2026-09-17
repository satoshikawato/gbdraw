"""Bind S09 review to inherited bytes, public witnesses and new browser evidence."""
from __future__ import annotations

import ast
import base64
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


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return digest(path.read_bytes())


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args])


def read(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def category(name):
    if name.startswith('gbdraw/'):
        return 'production'
    if name.startswith(('tests/', 'tools/')):
        return 'tests-and-tools'
    if '/results/data/' in name:
        return 'evidence'
    return 'documentation'


def categories(names):
    result = {kind: [] for kind in ('production', 'tests-and-tools', 'documentation', 'evidence')}
    for name in sorted(set(names)):
        result[category(name)].append(name)
    return result


def main():
    inherited = read(OUT / 'inheritance.json')
    source = Path(inherited['source'])
    head, base = inherited['inheritedHead'], inherited['base']
    assert git(ROOT, 'rev-parse', 'HEAD').decode().strip() == head
    assert git(ROOT, 'rev-parse', 'origin/dev').decode().strip() == base
    assert git(ROOT, 'branch', '--show-current').decode().strip() == inherited['branch']
    assert not git(ROOT, 'for-each-ref', '--format=%(upstream)', 'refs/heads/' + inherited['branch']).strip()
    hashes = {name: sha(ROOT / name) for name in inherited['sourceSha256']}
    changed = [name for name in hashes if hashes[name] != inherited['sourceSha256'][name]]
    assert changed == ['gbdraw/web/js/services/config.js'], changed
    for name, expected in inherited['rw04Inventory'].items():
        assert sha(source / name) == expected['sha256'], name
        if name not in {str((PLAN / 'MASTER_PLAN.md').relative_to(ROOT)),
                        'tests/run_protein_local_browser_acceptance.py'}:
            assert sha(ROOT / name) == expected['sha256'], name
    previous = read(DATA / 's08-followup/audit.json')
    unchanged_python = [name for name in hashes if name.endswith('.py')]
    assert len(unchanged_python) == 208
    for name in unchanged_python:
        assert hashes[name] == previous['sourceSha256'][name], name

    # Old result documents/audits are immutable snapshots; do not execute them.
    historical = git(ROOT, 'ls-files', str(PLAN.relative_to(ROOT))).decode().splitlines()
    historical_changes = git(ROOT, 'diff', '--name-only', head, '--', str(PLAN.relative_to(ROOT))).decode().splitlines()
    assert set(historical_changes) <= set(inherited['rw04Inventory']), historical_changes

    authorities = ['AGENTS.md', 'CLAUDE.md', 'gbdraw/web/CLAUDE.md',
        'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
        'docs/internal/PRODUCT_IMPACT_RATCHET.md',
        'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md',
        'tools/web-product-impact-map.json', 'tools/web-product-decisions.json',
        'tools/web-architecture-rules.json', 'tools/check-web-change-budget.mjs']
    for name in authorities:
        assert (ROOT / name).read_bytes() == git(ROOT, 'show', f'{base}:{name}'), name

    # Confirm publication from Git objects, independently of old report prose.
    publication = read(DATA / 's02-publication.json.gz')
    first_parent = set(git(ROOT, 'rev-list', '--first-parent', publication['main']).decode().splitlines())
    witnesses = []
    for spec, observation in publication['historicalPositiveArtifacts'].items():
        revision, name = spec.split(':', 1)
        assert revision == '0.13.0' or revision in first_parent
        raw = git(ROOT, 'show', spec)
        assert digest(raw) == observation['encodedSha256'], spec
        document = json.loads(gzip.decompress(raw) if name.endswith('.gz') else raw)
        typed = {}
        for resource_id, expected in observation['typedResources'].items():
            resource = document['resources'][resource_id]
            payload = base64.b64decode(resource['data'])
            assert digest(payload) == expected['sha256']
            assert json.loads(payload)['schema'] == expected['schema']
            typed[resource_id] = expected
        witnesses.append({'gitObject': spec, 'sha256': digest(raw), 'typedResources': typed})
    api_source = git(ROOT, 'show', '0.13.0:gbdraw/api/__init__.py')
    options_source = git(ROOT, 'show', '0.13.0:gbdraw/api/options.py')
    analysis_source = git(ROOT, 'show', '0.13.0:gbdraw/analysis/protein_colinearity.py')
    assert b'OrthologPath' in api_source and b'OrthologEdge' in api_source
    assert b'orthogroups: OrthogroupResult | None' in options_source
    assert b'class OrthogroupResult:' in analysis_source

    core_names = {'OrthogroupResult', 'OrthogroupGraphResult', 'compact_ortholog_paths',
                  'materialize_ortholog_paths', '_build_ortholog_path_indexes'}
    pc = 'gbdraw/analysis/protein_colinearity.py'
    def definitions(text):
        return {node.name: ast.dump(node, include_attributes=False)
                for node in ast.parse(text).body if getattr(node, 'name', None) in core_names}
    assert definitions(git(ROOT, 'show', f'ec16110c:{pc}')) == definitions((ROOT / pc).read_text())
    assert len(definitions((ROOT / pc).read_text())) == len(core_names)
    unchanged_s06 = ['gbdraw/analysis/ortholog_paths.py', 'gbdraw/session_request_codec.py',
        'gbdraw/api/session_compat.py', 'gbdraw/api/request_render.py',
        'gbdraw/web_support/orthogroup_metadata.py', 'gbdraw/web/js/services/session-request.js']
    for name in unchanged_s06:
        assert (ROOT / name).read_bytes() == git(ROOT, 'show', f'ec16110c:{name}'), name

    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as archive:
        for name in unchanged_python:
            assert archive.read(name) == (ROOT / name).read_bytes(), name
    installed = ROOT / '.venv/s09-package'
    support = run_path(str(ROOT / 'gbdraw/_build_support.py'))
    packaged = {str(path.relative_to(ROOT))
                for pattern in support['get_package_data_patterns'](include_browser_wheel=True)
                for path in (ROOT / 'gbdraw').glob(pattern) if path.is_file()}
    tracked = git(ROOT, 'ls-files', 'gbdraw').decode().splitlines()
    installed_names = [n for n in tracked if n.endswith('.py') or n in packaged]
    assert set(installed_names) == {n for n in tracked if (installed / n).is_file()}
    for name in installed_names:
        assert sha(ROOT / name) == sha(installed / name), name
    assert sha(wheel) == sha(installed / 'gbdraw/web' / wheel.name)

    parity = []
    for filename in ['browser-source.json.gz', 'browser-package.json.gz']:
        report = read(OUT / filename)
        assert report['runnerSha256'] == sha(ROOT / 'tests/run_protein_local_browser_acceptance.py')
        assert report['wheelSha256'] == sha(wheel)
        assert len(report['runs']) == 4
        prior = read(DATA / 's08-alias-validation' / filename)
        for index in range(0, 4, 2):
            before, after = report['runs'][index:index + 2]
            assert before['viewport'] == after['viewport'] == prior['runs'][index // 2]['viewport']
            assert before['fixtureSha256'] == prior['runs'][index // 2]['fixtureSha256']
            for run in (before, after):
                assert not run['external'] and not run['pageErrors']
                assert run['workersAfterPreviewLoad'] == run['liveWorkersAfterDispose'] == 0
            assert len(before['writerChecks']) == 4
            assert all(row['preserved'] and row['downloads'] == 0 and row['status'] is None
                       and row['error'] == 'Save Session requires a valid protein identity manifest.'
                       for row in before['writerChecks'])
            assert before['savedRawCount'] == 13
            for step in [before['savedRaw'], after['regenerated']]:
                assert step['telemetry']['cacheHits'] == 13 and not step['searches']
                for field in ['svgSha256', 'geometrySha256', 'provenance']:
                    assert step[field] == prior['runs'][index // 2]['savedRaw'][field], field
            parity.append({'artifact': filename, 'viewport': before['viewport'],
                           'savedRawHits': 13, 'searches': 0, 'invalidSaveCases': 4,
                           'svgGeometryProvenanceEqualToRW04': True})

    committed = categories(git(ROOT, 'diff', '--name-only', base, head).decode().splitlines())
    dirty = git(ROOT, 'diff', '--name-only', head).decode().splitlines()
    dirty += git(ROOT, 'ls-files', '--others', '--exclude-standard').decode().splitlines()
    new = categories(n for n in dirty if not (source / n).exists() or sha(ROOT / n) != sha(source / n))
    assert new['production'] == changed
    for kind in ('production', 'tests-and-tools', 'documentation'):
        parts = []
        for name in new[kind]:
            before = (source / name).read_text() if (source / name).exists() else ''
            parts.extend(difflib.unified_diff(before.splitlines(True), (ROOT / name).read_text().splitlines(True),
                         fromfile='inherited-final/' + name, tofile='s09/' + name))
        (OUT / f'new-{kind}.patch').write_text(''.join(parts))
    for location, expected in inherited['preserved'].items():
        actual = {'head': git(location, 'rev-parse', 'HEAD').decode().strip(),
                  'status': git(location, 'status', '--porcelain').decode(),
                  'diffSha256': digest(git(location, 'diff', 'HEAD', '--binary'))}
        assert actual == expected, location

    result = {'base': base, 'head': head, 'branch': inherited['branch'], 'upstream': None,
        'sourceIsHeadPlusUncommittedRW04AndS09': True,
        'sourceSha256': hashes, 'sourceAggregateSha256': digest(json.dumps(hashes, sort_keys=True).encode()),
        's09ChangedSourceFiles': changed, 'unchangedPythonFiles': len(unchanged_python),
        'testsSha256': {n: sha(ROOT / n) for n in ['tests/run_protein_local_browser_acceptance.py',
            'tests/web/losat-cache.test.mjs', 'tests/web/session-export-validation.test.mjs']},
        'authoritiesUnchangedFromBase': {n: sha(ROOT / n) for n in authorities},
        'baseSessionContractSha256': digest(git(ROOT, 'show', f'{base}:docs/SESSION_COMPATIBILITY.md')),
        'committedDiff': committed, 'rw04UncommittedDiff': categories(inherited['rw04Inventory']),
        's09Diff': new, 'historicalTrackedFilesPreserved': len(historical),
        'preservedWorktrees': list(inherited['preserved']), 'publicWitnesses': witnesses,
        'releaseApiWitnesses': {'0.13.0:gbdraw/api/__init__.py': digest(api_source),
            '0.13.0:gbdraw/api/options.py': digest(options_source),
            '0.13.0:gbdraw/analysis/protein_colinearity.py': digest(analysis_source)},
        's06UnchangedPathDefinitions': sorted(core_names), 's06UnchangedBoundaryFiles': unchanged_s06,
        'wheelSha256': sha(wheel), 'installedMatchingFiles': len(installed_names),
        'browserParity': parity, 'timingSamples': 0, 'warmups': 0,
        'humanS06ExceptionAndS08WriterReview': 'PENDING; no approved exact-head commit exists'}
    (OUT / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in {
        'sourceSha256', 'committedDiff', 'rw04UncommittedDiff', 's09Diff'}}, indent=2))


if __name__ == '__main__':
    main()
