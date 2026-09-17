"""S07.6 real Collinear helper parity using the packaged wheel, offline.

Frozen/current outputs must agree completely within native and within Pyodide.
This is correctness evidence, not a Web performance measurement.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
from urllib.parse import urlparse
import zipfile

from playwright.sync_api import sync_playwright

from gbdraw.analysis import collinearity as cc
from gbdraw.analysis import protein_colinearity as pc
from tests import run_losat_cache_browser_acceptance as acceptance
from tests.prototypes import collinearity_units_frozen as frozen
from tests.test_protein_comparison_benchmark import benchmark

ROOT = Path(__file__).resolve().parents[1]


def science(value):
    assert not value.get('error'), value.get('error')
    return {k: v for k, v in value.items() if k != 'cache'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-root', type=Path, required=True)
    args = parser.parse_args()
    wheel, = (ROOT/'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as archive:
        source_hashes = {name: hashlib.sha256(archive.read(name)).hexdigest()
                         for name in archive.namelist() if name.startswith('gbdraw/') and name.endswith('.py')}
        for name, sha in source_hashes.items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha, name
    baseline_wheel, = (args.baseline_root/'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(baseline_wheel) as archive:
        for name in source_hashes:
            assert archive.read(name) == (args.baseline_root/name).read_bytes(), name
        changed = [name for name, sha in source_hashes.items()
                   if hashlib.sha256(archive.read(name)).hexdigest() != sha]
        assert changed == ['gbdraw/analysis/collinearity_units.py'], changed
    inputs = benchmark.path_browser_inputs(ROOT, pc, cc, names=('gallery-collinear', 'gallery-collinear-off'))
    # Explicit execution settings, independent of the misleading saved active draft.
    for case in inputs:
        case['parameters'].update(collinearUnitMode='auto', collinearAnchorMode='rbh',
            collinearSearchScope='adjacent', collinearMinAnchors=1, collinearMaxUnitGap=0,
            collinearMaxDiagonalDrift=0, collinearMaxConflictsInMergeGap=1,
            collinearMaxParalogLinksPerOrthogroup=2, collinearMergeOrientation='either',
            orthogroupMembershipMode='anchor_core_v1', maxHits=5)
    helper = benchmark.helpers(ROOT)['convert_losatp_blastp_pairs_to_genomic_payload']
    expected = {}
    with tempfile.TemporaryDirectory() as directory:
        pairs, raw = Path(directory)/'pairs.json', Path(directory)/'raw.tsv'
        for case in inputs:
            pairs.write_text(case['pairsText'])
            raw.write_text(case['rawText'])
            kw = dict(mode='collinear', bitscore=50, orthogroup_member_max_hits=5,
                      collinear_infer_orthogroups=case['parameters']['collinearInferOrthogroups'])
            with patch.object(cc, 'build_collinearity_unit_index', frozen.build_collinearity_unit_index):
                old = science(json.loads(helper(str(pairs), str(raw), **kw)))
            current = science(json.loads(helper(str(pairs), str(raw), **kw)))
            assert current == old, case['name']
            expected[case['name']] = current
    report = {'artifact': 'source SPA + generated local wheel', 'wheelSha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
              'sourceHashes': source_hashes, 'nativeFrozenCurrentParity': True,
              'baselineWheelSha256': hashlib.sha256(baseline_wheel.read_bytes()).hexdigest(),
              'changedPackagedPython': changed,
              'nativeOutputSha256': {name: benchmark.digest(benchmark.json_bytes(benchmark.canonical(value)))
                                     for name, value in expected.items()},
              'scope': 'real helper ON/OFF; no raw search, Generate timing, save or replay claim', 'runs': []}
    browser_expected = {}
    with acceptance._serve_repo() as base, sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            for side, width, height in (('baseline', 1280, 720), ('current', 1280, 720),
                                        ('baseline', 390, 844), ('current', 390, 844)):
                context = browser.new_context(viewport={'width': width, 'height': height})
                external, errors, paths = [], [], set()
                def route(r):
                    url = urlparse(r.request.url)
                    if url.netloc != urlparse(base).netloc:
                        external.append(r.request.url)
                        r.abort()
                    elif side == 'baseline' and url.path.endswith(wheel.name):
                        paths.add(url.path)
                        r.fulfill(status=200, content_type='application/octet-stream', body=baseline_wheel.read_bytes())
                    else:
                        paths.add(url.path)
                        r.continue_()
                context.route('**/*', route)
                page = context.new_page()
                page.set_default_timeout(120_000)
                page.on('pageerror', lambda error: errors.append(str(error)))
                page.goto(base+'/gbdraw/web/index.html', wait_until='domcontentloaded')
                page.wait_for_function('() => Boolean(window.__GBDRAW_APP__)')
                run = {'source': side, 'viewport': [width, height], 'browser': browser.version, 'cases': {}}
                for case in inputs:
                    values = page.evaluate('''async input => {
                      const api = await import('/gbdraw/web/js/services/diagram-generation.js');
                      const encoder = new TextEncoder();
                      const values = [];
                      for (let i=0; i<2; i++) {
                        const response = await api.runDiagramHelperOperation(
                          api.DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD,
                          {...input.parameters, files: [
                            {role:'pairs', bytes:encoder.encode(input.pairsText).buffer},
                            {role:'rawTsv', bytes:encoder.encode(input.rawText).buffer}]});
                        values.push(response.result);
                      }
                      return values;
                    }''', case)
                    resource = values[0]['collinearityResult']
                    assert resource['schema'] == 3 and resource['value']['type'] == 'CollinearityResult'
                    key = (width, case['name'])
                    if side == 'baseline':
                        browser_expected[key] = science(values[0])
                    for value in values:
                        assert science(value) == browser_expected[key], (side, case['name'])
                    assert not values[0]['cache']['convertedPayloadHit']
                    assert values[1]['cache']['convertedPayloadHit']
                    groups = resource['value']['fields']['orthogroups']
                    assert (groups is not None) == case['parameters']['collinearInferOrthogroups']
                    run['cases'][case['name']] = {'parameters': case['parameters'], 'fixture': case['inventory'],
                        'coldRepeatFrozenParity': True, 'nativeExactParity': science(values[0]) == expected[case['name']], 'schema': 3, 'fullOutputSha256': benchmark.digest(
                            benchmark.json_bytes(benchmark.canonical(science(values[0]))))}
                page.evaluate("async () => (await import('/gbdraw/web/js/services/diagram-generation.js')).disposeDiagramGenerationWorker()")
                assert not external and not errors, (external, errors)
                run.update(externalRequests=external, pageErrors=errors, localPaths=sorted(paths))
                report['runs'].append(run)
                context.close()
        finally:
            browser.close()
    data = (json.dumps(report, indent=2)+'\n').encode()
    args.output.write_bytes(gzip.compress(data, mtime=0))
    print(json.dumps({'runs': len(report['runs']), 'fullOutputParity': True}))


if __name__ == '__main__':
    main()
