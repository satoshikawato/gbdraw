"""S07.7 correctness only: frozen/current helpers and final-wheel session flows."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse
import zipfile

from playwright.sync_api import sync_playwright

from gbdraw.analysis import collinearity as cc, protein_colinearity as pc
from tests import run_losat_cache_browser_acceptance as acceptance
from tests.test_protein_comparison_benchmark import benchmark

ROOT = Path(__file__).resolve().parents[1]
OBSERVE = """window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
  window.__GBDRAW_LOSAT_EXECUTOR_CALLS__++;
  throw new Error('S07.7 requires reuse of the saved raw evidence');
};
window.__s077 = {requests: [], helpers: []};
const originalPost = Worker.prototype.postMessage;
Worker.prototype.postMessage = function(message, ...args) {
  if (message?.type === 'run') window.__s077.requests.push(structuredClone(message.payload.request));
  if (message?.type === 'helper' && message.operation === 'convertLosatpPairsToGenomicPayload') {
    const {files, ...parameters} = message.payload;
    window.__s077.helpers.push(structuredClone(parameters));
  }
  return originalPost.call(this, message, ...args);
};"""


def digest(data):
    return hashlib.sha256(data).hexdigest()


def science(value):
    assert not value.get('error'), value.get('error')
    return benchmark.json_bytes(benchmark.canonical({k: v for k, v in value.items() if k != 'cache'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='Reuse completed checks only with identical wheel/source hashes')
    args = parser.parse_args()
    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    baseline, = (args.baseline_root / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as current, zipfile.ZipFile(baseline) as old:
        names = [n for n in current.namelist() if n.startswith('gbdraw/') and n.endswith('.py')]
        source_hashes = {n: digest(current.read(n)) for n in names}
        for name in names:
            assert current.read(name) == (ROOT / name).read_bytes(), name
            assert old.read(name) == (args.baseline_root / name).read_bytes(), name
        assert [n for n in names if current.read(n) != old.read(n)] == ['gbdraw/analysis/protein_colinearity.py']
    artifacts = ROOT / '.venv/s07-7-browser-artifacts'
    artifacts.mkdir(parents=True, exist_ok=True)
    inputs = benchmark.path_browser_inputs(ROOT, pc, cc,
        names=('gallery-collinear', 'gallery-collinear-off', 'gallery-orthogroup'))
    unbounded = dict(inputs[-1], name='gallery-orthogroup-unbounded',
                     parameters=dict(inputs[-1]['parameters'], orthogroupMemberMaxHits=None))
    inputs.append(unbounded)
    report = {'wheelSha256': digest(wheel.read_bytes()), 'baselineWheelSha256': digest(baseline.read_bytes()),
              'sourceHashes': source_hashes, 'helpers': [], 'flows': [],
              'performance': 'Not measured; no timing/profile/memory collection'}
    if args.resume:
        previous = json.loads(gzip.decompress(args.output.read_bytes()))
        for key in ('wheelSha256', 'baselineWheelSha256', 'sourceHashes'):
            assert previous[key] == report[key], key
        cases_by_name = {case['name']: case for case in inputs}
        for run in previous['helpers']:
            case = cases_by_name[run['case']]
            assert run['parameters'] == case['parameters']
            assert run['inputSha256'] == digest((case['pairsText'] + case['rawText']).encode())
        for run in previous['flows']:
            mode = run['case'].split('-')[1]
            fixture = ROOT / f'gbdraw/web/gallery/sessions/hepatoplasmataceae_{mode}.gbdraw-session.json.gz'
            assert run['fixtureSha256'] == digest(fixture.read_bytes())
            assert run['sessionSha256'] == digest((ROOT / run['sessionPath']).read_bytes())
        report = previous
    expected = {}
    with acceptance._serve_repo() as base, sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        if args.resume:
            assert report['browser'] == browser.version
        report['browser'] = browser.version
        def context_for(side, viewport):
            context = browser.new_context(accept_downloads=True, viewport={'width': viewport[0], 'height': viewport[1]})
            external, errors, paths = [], [], set()
            def route(r):
                url = urlparse(r.request.url)
                if url.netloc != urlparse(base).netloc:
                    external.append(r.request.url)
                    r.abort()
                elif side == 'baseline' and url.path.endswith(wheel.name):
                    paths.add(url.path)
                    r.fulfill(status=200, content_type='application/octet-stream', body=baseline.read_bytes())
                else:
                    paths.add(url.path)
                    r.continue_()
            context.route('**/*', route)
            context.add_init_script(OBSERVE)
            page = context.new_page()
            page.set_default_timeout(120_000)
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(base + '/gbdraw/web/index.html', wait_until='domcontentloaded')
            page.wait_for_function('() => Boolean(window.__GBDRAW_APP__)')
            return context, page, external, errors, paths
        try:
            for side in ('baseline', 'current'):
                if len([run for run in report['helpers'] if run['side'] == side]) == len(inputs):
                    continue
                context, page, external, errors, paths = context_for(side, (1280, 720))
                for case in inputs:
                    values = page.evaluate('''async input => {
                      const api = await import('/gbdraw/web/js/services/diagram-generation.js');
                      const encoder = new TextEncoder();
                      const values = [];
                      for (let i = 0; i < 2; i++) {
                        const response = await api.runDiagramHelperOperation(
                          api.DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD,
                          {...input.parameters, files: [
                            {role: 'pairs', bytes: encoder.encode(input.pairsText).buffer},
                            {role: 'rawTsv', bytes: encoder.encode(input.rawText).buffer}]});
                        values.push(response.result);
                      }
                      return values;
                    }''', case)
                    if side == 'baseline':
                        expected[case['name']] = science(values[0])
                    assert all(science(value) == expected[case['name']] for value in values), case['name']
                    assert not values[0]['cache']['convertedPayloadHit'] and values[1]['cache']['convertedPayloadHit']
                    report['helpers'].append({'side': side, 'case': case['name'], 'parameters': case['parameters'],
                        'inputSha256': digest((case['pairsText'] + case['rawText']).encode()),
                        'outputSha256': digest(expected[case['name']]), 'coldRepeatParity': True})
                    print('helper', side, case['name'], flush=True)
                assert not external and not errors, (external, errors)
                context.close()
            for viewport in ((1280, 720), (390, 844)):
                for mode, infer, limit in [('collinear', True, 5), ('collinear', False, 5),
                                           ('orthogroup', True, 5), ('orthogroup', True, None)]:
                    label = f'{viewport[0]}-{mode}-{infer}-{limit}'
                    if any(flow['case'] == label for flow in report['flows']):
                        continue
                    fixture = ROOT / f'gbdraw/web/gallery/sessions/hepatoplasmataceae_{mode}.gbdraw-session.json.gz'
                    checks = acceptance.AcceptanceChecks()
                    context, page, external, errors, paths = context_for('current', viewport)
                    acceptance._import_session(page, fixture, checks)
                    async_drop_derived = '''async () => {
                      const {state} = await import('/gbdraw/web/js/state.js');
                      const rawEntries = state.losatCache.value.size;
                      const resolved = state.files.linearCanonicalComparisons.length;
                      await window.__GBDRAW_APP__.setLinearComparisonGlobalAction('losat');
                      const removed = state.losatDerivedCache.value.size;
                      state.losatDerivedCache.value.clear();
                      if (state.losatCache.value.size !== rawEntries) throw new Error('Raw cache changed');
                      if (state.files.linearCanonicalComparisons.length) throw new Error('Resolved artifacts remain');
                      return {removed, resolved, rawEntries, action: 'Run LOSAT'};
                    }'''
                    cache_preparation = page.evaluate(async_drop_derived)
                    page.evaluate('''async options => {
                      const app = window.__GBDRAW_APP__;
                      app.losat.blastp.mode = options.mode;
                      await Vue.nextTick();
                      Object.assign(app.losat.blastp, {candidateLimit: null,
                        orthogroupMemberMaxHits: options.limit,
                        collinearInferOrthogroups: options.infer});
                      await Vue.nextTick();
                    }''', {'mode': mode, 'infer': infer, 'limit': limit})
                    run = {'case': label, 'viewport': viewport, 'fixtureSha256': digest(fixture.read_bytes()),
                           'cachePreparation': cache_preparation, 'generations': []}
                    def generate(page):
                        page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length > 0')
                        observed = acceptance._generate(page)
                        assert observed['result'] == {'status': 'ok'}, observed
                        assert observed['executorCalls'] == 0
                        acceptance._settle_app_render(page)
                        exported = acceptance._download_svg(page, checks)
                        acceptance._assert_svg_geometry_parity(page, exported, checks)
                        observed_requests = page.evaluate('() => window.__s077')
                        request = observed_requests['requests'][-1]
                        settings = next(c['settings'] for c in request['comparisons'] if c['kind'] == 'generatedProteinComparison')
                        assert settings['orthogroupMemberMaxHits'] == limit
                        assert settings['proteinBlastpCandidateLimit'] is None
                        if mode == 'collinear':
                            assert settings['collinearInferOrthogroups'] == infer
                            assert settings['collinearitySearchScope'] == 'adjacent'
                        helpers = observed_requests['helpers']
                        for parameters in helpers:
                            assert parameters['mode'] == mode and parameters['orthogroupMemberMaxHits'] == limit
                            if mode == 'collinear':
                                assert parameters['collinearInferOrthogroups'] == infer
                        geometry = acceptance._inspect_layout(page)['geometrySignature']
                        run['generations'].append({'settings': settings, 'helpers': helpers,
                            'result': observed['result'], 'executorCalls': 0, 'geometrySha256': digest(geometry.encode()),
                            'svgSha256': digest(exported.read_bytes())})
                        return geometry, exported.read_bytes()
                    first, svg = generate(page)
                    assert run['generations'][0]['helpers'], 'Initial Generate did not exercise real Python helper'
                    assert generate(page)[0] == first
                    (artifacts / f'{label}.svg').write_bytes(svg)
                    saved = acceptance._save_session(page, checks)
                    saved_path = artifacts / f'{label}.gbdraw-session.json.gz'
                    saved_path.write_bytes(saved.read_bytes())
                    run['sessionSha256'] = digest(saved.read_bytes())
                    run['sessionPath'] = str(saved_path.relative_to(ROOT))
                    assert not external and not errors, (external, errors)
                    context.close()
                    # A new browser context has no Worker, module, storage or HTTP cache from Generate.
                    context, page, external, errors, paths = context_for('current', viewport)
                    acceptance._import_session(page, saved_path, checks)
                    run['freshLoadCachePreparation'] = page.evaluate(async_drop_derived)
                    assert page.evaluate('() => window.__GBDRAW_APP__.losat.blastp.mode') == mode
                    assert generate(page)[0] == first
                    assert not external and not errors, (external, errors)
                    run.update(externalRequests=external, pageErrors=errors, localPaths=sorted(paths), assertions=checks.count)
                    report['flows'].append(run)
                    args.output.write_bytes(gzip.compress((json.dumps(report, indent=2) + '\n').encode(), mtime=0))
                    print('flow', label, 'passed', flush=True)
                    context.close()
        finally:
            browser.close()
    args.output.write_bytes(gzip.compress((json.dumps(report, indent=2) + '\n').encode(), mtime=0))
    print(json.dumps({'helperRuns': len(report['helpers']), 'flows': len(report['flows']), 'exactParity': True}))


if __name__ == '__main__':
    main()
