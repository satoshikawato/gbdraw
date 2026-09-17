"""Frozen helper parity and integrated real-browser acceptance/measurement."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import statistics
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
  throw new Error('S07.8 requires reuse of the saved raw evidence');
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
    parser.add_argument('--baseline-root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='Reuse completed checks only with identical wheel/source hashes')
    parser.add_argument('--integration', choices=('timing', 'lifecycle', 'replacement', 'offline', 'manifest-merge'))
    parser.add_argument('--cases', nargs='+', choices=tuple(INTEGRATION_CASES), default=list(INTEGRATION_CASES))
    parser.add_argument('--samples', type=int, choices=(1, 2, 3), default=3)
    parser.add_argument('--server-root', type=Path, default=ROOT,
                        help='Source checkout or installed distribution site-packages')
    args = parser.parse_args()
    if args.integration:
        return run_integration(args)
    if args.baseline_root is None:
        parser.error('--baseline-root is required for differential acceptance')
    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    baseline, = (args.baseline_root / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as current, zipfile.ZipFile(baseline) as old:
        names = [n for n in current.namelist() if n.startswith('gbdraw/') and n.endswith('.py')]
        source_hashes = {n: digest(current.read(n)) for n in names}
        for name in names:
            assert current.read(name) == (ROOT / name).read_bytes(), name
            assert old.read(name) == (args.baseline_root / name).read_bytes(), name
        assert {n for n in names if current.read(n) != old.read(n)} == {
            'gbdraw/analysis/protein_colinearity.py', 'gbdraw/analysis/collinearity_units.py'}
    artifacts = ROOT / '.venv' / (args.output.name.split('.')[0] + '-artifacts')
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


INTEGRATION_CASES = {
    'hep-on': ('hepatoplasmataceae_collinear', 'collinear', True, None),
    'hep-off': ('hepatoplasmataceae_collinear', 'collinear', False, None),
    'hep-similarity': ('hepatoplasmataceae_orthogroup', 'orthogroup', True, None),
    'vibrio-on': ('vibrio-harveyi-group-collinear', 'collinear', True, 5),
}

# Observe the real transport. No input, algorithm, response or cache is substituted.
INTEGRATION_OBSERVER = r"""(() => {
  const trace = window.__integration = {workers: [], messages: [], searches: []};
  const NativeWorker = window.Worker;
  window.Worker = class extends NativeWorker {
    constructor(url, options) {
      super(url, options);
      this.observation = {url: String(url), created: performance.now(), terminated: false};
      trace.workers.push(this.observation);
      this.pending = new Map();
      this.addEventListener('message', ({data}) => {
        const key = data?.requestId ?? data?.id ?? data?.type;
        const row = this.pending.get(key);
        if (row && data.type === row.type && typeof data.ok === 'boolean') {
          row.end = performance.now(); row.ok = data.ok;
          this.pending.delete(key);
        }
      });
    }
    postMessage(message, ...rest) {
      const row = {type: message.type, operation: message.operation,
        worker: this.observation.url, start: performance.now()};
      if (message.type === 'helper' || message.type === 'run' || message.type === 'init') {
        trace.messages.push(row);
        this.pending.set(message.requestId ?? message.id ?? message.type, row);
      }
      return super.postMessage(message, ...rest);
    }
    terminate() {
      this.observation.terminated = true;
      return super.terminate();
    }
  };
  window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
    const {runLosatPairsParallel} = await import('/gbdraw/web/js/services/losat.js');
    const row = {start: performance.now(), jobs: jobs.map(job => ({
      cacheKey: job.cacheKey, pairs: job.recordPairs, args: job.args
    }))};
    trace.searches.push(row);
    const result = await runLosatPairsParallel(jobs, options);
    row.end = performance.now(); row.completed = result.length;
    return result;
  };
})();"""


def integration_import(page, fixture):
    with page.expect_event('dialog', timeout=120_000) as event:
        page.locator('input[accept^=".json,"]').first.set_input_files(str(fixture))
    dialog = event.value
    assert dialog.message == 'Session loaded successfully!', dialog.message
    dialog.accept()
    page.wait_for_function('() => window.__GBDRAW_APP__.results.length > 0')


def integration_prepare(page, case, raw_cold=False):
    _, mode, infer, raw_limit = INTEGRATION_CASES[case]
    return page.evaluate('''async options => {
      const app = window.__GBDRAW_APP__;
      const {state} = await import('/gbdraw/web/js/state.js');
      await app.setLinearComparisonGlobalAction('losat');
      app.setLinearComparisonLosatMode('blastp');
      app.setLinearComparisonLosatpMode(options.mode);
      await Vue.nextTick();
      Object.assign(app.losat.blastp, {candidateLimit: options.rawLimit,
        orthogroupMemberMaxHits: 5, collinearInferOrthogroups: options.infer,
        collinearSearchScope: 'adjacent', collinearUnitMode: 'auto', collinearAnchorMode: 'rbh'});
      app.losat.executionMode = 'serial';
      if (options.rawCold) app.clearLosatCache();
      state.losatDerivedCache.value.clear();
      const {disposeDiagramGenerationWorker} = await import('/gbdraw/web/js/services/diagram-generation.js');
      disposeDiagramGenerationWorker();
      await Vue.nextTick();
      return {rawEntries: state.losatCache.value.size,
        derivedEntries: state.losatDerivedCache.value.size,
        resolved: state.files.linearCanonicalComparisons.length,
        settings: JSON.parse(JSON.stringify(app.losat)), records: app.linearSeqs.length};
    }''', {'mode': mode, 'infer': infer, 'rawLimit': raw_limit, 'rawCold': raw_cold})


def integration_generate(page, timed=False):
    return page.evaluate('''async timed => {
      const app = window.__GBDRAW_APP__;
      const trace = window.__integration;
      const startMessage = trace.messages.length, startSearch = trace.searches.length;
      const start = performance.now();
      const result = await app.runAnalysis();
      await Vue.nextTick();
      await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      const end = performance.now();
      const preview = document.querySelector('.shadow-xl.origin-top > svg');
      if (result.status !== 'ok' || !preview) throw new Error(JSON.stringify({result, error: app.errorLog}));
      const {state} = await import('/gbdraw/web/js/state.js');
      const hash = async text => Array.from(new Uint8Array(await crypto.subtle.digest(
        'SHA-256', new TextEncoder().encode(text))), b => b.toString(16).padStart(2,'0')).join('');
      const svg = String(app.results[app.selectedResultIndex].content);
      const geometry = [...preview.querySelectorAll('path')].map(p => p.getAttribute('d')).join('\\n');
      const row = {result, svgSha256: await hash(svg), geometrySha256: await hash(geometry),
        previewPathCount: preview.querySelectorAll('path').length,
        telemetry: JSON.parse(JSON.stringify(globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ || null)),
        rawKeys: [...state.losatCache.value.keys()].sort(),
        provenance: JSON.parse(JSON.stringify([...state.losatDerivedCache.value.values()].at(-1)?.payload?.provenance || null)),
        searches: trace.searches.slice(startSearch),
        messages: trace.messages.slice(startMessage),
        workerCount: trace.workers.length, liveWorkers: trace.workers.filter(w => !w.terminated).length};
      if (timed) Object.assign(row, {start, end, seconds: (end-start)/1000});
      return row;
    }''', timed)


def run_integration(args):
    from functools import partial
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
    import threading

    if args.integration == 'timing' and args.output.exists():
        raise ValueError('Timing evidence already exists; do not repeat completed samples.')

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *_):
            pass

    wheel, = (args.server_root / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    report = {'artifactRoot': str(args.server_root), 'wheelSha256': digest(wheel.read_bytes()),
              'runnerSha256': digest(Path(__file__).read_bytes()),
              'sourceSha256': {str(p.relative_to(ROOT)): digest(p.read_bytes())
                               for p in (ROOT / 'gbdraw').rglob('*.py')},
              'kind': args.integration, 'warmups': 0, 'cases': {}, 'runs': []}
    if args.integration == 'replacement' and args.resume:
        previous = json.loads(gzip.decompress(args.output.read_bytes()))
        for key in ('wheelSha256', 'sourceSha256', 'kind'):
            assert previous[key] == report[key], key
        prefix = previous
        while 'verifiedClearEvidence' in prefix:
            prefix = prefix['verifiedClearEvidence']
        cleared = next(row for row in prefix['inProgress']['steps'] if row['label'] == 'clear-cache')
        assert cleared['result'] == {'status': 'ok'} and cleared['searches']
        report['verifiedClearEvidence'] = previous
        report['resumptionKind'] = 'Reuse passed Clear Cache evidence; restart remaining checks from a fresh Gallery context'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        content = (json.dumps(report, indent=2, default=list) + '\n').encode()
        args.output.write_bytes(gzip.compress(content, mtime=0) if args.output.suffix == '.gz' else content)

    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(Handler, directory=str(args.server_root)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f'http://127.0.0.1:{server.server_port}'
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            report['browser'] = browser.version
            def context_for(viewport=(1280, 720)):
                context = browser.new_context(accept_downloads=True, viewport={'width': viewport[0], 'height': viewport[1]})
                record = {'viewport': viewport, 'external': [], 'pageErrors': [], 'paths': set()}
                def route(r):
                    url = urlparse(r.request.url)
                    if url.netloc != urlparse(base).netloc:
                        record['external'].append(r.request.url)
                        r.abort()
                    else:
                        record['paths'].add(url.path)
                        r.continue_()
                context.route('**/*', route)
                observer = INTEGRATION_OBSERVER
                if args.integration == 'manifest-merge':
                    # Corrupt only a real extraction response, before the client sees
                    # it. Success runs use the original transport and response bytes.
                    observer = observer.replace('const key = data?.requestId', """
                      if (window.__invalidManifest && data?.ok && data.result?.identity_manifest) {
                        data.result.identity_manifest.schema = -1;
                        window.__invalidManifestResponses++;
                      }
                      const key = data?.requestId""")
                context.add_init_script(observer)
                page = context.new_page()
                page.set_default_timeout(180_000)
                page.on('pageerror', lambda error: record['pageErrors'].append(str(error)))
                page.goto(base + '/gbdraw/web/index.html', wait_until='domcontentloaded')
                page.wait_for_function('() => window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length')
                assert page.evaluate('() => window.__integration.workers.length') == 0
                return context, page, record

            def close_context(context, page, record):
                record['beforeDispose'] = page.evaluate('() => window.__integration.workers')
                page.evaluate('''async () => {
                  const api = await import('/gbdraw/web/js/services/diagram-generation.js');
                  api.disposeDiagramGenerationWorker();
                }''')
                record['liveWorkersAfterDispose'] = page.evaluate('() => window.__integration.workers.filter(w => !w.terminated).length')
                assert record['liveWorkersAfterDispose'] == 0
                assert not record['external'] and not record['pageErrors'], record
                record['paths'] = sorted(record['paths'])
                report['runs'].append(record)
                context.close()

            try:
                if args.integration == 'manifest-merge':
                    run_manifest_merge(context_for, close_context, save)
                elif args.integration in {'lifecycle', 'replacement'}:
                    run_integration_lifecycle(args, report, context_for, close_context, save)
                elif args.integration == 'offline':
                    for viewport in ((1280, 720), (390, 844)):
                        context, page, record = context_for(viewport)
                        fixture = ROOT / 'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz'
                        integration_import(page, fixture)
                        record['preparation'] = integration_prepare(page, 'hep-on')
                        record['linear'] = integration_generate(page)
                        artifacts = ROOT / '.venv' / (args.output.name.split('.')[0] + '-artifacts')
                        artifacts.mkdir(parents=True, exist_ok=True)
                        for mode in ('linear', 'circular'):
                            if mode == 'circular':
                                page.evaluate('''async text => {
                                  const app = window.__GBDRAW_APP__;
                                  app.mode = 'circular'; app.cInputType = 'gb';
                                  app.files.c_gb = new File([text], 'HmmtDNA.gbk', {type:'text/plain'});
                                  app.files.c_gff = null; app.files.c_fasta = null;
                                  await Vue.nextTick();
                                }''', (ROOT / 'tests/test_inputs/HmmtDNA.gbk').read_text())
                                record['circular'] = integration_generate(page)
                                record['circularRepeat'] = integration_generate(page)
                                assert record['circular']['workerCount'] == record['linear']['workerCount']
                                assert record['circularRepeat']['geometrySha256'] == record['circular']['geometrySha256']
                            screenshot = artifacts / f'{viewport[0]}-{mode}.png'
                            page.screenshot(path=str(screenshot))
                            record[mode]['screenshot'] = str(screenshot.relative_to(ROOT))
                        close_context(context, page, record)
                        save()
                        print(viewport[0], 'offline Linear/Circular reuse and release passed', flush=True)
                else:
                    for case in args.cases:
                        fixture = ROOT / f'gbdraw/web/gallery/sessions/{INTEGRATION_CASES[case][0]}.gbdraw-session.json.gz'
                        values = report['cases'][case] = {'fixtureSha256': digest(fixture.read_bytes()), 'raw-cold': [], 'saved-raw': [], 'derived-warm': []}
                        # One real raw search; independent fresh-context saved-raw samples.
                        for state, count in (('raw-cold', 1), ('saved-raw', args.samples)):
                            for sample in range(count):
                                context, page, record = context_for()
                                record.update(case=case, state=state, sample=sample)
                                integration_import(page, fixture)
                                record['preparation'] = integration_prepare(page, case, state == 'raw-cold')
                                assert record['preparation']['resolved'] == 0
                                if state == 'raw-cold':
                                    assert record['preparation']['rawEntries'] == 0
                                observed = integration_generate(page, True)
                                values[state].append(observed)
                                save()
                                assert bool(observed['searches']) == (state == 'raw-cold'), observed['telemetry']
                                assert observed['telemetry']['proteinDerivedPayloadCacheMisses'] == 1
                                assert sum(len(row['jobs']) for row in observed['searches']) == observed['telemetry']['uniqueJobs']
                                if state == 'saved-raw':
                                    page.evaluate("() => window.__GBDRAW_APP__.setLinearComparisonGlobalAction('losat')")
                                    warm = integration_generate(page, True)
                                    values['derived-warm'].append(warm)
                                    save()
                                    assert not warm['searches']
                                    assert warm['telemetry']['proteinDerivedPayloadCacheHits'] == 1
                                    assert warm['geometrySha256'] == observed['geometrySha256']
                                close_context(context, page, record)
                                save()
                                print(case, state, sample + 1, 'passed', flush=True)
                        values['seconds'] = {}
                        for state in ('raw-cold', 'saved-raw', 'derived-warm'):
                            samples = [row['seconds'] for row in values[state]]
                            median = statistics.median(samples)
                            values['seconds'][state] = {'samples': samples, 'median': median,
                                'mad': statistics.median(abs(value-median) for value in samples) if len(samples)>1 else None}
                        save()
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
        save()


def run_manifest_merge(context_for, close_context, save):
    """Saved-raw merge and malformed extraction failure through the real Worker."""
    fixture = ROOT / 'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz'
    for viewport in ((1280, 720), (390, 844)):
        context, page, record = context_for(viewport)
        integration_import(page, fixture)
        record['fixtureSha256'] = digest(fixture.read_bytes())
        record['workersAfterPreviewLoad'] = page.evaluate('() => window.__integration.workers.length')
        assert record['workersAfterPreviewLoad'] == 0
        record['preparation'] = integration_prepare(page, 'hep-on')
        initial = record['savedRaw'] = integration_generate(page)
        assert not initial['searches']
        assert initial['telemetry']['proteinDerivedPayloadCacheMisses'] == 1
        assert initial['telemetry']['cacheHits'] == 13
        repeat = record['resolvedRepeat'] = integration_generate(page)
        assert repeat['svgSha256'] == initial['svgSha256']
        assert repeat['workerCount'] == initial['workerCount'] == 1
        assert not repeat['searches']
        # Give identical source bytes a fresh File owner to force extraction,
        # then corrupt its manifest. This does not replace the Worker or manufacture a Result.
        page.evaluate('''async () => {
          const app = window.__GBDRAW_APP__;
          await app.setLinearComparisonGlobalAction('losat');
          const file = app.linearSeqs[0].gb;
          window.__mergeOriginalFile = file;
          const {readFileText} = await import('/gbdraw/web/js/services/file-content-cache.js');
          app.setLinearSeqPrimaryFile(0, 'gb', new File([await readFileText(file)], file.name,
            {type: file.type, lastModified: file.lastModified}));
          await Vue.nextTick();
          window.__invalidManifestResponses = 0;
          window.__invalidManifest = true;
        }''')
        failed = record['invalidManifest'] = page.evaluate('''async () => {
          const app = window.__GBDRAW_APP__, history = window.__GBDRAW_HISTORY__;
          const {state} = await import('/gbdraw/web/js/state.js');
          const snapshot = () => JSON.stringify({results: app.results,
            selected: app.selectedResultIndex, undo: history.getUndoCount(), redo: history.getRedoCount(),
            raw: [...state.losatCache.value], derived: [...state.losatDerivedCache.value],
            manifest: state.proteinIdentityManifest.value});
          const before = snapshot();
          const trace = window.__integration;
          const messageStart = trace.messages.length, searchStart = trace.searches.length;
          const result = await app.runAnalysis();
          await Vue.nextTick();
          return {result, error: JSON.parse(JSON.stringify(app.errorLog)),
            preserved: before === snapshot(), undoCount: history.getUndoCount(),
            corruptResponses: window.__invalidManifestResponses,
            searches: trace.searches.slice(searchStart), messages: trace.messages.slice(messageStart),
            workerCount: trace.workers.length, liveWorkers: trace.workers.filter(w => !w.terminated).length,
            processing: app.processing};
        }''')
        assert failed['result'] == {'status': 'error'} and failed['preserved'], failed
        assert failed['error']['summary'] == (
            'Protein comparison metadata could not be validated. Reload the page and try again.'), failed
        assert failed['corruptResponses'] == 1 and failed['undoCount'] > 0, failed
        assert not failed['searches'] and not failed['processing'], failed
        assert all(message['type'] == 'helper' for message in failed['messages']), failed
        assert failed['workerCount'] == failed['liveWorkers'] == 1, failed
        page.evaluate('''async () => {
          window.__invalidManifest = false;
          const app = window.__GBDRAW_APP__;
          app.setLinearSeqPrimaryFile(0, 'gb', window.__mergeOriginalFile);
          await app.setLinearComparisonGlobalAction('losat');
          const {state} = await import('/gbdraw/web/js/state.js');
          state.losatDerivedCache.value.clear();
        }''')
        retry = record['savedRawRetry'] = integration_generate(page)
        for key in ('svgSha256', 'geometrySha256', 'provenance'):
            assert retry[key] == initial[key], key
        assert not retry['searches'] and retry['workerCount'] == 1
        close_context(context, page, record)
        save()
        print(viewport[0], 'saved raw, repeat, invalid manifest isolation and retry passed', flush=True)


def run_integration_lifecycle(args, report, context_for, close_context, save):
    """Real raw search, Python analysis/render, cancel/stale, retry and replacement."""
    artifacts = ROOT / '.venv' / (args.output.name.split('.')[0] + '-artifacts')
    artifacts.mkdir(parents=True, exist_ok=True)
    fixture = ROOT / 'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz'
    # Mobile regeneration/release is covered by offline mode and archived parity
    # flows. Exercise the longer real-search mutation sequence once on desktop.
    for viewport in ((1280, 720),):
        context, page, record = context_for(viewport)
        report['inProgress'] = record
        integration_import(page, fixture)
        record['fixtureSha256'] = digest(fixture.read_bytes())
        record['workersAfterPreviewLoad'] = page.evaluate('() => window.__integration.workers.length')
        assert record['workersAfterPreviewLoad'] == 0
        record['preparation'] = integration_prepare(page, 'hep-on')
        steps = record['steps'] = []
        def run(label):
            observed = integration_generate(page)
            observed['label'] = label
            steps.append(observed)
            save()
            print(viewport[0], label, 'passed', flush=True)
            return observed
        initial = run('saved-raw')
        assert not initial['searches']
        if args.integration != 'replacement':
            repeat = run('resolved-warm')
            assert repeat['geometrySha256'] == initial['geometrySha256']
            assert not repeat['searches']
            # A plain repeat consumes committed typed artifacts. The existing Run
            # LOSAT action removes those artifacts while retaining the derived cache.
            page.evaluate("() => window.__GBDRAW_APP__.setLinearComparisonGlobalAction('losat')")
            repeat = run('derived-warm')
            assert repeat['geometrySha256'] == initial['geometrySha256']
            assert repeat['telemetry']['proteinDerivedPayloadCacheHits'] == 1
            for label, expression in (
                ('color', "app.losat.blastp.collinearColorMode = 'orientation'"),
                ('block', 'app.losat.blastp.collinearMaxUnitGap = 1'),
                ('member', 'app.losat.blastp.orthogroupMemberMaxHits = 1'),
                ('filter', 'app.adv.min_bitscore = Number(app.adv.min_bitscore) + 1'),
                ('reverse', 'app.linearSeqs[0].region_reverse = !app.linearSeqs[0].region_reverse'),
            ):
                page.evaluate(f'() => {{ const app = window.__GBDRAW_APP__; {expression}; }}')
                observed = run(label)
                assert not observed['searches'], label
                assert observed['rawKeys'] == initial['rawKeys'], label
                assert observed['telemetry']['proteinDerivedPayloadCacheMisses'] == 1, label

            # All-record evidence permits reordered adjacent views to reuse raw.
            page.evaluate("() => { window.__GBDRAW_APP__.losat.blastp.collinearSearchScope = 'all'; }")
            run('all-record-evidence')
            page.evaluate('() => window.__GBDRAW_APP__.moveLinearSeqUp(1)')
            reordered = run('record-order')
            assert not reordered['searches']

            # Observe a completed real render, hold its response at the existing hook,
            # then cancel. Releasing that response must not replace the prior Result.
            page.evaluate('''() => {
              const app = window.__GBDRAW_APP__;
              window.__beforeCancel = app.results[app.selectedResultIndex].content;
              window.__beforeUndo = window.__GBDRAW_HISTORY__.getUndoCount();
              app.losat.blastp.collinearSearchScope = 'adjacent';
              app.losat.blastp.candidateLimit = 6;
              window.__responseHeld = false;
              const gate = new Promise(resolve => { window.__releaseResponse = resolve; });
              window.__GBDRAW_TEST_HOOKS__ = {beforeDiagramGenerationResponse() {
                window.__responseHeld = true; return gate;
              }};
              window.__cancelRunStatus = 'running';
              window.__cancelPromise = app.runAnalysis().then(result => {
                window.__cancelRunStatus = result; return result;
              });
            }''')
            try:
                page.wait_for_function("() => window.__responseHeld || window.__cancelRunStatus !== 'running'", timeout=240_000)
            finally:
                record['cancelBoundary'] = page.evaluate('''() => ({held: window.__responseHeld,
                  status: window.__cancelRunStatus, processing: window.__GBDRAW_APP__.processing,
                  error: window.__GBDRAW_APP__.errorLog, searches: window.__integration.searches,
                  messages: window.__integration.messages})''')
                save()
            assert record['cancelBoundary']['held'], record['cancelBoundary']
            before_cancel = page.evaluate('() => ({searches: window.__integration.searches, workers: window.__integration.workers})')
            assert before_cancel['searches'] and all(row.get('completed') for row in before_cancel['searches'])
            await_cancel = page.get_by_role('button', name=__import__('re').compile(r'Cancel$'))
            await_cancel.click()
            canceled = page.evaluate('''async () => {
              const result = await window.__cancelPromise;
              window.__releaseResponse();
              delete window.__GBDRAW_TEST_HOOKS__;
              await Vue.nextTick();
              await new Promise(resolve => requestAnimationFrame(resolve));
              return {result, preserved: window.__beforeCancel === window.__GBDRAW_APP__.results[window.__GBDRAW_APP__.selectedResultIndex].content,
                historyPreserved: window.__beforeUndo === window.__GBDRAW_HISTORY__.getUndoCount(),
                workers: window.__integration.workers};
            }''')
            assert canceled['result'] == {'status': 'canceled'} and canceled['preserved'] and canceled['historyPreserved']
            assert all(w['terminated'] for w in canceled['workers'] if 'diagram-generation-worker' in w['url'])
            record['cancelAfterRealRawAndRender'] = canceled
            page.evaluate('() => { window.__GBDRAW_APP__.losat.blastp.orthogroupMemberMaxHits = 5; }')
            retry = run('member-change-retry')
            assert not retry['searches']
            assert retry['workerCount'] > reordered['workerCount']
            assert retry['provenance']['orthogroup']['memberMaxHits'] == 5
        cleared = initial
        if not (args.integration == 'replacement' and args.resume):
            page.evaluate('''async () => {
              const app = window.__GBDRAW_APP__;
              app.clearLosatCache();
              await app.setLinearComparisonGlobalAction('losat');
              const {state} = await import('/gbdraw/web/js/state.js');
              if (state.losatCache.value.size || state.losatDerivedCache.value.size) throw new Error('Cache survived clear');
            }''')
            cleared = run('clear-cache')
            assert cleared['searches']

        # Change real biological source content, preserving its File owner.
        changed = page.evaluate(r'''async () => {
          const app = window.__GBDRAW_APP__;
          const file = app.linearSeqs[0].gb;
          const {readFileText} = await import('/gbdraw/web/js/services/file-content-cache.js');
          const text = await readFileText(file);
          const next = text.replace(/(\/translation="[A-Z])([A-Z])/, (full, start, amino) => start + (amino === 'A' ? 'G' : 'A'));
          if (text === next) throw new Error('No translation changed');
          app.setLinearSeqPrimaryFile(0, 'gb', new File([next], file.name, {type: file.type, lastModified: file.lastModified + 1}));
          return true;
        }''')
        assert changed
        source = run('changed-source')
        assert source['searches']
        assert source['rawKeys'] != cleared['rawKeys']

        history = page.evaluate('''async () => {
          const app = window.__GBDRAW_APP__, history = window.__GBDRAW_HISTORY__;
          const before = app.results[app.selectedResultIndex].content;
          const undone = await history.undo();
          const redo = await history.redo();
          return {undone, redo, restored: before === app.results[app.selectedResultIndex].content};
        }''')
        assert history == {'undone': True, 'redo': True, 'restored': True}
        record['history'] = history
        saved = acceptance._save_session(page, acceptance.AcceptanceChecks())
        target = artifacts / f'{viewport[0]}.gbdraw-session.json.gz'
        target.write_bytes(saved.read_bytes())
        record['savedSha256'] = digest(target.read_bytes())
        record['savedPath'] = str(target.relative_to(ROOT))
        close_context(context, page, record)
        save()
        context, page, restored_record = context_for(viewport)
        integration_import(page, target)
        restored_record['workersAfterPreviewLoad'] = page.evaluate('() => window.__integration.workers.length')
        assert restored_record['workersAfterPreviewLoad'] == 0
        # Use the saved mode/options, remove resolved/derived analysis only.
        page.evaluate('''async () => {
          await window.__GBDRAW_APP__.setLinearComparisonGlobalAction('losat');
          const {state} = await import('/gbdraw/web/js/state.js');
          state.losatDerivedCache.value.clear();
        }''')
        regenerated = integration_generate(page)
        assert not regenerated['searches']
        assert regenerated['geometrySha256'] == source['geometrySha256']
        restored_record['regenerated'] = regenerated
        close_context(context, page, restored_record)
        report.pop('inProgress', None)
        save()


if __name__ == '__main__':
    main()
