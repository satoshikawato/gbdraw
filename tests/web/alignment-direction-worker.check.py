"""Run the S01 source projection through the packaged, offline browser Worker.

Prepare the browser wheel, then run this file with Python Playwright installed.
"""
from __future__ import annotations

import base64
import functools
import json
import sys
import tempfile
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tests.test_alignment_direction_projection import build_fixture  # noqa: E402


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass


def main():
    with tempfile.TemporaryDirectory(prefix='gbdraw-s01-browser-') as directory:
        _request, helper, context, paths, _response, _sources = build_fixture(Path(directory))
        resources = {}
        for key, source_path in paths.items():
            data = Path(source_path).read_bytes()
            resources[key] = dict(encoding='base64', data=base64.b64encode(data).decode(), size=len(data),
                                  name=Path(source_path).name, kind='genbank')
        server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(QuietHandler, directory=str(ROOT)))
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(headless=True)
                page = browser.new_page()
                failures, external_requests = [], []
                page.on('pageerror', lambda error: failures.append(str(error)))
                def route(request_route):
                    if request_route.request.url.startswith(f'http://127.0.0.1:{server.server_port}/'):
                        request_route.continue_()
                    else:
                        external_requests.append(request_route.request.url)
                        request_route.abort()
                page.route('**/*', route)
                page.add_init_script("""
window.__S01_METRICS__ = [];
window.__GBDRAW_TEST_HOOKS__ = {onStructuralMetric:event=>window.__S01_METRICS__.push(event)};
""")
                page.goto(f'http://127.0.0.1:{server.server_port}/gbdraw/web/index.html', wait_until='domcontentloaded')
                page.wait_for_function('window.__GBDRAW_APP__')
                result = page.evaluate("""async ({helper,projection,resources}) => {
const {runDiagramHelperOperation,runDiagramGeneration,DIAGRAM_HELPER_OPERATIONS,cancelDiagramGeneration,isDiagramGenerationCanceled} = await import('./js/services/diagram-generation.js');
const {validateSimilarityAlignmentResolution,projectSimilarityAlignmentDirections} = await import('./js/app/similarity-alignment.js');
const {state} = await import('./js/state.js');
const {getCommittedCanonicalRenderRequest} = await import('./js/services/config.js');
const prior = JSON.stringify({results:state.results.value,request:getCommittedCanonicalRenderRequest()});
const operation = DIAGRAM_HELPER_OPERATIONS.RESOLVE_SIMILARITY_ALIGNMENT;
const first = (await runDiagramHelperOperation(operation,{request:helper,projection,resources})).result;
const preview = projectSimilarityAlignmentDirections({resolution:validateSimilarityAlignmentResolution(first,helper),intent:{mode:'right'},expectedBinding:first.projection.binding});
projection.orientations = Object.fromEntries(preview.records.map(row=>[row.recordKey,row.afterReverseComplement]));
const second = (await runDiagramHelperOperation(operation,{request:helper,projection,resources})).result;
const final = projectSimilarityAlignmentDirections({resolution:validateSimilarityAlignmentResolution(second,helper),intent:{mode:'right'},expectedBinding:preview.binding});
const plain = (await runDiagramHelperOperation(operation,{request:helper})).result;
let rejected = false;
try { await runDiagramHelperOperation(operation,{request:helper,projection:{...projection,unexpected:true},resources}); }
catch(error) { rejected = Boolean(error.message); }
const retry = (await runDiagramHelperOperation(operation,{request:helper,projection,resources})).result;
const generated = await runDiagramGeneration({request:projection.canonicalRequest,resources});
const output = {preview,final,generationResults:generated.results.length,plainProjection:plain.projection,rejected,retryBinding:retry.projection.binding,
 unchanged:prior===JSON.stringify({results:state.results.value,request:getCommittedCanonicalRenderRequest()}),
 workerCount:window.__S01_METRICS__.filter(event=>event.name==='workerConstructionCount').length,
 cacheHits:window.__S01_METRICS__.filter(event=>event.name==='workerResourceCacheHitCount').length};
const metricHook = window.__GBDRAW_TEST_HOOKS__.onStructuralMetric;
let cancelNextPreparation = true;
window.__GBDRAW_TEST_HOOKS__.onStructuralMetric = event => {
  metricHook(event);
  if (cancelNextPreparation && event.name === 'workerResourceCacheHitCount') {
    cancelNextPreparation = false;
    cancelDiagramGeneration();
  }
};
try { await runDiagramHelperOperation(operation,{request:helper,projection,resources}); }
catch(error) { output.canceledPreparation = isDiagramGenerationCanceled(error); }
await new Promise(resolve=>setTimeout(resolve,0));
output.workersAfterCancellation = window.__S01_METRICS__.filter(event=>event.name==='workerConstructionCount').length;
const afterCancel = (await runDiagramHelperOperation(operation,{request:helper,projection,resources})).result;
output.retryAfterCancelBinding = afterCancel.projection.binding;
output.workersAfterRetry = window.__S01_METRICS__.filter(event=>event.name==='workerConstructionCount').length;
cancelDiagramGeneration();
return output;
}""", dict(helper=helper, projection=context, resources=resources))
                assert result['final']['geometryValidated'] is True, result
                assert result['final']['signature'] == result['preview']['signature'], result
                assert result['final']['reference']['beforeX'] == result['final']['reference']['afterX'], result
                assert result['final']['binding'] == result['retryBinding'], result
                assert result['generationResults'] > 0, result
                assert result['plainProjection'] is None and result['unchanged'] and result['rejected'], result
                assert result['workerCount'] == 1 and result['cacheHits'] >= 6, result
                assert result['canceledPreparation'] and result['workersAfterCancellation'] == 1, result
                assert result['workersAfterRetry'] == 2 and result['retryAfterCancelBinding'] == result['final']['binding'], result
                assert not failures and not external_requests, (failures, external_requests)
                browser.close()
                print(json.dumps(dict(status='PASS', arrows=[row['afterArrow'] for row in result['final']['records']],
                                      reference=result['final']['reference'], workerCount=result['workerCount'],
                                      cacheHits=result['cacheHits'], generationResults=result['generationResults'], canceledPreparation=result['canceledPreparation'],
                                      successfulRetryAfterCancel=result['workersAfterRetry'] == 2, unchangedArtifact=result['unchanged'],
                                      externalRequests=external_requests), indent=2))
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


if __name__ == '__main__':
    main()
