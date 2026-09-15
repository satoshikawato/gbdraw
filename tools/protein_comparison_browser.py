"""Browser measurement adapter for benchmark_protein_comparison.py (one CLI)."""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import threading
from urllib.parse import urlparse


def run_browser(root, payload, expected_keys, samples):
    from playwright.sync_api import sync_playwright

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def end_headers(self):
            self.send_header("Cross-Origin-Opener-Policy", "same-origin")
            self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
            super().end_headers()

    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(Handler, directory=str(root)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    requests, blocked, errors = [], [], []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            try:
                context = browser.new_context(viewport={"width": 1440, "height": 1000})
                def route(request_route):
                    url = request_route.request.url
                    requests.append(url)
                    if urlparse(url).netloc != urlparse(base).netloc:
                        blocked.append(url)
                        request_route.abort()
                    else:
                        request_route.continue_()
                context.route("**/*", route)
                page = context.new_page()
                page.set_default_timeout(120_000)
                page.on("pageerror", lambda e: errors.append(str(e)))
                # Observe the production worker transport; no runtime algorithms,
                # cancellation results, helpers or payloads are substituted.
                page.add_init_script("""(() => {
                  window.__proteinBenchmarkMessages = [];
                  const original = Worker.prototype.postMessage;
                  const workers = new WeakMap();
                  Worker.prototype.postMessage = function(payload, ...rest) {
                    if (!workers.has(this)) {
                      const pending = new Map(); workers.set(this, pending);
                      this.addEventListener('message', event => {
                        const sample = pending.get(event.data?.requestId);
                        if (sample && event.data?.type === 'helper') {
                          sample.roundTripMs = performance.now() - sample.started;
                          sample.ok = event.data.ok;
                          pending.delete(event.data.requestId);
                        }
                      });
                    }
                    const start = performance.now();
                    const value = original.call(this, payload, ...rest);
                    if (payload?.type === 'helper') {
                      const sample = { operation: payload.operation,
                        started: start, postMessageMs: performance.now() - start };
                      workers.get(this).set(payload.requestId, sample);
                      window.__proteinBenchmarkMessages.push(sample);
                    }
                    return value;
                  };
                })();""")
                page.goto(base + "/gbdraw/web/index.html", wait_until="domcontentloaded")
                page.wait_for_function("() => Boolean(window.__GBDRAW_APP__)")
                page.evaluate("payload => { window.__proteinBenchmarkInput = payload; }", payload)
                result = page.evaluate("""async ({expectedKeys, samples}) => {
                  const api = await import('/gbdraw/web/js/services/diagram-generation.js');
                  const input = window.__proteinBenchmarkInput;
                  const serializedInputBytes = new TextEncoder().encode(JSON.stringify(input)).length;
                  const runs = [];
                  for (let i = 0; i <= samples; i++) {
                    const start = performance.now();
                    const {result: value} = await api.runDiagramHelperOperation(
                      api.DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEYS, input);
                    const elapsedMs = performance.now() - start;
                    if (JSON.stringify(value.keys) !== JSON.stringify(expectedKeys)) {
                      throw new Error('Browser / native ordered cache keys differ');
                    }
                    runs.push({kind: i === 0 ? 'cold-including-worker-init' : 'warm', elapsedMs});
                  }
                  const {prepareLosatSourceBatches} = await import('/gbdraw/web/js/app/linear-sources.js');
                  const hashText = async value => Array.from(new Uint8Array(await crypto.subtle.digest(
                    'SHA-256', new TextEncoder().encode(value))), b => b.toString(16).padStart(2, '0')).join('');
                  const files = [new File(['source A'], 'A.gb'), new File(['source B'], 'B.gb')];
                  const sequences = Array.from({length: 8}, (_, i) => ({uid: `record-${i}`,
                    gb: files[i < 6 ? 0 : 1], gff: null, fasta: null}));
                  const all = sequences.flatMap((_, queryIndex) => sequences.map((_, subjectIndex) => ({queryIndex, subjectIndex})));
                  const plans = [];
                  for (const infer of [true, false]) {
                    const specs = all.filter(x => infer || x.queryIndex !== x.subjectIndex);
                    const plan = await prepareLosatSourceBatches({sequences, specs,
                      getEntry: async i => ({fasta: `>protein-${i}\\nMKK\\n`}),
                      buildArgs: () => [], hashText, protein: true, excludeSelfComparisons: !infer});
                    if (!infer && plan.batches.some(b => b.query.indexes.some(i => b.subject.indexes.includes(i)))) {
                      throw new Error('OFF source batch searches a record against itself');
                    }
                    plans.push({infer, sources: 2, records: 8, directionalTables: specs.length,
                      sourceJobs: plan.batches.length, executedSearchInvocations: 0,
                      searchContexts: plan.batches.map(b => b.searchContext)});
                  }
                  if (plans[0].sourceJobs !== 4 || plans[1].sourceJobs !== 34) throw new Error('Source batching changed');
                  const messages = window.__proteinBenchmarkMessages;
                  if (messages.length !== samples + 1 || messages.some(m => !m.ok || !Number.isFinite(m.roundTripMs))) {
                    throw new Error('Incomplete real Worker helper observations');
                  }
                  return {runs, messages, serializedInputBytes, plans,
                    keysEqualNative: true, keyCount: expectedKeys.length,
                    crossOriginIsolated, heap: performance.memory ? {
                      usedJSHeapSize: performance.memory.usedJSHeapSize,
                      totalJSHeapSize: performance.memory.totalJSHeapSize,
                      scope: 'Chromium main-thread heap snapshot; not Pyodide/Wasm peak'
                    } : null};
                }""", {"expectedKeys": expected_keys, "samples": samples})
                result.update({"browser": browser.version, "viewport": {"width": 1440, "height": 1000},
                               "network": "fresh context, external requests blocked before navigation",
                               "requestedPaths": [u.replace(base, "<local>") for u in requests],
                               "blockedExternalRequests": blocked, "pageErrors": errors,
                               "artifact": "source SPA plus locally built browser wheel",
                               "limits": ["helper round trip includes clone, queue, Worker JSON and Pyodide work",
                                          "postMessage time measures synchronous enqueue/clone only",
                                          "no isolated network-transfer or pure Pyodide timing",
                                          "no LOSAT search, full Generate, mobile, or Wasm peak measurement"]})
                if blocked or errors:
                    raise ValueError(f"Browser boundary errors: blocked={blocked}, errors={errors}")
                return result
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
