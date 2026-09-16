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


# This diagnostic observer is appended only by the memory-run HTTP route. It is
# never packaged and does not replace a production algorithm or helper response.
_PATH_WORKER_OBSERVER = r'''
const _s06ProductionMessage = self.onmessage;
self.onmessage = event => {
  if (event.data?.type !== '__s06_probe') return _s06ProductionMessage(event);
  const {requestId, action} = event.data;
  try {
    const py = runtime.pyodide;
    let result;
    if (action === 'arm') {
      result = JSON.parse(String(py.runPython(`
import json as _s06_json, tracemalloc as _s06_trace
from gbdraw.analysis import protein_colinearity as _s06_pc
from gbdraw.analysis.ortholog_paths import OrthologPathCollection as _s06_collection
_s06_counts = {"pathObjects": 0, "iteratorCalls": 0, "adapterCalls": 0}
def _s06_wrap(function, name):
    def counted(*args, **kwargs):
        _s06_counts[name] += 1
        return function(*args, **kwargs)
    return counted
_s06_collection._path = _s06_wrap(_s06_collection._path, "pathObjects")
_s06_collection.iter_paths = _s06_wrap(_s06_collection.iter_paths, "iteratorCalls")
_s06_pc.materialize_ortholog_paths = _s06_wrap(_s06_pc.materialize_ortholog_paths, "adapterCalls")
_s06_trace.start()
_s06_json.dumps({"armed": True})
`)));
    } else if (action === 'memory') {
      result = JSON.parse(String(py.runPython(`
_s06_retained, _s06_peak = _s06_trace.get_traced_memory()
_s06_json.dumps({"pythonRetainedBytes": _s06_retained, "pythonPeakBytes": _s06_peak,
                 "operations": _s06_counts,
                 "filteredEntries": len(_WEB_LOSATP_FILTERED_HIT_CACHE),
                 "convertedEntries": len(_WEB_LOSATP_CONVERTED_PAYLOAD_CACHE),
                 "convertedUtf8Bytes": sum(len(v.encode()) for v in _WEB_LOSATP_CONVERTED_PAYLOAD_CACHE.values())})
`)));
    } else if (action === 'large') {
      result = JSON.parse(String(py.runPython(`
_s06_trace.stop()
from dataclasses import replace as _s06_replace
from gbdraw.session_request_codec import encode_canonical_typed_resource as _s06_encode
_s06_pm = {f"p{i}": _s06_pc.CdsProtein(protein_id=f"p{i}", record_index=i, feature_index=0,
           record_id=f"r{i}", start=0, end=300, strand=1, label=f"p{i}", protein_length=100, sequence="M"*100)
           for i in range(58)}
def _s06_edge(i, j):
    return _s06_pc.OrthologEdge("og_1", "og_1", "og_1", f"p{i}", f"p{j}", i, j,
        "rbh", "block_anchor", None, 90., 1e-30, 200., 100)
_s06_edges = [_s06_edge(i, j) for i in range(56) for j in range(i+1, 56)] + [_s06_edge(56, 57)]
_s06_updated, _s06_indexes = _s06_pc._build_ortholog_path_indexes({"og_1": _s06_edges}, _s06_pm)
_s06_result = _s06_pc.OrthogroupGraphResult({}, {}, ortholog_edges_by_orthogroup_id=_s06_updated,
                                         path_indexes_by_orthogroup_id=_s06_indexes)
_s06_index = _s06_indexes["og_1"]
_s06_wire = _s06_json.loads(_s06_encode("orthogroupResult", _s06_result))
_s06_before = dict(_s06_counts)
_s06_ranks = [1, 2, _s06_index.count//2, _s06_index.count-1, _s06_index.count]
assert all(_s06_index.rank_of(_s06_index.path_at(rank).protein_ids) == rank for rank in _s06_ranks)
_s06_json.dumps({"resource": _s06_wire, "count": str(_s06_index.count),
                 "ranks": list(map(str, _s06_ranks)), "normalOperations": _s06_before})
`)));
    } else throw new Error('Unknown S06 probe');
    self.postMessage({type: '__s06_probe', requestId, ok: true, result,
                      wasmCapacityBytes: py._module.HEAP8.buffer.byteLength});
  } catch (error) {
    self.postMessage({type: '__s06_probe', requestId, ok: false, error: String(error)});
  }
};
'''

_PATH_PAGE_OBSERVER = r'''(() => {
  const OriginalWorker = window.Worker;
  window.__s06 = {workers: [], messages: [], terminated: 0};
  window.Worker = class extends OriginalWorker {
    constructor(...args) {
      super(...args);
      if (String(args[0]).includes('diagram-generation-worker')) window.__s06.workers.push(this);
      this.addEventListener('error', event => console.error('Worker error: ' + event.message));
      this.addEventListener('message', event => {
        if (window.__s06MeasureTransfer && event.data?.type === 'helper') {
          window.__s06.messages.push({direction: 'response', requestId: event.data.requestId,
            // JSON transport is an object clone, not an ArrayBuffer. This is
            // its exact UTF-8 JSON size, not undocumented V8 clone storage.
            jsonUtf8Bytes: new TextEncoder().encode(JSON.stringify(event.data)).byteLength});
        }
      });
    }
    postMessage(message, ...rest) {
      if (window.__s06MeasureTransfer && message?.type === 'helper') {
        window.__s06.messages.push({direction: 'request', requestId: message.requestId,
          transferredArrayBufferBytes: (message.payload?.files || []).reduce((sum, f) => sum + f.bytes.byteLength, 0)});
      }
      return super.postMessage(message, ...rest);
    }
    terminate() {
      window.__s06.terminated++;
      window.__s06.workers = window.__s06.workers.filter(worker => worker !== this);
      return super.terminate();
    }
  };
})();'''


def run_path_browser(root, inputs, samples, measurement):
    """Real production helper/Worker, offline, cold/repeat/changed settings.

    Timing excludes instrumentation and separates runtime initialization. Memory
    uses a separate run and a diagnostic observer in the actual worker. Search,
    rendering and session acceptance are measured/checked by their own runners.
    """
    import statistics
    from playwright.sync_api import sync_playwright

    if samples < 1:
        raise ValueError("samples must be positive")

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
    runs = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-precise-memory-info"])
            try:
                for width, height in ((1280, 720), (390, 844)):
                    observed = {"viewport": {"width": width, "height": height}, "blockedExternal": [],
                                "requests": set(), "errors": [], "cases": {}}
                    context = browser.new_context(viewport=observed["viewport"])

                    def route(request_route):
                        url = urlparse(request_route.request.url)
                        if url.netloc != urlparse(base).netloc:
                            observed["blockedExternal"].append(request_route.request.url)
                            request_route.abort()
                        elif measurement == "memory" and url.path.endswith("/workers/diagram-generation-worker.js"):
                            observed["requests"].add(url.path)
                            source = (root / url.path.lstrip("/")).read_text()
                            request_route.fulfill(status=200, content_type="text/javascript", headers={"Cross-Origin-Embedder-Policy": "require-corp", "Cross-Origin-Opener-Policy": "same-origin"}, body=source + _PATH_WORKER_OBSERVER)
                        else:
                            observed["requests"].add(url.path)
                            request_route.continue_()

                    context.route("**/*", route)
                    page = context.new_page()
                    page.set_default_timeout(120_000)
                    page.on("pageerror", lambda error: observed["errors"].append(str(error)))
                    page.on("console", lambda message: print(message.text, flush=True) if message.type == "error" else None)
                    page.add_init_script("window.__s06MeasureTransfer = " + str(measurement == "memory").lower() + ";\n" + _PATH_PAGE_OBSERVER)
                    page.goto(base + "/gbdraw/web/index.html", wait_until="domcontentloaded")
                    page.wait_for_function("() => Boolean(window.__GBDRAW_APP__)")
                    for case in inputs:
                        print(f"path-browser {measurement}: {width} {case['name']}", flush=True)
                        result = page.evaluate(r'''async ({input, samples, measurement}) => {
                          const api = await import('/gbdraw/web/js/services/diagram-generation.js');
                          const encoder = new TextEncoder();
                          const messages = window.__s06.messages;
                          const heap = () => performance.memory ? {
                            used: performance.memory.usedJSHeapSize, capacity: performance.memory.totalJSHeapSize
                          } : null;
                          let probeId = 0;
                          const probe = action => new Promise((resolve, reject) => {
                            const worker = window.__s06.workers.at(-1);
                            const id = ++probeId;
                            const listener = event => {
                              if (event.data?.type !== '__s06_probe' || event.data.requestId !== id) return;
                              worker.removeEventListener('message', listener);
                              if (!event.data.ok) reject(new Error(JSON.stringify(event.data.error))); else resolve(event.data);
                            };
                            worker.addEventListener('message', listener);
                            worker.postMessage({type: '__s06_probe', action, requestId: id});
                          });
                          const payload = parameters => ({...parameters, files: [
                            {role: 'pairs', bytes: encoder.encode(input.pairsText).buffer},
                            {role: 'rawTsv', bytes: encoder.encode(input.rawText).buffer}
                          ]});
                          const summary = value => {
                            if (value.error) throw new Error(JSON.stringify(value.error));
                            const resource = value.orthogroupResult || value.collinearityResult;
                            if (resource?.schema !== 3) throw new Error('Not a current typed resource');
                            const groups = resource.value.type === 'CollinearityResult'
                              ? resource.value.fields.orthogroups : resource.value;
                            const indexes = groups?.fields?.pathIndexesByOrthogroupId || {};
                            if (input.parameters.collinearInferOrthogroups && groups?.type !== 'OrthogroupGraphResult') {
                              throw new Error('Default inference is not compact');
                            }
                            if (groups?.fields?.orthologPathsByOrthogroupId) throw new Error('Exhaustive default resource');
                            const counts = Object.values(indexes).map(v => v.fields.count);
                            if (counts.some(v => typeof v !== 'string' || !/^(0|[1-9][0-9]*)$/.test(v))) throw new Error('Inexact path count');
                            const total = counts.reduce((sum, v) => sum + BigInt(v), 0n).toString();
                            if (input.name === 'path-24' && total !== '4194304') throw new Error('R24 count mismatch');
                            return {pathCount: total, groups: Object.keys(indexes).length,
                              nodes: Object.values(indexes).reduce((n,v) => n + (v.fields.nodes?.length || 0), 0),
                              edges: Object.values(indexes).reduce((n,v) => n + (v.fields.transitions?.length || 0), 0),
                              serializedBytes: encoder.encode(JSON.stringify(value)).byteLength};
                          };
                          const rows = [];
                          const repetitions = measurement === 'timing' ? samples + 1 : 1;
                          for (let i = 0; i < repetitions; i++) {
                            api.disposeDiagramGenerationWorker();
                            const initStart = performance.now();
                            await api.runDiagramHelperOperation(api.DIAGRAM_HELPER_OPERATIONS.VALIDATE_CONFIG_OVERRIDES,
                              {mode: 'linear', configOverrides: {}});
                            const initializationMs = performance.now() - initStart;
                            const beforeHeap = heap();
                            if (measurement === 'memory') await probe('arm');
                            const row = {warmup: measurement === 'timing' && i === 0, initializationMs, beforeHeap, phases: {}};
                            let prior = null;
                            for (const phase of ['cold', 'repeat', 'changed']) {
                              const parameters = {...input.parameters};
                              if (phase === 'changed') parameters.bitscore = Number(parameters.bitscore) + 1;
                              const request = payload(parameters);
                              const start = performance.now();
                              const {result: value} = await api.runDiagramHelperOperation(
                                api.DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD, request);
                              const elapsedMs = performance.now() - start;
                              const info = summary(value);
                              const science = payload => Object.fromEntries(Object.entries(payload).filter(([key]) => key !== 'cache'));
                              if (phase === 'repeat' && JSON.stringify(science(value)) !== JSON.stringify(science(prior))) throw new Error('Repeated science changed');
                              if (Boolean(value.cache?.convertedPayloadHit) !== (phase === 'repeat')) throw new Error('Unexpected converted cache state');
                              const memory = measurement === 'memory' ? await probe('memory') : null;
                              if (memory && Object.values(memory.result.operations).some(n => n !== 0)) throw new Error('Normal consumer enumerated paths');
                              row.phases[phase] = {elapsedMs, ...info, cache: value.cache,
                                simultaneouslyRetainedResultJsonBytes: info.serializedBytes + (prior ? encoder.encode(JSON.stringify(prior)).byteLength : 0),
                                heap: heap(), memory, transfer: messages.slice(-2)};
                              prior = value;
                            }
                            if (measurement === 'memory' && input.name === 'path-24') {
                              const large = await probe('large');
                              const exact = '18014398509481985';
                              if (large.result.count !== exact || BigInt(large.result.count).toString() !== exact) throw new Error('Big integer transfer changed');
                              const count = large.result.resource.value.fields.pathIndexesByOrthogroupId.og_1.fields.count;
                              if (count !== exact || JSON.parse(JSON.stringify(large.result)).count !== exact) throw new Error('Typed bigint count changed');
                              if (Object.values(large.result.normalOperations).some(n => n !== 0)) throw new Error('Large normal construction expanded paths');
                              row.large = {count, ranks: large.result.ranks,
                                bytes: encoder.encode(JSON.stringify(large.result.resource)).byteLength,
                                normalOperations: large.result.normalOperations};
                            }
                            rows.push(row);
                          }
                          api.disposeDiagramGenerationWorker();
                          if (window.__s06.workers.length) throw new Error('Analysis Worker survived disposal');
                          return {rows, liveWorkers: window.__s06.workers.length, terminatedWorkers: window.__s06.terminated, rawSearchInvocations: 0};
                        }''', {"input": case, "samples": samples, "measurement": measurement})
                        if measurement == "timing":
                            result["timing"] = {}
                            for phase in ("cold", "repeat", "changed"):
                                values = [row["phases"][phase]["elapsedMs"] for row in result["rows"] if not row["warmup"]]
                                median = statistics.median(values)
                                mad = statistics.median(abs(v - median) for v in values)
                                result["timing"][phase] = {"samples": values, "medianMs": median, "madMs": mad,
                                                           "noisePct": 100 * mad / median}
                        observed["cases"][case["name"]] = {"input": case["inventory"], "parameters": case["parameters"], **result}
                    observed["requests"] = sorted(observed["requests"])
                    if observed["blockedExternal"] or observed["errors"]:
                        raise ValueError(f"Browser boundary failure: {observed}")
                    runs.append(observed)
                    context.close()
                return {"browser": browser.version, "measurement": measurement, "runs": runs,
                        "artifact": "source SPA and generated local S06 wheel",
                        "scope": "production conversion helper; runtime init separated; no fresh raw search or full Generate timing",
                        "memoryScopes": "Python allocations since arming; Wasm capacity, not live Python bytes; main-thread JS heap snapshots, not peaks; actual transferred buffers plus response JSON UTF-8 size, not V8 clone allocation",
                        "observer": "memory run only: diagnostic listener appended to actual Worker; no production algorithms replaced"}
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
