import assert from 'node:assert/strict';
import test from 'node:test';

// OV-3027, Owner D-42: when the diagram engine (Python in the Worker) runs out of memory, the page replaces it at
// once. The Worker's heap never shrinks, so the Worker that raised the MemoryError is terminated by the one owner of
// the Worker lifecycle (terminateWorker in services/diagram-generation.js), and the next request builds a new one
// (getWorker), which receives every input again. Generate, a committed live redraw and a label reflow all run through
// runDiagramGeneration; a helper operation through runAuxiliaryWorkerRequest.
globalThis.location = { href: 'https://example.test/gbdraw/web/' };
const metrics = {};
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric({ name, value }) { metrics[name] = (metrics[name] || 0) + Number(value || 0); }
};
const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);
class EngineWorker {
  static instances = [];
  constructor() {
    this.listeners = new Map();
    this.messages = [];
    this.terminated = false;
    EngineWorker.instances.push(this);
  }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, new Set());
    this.listeners.get(type).add(listener);
  }
  removeEventListener(type, listener) { this.listeners.get(type)?.delete(listener); }
  postMessage(message) {
    this.messages.push(message);
    if (message.type === 'init') {
      queueMicrotask(() => this.dispatch('message', { type: 'init', id: message.id, ok: true,
        capabilities: EXPECTED_WEB_RUNTIME_CAPABILITIES }));
    }
  }
  dispatch(type, data) {
    for (const listener of [...(this.listeners.get(type) || [])]) listener({ data });
  }
  terminate() { this.terminated = true; }
}
globalThis.Worker = EngineWorker;
const { runDiagramGeneration, runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS, disposeDiagramGenerationWorker } = await import(
  '../../gbdraw/web/js/services/diagram-generation.js'
);
const { normalizeUserFacingError, OUT_OF_MEMORY_CODE } = await import('../../gbdraw/web/js/utils/error-normalization.js');

const tick = () => new Promise((resolve) => setImmediate(resolve));
const SOURCE_BYTES = Buffer.from('ACGT').toString('base64');
const payload = () => ({
  request: { records: [{ source: { kind: 'genbank', resourceId: 'source' } }] },
  resources: { source: { name: 'source.gb', kind: 'genbank', size: 4, encoding: 'base64', data: SOURCE_BYTES } }
});
/** The newest message of `type` that the newest Worker received, once it is not `seen`. */
const sent = async (type, seen = null) => {
  for (let i = 0; i < 10_000; i += 1) {
    const worker = EngineWorker.instances.at(-1);
    const message = worker?.messages.findLast((entry) => entry.type === type);
    if (message && message !== seen) return { worker, message };
    await tick();
  }
  assert.fail(`no ${type} message was sent`);
};
const engineError = (operation, code) => ({ code, operation, stage: operation === 'generate' ? 'render' : 'helper', context: {} });
// The inputs a run message carries: each staged resource by its id, with its byte count.
const staged = (message) => (message.payload.stagedResources ?? []).map(({ resourceId, bytes }) => [resourceId, bytes.byteLength]);

test.afterEach(() => disposeDiagramGenerationWorker());

test('the wording owner words the code the service replaces the engine on', () => {
  const model = normalizeUserFacingError(engineError('generate', OUT_OF_MEMORY_CODE));
  assert.equal(model.code, OUT_OF_MEMORY_CODE);
  assert.deepEqual(model.actions, ['retry', 'save-session']);
});

test('a Generate that runs out of memory replaces the engine, and the next Generate runs on a new one', async () => {
  const constructed = metrics.workerConstructionCount || 0;
  const first = runDiagramGeneration(payload());
  const { worker, message } = await sent('run');
  assert.deepEqual(staged(message), [['source', 4]], 'the first engine receives the input');
  worker.dispatch('message', { type: 'run', requestId: message.requestId, ok: true,
    results: { error: engineError('generate', OUT_OF_MEMORY_CODE) } });
  const failed = await first;
  assert.equal(failed.results.error.code, OUT_OF_MEMORY_CODE, 'the caller words the failure');
  assert.equal(worker.terminated, true, 'the engine that ran out of memory is terminated');
  assert.equal([...worker.listeners.values()].reduce((total, set) => total + set.size, 0), 0);

  const retry = runDiagramGeneration(payload());
  const again = await sent('run');
  assert.notEqual(again.worker, worker, 'the next Generate builds a new Worker');
  assert.equal(metrics.workerConstructionCount, constructed + 2);
  assert.ok(again.worker.messages.some(({ type }) => type === 'init'), 'the new engine boots');
  assert.deepEqual(staged(again.message), [['source', 4]], 'the new engine receives the input again');
  again.worker.dispatch('message', { type: 'run', requestId: again.message.requestId, ok: true, results: [] });
  assert.deepEqual((await retry).results, []);
  assert.equal(again.worker.terminated, false);
});

test('another engine failure keeps the engine', async () => {
  const constructed = metrics.workerConstructionCount || 0;
  const first = runDiagramGeneration(payload());
  const { worker, message } = await sent('run');
  worker.dispatch('message', { type: 'run', requestId: message.requestId, ok: true,
    results: { error: { ...engineError('generate', 'RENDER_FAILED'), context: { exceptionType: 'ValidationError' } } } });
  assert.equal((await first).results.error.code, 'RENDER_FAILED');
  assert.equal(worker.terminated, false);
  const next = runDiagramGeneration(payload());
  const again = await sent('run', message);
  assert.equal(again.worker, worker, 'the same engine serves the next Generate');
  again.worker.dispatch('message', { type: 'run', requestId: again.message.requestId, ok: true, results: [] });
  await next;
  assert.equal(metrics.workerConstructionCount, constructed + 1);
});

test('a helper operation that runs out of memory replaces the engine and ends the Generate waiting on it', async () => {
  const constructed = metrics.workerConstructionCount || 0;
  const generate = runDiagramGeneration(payload());
  const { worker } = await sent('run');
  const helper = runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, { kind: 'color', features: [], rules: [] });
  const { message } = await sent('helper');
  worker.dispatch('message', { type: 'helper', requestId: message.requestId, ok: false,
    error: engineError(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, OUT_OF_MEMORY_CODE) });
  const helperError = await helper.then(() => assert.fail('the helper succeeded'), (error) => error);
  assert.equal(helperError.code, OUT_OF_MEMORY_CODE);
  assert.equal(worker.terminated, true);
  // The Generate posted to the replaced engine would never be answered: it ends with the same failure.
  const generateError = await generate.then(() => assert.fail('the Generate succeeded'), (error) => error);
  assert.deepEqual([generateError.code, generateError.operation], [OUT_OF_MEMORY_CODE, 'generate']);
  assert.deepEqual(generateError.actions, ['retry', 'save-session']);

  const retry = runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, { kind: 'color', features: [], rules: [] });
  const again = await sent('helper');
  assert.notEqual(again.worker, worker);
  assert.equal(metrics.workerConstructionCount, constructed + 2);
  again.worker.dispatch('message', { type: 'helper', requestId: again.message.requestId, ok: true, result: { winners: [] } });
  assert.deepEqual((await retry).result, { winners: [] });
});
