import assert from 'node:assert/strict';
import test from 'node:test';

// D7: on a cross-origin-isolated page, a Cancel stops Python through the
// Worker's interrupt buffer and keeps the Worker; without an answer in time it
// terminates the Worker, as a Cancel did before.
globalThis.location = { href: 'https://example.test/gbdraw/web/' };
globalThis.crossOriginIsolated = true;
const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);
class ControlledWorker {
  static instances = [];
  constructor() {
    this.listeners = new Map();
    this.messages = [];
    ControlledWorker.instances.push(this);
  }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, new Set());
    this.listeners.get(type).add(listener);
  }
  removeEventListener(type, listener) { this.listeners.get(type)?.delete(listener); }
  postMessage(message) { this.messages.push(message); }
  emit(data) {
    for (const listener of [...(this.listeners.get('message') || [])]) listener({ data });
  }
  initialize() {
    const { id } = this.messages.find(({ type }) => type === 'init');
    this.emit({ type: 'init', id, ok: true, capabilities: EXPECTED_WEB_RUNTIME_CAPABILITIES });
  }
  sent(type) { return this.messages.filter((message) => message.type === type); }
  drain() {
    const { id } = this.sent('cancel').at(-1);
    this.emit({ type: 'cancel', id, ok: true });
  }
  terminate() { this.terminated = true; }
}
globalThis.Worker = ControlledWorker;
const {
  DIAGRAM_WORKER_CANCEL_TIMEOUT_MS, runDiagramGeneration, runDiagramHelperOperation, cancelDiagramGeneration,
  disposeDiagramGenerationWorker
} = await import('../../gbdraw/web/js/services/diagram-generation.js');

const tick = () => new Promise((resolve) => setImmediate(resolve));
const flush = async () => { for (let i = 0; i < 20; i += 1) await tick(); };
// Waits for the page's asynchronous resource preparation (sha256) to post the next run.
const runsSent = async (worker, count) => {
  for (let i = 0; i < 10_000 && worker.sent('run').length < count; i += 1) await tick();
  assert.equal(worker.sent('run').length, count);
};
const payload = (text) => ({
  request: { records: [{ source: { kind: 'genbank', resourceId: 'source' } }] },
  resources: { source: { name: 'source.gb', kind: 'genbank', size: text.length, encoding: 'base64', data: Buffer.from(text).toString('base64') } }
});
const latestRun = (worker) => worker.sent('run').at(-1);
const transferred = (worker) => latestRun(worker).payload.stagedResources.map(({ resourceId }) => resourceId);
const listenerCount = (worker) => [...worker.listeners.values()].reduce((total, listeners) => total + listeners.size, 0);
const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
const drainOutcomes = () => metrics.filter(({ name }) => name === 'workerCancelDrain').map(({ outcome }) => outcome);

/** A booted Worker that has rendered once and holds the resource. */
const warmWorker = async () => {
  const first = runDiagramGeneration(payload('AAAA'));
  const worker = ControlledWorker.instances.at(-1);
  worker.initialize();
  await runsSent(worker, 1);
  assert.deepEqual(transferred(worker), ['source']);
  worker.emit({ type: 'run', requestId: latestRun(worker).requestId, ok: true, results: [] });
  await first;
  return worker;
};

test.afterEach(() => {
  disposeDiagramGenerationWorker();
  metrics.length = 0;
});

test('the init message carries an interrupt buffer in memory shared with the Worker', async () => {
  const worker = await warmWorker();
  const [{ interruptBuffer }] = worker.sent('init');
  assert.ok(interruptBuffer instanceof Uint8Array);
  assert.ok(interruptBuffer.buffer instanceof SharedArrayBuffer);
  assert.equal(interruptBuffer.length, 1);
  assert.equal(interruptBuffer[0], 0);
});

test('a Cancel during a render interrupts Python, keeps the Worker, and the next render waits for its answer', async () => {
  const worker = await warmWorker();
  const [{ interruptBuffer }] = worker.sent('init');
  const canceled = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const { requestId } = latestRun(worker);
  const rejected = assert.rejects(canceled, /canceled/);
  assert.equal(cancelDiagramGeneration(), true);
  await rejected;
  assert.equal(worker.terminated, undefined);
  assert.equal(interruptBuffer[0], 2);
  const [cancel] = worker.sent('cancel');
  assert.equal(cancel.through, requestId);
  assert.equal(typeof cancel.id, 'string');

  const next = runDiagramGeneration(payload('AAAA'));
  await flush();
  assert.equal(worker.sent('run').length, 2, 'nothing is sent before the Worker answers');
  worker.emit({ type: 'cancel', id: 'another-cancel', ok: true });
  await flush();
  assert.equal(worker.sent('run').length, 2, 'only the answer to this Cancel counts');
  worker.drain();
  await runsSent(worker, 3);
  assert.equal(ControlledWorker.instances.at(-1), worker);
  assert.equal(worker.sent('init').length, 1);
  assert.deepEqual(transferred(worker), ['source'], 'the Worker pruned for the canceled request: the bytes go again');
  worker.emit({ type: 'run', requestId: latestRun(worker).requestId, ok: true, results: [] });
  await next;
  assert.deepEqual(drainOutcomes(), ['drained']);
  assert.equal(listenerCount(worker), 0);
});

test('without an answer in time the Worker is terminated and the waiting render boots a new one', async (t) => {
  t.mock.timers.enable({ apis: ['setTimeout'] });
  const worker = await warmWorker();
  const canceled = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const rejected = assert.rejects(canceled, /canceled/);
  cancelDiagramGeneration();
  await rejected;
  const next = runDiagramGeneration(payload('AAAA'));
  t.mock.timers.tick(DIAGRAM_WORKER_CANCEL_TIMEOUT_MS - 1);
  await flush();
  assert.equal(worker.terminated, undefined);
  t.mock.timers.tick(1);
  await flush();
  assert.equal(worker.terminated, true);
  assert.deepEqual(drainOutcomes(), ['terminated']);
  const replacement = ControlledWorker.instances.at(-1);
  assert.notEqual(replacement, worker);
  replacement.initialize();
  await runsSent(replacement, 1);
  assert.deepEqual(transferred(replacement), ['source']);
  replacement.emit({ type: 'run', requestId: latestRun(replacement).requestId, ok: true, results: [] });
  await next;
  assert.equal(listenerCount(replacement), 0);
});

test('a Worker error during the drain terminates it', async () => {
  const worker = await warmWorker();
  const canceled = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const rejected = assert.rejects(canceled, /canceled/);
  cancelDiagramGeneration();
  await rejected;
  for (const listener of [...(worker.listeners.get('error') || [])]) listener({});
  assert.equal(worker.terminated, true);
  assert.deepEqual(drainOutcomes(), ['terminated']);
});

test('a Cancel during boot terminates the Worker as before', async () => {
  const pending = runDiagramGeneration(payload('AAAA'));
  const worker = ControlledWorker.instances.at(-1);
  const rejected = assert.rejects(pending, /canceled/);
  assert.equal(cancelDiagramGeneration(), true);
  await rejected;
  assert.equal(worker.terminated, true);
  assert.equal(worker.sent('cancel').length, 0);
});

test('a Cancel rejects auxiliary requests too, and a second Cancel during the drain writes and drains nothing new', async () => {
  const worker = await warmWorker();
  const [{ interruptBuffer }] = worker.sent('init');
  const helper = runDiagramHelperOperation('evaluateRules', { rules: [] });
  const render = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const rejected = Promise.all([assert.rejects(helper, /canceled/), assert.rejects(render, /canceled/)]);
  cancelDiagramGeneration();
  await rejected;
  assert.equal(worker.sent('cancel').at(-1).through, latestRun(worker).requestId);
  interruptBuffer[0] = 0;
  const during = runDiagramGeneration(payload('AAAA'));
  const duringRejected = assert.rejects(during, /canceled/);
  assert.equal(cancelDiagramGeneration(), true);
  await duringRejected;
  assert.equal(interruptBuffer[0], 0);
  assert.equal(worker.sent('cancel').length, 1);
  worker.drain();
  await flush();
  assert.equal(worker.sent('run').length, 2, 'the request canceled during the drain never reached the Worker');
  assert.equal(worker.terminated, undefined);
  assert.equal(listenerCount(worker), 0);
});

test('a Generate that waits for the drain and is disposed boots no Worker', async () => {
  const worker = await warmWorker();
  const canceled = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const rejected = assert.rejects(canceled, /canceled/);
  cancelDiagramGeneration();
  await rejected;
  const waiting = runDiagramGeneration(payload('AAAA'));
  const waitingRejected = assert.rejects(waiting, /disposed/);
  await flush();
  const helper = runDiagramHelperOperation('evaluateRules', { rules: [] });
  const helperRejected = assert.rejects(helper, /disposed/);
  await flush();
  const constructed = ControlledWorker.instances.length;
  disposeDiagramGenerationWorker();
  await Promise.all([waitingRejected, helperRejected]);
  await flush();
  assert.equal(worker.terminated, true);
  assert.equal(ControlledWorker.instances.length, constructed, 'no Worker boots for a request that no longer exists');
});

test('a second Cancel during the drain, then the timeout, boots no Worker for the canceled requests', async (t) => {
  t.mock.timers.enable({ apis: ['setTimeout'] });
  const worker = await warmWorker();
  const canceled = runDiagramGeneration(payload('AAAA'));
  await runsSent(worker, 2);
  const rejected = assert.rejects(canceled, /canceled/);
  cancelDiagramGeneration();
  await rejected;
  const waiting = runDiagramGeneration(payload('AAAA'));
  const helper = runDiagramHelperOperation('evaluateRules', { rules: [] });
  const waitingRejected = Promise.all([assert.rejects(waiting, /canceled/), assert.rejects(helper, /canceled/)]);
  await flush();
  assert.equal(cancelDiagramGeneration(), true);
  await waitingRejected;
  const constructed = ControlledWorker.instances.length;
  t.mock.timers.tick(DIAGRAM_WORKER_CANCEL_TIMEOUT_MS);
  await flush();
  assert.equal(worker.terminated, true);
  assert.deepEqual(drainOutcomes(), ['terminated']);
  assert.equal(ControlledWorker.instances.length, constructed, 'no Worker boots for a request that no longer exists');
});
