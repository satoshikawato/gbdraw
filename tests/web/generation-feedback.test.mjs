import assert from 'node:assert/strict';
import test from 'node:test';
import { assertKnownDefect } from './helpers/known-defect.mjs';

// Hold real client boundaries; no clock drives progress.
globalThis.location = { href: 'https://example.test/gbdraw/web/' };
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
    for (const listener of this.listeners.get('message') || []) listener({ data });
  }
  initialize() {
    const { id } = this.messages.find(({ type }) => type === 'init');
    this.emit({ type: 'init', id, ok: true, capabilities: EXPECTED_WEB_RUNTIME_CAPABILITIES });
  }
  terminate() { this.terminated = true; }
}
globalThis.Worker = ControlledWorker;
const { runDiagramGeneration, runDiagramHelperOperation, cancelDiagramGeneration, disposeDiagramGenerationWorker } = await import(
  '../../gbdraw/web/js/services/diagram-generation.js'
);
const flush = async () => { for (let i = 0; i < 12; i += 1) await Promise.resolve(); };
const payload = () => ({ request: {}, resources: {} });
const latestRequest = (worker) => worker.messages.filter(({ type }) => type === 'run').at(-1).requestId;
const assertClean = (worker) => {
  for (const listeners of worker.listeners.values()) assert.equal(listeners.size, 0);
};

test('Generate feedback follows cold/warm work, rejects duplicates, and ignores stale and settled progress', async () => {
  const events = [];
  const first = runDiagramGeneration(payload(), { onProgress: event => events.push(event) });
  const worker = ControlledWorker.instances.at(-1);
  assert.deepEqual(events.map(e => e.stage), ['preparing-runtime']);
  await flush();
  assert.equal(events.length, 1, 'pending initialization must not advance by itself');
  await assert.rejects(runDiagramGeneration(payload()), /already running/);
  worker.initialize();
  await flush();
  const requestId = latestRequest(worker);
  assert.deepEqual(events.map(e => e.stage), ['preparing-runtime', 'preparing-resources']);
  worker.emit({ type: 'progress', requestId: requestId + 1, stage: 'rendering' });
  assert.equal(events.length, 2);
  worker.emit({ type: 'progress', requestId, stage: 'rendering' });
  worker.emit({ type: 'progress', requestId, stage: 'finalizing' });
  assert.deepEqual(events.map(e => e.stage), ['preparing-runtime', 'preparing-resources', 'rendering', 'finalizing']);
  assert.ok(events.every(e => e.requestId === requestId));
  worker.emit({ type: 'run', requestId, ok: true, results: [{ name: 'out.svg', content: '<svg/>' }] });
  assert.equal((await first).results[0].name, 'out.svg');
  assertClean(worker);
  worker.emit({ type: 'progress', requestId, stage: 'rendering' });
  assert.equal(events.length, 4);

  const warmEvents = [];
  const second = runDiagramGeneration(payload(), { onProgress: event => warmEvents.push(event) });
  await flush();
  assert.deepEqual(warmEvents.map(e => e.stage), ['preparing-resources']);
  const secondId = latestRequest(worker);
  worker.emit({ type: 'progress', requestId, stage: 'finalizing' });
  assert.equal(warmEvents.length, 1);
  worker.emit({ type: 'progress', requestId: secondId, stage: 'rendering' });
  worker.emit({ type: 'run', requestId: secondId, ok: true, results: [] });
  await second;
  assert.equal(ControlledWorker.instances.length, 1);
  assertClean(worker);
  disposeDiagramGenerationWorker();
});

test('Cancel and failure clean progress listeners and permit a fresh run', async () => {
  for (const terminal of ['cancel-init', 'cancel-render', 'error']) {
    const events = [];
    const pending = runDiagramGeneration(payload(), { onProgress: event => events.push(event) });
    const rejection = assert.rejects(pending, terminal === 'error' ? error => error.code === 'UNKNOWN' && error.stage === 'render' : /cancel/i);
    const worker = ControlledWorker.instances.at(-1);
    if (terminal !== 'cancel-init') {
      worker.initialize();
      await flush();
      const requestId = latestRequest(worker);
      worker.emit({ type: 'progress', requestId, stage: 'rendering' });
      if (terminal === 'error') {
        worker.emit({ type: 'run', requestId, ok: false, error: { code: 'UNKNOWN', operation: 'generate', stage: 'render' } });
      }
    }
    if (terminal !== 'error') await cancelDiagramGeneration();
    await rejection;
    const count = events.length;
    worker.emit({ type: 'progress', requestId: 1, stage: 'finalizing' });
    assert.equal(events.length, count);
    assertClean(worker);
    assert.equal(worker.terminated, true);
  }
  const retry = runDiagramGeneration(payload());
  const worker = ControlledWorker.instances.at(-1);
  worker.initialize();
  await flush();
  worker.emit({ type: 'run', requestId: latestRequest(worker), ok: true, results: [] });
  await retry;
  assertClean(worker);
  disposeDiagramGenerationWorker();
});


test('Generate helper preparation reports real cold initialization and shares the warm render Worker', async () => {
  const events = [];
  const helper = runDiagramHelperOperation('evaluateRules', { rules: [] }, {
    onProgress: event => events.push(event)
  });
  const worker = ControlledWorker.instances.at(-1);
  const constructions = ControlledWorker.instances.length;
  assert.deepEqual(events.map(event => event.stage), ['preparing-runtime']);
  await flush();
  assert.equal(events.length, 1);
  worker.initialize();
  await flush();
  const requestId = worker.messages.find(message => message.type === 'helper').requestId;
  assert.equal(events[0].requestId, requestId);
  worker.emit({ type: 'helper', requestId, ok: true, result: { rules: [] } });
  await helper;
  assertClean(worker);
  const renderEvents = [];
  const render = runDiagramGeneration(payload(), { onProgress: event => renderEvents.push(event) });
  await flush();
  assert.deepEqual(renderEvents.map(event => event.stage), ['preparing-resources']);
  worker.emit({ type: 'run', requestId: latestRequest(worker), ok: true, results: [] });
  await render;
  assert.equal(ControlledWorker.instances.length, constructions);
  assertClean(worker);
  disposeDiagramGenerationWorker();

  const canceledEvents = [];
  const canceled = runDiagramHelperOperation('evaluateRules', { rules: [] }, {
    onProgress: event => canceledEvents.push(event)
  });
  const rejection = assert.rejects(canceled, /cancel/i);
  const canceledWorker = ControlledWorker.instances.at(-1);
  assert.deepEqual(canceledEvents.map(event => event.stage), ['preparing-runtime']);
  cancelDiagramGeneration();
  await rejection;
  assert.equal(canceledWorker.terminated, true);
  assertClean(canceledWorker);
  assert.equal(canceledEvents.length, 1);
});

// GE-09 (Web GUI audit 2026-09-30): a Cancel while no Worker-side request is
// active (for example, during JS-side preparation) keeps the warm Worker, so the
// next run is not a cold start.
test('Cancel without Worker-side work keeps the warm Worker (GE-09 known defect)', async () => {
  const helper = runDiagramHelperOperation('evaluateRules', { rules: [] });
  const worker = ControlledWorker.instances.at(-1);
  worker.initialize();
  await flush();
  const requestId = worker.messages.find(message => message.type === 'helper').requestId;
  worker.emit({ type: 'helper', requestId, ok: true, result: { rules: [] } });
  await helper;
  assertClean(worker);
  const constructions = ControlledWorker.instances.length;
  const canceled = cancelDiagramGeneration();
  await assertKnownDefect('GE-09', () => {
    assert.deepEqual(
      { canceled, terminated: Boolean(worker.terminated), constructions: ControlledWorker.instances.length },
      { canceled: false, terminated: false, constructions }
    );
  });
  disposeDiagramGenerationWorker();
});
