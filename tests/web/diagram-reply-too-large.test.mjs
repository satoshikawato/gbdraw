import assert from 'node:assert/strict';
import test from 'node:test';

// PD-OI-118 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md: the Worker drew the diagram, but the page
// cannot receive its reply (the run's `messageerror`). Generate, a committed live redraw and a label reflow all
// run through runDiagramGeneration, so they all get this error; Retry would fail the same way.
globalThis.location = { href: 'https://example.test/gbdraw/web/' };
const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);
class ReplyWorker {
  static instances = [];
  constructor() {
    this.listeners = new Map();
    this.messages = [];
    ReplyWorker.instances.push(this);
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
globalThis.Worker = ReplyWorker;
const { runDiagramGeneration, disposeDiagramGenerationWorker } = await import(
  '../../gbdraw/web/js/services/diagram-generation.js'
);

const REPLY_TOO_LARGE_TEXT = 'The diagram was drawn, but it is too large for this browser page to receive. '
  + 'Draw a smaller region, or turn off tracks or labels you do not need, or draw it with the gbdraw command line '
  + '(Save Session, then `gbdraw circular --session` or `gbdraw linear --session`).';
const RESULT_INVALID_TEXT = 'The generated result could not be accepted. Retry the operation.';

const tick = () => new Promise((resolve) => setImmediate(resolve));
const payload = () => ({
  request: { records: [{ source: { kind: 'genbank', resourceId: 'source' } }] },
  resources: { source: { name: 'source.gb', kind: 'genbank', size: 4, encoding: 'base64', data: Buffer.from('ACGT').toString('base64') } }
});
// A Generate whose run message the Worker has received.
const sentRun = async () => {
  const pending = runDiagramGeneration(payload());
  const worker = ReplyWorker.instances.at(-1);
  for (let i = 0; i < 10_000 && !worker.messages.some(({ type }) => type === 'run'); i += 1) await tick();
  const run = worker.messages.find(({ type }) => type === 'run');
  assert.ok(run, 'the run message was sent');
  return { pending, worker, requestId: run.requestId };
};
const listenerCount = (worker) => [...worker.listeners.values()].reduce((total, listeners) => total + listeners.size, 0);

test.afterEach(() => disposeDiagramGenerationWorker());

test('a run reply the page cannot receive is REPLY_TOO_LARGE with Save Session only', async () => {
  const { pending, worker } = await sentRun();
  worker.dispatch('messageerror', null);
  const error = await pending.then(() => assert.fail('the Generate was admitted'), (failure) => failure);
  assert.equal(error.code, 'RESULT_INVALID');
  assert.equal(error.operation, 'generate');
  assert.equal(error.stage, 'result-admission');
  assert.deepEqual(error.context, { reason: 'REPLY_TOO_LARGE' });
  assert.equal(error.summary, REPLY_TOO_LARGE_TEXT);
  assert.equal(error.message, REPLY_TOO_LARGE_TEXT);
  assert.deepEqual(error.actions, ['save-session']);
  assert.equal(worker.terminated, true, 'the Worker that holds the undelivered reply is released');
  assert.equal(listenerCount(worker), 0);
});

test('a run reply that fails admission for another cause keeps its text and Retry', async () => {
  const { pending, worker, requestId } = await sentRun();
  // A Result whose children bytes are not bytes: admission refuses it.
  worker.dispatch('message', { type: 'run', requestId, ok: true,
    results: { results: [{ name: 'out.svg', document: { head: '<svg>', tail: '</svg>', pieces: [] }, segments: 'not bytes' }] } });
  const error = await pending.then(() => assert.fail('the Generate was admitted'), (failure) => failure);
  assert.equal(error.code, 'RESULT_INVALID');
  assert.equal(error.stage, 'result-admission');
  assert.deepEqual(error.context, {});
  assert.equal(error.summary, RESULT_INVALID_TEXT);
  assert.deepEqual(error.actions, ['retry', 'save-session']);
});
