import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { normalizeUserFacingError } from '../../gbdraw/web/js/services/error-normalization.js';

const ownerUrl = new URL('../../gbdraw/web/js/services/losat.js', import.meta.url);
const ownerSource = (await readFile(ownerUrl, 'utf8'))
  .replace(/from '([^']+)'/g, (_match, path) => `from '${new URL(path, ownerUrl)}'`)
  .replaceAll('import.meta.url', JSON.stringify(ownerUrl.href));
const wasmBytes = await readFile(new URL('../../gbdraw/web/wasm/losat/losat-threaded.wasm', import.meta.url));
const privateText = '/private/specimen.fasta >patient_sequence\nACGT PRIVATE_LOSAT_STDERR';
const job = { pairIndex: 0, program: 'blastp', querySequenceKey: 'q', subjectSequenceKey: 's' };
const options = {
  totalThreadBudget: 2, concurrency: 1, threadsPerJob: 2,
  sequences: { q: '>q\n' + 'A'.repeat(300000), s: '>s\n' + 'A'.repeat(300000) }
};
let instance = 0;

const setup = async (t, { supportError, threadedError, holdThreaded = false } = {}) => {
  const original = Object.fromEntries(['window', 'crossOriginIsolated', 'fetch', 'Worker'].map(key => [key, Object.getOwnPropertyDescriptor(globalThis, key)]));
  const logs = [];
  for (const level of ['info', 'warn']) {
    const previous = console[level];
    console[level] = (...args) => logs.push(args);
    t.after(() => { console[level] = previous; });
  }
  t.after(() => {
    for (const [key, descriptor] of Object.entries(original)) {
      if (descriptor) Object.defineProperty(globalThis, key, descriptor);
      else delete globalThis[key];
    }
  });
  const workers = [];
  globalThis.window = { location: { href: 'https://example.test/gbdraw/web/' } };
  globalThis.crossOriginIsolated = true;
  globalThis.fetch = async () => {
    if (supportError) throw supportError;
    return { ok: true, arrayBuffer: async () => wasmBytes.buffer.slice(wasmBytes.byteOffset, wasmBytes.byteOffset + wasmBytes.byteLength) };
  };
  // Keep the real module compilation and orchestration; replace only Worker transport.
  globalThis.Worker = class {
    listeners = new Map();
    messages = [];
    constructor(url) { this.threaded = String(url).includes('losat-threaded-worker.js'); workers.push(this); }
    addEventListener(type, listener) {
      if (!this.listeners.has(type)) this.listeners.set(type, new Set());
      this.listeners.get(type).add(listener);
    }
    removeEventListener(type, listener) { this.listeners.get(type)?.delete(listener); }
    terminate() { this.terminated = true; }
    postMessage(message) {
      this.messages.push(message);
      if (holdThreaded && this.threaded && message.type === 'run') return;
      queueMicrotask(() => {
        const data = { id: message.id, type: message.type, ok: true, text: 'q\ts\t100\n' };
        if (this.threaded) Object.assign(data, threadedError
          ? { ok: false, error: threadedError }
          : { stderr: privateText });
        for (const listener of this.listeners.get('message') || []) listener({ data });
      });
    }
  };
  const owner = await import(`data:text/javascript;base64,${Buffer.from(`${ownerSource}\n// instance ${++instance}`).toString('base64')}`);
  return { owner, workers, logs };
};

test('threaded support failure preserves a known code without exposing the raw asset diagnostic', async t => {
  const error = Object.assign(new Error(privateText), { code: 'WORKER_INIT', stage: 'initialization' });
  const { owner, workers, logs } = await setup(t, { supportError: error });
  const status = await owner.getLosatThreadingSupport();
  assert.equal(status.state, 'unavailable');
  assert.equal(status.message, normalizeUserFacingError(error, { operation: 'generate', stage: 'initialization' }).summary);
  assert.doesNotMatch(JSON.stringify({ status, logs }), /private|patient_sequence|PRIVATE_LOSAT_STDERR/);
  assert.equal(workers.length, 0);
});

test('automatic threaded failure retains serial results and emits only bounded status and logs', async t => {
  const { owner, workers, logs } = await setup(t, { threadedError: privateText });
  const statuses = [];
  const result = await owner.runLosatPairsParallel([job], { ...options, onRuntimeStatus: status => statuses.push(status) });
  assert.deepEqual(result, [{ ...job, text: 'q\ts\t100\n' }]);
  assert.deepEqual(statuses.map(status => status.state), ['running', 'fallback']);
  assert.equal(statuses[1].mode, 'serial');
  assert.equal(statuses[1].fallbackReason, normalizeUserFacingError(new Error(privateText), { operation: 'generate', stage: 'helper' }).summary);
  assert.doesNotMatch(JSON.stringify({ statuses, logs }), /private|patient_sequence|PRIVATE_LOSAT_STDERR/);
  assert.deepEqual(workers.map(worker => worker.threaded), [true, false]);
  assert.ok(workers.every(worker => worker.terminated));
});

test('successful threaded results and progress do not automatically publish Worker stderr', async t => {
  const { owner, workers, logs } = await setup(t);
  const statuses = [];
  const progress = [];
  const result = await owner.runLosatPairsParallel([job], {
    ...options, executionMode: 'threaded', onRuntimeStatus: status => statuses.push(status),
    onProgress: value => progress.push(value)
  });
  assert.deepEqual(result, [{ ...job, text: 'q\ts\t100\n' }]);
  assert.deepEqual(statuses.map(status => status.state), ['running', 'available']);
  assert.deepEqual(progress, [{ completed: 1, total: 1, job, index: 0, threaded: true, spawnCount: 0 }]);
  assert.doesNotMatch(JSON.stringify({ statuses, logs }), /private|patient_sequence|PRIVATE_LOSAT_STDERR/);
  assert.equal(workers.length, 1);
  assert.ok(workers[0].terminated);
});

test('explicit threaded failure is still rejected without serial fallback', async t => {
  const { owner, workers, logs } = await setup(t, { threadedError: privateText });
  await assert.rejects(owner.runLosatPairsParallel([job], { ...options, executionMode: 'threaded' }), error => error.message.includes(privateText));
  assert.deepEqual(workers.map(worker => worker.threaded), [true]);
  assert.ok(workers[0].terminated);
  assert.doesNotMatch(JSON.stringify(logs), /private|patient_sequence|PRIVATE_LOSAT_STDERR/);
});

test('cancellation still rejects the initiating reason and never starts serial fallback', async t => {
  const { owner, workers, logs } = await setup(t, { holdThreaded: true });
  const controller = new AbortController();
  const result = owner.runLosatPairsParallel([job], { ...options, signal: controller.signal });
  while (!workers.some(worker => worker.messages.some(message => message.type === 'run'))) {
    await new Promise(resolve => setImmediate(resolve));
  }
  const reason = new Error('initiating cancellation');
  controller.abort(reason);
  await assert.rejects(result, error => error === reason);
  assert.deepEqual(workers.map(worker => worker.threaded), [true]);
  assert.ok(workers[0].terminated);
  assert.doesNotMatch(JSON.stringify(logs), /private|patient_sequence|PRIVATE_LOSAT_STDERR/);
});
