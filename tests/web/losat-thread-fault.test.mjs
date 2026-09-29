import assert from 'node:assert/strict';
import test from 'node:test';

// A WASI thread that traps leaves its job worker blocked in wasm, so the trap
// must reach the page that owns the job. These tests use synthetic modules:
// Chromium traps memory.copy/memory.fill in a thread whose memory size is stale
// after another thread grew shared memory, and LOSAT threads hit that trap.
const webRoot = new URL('../../gbdraw/web/', import.meta.url);
const threadWorkerUrl = new URL('js/workers/losat-wasi-thread-worker.js', webRoot);
const losatUrl = new URL('js/services/losat.js', webRoot);
const wasiShimUrl = new URL('vendor/browser_wasi_shim/dist/index.js', webRoot).href;

const leb = (value) => {
  const bytes = [];
  do {
    let byte = value & 0x7f;
    value >>>= 7;
    if (value) byte |= 0x80;
    bytes.push(byte);
  } while (value);
  return bytes;
};
const name = (text) => [...leb(text.length), ...Buffer.from(text)];
const section = (id, entries) => {
  const body = [...leb(entries.length), ...entries.flat()];
  return [id, ...leb(body.length), ...body];
};
const moduleBytes = (...sections) =>
  new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00, ...sections.flat()]);
const I32 = 0x7f;
const sharedMemory = (min, max) => [0x02, 0x03, ...leb(min), ...leb(max)];
const trappingThreadStart = [...leb(3), 0x00, 0x00, 0x0b]; // no locals; unreachable; end

// Thread side only: wasi_thread_start traps at once.
const threadTrapModule = moduleBytes(
  section(1, [[0x60, 2, I32, I32, 0]]),
  section(2, [[...name('env'), ...name('memory'), ...sharedMemory(1, 1)]]),
  section(3, [[0]]),
  section(7, [[...name('wasi_thread_start'), 0x00, 0]]),
  section(10, [trappingThreadStart])
);

// Threaded-LOSAT shape: _start spawns one thread, then waits on a futex that
// only that thread could notify; the thread traps.
const threadedJobModule = moduleBytes(
  section(1, [[0x60, 1, I32, 1, I32], [0x60, 0, 0], [0x60, 2, I32, I32, 0]]),
  section(2, [
    [...name('wasi'), ...name('thread-spawn'), 0x00, 0],
    [...name('env'), ...name('memory'), ...sharedMemory(1, 16384)]
  ]),
  section(3, [[1], [2]]),
  section(7, [[...name('_start'), 0x00, 1], [...name('wasi_thread_start'), 0x00, 2]]),
  section(10, [
    [...leb(18), 0x00, 0x41, 0x00, 0x10, 0x00, 0x1a,
      0x41, 0x00, 0x41, 0x00, 0x42, 0x7f, 0xfe, 0x01, 0x02, 0x00, 0x1a, 0x0b],
    trappingThreadStart
  ])
);

const settleWithin = (promise, ms, label) => Promise.race([
  promise,
  new Promise((_, reject) => setTimeout(() => reject(new Error(`${label} did not settle within ${ms} ms`)), ms).unref())
]);

const restoreGlobalsAfter = (t, keys) => {
  const original = Object.fromEntries(keys.map((key) => [key, Object.getOwnPropertyDescriptor(globalThis, key)]));
  t.after(() => {
    for (const [key, descriptor] of Object.entries(original)) {
      if (descriptor) Object.defineProperty(globalThis, key, descriptor);
      else delete globalThis[key];
    }
  });
};

test('a WASI thread that traps reports the trap on its job fault channel', async (t) => {
  assert.ok(WebAssembly.validate(threadTrapModule));
  restoreGlobalsAfter(t, ['self']);
  const posted = [];
  globalThis.self = { postMessage: (message) => posted.push(message), close() {} };
  await import(`${threadWorkerUrl.href}?instance=${Date.now()}`);
  const onmessage = globalThis.self.onmessage;

  const faultChannel = `gbdraw-test-thread-fault:${process.pid}:${Date.now()}`;
  const receiver = new BroadcastChannel(faultChannel);
  t.after(() => receiver.close());
  const fault = new Promise((resolve) => { receiver.onmessage = (event) => resolve(event.data); });

  await onmessage({
    data: {
      type: 'prepare', id: 'prepare-1', readyControl: new SharedArrayBuffer(4),
      module: await WebAssembly.compile(threadTrapModule),
      memory: new WebAssembly.Memory({ initial: 1, maximum: 1, shared: true }),
      args: [], env: [], wasiShimUrl, faultChannel
    }
  });
  assert.deepEqual(posted, [{ id: 'prepare-1', type: 'prepare', ok: true }]);

  const control = new SharedArrayBuffer(4);
  await onmessage({ data: { type: 'start', id: 'start-3', tid: 3, startArg: 0, control } });
  const report = await settleWithin(fault, 5000, 'thread fault report');
  assert.equal(report.type, 'thread-fault');
  assert.equal(report.tid, 3);
  assert.match(report.error, /unreachable/);
  assert.equal(report.trap, true);
  assert.equal(Atomics.load(new Int32Array(control), 0), -1);
});

const setupPage = async (t, { trappedAttempts = Infinity, trap = true } = {}) => {
  assert.ok(WebAssembly.validate(threadedJobModule));
  restoreGlobalsAfter(t, ['window', 'crossOriginIsolated', 'fetch', 'Worker']);
  for (const level of ['info', 'warn']) {
    const previous = console[level];
    console[level] = () => {};
    t.after(() => { console[level] = previous; });
  }
  const workers = [];
  globalThis.window = { location: { href: 'https://example.test/gbdraw/web/' } };
  globalThis.crossOriginIsolated = true;
  globalThis.fetch = async () => ({ ok: true, arrayBuffer: async () => threadedJobModule.slice().buffer });
  // Real orchestration; the Worker transport imitates threaded jobs. The first
  // trappedAttempts jobs have a thread that trapped: the job worker never
  // answers, and the thread reports on the fault channel named in the run
  // message. Later threaded jobs answer normally.
  let threadedAttempts = 0;
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
      if (this.threaded && (threadedAttempts += 1) <= trappedAttempts) {
        const channel = new BroadcastChannel(message.faultChannel);
        channel.postMessage({ type: 'thread-fault', tid: 2, error: 'memory access out of bounds', trap, stderr: '' });
        channel.close();
        return;
      }
      queueMicrotask(() => {
        const data = { id: message.id, type: message.type, ok: true, text: 'q\ts\t100\n' };
        for (const listener of this.listeners.get('message') || []) listener({ data });
      });
    }
  };
  const owner = await import(`${losatUrl.href}?instance=${Date.now()}-${Math.random()}`);
  return { owner, workers };
};

const job = { pairIndex: 0, program: 'blastp', querySequenceKey: 'q', subjectSequenceKey: 's' };
const options = {
  totalThreadBudget: 2, concurrency: 1, threadsPerJob: 2,
  // Large enough for Auto mode to choose threaded execution.
  sequences: { q: '>q\n' + 'M'.repeat(300000), s: '>s\n' + 'M'.repeat(300000) }
};

test('explicit threaded LOSAT runs a trapped job again and keeps its result', async (t) => {
  const { owner, workers } = await setupPage(t, { trappedAttempts: 1 });
  const result = await settleWithin(
    owner.runLosatPairsParallel([job], { ...options, executionMode: 'threaded' }),
    5000,
    'threaded LOSAT run'
  );
  assert.deepEqual(result, [{ ...job, text: 'q\ts\t100\n' }]);
  assert.deepEqual(workers.map((worker) => worker.threaded), [true, true]);
  assert.ok(workers.every((worker) => worker.terminated));
});

test('explicit threaded LOSAT fails the job after its bounded trap retries', async (t) => {
  const { owner, workers } = await setupPage(t);
  const run = owner.runLosatPairsParallel([job], { ...options, executionMode: 'threaded' });
  await assert.rejects(
    settleWithin(run, 5000, 'threaded LOSAT run'),
    { message: 'LOSAT pair #1: LOSAT thread 2 failed: memory access out of bounds' }
  );
  // One attempt and two retries.
  assert.deepEqual(workers.map((worker) => worker.threaded), [true, true, true]);
  assert.ok(workers.every((worker) => worker.terminated));
  assert.match(workers[0].messages[0].faultChannel, /^gbdraw-losat-thread-fault:/);
});

test('explicit threaded LOSAT does not retry a thread failure that is not a trap', async (t) => {
  const { owner, workers } = await setupPage(t, { trap: false });
  const run = owner.runLosatPairsParallel([job], { ...options, executionMode: 'threaded' });
  await assert.rejects(
    settleWithin(run, 5000, 'threaded LOSAT run'),
    { message: 'LOSAT pair #1: LOSAT thread 2 failed: memory access out of bounds' }
  );
  assert.deepEqual(workers.map((worker) => worker.threaded), [true]);
});

test('automatic threaded LOSAT falls back to serial execution after its trap retries', async (t) => {
  const { owner, workers } = await setupPage(t);
  const statuses = [];
  const result = await settleWithin(
    owner.runLosatPairsParallel([job], { ...options, onRuntimeStatus: (status) => statuses.push(status) }),
    5000,
    'automatic LOSAT run'
  );
  assert.deepEqual(result, [{ ...job, text: 'q\ts\t100\n' }]);
  assert.deepEqual(statuses.map((status) => status.state), ['running', 'fallback']);
  assert.deepEqual(workers.map((worker) => worker.threaded), [true, true, true, false]);
  assert.ok(workers.every((worker) => worker.terminated));
});
