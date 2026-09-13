import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

globalThis.location = { href: 'https://example.test/gbdraw/web/' };
globalThis.self = {};
const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);
// Exercise the private Worker staging owner without adding a production export.
const workerUrl = new URL('../../gbdraw/web/js/workers/diagram-generation-worker.js', import.meta.url);
const workerSource = (await readFile(workerUrl, 'utf8')).replace(
  /from '([^']+)'/g, (_match, path) => `from '${new URL(path, workerUrl)}'`
);
let workerNumber = 0;
const loadStagingOwner = () => import(`data:text/javascript;base64,${Buffer.from(
  `${workerSource}\nexport { stageRenderResources };\n// instance ${++workerNumber}`
).toString('base64')}`);

class MemoryFS {
  entries = new Map();
  analyzePath(path) { return { exists: this.entries.has(path) }; }
  mkdir(path) { this.entries.set(path, null); }
  writeFile(path, bytes) { this.entries.set(path, bytes); }
  createDataFile(directory, name, bytes) { this.writeFile(`${directory}/${name}`, bytes); }
  unlink(path) { assert.ok(this.entries.delete(path), path); }
  symlink(target, path) { this.entries.set(path, target); }
  read(path) {
    const value = this.entries.get(path);
    return typeof value === 'string' ? this.read(value) : value;
  }
  cleanWorkspace(workspace) {
    for (const path of this.entries.keys()) if (path.startsWith(`${workspace}/`) || path === workspace) this.entries.delete(path);
  }
  cacheFiles() { return [...this.entries].filter(([path, bytes]) => path.endsWith('.bin') && bytes instanceof Uint8Array); }
}

class StagingWorker {
  static instances = [];
  listeners = new Map();
  messages = [];
  fs = new MemoryFS();
  owner = loadStagingOwner();
  mode = 'success';
  constructor() { StagingWorker.instances.push(this); }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, new Set());
    this.listeners.get(type).add(listener);
  }
  removeEventListener(type, listener) { this.listeners.get(type)?.delete(listener); }
  emit(data) { for (const listener of [...(this.listeners.get('message') || [])]) listener({ data }); }
  terminate() { this.terminated = true; this.fs.entries.clear(); }
  postMessage(message) {
    const copy = structuredClone(message);
    this.messages.push(copy);
    if (message.type === 'init') {
      queueMicrotask(() => this.emit({ type: 'init', id: message.id, ok: true, capabilities: EXPECTED_WEB_RUNTIME_CAPABILITIES }));
    } else if (['helper', 'feature-extraction'].includes(message.type)) {
      queueMicrotask(() => this.emit({ type: message.type, requestId: message.requestId, ok: false, error: { message: 'Invalid helper source' } }));
    } else {
      this.operation = this.run(structuredClone(message));
    }
  }
  async run({ requestId, payload }) {
    const workspace = `/request-${requestId}`;
    try {
      const { stageRenderResources } = await this.owner;
      const paths = await stageRenderResources({ FS: this.fs }, workspace, payload.resourceManifest, payload.stagedResources, false, requestId);
      const content = Object.fromEntries(Object.entries(paths).map(([id, path]) => [id, Buffer.from(this.fs.read(path)).toString()]));
      this.fs.cleanWorkspace(workspace);
      this.staged?.();
      if (this.mode === 'hold') return;
      if (this.mode === 'protocol-error') throw new Error('Worker protocol failure after staging');
      this.emit({ type: 'run', requestId, ok: true, results: this.mode === 'python-error'
        ? { error: { type: 'ParseError', message: 'Rejected biological source' } }
        : [{ name: 'out.svg', content: JSON.stringify(content) }] });
    } catch (error) {
      this.fs.cleanWorkspace(workspace);
      this.emit({ type: 'run', requestId, ok: false, error: { message: error.message } });
    }
  }
  latestRun() { return this.messages.filter(m => m.type === 'run').at(-1); }
}
globalThis.Worker = StagingWorker;
const { runDiagramGeneration, runDiagramHelperOperation, runFeatureExtraction, cancelDiagramGeneration, disposeDiagramGenerationWorker } = await import(
  '../../gbdraw/web/js/services/diagram-generation.js'
);
const payload = values => ({
  request: { records: Object.keys(values).map(resourceId => ({ source: { kind: 'genbank', resourceId } })) },
  resources: Object.fromEntries(Object.entries(values).map(([id, value]) => [id, {
    name: `${id}.gb`, kind: 'genbank', size: Buffer.byteLength(value), encoding: 'base64', data: Buffer.from(value).toString('base64')
  }]))
});
const run = values => runDiagramGeneration(payload(values));
const transferred = worker => worker.latestRun().payload.stagedResources.map(r => r.resourceId);
const content = response => JSON.parse(response.results[0].content);
const current = () => StagingWorker.instances.at(-1);
test.afterEach(() => { disposeDiagramGenerationWorker(); delete globalThis.__GBDRAW_TEST_HOOKS__; });

test('W2/W4/W5: structured Python rejection acknowledges changed bytes, and A recovers once', async () => {
  const a = await run({ source: 'AAAA' });
  const worker = current();
  assert.deepEqual(transferred(worker), ['source']);
  assert.deepEqual(content(await run({ source: 'AAAA' })), content(a));
  assert.deepEqual(transferred(worker), []);
  const aToken = worker.latestRun().payload.resourceManifest[0].cacheToken;
  worker.mode = 'python-error';
  assert.ok((await run({ source: 'BBBB' })).results.error);
  assert.deepEqual(transferred(worker), ['source']);
  assert.notEqual(worker.latestRun().payload.resourceManifest[0].cacheToken, aToken);
  assert.equal(worker.terminated, undefined);
  await run({ source: 'BBBB' });
  assert.deepEqual(transferred(worker), [], 'rejected render still acknowledged raw resource ownership');
  worker.mode = 'success';
  assert.deepEqual(content(await run({ source: 'AAAA' })), content(a));
  assert.deepEqual(transferred(worker), ['source']);
  assert.equal(current(), worker);
  assert.equal(worker.fs.cacheFiles().length, 1);
});

test('W3: accepted A-B-A uses exact bytes with replacement under the same ID', async () => {
  for (const text of ['AAAA', 'BBBB', 'AAAA']) {
    assert.deepEqual(content(await run({ source: text })), { source: text });
    assert.deepEqual(transferred(current()), ['source']);
  }
});

test('W6: rejected A+C prunes B; restored A+B transfers B and reuses only A', async () => {
  const a = await run({ shared: 'A', original: 'B' });
  const worker = current();
  worker.mode = 'python-error';
  await run({ shared: 'A', candidate: 'C' });
  assert.deepEqual(transferred(worker), ['candidate']);
  assert.deepEqual(worker.fs.cacheFiles().map(([, bytes]) => Buffer.from(bytes).toString()).sort(), ['A', 'C']);
  worker.mode = 'success';
  assert.deepEqual(content(await run({ shared: 'A', original: 'B' })), content(a));
  assert.deepEqual(transferred(worker), ['original']);
  assert.deepEqual(worker.fs.cacheFiles().map(([, bytes]) => Buffer.from(bytes).toString()).sort(), ['A', 'B']);
});

test('W7: protocol failure after staging terminates the Worker and resets main knowledge', async () => {
  await run({ source: 'A' });
  const worker = current();
  worker.mode = 'protocol-error';
  await assert.rejects(run({ source: 'B' }), /Worker protocol failure/);
  assert.equal(worker.terminated, true);
  assert.equal(worker.fs.entries.size, 0);
  assert.deepEqual(content(await run({ source: 'A' })), { source: 'A' });
  assert.notEqual(current(), worker);
  assert.deepEqual(transferred(current()), ['source']);
});

test('partial staging failure releases replaced files and cannot advertise the old cache', async () => {
  await run({ first: 'A', second: 'B' });
  const worker = current();
  const create = worker.fs.createDataFile.bind(worker.fs);
  let writes = 0;
  worker.fs.createDataFile = (...args) => {
    if (++writes === 2) throw new Error('Injected cache file creation failure');
    create(...args);
  };
  await assert.rejects(run({ first: 'C', second: 'D' }), /cache file creation failure/);
  assert.equal(writes, 2);
  assert.equal(worker.terminated, true);
  assert.equal(worker.fs.entries.size, 0);
  assert.deepEqual(content(await run({ first: 'A', second: 'B' })), { first: 'A', second: 'B' });
  assert.deepEqual(transferred(current()), ['first', 'second']);
});

test('pre-dispatch validation and auxiliary input errors leave the render resource cache unchanged', async () => {
  await run({ source: 'A' });
  const worker = current();
  const invalid = payload({ source: 'B' });
  invalid.resources.source.size = -1;
  await assert.rejects(runDiagramGeneration(invalid), /invalid byte size/);
  assert.equal(worker.messages.filter(m => m.type === 'run').length, 1);
  await assert.rejects(runDiagramHelperOperation('listSequenceRecords', {}), /Invalid helper source/);
  await assert.rejects(runFeatureExtraction({}), /Invalid helper source/);
  assert.deepEqual(content(await run({ source: 'A' })), { source: 'A' });
  assert.deepEqual(transferred(worker), []);
  assert.equal(current(), worker);
  assert.equal(worker.terminated, undefined);
});

test('W8: cancel after staging releases Worker ownership and transfers A to a new Worker', async () => {
  await run({ source: 'A' });
  const worker = current();
  worker.mode = 'hold';
  const staged = new Promise(resolve => { worker.staged = resolve; });
  const pending = run({ source: 'B' });
  const rejected = assert.rejects(pending, /canceled/);
  await staged;
  cancelDiagramGeneration();
  await rejected;
  assert.equal(worker.terminated, true);
  assert.deepEqual(content(await run({ source: 'A' })), { source: 'A' });
  assert.deepEqual(transferred(current()), ['source']);
});

test('W9: unadmitted B has no delayed promotion; A recovers with its bytes', async () => {
  await run({ source: 'A' });
  const worker = current();
  const candidate = await run({ source: 'B' });
  assert.equal('finalizeResourcePromotion' in candidate, false);
  // Caller rejects the Result: no artifact publication or cache callback.
  assert.deepEqual(content(await run({ source: 'A' })), { source: 'A' });
  assert.deepEqual(transferred(worker), ['source']);
  assert.equal(current(), worker);
});

test('cancel during response handling cannot republish disposed Worker cache knowledge', async () => {
  await run({ source: 'A' });
  let release;
  let entered;
  const held = new Promise(resolve => { release = resolve; });
  const response = new Promise(resolve => { entered = resolve; });
  globalThis.__GBDRAW_TEST_HOOKS__ = { beforeDiagramGenerationResponse: () => { entered(); return held; } };
  const pending = run({ source: 'B' });
  const rejected = assert.rejects(pending, /canceled/);
  await response;
  cancelDiagramGeneration();
  await rejected;
  delete globalThis.__GBDRAW_TEST_HOOKS__;
  await run({ source: 'A' });
  assert.deepEqual(transferred(current()), ['source']);
  release();
  await Promise.resolve();
  await run({ source: 'A' });
  assert.deepEqual(transferred(current()), []);
});
