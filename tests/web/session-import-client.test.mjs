import assert from 'node:assert/strict';
import test from 'node:test';
import { importSessionFile } from '../../gbdraw/web/js/services/session-import-client.js';

class ImportWorker {
  static instances = [];
  constructor(url, options) {
    this.url = String(url);
    this.options = options;
    this.listeners = new Map();
    this.terminations = 0;
    ImportWorker.instances.push(this);
  }
  addEventListener(name, fn) { this.listeners.set(name, fn); }
  removeEventListener(name, fn) { assert.equal(this.listeners.get(name), fn); this.listeners.delete(name); }
  postMessage(message) { this.message = message; }
  terminate() { this.terminations += 1; }
  reply(data) { this.listeners.get('message')?.({ data }); }
}
const blob = new Blob(['{"private":"candidate"}']);
const clean = worker => {
  assert.equal(worker.terminations, 1);
  assert.equal(worker.listeners.size, 0);
};

test('whole-object success terminates before candidate use and ignores stale/duplicate replies', async () => {
  globalThis.Worker = ImportWorker;
  const pending = importSessionFile(blob);
  const worker = ImportWorker.instances.at(-1);
  assert.match(worker.url, /workers\/session-import-worker.js$/);
  assert.equal(worker.options.type, 'module');
  assert.equal(worker.message.file, blob, 'File/Blob itself reaches Worker; no main text');
  const operationId = worker.message.operationId;
  const receive = worker.listeners.get('message');
  worker.reply({ operationId: operationId - 1, status: 'ok', data: {} });
  assert.equal(worker.terminations, 0);
  const data = { nested: ['kept', null, 1] };
  worker.reply({ operationId, status: 'ok', data, characters: 42 });
  clean(worker);
  receive({ data: { operationId, status: 'error', error: { message: 'late' } } });
  const candidate = await pending;
  assert.equal(candidate.data, data, 'transport does not make an extra clone');
  clean(worker);
});

test('structured errors, crashes, unreadable and malformed replies settle and release listeners', async () => {
  globalThis.Worker = ImportWorker;
  for (const kind of ['structured', 'crash', 'messageerror', 'malformed']) {
    const pending = importSessionFile(blob);
    const worker = ImportWorker.instances.at(-1);
    if (kind === 'structured') {
      worker.reply({ operationId: worker.message.operationId, status: 'error',
        error: { name: 'SyntaxError', code: 'SESSION_IMPORT_PARSE_FAILED', stage: 'parse', message: 'Invalid JSON structure.' } });
    } else if (kind === 'crash') worker.listeners.get('error')({ preventDefault() {} });
    else if (kind === 'messageerror') worker.listeners.get('messageerror')();
    else worker.reply({ operationId: worker.message.operationId, status: 'ok' });
    await assert.rejects(pending, error => {
      assert.match(error.code, /^SESSION_IMPORT_/);
      if (kind === 'structured') {
        assert.equal(error.stage, 'parse');
        assert.equal(error.name, 'SyntaxError');
      }
      return true;
    });
    clean(worker);
  }
});

test('cancel/teardown terminates immediately, rejects late reply, and removes the signal listener', async () => {
  globalThis.Worker = ImportWorker;
  for (let i = 0; i < 20; i += 1) {
    const controller = new AbortController();
    let added = 0, removed = 0;
    const add = controller.signal.addEventListener.bind(controller.signal);
    const remove = controller.signal.removeEventListener.bind(controller.signal);
    controller.signal.addEventListener = (...args) => { added += 1; return add(...args); };
    controller.signal.removeEventListener = (...args) => { removed += 1; return remove(...args); };
    const pending = importSessionFile(blob, { signal: controller.signal });
    const worker = ImportWorker.instances.at(-1);
    const receive = worker.listeners.get('message');
    controller.abort();
    controller.abort();
    receive({ data: { operationId: worker.message.operationId, status: 'ok', data: {}, characters: 0 } });
    await assert.rejects(pending, { name: 'AbortError', code: 'SESSION_IMPORT_CANCELED' });
    clean(worker);
    assert.equal(added, 1);
    assert.equal(removed, 1);
  }
  const controller = new AbortController();
  controller.abort();
  const count = ImportWorker.instances.length;
  await assert.rejects(importSessionFile(blob, { signal: controller.signal }), { name: 'AbortError' });
  assert.equal(ImportWorker.instances.length, count);
});

test('Worker unavailable, constructor and postMessage failures have no main parser fallback', async () => {
  delete globalThis.Worker;
  await assert.rejects(importSessionFile(blob), { code: 'SESSION_IMPORT_UNAVAILABLE' });
  globalThis.Worker = class { constructor() { throw new Error('blocked'); } };
  await assert.rejects(importSessionFile(blob), { code: 'SESSION_IMPORT_START_FAILED' });
  globalThis.Worker = class extends ImportWorker { postMessage() { throw new Error('clone failed'); } };
  await assert.rejects(importSessionFile(blob), { code: 'SESSION_IMPORT_START_FAILED' });
  clean(ImportWorker.instances.at(-1));
  delete globalThis.Worker;
});
