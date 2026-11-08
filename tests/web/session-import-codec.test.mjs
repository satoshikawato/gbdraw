import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { importSessionFile } from '../../gbdraw/web/js/services/session-import-client.js';
import { gzipSync } from 'node:zlib';
import { readSessionText } from '../../gbdraw/web/js/services/session-file.js';
import { assertSafeObjectKeys, isUnsafeObjectKey } from '../../gbdraw/web/js/services/safe-object-keys.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';
import { writeSessionContainer } from '../../gbdraw/web/js/services/session-container.js';

installSessionImportWorker();
const decode = async blob => {
  try { return { status: 'ok', ...await importSessionFile(blob) }; }
  catch (error) { return { status: 'error', error }; }
};

test('production Worker reads plain/gzip magic and preserves whole JSON semantics', async () => {
  const data = { unicode: '日本語 😀', surrogate: '\ud800', rows: [null, true, 1.25], nested: { key: 'kept' } };
  const text = JSON.stringify(data);
  for (const file of [new Blob([text]), new Blob([gzipSync(text)], { type: 'text/plain' })]) {
    const reply = await decode(file);
    assert.equal(reply.status, 'ok');
    assert.equal(reply.characters, text.length);
    assert.deepEqual(reply.data, data);
    assert.ok(reply.timings.readMs >= 0 && reply.timings.parseMs >= 0);
  }
});

test('malformed JSON/gzip and fatal UTF-8 stay structured failures', async () => {
  for (const [bytes, stage] of [
    ['{private broken data', 'parse'],
    [new Uint8Array([0x1f, 0x8b, 1, 2, 3]), 'read'],
    [new Uint8Array([0xff]), 'read'],
    [gzipSync(new Uint8Array([0xc3])), 'read']
  ]) {
    const reply = await decode(new Blob([bytes]));
    assert.equal(reply.status, 'error');
    assert.equal(reply.error.stage, stage);
    // X-01 (SE-09): the Worker's final diagnostic code reaches the normalizer.
    assert.equal(reply.error.code, stage === 'parse' ? 'INPUT_INVALID' : 'INPUT_UNREADABLE');
    assert.equal(normalizeUserFacingError(reply.error).code, reply.error.code);
    assert.ok(!JSON.stringify(normalizeUserFacingError(reply.error)).includes('private broken data'));
    assert.ok(!reply.error.message.includes('private broken data'));
  }
});

// Receipt session.file-container: content decides. The Worker reads a container's directory and document and
// replies with the document and the entry of each declared resource (the page streams them into the copies);
// any other file stays a JSON Session.
test('the production Worker reads a Session file container by its content and replies with its entries', async () => {
  const document = { format: 'gbdraw-session', version: 51, resources: { 'resource-0001': { kind: 'web-file', name: 'a.gb', type: 'text/plain', size: 5, lastModified: 0 } } };
  const blob = await writeSessionContainer({
    document,
    documentBytes: () => new Blob([JSON.stringify(document)]).stream(),
    readResource: () => new Blob(['ACGT\n']).stream()
  });
  const reply = await decode(new File([blob], 'named.json'));
  assert.equal(reply.status, 'ok');
  assert.deepEqual(reply.data, document);
  assert.deepEqual(Object.keys(reply.container.resources), ['resource-0001']);
  assert.equal(reply.container.resources['resource-0001'].size, 5);
  assert.equal(reply.container.inputBytes, 5);
  const json = await decode(new File([JSON.stringify(document)], 'named.gbdraw'));
  assert.equal(json.status, 'ok');
  assert.equal(json.container, null);
});

// OV-2694 (D-30): a JSON Session with one key twice is refused as not JSON, as Python refuses it.
test('a JSON Session that holds one key twice is refused as not valid JSON', async () => {
  for (const text of ['{"version": 46, "version": 46}', '{"a": {"b": 1, "\\u0062": 2}}']) {
    for (const bytes of [text, gzipSync(text)]) {
      const reply = await decode(new Blob([bytes]));
      assert.equal(reply.status, 'error');
      assert.deepEqual([reply.error.code, reply.error.stage, reply.error.context.reason], ['INPUT_INVALID', 'parse', 'JSON_FORMAT']);
      assert.equal(normalizeUserFacingError(reply.error, { operation: 'session-load' }).summary,
        'This file is not a gbdraw Session: it is not valid JSON.');
    }
  }
  assert.equal((await decode(new Blob(['{"a": [{"b": 1}, {"b": 2}], "c": "{\\"a\\": 1}"}']))).status, 'ok');
});

// The Option Integrity Product Contract cites this title; the production Worker now refuses the key itself.
test('unsafe own keys remain untrusted after the production whole-object reply', async () => {
  for (const key of ['__proto__', 'constructor', 'prototype']) {
    const reply = await decode(new Blob([`{"nested":[{"kept":1},{"${key}":{}}]}`]));
    assert.equal(reply.status, 'error');
    // OV-952: the diagnostic the main thread gave before the walk moved into the Worker, never the key.
    assert.equal(reply.error.message, 'INPUT_INVALID/FIELDS');
    assert.equal(reply.error.code, 'INPUT_INVALID');
    assert.equal(reply.error.stage, 'request-validation');
    assert.deepEqual(reply.error.context, { field: 'schema', reason: 'FIELDS' });
    assert.ok(!JSON.stringify(normalizeUserFacingError(reply.error)).includes(key));
  }
  // The same walk still serves the routes that check a value on the main thread.
  assert.throws(() => assertSafeObjectKeys(JSON.parse('{"a":[{"constructor":1}]}')), (error) => (
    error.code === 'INPUT_INVALID' && error.context.field === 'schema' && error.context.reason === 'FIELDS'
  ));
  assert.doesNotThrow(() => assertSafeObjectKeys({ a: [{ b: 'constructor' }], proto: 1 }));
});

test('codec rejects exact file and expanded caps before reading/decoding oversize data', async () => {
  await assert.rejects(readSessionText({ size: 200 * 1024 * 1024 + 1,
    stream() { throw new Error('must not read'); } }), { code: 'SESSION_SIZE_LIMIT', stage: 'read' });
  const native = globalThis.DecompressionStream;
  globalThis.DecompressionStream = class {
    constructor() {
      const stream = new TransformStream({ transform(_chunk, controller) {
        controller.enqueue({ byteLength: 512 * 1024 * 1024 + 1 });
      } });
      this.readable = stream.readable; this.writable = stream.writable;
    }
  };
  try {
    await assert.rejects(readSessionText(new Blob([new Uint8Array([0x1f, 0x8b])])), { code: 'SESSION_SIZE_LIMIT', stage: 'read' });
  } finally { globalThis.DecompressionStream = native; }
});

// services/safe-object-keys.js is the one list of the keys that name an
// object's prototype or constructor: Load refuses them, and a merge skips them.
test('the unsafe object keys are listed once, in services/safe-object-keys.js', () => {
  assert.deepEqual(['__proto__', 'constructor', 'prototype', 'toString', 'proto'].map(isUnsafeObjectKey),
    [true, true, true, false, false]);
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const lists = readdirSync(root, { recursive: true }).map(String).filter((path) => path.endsWith('.js')
    && /['"]__proto__['"]\s*,\s*['"]constructor['"]/.test(readFileSync(new URL(path, root), 'utf8')));
  assert.deepEqual(lists.map((path) => path.replaceAll('\\', '/')), ['services/safe-object-keys.js']);
});
