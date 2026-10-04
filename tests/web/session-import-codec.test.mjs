import assert from 'node:assert/strict';
import test from 'node:test';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { importSessionFile } from '../../gbdraw/web/js/services/session-import-client.js';
import { gzipSync } from 'node:zlib';
import { readSessionText } from '../../gbdraw/web/js/services/session-file.js';
import { assertSafeObjectKeys } from '../../gbdraw/web/js/services/safe-object-keys.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/services/error-normalization.js';

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

test('unsafe own keys remain untrusted after the production whole-object reply', async () => {
  for (const key of ['__proto__', 'constructor', 'prototype']) {
    const reply = await decode(new Blob([`{"nested":{"${key}":{}}}`]));
    assert.equal(reply.status, 'ok');
    assert.equal(Object.hasOwn(reply.data.nested, key), true);
    assert.throws(() => assertSafeObjectKeys(reply.data, 'Session'), /unsafe key/);
  }
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
