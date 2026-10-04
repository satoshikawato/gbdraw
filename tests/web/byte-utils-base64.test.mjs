import assert from 'node:assert/strict';
import { cp, mkdtemp, rm, writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

const root = await mkdtemp(join(process.cwd(), '.gbdraw-base64-test-'));
let codec;
try {
  await cp(join(process.cwd(), 'gbdraw/web/js/services/byte-utils.js'), join(root, 'byte-utils.js'));
  await writeFile(join(root, 'package.json'), '{"type":"module"}');
  codec = await import(pathToFileURL(join(root, 'byte-utils.js')));
} finally {
  await rm(root, { recursive: true, force: true });
}
const { bytesToBase64, base64ToBytes, base64ToBytesInTasks, textToBase64 } = codec;

test('base64 preserves every byte and padding across chunk boundaries', () => {
  for (const size of [0, 1, 2, 3, 255, 256, 24575, 24576, 24577, 24578, 49151, 49152, 49153, 1048577]) {
    const bytes = Uint8Array.from({ length: size }, (_, index) => (index * 73 + 19) % 256);
    const encoded = bytesToBase64(bytes);
    assert.equal(encoded, Buffer.from(bytes).toString('base64'), `size ${size}`);
    assert.deepEqual(base64ToBytes(encoded), bytes, `round trip ${size}`);
  }
});

test('base64 bounds the temporary binary string while keeping the full resource', () => {
  const originalBtoa = globalThis.btoa;
  let maximumInput = 0;
  let calls = 0;
  globalThis.btoa = value => {
    maximumInput = Math.max(maximumInput, value.length);
    calls += 1;
    return originalBtoa(value);
  };
  const bytes = new Uint8Array(8 * 1024 * 1024 + 2).fill(255);
  let encoded;
  try { encoded = bytesToBase64(bytes); } finally { globalThis.btoa = originalBtoa; }
  assert.equal(encoded, Buffer.from(bytes).toString('base64'));
  assert(maximumInput <= 32768, 'Encoding must not materialize a full-resource binary string');
  assert(calls > 1);
});

test('text base64 retains native UTF-8 semantics for Unicode and lone surrogates', () => {
  const text = ('α界🧬\ud800\udfff\0\n').repeat(20000);
  assert.equal(textToBase64(text), Buffer.from(text, 'utf8').toString('base64'));
});


test('adopted-resource decoding writes bounded native chunks with exact bytes', { skip: typeof Uint8Array.prototype.setFromBase64 !== 'function' }, async () => {
  const bytes = Uint8Array.from({ length: 8 * 1024 * 1024 + 2 }, (_, index) => (index * 37 + 11) % 256);
  const encoded = Buffer.from(bytes).toString('base64');
  const nativeSet = Uint8Array.prototype.setFromBase64;
  let largestInput = 0;
  let calls = 0;
  Uint8Array.prototype.setFromBase64 = function (value, ...options) {
    largestInput = Math.max(largestInput, value.length);
    calls += 1;
    return nativeSet.call(this, value, ...options);
  };
  let decoded;
  try { decoded = await base64ToBytesInTasks(encoded); }
  finally { Uint8Array.prototype.setFromBase64 = nativeSet; }
  assert.deepEqual(decoded, bytes);
  assert(largestInput <= 0x10000);
  assert(calls > 1);
  assert.deepEqual(await base64ToBytesInTasks('T W\nE=\t'), Uint8Array.from(Buffer.from('Ma')));
  for (const invalid of ['A', 'TQ=', 'TQ==AAAA', 'TQ!']) {
    await assert.rejects(base64ToBytesInTasks(invalid));
  }
});

test('adopted-resource decoding preserves atob behavior on older engines', async () => {
  const nativeSet = Object.getOwnPropertyDescriptor(Uint8Array.prototype, 'setFromBase64');
  const originalAtob = globalThis.atob;
  let calls = 0;
  if (nativeSet) Object.defineProperty(Uint8Array.prototype, 'setFromBase64', { ...nativeSet, value: undefined });
  globalThis.atob = value => { calls += 1; return originalAtob(value); };
  try {
    const bytes = Uint8Array.from({ length: 65538 }, (_, index) => (index * 17 + 3) % 256);
    assert.deepEqual(await base64ToBytesInTasks(Buffer.from(bytes).toString('base64')), bytes);
    assert.deepEqual(await base64ToBytesInTasks('T W\nE=\t'), Uint8Array.from(Buffer.from('Ma')));
    await assert.rejects(base64ToBytesInTasks('TQ!'));
    assert(calls >= 3);
  } finally {
    if (nativeSet) Object.defineProperty(Uint8Array.prototype, 'setFromBase64', nativeSet);
    globalThis.atob = originalAtob;
  }
});
