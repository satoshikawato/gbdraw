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

const randomBytes = (size, seed) => {
  let state = seed >>> 0;
  return Uint8Array.from({ length: size }, () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    return state >>> 24;
  });
};
const ENCODER_SIZES = [0, 1, 2, 3, 4, 5, 24575, 24576, 24577, 65535, 65536, 65537, 1048577];

const withoutNativeToBase64 = (fn) => {
  const native = Object.getOwnPropertyDescriptor(Uint8Array.prototype, 'toBase64');
  if (native) Object.defineProperty(Uint8Array.prototype, 'toBase64', { ...native, value: undefined });
  try { return fn(); } finally {
    if (native) Object.defineProperty(Uint8Array.prototype, 'toBase64', native);
  }
};

test('base64 encoding uses the native encoder when the engine has one', () => {
  const native = Uint8Array.prototype.toBase64;
  let calls = 0;
  Uint8Array.prototype.toBase64 = function (...options) {
    calls += 1;
    return typeof native === 'function'
      ? native.apply(this, options)
      : Buffer.from(this).toString('base64');
  };
  let encoded;
  try { encoded = bytesToBase64(randomBytes(70000, 7)); }
  finally {
    if (native) Uint8Array.prototype.toBase64 = native;
    else delete Uint8Array.prototype.toBase64;
  }
  assert.equal(calls, 1);
  assert.equal(encoded, Buffer.from(randomBytes(70000, 7)).toString('base64'));
});

test('the native and portable base64 encoders give identical text', () => {
  for (const size of ENCODER_SIZES) {
    const bytes = randomBytes(size, size + 1);
    const expected = Buffer.from(bytes).toString('base64');
    assert.equal(bytesToBase64(bytes), expected, `native size ${size}`);
    assert.equal(withoutNativeToBase64(() => bytesToBase64(bytes)), expected, `portable size ${size}`);
    assert.equal(withoutNativeToBase64(() => bytesToBase64(bytes.buffer)), expected, `portable buffer ${size}`);
  }
});

test('the portable base64 encoder builds no binary string', () => {
  const originalBtoa = globalThis.btoa;
  let calls = 0;
  globalThis.btoa = (value) => { calls += 1; return originalBtoa(value); };
  const bytes = new Uint8Array(8 * 1024 * 1024 + 2).fill(255);
  let encoded;
  try { encoded = withoutNativeToBase64(() => bytesToBase64(bytes)); }
  finally { globalThis.btoa = originalBtoa; }
  assert.equal(encoded, Buffer.from(bytes).toString('base64'));
  assert.equal(calls, 0);
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

// OV-1566: `sha256Hex` is the one SHA-256-to-hex owner of the Web modules; the
// font identity, the Worker's artifact and font identities, and the Pyodide
// snapshot key call it instead of hex-encoding their own digests.
test('sha256Hex is the only Web module that hex-encodes a SHA-256 digest', async () => {
  const { readdir, readFile } = await import('node:fs/promises');
  const webRoot = join(process.cwd(), 'gbdraw/web/js');
  const files = (await readdir(webRoot, { recursive: true })).filter((path) => path.endsWith('.js'));
  const digesting = [];
  for (const path of files) {
    if ((await readFile(join(webRoot, path), 'utf8')).includes('subtle.digest(')) digesting.push(path);
  }
  assert.deepEqual(digesting.sort(), ['services/byte-utils.js']);
  const { createHash } = await import('node:crypto');
  const bytes = Uint8Array.from({ length: 300 }, (_, index) => index % 256);
  const expected = createHash('sha256').update(bytes.subarray(7, 200)).digest('hex');
  // A view of part of a buffer digests its own bytes, as the Worker's font and artifact readers pass them.
  assert.equal(await codec.sha256Hex(bytes.subarray(7, 200)), expected);
  assert.equal(await codec.sha256Hex(bytes.slice(7, 200).buffer), expected);
});
