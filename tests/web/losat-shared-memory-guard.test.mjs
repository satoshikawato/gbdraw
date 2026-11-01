import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { guardSharedMemoryBulkOps } from '../../gbdraw/web/js/services/losat-shared-memory-guard.js';

// In Chromium, memory.copy and memory.fill in a LOSAT thread trap "memory
// access out of bounds" when another thread has just grown the shared memory
// (the worker's view is older than the growth). The guard routes them through
// a helper that refreshes the view with memory.grow(0) first. A stale view is
// a race between isolates and is not reproduced here: these tests pin what
// the guard writes, that guarded operations behave as before, which modules it
// guards, and that compileLosatWasm compiles what it returns.
const webRoot = new URL('../../gbdraw/web/', import.meta.url);

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
const body = (...code) => [...leb(code.length + 1), 0x00, ...code];
const args3 = [0x20, 0, 0x20, 1, 0x20, 2];

// copy(dst, src, len) and fill(dst, value, len) on an imported memory.
const bulkModule = (memoryLimits) => moduleBytes(
  section(1, [[0x60, 3, I32, I32, I32, 0]]),
  section(2, [[...name('env'), ...name('memory'), 0x02, ...memoryLimits]]),
  section(3, [[0], [0]]),
  section(7, [[...name('copy'), 0x00, 0], [...name('fill'), 0x00, 1]]),
  section(10, [body(...args3, 0xfc, 0x0a, 0x00, 0x00, 0x0b), body(...args3, 0xfc, 0x0b, 0x00, 0x0b)])
);
const sharedBulkModule = bulkModule([0x03, 1, 2]);

test('the guard routes each bulk operation through a helper that behaves as the operation', async () => {
  const result = guardSharedMemoryBulkOps(sharedBulkModule);
  assert.ok(result);
  assert.equal(result.copy, 1);
  assert.equal(result.fill, 1);
  assert.ok(WebAssembly.validate(result.bytes));

  for (const bytes of [sharedBulkModule, result.bytes]) {
    const memory = new WebAssembly.Memory({ initial: 1, maximum: 2, shared: true });
    const { instance } = await WebAssembly.instantiate(bytes, { env: { memory } });
    const view = new Uint8Array(memory.buffer);
    instance.exports.fill(10, 7, 4);
    instance.exports.copy(20, 9, 4);
    assert.deepEqual([...view.subarray(9, 15)], [0, 7, 7, 7, 7, 0]);
    assert.deepEqual([...view.subarray(19, 25)], [0, 0, 7, 7, 7, 0]);
    // A range that lies outside the memory still traps, and grow(0) allocates nothing.
    assert.throws(() => instance.exports.copy(65530, 0, 16), WebAssembly.RuntimeError);
    assert.throws(() => instance.exports.fill(0, 1, 65537), WebAssembly.RuntimeError);
    assert.equal(memory.buffer.byteLength, 65536);
  }
});

test('the guard leaves a module without a shared memory or without bulk operations as it is', () => {
  const unshared = bulkModule([0x01, 1, 2]);
  assert.deepEqual(guardSharedMemoryBulkOps(unshared), { bytes: unshared, copy: 0, fill: 0 });
  const noBulk = moduleBytes(
    section(1, [[0x60, 0, 0]]),
    section(2, [[...name('env'), ...name('memory'), 0x02, 0x03, 1, 2]]),
    section(3, [[0]]),
    section(10, [body(0x0b)])
  );
  assert.deepEqual(guardSharedMemoryBulkOps(noBulk), { bytes: noBulk, copy: 0, fill: 0 });
});

// A shared-memory module with an instruction the reader does not know (an
// exception `throw`, with its tag section).
const unknownInstructionModule = moduleBytes(
  section(1, [[0x60, 0, 0]]),
  section(2, [[...name('env'), ...name('memory'), 0x02, 0x03, 1, 2]]),
  section(3, [[0]]),
  section(13, [[0x00, 0]]),
  section(10, [body(0x08, 0x00, 0x0b)])
);

test('the guard declines a module with an instruction it does not read', () => {
  assert.ok(WebAssembly.validate(unknownInstructionModule), 'the engine accepts the module');
  assert.equal(guardSharedMemoryBulkOps(unknownInstructionModule), null);
});

const stubCompile = (t, bytes) => {
  const original = { fetch: globalThis.fetch, compile: WebAssembly.compile, warn: console.warn };
  const compiled = [];
  const warnings = [];
  globalThis.fetch = async () => ({ ok: true, arrayBuffer: async () => bytes.slice().buffer });
  WebAssembly.compile = (input) => {
    compiled.push(new Uint8Array(input));
    return original.compile(input);
  };
  console.warn = (message) => warnings.push(String(message));
  t.after(() => {
    globalThis.fetch = original.fetch;
    WebAssembly.compile = original.compile;
    console.warn = original.warn;
  });
  return { compiled, warnings };
};

test('compileLosatWasm compiles the guarded module, and a declined module as it is', async (t) => {
  const { compileLosatWasm } = await import('../../gbdraw/web/js/services/losat-runtime.js');
  const { compiled, warnings } = stubCompile(t, sharedBulkModule);
  await compileLosatWasm('https://example.test/losat-threaded.wasm', 'threaded LOSAT wasm');
  assert.deepEqual(compiled, [guardSharedMemoryBulkOps(sharedBulkModule).bytes]);
  assert.deepEqual(warnings, []);

  const declined = stubCompile(t, unknownInstructionModule);
  await compileLosatWasm('https://example.test/other.wasm', 'threaded LOSAT wasm');
  assert.deepEqual(declined.compiled, [unknownInstructionModule]);
  assert.match(declined.warnings[0], /could not be guarded; it runs unguarded/);
});

test('the packaged threaded LOSAT build is guarded and the serial build is left as it is', () => {
  const threaded = new Uint8Array(readFileSync(new URL('wasm/losat/losat-threaded.wasm', webRoot)));
  const result = guardSharedMemoryBulkOps(threaded);
  assert.ok(result, 'the guard reads every instruction of the packaged threaded build');
  assert.ok(result.copy > 0, 'the threaded build has memory.copy operations to guard');
  assert.ok(WebAssembly.validate(result.bytes));
  const serial = new Uint8Array(readFileSync(new URL('wasm/losat/losat.wasm', webRoot)));
  assert.equal(guardSharedMemoryBulkOps(serial)?.bytes, serial);
});
