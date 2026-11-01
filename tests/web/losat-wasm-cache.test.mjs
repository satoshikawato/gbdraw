// LOSAT wasm loads revalidate instead of downloading again: every loader
// fetches through `compileLosatWasm` (services/losat-runtime.js) with
// `cache: 'no-cache'`, so an unchanged file costs a 304 and a replaced file
// at the same URL is never served stale.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

const { compileLosatWasm } = await import('../../gbdraw/web/js/services/losat-runtime.js');

const LOADERS = [
  'gbdraw/web/js/services/losat.js',
  'gbdraw/web/js/workers/losat-worker.js',
  'gbdraw/web/js/workers/losat-threaded-worker.js',
  'gbdraw/web/js/workers/losat-wasi-thread-worker.js',
];

test('compileLosatWasm revalidates the wasm and compiles it', async () => {
  const calls = [];
  const original = globalThis.fetch;
  globalThis.fetch = async (url, init) => {
    calls.push([url, init]);
    return new Response(new Uint8Array([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0]), { status: 200 });
  };
  try {
    const module = await compileLosatWasm('https://example.test/wasm/losat/losat.wasm', 'LOSAT wasm');
    assert.ok(module instanceof WebAssembly.Module);
  } finally {
    globalThis.fetch = original;
  }
  assert.deepEqual(calls, [['https://example.test/wasm/losat/losat.wasm', { cache: 'no-cache' }]]);
});

test('compileLosatWasm compiles the bytes as served: LOSAT guards the threaded file when it builds it', async (t) => {
  // A shared memory and a memory.fill, the shape a host-side rewrite would change.
  const served = new Uint8Array([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    1, 7, 1, 0x60, 3, 0x7f, 0x7f, 0x7f, 0,
    2, 11, 1, 3, 0x65, 0x6e, 0x76, 1, 0x6d, 0x02, 0x03, 1, 2,
    3, 2, 1, 0,
    10, 13, 1, 11, 0, 0x20, 0, 0x20, 1, 0x20, 2, 0xfc, 0x0b, 0, 0x0b
  ]);
  assert.ok(WebAssembly.validate(served));
  const originalFetch = globalThis.fetch;
  const originalCompile = WebAssembly.compile;
  const compiled = [];
  globalThis.fetch = async () => new Response(served, { status: 200 });
  WebAssembly.compile = (bytes) => {
    compiled.push(new Uint8Array(bytes));
    return originalCompile(bytes);
  };
  t.after(() => {
    globalThis.fetch = originalFetch;
    WebAssembly.compile = originalCompile;
  });
  await compileLosatWasm('https://example.test/wasm/losat/losat-threaded.wasm', 'threaded LOSAT wasm');
  assert.deepEqual(compiled, [served]);
});

test('compileLosatWasm names a missing wasm', async () => {
  const original = globalThis.fetch;
  globalThis.fetch = async () => new Response('', { status: 404 });
  try {
    await assert.rejects(
      compileLosatWasm('https://example.test/missing.wasm', 'threaded LOSAT wasm'),
      /Missing packaged asset: threaded LOSAT wasm not found at https:\/\/example\.test\/missing\.wasm/,
    );
  } finally {
    globalThis.fetch = original;
  }
});

test('every LOSAT wasm loader fetches through compileLosatWasm', () => {
  for (const path of LOADERS) {
    const source = readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');
    assert.doesNotMatch(source, /no-store/, path);
    assert.match(source, /compileLosatWasm\(/, path);
    assert.doesNotMatch(source, /fetch\(\s*(?:wasmUrl|resolvedWasmPath)/, path);
  }
});
