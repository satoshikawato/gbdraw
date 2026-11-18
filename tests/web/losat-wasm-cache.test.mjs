// LOSAT wasm loads revalidate instead of downloading again: every loader
// fetches through `compileLosatWasm` (services/losat-runtime.js) with
// `cache: 'no-cache'`, so an unchanged file costs a 304 and a replaced file
// at the same URL is never served stale. A page runs only the LOSAT build it
// records (OV-2440): it compiles a file only when its sha256 is the recorded
// one, and each LOSAT worker refuses a message from code of another build.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

import { fetchShippedLosatWasm, shippedLosatWasm } from './helpers/shipped-losat-wasm.mjs';

const web = new URL('../../gbdraw/web/', import.meta.url);
const { compileLosatWasm, losatBuildChangedFailure } = await import(new URL('js/services/losat-runtime.js', web).href);
const { LOSAT_WASM_BUILD } = await import(new URL('js/services/losat-build.generated.js', web).href);
const { normalizeUserFacingError } = await import(new URL('js/utils/error-normalization.js', web).href);
const wasiShimUrl = new URL('vendor/browser_wasi_shim/dist/index.js', web).href;

const LOADERS = [
  'gbdraw/web/js/services/losat.js',
  'gbdraw/web/js/workers/losat-worker.js',
  'gbdraw/web/js/workers/losat-threaded-worker.js',
  'gbdraw/web/js/workers/losat-wasi-thread-worker.js',
];
// A valid module whose sha256 the build does not record: a file of another build.
const UNRECORDED = new Uint8Array([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0]);
const OTHER_BUILD = '0.1.0+000000000000';

/** @param {unknown} error */
const isBuildChanged = (error) => error instanceof Error && error.code === 'LOSAT_BUILD_CHANGED' && error.stage === 'losat';

const stubFetch = (t, respond) => {
  const original = globalThis.fetch;
  globalThis.fetch = respond;
  t.after(() => { globalThis.fetch = original; });
};

test('compileLosatWasm revalidates the wasm and compiles it', async (t) => {
  const calls = [];
  stubFetch(t, async (url, init) => {
    calls.push([url, init]);
    return fetchShippedLosatWasm(url);
  });
  const module = await compileLosatWasm('https://example.test/wasm/losat/losat.wasm', 'losat.wasm');
  assert.ok(module instanceof WebAssembly.Module);
  assert.deepEqual(calls, [['https://example.test/wasm/losat/losat.wasm', { cache: 'no-cache' }]]);
});

test('compileLosatWasm compiles the bytes as served: LOSAT guards the threaded file when it builds it', async (t) => {
  const served = shippedLosatWasm('losat-threaded.wasm');
  const originalCompile = WebAssembly.compile;
  const compiled = [];
  stubFetch(t, fetchShippedLosatWasm);
  WebAssembly.compile = (bytes) => {
    compiled.push(new Uint8Array(bytes));
    return originalCompile(bytes);
  };
  t.after(() => { WebAssembly.compile = originalCompile; });
  await compileLosatWasm('https://example.test/wasm/losat/losat-threaded.wasm', 'losat-threaded.wasm');
  assert.equal(compiled.length, 1);
  assert.equal(Buffer.compare(compiled[0], served), 0);
});

test('compileLosatWasm names a missing wasm', async (t) => {
  stubFetch(t, async () => new Response('', { status: 404 }));
  await assert.rejects(
    compileLosatWasm('https://example.test/missing.wasm', 'losat-threaded.wasm'),
    /Missing packaged asset: losat-threaded\.wasm not found at https:\/\/example\.test\/missing\.wasm/,
  );
});

// OV-2440: a page open across an update fetched the new build's file and ran it under its own record of the build.
test('compileLosatWasm refuses bytes whose sha256 the build does not record, before compiling them', async (t) => {
  const originalCompile = WebAssembly.compile;
  let compiled = 0;
  WebAssembly.compile = (bytes) => { compiled += 1; return originalCompile(bytes); };
  t.after(() => { WebAssembly.compile = originalCompile; });
  stubFetch(t, async () => new Response(UNRECORDED, { status: 200 }));
  for (const file of ['losat.wasm', 'losat-threaded.wasm']) {
    await assert.rejects(compileLosatWasm(`https://example.test/wasm/losat/${file}`, file), isBuildChanged);
  }
  // The other shipped file is not this one either.
  stubFetch(t, async () => fetchShippedLosatWasm('losat.wasm'));
  await assert.rejects(
    compileLosatWasm('https://example.test/wasm/losat/losat-threaded.wasm', 'losat-threaded.wasm'), isBuildChanged
  );
  assert.equal(compiled, 0);
});

test('a changed LOSAT build tells the user to reload the page', () => {
  const model = normalizeUserFacingError(losatBuildChangedFailure('detail for the console'), { operation: 'generate' });
  assert.equal(model.code, 'LOSAT_BUILD_CHANGED');
  assert.equal(model.summary, 'gbdraw was updated while this page was open, so the comparison search did not run. Save a Session, reload the page, then run the search again.');
  assert.deepEqual(model.actions, ['save-session', 'reload']);
  assert.doesNotMatch(JSON.stringify(model), /detail for the console/);
});

/**
 * A LOSAT worker module loaded fresh, with what it posts.
 * @param {import('node:test').TestContext} t
 * @param {string} path
 */
const loadWorker = async (t, path) => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'self');
  t.after(() => { if (previous) Object.defineProperty(globalThis, 'self', previous); else delete globalThis.self; });
  const posted = [];
  globalThis.self = { postMessage: (message) => posted.push(message), close() {} };
  await import(`${new URL(path, web).href}?instance=${Date.now()}-${Math.random()}`);
  return { onmessage: globalThis.self.onmessage, posted };
};

// OV-2440: workers are made per run, so a page open across an update runs the new worker code with its old module.
test('the serial pool worker refuses a page of another LOSAT build', async (t) => {
  const { onmessage, posted } = await loadWorker(t, 'js/workers/losat-worker.js');
  const wasmModule = await WebAssembly.compile(shippedLosatWasm('losat.wasm'));
  await onmessage({ data: { type: 'init', id: 'old', losatBuild: OTHER_BUILD, wasmModule, wasiShimUrl } });
  await onmessage({ data: { type: 'init', id: 'none', wasmModule, wasiShimUrl } });
  await onmessage({ data: { type: 'init', id: 'same', losatBuild: LOSAT_WASM_BUILD.version, wasmModule, wasiShimUrl } });
  await onmessage({ data: { type: 'run', id: 'run', losatBuild: OTHER_BUILD, job: {} } });
  assert.deepEqual(posted.map(({ id, ok, buildChanged }) => ({ id, ok, buildChanged })), [
    { id: 'old', ok: false, buildChanged: true },
    { id: 'none', ok: false, buildChanged: true },
    { id: 'same', ok: true, buildChanged: undefined },
    { id: 'run', ok: false, buildChanged: true }
  ]);
});

test('the threaded job worker refuses a page of another LOSAT build without running the job', async (t) => {
  const { onmessage, posted } = await loadWorker(t, 'js/workers/losat-threaded-worker.js');
  const previousWorker = Object.getOwnPropertyDescriptor(globalThis, 'Worker');
  t.after(() => { if (previousWorker) Object.defineProperty(globalThis, 'Worker', previousWorker); else delete globalThis.Worker; });
  let made = 0;
  globalThis.Worker = class { constructor() { made += 1; throw new Error('not reached'); } };
  await onmessage({
    data: {
      type: 'run', id: 'job', losatBuild: OTHER_BUILD, module: await WebAssembly.compile(shippedLosatWasm('losat-threaded.wasm')),
      wasiShimUrl, threadsPerJob: 2, faultChannel: '', job: { program: 'blastn', extraArgs: [] }, queryFasta: '>q\nACGT\n', subjectFasta: '>s\nACGT\n'
    }
  });
  assert.deepEqual(posted.map(({ id, ok, buildChanged, hostFailure }) => ({ id, ok, buildChanged, hostFailure })),
    [{ id: 'job', ok: false, buildChanged: true, hostFailure: false }]);
  assert.equal(made, 0);
});

test('a WASI thread worker refuses a job worker of another LOSAT build', async (t) => {
  const { onmessage, posted } = await loadWorker(t, 'js/workers/losat-wasi-thread-worker.js');
  const readyControl = new SharedArrayBuffer(4);
  await onmessage({ data: { type: 'prepare', id: 'prepare', losatBuild: OTHER_BUILD, readyControl, wasiShimUrl } });
  const control = new SharedArrayBuffer(4);
  await onmessage({ data: { type: 'run', losatBuild: OTHER_BUILD, tid: 1, control, wasiShimUrl } });
  assert.deepEqual(posted.map(({ id, type, ok, buildChanged }) => ({ id, type, ok, buildChanged })), [
    { id: 'prepare', type: 'prepare', ok: false, buildChanged: true },
    { id: undefined, type: 'error', ok: undefined, buildChanged: true }
  ]);
  assert.equal(Atomics.load(new Int32Array(readyControl), 0), -1);
  assert.equal(Atomics.load(new Int32Array(control), 0), -1);
});

test('every LOSAT wasm loader fetches through compileLosatWasm and names the file it loads', () => {
  for (const path of LOADERS) {
    const source = readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');
    assert.doesNotMatch(source, /no-store/, path);
    assert.match(source, /compileLosatWasm\(/, path);
    assert.doesNotMatch(source, /fetch\(\s*(?:wasmUrl|resolvedWasmPath)/, path);
    assert.doesNotMatch(source, /compileLosatWasm\([^)]*'(?!losat(?:-threaded)?\.wasm')[^']*'\)/, path);
  }
});
