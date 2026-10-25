// The checks a Pyodide snapshot passes before the diagram Worker restores it
// (gbdraw/web/js/workers/pyodide-snapshot.js): the validity key covers every
// wheel's bytes, and a snapshot with a table slot or a JS reference it cannot
// restore is refused. A refused snapshot makes the Worker boot cold.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

const {
  PYTHON_HELPER_FILE, fetchPyodideSnapshot, installRuntimePackages, pyodideSnapshotKey, runtimeWheelDownloads
} = await import('../../gbdraw/web/js/workers/pyodide-snapshot.js');

const bytes = (text) => new TextEncoder().encode(text);
const FOLDER = 'http://test.invalid/snapshot/';
const GBDRAW_WHEEL = 'http://test.invalid/gbdraw-0.0.0-py3-none-any.whl';
const LOCAL_WHEELS = ['http://test.invalid/a-1.0-py3-none-any.whl', 'http://test.invalid/b-1.0-py3-none-any.whl'];
const helperFile = readFileSync(PYTHON_HELPER_FILE);
const wheels = { [LOCAL_WHEELS[0]]: bytes('wheel a'), [LOCAL_WHEELS[1]]: bytes('wheel b'), [GBDRAW_WHEEL]: bytes('gbdraw') };

const stateFor = async (changes = {}) => ({
  key: await pyodideSnapshotKey({
    localWheels: LOCAL_WHEELS.map((url) => [url, wheels[url]]), gbdrawWheel: wheels[GBDRAW_WHEEL], helperFile
  }),
  libs: [], extraSlots: [], unknownSlots: [], got: {}, tableLength: 0, firstLibTable: 0,
  loadedPackages: {}, unknownHiwire: [],
  ...changes
});

// Serves `files` (URL to bytes) while `run` runs, records each request, and
// answers 404 for any other URL.
const withServedFiles = async (files, run) => {
  const original = globalThis.fetch;
  const requests = [];
  globalThis.fetch = async (url, options = {}) => {
    requests.push([options.method || 'GET', String(url), options.cache || 'default']);
    return String(url) in files ? new Response(files[String(url)]) : new Response(null, { status: 404 });
  };
  try {
    return { value: await run(), requests };
  } finally {
    globalThis.fetch = original;
  }
};
const URLS = { localWheelUrls: LOCAL_WHEELS, gbdrawWheelUrl: GBDRAW_WHEEL };

// Serves the snapshot files and the wheels; `served` replaces some of them.
const fetchSnapshot = async (state, served = {}) => {
  const files = {
    [`${FOLDER}state.json`]: bytes(JSON.stringify(state)),
    [`${FOLDER}snapshot.bin`]: new Uint8Array(16),
    [`${FOLDER}fs.bin`]: new Uint8Array(4),
    [PYTHON_HELPER_FILE.href]: helperFile,
    ...wheels,
    ...served
  };
  return (await withServedFiles(files, () => fetchPyodideSnapshot({ folderUrl: FOLDER, wheels: runtimeWheelDownloads(URLS) }))).value;
};

const refusedWith = (pattern) => (error) => {
  assert.equal(error.code, 'WORKER_INIT');
  assert.equal(error.context.reason, 'pyodide-snapshot');
  assert.match(error.context.detail, pattern);
  return true;
};

test('a snapshot of the served wheels and helper source is accepted', async () => {
  const snapshot = await fetchSnapshot(await stateFor());
  assert.equal(snapshot.memory.length, 16);
});

test('a dependency wheel with other bytes under the same name refuses the snapshot', async () => {
  await assert.rejects(
    fetchSnapshot(await stateFor(), { [LOCAL_WHEELS[1]]: bytes('wheel b, rebuilt') }),
    refusedWith(/another recipe, wheel or helper source/)
  );
});

test('a snapshot with a table slot or a JS reference it cannot restore is refused', async () => {
  await assert.rejects(fetchSnapshot(await stateFor({ unknownSlots: [7] })), refusedWith(/cannot restore/));
  await assert.rejects(fetchSnapshot(await stateFor({ unknownHiwire: ['object [object Foo]'] })), refusedWith(/cannot restore/));
});

// The cold boot downloads each wheel once (a GET, no HEAD check) and
// installs them all in one micropip call from the Emscripten file system,
// after loading micropip and tzdata from the Pyodide distribution.
test('the recipe installs the downloaded wheels in one micropip call', async () => {
  const calls = [];
  const files = new Map();
  const pyodide = {
    loadPackage: async (names) => { calls.push(['loadPackage', names]); },
    pyimport: () => ({
      install: async (requirements) => {
        calls.push(['install', requirements, requirements.map((path) => new TextDecoder().decode(files.get(path.replace(/^emfs:/, ''))))]);
      }
    }),
    FS: {
      mkdirTree: () => {}, writeFile: (path, data) => files.set(path, data), unlink: (path) => files.delete(path), rmdir: () => {}
    }
  };
  const { requests } = await withServedFiles(wheels, () => installRuntimePackages(pyodide, {
    wheels: runtimeWheelDownloads(URLS), beforeInstall: async () => { calls.push(['beforeInstall']); }
  }));
  assert.deepEqual(requests, [...LOCAL_WHEELS, GBDRAW_WHEEL].map((url) => ['GET', url, 'default']));
  const folder = '/tmp/gbdraw-runtime-wheels';
  assert.deepEqual(calls, [
    ['loadPackage', ['micropip', 'tzdata']],
    ['beforeInstall'],
    ['install', [`emfs:${folder}/a-1.0-py3-none-any.whl`, `emfs:${folder}/b-1.0-py3-none-any.whl`,
      `emfs:${folder}/gbdraw-0.0.0-py3-none-any.whl`], ['wheel a', 'wheel b', 'gbdraw']]
  ]);
  assert.equal(files.size, 0);
});

test('a wheel that is not served fails the download with its name', async () => {
  const { value: downloads } = await withServedFiles({ [LOCAL_WHEELS[0]]: wheels[LOCAL_WHEELS[0]] }, () => runtimeWheelDownloads(URLS));
  await assert.rejects(downloads[1].bytes, { message: `Missing packaged asset: Pyodide dependency wheel #2 (404) at ${LOCAL_WHEELS[1]}` });
  await assert.rejects(downloads[2].bytes, { message: `Missing packaged asset: gbdraw browser wheel (404) at ${GBDRAW_WHEEL}` });
});
