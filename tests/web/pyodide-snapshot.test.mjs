// The checks a Pyodide snapshot passes before the diagram Worker restores it
// (gbdraw/web/js/workers/pyodide-snapshot.js): the validity key covers every
// wheel's bytes, and a snapshot with a table slot or a JS reference it cannot
// restore is refused. A refused snapshot makes the Worker boot cold.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

const { PYTHON_HELPER_FILE, fetchPyodideSnapshot, pyodideSnapshotKey } = await import(
  '../../gbdraw/web/js/workers/pyodide-snapshot.js'
);

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
  const original = globalThis.fetch;
  globalThis.fetch = async (url) => (String(url) in files
    ? new Response(files[String(url)])
    : new Response(null, { status: 404 }));
  try {
    return await fetchPyodideSnapshot({ folderUrl: FOLDER, gbdrawWheelUrl: GBDRAW_WHEEL, localWheelUrls: LOCAL_WHEELS });
  } finally {
    globalThis.fetch = original;
  }
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
