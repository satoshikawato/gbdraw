// OV-2750: the History byte budget counts only the bytes History itself keeps
// alive for the files its steps name (`residentBytesOf`,
// services/file-content-cache.js). A File's bytes stay on disk (the chosen
// file, or its copy, PD-OI-121), and a gzip File's expanded bytes belong to the
// content cache, which lets go of them on its own; dropping Undo steps would
// free neither. So a bound input larger than the budget no longer discards the
// older Undo steps; a Session resource view, whose bytes a step keeps, counts.
import assert from 'node:assert/strict';
import test from 'node:test';
import { webcrypto } from 'node:crypto';
import { gzipSync } from 'node:zlib';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const { createHistoryRouter } = await import('../../gbdraw/web/js/services/history.js');
const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { readFileBytes, releaseFileContentCache, residentBytesOf } = await import(
  '../../gbdraw/web/js/services/file-content-cache.js'
);

const MIB = 1024 * 1024;
/** A File the user chose, of `size` bytes on disk; the test never reads it. */
const diskFile = (name, size) => new (class extends File {
  get size() { return size; }
})([], name);

test('History keeps no bytes alive for a File, gzip or not; a Session view counts its size', async () => {
  assert.equal(residentBytesOf(diskFile('depth.tsv', 300 * 1000 * 1000)), 0);
  assert.equal(residentBytesOf({ name: 'view.gb', size: 1234 }), 1234, 'a Session resource view counts its size');
  const plain = new TextEncoder().encode('chr\t1\t7\n'.repeat(1000));
  const gzip = new File([gzipSync(plain)], 'depth.tsv.gz');
  await readFileBytes(gzip);
  assert.equal(residentBytesOf(gzip), 0, 'its expanded bytes, which the content cache holds while no copy does, are not History\'s');
  releaseFileContentCache(gzip);
});

test('a bound File larger than the budget leaves every undo step (OV-2750)', async () => {
  const fileStore = createInputFileStore({ residentBytesOf });
  const notices = [];
  let value = '';
  const input = diskFile('depth_U.tsv.gz', 376 * MIB);
  const buildState = () => ({ value, file: fileStore.registerFile(input) });
  const applyState = (snapshot) => { value = snapshot.value; };
  const router = createHistoryRouter({
    drawingIds: ['circular'],
    activeDrawingId: () => 'circular',
    maxBytes: 200 * MIB,
    onHistoryLimit: (message) => { notices.push(message); },
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    estimateFileBytes: fileStore.estimateBytes,
    retainOnlyFiles: fileStore.retainOnly,
    managerOptions: () => ({ buildIntent: buildState, applyIntent: applyState, buildCheckpoint: buildState, applyCheckpoint: applyState })
  });
  await router.captureBaseline();
  for (const next of ['a', 'b', 'c']) await router.runUndoable(`Edit ${next}`, () => { value = next; });
  assert.deepEqual(notices, [], 'no step is discarded');
  await router.undo();
  await router.undo();
  assert.equal(value, 'a');
  assert.equal(router.canUndo(), true, 'Undo reaches the first step');
});

test('a gzip File that expands past the budget before its copy is written leaves every undo step (OV-2750)', async () => {
  const plain = new TextEncoder().encode('chr\t1\t7\n'.repeat(40000));
  const gzip = new File([gzipSync(plain)], 'depth.tsv.gz');
  await readFileBytes(gzip);
  const fileStore = createInputFileStore({ residentBytesOf });
  const notices = [];
  let value = '';
  const buildState = () => ({ value, file: fileStore.registerFile(gzip) });
  const applyState = (snapshot) => { value = snapshot.value; };
  const router = createHistoryRouter({
    drawingIds: ['circular'],
    activeDrawingId: () => 'circular',
    maxBytes: Math.floor(plain.byteLength / 2),
    onHistoryLimit: (message) => { notices.push(message); },
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    estimateFileBytes: fileStore.estimateBytes,
    retainOnlyFiles: fileStore.retainOnly,
    managerOptions: () => ({ buildIntent: buildState, applyIntent: applyState, buildCheckpoint: buildState, applyCheckpoint: applyState })
  });
  try {
    await router.captureBaseline();
    for (const next of ['a', 'b', 'c']) await router.runUndoable(`Edit ${next}`, () => { value = next; });
    assert.deepEqual(notices, [], 'no step is discarded');
    await router.undo();
    await router.undo();
    assert.equal(value, 'a');
  } finally {
    releaseFileContentCache(gzip);
  }
});

test('the project input store counts resident bytes', async () => {
  globalThis.window ??= { Vue: { ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {} } };
  const { state } = await import('../../gbdraw/web/js/state.js');
  const ref = state.inputFiles.registerFile(diskFile('big.gb', 300 * 1000 * 1000));
  assert.equal(state.inputFiles.estimateBytes([ref.fileId]), 0);
  state.inputFiles.retainOnly(new Set());
});
