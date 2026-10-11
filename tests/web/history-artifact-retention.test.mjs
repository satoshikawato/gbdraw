import assert from 'node:assert/strict';
import test from 'node:test';

import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';

// OV-112: adjacent Generate entries name the same artifact (entry k's `after`
// is entry k+1's `before`). History counts each artifact's retained bytes
// once while any live entry names it, so the byte limit keeps the undo depth
// that the unique artifacts fit.
const ARTIFACT_BYTES = 1000;
const LIMIT_MESSAGE = 'Older undo history was discarded to stay within the history limit.';

const artifactFor = (id, { fingerprinted = true, compactSignature = 's' } = {}) => ({
  id,
  identity: {
    fingerprint: fingerprinted ? id.repeat(64).slice(0, 64) : '',
    compactSignature
  },
  retainedBytes: ARTIFACT_BYTES,
  fileIds: []
});

// Each capture is a fresh handle, as `captureGeneratedArtifactHandle` in
// `services/history-snapshot.js` returns one.
const createGenerateHistory = ({ maxBytes = Number.MAX_SAFE_INTEGER, fingerprinted = true } = {}) => {
  let artifact = artifactFor('loaded', { fingerprinted });
  const history = createHistoryManager({
    buildIntent: async () => ({}),
    applyIntent: async () => {},
    buildCheckpoint: () => {
      throw new Error('Generate replacement must not build a checkpoint.');
    },
    applyCheckpoint: async () => {},
    captureGeneratedArtifactHandle: () => ({ ...artifact }),
    restoreGeneratedArtifactHandle: async (handle) => {
      artifact = handle;
    },
    maxBytes
  });
  const generate = (id, options = { fingerprinted }) => history.runUndoableArtifactReplacement(
    `Generate ${id}`,
    async () => {
      artifact = artifactFor(id, options);
      return { status: 'ok' };
    },
    { shouldCommit: (result) => result.status === 'ok' }
  );
  return {
    history,
    generate,
    currentId: () => artifact.id,
    editArtifact: (compactSignature) => {
      artifact = { ...artifact, identity: { ...artifact.identity, compactSignature } };
    },
    retainedBytes: () => history.getDiagnostics().retainedEntryBytes
  };
};

const GENERATED_IDS = ['a', 'b', 'c', 'd', 'e', 'f'];

const undoAll = async (history) => {
  let steps = 0;
  while (history.canUndo()) {
    assert.equal(await history.undo(), true);
    steps += 1;
  }
  return steps;
};

test('six Generates keep every entry when the unique artifacts fit the limit', async () => {
  const { history, generate, currentId, retainedBytes } = createGenerateHistory({
    maxBytes: 7500
  });
  await history.initializeIntentBaseline('Loaded session');
  for (const id of GENERATED_IDS) await generate(id);

  assert.equal(history.getUndoCount(), 6);
  assert.equal(retainedBytes(), 7 * ARTIFACT_BYTES, 'loaded + six generated artifacts');
  assert.equal(history.historyLimitMessage.value, '');
  assert.equal(await undoAll(history), 6);
  assert.equal(currentId(), 'loaded');
  assert.equal(retainedBytes(), 7 * ARTIFACT_BYTES, 'Undo moves entries to Redo; they stay retained');
});

test('evicting the oldest entry keeps the artifact that the next entry still names', async () => {
  const { history, generate, currentId, retainedBytes } = createGenerateHistory({
    maxBytes: 5500
  });
  await history.initializeIntentBaseline('Loaded session');
  for (const id of GENERATED_IDS.slice(0, 4)) await generate(id);
  assert.equal(history.getUndoCount(), 4);
  assert.equal(retainedBytes(), 5 * ARTIFACT_BYTES);
  assert.equal(history.historyLimitMessage.value, '');

  // The fifth entry exceeds the limit: evicting "Generate a" releases only the
  // loaded artifact, because "Generate b" still names artifact a as `before`.
  await generate('e');
  assert.equal(history.getUndoCount(), 4);
  assert.equal(retainedBytes(), 5 * ARTIFACT_BYTES);
  assert.equal(history.historyLimitMessage.value, LIMIT_MESSAGE);

  await generate('f');
  assert.equal(history.getUndoCount(), 4);
  assert.equal(retainedBytes(), 5 * ARTIFACT_BYTES);
  assert.equal(history.undoLabel(), 'Generate f');
  assert.equal(await undoAll(history), 4);
  assert.equal(currentId(), 'b', 'the oldest kept entry restores the artifact it names');
});

test('a Generate after Undo releases the truncated Redo artifacts once', async () => {
  const { history, generate, currentId, retainedBytes } = createGenerateHistory();
  await history.initializeIntentBaseline('Loaded session');
  for (const id of GENERATED_IDS) await generate(id);
  for (let step = 0; step < 3; step += 1) await history.undo();
  assert.equal(currentId(), 'c');
  assert.equal(history.getRedoCount(), 3);
  assert.equal(retainedBytes(), 7 * ARTIFACT_BYTES);

  // Artifact c stays named by "Generate c" and the new entry; d, e, and f
  // were named only by the truncated Redo entries.
  await generate('x');
  assert.equal(history.getRedoCount(), 0);
  assert.equal(history.getUndoCount(), 4);
  assert.equal(retainedBytes(), 5 * ARTIFACT_BYTES, 'loaded, a, b, c, and x');
  assert.equal(await undoAll(history), 4);
  assert.equal(currentId(), 'loaded');

  await history.initializeIntentBaseline('Session replaced');
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(retainedBytes(), 0, 'clearing both stacks releases every artifact');
});

test('artifacts without a shared identity still count per handle', async () => {
  // Control: a handle without a fingerprint names only itself, so the
  // accounting and eviction match the per-entry sum.
  const unfingerprinted = createGenerateHistory({ maxBytes: 7500, fingerprinted: false });
  await unfingerprinted.history.initializeIntentBaseline('Loaded session');
  for (const id of GENERATED_IDS) await unfingerprinted.generate(id);
  assert.equal(unfingerprinted.history.getUndoCount(), 3);
  assert.equal(unfingerprinted.retainedBytes(), 6 * ARTIFACT_BYTES);
  assert.equal(unfingerprinted.history.historyLimitMessage.value, LIMIT_MESSAGE);
  assert.equal(await undoAll(unfingerprinted.history), 3);
  assert.equal(unfingerprinted.currentId(), 'c');

  // An edit between Generates changes the compact signature, so the next
  // `before` is another artifact than the previous `after`.
  const edited = createGenerateHistory();
  await edited.history.initializeIntentBaseline('Loaded session');
  await edited.generate('a');
  edited.editArtifact('edited');
  await edited.generate('b');
  assert.equal(edited.history.getUndoCount(), 2);
  assert.equal(edited.retainedBytes(), 4 * ARTIFACT_BYTES, 'loaded, a, a edited, and b');
});
