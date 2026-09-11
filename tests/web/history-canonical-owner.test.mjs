import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

globalThis.window = { Vue: {
  ref: value => ({ value }), reactive: value => value,
  computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
const {
  adoptCanonicalRenderArtifacts, canonicalRenderArtifactOwner,
  getCommittedCanonicalRenderRequest, buildConfigData, applyConfigData
} = await import('../../gbdraw/web/js/services/config.js');
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const { createHistoryFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');

const empty = canonicalRenderArtifactOwner.capture();
const canonical = JSON.parse(await readFile(
  'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', 'utf8'
));
adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
const a = canonicalRenderArtifactOwner.capture();
assert(Object.isFrozen(a));
assert.strictEqual(a.committedCanonicalSession.renderRequest, getCommittedCanonicalRenderRequest());
state.form.labels_mode = 'out';
state.results.value = [{ name: 'a.svg', content: '<svg id="a"/>' }];
const snapshots = createHistorySnapshotService({
  state, fileStore: createHistoryFileStore(), buildConfigData, applyConfigData
});
snapshots.setGeneratedArtifactRuntimeOwner(canonicalRenderArtifactOwner);
const history = createHistoryManager({
  buildIntent: snapshots.buildHistoryIntent, applyIntent: snapshots.applyHistoryIntent,
  buildCheckpoint: () => assert.fail('Generate must not clone a full checkpoint'),
  applyCheckpoint: () => assert.fail('Generate must restore its owner handle'),
  captureGeneratedArtifactHandle: snapshots.captureGeneratedArtifactHandle,
  restoreGeneratedArtifactHandle: snapshots.restoreGeneratedArtifactHandle,
  compareGeneratedArtifactHandles: snapshots.compareGeneratedArtifactHandles
});
await history.initializeIntentBaseline();
await history.runUndoable('Change labels', () => { state.form.labels_mode = 'none'; });
const next = structuredClone(canonical);
next.renderRequest.diagramOptions.configOverrides['labels.circular.scope'] = 'none';
await history.runUndoableArtifactReplacement('Generate diagram', () => {
  adoptCanonicalRenderArtifacts(next, { adoptOwnedRequest: true });
  state.results.value = [{ name: 'b.svg', content: '<svg id="b"/>' }];
});
const b = canonicalRenderArtifactOwner.capture();
assert.notStrictEqual(b.activeSessionResourceTable, a.activeSessionResourceTable);
await history.undo();
assert.strictEqual(getCommittedCanonicalRenderRequest(), a.committedCanonicalSession.renderRequest);
assert.strictEqual(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, a.activeSessionResourceTable);
assert.equal(state.results.value[0].name, 'a.svg');
assert.equal(state.form.labels_mode, 'none');
await history.redo();
assert.strictEqual(getCommittedCanonicalRenderRequest(), b.committedCanonicalSession.renderRequest);
assert.strictEqual(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, b.activeSessionResourceTable);
assert.equal(state.results.value[0].name, 'b.svg');
await history.undo();
await history.undo();
assert.equal(state.form.labels_mode, 'out');

await assert.rejects(history.runUndoableArtifactReplacement('Failed Generate', () => {
  adoptCanonicalRenderArtifacts(next, { adoptOwnedRequest: true });
  throw new Error('finalization failed');
}), /finalization failed/);
assert.strictEqual(getCommittedCanonicalRenderRequest(), a.committedCanonicalSession.renderRequest);
assert.strictEqual(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, a.activeSessionResourceTable);
canonicalRenderArtifactOwner.restore(empty);
assert.equal(getCommittedCanonicalRenderRequest(), null);
assert.equal(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, null);
state.results.value = [];
await history.initializeIntentBaseline();
await history.runUndoableArtifactReplacement('First Generate', () => {
  adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
  state.results.value = [{ name: 'first.svg', content: '<svg/>' }];
});
await history.undo();
assert.equal(getCommittedCanonicalRenderRequest(), null);
assert.equal(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, null);
assert.equal(state.results.value.length, 0);
await history.redo();
assert.deepEqual(getCommittedCanonicalRenderRequest(), canonical.renderRequest);
assert.equal(state.results.value[0].name, 'first.svg');
console.log('Generated History restores canonical request/resource owners without consuming draft History.');
