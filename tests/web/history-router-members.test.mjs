// The History router owns which drawings have a History (PD-OI-090, R16-3):
// `addMember` gives a new drawing an empty History, `removeMember` hands a
// drawing's History to the project step that removes it, which holds it in
// the shared byte budget until that step is evicted, and `replaceMembers`
// (a Load) clears the project stream and every History, then installs the
// new list.
import assert from 'node:assert/strict';
import test from 'node:test';

import { createHistoryRouter } from '../../gbdraw/web/js/services/history.js';
import { createInputFileStore } from '../../gbdraw/web/js/services/history-files.js';

const ref = (value) => ({ value });
const makeFile = (name, size = 10) => ({ name, size, type: 'text/plain', lastModified: 1 });

const makeRouter = ({ drawingIds, shownId, maxBytes = Number.MAX_SAFE_INTEGER }) => {
  const fileStore = createInputFileStore();
  const shown = ref(shownId);
  const values = {};
  const files = {};
  const router = createHistoryRouter({
    drawingIds,
    activeDrawingId: () => shown.value,
    maxBytes,
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    estimateFileBytes: fileStore.estimateBytes,
    retainOnlyFiles: fileStore.retainOnly,
    managerOptions: (drawingId) => {
      values[drawingId] ??= '';
      files[drawingId] ??= null;
      const buildState = () => ({
        value: values[drawingId],
        file: files[drawingId] ? fileStore.registerFile(files[drawingId]) : null
      });
      const applyState = (snapshot) => {
        values[drawingId] = snapshot.value;
        files[drawingId] = snapshot.file ? fileStore.fileOf(snapshot.file) : null;
      };
      return { buildIntent: buildState, applyIntent: applyState, buildCheckpoint: buildState, applyCheckpoint: applyState };
    }
  });
  return { router, shown, values, files, fileStore };
};

/** @param {string} reason */
const refusal = (reason) => (error) => (
  error.code === 'INPUT_INVALID' && error.stage === 'history-router' && error.context.reason === reason);

test('a drawing added to the router gets an empty History that records and undoes its steps', async () => {
  const { router, shown, values } = makeRouter({ drawingIds: ['circular', 'linear'], shownId: 'linear' });
  await router.captureBaseline();
  await router.runUndoable('Linear edit', () => { values.linear = 'a'; });

  router.addMember('linear-2');
  shown.value = 'linear-2';
  assert.equal(router.canUndo(), false, 'a new drawing starts with an empty History');
  await router.runUndoable('Linear 2 edit', () => { values['linear-2'] = 'b'; });
  assert.equal(router.undoLabel(), 'Linear 2 edit');
  await router.undo();
  assert.equal(values['linear-2'], '');
  assert.equal(values.linear, 'a', 'the Undo in the new drawing leaves the other drawing');
  assert.equal(router.canRedo(), true);

  assert.throws(() => router.addMember('linear'), refusal('DUPLICATE_DRAWING'));
  assert.throws(() => router.managerOf('linear-3'), refusal('UNKNOWN_DRAWING'));
});

test('a removed History stays in the budget, held by the project step that removed it, until that step is evicted', async () => {
  const probe = makeRouter({ drawingIds: ['circular'], shownId: 'circular' });
  await probe.router.captureBaseline();
  await probe.router.runUndoable('Probe', () => { probe.values.circular = 'x'.repeat(400); });
  const stepBytes = probe.router.getDiagnostics().retainedEntryBytes;
  assert.ok(stepBytes > 0);

  const { router, shown, values, files, fileStore } = makeRouter({
    drawingIds: ['circular', 'linear', 'linear-2'], shownId: 'linear-2', maxBytes: Math.floor(stepBytes * 1.5)
  });
  files['linear-2'] = makeFile('a.gb');
  await router.captureBaseline();
  const fileAId = router.managerOf('linear-2').getCurrentIntent().file.fileId;
  await router.runUndoable('Linear 2 edit', () => { values['linear-2'] = 'x'.repeat(400); files['linear-2'] = null; });
  shown.value = 'circular';

  assert.throws(() => router.removeMember('linear-2'), refusal('NO_PROJECT_STEP'), 'a removal is a project step');
  await router.runProjectStep('Delete Linear 2', () => { router.removeMember('linear-2'); });
  assert.throws(() => router.managerOf('linear-2'), refusal('UNKNOWN_DRAWING'));
  await assert.rejects(
    router.runProjectStep('Delete again', () => { router.removeMember('linear-2'); }),
    refusal('UNKNOWN_DRAWING')
  );
  assert.equal(router.undoLabel(), 'Delete Linear 2');
  assert.equal(fileStore.has(fileAId), true, 'the held History still names file A');

  // The held History's bytes count: one more step is over the budget, and the
  // project step that holds it goes before the shown drawing's step.
  await router.runUndoable('Circular edit', () => { values.circular = 'y'.repeat(400); });
  assert.equal(router.getUndoCount(), 1, 'the holding project step was evicted');
  assert.equal(router.undoLabel(), 'Circular edit');
  assert.equal(fileStore.has(fileAId), false, 'evicting the holder releases the held History and its files');
});

test('a project step that is not recorded gives a removed History back', async () => {
  const { router, shown, values } = makeRouter({ drawingIds: ['circular', 'linear'], shownId: 'linear' });
  await router.captureBaseline();
  await router.runUndoable('Linear edit', () => { values.linear = 'a'; });
  shown.value = 'circular';
  await assert.rejects(router.runProjectStep('Delete Linear', () => {
    router.removeMember('linear');
    throw new Error('the deletion failed');
  }), /the deletion failed/);
  assert.equal(router.managerOf('linear').canUndo(), true, 'the History is a member again with its steps');
  assert.equal(router.getUndoCount(), 0);
});

test('replacing the members clears the project stream and every History', async () => {
  const { router, shown, values } = makeRouter({ drawingIds: ['circular', 'linear'], shownId: 'circular' });
  await router.captureBaseline();
  await router.runProjectStep('Reset all drawings', () => { values.circular = 'r'; values.linear = 'r'; });
  await router.runUndoable('Circular edit', () => { values.circular = 'c'; });
  assert.equal(router.getUndoCount(), 2);

  assert.throws(() => router.replaceMembers(['circular', 'circular']), refusal('DUPLICATE_DRAWING'));
  assert.equal(router.getUndoCount(), 2, 'a refused list changes nothing');

  router.replaceMembers(['circular', 'linear-2']);
  assert.equal(router.getUndoCount(), 0, 'the shown drawing and the project stream hold no step');
  assert.equal(router.getRedoCount(), 0);
  assert.equal(router.canUndo(), false);
  assert.throws(() => router.managerOf('linear'), refusal('UNKNOWN_DRAWING'));
  shown.value = 'linear-2';
  assert.equal(router.canUndo(), false);
  await router.runUndoable('Linear 2 edit', () => { values['linear-2'] = 'l'; });
  await router.undo();
  assert.equal(values['linear-2'], '');
});

// OV-1769: a Load replaces the members before its History baseline; when the
// Load fails or is canceled after that, its rollback puts the members it found
// back, each History with its steps, and the project stream with its steps.
test('restoring captured members gives back every History and the project stream a replacement cleared', async () => {
  const { router, shown, values } = makeRouter({ drawingIds: ['circular', 'linear'], shownId: 'circular' });
  await router.captureBaseline();
  await router.runProjectStep('Reset all drawings', () => { values.circular = 'r'; values.linear = 'r'; });
  await router.runUndoable('Circular edit', () => { values.circular = 'c'; });
  await router.undo();
  assert.deepEqual([router.getUndoCount(), router.getRedoCount()], [1, 1]);

  const members = router.captureMembers();
  router.replaceMembers(['circular-2', 'linear-2']);
  shown.value = 'circular-2';
  assert.deepEqual([router.getUndoCount(), router.getRedoCount()], [0, 0]);

  shown.value = 'circular';
  router.restoreMembers(members);
  assert.throws(() => router.managerOf('circular-2'), refusal('UNKNOWN_DRAWING'));
  assert.deepEqual([router.getUndoCount(), router.getRedoCount()], [1, 1], 'the shown History and the project stream are back');
  assert.equal(router.redoLabel(), 'Circular edit');
  await router.redo();
  assert.equal(values.circular, 'c');
  await router.undo();
  await router.undo();
  assert.deepEqual([values.circular, values.linear], ['', ''], 'the project step undoes in both drawings');
});
