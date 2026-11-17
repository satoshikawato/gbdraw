// OV-3532 (OV-3504, OV-3505): a Linear drawing's Similarity alignment plan and
// record translations are its own (PD-OI-090: each drawing keeps its own
// record display and edits, History and Result; an edit, Generate, Undo or Load
// of one drawing never changes another's). They are artifact slot values: the
// shown drawing's are in `state`, another drawing's wait in its slot (R15), so
// a switch, Duplicate, New, Delete, Reset, Undo and Redo carry them per drawing.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} }
};
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');

/** @param {string} groupId */
const planOf = (groupId) => {
  const anchor = { recordKey: 'r1', biologicalFeatureId: `${groupId}-f1`, sourceFeatureIndex: 0, stableFeatureSvgId: null };
  return { schema: 2, groupId, reference: anchor, records: [{ recordKey: 'r1', status: 'reference', rationale: 'reference', anchor }] };
};
/** @param {number} x */
const translationsOf = (x) => [{ recordKey: 'r1', x, y: 0 }];
/** @param {string} name */
const resultOf = (name) => ({ name, content: `<svg id="${name}"/>` });
const shownAlignment = () => ({
  plan: state.similarityAlignmentPlan.value?.groupId ?? null,
  translations: state.linearRecordTranslations.value.map(({ x }) => x)
});
/** @param {ReturnType<typeof createAppSetup>} app @param {string} id */
const show = (app, id) => {
  if (state.activeDrawingId.value !== id) assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
};
/** The Align of a Similarity group: a Generate step that commits the plan and translations with its Result. */
const align = (history, groupId, x) => history.runUndoableArtifactReplacement('Align', () => {
  state.results.value = [resultOf(`${groupId}.svg`)];
  state.orthogroups.value = [{ id: groupId }];
  state.similarityAlignmentPlan.value = planOf(groupId);
  state.linearRecordTranslations.value = translationsOf(x);
  return true;
});
/** A manual record move: a step that writes the translations only. */
const drag = (history, x) => history.runUndoable('Move record', () => {
  state.linearRecordTranslations.value = translationsOf(x);
});

// OV-3532 (a): Align in a copy shows nothing on its source, and the copy's Undo
// changes only the copy.
test('an Align in the copy and its Undo leave the source drawing\'s alignment', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  await history.runUndoableArtifactReplacement('Generate diagram', () => {
    state.results.value = [resultOf('source.svg')];
    return true;
  });
  const { drawingId: copy } = app.drawingList.duplicateDrawing(source);
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'the copy has the source\'s (none)');

  await align(history, 'og-copy', 110);
  show(app, source);
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'the source shows no alignment');
  show(app, copy);
  assert.deepEqual(shownAlignment(), { plan: 'og-copy', translations: [110] });

  // The source aligns too; the copy's Undo takes back only the copy's Align.
  show(app, source);
  await align(history, 'og-source', 40);
  show(app, copy);
  assert.equal(await app.undoHistory(), true);
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'the copy\'s Undo restores the copy');
  show(app, source);
  assert.deepEqual(shownAlignment(), { plan: 'og-source', translations: [40] }, 'and leaves the source');
  show(app, copy);
  assert.equal(await app.redoHistory(), true);
  assert.deepEqual(shownAlignment(), { plan: 'og-copy', translations: [110] }, 'Redo gives the copy its Align back');
  show(app, source);
  assert.deepEqual(shownAlignment(), { plan: 'og-source', translations: [40] });
});

// OV-3505: Duplicate copies the record translations; a move in the copy stays in the copy.
test('Duplicate copies the plan and translations; a record move in the copy leaves the source', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  await align(history, 'og-1', 12);
  await drag(history, 15);
  const { drawingId: copy } = app.drawingList.duplicateDrawing(source);
  assert.deepEqual(shownAlignment(), { plan: 'og-1', translations: [15] }, 'Duplicate copies them');
  await drag(history, 99);
  show(app, source);
  assert.deepEqual(shownAlignment(), { plan: 'og-1', translations: [15] }, 'the source keeps its translations');
  show(app, copy);
  assert.deepEqual(shownAlignment(), { plan: 'og-1', translations: [99] });
  assert.equal(await app.undoHistory(), true);
  assert.deepEqual(shownAlignment(), { plan: 'og-1', translations: [15] }, 'Undo of the move is the copy\'s');
});

// OV-3532 (c), OV-3504: an edit that clears the plan clears only the drawing it edits.
test('a File edit in a new Linear drawing leaves another drawing\'s plan', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: aligned } = app.drawingList.newDrawing('linear');
  await align(history, 'og-2', 30);
  app.drawingList.newDrawing('linear');
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'a new drawing starts with none');
  await history.runUndoable('Add File', () => app.addLinearSeq());
  show(app, aligned);
  assert.deepEqual(shownAlignment(), { plan: 'og-2', translations: [30] });
});

test('New drawing with this style starts without an alignment and leaves the source\'s', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  await align(history, 'og-3', 7);
  assert.equal(app.drawingList.newDrawingWithStyle(source, 'linear').status, 'ok');
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] });
  show(app, source);
  assert.deepEqual(shownAlignment(), { plan: 'og-3', translations: [7] });
});

test('Delete holds a drawing\'s alignment and its Undo brings it back', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: kept } = app.drawingList.newDrawing('linear');
  const { drawingId: deleted } = app.drawingList.newDrawing('linear');
  await align(history, 'og-4', 21);
  assert.equal(app.drawingList.deleteDrawing(deleted).status, 'ok');
  show(app, kept);
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'the drawing shown after has none');
  await history.undo();
  assert.equal(state.activeDrawingId.value, deleted);
  assert.deepEqual(shownAlignment(), { plan: 'og-4', translations: [21] });
});

// Reset this drawing clears the shown drawing's plan; Reset all drawings each
// drawing's, a hidden one's in its slot, and one Undo restores each.
test('Reset clears the plan of each drawing it resets and its Undo restores each', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: first } = app.drawingList.newDrawing('linear');
  await align(history, 'og-5', 3);
  const { drawingId: second } = app.drawingList.newDrawing('linear');
  await align(history, 'og-6', 4);

  await app.resetShownDrawing();
  assert.equal(shownAlignment().plan, null, 'Reset this drawing clears the shown plan');
  show(app, first);
  assert.deepEqual(shownAlignment(), { plan: 'og-5', translations: [3] }, 'and leaves the other drawing\'s');
  show(app, second);
  assert.equal(await app.undoHistory(), true);
  assert.equal(shownAlignment().plan, 'og-6');

  show(app, first);
  await app.resetSettings();
  assert.equal(shownAlignment().plan, null);
  show(app, second);
  assert.equal(shownAlignment().plan, null, 'Reset all drawings clears the hidden drawing\'s plan in its slot');
  assert.equal(await app.undoHistory(), true);
  assert.equal(shownAlignment().plan, 'og-6', 'one Undo restores the hidden drawing\'s');
  show(app, first);
  assert.equal(shownAlignment().plan, 'og-5', 'and the shown one\'s');
});

// OV-381 on a switch: the arriving drawing's plan selects its group in the
// Orthogroups drawer, by the switch's own step (no watcher writes state on a
// switch), also between a drawing and its Duplicate, whose plan is the same object.
test('a switch selects the arriving plan\'s group, also between a drawing and its Duplicate', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  await align(history, 'og-7', 5);
  const { drawingId: copy } = app.drawingList.duplicateDrawing(source);
  assert.equal(state.selectedOrthogroupId.value, 'og-7', 'the Duplicate shows its plan\'s group');
  show(app, source);
  assert.equal(state.selectedOrthogroupId.value, 'og-7', 'and so does its source');
  app.drawingList.newDrawing('linear');
  assert.equal(state.selectedOrthogroupId.value, '', 'a drawing without a plan selects none');
  show(app, copy);
  assert.equal(state.selectedOrthogroupId.value, 'og-7');
});

// D-38 (PD-OI-039): an open Similarity alignment review belongs to the drawing
// it was opened in and stays open across a switch with its choices; the
// "Similarity alignment applied" summary was the departing drawing's.
test('a switch keeps an open Similarity alignment review and drops the summary', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: aligned } = app.drawingList.newDrawing('linear');
  await align(history, 'og-8', 6);
  const { drawingId: other } = app.drawingList.newDrawing('linear');
  show(app, aligned);
  const draft = { groupId: 'og-8', choices: [] };
  app.similarityAlignmentStatus.value = 'reviewing';
  app.similarityAlignmentDraft.value = draft;
  app.similarityAlignmentSummary.value = { text: 'Similarity alignment applied' };
  show(app, other);
  assert.equal(app.similarityAlignmentStatus.value, 'reviewing');
  assert.equal(app.similarityAlignmentDraft.value, draft, 'the review keeps its choices');
  assert.equal(app.similarityAlignmentSummary.value, null);
  assert.deepEqual(shownAlignment(), { plan: null, translations: [] }, 'the shown drawing keeps its own (none)');
  show(app, aligned);
  assert.equal(app.similarityAlignmentDraft.value, draft);
  assert.deepEqual(shownAlignment(), { plan: 'og-8', translations: [6] }, 'the committed plan stays the drawing\'s');
});

// D-40: a switch, and a step of the drawing list, wait for the Similarity alignment work in flight, through the
// one switch-wait owner, with the reason the Apply line gives.
test('a switch waits while a Similarity alignment resolves or applies', async () => {
  const app = createAppSetup();
  const { drawingId: first } = app.drawingList.newDrawing('linear');
  app.drawingList.newDrawing('linear');
  app.similarityAlignmentStatus.value = 'resolving';
  assert.deepEqual(app.selectDrawing(first), { status: 'busy', reason: 'Resolving the current alignment.' });
  app.similarityAlignmentStatus.value = 'applying';
  assert.deepEqual(app.selectDrawing(first), { status: 'busy', reason: 'Applying the alignment plan.' });
  assert.equal(app.drawingList.newDrawing('linear').status, 'busy');
  app.similarityAlignmentStatus.value = 'idle';
  assert.deepEqual(app.selectDrawing(first), { status: 'ok' });
});

// D-41 (S1): an Undo or Redo of a drawing list step switches the shown drawing inside History's traversal, so it
// waits for the Similarity alignment work in flight through the switch's own wait, with its reason; so does a pick in
// the History list that walks through such a step. Once the work ends, the step is taken.
test('an Undo or Redo of a drawing list step waits while a Similarity alignment resolves', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const before = state.activeDrawingId.value;
  const { drawingId: added } = app.drawingList.newDrawing('linear');
  const busy = { status: 'busy', reason: 'Resolving the current alignment.' };
  app.similarityAlignmentStatus.value = 'resolving';
  assert.deepEqual(await history.undo(), busy, 'the Undo of New waits');
  assert.deepEqual(await history.travelTo(null), busy, 'and so does a History list pick through it');
  assert.equal(state.activeDrawingId.value, added, 'the shown drawing stays');
  assert.ok(state.drawingTabs.value.some((tab) => tab.id === added), 'and so does the drawing');
  app.similarityAlignmentStatus.value = 'idle';
  assert.equal(await history.undo(), true);
  assert.equal(state.activeDrawingId.value, before);
  app.similarityAlignmentStatus.value = 'applying';
  assert.deepEqual(await history.redo(), { status: 'busy', reason: 'Applying the alignment plan.' }, 'the Redo waits');
  assert.equal(state.activeDrawingId.value, before);
  app.similarityAlignmentStatus.value = 'idle';
  assert.equal(await history.redo(), true);
  assert.equal(state.activeDrawingId.value, added);
});

// D-44 (S-1, R11): Undo and Redo and their enabled state read one check of the step they would take. The
// reviewer's sequence: Duplicate, then a review resolving in the copy: the next Undo is the Duplicate, a drawing list
// step that waits for the alignment work, so the header's Undo is disabled with the work's reason, and so are the
// History list's rows that walk through it. Once the work ends, both are enabled again.
test('Undo and Redo are disabled with the reason while the drawing list step they would take waits', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  app.drawingList.duplicateDrawing(source);
  const undoRows = () => app.historyList.rows.value.filter((row) => row.direction === 'undo' && row.steps > 0);
  app.similarityAlignmentStatus.value = 'resolving';
  assert.equal(app.canUndoHistory.value, false, 'the Undo of the Duplicate waits');
  assert.equal(app.undoHistoryTitle.value, 'Resolving the current alignment.');
  assert.ok(undoRows().length > 0 && undoRows().every((row) => app.historyList.rowDisabled(row)), 'and so do the rows');
  app.similarityAlignmentStatus.value = 'idle';
  assert.equal(app.canUndoHistory.value, true);
  assert.match(app.undoHistoryTitle.value, /^Undo /);
  assert.equal(await history.undo(), true);
  app.similarityAlignmentStatus.value = 'applying';
  assert.equal(app.canRedoHistory.value, false, 'the Redo of the Duplicate waits');
  assert.equal(app.redoHistoryTitle.value, 'Applying the alignment plan.');
  app.similarityAlignmentStatus.value = 'idle';
  assert.equal(app.canRedoHistory.value, true);
});

// D-57 (N-a): a History list row is disabled, with the reason as the list's hint, when any step of its walk would
// wait, read from the same one check the router uses; a pick of it takes no step, so no walk stops half way. The
// reviewer's sequence: Duplicate, a setting changed in the copy, then a review resolving: the "Duplicate drawing" row
// walks the setting step and then the Duplicate, which waits.
test('a History list row is disabled when any step of its walk waits, and its pick undoes nothing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  app.drawingList.duplicateDrawing(source);
  await drag(history, 7);
  const duplicateRow = app.historyList.rows.value.find((row) => row.direction === 'undo' && row.steps === 2);
  assert.ok(duplicateRow, 'the Duplicate is the second step back');
  const settingRow = app.historyList.rows.value.find((row) => row.direction === 'undo' && row.steps === 1);
  app.similarityAlignmentStatus.value = 'resolving';
  assert.equal(app.historyList.rowDisabled(settingRow), false, 'the setting step alone waits for nothing');
  assert.equal(app.historyList.rowDisabled(duplicateRow), true, 'the walk through the Duplicate waits');
  app.historyList.pointed.value = duplicateRow.key;
  assert.equal(app.historyList.hint.value, 'Resolving the current alignment.');
  assert.deepEqual(await history.travelTo(duplicateRow.target), { status: 'busy', reason: 'Resolving the current alignment.' });
  assert.deepEqual(shownAlignment().translations, [7], 'the pick undid nothing');
  app.similarityAlignmentStatus.value = 'idle';
  assert.equal(app.historyList.rowDisabled(duplicateRow), false);
  assert.match(app.historyList.hint.value, /^Undo 2 steps$/);
});

// D-44 (S-3, N-1): the Load rollback port's commit member replaces the History router's members only; the alignment
// review of a drawing the Load did not keep ends through its own member, which the Load calls once it succeeds.
test('the Load rollback port ends the review only on a Load that succeeded', async () => {
  const { createSessionImportRollbackState } = await import('../../gbdraw/web/js/app/app-setup.js');
  const replaced = [];
  const ended = [];
  const port = createSessionImportRollbackState({
    artifactSlots: /** @type {any} */ ({}), captureHistoryMembers: () => ({}), restoreHistoryMembers: () => {},
    replaceHistoryMembers: (ids) => replaced.push([...ids]),
    captureDisplayedArtifact: () => /** @type {any} */ ({}), installDisplayedArtifact: () => {},
    featureListScrollTop: { value: 0 }, featureListScrollRef: { value: null }, selectedPairwiseBlockOrthogroupId: { value: '' },
    alignment: { endReviewOutside: (ids) => ended.push([...ids]) }
  });
  port.replaceHistoryMembers(['linear-2']);
  assert.deepEqual({ replaced, ended }, { replaced: [['linear-2']], ended: [] }, 'the commit ends no review');
  port.endWorkOutside(['linear-2']);
  assert.deepEqual(ended, [['linear-2']]);
});

// D-44 (N-2): whether the review is the shown drawing's has one derivation, the alignment owner's (reviewIsShown);
// the review's hold (the root) and its marks (the review owner) read it, and compare no drawing IDs again.
test('only the alignment owner compares the review\'s drawing with the shown drawing', async () => {
  const { readFile } = await import('node:fs/promises');
  for (const file of ['app-setup.js', 'similarity-alignment-review.js']) {
    const source = await readFile(new URL(`../../gbdraw/web/js/app/${file}`, import.meta.url), 'utf8');
    assert.doesNotMatch(source, /[Rr]eviewDrawingId\(\)\s*[!=]==\s*state\.activeDrawingId\.value/, file);
  }
});

// D-41 (B1, OV-3618): a text that names a drawing reads its tab's name without throwing for an ID that has no tab.
test('the drawing bar names a drawing by its tab and has no name for an ID without one', () => {
  const app = createAppSetup();
  const { drawingId } = app.drawingList.newDrawing('linear');
  assert.equal(typeof app.drawingBar.nameOf(drawingId), 'string');
  assert.equal(app.drawingBar.nameOf('linear-99'), null);
});
