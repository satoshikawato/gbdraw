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
