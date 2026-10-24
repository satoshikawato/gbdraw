// The drawings of a project beyond the two mode drawings (PD-OI-090, B4).
// Reset resets each drawing by its own mode (R17-5): a drawing's reset leaves
// the other drawings and the one-per-Session values alone; the project part
// resets those values, as 0.14.0 Reset Settings did.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null };
installFakeSvgDom();
globalThis.File = File;

const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { resetProject, resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
const { createDefaultLosatExecution } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { createDefaultEditorDraftState, createDefaultPaletteDraftState } = await import('../../gbdraw/web/js/state.js');

// The one-per-Session values (R17-5): Reset this drawing keeps them.
const sessionValues = () => JSON.parse(JSON.stringify({
  losatExecution: state.losatExecution,
  richFeaturePopup: state.richFeaturePopup.value,
  autoLabelReflowEnabled: state.autoLabelReflowEnabled.value,
  downloadDpi: state.downloadDpi.value,
  paletteInstantPreviewEnabled: state.paletteInstantPreviewEnabled.value
}));
const editSessionValues = () => {
  state.losatExecution.threadsPerJob = 2;
  state.richFeaturePopup.value = false;
  state.autoLabelReflowEnabled.value = true;
  state.downloadDpi.value = 72;
  state.paletteInstantPreviewEnabled.value = true;
};
const DEFAULT_SESSION_VALUES = {
  losatExecution: createDefaultLosatExecution(),
  richFeaturePopup: true,
  autoLabelReflowEnabled: false,
  downloadDpi: createDefaultEditorDraftState().downloadDpi,
  paletteInstantPreviewEnabled: createDefaultPaletteDraftState().paletteInstantPreviewEnabled
};

test('a drawing\'s reset leaves the other drawing and the Session values; the project reset restores those', () => {
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  circular.form.plot_title = 'Circular title';
  linear.form.plot_title = 'Linear title';
  editSessionValues();
  const edited = sessionValues();
  assert.notDeepEqual(edited, JSON.parse(JSON.stringify(DEFAULT_SESSION_VALUES)));

  const artifact = resetSettings(state, circular);
  assert.equal(circular.form.plot_title, '');
  assert.equal(linear.form.plot_title, 'Linear title', 'the other drawing keeps its settings');
  assert.deepEqual(sessionValues(), edited, 'the one-per-Session values stay');
  assert.equal(artifact.appliedPaletteName, createDefaultPaletteDraftState().appliedPaletteName);

  resetProject(state);
  assert.deepEqual(sessionValues(), JSON.parse(JSON.stringify(DEFAULT_SESSION_VALUES)));
  assert.equal(linear.form.plot_title, 'Linear title', 'the project reset writes no drawing');
  linear.form.plot_title = '';
});

// The follow-ups run for every drawing by its mode, not for the drawing of a
// mode: a third drawing gets its Record Order rows and its comparison rings.
test('Reset all gives a second Circular drawing its Record Order rows and comparison rings', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const second = createDrawingState('circular-2', 'circular');
  state.addDrawing(second, 'Circular 2', 2);
  history.addMember('circular-2');
  second.circularRecordList.value = [{ selector: '#1', record_id: 'NC_1' }, { selector: '#2', record_id: 'NC_2' }];
  second.adv.multi_record_positions.splice(0, Infinity, { selector: '#2', row: 3 });
  bindDrawingInputFileTree(state.inputFiles, second.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'ring.fa')]);

  await app.resetSettings();
  assert.deepEqual(second.adv.multi_record_positions.map((position) => position.selector), ['#1', '#2'],
    'one Record Order row per loaded record, at its default');
  assert.equal(second.circularConservation.enabled, true, 'the kept ring file gives a ring');
  assert.equal(second.circularConservation.series.length, 1);
  assert.equal(second.files.c_conservation_fastas.length, 1, 'the ring file stays');
});

// A search belongs to the drawing it was typed in, not to its mode: a second
// Linear drawing starts without the first one's search.
test('a search typed in one Linear drawing is not shown in a second Linear drawing', () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const second = createDrawingState('linear-2', 'linear');
  state.addDrawing(second, 'Linear 2', 3);
  history.addMember('linear-2');
  app.selectDrawing('linear');
  app.previewFeatureSearchInput.value = 'dnaA';
  app.applyPreviewFeatureSearch();
  assert.equal(app.previewFeatureSearchQuery.value, 'dnaA');

  app.selectDrawing('linear-2');
  assert.equal(app.previewFeatureSearchQuery.value, '', 'the second drawing has no search');
  assert.equal(app.previewFeatureSearchInput.value, '');

  app.selectDrawing('linear');
  assert.equal(app.previewFeatureSearchQuery.value, 'dnaA', 'the first drawing keeps its own search');
});

// The operations on the drawing list (R11, R16-2): each is one project step
// of `app/drawing-list.js`, and a new drawing starts with an empty History.
const { allocateDrawingId, defaultDrawingName } = await import('../../gbdraw/web/js/services/drawing-view.js');
const tabIds = () => state.drawingTabs.value.map((tab) => tab.id);
const tabNameOf = (drawingId) => state.drawingTabs.value.find((tab) => tab.id === drawingId)?.name;
/** @param {string} reason */
const listRefusal = (reason) => (error) => (
  error.code === 'INPUT_INVALID' && error.stage === 'drawing-list' && error.context.reason === reason);

test('New is one project step: its Undo shows the drawing shown before, its Redo the new drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  app.selectDrawing('linear');
  const before = tabIds();
  const outcome = app.drawingList.newDrawing('circular');
  assert.equal(outcome.status, 'ok');
  const { drawingId } = outcome;
  assert.equal(drawingId, allocateDrawingId('circular', before), 'the ID of its kind over the live drawings');
  assert.deepEqual(tabIds(), [...before, drawingId], 'its tab comes last');
  assert.equal(tabNameOf(drawingId), defaultDrawingName('circular', drawingId));
  assert.equal(state.activeDrawingId.value, drawingId, 'the new drawing is shown');
  const drawing = state.drawingOf(drawingId);
  assert.equal(history.undoLabel(), 'New drawing');

  await history.undo();
  assert.equal(state.activeDrawingId.value, 'linear', 'Undo shows the drawing shown before');
  assert.deepEqual(tabIds(), before);
  assert.throws(() => history.managerOf(drawingId), 'its History left the router with it');

  await history.redo();
  assert.equal(state.activeDrawingId.value, drawingId, 'Redo shows the new drawing');
  assert.equal(state.drawingOf(drawingId), drawing, 'the same drawing comes back');
  assert.deepEqual(tabIds(), [...before, drawingId]);
  assert.equal(history.managerOf(drawingId).getUndoCount(), 0);
});

test('a new drawing starts with an empty History', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  app.selectDrawing('circular');
  await history.runUndoable('Plot title', () => { modeDrawing(state, 'circular').form.plot_title = 'Edited'; });
  assert.equal(history.getUndoCount(), 1, 'the shown drawing has a step');
  const outcome = app.drawingList.newDrawing('linear');
  assert.equal(outcome.status, 'ok');
  const { drawingId } = outcome;
  const own = history.managerOf(drawingId);
  assert.equal(own.getUndoCount(), 0);
  assert.equal(own.getRedoCount(), 0);
  assert.equal(history.undoLabel(), 'New drawing', 'Undo from the new drawing takes the New step');
  modeDrawing(state, 'circular').form.plot_title = '';
});

test('Rename refuses a blank or repeated name; a rename and a move are one step each', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const shown = state.activeDrawingId.value;
  const name = tabNameOf('linear');
  const other = tabNameOf('circular');
  assert.throws(() => app.drawingList.renameDrawing('linear', ''), listRefusal('BLANK_NAME'));
  assert.throws(() => app.drawingList.renameDrawing('linear', '   '), listRefusal('BLANK_NAME'));
  assert.throws(() => app.drawingList.renameDrawing('linear', other.toUpperCase()), listRefusal('REPEATED_NAME'));
  assert.throws(() => app.drawingList.renameDrawing('no-such-drawing', 'X'), listRefusal('UNKNOWN_DRAWING'));
  assert.equal(tabNameOf('linear'), name, 'a refused rename changes nothing');
  assert.equal(history.canUndo(), false, 'and records nothing');

  assert.equal(app.drawingList.renameDrawing('linear', 'Phage comparison').status, 'ok');
  assert.equal(tabNameOf('linear'), 'Phage comparison');
  assert.equal(history.undoLabel(), 'Rename drawing');
  await history.undo();
  assert.equal(tabNameOf('linear'), name, 'Undo restores the name');
  await history.redo();
  assert.equal(tabNameOf('linear'), 'Phage comparison');

  const order = tabIds();
  app.drawingList.moveDrawing('linear', 0);
  assert.equal(tabIds()[0], 'linear');
  assert.equal(history.undoLabel(), 'Move drawing');
  await history.undo();
  assert.deepEqual(tabIds(), order, 'Undo restores the tab order');
  assert.equal(history.undoLabel(), 'Rename drawing', 'the move was one step');
  await history.undo();
  assert.equal(state.activeDrawingId.value, shown, 'neither step switches the shown drawing');
});
