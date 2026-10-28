// The drawings of a project beyond the two mode drawings (PD-OI-090, B4).
// Reset resets each drawing by its own mode (R17-5): a drawing's reset leaves
// the other drawings and the one-per-Session values alone; the project part
// resets those values, as 0.14.0 Reset Settings did.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { withBootPalettes } from './helpers/boot-palettes.mjs';

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

const { state, createDrawingState, drawingHasEdits } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { PALETTE_SLICE } = await import('../../gbdraw/web/js/services/drawing/palette.js');
const { resetProject, resetSettings, withNewDrawingDefaults } = await import('../../gbdraw/web/js/services/reset.js');
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

test('New takes the defaults of its mode, the drawing Load and the generator make (OV-1761)', () => withBootPalettes(state, async () => {
  const app = createAppSetup();
  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const outcome = app.drawingList.newDrawing(mode);
    assert.equal(outcome.status, 'ok');
    const made = state.drawingOf(outcome.drawingId);
    const expected = withNewDrawingDefaults(state, createDrawingState(outcome.drawingId, mode));
    const colors = PALETTE_SLICE.capture(expected.slices.palette).currentColors;
    assert.ok(Object.keys(colors).length > 0, 'the default palette has colors');
    assert.deepEqual(PALETTE_SLICE.capture(made.slices.palette), PALETTE_SLICE.capture(expected.slices.palette),
      `a new ${mode} drawing takes its default palette's colors`);
    assert.equal(drawingHasEdits(made), false, 'it carries no edits');
  }
}));

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

// Delete (R16-2, R16-3, R17-6): one project step that holds the drawing, its
// Result slot, its History, and the files only it binds until it is evicted.
const { FEATURE_EDITS_SLICE } = await import('../../gbdraw/web/js/services/drawing/feature-edits.js');
/** @param {string} name */
const resultOf = (name) => ({ name, content: `<svg id="${name}"/>` });

test('Delete shows the next tab, else the previous; its Undo restores the drawing and shows it', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: first } = app.drawingList.newDrawing('circular');
  const { drawingId: second } = app.drawingList.newDrawing('linear');
  assert.deepEqual(app.selectDrawing(first), { status: 'ok' });
  const drawing = state.drawingOf(first);
  await history.runUndoable('Plot title', () => { drawing.form.plot_title = 'Deleted title'; });
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featureOverrides: { 'feature-1': { visible: false } } });
  state.results.value = [resultOf('first.svg')];
  const order = tabIds();
  const name = tabNameOf(first);

  assert.deepEqual(app.drawingList.deleteDrawing(first), { status: 'ok', drawingId: first });
  assert.equal(state.activeDrawingId.value, second, 'deleting the shown drawing shows the next tab');
  assert.deepEqual(tabIds(), order.filter((id) => id !== first));
  assert.throws(() => state.drawingOf(first));
  assert.throws(() => history.managerOf(first), 'its History is held by the step');
  assert.equal(history.undoLabel(), 'Delete drawing');

  await history.undo();
  assert.deepEqual(tabIds(), order, 'Undo puts its tab back where it was');
  assert.equal(tabNameOf(first), name);
  assert.equal(state.activeDrawingId.value, first, 'Undo shows the deleted drawing');
  assert.equal(state.drawingOf(first), drawing, 'with its settings and edits');
  assert.equal(drawing.form.plot_title, 'Deleted title');
  assert.deepEqual(Object.keys(drawing.featureOverrides), ['feature-1']);
  assert.deepEqual(state.results.value.map((result) => result.name), ['first.svg'], 'and its Result');
  assert.equal(history.managerOf(first).getUndoCount(), 1, 'and its History');

  await history.redo();
  assert.equal(state.activeDrawingId.value, second, 'Redo deletes it again and shows the drawing shown after');
  assert.deepEqual(tabIds(), order.filter((id) => id !== first));

  // The last tab has no next one: the previous tab is shown.
  const previous = tabIds()[tabIds().indexOf(second) - 1];
  app.drawingList.deleteDrawing(second);
  assert.equal(state.activeDrawingId.value, previous, 'deleting the last tab shows the previous one');
  await history.undo();
  assert.equal(state.activeDrawingId.value, second);
  await history.undo();
  assert.equal(state.activeDrawingId.value, first, 'each Undo restores its own drawing');
  assert.deepEqual(state.results.value.map((result) => result.name), ['first.svg']);
});

test('evicting the Delete step ends its Undo and releases the files only the deleted drawing binds', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(drawingId);
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'only-here.fa')]);
  const { fileId } = drawing.files.c_conservation_fastas[0];
  app.drawingList.deleteDrawing(drawingId);
  const kept = state.activeDrawingId.value;
  assert.equal(state.inputFiles.has(fileId), true, 'the held drawing keeps its file');

  // The step-count window evicts the Delete step after 30 newer project steps.
  for (let step = 0; step < 29; step += 1) app.drawingList.renameDrawing(kept, `Kept ${step}`);
  assert.equal(state.inputFiles.has(fileId), true, 'the file stays while the step is held');
  app.drawingList.renameDrawing(kept, 'Kept last');
  assert.equal(state.inputFiles.has(fileId), false, 'evicting the step releases the file');
  while (history.canUndo()) {
    assert.notEqual(history.undoLabel(), 'Delete drawing');
    await history.undo();
  }
  assert.equal(tabIds().includes(drawingId), false, 'no Undo brings the deleted drawing back');
});

test('the last drawing cannot be deleted; Delete asks first only for a drawing with a Result or edits', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(drawingId);
  assert.equal(app.drawingList.needsDeleteConfirmation(drawingId), false, 'a new drawing');
  drawing.form.plot_title = 'Changed setting';
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'ring.fa')]);
  assert.equal(app.drawingList.needsDeleteConfirmation(drawingId), false, 'files and changed settings only');
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featureOverrides: { 'feature-1': { visible: false } } });
  assert.equal(app.drawingList.needsDeleteConfirmation(drawingId), true, 'an edit');
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featureOverrides: {} });
  state.results.value = [resultOf('shown.svg')];
  assert.equal(app.drawingList.needsDeleteConfirmation(drawingId), true, 'a Result, shown');
  app.selectDrawing('circular');
  assert.equal(app.drawingList.needsDeleteConfirmation(drawingId), true, 'a Result, kept in its slot');
  assert.throws(() => app.drawingList.needsDeleteConfirmation('no-such-drawing'), listRefusal('UNKNOWN_DRAWING'));

  const tabs = state.drawingTabs.value;
  const drawings = state.drawingList();
  state.replaceDrawings([state.drawingOf('circular')], [{ id: 'circular', name: 'Only' }]);
  history.replaceMembers(['circular']);
  assert.throws(() => app.drawingList.deleteDrawing('circular'), listRefusal('LAST_DRAWING'));
  assert.deepEqual(tabIds(), ['circular'], 'a refused Delete changes nothing');
  assert.equal(history.canUndo(), false, 'and records nothing');
  state.replaceDrawings(drawings, tabs);
  history.replaceMembers(drawings.map((each) => each.id));
});

// The drawing bar (B4a-5): one tab per drawing in tab order, a default name in
// the UI language, and Rename and Move through the list operations. A refused
// name asks again with the error owner's wording (R6), never a matched message.
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { focusedTabIndex } = await import('../../gbdraw/web/js/app/drawing-list.js');

test('the drawing bar lists the tabs; Rename asks again on a refused name and Undo restores the name', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const shown = state.activeDrawingId.value;
  assert.deepEqual(app.drawingBar.tabs().map(({ id, label, mode, selected }) => ({ id, label, mode, selected })),
    state.drawingTabs.value.map((tab) => ({ id: tab.id, label: tab.name, mode: state.drawingOf(tab.id).mode,
      selected: tab.id === shown })));
  assert.equal(app.drawingBar.tabs()[0].label, defaultDrawingName(state.drawingOf(tabIds()[0]).mode, tabIds()[0]));

  const name = tabNameOf('linear');
  const other = tabNameOf('circular');
  const asked = () => ({ show: app.choicePrompt.show, value: app.choicePrompt.field?.value, description: app.choicePrompt.description });
  const answer = async (value) => {
    await new Promise((resolve) => setTimeout(resolve, 0));
    app.choicePrompt.field.value = value;
    app.resolveChoicePrompt('ok');
    await new Promise((resolve) => setTimeout(resolve, 0));
  };
  const renaming = app.drawingBar.renameDrawing('linear');
  await new Promise((resolve) => setTimeout(resolve, 0));
  assert.deepEqual(asked(), { show: true, value: name, description: '' }, 'the dialog starts from the stored name');
  await answer(other.toUpperCase());
  const refusal = normalizeUserFacingError(Object.assign(new Error('REPEATED_NAME'), {
    code: 'INPUT_INVALID', stage: 'drawing-list', context: { field: 'drawing', reason: 'REPEATED_NAME' } }));
  assert.equal(refusal.stage, 'drawing-list');
  assert.match(refusal.summary, /Another drawing has this name/);
  assert.deepEqual(asked(), { show: true, value: other.toUpperCase(), description: refusal.summary },
    'a repeated name asks again with the reason');
  await answer('   ');
  assert.match(app.choicePrompt.description, /Enter a name for the drawing/);
  await answer('Phage comparison');
  assert.equal((await renaming).status, 'ok');
  assert.equal(tabNameOf('linear'), 'Phage comparison');
  assert.equal(app.drawingBar.tabs().find((tab) => tab.id === 'linear').label, 'Phage comparison');
  assert.equal(history.undoLabel(), 'Rename drawing');
  await history.undo();
  assert.equal(tabNameOf('linear'), name, 'Undo restores the name');

  const cancelled = app.drawingBar.renameDrawing('linear');
  await new Promise((resolve) => setTimeout(resolve, 0));
  app.resolveChoicePrompt('cancel');
  assert.equal(await cancelled, null, 'Cancel renames nothing');

  const order = tabIds();
  const [first, second] = app.drawingBar.tabs();
  assert.deepEqual([first.canMoveLeft, first.canMoveRight, second.canMoveRight], [false, true, order.length > 2]);
  assert.equal(app.drawingBar.moveDrawing(second.id, -1).status, 'ok');
  assert.deepEqual(tabIds(), [order[1], order[0], ...order.slice(2)]);
  await history.undo();
  assert.deepEqual(tabIds(), order, 'Undo restores the order');
});

test('bar Delete asks once for a drawing with a Result and Cancel records no step; a bare drawing goes without asking', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const settle = () => new Promise((resolve) => setTimeout(resolve, 0));
  const { drawingId: bare } = app.drawingList.newDrawing('linear');
  const { drawingId: kept } = app.drawingList.newDrawing('circular');
  state.results.value = [resultOf('kept.svg')];
  const order = tabIds();
  const undoCount = history.getUndoCount();
  assert.equal(app.drawingBar.tabs().every((tab) => tab.canDelete), true, 'with several drawings every tab can be deleted');

  const cancelled = app.drawingBar.deleteDrawing(kept);
  await settle();
  assert.equal(app.choicePrompt.show, true, 'a drawing with a Result asks first');
  assert.match(app.choicePrompt.title, /^Delete /);
  app.resolveChoicePrompt('cancel');
  assert.equal(await cancelled, null);
  assert.deepEqual(tabIds(), order, 'Cancel deletes nothing');
  assert.equal(history.getUndoCount(), undoCount, 'and records no step');

  const confirmed = app.drawingBar.deleteDrawing(kept);
  await settle();
  app.resolveChoicePrompt('delete');
  assert.equal((await confirmed).status, 'ok');
  assert.deepEqual(tabIds(), order.filter((id) => id !== kept));
  await history.undo();
  assert.equal(state.activeDrawingId.value, kept, 'Undo shows the deleted drawing');
  assert.deepEqual(state.results.value.map((result) => result.name), ['kept.svg'], 'with its Result');

  state.results.value = [];
  app.selectDrawing(bare);
  const outcome = await app.drawingBar.deleteDrawing(bare);
  assert.equal(outcome.status, 'ok', 'a drawing without Result or edits is deleted without asking');
  assert.equal(app.choicePrompt.show, false);
  await history.undo();
});

test('arrow keys move the focus over the tabs and wrap; Home and End go to the ends', () => {
  assert.equal(focusedTabIndex('ArrowRight', 0, 3), 1);
  assert.equal(focusedTabIndex('ArrowRight', 2, 3), 0);
  assert.equal(focusedTabIndex('ArrowLeft', 0, 3), 2);
  assert.equal(focusedTabIndex('Home', 2, 3), 0);
  assert.equal(focusedTabIndex('End', 0, 3), 2);
  assert.equal(focusedTabIndex('Enter', 1, 3), null, 'Enter and Space select through the button (manual activation)');
});
