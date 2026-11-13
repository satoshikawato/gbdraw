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
globalThis.document = {
  createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} }
};
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

// A new Session starts with one Circular drawing; the Linear offer tab makes
// the Linear one (PD-OI-090 〔drawing bar と切り替え〕, T-B4-5). The tests below
// work on that pair.
test('a new Session lists one Circular drawing; the Linear offer tab makes the Linear drawing as one step', async () => {
  assert.deepEqual(state.drawingTabs.value.map(({ id, name }) => ({ id, name })), [{ id: 'circular', name: 'Circular' }]);
  assert.equal(state.activeDrawingId.value, 'circular');
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const offerOf = () => app.drawingBar.tabs().filter((tab) => tab.kind === 'offer');
  assert.deepEqual(offerOf().map(({ mode, label }) => ({ mode, label })), [{ mode: 'linear', label: 'Linear' }],
    'the offer tab shows while the project has no Linear drawing');
  assert.deepEqual(app.drawingBar.tabs().map((tab) => tab.kind), ['drawing', 'offer'], 'after the drawing tabs');

  const outcome = app.drawingBar.newDrawing('linear');
  assert.deepEqual(outcome, { status: 'ok', drawingId: 'linear' });
  assert.deepEqual(state.drawingTabs.value.map(({ id, name }) => ({ id, name })),
    [{ id: 'circular', name: 'Circular' }, { id: 'linear', name: 'Linear' }]);
  assert.equal(state.activeDrawingId.value, 'linear', 'the new Linear drawing is shown');
  assert.equal(state.drawingOf('linear').mode, 'linear');
  assert.deepEqual(offerOf(), [], 'no offer tab once a Linear drawing exists');
  assert.equal(history.timeline().undo[0]?.label, 'New drawing');
  assert.equal(history.managerOf('linear').getUndoCount(), 0, 'its History starts empty');

  await history.undo();
  assert.deepEqual(state.drawingTabs.value.map((tab) => tab.id), ['circular'], 'Undo removes the Linear drawing');
  assert.equal(state.activeDrawingId.value, 'circular', 'and shows the drawing shown before');
  assert.equal(offerOf().length, 1, 'the offer tab is back');
  assert.equal(history.timeline().undo.length, 0, 'the step was the only one');

  await history.redo();
  assert.deepEqual(state.drawingTabs.value.map((tab) => tab.id), ['circular', 'linear']);
  assert.equal(state.activeDrawingId.value, 'linear');
  app.selectDrawing('circular');
});

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
  second.circularRecordList.value = [
    { uid: 'record-1', k: 1, identity: 'record-1:1', selector: '#1', record_id: 'NC_1' },
    { uid: 'record-1', k: 2, identity: 'record-1:2', selector: '#2', record_id: 'NC_2' }
  ];
  second.adv.multi_record_positions.splice(0, Infinity, { selector: 'record-1:2', row: 3 });
  bindDrawingInputFileTree(state.inputFiles, second.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'ring.fa')]);

  await app.resetSettings();
  assert.deepEqual(second.adv.multi_record_positions.map((position) => position.selector), ['record-1:1', 'record-1:2'],
    'one Record Order row per loaded record, at its default');
  // OV-2818: the kept ring file keeps its ring row; Show Comparison Rings is off, as after any Reset.
  assert.equal(second.circularConservation.enabled, false);
  assert.equal(second.circularConservation.series.length, 1, 'the kept ring file gives a ring row');
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

// OV-2151: a switch names an existing drawing; an unknown ID is a programming
// error, not a switch to nothing (a spec that "switched to Linear" tested Circular).
test('selecting an unknown drawing ID throws and changes nothing; the shown drawing returns false', () => {
  const app = createAppSetup();
  const shown = state.activeDrawingId.value;
  const order = tabIds();
  assert.throws(() => app.selectDrawing('no-such-drawing'), /Unknown drawing: no-such-drawing/);
  assert.equal(state.activeDrawingId.value, shown);
  assert.deepEqual(tabIds(), order);
  assert.equal(app.selectDrawing(shown), false);
});

// R17-29 (OV-2156): an error that names a drawing leaves with it on Delete,
// as it leaves on a switch away (OV-372); Undo of the Delete does not bring it
// back. E1's Generate (PD-OI-079) therefore always finds the drawing it names.
test('deleting a drawing clears an error that names it, and only that error', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: gone } = app.drawingList.newDrawing('linear');
  const { drawingId: other } = app.drawingList.newDrawing('linear');
  const { drawingId: shown } = app.drawingList.newDrawing('circular');
  const errorNaming = (drawingId) => ({ code: 'SESSION_SAVE_REQUIRES_GENERATE', actions: ['generate'], context: { drawingId } });
  app.errorLog.value = errorNaming(other);
  app.drawingList.deleteDrawing(gone);
  assert.deepEqual(app.errorLog.value, errorNaming(other), 'an error naming another drawing stays');
  await history.undo();
  assert.ok(tabIds().includes(gone));
  app.errorLog.value = errorNaming(gone);
  app.drawingList.deleteDrawing(gone);
  assert.equal(state.activeDrawingId.value, shown);
  assert.equal(app.errorLog.value, null, 'the error naming the deleted drawing is gone');
  await history.undo();
  assert.ok(tabIds().includes(gone));
  assert.equal(app.errorLog.value, null, 'Undo of the Delete does not bring the error back');
});

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
  assert.equal(history.timeline().undo[0]?.label, 'New drawing');

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
  assert.equal(history.timeline().undo[0]?.label, 'New drawing', 'Undo from the new drawing takes the New step');
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
  assert.equal(history.timeline().undo[0]?.label, 'Rename drawing');
  await history.undo();
  assert.equal(tabNameOf('linear'), name, 'Undo restores the name');
  await history.redo();
  assert.equal(tabNameOf('linear'), 'Phage comparison');

  const order = tabIds();
  app.drawingList.moveDrawing('linear', 0);
  assert.equal(tabIds()[0], 'linear');
  assert.equal(history.timeline().undo[0]?.label, 'Move drawing');
  await history.undo();
  assert.deepEqual(tabIds(), order, 'Undo restores the tab order');
  assert.equal(history.timeline().undo[0]?.label, 'Rename drawing', 'the move was one step');
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
  assert.equal(history.timeline().undo[0]?.label, 'Delete drawing');

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

// Duplicate (PD-OI-090 〔作成〕, R16-2): the copy takes the source's files,
// settings, edits, and Result, with an empty History, and shares no member.
test('Duplicate copies a drawing as one step; its Undo shows the source and its Redo the copy', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('linear');
  const drawing = state.drawingOf(source);
  await history.runUndoable('Plot title', () => { drawing.form.plot_title = 'Source title'; });
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featureOverrides: { 'feature-1': { visible: false } } });
  state.results.value = [resultOf('source.svg')];
  const before = tabIds();
  const name = tabNameOf(source);

  const outcome = app.drawingList.duplicateDrawing(source);
  assert.equal(outcome.status, 'ok');
  const { drawingId: copyId } = outcome;
  assert.equal(copyId, allocateDrawingId('linear', before), 'the ID of its kind over the live drawings');
  assert.deepEqual(tabIds(), [...before, copyId], 'its tab comes last');
  assert.equal(tabNameOf(copyId), `${name} copy`);
  assert.equal(state.activeDrawingId.value, copyId, 'the copy is shown');
  const copy = state.drawingOf(copyId);
  assert.equal(copy.mode, 'linear');
  assert.equal(copy.form.plot_title, 'Source title', 'with its settings');
  assert.deepEqual(copy.featureOverrides, drawing.featureOverrides, 'and its edits');
  assert.deepEqual(state.results.value.map((result) => result.name), ['source.svg'], 'and its Result, without Generate');
  assert.equal(history.managerOf(copyId).getUndoCount(), 0, 'with an empty History');
  for (const key of Object.keys(drawing)) {
    const member = /** @type {Record<string, unknown>} */ (drawing)[key];
    if (member && typeof member === 'object') {
      assert.notEqual(/** @type {Record<string, unknown>} */ (copy)[key], member, `the copy shares no ${key}`);
    }
  }
  copy.form.plot_title = 'Copy title';
  assert.equal(drawing.form.plot_title, 'Source title', 'a write to the copy leaves the source');
  copy.form.plot_title = 'Source title';
  assert.equal(history.timeline().undo[0]?.label, 'Duplicate drawing');

  await history.undo();
  assert.equal(state.activeDrawingId.value, source, 'Undo shows the source');
  assert.deepEqual(tabIds(), before);
  assert.deepEqual(state.results.value.map((result) => result.name), ['source.svg']);

  await history.redo();
  assert.equal(state.activeDrawingId.value, copyId, 'Redo shows the copy');
  assert.equal(state.drawingOf(copyId), copy, 'the same copy comes back');
  assert.deepEqual(state.results.value.map((result) => result.name), ['source.svg'], 'with its Result');
  assert.equal(app.drawingList.duplicateDrawing(source).status, 'ok');
  assert.equal(tabNameOf(tabIds().at(-1)), `${name} copy 2`, 'a second copy takes the next free name');
});

// OV-2477 (Owner D-35 A): a name the app makes is made in the UI language and
// stored as is. A copy takes the source tab's shown label and the translated
// copy phrase; Rename starts from the shown label. A default name stays stored
// in English and shown translated.
const { setUiLocale } = await import('../../gbdraw/web/js/utils/i18n.js');

test('in a non-English UI, Duplicate names the copy from the shown label and Rename starts from it', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('circular');
  const stored = tabNameOf(source);
  assert.equal(stored, defaultDrawingName('circular', source), 'a default name is stored in English');
  const catalog = { Circular: '環状', '{name} copy': '{name} のコピー', '{name} copy {number}': '{name} のコピー {number}' };
  await setUiLocale('ja', async () => ({ default: catalog }));
  try {
    const shown = defaultDrawingName('circular', source, '環状');
    const labelOf = (drawingId) => app.drawingBar.tabs().find((tab) => tab.kind === 'drawing' && tab.id === drawingId).label;
    assert.equal(labelOf(source), shown, 'and shown in the UI language');
    assert.equal(tabNameOf(source), stored);

    const first = app.drawingBar.duplicateDrawing(source);
    assert.equal(first.status, 'ok');
    assert.equal(tabNameOf(first.drawingId), `${shown} のコピー`, 'the copy is named from the shown label');
    assert.equal(labelOf(first.drawingId), `${shown} のコピー`, 'and shown as stored');
    const second = app.drawingBar.duplicateDrawing(source);
    assert.equal(tabNameOf(second.drawingId), `${shown} のコピー 2`, 'a second copy takes the numbered phrase');

    const renaming = app.drawingBar.renameDrawing(source);
    await new Promise((resolve) => setTimeout(resolve, 0));
    assert.equal(app.choicePrompt.field?.value, shown, 'the Rename field starts from the shown label');
    app.resolveChoicePrompt('cancel');
    assert.equal(await renaming, null);
  } finally {
    await setUiLocale('en');
    for (let step = 0; step < 3; step += 1) await history.undo();
  }
  assert.ok(!tabIds().includes(source), 'Undo removes the drawings this test made');
});

// New drawing with this style (PD-OI-090 〔作成〕(3), R17-4): a new drawing of
// the chosen mode at its defaults, with the source's settings that its mode
// has and the table rows the user entered; no files, edits, or card rows.
test('New drawing with this style copies the settings its mode has and the user rows, as one step', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId: source } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(source);
  const linearDefaults = withNewDrawingDefaults(state, createDrawingState('linear-9', 'linear'));
  drawing.form.show_scale = !linearDefaults.form.show_scale;
  drawing.adv.legend_font_size = 31;
  drawing.adv.feature_width_circular = linearDefaults.adv.feature_width_circular + 7;
  drawing.manualPriorityRules.splice(0, drawing.manualPriorityRules.length,
    { feat: 'CDS', order: 'gene' }, { feat: 'tRNA', order: 'product', importedFrom: 'qualifier_priority' });
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featureOverrides: { 'feature-1': { visible: false } } });
  const before = tabIds();

  const outcome = app.drawingList.newDrawingWithStyle(source, 'linear');
  assert.equal(outcome.status, 'ok');
  const { drawingId } = outcome;
  assert.equal(drawingId, allocateDrawingId('linear', before));
  assert.deepEqual(tabIds(), [...before, drawingId], 'its tab comes last');
  assert.equal(state.activeDrawingId.value, drawingId, 'it is shown');
  const styled = state.drawingOf(drawingId);
  assert.equal(styled.mode, 'linear');
  assert.equal(styled.form.show_scale, !linearDefaults.form.show_scale, 'a setting both modes have');
  assert.equal(styled.adv.legend_font_size, 31, 'a style both modes have');
  assert.equal(styled.adv.feature_width_circular, linearDefaults.adv.feature_width_circular,
    'a Circular-only setting keeps the Linear default');
  assert.deepEqual(JSON.parse(JSON.stringify(styled.manualPriorityRules)), [{ feat: 'CDS', order: 'gene' }],
    'the Priority row the user entered, not the one a card file loaded');
  assert.deepEqual(JSON.parse(JSON.stringify(styled.featureOverrides)), {}, 'no feature edit');
  assert.equal(history.managerOf(drawingId).getUndoCount(), 0, 'with an empty History');
  assert.equal(drawing.manualPriorityRules.length, 2, 'the source keeps its rows');
  assert.equal(history.timeline().undo[0]?.label, 'New drawing with this style');

  await history.undo();
  assert.equal(state.activeDrawingId.value, source, 'Undo shows the source');
  assert.deepEqual(tabIds(), before);
  await history.redo();
  assert.equal(state.activeDrawingId.value, drawingId, 'Redo shows the new drawing');
  assert.equal(state.drawingOf(drawingId), styled);
});

// OV-2215 (D-15-182): across modes the copy follows the registry's `migrate`
// rule, the owner Python's mode split uses: a row both modes share is copied,
// and a `profile` row or a row only the target mode has takes the target
// mode's profile.
test('New drawing with this style in the other mode takes that mode\'s profile rows and own rows', () => {
  const app = createAppSetup();
  const { drawingId: source } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(source);
  const linearDefaults = withNewDrawingDefaults(state, createDrawingState('linear-9', 'linear'));
  drawing.adv.label_font_size = linearDefaults.adv.label_font_size + 3;
  drawing.form.plot_title = 'Source title';
  drawing.adv.def_font_size = linearDefaults.adv.def_font_size + 5;

  const { drawingId } = app.drawingList.newDrawingWithStyle(source, 'linear');
  const styled = state.drawingOf(drawingId);
  assert.equal(styled.adv.label_font_size, linearDefaults.adv.label_font_size + 3, 'a `copy` row both modes share');
  assert.equal(styled.adv.pairwise_match_style, 'curve', 'a Linear-only row takes the Linear default');
  for (const field of ['pairwise_match_style', 'min_bitscore', 'evalue', 'identity', 'alignment_length', 'def_font_size',
    'axis_stroke_color', 'plot_title_font_size']) {
    assert.deepEqual(styled.adv[field], linearDefaults.adv[field], `${field} takes the Linear profile`);
  }
  assert.equal(styled.form.plot_title, linearDefaults.form.plot_title, 'plot_title takes the Linear profile');
});

// Figure size (issue #1427, PD-OI-114): the size is a `copy` row of both modes, so Duplicate and New drawing with
// this style, in either mode, carry both sides as typed; Reset Settings returns them to Auto (null).
test('Duplicate and New drawing with this style copy the Figure size; Reset Settings returns it to Auto', () => {
  const app = createAppSetup();
  const { drawingId: source } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(source);
  drawing.form.figure_width = '180mm';
  drawing.form.figure_height = '5in';
  const sizeOf = (id) => [state.drawingOf(id).form.figure_width, state.drawingOf(id).form.figure_height];
  const copies = [
    app.drawingList.duplicateDrawing(source).drawingId,
    app.drawingList.newDrawingWithStyle(source, 'circular').drawingId,
    app.drawingList.newDrawingWithStyle(source, 'linear').drawingId
  ];
  for (const id of copies) assert.deepEqual(sizeOf(id), ['180mm', '5in'], id);
  resetSettings(state, drawing);
  assert.deepEqual(sizeOf(source), [null, null], 'Reset Settings: Auto');
  assert.deepEqual(sizeOf(copies[0]), ['180mm', '5in'], 'the copies keep theirs');
});

// OV-2215 (D-15-182, R17-4): the records a drawing draws, its crop, and their
// display name records of files the copy does not bind, so a style copy in
// any mode leaves them at the new drawing's defaults.
test('New drawing with this style copies no record selection, crop, or record display', () => {
  const app = createAppSetup();
  const { drawingId: source } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(source);
  const circularDefaults = withNewDrawingDefaults(state, createDrawingState('circular-9', 'circular'));
  drawing.form.circular_record_selector = 'NC_000001';
  drawing.form.circular_region_start = 1000;
  drawing.form.circular_region_end = 5000;
  drawing.form.circular_record_label = 'Chr1';
  drawing.form.circular_record_subtitle = 'Main chromosome';
  drawing.form.circular_reverse = true;
  drawing.adv.label_font_size = circularDefaults.adv.label_font_size + 3;
  drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, { selector: 'NC_000001', row: 2 });
  drawing.recordsOff.splice(0, drawing.recordsOff.length, 'circular-1');
  drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length, {
    scope: 'circular', sourceUid: 'circular', selector: '#1', recordId: 'NC_000001', topologyOverride: null,
    startCoordinate: 25, reverseComplementOverride: null, anchorIntent: null
  });

  const { drawingId } = app.drawingList.newDrawingWithStyle(source, 'circular');
  const styled = state.drawingOf(drawingId);
  assert.equal(styled.adv.label_font_size, circularDefaults.adv.label_font_size + 3, 'a style row is copied');
  for (const field of ['circular_record_selector', 'circular_region_start', 'circular_region_end', 'circular_record_label',
    'circular_record_subtitle', 'circular_reverse']) {
    assert.equal(styled.form[field], circularDefaults.form[field], `${field} keeps the new drawing's default`);
  }
  const plain = (/** @type {unknown} */ value) => JSON.parse(JSON.stringify(value));
  assert.deepEqual(plain(styled.adv.multi_record_positions), plain(circularDefaults.adv.multi_record_positions));
  assert.deepEqual(plain(styled.recordsOff), plain(circularDefaults.recordsOff));
  assert.deepEqual(plain(styled.recordDisplayDrafts), plain(circularDefaults.recordDisplayDrafts));
});

// D-15-185 (R17-4): a series label taken from the source drawing's Depth file
// is that file's, so the copy, which binds no file, has the default series
// name and Show Depth off (no Depth source); a typed label is style and stays.
// The source's files decide it, not the shown drawing's.
test('New drawing with this style leaves a series label taken from a Depth file and Show Depth behind', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  /** @param {string} name */
  const depthFile = (name) => new File(['chr\t1\t3\n'], name);
  /** @param {string} label @param {() => unknown} edit */
  const step = (label, edit) => history.runUndoable(label, edit);
  const plain = (/** @type {unknown} */ value) => JSON.parse(JSON.stringify(value));
  /** @param {DrawingState} drawing */
  const slotLabels = (drawing) => [...drawing.adv.circular_track_slots, ...drawing.adv.linear_track_slots]
    .filter((slot) => slot.renderer === 'depth').map((slot) => slot.params?.legend_label ?? '');
  /** @param {string} sourceId @param {'circular' | 'linear'} mode */
  const copyOf = (sourceId, mode) => {
    // Another drawing is shown, so the copy reads the source's files.
    app.drawingList.newDrawing(mode);
    const outcome = app.drawingList.newDrawingWithStyle(sourceId, mode);
    assert.equal(outcome.status, 'ok');
    return state.drawingOf(outcome.drawingId);
  };
  /** @param {DrawingState} copy @param {string} direction */
  const assertLeftBehind = (copy, direction) => {
    assert.deepEqual(plain(copy.adv.depth_tracks.map((track) => track.label)), ['', 'Coverage'],
      `${direction}: the file label becomes the default series name; the typed label and the series count stay`);
    assert.equal(copy.form.show_depth, false, `${direction}: Show Depth follows the copy's Depth sources (none)`);
    assert.ok(!slotLabels(copy).includes('sampleA'), `${direction}: no Depth row keeps the file label`);
  };

  app.selectDrawing('circular');
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('sampleA.tsv')));
  await step('Depth track', () => app.addCircularDepthTrack());
  await step('Depth label', () => app.setDepthTrackLabel(1, 'Coverage'));
  const circular = state.activeDrawing();
  assert.equal(app.getDepthTrackLabel(0), 'sampleA');
  assert.equal(circular.form.show_depth, true);

  const sameMode = copyOf(circular.id, 'circular');
  assertLeftBehind(sameMode, 'Circular to Circular');
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('sampleB.tsv')));
  assert.equal(app.getDepthTrackLabel(0), 'sampleB', 'a file uploaded into the copy names its series');
  assertLeftBehind(copyOf(circular.id, 'linear'), 'Circular to Linear');

  app.selectDrawing('linear');
  const linear = state.activeDrawing();
  await step('Depth file', () => app.setLinearDepthFile(linear.linearSeqs[0], 0, depthFile('sampleA.tsv')));
  await step('Depth track', () => app.addLinearDepthTrack());
  await step('Depth label', () => app.setDepthTrackLabel(1, 'Coverage'));
  assert.equal(app.getDepthTrackLabel(0), 'sampleA');
  assert.equal(linear.form.show_depth, true);
  assertLeftBehind(copyOf(linear.id, 'circular'), 'Linear to Circular');
  assertLeftBehind(copyOf(linear.id, 'linear'), 'Linear to Linear');
});

test('evicting the Delete step ends its Undo and releases the files only the deleted drawing binds', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { drawingId } = app.drawingList.newDrawing('circular');
  const drawing = state.drawingOf(drawingId);
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'only-here.fa')]);
  const { fileId } = drawing.files.c_conservation_fastas[0][0];
  app.drawingList.deleteDrawing(drawingId);
  const kept = state.activeDrawingId.value;
  assert.equal(state.inputFiles.has(fileId), true, 'the held drawing keeps its file');

  // The step-count window evicts the Delete step after 30 newer project steps.
  for (let step = 0; step < 29; step += 1) app.drawingList.renameDrawing(kept, `Kept ${step}`);
  assert.equal(state.inputFiles.has(fileId), true, 'the file stays while the step is held');
  app.drawingList.renameDrawing(kept, 'Kept last');
  assert.equal(state.inputFiles.has(fileId), false, 'evicting the step releases the file');
  while (history.canUndo()) {
    assert.notEqual(history.timeline().undo[0]?.label, 'Delete drawing');
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
  const shownLabel = app.drawingBar.tabs().find((tab) => tab.kind === 'drawing' && tab.id === 'linear').label;
  assert.deepEqual(asked(), { show: true, value: shownLabel, description: '' }, 'the dialog starts from the shown label');
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
  assert.equal(history.timeline().undo[0]?.label, 'Rename drawing');
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

// The auto label of a Depth series belongs to its drawing (OV-2202): the text
// one drawing took from its file is not auto in another drawing, and an Undo
// or Reset this drawing in one drawing leaves the other's labels.
test('a Depth series auto label belongs to its drawing; Undo and Reset in one leave the other (OV-2202)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  /** @param {string} name */
  const depthFile = (name) => new File(['chr\t1\t3\n'], name);
  const shownLabel = () => app.getDepthTrackLabel(0);
  /** @param {string} label @param {() => unknown} edit */
  const step = (label, edit) => history.runUndoable(label, edit);
  app.selectDrawing('circular');
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('alpha.tsv')));
  assert.equal(shownLabel(), 'alpha', 'a file names an auto-labeled series');

  const { drawingId: second } = app.drawingBar.newDrawing('circular');
  await step('Depth label', () => app.setDepthTrackLabel(0, 'alpha'));
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('beta.tsv')));
  assert.equal(shownLabel(), 'alpha', 'a label typed in drawing B stays: the text is drawing A\'s auto label, not B\'s');

  await step('Depth label', () => app.setDepthTrackLabel(0, ''));
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('gamma.tsv')));
  assert.equal(shownLabel(), 'gamma');
  await history.undo();
  await app.resetShownDrawing();
  app.selectDrawing('circular');
  assert.equal(shownLabel(), 'alpha', 'B\'s Undo and Reset this drawing leave A\'s label');
  await step('Depth label', () => app.setDepthTrackLabel(0, 'gamma'));
  await step('Depth file', () => app.setCircularDepthFile(0, depthFile('delta.tsv')));
  assert.equal(shownLabel(), 'gamma', 'a label typed in A stays: B\'s undone auto label is not A\'s');
  app.drawingList.deleteDrawing(second);
});

// D-32: removing a Circular drawing's source retires its record presentation:
// every record-bound form field returns to its value in the one list
// (`CIRCULAR_RECORD_FORM_DEFAULTS`), which a style copy also leaves out. The
// retire writes the fields by name, so this pins that it names each of them.
test('removing the Circular source returns every record-bound form field to its default', async () => {
  const { CIRCULAR_RECORD_FORM_DEFAULTS } = await import('../../gbdraw/web/js/services/config.js');
  const app = createAppSetup();
  app.selectDrawing('circular');
  const drawing = state.activeDrawing();
  const { replaceCircularSources, createCircularSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  replaceCircularSources(drawing.slices.inputs, [createCircularSource({ uid: 'card-1' })]);
  app.setCircularSourceFile('card-1', 'gb', new File(['LOCUS'], 'record.gb'));
  const edited = Object.fromEntries(Object.entries(CIRCULAR_RECORD_FORM_DEFAULTS).map(([field, value]) => [
    field, typeof value === 'boolean' ? !value : typeof value === 'string' ? 'edited'
      : Array.isArray(value) ? [{ start: 90, end: 10, reverseComplement: true, column: null }] : 7
  ]));
  // The Record select names a record by its identity; the presentation retires with the card the select names.
  Object.assign(drawing.form, edited, { circular_record_selector: 'card-1' });
  app.setCircularSourceFile('card-1', 'gb', null);
  assert.deepEqual(Object.fromEntries(Object.keys(CIRCULAR_RECORD_FORM_DEFAULTS).map((field) => [field, drawing.form[field]])),
    { ...CIRCULAR_RECORD_FORM_DEFAULTS });
});

// PD-OI-044 scenario revision 5: a grid or batch crop goes with its card's file when the card is removed, or when
// its file is replaced and held more than one record; a replaced file of one record keeps it, as a one-record crop.
test('removing or replacing a Circular file drops or keeps the grid and batch crops of its records', async () => {
  const app = createAppSetup();
  app.selectDrawing('circular');
  const drawing = state.activeDrawing();
  const { replaceCircularSources, createCircularSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  replaceCircularSources(drawing.slices.inputs, ['card-1', 'card-2', 'card-3'].map((uid) => createCircularSource({ uid })));
  for (const uid of ['card-1', 'card-2', 'card-3']) app.setCircularSourceFile(uid, 'gb', new File(['LOCUS'], `${uid}.gb`));
  const range = (source, record) => ({ source, record, windows: [{ start: 11, end: 90, reverseComplement: false, column: null }] });
  drawing.form.circular_record_ranges = [range('card-1', ''), range('card-2', ''), range('card-3', 'B')];
  // The records each card was read with (a replacement starts a new read).
  const read = () => {
    drawing.circularRecordList.value = [
      { uid: 'card-2', k: 1, identity: 'card-2', record_id: 'A' },
      { uid: 'card-3', k: 1, identity: 'card-3:1', record_id: 'A' }, { uid: 'card-3', k: 2, identity: 'card-3:2', record_id: 'B' }
    ];
  };
  read();
  app.setCircularSourceFile('card-1', 'gb', null);
  read();
  app.setCircularSourceFile('card-2', 'gb', new File(['LOCUS'], 'card-2-new.gb'));
  assert.deepEqual(drawing.form.circular_record_ranges.map(({ source }) => source), ['card-2', 'card-3']);
  read();
  app.setCircularSourceFile('card-3', 'gb', new File(['LOCUS'], 'card-3-new.gb'));
  assert.deepEqual(drawing.form.circular_record_ranges.map(({ source }) => source), ['card-2']);
});

// A Circular grid or batch drawing that holds two records' crops, read from one card of two records, with Record select
// available (Multi-Record Canvas off, the card's records read).
const CARD_RECORDS = {
  'card-1': [
    { uid: 'card-1', k: 1, identity: 'card-1:1', selector: '#1', record_id: 'A', record_length: 100 },
    { uid: 'card-1', k: 2, identity: 'card-1:2', selector: '#2', record_id: 'B', record_length: 100 }
  ],
  'card-2': [{ uid: 'card-2', k: 1, identity: 'card-2', selector: '#1', record_id: 'C', record_length: 100 }]
};
// The drawing's cards are read: their records are listed and the discovery is current.
const readCards = async (drawing) => {
  const { circularSourceInputs } = await import('../../gbdraw/web/js/services/input-types.js');
  drawing.circularRecordList.value = drawing.circularSources.flatMap(({ uid }) => CARD_RECORDS[uid] ?? []);
  Object.assign(drawing.circularRecordDiscovery, { status: 'ready', error: '', inputType: drawing.inputType.value,
    sources: circularSourceInputs(drawing.inputType.value, drawing.circularSources).map(({ uid, primaryFile, pairedFile }) => (
      { uid, primaryFile, pairedFile, error: '' })) });
};
const rangedCircularDrawing = async (app, uids = ['card-1']) => {
  app.selectDrawing('circular');
  const drawing = state.activeDrawing();
  const { replaceCircularSources, createCircularSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  replaceCircularSources(drawing.slices.inputs, uids.map((uid) => createCircularSource({ uid })));
  for (const uid of uids) app.setCircularSourceFile(uid, 'gb', new File(['LOCUS'], `${uid}.gb`));
  // The stand-in files hold no record: their real read ends in an error. Wait for it before the cards are read as
  // CARD_RECORDS says, so it cannot land on them during a later step.
  for (let tick = 0; tick < 500 && drawing.circularRecordDiscovery.status !== 'error'; tick += 1) {
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
  assert.equal(drawing.circularRecordDiscovery.status, 'error');
  drawing.form.multi_record_canvas = false;
  drawing.form.circular_record_selector = '';
  await readCards(drawing);
  drawing.form.circular_record_ranges = [
    { source: 'card-1', record: 'A', windows: [{ start: 11, end: 90, reverseComplement: false, column: null, acrossOrigin: false }] },
    { source: 'card-1', record: 'B', windows: [{ start: 21, end: 60, reverseComplement: true, column: null, acrossOrigin: false }] }
  ];
  return drawing;
};
const cropFields = ({ form }) => [form.circular_region_start, form.circular_region_end, form.circular_reverse];

// PD-OI-044 scenario revision 5: a record chosen alone puts its crop in the one-record crop fields; a range changed
// there stays through All records (which draw the loaded crops) and draws when that record is chosen alone again;
// choosing another record puts that record's crop in the fields.
test('a grid or batch crop changed in the crop fields stays through All records and comes back with its record', async () => {
  const app = createAppSetup();
  const drawing = await rangedCircularDrawing(app);
  const loaded = JSON.parse(JSON.stringify(drawing.form.circular_record_ranges));
  assert.equal(app.setCircularRecordPresentationSelector('card-1:1').status, 'ok');
  assert.deepEqual(cropFields(drawing), [11, 90, false]);
  drawing.form.circular_region_end = 50;
  app.setCircularRecordPresentationSelector('');
  assert.deepEqual(cropFields(drawing), [11, 50, false], 'All records keeps the fields');
  assert.deepEqual(drawing.form.circular_record_ranges, loaded, 'and draws the loaded crops');
  app.setCircularRecordPresentationSelector('card-1:1');
  assert.deepEqual(cropFields(drawing), [11, 50, false], 'the same record alone again draws the changed crop');
  app.setCircularRecordPresentationSelector('card-1:2');
  assert.deepEqual(cropFields(drawing), [21, 60, true], 'another record puts its own crop in the fields');
  app.setCircularRecordPresentationSelector('');
  app.setCircularRecordPresentationSelector('card-1:1');
  assert.deepEqual(cropFields(drawing), [11, 90, false], 'and the first record then takes its loaded crop again');
  app.setCircularRecordPresentationSelector('');
});

// PD-OI-044 scenario revision 5: the crop fields filled from a record's crop go with the record's card when it is
// removed, or when its file of several records is replaced, so another record chosen alone draws no crop (OV-3140).
for (const change of ['remove', 'replace']) {
  test(`a ${change}d card takes the crop fields filled from its record's crop with it (OV-3140)`, async () => {
    const app = createAppSetup();
    const drawing = await rangedCircularDrawing(app, ['card-1', 'card-2']);
    app.setCircularRecordPresentationSelector('card-1:1');
    assert.deepEqual(cropFields(drawing), [11, 90, false]);
    app.setCircularRecordPresentationSelector('');
    app.setCircularSourceFile('card-1', 'gb', change === 'remove' ? null : new File(['LOCUS'], 'card-1-new.gb'));
    await readCards(drawing);
    assert.deepEqual(drawing.form.circular_record_ranges, []);
    assert.equal(app.setCircularRecordPresentationSelector('card-2').status, 'ok');
    assert.deepEqual([...cropFields(drawing), drawing.form.circular_windows, drawing.form.circular_range_fields_record],
      [null, null, false, [], ''], 'the other record alone draws whole');
    app.setCircularRecordPresentationSelector('');
  });
}

// PD-OI-044 scenario revision 5: a grid or batch drawing's crops are its settings. They stay through a record turned
// OFF and ON and in a drawing's copy, with its card uids; Reset Settings clears them, so each record draws whole, and
// Undo brings them back with the other settings.
test('grid or batch crops stay through record OFF and ON and a duplicate; Reset Settings clears them and Undo restores them', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const drawing = await rangedCircularDrawing(app);
  const { circularRequestDiagrams } = await import('../../gbdraw/web/js/services/record-options.js');
  const lengths = (target) => circularRequestDiagrams(target.circularRecordList.value, target).map((diagram) => diagram.lengths);
  const loaded = JSON.parse(JSON.stringify(drawing.form.circular_record_ranges));
  assert.deepEqual(lengths(drawing), [[80], [40]], 'each record draws its crop');
  assert.equal(await app.recordSelection.setRecordsDrawn('circular', 'card-1', ['card-1:1'], false), true);
  await readCards(drawing);
  assert.deepEqual(drawing.recordsOff, ['card-1:1']);
  assert.deepEqual(lengths(drawing), [[40]]);
  assert.equal(await app.recordSelection.setRecordsDrawn('circular', 'card-1', ['card-1:1'], true), true);
  await readCards(drawing);
  assert.deepEqual(lengths(drawing), [[80], [40]], 'a record turned ON again draws its crop');

  const outcome = app.drawingList.duplicateDrawing('circular');
  assert.equal(outcome.status, 'ok');
  const copy = state.drawingOf(outcome.drawingId);
  assert.deepEqual(copy.form.circular_record_ranges, loaded, 'the copy keeps the crops');
  assert.deepEqual(copy.circularSources.map(({ uid }) => uid), ['card-1'], 'of the same cards');
  assert.notEqual(copy.form.circular_record_ranges, drawing.form.circular_record_ranges);
  await history.undo();
  app.selectDrawing('circular');

  // The crop fields hold the first record's crop when Reset Settings runs.
  app.setCircularRecordPresentationSelector('card-1:1');
  app.setCircularRecordPresentationSelector('');
  assert.deepEqual(cropFields(drawing), [11, 90, false]);
  await app.resetSettings();
  assert.deepEqual(drawing.form.circular_record_ranges, [], 'Reset Settings clears the crops');
  assert.deepEqual([...cropFields(drawing), drawing.form.circular_windows, drawing.form.circular_range_fields_record],
    [null, null, false, [], ''], 'and the crop fields filled from them');
  drawing.form.multi_record_canvas = false;
  await readCards(drawing);
  assert.deepEqual(lengths(drawing), [[100], [100]], 'so each record draws whole');
  await history.undo();
  assert.deepEqual(drawing.form.circular_record_ranges, loaded, 'Undo restores the crops');
});

// PD-OI-044 scenario revision 5: Delete settings of an OFF record clears its grid or batch crop, and the one-record
// crop fields filled from it, so the record draws whole, alone or not; the other record keeps its crop.
test('Delete settings clears a record\'s grid or batch crop and the crop fields filled from it', async () => {
  const app = createAppSetup();
  const drawing = await rangedCircularDrawing(app);
  app.setCircularRecordPresentationSelector('card-1:1');
  drawing.form.circular_region_end = 50;
  app.setCircularRecordPresentationSelector('');
  assert.equal(await app.recordSelection.setRecordsDrawn('circular', 'card-1', ['card-1:1'], false), true);
  await readCards(drawing);
  assert.equal(await app.recordSelection.deleteSettings('circular', ['card-1:1']), true);
  assert.deepEqual(drawing.form.circular_record_ranges.map(({ record }) => record), ['B']);
  assert.deepEqual([...cropFields(drawing), drawing.form.circular_range_fields_record], [null, null, false, '']);
  assert.equal(await app.recordSelection.setRecordsDrawn('circular', 'card-1', ['card-1:1'], true), true);
  await readCards(drawing);
  app.setCircularRecordPresentationSelector('card-1:2');
  assert.deepEqual(cropFields(drawing), [21, 60, true], 'the other record keeps its crop');
  app.setCircularRecordPresentationSelector('');
});

// D-15-186: "Every drawing that uses this file" stays all-or-nothing. While the
// drawing switch waits it writes nothing, not even the clicked card's upload,
// and the busy reason is shown once through `onBusy`.
test('Every drawing while the drawing switch is busy writes nothing and shows the busy reason (D-15-186)', async () => {
  const { createInputReplacement } = await import('../../gbdraw/web/js/app/drawing-list.js');
  const busy = { status: 'busy', reason: 'Generating diagram. Retry after generation finishes.' };
  let waiting = busy;
  const calls = [];
  const shared = { fileId: 'f1', name: 'a.gb' };
  const drawings = { circular: createDrawingState('circular', 'circular'), 'circular-2': createDrawingState('circular-2', 'circular') };
  const { replaceCircularSources, createCircularSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  for (const drawing of Object.values(drawings)) replaceCircularSources(drawing.slices.inputs, [createCircularSource({ uid: 'card-1', gb: shared })]);
  const replacement = createInputReplacement({
    drawingTabs: () => Object.keys(drawings).map((id) => ({ id })),
    shownDrawingId: () => 'circular',
    switchBusy: () => waiting,
    drawingOf: (drawingId) => drawings[drawingId],
    runProjectStep: (label, apply) => { calls.push(['step', label]); apply(); },
    runDrawingStep: (label, apply) => { calls.push(['drawingStep', label]); apply(); },
    applyInputReplacementToDrawings: () => calls.push(['replace']),
    askScope: async () => 'all',
    onBusy: (reason) => calls.push(['busy', reason])
  });
  const outcome = await replacement.replaceFromShownDrawing(
    shared, { fileId: 'f2', name: 'b.gb' }, () => calls.push(['writeShown'])
  );
  assert.deepEqual(calls, [['busy', busy]], 'nothing is written, the reason is shown once');
  assert.equal(outcome.status, 'busy');
  waiting = null;
});

// #1777 round 13 (OV-3205): a copy draws its source's rings from the source's Result, so it holds the source's LOSAT
// results, one row per drawing as a Session Load gives them: Save Raw LOSAT TSV saves the same tables in the copy,
// Clear still counts each search once (the searches are shared, never copied), and Undo and Redo of the Duplicate
// take the copy's rows away and back with the copy.
test('Duplicate gives the copy its source\'s LOSAT results; its Undo and Redo take them away and back', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { webLosatRuntimeRecord } = await import('../../gbdraw/web/js/services/losat-cache.js');
  const source = modeDrawing(state, 'circular').id;
  app.selectDrawing(source);
  const keys = ['ring-a', 'ring-b', 'ring-c'];
  const files = keys.map((key) => `${key}.tsv`);
  state.losatCache.value = new Map(keys.map((key) => [key, {
    schema: 2, kind: 'raw-losat', identityKind: 'nucleotide', key, text: `${key}\n`, program: 'blastn', outfmt: '6',
    args: [], queryCanonicalHash: `q-${key}`, subjectCanonicalHash: `s-${key}`, runtime: webLosatRuntimeRecord('blastn')
  }]));
  state.losatCacheInfo.value = keys.map((key, index) => ({ key, filename: files[index], display: true, drawingId: source }));
  const createElement = document.createElement;
  const saved = async () => {
    const names = [];
    document.createElement = () => ({ click() { names.push(this.download); } });
    try { await app.downloadLosatCache(); } finally { document.createElement = createElement; }
    return names;
  };
  const rowsOf = (drawingId) => state.losatCacheInfo.value.filter((row) => row.drawingId === drawingId).map(({ key }) => key);
  try {
    assert.deepEqual(await saved(), files);
    const { drawingId: copyId } = app.drawingList.duplicateDrawing(source);
    assert.equal(state.activeDrawingId.value, copyId);
    assert.deepEqual(await saved(), files, 'the copy saves the tables its rings draw');
    assert.equal(app.clearableLosatResultCount.value, 3, 'Clear still counts 3 searches');
    assert.deepEqual([rowsOf(source), rowsOf(copyId)], [keys, keys], 'one row per drawing that holds a search');
    assert.equal(state.losatCache.value.size, 3, 'the searches are shared, never copied');

    await history.undo();
    assert.deepEqual([rowsOf(source), rowsOf(copyId)], [keys, []], 'Undo takes the copy\'s rows with the copy');
    await history.redo();
    assert.equal(state.activeDrawingId.value, copyId);
    assert.deepEqual(await saved(), files, 'Redo brings them back');
    assert.deepEqual([rowsOf(source), rowsOf(copyId)], [keys, keys]);
    await history.undo();
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
    app.selectDrawing(source);
  }
});

// #1777 round 14 (OV-3209, OV-3208): a drawing's LOSAT rows are in the state only while the drawing is. Delete holds
// them with its slot and its Undo writes them back, so New or Duplicate under the deleted drawing's ID neither takes
// nor erases them.
test('Delete holds a drawing\'s LOSAT rows; New and Duplicate under its ID neither take nor erase them', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const source = modeDrawing(state, 'circular').id;
  app.selectDrawing(source);
  const rowsOf = (drawingId) => state.losatCacheInfo.value.filter((row) => row.drawingId === drawingId).map(({ key }) => key);
  try {
    const { drawingId: second } = app.drawingList.newDrawing('circular');
    state.losatCacheInfo.value = [
      { key: 'k1', filename: 'k1.tsv', display: true, drawingId: source },
      { key: 'k2', filename: 'k2.tsv', display: true, drawingId: second }
    ];
    state.losatCache.value = new Map(['k1', 'k2'].map((key) => [key, { key }]));
    app.selectDrawing(source);
    app.drawingList.deleteDrawing(second);
    assert.deepEqual([rowsOf(source), rowsOf(second)], [['k1'], []], 'the deleted drawing\'s rows leave with it');
    const { drawingId: fresh } = app.drawingList.newDrawing('circular');
    assert.equal(fresh, second, 'New takes the free ID');
    assert.deepEqual(app.losatResultRows.value, [], 'and none of the deleted drawing\'s results (OV-3208)');
    await history.undo();
    const { drawingId: copy } = app.drawingList.duplicateDrawing(source);
    assert.equal(copy, second, 'Duplicate takes the free ID');
    assert.deepEqual([rowsOf(source), rowsOf(copy)], [['k1'], ['k1']], 'the copy holds its source\'s rows only');
    await history.undo();
    await history.undo();
    assert.ok(state.drawingTabs.value.some((tab) => tab.id === second), 'Undo of the Delete');
    assert.deepEqual([rowsOf(source), rowsOf(second)], [['k1'], ['k2']], 'brings the drawing back with its rows');
    await history.redo();
    assert.deepEqual([rowsOf(source), rowsOf(second)], [['k1'], []], 'and Redo takes them again');
    await history.undo();
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
    app.selectDrawing(source);
  }
});

// #1777 round 14 (OV-3210): Linear's per-pair readers read the shown drawing's rows, so a Linear drawing and its copy
// each list a window pair once, and a raw file name typed in the copy names only the copy's row.
test('a Linear drawing and its copy each list a window pair once, from their own rows', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const shownBefore = state.activeDrawingId.value;
  const catalog = state.featureCatalog.value;
  try {
    const { drawingId: linear } = app.drawingList.newDrawing('linear');
    state.featureCatalog.value = { items: [{ recordKeys: ['a#1', 'a#2', 'b'], windowOf: { 'a#1': 'a', 'a#2': 'a' } }] };
    state.losatCacheInfo.value = ['e1', 'e2'].map((edgeKey, index) => ({
      key: `w${index + 1}`, filename: `w${index + 1}.tsv`, display: true, drawingId: linear, edgeKey,
      queryUid: `a#${index + 1}`, subjectUid: 'b'
    }));
    state.losatCache.value = new Map(['w1', 'w2'].map((key) => [key, { key }]));
    const { drawingId: copy } = app.drawingList.duplicateDrawing(linear);
    const listed = () => Object.values(app.linearLosatWindowPairs.value).flat().map(({ edgeKey }) => edgeKey);
    const byEdge = () => Object.values(app.linearLosatCacheInfoByEdgeKey.value).map(({ edgeKey, drawingId }) => [edgeKey, drawingId]);
    assert.deepEqual(listed(), ['e1', 'e2'], 'the copy lists each window pair once');
    assert.deepEqual(byEdge(), [['e1', copy], ['e2', copy]], 'from its own rows');
    // A name typed on a window pair lives in the drawing's comparison plan (#1845), not in these rows.
    app.selectDrawing(linear);
    assert.deepEqual(listed(), ['e1', 'e2'], 'the source lists each pair once');
    assert.deepEqual(byEdge(), [['e1', linear], ['e2', linear]]);
    await history.undo();
    await history.undo();
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
    state.featureCatalog.value = catalog;
    app.selectDrawing(shownBefore);
  }
});
