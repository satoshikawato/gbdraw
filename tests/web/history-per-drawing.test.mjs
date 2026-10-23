// PD-OI-090 (OIC-030), A3: each drawing has its own History. A drawing switch
// writes no setting and is no Undo step; Undo and Redo walk the shown
// drawing's steps and the project steps, newest first, and a step in one
// drawing never changes the other (checked by a hash of each drawing per
// step). Reset all drawings is one project step (D-15-51, D-15-52), and an
// edit of a Session-level value is a project step (R15-23).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
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

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { inputTypeKeyOf } = await import('../../gbdraw/web/js/services/drawing-view.js');
const { DRAWING_SLICES } = await import('./helpers/drawing-state.mjs');

// A3 review Low 5: every slice of the drawing and its input type.
/** @param {'circular' | 'linear'} id */
const hashOf = (id) => {
  const drawing = modeDrawing(state, id);
  const slices = Object.fromEntries(Object.entries(DRAWING_SLICES)
    .map(([name, slice]) => [name, slice.capture(drawing.slices[name])]));
  return createHash('sha256').update(JSON.stringify({
    slices, inputType: state[inputTypeKeyOf(drawing)].value
  })).digest('hex');
};

test('Undo walks the shown drawing\'s History and never changes the other drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  circular.unmanagedConfigOverrides['objects.circular.a'] = 1;
  linear.unmanagedConfigOverrides['objects.linear.b'] = 2;

  await app.resetUnmanagedConfigOverride('objects.circular.a');
  assert.equal(app.canUndoHistory.value, true);
  const circularAfterEdit = hashOf('circular');

  // A switch writes no setting and records no step.
  const before = { circular: hashOf('circular'), linear: hashOf('linear') };
  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  assert.deepEqual({ circular: hashOf('circular'), linear: hashOf('linear') }, before);
  assert.equal(app.canUndoHistory.value, false, 'the Linear drawing has no step yet');

  await app.resetUnmanagedConfigOverride('objects.linear.b');
  const linearAfterEdit = hashOf('linear');
  assert.equal(history.getUndoCount(), 1);
  await app.undoHistory();
  assert.equal(linear.unmanagedConfigOverrides['objects.linear.b'], 2, 'Undo restores the Linear edit');
  assert.equal(hashOf('circular'), circularAfterEdit, 'Undo in Linear leaves Circular unchanged');
  assert.equal(app.canUndoHistory.value, false);
  await app.redoHistory();
  assert.equal(hashOf('linear'), linearAfterEdit);
  assert.equal(hashOf('circular'), circularAfterEdit, 'Redo in Linear leaves Circular unchanged');

  // A checkpoint step and a Generate step of Linear walk the same way.
  await history.runUndoableCheckpoint('Linear checkpoint', () => { state.lInputType.value = 'gff'; linear.form.prefix = 'linear-checkpoint'; return true; });
  const linearAfterCheckpoint = hashOf('linear');
  assert.notEqual(linearAfterCheckpoint, linearAfterEdit);
  await history.runUndoableArtifactReplacement('Generate diagram', () => {
    state.results.value = [{ name: 'linear.svg', content: '<svg id="linear-generate"/>' }];
    return true;
  });
  const steps = [linearAfterCheckpoint, linearAfterEdit];
  for (const expected of steps) {
    await app.undoHistory();
    assert.equal(hashOf('linear'), expected);
    assert.equal(hashOf('circular'), circularAfterEdit, 'a Linear Undo leaves Circular unchanged');
  }
  for (const expected of [linearAfterCheckpoint, linearAfterCheckpoint]) {
    await app.redoHistory();
    assert.equal(hashOf('linear'), expected);
    assert.equal(hashOf('circular'), circularAfterEdit, 'a Linear Redo leaves Circular unchanged');
  }

  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  assert.equal(history.getUndoCount(), 1, 'the Circular drawing keeps its own step');
  await app.undoHistory();
  assert.equal(circular.unmanagedConfigOverrides['objects.circular.a'], 1, 'Undo restores the Circular edit');
  assert.equal(hashOf('linear'), linearAfterCheckpoint, 'Undo in Circular leaves Linear unchanged');
  assert.equal(state.lInputType.value, 'gff', 'and its input type');
});

test('Reset all drawings is one project step: one Undo restores every drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  if (state.activeDrawingId.value !== 'linear') app.selectDrawing('linear');
  // The app's palette colors hold every palette key, as a palette load leaves them.
  const { PALETTE_SLICE } = await import('../../gbdraw/web/js/services/drawing/palette.js');
  const { normalizePaletteColors } = await import('../../gbdraw/web/js/utils/color-utils.js');
  [circular, linear].forEach((drawing) => PALETTE_SLICE.install(drawing.slices.palette, {
    currentColors: normalizePaletteColors(PALETTE_SLICE.capture(drawing.slices.palette).currentColors)
  }));
  circular.form.prefix = 'circular-prefix';
  linear.form.prefix = 'linear-prefix';
  const circularBefore = hashOf('circular');
  const linearBefore = hashOf('linear');
  await app.resetSettings();
  assert.notEqual(circular.form.prefix, 'circular-prefix');
  assert.notEqual(linear.form.prefix, 'linear-prefix');
  const circularReset = hashOf('circular');
  const linearReset = hashOf('linear');

  // A newer step of the shown drawing is undone first.
  await history.runUndoable('Linear prefix', () => { linear.form.prefix = 'after-reset'; });
  assert.equal(history.undoLabel(), 'Linear prefix');
  await app.undoHistory();
  assert.equal(hashOf('linear'), linearReset);
  assert.equal(history.undoLabel(), 'Reset settings');

  await app.undoHistory();
  assert.equal(hashOf('circular'), circularBefore, 'one Undo in Linear restores Circular');
  assert.equal(hashOf('linear'), linearBefore, 'and Linear');
  await app.redoHistory();
  assert.equal(history.redoLabel(), 'Linear prefix', 'Redo replays in the order the steps ended');
  assert.equal(hashOf('circular'), circularReset, 'Redo resets Circular again');
  assert.equal(hashOf('linear'), linearReset, 'and Linear');

  // From the other drawing the Reset is the newest step too.
  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  assert.equal(history.undoLabel(), 'Reset settings');
  await app.undoHistory();
  assert.equal(hashOf('circular'), circularBefore, 'one Undo in Circular restores Circular');
  assert.equal(hashOf('linear'), linearBefore, 'and Linear');
});

test('a switch settles the departing drawing\'s open step in that drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  const circularSteps = history.getUndoCount();
  const before = circular.form.prefix;
  const title = state.sessionTitle.value;
  // A text control began on focus; the switch comes before its own commit.
  const step = await history.begin('Output prefix', { source: 'test-control' });
  assert.ok(step, 'the control opened a step');
  circular.form.prefix = 'typed-before-switch';
  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  const linearSteps = history.getUndoCount();
  // The Session title is shared; it changes in a step while Linear is shown.
  await history.runUndoable('Rename Session', () => { state.sessionTitle.value = 'title set in Linear'; });
  const titleSteps = history.getUndoCount();
  await history.commit(step);
  assert.equal(titleSteps, linearSteps + 1, 'the title edit is a project step');
  assert.equal(history.getUndoCount(), titleSteps, 'the Linear walk gains nothing from the Circular step');
  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  assert.equal(history.getUndoCount(), circularSteps + 2, 'the step is the Circular drawing\'s');
  // The step ended at the switch, before the title edit.
  assert.equal(history.undoLabel(), 'Rename Session');
  await app.undoHistory();
  assert.equal(state.sessionTitle.value, title, 'the title is undone first');
  assert.equal(circular.form.prefix, 'typed-before-switch');
  await app.undoHistory();
  assert.equal(circular.form.prefix, before);
  assert.equal(state.sessionTitle.value, title, 'the Circular step leaves the title');
});

// R15-23 (Q-DR-8 A): an edit of a Session-level value is a project step,
// undone newest first with the shown drawing's steps from any drawing.
test('the Session title edits are undone newest first from either drawing', async () => {
  for (const undoIn of /** @type {const} */ (['linear', 'circular'])) {
    const app = createAppSetup();
    const history = window.__GBDRAW_HISTORY__;
    if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
    const original = state.sessionTitle.value;
    await history.runUndoable('Rename Session', () => { state.sessionTitle.value = `A ${undoIn}`; });
    assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
    await history.runUndoable('Rename Session', () => { state.sessionTitle.value = `B ${undoIn}`; });
    if (undoIn !== 'linear') assert.deepEqual(app.selectDrawing(undoIn), { status: 'ok' });
    await app.undoHistory();
    assert.equal(state.sessionTitle.value, `A ${undoIn}`, `Undo in ${undoIn} undoes B first`);
    await app.undoHistory();
    assert.equal(state.sessionTitle.value, original, `then A (${undoIn})`);
    await app.redoHistory();
    assert.equal(state.sessionTitle.value, `A ${undoIn}`);
  }
});

test('a title edit in Circular is undone in Linear', async () => {
  const app = createAppSetup();
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  const original = state.sessionTitle.value;
  await window.__GBDRAW_HISTORY__.runUndoable('Rename Session', () => { state.sessionTitle.value = 'renamed in Circular'; });
  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  assert.equal(app.canUndoHistory.value, true);
  await app.undoHistory();
  assert.equal(state.sessionTitle.value, original);
});

// A3 review M1, R15-23: a checkpoint Undo restores what its drawing owns. The
// Session title and Auto Reflow edits made later in Linear are project steps,
// undone before the older Circular step, and the Linear input type stays.
test('a Circular step older than the Session edits is undone after them and keeps Linear\'s input type', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  state.lInputType.value = 'gb';
  const before = circular.form.prefix;
  const title = state.sessionTitle.value;
  const reflowBefore = state.autoLabelReflowEnabled.value;
  await history.runUndoableCheckpoint('Circular checkpoint', () => { circular.form.prefix = 'checkpoint-prefix'; return true; });
  assert.equal(circular.form.prefix, 'checkpoint-prefix');

  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  await history.runUndoable('Linear input type', () => { state.lInputType.value = 'gff'; });
  await history.runUndoable('Rename Session', () => { state.sessionTitle.value = 'title set in Linear'; });
  await history.runUndoable('Auto Reflow', () => { state.autoLabelReflowEnabled.value = !reflowBefore; });

  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  await app.undoHistory();
  assert.equal(state.autoLabelReflowEnabled.value, reflowBefore, 'the newest step, Auto Reflow, is undone first');
  assert.equal(circular.form.prefix, 'checkpoint-prefix');
  await app.undoHistory();
  assert.equal(state.sessionTitle.value, title, 'then the title');
  assert.equal(circular.form.prefix, 'checkpoint-prefix');
  await app.undoHistory();
  assert.equal(circular.form.prefix, before, 'then the Circular step');
  assert.equal(state.lInputType.value, 'gff', 'the Linear input type stays');
  assert.equal(state.sessionTitle.value, title);
});

// OV-1166 (Q-DR-6 A): Reset Settings gives every drawing's multi-record
// layout the defaults a new drawing has, whatever the canvas was before and
// whichever drawing is shown.
test('Reset Settings leaves the multi-record layout of a new drawing', async () => {
  const { createDefaultLayoutPreferences } = await import('../../gbdraw/web/js/services/layout-preferences.js');
  const fresh = createDefaultLayoutPreferences();
  for (const [shown, canvas] of /** @type {const} */ ([['circular', false], ['circular', true], ['linear', false]])) {
    const app = createAppSetup();
    if (state.activeDrawingId.value !== shown) app.selectDrawing(shown);
    const circular = modeDrawing(state, 'circular');
    circular.form.multi_record_canvas = canvas;
    circular.layoutPreferences.circular.multi.legend = 'top';
    await app.resetSettings();
    assert.deepEqual(JSON.parse(JSON.stringify(circular.layoutPreferences)), fresh,
      `the Circular layout after a Reset with ${shown} shown and the canvas ${canvas ? 'on' : 'off'}`);
    assert.deepEqual(JSON.parse(JSON.stringify(modeDrawing(state, 'linear').layoutPreferences)), fresh);
  }
});

// OV-1167 (Q-DR-7 B): Reset Settings re-derives the Circular comparison rings
// from the kept ring files, as a new drawing with those files uploaded has them.
test('Reset Settings re-derives the comparison rings from the kept ring files', async () => {
  const app = createAppSetup();
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  const rings = modeDrawing(state, 'circular').circularConservation;
  const ringsOf = () => JSON.parse(JSON.stringify({ enabled: rings.enabled, series: rings.series }));
  state.files.c_conservation_fastas = [new File(['>r\nACGT\n'], 'ring.fa')];
  app.changeCircularConservationSource();
  const uploaded = ringsOf();
  assert.equal(uploaded.enabled, true);
  assert.equal(uploaded.series.length, 1);
  rings.series[0].label = 'renamed ring';
  await app.resetSettings();
  assert.equal(state.files.c_conservation_fastas.length, 1, 'the ring file stays');
  assert.deepEqual(ringsOf(), uploaded);
  state.files.c_conservation_fastas = [];
  app.changeCircularConservationSource();
});

// OV-1168: the applied palette is the shown Result's (an artifact slot
// value). An Undo of Reset made in Circular, from Linear, restores Circular's
// part without writing Circular's applied palette onto the shown Linear.
test('an Undo of Reset from the other drawing keeps the shown drawing\'s applied palette', async () => {
  const app = createAppSetup();
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  state.appliedPaletteName.value = 'circular-palette';
  await app.resetSettings();
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  const linearApplied = state.appliedPaletteName.value;
  await app.undoHistory();
  assert.equal(state.appliedPaletteName.value, linearApplied);
});

/** @param {ReturnType<typeof createAppSetup>} app @param {'circular' | 'linear'} id */
const show = (app, id) => {
  if (state.activeDrawingId.value !== id) assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
};

// D-15-53 1 (PD-OI-090, OIC-030): a step that changes a drawing value and a
// Session-level value (Add rule clears the rule draft) is that drawing's step
// whole. The other drawing's History keeps its Redo and never takes part of
// it; its Undo restores the draft only where still as the step left it.
test('a step that also clears a draft is its drawing\'s step and leaves the other drawing\'s History', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  show(app, 'circular');
  const draftBefore = state.newSpecRule.val;
  await history.runUndoable('Type pattern', () => { state.newSpecRule.val = 'abc'; });
  show(app, 'linear');
  await history.runUndoable('Linear prefix', () => { linear.form.prefix = 'lin-1'; });
  await app.undoHistory();
  assert.equal(history.canRedo(), true);
  show(app, 'circular');
  const prefixBefore = circular.form.prefix;
  await history.runUndoable('Add specific color rule', () => { circular.form.prefix = 'rule-added'; state.newSpecRule.val = ''; });

  show(app, 'linear');
  assert.equal(history.canRedo(), true, 'a Circular step keeps Linear\'s Redo');
  assert.equal(history.undoLabel(), 'Type pattern', 'Linear does not see the Circular step');
  await app.undoHistory();
  assert.equal(state.newSpecRule.val, draftBefore, 'the draft stays as the Circular step left it');
  assert.equal(circular.form.prefix, 'rule-added');

  show(app, 'circular');
  assert.equal(history.undoLabel(), 'Add specific color rule');
  await app.undoHistory();
  assert.deepEqual({ prefix: circular.form.prefix, draft: state.newSpecRule.val }, { prefix: prefixBefore, draft: 'abc' },
    'one Undo restores the rule and the draft');
  assert.equal(history.redoLabel(), 'Type pattern', 'Redo replays in the order the steps ended');
  await app.redoHistory();
  await app.redoHistory();
  assert.deepEqual({ prefix: circular.form.prefix, draft: state.newSpecRule.val }, { prefix: 'rule-added', draft: '' });
});

// D-15-53 2: a Session-level value changed outside any step (an annotation
// click) marks nothing unsaved and rides with the next step, in that step's
// own drawing.
test('a selection outside any step is unsaved only with the next step, in its drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'circular');
  const selectedBefore = state.selectedAnnotation.value;
  const prefixBefore = circular.form.prefix;
  history.markSavePoint();
  state.selectedAnnotation.value = { setId: 'set-a', id: 'annotation-a' };
  assert.equal(await history.hasChangesSinceSavePoint(), false, 'a click alone is not unsaved work');
  await history.runUndoable('Circular prefix', () => { circular.form.prefix = 'after-click'; });
  assert.equal(await history.hasChangesSinceSavePoint(), true);
  show(app, 'linear');
  assert.equal(history.canUndo(), false, 'Linear does not see the Circular step');
  show(app, 'circular');
  assert.equal(history.undoLabel(), 'Circular prefix');
  await app.undoHistory();
  assert.deepEqual({ prefix: circular.form.prefix, selected: state.selectedAnnotation.value },
    { prefix: prefixBefore, selected: selectedBefore });
  // A title typed in an open step is unsaved work before the step ends.
  history.markSavePoint();
  const step = await history.begin('Rename Session', { source: 'test-control' });
  state.sessionTitle.value = 'typed title';
  assert.equal(await history.hasChangesSinceSavePoint(), true);
  await history.commit(step);
});

// D-15-53 2: the switch resets the feature panel tab (no step); the router
// takes it as the recorded Session value.
test('the feature panel tab reset at a switch is no change and no step', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  show(app, 'circular');
  state.featurePanelTab.value = 'labels';
  await history.initializeIntentBaseline('Load');
  history.markSavePoint();
  show(app, 'linear');
  assert.equal(state.featurePanelTab.value, 'colors');
  assert.equal(await history.hasChangesSinceSavePoint(), false);
  await history.runUndoable('Linear prefix', () => { modeDrawing(state, 'linear').form.prefix = 'lin-after-switch'; });
  show(app, 'circular');
  assert.equal(history.canUndo(), false, 'Circular does not see the Linear step');
  show(app, 'linear');
  await app.undoHistory();
  assert.equal(state.featurePanelTab.value, 'colors', 'the Linear Undo leaves the tab the switch set');
});

const settleTurns = async () => { for (let turn = 0; turn < 20; turn += 1) await new Promise((r) => setTimeout(r, 0)); };

// A3 review 3 Low 2: the switch takes only the tab reset as recorded; a title
// typed in a step still open at the switch is that step's change.
test('a title typed in a step open at the switch is undone, without the tab reset', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  show(app, 'circular');
  state.featurePanelTab.value = 'labels';
  await history.captureBaseline();
  const title = state.sessionTitle.value;
  const step = await history.begin('Rename Session', { source: 'test-control' });
  state.sessionTitle.value = 'typed before switch';
  show(app, 'linear');
  await settleTurns();
  await history.commit(step);
  assert.equal(state.featurePanelTab.value, 'colors');
  assert.equal(history.undoLabel(), 'Rename Session', 'the title edit is a project step');
  show(app, 'circular');
  assert.equal(history.undoLabel(), 'Rename Session');
  await app.undoHistory();
  assert.deepEqual({ title: state.sessionTitle.value, tab: state.featurePanelTab.value }, { title, tab: 'colors' });
});

// A3 review 3 Low 2: one open step changes a drawing value and a draft; the
// switch ends it in its drawing with both.
test('a prefix and a draft changed in a step open at the switch are undone together', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'circular');
  await history.captureBaseline();
  const before = { prefix: circular.form.prefix, draft: state.newSpecRule.val };
  await history.begin('Edit', { source: 'test-control' });
  circular.form.prefix = 'open-at-switch';
  state.newSpecRule.val = 'draft-at-switch';
  show(app, 'linear');
  await settleTurns();
  assert.equal(history.canUndo(), false, 'Linear does not see the Circular step');
  show(app, 'circular');
  assert.deepEqual({ label: history.undoLabel(), count: history.getUndoCount() }, { label: 'Edit', count: 1 });
  await app.undoHistory();
  assert.deepEqual({ prefix: circular.form.prefix, draft: state.newSpecRule.val }, before);
});

// A3 review 3 Medium 1 (D-15-53 2): a step that changes nothing records
// nothing, even after a selection made outside any step. The selection is no
// unsaved work and stays for the next step, which keeps it.
test('a step that changes nothing after a click records nothing and keeps the other drawing\'s Redo', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'linear');
  await history.captureBaseline();
  await history.runUndoable('Linear prefix', () => { modeDrawing(state, 'linear').form.prefix = 'redo-kept'; });
  await app.undoHistory();
  show(app, 'circular');
  history.markSavePoint();
  const selectedBefore = state.selectedAnnotation.value;
  const prefixBefore = circular.form.prefix;
  state.selectedAnnotation.value = { setId: 'set-b', id: 'annotation-b' };
  const step = await history.begin('Output prefix', { source: 'test-control' });
  assert.equal(await history.hasChangesSinceSavePoint(), false, 'a focus after a click is not unsaved work');
  await history.commit(step);
  assert.equal(await history.hasChangesSinceSavePoint(), false, 'a blur without an edit is not unsaved work');
  assert.equal(history.canUndo(), false, 'the blur records no step');
  show(app, 'linear');
  assert.deepEqual({ canRedo: history.canRedo(), canUndo: history.canUndo() }, { canRedo: true, canUndo: false });
  show(app, 'circular');
  await history.runUndoable('Circular prefix', () => { circular.form.prefix = 'after-blur'; });
  await app.undoHistory();
  assert.deepEqual({ prefix: circular.form.prefix, selected: state.selectedAnnotation.value },
    { prefix: prefixBefore, selected: selectedBefore }, 'the next step keeps the click');
});

// A3 review 2 M1 (OV-1168 class): the alignment Reset receipt is an artifact
// slot value, the shown Result's. An Undo of Reset made in Linear, from
// Circular, neither validates nor installs a receipt for the hidden Linear.
test('an Undo of Reset from Circular does not read Linear\'s alignment receipt', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  show(app, 'linear');
  state.similarityAlignmentResetReceipt.value = { binding: 'fake', directions: [], referenceDeltaX: null };
  await app.resetSettings();
  show(app, 'circular');
  assert.equal(history.undoLabel(), 'Reset settings');
  assert.equal(await app.undoHistory(), true);
  assert.notEqual(history.undoLabel(), 'Reset settings');
  assert.equal(state.similarityAlignmentResetReceipt.value, null, 'the shown Circular slot takes no receipt');
});

// A3 review 2 Low 4: a drawing step drops no project step; only the oldest
// title step is gone after the limit.
test('a step past the project limit drops only the oldest project step', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  show(app, 'circular');
  await history.captureBaseline();
  for (let index = 0; index <= 30; index += 1) {
    await history.runUndoable('Rename Session', () => { state.sessionTitle.value = `T${index}`; });
  }
  await history.runUndoable('Add rule', () => { modeDrawing(state, 'circular').form.prefix = 'mixed'; state.newSpecRule.val = 'x'; });
  while (history.canUndo()) assert.equal(await app.undoHistory(), true);
  assert.equal(state.sessionTitle.value, 'T0', 'only the first title step is gone');
});

// OV-1169: a drawing that lost a step newer than a project step to the limit
// restores its part of the project step only where still as the step left it.
test('an Undo of Reset after the limit dropped a newer step keeps that step\'s change', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'circular');
  circular.form.prefix = 'pre-reset';
  await history.captureBaseline();
  await app.resetSettings();
  await history.runUndoable('Circular prefix', () => { circular.form.prefix = 'later'; });
  for (let index = 0; index < 30; index += 1) {
    await history.runUndoable('Species', () => { circular.form.species = `S${index}`; });
  }
  while (history.undoLabel() === 'Species') await app.undoHistory();
  assert.equal(history.undoLabel(), 'Reset settings', 'the limit dropped the prefix step');
  await app.undoHistory();
  assert.equal(circular.form.prefix, 'later', 'the dropped step\'s change stays');
  await app.redoHistory();
  assert.equal(circular.form.prefix, 'later', 'and stays on Redo');
});

// A3 review 2 P4: a second Undo pressed while the first restores a project
// step returns false and leaves both stacks consistent.
test('two Undo presses at once over project steps take one step', async () => {
  createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  await history.captureBaseline();
  for (const title of ['T1', 'T2', 'T3']) await history.runUndoable('Rename Session', () => { state.sessionTitle.value = title; });
  assert.deepEqual(await Promise.all([history.undo(), history.undo()]), [true, false]);
  assert.deepEqual({ title: state.sessionTitle.value, undo: history.getUndoCount(), redo: history.getRedoCount() },
    { title: 'T2', undo: 2, redo: 1 });
  await history.redo();
  assert.equal(state.sessionTitle.value, 'T3');
});

// A3 review 2 P6: a hidden drawing's part of Reset restores only the values
// still as Reset left them; a later step of that drawing is kept.
test('an Undo of Reset from the other drawing keeps a later step of the hidden drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'linear');
  circular.form.prefix = 'pre-reset';
  await history.captureBaseline();
  await app.resetSettings();
  show(app, 'circular');
  await history.runUndoable('Circular prefix', () => { circular.form.prefix = 'later'; });
  show(app, 'linear');
  assert.equal(history.undoLabel(), 'Reset settings');
  await app.undoHistory();
  assert.equal(circular.form.prefix, 'later');
});

// OV-1195: a hidden drawing's artifact slot values (here Linear's alignment
// Reset receipt) are captured and restored on its own slot.
test('an Undo of Reset from Circular gives Linear its alignment receipt back', async () => {
  const app = createAppSetup();
  show(app, 'linear');
  const receipt = { binding: 'fake', directions: [], referenceDeltaX: null };
  state.similarityAlignmentResetReceipt.value = receipt;
  await app.resetSettings();
  assert.equal(state.similarityAlignmentResetReceipt.value, null);
  show(app, 'circular');
  assert.equal(await app.undoHistory(), true);
  assert.equal(state.similarityAlignmentResetReceipt.value, null, 'the shown Circular slot takes no receipt');
  show(app, 'linear');
  assert.deepEqual(state.similarityAlignmentResetReceipt.value, receipt);
  state.similarityAlignmentResetReceipt.value = null;
});

// OV-1168: Reset all drawings resets every drawing's applied palette, the
// hidden one's in its slot, and its Undo gives each drawing its own back.
test('Reset all drawings resets each drawing\'s applied palette and its Undo restores each', async () => {
  const app = createAppSetup();
  show(app, 'linear');
  state.appliedPaletteName.value = 'linear-palette';
  show(app, 'circular');
  state.appliedPaletteName.value = 'circular-palette';
  await app.resetSettings();
  assert.equal(state.appliedPaletteName.value, 'default');
  show(app, 'linear');
  assert.equal(state.appliedPaletteName.value, 'default', 'the hidden drawing\'s applied palette is reset too');
  await app.undoHistory();
  assert.equal(state.appliedPaletteName.value, 'linear-palette');
  show(app, 'circular');
  assert.equal(state.appliedPaletteName.value, 'circular-palette');
});

// OV-1160 (Owner R15-23 Q-DR-5 A): turning Instant Preview on applies every
// drawing's queued palette in one project step; its Undo queues each again.
test('turning Instant Preview on applies every drawing\'s queued palette in one step', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  show(app, 'circular');
  app.setPaletteInstantPreview(false);
  state.appliedPaletteName.value = 'circular-applied';
  await history.captureBaseline();
  await history.runUndoable('Change setting', () => app.selectPalette('queued', 'palette'));
  assert.equal(circular.pendingPaletteName.value, 'queued');
  show(app, 'linear');
  const linearApplied = state.appliedPaletteName.value;
  await app.setPaletteInstantPreview(true);
  assert.equal(state.paletteInstantPreviewEnabled.value, true);
  assert.equal(circular.pendingPaletteName.value, '', 'the hidden drawing\'s queue is applied');
  assert.equal(state.appliedPaletteName.value, linearApplied, 'the shown Result keeps its palette');
  assert.equal(history.undoLabel(), 'Change setting');
  await app.undoHistory();
  assert.equal(state.paletteInstantPreviewEnabled.value, false);
  assert.equal(circular.pendingPaletteName.value, 'queued', 'the Undo queues it again');
  show(app, 'circular');
  assert.equal(state.appliedPaletteName.value, 'circular-applied');
  await app.redoHistory();
  assert.equal(state.paletteInstantPreviewEnabled.value, true);
  assert.equal(circular.pendingPaletteName.value, '');
  assert.equal(state.appliedPaletteName.value, 'queued', 'Redo from the other drawing applies it again');
  assert.equal(await app.undoHistory(), true);
  assert.equal(state.appliedPaletteName.value, 'circular-applied');
  assert.equal(circular.pendingPaletteName.value, 'queued');
});

// OV-1655: a file card's Remove from Linear takes a Legend row from hidden
// Circular in one project step. When Circular is shown, its figure is drawn
// again and lays the kept rows out anew (their xPos and yPos move); that is no
// later edit of the list, so the Undo from Linear brings the row back before
// Circular is shown again. A later edit of a row (a new caption) is kept, and
// the Undo leaves the list (OIC-030).
for (const later of ['redraw', 'rename']) test(`an Undo from the other drawing ${later === 'redraw' ? 'restores a Legend list a redraw only laid out anew' : 'keeps a Legend list a later step changed'} (OV-1655)`, async () => {
  const { replaceLegendRows } = await import('../../gbdraw/web/js/services/drawing/paint.js');
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const circular = modeDrawing(state, 'circular');
  const row = (caption, yPos) => ({ caption, originalCaption: caption, color: '#000000', xPos: 10, yPos, featureIds: [] });
  show(app, 'linear');
  replaceLegendRows(circular.slices.paint, [row('Hand', 0), row('Imported', 40), row('GC content', 80)]);
  await history.captureBaseline();
  await history.runProjectStep('Remove specific rules file', () => {
    replaceLegendRows(circular.slices.paint, circular.legendEntries.value.filter((entry) => entry.caption !== 'Imported'));
  });
  const captions = () => circular.legendEntries.value.map((entry) => entry.caption);
  assert.deepEqual(captions(), ['Hand', 'GC content']);
  show(app, 'circular');
  if (later === 'redraw') replaceLegendRows(circular.slices.paint, [row('Hand', 0), row('GC content', 40)]);
  else await history.runUndoable('Rename Legend row', () => replaceLegendRows(circular.slices.paint, [row('Hand', 0), row('GC', 80)]));
  show(app, 'linear');
  assert.equal(history.undoLabel(), 'Remove specific rules file');
  await app.undoHistory();
  assert.deepEqual(captions(), later === 'redraw' ? ['Hand', 'Imported', 'GC content'] : ['Hand', 'GC']);
});
