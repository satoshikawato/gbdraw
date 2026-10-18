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

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');

/** @param {'circular' | 'linear'} id */
const hashOf = (id) => {
  const drawing = state.drawings[id];
  return createHash('sha256').update(JSON.stringify({
    form: drawing.form, adv: drawing.adv, overrides: drawing.unmanagedConfigOverrides
  })).digest('hex');
};

test('Undo walks the shown drawing\'s History and never changes the other drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { circular, linear } = state.drawings;
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

  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  assert.equal(history.getUndoCount(), 1, 'the Circular drawing keeps its own step');
  await app.undoHistory();
  assert.equal(circular.unmanagedConfigOverrides['objects.circular.a'], 1, 'Undo restores the Circular edit');
  assert.equal(hashOf('linear'), linearAfterEdit, 'Undo in Circular leaves Linear unchanged');
});

test('Reset all drawings is one project step: one Undo restores every drawing', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { circular, linear } = state.drawings;
  if (state.activeDrawingId.value !== 'linear') app.selectDrawing('linear');
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
  const { circular } = state.drawings;
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  const circularSteps = history.getUndoCount();
  const before = circular.form.prefix;
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
  // The step ended at the switch, before the title edit: the title is undone first.
  assert.equal(history.undoLabel(), 'Rename Session');
  await app.undoHistory();
  assert.equal(circular.form.prefix, 'typed-before-switch');
  await app.undoHistory();
  assert.equal(circular.form.prefix, before);
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
  const { circular } = state.drawings;
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
    const { circular } = state.drawings;
    circular.form.multi_record_canvas = canvas;
    circular.layoutPreferences.circular.multi.legend = 'top';
    await app.resetSettings();
    assert.deepEqual(JSON.parse(JSON.stringify(circular.layoutPreferences)), fresh,
      `the Circular layout after a Reset with ${shown} shown and the canvas ${canvas ? 'on' : 'off'}`);
    assert.deepEqual(JSON.parse(JSON.stringify(state.drawings.linear.layoutPreferences)), fresh);
  }
});
