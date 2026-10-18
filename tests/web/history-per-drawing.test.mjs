// PD-OI-090 (OIC-030), A3: each drawing has its own History. A drawing switch
// writes no setting and is no Undo step; Undo and Redo walk the shown
// drawing's steps, and a step in one drawing never changes the other (checked
// by a hash of each drawing per step); Reset Settings is one step in each
// drawing's History (D-15-50).
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

test('Reset Settings is one step in each drawing\'s History', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const { circular, linear } = state.drawings;
  if (state.activeDrawingId.value !== 'circular') app.selectDrawing('circular');
  circular.form.prefix = 'circular-prefix';
  linear.form.prefix = 'linear-prefix';
  const circularBefore = hashOf('circular');
  const linearBefore = hashOf('linear');
  await app.resetSettings();
  assert.notEqual(circular.form.prefix, 'circular-prefix');
  assert.notEqual(linear.form.prefix, 'linear-prefix');
  const linearReset = hashOf('linear');

  await app.undoHistory();
  assert.equal(hashOf('circular'), circularBefore, 'Undo in Circular restores its reset');
  assert.equal(hashOf('linear'), linearReset, 'and leaves the Linear reset');
  assert.deepEqual(app.selectDrawing('linear'), { status: 'ok' });
  assert.equal(history.undoLabel(), 'Reset settings');
  await app.undoHistory();
  assert.equal(hashOf('linear'), linearBefore, 'Undo in Linear restores its reset');
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
  // The Session title is shared; it changes while Linear is shown.
  state.sessionTitle.value = 'title set in Linear';
  await history.commit(step);
  assert.equal(history.getUndoCount(), linearSteps, 'the Linear History records nothing');
  assert.deepEqual(app.selectDrawing('circular'), { status: 'ok' });
  assert.equal(history.getUndoCount(), circularSteps + 1, 'the step is the Circular drawing\'s');
  await app.undoHistory();
  assert.equal(circular.form.prefix, before);
  assert.equal(state.sessionTitle.value, 'title set in Linear',
    'the step ended at the switch, so its Undo keeps what changed after it');
});
