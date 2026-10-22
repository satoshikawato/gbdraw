// OV-922: a History-tracked control records one step with two intent
// captures. Its reactive History state changes once per begin and commit, so
// the app renders for the step, not for each capture inside it (R11 steps,
// labels, and Undo availability unchanged).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const { ref, computed, watchEffect, nextTick } = loadPageVue();
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');

const settle = async () => { for (let i = 0; i < 5; i += 1) await nextTick(); };

const createApp = (buildIntentOf) => {
  const checked = ref(false);
  const history = createHistoryManager({
    buildIntent: buildIntentOf(() => ({ checked: checked.value })),
    applyIntent: (intent) => { checked.value = intent.checked; },
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {},
    makeRef: ref
  });
  // What the app template reads of History: Save availability, and the Undo
  // and Redo buttons with their titles.
  const saveAvailable = computed(() => !history.mutationPending());
  const canUndo = computed(() => { void history.revision.value; return history.canUndo(); });
  const canRedo = computed(() => { void history.revision.value; return history.canRedo(); });
  const undoTitle = computed(() => { void history.revision.value; return `Undo ${history.undoLabel()}`; });
  let renders = 0;
  watchEffect(() => {
    renders += 1;
    void [checked.value, saveAvailable.value, canUndo.value, canRedo.value, undoTitle.value];
  });
  return { checked, history, renders: () => renders };
};

const toggle = async (app) => {
  const tx = await app.history.begin('Normalize record lengths', { owner: 'checkbox' });
  app.checked.value = !app.checked.value;
  await nextTick();
  await app.history.commit(tx);
  await settle();
};

test('a tracked control renders for its value and its step, not per capture (OV-922)', async () => {
  const app = createApp((read) => () => read());
  await app.history.initializeIntentBaseline();
  await settle();
  const start = app.renders();
  const builds = app.history.getDiagnostics().intentBuilds;

  await toggle(app);
  // The value, then Undo becoming available.
  assert.equal(app.renders() - start, 2);
  assert.equal(app.history.getUndoCount(), 1);
  assert.equal(app.history.undoLabel(), 'Normalize record lengths');
  assert.equal(app.history.getDiagnostics().intentBuilds - builds, 2, 'one capture before, one after');

  const second = app.renders();
  await toggle(app);
  assert.equal(app.renders() - second, 1, 'only the value changes on screen');
  assert.equal(app.history.getUndoCount(), 2);

  // A button that changes nothing records no step and renders nothing.
  const idle = app.renders();
  const tx = await app.history.begin('Features', { owner: 'tab' });
  await app.history.commit(tx);
  await settle();
  assert.equal(app.renders(), idle);
  assert.equal(app.history.getUndoCount(), 2);
});

test('an asynchronous intent capture still reports the edit as applying (OV-922)', async () => {
  let release = () => {};
  const app = createApp((read) => () => new Promise((resolve) => { release = () => resolve(read()); }));
  const baseline = app.history.initializeIntentBaseline();
  release();
  await baseline;
  const begun = app.history.begin('Normalize record lengths', { owner: 'checkbox' });
  await nextTick();
  assert.equal(app.history.mutationPending(), true, 'Save waits while a capture is pending');
  release();
  const tx = await begun;
  assert.equal(app.history.mutationPending(), false);
  app.checked.value = true;
  const committed = app.history.commit(tx);
  await nextTick();
  release();
  await committed;
  assert.equal(app.history.getUndoCount(), 1);
});

// OV-946: a synchronous capture opens the step in the turn that asks for it.
// A control that takes focus while an action opens its step (a dialog that
// closes returns focus to the control that opened it) joins that step; it
// does not open a second step from the state before the action, which an
// Undo would then record as a new edit.
test('a control focused while an action opens its step joins that step (OV-946)', async () => {
  const app = createApp((read) => () => read());
  await app.history.initializeIntentBaseline();
  /** @type {Promise<unknown> | null} */
  let focused = null;
  // The render that closes the dialog runs before the action awaits anything.
  queueMicrotask(() => { focused = app.history.begin('Change setting', { source: 'input-adapter', owner: 'select' }); });
  await app.history.runUndoable('Change setting and reset Feature placements', () => { app.checked.value = true; });
  await focused;
  await settle();
  assert.equal(app.history.getUndoCount(), 1);
  assert.equal(app.history.mutationPending(), false, 'no step stays open');

  assert.equal(await app.history.undo(), true);
  assert.equal(app.checked.value, false);
  assert.equal(app.history.getUndoCount(), 0, 'Undo records no step');
  assert.equal(app.history.redoLabel(), 'Change setting and reset Feature placements');
  assert.equal(await app.history.redo(), true);
  assert.equal(app.checked.value, true);
  assert.equal(app.history.getUndoCount(), 1);
});

// OV-922: buttons that change only the view record no History step. They are
// History navigation (`data-history-ignore`), as the Editor toggle is, so a
// click captures no intent.
test('the editor tabs, Close editor, and Layout edit are History navigation (OV-922)', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const openTag = (marker) => {
    const at = html.indexOf(marker);
    assert.ok(at > 0, marker);
    return html.slice(html.lastIndexOf('<button', at), html.indexOf('>', at) + 1);
  };
  for (const marker of [
    `@click="openRightDrawerTab('legend')"`, `@click="openRightDrawerTab('features')"`,
    `@click="openRightDrawerTab('orthogroups')"`, `:aria-label="t('Close editor')"`,
    '@click="toggleRightDrawer"', '@click="layoutRepositionMode = !layoutRepositionMode"'
  ]) assert.match(openTag(marker), /\sdata-history-ignore[\s>]/, marker);
});
