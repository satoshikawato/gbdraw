// The drawn-track visibility watcher (app/svg-styles.js): a Show Depth (or a
// GC track) toggle on the shown drawing hides or shows the drawn track on the
// displayed Result. A Load commits new drawing states under the same IDs
// (OV-1778), so the watcher follows the shown drawing in the store, not the
// drawing objects that existed when the app started (OV-1805). Showing another
// drawing, or a Load that replaces the shown one, applies nothing: the
// Result binding owns those.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { createSvgStyles } = await import('../../gbdraw/web/js/app/svg-styles.js');

/** A displayed Result with one drawn Depth group. */
const fakeResult = () => {
  const attrs = new Map();
  const depth = {
    getAttribute: (name) => (attrs.has(name) ? attrs.get(name) : null),
    setAttribute: (name, value) => { attrs.set(name, String(value)); },
    removeAttribute: (name) => { attrs.delete(name); }
  };
  const svg = { querySelectorAll: (selector) => (selector.includes('depth') ? [depth] : []) };
  return { depth, container: { querySelector: (selector) => (selector === 'svg' ? svg : null) } };
};

test('Show Depth on the shown drawing applies after a Load replaced the drawings', async () => {
  const originalList = state.drawingList();
  const originalTabs = state.drawingTabs.value;
  const originalShown = state.activeDrawingId.value;
  const { depth, container } = fakeResult();
  /** @type {string[]} */
  const commits = [];
  createSvgStyles({
    state: { ...state, svgContent: Vue.ref('<svg/>'), svgContainer: Vue.ref(container) },
    watch: Vue.watch,
    nextTick: Vue.nextTick,
    commitActiveResultEdit: (reason) => { commits.push(reason); return true; },
    projectPaletteAndRules: () => true,
    readDrawnTracks: () => []
  });
  try {
    // A Load: new drawing states under the same IDs, the circular one shown.
    const loaded = [createDrawingState('circular', 'circular'), createDrawingState('linear', 'linear')];
    loaded[0].form.show_depth = true;
    state.replaceDrawings(loaded, originalTabs);
    state.activeDrawingId.value = 'circular';
    await Vue.nextTick();
    assert.deepEqual(commits, [], 'a Load that replaces the shown drawing applies nothing');

    state.activeDrawing().form.show_depth = false;
    await Vue.nextTick();
    assert.equal(depth.getAttribute('display'), 'none', 'Show Depth off hides the drawn Depth track');
    assert.deepEqual(commits, ['track-visibility']);

    state.activeDrawing().form.show_depth = true;
    await Vue.nextTick();
    assert.equal(depth.getAttribute('display'), null, 'Show Depth on shows it again');

    // Showing another drawing applies nothing; a toggle on a drawing not shown neither.
    commits.length = 0;
    state.activeDrawingId.value = 'linear';
    await Vue.nextTick();
    state.drawingOf('circular').form.show_depth = false;
    await Vue.nextTick();
    assert.deepEqual(commits, []);
  } finally {
    state.replaceDrawings(originalList, originalTabs);
    state.activeDrawingId.value = originalShown;
  }
});
