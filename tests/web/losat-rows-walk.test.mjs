// Guard: the LOSAT cache-info rows through every drawing-list step (#1777 round 14, after reviews 8 to 11 each found
// a defect in them). A deterministic walk (fixed seeds, bounded steps) over a Generate that writes the shown drawing's
// rows (Circular ring rows, Linear window-pair rows), New, Duplicate, Delete, a drawing switch, Clear LOSAT results,
// Undo and Redo, in the
// real app (`createAppSetup`, its drawing list and History). After every step it checks that:
//   (a) every row names a drawing of the project (a drawing's rows are in the state only while it is, OV-3209);
//   (b) a drawing holds one row per cache key (a Duplicate copies its source's rows once, OV-3205);
//   (c) every row's key is a cached search (also after a Clear and the Undo of a Delete before it);
//   (d) Undo then Redo of the last step gives back the same rows;
//   (f) the shown Linear drawing's per-pair readers list each of its edges once, from its own rows (OV-3210).
// The walk's Generate writes as the Circular Generate commits (the other drawings' rows stay) and is not a History
// step: History's Undo of a Generate restores the project's rows as they were before it (OV-3207, same on `dev`), so
// the follow-up that makes it restore the drawing's rows only adds the Generate's Undo here. Save and Load (check (e))
// are not walked: a Save writes a drawing's entries only where its committed request names them, which needs a
// Generate this harness does not run; the two-holder Load cases of `mode-scoped-session.test.mjs` pin the rows a
// Load gives.
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
const { losatRowsOfOtherDrawings } = await import('../../gbdraw/web/js/services/losat-cache.js');

/** @param {number} seed */
const random = (seed) => {
  let value = seed >>> 0;
  return () => {
    value = (value + 0x6D2B79F5) >>> 0;
    let next = Math.imul(value ^ (value >>> 15), 1 | value);
    next = (next + Math.imul(next ^ (next >>> 7), 61 | next)) ^ next;
    return ((next ^ (next >>> 14)) >>> 0) / 4294967296;
  };
};

const SEEDS = [1777, 3205, 3209, 3210, 42, 7, 11, 2026];
const STEPS = 120;

// The drawn Result's units: record `a` in two windows, record `b` whole, so each Linear row is a window pair.
const CATALOG = { items: [{ recordKeys: ['a#1', 'a#2', 'b'], windowOf: { 'a#1': 'a', 'a#2': 'a' } }] };

test('a drawing\'s LOSAT rows follow it through Generate, New, Duplicate, Delete, a switch, Clear, Undo and Redo', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const tabIds = () => state.drawingTabs.value.map((tab) => tab.id);
  /** The rows as one comparable value: per drawing, its rows' key, edge and file name, sorted. */
  const rowsByDrawing = () => JSON.stringify(Object.fromEntries(tabIds().map((id) => [id, state.losatCacheInfo.value
    .filter((row) => row.drawingId === id)
    .map(({ key, edgeKey = null, filename }) => [key, edgeKey, filename])
    .sort((left, right) => JSON.stringify(left).localeCompare(JSON.stringify(right)))])));
  let search = 0;
  /** @param {string} where */
  const check = (where) => {
    const rows = state.losatCacheInfo.value;
    const ids = new Set(tabIds());
    const strays = rows.filter((row) => !ids.has(row.drawingId)).map((row) => `${row.drawingId}:${row.key}`);
    assert.deepEqual(strays, [], `${where} (a): rows of a drawing the project does not have`);
    const pairs = rows.map((row) => `${row.drawingId}:${row.key}`);
    assert.equal(new Set(pairs).size, pairs.length, `${where} (b): one row per drawing and key: ${pairs}`);
    const unsearched = rows.filter((row) => !state.losatCache.value.has(row.key)).map((row) => row.key);
    assert.deepEqual(unsearched, [], `${where} (c): rows whose search is not cached`);
    const shown = state.activeDrawing();
    if (shown.mode !== 'linear') return;
    const own = rows.filter((row) => row.drawingId === shown.id).map((row) => row.edgeKey).sort();
    // The units come from the shown Result's catalog, which a switch replaces: the walk's Results all draw CATALOG.
    state.featureCatalog.value = CATALOG;
    const listed = Object.values(app.linearLosatWindowPairs.value).flat().map(({ edgeKey }) => edgeKey).sort();
    assert.deepEqual(listed, own, `${where} (f): the window pairs the shown Linear drawing lists`);
    const byEdge = Object.values(app.linearLosatCacheInfoByEdgeKey.value);
    assert.deepEqual(byEdge.filter((row) => row.drawingId !== shown.id), [], `${where} (f): another drawing's row`);
  };
  /** A Generate of the shown drawing that searched `count` pairs: its rows replace its own, the others stay. */
  const generate = (count) => {
    const drawing = state.activeDrawing();
    const rows = Array.from({ length: count }, (_, index) => {
      search += 1;
      const key = `search-${search}`;
      return drawing.mode === 'linear'
        ? { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id, edgeKey: `edge-${search}`,
          queryUid: `a#${(index % 2) + 1}`, subjectUid: 'b' }
        : { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id };
    });
    state.losatCache.value = new Map([...state.losatCache.value, ...rows.map((row) => [row.key, { key: row.key }])]);
    state.losatCacheInfo.value = [...losatRowsOfOtherDrawings(state.losatCacheInfo.value, drawing.id), ...rows];
  };
  const started = Date.now();
  try {
    for (const seed of SEEDS) {
      const next = random(seed);
      /** @template T @param {readonly T[]} items @returns {T} */
      const pick = (items) => items[Math.floor(next() * items.length)];
      const trail = [];
      for (let step = 0; step < STEPS; step += 1) {
        const ids = tabIds();
        const kind = pick(['generate', 'generate', 'new', 'duplicate', 'delete', 'switch', 'undo', 'redo', 'clear']);
        const timeline = history.timeline();
        let listStep = false;
        if (kind === 'generate') generate(1 + Math.floor(next() * 3));
        else if (kind === 'new') listStep = app.drawingList.newDrawing(pick(['circular', 'linear'])).status === 'ok';
        else if (kind === 'duplicate') listStep = app.drawingList.duplicateDrawing(pick(ids)).status === 'ok';
        else if (kind === 'delete' && ids.length > 1) listStep = app.drawingList.deleteDrawing(pick(ids)).status === 'ok';
        else if (kind === 'switch') app.selectDrawing(pick(ids));
        else if (kind === 'clear') app.clearLosatCache();
        else if (kind === 'undo' && timeline.undo.length > 0) { await history.undo(); listStep = true; }
        else if (kind === 'redo' && timeline.redo.length > 0) { await history.redo(); listStep = true; }
        trail.push(`${kind}@${state.activeDrawingId.value}`);
        const where = `seed ${seed} step ${step} [${trail.slice(-6).join(' ')}]`;
        check(where);
        if (listStep && history.timeline().undo.length > 0) {
          const before = rowsByDrawing();
          await history.undo();
          check(`${where} undo`);
          await history.redo();
          assert.equal(rowsByDrawing(), before, `${where} (d): Undo then Redo gives back the same rows`);
          check(`${where} redo`);
        }
      }
    }
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
  }
  console.log(`LOSAT rows walk: ${SEEDS.length} seeds x ${STEPS} steps in ${Date.now() - started} ms`);
});
