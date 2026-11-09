// Guard: the LOSAT cache-info rows through every drawing-list step and History (#1777 rounds 14 and 15, after reviews
// 8 to 12 each found a defect in them). A deterministic walk (fixed seeds, bounded steps) over a Generate (one History
// step that writes the shown drawing's rows: Circular ring rows, Linear window-pair rows), New, Duplicate, Delete, a
// drawing switch, Clear LOSAT results, Undo and Redo, in the real app (`createAppSetup`, its drawing list and History).
// After every step it checks that:
//   (a) every row names a drawing of the project (a drawing's rows are in the state only while it is, OV-3209);
//   (b) a drawing holds one row per cache key (a Duplicate copies its source's rows once, OV-3205);
//   (c) every row's key is a cached search (also after a Clear and the Undo of a Delete before it);
//   (d) Undo then Redo of a step the walk recorded (a Generate, New, Duplicate or Delete) gives back the same rows;
//   (f) the shown Linear drawing's per-pair readers list each of its edges once, from its own rows (OV-3210);
//   (g) Clear LOSAT results counts the searches it removes, so it is enabled, and the "No LOSAT results yet" hint
//       hides, exactly while one is cached (OV-3211), and no fewer than the rows name.
//   (h) an Undo or a Redo changes the rows of no drawing but the one it shows (OV-3207, OV-3212);
//   (i) a drawing that comes back (an Undo of its Delete, a Redo of its New or Duplicate) has the rows it left with,
//       or none after a Clear (OV-3209, OV-3212: another drawing's Generate and its Undo in between take none away).
// An Undo or a Redo of the walk is not taken back: the next Undo or Redo may be another drawing's step.
// The walk's Generate writes as Generate commits (the other drawings' rows stay), so its Undo and Redo put back the
// drawing's rows only (OV-3207, OV-3212). Save and Load (check (e)) are not walked: a Save writes a drawing's entries
// only where its committed request names them, which needs a Generate this harness does not run; the two-holder Load
// cases of `mode-scoped-session.test.mjs` pin the rows a Load gives.
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

/**
 * A Generate of the shown drawing as one History step, as Generate commits: its rows replace its own, the other
 * drawings' stay, and its searches join the cache.
 * @param {import('../../gbdraw/web/js/services/config.js').LosatCacheInfoEntry[]} rows
 */
const generateRows = (rows) => window.__GBDRAW_HISTORY__.runUndoableArtifactReplacement('Generate diagram', () => {
  const drawing = state.activeDrawing();
  state.losatCache.value = new Map([...state.losatCache.value, ...rows.map((row) => [row.key, { key: row.key }])]);
  state.losatCacheInfo.value = [...losatRowsOfOtherDrawings(state.losatCacheInfo.value, drawing.id),
    ...rows.map((row) => ({ ...row, drawingId: drawing.id }))];
  return { status: 'ok' };
}, { shouldCommit: (outcome) => outcome?.status === 'ok' });
/** @param {string} key */
const ringRow = (key) => ({ key, filename: `${key}.tsv`, display: true, drawingId: '' });
/** @param {string} drawingId */
const keysOf = (drawingId) => state.losatCacheInfo.value.filter((row) => row.drawingId === drawingId).map(({ key }) => key);
const cachedKeys = () => [...state.losatCache.value.keys()].sort();
const resetLosat = () => { state.losatCacheInfo.value = []; state.losatCache.value = new Map(); };

// OV-3207 (review 9's order): History's Undo of a Generate put back the project's rows and searches as they were
// before it, so another drawing's later Generate lost its results. It puts back the drawing's own rows only, and the
// searches a row still names.
test('Undo of one drawing\'s Generate keeps another drawing\'s LOSAT results; Redo gives its own back (OV-3207)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(a); await generateRows([ringRow('k1')]);
    app.selectDrawing(b); await generateRows([ringRow('k2')]);
    app.selectDrawing(a); await history.undo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k2']],
      'B keeps its result; the search no row names any more goes (PD-OI-022)');
    await history.redo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [['k1'], ['k2'], ['k1', 'k2']]);
  } finally { resetLosat(); }
});

// OV-3212 (review 12's order): a drawing whose Delete is undone between another drawing's Generate and that
// Generate's Undo keeps its rows, and so does one still deleted while that Undo runs (its Delete holds them, and
// the Generate's step has their searches); History takes a newer Delete back before an older Generate.
test('a deleted drawing\'s LOSAT results survive another drawing\'s Generate and its Undo, in either order (OV-3212)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  for (const deleteUndoneFirst of [true, false]) {
    resetLosat();
    try {
      const a = state.activeDrawingId.value;
      const { drawingId: b } = app.drawingList.newDrawing('circular');
      const { drawingId: c } = app.drawingList.newDrawing('circular');
      app.selectDrawing(b); await generateRows([ringRow('k2')]);
      app.selectDrawing(a);
      app.drawingList.deleteDrawing(b);
      app.selectDrawing(a); await generateRows([ringRow('k1')]);
      if (deleteUndoneFirst) {
        app.selectDrawing(c); await history.undo();
        app.selectDrawing(a); await history.undo();
      } else {
        await history.undo();
        assert.deepEqual(cachedKeys(), ['k2'], 'the search the Delete holds stays cached');
        await history.undo();
      }
      assert.ok(state.drawingTabs.value.some((tab) => tab.id === b), 'B is back');
      assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k2']], `Delete undone first: ${deleteUndoneFirst}`);
      app.drawingList.deleteDrawing(c);
    } finally { resetLosat(); }
  }
  resetLosat();
  try {
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(a); await generateRows([ringRow('k1')]);
    app.selectDrawing(b); await generateRows([ringRow('k2')]);
    app.selectDrawing(a);
    app.drawingList.deleteDrawing(b);
    app.selectDrawing(a);
    // The rows a Delete holds name searches every step it lets History install has: the Undo in A takes the newer
    // Delete back first, so A's Generate is undone with B's rows in the project.
    await history.undo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [['k1'], ['k2'], ['k1', 'k2']], 'the Delete is undone first');
    await history.undo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k2']], 'then A\'s Generate');
  } finally { resetLosat(); }
});

// OV-3211: Clear LOSAT results removes every cached search, one a Delete holds for its Undo included, so its count
// and both panels' hints count those searches: Clear is enabled and the hint is true while one is cached.
test('Clear counts every search it removes, the one a deleted drawing holds included (OV-3211)', () => {
  const app = createAppSetup();
  resetLosat();
  try {
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    state.losatCacheInfo.value = [{ ...ringRow('k1'), drawingId: a }, { ...ringRow('k2'), drawingId: b }];
    state.losatCache.value = new Map([['k1', { key: 'k1' }], ['k2', { key: 'k2' }]]);
    app.drawingList.deleteDrawing(b);
    assert.equal(app.clearableLosatResultCount.value, 2, 'S2: Clear removes 2 searches');
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map([['k2', { key: 'k2' }]]);
    app.selectDrawing(a);
    assert.deepEqual([app.clearableLosatResultCount.value, app.losatResultRows.value.length], [1, 0],
      'S3: Clear is enabled, and the shown drawing has no results ("No LOSAT results in this drawing")');
    app.clearLosatCache();
    assert.equal(app.clearableLosatResultCount.value, 0, '"No LOSAT results yet", Clear disabled');
  } finally { resetLosat(); }
});

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
    assert.equal(app.clearableLosatResultCount.value, state.losatCache.value.size, `${where} (g): Clear's count`);
    assert.ok(app.clearableLosatResultCount.value >= new Set(rows.map(({ key }) => key)).size, `${where} (g): uncounted`);
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
  /** A Generate of the shown drawing that searched `count` pairs, one History step. */
  const generate = (count) => {
    const drawing = state.activeDrawing();
    return generateRows(Array.from({ length: count }, (_, index) => {
      search += 1;
      const key = `search-${search}`;
      return drawing.mode === 'linear'
        ? { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id, edgeKey: `edge-${search}`,
          queryUid: `a#${(index % 2) + 1}`, subjectUid: 'b' }
        : { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id };
    }));
  };
  // (h) and (i): a drawing that leaves the project (a Delete, the Undo of a New or Duplicate) leaves with its rows, and
  // the same drawing coming back (by Undo or Redo) brings them back, or none after a Clear (its searches are gone).
  /** @type {Map<object, { rows: unknown, clears: number }>} */
  const away = new Map();
  let clears = 0;
  let returns = 0;
  const lives = () => new Map(tabIds().map((id) => [state.drawingOf(id), id]));
  /** @param {string} where @param {() => unknown} move @param {{ undoOrRedo?: boolean }} [options] */
  const tracked = async (where, move, { undoOrRedo = false } = {}) => {
    const before = { rows: JSON.parse(rowsByDrawing()), lives: lives() };
    await move();
    const after = { rows: JSON.parse(rowsByDrawing()), lives: lives() };
    before.lives.forEach((id, drawing) => {
      if (!after.lives.has(drawing)) away.set(drawing, { rows: before.rows[id], clears });
    });
    after.lives.forEach((id, drawing) => {
      const left = away.get(drawing);
      if (before.lives.has(drawing) || !left) return;
      returns += 1;
      away.delete(drawing);
      assert.deepEqual(after.rows[id], left.clears === clears ? left.rows : [], `${where} (i): ${id} comes back with its rows`);
    });
    if (!undoOrRedo) return;
    const changed = Object.keys(after.rows).filter((id) => id !== state.activeDrawingId.value && id in before.rows
      && JSON.stringify(after.rows[id]) !== JSON.stringify(before.rows[id]));
    assert.deepEqual(changed, [], `${where} (h): an Undo or Redo changed another drawing's rows`);
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
        const at = `seed ${seed} step ${step} [${trail.slice(-5).join(' ')} ${kind}]`;
        // Whether the step recorded a History step, which (d) takes back and does again.
        let recorded = false;
        await tracked(at, async () => {
          if (kind === 'generate') recorded = (await generate(1 + Math.floor(next() * 3)))?.status === 'ok';
          else if (kind === 'new') recorded = app.drawingList.newDrawing(pick(['circular', 'linear'])).status === 'ok';
          else if (kind === 'duplicate') recorded = app.drawingList.duplicateDrawing(pick(ids)).status === 'ok';
          else if (kind === 'delete' && ids.length > 1) recorded = app.drawingList.deleteDrawing(pick(ids)).status === 'ok';
          else if (kind === 'switch') app.selectDrawing(pick(ids));
          else if (kind === 'clear') { app.clearLosatCache(); clears += 1; }
          else if ((kind === 'undo' || kind === 'redo') && timeline[kind].length > 0) await history[kind]();
        }, { undoOrRedo: kind === 'undo' || kind === 'redo' });
        trail.push(`${kind}@${state.activeDrawingId.value}`);
        const where = `seed ${seed} step ${step} [${trail.slice(-6).join(' ')}]`;
        check(where);
        if (recorded) {
          const before = rowsByDrawing();
          await tracked(`${where} undo`, () => history.undo(), { undoOrRedo: true });
          check(`${where} undo`);
          await tracked(`${where} redo`, () => history.redo(), { undoOrRedo: true });
          assert.equal(rowsByDrawing(), before, `${where} (d): Undo then Redo gives back the same rows`);
          check(`${where} redo`);
        }
      }
    }
    assert.ok(returns > 100, `the walk brought drawings back ${returns} times`);
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
  }
  console.log(`LOSAT rows walk: ${SEEDS.length} seeds x ${STEPS} steps, ${returns} returns, in ${Date.now() - started} ms`);
});
