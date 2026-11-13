// Guard: the LOSAT cache-info rows and evidence through every drawing-list step and History (#1777 rounds 14 to 17,
// after reviews 8 to 13 each found a defect in them). A deterministic walk (fixed seeds, bounded steps) over a Generate
// (one History step that writes the shown drawing's rows: Circular ring rows; Linear window-pair rows with their protein
// searches, one all-against-all search without a row, the converted payload, and the drawing's protein identity
// manifest, whose record instance keys every Linear Generate shares with other meanings), New, Duplicate, Delete, a
// drawing switch, Clear LOSAT results, Undo, Redo, and review 13's orders (New or Duplicate, a Generate in the new
// drawing, a switch back, Undo, Undo, Redo, Redo), in the real app (`createAppSetup`, its drawing list and History).
// After every step it checks that:
//   (a) every row names a drawing of the project (a drawing's rows are in the state only while it is, OV-3209);
//   (b) a drawing holds one row per cache key (a Duplicate copies its source's rows once, OV-3205);
//   (c) every row's key is a cached search (also after a Clear and the Undo of a Delete before it);
//   (d) Undo then Redo of a step the walk recorded (a Generate, New, Duplicate or Delete) gives back the same rows;
//   (f) the shown Linear drawing's per-pair readers list each of its edges once, from its own rows (OV-3210);
//   (g) Clear LOSAT results counts the shown drawing's rows, one per key; after a Clear the shown drawing has none and
//       the cache holds only the searches the other drawings' rows name (or their manifests resolve) (OV-3211);
//   (h) an Undo or a Redo changes the rows of no drawing but the one it shows (OV-3207, OV-3212);
//   (i) a drawing that comes back (an Undo of its Delete, a Redo of its New or Duplicate) has the rows it left with,
//       or none after a Clear (OV-3209, OV-3212: another drawing's Generate and its Undo in between take none away),
//       and after a Clear an empty protein identity manifest, as the drawings Clear found (OV-3274);
//   (j) the shown drawing's evidence resolves: each protein row's search names record instances of the drawing's own
//       manifest (OV-3213, OV-3214), and the all-against-all search of the Generate its rows come from is cached;
//   (k') no step but Clear removes a cached search or converted payload (D-200, OV-3260: a later Generate reuses
//       them), and every converted payload's searches are cached.
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
const { losatRowsOfOtherDrawings, validateProteinRawEntryReferences, emptyProteinIdentityManifest } =
  await import('../../gbdraw/web/js/services/losat-cache.js');
const { LOSAT_DERIVED_CACHE_SCHEMA, PROTEIN_IDENTITY_MANIFEST_SCHEMA, PROTEIN_LOSAT_CACHE_SCHEMA } =
  await import('../../gbdraw/web/js/contracts/schema.generated.js');

// A Linear protein Generate's evidence, as its extraction keys it: one record instance per card uid (a Duplicate
// copies the uids), whose bindings (`meaning`: the cards' regions, visibility and names when it ran) differ between
// Generates; a runtime handle follows the feature and the uid only, as Python's `build_protein_runtime_handle`.
/** @param {string} uid @param {number} feature */
const handleOf = (uid, feature) => `h_${[...`${uid}x${feature}`].map((c) => 'abcdefghijklmnopqrstuvwxyz'[c.charCodeAt(0) % 26]).join('').padEnd(26, 'a')}`;
/** @param {number} feature */
const featureId = (feature) => `f_${String(feature).padStart(64, '0')}`;
/** @param {string} meaning @param {readonly string[]} uids */
const proteinManifest = (meaning, uids) => ({
  schema: PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  proteinSets: Object.fromEntries(uids.map((uid) => [`set-${meaning}-${uid}`, { schema: 1, proteins: [1, 2].map((n) => ({ featureAnalysisId: featureId(n) })) }])),
  recordAnalyses: Object.fromEntries(uids.map((uid) => [`ra-${meaning}-${uid}`, { schema: 1, proteinSetHash: `set-${meaning}-${uid}` }])),
  recordInstances: Object.fromEntries(uids.map((uid) => [uid, {
    schema: 2, recordAnalysisId: `ra-${meaning}-${uid}`, runtimeBindingHash: `bind-${meaning}-${uid}`, displayBindingHash: `disp-${meaning}-${uid}`,
    runtimeIds: Object.fromEntries([1, 2].map((n) => [featureId(n), handleOf(uid, n)])),
    featureMetadata: Object.fromEntries([1, 2].map((n) => [featureId(n), { displayAlias: `p${n}` }]))
  }]))
});
/** A protein raw search of `query` against `subject` under the bindings `meaning` names. @param {string} meaning @param {string} query @param {string} subject */
const proteinSearch = (meaning, query, subject) => ({
  kind: 'raw-losat', schema: PROTEIN_LOSAT_CACHE_SCHEMA, identityKind: 'protein', idEncoding: 'runtime-handle-v1', program: 'blastp',
  text: `${handleOf(query, 1)}\t${handleOf(subject, 2)}\t90.0\t10\t1\t0\t1\t10\t1\t10\t1e-5\t50\n`,
  queryProteinSetHash: `set-${meaning}-${query}`, subjectProteinSetHash: `set-${meaning}-${subject}`,
  queryRuntimeBindingHash: `bind-${meaning}-${query}`, subjectRuntimeBindingHash: `bind-${meaning}-${subject}`,
  queryRecordInstanceKey: query, subjectRecordInstanceKey: subject
});
/** A converted payload of a protein Generate, keyed by the raw searches it read. @param {string} key @param {readonly string[]} rawKeys */
const derivedEvidence = (key, rawKeys) => ({
  schema: LOSAT_DERIVED_CACHE_SCHEMA, kind: 'derived-losatp-payload', idEncoding: 'runtime-handle-v1', key, mode: 'pairwise', payload: {}, rawKeys
});

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
const MAX_DRAWINGS = 5;

// The drawn Result's units: record `a` in two windows, record `b` whole, so each Linear row is a window pair.
const CATALOG = { items: [{ recordKeys: ['a#1', 'a#2', 'b'], windowOf: { 'a#1': 'a', 'a#2': 'a' } }] };

/**
 * A Generate of the shown drawing as one History step, as Generate commits: its rows replace its own, the other
 * drawings' stay, and its searches join the cache; a protein Generate also writes its drawing's manifest and its
 * converted payloads.
 * @param {import('../../gbdraw/web/js/services/config.js').LosatCacheInfoEntry[]} rows
 * @param {{ searches?: ReadonlyMap<string, object>, manifest?: object, derived?: readonly { key: string }[] }} [evidence]
 */
const generateRows = (rows, { searches = new Map(), manifest, derived = [] } = {}) => window.__GBDRAW_HISTORY__.runUndoableArtifactReplacement('Generate diagram', () => {
  const drawing = state.activeDrawing();
  state.losatCache.value = new Map([...state.losatCache.value, ...rows.map((row) => [row.key, { key: row.key }]), ...searches]);
  state.losatCacheInfo.value = [...losatRowsOfOtherDrawings(state.losatCacheInfo.value, drawing.id),
    ...rows.map((row) => ({ ...row, drawingId: drawing.id }))];
  if (manifest) state.proteinIdentityManifest.value = manifest;
  if (derived.length) state.losatDerivedCache.value = new Map([...state.losatDerivedCache.value, ...derived.map((entry) => [entry.key, entry])]);
  return { status: 'ok' };
}, { shouldCommit: (outcome) => outcome?.status === 'ok' });
/** @param {string} key */
const ringRow = (key) => ({ key, filename: `${key}.tsv`, display: true, drawingId: '' });
/** @param {string} drawingId */
const keysOf = (drawingId) => state.losatCacheInfo.value.filter((row) => row.drawingId === drawingId).map(({ key }) => key);
const cachedKeys = () => [...state.losatCache.value.keys()].sort();
const resetLosat = () => {
  state.losatCacheInfo.value = []; state.losatCache.value = new Map(); state.losatDerivedCache.value = new Map();
  state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
};
/** Every drawing but the shown one leaves (a case starts from one drawing). @param {ReturnType<typeof createAppSetup>} app */
const onlyShownDrawing = (app) => state.drawingTabs.value.map((tab) => tab.id)
  .filter((id) => id !== state.activeDrawingId.value).forEach((id) => app.drawingList.deleteDrawing(id));

// OV-3207 (review 9's order): History's Undo of a Generate put back the project's rows and searches as they were
// before it, so another drawing's later Generate lost its results. It puts back the drawing's own rows only, and no
// search goes (D-200: a later Generate reuses it, OV-3260).
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
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k1', 'k2']],
      'B keeps its result; A\'s search stays for its next Generate');
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
        assert.deepEqual(cachedKeys(), ['k1', 'k2'], 'the search the Delete holds stays cached');
        await history.undo();
      }
      assert.ok(state.drawingTabs.value.some((tab) => tab.id === b), 'B is back');
      assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k1', 'k2']], `Delete undone first: ${deleteUndoneFirst}`);
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
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['k2'], ['k1', 'k2']], 'then A\'s Generate');
  } finally { resetLosat(); }
});

// OV-3212 (review 13's P1 and P2): a drawing whose New or Duplicate is undone holds its rows; the Undo of another
// drawing's older Generate in between keeps every cached search (D-200), so the Redo of the New or Duplicate brings
// the drawing back with them.
test('a drawing held by an undone New or Duplicate keeps its searches through another drawing\'s older Generate\'s Undo and Redo (OV-3212)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  for (const make of ['new', 'duplicate']) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const a = state.activeDrawingId.value;
      await generateRows([ringRow('k1')]);
      const { drawingId: b } = make === 'new' ? app.drawingList.newDrawing('circular') : app.drawingList.duplicateDrawing(a);
      app.selectDrawing(b); await generateRows([ringRow('k2')]);
      app.selectDrawing(a);
      await history.undo();
      assert.ok(!state.drawingTabs.value.some((tab) => tab.id === b), `${make}: B is held`);
      await history.undo();
      assert.deepEqual([keysOf(a), cachedKeys()], [[], ['k1', 'k2']], `${make}: the search B's held rows name stays`);
      await history.redo();
      await history.redo();
      assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [['k1'], ['k2'], ['k1', 'k2']], `${make}: B comes back with its rows`);
    } finally { resetLosat(); }
  }
});

// OV-3215 (review 13's P3): after Clear, the Undo of one drawing's Generate brings back that drawing's rows from
// before it, with their searches, and no search that no row names (PD-OI-022); another drawing's results stay
// (scenario revision 2).
test('Clear then the Undo of a Generate brings back no search that no row names (OV-3215)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    onlyShownDrawing(app);
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(b); await generateRows([ringRow('m2')]);
    app.selectDrawing(a); await generateRows([ringRow('m1')]);
    app.clearLosatCache();
    await history.undo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys(), app.losatResultRows.value.length], [[], ['m2'], ['m2'], 0]);
  } finally { resetLosat(); }
});

// OV-3274 (review 14): Clear empties every drawing's protein identity manifest (D-197), also that of a drawing History
// holds out of the project (a Delete, the Undo of a New or Duplicate): it comes back with an empty one and no rows.
test('a drawing that comes back after Clear has an empty protein identity manifest and no LOSAT rows (OV-3274)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  /** @param {string} key */
  const pairRow = (key) => ({ key, filename: `${key}.tsv`, display: true, drawingId: '', edgeKey: `edge-${key}`, queryUid: 'u1', subjectUid: 'u2' });
  for (const leave of ['delete', 'undo-new', 'undo-duplicate']) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const { drawingId: a } = app.drawingList.newDrawing('linear');
      onlyShownDrawing(app);
      const made = leave === 'delete' ? app.drawingList.newDrawing('linear')
        : leave === 'undo-new' ? app.drawingList.newDrawing('linear') : app.drawingList.duplicateDrawing(a);
      const b = made.drawingId;
      app.selectDrawing(b);
      await generateRows([pairRow('kB')], { searches: new Map([['kB', proteinSearch('b', 'u1', 'u2')]]), manifest: proteinManifest('b', ['u1', 'u2']) });
      app.selectDrawing(a);
      if (leave === 'delete') app.drawingList.deleteDrawing(b);
      else await history.undo();
      assert.ok(!state.drawingTabs.value.some((tab) => tab.id === b), `${leave}: B is held`);
      app.clearLosatCache();
      if (leave === 'delete') await history.undo();
      else await history.redo();
      app.selectDrawing(b);
      assert.deepEqual([keysOf(b), state.proteinIdentityManifest.value], [[], emptyProteinIdentityManifest()],
        `${leave}: B comes back with no rows and an empty manifest`);
    } finally { resetLosat(); }
  }
});

// OV-3214, OV-3213: each Linear drawing's protein pairs resolve through its own protein identity manifest, also when
// another drawing's protein Generate gave one card uid (a Duplicate copies the uids) another meaning, and after that
// Generate's Undo and Redo: a listed pair always downloads (the hydrate reads the shown drawing's manifest).
test('each Linear drawing\'s protein pairs resolve in its own manifest after another drawing\'s protein Generate, its Undo and Redo (OV-3214)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  /** @param {string} key */
  const pairRow = (key) => ({ key, filename: `${key}.tsv`, display: true, drawingId: '', edgeKey: `edge-${key}`, queryUid: 'u1', subjectUid: 'u2' });
  /** @param {string} drawingId @param {string} key */
  const resolves = (drawingId, key) => {
    app.selectDrawing(drawingId);
    return validateProteinRawEntryReferences(state.losatCache.value.get(key), state.proteinIdentityManifest.value);
  };
  try {
    onlyShownDrawing(app);
    const { drawingId: a } = app.drawingList.newDrawing('linear');
    onlyShownDrawing(app);
    await generateRows([pairRow('kA')], { searches: new Map([['kA', proteinSearch('whole', 'u1', 'u2')]]), manifest: proteinManifest('whole', ['u1', 'u2']) });
    const { drawingId: d } = app.drawingList.duplicateDrawing(a);
    app.selectDrawing(d);
    await generateRows([pairRow('kD')], { searches: new Map([['kD', proteinSearch('crop', 'u1', 'u2')]]), manifest: proteinManifest('crop', ['u1', 'u2']) });
    assert.deepEqual([resolves(a, 'kA'), resolves(d, 'kD')], [true, true], 'after the copy\'s protein Generate');
    app.selectDrawing(d); await history.undo();
    assert.deepEqual([keysOf(d), resolves(d, 'kA'), resolves(a, 'kA')], [['kA'], true, true], 'after its Undo');
    app.selectDrawing(d); await history.redo();
    assert.deepEqual([resolves(d, 'kD'), resolves(a, 'kA')], [true, true], 'after its Redo');
  } finally { resetLosat(); }
});

// OV-3211, scenario revision 2: the count is the shown drawing's results, and the hint says whether another drawing
// has some. A drawing a Delete holds out of the project is no drawing of it: a Clear keeps none of its searches, and
// the Undo of its Delete brings it back without LOSAT results (OV-3274).
test('Clear counts the shown drawing\'s results; a deleted drawing comes back without the searches a Clear removed (OV-3211)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    onlyShownDrawing(app);
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    await generateRows([ringRow('k2')]);
    app.selectDrawing(a);
    assert.deepEqual([app.losatResultRows.value.length, app.projectHasLosatResults.value], [0, true],
      'Clear is disabled; the hint says "No LOSAT results in this drawing"');
    await generateRows([ringRow('k1')]);
    app.drawingList.deleteDrawing(b);
    app.selectDrawing(a);
    assert.equal(app.losatResultRows.value.length, 1, '"Clear 1 LOSAT result"');
    app.clearLosatCache();
    assert.deepEqual([cachedKeys(), app.losatResultRows.value.length, app.projectHasLosatResults.value], [[], 0, false],
      'no drawing of the project names k2: "No LOSAT results yet"');
    await history.undo();
    assert.deepEqual([state.drawingTabs.value.some((tab) => tab.id === b), keysOf(b)], [true, []], 'B comes back without results');
  } finally { resetLosat(); }
});

// PD-OI-022 scenario revision 2 (Owner 2026-11-13, Q2=B): Clear LOSAT results clears the shown drawing's results. Its
// count is the shown drawing's results, one per cache key, and Clear takes every one of them from its list; a search
// another drawing names (a row of it, or a protein search that resolves in its manifest) stays for that drawing,
// whose list is unchanged; a search no remaining drawing names leaves the cache, with the converted payloads read
// from it.
test('Clear LOSAT results clears the shown drawing\'s results; a search another drawing names stays (PD-OI-022 rev 2)', async () => {
  const app = createAppSetup();
  /** @param {string} key */
  const pairRow = (key) => ({ key, filename: `${key}.tsv`, display: true, drawingId: '', edgeKey: `edge-${key}`, queryUid: 'u1', subjectUid: 'u2' });
  for (const mode of ['circular', 'linear']) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const { drawingId: a } = app.drawingList.newDrawing(mode);
      onlyShownDrawing(app);
      const { drawingId: b } = app.drawingList.newDrawing(mode);
      const linear = mode === 'linear';
      const row = linear ? pairRow : ringRow;
      // B: the shared search kS and its own kB; a Linear B also an all-against-all search that only its manifest names.
      app.selectDrawing(b);
      await generateRows([row('kS'), row('kB')], linear ? {
        searches: new Map([['kB', proteinSearch('b', 'u1', 'u2')], ['allB', proteinSearch('b', 'u2', 'u1')]]),
        manifest: proteinManifest('b', ['u1', 'u2']), derived: [derivedEvidence('pB', ['kB', 'allB'])]
      } : {});
      // A (shown): kS, its own kA, a search no row names (a failed Generate's) and, Linear, its own all-against-all.
      app.selectDrawing(a);
      await generateRows([row('kS'), row('kA')], linear ? {
        searches: new Map([['kA', proteinSearch('a', 'u1', 'u2')], ['allA', proteinSearch('a', 'u2', 'u1')], ['orphan', { key: 'orphan' }]]),
        manifest: proteinManifest('a', ['u1', 'u2']), derived: [derivedEvidence('pA', ['kA', 'allA']), derivedEvidence('pS', ['kS'])]
      } : { searches: new Map([['orphan', { key: 'orphan' }]]) });
      const bRowsBefore = state.losatCacheInfo.value.filter((entry) => entry.drawingId === b);
      assert.equal(app.losatResultRows.value.length, 2, `${mode}: "Clear 2 LOSAT results" counts the shown drawing's results`);
      app.clearLosatCache();
      assert.deepEqual([keysOf(a), app.losatResultRows.value.length], [[], 0], `${mode}: the shown drawing's list is empty`);
      assert.deepEqual(state.losatCacheInfo.value.filter((entry) => entry.drawingId === b), bRowsBefore, `${mode}: B's rows are unchanged`);
      assert.deepEqual(cachedKeys(), linear ? ['allB', 'kB', 'kS'] : ['kB', 'kS'],
        `${mode}: the searches B names stay; A's own and the one no row names go`);
      assert.deepEqual([...state.losatDerivedCache.value.keys()].sort(), linear ? ['pB', 'pS'] : [], `${mode}: payloads of kept searches stay`);
      assert.deepEqual(state.proteinIdentityManifest.value, emptyProteinIdentityManifest(), `${mode}: A's manifest is emptied`);
      app.selectDrawing(b);
      assert.equal(app.losatResultRows.value.length, 2, `${mode}: B still lists 2 results`);
      if (linear) assert.deepEqual(state.proteinIdentityManifest.value, proteinManifest('b', ['u1', 'u2']), 'B keeps its manifest');
      // Clear in B, now the only drawing with results: nothing is named any more, the cache empties.
      app.clearLosatCache();
      assert.deepEqual([state.losatCacheInfo.value, cachedKeys(), state.losatDerivedCache.value.size], [[], [], 0], `${mode}: B's Clear`);
    } finally { resetLosat(); }
  }
});

// A project of one drawing: Clear empties the whole cache, as before revision 2 (a search no row names included).
test('in a project of one drawing, Clear LOSAT results empties every LOSAT cache', async () => {
  const app = createAppSetup();
  resetLosat();
  try {
    onlyShownDrawing(app);
    await generateRows([ringRow('k1')], { searches: new Map([['orphan', { key: 'orphan' }]]), derived: [derivedEvidence('p1', ['k1'])] });
    app.clearLosatCache();
    assert.deepEqual([state.losatCacheInfo.value, cachedKeys(), state.losatDerivedCache.value.size, app.losatResultRows.value.length], [[], [], 0, 0]);
  } finally { resetLosat(); }
});

// Clear is no History step: an Undo after it takes back the shown drawing's last step, and History's install puts back
// the rows that step's drawing state had, with their searches (D-200); its Redo brings back the cleared rows.
test('Undo and Redo after Clear bring back the cleared drawing\'s rows with their searches', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    onlyShownDrawing(app);
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(b); await generateRows([ringRow('kS')]);
    app.selectDrawing(a);
    await generateRows([ringRow('k1')]);
    await generateRows([ringRow('kS'), ringRow('k2')]);
    app.clearLosatCache();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [[], ['kS'], ['kS']]);
    await history.undo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [['k1'], ['kS'], ['k1', 'kS']], 'the step before comes back');
    await history.redo();
    assert.deepEqual([keysOf(a), keysOf(b), cachedKeys()], [['kS', 'k2'], ['kS'], ['k1', 'k2', 'kS']], 'the cleared rows come back');
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
    assert.equal(app.losatResultRows.value.length, new Set(rows.filter((row) => row.drawingId === state.activeDrawingId.value)
      .map(({ key }) => key)).size, `${where} (g): the count is the shown drawing's results`);
    const shown = state.activeDrawing();
    const ownRows = rows.filter((row) => row.drawingId === shown.id);
    const unresolved = ownRows.filter(({ key }) => {
      const entry = state.losatCache.value.get(key);
      return entry?.identityKind === 'protein' && (!validateProteinRawEntryReferences(entry, state.proteinIdentityManifest.value)
        || !state.losatCache.value.has(`all-${generationOf.get(key)}`));
    }).map(({ key }) => key);
    assert.deepEqual(unresolved, [], `${where} (j): protein rows whose evidence does not resolve in the drawing`);
    const orphans = [...state.losatDerivedCache.value.values()].filter(({ rawKeys }) => !rawKeys.every((key) => state.losatCache.value.has(key)));
    assert.deepEqual(orphans.map(({ key }) => key), [], `${where} (k'): converted payloads whose searches are gone`);
    if (shown.mode !== 'linear') return;
    const own = rows.filter((row) => row.drawingId === shown.id).map((row) => row.edgeKey).sort();
    // The units come from the shown Result's catalog, which a switch replaces: the walk's Results all draw CATALOG.
    state.featureCatalog.value = CATALOG;
    const listed = Object.values(app.linearLosatWindowPairs.value).flat().map(({ edgeKey }) => edgeKey).sort();
    assert.deepEqual(listed, own, `${where} (f): the window pairs the shown Linear drawing lists`);
    const byEdge = Object.values(app.linearLosatCacheInfoByEdgeKey.value);
    assert.deepEqual(byEdge.filter((row) => row.drawingId !== shown.id), [], `${where} (f): another drawing's row`);
  };
  /** The Generate a row's search came from: a Linear Generate's searches and payload are named after it. */
  const generationOf = new Map();
  /** An Undo or a Redo. @param {'undo' | 'redo'} direction */
  const traverse = (direction) => history[direction]();
  /** A Generate of the shown drawing that searched `count` pairs, one History step. @param {number} count */
  const generate = (count) => {
    const drawing = state.activeDrawing();
    const meaning = `g${search + 1}`;
    const rows = Array.from({ length: count }, (_, index) => {
      search += 1;
      const key = `search-${search}`;
      generationOf.set(key, meaning);
      return drawing.mode === 'linear'
        ? { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id, edgeKey: `edge-${search}`,
          queryUid: `a#${(index % 2) + 1}`, subjectUid: 'b' }
        : { key, filename: `${key}.tsv`, display: true, drawingId: drawing.id };
    });
    if (drawing.mode !== 'linear') return generateRows(rows);
    const searches = new Map([...rows.map(({ key }) => [key, proteinSearch(meaning, 'u1', 'u2')]),
      [`all-${meaning}`, proteinSearch(meaning, 'u2', 'u1')]]);
    return generateRows(rows, { searches, manifest: proteinManifest(meaning, ['u1', 'u2']),
      derived: [derivedEvidence(`payload-${meaning}`, [...searches.keys()])] });
  };
  // (h) and (i): a drawing that leaves the project (a Delete, the Undo of a New or Duplicate) leaves with its rows, and
  // the same drawing coming back (by Undo or Redo) brings them back, or none after a Clear (its searches are gone).
  /** @type {Map<object, { rows: unknown, clears: number }>} */
  const away = new Map();
  let clears = 0;
  let returns = 0;
  let clearedReturns = 0;
  const lives = () => new Map(tabIds().map((id) => [state.drawingOf(id), id]));
  /** The cached searches' and converted payloads' keys. */
  const cacheKeys = () => [...[...state.losatCache.value.keys()].map((key) => `search:${key}`),
    ...[...state.losatDerivedCache.value.keys()].map((key) => `payload:${key}`)];
  /** @param {string} where @param {() => unknown} move @param {{ undoOrRedo?: boolean }} [options] */
  const tracked = async (where, move, { undoOrRedo = false } = {}) => {
    const before = { rows: JSON.parse(rowsByDrawing()), lives: lives(), clears, caches: cacheKeys() };
    await move();
    const after = { rows: JSON.parse(rowsByDrawing()), lives: lives() };
    const kept = new Set(cacheKeys());
    const gone = before.clears === clears ? before.caches.filter((key) => !kept.has(key)) : [];
    assert.deepEqual(gone, [], `${where} (k'): a step that is not Clear removed cached evidence`);
    before.lives.forEach((id, drawing) => {
      if (!after.lives.has(drawing)) away.set(drawing, { rows: before.rows[id], clears });
    });
    after.lives.forEach((id, drawing) => {
      const left = away.get(drawing);
      if (before.lives.has(drawing) || !left) return;
      returns += 1;
      away.delete(drawing);
      assert.deepEqual(after.rows[id], left.clears === clears ? left.rows : [], `${where} (i): ${id} comes back with its rows`);
      if (left.clears === clears) return;
      clearedReturns += 1;
      // A switch is no History step: the walk shows the drawing to read its manifest and shows the other one again.
      const shownNow = state.activeDrawingId.value;
      app.selectDrawing(id);
      assert.deepEqual(state.proteinIdentityManifest.value, emptyProteinIdentityManifest(), `${where} (i): ${id} comes back with an empty manifest`);
      app.selectDrawing(shownNow);
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
        const drawn = pick(['generate', 'generate', 'new', 'duplicate', 'delete', 'switch', 'undo', 'redo', 'clear', 'orders']);
        // A walk with more drawings each step finds nothing new and makes History's budget work grow with them.
        const kind = ids.length >= MAX_DRAWINGS && ['new', 'duplicate', 'orders'].includes(drawn) ? 'delete' : drawn;
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
          else if ((kind === 'undo' || kind === 'redo') && timeline[kind].length > 0) await traverse(kind);
        }, { undoOrRedo: kind === 'undo' || kind === 'redo' });
        if (kind === 'clear') {
          const named = new Set(state.losatCacheInfo.value.map(({ key }) => key));
          const kept = [...state.losatCache.value].filter(([key, entry]) => !named.has(key) && entry.identityKind !== 'protein');
          assert.deepEqual([app.losatResultRows.value.length, kept.map(([key]) => key)], [0, []], `${at} (g): Clear left results or unnamed searches`);
        }
        if (kind === 'orders') {
          // Review 13's P1 and P2: a new drawing's Generate, then the older steps' Undo and Redo from the first drawing.
          const first = state.activeDrawingId.value;
          const made = pick(['new', 'duplicate']) === 'new'
            ? app.drawingList.newDrawing(pick(['circular', 'linear'])) : app.drawingList.duplicateDrawing(first);
          const moves = [() => generate(1 + Math.floor(next() * 2)), () => app.selectDrawing(first),
            () => traverse('undo'), () => traverse('undo'), () => traverse('redo'), () => traverse('redo')];
          for (const [index, move] of (made.status === 'ok' ? moves : []).entries()) {
            await tracked(`${at} order ${index}`, move, { undoOrRedo: index >= 2 });
            check(`${at} order ${index}`);
          }
        }
        trail.push(`${kind}@${state.activeDrawingId.value}`);
        const where = `seed ${seed} step ${step} [${trail.slice(-6).join(' ')}]`;
        check(where);
        if (recorded) {
          const before = rowsByDrawing();
          await tracked(`${where} undo`, () => traverse('undo'), { undoOrRedo: true });
          check(`${where} undo`);
          await tracked(`${where} redo`, () => traverse('redo'), { undoOrRedo: true });
          assert.equal(rowsByDrawing(), before, `${where} (d): Undo then Redo gives back the same rows`);
          check(`${where} redo`);
        }
      }
    }
    assert.ok(returns > 100, `the walk brought drawings back ${returns} times`);
    assert.ok(clearedReturns > 10, `the walk brought drawings back after a Clear ${clearedReturns} times`);
  } finally {
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
  }
  console.log(`LOSAT rows walk: ${SEEDS.length} seeds x ${STEPS} steps, ${returns} returns (${clearedReturns} after a Clear), in ${Date.now() - started} ms`);
});
