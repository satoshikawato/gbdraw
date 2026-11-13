// Guard: each drawing's LOSAT evidence through every drawing-list step and History (#1777 rounds 14 to 17, after reviews
// 8 to 13 each found a defect in the rows; D-221: the evidence is the drawing's artifact values, so it moves with the
// drawing and History by reference and no operation counts holders). A deterministic walk (fixed seeds, bounded steps)
// over a Generate (one History step that writes the shown drawing's evidence: Circular ring rows; Linear window-pair
// rows with their protein searches, one all-against-all search without a row, the converted payload, and the drawing's
// protein identity manifest, whose record instance keys every Linear Generate shares with other meanings), New,
// Duplicate, Delete, a drawing switch, Clear LOSAT results (one History step, C5 (a)), Undo, Redo, and review 13's
// orders (New or Duplicate, a Generate in the new drawing, a switch back, Undo, Undo, Redo, Redo), in the real app
// (`createAppSetup`, its drawing list and History). After every step it checks that:
//   (a') every drawing's rows (shown or kept) have their keys in its own raw searches;
//   (b) a drawing holds one row per cache key;
//   (d) Undo then Redo of a step the walk recorded (a Generate, Clear, New, Duplicate or Delete) gives back every
//       drawing's values, the same objects;
//   (f) the shown Linear drawing's per-pair readers list each of its edges once (OV-3210);
//   (g) Clear LOSAT results counts the shown drawing's rows, one per key (D-223); after a Clear the shown drawing holds
//       fresh empty values and every other drawing's values are the same objects;
//   (h) an Undo or a Redo leaves the values of every drawing but the one whose step it takes the same objects;
//   (i) a drawing that comes back (an Undo of its Delete, a Redo of its New or Duplicate) has the values it left with,
//       the same objects, whatever Clear or Generate ran meanwhile (OV-3474);
//   (j) the shown drawing's protein rows resolve in its own manifest (OV-3213, OV-3214), and the all-against-all search
//       of the Generate they come from is in its own raw searches;
//   (k) the hint's "another drawing has results" is whether a drawing of the project lists one: a drawing History holds
//       out of the project (a Delete, an undone New or Duplicate) is none (D-224: its searches are not reused).
// The Generate lookup order (own, the project drawings' current searches, the retry) is pinned in
// `run-analysis-simple-path.test.mjs`; Save and Load per drawing in `mode-scoped-session.test.mjs`.
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
const { LOSAT_EVIDENCE_KEYS, emptyLosatEvidence, validateProteinRawEntryReferences } =
  await import('../../gbdraw/web/js/services/losat-cache.js');
const { createLinearComparisonView } = await import('../../gbdraw/web/js/app/linear-comparison-panel.js');
const { LEGACY_PROTEIN_CANDIDATE_SCHEMA } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
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
/** A converted payload of a protein Generate. @param {string} key */
const derivedEvidence = (key) => ({
  schema: LOSAT_DERIVED_CACHE_SCHEMA, kind: 'derived-losatp-payload', idEncoding: 'runtime-handle-v1', key, mode: 'pairwise', payload: {}
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
 * A Generate of the shown drawing as one History step, as Generate commits (D-221): its rows replace its own, its
 * searches join its own (new Maps; it keeps the searches it did not use until Clear), and a protein Generate writes its
 * manifest and converted payloads.
 * @param {import('../../gbdraw/web/js/services/losat-cache.js').LosatCacheInfoEntry[]} rows
 * @param {{ searches?: ReadonlyMap<string, object>, manifest?: object, derived?: readonly { key: string }[], legacy?: object }} [evidence]
 */
const generateRows = (rows, { searches = new Map(), manifest, derived = [], legacy } = {}) => window.__GBDRAW_HISTORY__.runUndoableArtifactReplacement('Generate diagram', () => {
  state.losatCache.value = new Map([...state.losatCache.value, ...rows.map((row) => [row.key, { key: row.key }]), ...searches]);
  state.losatCacheInfo.value = rows.map((row) => ({ ...row }));
  if (manifest) state.proteinIdentityManifest.value = manifest;
  if (derived.length) state.losatDerivedCache.value = new Map([...state.losatDerivedCache.value, ...derived.map((entry) => [entry.key, entry])]);
  if (legacy) state.legacyProteinRawCandidates.value = legacy;
  return { status: 'ok' };
}, { shouldCommit: (outcome) => outcome?.status === 'ok' });
/** @param {string} key */
const ringRow = (key) => ({ key, filename: `${key}.tsv`, display: true });
/** @param {string} key */
const pairRow = (key) => ({ key, filename: `${key}.tsv`, display: true, edgeKey: `edge-${key}`, queryUid: 'u1', subjectUid: 'u2' });
/** The shown drawing's LOSAT evidence: its six values. */
const shownValues = () => Object.fromEntries(LOSAT_EVIDENCE_KEYS.map((key) => [key, state[key].value]));
/**
 * Every drawing's LOSAT evidence, by drawing ID: the shown drawing's from the state, another's from its kept slot (the
 * empty values for a drawing that has no slot yet). Reads only.
 * @param {ReturnType<typeof createAppSetup>} app
 */
const projectValues = (app) => {
  const kept = new Map(app.sessionImportRollbackState.capture().artifactSlots.map((slot) => [slot.drawingId, slot.values]));
  return Object.fromEntries(state.drawingTabs.value.map(({ id }) => [id, id === state.activeDrawingId.value
    ? shownValues()
    : Object.fromEntries(LOSAT_EVIDENCE_KEYS.map((key) => [key, kept.get(id)?.[key] ?? null]))]));
};
/** @param {ReturnType<typeof createAppSetup>} app @param {string} drawingId */
const valuesOf = (app, drawingId) => projectValues(app)[drawingId];
/** @param {ReturnType<typeof createAppSetup>} app @param {string} drawingId */
const keysOf = (app, drawingId) => (valuesOf(app, drawingId)?.losatCacheInfo ?? []).map(({ key }) => key);
/** @param {ReturnType<typeof createAppSetup>} app @param {string} drawingId */
const searchKeysOf = (app, drawingId) => [...(valuesOf(app, drawingId)?.losatCache ?? new Map()).keys()].sort();
/** Whether two drawings' values are the same objects, key by key. @param {object} left @param {object} right */
const sameValues = (left, right) => LOSAT_EVIDENCE_KEYS.every((key) => left?.[key] === right?.[key]);
const resetLosat = () => {
  Object.entries(emptyLosatEvidence()).forEach(([key, value]) => { state[key].value = value; });
};
/** Every drawing but the shown one leaves (a case starts from one drawing). @param {ReturnType<typeof createAppSetup>} app */
const onlyShownDrawing = (app) => state.drawingTabs.value.map((tab) => tab.id)
  .filter((id) => id !== state.activeDrawingId.value).forEach((id) => app.drawingList.deleteDrawing(id));

// OV-3207 (review 9's order): History's Undo of one drawing's Generate leaves another drawing's results; D-221: it puts
// back the drawing's own evidence as it was before (its search goes; a later Generate searches again, D-224).
test('Undo of one drawing\'s Generate keeps another drawing\'s LOSAT results; Redo gives its own back (OV-3207)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    onlyShownDrawing(app);
    const a = state.activeDrawingId.value;
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(a); await generateRows([ringRow('k1')]);
    const generated = valuesOf(app, a);
    app.selectDrawing(b); await generateRows([ringRow('k2')]);
    const bValues = valuesOf(app, b);
    app.selectDrawing(a); await history.undo();
    assert.deepEqual([keysOf(app, a), searchKeysOf(app, a), keysOf(app, b)], [[], [], ['k2']], 'B keeps its result');
    assert.ok(sameValues(valuesOf(app, b), bValues), 'B\'s values are the same objects');
    await history.redo();
    assert.ok(sameValues(valuesOf(app, a), generated), 'Redo gives A its values back');
  } finally { resetLosat(); }
});

// OV-3212, OV-3474: a deleted drawing holds its values in its slot; another drawing's Generate, its Undo and Clear do
// not reach it, in any order, and the Undo of the Delete brings it back with the same objects.
test('a deleted drawing\'s LOSAT values survive another drawing\'s Generate, its Undo and Clear (OV-3212, OV-3474)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  for (const deleteUndoneFirst of [true, false]) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const a = state.activeDrawingId.value;
      const { drawingId: b } = app.drawingList.newDrawing('circular');
      const { drawingId: c } = app.drawingList.newDrawing('circular');
      app.selectDrawing(b); await generateRows([ringRow('k2')]);
      const bValues = valuesOf(app, b);
      app.selectDrawing(a);
      app.drawingList.deleteDrawing(b);
      app.selectDrawing(a); await generateRows([ringRow('k1')]);
      await app.clearLosatCache();
      if (deleteUndoneFirst) {
        app.selectDrawing(c); await history.undo();
        app.selectDrawing(a); await history.undo(); await history.undo();
      } else {
        await history.undo(); await history.undo();
        await history.undo();
      }
      assert.ok(state.drawingTabs.value.some((tab) => tab.id === b), 'B is back');
      assert.ok(sameValues(valuesOf(app, b), bValues), `B's values are the same objects (Delete undone first: ${deleteUndoneFirst})`);
      assert.deepEqual(keysOf(app, a), [], 'A is back before its Generate');
      app.drawingList.deleteDrawing(c);
    } finally { resetLosat(); }
  }
});

// OV-3212 (review 13's P1 and P2), OV-3474: a drawing whose New or Duplicate is undone holds its values; the Undo and
// Redo of another drawing's older Generate in between leave them, so the Redo brings the drawing back with them. (A new
// step, such as a Clear, after the Undo ends the Redo: the drawing is dropped with its values.)
test('a drawing held by an undone New or Duplicate comes back with its values after another drawing\'s Undo (OV-3212, OV-3474)', async () => {
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
      const bValues = valuesOf(app, b);
      app.selectDrawing(a);
      await history.undo();
      assert.ok(!state.drawingTabs.value.some((tab) => tab.id === b), `${make}: B is held`);
      await history.undo();
      assert.deepEqual([keysOf(app, a), searchKeysOf(app, a)], [[], []], `${make}: A before its Generate`);
      await history.redo();
      await history.redo();
      assert.ok(sameValues(valuesOf(app, b), bValues), `${make}: B comes back with its values`);
    } finally { resetLosat(); }
  }
});

// PD-OI-022 scenario revision 2 (Owner 2026-11-13, Q2=B) and C5 (a) (Owner 2026-11-14), replacing OV-3215: Clear LOSAT
// results clears the shown drawing's results as one History step of it. Its Undo brings back exactly what it cleared,
// the same objects; its Redo clears again; another drawing's values stay.
test('Clear LOSAT results is one History step: Undo brings back what it cleared, Redo clears again (C5 (a))', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  for (const mode of ['circular', 'linear']) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const { drawingId: a } = app.drawingList.newDrawing(mode);
      onlyShownDrawing(app);
      const { drawingId: b } = app.drawingList.newDrawing(mode);
      const linear = mode === 'linear';
      const row = linear ? pairRow : ringRow;
      app.selectDrawing(b); await generateRows([row('kB')]);
      const bValues = valuesOf(app, b);
      app.selectDrawing(a);
      await generateRows([row('k0')]);
      await generateRows([row('kS'), row('k1')], linear ? {
        searches: new Map([['k1', proteinSearch('a', 'u1', 'u2')], ['allA', proteinSearch('a', 'u2', 'u1')]]),
        manifest: proteinManifest('a', ['u1', 'u2']), derived: [derivedEvidence('pA')],
        legacy: { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [{ state: 'pending' }] }
      } : {});
      const cleared = shownValues();
      assert.equal(app.losatResultRows.value.length, 2, `${mode}: "Clear 2 LOSAT results"`);
      const undoCount = history.getUndoCount();
      await app.clearLosatCache();
      assert.equal(history.getUndoCount(), undoCount + 1, `${mode}: one History step`);
      const empty = emptyLosatEvidence();
      assert.deepEqual(Object.fromEntries(LOSAT_EVIDENCE_KEYS.map((key) => [key, state[key].value])), empty, `${mode}: A's values are empty`);
      assert.equal(app.losatResultRows.value.length, 0);
      assert.ok(sameValues(valuesOf(app, b), bValues), `${mode}: B's values stay`);
      await history.undo();
      assert.ok(sameValues(shownValues(), cleared), `${mode}: Undo brings back exactly what Clear cleared`);
      assert.equal(app.losatResultRows.value.length, 2);
      await history.redo();
      assert.deepEqual(shownValues(), empty, `${mode}: Redo clears again`);
      assert.ok(sameValues(valuesOf(app, b), bValues), `${mode}: B's values stay through Undo and Redo`);
      await history.undo();
      await history.undo();
      assert.deepEqual(keysOf(app, a), ['k0'], `${mode}: the Undo before it takes back the Generate`);
    } finally { resetLosat(); }
  }
});

// OV-3474 (rewrites OV-3274): a drawing History holds out of the project keeps its values through a Clear in another
// drawing: it comes back with its rows and manifest, the same objects.
test('a drawing that comes back after a Clear in another drawing has its rows and manifest (OV-3474)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  for (const leave of ['delete', 'undo-new', 'undo-duplicate']) {
    resetLosat();
    try {
      onlyShownDrawing(app);
      const { drawingId: a } = app.drawingList.newDrawing('linear');
      onlyShownDrawing(app);
      await generateRows([pairRow('kA')]);
      const made = leave === 'undo-duplicate' ? app.drawingList.duplicateDrawing(a) : app.drawingList.newDrawing('linear');
      const b = made.drawingId;
      app.selectDrawing(b);
      await generateRows([pairRow('kB')], { searches: new Map([['kB', proteinSearch('b', 'u1', 'u2')], ['allB', proteinSearch('b', 'u2', 'u1')]]),
        manifest: proteinManifest('b', ['u1', 'u2']) });
      const bValues = valuesOf(app, b);
      app.selectDrawing(a);
      if (leave === 'delete') app.drawingList.deleteDrawing(b);
      else { await history.undo(); await history.undo(); }
      assert.ok(!state.drawingTabs.value.some((tab) => tab.id === b), `${leave}: B is held`);
      await app.clearLosatCache();
      if (leave === 'delete') { await history.undo(); await history.undo(); } else { await history.redo(); await history.redo(); }
      assert.ok(state.drawingTabs.value.some((tab) => tab.id === b), `${leave}: B is back`);
      assert.ok(sameValues(valuesOf(app, b), bValues), `${leave}: B comes back with its rows and manifest`);
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
  /** @param {string} drawingId @param {string} key */
  const resolves = (drawingId, key) => {
    const values = valuesOf(app, drawingId);
    return validateProteinRawEntryReferences(values.losatCache.get(key), values.proteinIdentityManifest);
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
    assert.deepEqual([keysOf(app, d), resolves(d, 'kA'), resolves(a, 'kA')], [['kA'], true, true], 'after its Undo');
    app.selectDrawing(d); await history.redo();
    assert.deepEqual([resolves(d, 'kD'), resolves(a, 'kA')], [true, true], 'after its Redo');
  } finally { resetLosat(); }
});

// OV-3211, D-223 (rewritten for D-221): the count is the shown drawing's results; the hint says whether another drawing
// of the project lists some. A drawing a Delete holds is none of them; its Undo brings it back with its results.
test('Clear counts the shown drawing\'s results; the hint reads the project\'s drawings only (OV-3211)', async () => {
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
    await app.clearLosatCache();
    assert.deepEqual([searchKeysOf(app, a), app.losatResultRows.value.length, app.projectHasLosatResults.value], [[], 0, false],
      'no drawing of the project lists a result: "No LOSAT results yet"');
    await history.undo();
    await history.undo();
    assert.deepEqual([state.drawingTabs.value.some((tab) => tab.id === b), keysOf(app, b), app.projectHasLosatResults.value],
      [true, ['k2'], true], 'B comes back with its result');
  } finally { resetLosat(); }
});

// OV-3471 (F3), D-224: the older Session's protein candidates Load gave a drawing are its own: an Undo in another drawing
// leaves them as their drawing left them.
test('an Undo in one drawing leaves another drawing\'s legacy protein candidates (OV-3471)', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  resetLosat();
  try {
    onlyShownDrawing(app);
    const a = state.activeDrawingId.value;
    const pending = { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [{ state: 'pending' }] };
    const promoted = { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [{ state: 'promoted' }] };
    await generateRows([ringRow('k0')], { legacy: pending });
    const { drawingId: b } = app.drawingList.newDrawing('circular');
    app.selectDrawing(b);
    await generateRows([ringRow('kB')]);
    app.selectDrawing(a);
    await generateRows([ringRow('k1')], { legacy: promoted });
    app.selectDrawing(b);
    await history.undo();
    assert.equal(valuesOf(app, a).legacyProteinRawCandidates, promoted);
  } finally { resetLosat(); }
});

// OV-3470 (F2): a Linear comparison edit (a reindex of the pairs, or an invalidation) rewrites the rows of the drawing
// it names, through the History snapshot service's `editLosatRows`, and leaves every other drawing's.
test('a Linear comparison edit in one drawing leaves another Linear drawing\'s rows (OV-3470)', async () => {
  const app = createAppSetup();
  resetLosat();
  try {
    onlyShownDrawing(app);
    const { drawingId: l1 } = app.drawingList.newDrawing('linear');
    onlyShownDrawing(app);
    const { drawingId: l2 } = app.drawingList.newDrawing('linear');
    app.selectDrawing(l2); await generateRows([pairRow('k2')]);
    app.selectDrawing(l1); await generateRows([pairRow('k1')]);
    const before = valuesOf(app, l2);
    /** @type {Array<[string, number]>} */
    const edits = [];
    const view = createLinearComparisonView({
      state, mode: state.mode, shownDrawingOfMode: () => state.activeDrawing(),
      editLosatRows: (drawingId, edit) => { edits.push([drawingId, edit(state.losatCacheInfo.value).length]); }
    });
    view.invalidateLinearComparisonArtifacts(state.drawingOf(l1), { preserveLosatCacheInfo: true });
    view.invalidateLinearComparisonArtifacts(state.drawingOf(l1));
    assert.deepEqual(edits.map(([id]) => id), [l1, l1], 'each edit names the drawing it rewrites');
    assert.ok(sameValues(valuesOf(app, l2), before), 'L2\'s values stay');
  } finally { resetLosat(); }
});

// A project of one drawing: Clear empties every LOSAT value, a search no row names included.
test('in a project of one drawing, Clear LOSAT results empties every LOSAT value', async () => {
  const app = createAppSetup();
  resetLosat();
  try {
    onlyShownDrawing(app);
    await generateRows([ringRow('k1')], { searches: new Map([['orphan', { key: 'orphan' }]]), derived: [derivedEvidence('p1')] });
    await app.clearLosatCache();
    assert.deepEqual([state.losatCacheInfo.value, [...state.losatCache.value], state.losatDerivedCache.value.size, app.losatResultRows.value.length], [[], [], 0, 0]);
  } finally { resetLosat(); }
});

test('a drawing\'s LOSAT values follow it through Generate, New, Duplicate, Delete, a switch, Clear, Undo and Redo', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const tabIds = () => state.drawingTabs.value.map((tab) => tab.id);
  let search = 0;
  /** The Generate a row's search came from: a Linear Generate's searches and payload are named after it. */
  const generationOf = new Map();
  /** @param {string} where */
  const check = (where) => {
    const project = projectValues(app);
    for (const [id, values] of Object.entries(project)) {
      const rows = values.losatCacheInfo ?? [];
      const unsearched = rows.filter((row) => !values.losatCache?.has(row.key)).map((row) => row.key);
      assert.deepEqual(unsearched, [], `${where} (a'): ${id}'s rows whose search it does not hold`);
      const keys = rows.map(({ key }) => key);
      assert.equal(new Set(keys).size, keys.length, `${where} (b): ${id} holds one row per key: ${keys}`);
    }
    const shown = state.activeDrawing();
    const rows = state.losatCacheInfo.value;
    assert.equal(app.losatResultRows.value.length, new Set(rows.map(({ key }) => key)).size,
      `${where} (g): the count is the shown drawing's results`);
    const unresolved = rows.filter(({ key }) => {
      const entry = state.losatCache.value.get(key);
      return entry?.identityKind === 'protein' && (!validateProteinRawEntryReferences(entry, state.proteinIdentityManifest.value)
        || !state.losatCache.value.has(`all-${generationOf.get(key)}`));
    }).map(({ key }) => key);
    assert.deepEqual(unresolved, [], `${where} (j): protein rows whose evidence does not resolve in the drawing`);
    assert.equal(app.projectHasLosatResults.value, Object.values(project).some((values) => (values.losatCacheInfo ?? []).length > 0),
      `${where} (k): the hint reads the project's drawings`);
    if (shown.mode !== 'linear') return;
    const own = rows.map((row) => row.edgeKey).sort();
    // The units come from the shown Result's catalog, which a switch replaces: the walk's Results all draw CATALOG.
    state.featureCatalog.value = CATALOG;
    const listed = Object.values(app.linearLosatWindowPairs.value).flat().map(({ edgeKey }) => edgeKey).sort();
    assert.deepEqual(listed, own, `${where} (f): the window pairs the shown Linear drawing lists`);
  };
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
        ? { key, filename: `${key}.tsv`, display: true, edgeKey: `edge-${search}`, queryUid: `a#${(index % 2) + 1}`, subjectUid: 'b' }
        : { key, filename: `${key}.tsv`, display: true };
    });
    if (drawing.mode !== 'linear') return generateRows(rows);
    const searches = new Map([...rows.map(({ key }) => [key, proteinSearch(meaning, 'u1', 'u2')]),
      [`all-${meaning}`, proteinSearch(meaning, 'u2', 'u1')]]);
    return generateRows(rows, { searches, manifest: proteinManifest(meaning, ['u1', 'u2']), derived: [derivedEvidence(`payload-${meaning}`)] });
  };
  // (h) and (i): a drawing that leaves the project (a Delete, the Undo of a New or Duplicate) leaves with its values,
  // and the same drawing coming back (by Undo or Redo) brings them back, whatever ran meanwhile.
  /** @type {Map<object, object>} */
  const away = new Map();
  let returns = 0;
  let clearedReturns = 0;
  let clears = 0;
  /** @type {Map<object, number>} The Clear count when each drawing left. */
  const leftAt = new Map();
  const lives = () => new Map(tabIds().map((id) => [state.drawingOf(id), id]));
  /** @param {string} where @param {() => unknown} move @param {{ undoOrRedo?: boolean, clear?: boolean }} [options] */
  const tracked = async (where, move, { undoOrRedo = false, clear = false } = {}) => {
    const shownBefore = state.activeDrawingId.value;
    const before = { values: projectValues(app), lives: lives() };
    await move();
    const after = { values: projectValues(app), lives: lives() };
    before.lives.forEach((id, drawing) => {
      if (after.lives.has(drawing)) return;
      away.set(drawing, before.values[id]);
      leftAt.set(drawing, clears);
    });
    after.lives.forEach((id, drawing) => {
      const left = away.get(drawing);
      if (before.lives.has(drawing) || !left) return;
      returns += 1;
      if (leftAt.get(drawing) !== clears) clearedReturns += 1;
      away.delete(drawing);
      assert.ok(sameValues(after.values[id], left), `${where} (i): ${id} comes back with the values it left with`);
    });
    if (!undoOrRedo && !clear) return;
    // A Clear, an Undo or a Redo changes the values of the drawing whose step it is only.
    const changed = Object.keys(after.values).filter((id) => id !== shownBefore && id in before.values
      && !sameValues(after.values[id], before.values[id]));
    assert.deepEqual(changed, [], `${where} (${clear ? 'g' : 'h'}): another drawing's values changed`);
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
        const drawn = pick(['generate', 'generate', 'new', 'duplicate', 'delete', 'switch', 'undo', 'redo', 'clear', 'orders', 'held']);
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
          else if (kind === 'clear') { recorded = (await app.clearLosatCache()) === true; clears += 1; }
          else if ((kind === 'undo' || kind === 'redo') && timeline[kind].length > 0) await traverse(kind);
        }, { undoOrRedo: kind === 'undo' || kind === 'redo', clear: kind === 'clear' });
        if (kind === 'clear') {
          assert.deepEqual([app.losatResultRows.value.length, shownValues()], [0, emptyLosatEvidence()], `${at} (g): Clear left values`);
        }
        if (kind === 'held' && ids.length > 1) {
          // OV-3474: a Clear while a Delete holds a drawing, then both Undos. A Clear is a History step (C5 (a)), so it
          // ends the Redo that would bring an undone New or Duplicate back; a Delete's Undo still brings its drawing back.
          const held = pick(ids.filter((id) => id !== state.activeDrawingId.value));
          const moves = [() => app.drawingList.deleteDrawing(held), async () => { await generate(1); await app.clearLosatCache(); clears += 1; },
            () => traverse('undo'), () => traverse('undo'), () => traverse('undo')];
          for (const [index, move] of moves.entries()) {
            await tracked(`${at} held ${index}`, move, { undoOrRedo: index >= 2, clear: index === 1 });
            check(`${at} held ${index}`);
          }
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
          const before = projectValues(app);
          await tracked(`${where} undo`, () => traverse('undo'), { undoOrRedo: true });
          check(`${where} undo`);
          await tracked(`${where} redo`, () => traverse('redo'), { undoOrRedo: true });
          const after = projectValues(app);
          const differs = Object.keys(before).filter((id) => !sameValues(after[id], before[id]));
          assert.deepEqual(differs, [], `${where} (d): Undo then Redo gives back the same values`);
          check(`${where} redo`);
        }
      }
    }
    assert.ok(returns > 100, `the walk brought drawings back ${returns} times`);
    assert.ok(clearedReturns > 10, `the walk brought drawings back after a Clear ${clearedReturns} times`);
  } finally {
    resetLosat();
  }
  console.log(`LOSAT values walk: ${SEEDS.length} seeds x ${STEPS} steps, ${returns} returns (${clearedReturns} after a Clear), in ${Date.now() - started} ms`);
});
