// The window pairs of a Linear card pair (OV-2742): a card drawn in several windows is searched one window pair at a
// time, and Generate keeps one LOSAT cache-info row per window pair, keyed by the pair of drawn units. The Comparison
// card lists each window pair of a card pair, with its own Save Raw LOSAT TSV and file name, read from the drawn
// Result's units (its feature catalog); a record reorder keeps the rows. Python plans the units here, as in Generate.
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}, watch: () => {}
  }
};

const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { linearComparisonEdgeKey } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { linearDrawingUnitPlaces } = await import('../../gbdraw/web/js/services/session-request.js');
const {
  catalogUnitPlaces, isDrawnWindowPairRow, linearWindowPairsByCardEdge, planLinearUnitSearch, searchUnitPlaces,
  windowLosatFilenamesAfterRemoval
} = await import('../../gbdraw/web/js/services/linear-units.js');
const { windowPairLosatFilename } = await import('../../gbdraw/web/js/services/losat-cache.js');
const { createLinearComparisonView, createLinearLayoutReaders } = await import('../../gbdraw/web/js/app/linear-comparison-panel.js');
const { pythonLinearUnitPlan } = await import('./helpers/linear-unit-plan.mjs');
const { createSliceHandle } = await import('../../gbdraw/web/js/services/drawing/slice-handle.js');

const span = (start, end) => ({ start, end, reverseComplement: false, column: null });

/** The comparison plan over `cards`: one LOSAT edge per record pair. */
const resolutionOf = (cards, pairs) => {
  const indexOf = new Map(cards.map((card, index) => [card.uid, index]));
  const edges = pairs.map(([queryUid, subjectUid], ordinal) => Object.freeze({
    id: `edge-${ordinal}`, queryUid, subjectUid, queryIndex: indexOf.get(queryUid), subjectIndex: indexOf.get(subjectUid),
    edgeKey: linearComparisonEdgeKey(queryUid, subjectUid), ordinal, source: 'losat', included: true,
    file: null, fileActive: false, losatFilename: '', losatFilenameActive: false
  }));
  return Object.freeze({
    mode: 'selected', defaultSource: 'losat', edges: Object.freeze(edges), errors: Object.freeze([]), error: '',
    valid: true, hasComparisonIntent: edges.length > 0, hasLosatIntent: edges.length > 0, hasUploadIntent: false
  });
};

/**
 * What one Generate leaves: Python's unit plan, the cache-info rows it writes (one per unit pair, as
 * app/run-analysis.js writes them), and the drawn Result's catalog item, whose units name their records as Python's
 * catalog does (`recordKeys` the unit keys in drawing order, `windowOf` each window's record).
 */
const generated = async (cards, pairs) => {
  const resolution = resolutionOf(cards, pairs);
  const search = await planLinearUnitSearch(cards, resolution, linearDrawingUnitPlaces(
    { linearRecordLayoutEnabled: { value: false }, linearRecordRows: [], linearSeqs: cards }, () => ({ length: null, id: '', count: 1 })),
  async (payload) => ({ result: pythonLinearUnitPlan(payload) }));
  const unitKeys = search.sequences.map((sequence) => sequence.uid);
  const windowOf = Object.fromEntries(unitKeys.flatMap((unitKey, unit) => (
    search.windowOf[unit] === null ? [] : [[unitKey, search.windowOf[unit]]])));
  const cacheInfo = search.resolution.edges.map((edge) => ({
    key: `cache-${edge.ordinal}`, filename: `pair-${edge.ordinal}.tsv`, display: true, edgeKey: edge.edgeKey,
    ordinal: edge.ordinal, queryUid: edge.queryUid, subjectUid: edge.subjectUid, queryIndex: edge.queryIndex,
    subjectIndex: edge.subjectIndex
  }));
  const catalog = {
    schema: 1,
    items: [{ resultIndex: 0, resultName: 'r', recordKeys: unitKeys, ...(Object.keys(windowOf).length ? { windowOf } : {}) }]
  };
  return { resolution, search, cacheInfo, catalog };
};

const windowNumbers = (pairs) => pairs.map(({ queryWindow, subjectWindow }) => [queryWindow, subjectWindow]);

test('each window pair of a card pair is listed under the card pair, with its window numbers', async () => {
  const cards = [
    createLinearSeq({ uid: 'a', windows: [span(1, 150), span(151, 300)] }),
    createLinearSeq({ uid: 'b', windows: [span(10, 90), span(200, 280)] }),
    createLinearSeq({ uid: 'c' })
  ];
  const { cacheInfo, catalog } = await generated(cards, [['a', 'b'], ['b', 'c']]);
  const places = catalogUnitPlaces(catalog);
  assert.deepEqual([...places.values()].map(({ recordKey, window }) => [recordKey, window]),
    [['a', 1], ['a', 2], ['b', 1], ['b', 2], ['c', null]]);
  const byCard = linearWindowPairsByCardEdge(cacheInfo, places);
  assert.deepEqual([...byCard.keys()], [linearComparisonEdgeKey('a', 'b'), linearComparisonEdgeKey('b', 'c')]);
  assert.deepEqual(windowNumbers(byCard.get(linearComparisonEdgeKey('a', 'b'))), [[1, 1], [1, 2], [2, 1], [2, 2]]);
  assert.deepEqual(windowNumbers(byCard.get(linearComparisonEdgeKey('b', 'c'))), [[1, null], [2, null]]);
  // Each window pair is its own cache-info row.
  assert.deepEqual(byCard.get(linearComparisonEdgeKey('b', 'c')).map(({ edgeKey }) => edgeKey),
    cacheInfo.slice(4).map(({ edgeKey }) => edgeKey));
});

// OV-2742: Generate places the units it searches as Python planned them, which are the units of the Result it draws,
// and the suffix rule over those places gives each window pair of one card pair a name of its own while a pair of
// cards drawn whole keeps its card pair's name. (Generate's own naming is pinned by the Playwright case.)
test('Generate places its searched units as the Result does, and the suffix rule names each window pair apart', async () => {
  const cards = [
    createLinearSeq({ uid: 'a', windows: [span(1, 150), span(151, 300)] }),
    createLinearSeq({ uid: 'b', windows: [span(10, 90), span(200, 280)] }),
    createLinearSeq({ uid: 'c' }),
    createLinearSeq({ uid: 'd' })
  ];
  const { search, cacheInfo, catalog } = await generated(cards, [['a', 'b'], ['b', 'c'], ['c', 'd']]);
  const places = searchUnitPlaces(search);
  // The units Generate searches are placed as the drawn Result's catalog places them.
  assert.deepEqual([...places], [...catalogUnitPlaces(catalog)]);
  const names = cacheInfo.map(({ queryUid, subjectUid }) => windowPairLosatFilename(
    'A.B.losatn.tsv', places.get(queryUid).window, places.get(subjectUid).window));
  assert.deepEqual(names, [
    'A.B.losatn.w1-w1.tsv', 'A.B.losatn.w1-w2.tsv', 'A.B.losatn.w2-w1.tsv', 'A.B.losatn.w2-w2.tsv',
    'A.B.losatn.w1.tsv', 'A.B.losatn.w2.tsv',
    'A.B.losatn.tsv'
  ]);
  // A name typed on the card pair is the stem of its window pairs' names, not each one's whole name.
  assert.equal(windowPairLosatFilename('mine.TSV', 2, null), 'mine.w2.tsv');
});

// OV-2742: a name typed on a window pair stays with its window when Remove window renumbers the windows after it,
// and goes with its window, or with every window of the record when one is left.
test('Remove window moves each window-pair name with its window and drops the removed window\'s', async () => {
  const cards = [
    createLinearSeq({ uid: 'a', windows: [span(1, 100), span(101, 200), span(201, 300)] }),
    createLinearSeq({ uid: 'b' })
  ];
  const { catalog } = await generated(cards, [['a', 'b']]);
  const places = catalogUnitPlaces(catalog);
  const [w1, w2, w3] = [...places.keys()].slice(0, 3);
  const pair = (unit) => linearComparisonEdgeKey(unit, 'b');
  const names = { [pair(w1)]: 'first.tsv', [pair(w2)]: 'second.tsv', [pair(w3)]: 'third.tsv' };
  // Remove window 1: window 2 is now window 1, window 3 window 2.
  assert.deepEqual(windowLosatFilenamesAfterRemoval(names, places, 'a', 1, 2),
    { [pair(w1)]: 'second.tsv', [pair(w2)]: 'third.tsv' });
  // Remove window 2: window 1 keeps its name.
  assert.deepEqual(windowLosatFilenamesAfterRemoval(names, places, 'a', 2, 2),
    { [pair(w1)]: 'first.tsv', [pair(w2)]: 'third.tsv' });
  // A name of units the Result does not hold cannot be placed (the Web forms no unit key), so it goes; with no drawn
  // Result every name goes, rather than land on another window's span.
  assert.deepEqual(windowLosatFilenamesAfterRemoval({ ...names, 'x->y': 'lost.tsv' }, places, 'a', 3, 2),
    { [pair(w1)]: 'first.tsv', [pair(w2)]: 'second.tsv' });
  assert.deepEqual(windowLosatFilenamesAfterRemoval(names, new Map(), 'a', 3, 2), {});
  // One window left: the record is drawn whole, so no window pair of it is left to name.
  assert.deepEqual(windowLosatFilenamesAfterRemoval({ [pair(w1)]: 'first.tsv' }, places, 'a', 2, 1), {});
  // Another record's windows are not this record's.
  assert.deepEqual(windowLosatFilenamesAfterRemoval(names, places, 'b', 1, 2), names);
});

// OV-2742: only the edited record's side of a name is re-keyed. The other record may be one the last Generate left
// out (its record OFF), so the drawn Result does not place it; its side keeps its own key.
test('Remove window re-keys the edited record\'s side of a name and keeps a side the Result does not place', () => {
  const places = new Map([['R~w1', { recordKey: 'R', window: 1 }], ['R~w2', { recordKey: 'R', window: 2 }],
    ['R~w3', { recordKey: 'R', window: 3 }]]);
  assert.deepEqual(windowLosatFilenamesAfterRemoval({ [linearComparisonEdgeKey('R~w2', 'S')]: 'two.tsv' }, places, 'R', 1, 2),
    { [linearComparisonEdgeKey('R~w1', 'S')]: 'two.tsv' });
  assert.deepEqual(windowLosatFilenamesAfterRemoval({ [linearComparisonEdgeKey('S~w4', 'R~w3')]: 'three.tsv' }, places, 'R', 2, 2),
    { [linearComparisonEdgeKey('S~w4', 'R~w2')]: 'three.tsv' });
  assert.deepEqual(windowLosatFilenamesAfterRemoval({ [linearComparisonEdgeKey('R~w1', 'S')]: 'one.tsv' }, places, 'R', 1, 2), {});
});

test('a card pair whose cards draw one unit each lists no window pair (its row is the card pair\'s own)', async () => {
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 150)] }), createLinearSeq({ uid: 'b' })];
  const { cacheInfo, catalog } = await generated(cards, [['a', 'b']]);
  assert.deepEqual(cacheInfo.map(({ edgeKey }) => edgeKey), [linearComparisonEdgeKey('a', 'b')]);
  assert.equal(linearWindowPairsByCardEdge(cacheInfo, catalogUnitPlaces(catalog)).size, 0);
  assert.equal(catalogUnitPlaces(null).size, 0);
});

/** The panel ports over one drawing whose drawn cards, plan and committed catalog are `cards`, `resolution`, `catalog`. */
const panelOf = (cards, resolution, cacheInfo, catalog) => {
  const drawing = {
    drawnLinearSeqs: { value: cards }, linearComparisonResolution: { value: resolution }, slices: { inputs: createSliceHandle('inputs', { files: {} }) }
  };
  const losatCacheInfo = { value: cacheInfo };
  const state = { activeDrawing: () => drawing, featureCatalog: { value: catalog }, inputFiles: {} };
  const view = createLinearComparisonView({
    state, mode: { value: 'linear' }, shownDrawingOfMode: () => drawing,
    // The drawing's own rows (D-221): the panel edits only the drawing it names.
    editLosatRows: (drawingId, edit) => {
      assert.equal(drawingId, drawing.id);
      losatCacheInfo.value = edit(losatCacheInfo.value);
    }
  });
  const readers = createLinearLayoutReaders({ state, shownDrawingOfMode: () => drawing, losatCacheInfo });
  return { drawing, losatCacheInfo, view, readers };
};

test('a record reorder keeps the window pairs reachable from their card pair', async () => {
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 150), span(151, 300)] }), createLinearSeq({ uid: 'b' })];
  const { resolution, cacheInfo, catalog } = await generated(cards, [['a', 'b']]);
  const { drawing, losatCacheInfo, view, readers } = panelOf(cards, resolution, cacheInfo, catalog);
  const cardEdge = linearComparisonEdgeKey('a', 'b');
  assert.deepEqual(windowNumbers(readers.linearLosatWindowPairs.value[cardEdge]), [[1, null], [2, null]]);
  // Move b above a: the cards change places, the pair keeps its records.
  drawing.drawnLinearSeqs.value = [cards[1], cards[0]];
  drawing.linearComparisonResolution.value = resolutionOf(drawing.drawnLinearSeqs.value, [['a', 'b']]);
  view.invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
  assert.deepEqual(losatCacheInfo.value.map(({ edgeKey }) => edgeKey), cacheInfo.map(({ edgeKey }) => edgeKey));
  assert.deepEqual(windowNumbers(readers.linearLosatWindowPairs.value[cardEdge]), [[1, null], [2, null]]);
  // A record that is no longer drawn takes its window pairs with it.
  drawing.drawnLinearSeqs.value = [cards[0]];
  drawing.linearComparisonResolution.value = resolutionOf(drawing.drawnLinearSeqs.value, []);
  view.invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
  assert.deepEqual(losatCacheInfo.value, []);
});

test('a reorder keeps a card pair\'s own row as before (cards drawn whole)', async () => {
  const cards = [createLinearSeq({ uid: 'a' }), createLinearSeq({ uid: 'b' })];
  const { resolution, cacheInfo, catalog } = await generated(cards, [['a', 'b']]);
  const { drawing, losatCacheInfo, view, readers } = panelOf(cards, resolution, cacheInfo, catalog);
  drawing.drawnLinearSeqs.value = [cards[1], cards[0]];
  drawing.linearComparisonResolution.value = resolutionOf(drawing.drawnLinearSeqs.value, [['a', 'b']]);
  view.invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
  assert.deepEqual(losatCacheInfo.value, [{ ...cacheInfo[0], queryIndex: 1, subjectIndex: 0 }]);
  assert.deepEqual(Object.keys(readers.linearLosatCacheInfoByEdgeKey.value), [linearComparisonEdgeKey('a', 'b')]);
  assert.deepEqual(readers.linearLosatWindowPairs.value, {});
});

// OV-2918: a Session Load restores a window pair's row with its saved identity when the shown Result's catalog names
// its two units and their records are drawn (the rule a reorder applies too); a card pair's own row and a row of a
// record no longer drawn are not window pairs.
test('a window pair\'s row is kept while the shown Result names its units and their records are drawn', async () => {
  const cards = [
    createLinearSeq({ uid: 'a', windows: [span(1, 150), span(151, 300)] }),
    createLinearSeq({ uid: 'b' }),
    createLinearSeq({ uid: 'c' })
  ];
  const { cacheInfo, catalog } = await generated(cards, [['a', 'b'], ['b', 'c']]);
  const places = catalogUnitPlaces(catalog);
  const drawn = (uids) => (uid) => uids.includes(uid);
  assert.deepEqual(cacheInfo.map((entry) => isDrawnWindowPairRow(entry, places, drawn(['a', 'b', 'c']))), [true, true, false]);
  // Without b drawn, the window pairs of a and b are not reachable; without the catalog (no Result), none is.
  assert.deepEqual(cacheInfo.map((entry) => isDrawnWindowPairRow(entry, places, drawn(['a', 'c']))), [false, false, false]);
  assert.deepEqual(cacheInfo.map((entry) => isDrawnWindowPairRow(entry, new Map(), drawn(['a', 'b', 'c']))), [false, false, false]);
});
