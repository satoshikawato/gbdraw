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
const { linearWindowsPayload } = await import('../../gbdraw/web/js/services/session-request.js');
const { catalogUnitPlaces, isDrawnWindowPairRow, linearWindowPairsByCardEdge, planLinearUnitSearch } = await import('../../gbdraw/web/js/services/linear-units.js');
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
  const search = await planLinearUnitSearch(cards, resolution, (card) => linearWindowsPayload(card, { arranged: false }),
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
  // Each window pair is its own cache-info row, with its own file name.
  assert.deepEqual(byCard.get(linearComparisonEdgeKey('b', 'c')).map(({ edgeKey, filename }) => [edgeKey, filename]),
    cacheInfo.slice(4).map(({ edgeKey, filename }) => [edgeKey, filename]));
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
    state, mode: { value: 'linear' }, losatCacheInfo, shownDrawingOfMode: () => drawing,
    writeLosatCacheInfo: (next) => { losatCacheInfo.value = next; }
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
