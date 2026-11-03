// Step 3c A2 (CW-02, CW-06 row M2): the record and feature names of the
// Similarity alignment drawer read the catalog's name index, which the
// admission owner (services/feature-catalog.js) builds once per admitted
// catalog. A Generate that reuses the catalog (3a) changes the committed
// request only, so the names read no catalog row again. Counts, not timings.
import assert from 'node:assert/strict';
import test from 'node:test';

const ref = (value) => ({ value });
// Vue's deep `ref` hands out reactive proxies; `toRaw` returns the object behind one.
const raws = new WeakMap();
const reactive = (target) => {
  const proxy = new Proxy(target, {
    get: (object, key, receiver) => {
      const value = Reflect.get(object, key, receiver);
      return value && typeof value === 'object' ? reactive(value) : value;
    }
  });
  raws.set(proxy, target);
  return proxy;
};
globalThis.window = {
  Vue: {
    ref, computed: (getter) => ({ get value() { return getter(); } }), watch: () => {},
    toRaw: (value) => raws.get(value) ?? value
  }
};
const { createSimilarityAlignmentActions } = await import('../../gbdraw/web/js/app/similarity-alignment.js');
const { FEATURE_CATALOG_SCHEMA, admitFeatureCatalog, catalogNameIndex } = await import(
  '../../gbdraw/web/js/services/feature-catalog.js'
);

const row = (recordKey, id, index, recordIdx) => ({
  recordKey, biologicalFeatureId: id, stableFeatureId: `stable-${id}`, record_idx: recordIdx, sourceFeatureIndex: index,
  record_id: `ACC-${recordKey}`, type: 'CDS', start: 10 + index * 10, end: 30 + index * 10, strand: 1,
  anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' },
  gene: `gene-${id}`, qualifiers: { gene: [`gene-${id}`] }
});
// The drawer's members carry no gene: their names come from the catalog rows.
const member = (recordKey, id, index) => ({
  recordKey, biologicalFeatureId: id, sourceFeatureIndex: index, stableFeatureSvgId: `stable-${id}`,
  recordIndex: recordKey === 'a' ? 0 : 1, featureIndex: index, start: 10 + index * 10, end: 30 + index * 10,
  strand: '+', role: 'member', representative: false
});
const members = [member('a', 'x', 1), member('b', 'y', 2)];
const results = [{ name: 'diagram.svg' }];
// Reads of each item's biological features: one per item is one catalog traversal.
const itemReads = { count: 0 };
const newCatalog = () => {
  const catalog = { schema: FEATURE_CATALOG_SCHEMA, items: [{
    resultIndex: 0, resultName: 'diagram.svg', recordKeys: ['a', 'b'],
    features: [], biologicalFeatures: [row('a', 'x', 1, 0), row('a', 'z', 3, 0), row('b', 'y', 2, 1)],
    orthogroups: [], annotations: [], comparisonMatches: [], sequenceSources: []
  }] };
  catalog.items.forEach((item) => {
    const rows = item.biologicalFeatures;
    Object.defineProperty(item, 'biologicalFeatures', {
      enumerable: true, get: () => { itemReads.count += 1; return rows; }
    });
  });
  return catalog;
};
const admitted = (catalog, mode = 'linear') => {
  admitFeatureCatalog(catalog, results, { adopt: true, mode });
  return catalog;
};
// Item reads of `run`, without the admission's own traversal.
const readsOf = (run) => {
  itemReads.count = 0;
  const value = run();
  return { reads: itemReads.count, value };
};

const fixture = (catalog) => {
  const group = { id: 'og-1', members, orthologEdges: [] };
  let request = { records: ['a', 'b'].map((recordKey) => ({ recordKey, presentation: { reverseComplement: false } })) };
  const linearSeqs = ['a', 'b'].map((uid) => ({ uid, region_reverse: false }));
  const state = {
    activeDrawing: () => ({ mode: 'linear', recordsOff: [], linearSeqs }),
    mode: ref('linear'), featureCatalog: ref(catalog), selectedOrthogroupId: ref('og-1'),
    similarityAlignmentPlan: ref(null), similarityAlignmentResetReceipt: ref(null), errorLog: ref(null),
    linearRecordTranslations: ref([]),
    results: ref(results)
  };
  const actions = createSimilarityAlignmentActions({
    state,
    getOrthogroupById: (id) => (id === group.id ? group : null),
    getEnrichedOrthogroupMembers: () => group.members,
    selectedOrthogroupDetail: { get value() { return { group, id: group.id, members: group.members }; } },
    getCommittedRequest: () => request,
    getCommittedSession: () => ({ renderRequest: request }),
    projectCommittedAlignment: () => assert.fail('no alignment runs'),
    runRecordAlignment: () => assert.fail('no alignment runs'),
    runHelperOperation: () => assert.fail('no helper runs'),
    resolveOperation: 'resolveSimilarityAlignment'
  });
  const labels = () => actions.drawerReferenceOptions('og-1').map(({ label }) => label);
  // A Generate that reuses the catalog commits a new request (new record rows).
  const commitRequest = () => { request = { records: request.records.map((record) => ({ ...record })) }; };
  return { state, labels, commitRequest };
};
const expected = ['ACC-a · gene-x', 'ACC-b · gene-y'];

test('a reused catalog gives the drawer names without reading its rows again', () => {
  const { state, labels, commitRequest } = fixture(admitted(newCatalog()));
  assert.deepEqual(readsOf(labels), { reads: 1, value: expected }, 'the first read builds the index: one traversal');
  commitRequest();
  assert.deepEqual(readsOf(labels), { reads: 0, value: expected }, 'a new committed request reads no catalog row');
  assert.deepEqual(readsOf(labels), { reads: 0, value: expected });

  // A new catalog object is indexed once, also when its rows are equal (identity, not content).
  state.featureCatalog.value = admitted(newCatalog());
  assert.deepEqual(readsOf(labels), { reads: 1, value: expected });
  assert.deepEqual(readsOf(labels), { reads: 0, value: expected });
  const copy = admitted(structuredClone(state.featureCatalog.value));
  state.featureCatalog.value = copy;
  assert.deepEqual(readsOf(labels).value, expected, 'an equal copy is indexed again');
  assert.equal(catalogNameIndex(copy).catalog, copy);
});

test('a catalog held in a Vue ref is indexed by its raw object', () => {
  const catalog = admitted(newCatalog());
  const { labels } = fixture(reactive(catalog));
  assert.deepEqual(readsOf(labels), { reads: 1, value: expected });
  assert.deepEqual(readsOf(labels), { reads: 0, value: expected }, 'the proxy finds the raw catalog\'s index');
  // The rows are the raw catalog rows: the Reset helper request spreads them as before.
  const index = catalogNameIndex(reactive(catalog));
  assert.equal(index.catalog, catalog);
  const rawRows = catalog.items[0].biologicalFeatures;
  assert.equal(index.rows.length, 3);
  index.rows.forEach((indexed, position) => assert.equal(indexed, rawRows[position]));
});

test('a catalog without an admission is indexed for each read and the miss is counted', () => {
  const metrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
  try {
    const { labels } = fixture(newCatalog());
    const misses = () => metrics.filter(({ name }) => name === 'featureCatalogNameIndexMissCount')
      .reduce((total, { value }) => total + value, 0);
    // The drawer names look the catalog up twice (record and feature names): two traversals per read.
    assert.deepEqual(readsOf(labels), { reads: 2, value: expected }, 'the same names, uncached');
    assert.deepEqual(readsOf(labels), { reads: 2, value: expected });
    assert.equal(misses(), 4, 'each uncached lookup is counted');
    assert.equal(catalogNameIndex(null).rows.length, 0);
    assert.equal(misses(), 4, 'no catalog is not a miss');
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
});

test('a catalog admitted again in another mode is indexed again, once', () => {
  const catalog = admitted(newCatalog());
  const { labels } = fixture(catalog);
  assert.deepEqual(readsOf(labels), { reads: 1, value: expected });
  admitted(catalog, 'circular');
  assert.deepEqual(readsOf(labels), { reads: 1, value: expected }, 'the new admission holds a new index');
  assert.deepEqual(readsOf(labels), { reads: 0, value: expected });
});

// A record whose File names its accession does not read the catalog for one:
// its rows are not walked, also when none of them names a usable record ID.
test('a record with its own accession does not walk its catalog rows for one', () => {
  const recordIdReads = { count: 0 };
  const catalog = newCatalog();
  const rows = catalog.items[0].biologicalFeatures;
  rows.filter(({ recordKey }) => recordKey === 'a').forEach((feature) => {
    Object.defineProperty(feature, 'record_id', { enumerable: true, get: () => { recordIdReads.count += 1; return ''; } });
  });
  const { state, labels } = fixture(admitted(catalog));
  state.activeDrawing().linearSeqs[0].accession = 'SEQ-A';
  recordIdReads.count = 0;
  assert.deepEqual(labels(), ['SEQ-A · gene-x', 'ACC-b · gene-y']);
  assert.deepEqual(labels(), ['SEQ-A · gene-x', 'ACC-b · gene-y']);
  assert.equal(recordIdReads.count, 0, 'the rows of record a are not read for an accession');
});

test('a member without facts and with a null start names no coordinates', () => {
  const bare = { ...member('a', 'q', 5), start: null, end: 20, type: 'CDS' };
  const { state, labels } = fixture(admitted(newCatalog()));
  const group = { id: 'og-1', members: [bare] };
  const actions = createSimilarityAlignmentActions({
    state,
    getOrthogroupById: () => group,
    getEnrichedOrthogroupMembers: () => group.members,
    selectedOrthogroupDetail: { get value() { return { group, id: group.id, members: group.members }; } },
    getCommittedRequest: () => ({ records: [] }),
    getCommittedSession: () => ({ renderRequest: { records: [] } }),
    projectCommittedAlignment: () => assert.fail('no alignment runs'),
    runRecordAlignment: () => assert.fail('no alignment runs'),
    runHelperOperation: () => assert.fail('no helper runs'),
    resolveOperation: 'resolveSimilarityAlignment'
  });
  assert.ok(labels().length > 0);
  const [{ label }] = actions.drawerReferenceOptions('og-1');
  assert.ok(label.endsWith(' · CDS'), label);
  assert.ok(!label.includes('1..20'), label);
});

test('inparalog reference labels with a missing bound show only the present bound', () => {
  const bare = (id, index, bounds) => ({ ...member('a', id, index), type: 'CDS', ...bounds });
  const { state } = fixture(admitted(newCatalog()));
  const group = { id: 'og-1', members: [bare('q', 5, { start: null, end: 20 }), bare('r', 6, { start: 10, end: null })] };
  const actions = createSimilarityAlignmentActions({
    state,
    getOrthogroupById: () => group,
    getEnrichedOrthogroupMembers: () => group.members,
    selectedOrthogroupDetail: { get value() { return { group, id: group.id, members: group.members }; } },
    getCommittedRequest: () => ({ records: [] }),
    getCommittedSession: () => ({ renderRequest: { records: [] } }),
    projectCommittedAlignment: () => assert.fail('no alignment runs'),
    runRecordAlignment: () => assert.fail('no alignment runs'),
    runHelperOperation: () => assert.fail('no helper runs'),
    resolveOperation: 'resolveSimilarityAlignment'
  });
  const labels = actions.drawerReferenceOptions('og-1').map(({ label }) => label);
  assert.equal(labels.length, 2);
  assert.ok(labels[0].endsWith(' · ..20 (+)'), labels[0]);
  assert.ok(labels[1].endsWith(' · 11.. (+)'), labels[1]);
  labels.forEach((label) => assert.ok(!label.includes('NaN'), label));
});
