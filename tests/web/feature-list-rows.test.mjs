import assert from 'node:assert/strict';
import test from 'node:test';

import { resultCatalogFeatures } from '../../gbdraw/web/js/services/feature-catalog.js';
import {
  featureDrawnContext,
  listFeatureRows,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/app/feature-visibility.js';

// R-5 (Owner decision 2026-10-05): the Features list and Search features list
// the displayed Result's biological features that are of a selected type, have
// their own Feature visibility, or are drawn, with whether each is drawn.
const anchorProfile = { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' };
const biological = (recordKey, id, type, start) => ({
  recordKey, biologicalFeatureId: id, record_id: recordKey, type, start, end: start + 30, strand: 1,
  anchorProfile, qualifiers: { locus_tag: [id] }
});
const rendered = (recordKey, id, svgId) => ({
  svgId, recordKey, biologicalFeatureId: id, fillColor: '#000000',
  drawnSelector: { hash: svgId, location: null, recordLocation: null }
});
const item = (resultIndex, recordKey, features, drawn) => ({
  resultIndex,
  resultName: `result-${resultIndex}.svg`,
  recordKeys: [recordKey],
  biologicalFeatures: features.map(([id, type], index) => biological(recordKey, id, type, index * 100)),
  features: drawn.map((id) => rendered(recordKey, id, `svg-${id}`)),
  orthogroups: [],
  annotations: [],
  comparisonMatches: []
});
// A batch of two Results. The first draws A (CDS) and H (a gene a color rule
// draws); B, K (CDS) and G, J (genes) are not drawn.
const catalog = {
  schema: 5,
  items: [
    item(0, 'REC1', [['A', 'CDS'], ['B', 'CDS'], ['G', 'gene'], ['H', 'gene'], ['J', 'gene'], ['K', 'CDS']], ['A', 'H']),
    item(1, 'REC2', [['Z', 'CDS'], ['Y', 'gene']], ['Z'])
  ]
};
const results = [{ name: 'result-0.svg', content: '<svg />' }, { name: 'result-1.svg', content: '<svg />' }];

const listOf = (resultIndex, featureOverrides = {}, rules = []) => {
  const state = {
    featureCatalog: { value: catalog },
    results: { value: results },
    selectedResultIndex: { value: resultIndex },
    generatedMode: { value: 'circular' },
    featureOverrides,
    featureVisibilityManualRules: rules,
    manualSpecificRules: []
  };
  const catalogFeatures = resultCatalogFeatures(state);
  const list = listFeatureRows(catalogFeatures, featureDrawnContext(state, {
    diagramOptions: { selectedFeaturesSet: ['CDS'] }
  }));
  return {
    catalogFeatures,
    list,
    rows: list.rows.map((row) => [row.biological_feature_id, list.drawn.get(row), list.rendered.has(row)])
  };
};
const featureOf = (catalogFeatures, id) => catalogFeatures.biological.find((feature) => feature.biological_feature_id === id);

test('the list holds the selected types, the drawn features, and features with their own visibility', () => {
  const featureOverrides = {};
  const { catalogFeatures } = listOf(0);
  setFeatureVisibilityOverride(featureOverrides, featureOf(catalogFeatures, 'B'), 'off');
  // J was turned On and then Off: it stays listed so it can be turned On again.
  setFeatureVisibilityOverride(featureOverrides, featureOf(catalogFeatures, 'J'), 'off');
  // A rule hides K; its match is unknown until Python evaluates it, so K reads
  // as Python drew it (not drawn).
  const rules = [{ recordId: '*', featureType: 'CDS', qualifier: 'hash', value: '^x$', action: 'off' }];
  assert.deepEqual(listOf(0, featureOverrides, rules).rows, [
    ['A', true, true],
    ['B', false, false],
    ['H', true, true],
    ['J', false, false],
    ['K', false, false]
  ]);
});

test('a per-feature On lists and draws a feature of an unselected type', () => {
  const featureOverrides = {};
  const { catalogFeatures } = listOf(0);
  setFeatureVisibilityOverride(featureOverrides, featureOf(catalogFeatures, 'G'), 'on');
  const { rows } = listOf(0, featureOverrides);
  assert.deepEqual(rows.find(([id]) => id === 'G'), ['G', true, false]);
});

test('a drawn feature is listed as the rendered feature of the displayed Result', () => {
  const { list, catalogFeatures } = listOf(0);
  const row = list.rows[0];
  assert.equal(row, catalogFeatures.renderedByIdentity.get('REC1\u0000A'));
  assert.equal(row.svg_id, 'svg-A');
});

test('a batch Result lists only its own records', () => {
  assert.deepEqual(listOf(1).rows, [['Z', true, true]]);
});

// R13: the Features list reads the selected feature types from the displayed
// Result's committed metadata, the types of the request that drew it; nothing
// reaches `state` through a function.
test('the Features list reads the selected types the displayed Result was drawn with', async () => {
  globalThis.window ??= { Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
  } };
  const { state } = await import('../../gbdraw/web/js/state.js');
  const { admitFeatureCatalog } = await import('../../gbdraw/web/js/services/feature-catalog.js');
  const {
    admitCurrentSessionResults,
    createCurrentSessionResultSource,
    createEmptySvgMutationPlan,
    getCommittedSvgResultMetadata
  } = await import('../../gbdraw/web/js/services/svg-result-ingestion.js');
  const admit = (selectedFeatureTypes) => admitCurrentSessionResults(
    createCurrentSessionResultSource(results, admitFeatureCatalog(catalog, results, { mode: 'circular' })),
    {
      mutationPlan: createEmptySvgMutationPlan(results.length),
      sanitizer: { sanitize: (value) => value },
      selectedFeatureTypes
    }
  );
  state.featureCatalog.value = catalog;
  state.generatedMode.value = 'circular';
  state.selectedResultIndex.value = 0;
  const listed = () => state.featureList.value.rows.map((row) => [
    row.biological_feature_id, state.featureListState(row).drawn
  ]);

  state.results.value = admit(['CDS']);
  assert.deepEqual(getCommittedSvgResultMetadata(state.results.value[0]).selectedFeatureTypes, ['CDS']);
  // The genes are of an unselected type with no rule or edit of their own.
  assert.deepEqual(listed(), [['A', true], ['B', true], ['K', true]]);
  state.results.value = admit(['CDS', 'gene']);
  assert.deepEqual(listed(), [['A', true], ['B', true], ['G', true], ['H', true], ['J', true], ['K', true]]);
  // A request that names no types: each feature reads as the Result draws it.
  state.results.value = admit(undefined);
  assert.equal(getCommittedSvgResultMetadata(state.results.value[0]).selectedFeatureTypes, null);
  assert.deepEqual(listed(), [['A', true], ['B', false], ['G', false], ['H', true], ['J', false], ['K', false]]);
  assert.equal('committedDiagramOptions' in state, false);
  state.results.value = [];
  state.featureCatalog.value = null;
});
