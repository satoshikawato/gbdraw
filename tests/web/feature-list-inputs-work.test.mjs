// S4 (perf 0.14.x, counts not timings): the Features list (state.js
// `featureList`, R-5) is derived again only when one of its inputs changes. A
// write of a Result's content (each editor edit commits one, R1) keeps the
// Result names and committed metadata the list reads, so it does not list the
// 12.6k-18.8k catalog features of a Gallery Session again; every input does.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue };
const { watch } = Vue;
const { state } = await import('../../gbdraw/web/js/state.js');
const { installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const { createRulePreparation } = await import('../../gbdraw/web/js/app/rule-matching.js');
const { setFeatureVisibilityOverride } = await import('../../gbdraw/web/js/services/feature-visibility.js');
const { admitFeatureCatalog, resultCatalogFeatures } = await import('../../gbdraw/web/js/services/feature-catalog.js');
const {
  admitCurrentSessionResults, createCurrentSessionResultSource, createEmptySvgMutationPlan
} = await import('../../gbdraw/web/js/services/svg-result-ingestion.js');

const anchorProfile = { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' };
const item = (resultIndex, recordKey, ids, drawn) => ({
  resultIndex,
  resultName: `result-${resultIndex}.svg`,
  recordKeys: [recordKey],
  biologicalFeatures: ids.map((id, index) => ({
    recordKey, biologicalFeatureId: id, record_id: recordKey, type: 'CDS', start: index * 100, end: index * 100 + 30,
    strand: 1, anchorProfile, qualifiers: { locus_tag: [id] }
  })),
  features: drawn.map((id) => ({
    svgId: `svg-${id}`, recordKey, biologicalFeatureId: id, fillColor: '#000000',
    drawnSelector: { hash: `svg-${id}`, location: null, recordLocation: null }
  })),
  orthogroups: [],
  annotations: [],
  comparisonMatches: []
});
const catalog = () => ({ schema: 5, items: [item(0, 'REC1', ['A', 'B', 'C'], ['A', 'B']), item(1, 'REC2', ['Z'], ['Z'])] });
// Committed Results, as Session load and Generate admit them.
const resultsNamed = (...names) => {
  const results = names.map((name) => ({ name, content: '<svg />' }));
  const admission = admitFeatureCatalog(catalog(), results, { mode: 'circular' });
  return admitCurrentSessionResults(createCurrentSessionResultSource(results, admission), {
    mutationPlan: createEmptySvgMutationPlan(results.length), sanitizer: { sanitize: (value) => value }
  });
};

const setup = () => {
  installShownDrawing('circular');
  state.generatedMode.value = 'circular';
  state.featureCatalog.value = catalog();
  state.results.value = resultsNamed('result-0.svg', 'result-1.svg');
  state.selectedResultIndex.value = 0;
  // The features the rule preparation evaluates, as Generate admission sets them.
  const admitted = resultCatalogFeatures(state);
  state.extractedFeatures.value = [...admitted.rendered.values()];
  state.biologicalFeatures.value = admitted.biological;
  const drawing = state.activeDrawing();
  Object.keys(drawing.featureOverrides).forEach((key) => delete drawing.featureOverrides[key]);
  drawing.featureVisibilityManualRules.splice(0);
  let derivations = 0;
  const stop = watch(() => state.featureList.value, () => { derivations += 1; }, { flush: 'sync' });
  const rows = () => state.featureList.value.rows.map((row) => [row.biological_feature_id, state.featureListState(row).drawn]);
  return { drawing, rows, derivations: () => derivations, stop };
};
// As preview-runtime.js `writeResultContent`: the Result is spread through its
// reactive proxy, so its committed state comes back as a reactive copy; the
// committed metadata is frozen, so it reads as the same object.
const writeContent = (index, content) => {
  const next = [...state.results.value];
  next[index] = { ...state.results.value[index], content };
  state.results.value = next;
};

test('a write of the Results\' content does not list the features again', () => {
  const { rows, derivations, stop } = setup();
  assert.deepEqual(rows(), [['A', true], ['B', true], ['C', false]]);
  writeContent(0, '<svg><g/></svg>');
  writeContent(0, '<svg><g/><g/></svg>');
  writeContent(1, '<svg><path/></svg>');
  assert.equal(derivations(), 0);
  stop();
});

test('each input of the list lists the features again', async () => {
  const { drawing, rows, derivations, stop } = setup();
  let seen = 0;
  const changed = (name) => {
    assert.ok(derivations() > seen, `${name} did not list the features again`);
    seen = derivations();
  };
  setFeatureVisibilityOverride(drawing.slices.featureEdits, state.featureList.value.rows[1], 'off');
  changed('a Feature visibility edit');
  assert.deepEqual(rows(), [['A', true], ['B', false], ['C', false]]);

  // A visibility rule: its matches are unknown until the rule preparation records them.
  drawing.featureVisibilityManualRules.push({ recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^A$', action: 'off' });
  changed('a visibility rule');
  const preparation = createRulePreparation({
    state,
    pending: state.ruleMatchingPending,
    evaluate: async ({ features }) => ({ matches: features.map((feature) => (feature.qualifiers.locus_tag[0] === 'A' ? [0] : [])) }),
    visibilityRules: () => drawing.featureVisibilityManualRules
  });
  assert.equal(await preparation.prepareDrawn(), true);
  changed('recorded rule matches');
  assert.deepEqual(rows(), [['A', false], ['B', false], ['C', false]]);

  state.selectedResultIndex.value = 1;
  changed('a Result selection');
  assert.deepEqual(rows().map(([id]) => id), ['Z']);
  state.selectedResultIndex.value = 0;
  changed('a Result selection');

  state.featureCatalog.value = catalog();
  changed('a new catalog');
  state.results.value = [{ name: 'result-0.svg', content: '<svg />' }];
  changed('a Result list of other names');
  stop();
});

// Step 3c A3 (CW-06 row M3): a Generate that reuses the catalog (3a) commits
// new Results whose committed metadata is a new object with the inputs the
// list reads unchanged: the same rendered identities and equal selected
// feature types. The list compares those inputs, not the metadata object.
const { toRaw } = Vue;
const generate = ({ catalogObject = toRaw(state.featureCatalog.value), selectedFeatureTypes = null } = {}) => {
  const results = ['result-0.svg', 'result-1.svg'].map((name) => ({ name, content: '<svg />' }));
  // The same catalog object is an admission hit, as Generate's reuse is.
  const admission = admitFeatureCatalog(catalogObject, results, { mode: 'circular' });
  return admitCurrentSessionResults(createCurrentSessionResultSource(results, admission), {
    mutationPlan: createEmptySvgMutationPlan(results.length), sanitizer: { sanitize: (value) => value }, selectedFeatureTypes
  });
};

test('a Generate that leaves the listed inputs equal does not list the features again', () => {
  const { rows, derivations, stop } = setup();
  let seen = derivations();
  const derived = (name, count) => {
    assert.equal(derivations() - seen, count, `${name}: ${count ? 'lists again' : 'keeps the list'}`);
    seen = derivations();
  };
  state.results.value = generate({ selectedFeatureTypes: ['CDS', 'tRNA'] });
  derived('Results of the shown catalog', 1);
  const listed = rows();
  state.results.value = generate({ selectedFeatureTypes: ['CDS', 'tRNA'] });
  derived('equal feature types and the same identities', 0);
  assert.deepEqual(rows(), listed);
  state.results.value = generate({ selectedFeatureTypes: ['CDS'] });
  derived('other feature types', 1);
  state.results.value = generate({ selectedFeatureTypes: null });
  derived('no feature types after some', 1);
  state.results.value = generate({ selectedFeatureTypes: null });
  derived('no feature types on both sides', 0);
  state.results.value = generate({ selectedFeatureTypes: ['CDS'] });
  derived('feature types after none', 1);
  // The same catalog content admitted as another object has other identities.
  state.results.value = generate({ catalogObject: catalog(), selectedFeatureTypes: ['CDS'] });
  derived('other rendered identities', 1);
  state.featureCatalog.value = catalog();
  derived('a new catalog', 1);
  state.results.value = generate({ selectedFeatureTypes: ['CDS'] });
  derived('Results of the new catalog', 1);
  state.results.value = state.results.value.map((result, index) => (index === 1 ? { ...result, name: 'other.svg' } : result));
  derived('a changed Result name', 1);
  stop();
});

test('without a catalog, a Generate that keeps the rendered identities does not list again', () => {
  const { derivations, stop } = setup();
  const shown = toRaw(state.featureCatalog.value);
  state.results.value = generate({ catalogObject: shown });
  state.featureCatalog.value = null;
  let seen = derivations();
  const listed = () => state.featureList.value.rows.map((row) => row.svg_id);
  assert.deepEqual(listed(), ['svg-A', 'svg-B']);
  state.results.value = generate({ catalogObject: shown });
  assert.equal(derivations() - seen, 0, 'the same identities keep the list');
  state.results.value = generate({ catalogObject: catalog() });
  assert.equal(derivations() - seen, 1, 'other identities list again');
  assert.deepEqual(listed(), ['svg-A', 'svg-B']);
  stop();
});

// A3, second consumer (`app/feature-search/preview-actions.js`): the list's
// watcher no longer runs after such a Generate. The mounted Result's ready
// hook refreshes an active search, which keeps its matches and its index.
test('an active search keeps its matches and its index across a Generate that reuses the catalog', async () => {
  const { derivations, stop } = setup();
  state.results.value = generate({ selectedFeatureTypes: ['CDS'] });
  const metrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
  globalThis.requestAnimationFrame ??= (callback) => { callback(); return 1; };
  try {
    const { computed, nextTick } = Vue;
    const { createPreviewFeatureSearch } = await import('../../gbdraw/web/js/app/feature-search/preview-actions.js');
    const search = createPreviewFeatureSearch({
      state, watch, nextTick, computed, isActiveResultReady: () => true,
      openFeatureEditorForFeature: () => {}, focusSearchQuery: () => {}, focusSearchToggle: () => {}
    });
    const builds = () => metrics.filter(({ name }) => name === 'featureSearchIndexBuildCount').length;
    state.previewFeatureSearchInput.value = 'B';
    search.applySearch();
    await nextTick();
    const matches = [...state.previewFeatureSearchMatches.value];
    assert.equal(matches.length, 1, 'the search finds the feature');
    assert.equal(builds(), 1);
    const seen = derivations();

    state.results.value = generate({ selectedFeatureTypes: ['CDS'] });
    await nextTick();
    search.handleMountedResultReady();
    await nextTick();
    assert.equal(derivations(), seen, 'the list is kept');
    assert.deepEqual(state.previewFeatureSearchMatches.value, matches, 'the search keeps its matches');
    assert.equal(builds(), 1, 'the search index is not built again');
    search.clearSearch();
    search.dispose();
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
    stop();
  }
});
