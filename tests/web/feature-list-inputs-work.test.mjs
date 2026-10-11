// S4 (perf 0.14.x, counts not timings): the Features list (state.js
// `featureList`, R-5) is derived again only when one of its inputs changes. A
// write of a Result's content (each editor edit commits one, R1) keeps the
// Result names and committed metadata the list reads, so it does not list the
// 12.6k-18.8k catalog features of a Gallery Session again; every input does.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
globalThis.window = { Vue: context.Vue };
const { watch } = context.Vue;
const { state } = await import('../../gbdraw/web/js/state.js');
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
  state.mode.value = 'circular';
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
  setFeatureVisibilityOverride(drawing.featureOverrides, state.featureList.value.rows[1], 'off');
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
