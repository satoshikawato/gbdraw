import assert from 'node:assert/strict';

const ref = (value) => ({ value });
const watchers = [];
const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    metrics.push(metric);
  }
};
globalThis.requestAnimationFrame = (callback) => {
  callback();
  return 1;
};
globalThis.cancelAnimationFrame = () => {};

const attributes = new Map([
  ['id', 'feature-1'],
  ['data-gbdraw-feature-id', 'feature-1']
]);
const classTokens = new Set();
const featureElement = {
  id: 'feature-1',
  classList: {
    get length() { return classTokens.size; },
    toggle(token, enabled) {
      if (enabled) classTokens.add(token);
      else classTokens.delete(token);
    },
    remove(token) { classTokens.delete(token); }
  },
  getAttribute: (name) => attributes.get(name) ?? null,
  setAttribute: (name, value) => attributes.set(name, String(value)),
  removeAttribute: (name) => attributes.delete(name)
};
const svgClassTokens = new Set();
const svg = {
  matches: (selector) => selector === 'svg',
  classList: {
    get length() { return svgClassTokens.size; },
    toggle(token, enabled) {
      if (enabled) svgClassTokens.add(token);
      else svgClassTokens.delete(token);
    },
    remove(token) { svgClassTokens.delete(token); }
  },
  getAttribute: () => null,
  setAttribute() {},
  removeAttribute() {},
  querySelectorAll(selector) {
    if (String(selector).includes('data-gbdraw-feature-id')) return [featureElement];
    return [];
  },
  getElementById: (id) => (id === 'feature-1' ? featureElement : null)
};

const feature = {
  svg_id: 'feature-1',
  type: 'CDS',
  record_id: 'record-1',
  start: 1,
  end: 9,
  strand: 1,
  qualifiers: { product: ['kinase'] }
};
const state = {
  svgContainer: ref(svg),
  canvasContainerRef: ref(null),
  canvasPan: { x: 0, y: 0 },
  selectedResultIndex: ref(0),
  svgContent: ref('<svg></svg>'),
  extractedFeatures: ref([feature]),
  featuresBySvgId: ref(new Map([['feature-1', feature]])),
  orthogroups: ref([]),
  adv: { rich_feature_popup: true },
  previewFeatureSearchInput: ref(''),
  previewFeatureSearchQuery: ref(''),
  previewFeatureSearchField: ref('all'),
  previewFeatureSearchQualifierKey: ref(''),
  previewFeatureSearchUseRegex: ref(false),
  previewFeatureSearchMatches: ref([]),
  previewFeatureSearchMatchDetails: ref({}),
  previewFeatureSearchActiveIndex: ref(-1),
  previewFeatureSearchError: ref(''),
  previewFeatureSearchRenderedCount: ref(0),
  clickedFeature: ref(null),
  showRightDrawer: ref(false)
};

const { createPreviewFeatureSearch } = await import(
  '../../gbdraw/web/js/app/feature-search/preview-actions.js'
);
const search = createPreviewFeatureSearch({
  state,
  watch(source, callback) {
    watchers.push({ source, callback });
  },
  nextTick(callback) {
    if (typeof callback === 'function') callback();
    return Promise.resolve();
  },
  computed(getter) {
    return { get value() { return getter(); } };
  },
  reactive: (value) => value,
  previewRuntime: { isActiveResultReady: () => true },
  openFeatureEditorForFeature() {}
});

const metricTotal = (name) => metrics
  .filter((metric) => metric.name === name)
  .reduce((total, metric) => total + Number(metric.value || 0), 0);

assert.equal(metricTotal('featureSearchIndexBuildCount'), 0);
assert.equal(metricTotal('featureDomFullScanCount'), 0);
search.refreshSearch();
assert.equal(metricTotal('featureSearchIndexBuildCount'), 0);
assert.equal(metricTotal('featureDomFullScanCount'), 0);

state.previewFeatureSearchInput.value = 'kinase';
search.applySearch();
assert.equal(metricTotal('featureSearchIndexBuildCount'), 1);
assert.equal(metricTotal('featureDomFullScanCount'), 1);
assert.deepEqual(state.previewFeatureSearchMatches.value, ['feature-1']);

search.applySearch();
assert.equal(metricTotal('featureSearchIndexBuildCount'), 1);
assert.equal(metricTotal('featureDomFullScanCount'), 1);

watchers[1].callback();
assert.equal(metricTotal('featureSearchIndexBuildCount'), 2);
assert.equal(metricTotal('featureDomFullScanCount'), 1);

watchers[2].callback();
search.refreshSearch();
assert.equal(metricTotal('featureSearchIndexBuildCount'), 2);
assert.equal(metricTotal('featureDomFullScanCount'), 2);

console.log('preview feature search tests passed');
