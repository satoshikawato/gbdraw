import assert from 'node:assert/strict';
import { withDrawings } from './helpers/drawing-state.mjs';

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

const hidden = {
  id: 'record-1\u0000hidden',
  svg_id: 'stable-hidden',
  type: 'CDS',
  record_id: 'record-1',
  start: 20,
  end: 29,
  strand: 1,
  qualifiers: { product: ['phosphatase'] }
};
const state = {
  svgContainer: ref(svg),
  canvasContainerRef: ref(null),
  canvasPan: { x: 0, y: 0 },
  selectedResultIndex: ref(0),
  svgContent: ref('<svg></svg>'),
  extractedFeatures: ref([]),
  featureList: ref({ rows: [hidden] }),
  featureListState: () => ({ drawn: false, renderedId: '' }),
  orthogroups: ref([]),
  richFeaturePopup: ref(true),
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
  clickedFeature: ref(null)
};

const { createPreviewFeatureSearch } = await import(
  '../../gbdraw/web/js/app/feature-search/preview-actions.js'
);
const search = createPreviewFeatureSearch({
  state: withDrawings(state),
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
  isActiveResultReady: () => true,
  openFeatureEditorForFeature() {}
});

// OV-577: a feature without a rendered label has no label text to look up, so a
// bulk override stored for the text `undefined` must not become its search text.
state.labelTextBulkOverrides = { undefined: 'zzmarkerzz' };
state.previewFeatureSearchInput.value = 'zzmarkerzz';
search.applySearch();
assert.deepEqual(state.previewFeatureSearchMatches.value, []);

// An empty source text and a source that is not text name no bulk override either (review U7).
hidden.label = '';
hidden.gene = 5;
state.labelTextBulkOverrides = { '': 'zzemptyzz', 5: 'zzfivezz' };
// The watchers of the rows and of the override tables drop the search index.
watchers.forEach(({ callback }) => callback());
for (const text of ['zzemptyzz', 'zzfivezz']) {
  state.previewFeatureSearchInput.value = text;
  search.applySearch();
  assert.deepEqual(state.previewFeatureSearchMatches.value, [], `no feature matches ${text}`);
}

console.log('preview feature search sources tests passed');
