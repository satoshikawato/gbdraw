// The Nucleotide and Amino acid fields of Feature search read the features' sequences from Python through the
// catalog-sequence door (D5): once per shown Result, never for All or another field, a stale reply writes nothing,
// a failure keeps its error with Retry, and the amino acid sequence comes only from the reply (no /translation
// fallback in the page).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { featureSequenceKey } from '../../gbdraw/web/js/services/catalog-sequences.js';
import { diagnosticError } from '../../gbdraw/web/js/utils/error-normalization.js';
import {
  buildFeatureSearchIndex,
  runFeatureSearch
} from '../../gbdraw/web/js/app/feature-search/search-core.js';

globalThis.requestAnimationFrame = (callback) => {
  callback();
  return 1;
};
globalThis.cancelAnimationFrame = () => {};

const ref = (value) => ({ value });
const featureRow = (id, start, end, extra = {}) => ({
  id: `rec-1\u0000${id}`,
  svg_id: `svg-${id}`,
  recordKey: 'rec-1',
  biologicalFeatureId: id,
  type: 'CDS',
  record_id: 'REC',
  start,
  end,
  strand: 1,
  qualifiers: { product: [`${id} product`] },
  ...extra
});
const kinase = featureRow('kinase', 0, 12, { qualifiers: { product: ['kinase'], translation: ['MSGYRA'] } });
const phosphatase = featureRow('phosphatase', 20, 29);
const reply = (row, nucleotideSequence, aminoAcidSequence) => [
  featureSequenceKey({ recordKey: row.recordKey, biologicalFeatureId: row.biologicalFeatureId }),
  { recordKey: row.recordKey, biologicalFeatureId: row.biologicalFeatureId, nucleotideSequence, aminoAcidSequence,
    sequenceWarnings: [] }
];
const sequences = new Map([
  reply(kinase, 'ATGTCCGGCTAC', 'MKKK'),
  reply(phosphatase, 'ATGCCCTGA', 'MP')
]);
const rendered = new Set([kinase.svg_id, phosphatase.svg_id]);

test('the sequence fields search only the sequences read from Python', () => {
  const search = (field, query, withSequences) => runFeatureSearch({
    features: [kinase, phosphatase], renderedFeatureIds: rendered, query, field,
    searchIndex: buildFeatureSearchIndex({ features: [kinase, phosphatase], sequences: withSequences ? sequences : null })
  }).matches;
  assert.deepEqual(search('nucleotide', 'TCCGGN', true), [kinase.svg_id]);
  assert.deepEqual(search('amino-acid', 'MKKK', true), [kinase.svg_id]);
  // The amino acid sequence is Python's: the page never falls back to the /translation qualifier.
  assert.deepEqual(search('amino-acid', 'MSGYRA', true), []);
  assert.deepEqual(search('amino-acid', 'MSGYRA', false), []);
  // A sequence a row still carried is not read either: the catalog row is not a sequence source.
  const carried = { ...phosphatase, nucleotide_sequence: 'GGGGGG', amino_acid_sequence: 'WWW' };
  assert.deepEqual(runFeatureSearch({
    features: [carried], renderedFeatureIds: rendered, query: 'GGGGGG', field: 'nucleotide'
  }).matches, []);
  // All excludes sequence content (PD-OI-056).
  assert.deepEqual(search('all', 'ATGTCC', true), []);
});

/** A deferred reply the test settles by hand. */
const deferred = () => {
  let resolve;
  let reject;
  const promise = new Promise((res, rej) => { resolve = res; reject = rej; });
  return { promise, resolve, reject };
};
const settle = () => new Promise((resolve) => setImmediate(resolve));

const createFixture = () => {
  const watchers = [];
  const svg = {
    matches: (selector) => selector === 'svg',
    classList: { length: 0, toggle() {}, remove() {} },
    getAttribute: () => null,
    setAttribute() {},
    removeAttribute() {},
    querySelectorAll: () => [],
    getElementById: () => null
  };
  const state = withDrawings({
    svgContainer: ref(svg),
    canvasContainerRef: ref(null),
    canvasPan: { x: 0, y: 0 },
    selectedResultIndex: ref(0),
    svgContent: ref('<svg></svg>'),
    extractedFeatures: ref([kinase, phosphatase]),
    featureList: ref({ rows: [kinase, phosphatase] }),
    featureListState: (row) => ({ drawn: true, renderedId: row.svg_id }),
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
    previewFeatureSearchVisible: ref(true),
    clickedFeature: ref(null)
  });
  const reads = [];
  let next = () => Promise.resolve(sequences);
  return { state, watchers, reads, setNext: (fn) => { next = fn; }, ports: {
    state,
    watch(source, callback) { watchers.push({ source, callback }); },
    nextTick(callback) {
      if (typeof callback === 'function') callback();
      return Promise.resolve();
    },
    computed(getter) { return { get value() { return getter(); } }; },
    isActiveResultReady: () => true,
    openFeatureEditorForFeature() {},
    focusSearchQuery() {},
    focusSearchToggle() {},
    readFeatureSequences(features) {
      reads.push(features);
      return next();
    }
  } };
};

const { createPreviewFeatureSearch } = await import('../../gbdraw/web/js/app/feature-search/preview-actions.js');
const searchFor = (fixture, search, field, query) => {
  fixture.state.previewFeatureSearchField.value = field;
  fixture.state.previewFeatureSearchInput.value = query;
  search.applySearch();
};

test('a sequence field reads every searchable feature once per shown Result', async () => {
  const fixture = createFixture();
  const search = createPreviewFeatureSearch(fixture.ports);
  const pending = deferred();
  fixture.setNext(() => pending.promise);

  searchFor(fixture, search, 'nucleotide', 'TCCGGN');
  assert.deepEqual(fixture.reads, [[
    { recordKey: 'rec-1', biologicalFeatureId: 'kinase', bases: 12 },
    { recordKey: 'rec-1', biologicalFeatureId: 'phosphatase', bases: 9 }
  ]]);
  assert.equal(search.previewFeatureSearchStatusText.value, 'Reading the sequence…');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, []);
  pending.resolve(sequences);
  await settle();
  assert.equal(search.previewFeatureSearchSequenceFailure.value, null);
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [kinase.svg_id]);
  assert.equal(search.previewFeatureSearchStatusText.value, '1 of 1 match');

  // Another sequence search on the same Result reads nothing again.
  searchFor(fixture, search, 'amino-acid', 'MP');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [phosphatase.svg_id]);
  assert.equal(fixture.reads.length, 1);

  // All and the other fields never read, and stay synchronous.
  searchFor(fixture, search, 'all', 'kinase');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [kinase.svg_id]);
  searchFor(fixture, search, 'label', 'phosphatase');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [phosphatase.svg_id]);
  assert.equal(fixture.reads.length, 1);

  // A new Result (the Features list changes) drops the read sequences, so the next sequence search reads again.
  fixture.state.featureList.value = { rows: [kinase, phosphatase] };
  fixture.watchers[1].callback();
  searchFor(fixture, search, 'nucleotide', 'ATGCCC');
  await settle();
  assert.equal(fixture.reads.length, 2);
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [phosphatase.svg_id]);
});

test('a field with no search reads nothing; All never reads', () => {
  const fixture = createFixture();
  const search = createPreviewFeatureSearch(fixture.ports);
  fixture.state.previewFeatureSearchField.value = 'nucleotide';
  search.refreshSearch();
  searchFor(fixture, search, 'all', 'ATGTCC');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, []);
  assert.equal(fixture.reads.length, 0);
});

test('a stale reply writes nothing', async () => {
  const fixture = createFixture();
  const search = createPreviewFeatureSearch(fixture.ports);
  fixture.setNext(() => Promise.resolve(null));
  searchFor(fixture, search, 'nucleotide', 'TCCGGN');
  await settle();
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, []);
  assert.equal(search.previewFeatureSearchSequenceFailure.value, null);
  assert.notEqual(search.previewFeatureSearchStatusText.value, 'Reading the sequence…');
  assert.notEqual(search.previewFeatureSearchStatusText.value, '1 of 1 match');
});

test('a reply that arrives after a new Result is dropped', async () => {
  const fixture = createFixture();
  const search = createPreviewFeatureSearch(fixture.ports);
  const pending = deferred();
  fixture.setNext(() => pending.promise);
  searchFor(fixture, search, 'nucleotide', 'TCCGGN');
  fixture.setNext(() => new Promise(() => {}));
  fixture.watchers[1].callback();
  pending.resolve(sequences);
  await settle();
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, []);
});

test('a failed read keeps its error and Retry reads again', async () => {
  const fixture = createFixture();
  const search = createPreviewFeatureSearch(fixture.ports);
  fixture.setNext(() => Promise.reject(diagnosticError('SEQUENCE_REQUIRES_GENERATE', {},
    { stage: 'helper', operation: 'readFeatureSequences' })));
  searchFor(fixture, search, 'nucleotide', 'TCCGGN');
  await settle();
  const failure = search.previewFeatureSearchSequenceFailure.value?.error;
  assert.ok(failure, 'the failure is kept');
  assert.equal(failure.code, 'SEQUENCE_REQUIRES_GENERATE');
  assert.equal(failure.operation, 'readFeatureSequences');
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, []);
  assert.notEqual(search.previewFeatureSearchStatusText.value, 'Reading the sequence…');

  fixture.setNext(() => Promise.resolve(sequences));
  search.previewFeatureSearchSequenceFailure.value.retry();
  assert.equal(search.previewFeatureSearchSequenceFailure.value, null);
  await settle();
  assert.equal(fixture.reads.length, 2);
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [kinase.svg_id]);

  // A search on another field ends the failed state.
  fixture.setNext(() => Promise.reject(new Error('worker lost')));
  fixture.watchers[1].callback();
  searchFor(fixture, search, 'nucleotide', 'TCCGGN');
  await settle();
  assert.ok(search.previewFeatureSearchSequenceFailure.value);
  searchFor(fixture, search, 'label', 'kinase');
  assert.equal(search.previewFeatureSearchSequenceFailure.value, null);
  assert.deepEqual(fixture.state.previewFeatureSearchMatches.value, [kinase.svg_id]);
});
