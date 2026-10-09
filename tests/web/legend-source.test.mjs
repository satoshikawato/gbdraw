import assert from 'node:assert/strict';
import test from 'node:test';

import { resultCatalogFeatures } from '../../gbdraw/web/js/services/feature-catalog.js';
import {
  featureDrawnContext,
  legendRowsShowable,
  resultLegendRowKeys,
  resultLegendSources,
  sameLegendSources,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/services/feature-visibility.js';
import { recordRuleMatches, ruleKey } from '../../gbdraw/web/js/services/rule-matchers.js';

// OV-42, OV-43 (Owner decision 2026-10-06, option A): a live edit that changes
// what a Result's Legend derives from asks for the automatic rerender, so the
// source must change with the drawn types, their first-drawn order, and the
// specific color rules a drawn feature uses, and with nothing else.
const anchorProfile = { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' };
const biological = (id, type, start) => ({
  recordKey: 'REC1', biologicalFeatureId: id, record_id: 'REC1', type, start, end: start + 30, strand: 1,
  anchorProfile, qualifiers: { locus_tag: [id] }
});
const rendered = (id) => ({
  svgId: `svg-${id}`, recordKey: 'REC1', biologicalFeatureId: id, fillColor: '#000000',
  drawnSelector: { hash: `svg-${id}`, location: null, recordLocation: null }
});
// Python draws the CDS A first, then the repeat_region R, then the CDS B.
const catalog = {
  schema: 5,
  items: [{
    resultIndex: 0,
    resultName: 'result-0.svg',
    recordKeys: ['REC1'],
    biologicalFeatures: [biological('A', 'CDS', 0), biological('R', 'repeat_region', 100), biological('B', 'CDS', 200)],
    features: ['A', 'R', 'B'].map(rendered),
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  }]
};
const stateOf = (featureOverrides = {}, rules = []) => ({
  featureCatalog: { value: catalog },
  results: { value: [{ name: 'result-0.svg', content: '<svg />' }] },
  selectedResultIndex: { value: 0 },
  generatedMode: { value: 'circular' },
  featureOverrides,
  featureVisibilityManualRules: [],
  manualSpecificRules: rules
});
const sourcesOf = (hidden = [], rules = []) => {
  const featureOverrides = {};
  const { biological: features } = resultCatalogFeatures(stateOf());
  hidden.forEach((id) => setFeatureVisibilityOverride(
    featureOverrides, features.find((feature) => feature.biological_feature_id === id), 'off'
  ));
  const state = stateOf(featureOverrides, rules);
  const context = {
    ...featureDrawnContext(state, { diagramOptions: { selectedFeaturesSet: ['CDS', 'repeat_region'] } }),
    colorRules: rules
  };
  return resultLegendSources(state, context);
};
const rule = (cap, color = '#ff0000') => ({ feat: 'CDS', qual: 'locus_tag', val: '^A$', color, cap });

test('hiding the first drawn feature of a type that Python draws first changes the source', () => {
  assert.equal(sameLegendSources(sourcesOf(), sourcesOf(['A'])), false, 'the CDS row moves behind repeat_region');
  assert.equal(sameLegendSources(sourcesOf(), sourcesOf(['A', 'B'])), false, 'the CDS row is gone');
});

test('hiding a later feature of a drawn type leaves the source as it is', () => {
  assert.equal(sameLegendSources(sourcesOf(), sourcesOf(['B'])), true);
});

test('a rule match that is not known yet reads as a changed source', () => {
  // The match of a rule is Python's, prepared before the edit commits; until
  // then no source equals another, so a missed change cannot hide here.
  const unknown = sourcesOf([], [rule('alpha')]);
  assert.deepEqual(unknown, [null]);
  assert.equal(sameLegendSources(unknown, unknown), false);
  assert.equal(sameLegendSources([], []), true);
  assert.equal(sameLegendSources(['a'], ['a', 'b']), false);
});

test('a type without a captioned rule needs no match, whatever its rules match', () => {
  // The rule has no caption, so Python draws the type's default row either way.
  const uncaptioned = { feat: 'CDS', qual: 'locus_tag', val: '^A$', color: '#ff0000', cap: '' };
  const [source] = sourcesOf([], [uncaptioned]);
  assert.notEqual(source, null);
  assert.equal(sameLegendSources(sourcesOf(), sourcesOf([], [uncaptioned])), true, 'an uncaptioned rule draws no row');
  assert.equal(sameLegendSources(sourcesOf([], [uncaptioned]), sourcesOf([], [{ ...uncaptioned, val: '^B$' }])), true);
});

// U3a 1d (R14-8): a rule change asks for the automatic rerender only when the
// Legend rows it regroups cannot be shown live. Each row of the displayed
// Result keeps its features together: it is relabeled (its old row hidden,
// the new one shown), merged into a row it joins, or recolored. A row that
// loses part of its features, a row built from parts of others, a merge into
// a new row, an unknown match, or a row of another batch Result that changes
// asks Python. `rows` are each Result's features, `[id, type]`; the rules
// match by locus_tag.
const batchStateOf = (rows, rules) => {
  const items = rows.map((features, resultIndex) => {
    const record = `REC${resultIndex + 1}`;
    const drawn = features.map(([id, type], index) => ({
      ...biological(id, type, index * 100), recordKey: record, record_id: record
    }));
    return {
      resultIndex,
      resultName: `result-${resultIndex}.svg`,
      recordKeys: [record],
      biologicalFeatures: drawn,
      features: features.map(([id]) => ({ ...rendered(id), recordKey: record })),
      orthogroups: [],
      annotations: [],
      comparisonMatches: []
    };
  });
  return {
    featureCatalog: { value: { schema: 5, items } },
    results: { value: rows.map((_, index) => ({ name: `result-${index}.svg`, content: '<svg />' })) },
    selectedResultIndex: { value: 0 },
    generatedMode: { value: 'circular' },
    featureOverrides: {},
    featureVisibilityManualRules: [],
    manualSpecificRules: rules
  };
};
const tagRule = (pattern, cap, color = '#ff0000') => ({ feat: 'CDS', qual: 'locus_tag', val: `^(${pattern})$`, color, cap });
// Records each rule's match on the rendered features of every Result, as the
// rule preparation does with Python's answers.
const prepareMatches = (state, rules) => {
  state.results.value.forEach((_, index) => {
    const features = [...resultCatalogFeatures(state, index).renderedByIdentity.values()];
    const keys = rules.map(ruleKey);
    recordRuleMatches(features, keys, (featureIndex) => {
      const tag = features[featureIndex].locus_tag || features[featureIndex].biological_feature_id;
      const matched = rules.flatMap((rule, ruleIndex) => (
        rule.feat === features[featureIndex].type && new RegExp(rule.val).test(tag) ? [ruleIndex] : []
      ));
      return { matched, priorities: matched.map(() => 0), declined: [] };
    });
  });
};
const rowKeysOf = (rows, rules, { prepared = true } = {}) => {
  const state = batchStateOf(rows, rules);
  if (prepared) prepareMatches(state, rules);
  const context = {
    ...featureDrawnContext(state, { diagramOptions: { selectedFeaturesSet: ['CDS', 'repeat_region'] } }),
    colorRules: rules
  };
  return resultLegendRowKeys(state, context);
};
/** @param {string[]} shown */
const showing = (shown) => (key) => shown.includes(key);
const ONE = [[['A', 'CDS'], ['R', 'repeat_region'], ['B', 'CDS']]];

test('a whole-row relabel, a merge into a shown row and a recolor are shown live', () => {
  const alpha = rowKeysOf(ONE, [tagRule('A|B', 'alpha')]);
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'beta')]),
    { displayed: 0, shows: showing(['beta', 'repeat_region']) }), true, 'relabel: the old row hidden, the new row shown');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'beta')]),
    { displayed: 0, shows: showing(['alpha', 'beta', 'repeat_region']) }), false, 'relabel while the old row stays shown');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'beta')]),
    { displayed: 0, shows: showing(['repeat_region']) }), false, 'relabel without the new row');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'beta')]),
    { displayed: 0, shows: showing(['beta', 'repeat_region']), placed: (key, next) => key === 'alpha' && next === 'beta' }),
  true, 'relabel whose new row takes the old row\'s place');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'beta')]),
    { displayed: 0, shows: showing(['beta', 'repeat_region']), placed: () => false }),
  false, 'relabel whose new row is appended: Python keeps it in the old row\'s place');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A|B', 'alpha', '#00ff00')]),
    { displayed: 0, shows: showing(['alpha', 'repeat_region']) }), true, 'recolor');
  const two = rowKeysOf(ONE, [tagRule('A', 'alpha'), tagRule('B', 'beta')]);
  assert.equal(legendRowsShowable(two, rowKeysOf(ONE, [tagRule('A', 'beta'), tagRule('B', 'beta')]),
    { displayed: 0, shows: showing(['beta', 'repeat_region']) }), true, 'merge into the shown row beta');
  // The palette row of a type with features is renamed by rules of its features (P-1).
  const plain = rowKeysOf(ONE, []);
  assert.equal(legendRowsShowable(plain, rowKeysOf(ONE, [tagRule('A|B', 'gamma')]),
    { displayed: 0, shows: showing(['gamma', 'repeat_region']) }), true, 'the CDS row relabeled');
});

test('a split, a merge into a new row, a swap, an unknown match or a generated caption asks Python', () => {
  const alpha = rowKeysOf(ONE, [tagRule('A|B', 'alpha')]);
  const all = showing(['alpha', 'beta', 'gamma', 'other proteins', 'repeat_region']);
  assert.equal(legendRowsShowable(alpha, rowKeysOf(ONE, [tagRule('A', 'alpha')]), { displayed: 0, shows: all }), false,
    'B leaves the row for other proteins');
  assert.equal(legendRowsShowable(rowKeysOf(ONE, []), rowKeysOf(ONE, [tagRule('A', 'gamma')]), { displayed: 0, shows: all }), false,
    'a row built from part of the CDS row');
  const two = rowKeysOf(ONE, [tagRule('A', 'alpha'), tagRule('B', 'beta')]);
  assert.equal(legendRowsShowable(two, rowKeysOf(ONE, [tagRule('A', 'gamma'), tagRule('B', 'gamma')]),
    { displayed: 0, shows: showing(['gamma', 'repeat_region']) }), false, 'a merge into a new row');
  assert.equal(legendRowsShowable(two, rowKeysOf(ONE, [tagRule('A', 'beta'), tagRule('B', 'alpha')]), { displayed: 0, shows: all }), false,
    'a swap of two rows');
  assert.deepEqual(rowKeysOf(ONE, [tagRule('A|B', 'alpha')], { prepared: false }), [null], 'a match not known yet');
  assert.equal(legendRowsShowable([null], alpha, { displayed: 0, shows: all }), false);
  assert.deepEqual(rowKeysOf(ONE, [tagRule('A', 'repeat_region')]), [null], 'Python suffixes a caption of a generated row');
});

test('another batch Result changes no row live; its colors may change', () => {
  const BATCH = [...ONE, [['C', 'CDS'], ['D', 'CDS']]];
  const alpha = rowKeysOf(BATCH, [tagRule('A|B|C|D', 'alpha')]);
  const shows = showing(['beta', 'repeat_region']);
  assert.equal(legendRowsShowable(alpha, rowKeysOf(BATCH, [tagRule('A|B|C|D', 'beta')]), { displayed: 0, shows }), false,
    'Result 2 draws the relabeled row');
  assert.equal(legendRowsShowable(alpha, rowKeysOf(BATCH, [tagRule('A|B|C|D', 'alpha', '#00ff00')]),
    { displayed: 0, shows: showing(['alpha', 'repeat_region']) }), true, 'a recolor');
  const own = rowKeysOf(BATCH, [tagRule('A|B', 'alpha')]);
  assert.equal(legendRowsShowable(own, rowKeysOf(BATCH, [tagRule('A|B', 'beta')]), { displayed: 0, shows }), true,
    'a row only the displayed Result draws');
});
