import assert from 'node:assert/strict';
import test from 'node:test';

import { resultCatalogFeatures } from '../../gbdraw/web/js/services/feature-catalog.js';
import {
  featureDrawnContext,
  resultLegendSources,
  sameLegendSources,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/services/feature-visibility.js';

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
