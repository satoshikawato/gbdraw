import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import {
  featureDrawnContext,
  requestFeatureVisibilityRules,
  resolveFeatureDrawn,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/app/feature-visibility.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';

// R4: tests/test_web_rule_matching.py runs the same cases through
// gbdraw/features/visibility.py::should_render_feature, which Generate uses.
const DRAWN = JSON.parse(readFileSync(new URL('../fixtures/feature_drawn_cases.json', import.meta.url), 'utf8'));

// A rendered feature of the catalog, as the Web holds it (admission stamps
// the Result's mode as `scope`, which the identity key of an edit names).
const catalogFeature = (name, index) => {
  const spec = DRAWN.features[name];
  return {
    scope: 'circular',
    svg_id: spec.drawnSelector.hash,
    type: spec.type,
    record_id: spec.recordId,
    record_key: spec.recordId,
    biological_feature_id: `${name}-${index}`,
    qualifiers: spec.qualifiers,
    drawnSelector: spec.drawnSelector
  };
};

const preparationFor = (state, evaluate = evaluatePythonRules) => createRulePreparation({
  state,
  evaluate,
  visibilityRules: () => requestFeatureVisibilityRules(state.featureVisibilityManualRules)
});

test('resolveFeatureDrawn answers every shared case as should_render_feature does', async () => {
  const features = DRAWN.cases.map((vector, index) => catalogFeature(vector.feature, index));
  // One preparation evaluates every case's rules; a match is a fact of one
  // feature and one rule, so each case reads only its own rules.
  const state = {
    extractedFeatures: { value: features },
    featureOverrides: {},
    featureVisibilityManualRules: DRAWN.cases.flatMap((vector) => vector.visibilityRules || []),
    manualSpecificRules: DRAWN.cases.flatMap((vector) => vector.colorRules || [])
  };
  assert.equal(await preparationFor(state).prepareDrawn(), true);
  DRAWN.cases.forEach((vector, index) => {
    const featureOverrides = {};
    if (vector.override) setFeatureVisibilityOverride(featureOverrides, features[index], vector.override);
    const context = featureDrawnContext({
      featureOverrides,
      featureVisibilityManualRules: vector.visibilityRules || [],
      manualSpecificRules: vector.colorRules || []
    }, { diagramOptions: { selectedFeaturesSet: vector.selectedFeatures } });
    assert.equal(resolveFeatureDrawn(features[index], context), vector.drawn, vector.name);
  });
});

const offRule = (qualifier, value) => ({ recordId: '*', featureType: '*', qualifier, value, action: 'off' });
const contextFor = (state) => featureDrawnContext(state, { diagramOptions: { selectedFeaturesSet: ['CDS'] } });

test('resolveFeatureDrawn is unknown until Python has matched the rule it reaches', async () => {
  const feature = catalogFeature('cds_fl1', 0);
  const state = {
    extractedFeatures: { value: [feature] },
    featureOverrides: {},
    featureVisibilityManualRules: [offRule('locus_tag', '^fl1$')],
    manualSpecificRules: []
  };
  assert.equal(resolveFeatureDrawn(feature, contextFor(state)), null);
  // An override decides before any rule, so it needs no match.
  setFeatureVisibilityOverride(state.featureOverrides, feature, 'on');
  assert.equal(resolveFeatureDrawn(feature, contextFor(state)), true);
  setFeatureVisibilityOverride(state.featureOverrides, feature, 'default');
  let calls = 0;
  const preparation = preparationFor(state, async (payload) => { calls += 1; return evaluatePythonRules(payload); });
  assert.equal(await preparation.prepareDrawn(), true);
  assert.equal(resolveFeatureDrawn(feature, contextFor(state)), false);
  assert.equal(preparation.prepareDrawn(), true, 'prepared matches are reused synchronously');
  assert.equal(calls, 1);
  // A rule Generate rejects leaves its matches unknown, and the preparation
  // resolves to Generate's error, which names the table row (OV-19).
  state.featureVisibilityManualRules.push(offRule('product', '('));
  const rejected = await preparation.prepareDrawn();
  assert.match(String(rejected.error?.message), /Invalid regex in feature visibility table at row 2:/);
  state.featureVisibilityManualRules.reverse();
  assert.equal(resolveFeatureDrawn(feature, contextFor(state)), null);
});

// A catalog feature that Python did not render has no drawn selector values,
// so a rule on them is not matched live; a qualifier rule still is.
test('resolveFeatureDrawn declines drawn-coordinate rules for a feature Python did not render', async () => {
  const { drawnSelector, ...biological } = catalogFeature('cds_fl1', 0);
  const state = {
    extractedFeatures: { value: [] },
    biologicalFeatures: { value: [biological] },
    featureOverrides: {},
    featureVisibilityManualRules: [offRule('location', '^100\\.\\.700$')],
    manualSpecificRules: []
  };
  assert.equal(await preparationFor(state).prepareDrawn(), true);
  assert.equal(resolveFeatureDrawn(biological, contextFor(state)), null);
  state.featureVisibilityManualRules.splice(0, 1, offRule('locus_tag', 'FL1'));
  assert.equal(await preparationFor(state).prepareDrawn(), true);
  assert.equal(resolveFeatureDrawn(biological, contextFor(state)), false);
});
