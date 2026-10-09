// OV-253: Reset fill color on a feature whose caption no other feature
// shares commits without a dialog. The reset is one undoable step: its
// History step ends after the reset's rule commit, not before it.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation, runWhenPrepared } from '../../gbdraw/web/js/app/rule-matching.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });

const features = ['p0', 'p1'].map((product, index) => ({
  id: `feature-${index}`, svg_id: `f${index}`, type: 'CDS', product, qualifiers: { product: [product] }
}));
const savedRule = { feat: 'CDS', qual: 'product', val: '^p1$', color: '#111111', cap: 'p1' };
const hashRule = { feat: 'CDS', qual: 'hash', val: 'f0', color: '#222222', cap: 'p0' };

// Python's answer for `hash` and `product` rules (gbdraw/web_support/rule_matching.py).
const evaluateRules = async ({ features: payloads, rules }) => {
  const matches = payloads.map((payload) => rules.flatMap((rule, index) => {
    const values = rule.qual === 'hash' ? [payload.selector.hash] : payload.qualifiers[rule.qual] || [];
    return values.some((value) => new RegExp(rule.val).test(value)) ? [index] : [];
  }));
  return { matches, priorities: matches.map((row) => row.map(() => 0)) };
};

test('a popup fill reset without a dialog commits inside its History step', async () => {
  /** @type {string[]} */
  const stages = [];
  const manualSpecificRules = [savedRule, hashRule];
  const resetColorDialog = { show: false };
  const state = withDrawings({
    extractedFeatures: ref(features), biologicalFeatures: ref(features), manualSpecificRules,
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#cccccc' }]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()), results: ref([]), selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => null }),
    clickedFeature: ref({ feat: features[0], svg_id: features[0].svg_id, legendName: '' }),
    featureStyleScopeDialog: { show: false }, resetColorDialog, legendRenameDialog: {},
    originalLegendOrder: ref([]), originalLegendColors: ref({}), originalSvgStroke: ref({ color: null, width: null }),
    appliedPaletteColors: ref({ CDS: '#cccccc' }), skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false)
  });
  const preparation = createRulePreparation({
    state,
    evaluate: (payload) => {
      stages.push(`worker:evaluateRules:${payload.kind}:${payload.rules.map((rule) => rule.qual).join('+')}`);
      return evaluateRules(payload);
    }
  });
  const history = createHistoryManager({
    buildIntent: () => { stages.push('history:buildIntent'); return { rules: JSON.parse(JSON.stringify(manualSpecificRules)) }; },
    signatureFor: (value) => { stages.push('history:signature'); return JSON.stringify(value); },
    applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  const actions = createFeatureColorActions({
    state, nextTick: async () => {}, onLegendGeometryChanged: () => {}, extractLegendEntries: () => {},
    getFeatureElements: () => [], getFeatureFillElements: () => [],
    ruleActions: {
      runWithRuleMatches: (rules, commit) => runWhenPrepared(state, () => [preparation.prepare(rules)], commit),
      // The rule owner's commit lands after its own preparation (a later task).
      commitSpecificRules: async (rules) => {
        stages.push('commitSpecificRules');
        await new Promise((resolve) => setTimeout(resolve, 0));
        manualSpecificRules.splice(0, manualSpecificRules.length, ...rules);
        return true;
      },
      countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null,
      findFeaturesWithSameDisplayedLabel: () => [], findFeaturesWithSameIndividualLabel: () => [],
      findFeaturesWithSameLegendItem: () => [], findMatchingRegexRule: () => null, getDisplayedFeatureLabel: () => '',
      getIndividualFeatureLabel: (feature) => feature.qualifiers.product[0],
      effectiveLegendCaptions: () => () => 'CDS',
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      getLabelSpecificRule: (feature, label) => (label ? { feat: feature.type, qual: 'product', val: `^${label}$` } : null),
      getLegendRowRules: () => []
    }
  });
  assert.equal(await preparation.prepare([savedRule, hashRule]), true);
  stages.length = 0;
  // app-setup.js: `undoableAction('Reset feature color', resetClickedFeatureFillColor)`.
  await history.runUndoable('Reset feature color', () => actions.resetClickedFeatureFillColor());
  assert.equal(resetColorDialog.show, false);
  assert.deepEqual(manualSpecificRules, [savedRule]);
  assert.deepEqual(stages, [
    'history:buildIntent', 'history:signature', 'worker:evaluateRules:color:product', 'commitSpecificRules',
    'history:buildIntent', 'history:signature'
  ]);
  assert.equal(history.getUndoCount(), 1);
});
