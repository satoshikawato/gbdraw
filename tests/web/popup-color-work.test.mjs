// OV-225 (R14-3, a work allowlist, not timings): from a popup color pick to
// the scope dialog, the app runs exactly the listed History stages and Worker
// requests. The dialog reads only the saved rules, so the rules the edit may
// add (the clicked feature's hash and label rules) are prepared when it
// commits. On a Session just loaded, that preparation started the Worker's
// Python runtime before the dialog opened (about 6 s on Vnig_TUMSAT-TG-2018).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation, runWhenPrepared } from '../../gbdraw/web/js/app/rule-matching.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });

const features = Array.from({ length: 4 }, (_, index) => ({
  id: `feature-${index}`, svg_id: `f${index}`, type: 'CDS', qualifiers: { product: [`p${index}`] }
}));
const savedRule = { feat: 'CDS', qual: 'product', val: '^p3$', color: '#111111', cap: 'p3' };

// Python's answer for `hash` and `product` rules (gbdraw/web_support/rule_matching.py).
const evaluateRules = async ({ features: payloads, rules }) => {
  const matches = payloads.map((payload) => rules.flatMap((rule, index) => {
    const values = rule.qual === 'hash' ? [payload.selector.hash] : payload.qualifiers[rule.qual] || [];
    return values.some((value) => new RegExp(rule.val).test(value)) ? [index] : [];
  }));
  return { matches, priorities: matches.map((row) => row.map(() => 0)) };
};

const setup = ({ siblings }) => {
  /** @type {string[]} */
  const stages = [];
  /** @type {string[] | null} */
  let atDialog = null;
  const featureStyleScopeDialog = {
    _show: false,
    get show() { return this._show; },
    set show(value) { if (value && !this._show) atDialog = [...stages]; this._show = value; }
  };
  const manualSpecificRules = [savedRule];
  const state = withDrawings({
    extractedFeatures: ref(features), biologicalFeatures: ref(features), manualSpecificRules,
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#cccccc' }]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()), results: ref([]), selectedResultIndex: ref(0), svgContainer: ref(null),
    clickedFeature: ref({ feat: features[0], svg_id: features[0].svg_id, legendName: '' }), featureStyleScopeDialog,
    resetColorDialog: {}, legendRenameDialog: {}, originalLegendOrder: ref([]), originalLegendColors: ref({}),
    originalSvgStroke: ref({ color: null, width: null }), appliedPaletteColors: ref({ CDS: '#cccccc' }),
    skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false)
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
      commitSpecificRules: async () => { stages.push('commitSpecificRules'); return true; },
      countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null,
      findFeaturesWithSameDisplayedLabel: () => [], findFeaturesWithSameIndividualLabel: () => [],
      findFeaturesWithSameLegendItem: (feature) => (siblings ? features.filter((other) => other !== feature) : []),
      findMatchingRegexRule: () => null, getDisplayedFeatureLabel: () => '',
      getIndividualFeatureLabel: (feature) => feature.qualifiers.product[0],
      effectiveLegendCaptions: () => () => 'CDS',
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      getLabelSpecificRule: (feature, label) => (label ? { feat: feature.type, qual: 'product', val: `^${label}$` } : null),
      getLegendRowRules: () => []
    }
  });
  // app-setup.js: `undoableAction('Change feature color', updateClickedFeatureColor)`.
  const pick = (color) => history.runUndoable('Change feature color', () => actions.updateClickedFeatureColor(color));
  return { stages, preparation, pick, dialogStages: () => atDialog, featureStyleScopeDialog };
};

test('a popup color pick opens the scope dialog after the History capture alone', async () => {
  const { stages, preparation, pick, dialogStages, featureStyleScopeDialog } = setup({ siblings: true });
  assert.equal(await preparation.prepare([savedRule]), true);
  stages.length = 0;
  await pick('#123456');
  assert.equal(featureStyleScopeDialog.show, true);
  assert.deepEqual(dialogStages(), ['history:buildIntent', 'history:signature']);
});

test('a popup color pick without a dialog prepares the rules it commits', async () => {
  const { stages, preparation, pick, dialogStages } = setup({ siblings: false });
  assert.equal(await preparation.prepare([savedRule]), true);
  stages.length = 0;
  await pick('#123456');
  assert.equal(dialogStages(), null);
  assert.deepEqual(stages, [
    'history:buildIntent', 'history:signature', 'worker:evaluateRules:color:hash+product', 'commitSpecificRules',
    'history:buildIntent', 'history:signature'
  ]);
});
