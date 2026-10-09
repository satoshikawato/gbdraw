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

const featuresOf = (products) => products.map((product, index) => ({
  id: `feature-${index}`, svg_id: `f${index}`, type: 'CDS', product, qualifiers: { product: [`p${index}`] }
}));
const savedRule = { feat: 'CDS', qual: 'product', val: '^p3$', color: '#111111', cap: 'p3' };

// Python's answer for `hash` and `product` rules (gbdraw/web_support/rule_matching.py).
const evaluateRules = async ({ features: payloads, rules, kind }) => {
  if (kind === 'color-captions') return { rules };
  const matches = payloads.map((payload) => rules.flatMap((rule, index) => {
    const values = rule.qual === 'hash' ? [payload.selector.hash] : payload.qualifiers[rule.qual] || [];
    return values.some((value) => new RegExp(rule.val).test(value)) ? [index] : [];
  }));
  return { matches, priorities: matches.map((row) => row.map(() => 0)) };
};

const setup = ({ siblings, features = featuresOf(['p', 'p', 'p', 'p']), savedRules = [savedRule], commit = async () => true }) => {
  /** @type {string[]} */
  const stages = [];
  /** @type {string[] | null} */
  let atDialog = null;
  // A dialog records the stages run before it opened.
  const dialog = (fields = {}) => ({
    ...fields,
    _show: false,
    get show() { return this._show; },
    set show(value) { if (value && !this._show) atDialog = [...stages]; this._show = value; }
  });
  const featureStyleScopeDialog = dialog();
  const manualSpecificRules = [...savedRules];
  const state = withDrawings({
    extractedFeatures: ref(features), biologicalFeatures: ref(features), manualSpecificRules,
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#cccccc' }]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()), results: ref([]), selectedResultIndex: ref(0), svgContainer: ref({ querySelector: () => null }),
    clickedFeature: ref({ feat: features[0], svg_id: features[0].svg_id, legendName: '' }), featureStyleScopeDialog,
    resetColorDialog: dialog(), legendRenameDialog: dialog(), originalLegendOrder: ref([]), originalLegendColors: ref({}),
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
      commitSpecificRules: async (rules) => {
        stages.push('commitSpecificRules');
        if (!await commit()) return false;
        manualSpecificRules.splice(0, manualSpecificRules.length, ...rules);
        return true;
      },
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
  // app-setup.js: `undoableAction(label, action)`.
  const pick = (color) => history.runUndoable('Change feature color', () => actions.updateClickedFeatureColor(color));
  const rename = (caption) => {
    state.clickedFeature.value.legendName = caption;
    return history.runUndoable('Rename legend item', () => actions.handleLegendNameCommit());
  };
  const reset = () => history.runUndoable('Reset feature color', () => actions.resetClickedFeatureFillColor());
  return {
    stages, preparation, pick, rename, reset, history, manualSpecificRules, dialogStages: () => atDialog,
    featureStyleScopeDialog, legendRenameDialog: state.legendRenameDialog, resetColorDialog: state.resetColorDialog
  };
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

// The popup's Legend name and its fill reset open their dialogs from the saved
// rules too, and prepare the rules they may add when the choice commits.
for (const [name, run, dialogOf] of [
  ['a popup Legend rename', (setup_) => setup_.rename('Renamed'), (setup_) => setup_.legendRenameDialog],
  ['a popup fill reset', (setup_) => setup_.reset(), (setup_) => setup_.resetColorDialog]
]) {
  test(`${name} opens its dialog after the History capture alone`, async () => {
    const setup_ = setup({ siblings: true });
    assert.equal(await setup_.preparation.prepare([savedRule]), true);
    setup_.stages.length = 0;
    await run(setup_);
    assert.equal(dialogOf(setup_).show, true);
    assert.deepEqual(setup_.dialogStages(), ['history:buildIntent', 'history:signature']);
  });
}
