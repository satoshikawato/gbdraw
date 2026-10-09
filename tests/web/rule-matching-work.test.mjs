// OV-193 (R14-3, counts not timings): with K features and K `hash` rules, a
// whole-feature pass over the rule matches and the "All features with legend
// item" commit on the row those rules draw read each feature and build each
// rule key a bounded number of times, so their work grows linearly with K.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { createSvgStyles } from '../../gbdraw/web/js/app/svg-styles.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });
const hashOf = (index) => `f${index.toString(16).padStart(8, '0')}`;

// Property reads of the features and rules, and JSON.stringify calls (rule keys).
const counting = () => {
  const counts = { reads: 0, keys: 0 };
  const stringify = JSON.stringify;
  return {
    counts,
    feature: (target) => new Proxy(target, { get: (object, key) => { counts.reads += 1; return Reflect.get(object, key); } }),
    measure: async (run) => {
      counts.reads = 0;
      counts.keys = 0;
      JSON.stringify = (...args) => { counts.keys += 1; return stringify(...args); };
      try { await run(); } finally { JSON.stringify = stringify; }
      return { ...counts };
    }
  };
};

// Python's answer for `hash` rules (gbdraw/web_support/rule_matching.py):
// a feature matches the rules whose value is its hash.
const evaluateHashRules = async ({ features, rules }) => {
  const byValue = new Map();
  rules.forEach((rule, index) => byValue.set(rule.val, [...(byValue.get(rule.val) || []), index]));
  const matches = features.map((feature) => byValue.get(feature.selector.hash) || []);
  // Dense rows of equal ranks read the same through the aligned reader, so the
  // fake also answers the dense reader of the parent commit (failing-first by counts).
  return { matches, priorities: matches.map(() => rules.map(() => 0)), winners: [] };
};

const setup = (count, count_) => {
  const features = Array.from({ length: count }, (_, index) => count_.feature({
    id: `feature-${index}`, svg_id: hashOf(index), type: 'CDS', qualifiers: { product: [`p${index}`] }
  }));
  const rules = features.map((_, index) => count_.feature({ feat: 'CDS', qual: 'hash', val: hashOf(index), color: '#111111', cap: 'CDS' }));
  const state = withDrawings({
    extractedFeatures: ref(features), biologicalFeatures: ref(features), manualSpecificRules: rules,
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#111111' }]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set())
  });
  return { features, rules, state };
};

const linear = (small, large) => {
  for (const key of ['reads', 'keys']) {
    assert.ok(large[key] <= small[key] * 4.5, `${key}: ${small[key]} at K=50, ${large[key]} at K=200`);
  }
};

test('a whole-feature pass over K features and K hash rules does O(K) work', async () => {
  const measured = [];
  for (const count of [50, 200]) {
    const count_ = counting();
    const { features, rules, state } = setup(count, count_);
    const preparation = createRulePreparation({ state, evaluate: evaluateHashRules });
    assert.equal(await preparation.prepare(rules), true);
    const elements = features.map((feature) => {
      const attributes = {
        id: feature.svg_id, 'data-gbdraw-feature-id': feature.svg_id, 'data-gbdraw-feature-part': 'block', fill: '#000000'
      };
      return { getAttribute: (key) => attributes[key] ?? null, setAttribute: (key, value) => { attributes[key] = value; } };
    });
    const svg = { querySelectorAll: (selector) => (selector.includes('data-gbdraw-feature-id') ? elements : []) };
    Object.assign(state, {
      svgContent: ref('<svg/>'), svgContainer: ref({ querySelector: () => svg }), appliedPaletteColors: ref({ CDS: '#cccccc' }),
      featuresBySvgId: ref(new Map())
    });
    const styles = createSvgStyles({ state, watch() {}, nextTick: (fn) => fn?.(), commitActiveResultEdit: () => true });
    measured.push(await count_.measure(() => {
      assert.equal(preparation.isPrepared(rules), true);
      styles.applySpecificRulesToSvg();
    }));
    assert.equal(elements.at(-1).getAttribute('fill'), '#111111');
  }
  linear(...measured);
});

test('"All features with legend item" over K features with K hash rules does O(K) work', async () => {
  const measured = [];
  for (const count of [50, 200]) {
    const count_ = counting();
    const { features, rules, state } = setup(count, count_);
    const preparation = createRulePreparation({ state, evaluate: evaluateHashRules });
    assert.equal(await preparation.prepare(rules), true);
    let committed = null;
    const featureStyleScopeDialog = { show: true, feat: features[0], color: '#abcdef', legendName: 'CDS' };
    Object.assign(state, {
      results: ref([]), selectedResultIndex: ref(0), svgContainer: ref(null), clickedFeature: ref(null),
      featureStyleScopeDialog, resetColorDialog: {}, legendRenameDialog: {}, originalLegendOrder: ref([]),
      originalLegendColors: ref({}), originalSvgStroke: ref({ color: null, width: null }), appliedPaletteColors: ref({}),
      skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false), hasPendingPaletteDraft: ref(false)
    });
    const caption = () => 'CDS';
    const actions = createFeatureColorActions({
      state, nextTick: async () => {}, onLegendGeometryChanged: () => {}, extractLegendEntries: () => {},
      getFeatureElements: () => [], getFeatureFillElements: () => [],
      readUserDefaultColor: () => null, setDefaultColor: () => assert.fail('the row is drawn by rules'),
      ruleActions: {
        runWithRuleMatches: (_, commit) => commit(),
        commitSpecificRules: async (next) => { committed = next; return true; },
        countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null,
        findFeaturesWithSameDisplayedLabel: () => [], findFeaturesWithSameIndividualLabel: () => [],
        findFeaturesWithSameLegendItem: (feature) => features.filter((other) => other !== feature),
        findMatchingRegexRule: () => null, getDisplayedFeatureLabel: () => '', getIndividualFeatureLabel: () => '',
        getEffectiveLegendCaption: caption, effectiveLegendCaptions: () => caption,
        getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
        getLabelSpecificRule: () => null, getLegendRowRules: () => rules
      }
    });
    measured.push(await count_.measure(() => actions.handleColorScopeChoice('caption')));
    assert.equal(committed.length, count);
    assert.ok(committed.every((rule, index) => rule.val === hashOf(index) && rule.color === '#abcdef'));
  }
  linear(...measured);
});

test('the Legend item of K features, whose K hash rules name a row of another color, is read in O(K)', async () => {
  const measured = [];
  for (const count of [50, 200]) {
    const count_ = counting();
    const { features, rules, state } = setup(count, count_);
    rules.forEach((rule) => { rule.color = '#abcdef'; });
    const preparation = createRulePreparation({ state, evaluate: evaluateHashRules });
    assert.equal(await preparation.prepare(rules), true);
    Object.assign(state, { results: ref([]), files: { t_color: null }, fileLegendCaptions: ref(new Set()), originalLegendOrder: ref(['CDS']) });
    const actions = createFeatureRuleActions({
      ref, computed: (get) => ({ get value() { return get(); } }), watch() {}, state, rulePreparation: preparation,
      runUndoable: async (_, commit) => commit(), runUndoableCheckpoint: async (_, commit) => commit(),
      prepareFileLegendEntries: async () => false, projectPaletteAndRules: () => true,
      ports: { requestAutomaticRerender: () => true }, getCommittedRequest: () => null, isPatternEditAvailable: () => true
    });
    /** @type {Record<string, any>[]} */
    let siblings = [];
    measured.push(await count_.measure(() => { siblings = actions.findFeaturesWithSameLegendItem(features[0]); }));
    assert.equal(siblings.length, count - 1);
  }
  linear(...measured);
});
