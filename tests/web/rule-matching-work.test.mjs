// OV-193 (R14-3, counts not timings): with K features and K `hash` rules, a
// whole-feature pass over the rule matches and the "All features with legend
// item" commit read each feature and build each rule key a bounded number of
// times, so their work grows linearly with K.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { displayedFeatureAddressing } from '../../gbdraw/web/js/services/feature-override-identity.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { withDrawings } from './helpers/drawing-state.mjs';

// The Legend owner's writer of a row's own color (app/feature-editor.js wires it); these tests record its calls.
const legendRowColorWrites = [];
/** @param {unknown[]} write */
const recordLegendRowColor = (...write) => { legendRowColorWrites.push(write); };

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
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#111111' }]), deletedLegendEntries: ref([]),
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
    // The live rule paint of the displayed Result (app/candidate-render.js, the
    // compile's `rules` stage).
    const catalogAdmission = displayedFeatureAddressing(features, ['diagram.svg'], 0);
    /** @type {{ renderedId: string, color: string | null }[]} */
    let fills = [];
    measured.push(await count_.measure(() => {
      assert.equal(preparation.isPrepared(rules), true);
      fills = compileDirectEditorMutationPlan({
        catalogAdmission, manualSpecificRules: rules,
        livePreview: { drawnTracks: [], domains: ['featureFills'], paletteColors: { CDS: '#cccccc' }, drawnContext: null }
      }).operationsByResult[0].featureFills;
    }));
    assert.deepEqual(fills.at(-1), { renderedId: hashOf(count - 1), color: '#111111' });
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
    // Not a palette row (D-15): the choice writes one rule per feature.
    const featureStyleScopeDialog = {
      show: true, feat: features[0], color: '#abcdef', legendName: 'CDS', defaultColorType: null
    };
    Object.assign(state, {
      results: ref([]), selectedResultIndex: ref(0), svgContainer: ref(null), clickedFeature: ref(null),
      featureStyleScopeDialog, resetColorDialog: {}, legendRenameDialog: {}, originalLegendOrder: ref([]),
      originalLegendColors: ref({}), originalSvgStroke: ref({ color: null, width: null }), appliedPaletteColors: ref({}),
      skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false), hasPendingPaletteDraft: ref(false)
    });
    const caption = () => 'CDS';
    const actions = createFeatureColorActions({ readDrawnTracks: () => [], writeLegendRowColor: recordLegendRowColor,
      state, nextTick: async () => {}, onLegendGeometryChanged: () => {}, extractLegendEntries: () => {},
      getFeatureElements: () => [], getFeatureFillElements: () => [],
      readUserDefaultColor: () => null, setDefaultColor: () => assert.fail('the row is drawn by rules'),
      ruleActions: {
        runWithRuleMatches: (_, commit) => commit(),
        commitSpecificRules: async (build) => { committed = build(rules, (row, patch) => ({ ...row, ...patch })); return true; },
        countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null,
        findFeaturesWithSameDisplayedLabel: () => [], findFeaturesWithSameIndividualLabel: () => [],
        findFeaturesWithSameLegendItem: (feature) => features.filter((other) => other !== feature),
        findMatchingRegexRule: () => null, getDisplayedFeatureLabel: () => '', getIndividualFeatureLabel: () => '',
        getEffectiveLegendCaption: caption, effectiveLegendCaptions: () => caption,
        getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
        getLabelSpecificRule: () => null, getLegendRowRules: () => []
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
      ports: { requestAutomaticRerender: () => true, isRedrawOwed: () => false }, getCommittedRequest: () => null, isPatternEditAvailable: () => true
    });
    /** @type {Record<string, any>[]} */
    let siblings = [];
    measured.push(await count_.measure(() => { siblings = actions.findFeaturesWithSameLegendItem(features[0]); }));
    assert.equal(siblings.length, count - 1);
  }
  linear(...measured);
});

// Step 3c A1 (CW-06 row M1): the rule preparation builds the union of the
// catalog and biological features once per pair of array identities. A
// Generate that reuses the catalog (3a) keeps both arrays, so `prepare`,
// `isPrepared`, and the visibility preparation read the union already built.
test('the rule preparation builds the feature union once per pair of array identities', async () => {
  const iterations = new Map();
  // An array whose iterations are counted by name: the union spreads it.
  const counted = (name, rows) => {
    const array = [...rows];
    iterations.set(name, 0);
    Object.defineProperty(array, Symbol.iterator, {
      value() { iterations.set(name, iterations.get(name) + 1); return Array.prototype[Symbol.iterator].call(this); }
    });
    return array;
  };
  const row = (index) => ({ id: `feature-${index}`, svg_id: hashOf(index), type: 'CDS', qualifiers: { product: [`p${index}`] } });
  const rows = [row(0), row(1), row(2)];
  const rules = rows.map((_, index) => ({ feat: 'CDS', qual: 'hash', val: hashOf(index), color: '#111111', cap: 'CDS' }));
  const state = withDrawings({
    extractedFeatures: ref(counted('catalog', rows)), biologicalFeatures: ref(counted('biology', rows.slice(0, 2))),
    manualSpecificRules: rules, svgResultIdentity: ref('one'), legendEntries: ref([]), deletedLegendEntries: ref([]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set())
  });
  /** @type {string[][]} */
  const evaluated = [];
  const evaluate = async (payload) => {
    evaluated.push(payload.features.map((feature) => feature.selector.hash));
    return evaluateHashRules(payload);
  };
  const preparation = createRulePreparation({ state, evaluate });
  const reads = () => Object.fromEntries(iterations);

  assert.equal(await preparation.prepare(rules), true);
  assert.equal(await preparation.prepare(rules), true);
  assert.equal(preparation.isPrepared(rules), true);
  assert.deepEqual(reads(), { catalog: 1, biology: 1 }, 'unchanged arrays: one union for every call');

  // A new catalog array with equal rows (identity, not content) builds again.
  state.extractedFeatures.value = counted('catalog', rows);
  assert.equal(preparation.isPrepared(rules), true);
  assert.equal(preparation.isPrepared(rules), true);
  assert.deepEqual(reads(), { catalog: 1, biology: 2 }, 'a replaced catalog array builds the union again, once');

  // The biological array alone replaced builds again.
  state.biologicalFeatures.value = counted('biology', rows.slice(0, 1));
  assert.equal(preparation.isPrepared(rules), true);
  assert.deepEqual(reads(), { catalog: 2, biology: 1 }, 'a replaced biological array builds the union again');

  // A Session Load replaces both arrays: the next call reads only the new rows.
  const loaded = [row(7), row(8)];
  state.extractedFeatures.value = counted('catalog', loaded);
  state.biologicalFeatures.value = counted('biology', []);
  evaluated.length = 0;
  assert.equal(await preparation.prepare([{ feat: 'CDS', qual: 'hash', val: hashOf(7), color: '#222222', cap: 'CDS' }]), true);
  assert.deepEqual(evaluated, [[hashOf(7), hashOf(8)]], 'the next preparation evaluates the loaded features only');
  assert.deepEqual(reads(), { catalog: 1, biology: 1 });
});
