import { evaluatePythonRules, groupColorRulesOver, synthesizePythonColorRules } from './helpers/python-evaluator.mjs';
import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';

// The Legend owner's writer of a row's own color (app/feature-editor.js wires it); these tests record its calls.
const legendRowColorWrites = [];
/** @param {unknown[]} write */
const recordLegendRowColor = (...write) => { legendRowColorWrites.push(write); };

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-color-actions-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await cp(sourceDir, tempDir, { recursive: true });
const colorActionsSource = await readFile(join(sourceDir, 'app', 'feature-editor', 'color-actions.js'), 'utf8');

const { LIVE_EDIT_DOMAINS } = await import(pathToFileURL(join(tempDir, 'app', 'candidate-render.js')));
const { createFeatureColorActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'color-actions.js'))
);
const { getFeatureColorRuleHash } = await import(pathToFileURL(join(tempDir, 'services', 'feature-utils.js')));
const { buildLegendIntents, colorRuleLegendGroups, legendRowRules } = await import(pathToFileURL(join(tempDir, 'services', 'specific-color-rules.js')));

assert.doesNotMatch(colorActionsSource, /serializeCleanSvg|results\.value\[[^\]]+\]\s*=/);

const ref = (value) => ({ value });
// The mounted SVG as `getAllFeatureLegendGroups` (services/legend-svg.js) reads it:
// `#legend` holding one `#feature_legend` group, or no legend.
const legendSvg = (featureLegend = null) => ({
  getElementById: (id) => (id === 'legend' && featureLegend
    ? { querySelector: (selector) => (selector === '#feature_legend' ? featureLegend : null) }
    : null)
});

// Each catalog feature carries its catalog identity (`record_key`,
// `biological_feature_id`), from which Python composes its feature instance.
const featureA = {
  id: 'feature-a',
  svg_id: 'hash-a',
  record_key: 'record-1',
  biological_feature_id: 'hash-a',
  type: 'CDS',
  qualifiers: { gene_kind: 'core biosynthetic genes' },
  start: 1,
  end: 10
};
const featureB = {
  id: 'feature-b',
  svg_id: 'hash-b',
  record_key: 'record-1',
  biological_feature_id: 'hash-b',
  type: 'CDS',
  qualifiers: { gene_kind: 'core biosynthetic genes' },
  start: 11,
  end: 20
};
const hashOnlyFeature = {
  id: 'feature-c',
  svg_id: 'hash-c',
  record_key: 'record-1',
  biological_feature_id: 'hash-c',
  type: 'CDS',
  qualifiers: { gene_kind: 'transport' },
  start: 21,
  end: 30
};

const specificRule = {
  feat: 'CDS',
  qual: 'gene_kind',
  val: '^core biosynthetic genes$',
  color: '#111111',
  cap: 'Core'
};
const manualSpecificRules = [
  specificRule,
  { feat: 'CDS', qual: 'hash', val: 'hash-a', color: '#222222', cap: 'Core' },
  { feat: 'CDS', qual: 'hash', val: 'hash-b', color: '#222222', cap: 'Core' },
  { feat: 'CDS', qual: 'hash', val: 'hash-c', color: '#333333', cap: 'Core' },
  { feat: 'CDS', qual: 'hash', val: 'hash-z', color: '#444444', cap: 'Other' }
];
const featureColorOverrides = {};
const legendColorOverrides = {};
const committedLegendIntents = [];
const extractedFeatures = ref([featureA, featureB, hashOnlyFeature]);
const biologicalFeatures = ref([featureA, featureB, hashOnlyFeature]);
const legendEntries = ref([{ caption: 'Core', color: '#111111', featureIds: ['hash-a', 'hash-b', 'hash-c'] }]);
const featureStyleScopeDialog = {
  show: true,
  kind: 'fill',
  feat: featureA,
  color: '#abcdef',
  strokeColor: null,
  strokeWidth: null,
  matchingRule: specificRule,
  ruleMatchCount: 2,
  legendName: 'Core',
  siblingCount: 2,
  displayLabel: null,
  displayLabelSiblingCount: 0,
  annotationLabel: null,
  annotationLabelSiblingCount: 0,
  existingCaptionRule: null,
  existingCaptionColor: null,
  resolve: null
};

let applySpecificRulesCount = 0;
/** @type {Array<{ domains: readonly string[] }>} */
const shownIntents = [];
const previewFillColors = new Map();
const featureElementsById = new Map();
let previewFillApplyCount = 0;
const svgContainer = ref(null);
const clickedFeature = ref(null);
const featureStrokeOverrides = {};
const legendStrokeOverrides = {};
// The rule commit shows its fills itself (svg-styles); the color action
// commits only its own DOM edits.
const applyRulePreviewFill = (featureId, color) => {
  if (previewFillColors.get(featureId) === color) return;
  previewFillColors.set(featureId, color);
  previewFillApplyCount += 1;
};

const { featureOverrideKey } = await import(pathToFileURL(join(tempDir, 'services', 'feature-override-identity.js')));
const { createRulePreparation, ruleFeaturePayload, runWhenPrepared } = await import(pathToFileURL(join(tempDir, 'app', 'rule-matching.js')));
const { ruleMatcher } = await import(pathToFileURL(join(tempDir, 'services', 'rule-matchers.js')));
const firstMatchingRule = (feature, rules) => ruleMatcher(rules).first(feature);
// The rule owner's `runWithRuleMatches` for these fakes: the color matches of the rules, prepared once.
const runWithRuleMatchesOf = (preparation, state) => (rules, commit) => runWhenPrepared(state, () => [preparation.prepare(rules)], commit);
const preparationState = { extractedFeatures, biologicalFeatures, manualSpecificRules };
let legendSiblings = [featureB, hashOnlyFeature];
let matchingRegexRule = specificRule;
const preparedRuleSets = [];
const actionsPreparation = createRulePreparation({ state: withDrawings(preparationState), evaluate: evaluatePythonRules, synthesize: synthesizePythonColorRules });
const runWithRuleMatchesForActions = runWithRuleMatchesOf(actionsPreparation, preparationState);
const actions = createFeatureColorActions({ readDrawnTracks: () => [], writeLegendRowColor: recordLegendRowColor,
  state: withDrawings({
    results: ref([]),
    selectedResultIndex: ref(0),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    manualSpecificRules,
    extractedFeatures,
    biologicalFeatures,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    featureStyleScopeDialog,
    resetColorDialog: {},
    legendRenameDialog: {},
    legendEntries,
    deletedLegendEntries: ref([]),
    dormantLegendEntries: ref([]),
    legendStrokeOverrides,
    legendColorOverrides,
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    featureStrokeOverrides,
    skipCaptureBaseConfig: ref(false),
    skipExtractOnSvgChange: ref(false),
    addedLegendCaptions: ref(new Set())
  }),
  nextTick: async () => {},
  showEditorIntent: (options) => { shownIntents.push(options); },
  ruleActions: {
    runWithRuleMatches: (rules, commit) => {
      preparedRuleSets.push(rules.map((rule) => `${rule.feat}|${rule.qual}|${rule.val}`));
      return runWithRuleMatchesForActions(rules, commit);
    },
    commitSpecificRules: async (build, _label, {afterCommit = () => {}, previousLegendIntents = []} = {}) => {
      const rules = await build(manualSpecificRules, (row, patch) => ({ ...row, ...patch }));
      if (!rules) return false;
      committedLegendIntents.push(previousLegendIntents);
      const preparation = createRulePreparation({state:withDrawings({extractedFeatures,biologicalFeatures,manualSpecificRules}),evaluate:evaluatePythonRules});
      const candidate = await preparation.prepareCandidate(rules);
      if (!candidate) return false;
      const drawn = new Set(extractedFeatures.value.map(feature => firstMatchingRule(feature, candidate.rules)));
      const { intents } = buildLegendIntents(candidate.rules, [], rule => drawn.has(rule));
      manualSpecificRules.splice(0,manualSpecificRules.length,...candidate.rules);
      legendEntries.value = intents.map(entry=>({...entry}));
      for (const feature of extractedFeatures.value) {
        const rule=firstMatchingRule(feature,manualSpecificRules);
        if(rule) featureColorOverrides[featureOverrideKey(feature)]={color:rule.color,caption:rule.cap};
        else delete featureColorOverrides[featureOverrideKey(feature)];
      }
      afterCommit(intents);
      // As `commitOnce`: the fills show after `afterCommit`, at the end of the step.
      for (const feature of extractedFeatures.value) {
        const rule=firstMatchingRule(feature,manualSpecificRules);
        if(rule) applyRulePreviewFill(feature.svg_id, rule.color);
      }
      applySpecificRulesCount++;
      return true;
    },
    countFeaturesMatchingRule: () => 0,
    findExistingColorForCaption: () => null,
    findFeaturesWithSameDisplayedLabel: (currentFeature, label) => extractedFeatures.value.filter(
      (feature) =>
        feature.svg_id !== currentFeature.svg_id &&
        (feature.displayLabel || feature.product) === label
    ),
    findFeaturesWithSameIndividualLabel: () => [],
    findFeaturesWithSameLegendItem: () => legendSiblings,
    findMatchingRegexRule: () => matchingRegexRule,
    getDisplayedFeatureLabel: (feature) => feature.displayLabel || feature.product || '',
    effectiveLegendCaptions: () => () => 'Core',
    // The port's answer for this table, whose row `Core` holds every rule captioned `Core` (the key
    // allocation, which gives one caption in several colors a row per color, has its own fixture).
    getLegendRowRules: (caption) => manualSpecificRules.filter((rule) => rule.cap === caption),
    getIndividualFeatureLabel: (feature) => feature.product || '',
    // FE-09 (rule-actions.js): a feature's single-feature condition is its stable hash.
    getFeatureQualifier: (feature) => ({ qual: 'hash', val: getFeatureColorRuleHash(feature) }),
    groupColorRules: groupColorRulesOver(actionsPreparation, (feature) => [feature.displayLabel || feature.product || '', feature.product || ''])
  },
  getFeatureElements: (_svg, featureId) => featureElementsById.get(featureId) || [],
  getFeatureFillElements: (_svg, featureId) => featureElementsById.get(featureId) || []
});

await actions.handleColorScopeChoice('caption');

assert.equal(applySpecificRulesCount, 1);
assert.equal(manualSpecificRules.find(rule => rule.qual === 'gene_kind').color, '#abcdef');
assert.equal(legendEntries.value[0].color, '#abcdef');
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-a'), false);
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-b'), false);
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-c')?.color, '#abcdef');
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-z')?.color, '#444444');
assert.deepEqual(featureColorOverrides[featureOverrideKey(featureA)], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides[featureOverrideKey(featureB)], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides[featureOverrideKey(hashOnlyFeature)], { color: '#abcdef', caption: 'Core' });

// Replacing an inherited caption requires every existing contributor.
manualSpecificRules.splice(0);
legendEntries.value = [{ caption: 'Core', color: '#111111' }];
Object.assign(featureStyleScopeDialog, { show: true, feat: featureA, color: '#abcdef', legendName: 'Core' });
await actions.handleColorScopeChoice('caption');
assert.deepEqual(committedLegendIntents.at(-1), [{ caption: 'Core', color: '#111111' }]);
assert.equal(legendColorOverrides.Core, '#abcdef');

// A partial group or a new caption cannot adopt an unrelated existing row.
manualSpecificRules.splice(0);
legendEntries.value = [{ caption: 'Core', color: '#111111' }];
delete legendColorOverrides.Core;
await actions.setFeatureColor(featureA, '#abcdef', 'Core');
assert.deepEqual(committedLegendIntents.at(-1), []);
assert.equal(legendColorOverrides.Core, undefined);
manualSpecificRules.splice(0);
legendEntries.value = [{ caption: 'Manual row', color: '#111111' }];
await actions.setFeatureColor(featureA, '#abcdef', 'Manual row');
assert.deepEqual(committedLegendIntents.at(-1), []);
assert.equal(legendColorOverrides['Manual row'], undefined);

// A popup color change prepares only the clicked feature's rules: the sibling
// features of its legend item add no rules, so the request stays constant
// (it was siblings x features on a large Session).
{
  manualSpecificRules.splice(0);
  legendSiblings = Array.from({ length: 50 }, (_, index) => ({
    id: `sibling-${index}`,
    svg_id: `hash-sibling-${index}`,
    type: 'CDS',
    product: `product ${index}`,
    qualifiers: { product: [`product ${index}`] },
    start: 100 + index,
    end: 110 + index
  }));
  preparedRuleSets.length = 0;
  clickedFeature.value = { feat: featureA, legendName: 'Core' };
  await actions.updateClickedFeatureColor('#123456');
  const distinctRuleKeys = new Set(preparedRuleSets.flat());
  assert.ok(distinctRuleKeys.size <= 3, `prepared ${distinctRuleKeys.size} distinct rules for 50 siblings`);
  legendSiblings = [featureB, hashOnlyFeature];
  clickedFeature.value = null;
  Object.assign(featureStyleScopeDialog, { show: false });
  manualSpecificRules.splice(0);
}

// Each popup action prepares an allowlisted set of rules (OV-198): a stroke
// edit reads only the saved rules, a color edit also the clicked feature's hash
// rule, once it commits; the scope dialog reads only the saved rules (OV-225).
// A group edit's other rules are Python's (OV-244) and prepared at its commit.
{
  const popupFeature = {
    id: 'popup-feature', svg_id: 'hash-popup', type: 'CDS', product: 'popup product',
    qualifiers: { product: ['popup product'] }, start: 40, end: 50
  };
  const ruleKey = (rule) => `${rule.feat}|${rule.qual}|${rule.val}`;
  const savedRules = [specificRule, { feat: 'CDS', qual: 'hash', val: 'hash-z', color: '#444444', cap: 'Other' }];
  const savedKeys = savedRules.map(ruleKey);
  const clickedKeys = [...savedKeys, `CDS|hash|${getFeatureColorRuleHash(popupFeature)}`];
  const strokeScope = (kind) => () => {
    Object.assign(featureStyleScopeDialog, { show: true, kind, feat: popupFeature, color: '#123456', strokeColor: '#123456' });
    return actions.handleFeatureStyleScopeChoice('cancel');
  };
  const allowlist = [
    ['applyStrokeToSelectedFeatures', () => actions.applyStrokeToSelectedFeatures([popupFeature], '#112233', 2), savedKeys],
    ['resetClickedFeatureStroke', () => actions.resetClickedFeatureStroke(), savedKeys],
    ['setClickedFeatureStrokeColorValue', () => actions.setClickedFeatureStrokeColorValue('#112233'), savedKeys],
    ['setClickedFeatureStrokeWidthValue', () => actions.setClickedFeatureStrokeWidthValue(3), savedKeys],
    ['updateClickedFeatureStroke', () => actions.updateClickedFeatureStroke('#112233', 2), savedKeys],
    ['handleFeatureStyleScopeChoice (stroke)', strokeScope('stroke'), savedKeys],
    ['handleFeatureStyleScopeChoice (fill)', strokeScope('fill'), clickedKeys],
    ['updateClickedFeatureColor (scope dialog)', () => actions.updateClickedFeatureColor('#123456'), savedKeys],
    ['updateClickedFeatureColor (no dialog)', async () => {
      matchingRegexRule = null;
      legendSiblings = [];
      try { await actions.updateClickedFeatureColor('#123456'); } finally {
        matchingRegexRule = specificRule;
        legendSiblings = [featureB, hashOnlyFeature];
      }
    }, clickedKeys],
    // Review r1: Selected features writes Python's table, which Python verified,
    // so it prepares only the saved rules before the synthesis, not a hash
    // rule per selected feature that it no longer writes.
    ['applyColorToSelectedFeatures', async () => {
      await actions.applyColorToSelectedFeatures([featureA, featureB], '#13579b', 'Picked');
      assert.deepEqual([featureA, featureB, hashOnlyFeature].map((feature) => firstMatchingRule(feature, manualSpecificRules)?.color),
        ['#13579b', '#13579b', undefined], 'the selected features take the color, and no other feature does');
    }, savedKeys]
  ];
  for (const [name, run, expected] of allowlist) {
    manualSpecificRules.splice(0, manualSpecificRules.length, ...savedRules);
    clickedFeature.value = { feat: popupFeature, svg_id: popupFeature.svg_id, legendName: 'Core', strokeColor: '#000000', strokeWidth: 1 };
    Object.assign(featureStyleScopeDialog, { show: false, kind: 'fill', feat: null });
    preparedRuleSets.length = 0;
    await run();
    assert.deepEqual([...new Set(preparedRuleSets.flat())].sort(), [...expected].sort(), `${name} prepares exactly its allowlisted rules`);
  }
  clickedFeature.value = null;
  Object.assign(featureStyleScopeDialog, { show: false, kind: 'fill', feat: null });
  manualSpecificRules.splice(0);
}

const labelFeatureA = {
  id: 'label-feature-a',
  svg_id: 'f11111111_record_1',
  record_key: 'record-1',
  biological_feature_id: 'f11111111',
  rendered_feature_svg_id: 'f11111111_record_1',
  stable_svg_id: 'faaaaaaaa',
  type: 'CDS',
  product: 'wsv360-like protein',
  qualifiers: { product: ['wsv360-like protein'] },
  selector: { hash: 'f11111111', qualifiers: { product: ['wsv360-like protein'] } }
};
const labelFeatureB = {
  id: 'label-feature-b',
  svg_id: 'f22222222_record_2',
  record_key: 'record-2',
  biological_feature_id: 'f22222222',
  rendered_feature_svg_id: 'f22222222_record_2',
  stable_svg_id: 'fbbbbbbbb',
  type: 'CDS',
  product: 'wsv360-like protein',
  qualifiers: { product: ['wsv360-like protein'] },
  selector: { hash: 'f22222222', qualifiers: { product: ['wsv360-like protein'] } }
};

manualSpecificRules.splice(
  0,
  manualSpecificRules.length,
  {
    feat: 'CDS',
    qual: 'product',
    val: '^wsv.*$',
    color: '#999999',
    cap: 'broad product rule'
  },
  {
    feat: 'CDS',
    qual: 'product',
    val: '^wsv360-like protein$',
    color: '#000000',
    cap: 'old caption',
    fromFile: true
  }
);
legendEntries.value = [];
Object.keys(featureColorOverrides).forEach((key) => delete featureColorOverrides[key]);
extractedFeatures.value = [labelFeatureA, labelFeatureB];
biologicalFeatures.value = [labelFeatureA, labelFeatureB];
Object.assign(featureStyleScopeDialog, {
  show: true,
  feat: labelFeatureA,
  color: '#8cf04f',
  matchingRule: null,
  ruleMatchCount: 0,
  legendName: 'CDS',
  siblingCount: 0,
  displayLabel: 'wsv360-like protein',
  displayLabelSiblingCount: 1,
  annotationLabel: 'wsv360-like protein',
  annotationLabelSiblingCount: 1,
  existingCaptionColor: null
});

await actions.handleColorScopeChoice('displayLabel');

assert.deepEqual(manualSpecificRules, [
  {
    feat: 'CDS',
    qual: 'product',
    val: '^wsv360-like protein$',
    color: '#8cf04f',
    cap: 'wsv360-like protein'
  },
  {
    feat: 'CDS',
    qual: 'product',
    val: '^wsv.*$',
    color: '#999999',
    cap: 'broad product rule'
  }
]);
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash'), false);
assert.equal(Object.hasOwn(manualSpecificRules[0], 'fromFile'), false);
assert.deepEqual(featureColorOverrides[featureOverrideKey(labelFeatureA)], {
  color: '#8cf04f',
  caption: 'wsv360-like protein'
});
assert.deepEqual(featureColorOverrides[featureOverrideKey(labelFeatureB)], {
  color: '#8cf04f',
  caption: 'wsv360-like protein'
});

manualSpecificRules.splice(0);
Object.keys(featureColorOverrides).forEach((key) => delete featureColorOverrides[key]);
extractedFeatures.value = [labelFeatureA];
biologicalFeatures.value = [labelFeatureA];
await actions.setFeatureColor(labelFeatureA, '#123456', 'single feature');
assert.deepEqual(manualSpecificRules, [{
  feat: 'CDS',
  qual: 'hash',
  // "This feature only" names the feature instance (PD-OI-069 revision 2).
  val: '^f11111111@record-1$',
  color: '#123456',
  cap: 'single feature'
}]);
assert.equal(getFeatureColorRuleHash(labelFeatureA), 'f11111111');

legendEntries.value = [{ caption: 'single feature', color: '#123456', featureIds: ['f11111111_record_1'] }];
const noOpFillCount = previewFillApplyCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#123456', 'single feature'), false);
assert.equal(previewFillApplyCount, noOpFillCount);
// Review L-3 (D-37): the row a name names is found through the Legend owner's
// key port (the row of that name in the color first), never by reading the
// name as a key: a row whose key is the name in another color is not it.
legendEntries.value = [
  { caption: 'single feature', originalCaption: 'single feature', color: '#ff0000', featureIds: [] },
  { caption: 'single feature', originalCaption: 'single feature [#123456]', color: '#123456', featureIds: ['f11111111_record_1'] }
];
assert.equal(await actions.setFeatureColor(labelFeatureA, '#123456', 'single feature'), false, 'the same color again is no edit');
assert.equal(previewFillApplyCount, noOpFillCount);
legendEntries.value = [{ caption: 'single feature', color: '#123456', featureIds: ['f11111111_record_1'] }];

const compoundFillCount = previewFillApplyCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#654321', 'renamed feature'), true);
assert.ok(previewFillApplyCount > compoundFillCount);
assert.equal(manualSpecificRules[0].cap, 'renamed feature');
assert.equal(legendEntries.value[0].caption, 'renamed feature');
legendEntries.value = [];

const outsideLabelGroup = {
  ...labelFeatureB,
  id: 'label-feature-outside',
  svg_id: 'f33333333_record_2',
  rendered_feature_svg_id: 'f33333333_record_2',
  biological_feature_id: 'f33333333',
  // Its own feature: a feature that shares a hash is one target for Python (OV-244).
  selector: { ...labelFeatureB.selector, hash: 'f33333333' },
  displayLabel: 'a different edited label'
};
manualSpecificRules.splice(0);
Object.keys(featureColorOverrides).forEach((key) => delete featureColorOverrides[key]);
extractedFeatures.value = [labelFeatureA, labelFeatureB];
biologicalFeatures.value = [labelFeatureA, labelFeatureB, outsideLabelGroup];
Object.assign(featureStyleScopeDialog, {
  show: true,
  feat: labelFeatureA,
  color: '#654321',
  displayLabel: 'wsv360-like protein',
  displayLabelSiblingCount: 1
});

await actions.handleColorScopeChoice('displayLabel');

assert.deepEqual(
  manualSpecificRules.map(({ feat, qual, val, color }) => ({ feat, qual, val, color })),
  [
    { feat: 'CDS', qual: 'hash', val: '^f11111111@record-1$', color: '#654321' },
    { feat: 'CDS', qual: 'hash', val: '^f22222222@record-2$', color: '#654321' }
  ]
);
assert.equal(featureColorOverrides[featureOverrideKey(outsideLabelGroup)], undefined);

const conflictingFeatureA = {
  ...labelFeatureA,
  qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] },
  selector: {
    hash: 'f11111111',
    record_location: 'RecA:0..90:+',
    qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] }
  }
};
const conflictingFeatureB = {
  ...labelFeatureB,
  qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] },
  selector: {
    hash: 'f22222222',
    record_location: 'RecB:0..90:+',
    qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] }
  }
};
manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'CDS',
  qual: 'gene',
  val: '^wsv360$',
  color: '#101010',
  cap: 'existing gene rule'
});
extractedFeatures.value = [conflictingFeatureA, conflictingFeatureB];
biologicalFeatures.value = [conflictingFeatureA, conflictingFeatureB];
Object.assign(featureStyleScopeDialog, {
  show: true,
  feat: conflictingFeatureA,
  color: '#abcdef',
  displayLabel: 'wsv360-like protein',
  displayLabelSiblingCount: 1
});
await actions.handleColorScopeChoice('displayLabel');
// Python's matcher reads the product before the gene of these features, so the
// product rule colors both and no other feature changes (OV-244, Owner R-07).
assert.deepEqual(manualSpecificRules.map(({ qual, color }) => [qual, color]), [['gene', '#101010'], ['product', '#abcdef']]);

manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'CDS',
  qual: 'record_location',
  val: '^RecA:0\\.\\.90:\\+$',
  color: '#202020',
  cap: 'existing record rule'
});
Object.assign(featureStyleScopeDialog, { show: true, feat: conflictingFeatureA, color: '#aabbcc' });
await actions.handleColorScopeChoice('displayLabel');
// The record rule wins over a product rule for feature A, which takes a hash rule.
assert.deepEqual(manualSpecificRules.map(({ qual, val, color }) => [qual, qual === 'hash' ? val : '', color]), [
  ['record_location', '', '#202020'], ['product', '', '#aabbcc'], ['hash', '^f11111111@record-1$', '#aabbcc']
]);

const geneLabelFeature = {
  ...labelFeatureB,
  product: '',
  gene: 'wsv360-like protein',
  displayLabel: 'wsv360-like protein',
  qualifiers: { gene: ['wsv360-like protein'] },
  selector: { hash: 'f22222222', qualifiers: { gene: ['wsv360-like protein'] } }
};
manualSpecificRules.splice(0);
extractedFeatures.value = [labelFeatureA, geneLabelFeature];
biologicalFeatures.value = [labelFeatureA, geneLabelFeature];
Object.assign(featureStyleScopeDialog, {
  show: true,
  feat: labelFeatureA,
  color: '#fedcba',
  displayLabel: 'wsv360-like protein',
  displayLabelSiblingCount: 1
});
await actions.handleColorScopeChoice('displayLabel');
// The label comes from another qualifier on each: one exact rule each (OV-244).
assert.deepEqual(manualSpecificRules.map(({ qual, val }) => [qual, val]).sort(),
  [['gene', '^wsv360-like protein$'], ['product', '^wsv360-like protein$']]);

const duplicateFeatureA = {
  ...labelFeatureA,
  id: 'duplicate-a',
  svg_id: 'f44444444_record_1',
  rendered_feature_svg_id: 'f44444444_record_1',
  product: '',
  qualifiers: {},
  // The catalog identity of the copy in the first Linear card.
  record_key: 'linear-seq-a',
  biological_feature_id: 'f44444444',
  selector: { hash: 'f44444444', qualifiers: {} }
};
const duplicateFeatureB = {
  ...duplicateFeatureA,
  id: 'duplicate-b',
  svg_id: 'f44444444_record_2',
  rendered_feature_svg_id: 'f44444444_record_2',
  record_key: 'linear-seq-b'
};
const sharedDuplicateRule = { feat: 'CDS', qual: 'hash', val: 'f44444444', color: '#999999', cap: 'shared duplicate rule' };
const duplicateRules = () => manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap }));
const duplicateWinners = async () => (await evaluatePythonRules({
  features: [duplicateFeatureA, duplicateFeatureB].map((feature) => ruleFeaturePayload(feature)),
  rules: manualSpecificRules,
  kind: 'color'
})).winners;
extractedFeatures.value = [duplicateFeatureA, duplicateFeatureB];
biologicalFeatures.value = [duplicateFeatureA, duplicateFeatureB];
// The same file drawn twice, with a rule an older "This feature only" wrote: it
// names the stable hash both copies share. Resetting one copy keeps the other
// copy's color and Legend row: the rule names that copy's instance instead
// (PD-OI-069 revision 2).
manualSpecificRules.splice(0, manualSpecificRules.length, { ...sharedDuplicateRule });
await actions.setFeatureColorValue(duplicateFeatureA, null);
assert.deepEqual(duplicateRules(), [{ ...sharedDuplicateRule, val: '^f44444444@linear-seq-b$' }]);
assert.deepEqual(await duplicateWinners(), [-1, 0]);
// "This feature only" names the clicked copy's instance, before the
// stable-hash rule both copies share, so the other copy keeps that rule's color.
manualSpecificRules.splice(0, manualSpecificRules.length, { ...sharedDuplicateRule });
await actions.setFeatureColor(duplicateFeatureA, '#112233', 'one duplicate');
const ownDuplicateRule = { feat: 'CDS', qual: 'hash', val: '^f44444444@linear-seq-a$', color: '#112233', cap: 'one duplicate' };
assert.deepEqual(duplicateRules(), [ownDuplicateRule, sharedDuplicateRule]);
assert.deepEqual(await duplicateWinners(), [0, 1]);
// The other copy's reset takes back what names it. The shared rule colors no
// other copy (the clicked copy's own rule wins), so it leaves no row behind.
await actions.setFeatureColorValue(duplicateFeatureB, null);
assert.deepEqual(duplicateRules(), [ownDuplicateRule]);
assert.deepEqual(await duplicateWinners(), [0, -1]);

const sharedHashCds = {
  ...labelFeatureA,
  id: 'shared-hash-cds',
  svg_id: 'f55555555',
  rendered_feature_svg_id: 'f55555555',
  product: ''
};
manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'tRNA',
  qual: 'hash',
  val: 'f55555555',
  color: '#aaaaaa',
  cap: 'tRNA rule'
});
extractedFeatures.value = [sharedHashCds];
biologicalFeatures.value = [sharedHashCds];
await actions.setFeatureColor(sharedHashCds, '#445566', 'CDS rule');
assert.equal(manualSpecificRules.find((rule) => rule.feat === 'tRNA')?.color, '#aaaaaa');
assert.equal(manualSpecificRules.find((rule) => rule.feat === 'CDS')?.color, '#445566');

manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'CDS',
  qual: 'hash',
  val: '^f.*$',
  color: '#999999',
  cap: 'broad hash rule'
});
extractedFeatures.value = [labelFeatureA];
biologicalFeatures.value = [labelFeatureA];
await actions.setFeatureColor(labelFeatureA, '#778899', 'exact hash rule');
assert.equal(manualSpecificRules[0].val, '^f11111111@record-1$');
assert.equal(manualSpecificRules[1].val, '^f.*$');

const stableColorFeature = {
  ...labelFeatureA,
  id: 'legacy-color-id',
  record_key: 'record-key-a',
  biological_feature_id: 'biological-a',
  recordKey: 'record-key-a',
  biologicalFeatureId: 'biological-a'
};
const stableColorKey = 'record-key-a\u0000biological-a';
extractedFeatures.value = [stableColorFeature];
biologicalFeatures.value = [stableColorFeature];
manualSpecificRules.splice(0);
featureColorOverrides[stableColorKey] = {
  color: '#123456',
  caption: 'Stable feature'
};
await actions.setFeatureColorValue(stableColorFeature, null);
assert.equal(featureColorOverrides[stableColorKey], undefined);

await actions.setFeatureColorValue(stableColorFeature, 'none', 'No fill');
assert.deepEqual(featureColorOverrides[stableColorKey], {
  color: 'none',
  caption: 'No fill'
});
assert.equal(manualSpecificRules[0].color, 'none');

const strokeAttributes = new Map([
  ['stroke', '#111111'],
  ['stroke-width', '1']
]);
let strokeMutationCount = 0;
const strokeElement = {
  getAttribute: (name) => strokeAttributes.has(name) ? strokeAttributes.get(name) : null,
  setAttribute: (name, value) => {
    strokeAttributes.set(name, String(value));
    strokeMutationCount += 1;
  },
  removeAttribute: (name) => {
    strokeAttributes.delete(name);
    strokeMutationCount += 1;
  }
};
const strokeFeature = {
  id: 'stroke-feature',
  svg_id: 'stroke-feature-svg',
  type: 'CDS',
  product: 'Stroke feature'
};
const strokeSvg = legendSvg();
svgContainer.value = { querySelector: (selector) => selector === 'svg' ? strokeSvg : null };
featureElementsById.set(strokeFeature.svg_id, [strokeElement]);
clickedFeature.value = {
  svg_id: strokeFeature.svg_id,
  feat: strokeFeature,
  color: '#cccccc',
  strokeColor: '#111111',
  strokeWidth: 1
};
// A stroke edit writes the intent only, with the stroke Python drew; the
// composition root shows it through the executor (EU U2a).
const strokeKey = featureOverrideKey(strokeFeature) || strokeFeature.svg_id;
assert.equal(await actions.updateClickedFeatureStroke('#111111', 1), false);
assert.equal(await actions.updateClickedFeatureStroke('#222222', 2), true);
assert.deepEqual(featureStrokeOverrides[strokeKey], {
  originalStrokeColor: '#111111', originalStrokeWidth: 1, strokeColor: '#222222', strokeWidth: 2
});
assert.equal(strokeMutationCount, 0);
assert.equal(await actions.resetClickedFeatureStroke(), true);
assert.equal(featureStrokeOverrides[strokeKey], undefined);
assert.deepEqual([clickedFeature.value.strokeColor, clickedFeature.value.strokeWidth], ['#111111', 1],
  'the popup shows the stroke Python drew');
assert.equal(await actions.resetClickedFeatureStroke(), false);
assert.equal(await actions.applyStrokeToSelectedFeatures([strokeFeature], '#111111', 1), false);
[featureB, hashOnlyFeature].forEach((feature) => {
  const attributes = new Map([
    ['stroke', '#111111'],
    ['stroke-width', '1']
  ]);
  featureElementsById.set(feature.svg_id, [{
    getAttribute: (name) => attributes.has(name) ? attributes.get(name) : null,
    setAttribute: (name, value) => attributes.set(name, String(value)),
    removeAttribute: (name) => attributes.delete(name)
  }]);
});
legendEntries.value = [{
  caption: 'Core',
  color: '#cccccc',
  featureIds: [strokeFeature.svg_id, featureB.svg_id, hashOnlyFeature.svg_id]
}];

featureStyleScopeDialog.show = false;
assert.equal(await actions.setClickedFeatureStrokeWidthValue(1), false);
assert.equal(featureStyleScopeDialog.show, false);
assert.equal(await actions.setClickedFeatureStrokeWidthValue(''), false);
assert.equal(featureStyleScopeDialog.show, false);
assert.equal(await actions.setClickedFeatureStrokeWidthValue(2.5), false);
assert.equal(featureStyleScopeDialog.show, true);
assert.equal(featureStyleScopeDialog.kind, 'stroke');
assert.equal(featureStyleScopeDialog.strokeColor, '#111111');
assert.equal(featureStyleScopeDialog.strokeWidth, 2.5);
assert.equal(clickedFeature.value?.svg_id, strokeFeature.svg_id, 'OV-1278: the popup stays open while the dialog asks');
assert.equal(strokeAttributes.get('stroke'), '#111111');
assert.equal(strokeAttributes.get('stroke-width'), '1');
assert.equal(await actions.handleFeatureStyleScopeChoice('cancel'), false);
assert.equal(featureStyleScopeDialog.show, false);
assert.deepEqual([clickedFeature.value?.svg_id, clickedFeature.value?.strokeWidth], [strokeFeature.svg_id, 1],
  'OV-1278: Cancel leaves the popup open with the stroke the feature has');

clickedFeature.value = {
  svg_id: strokeFeature.svg_id,
  feat: strokeFeature,
  color: '#cccccc',
  legendName: 'Core',
  strokeColor: '#111111',
  strokeWidth: 1
};
assert.equal(await actions.setClickedFeatureStrokeColorValue('#445566'), false);
assert.equal(await actions.handleFeatureStyleScopeChoice('caption'), true);
assert.equal(strokeAttributes.get('stroke'), '#111111');
assert.deepEqual(legendStrokeOverrides.Core, {
  originalStrokeColor: null,
  originalStrokeWidth: null,
  strokeColor: '#445566',
  strokeWidth: 1
});
assert.equal(strokeMutationCount, 0);

// Rule scope uses Python's wildcard feature type and inline flags too.
const wildcardRule = { feat: '*', qual: 'gene_kind', val: '(?i)core', color: '#111111', cap: '' };
manualSpecificRules.splice(0, manualSpecificRules.length, wildcardRule);
extractedFeatures.value = [featureB, hashOnlyFeature];
biologicalFeatures.value = [featureB, hashOnlyFeature];
Object.assign(featureStyleScopeDialog, { show: true, kind: 'stroke', feat: featureB,
  matchingRule: wildcardRule, strokeColor: '#abcdef', strokeWidth: 2 });
assert.equal(await actions.handleFeatureStyleScopeChoice('rule'), true);
assert.deepEqual([featureB, hashOnlyFeature].map((feature) => (
  featureStrokeOverrides[featureOverrideKey(feature) || feature.svg_id]?.strokeColor ?? null
)), ['#abcdef', '#445566'], 'the rule scope restrokes only the features its rule matches');

clickedFeature.value = null;
featureElementsById.clear();

globalThis.CSS = { escape: (value) => String(value) };
const legendText = { textContent: 'Short caption' };
const legendPath = {
  getAttribute: (name) => name === 'fill' ? '#123456' : null,
  setAttribute: () => {}
};
const legendAttributes = new Map([['data-legend-key', 'Short caption']]);
const legendEntryGroup = {
  getAttribute: (name) => legendAttributes.get(name) || null,
  setAttribute: (name, value) => legendAttributes.set(name, value),
  querySelector: (selector) => selector === 'text' ? legendText : null,
  querySelectorAll: (selector) => selector === 'path' ? [legendPath] : []
};
const legendFeatureGroup = {
  querySelector: (selector) => selector.includes('Short caption') ? legendEntryGroup : null
};
const svg = legendSvg(legendFeatureGroup);
svgContainer.value = { querySelector: (selector) => selector === 'svg' ? svg : null };
legendEntries.value = [{ caption: 'Short caption', color: '#123456', featureIds: [] }];
extractedFeatures.value = [];
biologicalFeatures.value = [];

await actions.renameLegendEntry(0, 'Oxidative phosphorylation');

// U3a A2a (gap 3): the rename of an editor row writes the intent, its name;
// the row keeps its key (D-25), and the port shows the row with its styles in
// one compile.
assert.deepEqual(legendEntries.value[0], {
  caption: 'Oxidative phosphorylation', originalCaption: 'Short caption', color: '#123456', featureIds: []
});
assert.deepEqual(shownIntents.map(({ domains }) => [...domains]), [[
  ...LIVE_EDIT_DOMAINS.legendStructure, ...LIVE_EDIT_DOMAINS.legendFills, ...LIVE_EDIT_DOMAINS.strokes
]]);
assert.equal(legendText.textContent, 'Short caption', 'the executor renames the row');
assert.equal(legendAttributes.get('data-legend-key'), 'Short caption');

// FE-10: Reset fill uses the palette default of the feature being reset. A
// canceled or completed Reset dialog of another feature type must not leave a
// color behind for a later Reset that has no dialog.
{
  const feature = (svgId, type, product, start, end) => ({
    id: svgId, svg_id: svgId, record_key: 'record-1', biological_feature_id: svgId,
    type, product, qualifiers: { product: [product] }, start, end
  });
  const trnaA = feature('trna-a', 'tRNA', 'tRNA-Leu', 1, 70);
  const trnaB = feature('trna-b', 'tRNA', 'tRNA-Leu', 80, 150);
  const rrna = feature('rrna-s', 'rRNA', 's-rRNA', 200, 900);
  const resetRules = [{ feat: 'rRNA', qual: 'product', val: '^s-rRNA$', color: '#ff00ff', cap: 'small rRNA' }];
  const committed = [];
  const resetDialog = { show: false, caption: '', siblingCount: 0 };
  const resetClicked = ref(null);
  const resetFeatures = ref([trnaA, trnaB, rrna]);
  const resetPreparation = createRulePreparation({
    state: withDrawings({ extractedFeatures: resetFeatures, biologicalFeatures: resetFeatures, manualSpecificRules: resetRules }),
    evaluate: evaluatePythonRules,
    synthesize: synthesizePythonColorRules
  });
  const resetActions = createFeatureColorActions({ readDrawnTracks: () => [], writeLegendRowColor: recordLegendRowColor,
    state: withDrawings({
      results: ref([]),
      selectedResultIndex: ref(0),
      appliedPaletteColors: ref({ tRNA: '#e8b441', rRNA: '#71ee7d' }),
      manualSpecificRules: resetRules,
      extractedFeatures: resetFeatures,
      biologicalFeatures: resetFeatures,
      featureColorOverrides: {},
      svgContainer: ref({ querySelector: () => null }),
      clickedFeature: resetClicked,
      featureStyleScopeDialog: {},
      resetColorDialog: resetDialog,
      legendRenameDialog: {},
      legendEntries: ref([]),
      legendStrokeOverrides: {},
      legendColorOverrides: {},
      originalLegendOrder: ref([]),
      originalLegendColors: ref({}),
      featureStrokeOverrides: {},
      skipCaptureBaseConfig: ref(false),
      skipExtractOnSvgChange: ref(false),
      addedLegendCaptions: ref(new Set())
    }),
    nextTick: async () => {},
    ruleActions: {
      runWithRuleMatches: runWithRuleMatchesOf(resetPreparation, { extractedFeatures: resetFeatures, biologicalFeatures: resetFeatures, manualSpecificRules: resetRules }),
      commitSpecificRules: async (next) => {
        committed.push((await next(resetRules, (row, patch) => ({ ...row, ...patch }))).map((rule) => ({ ...rule })));
        return true;
      },
      effectiveLegendCaptions: () => (feature) => feature.type,
      getLegendRowRules: (caption) => legendRowRules(caption, { rules: resetRules, pythonRows: new Map() }),
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      findFeaturesWithSameLegendItem: () => [],
      findFeaturesWithSameDisplayedLabel: () => [],
      findFeaturesWithSameIndividualLabel: () => [],
      getDisplayedFeatureLabel: (feature) => feature.product,
      getIndividualFeatureLabel: (feature) => feature.product,
      groupColorRules: groupColorRulesOver(resetPreparation, (feature) => [feature.product])
    },
    getFeatureElements: () => [],
    getFeatureFillElements: () => []
  });
  for (const choice of ['cancel', 'this']) {
    committed.length = 0;
    resetClicked.value = { svg_id: trnaA.svg_id, feat: trnaA };
    await resetActions.resetClickedFeatureFillColor();
    assert.equal(resetDialog.show, true, 'a tRNA with a sibling opens the Reset dialog');
    assert.equal(Object.hasOwn(resetDialog, 'defaultColor'), false, 'the dialog holds display values only');
    await resetActions.handleResetColorChoice(choice);
    resetClicked.value = { svg_id: rrna.svg_id, feat: rrna };
    await resetActions.resetClickedFeatureFillColor();
    // Python takes the reset feature's rules back and names its instance (PD-OI-069 revision 2).
    const rrnaReset = committed.at(-1).find((rule) => rule.qual === 'hash' && rule.val === `^${rrna.svg_id}@record-1$`);
    assert.equal(rrnaReset?.color, '#71ee7d', `Reset after a ${choice} dialog uses the rRNA default`);
  }
  // OV-548: a type is file text. `constructor` is no palette key, so the Reset
  // button has no default to offer: no dialog, no rule write.
  const ctorA = feature('ctor-a', 'constructor', 'x', 1000, 1100);
  const ctorB = feature('ctor-b', 'constructor', 'x', 1200, 1300);
  resetFeatures.value = [trnaA, trnaB, rrna, ctorA, ctorB];
  resetDialog.show = false;
  committed.length = 0;
  resetClicked.value = { svg_id: ctorA.svg_id, feat: ctorA };
  await resetActions.resetClickedFeatureFillColor();
  assert.equal(resetDialog.show, false, 'a type that names an inherited member opens no Reset dialog');
  assert.deepEqual(committed, [], 'and Reset writes no rule for it');
}

// PV-02 and PV-04 (D-06, PD-OI-061): Legend editor renames.
{
  const { compileDirectEditorMutationPlan } = await import(pathToFileURL(join(tempDir, 'app', 'candidate-render.js')));
  const trna = { id: 't1', svg_id: 't1', type: 'tRNA', product: 'tRNA-Leu', qualifiers: { product: ['tRNA-Leu'] }, start: 1, end: 70 };
  const fakeEntry = (caption) => {
    const attributes = new Map([['data-legend-key', caption]]);
    const text = { textContent: caption };
    return {
      getAttribute: (name) => attributes.get(name) || null,
      setAttribute: (name, value) => attributes.set(name, value),
      querySelector: (selector) => selector === 'text' ? text : null,
      querySelectorAll: () => []
    };
  };
  const build = ({ entries, order, rules = [], features = [], legendColorOverrides = {}, deleted = [], scopeDialog = {}, clicked = null, legendRowIsNamed,
    legendRowName, findLegendRowNamed, keepBothNamesOf, ruleOverrides = {} }) => {
    // Each derivation of the Legend row colors reads the displayed Result once.
    const derivations = { count: 0 };
    const groupEntries = new Map(entries.map((entry) => [entry.caption, fakeEntry(entry.caption)]));
    const legendGroup = {
      querySelector: (selector) => [...groupEntries].find(([caption]) => selector.includes(`"${caption}"`))?.[1] || null,
      querySelectorAll: () => []
    };
    const svgRoot = legendSvg(legendGroup);
    const committed = [];
    const commitOptions = [];
    const errorLog = ref(null);
    const legendRenameDialog = {};
    const stateLegendEntries = ref(entries.map((entry) => ({ featureIds: [], originalCaption: entry.caption, ...entry })));
    const originalOrder = ref([...order]);
    // Catalog features: each has its catalog identity, here its rendered id in one record.
    const featureList = ref(features.map((feature) => ({ record_key: 'record-1', biological_feature_id: feature.svg_id, ...feature })));
    const renameState = { extractedFeatures: featureList, biologicalFeatures: featureList, manualSpecificRules: rules };
    const renamePreparation = createRulePreparation({ state: withDrawings(renameState), evaluate: evaluatePythonRules, synthesize: synthesizePythonColorRules });
    const renameActions = createFeatureColorActions({ readDrawnTracks: () => [], writeLegendRowColor: recordLegendRowColor, legendRowIsNamed, legendRowName, findLegendRowNamed, keepBothNamesOf,
      state: withDrawings({
        results: ref([]), selectedResultIndex: ref(0), appliedPaletteColors: ref({ tRNA: '#e8b441' }),
        manualSpecificRules: rules, extractedFeatures: featureList, biologicalFeatures: featureList,
        featureColorOverrides: {}, svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svgRoot : null }),
        clickedFeature: ref(clicked), featureStyleScopeDialog: scopeDialog, resetColorDialog: {}, legendRenameDialog,
        legendEntries: stateLegendEntries, deletedLegendEntries: ref(deleted), dormantLegendEntries: ref([]),
        legendStrokeOverrides: {}, legendColorOverrides,
        originalLegendOrder: originalOrder, originalLegendColors: ref({}),
        featureStrokeOverrides: {}, skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false),
        addedLegendCaptions: ref(new Set()), errorLog,
        svgContent: { get value() { derivations.count += 1; return null; } }
      }),
      nextTick: async () => {},
      ruleActions: {
        runWithRuleMatches: runWithRuleMatchesOf(renamePreparation, renameState),
        commitSpecificRules: async (next, _label, options = {}) => {
          committed.push((await next(rules, (row, patch) => ({ ...row, ...patch }))).map((rule) => ({ ...rule })));
          commitOptions.push(options);
          return true;
        },
        effectiveLegendCaptions: () => (feature) => feature.legendCaption || rules.find((rule) => rule.feat === feature.type)?.cap || feature.type,
        // Python's rows of the displayed Result: the listed rows as Python drew them.
        getLegendRowRules: (caption) => legendRowRules(caption, {
          rules, originalLegendOrder: originalOrder.value, pythonRows: new Map(stateLegendEntries.value.map((entry) => [
            entry.originalCaption, { key: entry.originalCaption, color: entry.color }
          ]))
        }),
        getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
        findFeaturesWithSameLegendItem: () => [], findFeaturesWithSameDisplayedLabel: () => [],
        findFeaturesWithSameIndividualLabel: () => [], getDisplayedFeatureLabel: (feature) => feature.product,
        getIndividualFeatureLabel: (feature) => feature.product,
        groupColorRules: groupColorRulesOver(renamePreparation, (feature) => [feature.product]),
        ...ruleOverrides
      },
      getFeatureElements: () => [],
      getFeatureFillElements: () => []
    });
    return { renameActions, committed, commitOptions, errorLog, scopeDialog, legendRenameDialog, stateLegendEntries, originalOrder, legendColorOverrides, deleted, derivations };
  };

  // A row's name field shows the name the Legend draws (`Shared (1)`), not the
  // row's key: typing that, or the row's own name, is no rename (PD-OI-042
  // scenario revision 2); another name renames the row.
  const sharedRows = {
    entries: [{ caption: 'Shared [#112233]', color: '#112233' }, { caption: 'Shared [#445566]', color: '#445566' }],
    order: ['Shared [#112233]', 'Shared [#445566]'],
    legendRowIsNamed: (caption, typed) => typed === caption || (caption === 'Shared [#112233]' && ['Shared', 'Shared (1)'].includes(typed))
  };
  for (const typed of ['Shared (1)', 'Shared', 'Shared [#112233]']) {
    const own = build(sharedRows);
    await own.renameActions.renameLegendEntry(0, typed);
    assert.equal(own.stateLegendEntries.value[0].caption, 'Shared [#112233]', `"${typed}" is the row itself`);
    assert.notEqual(own.legendRenameDialog.show, true, `"${typed}" asks nothing`);
  }
  const otherName = build(sharedRows);
  await otherName.renameActions.renameLegendEntry(0, 'Packaging');
  assert.equal(otherName.stateLegendEntries.value[0].caption, 'Packaging');
  // OV-1216: a renamed row's key names it back; its name is no rename.
  for (const [typed, caption] of [['depth', 'depth'], ['Coverage', 'Coverage']]) {
    const renamed = build({
      entries: [{ caption: 'Coverage', originalCaption: 'depth', color: '#555555' }],
      order: ['depth'],
      legendRowIsNamed: (key, text) => text === key || (key === 'depth' && text === 'Coverage')
    });
    await renamed.renameActions.renameLegendEntry(0, typed);
    assert.equal(renamed.stateLegendEntries.value[0].caption, caption, `"${typed}" on a row renamed Coverage`);
  }
  // PD-OI-061 scenario revision 5: a shown numbered name (`Shared (2)`) names
  // that row, which the Legend owner finds by its drawn name; a Merge adopts
  // the row's name, never its key (OV-640).
  const numbered = build({
    entries: [...sharedRows.entries, { caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }],
    order: [...sharedRows.order, 'tRNA'],
    features: [{ ...trna, legendCaption: 'tRNA' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Shared [#445566]' }],
    findLegendRowNamed: (typed) => ({ 'Shared (1)': 'Shared [#112233]', 'Shared (2)': 'Shared [#445566]' })[typed] ?? null,
    legendRowName: (key) => (/^Shared \[/.test(key) ? 'Shared' : key)
  });
  await numbered.renameActions.renameLegendEntry(2, 'Shared (2)');
  assert.deepEqual([numbered.legendRenameDialog.show, numbered.legendRenameDialog.mode, numbered.legendRenameDialog.mergeAvailable,
    numbered.legendRenameDialog.targetColor], [true, 'target', true, '#445566'], 'a shown numbered name targets its row');
  await numbered.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(numbered.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Shared', '#445566']], 'Merge adopts the name, not the key');
  // Review M-3 (D-37): with that row's color and feature type no dialog asks,
  // and the row joins it under its name, never the typed numbered name.
  const sameRow = build({
    entries: [...sharedRows.entries, { caption: 'tRNA', color: '#445566', featureIds: ['t1'] }],
    order: [...sharedRows.order, 'tRNA'],
    features: [{ ...trna, legendCaption: 'tRNA' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Shared [#445566]' }],
    findLegendRowNamed: (typed) => ({ 'Shared (1)': 'Shared [#112233]', 'Shared (2)': 'Shared [#445566]' })[typed] ?? null,
    legendRowName: (key) => (/^Shared \[/.test(key) ? 'Shared' : key)
  });
  await sameRow.renameActions.renameLegendEntry(2, 'Shared (2)');
  assert.notEqual(sameRow.legendRenameDialog.show, true, 'one color and one type: no dialog');
  assert.deepEqual(sameRow.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Shared', '#445566']],
    'a typed shown numbered name joins that row under its name');

  // PV-02: a renamed renderer-generated row keeps its generated identity, so
  // Generate replays the rename.
  const gc = build({
    entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'GC content', color: '#a1a1a1' }],
    order: ['CDS', 'GC content']
  });
  await gc.renameActions.renameLegendEntry(1, 'GC percent');
  // Review L4: a rename onto a caption no row shows derives the row colors
  // once, for the source row's color, and none for the absent target.
  assert.equal(gc.derivations.count, 1, 'one Legend row color derivation');
  assert.equal(gc.stateLegendEntries.value[1].caption, 'GC percent');
  assert.equal(gc.stateLegendEntries.value[1].originalCaption, 'GC content');
  assert.deepEqual(gc.originalOrder.value, ['CDS', 'GC content']);
  const gcPlan = compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: gc.stateLegendEntries.value,
    originalLegendOrder: gc.originalOrder.value
  });
  assert.deepEqual(gcPlan.operationsByResult[0].legendRenames.map(({ from, to }) => [from, to]), [['GC content', 'GC percent']]);

  // PR 3 (D-25): a Legend row's styles are stored under its key, so a rename
  // moves none: the renamed row keeps its own, and a style stored under the
  // key of another row (OV-60's earlier row) stays that row's. Generate
  // replays the rename also once a color rule without a drawn row carries the
  // name (Add rule).
  const unusedRule = { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: '#c83366', cap: 'Zeta' };
  const gcRows = { entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'GC content', color: '#a1a1a1' }], order: ['CDS', 'GC content'] };
  const staleStyle = build({ ...gcRows, legendColorOverrides: { Zeta: '#7b2cbf', 'GC content': '#123123' } });
  await staleStyle.renameActions.renameLegendEntry(1, 'Zeta');
  assert.notEqual(staleStyle.legendRenameDialog.show, true);
  assert.deepEqual(staleStyle.legendColorOverrides, { Zeta: '#7b2cbf', 'GC content': '#123123' }, 'a rename moves no style');
  const onRulePlan = compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: staleStyle.stateLegendEntries.value,
    originalLegendOrder: staleStyle.originalOrder.value,
    legendColorOverrides: staleStyle.legendColorOverrides,
    manualSpecificRules: [unusedRule]
  });
  assert.deepEqual(onRulePlan.operationsByResult[0].legendRenames.map(({ from, to, allowMissing }) => [from, to, allowMissing]),
    [['GC content', 'Zeta', true]]);
  assert.deepEqual(onRulePlan.operationsByResult[0].legendFills.map(({ caption, color }) => [caption, color]),
    [['Zeta', '#7b2cbf'], ['GC content', '#123123']], 'each style fills the row of its key');
  // PD-OI-061 revision 4: the name of a color rule is a rename target also
  // when the rule draws no row. A row without features cannot merge into it:
  // Keep both or Cancel, and no rule commit.
  const onRule = build({ ...gcRows, rules: [unusedRule] });
  await onRule.renameActions.renameLegendEntry(1, 'Zeta');
  assert.equal(onRule.legendRenameDialog.show, true, 'a rename onto the name of a rule without a row asks');
  assert.deepEqual([onRule.legendRenameDialog.mergeAvailable, onRule.legendRenameDialog.targetCaption, onRule.legendRenameDialog.targetColor,
    onRule.legendRenameDialog.keepBothName], [false, 'Zeta', '#c83366', 'Zeta']);
  await onRule.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(onRule.committed, [], 'no rule commit');
  assert.deepEqual(onRule.stateLegendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['CDS', 'CDS'], ['Zeta', 'GC content']], 'Keep both keeps the typed name (PD-OI-061 scenario revision 5)');
  // OV-294 residual: a renamed row's rules are read by its Python key, not by
  // the caption it shows, so renaming a row that shows a rule's caption
  // renames the row and leaves that rule alone.
  const showsRuleCaption = build({
    entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'Zeta', originalCaption: 'GC content', color: '#a1a1a1' }],
    order: ['CDS', 'GC content'], rules: [unusedRule]
  });
  await showsRuleCaption.renameActions.renameLegendEntry(1, 'GC percent');
  assert.deepEqual(showsRuleCaption.committed, [], 'no rule commit');
  assert.deepEqual(showsRuleCaption.stateLegendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['CDS', 'CDS'], ['GC percent', 'GC content']]);
  // Review M1 of 2b96350f: the popup and Apply to all hold the rule's row key
  // (the feature's effective caption), so a renamed row that shows the same
  // caption does not take the rule's place: Apply to all recolors the rule and
  // a group rename recaptions it, as on 579ea705.
  const locusFeature = (tag) => ({ id: tag, svg_id: tag, type: 'CDS', product: 'p', qualifiers: { locus_tag: [tag] }, start: 1, end: 90 });
  const [fl1, fl2] = [locusFeature('FL1'), locusFeature('FL2')];
  const drawnRule = { feat: 'CDS', qual: 'locus_tag', val: '^FL[12]$', color: '#c83366', cap: 'Zeta' };
  const shownByRenamedRow = (options) => build({
    entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'Zeta', originalCaption: 'GC content', color: '#a1a1a1' }],
    order: ['CDS', 'GC content'], rules: [drawnRule], features: [fl1, fl2], ...options
  });
  const applyToAll = shownByRenamedRow({ scopeDialog: { feat: fl1, color: '#123456', legendName: 'Zeta' } });
  await applyToAll.renameActions.handleColorScopeChoice('caption');
  assert.deepEqual(applyToAll.committed.map((rules) => rules.map(({ val, color, cap }) => [val, color, cap])),
    [[['^FL[12]$', '#123456', 'Zeta']]], 'Apply to all recolors the rule');
  const popupRename = shownByRenamedRow({ clicked: { svg_id: 'FL1', feat: fl1, legendName: 'Eta', appliedLegendName: 'Zeta' } });
  await popupRename.renameActions.handleLegendNameCommit();
  await popupRename.renameActions.handleLegendRenameChoice('group');
  assert.deepEqual(popupRename.committed.map((rules) => rules.map(({ val, cap }) => [val, cap])),
    [[['^FL[12]$', 'Eta']]], 'a group rename recaptions the rule');

  // PV-04 (PD-OI-061 amended, OV-62): rows that draw features of one same type
  // ask Merge, Suffix, or Cancel before any rule commit. Any other pair, with
  // or without features, asks Suffix or Cancel only.
  const collide = () => build({
    entries: [{ caption: 'Leu A', color: '#e8b441', featureIds: ['t1'] }, { caption: 'Leu B', color: '#71ee7d', featureIds: ['t2'] }],
    order: ['Leu A', 'Leu B'],
    features: [{ ...trna, legendCaption: 'Leu A' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Leu B' }]
  });
  const merge = collide();
  await merge.renameActions.renameLegendEntry(0, 'Leu B');
  assert.equal(merge.legendRenameDialog.show, true);
  assert.equal(merge.legendRenameDialog.mode, 'target');
  assert.equal(merge.legendRenameDialog.mergeAvailable, true);
  assert.equal(merge.committed.length, 0);
  await merge.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(merge.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Leu B', '#71ee7d']]);
  const keepBoth = collide();
  await keepBoth.renameActions.renameLegendEntry(0, 'Leu B');
  await keepBoth.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(keepBoth.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Leu B', '#e8b441']]);
  const cancel = collide();
  await cancel.renameActions.renameLegendEntry(0, 'Leu B');
  await cancel.renameActions.handleLegendRenameChoice('cancel');
  assert.equal(cancel.committed.length, 0);
  assert.equal(cancel.legendRenameDialog.show, false);

  // OV-648: the Color Change Scope dialog leaves out a label scope that selects
  // the features the Legend item scope selects (or, for the annotation label,
  // the display label scope), whatever the names: a row key such as
  // `Shared [#112233]` is no label, and a label that repeats a row's name can
  // select other features.
  const scoped = (legendName, { legend, display, annotation }) => {
    const others = (ids) => () => ids.map((id) => ({ ...trna, id, svg_id: id, product: 'Shared' }));
    const fixture = build({ entries: [{ caption: legendName, color: '#e8b441', featureIds: ['t1'] }], order: [legendName],
      features: [{ ...trna, product: 'Shared', legendCaption: legendName }],
      scopeDialog: {}, ruleOverrides: { findFeaturesWithSameLegendItem: others(legend), findFeaturesWithSameDisplayedLabel: others(display),
        findFeaturesWithSameIndividualLabel: others(annotation), getDisplayedFeatureLabel: () => 'Shared',
        findMatchingRegexRule: () => null, countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null } });
    return fixture;
  };
  for (const [name, legendName, sets, offered] of [
    ['a keyed row selecting the labeled features', 'Shared [#112233]', { legend: ['t2'], display: ['t2'], annotation: ['t2'] }, [false, false]],
    ['a row named like the label selecting other features', 'Shared', { legend: ['t2'], display: ['t2', 't3'], annotation: ['t2', 't3'] }, [true, false]],
    ['an annotation label selecting other features than both', 'Shared', { legend: ['t2'], display: ['t2'], annotation: ['t2', 't4'] }, [false, true]]
  ]) {
    const fixture = scoped(legendName, sets);
    await fixture.renameActions.requestFeatureColorChange({ ...trna, product: 'Shared', legendCaption: legendName }, '#123456');
    assert.equal(fixture.errorLog.value, null, name);
    assert.deepEqual([fixture.scopeDialog.displayLabelOffered, fixture.scopeDialog.annotationLabelOffered], offered, name);
  }

  // Review r1 (OV-244): a group edit whose rule synthesis Python refuses shows
  // the failure, on the Legend panel rename and on Selected features alike.
  for (const [name, run] of [
    ['the Legend panel rename', (fixture) => fixture.renameActions.renameLegendEntry(0, 'Leucine')],
    ['Selected features', (fixture) => fixture.renameActions.applyColorToSelectedFeatures([{ ...trna, legendCaption: 'tRNA' }], '#123456', 'Leucine')]
  ]) {
    const synthesized = [];
    const fixture = build({ entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }], order: ['tRNA'],
      features: [{ ...trna, legendCaption: 'tRNA' }],
      ruleOverrides: { groupColorRules: () => async () => {
        synthesized.push(name);
        throw new Error('The rule synthesis could not color the selected features.');
      } } });
    await run(fixture);
    assert.deepEqual(synthesized, [name], `${name} asks Python for its rules`);
    assert.notEqual(fixture.errorLog.value, null, `${name} reports the refused synthesis`);
    assert.deepEqual(fixture.committed, [], `${name} writes no rules`);
  }

  // OV-2189 (Owner 2026-10-29): a popup Legend name on the only feature of a
  // rule row asks Python for its rules, which keep the renamed rule in its
  // place in the table (the old rule colored the feature there), so the row
  // keeps its place at Generate (OV-158).
  {
    const alphaRule = { feat: 'tRNA', qual: 'product', val: '^alpha$', color: '#e63946', cap: 'alpha' };
    const betaRule = { feat: 'tRNA', qual: 'product', val: '^beta$', color: '#2a9d8f', cap: 'beta' };
    const alpha = { ...trna, product: 'alpha', qualifiers: { product: ['alpha'] }, legendCaption: 'alpha' };
    const beta = { ...trna, id: 't2', svg_id: 't2', product: 'beta', qualifiers: { product: ['beta'] }, start: 101, end: 170, legendCaption: 'beta' };
    const fixture = build({ entries: [{ caption: 'alpha', color: '#e63946', featureIds: ['t1'] }, { caption: 'beta', color: '#2a9d8f', featureIds: ['t2'] }],
      order: ['alpha', 'beta'], rules: [alphaRule, betaRule], features: [alpha, beta],
      clicked: { svg_id: 't2', feat: beta, color: '#2a9d8f', legendName: 'gamma', appliedLegendName: 'beta' } });
    await fixture.renameActions.handleLegendNameCommit();
    assert.equal(fixture.errorLog.value, null);
    assert.deepEqual(fixture.committed.at(-1).map(({ val, cap, color }) => [val, cap, color]),
      [['^alpha$', 'alpha', '#e63946'], ['^beta$', 'gamma', '#2a9d8f']], 'the renamed rule keeps its place');
    assert.deepEqual(fixture.commitOptions.at(-1).legendPlacement, { name: 'gamma', color: '#2a9d8f', at: 'beta', replaces: true });
  }

  // PD-OI-061 scenario revision 5: the dialog names what Keep both draws (the
  // Legend owner's `keepBothNamesOf`). A row renamed whole keeps its place and
  // its key's records; a feature split off its row is placed right after the
  // target row (D-29).
  const names = [];
  const keepBothNamesOf = (rename) => {
    names.push(rename);
    return { keepBothName: 'Leu B (1)', targetNameAfter: 'Leu B (2)', restores: false };
  };
  const whole = build({ ...{ entries: [{ caption: 'Leu A', color: '#e8b441', featureIds: ['t1'] }, { caption: 'Leu B', color: '#71ee7d', featureIds: ['t2'] }],
    order: ['Leu A', 'Leu B'], features: [{ ...trna, legendCaption: 'Leu A' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Leu B' }] }, keepBothNamesOf });
  await whole.renameActions.renameLegendEntry(0, 'Leu B');
  assert.deepEqual(names.at(-1), { name: 'Leu B', renamedKey: 'Leu A', targetKey: 'Leu B', color: '#e8b441', deletedKey: null });
  assert.deepEqual([whole.legendRenameDialog.keepBothName, whole.legendRenameDialog.targetNameAfter], ['Leu B (1)', 'Leu B (2)']);
  await whole.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(whole.committed.at(-1).map(({ cap }) => cap), ['Leu B'], 'the typed name, not the drawn one');
  assert.deepEqual(whole.commitOptions.at(-1).legendPlacement, { name: 'Leu B', color: '#e8b441', at: 'Leu A', replaces: true });
  const leuA2 = { ...trna, id: 't3', svg_id: 't3', legendCaption: 'Leu A' };
  const split = build({ entries: [{ caption: 'Leu A', color: '#e8b441', featureIds: ['t1', 't3'] }, { caption: 'Leu B', color: '#71ee7d', featureIds: ['t2'] }],
    order: ['Leu A', 'Leu B'], features: [{ ...trna, legendCaption: 'Leu A' }, leuA2, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Leu B' }],
    clicked: { svg_id: 't1', feat: { ...trna, legendCaption: 'Leu A' }, legendName: 'Leu B', appliedLegendName: 'Leu A' }, keepBothNamesOf });
  await split.renameActions.handleLegendNameCommit();
  await split.renameActions.handleLegendRenameChoice('single');
  assert.equal(split.legendRenameDialog.mode, 'target');
  assert.deepEqual(names.at(-1), { name: 'Leu B', renamedKey: null, targetKey: 'Leu B', color: '#e8b441', deletedKey: null }, 'a split row is named as placed after the target');
  await split.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(split.commitOptions.at(-1).legendPlacement, { name: 'Leu B', color: '#e8b441', at: 'Leu B' }, 'D-29: right after the target row');

  // Different types, a row without features, and a target without live features
  // offer no Merge; a programmatic Merge is refused like Cancel (OV-62).
  // `source` is the renamed row's index, `target` the caption it takes.
  const gcSkew = (color) => build({
    entries: [{ caption: 'CDS', color: '#54bcf8', featureIds: ['t1'] },
      { caption: 'GC skew (+)', color: '#6dded3' }, { caption: 'GC skew (-)', color }],
    order: ['CDS', 'GC skew (+)', 'GC skew (-)'],
    features: [{ ...trna, type: 'CDS', legendCaption: 'CDS' }]
  });
  const noMerge = {
    differentTypes: { source: 0, target: 'rRNA', make: () => build({
      entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'rRNA', color: '#71ee7d', featureIds: ['r1'] }],
      order: ['tRNA', 'rRNA'],
      features: [trna, { ...trna, id: 'r1', svg_id: 'r1', type: 'rRNA' }]
    }) },
    // The colors match: the row still asks, as a silent rename would draw two rows of one caption.
    differentTypesSameColor: { source: 0, target: 'rRNA', make: () => build({
      entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'rRNA', color: '#e8b441', featureIds: ['r1'] }],
      order: ['tRNA', 'rRNA'],
      features: [trna, { ...trna, id: 'r1', svg_id: 'r1', type: 'rRNA' }]
    }) },
    featureless: { source: 1, target: 'tRNA', make: () => build({
      entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'GC content', color: '#a1a1a1' }],
      order: ['tRNA', 'GC content'],
      features: [trna]
    }) },
    // OV-151 (GUI audit FL-04, PD-OI-061 revision 2): two rows without features
    // on one track, GC skew (+) onto GC skew (-), with other and with equal colors.
    featurelessSameTrack: { source: 1, target: 'GC skew (-)', renamed: true, make: () => gcSkew('#ad72e3') },
    featurelessSameTrackSameColor: { source: 1, target: 'GC skew (-)', renamed: true, make: () => gcSkew('#6dded3') }
  };
  for (const [name, { source, target, renamed, make }] of Object.entries(noMerge)) {
    const refused = make();
    const captions = refused.stateLegendEntries.value.map((entry) => entry.caption);
    await refused.renameActions.renameLegendEntry(source, target);
    assert.equal(refused.legendRenameDialog.show, true, name);
    assert.equal(refused.legendRenameDialog.mergeAvailable, false, name);
    await refused.renameActions.handleLegendRenameChoice('merge');
    assert.equal(refused.legendRenameDialog.show, false, name);
    assert.equal(refused.committed.length, 0, name);
    assert.deepEqual(refused.stateLegendEntries.value.map((entry) => entry.caption), captions, name);
    const keptBoth = make();
    await keptBoth.renameActions.renameLegendEntry(source, target);
    await keptBoth.renameActions.handleLegendRenameChoice('keep-both');
    const keptCaption = keptBoth.committed.length
      ? keptBoth.committed.at(-1)[0].cap : keptBoth.stateLegendEntries.value[source].caption;
    assert.equal(keptCaption, target, `${name}: Keep both keeps the typed name (PD-OI-061 scenario revision 5)`);
    if (!renamed) continue;
    // A row without features renames in the Legend only: no rule commit, one
    // rename that Generate replays, and nothing that reaches the target row.
    assert.equal(keptBoth.committed.length, 0, name);
    assert.deepEqual(keptBoth.stateLegendEntries.value.map((entry) => entry.caption),
      captions.map((caption, index) => (index === source ? target : caption)), `${name}: the renamed row keeps its place`);
    const plan = compileDirectEditorMutationPlan({
      catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
      legendEntries: keptBoth.stateLegendEntries.value,
      originalLegendOrder: keptBoth.originalOrder.value
    });
    const operations = plan.operationsByResult[0];
    assert.deepEqual(operations.legendRenames.map(({ from, to }) => [from, to]), [[captions[source], target]], name);
    // Operations address rows by key; a rename's `to` is a name.
    const touchingTarget = Object.entries(operations).flatMap(([kind, list]) => (Array.isArray(list) ? list : [])
      .filter((operation) => [operation?.caption, operation?.from].includes(target))
      .map((operation) => [kind, operation]));
    assert.deepEqual(touchingTarget, [], `${name}: no operation reaches ${target}`);
  }

  // PD-OI-061 revision 4 (OV-580 to OV-582): the name of a color rule is a
  // rename target like any other row's, drawn or not, and a pick from the
  // popup's list joins the picked group with its name and its color.
  const specialRule = { feat: 'tRNA', qual: 'product', val: '^NOMATCH$', color: '#ff0000', cap: 'Special' };
  const connectorRule = { feat: 'tRNA', qual: 'product', val: '^NOMATCH2$', color: '#5a5a5a', cap: 'Connector' };
  const cdsRule = { feat: 'CDS', qual: 'product', val: '^NOMATCH3$', color: '#00aa00', cap: 'Tail' };
  // The groups are the rule table's (specific-color-rules.js): one per name
  // and color, with the feature types of its rules. The popup lists those of
  // the clicked feature's type.
  assert.deepEqual(colorRuleLegendGroups([specialRule, connectorRule, cdsRule, { ...specialRule, val: '^NOMATCH4$' },
    { ...cdsRule, feat: 'tRNA' }, { ...cdsRule, cap: '' }]), [
    { name: 'Special', color: '#ff0000', types: ['tRNA'] },
    { name: 'Connector', color: '#5a5a5a', types: ['tRNA'] },
    { name: 'Tail', color: '#00aa00', types: ['CDS', 'tRNA'] }
  ]);
  assert.deepEqual(colorRuleLegendGroups([specialRule, connectorRule, cdsRule], 'tRNA').map(({ name }) => name), ['Special', 'Connector'],
    'a rule of another feature type is not offered');
  assert.deepEqual(colorRuleLegendGroups([specialRule, { ...specialRule, color: '#F00' }, { ...specialRule, color: '#0000ff' }], 'tRNA'),
    [{ name: 'Special', color: '#ff0000', types: ['tRNA'] }, { name: 'Special', color: '#0000ff', types: ['tRNA'] }],
    'one group per name and color');
  const [special, connector] = colorRuleLegendGroups([specialRule, connectorRule], 'tRNA');
  const [tail] = colorRuleLegendGroups([cdsRule], 'CDS');
  const t1 = { ...trna, legendCaption: 'tRNA' };
  const t2 = { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Special' };
  const t3 = { ...trna, id: 't3', svg_id: 't3', legendCaption: 'tRNA' };
  const ruleRows = (options = {}) => build({
    entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'Special', color: '#ff0000', featureIds: ['t2'] }],
    order: ['tRNA', 'Special'],
    rules: [specialRule, connectorRule, cdsRule],
    features: [t1, t2],
    clicked: { svg_id: 't1', feat: t1, legendName: 'tRNA', appliedLegendName: 'tRNA' },
    ...options
  });
  // A feature's own rule names its instance (`^<id>@record-1$`); read it by the feature's id.
  const ownId = (val) => /^\^(.*)@record-1\$$/.exec(val)?.[1] ?? val;
  const committedFor = (fixture, ...ids) => (fixture.committed.at(-1) || []).filter(({ val }) => ids.includes(ownId(val)))
    .map(({ val, cap, color }) => [ownId(val), cap, color]);

  // A typed rename onto a rule's listed row asks; Merge adopts its color.
  const ruleOwned = ruleRows();
  await ruleOwned.renameActions.renameLegendEntry(0, 'Special');
  assert.equal(ruleOwned.legendRenameDialog.show, true, 'a rename onto a color rule\'s row asks');
  assert.deepEqual([ruleOwned.legendRenameDialog.mode, ruleOwned.legendRenameDialog.mergeAvailable, ruleOwned.legendRenameDialog.reason],
    ['target', true, '"Special" already exists with another color.']);
  assert.equal(ruleOwned.committed.length, 0);
  await ruleOwned.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(committedFor(ruleOwned, 't1'), [['t1', 'Special', '#ff0000']]);
  // A typed rename onto the name of a rule that draws no row asks too, from
  // the Legend editor and from the popup's field.
  for (const typed of [
    (fixture) => fixture.renameActions.renameLegendEntry(0, 'Connector'),
    (fixture) => fixture.renameActions.handleLegendNameCommit()
  ]) {
    const asked = ruleRows({ clicked: { svg_id: 't1', feat: t1, legendName: 'Connector', appliedLegendName: 'tRNA' } });
    await typed(asked);
    assert.equal(asked.legendRenameDialog.show, true, 'a rename onto the name of a rule without a row asks');
    assert.deepEqual([asked.legendRenameDialog.mode, asked.legendRenameDialog.mergeAvailable, asked.legendRenameDialog.targetCaption,
      asked.legendRenameDialog.targetColor, asked.legendRenameDialog.keepBothName], ['target', true, 'Connector', '#5a5a5a', 'Connector']);
    assert.equal(asked.legendRenameDialog.reason, '"Connector" is the name of a color rule with another color.',
      'OV-586: the dialog says the name is a color rule\'s, not a Legend item\'s');
    assert.equal(asked.committed.length, 0);
    await asked.renameActions.handleLegendRenameChoice('merge');
    assert.deepEqual(committedFor(asked, 't1'), [['t1', 'Connector', '#5a5a5a']]);
    const kept = ruleRows({ clicked: { svg_id: 't1', feat: t1, legendName: 'Connector', appliedLegendName: 'tRNA' } });
    await typed(kept);
    await kept.renameActions.handleLegendRenameChoice('keep-both');
    assert.deepEqual(committedFor(kept, 't1'), [['t1', 'Connector', '#e8b441']]);
  }
  // The name of a rule of another feature type offers no Merge.
  const otherType = ruleRows();
  await otherType.renameActions.renameLegendEntry(0, 'Tail');
  assert.deepEqual([otherType.legendRenameDialog.show, otherType.legendRenameDialog.mergeAvailable], [true, false]);
  assert.equal(otherType.legendRenameDialog.reason, '"Tail" is the name of a color rule with another color.');
  await otherType.renameActions.handleLegendRenameChoice('merge');
  assert.equal(otherType.committed.length, 0, 'a Merge that is not offered is refused');
  // OV-649 (Owner, Q-01 a): in the renamed row's own color the reason is the
  // feature type, not the color.
  const sameColorOtherType = ruleRows({ rules: [specialRule, connectorRule, { ...cdsRule, color: '#e8b441' }] });
  await sameColorOtherType.renameActions.renameLegendEntry(0, 'Tail');
  assert.deepEqual([sameColorOtherType.legendRenameDialog.show, sameColorOtherType.legendRenameDialog.mergeAvailable,
    sameColorOtherType.legendRenameDialog.reason], [true, false, '"Tail" is the name of a color rule for another feature type.']);

  // A pick joins the group without the conflict dialog: a group that draws a
  // row, and one that draws none.
  for (const [group, expected] of [[special, ['t1', 'Special', '#ff0000']], [connector, ['t1', 'Connector', '#5a5a5a']]]) {
    const picked = ruleRows();
    await picked.renameActions.joinLegendGroup(group);
    assert.notEqual(picked.legendRenameDialog.show, true, `a pick of ${group.name} asks nothing`);
    assert.deepEqual(committedFor(picked, 't1'), [expected]);
  }
  // The pick carries its group: the color is the picked group's, not the first
  // rule's of that name.
  const twins = [{ ...connectorRule, cap: 'Twin', color: '#111111' }, { ...connectorRule, cap: 'Twin', color: '#222222' }];
  const twinPick = ruleRows({ rules: twins });
  await twinPick.renameActions.joinLegendGroup(colorRuleLegendGroups(twins, 'tRNA')[1]);
  assert.deepEqual(committedFor(twinPick, 't1'), [['t1', 'Twin', '#222222']]);
  // A typed name that two groups carry is the group in the renamed row's
  // color (one name and one color share a Legend row), not the first group of
  // the name: nothing to ask. In a third color it asks, with the first group.
  const twinTyped = ruleRows({
    rules: twins, entries: [{ caption: 'tRNA', color: '#222222', featureIds: ['t1'] }], order: ['tRNA'], features: [t1]
  });
  await twinTyped.renameActions.renameLegendEntry(0, 'Twin');
  assert.notEqual(twinTyped.legendRenameDialog.show, true, 'the typed name and the row\'s color are a group\'s: nothing to ask');
  // OV-244: the row's features are a group; Python writes the exact product rule that colors only them.
  assert.deepEqual(committedFor(twinTyped, '^tRNA-Leu$'), [['^tRNA-Leu$', 'Twin', '#222222']]);
  const twinAsked = ruleRows({ rules: twins });
  await twinAsked.renameActions.renameLegendEntry(0, 'Twin');
  assert.deepEqual([twinAsked.legendRenameDialog.show, twinAsked.legendRenameDialog.targetColor, twinAsked.legendRenameDialog.reason],
    [true, '#111111', '"Twin" is the name of a color rule with another color.'], 'in another color the first group of the name is the target');
  // Other features share the feature's row: the scope question stays, and
  // neither answer is followed by the conflict dialog.
  for (const [scope, ids] of [['single', ['t1']], ['group', ['t1', 't3']]]) {
    const shared = ruleRows({ features: [t1, t2, t3] });
    await shared.renameActions.joinLegendGroup(special);
    assert.deepEqual([shared.legendRenameDialog.show, shared.legendRenameDialog.mode, shared.legendRenameDialog.siblingCount],
      [true, 'scope', 1]);
    assert.equal(shared.committed.length, 0);
    await shared.renameActions.handleLegendRenameChoice(scope);
    assert.equal(shared.legendRenameDialog.show, false, `${scope}: no conflict dialog`);
    assert.deepEqual(committedFor(shared, 't1', 't3'), ids.map((id) => [id, 'Special', '#ff0000']));
  }
  // Where Merge is not allowed the pick asks as a typed rename does: a group
  // of another feature type, and a group whose row is deleted (OV-285).
  const pickedOtherType = ruleRows();
  await pickedOtherType.renameActions.joinLegendGroup(tail);
  assert.deepEqual([pickedOtherType.legendRenameDialog.show, pickedOtherType.legendRenameDialog.mode,
    pickedOtherType.legendRenameDialog.mergeAvailable], [true, 'target', false]);
  assert.equal(pickedOtherType.committed.length, 0);
  // D-07: a rule for every feature type (`*`) is a rule for each type, so the
  // popup lists its group for any feature. A pick joins every group the list
  // offers: a group of `*` rules only, a group with rules of several types,
  // and a group whose row draws only features of another type.
  const everyRule = { feat: '*', qual: 'note', val: '^NOMATCH5$', color: '#123456', cap: 'Everything' };
  const mixedRules = [{ ...cdsRule, cap: 'Mixed' }, { ...cdsRule, cap: 'Mixed', feat: 'tRNA', val: '^NOMATCH6$' }];
  assert.deepEqual(colorRuleLegendGroups([everyRule, cdsRule, ...mixedRules], 'tRNA').map(({ name, types }) => [name, types]),
    [['Everything', ['*']], ['Mixed', ['CDS', 'tRNA']]], 'the list offers the group of a rule for every type');
  const c1 = { id: 'c1', svg_id: 'c1', type: 'CDS', product: 'p', qualifiers: { product: ['p'] }, start: 100, end: 190, legendCaption: 'Mixed' };
  // OV-244: where t1 is the only tRNA, Python writes the exact product rule that colors only it.
  for (const [rules, name, color, options, ruleValue = 't1'] of [
    [[everyRule], 'Everything', '#123456', {}],
    [mixedRules, 'Mixed', '#00aa00', {}],
    [mixedRules, 'Mixed', '#00aa00', {
      entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'Mixed', color: '#00aa00', featureIds: ['c1'] }],
      order: ['tRNA', 'Mixed'], features: [t1, c1]
    }, '^tRNA-Leu$']
  ]) {
    const picked = ruleRows({ rules, ...options });
    await picked.renameActions.joinLegendGroup(colorRuleLegendGroups(rules, 'tRNA')[0]);
    assert.notEqual(picked.legendRenameDialog.show, true, `a pick of ${name} asks nothing`);
    assert.deepEqual(committedFor(picked, ruleValue), [[ruleValue, name, color]]);
    // A typed rename keeps the Merge rule of rows of one same feature type.
    const typed = ruleRows({ rules, ...options });
    await typed.renameActions.renameLegendEntry(0, name);
    assert.deepEqual([typed.legendRenameDialog.show, typed.legendRenameDialog.mergeAvailable], [true, false],
      `a typed rename onto ${name} offers no Merge`);
  }
  // D-11 (OV-587): the two results above are told apart by where the name
  // comes from, so the Legend name field suggests no group: a group's name is
  // a typed rename only when the user types it, and the list under the field
  // is the only place that offers a group.
  const template = await readFile(join(repoRoot, 'gbdraw', 'web', 'index.html'), 'utf8');
  const nameField = template.match(/<input\b[^>]*aria-label="t\('Feature legend name'\)"[^>]*>/)?.[0] || '';
  assert.match(nameField, /@change="handleLegendNameCommit"/, 'the popup has its Legend name field');
  assert.doesNotMatch(nameField, /\blist=/, 'the Legend name field suggests no group');
  assert.doesNotMatch(template, /<datalist\b[^>]*legend-name/);
  const pickedDeleted = ruleRows({
    entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }],
    deleted: [{ caption: 'Special', originalCaption: 'Special', color: '#ff0000', featureIds: ['t2'] }]
  });
  await pickedDeleted.renameActions.joinLegendGroup(special);
  assert.deepEqual([pickedDeleted.legendRenameDialog.show, pickedDeleted.legendRenameDialog.mode,
    pickedDeleted.legendRenameDialog.deletedTargetKey], [true, 'target', 'Special']);
  assert.equal(pickedDeleted.committed.length, 0);
  pickedDeleted.stateLegendEntries.value = [...pickedDeleted.stateLegendEntries.value, ...pickedDeleted.deleted.splice(0)];
  await pickedDeleted.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(committedFor(pickedDeleted, 't1'), [['t1', 'Special', '#ff0000']]);

  // R15-3 (OV-285): a rename onto a deleted row's caption asks as onto a
  // listed row. The dialog names the deleted row by its Python key; the merge
  // goes into that row once the composition root restored it, so this owner
  // alone refuses it like a Merge it does not offer.
  // The Legend owner answers whether the deleted row has the renamed row's
  // swatch (review L-1); here, whether it has its color.
  /** @type {unknown[]} */
  const asks = [];
  const ontoDeleted = (deletedColor) => build({
    entries: [{ caption: 'Leu A', color: '#e8b441', featureIds: ['t1'] }],
    deleted: [{ caption: 'Leu B', originalCaption: 'tRNA', color: deletedColor, featureIds: ['t2'] }],
    order: ['tRNA'],
    // A feature's Legend item is its row's key: the deleted row's is `tRNA`.
    features: [{ ...trna, legendCaption: 'Leu A' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'tRNA' }],
    keepBothNamesOf: (rename) => {
      asks.push(rename.deletedKey);
      return { keepBothName: rename.name, targetNameAfter: null, restores: deletedColor === '#e8b441' };
    }
  });
  for (const deletedColor of ['#71ee7d', '#e8b441']) {
    const asked = ontoDeleted(deletedColor);
    await asked.renameActions.renameLegendEntry(0, 'Leu B');
    assert.equal(asked.legendRenameDialog.show, true, `a deleted caption asks (${deletedColor})`);
    assert.equal(asked.legendRenameDialog.mode, 'target');
    assert.equal(asked.legendRenameDialog.mergeAvailable, true, 'one feature type: Restore and merge is offered');
    assert.equal(asked.legendRenameDialog.deletedTargetKey, 'tRNA', 'the deleted row is named by its Python key');
    assert.equal(asked.committed.length, 0);
    await asked.renameActions.handleLegendRenameChoice('merge');
    assert.equal(asked.legendRenameDialog.show, false, 'a merge into a row still deleted is refused');
    assert.equal(asked.committed.length, 0);
    assert.deepEqual(asked.stateLegendEntries.value.map((entry) => entry.caption), ['Leu A']);
    // D-28: onto a deleted row in the renamed row's color, Keep both restores
    // it first (the composition root), as the two are one row.
    assert.equal(asked.legendRenameDialog.keepBothRestores, false, 'the dialog is closed');
  }
  for (const [deletedColor, restores] of [['#71ee7d', false], ['#e8b441', true]]) {
    const asked = ontoDeleted(deletedColor);
    await asked.renameActions.renameLegendEntry(0, 'Leu B');
    assert.equal(asked.legendRenameDialog.keepBothRestores, restores, `Keep both restores the deleted row of one swatch only (${deletedColor})`);
    assert.equal(asks.at(-1), 'tRNA', 'the Legend owner is asked about the deleted row by its key');
    assert.equal(asked.legendRenameDialog.targetNameAfter, '', 'a deleted target gets no sentence');
  }
  const restored = ontoDeleted('#71ee7d');
  await restored.renameActions.renameLegendEntry(0, 'Leu B');
  // The composition root's Restore returns the row to the list before the merge.
  restored.stateLegendEntries.value = [...restored.stateLegendEntries.value, ...restored.deleted.splice(0)];
  await restored.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(restored.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Leu B', '#71ee7d']]);
  // A deleted row a color rule draws: its merge adopts its color, as any row's does.
  const ruleOwnedDeleted = build({
    entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }],
    deleted: [{ caption: 'Special', originalCaption: 'Special', color: '#ff0000', featureIds: ['t2'] }],
    order: ['tRNA', 'Special'],
    rules: [{ feat: 'tRNA', qual: 'product', val: '^NOMATCH$', color: '#ff0000', cap: 'Special' }],
    features: [{ ...trna, legendCaption: 'tRNA' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Special' }]
  });
  await ruleOwnedDeleted.renameActions.renameLegendEntry(0, 'Special');
  assert.equal(ruleOwnedDeleted.legendRenameDialog.show, true, 'a rule-owned deleted caption asks');
  assert.equal(ruleOwnedDeleted.legendRenameDialog.deletedTargetKey, 'Special');
  ruleOwnedDeleted.stateLegendEntries.value = [...ruleOwnedDeleted.stateLegendEntries.value, ...ruleOwnedDeleted.deleted.splice(0)];
  await ruleOwnedDeleted.renameActions.handleLegendRenameChoice('merge');
  assert.deepEqual(ruleOwnedDeleted.committed.at(-1).filter(({ val }) => ownId(val) === 't1').map(({ cap, color }) => [cap, color]),
    [['Special', '#ff0000']]);
  const keptOntoDeleted = ontoDeleted('#71ee7d');
  await keptOntoDeleted.renameActions.renameLegendEntry(0, 'Leu B');
  await keptOntoDeleted.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(keptOntoDeleted.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Leu B', '#e8b441']]);
  const cancelledOntoDeleted = ontoDeleted('#71ee7d');
  await cancelledOntoDeleted.renameActions.renameLegendEntry(0, 'Leu B');
  await cancelledOntoDeleted.renameActions.handleLegendRenameChoice('cancel');
  assert.equal(cancelledOntoDeleted.committed.length, 0);
  assert.equal(cancelledOntoDeleted.legendRenameDialog.show, false);
  assert.equal(cancelledOntoDeleted.legendRenameDialog.deletedTargetKey, '');
  assert.deepEqual(cancelledOntoDeleted.deleted.map((entry) => entry.caption), ['Leu B']);
  // A row without features (GC skew) onto a deleted GC skew row: Keep both
  // or Cancel only; Keep both keeps the typed name beside the deleted row.
  const gcOntoDeleted = () => build({
    entries: [{ caption: 'CDS', color: '#54bcf8', featureIds: ['t1'] }, { caption: 'GC skew (+)', color: '#6dded3' }],
    deleted: [{ caption: 'GC skew (-)', originalCaption: 'GC skew (-)', color: '#ad72e3' }],
    order: ['CDS', 'GC skew (+)', 'GC skew (-)'],
    features: [{ ...trna, type: 'CDS', legendCaption: 'CDS' }]
  });
  const gcAsked = gcOntoDeleted();
  await gcAsked.renameActions.renameLegendEntry(1, 'GC skew (-)');
  assert.equal(gcAsked.legendRenameDialog.show, true, 'a featureless rename onto a deleted caption asks');
  assert.equal(gcAsked.legendRenameDialog.mergeAvailable, false);
  assert.equal(gcAsked.legendRenameDialog.deletedTargetKey, 'GC skew (-)');
  assert.deepEqual(gcAsked.stateLegendEntries.value.map((entry) => entry.caption), ['CDS', 'GC skew (+)']);
  await gcAsked.renameActions.handleLegendRenameChoice('keep-both');
  assert.deepEqual(gcAsked.stateLegendEntries.value.map((entry) => entry.caption), ['CDS', 'GC skew (-)']);
  assert.equal(gcAsked.committed.length, 0);
}

// Review M1 of OV-288: the popup finds the stroked Legend row a feature shows
// through the rule the compile reads (`draftLegendRows`), so a deleted row's
// stroke does not reach it. CDS is stroked and deleted; FL1 is given a stroke
// of its own and then Reset: it shows the stroke Python drew, and a later
// width edit saves that color, not the deleted row's.
{
  const fl1 = { id: 'fl1', svg_id: 'fl1-svg', type: 'CDS', product: 'FL1', start: 1, end: 10 };
  const fl1Attributes = new Map([['fill', '#cccccc'], ['stroke', '#000000'], ['stroke-width', '1']]);
  const fl1Element = {
    getAttribute: (name) => fl1Attributes.get(name) ?? null, setAttribute: () => {}, removeAttribute: () => {}
  };
  const fl1Svg = legendSvg();
  const fl1StrokeOverrides = {};
  const fl1Clicked = ref(null);
  const fl1Features = ref([fl1]);
  const fl1Rules = [];
  const fl1State = { extractedFeatures: fl1Features, biologicalFeatures: fl1Features, manualSpecificRules: fl1Rules };
  const fl1Actions = createFeatureColorActions({ readDrawnTracks: () => [], writeLegendRowColor: recordLegendRowColor,
    state: withDrawings({
      results: ref([]), selectedResultIndex: ref(0), appliedPaletteColors: ref({ CDS: '#cccccc' }),
      manualSpecificRules: fl1Rules, extractedFeatures: fl1Features, biologicalFeatures: fl1Features,
      featureColorOverrides: {}, svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? fl1Svg : null }),
      clickedFeature: fl1Clicked, featureStyleScopeDialog: {}, resetColorDialog: {}, legendRenameDialog: {},
      legendEntries: ref([{ caption: 'GC content', originalCaption: 'GC content', color: '#a1a1a1', featureIds: [] }]),
      deletedLegendEntries: ref([{ caption: 'CDS', originalCaption: 'CDS', color: '#cccccc', featureIds: [] }]),
      dormantLegendEntries: ref([]),
      legendStrokeOverrides: { CDS: { strokeColor: '#e63946', strokeWidth: 3 } }, legendColorOverrides: {},
      originalLegendOrder: ref(['CDS', 'GC content']), originalLegendColors: ref({}),
      featureStrokeOverrides: fl1StrokeOverrides, skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false),
      addedLegendCaptions: ref(new Set())
    }),
    nextTick: async () => {},
    showEditorIntent: () => {},
    ruleActions: {
      runWithRuleMatches: runWithRuleMatchesOf(createRulePreparation({ state: withDrawings(fl1State), evaluate: evaluatePythonRules }), fl1State),
      commitSpecificRules: async () => true,
      countFeaturesMatchingRule: () => 0,
      effectiveLegendCaptions: () => () => 'CDS',
      getLegendRowRules: () => [],
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      findMatchingRegexRule: () => null,
      findFeaturesWithSameLegendItem: () => [], findFeaturesWithSameDisplayedLabel: () => [],
      findFeaturesWithSameIndividualLabel: () => [], getDisplayedFeatureLabel: () => '',
      getIndividualFeatureLabel: () => ''
    },
    getFeatureElements: () => [fl1Element],
    getFeatureFillElements: () => [fl1Element]
  });
  const fl1Key = featureOverrideKey(fl1) || fl1.svg_id;
  fl1Clicked.value = { svg_id: fl1.svg_id, feat: fl1, color: '#cccccc', strokeColor: '#000000', strokeWidth: 1 };
  assert.equal(await fl1Actions.updateClickedFeatureStroke('#2a9d8f', 1), true);
  assert.equal(fl1StrokeOverrides[fl1Key]?.strokeColor, '#2a9d8f');
  assert.equal(await fl1Actions.resetClickedFeatureStroke(), true);
  assert.equal(fl1StrokeOverrides[fl1Key], undefined);
  assert.equal(fl1Clicked.value.strokeColor, '#000000', 'Reset shows the stroke Python drew, not the deleted row\'s');
  await fl1Actions.setClickedFeatureStrokeWidthValue(4);
  assert.deepEqual(
    [fl1StrokeOverrides[fl1Key]?.strokeColor, fl1StrokeOverrides[fl1Key]?.strokeWidth], ['#000000', 4],
    'the width edit saves the stroke Python drew, not the deleted row\'s'
  );
}

// A feature type is file text: a type named like an Object.prototype member
// (`constructor`) is not a palette key, so Reset has no palette color to write.
{
  const ctorFeature = { id: 'feature-ctor', svg_id: 'hash-ctor', type: 'constructor', qualifiers: {}, start: 1, end: 2 };
  manualSpecificRules.splice(0);
  clickedFeature.value = { feat: ctorFeature, svg_id: ctorFeature.svg_id };
  await actions.handleResetColorChoice('this_with_legend');
  assert.deepEqual(manualSpecificRules, [], 'a type that is not a palette key writes no rule');
  clickedFeature.value = null;
}
