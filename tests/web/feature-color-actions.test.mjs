import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';
import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-color-actions-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await cp(sourceDir, tempDir, { recursive: true });
const colorActionsSource = await readFile(join(sourceDir, 'app', 'feature-editor', 'color-actions.js'), 'utf8');

const { createFeatureColorActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'color-actions.js'))
);
const { getFeatureColorRuleHash } = await import(pathToFileURL(join(tempDir, 'services', 'feature-utils.js')));
const { resolveFeatureLabelSelector } = await import(pathToFileURL(join(tempDir, 'services', 'feature-selector.js')));
const { buildLegendIntents, legendRowRules } = await import(pathToFileURL(join(tempDir, 'services', 'specific-color-rules.js')));

assert.doesNotMatch(colorActionsSource, /serializeCleanSvg|results\.value\[[^\]]+\]\s*=/);

const ref = (value) => ({ value });
// The mounted SVG as `getAllFeatureLegendGroups` (services/legend-svg.js) reads it:
// `#legend` holding one `#feature_legend` group, or no legend.
const legendSvg = (featureLegend = null) => ({
  getElementById: (id) => (id === 'legend' && featureLegend
    ? { querySelector: (selector) => (selector === '#feature_legend' ? featureLegend : null) }
    : null)
});

const featureA = {
  id: 'feature-a',
  svg_id: 'hash-a',
  type: 'CDS',
  qualifiers: { gene_kind: 'core biosynthetic genes' },
  start: 1,
  end: 10
};
const featureB = {
  id: 'feature-b',
  svg_id: 'hash-b',
  type: 'CDS',
  qualifiers: { gene_kind: 'core biosynthetic genes' },
  start: 11,
  end: 20
};
const hashOnlyFeature = {
  id: 'feature-c',
  svg_id: 'hash-c',
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
let legendGeometryChangedCount = 0;
const previewFillColors = new Map();
const featureElementsById = new Map();
let previewFillApplyCount = 0;
let previewCommitCount = 0;
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
const commitActiveResultEdit = () => {
  previewCommitCount += 1;
  return true;
};

const { featureOverrideKey } = await import(pathToFileURL(join(tempDir, 'services', 'feature-override-identity.js')));
const { createRulePreparation, runWhenPrepared } = await import(pathToFileURL(join(tempDir, 'app', 'rule-matching.js')));
const { ruleMatcher } = await import(pathToFileURL(join(tempDir, 'services', 'rule-matchers.js')));
const firstMatchingRule = (feature, rules) => ruleMatcher(rules).first(feature);
// The rule owner's `runWithRuleMatches` for these fakes: the color matches of the rules, prepared once.
const runWithRuleMatchesOf = (preparation, state) => (rules, commit) => runWhenPrepared(state, () => [preparation.prepare(rules)], commit);
const preparationState = { extractedFeatures, biologicalFeatures, manualSpecificRules };
let legendSiblings = [featureB, hashOnlyFeature];
let matchingRegexRule = specificRule;
const preparedRuleSets = [];
const runWithRuleMatchesForActions = runWithRuleMatchesOf(createRulePreparation({ state: withDrawings(preparationState), evaluate: evaluatePythonRules }), preparationState);
const actions = createFeatureColorActions({
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
  onLegendGeometryChanged: () => {
    legendGeometryChangedCount += 1;
  },
  extractLegendEntries: () => {},
  ruleActions: {
    runWithRuleMatches: (rules, commit) => {
      preparedRuleSets.push(rules.map((rule) => `${rule.feat}|${rule.qual}|${rule.val}`));
      return runWithRuleMatchesForActions(rules, commit);
    },
    commitSpecificRules: async (rules, _label, {afterCommit = () => {}, previousLegendIntents = []} = {}) => {
      committedLegendIntents.push(previousLegendIntents);
      const preparation = createRulePreparation({state:withDrawings({extractedFeatures,biologicalFeatures,manualSpecificRules}),evaluate:evaluatePythonRules});
      const candidate = await preparation.prepareCandidate(rules);
      if (!candidate) return false;
      const drawn = new Set(extractedFeatures.value.map(feature => firstMatchingRule(feature, candidate.rules)));
      const { intents } = buildLegendIntents(candidate.rules.filter(rule => drawn.has(rule)));
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
    getLegendRowRules: (caption) => legendRowRules(caption, { rules: manualSpecificRules, legendEntries: legendEntries.value }),
    getIndividualFeatureLabel: (feature) => feature.product || '',
    // FE-09 (D-14): "This feature only" always writes the stable hash.
    getFeatureQualifier: (feature) => ({ qual: 'hash', val: getFeatureColorRuleHash(feature) }),
    getLabelSpecificRule: (feature, label) => {
      const selector = resolveFeatureLabelSelector(feature, label);
      return selector
        ? { feat: feature.type, qual: selector.qualifier, val: selector.pattern }
        : null;
    }
  },
  getFeatureElements: (_svg, featureId) => featureElementsById.get(featureId) || [],
  getFeatureFillElements: (_svg, featureId) => featureElementsById.get(featureId) || [],
  commitActiveResultEdit
});

await actions.handleColorScopeChoice('caption');

assert.equal(applySpecificRulesCount, 1);
assert.equal(manualSpecificRules.find(rule => rule.qual === 'gene_kind').color, '#abcdef');
assert.equal(legendEntries.value[0].color, '#abcdef');
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-a'), false);
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-b'), false);
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-c')?.color, '#abcdef');
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-z')?.color, '#444444');
assert.deepEqual(featureColorOverrides['feature-a'], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides['feature-b'], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides['feature-c'], { color: '#abcdef', caption: 'Core' });

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
// rule and its label rules, once it commits; the scope dialog reads only the
// saved rules (OV-225).
{
  const popupFeature = {
    id: 'popup-feature', svg_id: 'hash-popup', type: 'CDS', product: 'popup product',
    qualifiers: { product: ['popup product'] }, start: 40, end: 50
  };
  const ruleKey = (rule) => `${rule.feat}|${rule.qual}|${rule.val}`;
  const savedRules = [specificRule, { feat: 'CDS', qual: 'hash', val: 'hash-z', color: '#444444', cap: 'Other' }];
  const savedKeys = savedRules.map(ruleKey);
  const label = resolveFeatureLabelSelector(popupFeature, 'popup product');
  assert.ok(label, 'the popup feature has a label rule');
  const clickedKeys = [
    ...savedKeys,
    `CDS|hash|${getFeatureColorRuleHash(popupFeature)}`,
    `CDS|${label.qualifier}|${label.pattern}`
  ];
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
    }, clickedKeys]
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
assert.deepEqual(featureColorOverrides['label-feature-a'], {
  color: '#8cf04f',
  caption: 'wsv360-like protein'
});
assert.deepEqual(featureColorOverrides['label-feature-b'], {
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
  val: 'f11111111',
  color: '#123456',
  cap: 'single feature'
}]);
assert.equal(getFeatureColorRuleHash(labelFeatureA), 'f11111111');

legendEntries.value = [{ caption: 'single feature', color: '#123456', featureIds: ['f11111111_record_1'] }];
const noOpFillCount = previewFillApplyCount;
const noOpCommitCount = previewCommitCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#123456', 'single feature'), false);
assert.equal(previewFillApplyCount, noOpFillCount);
assert.equal(previewCommitCount, noOpCommitCount);

const compoundCommitCount = previewCommitCount;
const compoundFillCount = previewFillApplyCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#654321', 'renamed feature'), true);
assert.ok(previewFillApplyCount > compoundFillCount);
assert.equal(previewCommitCount, compoundCommitCount);
assert.equal(manualSpecificRules[0].cap, 'renamed feature');
assert.equal(legendEntries.value[0].caption, 'renamed feature');
legendEntries.value = [];

const outsideLabelGroup = {
  ...labelFeatureB,
  id: 'label-feature-outside',
  svg_id: 'f33333333_record_2',
  rendered_feature_svg_id: 'f33333333_record_2',
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
    { feat: 'CDS', qual: 'hash', val: 'f11111111', color: '#654321' },
    { feat: 'CDS', qual: 'hash', val: 'f22222222', color: '#654321' }
  ]
);
assert.equal(featureColorOverrides['label-feature-outside'], undefined);

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
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'product'), false);
assert.equal(manualSpecificRules.filter((rule) => rule.qual === 'hash').length, 2);

manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'CDS',
  qual: 'record_location',
  val: '^RecA:0\\.\\.90:\\+$',
  color: '#202020',
  cap: 'existing record rule'
});
Object.assign(featureStyleScopeDialog, { show: true, feat: conflictingFeatureA, color: '#aabbcc' });
await actions.handleColorScopeChoice('displayLabel');
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'product'), false);
assert.equal(manualSpecificRules.filter((rule) => rule.qual === 'hash').length, 2);

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
assert.equal(manualSpecificRules.every((rule) => rule.qual === 'hash'), true);
assert.equal(manualSpecificRules.length, 2);

const duplicateFeatureA = {
  ...labelFeatureA,
  id: 'duplicate-a',
  svg_id: 'f44444444_record_1',
  rendered_feature_svg_id: 'f44444444_record_1',
  product: '',
  qualifiers: {},
  selector: { hash: 'f44444444', qualifiers: {} }
};
const duplicateFeatureB = {
  ...duplicateFeatureA,
  id: 'duplicate-b',
  svg_id: 'f44444444_record_2',
  rendered_feature_svg_id: 'f44444444_record_2'
};
manualSpecificRules.splice(0, manualSpecificRules.length, {
  feat: 'CDS',
  qual: 'hash',
  val: 'f44444444',
  color: '#999999',
  cap: 'shared duplicate rule'
});
extractedFeatures.value = [duplicateFeatureA, duplicateFeatureB];
biologicalFeatures.value = [duplicateFeatureA, duplicateFeatureB];
await actions.setFeatureColor(duplicateFeatureA, '#112233', 'one duplicate');
// Duplicate records share the stable hash Python matches, so the edit replaces
// the shared rule and colors both copies (D-14 accepted residual risk).
assert.deepEqual(
  manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
  [{ feat: 'CDS', qual: 'hash', val: 'f44444444', color: '#112233', cap: 'one duplicate' }]
);
await actions.setFeatureColorValue(duplicateFeatureB, null);
assert.deepEqual(manualSpecificRules, []);

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
assert.equal(manualSpecificRules[0].val, 'f11111111');
assert.equal(manualSpecificRules[1].val, '^f.*$');

const stableColorFeature = {
  ...labelFeatureA,
  id: 'legacy-color-id',
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
const resetStrokeCommitCount = previewCommitCount;
assert.equal(await actions.updateClickedFeatureStroke('#111111', 1), false);
assert.equal(await actions.updateClickedFeatureStroke('#222222', 2), true);
assert.deepEqual(featureStrokeOverrides[strokeKey], {
  originalStrokeColor: '#111111', originalStrokeWidth: 1, strokeColor: '#222222', strokeWidth: 2
});
assert.equal(strokeMutationCount, 0);
assert.equal(previewCommitCount, resetStrokeCommitCount);
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
assert.equal(clickedFeature.value, null);
assert.equal(strokeAttributes.get('stroke'), '#111111');
assert.equal(strokeAttributes.get('stroke-width'), '1');
assert.equal(previewCommitCount, resetStrokeCommitCount);
assert.equal(await actions.handleFeatureStyleScopeChoice('cancel'), false);
assert.equal(featureStyleScopeDialog.show, false);

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
assert.equal(previewCommitCount, resetStrokeCommitCount);

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

assert.equal(legendGeometryChangedCount, 1);
assert.equal(legendText.textContent, 'Oxidative phosphorylation');
assert.equal(legendAttributes.get('data-legend-key'), 'Oxidative phosphorylation');

// FE-10: Reset fill uses the palette default of the feature being reset. A
// canceled or completed Reset dialog of another feature type must not leave a
// color behind for a later Reset that has no dialog.
{
  const feature = (svgId, type, product, start, end) => ({
    id: svgId, svg_id: svgId, type, product, qualifiers: { product: [product] }, start, end
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
    evaluate: evaluatePythonRules
  });
  const resetActions = createFeatureColorActions({
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
      commitSpecificRules: async (rules) => {
        committed.push(rules.map((rule) => ({ ...rule })));
        return true;
      },
      effectiveLegendCaptions: () => (feature) => feature.type,
      getLegendRowRules: (caption) => legendRowRules(caption, { rules: resetRules }),
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      findFeaturesWithSameLegendItem: () => [],
      findFeaturesWithSameDisplayedLabel: () => [],
      findFeaturesWithSameIndividualLabel: () => [],
      getDisplayedFeatureLabel: (feature) => feature.product,
      getIndividualFeatureLabel: (feature) => feature.product,
      getLabelSpecificRule: () => null
    },
    getFeatureElements: () => [],
    getFeatureFillElements: () => [],
    commitActiveResultEdit: null
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
    const rrnaReset = committed.at(-1).find((rule) => rule.qual === 'hash' && rule.val === rrna.svg_id);
    assert.equal(rrnaReset?.color, '#71ee7d', `Reset after a ${choice} dialog uses the rRNA default`);
  }
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
  const build = ({ entries, order, rules = [], features = [], legendColorOverrides = {} }) => {
    const groupEntries = new Map(entries.map((entry) => [entry.caption, fakeEntry(entry.caption)]));
    const legendGroup = {
      querySelector: (selector) => [...groupEntries].find(([caption]) => selector.includes(`"${caption}"`))?.[1] || null
    };
    const svgRoot = legendSvg(legendGroup);
    const committed = [];
    const legendRenameDialog = {};
    const stateLegendEntries = ref(entries.map((entry) => ({ featureIds: [], originalCaption: entry.caption, ...entry })));
    const originalOrder = ref([...order]);
    const featureList = ref(features);
    const renameState = { extractedFeatures: featureList, biologicalFeatures: featureList, manualSpecificRules: rules };
    const renameActions = createFeatureColorActions({
      state: withDrawings({
        results: ref([]), selectedResultIndex: ref(0), appliedPaletteColors: ref({ tRNA: '#e8b441' }),
        manualSpecificRules: rules, extractedFeatures: featureList, biologicalFeatures: featureList,
        featureColorOverrides: {}, svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svgRoot : null }),
        clickedFeature: ref(null), featureStyleScopeDialog: {}, resetColorDialog: {}, legendRenameDialog,
        legendEntries: stateLegendEntries, legendStrokeOverrides: {}, legendColorOverrides,
        originalLegendOrder: originalOrder, originalLegendColors: ref({}),
        featureStrokeOverrides: {}, skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false),
        addedLegendCaptions: ref(new Set())
      }),
      nextTick: async () => {},
      onLegendGeometryChanged: () => {}, extractLegendEntries: () => {},
      ruleActions: {
        runWithRuleMatches: runWithRuleMatchesOf(createRulePreparation({ state: withDrawings(renameState), evaluate: evaluatePythonRules }), renameState),
        commitSpecificRules: async (nextRules) => { committed.push(nextRules.map((rule) => ({ ...rule }))); return true; },
        effectiveLegendCaptions: () => (feature) => feature.legendCaption || rules.find((rule) => rule.feat === feature.type)?.cap || feature.type,
        getLegendRowRules: (caption) => legendRowRules(caption, {
          rules, legendEntries: stateLegendEntries.value, originalLegendOrder: originalOrder.value
        }),
        getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
        findFeaturesWithSameLegendItem: () => [], findFeaturesWithSameDisplayedLabel: () => [],
        findFeaturesWithSameIndividualLabel: () => [], getDisplayedFeatureLabel: (feature) => feature.product,
        getIndividualFeatureLabel: (feature) => feature.product, getLabelSpecificRule: () => null
      },
      getFeatureElements: () => [],
      getFeatureFillElements: () => [],
      commitActiveResultEdit: null
    });
    return { renameActions, committed, legendRenameDialog, stateLegendEntries, originalOrder, legendColorOverrides };
  };

  // PV-02: a renamed renderer-generated row keeps its generated identity, so
  // Generate replays the rename.
  const gc = build({
    entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'GC content', color: '#a1a1a1' }],
    order: ['CDS', 'GC content']
  });
  await gc.renameActions.renameLegendEntry(1, 'GC percent');
  assert.equal(gc.stateLegendEntries.value[1].caption, 'GC percent');
  assert.equal(gc.stateLegendEntries.value[1].originalCaption, 'GC content');
  assert.deepEqual(gc.originalOrder.value, ['CDS', 'GC content']);
  const gcPlan = compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: gc.stateLegendEntries.value,
    originalLegendOrder: gc.originalOrder.value
  });
  assert.deepEqual(gcPlan.operationsByResult[0].legendRenames.map(({ from, to }) => [from, to]), [['GC content', 'GC percent']]);

  // OV-60: a renamed row onto a color rule's caption (a rule without a drawn
  // row) is replayed by Generate, and a style an earlier row left under the
  // caption does not follow the renamed row.
  const unusedRule = { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: '#c83366', cap: 'Zeta' };
  const onRule = build({
    entries: [{ caption: 'CDS', color: '#54bcf8' }, { caption: 'GC content', color: '#a1a1a1' }],
    order: ['CDS', 'GC content'],
    rules: [unusedRule],
    legendColorOverrides: { Zeta: '#7b2cbf' }
  });
  await onRule.renameActions.renameLegendEntry(1, 'Zeta');
  assert.notEqual(onRule.legendRenameDialog.show, true);
  assert.deepEqual(onRule.legendColorOverrides, {}, 'the stale style of the caption is dropped');
  const onRulePlan = compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: onRule.stateLegendEntries.value,
    originalLegendOrder: onRule.originalOrder.value,
    manualSpecificRules: [unusedRule]
  });
  assert.deepEqual(onRulePlan.operationsByResult[0].legendRenames.map(({ from, to, allowMissing }) => [from, to, allowMissing]),
    [['GC content', 'Zeta', true]]);
  assert.deepEqual(onRulePlan.operationsByResult[0].legendFills, []);

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
  const suffix = collide();
  await suffix.renameActions.renameLegendEntry(0, 'Leu B');
  await suffix.renameActions.handleLegendRenameChoice('suffix');
  assert.deepEqual(suffix.committed.at(-1).map(({ cap, color }) => [cap, color]), [['Leu B (1)', '#e8b441']]);
  const cancel = collide();
  await cancel.renameActions.renameLegendEntry(0, 'Leu B');
  await cancel.renameActions.handleLegendRenameChoice('cancel');
  assert.equal(cancel.committed.length, 0);
  assert.equal(cancel.legendRenameDialog.show, false);

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
    featurelessSameTrack: { source: 1, target: 'GC skew (-)', suffixed: 'GC skew (-) (1)', make: () => gcSkew('#ad72e3') },
    featurelessSameTrackSameColor: { source: 1, target: 'GC skew (-)', suffixed: 'GC skew (-) (1)', make: () => gcSkew('#6dded3') }
  };
  for (const [name, { source, target, suffixed: expectedSuffix, make }] of Object.entries(noMerge)) {
    const refused = make();
    const captions = refused.stateLegendEntries.value.map((entry) => entry.caption);
    await refused.renameActions.renameLegendEntry(source, target);
    assert.equal(refused.legendRenameDialog.show, true, name);
    assert.equal(refused.legendRenameDialog.mergeAvailable, false, name);
    await refused.renameActions.handleLegendRenameChoice('merge');
    assert.equal(refused.legendRenameDialog.show, false, name);
    assert.equal(refused.committed.length, 0, name);
    assert.deepEqual(refused.stateLegendEntries.value.map((entry) => entry.caption), captions, name);
    const suffixed = make();
    await suffixed.renameActions.renameLegendEntry(source, target);
    await suffixed.renameActions.handleLegendRenameChoice('suffix');
    const suffixedCaption = suffixed.committed.length
      ? suffixed.committed.at(-1)[0].cap : suffixed.stateLegendEntries.value[source].caption;
    assert.notEqual(suffixedCaption, target, name);
    assert.ok(suffixedCaption.startsWith(target), name);
    if (!expectedSuffix) continue;
    // A row without features renames in the Legend only: no rule commit, one
    // rename that Generate replays, and nothing that reaches the target row.
    assert.equal(suffixedCaption, expectedSuffix, name);
    assert.equal(suffixed.committed.length, 0, name);
    assert.deepEqual(suffixed.stateLegendEntries.value.map((entry) => entry.caption),
      captions.map((caption, index) => (index === source ? expectedSuffix : caption)), `${name}: one row per caption, order kept`);
    const plan = compileDirectEditorMutationPlan({
      catalogAdmission: { resultNames: ['a.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
      legendEntries: suffixed.stateLegendEntries.value,
      originalLegendOrder: suffixed.originalOrder.value
    });
    const operations = plan.operationsByResult[0];
    assert.deepEqual(operations.legendRenames.map(({ from, to }) => [from, to]), [[captions[source], expectedSuffix]], name);
    const touchingTarget = Object.entries(operations).flatMap(([kind, list]) => (Array.isArray(list) ? list : [])
      .filter((operation) => [operation?.caption, operation?.from, operation?.to].includes(target))
      .map((operation) => [kind, operation]));
    assert.deepEqual(touchingTarget, [], `${name}: no operation reaches ${target}`);
  }

  // A target owned by a specific-color rule keeps PD-OI-042 disambiguation
  // for a row of the same type.
  const ruleOwned = build({
    entries: [{ caption: 'tRNA', color: '#e8b441', featureIds: ['t1'] }, { caption: 'Special', color: '#ff0000', featureIds: ['t2'] }],
    order: ['tRNA', 'Special'],
    rules: [{ feat: 'tRNA', qual: 'product', val: '^NOMATCH$', color: '#ff0000', cap: 'Special' }],
    features: [{ ...trna, legendCaption: 'tRNA' }, { ...trna, id: 't2', svg_id: 't2', legendCaption: 'Special' }]
  });
  await ruleOwned.renameActions.renameLegendEntry(0, 'Special');
  assert.notEqual(ruleOwned.legendRenameDialog.show, true);
  assert.equal(ruleOwned.committed.length, 1);
}
