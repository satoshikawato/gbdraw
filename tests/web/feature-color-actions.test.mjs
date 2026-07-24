import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-color-actions-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app', 'feature-editor'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(
  join(tempDir, 'app', 'feature-editor', 'color-actions.js'),
  await readFile(join(sourceDir, 'app', 'feature-editor', 'color-actions.js'), 'utf8'),
  'utf8'
);
await writeFile(join(tempDir, 'app', 'feature-utils.js'), await readFile(join(sourceDir, 'app', 'feature-utils.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-selector.js'), await readFile(join(sourceDir, 'app', 'feature-selector.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'color-utils.js'), await readFile(join(sourceDir, 'app', 'color-utils.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'svg-serialization.js'), await readFile(join(sourceDir, 'services', 'svg-serialization.js'), 'utf8'), 'utf8');

const { createFeatureColorActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'color-actions.js'))
);
const { getFeatureGenerationHash } = await import(pathToFileURL(join(tempDir, 'app', 'feature-utils.js')));
const { resolveFeatureLabelSelector } = await import(pathToFileURL(join(tempDir, 'app', 'feature-selector.js')));

const ref = (value) => ({ value });

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
const extractedFeatures = ref([featureA, featureB, hashOnlyFeature]);
const biologicalFeatures = ref([featureA, featureB, hashOnlyFeature]);
const legendEntries = ref([{ caption: 'Core', color: '#111111', featureIds: ['hash-a', 'hash-b', 'hash-c'] }]);
const colorScopeDialog = {
  show: true,
  feat: featureA,
  color: '#abcdef',
  matchingRule: specificRule,
  ruleMatchCount: 2,
  legendName: 'Core',
  siblingCount: 2,
  displayLabel: null,
  displayLabelSiblingCount: 0,
  annotationLabel: null,
  annotationLabelSiblingCount: 0,
  individualLabel: null,
  individualLabelSiblingCount: 0,
  existingCaptionRule: null,
  existingCaptionColor: null,
  resolve: null
};

let addLegendEntryCount = 0;
let applySpecificRulesCount = 0;

const actions = createFeatureColorActions({
  state: {
    pyodideReady: ref(true),
    results: ref([]),
    selectedResultIndex: ref(0),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    manualSpecificRules,
    extractedFeatures,
    biologicalFeatures,
    featureColorOverrides,
    svgContainer: ref(null),
    clickedFeature: ref(null),
    colorScopeDialog,
    resetColorDialog: {},
    legendRenameDialog: {},
    legendEntries,
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    originalSvgStroke: ref({ color: null, width: null }),
    featureStrokeOverrides: {},
    skipCaptureBaseConfig: ref(false),
    skipExtractOnSvgChange: ref(false),
    addedLegendCaptions: ref(new Set())
  },
  nextTick: async () => {},
  legendActions: {
    addLegendEntry: async () => {
      addLegendEntryCount += 1;
      return '';
    },
    removeLegendEntry: () => {},
    updateLegendEntryColorByCaption: (caption, color) => {
      const entry = legendEntries.value.find((candidate) => candidate.caption === caption);
      if (entry) entry.color = color;
    },
    compactLegendEntries: () => {},
    recenterCurrentLegendRoot: () => {},
    extractLegendEntries: () => {},
    getAllFeatureLegendGroups: () => []
  },
  svgActions: {
    applySpecificRulesToSvg: () => {
      applySpecificRulesCount += 1;
    }
  },
  ruleActions: {
    countFeaturesMatchingRule: () => 0,
    findExistingColorForCaption: () => null,
    findFeaturesWithSameDisplayedLabel: (currentFeature, label) => extractedFeatures.value.filter(
      (feature) =>
        feature.svg_id !== currentFeature.svg_id &&
        (feature.displayLabel || feature.product) === label
    ),
    findFeaturesWithSameIndividualLabel: () => [],
    findFeaturesWithSameLegendItem: () => [featureB, hashOnlyFeature],
    findMatchingRegexRule: () => specificRule,
    getDisplayedFeatureLabel: (feature) => feature.displayLabel || feature.product || '',
    getEffectiveLegendCaption: () => 'Core',
    getIndividualFeatureLabel: (feature) => feature.product || '',
    getFeatureQualifier: (feature) => {
      const generationHash = getFeatureGenerationHash(feature);
      const collisionCount = extractedFeatures.value.filter(
        (candidate) => candidate.type === feature.type && getFeatureGenerationHash(candidate) === generationHash
      ).length;
      return {
        qual: 'hash',
        val: collisionCount > 1 ? feature.svg_id : generationHash
      };
    },
    getLabelSpecificRule: (feature, label) => {
      const selector = resolveFeatureLabelSelector(feature, label);
      return selector
        ? { feat: feature.type, qual: selector.qualifier, val: selector.pattern }
        : null;
    }
  },
  featureSvgActions: {
    applyInstantPreview: () => {},
    getFeatureElements: () => [],
    getFeatureFillElements: () => []
  }
});

await actions.handleColorScopeChoice('caption');

assert.equal(addLegendEntryCount, 0);
assert.equal(applySpecificRulesCount, 1);
assert.equal(specificRule.color, '#abcdef');
assert.equal(legendEntries.value[0].color, '#abcdef');
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-a'), false);
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'hash-b'), false);
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-c')?.color, '#abcdef');
assert.equal(manualSpecificRules.find((rule) => rule.qual === 'hash' && rule.val === 'hash-z')?.color, '#444444');
assert.deepEqual(featureColorOverrides['feature-a'], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides['feature-b'], { color: '#abcdef', caption: 'Core' });
assert.deepEqual(featureColorOverrides['feature-c'], { color: '#abcdef', caption: 'Core' });

const labelFeatureA = {
  id: 'label-feature-a',
  svg_id: 'f11111111_record_1',
  rendered_feature_svg_id: 'f11111111_record_1',
  stable_svg_id: 'faaaaaaaa',
  type: 'CDS',
  product: 'wsv360-like protein',
  qualifiers: { product: ['wsv360-like protein'] },
  selector: { hash: 'faaaaaaaa', qualifiers: { product: ['wsv360-like protein'] } }
};
const labelFeatureB = {
  id: 'label-feature-b',
  svg_id: 'f22222222_record_2',
  rendered_feature_svg_id: 'f22222222_record_2',
  stable_svg_id: 'fbbbbbbbb',
  type: 'CDS',
  product: 'wsv360-like protein',
  qualifiers: { product: ['wsv360-like protein'] },
  selector: { hash: 'fbbbbbbbb', qualifiers: { product: ['wsv360-like protein'] } }
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
Object.assign(colorScopeDialog, {
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
assert.equal(getFeatureGenerationHash(labelFeatureA), 'f11111111');

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
Object.assign(colorScopeDialog, {
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
    hash: 'faaaaaaaa',
    record_location: 'RecA:0..90:+',
    qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] }
  }
};
const conflictingFeatureB = {
  ...labelFeatureB,
  qualifiers: { product: ['wsv360-like protein'], gene: ['wsv360'] },
  selector: {
    hash: 'fbbbbbbbb',
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
Object.assign(colorScopeDialog, {
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
Object.assign(colorScopeDialog, { show: true, feat: conflictingFeatureA, color: '#aabbcc' });
await actions.handleColorScopeChoice('displayLabel');
assert.equal(manualSpecificRules.some((rule) => rule.qual === 'product'), false);
assert.equal(manualSpecificRules.filter((rule) => rule.qual === 'hash').length, 2);

const geneLabelFeature = {
  ...labelFeatureB,
  product: '',
  gene: 'wsv360-like protein',
  displayLabel: 'wsv360-like protein',
  qualifiers: { gene: ['wsv360-like protein'] },
  selector: { hash: 'fbbbbbbbb', qualifiers: { gene: ['wsv360-like protein'] } }
};
manualSpecificRules.splice(0);
extractedFeatures.value = [labelFeatureA, geneLabelFeature];
biologicalFeatures.value = [labelFeatureA, geneLabelFeature];
Object.assign(colorScopeDialog, {
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
assert.equal(manualSpecificRules[0].val, 'f44444444_record_1');
assert.equal(manualSpecificRules[1].val, 'f44444444');

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
