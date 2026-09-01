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
const colorActionsSource = await readFile(join(sourceDir, 'app', 'feature-editor', 'color-actions.js'), 'utf8');
await writeFile(
  join(tempDir, 'app', 'feature-editor', 'color-actions.js'),
  colorActionsSource,
  'utf8'
);
await writeFile(join(tempDir, 'app', 'feature-utils.js'), await readFile(join(sourceDir, 'app', 'feature-utils.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-selector.js'), await readFile(join(sourceDir, 'app', 'feature-selector.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'color-utils.js'), await readFile(join(sourceDir, 'app', 'color-utils.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'losat-normalization.js'), await readFile(join(sourceDir, 'app', 'losat-normalization.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'svg-serialization.js'), await readFile(join(sourceDir, 'services', 'svg-serialization.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'feature-catalog.js'), await readFile(join(sourceDir, 'services', 'feature-catalog.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'feature-identity.js'), await readFile(join(sourceDir, 'services', 'feature-identity.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'orthogroup-feature-metadata.js'), await readFile(join(sourceDir, 'services', 'orthogroup-feature-metadata.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'runtime-test-hooks.js'), await readFile(join(sourceDir, 'services', 'runtime-test-hooks.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'feature-override-identity.js'), await readFile(join(sourceDir, 'services', 'feature-override-identity.js'), 'utf8'), 'utf8');

const { createFeatureColorActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'color-actions.js'))
);
const { getFeatureGenerationHash } = await import(pathToFileURL(join(tempDir, 'app', 'feature-utils.js')));
const { resolveFeatureLabelSelector } = await import(pathToFileURL(join(tempDir, 'app', 'feature-selector.js')));

assert.doesNotMatch(colorActionsSource, /serializeCleanSvg|results\.value\[[^\]]+\]\s*=/);

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

let addLegendEntryCount = 0;
let applySpecificRulesCount = 0;
let legendGeometryChangedCount = 0;
const addLegendEntryOptions = [];
const removeLegendEntryOptions = [];
const updateLegendEntryColorOptions = [];
const previewFillColors = new Map();
const featureElementsById = new Map();
let previewFillApplyCount = 0;
let previewFlushCount = 0;
let previewDirty = false;
const svgContainer = ref(null);
const clickedFeature = ref(null);
const originalSvgStroke = ref({ color: null, width: null });
const featureStrokeOverrides = {};
const legendStrokeOverrides = {};
const previewRuntime = {
  applyFeatureFillChanges: (changes) => {
    let updated = false;
    for (const change of changes) {
      if (previewFillColors.get(change.featureId) === change.color) continue;
      previewFillColors.set(change.featureId, change.color);
      previewFillApplyCount += 1;
      updated = true;
    }
    if (updated) previewDirty = true;
    return updated;
  },
  flushActiveResult: () => {
    if (!previewDirty) return false;
    previewDirty = false;
    previewFlushCount += 1;
    return true;
  },
  getActiveRuntime: () => ({ dirty: previewDirty }),
  markActiveResultDirty: () => {
    previewDirty = true;
    return true;
  }
};

const actions = createFeatureColorActions({
  state: {
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
    legendColorOverrides: {},
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    originalSvgStroke,
    featureStrokeOverrides,
    skipCaptureBaseConfig: ref(false),
    skipExtractOnSvgChange: ref(false),
    addedLegendCaptions: ref(new Set())
  },
  nextTick: async () => {},
  legendActions: {
    addLegendEntry: async (_caption, _color, options = {}) => {
      addLegendEntryCount += 1;
      addLegendEntryOptions.push(options);
      return '';
    },
    removeLegendEntry: (_caption, options = {}) => {
      removeLegendEntryOptions.push(options);
      return false;
    },
    updateLegendEntryColorByCaption: (caption, color, options = {}) => {
      updateLegendEntryColorOptions.push(options);
      const entry = legendEntries.value.find((candidate) => candidate.caption === caption);
      if (!entry || entry.color === color) return false;
      entry.color = color;
      return true;
    },
    compactLegendEntries: () => {},
    onLegendGeometryChanged: () => {
      legendGeometryChangedCount += 1;
    },
    extractLegendEntries: () => {},
    getAllFeatureLegendGroups: (svg) => svg?.legendGroups || []
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
    getFeatureElements: (_svg, featureId) => featureElementsById.get(featureId) || [],
    getFeatureFillElements: (_svg, featureId) => featureElementsById.get(featureId) || []
  },
  previewRuntime
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
assert.equal(getFeatureGenerationHash(labelFeatureA), 'f11111111');

legendEntries.value = [{ caption: 'single feature', color: '#123456', featureIds: ['f11111111_record_1'] }];
const noOpFillCount = previewFillApplyCount;
const noOpFlushCount = previewFlushCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#123456', 'single feature'), false);
assert.equal(previewFillApplyCount, noOpFillCount);
assert.equal(previewFlushCount, noOpFlushCount);

const compoundFlushCount = previewFlushCount;
assert.equal(await actions.setFeatureColor(labelFeatureA, '#654321', 'renamed feature'), true);
assert.equal(previewFlushCount - compoundFlushCount, 1);
assert.equal(addLegendEntryOptions.at(-1)?.commit, false);
assert.equal(removeLegendEntryOptions.at(-1)?.commit, false);
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
  selector: { hash: 'fbbbbbbbb', qualifiers: { gene: ['wsv360-like protein'] } }
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

const stableColorFeature = {
  ...labelFeatureA,
  id: 'legacy-color-id',
  recordKey: 'record-key-a',
  biologicalFeatureId: 'biological-a'
};
const stableColorKey = 'record-key-a\u0000biological-a';
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
const strokeSvg = { legendGroups: [] };
svgContainer.value = { querySelector: (selector) => selector === 'svg' ? strokeSvg : null };
featureElementsById.set(strokeFeature.svg_id, [strokeElement]);
clickedFeature.value = {
  svg_id: strokeFeature.svg_id,
  feat: strokeFeature,
  color: '#cccccc',
  strokeColor: '#111111',
  strokeWidth: 1,
  originalStrokeColor: '#111111',
  originalStrokeWidth: 1
};
originalSvgStroke.value = { color: '#111111', width: 1 };

const sameStrokeFlushCount = previewFlushCount;
assert.equal(await actions.updateClickedFeatureStroke('#111111', 1), false);
assert.equal(strokeMutationCount, 0);
assert.equal(previewFlushCount, sameStrokeFlushCount);

assert.equal(await actions.updateClickedFeatureStroke('#222222', 2), true);
assert.equal(strokeMutationCount, 2);
assert.equal(previewFlushCount - sameStrokeFlushCount, 1);
const changedStrokeFlushCount = previewFlushCount;
assert.equal(await actions.updateClickedFeatureStroke('#222222', 2), false);
assert.equal(strokeMutationCount, 2);
assert.equal(previewFlushCount, changedStrokeFlushCount);

assert.equal(await actions.resetClickedFeatureStroke(), true);
const resetStrokeMutationCount = strokeMutationCount;
const resetStrokeFlushCount = previewFlushCount;
assert.equal(await actions.resetClickedFeatureStroke(), false);
assert.equal(strokeMutationCount, resetStrokeMutationCount);
assert.equal(previewFlushCount, resetStrokeFlushCount);
assert.equal(await actions.applyStrokeToSelectedFeatures([strokeFeature], '#111111', 1), false);
const siblingStrokeAttributes = [featureB, hashOnlyFeature].map((feature) => {
  const attributes = new Map([
    ['stroke', '#111111'],
    ['stroke-width', '1']
  ]);
  featureElementsById.set(feature.svg_id, [{
    getAttribute: (name) => attributes.has(name) ? attributes.get(name) : null,
    setAttribute: (name, value) => attributes.set(name, String(value)),
    removeAttribute: (name) => attributes.delete(name)
  }]);
  return attributes;
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
assert.equal(previewFlushCount, resetStrokeFlushCount);
assert.equal(await actions.handleFeatureStyleScopeChoice('cancel'), false);
assert.equal(featureStyleScopeDialog.show, false);

clickedFeature.value = {
  svg_id: strokeFeature.svg_id,
  feat: strokeFeature,
  color: '#cccccc',
  legendName: 'Core',
  strokeColor: '#111111',
  strokeWidth: 1,
  originalStrokeColor: '#111111',
  originalStrokeWidth: 1
};
assert.equal(await actions.setClickedFeatureStrokeColorValue('#445566'), false);
assert.equal(await actions.handleFeatureStyleScopeChoice('caption'), true);
assert.equal(strokeAttributes.get('stroke'), '#445566');
assert.equal(strokeAttributes.get('stroke-width'), '1');
siblingStrokeAttributes.forEach((attributes) => {
  assert.equal(attributes.get('stroke'), '#445566');
  assert.equal(attributes.get('stroke-width'), '1');
});
assert.deepEqual(legendStrokeOverrides.Core, {
  originalStrokeColor: null,
  originalStrokeWidth: null,
  strokeColor: '#445566',
  strokeWidth: 1
});
assert.equal(previewFlushCount, resetStrokeFlushCount + 1);

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
const svg = { legendGroups: [legendFeatureGroup] };
svgContainer.value = { querySelector: (selector) => selector === 'svg' ? svg : null };
legendEntries.value = [{ caption: 'Short caption', color: '#123456', featureIds: [] }];
extractedFeatures.value = [];
biologicalFeatures.value = [];

await actions.renameLegendEntry(0, 'Oxidative phosphorylation');

assert.equal(legendGeometryChangedCount, 1);
assert.equal(legendText.textContent, 'Oxidative phosphorylation');
assert.equal(legendAttributes.get('data-legend-key'), 'Oxidative phosphorylation');
assert.equal(addLegendEntryOptions.length > 0, true);
assert.equal(removeLegendEntryOptions.length > 0, true);
assert.equal(updateLegendEntryColorOptions.length > 0, true);
assert.equal(addLegendEntryOptions.every((options) => options.commit === false), true);
assert.equal(removeLegendEntryOptions.every((options) => options.commit === false), true);
assert.equal(updateLegendEntryColorOptions.every((options) => options.commit === false), true);
