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
await writeFile(join(tempDir, 'app', 'color-utils.js'), await readFile(join(sourceDir, 'app', 'color-utils.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'services', 'svg-serialization.js'), await readFile(join(sourceDir, 'services', 'svg-serialization.js'), 'utf8'), 'utf8');

const { createFeatureColorActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'color-actions.js'))
);

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
    extractedFeatures: ref([featureA, featureB, hashOnlyFeature]),
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
    findFeaturesWithSameDisplayedLabel: () => [],
    findFeaturesWithSameIndividualLabel: () => [],
    findFeaturesWithSameLegendItem: () => [featureB, hashOnlyFeature],
    findMatchingRegexRule: () => specificRule,
    getDisplayedFeatureLabel: () => '',
    getEffectiveLegendCaption: () => 'Core',
    getIndividualFeatureLabel: () => '',
    getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id })
  },
  featureSvgActions: {
    applyInstantPreview: () => {},
    getFeatureElements: () => []
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
