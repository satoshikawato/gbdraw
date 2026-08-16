import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

// Run this test against the detached PR base with:
// GBDRAW_CHARACTERIZATION_ROOT=/path/to/base node tests/web/live-svg-edit-base-characterization.test.mjs
const sourceRoot = process.env.GBDRAW_CHARACTERIZATION_ROOT || process.cwd();
const sourceDir = join(sourceRoot, 'gbdraw', 'web', 'js', 'app');
const visibilityActionsSource = await readFile(
  join(sourceDir, 'feature-editor', 'visibility-actions.js'),
  'utf8'
);
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-live-svg-base-oracle-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app', 'feature-editor'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });

for (const relativePath of [
  ['app', 'feature-editor', 'visibility-actions.js'],
  ['app', 'feature-visibility.js'],
  ['app', 'feature-selector.js'],
  ['app', 'feature-utils.js'],
  ['services', 'text-download.js'],
  ['services', 'feature-identity.js']
]) {
  const [owner, ...parts] = relativePath;
  const source = owner === 'app'
    ? join(sourceDir, ...parts)
    : join(sourceDir, '..', owner, ...parts);
  await writeFile(join(tempDir, owner, ...parts), await readFile(source, 'utf8'), 'utf8');
}

const { createFeatureVisibilityActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'visibility-actions.js'))
);

const ref = (value) => ({ value });
const featureA = {
  svg_id: 'feature-a',
  type: 'CDS',
  product: 'shared product',
  qualifiers: { protein_id: ['WP_A.1'] }
};
const featureB = {
  svg_id: 'feature-b',
  type: 'CDS',
  product: 'shared product',
  qualifiers: { protein_id: ['WP_B.1'] }
};
const manualRules = [];
const overrides = {};
const preview = [];
const recordChanges = (changes) => {
  (Array.isArray(changes) ? changes : []).forEach((change) => {
    preview.push([String(change.featureId), String(change.mode)]);
  });
  return true;
};

const actions = createFeatureVisibilityActions({
  state: {
    clickedFeature: ref({ svg_id: 'feature-a', featureVisibility: 'default', feat: featureA }),
    extractedFeatures: ref([featureA, featureB]),
    orthogroups: ref([]),
    featureVisibilityManualRules: manualRules,
    featureVisibilityRules: ref([]),
    featureVisibilityOverrides: overrides,
    featureVisibilitySelectorCache: {},
    featureVisibilityScopeDialog: {},
    labelLayoutDirtyReason: ref(''),
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => ({}) })
  },
  featureSvgActions: {
    applyVisibilityPreviewBySvgId: (featureId, mode) => recordChanges([{ featureId, mode }]),
    applyVisibilityPreviewChanges: recordChanges
  },
  previewRuntime: {
    flushActiveResult: () => true,
    runDomEditSync: ({ action }) => action()
  }
});

const productScope = {
  id: 'product',
  featureType: 'CDS',
  qualifier: 'product',
  value: 'shared product',
  label: 'Exact product: shared product'
};
actions.setFeatureVisibility(featureA, 'off', { triggerReflow: false, scope: productScope });
assert.deepEqual(preview.splice(0), [['feature-a', 'off'], ['feature-b', 'off']]);
assert.equal(manualRules.at(-1).qualifier, 'product');

overrides['feature-b'] = 'on';
actions.setFeatureVisibility(featureA, 'default', { triggerReflow: false, scope: productScope });
const noDecisionMode = visibilityActionsSource.includes('runDomEditSync') ? 'default' : 'on';
assert.deepEqual(preview.splice(0), [['feature-a', noDecisionMode], ['feature-b', 'on']]);
assert.equal(manualRules.length, 0);

const proteinScope = {
  id: 'protein_id',
  featureType: 'CDS',
  qualifier: 'protein_id',
  value: 'WP_A.1',
  label: 'Exact protein ID: WP_A.1'
};
actions.setFeatureVisibility(featureA, 'off', { triggerReflow: false, scope: proteinScope });
assert.deepEqual(preview.splice(0), [['feature-a', 'off']]);
assert.equal(manualRules.at(-1).qualifier, 'protein_id');

console.log(`live SVG base characterization passed for ${sourceRoot}`);
